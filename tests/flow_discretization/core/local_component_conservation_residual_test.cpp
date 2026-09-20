#include <mpmc/flow_discretization/local_component_conservation_residual.hpp>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <iostream>
#include <limits>
#include <source_location>
#include <span>
#include <stdexcept>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

bool local_component_conservation_residual_header();

namespace {

namespace fd = mpmc::flow_discretization;
namespace flow = mpmc::flow;
namespace mesh = mpmc::mesh;

void require(
    bool condition,
    std::string_view message,
    std::source_location where =
        std::source_location::current()) {
    if (!condition) {
        throw std::runtime_error(
            std::string{where.file_name()} +
            ":" +
            std::to_string(where.line()) +
            ": " +
            std::string{message});
    }
}

void near(
    double actual,
    double expected,
    double relative = 3.0e-13,
    double absolute = 3.0e-13,
    std::source_location where =
        std::source_location::current()) {
    if (!std::isfinite(actual) ||
        !std::isfinite(expected) ||
        std::abs(actual - expected) >
            absolute +
                relative *
                    std::max(
                        std::abs(actual),
                        std::abs(expected))) {
        std::cerr
            << "actual=" << actual
            << " expected=" << expected
            << '\n';
        require(false, "numeric mismatch", where);
    }
}

template <class Function>
void expect_invalid(
    Function&& function,
    std::string_view fragment) {
    try {
        function();
    } catch (const std::invalid_argument& error) {
        require(
            std::string_view{error.what()}.find(fragment) !=
                std::string_view::npos,
            "invalid_argument diagnostic changed");
        return;
    }
    throw std::runtime_error(
        "expected std::invalid_argument");
}

flow::NaturalVariableStateIdentity3P identity(
    std::size_t kind) {
    if (kind > 2U) {
        throw std::invalid_argument(
            "identity fixture kind out of range");
    }

    const std::array<
        std::array<std::size_t, 3>,
        3>
        pivots{{
            {1U, 0U, 2U},
            {0U, 2U, 1U},
            {2U, 1U, 0U}}};

    const std::array<double, 3>
        pressures{1.00e6, 0.99e6, 1.02e6};
    const std::array<double, 3>
        temperatures{350.0, 351.0, 349.0};

    const std::array<
        std::array<double, 3>,
        3>
        saturations{{
            {0.20, 0.30, 0.50},
            {0.25, 0.35, 0.40},
            {0.30, 0.20, 0.50}}};

    const std::array<
        std::array<std::vector<double>, 3>,
        3>
        compositions{{
            {
                std::vector<double>{0.10, 0.70, 0.20},
                std::vector<double>{0.60, 0.20, 0.20},
                std::vector<double>{0.20, 0.30, 0.50}},
            {
                std::vector<double>{0.55, 0.25, 0.20},
                std::vector<double>{0.25, 0.25, 0.50},
                std::vector<double>{0.20, 0.60, 0.20}},
            {
                std::vector<double>{0.20, 0.20, 0.60},
                std::vector<double>{0.30, 0.50, 0.20},
                std::vector<double>{0.65, 0.15, 0.20}}
        }};

    const auto pivot =
        flow::NaturalVariableCompositionPivot3P::
            from_dependent_components(
                3U,
                pivots[kind]);

    return {
        flow::NaturalVariableLayout3P{pivot},
        {"A", "B", "C"},
        pressures[kind],
        temperatures[kind],
        saturations[kind],
        compositions[kind]};
}

std::vector<double> component_jacobian(
    std::size_t component_count,
    std::size_t column_count,
    double scale) {
    std::vector<double> result(
        component_count * column_count,
        0.0);
    for (std::size_t component = 0U;
         component < component_count;
         ++component) {
        for (std::size_t column = 0U;
             column < column_count;
             ++column) {
            result[
                component * column_count +
                column] =
                scale *
                static_cast<double>(
                    1U +
                    component * 10U +
                    column);
        }
    }
    return result;
}

std::vector<double> total_gradient(
    const std::vector<double>& component_values,
    std::size_t component_count,
    std::size_t column_count) {
    std::vector<double> result(
        column_count,
        0.0);
    for (std::size_t column = 0U;
         column < column_count;
         ++column) {
        for (std::size_t component = 0U;
             component < component_count;
             ++component) {
            result[column] +=
                component_values[
                    component * column_count +
                    column];
        }
    }
    return result;
}

fd::ConservativeComponentFaceRateScatterLinearization3D
scatter(
    std::size_t face_id,
    flow::NaturalVariableStateIdentity3P owner_identity,
    flow::NaturalVariableStateIdentity3P neighbour_identity,
    std::vector<double> owner_rate,
    double owner_column_scale,
    double neighbour_column_scale) {
    const std::size_t n = owner_rate.size();
    const std::size_t owner_q =
        owner_identity.layout.unknown_count();
    const std::size_t neighbour_q =
        neighbour_identity.layout.unknown_count();

    auto neighbour_rate = owner_rate;
    for (double& value : neighbour_rate) {
        value = -value;
    }

    const auto owner_owner =
        component_jacobian(
            n,
            owner_q,
            owner_column_scale);
    const auto owner_neighbour =
        component_jacobian(
            n,
            neighbour_q,
            neighbour_column_scale);

    std::vector<double> neighbour_owner(
        owner_owner.size());
    std::vector<double> neighbour_neighbour(
        owner_neighbour.size());
    for (std::size_t index = 0U;
         index < owner_owner.size();
         ++index) {
        neighbour_owner[index] =
            -owner_owner[index];
    }
    for (std::size_t index = 0U;
         index < owner_neighbour.size();
         ++index) {
        neighbour_neighbour[index] =
            -owner_neighbour[index];
    }

    double owner_total = 0.0;
    for (double value : owner_rate) {
        owner_total += value;
    }

    const auto owner_total_owner =
        total_gradient(
            owner_owner,
            n,
            owner_q);
    const auto owner_total_neighbour =
        total_gradient(
            owner_neighbour,
            n,
            neighbour_q);

    std::vector<double> neighbour_total_owner(
        owner_q);
    std::vector<double> neighbour_total_neighbour(
        neighbour_q);
    for (std::size_t column = 0U;
         column < owner_q;
         ++column) {
        neighbour_total_owner[column] =
            -owner_total_owner[column];
    }
    for (std::size_t column = 0U;
         column < neighbour_q;
         ++column) {
        neighbour_total_neighbour[column] =
            -owner_total_neighbour[column];
    }

    return {
        mesh::LocalIndex{
            static_cast<mesh::LocalIndex::value_type>(
                face_id)},
        {"A", "B", "C"},
        std::move(owner_identity),
        std::move(neighbour_identity),
        std::move(owner_rate),
        std::move(neighbour_rate),
        owner_owner,
        owner_neighbour,
        std::move(neighbour_owner),
        std::move(neighbour_neighbour),
        owner_total,
        -owner_total,
        owner_total_owner,
        owner_total_neighbour,
        std::move(neighbour_total_owner),
        std::move(neighbour_total_neighbour)};
}

fd::NormalizedComponentFaceContributionLinearization3D
owner_incident_face() {
    return fd::normalize_component_face_rate_by_bulk_volume(
        scatter(
            13U,
            identity(0U),
            identity(1U),
            {2.0, -1.0, 3.0},
            0.01,
            -0.02),
        {2.0, 5.0});
}

fd::NormalizedComponentFaceContributionLinearization3D
neighbour_incident_face() {
    return fd::normalize_component_face_rate_by_bulk_volume(
        scatter(
            17U,
            identity(2U),
            identity(0U),
            {-4.0, 2.0, 1.0},
            0.03,
            -0.04),
        {7.0, 2.0});
}

flow::BackwardEulerComponentAccumulationResidual3P
accumulation() {
    auto local = identity(0U);
    const std::size_t q =
        local.layout.unknown_count();

    const std::vector<double>
        component_residual{0.4, -0.1, 0.2};
    const auto jacobian =
        component_jacobian(
            3U,
            q,
            0.005);
    const auto gradient =
        total_gradient(
            jacobian,
            3U,
            q);

    double total = 0.0;
    for (double value : component_residual) {
        total += value;
    }

    return {
        local.layout,
        0.22,
        100.0,
        {"A", "B", "C"},
        component_residual,
        total,
        q,
        jacobian,
        gradient};
}

fd::LocalComponentConservationResidualLinearization3D
assemble(
    const flow::BackwardEulerComponentAccumulationResidual3P&
        acc,
    const fd::NormalizedComponentFaceContributionLinearization3D&
        first,
    const fd::NormalizedComponentFaceContributionLinearization3D&
        second) {
    const std::array<
        fd::LocalCellNormalizedFaceContributionBinding3D,
        2>
        bindings{{
            {
                fd::IncidentFaceLocalSide3D::owner,
                &first},
            {
                fd::IncidentFaceLocalSide3D::neighbour,
                &second}
        }};

    return fd::build_local_component_conservation_residual(
        identity(0U),
        2.0,
        acc,
        std::span<
            const fd::
                LocalCellNormalizedFaceContributionBinding3D>{
                    bindings});
}

void multi_face_residual_and_blocks() {
    const auto acc = accumulation();
    const auto first = owner_incident_face();
    const auto second = neighbour_incident_face();
    const auto result =
        assemble(acc, first, second);

    require(
        result.component_ids ==
                std::vector<std::string>{"A", "B", "C"} &&
            result.neighbour_blocks.size() == 2U,
        "local conservation lost component identity or neighbour blocks");

    near(result.cell_bulk_volume_m3, 2.0);
    near(result.porosity, acc.porosity);
    near(result.time_step_seconds, acc.time_step_seconds);

    const std::size_t local_q =
        result.cell_state_identity.layout
            .unknown_count();
    for (std::size_t component = 0U;
         component < result.component_count();
         ++component) {
        near(
            result.residual(component),
            acc.residual(component) +
                first.owner_contribution(component) +
                second.neighbour_contribution(component));

        for (std::size_t column = 0U;
             column < local_q;
             ++column) {
            near(
                result.d_local(
                    component,
                    column),
                acc.d_residual(
                    component,
                    column) +
                    first.d_owner_contribution_wrt_owner(
                        component,
                        column) +
                    second.d_neighbour_contribution_wrt_neighbour(
                        component,
                        column));
        }
    }

    near(
        result.total_residual_mol_per_bulk_m3_s,
        acc.total_residual_mol_per_bulk_m3_s +
            first.owner_total_contribution_mol_per_bulk_m3_s +
            second.neighbour_total_contribution_mol_per_bulk_m3_s);

    const auto& first_block =
        result.neighbour_blocks[0];
    require(
        first_block.face ==
                mesh::LocalIndex{
                    mesh::LocalIndex::value_type{13}} &&
            first_block.local_side ==
                fd::IncidentFaceLocalSide3D::owner &&
            first_block.neighbour_state_identity.component_ids ==
                std::vector<std::string>{"A", "B", "C"},
        "owner-oriented incident face neighbour block identity changed");

    for (std::size_t component = 0U;
         component < result.component_count();
         ++component) {
        for (std::size_t column = 0U;
             column <
             first.neighbour_state_identity.layout
                 .unknown_count();
             ++column) {
            near(
                first_block.d_component(
                    component,
                    column),
                first
                    .d_owner_contribution_wrt_neighbour(
                        component,
                        column));
        }
    }

    const auto& second_block =
        result.neighbour_blocks[1];
    require(
        second_block.face ==
                mesh::LocalIndex{
                    mesh::LocalIndex::value_type{17}} &&
            second_block.local_side ==
                fd::IncidentFaceLocalSide3D::neighbour,
        "neighbour-oriented incident face block identity changed");

    for (std::size_t component = 0U;
         component < result.component_count();
         ++component) {
        for (std::size_t column = 0U;
             column <
             second.owner_state_identity.layout
                 .unknown_count();
             ++column) {
            near(
                second_block.d_component(
                    component,
                    column),
                second
                    .d_neighbour_contribution_wrt_owner(
                        component,
                        column));
        }
    }
}

flow::BackwardEulerComponentAccumulationResidual3P
perturb_accumulation(
    flow::BackwardEulerComponentAccumulationResidual3P value,
    std::size_t column,
    double delta) {
    const std::size_t n = value.component_count();
    const std::size_t q = value.input_count;
    value.total_residual_mol_per_bulk_m3_s +=
        value.total_residual_gradient[column] *
        delta;
    for (std::size_t component = 0U;
         component < n;
         ++component) {
        value
            .component_residual_mol_per_bulk_m3_s[
                component] +=
            value.component_jacobian[
                component * q +
                column] *
            delta;
    }
    return value;
}

fd::NormalizedComponentFaceContributionLinearization3D
perturb_face(
    fd::NormalizedComponentFaceContributionLinearization3D
        value,
    bool owner_column,
    std::size_t column,
    double delta) {
    const std::size_t n =
        value.component_count();
    const std::size_t q =
        owner_column
            ? value.owner_state_identity.layout
                  .unknown_count()
            : value.neighbour_state_identity.layout
                  .unknown_count();

    const auto& owner_block =
        owner_column
            ? value.owner_row_owner_column_jacobian
            : value.owner_row_neighbour_column_jacobian;
    const auto& neighbour_block =
        owner_column
            ? value.neighbour_row_owner_column_jacobian
            : value.neighbour_row_neighbour_column_jacobian;
    const auto& owner_total_gradient =
        owner_column
            ? value.owner_total_row_owner_column_gradient
            : value.owner_total_row_neighbour_column_gradient;
    const auto& neighbour_total_gradient =
        owner_column
            ? value.neighbour_total_row_owner_column_gradient
            : value.neighbour_total_row_neighbour_column_gradient;

    value.owner_total_contribution_mol_per_bulk_m3_s +=
        owner_total_gradient[column] *
        delta;
    value.neighbour_total_contribution_mol_per_bulk_m3_s +=
        neighbour_total_gradient[column] *
        delta;

    for (std::size_t component = 0U;
         component < n;
         ++component) {
        const std::size_t index =
            component * q +
            column;
        value
            .owner_component_contribution_mol_per_bulk_m3_s[
                component] +=
            owner_block[index] *
            delta;
        value
            .neighbour_component_contribution_mol_per_bulk_m3_s[
                component] +=
            neighbour_block[index] *
            delta;
    }

    return value;
}

void jacobian_fresh_source_perturbation() {
    constexpr double step = 1.0e-6;

    const auto acc = accumulation();
    const auto first = owner_incident_face();
    const auto second = neighbour_incident_face();
    const auto base =
        assemble(acc, first, second);

    const std::size_t local_q =
        base.cell_state_identity.layout
            .unknown_count();

    for (std::size_t column = 0U;
         column < local_q;
         ++column) {
        const auto plus =
            assemble(
                perturb_accumulation(
                    acc,
                    column,
                    step),
                perturb_face(
                    first,
                    true,
                    column,
                    step),
                perturb_face(
                    second,
                    false,
                    column,
                    step));
        const auto minus =
            assemble(
                perturb_accumulation(
                    acc,
                    column,
                    -step),
                perturb_face(
                    first,
                    true,
                    column,
                    -step),
                perturb_face(
                    second,
                    false,
                    column,
                    -step));

        for (std::size_t component = 0U;
             component < base.component_count();
             ++component) {
            near(
                base.d_local(
                    component,
                    column),
                (plus.residual(component) -
                 minus.residual(component)) /
                    (2.0 * step),
                2.0e-9,
                2.0e-9);
        }
        near(
            base.d_total_local(column),
            (plus.total_residual_mol_per_bulk_m3_s -
             minus.total_residual_mol_per_bulk_m3_s) /
                (2.0 * step),
            2.0e-9,
            2.0e-9);
    }

    for (std::size_t block_index = 0U;
         block_index < base.neighbour_blocks.size();
         ++block_index) {
        const bool first_face =
            block_index == 0U;
        const auto& block =
            base.neighbour_blocks[block_index];
        const std::size_t neighbour_q =
            block.neighbour_state_identity.layout
                .unknown_count();

        for (std::size_t column = 0U;
             column < neighbour_q;
             ++column) {
            const auto plus =
                first_face
                    ? assemble(
                          acc,
                          perturb_face(
                              first,
                              false,
                              column,
                              step),
                          second)
                    : assemble(
                          acc,
                          first,
                          perturb_face(
                              second,
                              true,
                              column,
                              step));
            const auto minus =
                first_face
                    ? assemble(
                          acc,
                          perturb_face(
                              first,
                              false,
                              column,
                              -step),
                          second)
                    : assemble(
                          acc,
                          first,
                          perturb_face(
                              second,
                              true,
                              column,
                              -step));

            for (std::size_t component = 0U;
                 component < base.component_count();
                 ++component) {
                near(
                    block.d_component(
                        component,
                        column),
                    (plus.residual(component) -
                     minus.residual(component)) /
                        (2.0 * step),
                    2.0e-9,
                    2.0e-9);
            }
            near(
                block.d_total(column),
                (plus.total_residual_mol_per_bulk_m3_s -
                 minus.total_residual_mol_per_bulk_m3_s) /
                    (2.0 * step),
                2.0e-9,
                2.0e-9);
        }
    }
}

void zero_primal_nonzero_jacobian() {
    auto acc = accumulation();
    std::fill(
        acc.component_residual_mol_per_bulk_m3_s.begin(),
        acc.component_residual_mol_per_bulk_m3_s.end(),
        0.0);
    acc.total_residual_mol_per_bulk_m3_s = 0.0;

    auto first = owner_incident_face();
    auto second = neighbour_incident_face();

    for (auto* face : {&first, &second}) {
        std::fill(
            face->owner_component_contribution_mol_per_bulk_m3_s.begin(),
            face->owner_component_contribution_mol_per_bulk_m3_s.end(),
            0.0);
        std::fill(
            face->neighbour_component_contribution_mol_per_bulk_m3_s.begin(),
            face->neighbour_component_contribution_mol_per_bulk_m3_s.end(),
            0.0);
        face->owner_total_contribution_mol_per_bulk_m3_s =
            0.0;
        face->neighbour_total_contribution_mol_per_bulk_m3_s =
            0.0;
    }

    const auto result =
        assemble(acc, first, second);

    require(
        std::all_of(
            result
                .component_residual_mol_per_bulk_m3_s
                .begin(),
            result
                .component_residual_mol_per_bulk_m3_s
                .end(),
            [](double value) {
                return value == 0.0;
            }) &&
            result.total_residual_mol_per_bulk_m3_s ==
                0.0,
        "zero primal accumulation/spatial contributions did not produce zero local residual");

    require(
        std::any_of(
            result.local_component_jacobian.begin(),
            result.local_component_jacobian.end(),
            [](double value) {
                return value != 0.0;
            }) &&
            std::any_of(
                result.neighbour_blocks[0]
                    .component_jacobian.begin(),
                result.neighbour_blocks[0]
                    .component_jacobian.end(),
                [](double value) {
                    return value != 0.0;
                }),
        "zero local residual incorrectly erased structural Jacobian");
}

void invalid_inputs() {
    const auto acc = accumulation();
    const auto first = owner_incident_face();
    const auto second = neighbour_incident_face();

    {
        const std::array<
            fd::LocalCellNormalizedFaceContributionBinding3D,
            2>
            duplicate{{
                {
                    fd::IncidentFaceLocalSide3D::owner,
                    &first},
                {
                    fd::IncidentFaceLocalSide3D::owner,
                    &first}
            }};
        expect_invalid(
            [&] {
                (void)fd::
                    build_local_component_conservation_residual(
                        identity(0U),
                        2.0,
                        acc,
                        duplicate);
            },
            "duplicate incident face");
    }

    {
        const std::array<
            fd::LocalCellNormalizedFaceContributionBinding3D,
            1>
            null_binding{{
                {
                    fd::IncidentFaceLocalSide3D::owner,
                    nullptr}
            }};
        expect_invalid(
            [&] {
                (void)fd::
                    build_local_component_conservation_residual(
                        identity(0U),
                        2.0,
                        acc,
                        null_binding);
            },
            "must not be null");
    }

    {
        const std::array<
            fd::LocalCellNormalizedFaceContributionBinding3D,
            1>
            wrong_side{{
                {
                    fd::IncidentFaceLocalSide3D::neighbour,
                    &first}
            }};
        expect_invalid(
            [&] {
                (void)fd::
                    build_local_component_conservation_residual(
                        identity(0U),
                        2.0,
                        acc,
                        wrong_side);
            },
            "selected side");
    }

    {
        const std::array<
            fd::LocalCellNormalizedFaceContributionBinding3D,
            1>
            binding{{
                {
                    fd::IncidentFaceLocalSide3D::owner,
                    &first}
            }};
        expect_invalid(
            [&] {
                (void)fd::
                    build_local_component_conservation_residual(
                        identity(0U),
                        2.1,
                        acc,
                        binding);
            },
            "bulk volume");
    }

    {
        auto wrong_acc = acc;
        wrong_acc.component_ids[1] = "X";
        expect_invalid(
            [&] {
                const std::array<
                    fd::LocalCellNormalizedFaceContributionBinding3D,
                    0>
                    none{};
                (void)fd::
                    build_local_component_conservation_residual(
                        identity(0U),
                        2.0,
                        wrong_acc,
                        none);
            },
            "does not match");
    }

    {
        auto wrong_acc = acc;
        wrong_acc.total_residual_mol_per_bulk_m3_s +=
            1.0;
        expect_invalid(
            [&] {
                const std::array<
                    fd::LocalCellNormalizedFaceContributionBinding3D,
                    0>
                    none{};
                (void)fd::
                    build_local_component_conservation_residual(
                        identity(0U),
                        2.0,
                        wrong_acc,
                        none);
            },
            "do not close");
    }

    {
        auto wrong_face = first;
        wrong_face
            .owner_component_contribution_mol_per_bulk_m3_s[0] +=
            1.0;
        const std::array<
            fd::LocalCellNormalizedFaceContributionBinding3D,
            1>
            binding{{
                {
                    fd::IncidentFaceLocalSide3D::owner,
                    &wrong_face}
            }};
        expect_invalid(
            [&] {
                (void)fd::
                    build_local_component_conservation_residual(
                        identity(0U),
                        2.0,
                        acc,
                        binding);
            },
            "weighted conservation");
    }

    {
        const std::array<
            fd::LocalCellNormalizedFaceContributionBinding3D,
            1>
            invalid_enum{{
                {
                    static_cast<
                        fd::IncidentFaceLocalSide3D>(
                            99),
                    &first}
            }};
        expect_invalid(
            [&] {
                (void)fd::
                    build_local_component_conservation_residual(
                        identity(0U),
                        2.0,
                        acc,
                        invalid_enum);
            },
            "unsupported incident-face local side");
    }

    (void)second;
}

void headers() {
    require(
        local_component_conservation_residual_header(),
        "local component conservation header probe failed");
}

using Test =
    std::pair<std::string_view, void (*)()>;

constexpr Test tests[]{
    {"multi_face_residual_and_blocks",
     multi_face_residual_and_blocks},
    {"jacobian_fresh_source_perturbation",
     jacobian_fresh_source_perturbation},
    {"zero_primal_nonzero_jacobian",
     zero_primal_nonzero_jacobian},
    {"invalid_inputs", invalid_inputs},
    {"headers", headers}};

} // namespace

int main(int argc, char** argv) {
    try {
        if (argc != 2) {
            throw std::invalid_argument(
                "one test name required");
        }
        for (const auto& [name, run] : tests) {
            if (name == argv[1]) {
                run();
                std::cout
                    << "[PASS] "
                    << name
                    << '\n';
                return 0;
            }
        }
        throw std::invalid_argument(
            "unknown test");
    } catch (const std::exception& error) {
        std::cerr
            << "[FAIL] "
            << error.what()
            << '\n';
        return 1;
    }
}
