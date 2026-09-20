#include <mpmc/flow_discretization/normalized_component_face_contribution.hpp>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <iostream>
#include <limits>
#include <source_location>
#include <stdexcept>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

bool normalized_component_face_contribution_header();

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
    bool owner) {
    const auto pivot =
        flow::NaturalVariableCompositionPivot3P::
            from_dependent_components(
                3U,
                owner
                    ? std::array<std::size_t, 3>{1U, 0U, 2U}
                    : std::array<std::size_t, 3>{0U, 2U, 1U});

    return {
        flow::NaturalVariableLayout3P{pivot},
        {"A", "B", "C"},
        owner ? 1.0e6 : 0.99e6,
        owner ? 350.0 : 351.0,
        owner
            ? std::array<double, 3>{0.20, 0.30, 0.50}
            : std::array<double, 3>{0.25, 0.35, 0.40},
        owner
            ? std::array<std::vector<double>, 3>{
                  std::vector<double>{0.10, 0.70, 0.20},
                  std::vector<double>{0.60, 0.20, 0.20},
                  std::vector<double>{0.20, 0.30, 0.50}}
            : std::array<std::vector<double>, 3>{
                  std::vector<double>{0.55, 0.25, 0.20},
                  std::vector<double>{0.25, 0.25, 0.50},
                  std::vector<double>{0.20, 0.60, 0.20}}};
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
source_scatter() {
    auto owner_identity = identity(true);
    auto neighbour_identity = identity(false);
    const std::size_t owner_q =
        owner_identity.layout.unknown_count();
    const std::size_t neighbour_q =
        neighbour_identity.layout.unknown_count();

    const std::vector<double>
        owner_rate{2.0, -1.0, 3.0};
    std::vector<double> neighbour_rate{
        -2.0, 1.0, -3.0};

    const auto owner_owner =
        component_jacobian(
            3U,
            owner_q,
            0.01);
    const auto owner_neighbour =
        component_jacobian(
            3U,
            neighbour_q,
            -0.02);

    std::vector<double> neighbour_owner(
        owner_owner.size());
    std::vector<double> neighbour_neighbour(
        owner_neighbour.size());
    for (std::size_t i = 0U;
         i < owner_owner.size();
         ++i) {
        neighbour_owner[i] =
            -owner_owner[i];
    }
    for (std::size_t i = 0U;
         i < owner_neighbour.size();
         ++i) {
        neighbour_neighbour[i] =
            -owner_neighbour[i];
    }

    double owner_total = 0.0;
    for (double value : owner_rate) {
        owner_total += value;
    }

    const auto owner_total_owner =
        total_gradient(
            owner_owner,
            3U,
            owner_q);
    const auto owner_total_neighbour =
        total_gradient(
            owner_neighbour,
            3U,
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
            mesh::LocalIndex::value_type{13}},
        {"A", "B", "C"},
        std::move(owner_identity),
        std::move(neighbour_identity),
        owner_rate,
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

void normalized_values_and_weighted_conservation() {
    constexpr double owner_volume = 2.0;
    constexpr double neighbour_volume = 5.0;

    const auto source = source_scatter();
    const auto result =
        fd::normalize_component_face_rate_by_bulk_volume(
            source,
            {owner_volume, neighbour_volume});

    require(
        result.face ==
                mesh::LocalIndex{
                    mesh::LocalIndex::value_type{13}} &&
            result.component_ids ==
                std::vector<std::string>{"A", "B", "C"} &&
            result.owner_state_identity
                    .layout
                    .composition_pivot()
                    .dependent_components() ==
                source.owner_state_identity
                    .layout
                    .composition_pivot()
                    .dependent_components() &&
            result.neighbour_state_identity
                    .layout
                    .composition_pivot()
                    .dependent_components() ==
                source.neighbour_state_identity
                    .layout
                    .composition_pivot()
                    .dependent_components(),
        "normalized contribution lost face/component/chart identity");

    near(
        result.bulk_volume.owner_bulk_volume_m3,
        owner_volume);
    near(
        result.bulk_volume.neighbour_bulk_volume_m3,
        neighbour_volume);

    bool observed_non_antisymmetric_normalized_value =
        false;
    for (std::size_t component = 0U;
         component < source.component_count();
         ++component) {
        near(
            result.owner_contribution(component),
            source.owner_rate(component) /
                owner_volume);
        near(
            result.neighbour_contribution(component),
            source.neighbour_rate(component) /
                neighbour_volume);

        near(
            owner_volume *
                    result.owner_contribution(component) +
                neighbour_volume *
                    result.neighbour_contribution(component),
            0.0,
            0.0,
            2.0e-15);

        if (result.owner_contribution(component) +
                result.neighbour_contribution(component) !=
            0.0) {
            observed_non_antisymmetric_normalized_value =
                true;
        }
    }
    require(
        observed_non_antisymmetric_normalized_value,
        "unequal cell volumes incorrectly preserved direct normalized antisymmetry");

    near(
        result.owner_total_contribution_mol_per_bulk_m3_s,
        source.owner_total_face_rate_mol_per_s /
            owner_volume);
    near(
        result.neighbour_total_contribution_mol_per_bulk_m3_s,
        source.neighbour_total_face_rate_mol_per_s /
            neighbour_volume);
    near(
        owner_volume *
                result.owner_total_contribution_mol_per_bulk_m3_s +
            neighbour_volume *
                result.neighbour_total_contribution_mol_per_bulk_m3_s,
        0.0,
        0.0,
        2.0e-15);
}

void normalized_jacobian() {
    constexpr double owner_volume = 2.0;
    constexpr double neighbour_volume = 5.0;

    const auto source = source_scatter();
    const auto result =
        fd::normalize_component_face_rate_by_bulk_volume(
            source,
            {owner_volume, neighbour_volume});

    const std::size_t owner_q =
        source.owner_state_identity.layout
            .unknown_count();
    const std::size_t neighbour_q =
        source.neighbour_state_identity.layout
            .unknown_count();

    for (std::size_t component = 0U;
         component < source.component_count();
         ++component) {
        for (std::size_t column = 0U;
             column < owner_q;
             ++column) {
            near(
                result
                    .d_owner_contribution_wrt_owner(
                        component,
                        column),
                source
                    .d_owner_rate_wrt_owner(
                        component,
                        column) /
                    owner_volume);
            near(
                result
                    .d_neighbour_contribution_wrt_owner(
                        component,
                        column),
                source
                    .d_neighbour_rate_wrt_owner(
                        component,
                        column) /
                    neighbour_volume);
            near(
                owner_volume *
                        result
                            .d_owner_contribution_wrt_owner(
                                component,
                                column) +
                    neighbour_volume *
                        result
                            .d_neighbour_contribution_wrt_owner(
                                component,
                                column),
                0.0,
                0.0,
                2.0e-15);
        }

        for (std::size_t column = 0U;
             column < neighbour_q;
             ++column) {
            near(
                result
                    .d_owner_contribution_wrt_neighbour(
                        component,
                        column),
                source
                    .d_owner_rate_wrt_neighbour(
                        component,
                        column) /
                    owner_volume);
            near(
                result
                    .d_neighbour_contribution_wrt_neighbour(
                        component,
                        column),
                source
                    .d_neighbour_rate_wrt_neighbour(
                        component,
                        column) /
                    neighbour_volume);
            near(
                owner_volume *
                        result
                            .d_owner_contribution_wrt_neighbour(
                                component,
                                column) +
                    neighbour_volume *
                        result
                            .d_neighbour_contribution_wrt_neighbour(
                                component,
                                column),
                0.0,
                0.0,
                2.0e-15);
        }
    }
}

fd::ConservativeComponentFaceRateScatterLinearization3D
perturb_source(
    fd::ConservativeComponentFaceRateScatterLinearization3D
        source,
    bool owner_column,
    std::size_t column,
    double delta) {
    const std::size_t n =
        source.component_count();
    const std::size_t columns =
        owner_column
            ? source.owner_state_identity.layout
                  .unknown_count()
            : source.neighbour_state_identity.layout
                  .unknown_count();

    source.owner_total_face_rate_mol_per_s = 0.0;
    source.neighbour_total_face_rate_mol_per_s = 0.0;

    for (std::size_t component = 0U;
         component < n;
         ++component) {
        const std::size_t index =
            component * columns +
            column;
        const double owner_derivative =
            owner_column
                ? source
                      .owner_row_owner_column_jacobian[
                          index]
                : source
                      .owner_row_neighbour_column_jacobian[
                          index];
        const double neighbour_derivative =
            owner_column
                ? source
                      .neighbour_row_owner_column_jacobian[
                          index]
                : source
                      .neighbour_row_neighbour_column_jacobian[
                          index];

        source.owner_component_face_rate_mol_per_s[
            component] +=
            owner_derivative *
            delta;
        source.neighbour_component_face_rate_mol_per_s[
            component] +=
            neighbour_derivative *
            delta;

        source.owner_total_face_rate_mol_per_s +=
            source.owner_component_face_rate_mol_per_s[
                component];
        source.neighbour_total_face_rate_mol_per_s +=
            source.neighbour_component_face_rate_mol_per_s[
                component];
    }

    return source;
}

void jacobian_fresh_source_perturbation() {
    constexpr fd::TwoCellBulkVolume3D volume{
        2.0, 5.0};
    constexpr double step = 1.0e-6;

    const auto source = source_scatter();
    const auto base =
        fd::normalize_component_face_rate_by_bulk_volume(
            source,
            volume);

    for (std::size_t column = 0U;
         column <
         source.owner_state_identity.layout
             .unknown_count();
         ++column) {
        const auto plus =
            fd::normalize_component_face_rate_by_bulk_volume(
                perturb_source(
                    source,
                    true,
                    column,
                    step),
                volume);
        const auto minus =
            fd::normalize_component_face_rate_by_bulk_volume(
                perturb_source(
                    source,
                    true,
                    column,
                    -step),
                volume);

        for (std::size_t component = 0U;
             component < source.component_count();
             ++component) {
            near(
                base.d_owner_contribution_wrt_owner(
                    component,
                    column),
                (plus.owner_contribution(component) -
                 minus.owner_contribution(component)) /
                    (2.0 * step),
                2.0e-9,
                2.0e-9);
            near(
                base.d_neighbour_contribution_wrt_owner(
                    component,
                    column),
                (plus.neighbour_contribution(component) -
                 minus.neighbour_contribution(component)) /
                    (2.0 * step),
                2.0e-9,
                2.0e-9);
        }
    }

    for (std::size_t column = 0U;
         column <
         source.neighbour_state_identity.layout
             .unknown_count();
         ++column) {
        const auto plus =
            fd::normalize_component_face_rate_by_bulk_volume(
                perturb_source(
                    source,
                    false,
                    column,
                    step),
                volume);
        const auto minus =
            fd::normalize_component_face_rate_by_bulk_volume(
                perturb_source(
                    source,
                    false,
                    column,
                    -step),
                volume);

        for (std::size_t component = 0U;
             component < source.component_count();
             ++component) {
            near(
                base.d_owner_contribution_wrt_neighbour(
                    component,
                    column),
                (plus.owner_contribution(component) -
                 minus.owner_contribution(component)) /
                    (2.0 * step),
                2.0e-9,
                2.0e-9);
            near(
                base.d_neighbour_contribution_wrt_neighbour(
                    component,
                    column),
                (plus.neighbour_contribution(component) -
                 minus.neighbour_contribution(component)) /
                    (2.0 * step),
                2.0e-9,
                2.0e-9);
        }
    }
}

void zero_primal_nonzero_jacobian() {
    auto source = source_scatter();
    std::fill(
        source.owner_component_face_rate_mol_per_s.begin(),
        source.owner_component_face_rate_mol_per_s.end(),
        0.0);
    std::fill(
        source.neighbour_component_face_rate_mol_per_s.begin(),
        source.neighbour_component_face_rate_mol_per_s.end(),
        0.0);
    source.owner_total_face_rate_mol_per_s = 0.0;
    source.neighbour_total_face_rate_mol_per_s = 0.0;

    const auto result =
        fd::normalize_component_face_rate_by_bulk_volume(
            source,
            {2.0, 5.0});

    require(
        std::all_of(
            result
                .owner_component_contribution_mol_per_bulk_m3_s
                .begin(),
            result
                .owner_component_contribution_mol_per_bulk_m3_s
                .end(),
            [](double value) {
                return value == 0.0;
            }) &&
            std::all_of(
                result
                    .neighbour_component_contribution_mol_per_bulk_m3_s
                    .begin(),
                result
                    .neighbour_component_contribution_mol_per_bulk_m3_s
                    .end(),
                [](double value) {
                    return value == 0.0;
                }),
        "zero primal face rate did not normalize to zero");

    require(
        std::any_of(
            result.owner_row_owner_column_jacobian.begin(),
            result.owner_row_owner_column_jacobian.end(),
            [](double value) {
                return value != 0.0;
            }) &&
            std::any_of(
                result.neighbour_row_neighbour_column_jacobian.begin(),
                result.neighbour_row_neighbour_column_jacobian.end(),
                [](double value) {
                    return value != 0.0;
                }),
        "zero primal face rate incorrectly erased normalized Jacobian");
}

void invalid_inputs() {
    const auto source = source_scatter();

    expect_invalid(
        [&] {
            (void)fd::
                normalize_component_face_rate_by_bulk_volume(
                    source,
                    {0.0, 5.0});
        },
        "owner cell bulk volume");

    expect_invalid(
        [&] {
            (void)fd::
                normalize_component_face_rate_by_bulk_volume(
                    source,
                    {2.0,
                     std::numeric_limits<double>::
                         quiet_NaN()});
        },
        "neighbour cell bulk volume");

    auto nonconservative = source;
    nonconservative
        .neighbour_component_face_rate_mol_per_s[0] +=
        1.0;
    expect_invalid(
        [&] {
            (void)fd::
                normalize_component_face_rate_by_bulk_volume(
                    nonconservative,
                    {2.0, 5.0});
        },
        "not exactly conservative");

    auto wrong_total = source;
    wrong_total.owner_total_face_rate_mol_per_s +=
        1.0;
    wrong_total.neighbour_total_face_rate_mol_per_s -=
        1.0;
    expect_invalid(
        [&] {
            (void)fd::
                normalize_component_face_rate_by_bulk_volume(
                    wrong_total,
                    {2.0, 5.0});
        },
        "does not close");

    auto wrong_total_derivative = source;
    wrong_total_derivative
        .owner_total_row_owner_column_gradient[0] +=
        1.0;
    wrong_total_derivative
        .neighbour_total_row_owner_column_gradient[0] -=
        1.0;
    expect_invalid(
        [&] {
            (void)fd::
                normalize_component_face_rate_by_bulk_volume(
                    wrong_total_derivative,
                    {2.0, 5.0});
        },
        "does not close");
}

void headers() {
    require(
        normalized_component_face_contribution_header(),
        "normalized component face contribution header probe failed");
}

using Test =
    std::pair<std::string_view, void (*)()>;

constexpr Test tests[]{
    {"normalized_values_and_weighted_conservation",
     normalized_values_and_weighted_conservation},
    {"normalized_jacobian", normalized_jacobian},
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
