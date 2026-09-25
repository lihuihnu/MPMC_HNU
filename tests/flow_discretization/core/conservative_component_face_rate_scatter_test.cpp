#include <mpmc/flow_discretization/conservative_component_face_rate_scatter.hpp>

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

bool conservative_component_face_rate_scatter_header();

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
    double relative = 2.0e-13,
    double absolute = 2.0e-13,
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

std::vector<double> jacobian(
    std::size_t components,
    std::size_t columns,
    double scale) {
    std::vector<double> result(
        components * columns,
        0.0);
    for (std::size_t component = 0U;
         component < components;
         ++component) {
        for (std::size_t column = 0U;
             column < columns;
             ++column) {
            result[
                component * columns +
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
    const std::vector<double>& component_jacobian,
    std::size_t components,
    std::size_t columns) {
    std::vector<double> result(
        columns,
        0.0);
    for (std::size_t column = 0U;
         column < columns;
         ++column) {
        for (std::size_t component = 0U;
             component < components;
             ++component) {
            result[column] +=
                component_jacobian[
                    component * columns +
                    column];
        }
    }
    return result;
}

fd::MaterializedTpfaInternalFaceComponentMolarFluxLinearization3D
source_flux() {
    auto owner = identity(true);
    auto neighbour = identity(false);
    const std::size_t owner_q =
        owner.layout.unknown_count();
    const std::size_t neighbour_q =
        neighbour.layout.unknown_count();

    const std::vector<double>
        component_flux{1.25, -0.50, 2.00};
    const auto owner_jacobian =
        jacobian(3U, owner_q, 0.01);
    const auto neighbour_jacobian =
        jacobian(3U, neighbour_q, -0.02);

    double total = 0.0;
    for (double value : component_flux) {
        total += value;
    }

    return {
        mesh::LocalIndex{
            mesh::LocalIndex::value_type{11}},
        {"A", "B", "C"},
        std::move(owner),
        std::move(neighbour),
        {},
        component_flux,
        owner_jacobian,
        neighbour_jacobian,
        total,
        total_gradient(
            owner_jacobian,
            3U,
            owner_q),
        total_gradient(
            neighbour_jacobian,
            3U,
            neighbour_q)};
}

void conservative_values() {
    const auto source = source_flux();
    const auto scatter =
        fd::scatter_component_molar_face_flux_conservatively(
            source);

    require(
        scatter.face ==
                mesh::LocalIndex{
                    mesh::LocalIndex::value_type{11}} &&
            scatter.component_ids ==
                std::vector<std::string>{"A", "B", "C"} &&
            scatter.owner_state_identity
                    .layout
                    .composition_pivot()
                    .dependent_components() ==
                source.owner_state_identity
                    .layout
                    .composition_pivot()
                    .dependent_components() &&
            scatter.neighbour_state_identity
                    .layout
                    .composition_pivot()
                    .dependent_components() ==
                source.neighbour_state_identity
                    .layout
                    .composition_pivot()
                    .dependent_components(),
        "scatter lost face/component/chart identity");

    for (std::size_t component = 0U;
         component < source.component_count();
         ++component) {
        near(
            scatter.owner_rate(component),
            source.component_flux(component));
        near(
            scatter.neighbour_rate(component),
            -source.component_flux(component));
        near(
            scatter.owner_rate(component) +
                scatter.neighbour_rate(component),
            0.0,
            0.0,
            0.0);
    }

    near(
        scatter.owner_total_face_rate_mol_per_s,
        source.total_molar_flux_mol_per_s);
    near(
        scatter.neighbour_total_face_rate_mol_per_s,
        -source.total_molar_flux_mol_per_s);
    near(
        scatter.owner_total_face_rate_mol_per_s +
            scatter.neighbour_total_face_rate_mol_per_s,
        0.0,
        0.0,
        0.0);
}

void strict_jacobian_antisymmetry() {
    const auto source = source_flux();
    const auto scatter =
        fd::scatter_component_molar_face_flux_conservatively(
            source);
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
            const double source_value =
                source.d_component_flux_owner(
                    component,
                    column);
            near(
                scatter.d_owner_rate_wrt_owner(
                    component,
                    column),
                source_value);
            near(
                scatter.d_neighbour_rate_wrt_owner(
                    component,
                    column),
                -source_value);
            near(
                scatter.d_owner_rate_wrt_owner(
                    component,
                    column) +
                    scatter.d_neighbour_rate_wrt_owner(
                        component,
                        column),
                0.0,
                0.0,
                0.0);
        }

        for (std::size_t column = 0U;
             column < neighbour_q;
             ++column) {
            const double source_value =
                source.d_component_flux_neighbour(
                    component,
                    column);
            near(
                scatter.d_owner_rate_wrt_neighbour(
                    component,
                    column),
                source_value);
            near(
                scatter.d_neighbour_rate_wrt_neighbour(
                    component,
                    column),
                -source_value);
            near(
                scatter.d_owner_rate_wrt_neighbour(
                    component,
                    column) +
                    scatter.d_neighbour_rate_wrt_neighbour(
                        component,
                        column),
                0.0,
                0.0,
                0.0);
        }
    }

    for (std::size_t column = 0U;
         column < owner_q;
         ++column) {
        near(
            scatter
                .owner_total_row_owner_column_gradient[
                    column],
            source.owner_total_molar_flux_gradient[
                column]);
        near(
            scatter
                    .owner_total_row_owner_column_gradient[
                        column] +
                scatter
                    .neighbour_total_row_owner_column_gradient[
                        column],
            0.0,
            0.0,
            0.0);
    }

    for (std::size_t column = 0U;
         column < neighbour_q;
         ++column) {
        near(
            scatter
                .owner_total_row_neighbour_column_gradient[
                    column],
            source.neighbour_total_molar_flux_gradient[
                column]);
        near(
            scatter
                    .owner_total_row_neighbour_column_gradient[
                        column] +
                scatter
                    .neighbour_total_row_neighbour_column_gradient[
                        column],
            0.0,
            0.0,
            0.0);
    }
}

fd::MaterializedTpfaInternalFaceComponentMolarFluxLinearization3D
perturb_source(
    fd::MaterializedTpfaInternalFaceComponentMolarFluxLinearization3D source,
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
    const auto& jacobian_values =
        owner_column
            ? source.owner_component_flux_jacobian
            : source.neighbour_component_flux_jacobian;

    source.total_molar_flux_mol_per_s = 0.0;
    for (std::size_t component = 0U;
         component < n;
         ++component) {
        source.component_molar_flux_mol_per_s[
            component] +=
            jacobian_values[
                component * columns +
                column] *
            delta;
        source.total_molar_flux_mol_per_s +=
            source.component_molar_flux_mol_per_s[
                component];
    }
    return source;
}

void jacobian_fresh_source_perturbation() {
    const auto source = source_flux();
    const auto base =
        fd::scatter_component_molar_face_flux_conservatively(
            source);
    constexpr double step = 1.0e-6;

    for (std::size_t column = 0U;
         column <
         source.owner_state_identity.layout
             .unknown_count();
         ++column) {
        const auto plus =
            fd::scatter_component_molar_face_flux_conservatively(
                perturb_source(
                    source,
                    true,
                    column,
                    step));
        const auto minus =
            fd::scatter_component_molar_face_flux_conservatively(
                perturb_source(
                    source,
                    true,
                    column,
                    -step));

        for (std::size_t component = 0U;
             component < source.component_count();
             ++component) {
            const double fd_owner =
                (plus.owner_rate(component) -
                 minus.owner_rate(component)) /
                (2.0 * step);
            const double fd_neighbour =
                (plus.neighbour_rate(component) -
                 minus.neighbour_rate(component)) /
                (2.0 * step);
            near(
                base.d_owner_rate_wrt_owner(
                    component,
                    column),
                fd_owner,
                2.0e-9,
                2.0e-9);
            near(
                base.d_neighbour_rate_wrt_owner(
                    component,
                    column),
                fd_neighbour,
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
            fd::scatter_component_molar_face_flux_conservatively(
                perturb_source(
                    source,
                    false,
                    column,
                    step));
        const auto minus =
            fd::scatter_component_molar_face_flux_conservatively(
                perturb_source(
                    source,
                    false,
                    column,
                    -step));

        for (std::size_t component = 0U;
             component < source.component_count();
             ++component) {
            const double fd_owner =
                (plus.owner_rate(component) -
                 minus.owner_rate(component)) /
                (2.0 * step);
            const double fd_neighbour =
                (plus.neighbour_rate(component) -
                 minus.neighbour_rate(component)) /
                (2.0 * step);
            near(
                base.d_owner_rate_wrt_neighbour(
                    component,
                    column),
                fd_owner,
                2.0e-9,
                2.0e-9);
            near(
                base.d_neighbour_rate_wrt_neighbour(
                    component,
                    column),
                fd_neighbour,
                2.0e-9,
                2.0e-9);
        }
    }
}

void zero_primal_nonzero_jacobian() {
    auto source = source_flux();
    std::fill(
        source.component_molar_flux_mol_per_s.begin(),
        source.component_molar_flux_mol_per_s.end(),
        0.0);
    source.total_molar_flux_mol_per_s = 0.0;

    const auto scatter =
        fd::scatter_component_molar_face_flux_conservatively(
            source);

    require(
        std::all_of(
            scatter.owner_component_face_rate_mol_per_s.begin(),
            scatter.owner_component_face_rate_mol_per_s.end(),
            [](double value) {
                return value == 0.0;
            }) &&
            std::all_of(
                scatter.neighbour_component_face_rate_mol_per_s.begin(),
                scatter.neighbour_component_face_rate_mol_per_s.end(),
                [](double value) {
                    return value == 0.0;
                }),
        "zero primal component face rate did not scatter to zero");

    require(
        std::any_of(
            scatter.owner_row_owner_column_jacobian.begin(),
            scatter.owner_row_owner_column_jacobian.end(),
            [](double value) {
                return value != 0.0;
            }) &&
            std::any_of(
                scatter.neighbour_row_neighbour_column_jacobian.begin(),
                scatter.neighbour_row_neighbour_column_jacobian.end(),
                [](double value) {
                    return value != 0.0;
                }),
        "zero primal component face rate incorrectly erased Jacobian");
}

void invalid_source() {
    auto wrong_identity = source_flux();
    wrong_identity.neighbour_state_identity
        .component_ids[1] = "X";
    expect_invalid(
        [&] {
            (void)fd::
                scatter_component_molar_face_flux_conservatively(
                    wrong_identity);
        },
        "identity/chart");

    auto wrong_shape = source_flux();
    wrong_shape.owner_component_flux_jacobian.pop_back();
    expect_invalid(
        [&] {
            (void)fd::
                scatter_component_molar_face_flux_conservatively(
                    wrong_shape);
        },
        "shape");

    auto nonfinite = source_flux();
    nonfinite
        .neighbour_component_flux_jacobian[0] =
        std::numeric_limits<double>::infinity();
    expect_invalid(
        [&] {
            (void)fd::
                scatter_component_molar_face_flux_conservatively(
                    nonfinite);
        },
        "non-finite");

    auto wrong_total = source_flux();
    wrong_total.total_molar_flux_mol_per_s +=
        1.0;
    expect_invalid(
        [&] {
            (void)fd::
                scatter_component_molar_face_flux_conservatively(
                    wrong_total);
        },
        "does not close");

    auto wrong_derivative = source_flux();
    wrong_derivative
        .owner_total_molar_flux_gradient[0] +=
        1.0;
    expect_invalid(
        [&] {
            (void)fd::
                scatter_component_molar_face_flux_conservatively(
                    wrong_derivative);
        },
        "Jacobian does not close");
}

void headers() {
    require(
        conservative_component_face_rate_scatter_header(),
        "conservative component face-rate scatter header probe failed");
}

using Test =
    std::pair<std::string_view, void (*)()>;

constexpr Test tests[]{
    {"conservative_values", conservative_values},
    {"strict_jacobian_antisymmetry", strict_jacobian_antisymmetry},
    {"jacobian_fresh_source_perturbation", jacobian_fresh_source_perturbation},
    {"zero_primal_nonzero_jacobian", zero_primal_nonzero_jacobian},
    {"invalid_source", invalid_source},
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
