#include <mpmc/flow_discretization/tpfa_component_molar_flux.hpp>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <iostream>
#include <limits>
#include <optional>
#include <source_location>
#include <span>
#include <stdexcept>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

bool tpfa_component_molar_flux_header();

namespace {

namespace disc = mpmc::discretization;
namespace flow = mpmc::flow;
namespace fd = mpmc::flow_discretization;
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
    double relative = 3.0e-10,
    double absolute = 3.0e-10,
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

struct SideChart {
    flow::NaturalVariableLayout3P layout;
    std::vector<double> base_inputs;
};

std::vector<double> independent_values(
    std::span<const double> full,
    std::size_t dependent) {
    std::vector<double> result;
    result.reserve(full.size() - 1U);
    for (std::size_t component = 0U;
         component < full.size();
         ++component) {
        if (component != dependent) {
            result.push_back(full[component]);
        }
    }
    return result;
}

std::array<std::vector<double>, 3>
base_compositions(bool owner) {
    if (owner) {
        return {
            std::vector<double>{0.10, 0.70, 0.20},
            std::vector<double>{0.60, 0.20, 0.20},
            std::vector<double>{0.20, 0.30, 0.50}};
    }
    return {
        std::vector<double>{0.55, 0.25, 0.20},
        std::vector<double>{0.25, 0.25, 0.50},
        std::vector<double>{0.20, 0.60, 0.20}};
}

SideChart side_chart(bool owner) {
    const auto compositions =
        base_compositions(owner);
    const auto pivot =
        flow::NaturalVariableCompositionPivot3P::
            select(compositions);
    flow::NaturalVariableLayout3P layout{pivot};

    std::vector<double> inputs(
        layout.unknown_count(),
        0.0);
    inputs[layout.pressure_unknown_index()] =
        owner ? 1.0e6 : 0.99e6;
    inputs[layout.temperature_unknown_index()] =
        owner ? 350.0 : 351.0;
    inputs[*layout.independent_saturation_unknown_index(
        flow::PhaseSlot3::phase0)] =
        owner ? 0.20 : 0.25;
    inputs[*layout.independent_saturation_unknown_index(
        flow::PhaseSlot3::phase1)] =
        owner ? 0.30 : 0.35;

    for (std::size_t phase = 0U;
         phase < 3U;
         ++phase) {
        const auto slot =
            static_cast<flow::PhaseSlot3>(
                phase);
        for (std::size_t component = 0U;
             component < 3U;
             ++component) {
            const auto column =
                layout
                    .independent_composition_unknown_index(
                        slot,
                        component);
            if (column) {
                inputs[*column] =
                    compositions[phase][component];
            }
        }
    }
    return {std::move(layout), std::move(inputs)};
}

std::array<std::vector<double>, 3>
reconstruct_compositions(
    const flow::NaturalVariableLayout3P& layout,
    std::span<const double> inputs) {
    const std::size_t n =
        layout.component_count();
    std::array<std::vector<double>, 3> result;

    for (std::size_t phase = 0U;
         phase < 3U;
         ++phase) {
        const auto slot =
            static_cast<flow::PhaseSlot3>(
                phase);
        const std::size_t dependent =
            layout.dependent_composition_component(
                slot);
        result[phase].assign(n, 0.0);
        double sum = 0.0;
        for (std::size_t component = 0U;
             component < n;
             ++component) {
            const auto column =
                layout
                    .independent_composition_unknown_index(
                        slot,
                        component);
            if (!column) {
                continue;
            }
            result[phase][component] =
                inputs[*column];
            sum += inputs[*column];
        }
        result[phase][dependent] =
            1.0 - sum;
    }
    return result;
}

std::array<std::vector<double>, 3>
density_gradients(
    const flow::NaturalVariableLayout3P& layout,
    bool owner) {
    std::array<std::vector<double>, 3> result;
    const double side_scale =
        owner ? 1.0 : 1.3;
    for (std::size_t phase = 0U;
         phase < 3U;
         ++phase) {
        result[phase].resize(
            layout.unknown_count());
        for (std::size_t column = 0U;
             column < layout.unknown_count();
             ++column) {
            result[phase][column] =
                side_scale *
                1.0e-3 *
                static_cast<double>(
                    1U +
                    phase * 10U +
                    column);
        }
    }
    return result;
}

std::array<double, 3> molar_density(
    bool owner,
    const SideChart& chart,
    std::span<const double> inputs) {
    const auto gradient =
        density_gradients(
            chart.layout,
            owner);
    const std::array<double, 3> base{
        owner ? 1000.0 : 1100.0,
        owner ? 2000.0 : 2100.0,
        owner ? 3000.0 : 3100.0};

    std::array<double, 3> result =
        base;
    for (std::size_t phase = 0U;
         phase < 3U;
         ++phase) {
        for (std::size_t column = 0U;
             column < inputs.size();
             ++column) {
            result[phase] +=
                gradient[phase][column] *
                (inputs[column] -
                 chart.base_inputs[column]);
        }
    }
    return result;
}

flow::NaturalVariableCellState3P make_state(
    bool owner,
    const SideChart& chart,
    std::span<const double> inputs) {
    const auto compositions =
        reconstruct_compositions(
            chart.layout,
            inputs);
    const auto density =
        molar_density(
            owner,
            chart,
            inputs);

    flow::NaturalVariableCellStateInput3P input;
    input.component_ids = {"A", "B", "C"};
    input.reference_pressure_pa =
        inputs[
            chart.layout
                .pressure_unknown_index()];
    input.temperature_k =
        inputs[
            chart.layout
                .temperature_unknown_index()];
    input.independent_saturations = {
        inputs[
            *chart.layout
                 .independent_saturation_unknown_index(
                     flow::PhaseSlot3::phase0)],
        inputs[
            *chart.layout
                 .independent_saturation_unknown_index(
                     flow::PhaseSlot3::phase1)]};
    input.composition_pivot =
        chart.layout.composition_pivot();

    for (std::size_t phase = 0U;
         phase < 3U;
         ++phase) {
        const auto slot =
            static_cast<flow::PhaseSlot3>(
                phase);
        input.independent_phase_compositions[phase] =
            independent_values(
                compositions[phase],
                chart.layout
                    .dependent_composition_component(
                        slot));
        input.phase_properties[phase] = {
            density[phase],
            700.0 +
                20.0 *
                    static_cast<double>(phase),
            1.0e-3 +
                2.0e-4 *
                    static_cast<double>(phase),
            2.0e5,
            1.5e5};
    }

    return flow::NaturalVariableCellState3P::create(
        std::move(input));
}

flow::PhaseMolarDensityNaturalVariableLinearization3P
density_linearization(
    bool owner,
    const SideChart& chart,
    std::span<const double> inputs) {
    return {
        chart.layout,
        molar_density(
            owner,
            chart,
            inputs),
        density_gradients(
            chart.layout,
            owner)};
}

flow::NaturalVariableStateIdentity3P identity(
    const flow::NaturalVariableCellState3P& state) {
    flow::NaturalVariableStateIdentity3P result{
        state.layout(),
        std::vector<std::string>{
            state.component_ids().begin(),
            state.component_ids().end()},
        state.reference_pressure_pa(),
        state.temperature_k(),
        {},
        {}};

    for (std::size_t phase = 0U;
         phase < 3U;
         ++phase) {
        const auto slot =
            static_cast<flow::PhaseSlot3>(
                phase);
        result.saturation[phase] =
            state.phase_saturation(slot);
        const auto x =
            state.phase_composition(slot);
        result.phase_composition[phase].assign(
            x.begin(),
            x.end());
    }
    return result;
}

disc::CombinedTransmissibilityAdmissibility3D
direct_admissibility() {
    return {
        disc::
            CombinedTransmissibilityAdmissibilityDisposition3D::
                direct_normal_projection_k_orthogonal_candidate,
        {
            disc::
                TransmissibilityGeometryDisposition3D::
                    direct_normal_projection_allowed,
            0.0,
            0.05},
        {
            disc::
                KOrthogonalityDisposition3D::
                    k_orthogonal_within_policy,
            std::nullopt,
            std::nullopt,
            {std::nullopt, std::nullopt},
            0.05}};
}

disc::TpfaInternalFaceTransmissibilityEntry3D
materialized_entry() {
    return {
        mesh::LocalIndex{
            mesh::LocalIndex::value_type{9}},
        disc::
            TpfaInternalFaceTransmissibilityDisposition3D::
                materialized,
        direct_admissibility(),
        disc::TpfaStaticFaceTransmissibility3D{
            disc::
                TpfaStaticFaceTransmissibilityDisposition3D::
                    positive_harmonic_combination,
            1.0,
            2.0e-12}};
}

struct LocalFaceFixture {
    flow::NaturalVariableCellState3P owner_state;
    flow::PhaseMolarDensityNaturalVariableLinearization3P
        owner_density;
    flow::NaturalVariableCellState3P neighbour_state;
    flow::PhaseMolarDensityNaturalVariableLinearization3P
        neighbour_density;
    fd::MaterializedTpfaInternalFacePhaseDarcyFluxLinearization3D
        phase_flux;
};

LocalFaceFixture build_fixture(
    const SideChart& owner_chart,
    std::span<const double> owner_inputs,
    const SideChart& neighbour_chart,
    std::span<const double> neighbour_inputs) {
    auto owner_state =
        make_state(
            true,
            owner_chart,
            owner_inputs);
    auto neighbour_state =
        make_state(
            false,
            neighbour_chart,
            neighbour_inputs);

    const std::size_t owner_q =
        owner_chart.layout.unknown_count();
    const std::size_t neighbour_q =
        neighbour_chart.layout.unknown_count();

    flow::TwoCellPhasePotentialUpwindLinearization3P
        potential{
            flow::FacePhaseDensityPolicy3P::
                arithmetic_mean_owner_neighbour,
            {0.0, 0.0, 0.0},
            {1.0, 0.0, 0.0},
            identity(owner_state),
            identity(neighbour_state),
            {}};

    const auto owner_delta =
        [&] {
            std::vector<double> value(owner_q);
            for (std::size_t column = 0U;
                 column < owner_q;
                 ++column) {
                value[column] =
                    owner_inputs[column] -
                    owner_chart.base_inputs[column];
            }
            return value;
        }();
    const auto neighbour_delta =
        [&] {
            std::vector<double> value(neighbour_q);
            for (std::size_t column = 0U;
                 column < neighbour_q;
                 ++column) {
                value[column] =
                    neighbour_inputs[column] -
                    neighbour_chart.base_inputs[column];
            }
            return value;
        }();

    for (std::size_t phase = 0U;
         phase < 3U;
         ++phase) {
        auto& entry =
            potential.phase[phase];
        entry.face_mass_density_kg_per_m3 =
            1.0;
        entry.owner_face_density_gradient.assign(
            owner_q,
            0.0);
        entry.neighbour_face_density_gradient.assign(
            neighbour_q,
            0.0);
        entry.gravity_projection_m2_per_s2 =
            0.0;
        entry.gravity_pressure_difference_pa =
            0.0;

        entry.owner_phase_potential_gradient.assign(
            owner_q,
            0.0);
        entry.neighbour_phase_potential_gradient.assign(
            neighbour_q,
            0.0);
        entry.owner_upwind_mobility_gradient.assign(
            owner_q,
            0.0);
        entry.neighbour_upwind_mobility_gradient.assign(
            neighbour_q,
            0.0);

        if (phase == 0U) {
            entry.phase_potential_difference_pa =
                10000.0;
            for (std::size_t column = 0U;
                 column < owner_q;
                 ++column) {
                entry.owner_phase_potential_gradient[column] =
                    0.5 *
                    static_cast<double>(column + 1U);
                entry.phase_potential_difference_pa +=
                    entry.owner_phase_potential_gradient[column] *
                    owner_delta[column];
            }
            for (std::size_t column = 0U;
                 column < neighbour_q;
                 ++column) {
                entry.neighbour_phase_potential_gradient[column] =
                    0.75 *
                    static_cast<double>(column + 1U);
                entry.phase_potential_difference_pa +=
                    entry.neighbour_phase_potential_gradient[column] *
                    neighbour_delta[column];
                entry.neighbour_upwind_mobility_gradient[column] =
                    1.0e-5 *
                    static_cast<double>(column + 1U);
            }
            entry.upwind_selection =
                flow::UpwindCellSelection3P::
                    neighbour_positive_phase_potential;
            entry.upwind_mobility_per_pa_s =
                4.0;
            for (std::size_t column = 0U;
                 column < neighbour_q;
                 ++column) {
                entry.upwind_mobility_per_pa_s +=
                    entry.neighbour_upwind_mobility_gradient[column] *
                    neighbour_delta[column];
            }
        } else if (phase == 1U) {
            entry.phase_potential_difference_pa =
                -5000.0;
            for (std::size_t column = 0U;
                 column < owner_q;
                 ++column) {
                entry.owner_phase_potential_gradient[column] =
                    -0.4 *
                    static_cast<double>(column + 1U);
                entry.phase_potential_difference_pa +=
                    entry.owner_phase_potential_gradient[column] *
                    owner_delta[column];
                entry.owner_upwind_mobility_gradient[column] =
                    1.5e-5 *
                    static_cast<double>(column + 1U);
            }
            for (std::size_t column = 0U;
                 column < neighbour_q;
                 ++column) {
                entry.neighbour_phase_potential_gradient[column] =
                    0.25 *
                    static_cast<double>(column + 1U);
                entry.phase_potential_difference_pa +=
                    entry.neighbour_phase_potential_gradient[column] *
                    neighbour_delta[column];
            }
            entry.upwind_selection =
                flow::UpwindCellSelection3P::
                    owner_negative_phase_potential;
            entry.upwind_mobility_per_pa_s =
                2.0;
            for (std::size_t column = 0U;
                 column < owner_q;
                 ++column) {
                entry.upwind_mobility_per_pa_s +=
                    entry.owner_upwind_mobility_gradient[column] *
                    owner_delta[column];
            }
        } else {
            entry.phase_potential_difference_pa =
                0.0;
            entry.upwind_selection =
                flow::UpwindCellSelection3P::
                    owner_exact_zero_tie;
            entry.upwind_mobility_per_pa_s =
                3.0;
            for (std::size_t column = 0U;
                 column < owner_q;
                 ++column) {
                entry.owner_upwind_mobility_gradient[column] =
                    1.0e-5 *
                    static_cast<double>(column + 1U);
                entry.upwind_mobility_per_pa_s +=
                    entry.owner_upwind_mobility_gradient[column] *
                    owner_delta[column];
            }
        }
    }

    auto phase_flux =
        fd::build_materialized_tpfa_internal_face_phase_darcy_flux(
            materialized_entry(),
            potential);

    return {
        std::move(owner_state),
        density_linearization(
            true,
            owner_chart,
            owner_inputs),
        std::move(neighbour_state),
        density_linearization(
            false,
            neighbour_chart,
            neighbour_inputs),
        std::move(phase_flux)};
}

fd::MaterializedTpfaInternalFaceComponentMolarFluxLinearization3D
build_component_flux(
    const SideChart& owner_chart,
    std::span<const double> owner_inputs,
    const SideChart& neighbour_chart,
    std::span<const double> neighbour_inputs) {
    auto fixture =
        build_fixture(
            owner_chart,
            owner_inputs,
            neighbour_chart,
            neighbour_inputs);
    return fd::
        build_materialized_tpfa_internal_face_component_molar_flux(
            fixture.phase_flux,
            fixture.owner_state,
            fixture.owner_density,
            fixture.neighbour_state,
            fixture.neighbour_density);
}

void upwind_phase_molar_content() {
    const auto owner_chart =
        side_chart(true);
    const auto neighbour_chart =
        side_chart(false);
    auto fixture =
        build_fixture(
            owner_chart,
            owner_chart.base_inputs,
            neighbour_chart,
            neighbour_chart.base_inputs);

    const auto result =
        fd::
            build_materialized_tpfa_internal_face_component_molar_flux(
                fixture.phase_flux,
                fixture.owner_state,
                fixture.owner_density,
                fixture.neighbour_state,
                fixture.neighbour_density);

    require(
        result.phase_molar_content[0].upwind_selection ==
                flow::UpwindCellSelection3P::
                    neighbour_positive_phase_potential &&
            result.phase_molar_content[1].upwind_selection ==
                flow::UpwindCellSelection3P::
                    owner_negative_phase_potential &&
            result.phase_molar_content[2].upwind_selection ==
                flow::UpwindCellSelection3P::
                    owner_exact_zero_tie,
        "component molar content did not reuse phase upwind selection");

    const std::array<bool, 3>
        owner_upwind{false, true, true};
    for (std::size_t phase = 0U;
         phase < 3U;
         ++phase) {
        const auto& upstream_state =
            owner_upwind[phase]
                ? fixture.owner_state
                : fixture.neighbour_state;
        const auto& upstream_density =
            owner_upwind[phase]
                ? fixture.owner_density
                : fixture.neighbour_density;
        const auto x =
            upstream_state.phase_composition(
                static_cast<flow::PhaseSlot3>(
                    phase));
        const auto& content =
            result.phase_molar_content[phase];

        near(
            content.upwind_molar_density_mol_per_m3,
            upstream_density
                .molar_density_mol_per_m3[phase]);
        double sum = 0.0;
        for (std::size_t component = 0U;
             component < 3U;
             ++component) {
            near(
                content.upwind_mole_fraction[component],
                x[component]);
            near(
                content
                    .component_molar_content_mol_per_m3
                    [component],
                upstream_density
                    .molar_density_mol_per_m3[phase] *
                    x[component]);
            sum +=
                content
                    .component_molar_content_mol_per_m3
                    [component];
        }
        near(
            sum,
            content.upwind_molar_density_mol_per_m3);

        if (owner_upwind[phase]) {
            require(
                std::all_of(
                    content
                        .neighbour_component_molar_content_jacobian
                        .begin(),
                    content
                        .neighbour_component_molar_content_jacobian
                        .end(),
                    [](double value) {
                        return value == 0.0;
                    }),
                "owner-upwind phase acquired neighbour molar-content derivatives");
        } else {
            require(
                std::all_of(
                    content
                        .owner_component_molar_content_jacobian
                        .begin(),
                    content
                        .owner_component_molar_content_jacobian
                        .end(),
                    [](double value) {
                        return value == 0.0;
                    }),
                "neighbour-upwind phase acquired owner molar-content derivatives");
        }
    }
}

void component_flux_and_closure() {
    const auto owner_chart =
        side_chart(true);
    const auto neighbour_chart =
        side_chart(false);
    auto fixture =
        build_fixture(
            owner_chart,
            owner_chart.base_inputs,
            neighbour_chart,
            neighbour_chart.base_inputs);

    const auto result =
        fd::
            build_materialized_tpfa_internal_face_component_molar_flux(
                fixture.phase_flux,
                fixture.owner_state,
                fixture.owner_density,
                fixture.neighbour_state,
                fixture.neighbour_density);

    std::array<double, 3> expected{};
    double expected_total = 0.0;
    for (std::size_t phase = 0U;
         phase < 3U;
         ++phase) {
        const auto& content =
            result.phase_molar_content[phase];
        const double phase_flux =
            fixture.phase_flux.phase[phase]
                .volumetric_flux_m3_per_s;
        expected_total +=
            content.upwind_molar_density_mol_per_m3 *
            phase_flux;
        for (std::size_t component = 0U;
             component < 3U;
             ++component) {
            expected[component] +=
                content
                    .component_molar_content_mol_per_m3
                    [component] *
                phase_flux;
        }
    }

    double component_sum = 0.0;
    for (std::size_t component = 0U;
         component < 3U;
         ++component) {
        near(
            result.component_flux(component),
            expected[component],
            3.0e-12,
            3.0e-12);
        component_sum +=
            result.component_flux(component);
    }
    near(
        component_sum,
        expected_total,
        3.0e-12,
        3.0e-12);
    near(
        result.total_molar_flux_mol_per_s,
        expected_total,
        3.0e-12,
        3.0e-12);
}

void jacobian_fresh_perturbation() {
    const auto owner_chart =
        side_chart(true);
    const auto neighbour_chart =
        side_chart(false);
    const auto base =
        build_component_flux(
            owner_chart,
            owner_chart.base_inputs,
            neighbour_chart,
            neighbour_chart.base_inputs);

    for (std::size_t column = 0U;
         column <
         owner_chart.layout.unknown_count();
         ++column) {
        double step = 1.0e-6;
        if (column ==
            owner_chart.layout.pressure_unknown_index()) {
            step = 10.0;
        } else if (
            column ==
            owner_chart.layout.temperature_unknown_index()) {
            step = 1.0e-4;
        }

        auto plus = owner_chart.base_inputs;
        auto minus = owner_chart.base_inputs;
        plus[column] += step;
        minus[column] -= step;

        const auto plus_flux =
            build_component_flux(
                owner_chart,
                plus,
                neighbour_chart,
                neighbour_chart.base_inputs);
        const auto minus_flux =
            build_component_flux(
                owner_chart,
                minus,
                neighbour_chart,
                neighbour_chart.base_inputs);

        for (std::size_t component = 0U;
             component < 3U;
             ++component) {
            const double fd_value =
                (plus_flux.component_flux(component) -
                 minus_flux.component_flux(component)) /
                (2.0 * step);
            near(
                base.d_component_flux_owner(
                    component,
                    column),
                fd_value,
                3.0e-6,
                3.0e-8);
        }
    }

    for (std::size_t column = 0U;
         column <
         neighbour_chart.layout.unknown_count();
         ++column) {
        double step = 1.0e-6;
        if (column ==
            neighbour_chart.layout.pressure_unknown_index()) {
            step = 10.0;
        } else if (
            column ==
            neighbour_chart.layout.temperature_unknown_index()) {
            step = 1.0e-4;
        }

        auto plus = neighbour_chart.base_inputs;
        auto minus = neighbour_chart.base_inputs;
        plus[column] += step;
        minus[column] -= step;

        const auto plus_flux =
            build_component_flux(
                owner_chart,
                owner_chart.base_inputs,
                neighbour_chart,
                plus);
        const auto minus_flux =
            build_component_flux(
                owner_chart,
                owner_chart.base_inputs,
                neighbour_chart,
                minus);

        for (std::size_t component = 0U;
             component < 3U;
             ++component) {
            const double fd_value =
                (plus_flux.component_flux(component) -
                 minus_flux.component_flux(component)) /
                (2.0 * step);
            near(
                base.d_component_flux_neighbour(
                    component,
                    column),
                fd_value,
                3.0e-6,
                3.0e-8);
        }
    }
}

void exact_zero_tie_uses_owner_content() {
    const auto owner_chart =
        side_chart(true);
    const auto neighbour_chart =
        side_chart(false);
    auto fixture =
        build_fixture(
            owner_chart,
            owner_chart.base_inputs,
            neighbour_chart,
            neighbour_chart.base_inputs);
    const auto result =
        fd::
            build_materialized_tpfa_internal_face_component_molar_flux(
                fixture.phase_flux,
                fixture.owner_state,
                fixture.owner_density,
                fixture.neighbour_state,
                fixture.neighbour_density);

    const auto owner_x =
        fixture.owner_state.phase_composition(
            flow::PhaseSlot3::phase2);
    const auto& tie =
        result.phase_molar_content[2];

    require(
        tie.upwind_selection ==
            flow::UpwindCellSelection3P::
                owner_exact_zero_tie,
        "phase2 is no longer exact-zero owner tie");
    near(
        tie.upwind_molar_density_mol_per_m3,
        fixture.owner_density
            .molar_density_mol_per_m3[2]);
    for (std::size_t component = 0U;
         component < 3U;
         ++component) {
        near(
            tie.upwind_mole_fraction[component],
            owner_x[component]);
    }
    require(
        fixture.phase_flux.phase[2]
                .volumetric_flux_m3_per_s ==
            0.0,
        "exact-zero tie phase flux changed");
}

void invalid_inputs() {
    const auto owner_chart =
        side_chart(true);
    const auto neighbour_chart =
        side_chart(false);
    auto fixture =
        build_fixture(
            owner_chart,
            owner_chart.base_inputs,
            neighbour_chart,
            neighbour_chart.base_inputs);

    auto wrong_density =
        fixture.owner_density;
    wrong_density
        .molar_density_mol_per_m3[0] +=
        100.0;
    expect_invalid(
        [&] {
            (void)fd::
                build_materialized_tpfa_internal_face_component_molar_flux(
                    fixture.phase_flux,
                    fixture.owner_state,
                    wrong_density,
                    fixture.neighbour_state,
                    fixture.neighbour_density);
        },
        "molar-density primal");

    auto wrong_shape =
        fixture.neighbour_density;
    wrong_shape.gradient[1].pop_back();
    expect_invalid(
        [&] {
            (void)fd::
                build_materialized_tpfa_internal_face_component_molar_flux(
                    fixture.phase_flux,
                    fixture.owner_state,
                    fixture.owner_density,
                    fixture.neighbour_state,
                    wrong_shape);
        },
        "gradient shape");

    auto wrong_identity =
        fixture.phase_flux;
    wrong_identity
        .owner_state_identity
        .temperature_k +=
        1.0;
    expect_invalid(
        [&] {
            (void)fd::
                build_materialized_tpfa_internal_face_component_molar_flux(
                    wrong_identity,
                    fixture.owner_state,
                    fixture.owner_density,
                    fixture.neighbour_state,
                    fixture.neighbour_density);
        },
        "state identity");

    auto wrong_selection =
        fixture.phase_flux;
    wrong_selection.phase[0]
        .upwind_selection =
        flow::UpwindCellSelection3P::
            owner_negative_phase_potential;
    expect_invalid(
        [&] {
            (void)fd::
                build_materialized_tpfa_internal_face_component_molar_flux(
                    wrong_selection,
                    fixture.owner_state,
                    fixture.owner_density,
                    fixture.neighbour_state,
                    fixture.neighbour_density);
        },
        "owner-upwind phase");

    auto wrong_flux =
        fixture.phase_flux;
    wrong_flux.phase[1]
        .volumetric_flux_m3_per_s *=
        2.0;
    expect_invalid(
        [&] {
            (void)fd::
                build_materialized_tpfa_internal_face_component_molar_flux(
                    wrong_flux,
                    fixture.owner_state,
                    fixture.owner_density,
                    fixture.neighbour_state,
                    fixture.neighbour_density);
        },
        "inconsistent");
}

void headers() {
    require(
        tpfa_component_molar_flux_header(),
        "component molar-flux header probe failed");
}

using Test =
    std::pair<std::string_view, void (*)()>;

constexpr Test tests[]{
    {"upwind_phase_molar_content", upwind_phase_molar_content},
    {"component_flux_and_closure", component_flux_and_closure},
    {"jacobian_fresh_perturbation", jacobian_fresh_perturbation},
    {"exact_zero_tie_uses_owner_content", exact_zero_tie_uses_owner_content},
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
