#include <mpmc/flow/two_phase_natural_variable.hpp>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <stdexcept>
#include <string>
#include <vector>

namespace flow = mpmc::flow;

namespace {

void require(bool condition, const char* message) {
    if (!condition) {
        throw std::runtime_error(message);
    }
}

void near(double actual, double expected, double tol = 1.0e-12) {
    if (!std::isfinite(actual) ||
        !std::isfinite(expected) ||
        std::abs(actual - expected) >
            tol * std::max(
                1.0,
                std::max(
                    std::abs(actual),
                    std::abs(expected)))) {
        throw std::runtime_error(
            "two-phase numeric mismatch");
    }
}

flow::NaturalVariableCellState2P state() {
    flow::NaturalVariableCellStateInput2P input;
    input.component_ids = {"A", "B", "C"};
    input.reference_pressure_pa = 2.0e6;
    input.temperature_k = 360.0;
    input.independent_saturation = 0.35;
    input.composition_pivot =
        flow::NaturalVariableCompositionPivot2P::
            from_dependent_components(
                3U,
                {1U, 2U});
    input.independent_phase_compositions[0] =
        {0.20, 0.30};
    input.independent_phase_compositions[1] =
        {0.40, 0.35};
    input.phase_properties[0] =
        flow::PhasePropertyPrerequisiteInput{
            10.0,
            800.0,
            2.0e-3,
            1000.0,
            700.0};
    input.phase_properties[1] =
        flow::PhasePropertyPrerequisiteInput{
            5.0,
            120.0,
            1.0e-3,
            1500.0,
            1100.0};
    return flow::NaturalVariableCellState2P::
        create(std::move(input));
}

} // namespace

void two_phase_natural_variable_contract() {
    const auto current = state();
    const auto& layout = current.layout();

    require(
        layout.component_count() == 3U &&
            layout.phase_count() == 2U &&
            layout.unknown_count() == 7U &&
            layout.equation_count() == 7U,
        "two-phase cardinality is not 2*Nc+1");
    require(
        layout.pressure_unknown_index() == 0U &&
            layout.temperature_unknown_index() == 1U &&
            layout.independent_saturation_unknown_index() == 2U &&
            layout.energy_equation_index() == 3U,
        "two-phase primary row/column indexing mismatch");
    require(
        layout.fugacity_equilibrium_equation_index(0U) == 4U &&
            layout.fugacity_equilibrium_equation_index(2U) == 6U,
        "two-phase fugacity equation block mismatch");

    require(
        layout
                .independent_composition_unknown_index(
                    0U,
                    0U)
                .value() ==
            3U &&
        !layout
             .independent_composition_unknown_index(
                 0U,
                 1U)
             .has_value() &&
        layout
                .independent_composition_unknown_index(
                    0U,
                    2U)
                .value() ==
            4U &&
        layout
                .independent_composition_unknown_index(
                    1U,
                    0U)
                .value() ==
            5U &&
        layout
                .independent_composition_unknown_index(
                    1U,
                    1U)
                .value() ==
            6U &&
        !layout
             .independent_composition_unknown_index(
                 1U,
                 2U)
             .has_value(),
        "two-phase composition chart indexing mismatch");

    near(current.phase_saturation(0U), 0.35);
    near(current.phase_saturation(1U), 0.65);

    const auto x0 =
        current.phase_composition(0U);
    const auto x1 =
        current.phase_composition(1U);
    near(x0[0], 0.20);
    near(x0[1], 0.50);
    near(x0[2], 0.30);
    near(x1[0], 0.40);
    near(x1[1], 0.35);
    near(x1[2], 0.25);

    const auto identity =
        flow::two_phase_detail::
            make_state_identity(current);
    require(
        identity.layout.phase_count() == 2U &&
            identity.layout.unknown_count() == 7U &&
            identity.saturation[0] == 0.35 &&
            identity.saturation[1] == 0.65 &&
            identity.saturation[2] == 0.0 &&
            identity.phase_composition[0].size() == 3U &&
            identity.phase_composition[1].size() == 3U &&
            identity.phase_composition[2].empty(),
        "two-phase identity retained invalid inactive phase data");

    const std::size_t q =
        layout.unknown_count();
    std::array<std::vector<double>, 2> dc{
        std::vector<double>(q, 0.0),
        std::vector<double>(q, 0.0)};
    dc[0][0] = 2.0e-6;
    dc[1][0] = 1.0e-6;
    const auto molar =
        flow::make_two_phase_molar_density_linearization(
            current,
            dc);

    std::array<std::vector<double>, 2> zero{
        std::vector<double>(q, 0.0),
        std::vector<double>(q, 0.0)};
    const flow::TransportPropertyProvenance provenance{
        "two-phase-test",
        "controlled",
        "v1"};

    const auto transport =
        flow::make_two_phase_transport_linearization(
            current,
            zero,
            zero,
            {0.7, 0.3},
            zero,
            provenance,
            provenance);

    std::array<std::vector<double>, 2> dh{
        std::vector<double>(q, 0.0),
        std::vector<double>(q, 0.0)};
    std::array<std::vector<double>, 2> du{
        std::vector<double>(q, 0.0),
        std::vector<double>(q, 0.0)};
    dh[0][1] = 3.0;
    dh[1][1] = 4.0;
    du[0][1] = 2.0;
    du[1][1] = 2.5;
    const auto caloric =
        flow::make_two_phase_caloric_linearization(
            current,
            dh,
            du,
            provenance,
            provenance);

    std::vector<double> de_rock(q, 0.0);
    de_rock[1] = 5.0;
    const auto rock =
        flow::make_two_phase_rock_thermal_storage_linearization(
            current,
            2000.0,
            de_rock,
            provenance);

    const auto accumulation =
        flow::build_two_phase_component_accumulation(
            current,
            0.25);
    const auto accumulation_linearization =
        flow::build_two_phase_component_accumulation_linearization(
            current,
            0.25,
            molar);

    double component_sum = 0.0;
    for (double value :
         accumulation
             .component_accumulation_mol_per_bulk_m3) {
        component_sum += value;
    }
    near(
        component_sum,
        accumulation.total_accumulation_mol_per_bulk_m3);

    for (std::size_t column = 0U;
         column < q;
         ++column) {
        double sum = 0.0;
        for (std::size_t component = 0U;
             component < 3U;
             ++component) {
            sum +=
                accumulation_linearization.d_component(
                    component,
                    column);
        }
        near(
            sum,
            accumulation_linearization.d_total(
                column));
    }

    const auto energy =
        flow::build_two_phase_energy_accumulation_snapshot(
            current,
            0.25,
            transport,
            caloric,
            rock);
    const auto energy_linearization =
        flow::build_two_phase_energy_accumulation_linearization(
            current,
            0.25,
            transport,
            caloric,
            rock);
    require(
        energy.state_identity.layout.phase_count() ==
            2U &&
            energy_linearization.input_count == q &&
            energy_linearization
                    .total_internal_energy_gradient
                    .size() ==
                q,
        "two-phase energy reduction shape mismatch");

    std::vector<double> fugacity_value{
        0.1, -0.2, 0.3};
    std::vector<double> fugacity_jacobian(
        3U * q,
        0.0);
    fugacity_jacobian[0U * q + 3U] = 1.0;
    fugacity_jacobian[1U * q + 2U] = 1.0;
    fugacity_jacobian[2U * q + 6U] = 1.0;
    const auto fugacity =
        flow::make_two_phase_fugacity_equilibrium_linearization(
            current,
            fugacity_value,
            fugacity_jacobian);
    require(
        fugacity.layout.phase_count() == 2U &&
            fugacity.component_count() == 3U &&
            fugacity.residual_count() == 3U &&
            fugacity.input_count == 7U,
        "two-phase fugacity carrier shape mismatch");

    const auto pair =
        flow::make_pore_volume_component_accumulation_pair(
            accumulation,
            accumulation);
    const auto component_residual =
        flow::build_backward_euler_component_accumulation_residual(
            pair,
            accumulation_linearization,
            10.0);
    const auto energy_residual =
        flow::build_backward_euler_energy_accumulation_residual(
            energy,
            energy_linearization,
            energy,
            10.0);

    require(
        component_residual.current_layout.phase_count() ==
            2U &&
            component_residual.input_count == q &&
            energy_residual.current_state_identity.layout
                    .phase_count() ==
                2U &&
            energy_residual.input_count == q,
        "two-phase backward-Euler reduction did not preserve 2*Nc+1 layout");
}
