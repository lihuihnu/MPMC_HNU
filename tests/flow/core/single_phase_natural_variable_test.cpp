#include <mpmc/flow/single_phase_natural_variable.hpp>

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
            "single-phase numeric mismatch");
    }
}

flow::NaturalVariableCellState1P state() {
    flow::NaturalVariableCellStateInput1P input;
    input.component_ids = {"A", "B", "C"};
    input.reference_pressure_pa = 2.0e6;
    input.temperature_k = 360.0;
    input.composition_pivot =
        flow::NaturalVariableCompositionPivot1P::
            from_dependent_component(
                3U,
                1U);
    input.independent_composition = {
        0.20,
        0.30};
    input.phase_properties =
        flow::PhasePropertyPrerequisiteInput{
            10.0,
            800.0,
            2.0e-3,
            1000.0,
            700.0};
    return flow::NaturalVariableCellState1P::
        create(std::move(input));
}

} // namespace

void single_phase_natural_variable_contract() {
    const auto current = state();
    const auto& layout = current.layout();

    require(
        layout.component_count() == 3U &&
            layout.phase_count() == 1U &&
            layout.unknown_count() == 4U &&
            layout.equation_count() == 4U,
        "single-phase cardinality is not Nc+1");
    require(
        layout.pressure_unknown_index() == 0U &&
            layout.temperature_unknown_index() == 1U &&
            layout.energy_equation_index() == 3U,
        "single-phase row/column indexing mismatch");
    require(
        !layout.independent_composition_unknown_index(1U)
             .has_value() &&
            layout
                    .independent_composition_unknown_index(
                        0U)
                    .value() ==
                2U &&
            layout
                    .independent_composition_unknown_index(
                        2U)
                    .value() ==
                3U,
        "single-phase pivoted composition indexing mismatch");

    const auto x =
        current.phase_composition();
    require(
        x.size() == 3U,
        "single-phase composition size mismatch");
    near(x[0], 0.20);
    near(x[1], 0.50);
    near(x[2], 0.30);
    near(current.phase_saturation(), 1.0);

    const auto identity =
        flow::single_phase_detail::
            make_state_identity(current);
    require(
        identity.layout.phase_count() == 1U &&
            identity.layout.unknown_count() == 4U &&
            identity.saturation[0] == 1.0 &&
            identity.saturation[1] == 0.0 &&
            identity.saturation[2] == 0.0 &&
            identity.phase_composition[0].size() == 3U &&
            identity.phase_composition[1].empty() &&
            identity.phase_composition[2].empty(),
        "single-phase identity retained stale inactive phase data");

    const std::size_t q =
        layout.unknown_count();
    std::vector<double> dc(q, 0.0);
    dc[0] = 2.0e-6;
    dc[1] = -1.0e-3;
    const auto molar =
        flow::make_single_phase_molar_density_linearization(
            current,
            dc);

    std::vector<double> zero(q, 0.0);
    const flow::TransportPropertyProvenance provenance{
        "single-phase-test",
        "controlled",
        "v1"};

    const auto transport =
        flow::make_single_phase_transport_linearization(
            current,
            zero,
            zero,
            1.0,
            zero,
            provenance,
            provenance);

    std::vector<double> dh(q, 0.0);
    std::vector<double> du(q, 0.0);
    dh[1] = 3.0;
    du[1] = 2.0;
    const auto caloric =
        flow::make_single_phase_caloric_linearization(
            current,
            dh,
            du,
            provenance,
            provenance);

    std::vector<double> de_rock(q, 0.0);
    de_rock[1] = 5.0;
    const auto rock =
        flow::make_single_phase_rock_thermal_storage_linearization(
            current,
            2000.0,
            de_rock,
            provenance);

    const auto accumulation =
        flow::build_single_phase_component_accumulation(
            current,
            0.25);
    const auto accumulation_linearization =
        flow::build_single_phase_component_accumulation_linearization(
            current,
            0.25,
            molar);

    near(
        accumulation.total_accumulation_mol_per_bulk_m3,
        2.5);
    near(
        accumulation.component_accumulation_mol_per_bulk_m3[0],
        0.5);
    near(
        accumulation.component_accumulation_mol_per_bulk_m3[1],
        1.25);
    near(
        accumulation.component_accumulation_mol_per_bulk_m3[2],
        0.75);

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
        flow::build_single_phase_energy_accumulation_snapshot(
            current,
            0.25,
            transport,
            caloric,
            rock);
    const auto energy_linearization =
        flow::build_single_phase_energy_accumulation_linearization(
            current,
            0.25,
            transport,
            caloric,
            rock);
    require(
        energy.state_identity.layout.phase_count() ==
            1U &&
            energy_linearization.input_count == q &&
            energy_linearization
                    .total_internal_energy_gradient
                    .size() ==
                q,
        "single-phase energy reduction shape mismatch");

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
            1U &&
            component_residual.input_count == q &&
            energy_residual.current_state_identity.layout
                    .phase_count() ==
                1U &&
            energy_residual.input_count == q,
        "single-phase backward-Euler reduction did not preserve Nc+1 layout");

    for (double value :
         component_residual
             .component_residual_mol_per_bulk_m3_s) {
        near(value, 0.0);
    }
    near(
        energy_residual.residual_w_per_bulk_m3,
        0.0);
}
