#ifndef MPMC_FLOW_TWO_PHASE_ACCUMULATION_HPP
#define MPMC_FLOW_TWO_PHASE_ACCUMULATION_HPP

#include <mpmc/flow/two_phase_properties.hpp>
#include <mpmc/flow/component_accumulation_time.hpp>
#include <mpmc/flow/energy_accumulation.hpp>

#include <cmath>
#include <cstddef>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

// Component and energy storage, with the existing identity and conservation checks.
namespace mpmc::flow {

[[nodiscard]] inline
PoreVolumeComponentAccumulationSnapshot3P
build_two_phase_component_accumulation(
    const NaturalVariableCellState2P& state,
    double porosity) {
    component_accumulation_detail::
        validate_porosity(porosity);

    const std::size_t n =
        state.layout().component_count();
    PoreVolumeComponentAccumulationSnapshot3P
        result;
    result.porosity = porosity;
    result.component_ids.assign(
        state.component_ids().begin(),
        state.component_ids().end());
    result.component_accumulation_mol_per_bulk_m3
        .assign(n, 0.0);

    for (std::size_t phase = 0U;
         phase < 2U;
         ++phase) {
        const double saturation =
            state.phase_saturation(phase);
        const double density =
            state.phase_properties(phase)
                .molar_density_mol_per_m3;
        const auto composition =
            state.phase_composition(phase);
        result.total_accumulation_mol_per_bulk_m3 +=
            porosity *
            saturation *
            density;
        for (std::size_t component = 0U;
             component < n;
             ++component) {
            result
                .component_accumulation_mol_per_bulk_m3[
                    component] +=
                porosity *
                saturation *
                density *
                composition[component];
        }
    }
    component_accumulation_detail::
        validate_snapshot_identity(result);
    return result;
}

[[nodiscard]] inline
PoreVolumeComponentAccumulationLinearization3P
build_two_phase_component_accumulation_linearization(
    const NaturalVariableCellState2P& state,
    double porosity,
    const TwoPhaseMolarDensityNaturalVariableLinearization&
        density) {
    component_accumulation_detail::
        validate_porosity(porosity);
    if (!two_phase_detail::same_layout(
            state.layout().descriptor(),
            density.layout)) {
        throw std::invalid_argument(
            "mpmc::flow: two-phase molar-density layout mismatch");
    }

    const std::size_t n =
        state.layout().component_count();
    const std::size_t q =
        state.layout().unknown_count();
    for (std::size_t phase = 0U;
         phase < 2U;
         ++phase) {
        two_phase_detail::validate_gradient(
            density.gradient[phase],
            q,
            "molar-density");
        if (!component_accumulation_detail::
                near_roundoff(
                    density
                        .molar_density_mol_per_m3[
                            phase],
                    state.phase_properties(phase)
                        .molar_density_mol_per_m3)) {
            throw std::invalid_argument(
                "mpmc::flow: two-phase molar-density primal mismatch");
        }
    }

    PoreVolumeComponentAccumulationLinearization3P
        result{
            state.layout().descriptor(),
            porosity,
            std::vector<std::string>{
                state.component_ids().begin(),
                state.component_ids().end()},
            q,
            std::vector<double>(
                n * q,
                0.0),
            std::vector<double>(
                q,
                0.0)};

    for (std::size_t column = 0U;
         column < q;
         ++column) {
        for (std::size_t phase = 0U;
             phase < 2U;
             ++phase) {
            const double saturation =
                state.phase_saturation(phase);
            const double ds =
                two_phase_detail::d_saturation(
                    phase,
                    column);
            const double c =
                density
                    .molar_density_mol_per_m3[
                        phase];
            const double dc =
                density.gradient[phase][column];
            const auto x =
                state.phase_composition(phase);

            result.total_accumulation_gradient[
                column] +=
                porosity *
                (ds * c +
                 saturation * dc);

            for (std::size_t component = 0U;
                 component < n;
                 ++component) {
                const double dx =
                    two_phase_detail::
                        d_composition(
                            state.layout(),
                            phase,
                            component,
                            column);
                result.component_jacobian[
                    component * q +
                    column] +=
                    porosity *
                    (ds * c *
                         x[component] +
                     saturation *
                         dc *
                         x[component] +
                     saturation *
                         c *
                         dx);
            }
        }
    }

    component_accumulation_time_detail::
        validate_current_linearization(
            build_two_phase_component_accumulation(
                state,
                porosity),
            result);
    return result;
}

[[nodiscard]] inline
PoreVolumeEnergyAccumulationSnapshot3P
build_two_phase_energy_accumulation_snapshot(
    const NaturalVariableCellState2P& state,
    double porosity,
    const TwoPhaseTransportNaturalVariableLinearization&
        transport,
    const TwoPhaseCaloricNaturalVariableLinearization&
        caloric,
    const TwoPhaseRockThermalStorageLinearization&
        rock) {
    energy_accumulation_detail::
        validate_porosity(porosity);
    const auto identity =
        two_phase_detail::make_state_identity(
            state);
    if (!energy_accumulation_detail::
            same_state_identity(
                identity,
                transport.state_identity) ||
        !energy_accumulation_detail::
            same_state_identity(
                identity,
                caloric.state_identity) ||
        !energy_accumulation_detail::
            same_state_identity(
                identity,
                rock.state_identity)) {
        throw std::invalid_argument(
            "mpmc::flow: two-phase energy property identity mismatch");
    }

    double fluid = 0.0;
    for (std::size_t phase = 0U;
         phase < 2U;
         ++phase) {
        fluid +=
            porosity *
            state.phase_saturation(phase) *
            transport.mass_density_kg_per_m3[
                phase] *
            caloric.specific_internal_energy_j_per_kg[
                phase];
    }
    const double rock_bulk =
        (1.0 - porosity) *
        rock.volumetric_internal_energy_j_per_rock_m3;
    const double total =
        fluid + rock_bulk;
    if (!std::isfinite(fluid) ||
        !std::isfinite(rock_bulk) ||
        !std::isfinite(total)) {
        throw std::range_error(
            "mpmc::flow: two-phase energy accumulation is non-finite");
    }
    return {
        identity,
        porosity,
        fluid,
        rock_bulk,
        total};
}

[[nodiscard]] inline
PoreVolumeEnergyAccumulationLinearization3P
build_two_phase_energy_accumulation_linearization(
    const NaturalVariableCellState2P& state,
    double porosity,
    const TwoPhaseTransportNaturalVariableLinearization&
        transport,
    const TwoPhaseCaloricNaturalVariableLinearization&
        caloric,
    const TwoPhaseRockThermalStorageLinearization&
        rock) {
    const auto primal =
        build_two_phase_energy_accumulation_snapshot(
            state,
            porosity,
            transport,
            caloric,
            rock);
    const std::size_t q =
        state.layout().unknown_count();

    std::vector<double> gradient(
        q,
        0.0);
    for (std::size_t column = 0U;
         column < q;
         ++column) {
        gradient[column] =
            (1.0 - porosity) *
            rock.volumetric_internal_energy_gradient[
                column];

        for (std::size_t phase = 0U;
             phase < 2U;
             ++phase) {
            two_phase_detail::validate_gradient(
                transport
                    .mass_density_gradient[phase],
                q,
                "mass-density");
            two_phase_detail::validate_gradient(
                caloric
                    .specific_internal_energy_gradient[
                        phase],
                q,
                "internal-energy");
            const double saturation =
                state.phase_saturation(phase);
            const double ds =
                two_phase_detail::d_saturation(
                    phase,
                    column);
            const double rho =
                transport.mass_density_kg_per_m3[
                    phase];
            const double drho =
                transport
                    .mass_density_gradient[phase][
                        column];
            const double u =
                caloric
                    .specific_internal_energy_j_per_kg[
                        phase];
            const double du =
                caloric
                    .specific_internal_energy_gradient[
                        phase][column];

            gradient[column] +=
                porosity *
                (ds * rho * u +
                 saturation * drho * u +
                 saturation * rho * du);
        }
        if (!std::isfinite(
                gradient[column])) {
            throw std::range_error(
                "mpmc::flow: two-phase energy accumulation derivative is non-finite");
        }
    }

    return {
        primal.state_identity,
        porosity,
        primal.total_internal_energy_j_per_bulk_m3,
        q,
        std::move(gradient)};
}

} // namespace mpmc::flow

#endif // MPMC_FLOW_TWO_PHASE_ACCUMULATION_HPP
