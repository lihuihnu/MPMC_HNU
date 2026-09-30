#ifndef MPMC_FLOW_SINGLE_PHASE_ACCUMULATION_HPP
#define MPMC_FLOW_SINGLE_PHASE_ACCUMULATION_HPP

#include <mpmc/flow/single_phase_properties.hpp>
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
build_single_phase_component_accumulation(
    const NaturalVariableCellState1P& state,
    double porosity) {
    component_accumulation_detail::
        validate_porosity(porosity);

    PoreVolumeComponentAccumulationSnapshot3P
        result;
    result.porosity = porosity;
    result.component_ids.assign(
        state.component_ids().begin(),
        state.component_ids().end());
    result.component_accumulation_mol_per_bulk_m3
        .assign(
            result.component_ids.size(),
            0.0);

    const double density =
        state.phase_properties()
            .molar_density_mol_per_m3;
    const auto composition =
        state.phase_composition();
    for (std::size_t component = 0U;
         component < composition.size();
         ++component) {
        result
            .component_accumulation_mol_per_bulk_m3[
                component] =
            porosity *
            density *
            composition[component];
    }
    result.total_accumulation_mol_per_bulk_m3 =
        porosity * density;
    component_accumulation_detail::
        validate_snapshot_identity(result);
    return result;
}

[[nodiscard]] inline
PoreVolumeComponentAccumulationLinearization3P
build_single_phase_component_accumulation_linearization(
    const NaturalVariableCellState1P& state,
    double porosity,
    const SinglePhaseMolarDensityNaturalVariableLinearization&
        density) {
    component_accumulation_detail::
        validate_porosity(porosity);
    const auto descriptor =
        state.layout().descriptor();
    if (!single_phase_detail::same_layout(
            descriptor,
            density.layout) ||
        !component_accumulation_detail::
            near_roundoff(
                density.molar_density_mol_per_m3,
                state.phase_properties()
                    .molar_density_mol_per_m3)) {
        throw std::invalid_argument(
            "mpmc::flow: single-phase molar-density linearization identity mismatch");
    }
    const std::size_t n =
        state.layout().component_count();
    const std::size_t q =
        state.layout().unknown_count();
    single_phase_detail::validate_gradient(
        density.gradient,
        q,
        "molar-density");

    PoreVolumeComponentAccumulationLinearization3P
        result{
            descriptor,
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

    const double c =
        density.molar_density_mol_per_m3;
    const auto x =
        state.phase_composition();
    for (std::size_t column = 0U;
         column < q;
         ++column) {
        result.total_accumulation_gradient[
            column] =
            porosity *
            density.gradient[column];
        for (std::size_t component = 0U;
             component < n;
             ++component) {
            const double dx =
                single_phase_detail::
                    d_composition(
                        state.layout(),
                        component,
                        column);
            result.component_jacobian[
                component * q + column] =
                porosity *
                (density.gradient[column] *
                     x[component] +
                 c * dx);
        }
    }
    component_accumulation_time_detail::
        validate_current_linearization(
            build_single_phase_component_accumulation(
                state,
                porosity),
            result);
    return result;
}

[[nodiscard]] inline
PoreVolumeEnergyAccumulationSnapshot3P
build_single_phase_energy_accumulation_snapshot(
    const NaturalVariableCellState1P& state,
    double porosity,
    const SinglePhaseTransportNaturalVariableLinearization&
        transport,
    const SinglePhaseCaloricNaturalVariableLinearization&
        caloric,
    const SinglePhaseRockThermalStorageLinearization&
        rock) {
    energy_accumulation_detail::
        validate_porosity(porosity);
    const auto identity =
        single_phase_detail::
            make_state_identity(state);
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
            "mpmc::flow: single-phase energy property identity mismatch");
    }
    const double fluid =
        porosity *
        transport.mass_density_kg_per_m3 *
        caloric.specific_internal_energy_j_per_kg;
    const double rock_bulk =
        (1.0 - porosity) *
        rock.volumetric_internal_energy_j_per_rock_m3;
    const double total =
        fluid + rock_bulk;
    if (!std::isfinite(fluid) ||
        !std::isfinite(rock_bulk) ||
        !std::isfinite(total)) {
        throw std::range_error(
            "mpmc::flow: single-phase energy accumulation is non-finite");
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
build_single_phase_energy_accumulation_linearization(
    const NaturalVariableCellState1P& state,
    double porosity,
    const SinglePhaseTransportNaturalVariableLinearization&
        transport,
    const SinglePhaseCaloricNaturalVariableLinearization&
        caloric,
    const SinglePhaseRockThermalStorageLinearization&
        rock) {
    const auto primal =
        build_single_phase_energy_accumulation_snapshot(
            state,
            porosity,
            transport,
            caloric,
            rock);
    const std::size_t q =
        state.layout().unknown_count();
    single_phase_detail::validate_gradient(
        transport.mass_density_gradient,
        q,
        "mass-density");
    single_phase_detail::validate_gradient(
        caloric.specific_internal_energy_gradient,
        q,
        "internal-energy");
    single_phase_detail::validate_gradient(
        rock.volumetric_internal_energy_gradient,
        q,
        "rock-internal-energy");

    std::vector<double> gradient(
        q,
        0.0);
    const double rho =
        transport.mass_density_kg_per_m3;
    const double u =
        caloric.specific_internal_energy_j_per_kg;
    for (std::size_t column = 0U;
         column < q;
         ++column) {
        gradient[column] =
            porosity *
                (transport
                     .mass_density_gradient[column] *
                     u +
                 rho *
                     caloric
                         .specific_internal_energy_gradient[
                             column]) +
            (1.0 - porosity) *
                rock.volumetric_internal_energy_gradient[
                    column];
        if (!std::isfinite(
                gradient[column])) {
            throw std::range_error(
                "mpmc::flow: single-phase energy accumulation derivative is non-finite");
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

#endif // MPMC_FLOW_SINGLE_PHASE_ACCUMULATION_HPP
