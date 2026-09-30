#ifndef MPMC_FLOW_SINGLE_PHASE_PROPERTIES_HPP
#define MPMC_FLOW_SINGLE_PHASE_PROPERTIES_HPP

#include <mpmc/flow/single_phase_cell_state.hpp>
#include <mpmc/flow/detail/validation.hpp>
#include <mpmc/flow/detail/composition_coordinates.hpp>
#include <mpmc/flow/phase_transport.hpp>

#include <array>
#include <cmath>
#include <cstddef>
#include <limits>
#include <span>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

// Owned property values and gradients on one frozen state; no storage or face assembly.
namespace mpmc::flow {

namespace single_phase_detail {

[[nodiscard]] inline NaturalVariableStateIdentity
make_state_identity(
    const NaturalVariableCellState1P& state) {
    return {
        state.layout().descriptor(),
        std::vector<std::string>{
            state.component_ids().begin(),
            state.component_ids().end()},
        state.reference_pressure_pa(),
        state.temperature_k(),
        {1.0, 0.0, 0.0},
        {
            std::vector<double>{
                state.phase_composition().begin(),
                state.phase_composition().end()},
            std::vector<double>{},
            std::vector<double>{}}};
}

[[nodiscard]] inline bool
same_layout(
    const NaturalVariableLayoutDescriptor& first,
    const NaturalVariableLayoutDescriptor& second) {
    return mpmc::flow::validation_detail::same_layout(first, second);
}

inline void validate_gradient(
    std::span<const double> gradient,
    std::size_t q,
    const char* name) {
    if (gradient.size() != q) {
        throw std::invalid_argument(
            std::string{"mpmc::flow: single-phase "} +
            name +
            " gradient shape mismatch");
    }
    for (double value : gradient) {
        if (!std::isfinite(value)) {
            throw std::invalid_argument(
                std::string{"mpmc::flow: single-phase "} +
                name +
                " gradient contains non-finite derivative");
        }
    }
}

[[nodiscard]] inline double
d_composition(
    const NaturalVariableLayout1P& layout,
    std::size_t component,
    std::size_t column) {
    return composition_coordinate_detail::block_derivative(
        layout.component_count(), layout.dependent_composition_component(),
        2U, component, column);
}

} // namespace single_phase_detail

struct SinglePhaseMolarDensityNaturalVariableLinearization {
    NaturalVariableLayoutDescriptor layout{
        std::size_t{2U},
        std::size_t{1U},
        std::vector<std::size_t>{1U}};
    double molar_density_mol_per_m3{};
    std::vector<double> gradient;
};

struct SinglePhaseTransportNaturalVariableLinearization {
    NaturalVariableStateIdentity state_identity;
    TransportPropertyProvenance mass_density_provenance;
    TransportPropertyProvenance viscosity_provenance;
    double mass_density_kg_per_m3{};
    double dynamic_viscosity_pa_s{};
    double relative_permeability{};
    std::vector<double> mass_density_gradient;
    std::vector<double> dynamic_viscosity_gradient;
    std::vector<double> relative_permeability_gradient;
};

struct SinglePhaseCaloricNaturalVariableLinearization {
    NaturalVariableStateIdentity state_identity;
    TransportPropertyProvenance enthalpy_provenance;
    TransportPropertyProvenance internal_energy_provenance;
    double specific_enthalpy_j_per_kg{};
    double specific_internal_energy_j_per_kg{};
    std::vector<double> specific_enthalpy_gradient;
    std::vector<double> specific_internal_energy_gradient;
};

struct SinglePhaseRockThermalStorageLinearization {
    NaturalVariableStateIdentity state_identity;
    TransportPropertyProvenance provenance;
    double volumetric_internal_energy_j_per_rock_m3{};
    std::vector<double>
        volumetric_internal_energy_gradient;
};


[[nodiscard]] inline
SinglePhaseMolarDensityNaturalVariableLinearization
make_single_phase_molar_density_linearization(
    const NaturalVariableCellState1P& state,
    std::vector<double> gradient) {
    single_phase_detail::validate_gradient(
        gradient,
        state.layout().unknown_count(),
        "molar-density");
    return {
        state.layout().descriptor(),
        state.phase_properties()
            .molar_density_mol_per_m3,
        std::move(gradient)};
}

[[nodiscard]] inline
SinglePhaseTransportNaturalVariableLinearization
make_single_phase_transport_linearization(
    const NaturalVariableCellState1P& state,
    std::vector<double> mass_density_gradient,
    std::vector<double> viscosity_gradient,
    double relative_permeability,
    std::vector<double>
        relative_permeability_gradient,
    TransportPropertyProvenance
        mass_density_provenance,
    TransportPropertyProvenance
        viscosity_provenance) {
    const std::size_t q =
        state.layout().unknown_count();
    single_phase_detail::validate_gradient(
        mass_density_gradient,
        q,
        "mass-density");
    single_phase_detail::validate_gradient(
        viscosity_gradient,
        q,
        "viscosity");
    single_phase_detail::validate_gradient(
        relative_permeability_gradient,
        q,
        "relative-permeability");
    phase_transport_detail::validate_provenance(
        mass_density_provenance,
        "single-phase-mass-density");
    phase_transport_detail::validate_provenance(
        viscosity_provenance,
        "single-phase-viscosity");
    if (!std::isfinite(relative_permeability) ||
        relative_permeability < 0.0) {
        throw std::invalid_argument(
            "mpmc::flow: single-phase relative permeability must be finite and nonnegative");
    }

    const auto& property =
        state.phase_properties();
    return {
        single_phase_detail::make_state_identity(
            state),
        std::move(mass_density_provenance),
        std::move(viscosity_provenance),
        property.mass_density_kg_per_m3,
        property.dynamic_viscosity_pa_s,
        relative_permeability,
        std::move(mass_density_gradient),
        std::move(viscosity_gradient),
        std::move(
            relative_permeability_gradient)};
}

[[nodiscard]] inline
SinglePhaseCaloricNaturalVariableLinearization
make_single_phase_caloric_linearization(
    const NaturalVariableCellState1P& state,
    std::vector<double> enthalpy_gradient,
    std::vector<double> internal_energy_gradient,
    TransportPropertyProvenance
        enthalpy_provenance,
    TransportPropertyProvenance
        internal_energy_provenance) {
    const std::size_t q =
        state.layout().unknown_count();
    single_phase_detail::validate_gradient(
        enthalpy_gradient,
        q,
        "enthalpy");
    single_phase_detail::validate_gradient(
        internal_energy_gradient,
        q,
        "internal-energy");
    phase_transport_detail::validate_provenance(
        enthalpy_provenance,
        "single-phase-enthalpy");
    phase_transport_detail::validate_provenance(
        internal_energy_provenance,
        "single-phase-internal-energy");
    const auto& property =
        state.phase_properties();
    return {
        single_phase_detail::make_state_identity(
            state),
        std::move(enthalpy_provenance),
        std::move(internal_energy_provenance),
        property.specific_enthalpy_j_per_kg,
        property.specific_internal_energy_j_per_kg,
        std::move(enthalpy_gradient),
        std::move(
            internal_energy_gradient)};
}

[[nodiscard]] inline
SinglePhaseRockThermalStorageLinearization
make_single_phase_rock_thermal_storage_linearization(
    const NaturalVariableCellState1P& state,
    double volumetric_internal_energy_j_per_rock_m3,
    std::vector<double> gradient,
    TransportPropertyProvenance provenance) {
    single_phase_detail::validate_gradient(
        gradient,
        state.layout().unknown_count(),
        "rock-internal-energy");
    phase_transport_detail::validate_provenance(
        provenance,
        "single-phase-rock-internal-energy");
    if (!std::isfinite(
            volumetric_internal_energy_j_per_rock_m3)) {
        throw std::invalid_argument(
            "mpmc::flow: single-phase rock internal energy must be finite");
    }
    return {
        single_phase_detail::make_state_identity(
            state),
        std::move(provenance),
        volumetric_internal_energy_j_per_rock_m3,
        std::move(gradient)};
}

} // namespace mpmc::flow

#endif // MPMC_FLOW_SINGLE_PHASE_PROPERTIES_HPP
