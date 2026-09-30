#ifndef MPMC_FLOW_TWO_PHASE_PROPERTIES_HPP
#define MPMC_FLOW_TWO_PHASE_PROPERTIES_HPP

#include <mpmc/flow/two_phase_cell_state.hpp>
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

namespace two_phase_detail {

[[nodiscard]] inline NaturalVariableStateIdentity
make_state_identity(
    const NaturalVariableCellState2P& state) {
    return {
        state.layout().descriptor(),
        std::vector<std::string>{
            state.component_ids().begin(),
            state.component_ids().end()},
        state.reference_pressure_pa(),
        state.temperature_k(),
        {
            state.phase_saturation(0U),
            state.phase_saturation(1U),
            0.0},
        {
            std::vector<double>{
                state.phase_composition(0U).begin(),
                state.phase_composition(0U).end()},
            std::vector<double>{
                state.phase_composition(1U).begin(),
                state.phase_composition(1U).end()},
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
            std::string{"mpmc::flow: two-phase "} +
            name +
            " gradient shape mismatch");
    }
    for (double value : gradient) {
        if (!std::isfinite(value)) {
            throw std::invalid_argument(
                std::string{"mpmc::flow: two-phase "} +
                name +
                " gradient contains non-finite derivative");
        }
    }
}

[[nodiscard]] inline double
d_saturation(
    std::size_t phase,
    std::size_t column) {
    if (phase >= 2U) {
        throw std::out_of_range(
            "mpmc::flow: two-phase saturation phase out of range");
    }
    if (column !=
        NaturalVariableLayout2P::
            independent_saturation_unknown_index()) {
        return 0.0;
    }
    return phase == 0U
        ? 1.0
        : -1.0;
}

[[nodiscard]] inline double
d_composition(
    const NaturalVariableLayout2P& layout,
    std::size_t phase,
    std::size_t component,
    std::size_t column) {
    if (phase >= 2U || component >= layout.component_count()) {
        throw std::out_of_range(
            "mpmc::flow: two-phase composition derivative index out of range");
    }
    const std::size_t first_column =
        3U + phase * (layout.component_count() - 1U);
    return composition_coordinate_detail::block_derivative(
        layout.component_count(), layout.dependent_composition_component(phase),
        first_column, component, column);
}

} // namespace two_phase_detail

struct TwoPhaseMolarDensityNaturalVariableLinearization {
    NaturalVariableLayoutDescriptor layout{
        std::size_t{2U},
        std::size_t{2U},
        std::vector<std::size_t>{1U, 1U}};
    std::array<double, 2>
        molar_density_mol_per_m3{};
    std::array<std::vector<double>, 2>
        gradient;
};

struct TwoPhaseTransportNaturalVariableLinearization {
    NaturalVariableStateIdentity state_identity;
    TransportPropertyProvenance mass_density_provenance;
    TransportPropertyProvenance viscosity_provenance;
    std::array<double, 2>
        mass_density_kg_per_m3{};
    std::array<double, 2>
        dynamic_viscosity_pa_s{};
    std::array<double, 2>
        relative_permeability{};
    std::array<std::vector<double>, 2>
        mass_density_gradient;
    std::array<std::vector<double>, 2>
        dynamic_viscosity_gradient;
    std::array<std::vector<double>, 2>
        relative_permeability_gradient;
};

struct TwoPhaseCaloricNaturalVariableLinearization {
    NaturalVariableStateIdentity state_identity;
    TransportPropertyProvenance enthalpy_provenance;
    TransportPropertyProvenance internal_energy_provenance;
    std::array<double, 2>
        specific_enthalpy_j_per_kg{};
    std::array<double, 2>
        specific_internal_energy_j_per_kg{};
    std::array<std::vector<double>, 2>
        specific_enthalpy_gradient;
    std::array<std::vector<double>, 2>
        specific_internal_energy_gradient;
};

struct TwoPhaseRockThermalStorageLinearization {
    NaturalVariableStateIdentity state_identity;
    TransportPropertyProvenance provenance;
    double volumetric_internal_energy_j_per_rock_m3{};
    std::vector<double>
        volumetric_internal_energy_gradient;
};


struct TwoPhaseFugacityEquilibriumLinearization {
    NaturalVariableLayoutDescriptor layout{
        std::size_t{2U},
        std::size_t{2U},
        std::vector<std::size_t>{1U, 1U}};
    std::vector<std::string> component_ids;
    std::vector<double> residual;
    std::size_t input_count{};
    std::vector<double> jacobian;

    [[nodiscard]] std::size_t
    component_count() const noexcept {
        return component_ids.size();
    }

    [[nodiscard]] std::size_t
    residual_count() const noexcept {
        return residual.size();
    }

    [[nodiscard]] double d_residual(
        std::size_t component,
        std::size_t column) const {
        if (component >= component_count() ||
            column >= input_count) {
            throw std::out_of_range(
                "mpmc::flow: two-phase fugacity Jacobian index out of range");
        }
        return jacobian.at(
            component * input_count +
            column);
    }
};

[[nodiscard]] inline
TwoPhaseMolarDensityNaturalVariableLinearization
make_two_phase_molar_density_linearization(
    const NaturalVariableCellState2P& state,
    std::array<std::vector<double>, 2>
        gradient) {
    const std::size_t q =
        state.layout().unknown_count();
    TwoPhaseMolarDensityNaturalVariableLinearization
        result;
    result.layout =
        state.layout().descriptor();
    for (std::size_t phase = 0U;
         phase < 2U;
         ++phase) {
        two_phase_detail::validate_gradient(
            gradient[phase],
            q,
            "molar-density");
        result.molar_density_mol_per_m3[phase] =
            state.phase_properties(phase)
                .molar_density_mol_per_m3;
        result.gradient[phase] =
            std::move(gradient[phase]);
    }
    return result;
}

[[nodiscard]] inline
TwoPhaseTransportNaturalVariableLinearization
make_two_phase_transport_linearization(
    const NaturalVariableCellState2P& state,
    std::array<std::vector<double>, 2>
        mass_density_gradient,
    std::array<std::vector<double>, 2>
        viscosity_gradient,
    std::array<double, 2>
        relative_permeability,
    std::array<std::vector<double>, 2>
        relative_permeability_gradient,
    TransportPropertyProvenance
        mass_density_provenance,
    TransportPropertyProvenance
        viscosity_provenance) {
    phase_transport_detail::validate_provenance(
        mass_density_provenance,
        "two-phase-mass-density");
    phase_transport_detail::validate_provenance(
        viscosity_provenance,
        "two-phase-viscosity");

    TwoPhaseTransportNaturalVariableLinearization
        result;
    result.state_identity =
        two_phase_detail::make_state_identity(
            state);
    result.mass_density_provenance =
        std::move(mass_density_provenance);
    result.viscosity_provenance =
        std::move(viscosity_provenance);

    const std::size_t q =
        state.layout().unknown_count();
    for (std::size_t phase = 0U;
         phase < 2U;
         ++phase) {
        two_phase_detail::validate_gradient(
            mass_density_gradient[phase],
            q,
            "mass-density");
        two_phase_detail::validate_gradient(
            viscosity_gradient[phase],
            q,
            "viscosity");
        two_phase_detail::validate_gradient(
            relative_permeability_gradient[phase],
            q,
            "relative-permeability");
        if (!std::isfinite(
                relative_permeability[phase]) ||
            relative_permeability[phase] < 0.0) {
            throw std::invalid_argument(
                "mpmc::flow: two-phase relative permeability must be finite and nonnegative");
        }
        const auto& property =
            state.phase_properties(phase);
        result.mass_density_kg_per_m3[phase] =
            property.mass_density_kg_per_m3;
        result.dynamic_viscosity_pa_s[phase] =
            property.dynamic_viscosity_pa_s;
        result.relative_permeability[phase] =
            relative_permeability[phase];
        result.mass_density_gradient[phase] =
            std::move(
                mass_density_gradient[phase]);
        result.dynamic_viscosity_gradient[phase] =
            std::move(
                viscosity_gradient[phase]);
        result.relative_permeability_gradient[phase] =
            std::move(
                relative_permeability_gradient[phase]);
    }
    return result;
}

[[nodiscard]] inline
TwoPhaseCaloricNaturalVariableLinearization
make_two_phase_caloric_linearization(
    const NaturalVariableCellState2P& state,
    std::array<std::vector<double>, 2>
        enthalpy_gradient,
    std::array<std::vector<double>, 2>
        internal_energy_gradient,
    TransportPropertyProvenance
        enthalpy_provenance,
    TransportPropertyProvenance
        internal_energy_provenance) {
    phase_transport_detail::validate_provenance(
        enthalpy_provenance,
        "two-phase-enthalpy");
    phase_transport_detail::validate_provenance(
        internal_energy_provenance,
        "two-phase-internal-energy");

    TwoPhaseCaloricNaturalVariableLinearization
        result;
    result.state_identity =
        two_phase_detail::make_state_identity(
            state);
    result.enthalpy_provenance =
        std::move(enthalpy_provenance);
    result.internal_energy_provenance =
        std::move(internal_energy_provenance);

    const std::size_t q =
        state.layout().unknown_count();
    for (std::size_t phase = 0U;
         phase < 2U;
         ++phase) {
        two_phase_detail::validate_gradient(
            enthalpy_gradient[phase],
            q,
            "enthalpy");
        two_phase_detail::validate_gradient(
            internal_energy_gradient[phase],
            q,
            "internal-energy");
        const auto& property =
            state.phase_properties(phase);
        result.specific_enthalpy_j_per_kg[phase] =
            property.specific_enthalpy_j_per_kg;
        result.specific_internal_energy_j_per_kg[phase] =
            property.specific_internal_energy_j_per_kg;
        result.specific_enthalpy_gradient[phase] =
            std::move(
                enthalpy_gradient[phase]);
        result.specific_internal_energy_gradient[phase] =
            std::move(
                internal_energy_gradient[phase]);
    }
    return result;
}

[[nodiscard]] inline
TwoPhaseRockThermalStorageLinearization
make_two_phase_rock_thermal_storage_linearization(
    const NaturalVariableCellState2P& state,
    double volumetric_internal_energy_j_per_rock_m3,
    std::vector<double> gradient,
    TransportPropertyProvenance provenance) {
    two_phase_detail::validate_gradient(
        gradient,
        state.layout().unknown_count(),
        "rock-internal-energy");
    phase_transport_detail::validate_provenance(
        provenance,
        "two-phase-rock-internal-energy");
    if (!std::isfinite(
            volumetric_internal_energy_j_per_rock_m3)) {
        throw std::invalid_argument(
            "mpmc::flow: two-phase rock internal energy must be finite");
    }
    return {
        two_phase_detail::make_state_identity(
            state),
        std::move(provenance),
        volumetric_internal_energy_j_per_rock_m3,
        std::move(gradient)};
}

[[nodiscard]] inline
TwoPhaseFugacityEquilibriumLinearization
make_two_phase_fugacity_equilibrium_linearization(
    const NaturalVariableCellState2P& state,
    std::vector<double> residual,
    std::vector<double> jacobian) {
    const std::size_t n =
        state.layout().component_count();
    const std::size_t q =
        state.layout().unknown_count();
    if (residual.size() != n ||
        n > std::numeric_limits<std::size_t>::max() /
                q ||
        jacobian.size() != n * q) {
        throw std::invalid_argument(
            "mpmc::flow: two-phase fugacity residual/Jacobian shape mismatch");
    }
    for (double value : residual) {
        if (!std::isfinite(value)) {
            throw std::invalid_argument(
                "mpmc::flow: two-phase fugacity residual contains non-finite value");
        }
    }
    for (double value : jacobian) {
        if (!std::isfinite(value)) {
            throw std::invalid_argument(
                "mpmc::flow: two-phase fugacity Jacobian contains non-finite value");
        }
    }
    return {
        state.layout().descriptor(),
        std::vector<std::string>{
            state.component_ids().begin(),
            state.component_ids().end()},
        std::move(residual),
        q,
        std::move(jacobian)};
}

} // namespace mpmc::flow

#endif // MPMC_FLOW_TWO_PHASE_PROPERTIES_HPP
