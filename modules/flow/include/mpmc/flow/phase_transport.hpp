#ifndef MPMC_FLOW_PHASE_TRANSPORT_HPP
#define MPMC_FLOW_PHASE_TRANSPORT_HPP

#include <mpmc/flow/natural_variable_cell_state.hpp>
#include <mpmc/flow/saturation_constitutive.hpp>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <limits>
#include <span>
#include <stdexcept>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

namespace mpmc::flow {

inline constexpr std::string_view
    phase_transport_linearization_convention =
        "flow/phase-transport-property-linearization/fixed-three-phase/v1";

inline constexpr std::string_view
    local_phase_mobility_convention =
        "flow/local-phase-mobility/fixed-three-phase/v1";

struct TransportPropertyProvenance {
    std::string model;
    std::string dataset_id;
    std::string revision;
};

/// Owned exact-state identity used to prevent mixing derivatives from another
/// natural-variable state with the current local constitutive evaluation.
struct NaturalVariableStateIdentity {
    NaturalVariableLayoutDescriptor layout{
        std::size_t{2U},
        std::size_t{1U},
        std::vector<std::size_t>{1U}};
    std::vector<std::string> component_ids;
    double reference_pressure_pa{};
    double temperature_k{};
    std::vector<double> saturation;
    std::vector<std::vector<double>>
        phase_composition;
};

using NaturalVariableStateIdentity3P =
    NaturalVariableStateIdentity;

struct PhaseTransportPropertyNaturalVariableLinearization3P {
    static constexpr std::string_view convention =
        phase_transport_linearization_convention;

    NaturalVariableStateIdentity3P state_identity;
    TransportPropertyProvenance mass_density_provenance;
    TransportPropertyProvenance viscosity_provenance;

    std::array<double, 3>
        mass_density_kg_per_m3{};
    std::array<double, 3>
        dynamic_viscosity_pa_s{};

    std::array<std::vector<double>, 3>
        mass_density_gradient;
    std::array<std::vector<double>, 3>
        dynamic_viscosity_gradient;
};

/// Derivatives of the saturation constitutive laws in their native two-column
/// chart (S0,S1).
///
/// phase0 capillary offset is identically zero by convention and therefore its
/// derivative row must be exactly zero.
struct ThreePhaseSaturationCoordinateDerivatives3P {
    std::array<std::array<double, 2>, 3>
        relative_permeability;
    std::array<std::array<double, 2>, 3>
        capillary_pressure_offset_pa;
};

/// Saturation constitutive primal and Jacobian promoted from the native
/// (S0,S1) chart to the full current NaturalVariableLayout3P.
///
/// kr depends only on saturation in this v1 contract.
/// p_alpha has dp_alpha/dp_ref = 1 and saturation derivatives inherited from
/// the capillary-pressure law. T/composition columns are zero.
struct ThreePhaseSaturationConstitutiveNaturalVariableLinearization3P {
    NaturalVariableStateIdentity3P state_identity;

    std::array<double, 3>
        relative_permeability{};
    std::array<double, 3>
        capillary_pressure_offset_pa{};
    std::array<double, 3>
        phase_pressure_pa{};

    std::array<std::vector<double>, 3>
        relative_permeability_gradient;
    std::array<std::vector<double>, 3>
        phase_pressure_gradient;
};

/// Complete local Darcy-prerequisite property package, still without K,
/// transmissibility, gravity force or a face flux.
///
/// mobility = kr / mu has units [1 / (Pa s)].
struct LocalPhaseMobilityLinearization3P {
    static constexpr std::string_view convention =
        local_phase_mobility_convention;

    NaturalVariableStateIdentity3P state_identity;
    TransportPropertyProvenance mass_density_provenance;
    TransportPropertyProvenance viscosity_provenance;

    std::array<double, 3> phase_pressure_pa{};
    std::array<double, 3> mass_density_kg_per_m3{};
    std::array<double, 3> dynamic_viscosity_pa_s{};
    std::array<double, 3> relative_permeability{};
    std::array<double, 3> mobility_per_pa_s{};

    std::array<std::vector<double>, 3>
        phase_pressure_gradient;
    std::array<std::vector<double>, 3>
        mass_density_gradient;
    std::array<std::vector<double>, 3>
        dynamic_viscosity_gradient;
    std::array<std::vector<double>, 3>
        relative_permeability_gradient;
    std::array<std::vector<double>, 3>
        mobility_gradient;
};

namespace phase_transport_detail {

[[nodiscard]] inline bool near_roundoff(
    double first,
    double second) {
    if (!std::isfinite(first) ||
        !std::isfinite(second)) {
        return false;
    }
    const double scale =
        std::max(
            {1.0,
             std::abs(first),
             std::abs(second)});
    return std::abs(first - second) <=
        4096.0 *
            std::numeric_limits<double>::epsilon() *
            scale;
}

inline void validate_provenance(
    const TransportPropertyProvenance& provenance,
    const char* property_name) {
    if (provenance.model.empty() ||
        provenance.dataset_id.empty() ||
        provenance.revision.empty()) {
        throw std::invalid_argument(
            std::string{"mpmc::flow: "} +
            property_name +
            " provenance requires nonempty model/dataset/revision");
    }
}

[[nodiscard]] inline NaturalVariableStateIdentity3P
make_state_identity(
    const NaturalVariableCellState3P& state) {
    NaturalVariableStateIdentity3P result{
        NaturalVariableLayoutDescriptor{
            state.layout()},
        std::vector<std::string>{
            state.component_ids().begin(),
            state.component_ids().end()},
        state.reference_pressure_pa(),
        state.temperature_k(),
        std::vector<double>(
            fixed_three_phase_count,
            0.0),
        std::vector<std::vector<double>>(
            fixed_three_phase_count)};

    for (std::size_t phase = 0U;
         phase < fixed_three_phase_count;
         ++phase) {
        const auto slot =
            static_cast<PhaseSlot3>(phase);
        result.saturation[phase] =
            state.phase_saturation(slot);
        const auto composition =
            state.phase_composition(slot);
        result.phase_composition[phase].assign(
            composition.begin(),
            composition.end());
    }
    return result;
}

inline void validate_state_identity(
    const NaturalVariableCellState3P& state,
    const NaturalVariableStateIdentity3P& identity) {
    const auto& layout =
        state.layout();
    if (identity.component_ids !=
            std::vector<std::string>{
                state.component_ids().begin(),
                state.component_ids().end()} ||
        identity.layout.component_count() !=
            layout.component_count() ||
        identity.layout.phase_count() !=
            fixed_three_phase_count ||
        identity.layout.composition_pivot()
                .dependent_components() !=
            std::vector<std::size_t>{
                layout.composition_pivot()
                    .dependent_components()
                    .begin(),
                layout.composition_pivot()
                    .dependent_components()
                    .end()} ||
        identity.layout.unknown_count() !=
            layout.unknown_count() ||
        identity.saturation.size() !=
            fixed_three_phase_count ||
        identity.phase_composition.size() !=
            fixed_three_phase_count ||
        !near_roundoff(
            identity.reference_pressure_pa,
            state.reference_pressure_pa()) ||
        !near_roundoff(
            identity.temperature_k,
            state.temperature_k())) {
        throw std::invalid_argument(
            "mpmc::flow: phase transport linearization state identity/chart mismatch");
    }

    for (std::size_t phase = 0U;
         phase < fixed_three_phase_count;
         ++phase) {
        const auto slot =
            static_cast<PhaseSlot3>(phase);
        if (!near_roundoff(
                identity.saturation[phase],
                state.phase_saturation(slot))) {
            throw std::invalid_argument(
                "mpmc::flow: phase transport linearization saturation state mismatch");
        }
        const auto composition =
            state.phase_composition(slot);
        if (identity.phase_composition[phase].size() !=
            composition.size()) {
            throw std::invalid_argument(
                "mpmc::flow: phase transport linearization composition shape mismatch");
        }
        for (std::size_t component = 0U;
             component < composition.size();
             ++component) {
            if (!near_roundoff(
                    identity.phase_composition[phase][component],
                    composition[component])) {
                throw std::invalid_argument(
                    "mpmc::flow: phase transport linearization composition state mismatch");
            }
        }
    }
}

inline void validate_gradient_shape(
    const NaturalVariableLayout3P& layout,
    const std::array<std::vector<double>, 3>& gradient,
    const char* property_name) {
    const std::size_t q =
        layout.unknown_count();
    for (const auto& phase_gradient :
         gradient) {
        if (phase_gradient.size() != q) {
            throw std::invalid_argument(
                std::string{"mpmc::flow: "} +
                property_name +
                " gradient shape does not match natural-variable chart");
        }
        for (double value :
             phase_gradient) {
            if (!std::isfinite(value)) {
                throw std::invalid_argument(
                    std::string{"mpmc::flow: "} +
                    property_name +
                    " gradient contains non-finite derivative");
            }
        }
    }
}

inline void validate_transport_linearization(
    const NaturalVariableCellState3P& state,
    const PhaseTransportPropertyNaturalVariableLinearization3P&
        transport) {
    validate_state_identity(
        state,
        transport.state_identity);
    validate_provenance(
        transport.mass_density_provenance,
        "mass-density");
    validate_provenance(
        transport.viscosity_provenance,
        "viscosity");

    for (std::size_t phase = 0U;
         phase < fixed_three_phase_count;
         ++phase) {
        const auto slot =
            static_cast<PhaseSlot3>(phase);
        const auto& properties =
            state.phase_properties(slot);

        if (!std::isfinite(
                transport.mass_density_kg_per_m3[phase]) ||
            !(transport.mass_density_kg_per_m3[phase] > 0.0) ||
            !near_roundoff(
                transport.mass_density_kg_per_m3[phase],
                properties.mass_density_kg_per_m3)) {
            throw std::invalid_argument(
                "mpmc::flow: phase mass-density primal does not match exact current state");
        }

        if (!std::isfinite(
                transport.dynamic_viscosity_pa_s[phase]) ||
            !(transport.dynamic_viscosity_pa_s[phase] > 0.0) ||
            !near_roundoff(
                transport.dynamic_viscosity_pa_s[phase],
                properties.dynamic_viscosity_pa_s)) {
            throw std::invalid_argument(
                "mpmc::flow: phase viscosity primal does not match exact current state");
        }
    }

    validate_gradient_shape(
        state.layout(),
        transport.mass_density_gradient,
        "mass-density");
    validate_gradient_shape(
        state.layout(),
        transport.dynamic_viscosity_gradient,
        "viscosity");
}

inline void validate_saturation_linearization(
    const NaturalVariableCellState3P& state,
    const ThreePhaseSaturationConstitutiveNaturalVariableLinearization3P&
        constitutive) {
    validate_state_identity(
        state,
        constitutive.state_identity);
    validate_gradient_shape(
        state.layout(),
        constitutive.relative_permeability_gradient,
        "relative-permeability");
    validate_gradient_shape(
        state.layout(),
        constitutive.phase_pressure_gradient,
        "phase-pressure");

    for (std::size_t phase = 0U;
         phase < fixed_three_phase_count;
         ++phase) {
        if (!std::isfinite(
                constitutive.relative_permeability[phase]) ||
            constitutive.relative_permeability[phase] < 0.0 ||
            !std::isfinite(
                constitutive.phase_pressure_pa[phase]) ||
            !(constitutive.phase_pressure_pa[phase] > 0.0) ||
            !std::isfinite(
                constitutive.capillary_pressure_offset_pa[phase])) {
            throw std::invalid_argument(
                "mpmc::flow: invalid saturation-constitutive primal in transport linearization");
        }
    }

    if (constitutive.capillary_pressure_offset_pa[0] != 0.0) {
        throw std::invalid_argument(
            "mpmc::flow: reference-phase capillary offset must be exactly zero");
    }

    for (std::size_t phase = 0U;
         phase < fixed_three_phase_count;
         ++phase) {
        const double expected_pressure =
            state.reference_pressure_pa() +
            constitutive.capillary_pressure_offset_pa[
                phase];
        if (!near_roundoff(
                constitutive.phase_pressure_pa[phase],
                expected_pressure)) {
            throw std::invalid_argument(
                "mpmc::flow: resolved phase pressure is inconsistent with reference pressure plus capillary offset");
        }
    }
}

} // namespace phase_transport_detail

[[nodiscard]] inline
PhaseTransportPropertyNaturalVariableLinearization3P
make_phase_transport_property_linearization(
    const NaturalVariableCellState3P& state,
    std::array<std::vector<double>, 3>
        mass_density_gradient,
    std::array<std::vector<double>, 3>
        dynamic_viscosity_gradient,
    TransportPropertyProvenance
        mass_density_provenance,
    TransportPropertyProvenance
        viscosity_provenance) {
    std::array<double, 3>
        mass_density{};
    std::array<double, 3>
        viscosity{};

    for (std::size_t phase = 0U;
         phase < fixed_three_phase_count;
         ++phase) {
        const auto& properties =
            state.phase_properties(
                static_cast<PhaseSlot3>(phase));
        mass_density[phase] =
            properties.mass_density_kg_per_m3;
        viscosity[phase] =
            properties.dynamic_viscosity_pa_s;
    }

    PhaseTransportPropertyNaturalVariableLinearization3P
        result{
            phase_transport_detail::
                make_state_identity(state),
            std::move(mass_density_provenance),
            std::move(viscosity_provenance),
            mass_density,
            viscosity,
            std::move(mass_density_gradient),
            std::move(dynamic_viscosity_gradient)};

    phase_transport_detail::
        validate_transport_linearization(
            state,
            result);
    return result;
}

/// Promote a scalar-generic saturation-law derivative payload from native
/// (S0,S1) coordinates into the full current natural-variable chart.
[[nodiscard]] inline
ThreePhaseSaturationConstitutiveNaturalVariableLinearization3P
make_saturation_constitutive_natural_variable_linearization(
    const NaturalVariableCellState3P& state,
    const ThreePhaseSaturationConstitutiveEvaluation3P<double>&
        primal,
    const ThreePhaseSaturationCoordinateDerivatives3P&
        derivatives) {
    const auto& layout =
        state.layout();
    const std::size_t q =
        layout.unknown_count();

    for (std::size_t phase = 0U;
         phase < fixed_three_phase_count;
         ++phase) {
        const auto slot =
            static_cast<PhaseSlot3>(phase);
        if (!phase_transport_detail::near_roundoff(
                primal.saturation_state.saturation[phase],
                state.phase_saturation(slot))) {
            throw std::invalid_argument(
                "mpmc::flow: saturation constitutive primal does not match exact current saturation");
        }
        for (std::size_t direction = 0U;
             direction < 2U;
             ++direction) {
            if (!std::isfinite(
                    derivatives.relative_permeability[phase][direction]) ||
                !std::isfinite(
                    derivatives.capillary_pressure_offset_pa[phase][direction])) {
                throw std::invalid_argument(
                    "mpmc::flow: saturation constitutive coordinate derivative is non-finite");
            }
        }
    }

    if (!phase_transport_detail::near_roundoff(
            primal.phase_pressure_pa[0],
            state.reference_pressure_pa()) ||
        primal.capillary_pressure_offset_pa[0] != 0.0 ||
        derivatives.capillary_pressure_offset_pa[0][0] != 0.0 ||
        derivatives.capillary_pressure_offset_pa[0][1] != 0.0) {
        throw std::invalid_argument(
            "mpmc::flow: reference-phase pressure/capillary convention mismatch");
    }

    std::array<std::vector<double>, 3>
        kr_gradient;
    std::array<std::vector<double>, 3>
        pressure_gradient;

    const auto s0_column =
        layout.independent_saturation_unknown_index(
            PhaseSlot3::phase0);
    const auto s1_column =
        layout.independent_saturation_unknown_index(
            PhaseSlot3::phase1);
    if (!s0_column || !s1_column) {
        throw std::logic_error(
            "mpmc::flow: fixed-three-phase saturation chart is malformed");
    }

    for (std::size_t phase = 0U;
         phase < fixed_three_phase_count;
         ++phase) {
        kr_gradient[phase].assign(q, 0.0);
        pressure_gradient[phase].assign(q, 0.0);

        kr_gradient[phase][*s0_column] =
            derivatives.relative_permeability[phase][0];
        kr_gradient[phase][*s1_column] =
            derivatives.relative_permeability[phase][1];

        pressure_gradient[phase][
            layout.pressure_unknown_index()] =
            1.0;
        pressure_gradient[phase][*s0_column] =
            derivatives.capillary_pressure_offset_pa[phase][0];
        pressure_gradient[phase][*s1_column] =
            derivatives.capillary_pressure_offset_pa[phase][1];
    }

    ThreePhaseSaturationConstitutiveNaturalVariableLinearization3P
        result{
            phase_transport_detail::
                make_state_identity(state),
            primal.relative_permeability,
            primal.capillary_pressure_offset_pa,
            primal.phase_pressure_pa,
            std::move(kr_gradient),
            std::move(pressure_gradient)};

    phase_transport_detail::
        validate_saturation_linearization(
            state,
            result);
    return result;
}

[[nodiscard]] inline
LocalPhaseMobilityLinearization3P
build_local_phase_mobility_linearization(
    const NaturalVariableCellState3P& state,
    const PhaseTransportPropertyNaturalVariableLinearization3P&
        transport,
    const ThreePhaseSaturationConstitutiveNaturalVariableLinearization3P&
        constitutive) {
    phase_transport_detail::
        validate_transport_linearization(
            state,
            transport);
    phase_transport_detail::
        validate_saturation_linearization(
            state,
            constitutive);

    const std::size_t q =
        state.layout().unknown_count();

    LocalPhaseMobilityLinearization3P
        result{
            phase_transport_detail::
                make_state_identity(state),
            transport.mass_density_provenance,
            transport.viscosity_provenance,
            constitutive.phase_pressure_pa,
            transport.mass_density_kg_per_m3,
            transport.dynamic_viscosity_pa_s,
            constitutive.relative_permeability,
            {},
            constitutive.phase_pressure_gradient,
            transport.mass_density_gradient,
            transport.dynamic_viscosity_gradient,
            constitutive.relative_permeability_gradient,
            {}};

    for (std::size_t phase = 0U;
         phase < fixed_three_phase_count;
         ++phase) {
        const double kr =
            result.relative_permeability[phase];
        const double viscosity =
            result.dynamic_viscosity_pa_s[phase];
        const double mobility =
            kr / viscosity;
        if (!std::isfinite(mobility) ||
            mobility < 0.0) {
            throw std::range_error(
                "mpmc::flow: local phase mobility is non-finite or negative");
        }
        result.mobility_per_pa_s[phase] =
            mobility;
        result.mobility_gradient[phase].assign(
            q,
            0.0);

        const double viscosity_squared =
            viscosity * viscosity;
        if (!std::isfinite(viscosity_squared) ||
            !(viscosity_squared > 0.0)) {
            throw std::range_error(
                "mpmc::flow: phase viscosity square is invalid");
        }

        for (std::size_t column = 0U;
             column < q;
             ++column) {
            const double derivative =
                (result.relative_permeability_gradient[phase][column] *
                     viscosity -
                 kr *
                     result.dynamic_viscosity_gradient[phase][column]) /
                viscosity_squared;
            if (!std::isfinite(derivative)) {
                throw std::range_error(
                    "mpmc::flow: local phase mobility derivative is non-finite");
            }
            result.mobility_gradient[phase][column] =
                derivative;
        }
    }

    return result;
}

} // namespace mpmc::flow

#endif // MPMC_FLOW_PHASE_TRANSPORT_HPP
