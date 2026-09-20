#ifndef MPMC_FLOW_ENERGY_ACCUMULATION_HPP
#define MPMC_FLOW_ENERGY_ACCUMULATION_HPP

#include <mpmc/flow/phase_transport.hpp>

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
    phase_caloric_linearization_convention =
        "flow/phase-caloric-property-linearization/fixed-three-phase/v1";
inline constexpr std::string_view
    stationary_rock_thermal_storage_convention =
        "flow/stationary-rock-volumetric-internal-energy/v1";
inline constexpr std::string_view
    energy_accumulation_convention =
        "flow/pore-volume-energy-accumulation/fixed-three-phase/v1";
inline constexpr std::string_view
    backward_euler_energy_accumulation_convention =
        "flow/energy-accumulation/backward-euler-local/v1";

/// Caloric phase properties and derivatives in one frozen natural-variable chart.
///
/// Primal values are copied from NaturalVariableCellState3P, so this contract
/// cannot silently substitute a second enthalpy/internal-energy state.
struct PhaseCaloricPropertyNaturalVariableLinearization3P {
    static constexpr std::string_view convention =
        phase_caloric_linearization_convention;

    NaturalVariableStateIdentity3P state_identity;
    TransportPropertyProvenance enthalpy_provenance;
    TransportPropertyProvenance internal_energy_provenance;

    std::array<double, 3>
        specific_enthalpy_j_per_kg{};
    std::array<double, 3>
        specific_internal_energy_j_per_kg{};

    std::array<std::vector<double>, 3>
        specific_enthalpy_gradient;
    std::array<std::vector<double>, 3>
        specific_internal_energy_gradient;
};

/// Stationary-rock thermal storage.
///
/// e_r is rock-volume internal energy [J / rock-m^3]. The flow module does not
/// invent rock density, heat capacity or a reference temperature. A caller may
/// derive e_r from such a model and supply its natural-variable gradient.
///
/// "stationary" means the rock contributes accumulation only; it has no
/// advective mass/enthalpy flux.
struct StationaryRockThermalStorageLinearization3P {
    static constexpr std::string_view convention =
        stationary_rock_thermal_storage_convention;

    NaturalVariableStateIdentity3P state_identity;
    TransportPropertyProvenance provenance;
    double volumetric_internal_energy_j_per_rock_m3{};
    std::vector<double>
        volumetric_internal_energy_gradient;
};

struct PoreVolumeEnergyAccumulationSnapshot3P {
    static constexpr std::string_view convention =
        energy_accumulation_convention;

    NaturalVariableStateIdentity3P state_identity;
    double porosity{};
    double fluid_internal_energy_j_per_bulk_m3{};
    double rock_internal_energy_j_per_bulk_m3{};
    double total_internal_energy_j_per_bulk_m3{};
};

struct PoreVolumeEnergyAccumulationLinearization3P {
    static constexpr std::string_view convention =
        energy_accumulation_convention;

    NaturalVariableStateIdentity3P state_identity;
    double porosity{};
    double total_internal_energy_j_per_bulk_m3{};
    std::size_t input_count{};
    std::vector<double>
        total_internal_energy_gradient;
};

/// Backward-Euler energy accumulation residual.
///
/// Units:
///   W / bulk-m^3 = J / (bulk-m^3 s)
///
/// Previous-state storage and dt are frozen. Current porosity is also frozen in
/// this rigid/static porous-medium contract; d(phi)=0 and d(V_b)=0.
struct BackwardEulerEnergyAccumulationResidual3P {
    static constexpr std::string_view convention =
        backward_euler_energy_accumulation_convention;

    NaturalVariableStateIdentity3P
        current_state_identity;
    double porosity{};
    double time_step_seconds{};
    double residual_w_per_bulk_m3{};
    std::size_t input_count{};
    std::vector<double> gradient;

    [[nodiscard]] double d_residual(
        std::size_t column) const {
        return gradient.at(column);
    }
};

namespace energy_accumulation_detail {

[[nodiscard]] inline bool same_state_identity(
    const NaturalVariableStateIdentity3P& first,
    const NaturalVariableStateIdentity3P& second) {
    if (first.component_ids != second.component_ids ||
        first.layout.component_count() !=
            second.layout.component_count() ||
        first.layout.unknown_count() !=
            second.layout.unknown_count() ||
        first.layout.composition_pivot()
                .dependent_components() !=
            second.layout.composition_pivot()
                .dependent_components() ||
        !phase_transport_detail::near_roundoff(
            first.reference_pressure_pa,
            second.reference_pressure_pa) ||
        !phase_transport_detail::near_roundoff(
            first.temperature_k,
            second.temperature_k)) {
        return false;
    }

    for (std::size_t phase = 0U;
         phase < fixed_three_phase_count;
         ++phase) {
        if (!phase_transport_detail::near_roundoff(
                first.saturation[phase],
                second.saturation[phase]) ||
            first.phase_composition[phase].size() !=
                second.phase_composition[phase].size()) {
            return false;
        }
        for (std::size_t component = 0U;
             component <
             first.phase_composition[phase].size();
             ++component) {
            if (!phase_transport_detail::near_roundoff(
                    first.phase_composition[phase][component],
                    second.phase_composition[phase][component])) {
                return false;
            }
        }
    }
    return true;
}

inline void validate_porosity(
    double porosity) {
    if (!std::isfinite(porosity) ||
        !(porosity > 0.0) ||
        !(porosity < 1.0)) {
        throw std::invalid_argument(
            "mpmc::flow: energy accumulation porosity must be finite and strictly between zero and one");
    }
}

inline void validate_snapshot(
    const PoreVolumeEnergyAccumulationSnapshot3P& snapshot) {
    validate_porosity(snapshot.porosity);
    if (!std::isfinite(
            snapshot.fluid_internal_energy_j_per_bulk_m3) ||
        !std::isfinite(
            snapshot.rock_internal_energy_j_per_bulk_m3) ||
        !std::isfinite(
            snapshot.total_internal_energy_j_per_bulk_m3)) {
        throw std::invalid_argument(
            "mpmc::flow: energy accumulation snapshot contains non-finite storage");
    }
    const double sum =
        snapshot.fluid_internal_energy_j_per_bulk_m3 +
        snapshot.rock_internal_energy_j_per_bulk_m3;
    const double scale =
        std::abs(
            snapshot.fluid_internal_energy_j_per_bulk_m3) +
        std::abs(
            snapshot.rock_internal_energy_j_per_bulk_m3) +
        std::abs(
            snapshot.total_internal_energy_j_per_bulk_m3);
    const double reference =
        std::max(
            {1.0, scale});
    if (std::abs(
            sum -
            snapshot.total_internal_energy_j_per_bulk_m3) >
        8192.0 *
            std::numeric_limits<double>::epsilon() *
            reference) {
        throw std::invalid_argument(
            "mpmc::flow: energy accumulation snapshot fluid/rock storage does not close to total");
    }
}


inline void validate_caloric(
    const NaturalVariableCellState3P& state,
    const PhaseCaloricPropertyNaturalVariableLinearization3P&
        caloric) {
    phase_transport_detail::validate_state_identity(
        state,
        caloric.state_identity);
    phase_transport_detail::validate_provenance(
        caloric.enthalpy_provenance,
        "specific-enthalpy");
    phase_transport_detail::validate_provenance(
        caloric.internal_energy_provenance,
        "specific-internal-energy");
    phase_transport_detail::validate_gradient_shape(
        state.layout(),
        caloric.specific_enthalpy_gradient,
        "specific-enthalpy");
    phase_transport_detail::validate_gradient_shape(
        state.layout(),
        caloric.specific_internal_energy_gradient,
        "specific-internal-energy");

    for (std::size_t phase = 0U;
         phase < fixed_three_phase_count;
         ++phase) {
        const auto& properties =
            state.phase_properties(
                static_cast<PhaseSlot3>(
                    phase));
        if (!std::isfinite(
                caloric.specific_enthalpy_j_per_kg[phase]) ||
            !std::isfinite(
                caloric.specific_internal_energy_j_per_kg[phase]) ||
            !phase_transport_detail::near_roundoff(
                caloric.specific_enthalpy_j_per_kg[phase],
                properties.specific_enthalpy_j_per_kg) ||
            !phase_transport_detail::near_roundoff(
                caloric.specific_internal_energy_j_per_kg[phase],
                properties.specific_internal_energy_j_per_kg)) {
            throw std::invalid_argument(
                "mpmc::flow: caloric-property primal does not match exact current state");
        }
    }
}

inline void validate_rock(
    const NaturalVariableCellState3P& state,
    const StationaryRockThermalStorageLinearization3P&
        rock) {
    phase_transport_detail::validate_state_identity(
        state,
        rock.state_identity);
    phase_transport_detail::validate_provenance(
        rock.provenance,
        "rock-volumetric-internal-energy");
    if (!std::isfinite(
            rock.volumetric_internal_energy_j_per_rock_m3) ||
        rock.volumetric_internal_energy_gradient.size() !=
            state.layout().unknown_count()) {
        throw std::invalid_argument(
            "mpmc::flow: malformed stationary-rock thermal storage");
    }
    for (double derivative :
         rock.volumetric_internal_energy_gradient) {
        if (!std::isfinite(derivative)) {
            throw std::invalid_argument(
                "mpmc::flow: stationary-rock thermal storage gradient contains non-finite derivative");
        }
    }
}

[[nodiscard]] inline double saturation_derivative(
    const NaturalVariableLayout3P& layout,
    std::size_t phase,
    std::size_t column) {
    const auto s0 =
        layout.independent_saturation_unknown_index(
            PhaseSlot3::phase0);
    const auto s1 =
        layout.independent_saturation_unknown_index(
            PhaseSlot3::phase1);
    if (!s0.has_value() ||
        !s1.has_value()) {
        throw std::logic_error(
            "mpmc::flow: fixed-three-phase saturation chart is malformed");
    }
    if (phase == 0U) {
        return column == *s0 ? 1.0 : 0.0;
    }
    if (phase == 1U) {
        return column == *s1 ? 1.0 : 0.0;
    }
    if (phase == 2U) {
        return (column == *s0 ||
                column == *s1)
            ? -1.0
            : 0.0;
    }
    throw std::out_of_range(
        "mpmc::flow: phase index out of range in energy accumulation");
}

} // namespace energy_accumulation_detail

[[nodiscard]] inline
PhaseCaloricPropertyNaturalVariableLinearization3P
make_phase_caloric_property_linearization(
    const NaturalVariableCellState3P& state,
    std::array<std::vector<double>, 3>
        specific_enthalpy_gradient,
    std::array<std::vector<double>, 3>
        specific_internal_energy_gradient,
    TransportPropertyProvenance
        enthalpy_provenance,
    TransportPropertyProvenance
        internal_energy_provenance) {
    std::array<double, 3> enthalpy{};
    std::array<double, 3> internal_energy{};
    for (std::size_t phase = 0U;
         phase < fixed_three_phase_count;
         ++phase) {
        const auto& properties =
            state.phase_properties(
                static_cast<PhaseSlot3>(
                    phase));
        enthalpy[phase] =
            properties.specific_enthalpy_j_per_kg;
        internal_energy[phase] =
            properties.specific_internal_energy_j_per_kg;
    }

    PhaseCaloricPropertyNaturalVariableLinearization3P
        result{
            phase_transport_detail::
                make_state_identity(state),
            std::move(enthalpy_provenance),
            std::move(internal_energy_provenance),
            enthalpy,
            internal_energy,
            std::move(specific_enthalpy_gradient),
            std::move(specific_internal_energy_gradient)};

    energy_accumulation_detail::validate_caloric(
        state,
        result);
    return result;
}

[[nodiscard]] inline
StationaryRockThermalStorageLinearization3P
make_stationary_rock_thermal_storage_linearization(
    const NaturalVariableCellState3P& state,
    double volumetric_internal_energy_j_per_rock_m3,
    std::vector<double>
        volumetric_internal_energy_gradient,
    TransportPropertyProvenance provenance) {
    StationaryRockThermalStorageLinearization3P
        result{
            phase_transport_detail::
                make_state_identity(state),
            std::move(provenance),
            volumetric_internal_energy_j_per_rock_m3,
            std::move(
                volumetric_internal_energy_gradient)};
    energy_accumulation_detail::validate_rock(
        state,
        result);
    return result;
}

[[nodiscard]] inline
PoreVolumeEnergyAccumulationSnapshot3P
build_pore_volume_energy_accumulation_snapshot(
    const NaturalVariableCellState3P& state,
    double porosity,
    const PhaseTransportPropertyNaturalVariableLinearization3P&
        transport,
    const PhaseCaloricPropertyNaturalVariableLinearization3P&
        caloric,
    const StationaryRockThermalStorageLinearization3P&
        rock) {
    using namespace energy_accumulation_detail;
    validate_porosity(porosity);
    phase_transport_detail::validate_transport_linearization(
        state,
        transport);
    validate_caloric(
        state,
        caloric);
    validate_rock(
        state,
        rock);

    if (!same_state_identity(
            transport.state_identity,
            caloric.state_identity) ||
        !same_state_identity(
            transport.state_identity,
            rock.state_identity)) {
        throw std::invalid_argument(
            "mpmc::flow: energy accumulation state identities disagree");
    }

    double fluid =
        0.0;
    for (std::size_t phase = 0U;
         phase < fixed_three_phase_count;
         ++phase) {
        const double phase_energy =
            state.phase_saturation(
                static_cast<PhaseSlot3>(
                    phase)) *
            transport.mass_density_kg_per_m3[phase] *
            caloric.specific_internal_energy_j_per_kg[phase];
        if (!std::isfinite(phase_energy)) {
            throw std::range_error(
                "mpmc::flow: fluid internal-energy accumulation is non-finite");
        }
        fluid +=
            porosity *
            phase_energy;
    }

    const double rock_bulk =
        (1.0 - porosity) *
        rock.volumetric_internal_energy_j_per_rock_m3;
    const double total =
        fluid +
        rock_bulk;
    if (!std::isfinite(fluid) ||
        !std::isfinite(rock_bulk) ||
        !std::isfinite(total)) {
        throw std::range_error(
            "mpmc::flow: total energy accumulation is non-finite");
    }

    return {
        phase_transport_detail::
            make_state_identity(state),
        porosity,
        fluid,
        rock_bulk,
        total};
}

[[nodiscard]] inline
PoreVolumeEnergyAccumulationLinearization3P
build_pore_volume_energy_accumulation_linearization(
    const NaturalVariableCellState3P& state,
    double porosity,
    const PhaseTransportPropertyNaturalVariableLinearization3P&
        transport,
    const PhaseCaloricPropertyNaturalVariableLinearization3P&
        caloric,
    const StationaryRockThermalStorageLinearization3P&
        rock) {
    using namespace energy_accumulation_detail;

    const auto primal =
        build_pore_volume_energy_accumulation_snapshot(
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
        double derivative =
            (1.0 - porosity) *
            rock.volumetric_internal_energy_gradient[
                column];

        for (std::size_t phase = 0U;
             phase < fixed_three_phase_count;
             ++phase) {
            const double saturation =
                state.phase_saturation(
                    static_cast<PhaseSlot3>(
                        phase));
            const double ds =
                saturation_derivative(
                    state.layout(),
                    phase,
                    column);
            const double rho =
                transport.mass_density_kg_per_m3[
                    phase];
            const double u =
                caloric.specific_internal_energy_j_per_kg[
                    phase];
            const double drho =
                transport.mass_density_gradient[
                    phase][column];
            const double du =
                caloric.specific_internal_energy_gradient[
                    phase][column];

            derivative +=
                porosity *
                (ds * rho * u +
                 saturation *
                     (drho * u +
                      rho * du));
        }

        if (!std::isfinite(derivative)) {
            throw std::range_error(
                "mpmc::flow: energy accumulation derivative is non-finite");
        }
        gradient[column] =
            derivative;
    }

    return {
        primal.state_identity,
        porosity,
        primal.total_internal_energy_j_per_bulk_m3,
        q,
        std::move(gradient)};
}

[[nodiscard]] inline
BackwardEulerEnergyAccumulationResidual3P
build_backward_euler_energy_accumulation_residual(
    const PoreVolumeEnergyAccumulationSnapshot3P&
        current,
    const PoreVolumeEnergyAccumulationLinearization3P&
        current_linearization,
    const PoreVolumeEnergyAccumulationSnapshot3P&
        previous,
    double time_step_seconds) {
    using namespace energy_accumulation_detail;

    if (!std::isfinite(time_step_seconds) ||
        !(time_step_seconds > 0.0)) {
        throw std::invalid_argument(
            "mpmc::flow: energy backward-Euler time step must be finite and strictly positive [s]");
    }
    validate_snapshot(current);
    validate_snapshot(previous);
    for (double derivative :
         current_linearization.total_internal_energy_gradient) {
        if (!std::isfinite(derivative)) {
            throw std::invalid_argument(
                "mpmc::flow: current energy accumulation linearization contains non-finite derivative");
        }
    }
    if (!same_state_identity(
            current.state_identity,
            current_linearization.state_identity) ||
        current.state_identity.component_ids !=
            previous.state_identity.component_ids ||
        !phase_transport_detail::near_roundoff(
            current.porosity,
            current_linearization.porosity) ||
        !phase_transport_detail::near_roundoff(
            current.porosity,
            previous.porosity) ||
        current_linearization.input_count !=
            current.state_identity.layout.unknown_count() ||
        current_linearization.total_internal_energy_gradient.size() !=
            current_linearization.input_count ||
        !phase_transport_detail::near_roundoff(
            current.total_internal_energy_j_per_bulk_m3,
            current_linearization.total_internal_energy_j_per_bulk_m3)) {
        throw std::invalid_argument(
            "mpmc::flow: energy backward-Euler snapshot/linearization identity mismatch");
    }

    const double residual =
        (current.total_internal_energy_j_per_bulk_m3 -
         previous.total_internal_energy_j_per_bulk_m3) /
        time_step_seconds;
    if (!std::isfinite(residual)) {
        throw std::range_error(
            "mpmc::flow: energy backward-Euler residual is non-finite");
    }

    std::vector<double> gradient =
        current_linearization
            .total_internal_energy_gradient;
    for (double& derivative : gradient) {
        derivative /=
            time_step_seconds;
        if (!std::isfinite(derivative)) {
            throw std::range_error(
                "mpmc::flow: energy backward-Euler derivative is non-finite");
        }
    }

    return {
        current.state_identity,
        current.porosity,
        time_step_seconds,
        residual,
        current_linearization.input_count,
        std::move(gradient)};
}

} // namespace mpmc::flow

#endif // MPMC_FLOW_ENERGY_ACCUMULATION_HPP
