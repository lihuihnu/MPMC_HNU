#ifndef MPMC_WELL_DISCRETIZATION_ENERGY_RATE_HPP
#define MPMC_WELL_DISCRETIZATION_ENERGY_RATE_HPP

#include <mpmc/flow/energy_accumulation.hpp>
#include <mpmc/well_discretization/pressure_drawdown_rate.hpp>

#include <array>
#include <cmath>
#include <cstddef>
#include <stdexcept>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

namespace mpmc::well_discretization {

inline constexpr std::string_view
    connection_advective_energy_rate_convention =
        "well-discretization/connection-advective-energy-rate/directional-enthalpy/v1";

enum class WellConnectionEnergyEnthalpySource3P {
    reservoir_production_or_zero_tie,
    explicit_well_injection
};

/// Explicit well-side specific enthalpy for injection.
///
/// Values are [J/kg]. They are external well data and therefore have no
/// derivative with respect to the reservoir natural-variable chart in this
/// contract. Enthalpy may be negative relative to a chosen reference, so only
/// finiteness is required.
struct InjectionPhaseSpecificEnthalpy3P {
    std::string provenance;
    std::array<double, 3>
        specific_enthalpy_j_per_kg{};
};

/// Three-phase connection advective energy-rate linearization.
///
/// Phase volumetric rate remains reservoir-volume based and production-positive.
/// Therefore reservoir-side mass density is used for both flow directions.
///
/// For q_alpha >= 0 (including exact zero):
///   E_dot_alpha = q_alpha * rho_alpha,cell * h_alpha,cell.
///
/// For q_alpha < 0:
///   E_dot_alpha = q_alpha * rho_alpha,cell * h_alpha,inj.
///
/// Injection enthalpy is explicit well-side data. It is never substituted by
/// reservoir enthalpy and has zero reservoir-state derivative here.
struct WellConnectionAdvectiveEnergyRateLinearization3P {
    static constexpr std::string_view convention =
        connection_advective_energy_rate_convention;

    mpmc::flow::NaturalVariableStateIdentity3P
        state_identity;
    double bottom_hole_pressure_pa{};

    std::array<double, 3>
        phase_volumetric_rate_m3_per_s{};
    std::array<double, 3>
        phase_mass_density_kg_per_m3{};
    std::array<double, 3>
        selected_specific_enthalpy_j_per_kg{};
    std::array<
        WellConnectionEnergyEnthalpySource3P,
        3>
        enthalpy_source{};
    std::array<double, 3>
        phase_advective_energy_rate_w{};

    double total_advective_energy_rate_w{};
    std::vector<double>
        reservoir_natural_variable_gradient;
    double bottom_hole_pressure_derivative_w_per_pa{};

    [[nodiscard]] double
    d_energy_rate_d_reservoir_unknown(
        std::size_t column) const {
        return reservoir_natural_variable_gradient.at(
            column);
    }

    [[nodiscard]] double
    d_energy_rate_d_bottom_hole_pressure() const noexcept {
        return bottom_hole_pressure_derivative_w_per_pa;
    }
};

namespace energy_rate_detail {

[[nodiscard]] inline bool
same_state_identity(
    const mpmc::flow::NaturalVariableStateIdentity3P&
        first,
    const mpmc::flow::NaturalVariableStateIdentity3P&
        second) {
    return
        first.component_ids ==
            second.component_ids &&
        first.layout.component_count() ==
            second.layout.component_count() &&
        first.layout.phase_count() ==
            second.layout.phase_count() &&
        first.layout.unknown_count() ==
            second.layout.unknown_count() &&
        first.layout.composition_pivot()
                .dependent_components() ==
            second.layout.composition_pivot()
                .dependent_components() &&
        first.reference_pressure_pa ==
            second.reference_pressure_pa &&
        first.temperature_k ==
            second.temperature_k &&
        first.saturation ==
            second.saturation &&
        first.phase_composition ==
            second.phase_composition;
}

inline void
validate_injection_enthalpy(
    const InjectionPhaseSpecificEnthalpy3P&
        injection) {
    if (injection.provenance.empty()) {
        throw std::invalid_argument(
            "mpmc::well_discretization::make_connection_advective_energy_rate_linearization_3p: injection enthalpy provenance must be nonempty");
    }
    for (double value :
         injection.specific_enthalpy_j_per_kg) {
        if (!std::isfinite(value)) {
            throw std::invalid_argument(
                "mpmc::well_discretization::make_connection_advective_energy_rate_linearization_3p: injection enthalpy must be finite [J/kg]");
        }
    }
}

inline void validate_reservoir_energy_properties(
    const mpmc::flow::NaturalVariableStateIdentity3P&
        identity,
    const mpmc::flow::
        PhaseTransportPropertyNaturalVariableLinearization3P&
            transport,
    const mpmc::flow::
        PhaseCaloricPropertyNaturalVariableLinearization3P&
            caloric) {
    if (!same_state_identity(
            identity,
            transport.state_identity) ||
        !same_state_identity(
            identity,
            caloric.state_identity)) {
        throw std::invalid_argument(
            "mpmc::well_discretization::make_connection_advective_energy_rate_linearization_3p: mobility/transport/caloric state identities disagree");
    }

    const std::size_t q =
        identity.layout.unknown_count();
    if (identity.layout.phase_count() !=
            mpmc::flow::fixed_three_phase_count ||
        q == 0U ||
        identity.component_ids.size() !=
            identity.layout.component_count()) {
        throw std::invalid_argument(
            "mpmc::well_discretization: malformed three-phase energy state identity");
    }

    for (std::size_t phase = 0U;
         phase < mpmc::flow::fixed_three_phase_count;
         ++phase) {
        if (!std::isfinite(
                transport.mass_density_kg_per_m3[
                    phase]) ||
            !(transport.mass_density_kg_per_m3[
                  phase] >
              0.0) ||
            transport.mass_density_gradient[
                phase].size() != q ||
            !std::isfinite(
                caloric.specific_enthalpy_j_per_kg[
                    phase]) ||
            caloric.specific_enthalpy_gradient[
                phase].size() != q) {
            throw std::invalid_argument(
                "mpmc::well_discretization: malformed reservoir density/enthalpy energy property");
        }

        for (double derivative :
             transport.mass_density_gradient[
                 phase]) {
            if (!std::isfinite(derivative)) {
                throw std::invalid_argument(
                    "mpmc::well_discretization: reservoir mass-density gradient contains non-finite derivative");
            }
        }
        for (double derivative :
             caloric.specific_enthalpy_gradient[
                 phase]) {
            if (!std::isfinite(derivative)) {
                throw std::invalid_argument(
                    "mpmc::well_discretization: reservoir enthalpy gradient contains non-finite derivative");
            }
        }
    }
}

} // namespace energy_rate_detail

/// Build the three-phase advective energy rate for one well connection.
///
/// This is an energy-flux constitutive contract only. It does not include
/// conductive heat exchange with the wellbore, wellbore heat loss, a thermal
/// well-control equation, source normalization, or PETSc assembly.
[[nodiscard]] inline
WellConnectionAdvectiveEnergyRateLinearization3P
make_connection_advective_energy_rate_linearization_3p(
    const mpmc::well::PeacemanWellIndex3D&
        connection,
    const mpmc::flow::LocalPhaseMobilityLinearization3P&
        mobility,
    const mpmc::flow::
        PhaseTransportPropertyNaturalVariableLinearization3P&
            transport,
    const mpmc::flow::
        PhaseCaloricPropertyNaturalVariableLinearization3P&
            caloric,
    const InjectionPhaseSpecificEnthalpy3P&
        injection_enthalpy,
    double bottom_hole_pressure_pa) {
    using namespace energy_rate_detail;

    validate_injection_enthalpy(
        injection_enthalpy);
    validate_reservoir_energy_properties(
        mobility.state_identity,
        transport,
        caloric);

    const std::size_t q =
        mobility.state_identity.layout
            .unknown_count();

    WellConnectionAdvectiveEnergyRateLinearization3P
        result;
    result.state_identity =
        mobility.state_identity;
    result.bottom_hole_pressure_pa =
        bottom_hole_pressure_pa;
    result.phase_mass_density_kg_per_m3 =
        transport.mass_density_kg_per_m3;
    result.reservoir_natural_variable_gradient.assign(
        q,
        0.0);

    for (std::size_t phase = 0U;
         phase < mpmc::flow::fixed_three_phase_count;
         ++phase) {
        const auto slot =
            static_cast<
                mpmc::flow::PhaseSlot3>(
                    phase);
        const auto phase_rate =
            make_pressure_drawdown_phase_rate_linearization_3p(
                connection,
                mobility,
                slot,
                bottom_hole_pressure_pa);

        const double volumetric_rate =
            phase_rate
                .phase_volumetric_rate_m3_per_s;
        const double density =
            transport.mass_density_kg_per_m3[
                phase];
        const bool production_or_zero =
            volumetric_rate >= 0.0;
        const double enthalpy =
            production_or_zero
                ? caloric
                      .specific_enthalpy_j_per_kg[
                          phase]
                : injection_enthalpy
                      .specific_enthalpy_j_per_kg[
                          phase];

        result.phase_volumetric_rate_m3_per_s[
            phase] =
            volumetric_rate;
        result.selected_specific_enthalpy_j_per_kg[
            phase] =
            enthalpy;
        result.enthalpy_source[phase] =
            production_or_zero
                ? WellConnectionEnergyEnthalpySource3P::
                      reservoir_production_or_zero_tie
                : WellConnectionEnergyEnthalpySource3P::
                      explicit_well_injection;

        const double phase_energy_rate =
            volumetric_rate *
            density *
            enthalpy;
        if (!std::isfinite(phase_energy_rate)) {
            throw std::range_error(
                "mpmc::well_discretization: phase advective energy rate is outside representable range");
        }
        result.phase_advective_energy_rate_w[
            phase] =
            phase_energy_rate;
        result.total_advective_energy_rate_w +=
            phase_energy_rate;

        for (std::size_t column = 0U;
             column < q;
             ++column) {
            const double dq =
                phase_rate
                    .d_phase_rate_d_reservoir_unknown(
                        column);
            const double drho =
                transport.mass_density_gradient[
                    phase][column];
            const double dh =
                production_or_zero
                    ? caloric
                          .specific_enthalpy_gradient[
                              phase][column]
                    : 0.0;

            const double derivative =
                dq * density * enthalpy +
                volumetric_rate *
                    drho *
                    enthalpy +
                volumetric_rate *
                    density *
                    dh;
            if (!std::isfinite(derivative)) {
                throw std::range_error(
                    "mpmc::well_discretization: reservoir advective-energy derivative is outside representable range");
            }
            result.reservoir_natural_variable_gradient[
                column] +=
                derivative;
        }

        const double bhp_derivative =
            phase_rate
                .d_phase_rate_d_bottom_hole_pressure() *
            density *
            enthalpy;
        if (!std::isfinite(bhp_derivative)) {
            throw std::range_error(
                "mpmc::well_discretization: BHP advective-energy derivative is outside representable range");
        }
        result.bottom_hole_pressure_derivative_w_per_pa +=
            bhp_derivative;
    }

    if (!std::isfinite(
            result.total_advective_energy_rate_w) ||
        !std::isfinite(
            result.bottom_hole_pressure_derivative_w_per_pa)) {
        throw std::range_error(
            "mpmc::well_discretization: total advective-energy rate/derivative is outside representable range");
    }

    return result;
}

} // namespace mpmc::well_discretization

#endif // MPMC_WELL_DISCRETIZATION_ENERGY_RATE_HPP
