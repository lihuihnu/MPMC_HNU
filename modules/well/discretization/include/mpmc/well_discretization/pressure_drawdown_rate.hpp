#ifndef MPMC_WELL_DISCRETIZATION_PRESSURE_DRAWDOWN_RATE_HPP
#define MPMC_WELL_DISCRETIZATION_PRESSURE_DRAWDOWN_RATE_HPP

#include <mpmc/well_discretization/hydraulic_conductance.hpp>

#include <cmath>
#include <cstddef>
#include <stdexcept>
#include <string_view>
#include <utility>
#include <vector>

namespace mpmc::well_discretization {

inline constexpr std::string_view
    pressure_drawdown_phase_rate_convention =
        "well-discretization/pressure-drawdown-phase-rate/cell-minus-bhp/v1";

/// One connection's one-phase reservoir-volume rate linearization.
///
/// Sign convention:
///   q_alpha > 0  : reservoir cell -> well (production)
///   q_alpha < 0  : well -> reservoir cell (injection)
///
/// Contract:
///   Delta p_alpha = p_alpha,cell - p_bhp                 [Pa]
///   q_alpha       = C_alpha * Delta p_alpha              [m3/s]
///
/// Reservoir natural-variable derivative:
///   dq_alpha/dq_j =
///       dC_alpha/dq_j * Delta p_alpha
///       + C_alpha * dp_alpha,cell/dq_j
///
/// Explicit bottom-hole-pressure derivative:
///   dq_alpha/dp_bhp = -C_alpha                           [m3/(Pa s)]
///
/// p_bhp is an explicit scalar input in this slice. It is not a globally
/// numbered well unknown and no control equation is created here.
struct WellConnectionPressureDrawdownPhaseRateLinearization3P {
    static constexpr std::string_view convention =
        pressure_drawdown_phase_rate_convention;

    mpmc::flow::NaturalVariableStateIdentity3P
        state_identity;
    mpmc::flow::PhaseSlot3 phase{
        mpmc::flow::PhaseSlot3::phase0};

    double hydraulic_conductance_m3_per_pa_s{};
    double cell_phase_pressure_pa{};
    double bottom_hole_pressure_pa{};
    double pressure_drawdown_pa{};
    double phase_volumetric_rate_m3_per_s{};

    std::vector<double>
        reservoir_natural_variable_gradient;
    double bottom_hole_pressure_derivative_m3_per_pa_s{};

    [[nodiscard]] double
    d_phase_rate_d_reservoir_unknown(
        std::size_t column) const {
        return reservoir_natural_variable_gradient.at(
            column);
    }

    [[nodiscard]] double
    d_phase_rate_d_bottom_hole_pressure() const noexcept {
        return bottom_hole_pressure_derivative_m3_per_pa_s;
    }
};

/// Build one phase's connection rate from the already validated Peaceman
/// connection plus the existing local phase mobility/pressure linearization.
///
/// No component split, molar/mass conversion, connection aggregation,
/// hydrostatic/friction wellbore correction, source residual or PETSc object is
/// introduced by this function.
[[nodiscard]] inline
WellConnectionPressureDrawdownPhaseRateLinearization3P
make_pressure_drawdown_phase_rate_linearization_3p(
    const mpmc::well::PeacemanWellIndex3D&
        connection,
    const mpmc::flow::LocalPhaseMobilityLinearization3P&
        mobility,
    mpmc::flow::PhaseSlot3 phase,
    double bottom_hole_pressure_pa) {
    if (!std::isfinite(bottom_hole_pressure_pa) ||
        !(bottom_hole_pressure_pa > 0.0)) {
        throw std::invalid_argument(
            "mpmc::well_discretization::make_pressure_drawdown_phase_rate_linearization_3p: bottom-hole pressure must be finite and strictly positive [Pa]");
    }

    const auto conductance =
        make_hydraulic_conductance_linearization_3p(
            connection,
            mobility,
            phase);

    const std::size_t phase_index =
        static_cast<std::size_t>(phase);
    const std::size_t q =
        conductance.state_identity.layout
            .unknown_count();

    const double phase_pressure =
        mobility.phase_pressure_pa[
            phase_index];
    const auto& pressure_gradient =
        mobility.phase_pressure_gradient[
            phase_index];

    if (!std::isfinite(phase_pressure) ||
        !(phase_pressure > 0.0) ||
        pressure_gradient.size() != q) {
        throw std::invalid_argument(
            "mpmc::well_discretization::make_pressure_drawdown_phase_rate_linearization_3p: phase pressure value/gradient shape is invalid");
    }
    for (double derivative :
         pressure_gradient) {
        if (!std::isfinite(derivative)) {
            throw std::invalid_argument(
                "mpmc::well_discretization::make_pressure_drawdown_phase_rate_linearization_3p: phase pressure gradient contains non-finite derivative");
        }
    }

    const double drawdown =
        phase_pressure -
        bottom_hole_pressure_pa;
    const double phase_rate =
        conductance
            .hydraulic_conductance_m3_per_pa_s *
        drawdown;
    if (!std::isfinite(drawdown) ||
        !std::isfinite(phase_rate)) {
        throw std::range_error(
            "mpmc::well_discretization::make_pressure_drawdown_phase_rate_linearization_3p: pressure drawdown or phase rate is outside representable range");
    }

    std::vector<double> reservoir_gradient(
        q,
        0.0);
    for (std::size_t column = 0U;
         column < q;
         ++column) {
        reservoir_gradient[column] =
            conductance
                .d_hydraulic_conductance(
                    column) *
                drawdown +
            conductance
                .hydraulic_conductance_m3_per_pa_s *
                pressure_gradient[column];
        if (!std::isfinite(
                reservoir_gradient[column])) {
            throw std::range_error(
                "mpmc::well_discretization::make_pressure_drawdown_phase_rate_linearization_3p: reservoir phase-rate derivative is outside representable range");
        }
    }

    const double bhp_derivative =
        -conductance
             .hydraulic_conductance_m3_per_pa_s;
    if (!std::isfinite(bhp_derivative)) {
        throw std::range_error(
            "mpmc::well_discretization::make_pressure_drawdown_phase_rate_linearization_3p: BHP derivative is outside representable range");
    }

    return {
        conductance.state_identity,
        phase,
        conductance
            .hydraulic_conductance_m3_per_pa_s,
        phase_pressure,
        bottom_hole_pressure_pa,
        drawdown,
        phase_rate,
        std::move(reservoir_gradient),
        bhp_derivative};
}

} // namespace mpmc::well_discretization

#endif // MPMC_WELL_DISCRETIZATION_PRESSURE_DRAWDOWN_RATE_HPP
