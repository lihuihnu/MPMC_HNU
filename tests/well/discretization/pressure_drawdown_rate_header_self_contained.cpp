#include <mpmc/well_discretization/pressure_drawdown_rate.hpp>

#include <string_view>

static_assert(
    mpmc::well_discretization::
        WellConnectionPressureDrawdownPhaseRateLinearization3P::
            convention ==
    mpmc::well_discretization::
        pressure_drawdown_phase_rate_convention);

bool pressure_drawdown_rate_header_self_contained() {
    return
        mpmc::well_discretization::
            pressure_drawdown_phase_rate_convention ==
        std::string_view{
            "well-discretization/pressure-drawdown-phase-rate/cell-minus-bhp/v1"};
}
