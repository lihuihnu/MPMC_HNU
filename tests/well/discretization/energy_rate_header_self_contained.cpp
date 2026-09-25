#include <mpmc/well_discretization/energy_rate.hpp>

#include <string_view>

static_assert(
    mpmc::well_discretization::
        WellConnectionAdvectiveEnergyRateLinearization3P::
            convention ==
    mpmc::well_discretization::
        connection_advective_energy_rate_convention);

int main() {
    return
        mpmc::well_discretization::
            connection_advective_energy_rate_convention ==
        std::string_view{
            "well-discretization/connection-advective-energy-rate/directional-enthalpy/v1"}
        ? 0
        : 1;
}
