#include <mpmc/well_discretization/component_molar_rate.hpp>

#include <string_view>

static_assert(
    mpmc::well_discretization::
        WellConnectionComponentMolarRateLinearization3P::
            convention ==
    mpmc::well_discretization::
        connection_component_molar_rate_convention);

int main() {
    return
        mpmc::well_discretization::
            connection_component_molar_rate_convention ==
        std::string_view{
            "well-discretization/connection-component-molar-rate/three-phase-sum/v1"}
        ? 0
        : 1;
}
