#include <mpmc/well_discretization/hydraulic_conductance.hpp>

#include <string_view>

static_assert(
    mpmc::well_discretization::
        WellConnectionHydraulicConductanceLinearization3P::
            convention ==
    mpmc::well_discretization::
        hydraulic_conductance_convention);

static_assert(
    mpmc::well_discretization::
        WellConnectionHydraulicConductanceLinearization3P::
            well_index_derivative_is_zero);

int main() {
    return
        mpmc::well_discretization::
            hydraulic_conductance_convention ==
        std::string_view{
            "well-discretization/hydraulic-conductance/frozen-wi-times-local-phase-mobility/v1"}
        ? 0
        : 1;
}
