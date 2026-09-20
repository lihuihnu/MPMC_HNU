#include <mpmc/flow/fugacity_equilibrium_linearization.hpp>

#include <string_view>

static_assert(
    mpmc::flow::
        FugacityEquilibriumResidualLinearization3P::
            convention ==
    mpmc::flow::
        fugacity_equilibrium_linearization_convention);

bool fugacity_equilibrium_linearization_header() {
    return true;
}
