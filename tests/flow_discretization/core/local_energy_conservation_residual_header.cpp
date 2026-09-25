#include <mpmc/flow_discretization/local_energy_conservation_residual.hpp>

#include <string_view>

static_assert(
    mpmc::flow_discretization::
        LocalEnergyConservationResidualLinearization3D::
            convention ==
    mpmc::flow_discretization::
        local_energy_conservation_residual_convention);

bool local_energy_conservation_residual_header() {
    return true;
}
