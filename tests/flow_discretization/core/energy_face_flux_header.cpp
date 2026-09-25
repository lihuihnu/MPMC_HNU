#include <mpmc/flow_discretization/energy_face_flux.hpp>

#include <string_view>

static_assert(
    mpmc::flow_discretization::
        InternalEnergyFaceRateLinearization3D::
            convention ==
    mpmc::flow_discretization::
        internal_energy_face_rate_convention);

bool energy_face_flux_header() {
    return true;
}
