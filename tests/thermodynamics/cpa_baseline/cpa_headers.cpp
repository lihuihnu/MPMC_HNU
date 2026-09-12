#include <mpmc/thermodynamics/cpa_parameters.hpp>
#include <mpmc/thermodynamics/cpa_association.hpp>
#include <mpmc/thermodynamics/cpa_phase.hpp>

bool cpa_public_headers_self_contained() {
    return mpmc::thermodynamics::cpa_profile.size() > 0U;
}
