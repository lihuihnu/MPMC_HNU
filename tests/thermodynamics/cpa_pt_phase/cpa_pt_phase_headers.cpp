#include <mpmc/thermodynamics/cpa_pt_phase.hpp>

bool cpa_pt_phase_public_header_self_contained() {
    return mpmc::thermodynamics::cpa_pt_convention.size() > 0U;
}
