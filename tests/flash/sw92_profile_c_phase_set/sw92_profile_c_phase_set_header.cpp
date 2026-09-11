#include <mpmc/flash/sw92_profile_c_phase_set.hpp>

bool sw92_profile_c_phase_set_header() {
    mpmc::flash::Sw92ProfileCPtPhaseSetResult result;
    return !result.global_stability_proven &&
           !result.morphology_resolved &&
           !result.accepted_phase_set_published();
}
