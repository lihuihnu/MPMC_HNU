#include <mpmc/flash/sw92_phase_assigned_pt.hpp>

bool sw92_phase_assigned_pt_header() {
    mpmc::flash::Sw92PhaseAssignedPtResult result;
    return !result.global_stability_proven &&
           !result.accepted_phase_set_published &&
           !result.morphology_resolved &&
           result.maximum_phase_count == 3U;
}
