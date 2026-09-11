#include <mpmc/flash/sw92_phase_assigned_no_w.hpp>

bool sw92_phase_assigned_no_w_header() {
    mpmc::flash::Sw92PhaseAssignedNoWResult result;
    return !result.global_stability_proven &&
           !result.accepted_phase_set_published &&
           !result.morphology_resolved;
}
