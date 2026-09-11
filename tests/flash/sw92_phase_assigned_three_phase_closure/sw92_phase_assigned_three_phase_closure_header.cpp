#include <mpmc/flash/sw92_phase_assigned_three_phase_closure.hpp>

bool sw92_phase_assigned_three_phase_closure_header() {
    mpmc::flash::Sw92PhaseAssignedC2b2Result result;
    return !result.global_stability_proven &&
           !result.accepted_phase_set_published &&
           !result.morphology_resolved;
}
