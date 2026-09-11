#include <mpmc/flash/sw92_phase_assigned_three_phase.hpp>

bool sw92_phase_assigned_three_phase_header() {
    mpmc::flash::Sw92PhaseAssignedThreePhaseOptions options;
    return options.max_evaluations > 0U &&
           mpmc::flash::Sw92PhaseAssignedThreePhaseResult::maximum_phase_count == 3U &&
           !mpmc::flash::Sw92PhaseAssignedThreePhaseResult::global_stability_proven &&
           !mpmc::flash::Sw92PhaseAssignedThreePhaseResult::accepted_phase_set_published;
}
