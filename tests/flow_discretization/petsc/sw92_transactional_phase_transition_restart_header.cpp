#include <mpmc/flow_discretization_petsc/sw92_transactional_phase_transition_restart.hpp>

bool sw92_transactional_phase_transition_restart_header() {
    return mpmc::flow_discretization_petsc::
               sw92_transactional_phase_transition_restart_convention ==
           "flow_discretization_petsc/sw92-transactional-phase-transition-restart/no-cross-cardinality-face/v1";
}
