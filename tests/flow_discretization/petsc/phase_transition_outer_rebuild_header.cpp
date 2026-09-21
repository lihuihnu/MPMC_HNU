#include <mpmc/flow_discretization_petsc/phase_transition_outer_rebuild.hpp>

bool phase_transition_outer_rebuild_header() {
    return !mpmc::flow_discretization_petsc::
        phase_transition_outer_rebuild_convention.empty() &&
        !mpmc::flow_discretization_petsc::
        frozen_absent_phase_coordinate_registry_convention.empty();
}
