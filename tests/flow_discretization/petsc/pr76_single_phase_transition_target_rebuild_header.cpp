#include <mpmc/flow_discretization_petsc/pr76_single_phase_transition_target_rebuild.hpp>

#include <string_view>

bool pr76_single_phase_transition_target_rebuild_header() {
    return
        mpmc::flow_discretization_petsc::
            pr76_single_phase_transition_target_rebuild_convention ==
        std::string_view{
            "flow_discretization_petsc/pr76-single-phase-transition-target-rebuild/v1"};
}
