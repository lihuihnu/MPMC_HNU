#include <mpmc/flow_discretization_petsc/post_snes_phase_transition_controller.hpp>

bool post_snes_phase_transition_controller_header() {
    return !mpmc::flow_discretization_petsc::
        post_snes_phase_transition_controller_convention.empty();
}
