#include <mpmc/flow_discretization_petsc/post_snes_pt_flash_phase_transition_scanner.hpp>

bool post_snes_pt_flash_phase_transition_scanner_header() {
    return !mpmc::flow_discretization_petsc::
        post_snes_pt_flash_phase_transition_scanner_convention.empty();
}
