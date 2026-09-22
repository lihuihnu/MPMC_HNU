#include <mpmc/flow_discretization_petsc/post_snes_phase_transition_handoff_scanner.hpp>

#include <string_view>

bool post_snes_phase_transition_handoff_scanner_header() {
    return
        mpmc::flow_discretization_petsc::
            post_snes_phase_transition_handoff_scanner_convention ==
        std::string_view{
            "flow_discretization_petsc/post-snes-phase-transition-handoff-scanner/v1"};
}
