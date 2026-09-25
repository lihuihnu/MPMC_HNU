#include <mpmc/flow_discretization_petsc/post_snes_sw92_profile_c_phase_transition_scanner.hpp>

bool post_snes_sw92_profile_c_phase_transition_scanner_header() {
    return mpmc::flow_discretization_petsc::
               post_snes_sw92_profile_c_phase_transition_scanner_convention ==
           "flow_discretization_petsc/post-snes-sw92-profile-c-authoritative-target/v1";
}
