#include <mpmc/flow/phase_set_transition.hpp>
#include <mpmc/flow/phase_set_transition_flash_adapter.hpp>

bool phase_set_transition_header() {
    return
        mpmc::flow::
            phase_set_transition_convention ==
            "flow/natural-variable/phase-set-transition/v1" &&
        mpmc::flow::
            phase_set_transition_flash_adapter_convention ==
            "flow/flash/phase-set-transition-adapter/v1";
}
