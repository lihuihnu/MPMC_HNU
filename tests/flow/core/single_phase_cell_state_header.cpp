#include <mpmc/flow/single_phase_cell_state.hpp>

bool single_phase_cell_state_header() {
    const mpmc::flow::NaturalVariableLayout1P layout{3U};
    return layout.unknown_count() == 4U;
}
