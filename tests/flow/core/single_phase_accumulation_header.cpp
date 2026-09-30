#include <mpmc/flow/single_phase_accumulation.hpp>

bool single_phase_accumulation_header() {
    const mpmc::flow::NaturalVariableLayout1P layout{3U};
    return layout.unknown_count() == 4U;
}
