#include <mpmc/flow/two_phase_accumulation.hpp>

bool two_phase_accumulation_header() {
    const mpmc::flow::NaturalVariableLayout2P layout{3U};
    return layout.unknown_count() == 7U;
}
