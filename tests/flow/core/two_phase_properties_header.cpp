#include <mpmc/flow/two_phase_properties.hpp>

bool two_phase_properties_header() {
    const mpmc::flow::NaturalVariableLayout2P layout{3U};
    return layout.unknown_count() == 7U;
}
