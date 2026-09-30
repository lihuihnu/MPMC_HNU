#include <mpmc/flow/single_phase_properties.hpp>

bool single_phase_properties_header() {
    const mpmc::flow::NaturalVariableLayout1P layout{3U};
    return layout.unknown_count() == 4U;
}
