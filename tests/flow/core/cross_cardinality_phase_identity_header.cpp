#include <mpmc/flow/cross_cardinality_phase_identity.hpp>

#include <string_view>

bool cross_cardinality_phase_identity_header() {
    using namespace mpmc::flow;
    return !cross_cardinality_phase_identity_convention.empty() &&
        CrossCardinalityAbsentPhaseSemantics::
                mobility_per_pa_s ==
            0.0;
}
