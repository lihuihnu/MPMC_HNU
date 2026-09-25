#include <mpmc/flow_discretization_petsc/cross_cardinality_tpfa_bridge.hpp>

bool cross_cardinality_tpfa_bridge_header() {
    return !mpmc::flow_discretization_petsc::
        cross_cardinality_tpfa_bridge_convention.empty();
}
