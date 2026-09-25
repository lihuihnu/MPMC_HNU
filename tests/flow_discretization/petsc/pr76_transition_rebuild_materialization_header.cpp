#include <mpmc/flow_discretization_petsc/pr76_transition_rebuild_materialization.hpp>

#include <string_view>

bool pr76_transition_rebuild_materialization_header() {
    return
        mpmc::flow_discretization_petsc::
            pr76_transition_rebuild_materialization_convention ==
        std::string_view{
            "flow_discretization_petsc/pr76-transition-rebuild-materialization/v1"};
}
