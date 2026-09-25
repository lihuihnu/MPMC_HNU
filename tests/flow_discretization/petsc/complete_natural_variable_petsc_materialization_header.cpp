#include <mpmc/flow_discretization_petsc/complete_natural_variable_petsc_materialization.hpp>

#include <string_view>

static_assert(
    mpmc::flow_discretization_petsc::
        complete_natural_variable_petsc_materialization_convention ==
    std::string_view{
        "flow_discretization_petsc/complete-natural-variable-petsc-materialization/v1"});

bool complete_natural_variable_petsc_materialization_header() {
    return true;
}
