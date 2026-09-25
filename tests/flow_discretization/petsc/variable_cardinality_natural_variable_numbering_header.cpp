#include <mpmc/flow_discretization_petsc/variable_cardinality_natural_variable_numbering.hpp>

bool variable_cardinality_natural_variable_numbering_header() {
    return
        mpmc::flow_discretization_petsc::
            variable_cardinality_natural_variable_numbering_convention ==
        "flow_discretization_petsc/variable-cardinality-natural-variable-numbering/v1";
}
