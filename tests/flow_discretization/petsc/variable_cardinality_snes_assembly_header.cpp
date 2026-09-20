#include <mpmc/flow_discretization_petsc/variable_cardinality_snes_assembly.hpp>

bool variable_cardinality_snes_assembly_header() {
    return !mpmc::flow_discretization_petsc::
        variable_cardinality_snes_assembly_convention.empty();
}
