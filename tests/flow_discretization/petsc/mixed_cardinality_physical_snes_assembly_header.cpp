#include <mpmc/flow_discretization_petsc/mixed_cardinality_physical_snes_assembly.hpp>

bool mixed_cardinality_physical_snes_assembly_header() {
    return !mpmc::flow_discretization_petsc::
        mixed_cardinality_physical_snes_assembly_convention.empty();
}
