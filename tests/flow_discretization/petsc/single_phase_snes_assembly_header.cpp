#include <mpmc/flow_discretization_petsc/single_phase_snes_assembly.hpp>

bool single_phase_snes_assembly_header() {
    return
        mpmc::flow_discretization_petsc::
            single_phase_snes_assembly_convention ==
        "flow_discretization_petsc/single-phase-snes-assembly/v1";
}
