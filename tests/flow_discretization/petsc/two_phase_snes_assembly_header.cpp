#include <mpmc/flow_discretization_petsc/two_phase_snes_assembly.hpp>

bool two_phase_snes_assembly_header() {
    return
        mpmc::flow_discretization_petsc::
            two_phase_snes_assembly_convention ==
        "flow_discretization_petsc/two-phase-snes-assembly/v1";
}
