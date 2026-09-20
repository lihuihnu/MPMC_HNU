#include <mpmc/flow_discretization_petsc/fixed_three_phase_snes_assembly.hpp>

#include <string_view>

static_assert(
    mpmc::flow_discretization_petsc::
        FixedThreePhaseSnesAssemblyContext3D::
            convention ==
    mpmc::flow_discretization_petsc::
        fixed_three_phase_snes_assembly_convention);

bool fixed_three_phase_snes_assembly_header() {
    return true;
}
