#include <mpmc/flow_discretization_petsc/natural_variable_snes_solver.hpp>

#include <string_view>
#include <type_traits>
#include <utility>

static_assert(
    mpmc::flow_discretization_petsc::
        NaturalVariableSnesSolveReport3D::
            convention ==
    mpmc::flow_discretization_petsc::
        natural_variable_snes_solver_convention);

static_assert(
    std::is_same_v<
        decltype(
            (std::declval<
                const mpmc::flow_discretization_petsc::
                    NaturalVariableSnesSolutionEntry3D&>()
                 .petsc_global_scalar)),
        const PetscInt&>);

bool natural_variable_snes_solver_header() {
    return true;
}
