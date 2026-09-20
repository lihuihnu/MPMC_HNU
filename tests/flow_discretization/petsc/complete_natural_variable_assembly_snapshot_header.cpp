#include <mpmc/flow_discretization_petsc/complete_natural_variable_assembly_snapshot.hpp>

#include <string_view>
#include <type_traits>
#include <utility>

static_assert(
    mpmc::flow_discretization_petsc::
        CompleteNaturalVariableAssemblySnapshot3D::
            convention ==
    mpmc::flow_discretization_petsc::
        complete_natural_variable_assembly_snapshot_convention);

static_assert(
    std::is_same_v<
        decltype(
            (std::declval<
                const mpmc::flow_discretization_petsc::
                    CompleteNaturalVariableJacobianEntry3D&>()
                 .petsc_global_column)),
        const PetscInt&>);

bool complete_natural_variable_assembly_snapshot_header() {
    return true;
}
