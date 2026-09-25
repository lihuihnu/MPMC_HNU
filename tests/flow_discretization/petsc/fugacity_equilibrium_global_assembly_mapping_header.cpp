#include <mpmc/flow_discretization_petsc/fugacity_equilibrium_global_assembly_mapping.hpp>

#include <string_view>
#include <type_traits>
#include <utility>

static_assert(
    mpmc::flow_discretization_petsc::
        FugacityEquilibriumGlobalAssemblyEntries3D::
            convention ==
    mpmc::flow_discretization_petsc::
        fugacity_equilibrium_global_assembly_mapping_convention);

static_assert(
    std::is_same_v<
        decltype(
            (std::declval<
                const mpmc::flow_discretization_petsc::
                    AssemblyReadyFugacityJacobianEntry3D&>()
                 .petsc_global_column)),
        const PetscInt&>);

bool fugacity_equilibrium_global_assembly_mapping_header() {
    return true;
}
