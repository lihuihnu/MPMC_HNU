#include <mpmc/flow_discretization_petsc/distributed_energy_conservation.hpp>

#include <string_view>
#include <type_traits>
#include <utility>

static_assert(
    mpmc::flow_discretization_petsc::
        DistributedOwnedMultiCellEnergyConservationSnapshot3D::
            convention ==
    mpmc::flow_discretization_petsc::
        distributed_owned_energy_conservation_convention);

static_assert(
    std::is_same_v<
        decltype(
            (std::declval<
                const mpmc::flow_discretization_petsc::
                    OwnedCellPairEnergyJacobianBlock3D&>()
                 .column_cell_global)),
        const mpmc::mesh::GlobalEntityId&>);

bool distributed_energy_conservation_header() {
    return true;
}
