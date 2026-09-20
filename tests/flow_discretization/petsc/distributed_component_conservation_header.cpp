#include <mpmc/flow_discretization_petsc/distributed_component_conservation.hpp>

#include <string_view>
#include <type_traits>
#include <utility>

static_assert(
    mpmc::flow_discretization_petsc::
        DistributedOwnedMultiCellComponentConservationSnapshot3D::
            convention ==
    mpmc::flow_discretization_petsc::
        distributed_owned_component_conservation_convention);

static_assert(
    std::is_same_v<
        decltype(
            (std::declval<
                const mpmc::flow_discretization_petsc::
                    DistributedCellStateBinding3D&>()
                 .cell_global)),
        const mpmc::mesh::GlobalEntityId&>);

bool distributed_component_conservation_header() {
    return true;
}
