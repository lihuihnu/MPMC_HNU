#include <mpmc/flow_discretization/owned_multi_cell_component_conservation.hpp>

#include <string_view>
#include <type_traits>
#include <utility>

static_assert(
    mpmc::flow_discretization::
        SerialOwnedMultiCellComponentConservationSnapshot3D::
            convention ==
    mpmc::flow_discretization::
        serial_owned_multi_cell_component_conservation_convention);

static_assert(
    std::is_same_v<
        decltype(
            (std::declval<
                const mpmc::flow_discretization::
                    OwnedCellComponentConservationRow3D&>()
                 .cell_global)),
        const mpmc::mesh::GlobalEntityId&>);

bool owned_multi_cell_component_conservation_header() {
    return true;
}
