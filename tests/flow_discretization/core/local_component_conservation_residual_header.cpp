#include <mpmc/flow_discretization/local_component_conservation_residual.hpp>

#include <string_view>
#include <type_traits>
#include <utility>
#include <vector>

static_assert(
    mpmc::flow_discretization::
        LocalComponentConservationResidualLinearization3D::
            convention ==
    mpmc::flow_discretization::
        local_component_conservation_residual_convention);

static_assert(
    std::is_same_v<
        decltype(
            (std::declval<
                const mpmc::flow_discretization::
                    LocalComponentConservationResidualLinearization3D&>()
                 .neighbour_blocks)),
        const std::vector<
            mpmc::flow_discretization::
                LocalComponentConservationNeighbourBlock3D>&>);

bool local_component_conservation_residual_header() {
    return true;
}
