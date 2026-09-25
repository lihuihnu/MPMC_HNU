#include <mpmc/flow_discretization/conservative_component_face_rate_scatter.hpp>

#include <string_view>
#include <type_traits>
#include <utility>

static_assert(
    mpmc::flow_discretization::
        ConservativeComponentFaceRateScatterLinearization3D::
            convention ==
    mpmc::flow_discretization::
        conservative_component_face_rate_scatter_convention);

static_assert(
    std::is_same_v<
        decltype(
            (std::declval<
                const mpmc::flow_discretization::
                    ConservativeComponentFaceRateScatterLinearization3D&>()
                 .owner_total_face_rate_mol_per_s)),
        const double&>);

bool conservative_component_face_rate_scatter_header() {
    return true;
}
