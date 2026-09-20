#include <mpmc/flow_discretization/normalized_component_face_contribution.hpp>

#include <string_view>
#include <type_traits>
#include <utility>

static_assert(
    mpmc::flow_discretization::
        NormalizedComponentFaceContributionLinearization3D::
            convention ==
    mpmc::flow_discretization::
        normalized_component_face_contribution_convention);

static_assert(
    mpmc::flow_discretization::
        NormalizedComponentFaceContributionLinearization3D::
            bulk_volume_derivative_is_zero);

static_assert(
    std::is_same_v<
        decltype(
            (std::declval<
                const mpmc::flow_discretization::
                    NormalizedComponentFaceContributionLinearization3D&>()
                 .owner_total_contribution_mol_per_bulk_m3_s)),
        const double&>);

bool normalized_component_face_contribution_header() {
    return true;
}
