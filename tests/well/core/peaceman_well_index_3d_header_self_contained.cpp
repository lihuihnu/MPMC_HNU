#include <mpmc/well/peaceman_well_index_3d.hpp>

#include <string_view>

static_assert(
    mpmc::well::
        PeacemanWellIndex3D::convention ==
    mpmc::well::
        peaceman_well_index_3d_convention);

int main() {
    return
        mpmc::well::
            peaceman_well_index_3d_convention ==
        std::string_view{
            "well/peaceman-well-index/axis-aligned-cartesian-diagonal-k/v1"}
        ? 0
        : 1;
}
