#include <mpmc/discretization/peaceman_well_index_3d.hpp>

#include <string_view>

static_assert(
    mpmc::discretization::
        PeacemanWellIndex3D::convention ==
    mpmc::discretization::
        peaceman_well_index_3d_convention);

int main() {
    return
        mpmc::discretization::
            peaceman_well_index_3d_convention ==
        std::string_view{
            "discretization/peaceman-well-index/axis-aligned-cartesian-diagonal-k/v1"}
        ? 0
        : 1;
}
