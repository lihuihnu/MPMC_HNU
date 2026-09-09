#include <mpmc/flash/pt_stability.hpp>
#include <mpmc/flash/pt_stability.hpp>

bool stability_generic_header() {
    const std::vector<double> z{1};
    const mpmc::flash::StabilityPhase phase{{0},0,true};
    return mpmc::flash::tangent_plane_distance(z,z,phase,phase).value==0;
}
