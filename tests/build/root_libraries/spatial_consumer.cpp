#include <mpmc/mesh/entity.hpp>
#if MPMC_TEST_DISCRETIZATION
#include <mpmc/discretization/tpfa_half_transmissibility_3d.hpp>
#endif

// Isolated spatial-library usage probe. No AD, EOS, PETSc or scientific oracle.
int main() {
    const mpmc::mesh::GlobalEntityId first{7};
    const mpmc::mesh::GlobalEntityId second{8};
    if (!(first < second)) {
        return 1;
    }
#if MPMC_TEST_DISCRETIZATION
    // Synthetic zero-permeability half connection: zero coefficient [m],
    // positive squared distance [m2] and face area [m2] produce exact zero [m3].
    const mpmc::discretization::TpfaHalfConnectionCoefficient3D half{
        mpmc::discretization::TpfaHalfConnectionProjection3D::zero_projection,
        0.0, 1.0, 0.0};
    const auto scaled = mpmc::discretization::make_area_scaled_tpfa_half_transmissibility_3d(
        1.0, half);
    if (scaled.half_transmissibility_m3 != 0.0) {
        return 2;
    }
#endif
    return 0;
}
