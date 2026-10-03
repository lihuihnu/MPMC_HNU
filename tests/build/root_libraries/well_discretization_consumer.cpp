#include <mpmc/discretization/tpfa_half_transmissibility_3d.hpp>
#include <mpmc/well_discretization/cell_source_adapter.hpp>

// Link only well_discretization: its public adapter must bring the well, flow,
// flow_discretization and spatial headers without EOS/AD include directories.
// Synthetic type/zero-transmissibility probe, not a physical well calculation.
int main() {
    mpmc::well_discretization::WellConnectionCellSourceAdapterResult3P adapter;
    adapter.cell_source.component_ids = {"first", "second"};
    const mpmc::discretization::TpfaHalfConnectionCoefficient3D half{
        mpmc::discretization::TpfaHalfConnectionProjection3D::zero_projection,
        0.0, 1.0, 0.0};
    const auto scaled = mpmc::discretization::make_area_scaled_tpfa_half_transmissibility_3d(
        1.0, half);
    return adapter.cell_source.component_count() == 2 &&
                   scaled.half_transmissibility_m3 == 0.0
               ? 0
               : 1;
}
