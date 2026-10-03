#include <mpmc/discretization/tpfa_half_transmissibility_3d.hpp>
#include <mpmc/flow/natural_variable_cell_state.hpp>
#include <mpmc/flow_discretization/cell_source.hpp>
#include <mpmc/mesh/entity.hpp>

// Bridge usage probe: no thermodynamics/AD target is linked to this executable.
// Exercise each layer with synthetic values; no physical solver claim is made.
int main() {
    const mpmc::flow::NaturalVariableLayoutDescriptor layout{2, 1, {1}};
    mpmc::flow_discretization::CellSourceLinearization3D source;
    source.component_ids = {"first", "second"};
    const mpmc::mesh::GlobalEntityId cell{7};
    const mpmc::discretization::TpfaHalfConnectionCoefficient3D half{
        mpmc::discretization::TpfaHalfConnectionProjection3D::zero_projection,
        0.0, 1.0, 0.0};
    const auto scaled = mpmc::discretization::make_area_scaled_tpfa_half_transmissibility_3d(
        1.0, half);
    return layout.unknown_count() == 3 && source.component_count() == 2 &&
                   cell.value() == 7 && scaled.half_transmissibility_m3 == 0.0
               ? 0
               : 1;
}
