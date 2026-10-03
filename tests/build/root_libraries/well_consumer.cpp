#include <mpmc/mesh/entity.hpp>
#include <mpmc/well/peaceman_well_index_3d.hpp>

// Public type and mesh dependency probe, not a Peaceman formula validation.
int main() {
    mpmc::well::CellPeacemanWellIndex3D completion;
    completion.cell = mpmc::mesh::LocalIndex{7};
    completion.connection.cell_dimensions_m = {1.0, 2.0, 3.0};
    completion.connection.well_direction = mpmc::well::AxisAlignedWellDirection3D::z;
    return completion.cell.value() == 7 &&
                   completion.connection.cell_dimensions_m.dz_m == 3.0 &&
                   completion.connection.well_direction == mpmc::well::AxisAlignedWellDirection3D::z
               ? 0
               : 1;
}
