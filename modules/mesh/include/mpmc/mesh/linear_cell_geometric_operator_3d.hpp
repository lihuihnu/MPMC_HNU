#ifndef MPMC_MESH_LINEAR_CELL_GEOMETRIC_OPERATOR_3D_HPP
#define MPMC_MESH_LINEAR_CELL_GEOMETRIC_OPERATOR_3D_HPP

#include <mpmc/mesh/cell_face_geometric_operator_3d.hpp>
#include <mpmc/mesh/linear_cell_mesh_3d.hpp>

namespace mpmc::mesh {

/// Bridge for a validated linear tetrahedron/hexahedron/wedge/pyramid mesh,
/// including conforming mixtures. Uses the same vertex-mean reference points
/// as make_linear_mesh_3d's normal orientation; these are not volume centroids.
/// Strict positive distances remain required: supported I/O does not imply
/// admissibility for TPFA or any other numerical discretization.
[[nodiscard]] inline CellFaceGeometricOperator3D
make_cell_face_geometric_operator_3d(const LinearMesh3D& mesh) {
    return cell_face_geometric_operator_3d_detail::from_vertex_means(
        mesh.topology, mesh.vertex_coordinates_m, mesh.face_geometry, false);
}

} // namespace mpmc::mesh
#endif // MPMC_MESH_LINEAR_CELL_GEOMETRIC_OPERATOR_3D_HPP
