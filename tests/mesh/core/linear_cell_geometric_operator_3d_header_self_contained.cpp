#include <mpmc/mesh/linear_cell_geometric_operator_3d.hpp>

bool linear_cell_geometric_operator_3d_header_self_contained() {
    using namespace mpmc::mesh;
    auto* factory = static_cast<CellFaceGeometricOperator3D (*)(const LinearMesh3D&)>(
        &make_cell_face_geometric_operator_3d);
    return factory != nullptr;
}
