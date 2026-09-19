#include <mpmc/mesh/cell_face_geometric_operator_3d.hpp>

#include <type_traits>

static_assert(std::is_copy_constructible_v<
              mpmc::mesh::CellFaceGeometricOperator3D>);
static_assert(!std::is_default_constructible_v<
              mpmc::mesh::CellFaceGeometricOperator3D>);
