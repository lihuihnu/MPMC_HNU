#include <mpmc/mesh/face_geometry_3d.hpp>

#include <type_traits>

static_assert(std::is_copy_constructible_v<mpmc::mesh::FaceGeometry3D>);
static_assert(!std::is_default_constructible_v<mpmc::mesh::FaceGeometry3D>);
