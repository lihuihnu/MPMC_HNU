#include <mpmc/mesh/permeability_tensor_3d.hpp>

#include <type_traits>

static_assert(std::is_copy_constructible_v<
              mpmc::mesh::CellCartesianDiagonalPermeability3D>);
static_assert(!std::is_default_constructible_v<
              mpmc::mesh::CellCartesianDiagonalPermeability3D>);
