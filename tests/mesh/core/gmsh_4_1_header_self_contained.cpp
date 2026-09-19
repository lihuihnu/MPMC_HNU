#include <mpmc/mesh/gmsh_4_1.hpp>

#include <type_traits>

static_assert(std::is_move_constructible_v<mpmc::mesh::Gmsh41ImportResult>);
static_assert(!std::is_default_constructible_v<mpmc::mesh::Gmsh41ImportResult>);
