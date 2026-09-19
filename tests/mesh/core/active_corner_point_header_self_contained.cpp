#include <mpmc/mesh/active_corner_point.hpp>

#include <type_traits>

static_assert(std::is_move_constructible_v<mpmc::mesh::ActiveCornerPointGrid>);
static_assert(!std::is_default_constructible_v<mpmc::mesh::ActiveCornerPointGrid>);
