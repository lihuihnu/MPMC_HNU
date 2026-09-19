#include <mpmc/mesh/vtu.hpp>

#include <type_traits>

static_assert(std::is_move_constructible_v<mpmc::mesh::VtuImportResult>);
static_assert(!std::is_default_constructible_v<mpmc::mesh::VtuImportResult>);
