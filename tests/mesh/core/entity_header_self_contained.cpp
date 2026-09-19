#include <mpmc/mesh/entity.hpp>

#include <type_traits>

static_assert(std::is_trivially_copyable_v<mpmc::mesh::LocalIndex>);
static_assert(std::is_trivially_copyable_v<mpmc::mesh::GlobalEntityId>);
