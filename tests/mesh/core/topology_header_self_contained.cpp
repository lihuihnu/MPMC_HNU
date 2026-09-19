#include <mpmc/mesh/topology.hpp>

#include <span>
#include <type_traits>
#include <utility>

static_assert(!std::is_default_constructible_v<mpmc::mesh::Topology>);
static_assert(std::is_copy_constructible_v<mpmc::mesh::Topology>);
static_assert(std::is_move_constructible_v<mpmc::mesh::Topology>);
static_assert(!std::is_copy_assignable_v<mpmc::mesh::Topology>);
static_assert(!std::is_move_assignable_v<mpmc::mesh::Topology>);
static_assert(std::is_same_v<
              decltype(std::declval<const mpmc::mesh::Topology&>().global_ids(
                  mpmc::mesh::EntityKind::cell)),
              std::span<const mpmc::mesh::GlobalEntityId>>);
