#include <mpmc/mesh/partition_snapshot.hpp>

#include <cstdint>
#include <span>
#include <type_traits>
#include <utility>

static_assert(sizeof(mpmc::mesh::PartitionRank) == sizeof(std::uint32_t));
static_assert(sizeof(mpmc::mesh::EntityOwnership) == sizeof(std::uint8_t));
static_assert(std::is_trivially_copyable_v<mpmc::mesh::PartitionRank>);
static_assert(!std::is_default_constructible_v<mpmc::mesh::PartitionSnapshot>);
static_assert(std::is_copy_constructible_v<mpmc::mesh::PartitionSnapshot>);
static_assert(std::is_move_constructible_v<mpmc::mesh::PartitionSnapshot>);
static_assert(!std::is_copy_assignable_v<mpmc::mesh::PartitionSnapshot>);
static_assert(!std::is_move_assignable_v<mpmc::mesh::PartitionSnapshot>);
static_assert(std::is_same_v<
              decltype(std::declval<const mpmc::mesh::PartitionSnapshot&>().global_ids(
                  mpmc::mesh::EntityKind::cell)),
              std::span<const mpmc::mesh::GlobalEntityId>>);
