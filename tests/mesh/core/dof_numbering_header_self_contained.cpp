#include <mpmc/mesh/dof_numbering.hpp>

#include <cstdint>
#include <type_traits>

static_assert(sizeof(mpmc::mesh::GlobalDofIndex) == sizeof(std::uint64_t));
static_assert(sizeof(mpmc::mesh::GlobalEntityOrdinal) == sizeof(std::uint64_t));
static_assert(std::is_trivially_copyable_v<mpmc::mesh::GlobalDofIndex>);
static_assert(std::is_trivially_copyable_v<mpmc::mesh::GlobalEntityOrdinal>);
static_assert(!std::is_default_constructible_v<mpmc::mesh::DofNumberingSnapshot>);
static_assert(std::is_copy_constructible_v<mpmc::mesh::DofNumberingSnapshot>);
static_assert(std::is_move_constructible_v<mpmc::mesh::DofNumberingSnapshot>);
static_assert(!std::is_copy_assignable_v<mpmc::mesh::DofNumberingSnapshot>);
static_assert(!std::is_move_assignable_v<mpmc::mesh::DofNumberingSnapshot>);
