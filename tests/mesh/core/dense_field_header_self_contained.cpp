#include <mpmc/mesh/dense_field.hpp>

#include <span>
#include <type_traits>
#include <utility>

static_assert(!std::is_default_constructible_v<mpmc::mesh::DenseFieldSnapshot>);
static_assert(std::is_copy_constructible_v<mpmc::mesh::DenseFieldSnapshot>);
static_assert(std::is_move_constructible_v<mpmc::mesh::DenseFieldSnapshot>);
static_assert(!std::is_copy_assignable_v<mpmc::mesh::DenseFieldSnapshot>);
static_assert(!std::is_move_assignable_v<mpmc::mesh::DenseFieldSnapshot>);
static_assert(std::is_same_v<
              decltype(std::declval<const mpmc::mesh::DenseFieldSnapshot&>().values()),
              std::span<const double>>);
