#include <mpmc/mesh/shared_entity_plan.hpp>

#include <span>
#include <type_traits>
#include <utility>

static_assert(!std::is_default_constructible_v<mpmc::mesh::SharedEntityPlan>);
static_assert(std::is_copy_constructible_v<mpmc::mesh::SharedEntityPlan>);
static_assert(std::is_move_constructible_v<mpmc::mesh::SharedEntityPlan>);
static_assert(!std::is_copy_assignable_v<mpmc::mesh::SharedEntityPlan>);
static_assert(!std::is_move_assignable_v<mpmc::mesh::SharedEntityPlan>);
static_assert(std::is_same_v<
              decltype(std::declval<const mpmc::mesh::SharedEntityPlan&>().send_entities()),
              std::span<const mpmc::mesh::HaloEntityRef>>);
static_assert(std::is_same_v<
              decltype(std::declval<const mpmc::mesh::SharedEntityPlan&>().neighbors()),
              std::span<const mpmc::mesh::NeighborExchangeRange>>);
