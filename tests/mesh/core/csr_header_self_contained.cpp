#include <mpmc/mesh/csr_adjacency.hpp>

#include <span>
#include <type_traits>
#include <utility>

static_assert(
    std::is_same_v<decltype(std::declval<const mpmc::mesh::CsrAdjacency&>().indices()),
                   std::span<const mpmc::mesh::LocalIndex>>);
