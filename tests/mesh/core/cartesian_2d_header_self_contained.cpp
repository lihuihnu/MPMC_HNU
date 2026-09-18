#include <mpmc/mesh/cartesian_2d.hpp>

#include <type_traits>

static_assert(std::is_same_v<
              decltype(mpmc::mesh::make_cartesian_topology_2d(1U, 1U)),
              mpmc::mesh::Topology>);
