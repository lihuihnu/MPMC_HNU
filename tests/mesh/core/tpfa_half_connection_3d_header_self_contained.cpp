#include <mpmc/mesh/tpfa_half_connection_3d.hpp>

#include <type_traits>

static_assert(std::is_trivially_copyable_v<
              mpmc::mesh::TpfaHalfConnectionCoefficient3D>);
