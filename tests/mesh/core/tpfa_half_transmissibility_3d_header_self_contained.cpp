#include <mpmc/mesh/tpfa_half_transmissibility_3d.hpp>

#include <type_traits>

static_assert(std::is_trivially_copyable_v<
              mpmc::mesh::TpfaAreaScaledHalfTransmissibility3D>);
