#include <mpmc/discretization/tpfa_half_transmissibility_3d.hpp>

#include <type_traits>

static_assert(std::is_trivially_copyable_v<
              mpmc::discretization::TpfaAreaScaledHalfTransmissibility3D>);
