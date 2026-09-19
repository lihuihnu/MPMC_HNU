#include <mpmc/discretization/tpfa_internal_face_transmissibility_snapshot_3d.hpp>

#include <type_traits>

static_assert(std::is_copy_constructible_v<
              mpmc::discretization::TpfaInternalFaceTransmissibilitySnapshot3D>);
static_assert(!std::is_default_constructible_v<
              mpmc::discretization::TpfaInternalFaceTransmissibilitySnapshot3D>);
