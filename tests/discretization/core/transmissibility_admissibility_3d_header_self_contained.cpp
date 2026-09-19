#include <mpmc/discretization/transmissibility_admissibility_3d.hpp>

#include <type_traits>

static_assert(std::is_same_v<
              decltype(&mpmc::discretization::
                           classify_internal_face_transmissibility_admissibility),
              mpmc::discretization::CombinedTransmissibilityAdmissibility3D (*)(
                  const mpmc::mesh::CellFaceGeometricOperator3D&,
                  const mpmc::mesh::CellCartesianDiagonalPermeability3D&,
                  mpmc::mesh::LocalIndex,
                  mpmc::discretization::
                      TransmissibilityGeometryAdmissibilityPolicy3D,
                  mpmc::discretization::
                      KOrthogonalityAdmissibilityPolicy3D)>);
