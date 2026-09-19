#ifndef MPMC_MESH_TPFA_HALF_TRANSMISSIBILITY_3D_HPP
#define MPMC_MESH_TPFA_HALF_TRANSMISSIBILITY_3D_HPP

#include <mpmc/mesh/tpfa_half_connection_3d.hpp>

#include <cmath>
#include <limits>
#include <stdexcept>

namespace mpmc::mesh {

/// Area-scaled one-sided TPFA half transmissibility.
///
/// The underlying half-connection coefficient is retained verbatim for audit:
///
///   c_half = (n^T K d) / (d^T d)       [m]
///   T_half = A_face * c_half            [m3]
///
/// This remains a single cell-face half quantity. No owner/neighbour harmonic
/// combination, two-point face transmissibility, Darcy flux, or residual is
/// defined here.
struct TpfaAreaScaledHalfTransmissibility3D {
    TpfaHalfConnectionCoefficient3D half_connection;
    double face_area_m2;
    double half_transmissibility_m3;
};

namespace tpfa_half_transmissibility_3d_detail {

inline void validate_half_connection(
    const TpfaHalfConnectionCoefficient3D& half_connection) {
    if (!std::isfinite(
            half_connection
                .normal_permeability_displacement_m3) ||
        !std::isfinite(
            half_connection.squared_distance_m2) ||
        !std::isfinite(
            half_connection.coefficient_m) ||
        half_connection.squared_distance_m2 <= 0.0 ||
        half_connection
                .normal_permeability_displacement_m3 <
            0.0 ||
        half_connection.coefficient_m < 0.0) {
        throw std::invalid_argument(
            "mpmc::mesh::make_area_scaled_tpfa_half_transmissibility_3d: invalid half-connection coefficient");
    }

    switch (half_connection.projection) {
    case TpfaHalfConnectionProjection3D::
        positive_projection:
        if (half_connection
                    .normal_permeability_displacement_m3 <=
                0.0 ||
            half_connection.coefficient_m <= 0.0) {
            throw std::invalid_argument(
                "mpmc::mesh::make_area_scaled_tpfa_half_transmissibility_3d: positive projection requires positive numerator and coefficient");
        }
        break;
    case TpfaHalfConnectionProjection3D::
        zero_projection:
        if (half_connection
                    .normal_permeability_displacement_m3 !=
                0.0 ||
            half_connection.coefficient_m != 0.0) {
            throw std::invalid_argument(
                "mpmc::mesh::make_area_scaled_tpfa_half_transmissibility_3d: zero projection requires exact zero numerator and coefficient");
        }
        break;
    default:
        throw std::invalid_argument(
            "mpmc::mesh::make_area_scaled_tpfa_half_transmissibility_3d: invalid projection state");
    }

    const double reconstructed_coefficient_m =
        half_connection
            .normal_permeability_displacement_m3 /
        half_connection.squared_distance_m2;
    const double scale_m =
        std::max(
            std::abs(
                reconstructed_coefficient_m),
            std::abs(
                half_connection.coefficient_m));
    const double tolerance_m =
        4096.0 *
        std::numeric_limits<double>::epsilon() *
        scale_m;
    if (!std::isfinite(
            reconstructed_coefficient_m) ||
        std::abs(
            reconstructed_coefficient_m -
            half_connection.coefficient_m) >
            tolerance_m) {
        throw std::invalid_argument(
            "mpmc::mesh::make_area_scaled_tpfa_half_transmissibility_3d: half-connection numerator/denominator/coefficient are inconsistent");
    }
}

} // namespace tpfa_half_transmissibility_3d_detail

/// Scale a validated one-sided coefficient by a strictly positive face area.
///
/// Units:
///   face_area_m2                  [m2]
///   half_connection.coefficient_m [m]
///   half_transmissibility_m3      [m3]
///
/// A zero_projection is preserved exactly as zero m3. Negative or zero face
/// area is rejected because FaceGeometry3D already requires a finite positive
/// physical face measure.
[[nodiscard]] inline TpfaAreaScaledHalfTransmissibility3D
make_area_scaled_tpfa_half_transmissibility_3d(
    double face_area_m2,
    TpfaHalfConnectionCoefficient3D half_connection) {
    if (!std::isfinite(face_area_m2) ||
        face_area_m2 <= 0.0) {
        throw std::invalid_argument(
            "mpmc::mesh::make_area_scaled_tpfa_half_transmissibility_3d: face area must be finite and strictly positive");
    }

    tpfa_half_transmissibility_3d_detail::
        validate_half_connection(
            half_connection);

    const double half_transmissibility_m3 =
        face_area_m2 *
        half_connection.coefficient_m;
    if (!std::isfinite(
            half_transmissibility_m3) ||
        half_transmissibility_m3 < 0.0) {
        throw std::invalid_argument(
            "mpmc::mesh::make_area_scaled_tpfa_half_transmissibility_3d: area-scaled half transmissibility must be finite and non-negative");
    }

    if (half_connection.projection ==
            TpfaHalfConnectionProjection3D::
                zero_projection &&
        half_transmissibility_m3 != 0.0) {
        throw std::logic_error(
            "mpmc::mesh::make_area_scaled_tpfa_half_transmissibility_3d: zero projection did not remain zero after area scaling");
    }
    if (half_connection.projection ==
            TpfaHalfConnectionProjection3D::
                positive_projection &&
        half_transmissibility_m3 <= 0.0) {
        throw std::invalid_argument(
            "mpmc::mesh::make_area_scaled_tpfa_half_transmissibility_3d: positive projection must produce positive area-scaled half transmissibility");
    }

    return TpfaAreaScaledHalfTransmissibility3D{
        half_connection,
        face_area_m2,
        half_transmissibility_m3};
}

/// Owner-side area-scaled half transmissibility.
///
/// Boundary and internal owner half connections are both valid. The face area
/// is taken from the same FaceGeometry3D-derived geometric snapshot as the
/// owner half-connection inputs.
[[nodiscard]] inline TpfaAreaScaledHalfTransmissibility3D
make_owner_area_scaled_tpfa_half_transmissibility_3d(
    const CellFaceGeometricOperator3D& geometry,
    const CellCartesianDiagonalPermeability3D& permeability,
    LocalIndex face) {
    const auto half_connection =
        make_owner_tpfa_half_connection_coefficient_3d(
            geometry,
            permeability,
            face);
    const double face_area_m2 =
        geometry.transmissibility_geometry(face)
            .area_m2;
    return make_area_scaled_tpfa_half_transmissibility_3d(
        face_area_m2,
        half_connection);
}

/// Neighbour-side area-scaled half transmissibility for an internal face.
///
/// This independently scales the neighbour half coefficient by the common face
/// area. It does not combine owner and neighbour values.
[[nodiscard]] inline TpfaAreaScaledHalfTransmissibility3D
make_neighbour_area_scaled_tpfa_half_transmissibility_3d(
    const CellFaceGeometricOperator3D& geometry,
    const CellCartesianDiagonalPermeability3D& permeability,
    LocalIndex face) {
    const auto half_connection =
        make_neighbour_tpfa_half_connection_coefficient_3d(
            geometry,
            permeability,
            face);
    const double face_area_m2 =
        geometry.transmissibility_geometry(face)
            .area_m2;
    return make_area_scaled_tpfa_half_transmissibility_3d(
        face_area_m2,
        half_connection);
}

} // namespace mpmc::mesh

#endif // MPMC_MESH_TPFA_HALF_TRANSMISSIBILITY_3D_HPP
