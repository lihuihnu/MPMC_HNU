#ifndef MPMC_DISCRETIZATION_TPFA_HALF_CONNECTION_3D_HPP
#define MPMC_DISCRETIZATION_TPFA_HALF_CONNECTION_3D_HPP

#include <mpmc/mesh/cell_face_geometric_operator_3d.hpp>
#include <mpmc/mesh/permeability_tensor_3d.hpp>

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>

namespace mpmc::discretization {

/// Projection state for a single TPFA cell-face half connection.
///
/// positive_projection means n^T K d is strictly positive above the
/// permeability/length-scaled floating-point tolerance.
/// zero_projection preserves a physically valid zero-permeability projection
/// as a zero coefficient rather than rejecting it as invalid geometry.
enum class TpfaHalfConnectionProjection3D {
    positive_projection,
    zero_projection
};

/// One-sided TPFA coefficient for a single cell-face half connection.
///
/// Sign convention:
/// - d points from this cell centroid to the face centroid [m].
/// - n is this cell's outward unit face normal [-].
/// - K is this cell's permeability tensor [m2].
///
/// The core scalar is
///
///   c = (n^T K d) / (d^T d)
///
/// with numerator [m3], denominator [m2], and coefficient c [m].
///
/// Face area is deliberately absent. Multiplying by area and combining owner
/// and neighbour half connections belong to later contracts.
struct TpfaHalfConnectionCoefficient3D {
    TpfaHalfConnectionProjection3D projection;
    double normal_permeability_displacement_m3;
    double squared_distance_m2;
    double coefficient_m;
};

namespace tpfa_half_connection_3d_detail {

inline void require_tensor(
    mpmc::mesh::CartesianDiagonalPermeabilityTensor3D tensor) {
    if (!std::isfinite(tensor.kxx_m2) ||
        !std::isfinite(tensor.kyy_m2) ||
        !std::isfinite(tensor.kzz_m2) ||
        tensor.kxx_m2 < 0.0 ||
        tensor.kyy_m2 < 0.0 ||
        tensor.kzz_m2 < 0.0) {
        throw std::invalid_argument(
            "mpmc::discretization::make_tpfa_half_connection_coefficient_3d: permeability tensor must be finite and non-negative");
    }
}

inline void require_unit_normal(mpmc::mesh::UnitVector3D normal) {
    if (!std::isfinite(normal.x) ||
        !std::isfinite(normal.y) ||
        !std::isfinite(normal.z)) {
        throw std::invalid_argument(
            "mpmc::discretization::make_tpfa_half_connection_coefficient_3d: outward normal must be finite");
    }

    const double magnitude =
        std::sqrt(
            normal.x * normal.x +
            normal.y * normal.y +
            normal.z * normal.z);
    constexpr double tolerance =
        128.0 *
        std::numeric_limits<double>::epsilon();
    if (!std::isfinite(magnitude) ||
        std::abs(magnitude - 1.0) >
            tolerance) {
        throw std::invalid_argument(
            "mpmc::discretization::make_tpfa_half_connection_coefficient_3d: outward normal must have unit length");
    }
}

inline void require_displacement(mpmc::mesh::Displacement3D displacement) {
    if (!std::isfinite(displacement.x_m) ||
        !std::isfinite(displacement.y_m) ||
        !std::isfinite(displacement.z_m)) {
        throw std::invalid_argument(
            "mpmc::discretization::make_tpfa_half_connection_coefficient_3d: cell-to-face displacement must be finite");
    }
}

[[nodiscard]] inline double maximum_permeability(
    mpmc::mesh::CartesianDiagonalPermeabilityTensor3D tensor) {
    return std::max(
        {tensor.kxx_m2,
         tensor.kyy_m2,
         tensor.kzz_m2});
}

[[nodiscard]] inline double displacement_magnitude(
    mpmc::mesh::Displacement3D displacement) {
    return std::sqrt(
        displacement.x_m * displacement.x_m +
        displacement.y_m * displacement.y_m +
        displacement.z_m * displacement.z_m);
}

} // namespace tpfa_half_connection_3d_detail

/// Construct the one-sided coefficient from explicit single-cell inputs.
///
/// Degeneracy/error rules:
/// - non-finite or negative permeability components are rejected;
/// - n must be finite unit length;
/// - d must be finite and d^T d must be strictly positive;
/// - n^T K d < -tol is rejected because it violates the outward-normal /
///   cell-to-face sign convention;
/// - |n^T K d| <= tol is retained as an explicit zero_projection with c=0.
///
/// This function does not decide geometry/K admissibility. Callers may compute
/// the scalar for diagnostics, but direct TPFA use remains gated by the
/// admissibility contracts.
[[nodiscard]] inline TpfaHalfConnectionCoefficient3D
make_tpfa_half_connection_coefficient_3d(
    mpmc::mesh::CartesianDiagonalPermeabilityTensor3D permeability,
    mpmc::mesh::UnitVector3D outward_unit_normal,
    mpmc::mesh::Displacement3D cell_to_face_displacement_m) {
    using namespace tpfa_half_connection_3d_detail;

    require_tensor(permeability);
    require_unit_normal(outward_unit_normal);
    require_displacement(cell_to_face_displacement_m);

    const double squared_distance_m2 =
        cell_to_face_displacement_m.x_m *
            cell_to_face_displacement_m.x_m +
        cell_to_face_displacement_m.y_m *
            cell_to_face_displacement_m.y_m +
        cell_to_face_displacement_m.z_m *
            cell_to_face_displacement_m.z_m;
    if (!std::isfinite(squared_distance_m2) ||
        squared_distance_m2 <= 0.0) {
        throw std::invalid_argument(
            "mpmc::discretization::make_tpfa_half_connection_coefficient_3d: d^T d must be finite and strictly positive");
    }

    const double kdx =
        permeability.kxx_m2 *
        cell_to_face_displacement_m.x_m;
    const double kdy =
        permeability.kyy_m2 *
        cell_to_face_displacement_m.y_m;
    const double kdz =
        permeability.kzz_m2 *
        cell_to_face_displacement_m.z_m;
    const double numerator_m3 =
        outward_unit_normal.x * kdx +
        outward_unit_normal.y * kdy +
        outward_unit_normal.z * kdz;
    if (!std::isfinite(numerator_m3)) {
        throw std::invalid_argument(
            "mpmc::discretization::make_tpfa_half_connection_coefficient_3d: n^T K d must be finite");
    }

    const double projection_scale_m3 =
        maximum_permeability(permeability) *
        displacement_magnitude(
            cell_to_face_displacement_m);
    if (!std::isfinite(projection_scale_m3)) {
        throw std::invalid_argument(
            "mpmc::discretization::make_tpfa_half_connection_coefficient_3d: permeability/displacement scale must be finite");
    }
    const double projection_tolerance_m3 =
        4096.0 *
        std::numeric_limits<double>::epsilon() *
        projection_scale_m3;

    if (numerator_m3 <
        -projection_tolerance_m3) {
        throw std::invalid_argument(
            "mpmc::discretization::make_tpfa_half_connection_coefficient_3d: n^T K d is negative under the outward-normal/cell-to-face sign convention");
    }

    if (std::abs(numerator_m3) <=
        projection_tolerance_m3) {
        return TpfaHalfConnectionCoefficient3D{
            TpfaHalfConnectionProjection3D::
                zero_projection,
            0.0,
            squared_distance_m2,
            0.0};
    }

    const double coefficient_m =
        numerator_m3 /
        squared_distance_m2;
    if (!std::isfinite(coefficient_m) ||
        coefficient_m <= 0.0) {
        throw std::invalid_argument(
            "mpmc::discretization::make_tpfa_half_connection_coefficient_3d: positive projection must produce a finite positive coefficient");
    }

    return TpfaHalfConnectionCoefficient3D{
        TpfaHalfConnectionProjection3D::
            positive_projection,
        numerator_m3,
        squared_distance_m2,
        coefficient_m};
}

/// Owner-side half connection for a face.
///
/// Uses the canonical owner's outward normal and owner-centroid->face-centroid
/// displacement. Boundary and internal faces are both valid owner half
/// connections.
[[nodiscard]] inline TpfaHalfConnectionCoefficient3D
make_owner_tpfa_half_connection_coefficient_3d(
    const mpmc::mesh::CellFaceGeometricOperator3D& geometry,
    const mpmc::mesh::CellCartesianDiagonalPermeability3D& permeability,
    mpmc::mesh::LocalIndex face) {
    if (geometry.cell_count() !=
        permeability.cell_count()) {
        throw std::invalid_argument(
            "mpmc::discretization::make_owner_tpfa_half_connection_coefficient_3d: geometry/permeability cell counts do not match");
    }

    const mpmc::mesh::LocalIndex owner =
        geometry.face_owner(face);
    const auto face_geometry =
        geometry.face_connection_geometry(face);
    return make_tpfa_half_connection_coefficient_3d(
        permeability.tensor(owner),
        face_geometry.owner_unit_normal,
        geometry.owner_to_face_displacement_m(face));
}

/// Neighbour-side half connection for an internal face.
///
/// The neighbour outward normal is exactly -owner_normal. The displacement is
/// independently taken from neighbour centroid to the same face centroid.
/// No harmonic/two-point combination is performed here.
[[nodiscard]] inline TpfaHalfConnectionCoefficient3D
make_neighbour_tpfa_half_connection_coefficient_3d(
    const mpmc::mesh::CellFaceGeometricOperator3D& geometry,
    const mpmc::mesh::CellCartesianDiagonalPermeability3D& permeability,
    mpmc::mesh::LocalIndex face) {
    if (geometry.cell_count() !=
        permeability.cell_count()) {
        throw std::invalid_argument(
            "mpmc::discretization::make_neighbour_tpfa_half_connection_coefficient_3d: geometry/permeability cell counts do not match");
    }

    const auto neighbour =
        geometry.face_neighbour(face);
    const auto displacement =
        geometry.neighbour_to_face_displacement_m(face);
    if (!neighbour.has_value() ||
        !displacement.has_value()) {
        throw std::invalid_argument(
            "mpmc::discretization::make_neighbour_tpfa_half_connection_coefficient_3d: internal neighbour half connection is required");
    }

    const auto owner_normal =
        geometry.face_connection_geometry(face)
            .owner_unit_normal;
    const mpmc::mesh::UnitVector3D neighbour_normal{
        -owner_normal.x,
        -owner_normal.y,
        -owner_normal.z};

    return make_tpfa_half_connection_coefficient_3d(
        permeability.tensor(*neighbour),
        neighbour_normal,
        *displacement);
}

} // namespace mpmc::discretization

#endif // MPMC_DISCRETIZATION_TPFA_HALF_CONNECTION_3D_HPP
