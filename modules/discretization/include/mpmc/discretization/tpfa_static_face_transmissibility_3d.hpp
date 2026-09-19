#ifndef MPMC_DISCRETIZATION_TPFA_STATIC_FACE_TRANSMISSIBILITY_3D_HPP
#define MPMC_DISCRETIZATION_TPFA_STATIC_FACE_TRANSMISSIBILITY_3D_HPP

#include <mpmc/discretization/tpfa_half_transmissibility_3d.hpp>

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>

namespace mpmc::discretization {

/// Symmetric disposition for combining two internal-face half transmissibilities.
///
/// The state is deliberately invariant under exchanging owner and neighbour.
enum class TpfaStaticFaceTransmissibilityDisposition3D {
    positive_harmonic_combination,
    zero_due_to_one_half,
    zero_due_to_both_halves
};

/// Static two-sided TPFA face transmissibility.
///
/// For two positive finite one-sided values T_o and T_n [m3]:
///
///   T_f = T_o T_n / (T_o + T_n)        [m3]
///
/// A zero side blocks the two-sided static connection and yields T_f=0.
/// Pressure, mobility, Darcy flux, and residual terms are outside this contract.
struct TpfaStaticFaceTransmissibility3D {
    TpfaStaticFaceTransmissibilityDisposition3D disposition;
    double face_area_m2;
    double face_transmissibility_m3;
};

namespace tpfa_static_face_transmissibility_3d_detail {

inline void validate_area_scaled_half(
    const TpfaAreaScaledHalfTransmissibility3D& half) {
    if (!std::isfinite(half.face_area_m2) ||
        half.face_area_m2 <= 0.0 ||
        !std::isfinite(
            half.half_transmissibility_m3) ||
        half.half_transmissibility_m3 < 0.0) {
        throw std::invalid_argument(
            "mpmc::discretization::combine_internal_face_tpfa_half_transmissibilities_3d: half area/transmissibility must be finite with positive area and non-negative transmissibility");
    }

    tpfa_half_transmissibility_3d_detail::
        validate_half_connection(
            half.half_connection);

    const double reconstructed_m3 =
        half.face_area_m2 *
        half.half_connection.coefficient_m;
    if (!std::isfinite(reconstructed_m3)) {
        throw std::invalid_argument(
            "mpmc::discretization::combine_internal_face_tpfa_half_transmissibilities_3d: reconstructed half transmissibility is non-finite");
    }

    const double scale_m3 =
        std::max(
            std::abs(reconstructed_m3),
            std::abs(
                half.half_transmissibility_m3));
    const double tolerance_m3 =
        4096.0 *
        std::numeric_limits<double>::epsilon() *
        scale_m3;
    if (std::abs(
            reconstructed_m3 -
            half.half_transmissibility_m3) >
        tolerance_m3) {
        throw std::invalid_argument(
            "mpmc::discretization::combine_internal_face_tpfa_half_transmissibilities_3d: area/coefficient/transmissibility are inconsistent");
    }

    switch (half.half_connection.projection) {
    case TpfaHalfConnectionProjection3D::
        positive_projection:
        if (half.half_transmissibility_m3 <=
            0.0) {
            throw std::invalid_argument(
                "mpmc::discretization::combine_internal_face_tpfa_half_transmissibilities_3d: positive projection requires positive half transmissibility");
        }
        break;
    case TpfaHalfConnectionProjection3D::
        zero_projection:
        if (half.half_transmissibility_m3 !=
            0.0) {
            throw std::invalid_argument(
                "mpmc::discretization::combine_internal_face_tpfa_half_transmissibilities_3d: zero projection requires exact zero half transmissibility");
        }
        break;
    default:
        throw std::invalid_argument(
            "mpmc::discretization::combine_internal_face_tpfa_half_transmissibilities_3d: invalid half projection state");
    }
}

[[nodiscard]] inline double common_face_area_m2(
    const TpfaAreaScaledHalfTransmissibility3D& owner,
    const TpfaAreaScaledHalfTransmissibility3D& neighbour) {
    const double scale_m2 =
        std::max(
            owner.face_area_m2,
            neighbour.face_area_m2);
    const double tolerance_m2 =
        4096.0 *
        std::numeric_limits<double>::epsilon() *
        scale_m2;
    if (std::abs(
            owner.face_area_m2 -
            neighbour.face_area_m2) >
        tolerance_m2) {
        throw std::invalid_argument(
            "mpmc::discretization::combine_internal_face_tpfa_half_transmissibilities_3d: owner/neighbour halves must use the same physical face area");
    }
    return 0.5 * owner.face_area_m2 +
           0.5 * neighbour.face_area_m2;
}

} // namespace tpfa_static_face_transmissibility_3d_detail

/// Harmonic combination of two validated area-scaled internal-face halves.
///
/// Semantics:
/// - positive + positive: T_f = T_o*T_n/(T_o+T_n), computed with the
///   algebraically equivalent min/(1+min/max) form to avoid product overflow;
/// - exactly one zero_projection side: T_f=0, zero_due_to_one_half;
/// - both zero_projection sides: T_f=0, zero_due_to_both_halves;
/// - negative, non-finite, or internally inconsistent half values are rejected.
///
/// The two halves must represent the same physical face area. The result is
/// symmetric under exchanging owner and neighbour.
[[nodiscard]] inline TpfaStaticFaceTransmissibility3D
combine_internal_face_tpfa_half_transmissibilities_3d(
    const TpfaAreaScaledHalfTransmissibility3D& owner,
    const TpfaAreaScaledHalfTransmissibility3D& neighbour) {
    using namespace
        tpfa_static_face_transmissibility_3d_detail;

    validate_area_scaled_half(owner);
    validate_area_scaled_half(neighbour);
    const double face_area_m2 =
        common_face_area_m2(
            owner,
            neighbour);

    const double owner_m3 =
        owner.half_transmissibility_m3;
    const double neighbour_m3 =
        neighbour.half_transmissibility_m3;

    const bool owner_zero =
        owner_m3 == 0.0;
    const bool neighbour_zero =
        neighbour_m3 == 0.0;
    if (owner_zero || neighbour_zero) {
        return TpfaStaticFaceTransmissibility3D{
            owner_zero && neighbour_zero
                ? TpfaStaticFaceTransmissibilityDisposition3D::
                      zero_due_to_both_halves
                : TpfaStaticFaceTransmissibilityDisposition3D::
                      zero_due_to_one_half,
            face_area_m2,
            0.0};
    }

    const double smaller_m3 =
        std::min(owner_m3, neighbour_m3);
    const double larger_m3 =
        std::max(owner_m3, neighbour_m3);
    const double ratio =
        smaller_m3 / larger_m3;
    const double face_transmissibility_m3 =
        smaller_m3 /
        (1.0 + ratio);

    if (!std::isfinite(
            face_transmissibility_m3) ||
        face_transmissibility_m3 <= 0.0 ||
        face_transmissibility_m3 >
            smaller_m3) {
        throw std::invalid_argument(
            "mpmc::discretization::combine_internal_face_tpfa_half_transmissibilities_3d: positive harmonic combination must be finite, positive, and no larger than either half");
    }

    return TpfaStaticFaceTransmissibility3D{
        TpfaStaticFaceTransmissibilityDisposition3D::
            positive_harmonic_combination,
        face_area_m2,
        face_transmissibility_m3};
}

/// Construct and combine the owner/neighbour half transmissibilities for one
/// internal face. Boundary faces are rejected by the neighbour-half contract.
///
/// This wrapper still produces only a static transmissibility. It does not
/// apply pressure differences, mobility, gravity, Darcy flux, or residual terms.
[[nodiscard]] inline TpfaStaticFaceTransmissibility3D
make_internal_face_static_tpfa_transmissibility_3d(
    const mpmc::mesh::CellFaceGeometricOperator3D& geometry,
    const mpmc::mesh::CellCartesianDiagonalPermeability3D& permeability,
    mpmc::mesh::LocalIndex face) {
    if (!geometry.face_neighbour(face).has_value()) {
        throw std::invalid_argument(
            "mpmc::discretization::make_internal_face_static_tpfa_transmissibility_3d: internal face is required");
    }

    const auto owner =
        make_owner_area_scaled_tpfa_half_transmissibility_3d(
            geometry,
            permeability,
            face);
    const auto neighbour =
        make_neighbour_area_scaled_tpfa_half_transmissibility_3d(
            geometry,
            permeability,
            face);
    return combine_internal_face_tpfa_half_transmissibilities_3d(
        owner,
        neighbour);
}

} // namespace mpmc::discretization

#endif // MPMC_DISCRETIZATION_TPFA_STATIC_FACE_TRANSMISSIBILITY_3D_HPP
