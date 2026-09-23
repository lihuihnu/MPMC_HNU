#ifndef MPMC_DISCRETIZATION_PEACEMAN_WELL_INDEX_3D_HPP
#define MPMC_DISCRETIZATION_PEACEMAN_WELL_INDEX_3D_HPP

#include <mpmc/mesh/permeability_tensor_3d.hpp>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <numbers>
#include <stdexcept>
#include <string_view>

namespace mpmc::discretization {

inline constexpr std::string_view
    peaceman_well_index_3d_convention =
        "discretization/peaceman-well-index/axis-aligned-cartesian-diagonal-k/v1";

/// Axis of one straight well segment that fully penetrates one Cartesian cell.
///
/// This first contract deliberately excludes deviated wells and off-diagonal
/// permeability. The transverse plane is therefore unambiguous:
///   X-well -> YZ, Y-well -> XZ, Z-well -> XY.
enum class AxisAlignedWellDirection3D : std::uint8_t {
    x,
    y,
    z
};

/// Full Cartesian cell edge lengths in the mesh/world xyz basis [m].
struct CartesianCellDimensions3D {
    double dx_m;
    double dy_m;
    double dz_m;
};

/// Peaceman connection geometry/property result for one full-cell completion.
///
/// Scientific contract:
///   Peaceman, D.W., SPEJ 23(3), 531-543 (1983),
///   DOI 10.2118/10528-PA.
///
/// For transverse axes i,j:
///   k_h = sqrt(K_i K_j)
///   r_0 = 0.28 * sqrt(
///       dx_i^2 * sqrt(K_j/K_i) +
///       dx_j^2 * sqrt(K_i/K_j))
///       / ((K_j/K_i)^(1/4) + (K_i/K_j)^(1/4))
///   WI = 2*pi*k_h*h / (ln(r_0/r_w) + skin)
///
/// SI units:
///   K [m2], dx/h/r_0/r_w [m], WI [m3].
///
/// The residual/rate equation, mobility, BHP/rate controls and phase physics
/// are intentionally outside this geometry/property contract.
struct PeacemanWellIndex3D {
    static constexpr std::string_view convention =
        peaceman_well_index_3d_convention;

    AxisAlignedWellDirection3D well_direction{
        AxisAlignedWellDirection3D::z};
    CartesianCellDimensions3D cell_dimensions_m{};
    mpmc::mesh::CartesianDiagonalPermeabilityTensor3D
        permeability_m2{};

    double transverse_permeability_i_m2{};
    double transverse_permeability_j_m2{};
    double effective_radial_permeability_m2{};
    double completion_length_m{};
    double equivalent_radius_m{};
    double wellbore_radius_m{};
    double skin_factor{};
    double logarithmic_denominator{};
    double well_index_m3{};
};

/// Same connection result with the local cell identity made explicit.
struct CellPeacemanWellIndex3D {
    mpmc::mesh::LocalIndex cell{
        mpmc::mesh::LocalIndex::value_type{0}};
    PeacemanWellIndex3D connection;
};

namespace peaceman_well_index_3d_detail {

struct TransversePlaneData3D {
    double dimension_i_m;
    double dimension_j_m;
    double completion_length_m;
    double permeability_i_m2;
    double permeability_j_m2;
};

inline void validate_dimensions(
    const CartesianCellDimensions3D& dimensions) {
    if (!std::isfinite(dimensions.dx_m) ||
        !std::isfinite(dimensions.dy_m) ||
        !std::isfinite(dimensions.dz_m) ||
        !(dimensions.dx_m > 0.0) ||
        !(dimensions.dy_m > 0.0) ||
        !(dimensions.dz_m > 0.0)) {
        throw std::invalid_argument(
            "mpmc::discretization::make_peaceman_well_index_3d: Cartesian cell dimensions must be finite and strictly positive [m]");
    }
}

inline void validate_permeability(
    const mpmc::mesh::
        CartesianDiagonalPermeabilityTensor3D&
            permeability) {
    if (!std::isfinite(permeability.kxx_m2) ||
        !std::isfinite(permeability.kyy_m2) ||
        !std::isfinite(permeability.kzz_m2) ||
        permeability.kxx_m2 < 0.0 ||
        permeability.kyy_m2 < 0.0 ||
        permeability.kzz_m2 < 0.0) {
        throw std::invalid_argument(
            "mpmc::discretization::make_peaceman_well_index_3d: Cartesian diagonal permeability must be finite and non-negative [m2]");
    }
}

[[nodiscard]] inline TransversePlaneData3D
transverse_plane(
    const CartesianCellDimensions3D& dimensions,
    const mpmc::mesh::
        CartesianDiagonalPermeabilityTensor3D&
            permeability,
    AxisAlignedWellDirection3D direction) {
    switch (direction) {
    case AxisAlignedWellDirection3D::x:
        return {
            dimensions.dy_m,
            dimensions.dz_m,
            dimensions.dx_m,
            permeability.kyy_m2,
            permeability.kzz_m2};
    case AxisAlignedWellDirection3D::y:
        return {
            dimensions.dx_m,
            dimensions.dz_m,
            dimensions.dy_m,
            permeability.kxx_m2,
            permeability.kzz_m2};
    case AxisAlignedWellDirection3D::z:
        return {
            dimensions.dx_m,
            dimensions.dy_m,
            dimensions.dz_m,
            permeability.kxx_m2,
            permeability.kyy_m2};
    }
    throw std::invalid_argument(
        "mpmc::discretization::make_peaceman_well_index_3d: unsupported well direction");
}

} // namespace peaceman_well_index_3d_detail

/// Build the Peaceman connection factor for one full-cell, axis-aligned
/// Cartesian completion with diagonal anisotropic permeability.
///
/// The two transverse permeabilities must be strictly positive because the
/// anisotropic equivalent-radius mapping contains their ratio. The axial
/// permeability may be zero because it is not part of the classical
/// axis-aligned Peaceman connection factor.
///
/// Skin may be positive or negative, but the physical logarithmic denominator
/// must remain strictly positive. No absolute value is applied to turn an
/// invalid radius/skin combination into a positive WI.
[[nodiscard]] inline PeacemanWellIndex3D
make_peaceman_well_index_3d(
    CartesianCellDimensions3D dimensions,
    mpmc::mesh::CartesianDiagonalPermeabilityTensor3D
        permeability,
    AxisAlignedWellDirection3D direction,
    double wellbore_radius_m,
    double skin_factor = 0.0) {
    using namespace peaceman_well_index_3d_detail;

    validate_dimensions(dimensions);
    validate_permeability(permeability);
    if (!std::isfinite(wellbore_radius_m) ||
        !(wellbore_radius_m > 0.0)) {
        throw std::invalid_argument(
            "mpmc::discretization::make_peaceman_well_index_3d: wellbore radius must be finite and strictly positive [m]");
    }
    if (!std::isfinite(skin_factor)) {
        throw std::invalid_argument(
            "mpmc::discretization::make_peaceman_well_index_3d: skin factor must be finite");
    }

    const auto plane =
        transverse_plane(
            dimensions,
            permeability,
            direction);
    if (!std::isfinite(plane.permeability_i_m2) ||
        !std::isfinite(plane.permeability_j_m2) ||
        !(plane.permeability_i_m2 > 0.0) ||
        !(plane.permeability_j_m2 > 0.0)) {
        throw std::invalid_argument(
            "mpmc::discretization::make_peaceman_well_index_3d: both transverse permeabilities must be strictly positive [m2]");
    }

    const double sqrt_ki =
        std::sqrt(
            plane.permeability_i_m2);
    const double sqrt_kj =
        std::sqrt(
            plane.permeability_j_m2);
    const double effective_radial_permeability_m2 =
        sqrt_ki * sqrt_kj;
    if (!std::isfinite(
            effective_radial_permeability_m2) ||
        !(effective_radial_permeability_m2 > 0.0)) {
        throw std::invalid_argument(
            "mpmc::discretization::make_peaceman_well_index_3d: effective transverse permeability is outside representable positive range");
    }

    // Compute fourth-root anisotropy weights through sqrt(K) ratios. This is
    // algebraically identical to (Kj/Ki)^(1/4) and avoids forming K ratios
    // directly before the square root.
    const double fourth_root_kj_over_ki =
        std::sqrt(sqrt_kj / sqrt_ki);
    const double fourth_root_ki_over_kj =
        std::sqrt(sqrt_ki / sqrt_kj);
    if (!std::isfinite(
            fourth_root_kj_over_ki) ||
        !std::isfinite(
            fourth_root_ki_over_kj) ||
        !(fourth_root_kj_over_ki > 0.0) ||
        !(fourth_root_ki_over_kj > 0.0)) {
        throw std::invalid_argument(
            "mpmc::discretization::make_peaceman_well_index_3d: anisotropy ratio is outside representable positive range");
    }

    const double weighted_radius_numerator_m =
        std::hypot(
            plane.dimension_i_m *
                fourth_root_kj_over_ki,
            plane.dimension_j_m *
                fourth_root_ki_over_kj);
    const double anisotropy_denominator =
        fourth_root_kj_over_ki +
        fourth_root_ki_over_kj;
    const double equivalent_radius_m =
        0.28 *
        weighted_radius_numerator_m /
        anisotropy_denominator;
    if (!std::isfinite(equivalent_radius_m) ||
        !(equivalent_radius_m > 0.0)) {
        throw std::invalid_argument(
            "mpmc::discretization::make_peaceman_well_index_3d: equivalent well-block radius is not finite and positive");
    }

    const double logarithmic_denominator =
        std::log(equivalent_radius_m) -
        std::log(wellbore_radius_m) +
        skin_factor;
    if (!std::isfinite(logarithmic_denominator) ||
        !(logarithmic_denominator > 0.0)) {
        throw std::invalid_argument(
            "mpmc::discretization::make_peaceman_well_index_3d: ln(r0/rw)+skin must be finite and strictly positive");
    }

    const double well_index_m3 =
        2.0 *
        std::numbers::pi_v<double> *
        effective_radial_permeability_m2 *
        plane.completion_length_m /
        logarithmic_denominator;
    if (!std::isfinite(well_index_m3) ||
        !(well_index_m3 > 0.0)) {
        throw std::invalid_argument(
            "mpmc::discretization::make_peaceman_well_index_3d: well index must be finite and strictly positive [m3]");
    }

    return PeacemanWellIndex3D{
        direction,
        dimensions,
        permeability,
        plane.permeability_i_m2,
        plane.permeability_j_m2,
        effective_radial_permeability_m2,
        plane.completion_length_m,
        equivalent_radius_m,
        wellbore_radius_m,
        skin_factor,
        logarithmic_denominator,
        well_index_m3};
}

/// Cell-connected overload using the repository's typed diagonal permeability
/// snapshot. Geometry remains explicit because this slice intentionally does
/// not guess Cartesian dimensions from a general corner-point cell.
[[nodiscard]] inline CellPeacemanWellIndex3D
make_cell_peaceman_well_index_3d(
    const mpmc::mesh::CellCartesianDiagonalPermeability3D&
        permeability,
    mpmc::mesh::LocalIndex cell,
    CartesianCellDimensions3D dimensions,
    AxisAlignedWellDirection3D direction,
    double wellbore_radius_m,
    double skin_factor = 0.0) {
    return {
        cell,
        make_peaceman_well_index_3d(
            dimensions,
            permeability.tensor(cell),
            direction,
            wellbore_radius_m,
            skin_factor)};
}

} // namespace mpmc::discretization

#endif // MPMC_DISCRETIZATION_PEACEMAN_WELL_INDEX_3D_HPP
