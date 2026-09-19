#ifndef MPMC_MESH_CARTESIAN_SYMMETRIC_TENSOR_3D_HPP
#define MPMC_MESH_CARTESIAN_SYMMETRIC_TENSOR_3D_HPP

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <limits>

namespace mpmc::mesh {

/// Explicit basis for symmetric 3D material tensors attached to mesh entities.
enum class CartesianTensorBasis3D : std::uint8_t {
    mesh_world_xyz = 0,
};

/// Property-neutral six-component representation of a symmetric 3D tensor.
///
/// Component order is [xx, yy, zz, xy, xz, yz]. The mirrored components
/// yx/zx/zy are intentionally not duplicated.
struct CartesianSymmetricTensorComponents3D {
    double xx;
    double yy;
    double zz;
    double xy;
    double xz;
    double yz;
};

enum class PositiveSemidefiniteStatus3D : std::uint8_t {
    valid = 0,
    non_finite = 1,
    negative_diagonal = 2,
    indefinite = 3,
};

/// Scale-aware positive-semidefinite classification for a real symmetric 3x3
/// tensor. Zero eigen-directions are allowed.
///
/// The test uses the diagonal entries, all 2x2 principal minors and the 3x3
/// determinant after normalization by the largest absolute component. This is
/// the shared numerical validity rule for symmetric permeability/conductivity
/// contracts; it does not modify or clip the supplied tensor.
[[nodiscard]] inline PositiveSemidefiniteStatus3D
classify_positive_semidefinite(
    CartesianSymmetricTensorComponents3D tensor) noexcept {
    const double scale = std::max(
        {std::abs(tensor.xx),
         std::abs(tensor.yy),
         std::abs(tensor.zz),
         std::abs(tensor.xy),
         std::abs(tensor.xz),
         std::abs(tensor.yz)});

    if (!std::isfinite(scale)) {
        return PositiveSemidefiniteStatus3D::non_finite;
    }
    if (tensor.xx < 0.0 ||
        tensor.yy < 0.0 ||
        tensor.zz < 0.0) {
        return PositiveSemidefiniteStatus3D::negative_diagonal;
    }
    if (scale == 0.0) {
        return PositiveSemidefiniteStatus3D::valid;
    }

    const double a = tensor.xx / scale;
    const double b = tensor.yy / scale;
    const double c = tensor.zz / scale;
    const double d = tensor.xy / scale;
    const double e = tensor.xz / scale;
    const double f = tensor.yz / scale;

    const double minor_xy = a * b - d * d;
    const double minor_xz = a * c - e * e;
    const double minor_yz = b * c - f * f;
    const double determinant =
        a * b * c +
        2.0 * d * e * f -
        a * f * f -
        b * e * e -
        c * d * d;

    constexpr double tolerance =
        4096.0 *
        std::numeric_limits<double>::epsilon();
    if (!std::isfinite(minor_xy) ||
        !std::isfinite(minor_xz) ||
        !std::isfinite(minor_yz) ||
        !std::isfinite(determinant) ||
        minor_xy < -tolerance ||
        minor_xz < -tolerance ||
        minor_yz < -tolerance ||
        determinant < -tolerance) {
        return PositiveSemidefiniteStatus3D::indefinite;
    }
    return PositiveSemidefiniteStatus3D::valid;
}

} // namespace mpmc::mesh

#endif // MPMC_MESH_CARTESIAN_SYMMETRIC_TENSOR_3D_HPP
