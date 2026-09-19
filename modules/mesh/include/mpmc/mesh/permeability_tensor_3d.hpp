#ifndef MPMC_MESH_PERMEABILITY_TENSOR_3D_HPP
#define MPMC_MESH_PERMEABILITY_TENSOR_3D_HPP

#include <mpmc/mesh/cell_face_geometric_operator_3d.hpp>
#include <mpmc/mesh/dense_field.hpp>
#include <mpmc/mesh/topology.hpp>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <optional>
#include <stdexcept>
#include <string_view>
#include <utility>
#include <vector>

namespace mpmc::mesh {

/// Cell permeability tensor for the current baseline.
///
/// PERMX/PERMY/PERMZ are interpreted as Kxx/Kyy/Kzz in the existing mesh/world
/// Cartesian x/y/z basis. No local-axis rotation or off-diagonal permeability
/// components are inferred by this contract.
struct CartesianDiagonalPermeabilityTensor3D {
    double kxx_m2;
    double kyy_m2;
    double kzz_m2;
};

/// Immutable cell-aligned diagonal permeability tensor snapshot.
class CellCartesianDiagonalPermeability3D {
public:
    explicit CellCartesianDiagonalPermeability3D(
        std::vector<CartesianDiagonalPermeabilityTensor3D> tensors)
        : tensors_(std::move(tensors)) {
        for (const auto tensor : tensors_) {
            if (!std::isfinite(tensor.kxx_m2) ||
                !std::isfinite(tensor.kyy_m2) ||
                !std::isfinite(tensor.kzz_m2) ||
                tensor.kxx_m2 < 0.0 ||
                tensor.kyy_m2 < 0.0 ||
                tensor.kzz_m2 < 0.0) {
                throw std::invalid_argument(
                    "mpmc::mesh::CellCartesianDiagonalPermeability3D: tensor components must be finite and non-negative");
            }
        }
    }

    CellCartesianDiagonalPermeability3D(
        const CellCartesianDiagonalPermeability3D&) = default;
    CellCartesianDiagonalPermeability3D(
        CellCartesianDiagonalPermeability3D&&) noexcept = default;
    CellCartesianDiagonalPermeability3D& operator=(
        const CellCartesianDiagonalPermeability3D&) = delete;
    CellCartesianDiagonalPermeability3D& operator=(
        CellCartesianDiagonalPermeability3D&&) = delete;
    ~CellCartesianDiagonalPermeability3D() = default;

    [[nodiscard]] std::size_t cell_count() const noexcept {
        return tensors_.size();
    }

    [[nodiscard]] CartesianDiagonalPermeabilityTensor3D
    tensor(LocalIndex cell) const {
        return tensors_.at(
            static_cast<std::size_t>(cell.value()));
    }

private:
    std::vector<CartesianDiagonalPermeabilityTensor3D> tensors_;
};

namespace permeability_tensor_3d_detail {

inline void require_permeability_field(
    const DenseFieldSnapshot& field,
    const Topology& topology,
    std::string_view expected_id) {
    if (field.location() != EntityKind::cell ||
        field.component_count() != 1U ||
        field.entity_count() !=
            topology.entity_count(EntityKind::cell)) {
        throw std::invalid_argument(
            "mpmc::mesh::make_cell_cartesian_diagonal_permeability_3d: permeability fields must be scalar cell fields aligned to topology");
    }
    if (field.metadata().id != expected_id ||
        field.metadata().unit != "m2") {
        throw std::invalid_argument(
            "mpmc::mesh::make_cell_cartesian_diagonal_permeability_3d: expected PERMX/PERMY/PERMZ with SI unit m2");
    }
}

struct Vector3D {
    double x;
    double y;
    double z;
};

[[nodiscard]] inline Vector3D apply(
    CartesianDiagonalPermeabilityTensor3D tensor,
    Displacement3D vector) {
    return Vector3D{
        tensor.kxx_m2 * vector.x_m,
        tensor.kyy_m2 * vector.y_m,
        tensor.kzz_m2 * vector.z_m};
}

[[nodiscard]] inline Vector3D apply(
    CartesianDiagonalPermeabilityTensor3D tensor,
    UnitVector3D vector) {
    return Vector3D{
        tensor.kxx_m2 * vector.x,
        tensor.kyy_m2 * vector.y,
        tensor.kzz_m2 * vector.z};
}

[[nodiscard]] inline double dot(
    Vector3D left,
    Displacement3D right) {
    return left.x * right.x_m +
           left.y * right.y_m +
           left.z * right.z_m;
}

[[nodiscard]] inline double dot(
    Vector3D left,
    UnitVector3D right) {
    return left.x * right.x +
           left.y * right.y +
           left.z * right.z;
}

[[nodiscard]] inline double magnitude(Vector3D value) {
    return std::sqrt(
        value.x * value.x +
        value.y * value.y +
        value.z * value.z);
}

[[nodiscard]] inline double magnitude(Displacement3D value) {
    return std::sqrt(
        value.x_m * value.x_m +
        value.y_m * value.y_m +
        value.z_m * value.z_m);
}

[[nodiscard]] inline double maximum_component(
    CartesianDiagonalPermeabilityTensor3D tensor) {
    return std::max(
        {tensor.kxx_m2,
         tensor.kyy_m2,
         tensor.kzz_m2});
}

} // namespace permeability_tensor_3d_detail

[[nodiscard]] inline CellCartesianDiagonalPermeability3D
make_cell_cartesian_diagonal_permeability_3d(
    const Topology& topology,
    const DenseFieldSnapshot& permx,
    const DenseFieldSnapshot& permy,
    const DenseFieldSnapshot& permz) {
    permeability_tensor_3d_detail::require_permeability_field(
        permx, topology, "PERMX");
    permeability_tensor_3d_detail::require_permeability_field(
        permy, topology, "PERMY");
    permeability_tensor_3d_detail::require_permeability_field(
        permz, topology, "PERMZ");

    const std::size_t cell_count =
        topology.entity_count(EntityKind::cell);
    std::vector<CartesianDiagonalPermeabilityTensor3D>
        tensors;
    tensors.reserve(cell_count);
    for (std::size_t cell = 0U;
         cell < cell_count;
         ++cell) {
        const LocalIndex local{
            static_cast<LocalIndex::value_type>(cell)};
        tensors.push_back(
            CartesianDiagonalPermeabilityTensor3D{
                permx.value(local, 0U),
                permy.value(local, 0U),
                permz.value(local, 0U)});
    }

    return CellCartesianDiagonalPermeability3D{
        std::move(tensors)};
}

/// Direction alignment diagnostic. This is dimensionless and compares two
/// same-direction vectors after normalization.
struct DirectionAlignment3D {
    double alignment_cosine;
    double angle_rad;
};

/// Requested centre-line diagnostic: K*d_cc versus owner-relative face normal.
///
/// This is exposed for auditability but is not used as the mathematical
/// K-orthogonality criterion. For anisotropic K, K*d_cc || n is generally not
/// equivalent to the standard half-face co-normal condition.
struct CenterLinePermeabilityDirectionDiagnostic3D {
    std::optional<DirectionAlignment3D> owner;
    std::optional<DirectionAlignment3D> neighbour;
};

enum class KOrthogonalityDisposition3D {
    k_orthogonal_within_policy,
    requires_k_non_orthogonal_treatment,
    degenerate_permeability_direction
};

struct KOrthogonalityAdmissibilityPolicy3D {
    double max_half_face_co_normal_angle_rad;
};

/// Per-half-face K-orthogonality result using the standard co-normal condition:
/// K*n_out must align with the cell-centre-to-face displacement.
struct InternalFaceKOrthogonality3D {
    KOrthogonalityDisposition3D disposition;
    std::optional<DirectionAlignment3D> owner_half_face;
    std::optional<DirectionAlignment3D> neighbour_half_face;
    CenterLinePermeabilityDirectionDiagnostic3D
        center_line_diagnostic;
    double max_half_face_co_normal_angle_rad;
};

namespace permeability_tensor_3d_detail {

[[nodiscard]] inline std::optional<DirectionAlignment3D>
alignment(
    Vector3D left,
    Displacement3D right,
    double scale) {
    const double left_magnitude =
        magnitude(left);
    const double right_magnitude =
        magnitude(right);
    const double tolerance =
        4096.0 *
        std::numeric_limits<double>::epsilon() *
        scale;
    if (!std::isfinite(left_magnitude) ||
        !std::isfinite(right_magnitude) ||
        left_magnitude <= tolerance ||
        right_magnitude <= 0.0) {
        return std::nullopt;
    }

    double cosine =
        dot(left, right) /
        (left_magnitude * right_magnitude);
    if (!std::isfinite(cosine)) {
        return std::nullopt;
    }
    cosine = std::clamp(cosine, -1.0, 1.0);
    return DirectionAlignment3D{
        cosine,
        std::acos(cosine)};
}

[[nodiscard]] inline std::optional<DirectionAlignment3D>
alignment(
    Vector3D left,
    UnitVector3D right,
    double scale) {
    const double left_magnitude =
        magnitude(left);
    const double tolerance =
        4096.0 *
        std::numeric_limits<double>::epsilon() *
        scale;
    if (!std::isfinite(left_magnitude) ||
        left_magnitude <= tolerance) {
        return std::nullopt;
    }

    double cosine =
        dot(left, right) /
        left_magnitude;
    if (!std::isfinite(cosine)) {
        return std::nullopt;
    }
    cosine = std::clamp(cosine, -1.0, 1.0);
    return DirectionAlignment3D{
        cosine,
        std::acos(cosine)};
}

} // namespace permeability_tensor_3d_detail

[[nodiscard]] inline InternalFaceKOrthogonality3D
classify_internal_face_k_orthogonality(
    const CellFaceGeometricOperator3D& geometry,
    const CellCartesianDiagonalPermeability3D& permeability,
    LocalIndex face,
    KOrthogonalityAdmissibilityPolicy3D policy) {
    using namespace permeability_tensor_3d_detail;

    const double half_pi =
        0.5 * std::acos(-1.0);
    if (!std::isfinite(
            policy.max_half_face_co_normal_angle_rad) ||
        policy.max_half_face_co_normal_angle_rad < 0.0 ||
        policy.max_half_face_co_normal_angle_rad >=
            half_pi) {
        throw std::invalid_argument(
            "mpmc::mesh::classify_internal_face_k_orthogonality: maximum co-normal angle must be finite in [0, pi/2)");
    }
    if (permeability.cell_count() !=
        geometry.cell_count()) {
        throw std::invalid_argument(
            "mpmc::mesh::classify_internal_face_k_orthogonality: permeability/geometry cell counts do not match");
    }

    const auto neighbour =
        geometry.face_neighbour(face);
    const auto non_orthogonality =
        geometry.internal_non_orthogonality(face);
    const auto neighbour_to_face =
        geometry.neighbour_to_face_displacement_m(face);
    if (!neighbour.has_value() ||
        !non_orthogonality.has_value() ||
        !neighbour_to_face.has_value()) {
        throw std::invalid_argument(
            "mpmc::mesh::classify_internal_face_k_orthogonality: internal face geometry is required");
    }

    const LocalIndex owner =
        geometry.face_owner(face);
    const auto owner_tensor =
        permeability.tensor(owner);
    const auto neighbour_tensor =
        permeability.tensor(*neighbour);
    const auto owner_normal =
        geometry.transmissibility_geometry(face)
            .owner_unit_normal;
    const UnitVector3D neighbour_normal{
        -owner_normal.x,
        -owner_normal.y,
        -owner_normal.z};

    const auto owner_to_face =
        geometry.owner_to_face_displacement_m(face);
    const auto owner_conormal =
        apply(owner_tensor, owner_normal);
    const auto neighbour_conormal =
        apply(neighbour_tensor, neighbour_normal);

    // K*n has permeability units (m2), so its degeneracy tolerance must
    // be scaled only by a permeability magnitude. Do not multiply by length.
    const double owner_scale =
        maximum_component(owner_tensor);
    const double neighbour_scale =
        maximum_component(neighbour_tensor);

    const auto owner_half_face =
        alignment(
            owner_conormal,
            owner_to_face,
            owner_scale);
    const auto neighbour_half_face =
        alignment(
            neighbour_conormal,
            *neighbour_to_face,
            neighbour_scale);

    const auto center_vector =
        non_orthogonality
            ->owner_to_neighbour_displacement_m;
    const auto owner_k_center =
        apply(owner_tensor, center_vector);
    const auto neighbour_k_center =
        apply(neighbour_tensor, center_vector);
    const auto owner_center_alignment =
        alignment(
            owner_k_center,
            owner_normal,
            maximum_component(owner_tensor) *
                non_orthogonality->center_distance_m);
    const auto neighbour_center_alignment =
        alignment(
            neighbour_k_center,
            owner_normal,
            maximum_component(neighbour_tensor) *
                non_orthogonality->center_distance_m);

    KOrthogonalityDisposition3D disposition =
        KOrthogonalityDisposition3D::
            degenerate_permeability_direction;
    if (owner_half_face.has_value() &&
        neighbour_half_face.has_value()) {
        const bool aligned =
            owner_half_face->alignment_cosine > 0.0 &&
            neighbour_half_face->alignment_cosine > 0.0 &&
            owner_half_face->angle_rad <=
                policy.max_half_face_co_normal_angle_rad &&
            neighbour_half_face->angle_rad <=
                policy.max_half_face_co_normal_angle_rad;
        disposition =
            aligned
                ? KOrthogonalityDisposition3D::
                      k_orthogonal_within_policy
                : KOrthogonalityDisposition3D::
                      requires_k_non_orthogonal_treatment;
    }

    return InternalFaceKOrthogonality3D{
        disposition,
        owner_half_face,
        neighbour_half_face,
        CenterLinePermeabilityDirectionDiagnostic3D{
            owner_center_alignment,
            neighbour_center_alignment},
        policy.max_half_face_co_normal_angle_rad};
}

enum class CombinedTransmissibilityAdmissibilityDisposition3D {
    direct_normal_projection_k_orthogonal_candidate,
    requires_geometry_non_orthogonal_treatment,
    requires_k_non_orthogonal_treatment,
    requires_geometry_and_k_non_orthogonal_treatment,
    degenerate_permeability_direction
};

struct CombinedTransmissibilityAdmissibility3D {
    CombinedTransmissibilityAdmissibilityDisposition3D
        disposition;
    TransmissibilityGeometryAdmissibility3D geometry;
    InternalFaceKOrthogonality3D k_orthogonality;
};

[[nodiscard]] inline CombinedTransmissibilityAdmissibility3D
classify_internal_face_transmissibility_admissibility(
    const CellFaceGeometricOperator3D& geometry,
    const CellCartesianDiagonalPermeability3D& permeability,
    LocalIndex face,
    TransmissibilityGeometryAdmissibilityPolicy3D
        geometry_policy,
    KOrthogonalityAdmissibilityPolicy3D k_policy) {
    const auto geometry_result =
        classify_internal_face_transmissibility_geometry(
            geometry,
            face,
            geometry_policy);
    const auto k_result =
        classify_internal_face_k_orthogonality(
            geometry,
            permeability,
            face,
            k_policy);

    CombinedTransmissibilityAdmissibilityDisposition3D
        disposition =
            CombinedTransmissibilityAdmissibilityDisposition3D::
                degenerate_permeability_direction;

    if (k_result.disposition !=
        KOrthogonalityDisposition3D::
            degenerate_permeability_direction) {
        const bool geometry_direct =
            geometry_result.disposition ==
            TransmissibilityGeometryDisposition3D::
                direct_normal_projection_allowed;
        const bool k_direct =
            k_result.disposition ==
            KOrthogonalityDisposition3D::
                k_orthogonal_within_policy;

        if (geometry_direct && k_direct) {
            disposition =
                CombinedTransmissibilityAdmissibilityDisposition3D::
                    direct_normal_projection_k_orthogonal_candidate;
        } else if (!geometry_direct && k_direct) {
            disposition =
                CombinedTransmissibilityAdmissibilityDisposition3D::
                    requires_geometry_non_orthogonal_treatment;
        } else if (geometry_direct && !k_direct) {
            disposition =
                CombinedTransmissibilityAdmissibilityDisposition3D::
                    requires_k_non_orthogonal_treatment;
        } else {
            disposition =
                CombinedTransmissibilityAdmissibilityDisposition3D::
                    requires_geometry_and_k_non_orthogonal_treatment;
        }
    }

    return CombinedTransmissibilityAdmissibility3D{
        disposition,
        geometry_result,
        k_result};
}

} // namespace mpmc::mesh

#endif // MPMC_MESH_PERMEABILITY_TENSOR_3D_HPP
