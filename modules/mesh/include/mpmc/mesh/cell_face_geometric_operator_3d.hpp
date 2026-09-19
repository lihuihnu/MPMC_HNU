#ifndef MPMC_MESH_CELL_FACE_GEOMETRIC_OPERATOR_3D_HPP
#define MPMC_MESH_CELL_FACE_GEOMETRIC_OPERATOR_3D_HPP

#include <mpmc/mesh/corner_point_geometry_3d.hpp>
#include <mpmc/mesh/face_geometry_3d.hpp>
#include <mpmc/mesh/topology.hpp>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <optional>
#include <span>
#include <stdexcept>
#include <utility>
#include <vector>

namespace mpmc::mesh {

/// Vector from a cell centroid to a face centroid, in SI metres.
struct Displacement3D {
    double x_m;
    double y_m;
    double z_m;
};

/// Internal-face non-orthogonality measured between the owner-to-neighbour
/// cell-centre vector and the owner-relative face unit normal.
///
/// theta=0 is orthogonal. Increasing theta indicates increasing non-orthogonality.
struct InternalFaceNonOrthogonality3D {
    Displacement3D owner_to_neighbour_displacement_m;
    double center_distance_m;
    double normal_alignment_cosine;
    double angle_rad;
};

/// Geometry-only inputs that a later two-point transmissibility kernel may consume.
///
/// This contract deliberately contains no permeability, mobility, Darcy flux,
/// residual, or transmissibility value. The owner normal points away from the
/// owner cell. For an internal face both normal distances are strictly positive.
struct FaceTransmissibilityGeometry3D {
    LocalIndex owner;
    std::optional<LocalIndex> neighbour;
    double area_m2;
    UnitVector3D owner_unit_normal;
    double owner_normal_distance_m;
    std::optional<double> neighbour_normal_distance_m;
    std::optional<InternalFaceNonOrthogonality3D>
        internal_non_orthogonality;
};

/// Immutable 3D cell/face geometric-operator snapshot.
///
/// The snapshot is aligned to the supplied Topology local cell/face ordering.
/// The current corner-point baseline defines each cell centroid as the
/// arithmetic mean of its eight processed shared vertices. This is the same
/// cell reference point used by the active corner-point processor to orient
/// FaceGeometry3D owner normals; no alternative volumetric-centroid convention
/// is introduced here.
///
/// Owner displacement is (face centroid - owner cell centroid). For an internal
/// face, neighbour displacement is (face centroid - neighbour cell centroid).
/// With the owner-relative unit normal n:
///
///   d_owner     = dot(owner_to_face, n)
///   d_neighbour = -dot(neighbour_to_face, n)
///
/// Both distances must be finite and strictly positive above a scale-aware
/// floating-point tolerance. For an internal face the owner-to-neighbour
/// centre vector d_cc is also retained and the non-orthogonality angle is
///
///   theta = acos(dot(d_cc, n) / |d_cc|)
///
/// where theta=0 denotes an orthogonal face. Boundary faces have no neighbour
/// displacement, neighbour distance, or non-orthogonality quantity.
class CellFaceGeometricOperator3D {
public:
    CellFaceGeometricOperator3D(
        std::vector<Coordinate3D> cell_centroids_m,
        std::vector<LocalIndex> face_owners,
        std::vector<std::optional<LocalIndex>> face_neighbours,
        std::vector<Displacement3D> owner_to_face_displacements_m,
        std::vector<std::optional<Displacement3D>>
            neighbour_to_face_displacements_m,
        std::vector<double> owner_normal_distances_m,
        std::vector<std::optional<double>> neighbour_normal_distances_m,
        std::vector<std::optional<InternalFaceNonOrthogonality3D>>
            internal_non_orthogonality,
        std::vector<double> face_areas_m2,
        std::vector<UnitVector3D> face_owner_unit_normals)
        : cell_centroids_m_(std::move(cell_centroids_m)),
          face_owners_(std::move(face_owners)),
          face_neighbours_(std::move(face_neighbours)),
          owner_to_face_displacements_m_(
              std::move(owner_to_face_displacements_m)),
          neighbour_to_face_displacements_m_(
              std::move(neighbour_to_face_displacements_m)),
          owner_normal_distances_m_(
              std::move(owner_normal_distances_m)),
          neighbour_normal_distances_m_(
              std::move(neighbour_normal_distances_m)),
          internal_non_orthogonality_(
              std::move(internal_non_orthogonality)),
          face_areas_m2_(std::move(face_areas_m2)),
          face_owner_unit_normals_(
              std::move(face_owner_unit_normals)) {
        validate();
    }

    CellFaceGeometricOperator3D(
        const CellFaceGeometricOperator3D&) = default;
    CellFaceGeometricOperator3D(
        CellFaceGeometricOperator3D&&) noexcept = default;
    CellFaceGeometricOperator3D& operator=(
        const CellFaceGeometricOperator3D&) = delete;
    CellFaceGeometricOperator3D& operator=(
        CellFaceGeometricOperator3D&&) = delete;
    ~CellFaceGeometricOperator3D() = default;

    [[nodiscard]] std::size_t cell_count() const noexcept {
        return cell_centroids_m_.size();
    }

    [[nodiscard]] std::size_t face_count() const noexcept {
        return face_owners_.size();
    }

    [[nodiscard]] Coordinate3D
    cell_centroid_m(LocalIndex cell) const {
        return cell_centroids_m_.at(
            static_cast<std::size_t>(cell.value()));
    }

    [[nodiscard]] LocalIndex
    face_owner(LocalIndex face) const {
        return face_owners_.at(
            static_cast<std::size_t>(face.value()));
    }

    [[nodiscard]] std::optional<LocalIndex>
    face_neighbour(LocalIndex face) const {
        return face_neighbours_.at(
            static_cast<std::size_t>(face.value()));
    }

    [[nodiscard]] Displacement3D
    owner_to_face_displacement_m(LocalIndex face) const {
        return owner_to_face_displacements_m_.at(
            static_cast<std::size_t>(face.value()));
    }

    [[nodiscard]] std::optional<Displacement3D>
    neighbour_to_face_displacement_m(LocalIndex face) const {
        return neighbour_to_face_displacements_m_.at(
            static_cast<std::size_t>(face.value()));
    }

    [[nodiscard]] double
    owner_normal_distance_m(LocalIndex face) const {
        return owner_normal_distances_m_.at(
            static_cast<std::size_t>(face.value()));
    }

    [[nodiscard]] std::optional<double>
    neighbour_normal_distance_m(LocalIndex face) const {
        return neighbour_normal_distances_m_.at(
            static_cast<std::size_t>(face.value()));
    }

    [[nodiscard]] std::optional<InternalFaceNonOrthogonality3D>
    internal_non_orthogonality(LocalIndex face) const {
        return internal_non_orthogonality_.at(
            static_cast<std::size_t>(face.value()));
    }

    [[nodiscard]] FaceTransmissibilityGeometry3D
    transmissibility_geometry(LocalIndex face) const {
        const std::size_t local =
            static_cast<std::size_t>(face.value());
        return FaceTransmissibilityGeometry3D{
            face_owners_.at(local),
            face_neighbours_.at(local),
            face_areas_m2_.at(local),
            face_owner_unit_normals_.at(local),
            owner_normal_distances_m_.at(local),
            neighbour_normal_distances_m_.at(local),
            internal_non_orthogonality_.at(local)};
    }

private:
    static void require_finite_coordinate(
        Coordinate3D value,
        const char* message) {
        if (!std::isfinite(value.x_m) ||
            !std::isfinite(value.y_m) ||
            !std::isfinite(value.z_m)) {
            throw std::invalid_argument(message);
        }
    }

    static void require_finite_displacement(
        Displacement3D value,
        const char* message) {
        if (!std::isfinite(value.x_m) ||
            !std::isfinite(value.y_m) ||
            !std::isfinite(value.z_m)) {
            throw std::invalid_argument(message);
        }
    }

    void validate() const {
        const std::size_t faces = face_owners_.size();
        if (face_neighbours_.size() != faces ||
            owner_to_face_displacements_m_.size() != faces ||
            neighbour_to_face_displacements_m_.size() != faces ||
            owner_normal_distances_m_.size() != faces ||
            neighbour_normal_distances_m_.size() != faces ||
            internal_non_orthogonality_.size() != faces ||
            face_areas_m2_.size() != faces ||
            face_owner_unit_normals_.size() != faces) {
            throw std::invalid_argument(
                "mpmc::mesh::CellFaceGeometricOperator3D: face array sizes do not match");
        }

        for (const auto centroid : cell_centroids_m_) {
            require_finite_coordinate(
                centroid,
                "mpmc::mesh::CellFaceGeometricOperator3D: non-finite cell centroid");
        }

        constexpr double normal_tolerance =
            128.0 *
            std::numeric_limits<double>::epsilon();

        for (std::size_t face = 0U;
             face < faces;
             ++face) {
            const std::size_t owner =
                static_cast<std::size_t>(
                    face_owners_[face].value());
            if (owner >= cell_centroids_m_.size()) {
                throw std::out_of_range(
                    "mpmc::mesh::CellFaceGeometricOperator3D: face owner index out of range");
            }

            const bool has_neighbour =
                face_neighbours_[face].has_value();
            if (has_neighbour !=
                    neighbour_to_face_displacements_m_[face].has_value() ||
                has_neighbour !=
                    neighbour_normal_distances_m_[face].has_value() ||
                has_neighbour !=
                    internal_non_orthogonality_[face].has_value()) {
                throw std::invalid_argument(
                    "mpmc::mesh::CellFaceGeometricOperator3D: neighbour arrays are inconsistent");
            }
            if (has_neighbour) {
                const std::size_t neighbour =
                    static_cast<std::size_t>(
                        face_neighbours_[face]->value());
                if (neighbour >= cell_centroids_m_.size() ||
                    neighbour == owner) {
                    throw std::out_of_range(
                        "mpmc::mesh::CellFaceGeometricOperator3D: neighbour index is invalid");
                }
            }

            require_finite_displacement(
                owner_to_face_displacements_m_[face],
                "mpmc::mesh::CellFaceGeometricOperator3D: non-finite owner displacement");
            if (has_neighbour) {
                require_finite_displacement(
                    *neighbour_to_face_displacements_m_[face],
                    "mpmc::mesh::CellFaceGeometricOperator3D: non-finite neighbour displacement");
            }

            const double owner_distance =
                owner_normal_distances_m_[face];
            if (!std::isfinite(owner_distance) ||
                owner_distance <= 0.0) {
                throw std::invalid_argument(
                    "mpmc::mesh::CellFaceGeometricOperator3D: owner normal distance must be finite and positive");
            }
            if (has_neighbour) {
                const double neighbour_distance =
                    *neighbour_normal_distances_m_[face];
                if (!std::isfinite(neighbour_distance) ||
                    neighbour_distance <= 0.0) {
                    throw std::invalid_argument(
                        "mpmc::mesh::CellFaceGeometricOperator3D: neighbour normal distance must be finite and positive");
                }
            }

            if (has_neighbour) {
                const auto non_orthogonality =
                    *internal_non_orthogonality_[face];
                require_finite_displacement(
                    non_orthogonality
                        .owner_to_neighbour_displacement_m,
                    "mpmc::mesh::CellFaceGeometricOperator3D: non-finite owner-to-neighbour displacement");
                if (!std::isfinite(
                        non_orthogonality.center_distance_m) ||
                    non_orthogonality.center_distance_m <= 0.0 ||
                    !std::isfinite(
                        non_orthogonality.normal_alignment_cosine) ||
                    non_orthogonality.normal_alignment_cosine <= 0.0 ||
                    non_orthogonality.normal_alignment_cosine >
                        1.0 + normal_tolerance ||
                    !std::isfinite(
                        non_orthogonality.angle_rad) ||
                    non_orthogonality.angle_rad < 0.0 ||
                    non_orthogonality.angle_rad >
                        0.5 * std::acos(-1.0)) {
                    throw std::invalid_argument(
                        "mpmc::mesh::CellFaceGeometricOperator3D: invalid internal-face non-orthogonality geometry");
                }
            }

            const double area = face_areas_m2_[face];
            if (!std::isfinite(area) ||
                area <= 0.0) {
                throw std::invalid_argument(
                    "mpmc::mesh::CellFaceGeometricOperator3D: face area must be finite and positive");
            }

            const auto normal =
                face_owner_unit_normals_[face];
            const double normal_magnitude =
                std::sqrt(
                    normal.x * normal.x +
                    normal.y * normal.y +
                    normal.z * normal.z);
            if (!std::isfinite(normal.x) ||
                !std::isfinite(normal.y) ||
                !std::isfinite(normal.z) ||
                !std::isfinite(normal_magnitude) ||
                std::abs(normal_magnitude - 1.0) >
                    normal_tolerance) {
                throw std::invalid_argument(
                    "mpmc::mesh::CellFaceGeometricOperator3D: face normal is not finite unit length");
            }
        }
    }

    std::vector<Coordinate3D> cell_centroids_m_;
    std::vector<LocalIndex> face_owners_;
    std::vector<std::optional<LocalIndex>> face_neighbours_;
    std::vector<Displacement3D>
        owner_to_face_displacements_m_;
    std::vector<std::optional<Displacement3D>>
        neighbour_to_face_displacements_m_;
    std::vector<double> owner_normal_distances_m_;
    std::vector<std::optional<double>>
        neighbour_normal_distances_m_;
    std::vector<std::optional<InternalFaceNonOrthogonality3D>>
        internal_non_orthogonality_;
    std::vector<double> face_areas_m2_;
    std::vector<UnitVector3D>
        face_owner_unit_normals_;
};

namespace cell_face_geometric_operator_3d_detail {

[[nodiscard]] inline Displacement3D subtract(
    Coordinate3D left,
    Coordinate3D right) {
    return Displacement3D{
        left.x_m - right.x_m,
        left.y_m - right.y_m,
        left.z_m - right.z_m};
}

[[nodiscard]] inline double dot(
    Displacement3D value,
    UnitVector3D normal) {
    return value.x_m * normal.x +
           value.y_m * normal.y +
           value.z_m * normal.z;
}

[[nodiscard]] inline double magnitude(
    Displacement3D value) {
    return std::sqrt(
        value.x_m * value.x_m +
        value.y_m * value.y_m +
        value.z_m * value.z_m);
}

[[nodiscard]] inline Displacement3D
owner_to_neighbour(
    Displacement3D owner_to_face,
    Displacement3D neighbour_to_face) {
    return Displacement3D{
        owner_to_face.x_m - neighbour_to_face.x_m,
        owner_to_face.y_m - neighbour_to_face.y_m,
        owner_to_face.z_m - neighbour_to_face.z_m};
}

[[nodiscard]] inline Coordinate3D cell_vertex_mean(
    const Topology& topology,
    std::span<const Coordinate3D> vertex_coordinates_m,
    LocalIndex cell) {
    const auto vertices =
        topology.relation(
            EntityKind::cell,
            EntityKind::vertex)
            .adjacent(cell);
    if (vertices.size() != 8U) {
        throw std::invalid_argument(
            "mpmc::mesh::make_cell_face_geometric_operator_3d: each 3D cell must have exactly eight vertices");
    }

    std::vector<LocalIndex> unique(
        vertices.begin(),
        vertices.end());
    std::sort(
        unique.begin(),
        unique.end(),
        [](LocalIndex left, LocalIndex right) {
            return left.value() < right.value();
        });
    if (std::adjacent_find(
            unique.begin(),
            unique.end()) != unique.end()) {
        throw std::invalid_argument(
            "mpmc::mesh::make_cell_face_geometric_operator_3d: cell has duplicate vertices");
    }

    Coordinate3D centroid{0.0, 0.0, 0.0};
    for (const auto vertex : vertices) {
        const std::size_t local =
            static_cast<std::size_t>(
                vertex.value());
        if (local >= vertex_coordinates_m.size()) {
            throw std::out_of_range(
                "mpmc::mesh::make_cell_face_geometric_operator_3d: cell vertex index out of range");
        }
        const auto coordinate =
            vertex_coordinates_m[local];
        if (!std::isfinite(coordinate.x_m) ||
            !std::isfinite(coordinate.y_m) ||
            !std::isfinite(coordinate.z_m)) {
            throw std::invalid_argument(
                "mpmc::mesh::make_cell_face_geometric_operator_3d: non-finite vertex coordinate");
        }
        centroid.x_m += coordinate.x_m;
        centroid.y_m += coordinate.y_m;
        centroid.z_m += coordinate.z_m;
    }
    centroid.x_m /= 8.0;
    centroid.y_m /= 8.0;
    centroid.z_m /= 8.0;
    return centroid;
}

inline void require_positive_normal_distance(
    double distance,
    double scale,
    const char* message) {
    if (!std::isfinite(distance) ||
        !std::isfinite(scale) ||
        scale <= 0.0) {
        throw std::invalid_argument(message);
    }
    const double tolerance =
        4096.0 *
        std::numeric_limits<double>::epsilon() *
        scale;
    if (distance <= tolerance) {
        throw std::invalid_argument(message);
    }
}

} // namespace cell_face_geometric_operator_3d_detail

[[nodiscard]] inline CellFaceGeometricOperator3D
make_cell_face_geometric_operator_3d(
    const Topology& topology,
    std::span<const Coordinate3D> vertex_coordinates_m,
    const FaceGeometry3D& face_geometry) {
    using namespace
        cell_face_geometric_operator_3d_detail;

    const std::size_t cell_count =
        topology.entity_count(EntityKind::cell);
    const std::size_t face_count =
        topology.entity_count(EntityKind::face);
    const std::size_t vertex_count =
        topology.entity_count(EntityKind::vertex);

    if (vertex_coordinates_m.size() !=
        vertex_count) {
        throw std::invalid_argument(
            "mpmc::mesh::make_cell_face_geometric_operator_3d: vertex coordinate count does not match topology");
    }
    if (face_geometry.cell_count() !=
            cell_count ||
        face_geometry.face_count() !=
            face_count) {
        throw std::invalid_argument(
            "mpmc::mesh::make_cell_face_geometric_operator_3d: FaceGeometry3D counts do not match topology");
    }
    if (!topology.has_relation(
            EntityKind::cell,
            EntityKind::vertex) ||
        !topology.has_relation(
            EntityKind::face,
            EntityKind::cell)) {
        throw std::invalid_argument(
            "mpmc::mesh::make_cell_face_geometric_operator_3d: cell->vertex and face->cell relations are required");
    }

    std::vector<Coordinate3D>
        cell_centroids_m;
    cell_centroids_m.reserve(cell_count);
    for (std::size_t cell = 0U;
         cell < cell_count;
         ++cell) {
        cell_centroids_m.push_back(
            cell_vertex_mean(
                topology,
                vertex_coordinates_m,
                LocalIndex{
                    static_cast<
                        LocalIndex::value_type>(
                            cell)}));
    }

    const auto& face_cells =
        topology.relation(
            EntityKind::face,
            EntityKind::cell);

    std::vector<LocalIndex> face_owners;
    std::vector<std::optional<LocalIndex>>
        face_neighbours;
    std::vector<Displacement3D>
        owner_displacements;
    std::vector<std::optional<Displacement3D>>
        neighbour_displacements;
    std::vector<double>
        owner_distances;
    std::vector<std::optional<double>>
        neighbour_distances;
    std::vector<std::optional<InternalFaceNonOrthogonality3D>>
        internal_non_orthogonality;
    std::vector<double> face_areas;
    std::vector<UnitVector3D> face_normals;

    face_owners.reserve(face_count);
    face_neighbours.reserve(face_count);
    owner_displacements.reserve(face_count);
    neighbour_displacements.reserve(face_count);
    owner_distances.reserve(face_count);
    neighbour_distances.reserve(face_count);
    internal_non_orthogonality.reserve(face_count);
    face_areas.reserve(face_count);
    face_normals.reserve(face_count);

    for (std::size_t face = 0U;
         face < face_count;
         ++face) {
        const LocalIndex local_face{
            static_cast<
                LocalIndex::value_type>(
                    face)};
        const auto support =
            face_cells.adjacent(local_face);
        if (support.size() != 1U &&
            support.size() != 2U) {
            throw std::invalid_argument(
                "mpmc::mesh::make_cell_face_geometric_operator_3d: face support must contain one boundary cell or two internal cells");
        }

        const LocalIndex owner =
            face_geometry.face_owner(
                local_face);
        if (support.front() != owner) {
            throw std::invalid_argument(
                "mpmc::mesh::make_cell_face_geometric_operator_3d: FaceGeometry3D owner must match canonical face->cell owner");
        }

        const std::size_t owner_local =
            static_cast<std::size_t>(
                owner.value());
        if (owner_local >= cell_count) {
            throw std::out_of_range(
                "mpmc::mesh::make_cell_face_geometric_operator_3d: owner cell index out of range");
        }

        const Coordinate3D face_centroid =
            face_geometry.face_centroid_m(
                local_face);
        const double area =
            face_geometry.face_area_m2(
                local_face);
        const UnitVector3D normal =
            face_geometry.face_owner_unit_normal(
                local_face);

        const Displacement3D owner_to_face =
            subtract(
                face_centroid,
                cell_centroids_m[owner_local]);
        const double owner_distance =
            dot(owner_to_face, normal);
        const double owner_scale =
            std::max(
                magnitude(owner_to_face),
                std::sqrt(area));
        require_positive_normal_distance(
            owner_distance,
            owner_scale,
            "mpmc::mesh::make_cell_face_geometric_operator_3d: owner normal distance is non-positive or degenerate");

        std::optional<LocalIndex> neighbour;
        std::optional<Displacement3D>
            neighbour_to_face;
        std::optional<double>
            neighbour_distance;
        std::optional<InternalFaceNonOrthogonality3D>
            non_orthogonality;

        if (support.size() == 2U) {
            if (support[1] == owner) {
                throw std::invalid_argument(
                    "mpmc::mesh::make_cell_face_geometric_operator_3d: internal face has duplicate owner/neighbour cell");
            }
            neighbour = support[1];
            const std::size_t neighbour_local =
                static_cast<std::size_t>(
                    neighbour->value());
            if (neighbour_local >= cell_count) {
                throw std::out_of_range(
                    "mpmc::mesh::make_cell_face_geometric_operator_3d: neighbour cell index out of range");
            }

            neighbour_to_face =
                subtract(
                    face_centroid,
                    cell_centroids_m[
                        neighbour_local]);
            neighbour_distance =
                -dot(
                    *neighbour_to_face,
                    normal);
            const double neighbour_scale =
                std::max(
                    magnitude(*neighbour_to_face),
                    std::sqrt(area));
            require_positive_normal_distance(
                *neighbour_distance,
                neighbour_scale,
                "mpmc::mesh::make_cell_face_geometric_operator_3d: neighbour normal distance is non-positive or degenerate");

            const Displacement3D centre_vector =
                owner_to_neighbour(
                    owner_to_face,
                    *neighbour_to_face);
            const double centre_distance =
                magnitude(centre_vector);
            const double centre_scale =
                std::max(
                    centre_distance,
                    std::sqrt(area));
            require_positive_normal_distance(
                centre_distance,
                centre_scale,
                "mpmc::mesh::make_cell_face_geometric_operator_3d: owner-to-neighbour cell-centre distance is degenerate");

            const double normal_projection =
                dot(centre_vector, normal);
            const double expected_projection =
                owner_distance +
                *neighbour_distance;
            const double projection_tolerance =
                4096.0 *
                std::numeric_limits<double>::epsilon() *
                centre_scale;
            if (!std::isfinite(normal_projection) ||
                std::abs(
                    normal_projection -
                    expected_projection) >
                    projection_tolerance ||
                normal_projection <=
                    projection_tolerance) {
                throw std::invalid_argument(
                    "mpmc::mesh::make_cell_face_geometric_operator_3d: inconsistent owner-to-neighbour normal projection");
            }

            double alignment_cosine =
                normal_projection /
                centre_distance;
            if (!std::isfinite(alignment_cosine) ||
                alignment_cosine <= 0.0 ||
                alignment_cosine >
                    1.0 +
                        128.0 *
                        std::numeric_limits<double>::epsilon()) {
                throw std::invalid_argument(
                    "mpmc::mesh::make_cell_face_geometric_operator_3d: invalid internal-face normal alignment");
            }
            alignment_cosine =
                std::clamp(
                    alignment_cosine,
                    0.0,
                    1.0);
            non_orthogonality =
                InternalFaceNonOrthogonality3D{
                    centre_vector,
                    centre_distance,
                    alignment_cosine,
                    std::acos(alignment_cosine)};
        }

        face_owners.push_back(owner);
        face_neighbours.push_back(neighbour);
        owner_displacements.push_back(
            owner_to_face);
        neighbour_displacements.push_back(
            neighbour_to_face);
        owner_distances.push_back(
            owner_distance);
        neighbour_distances.push_back(
            neighbour_distance);
        internal_non_orthogonality.push_back(
            non_orthogonality);
        face_areas.push_back(area);
        face_normals.push_back(normal);
    }

    return CellFaceGeometricOperator3D{
        std::move(cell_centroids_m),
        std::move(face_owners),
        std::move(face_neighbours),
        std::move(owner_displacements),
        std::move(neighbour_displacements),
        std::move(owner_distances),
        std::move(neighbour_distances),
        std::move(internal_non_orthogonality),
        std::move(face_areas),
        std::move(face_normals)};
}

} // namespace mpmc::mesh

#endif // MPMC_MESH_CELL_FACE_GEOMETRIC_OPERATOR_3D_HPP
