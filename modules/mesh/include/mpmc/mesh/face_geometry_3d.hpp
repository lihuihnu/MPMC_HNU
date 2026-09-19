#ifndef MPMC_MESH_FACE_GEOMETRY_3D_HPP
#define MPMC_MESH_FACE_GEOMETRY_3D_HPP

#include <mpmc/mesh/corner_point_geometry_3d.hpp>
#include <mpmc/mesh/entity.hpp>

#include <cmath>
#include <cstddef>
#include <limits>
#include <span>
#include <stdexcept>
#include <utility>
#include <vector>

namespace mpmc::mesh {

struct UnitVector3D {
    double x;
    double y;
    double z;
};

/// Immutable 3D quad-face metric snapshot aligned to a Topology face ordering.
///
/// Face centroids use SI metres, areas use square metres, and normals are unit
/// vectors directed outward from the stored owner cell.
class FaceGeometry3D {
public:
    FaceGeometry3D(
        std::size_t cell_count,
        std::vector<Coordinate3D> face_centroids_m,
        std::vector<double> face_areas_m2,
        std::vector<LocalIndex> face_owners,
        std::vector<UnitVector3D> face_owner_unit_normals)
        : cell_count_(cell_count),
          face_centroids_m_(
              std::move(face_centroids_m)),
          face_areas_m2_(
              std::move(face_areas_m2)),
          face_owners_(
              std::move(face_owners)),
          face_owner_unit_normals_(
              std::move(face_owner_unit_normals)) {
        validate();
    }

    FaceGeometry3D(const FaceGeometry3D&) = default;
    FaceGeometry3D(FaceGeometry3D&&) noexcept = default;
    FaceGeometry3D& operator=(const FaceGeometry3D&) = delete;
    FaceGeometry3D& operator=(FaceGeometry3D&&) = delete;
    ~FaceGeometry3D() = default;

    [[nodiscard]] std::size_t cell_count() const noexcept {
        return cell_count_;
    }

    [[nodiscard]] std::size_t face_count() const noexcept {
        return face_centroids_m_.size();
    }

    [[nodiscard]] std::span<const Coordinate3D>
    face_centroids_m() const noexcept {
        return face_centroids_m_;
    }

    [[nodiscard]] std::span<const double>
    face_areas_m2() const noexcept {
        return face_areas_m2_;
    }

    [[nodiscard]] std::span<const LocalIndex>
    face_owners() const noexcept {
        return face_owners_;
    }

    [[nodiscard]] std::span<const UnitVector3D>
    face_owner_unit_normals() const noexcept {
        return face_owner_unit_normals_;
    }

    [[nodiscard]] Coordinate3D
    face_centroid_m(LocalIndex face) const {
        return face_centroids_m_.at(
            static_cast<std::size_t>(
                face.value()));
    }

    [[nodiscard]] double
    face_area_m2(LocalIndex face) const {
        return face_areas_m2_.at(
            static_cast<std::size_t>(
                face.value()));
    }

    [[nodiscard]] LocalIndex
    face_owner(LocalIndex face) const {
        return face_owners_.at(
            static_cast<std::size_t>(
                face.value()));
    }

    [[nodiscard]] UnitVector3D
    face_owner_unit_normal(LocalIndex face) const {
        return face_owner_unit_normals_.at(
            static_cast<std::size_t>(
                face.value()));
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

    void validate() const {
        if (face_centroids_m_.size() !=
                face_areas_m2_.size() ||
            face_centroids_m_.size() !=
                face_owners_.size() ||
            face_centroids_m_.size() !=
                face_owner_unit_normals_.size()) {
            throw std::invalid_argument(
                "mpmc::mesh::FaceGeometry3D: face geometry array sizes do not match");
        }

        constexpr double normal_tolerance =
            128.0 *
            std::numeric_limits<double>::epsilon();

        for (std::size_t face = 0U;
             face < face_centroids_m_.size();
             ++face) {
            require_finite_coordinate(
                face_centroids_m_[face],
                "mpmc::mesh::FaceGeometry3D: non-finite face centroid");

            const double area =
                face_areas_m2_[face];
            if (!std::isfinite(area) ||
                area <= 0.0) {
                throw std::invalid_argument(
                    "mpmc::mesh::FaceGeometry3D: face area must be finite and positive");
            }

            if (static_cast<std::size_t>(
                    face_owners_[face].value()) >=
                cell_count_) {
                throw std::out_of_range(
                    "mpmc::mesh::FaceGeometry3D: face owner index out of range");
            }

            const auto normal =
                face_owner_unit_normals_[face];
            if (!std::isfinite(normal.x) ||
                !std::isfinite(normal.y) ||
                !std::isfinite(normal.z)) {
                throw std::invalid_argument(
                    "mpmc::mesh::FaceGeometry3D: non-finite face unit normal");
            }
            const double magnitude =
                std::sqrt(
                    normal.x * normal.x +
                    normal.y * normal.y +
                    normal.z * normal.z);
            if (!std::isfinite(magnitude) ||
                std::abs(magnitude - 1.0) >
                    normal_tolerance) {
                throw std::invalid_argument(
                    "mpmc::mesh::FaceGeometry3D: face normal is not unit length");
            }
        }
    }

    std::size_t cell_count_;
    std::vector<Coordinate3D> face_centroids_m_;
    std::vector<double> face_areas_m2_;
    std::vector<LocalIndex> face_owners_;
    std::vector<UnitVector3D> face_owner_unit_normals_;
};

} // namespace mpmc::mesh

#endif // MPMC_MESH_FACE_GEOMETRY_3D_HPP
