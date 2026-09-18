#ifndef MPMC_MESH_GEOMETRY_2D_HPP
#define MPMC_MESH_GEOMETRY_2D_HPP

#include <mpmc/mesh/entity.hpp>
#include <mpmc/mesh/topology.hpp>

#include <array>
#include <cmath>
#include <cstddef>
#include <limits>
#include <numeric>
#include <span>
#include <stdexcept>
#include <utility>
#include <vector>

namespace mpmc::mesh {

struct Coordinate2D {
    double x_m;
    double y_m;
};

struct UnitVector2D {
    double x;
    double y;
};

/// Immutable two-dimensional metric snapshot aligned to a Topology local ordering.
///
/// Coordinates and lengths use SI metres; cell measures use square metres.
/// Face normals are unit vectors directed outward from the stored owner cell.
/// The snapshot owns all arrays and retains no reference to the source Topology.
class Geometry2D {
public:
    Geometry2D(std::vector<Coordinate2D> vertex_coordinates_m,
               std::vector<Coordinate2D> cell_centroids_m,
               std::vector<double> cell_areas_m2,
               std::vector<Coordinate2D> face_centroids_m,
               std::vector<double> face_lengths_m,
               std::vector<LocalIndex> face_owners,
               std::vector<UnitVector2D> face_owner_unit_normals)
        : vertex_coordinates_m_(std::move(vertex_coordinates_m)),
          cell_centroids_m_(std::move(cell_centroids_m)),
          cell_areas_m2_(std::move(cell_areas_m2)),
          face_centroids_m_(std::move(face_centroids_m)),
          face_lengths_m_(std::move(face_lengths_m)),
          face_owners_(std::move(face_owners)),
          face_owner_unit_normals_(std::move(face_owner_unit_normals)) {
        validate();
    }

    Geometry2D(const Geometry2D&) = default;
    Geometry2D(Geometry2D&&) noexcept = default;
    Geometry2D& operator=(const Geometry2D&) = delete;
    Geometry2D& operator=(Geometry2D&&) = delete;
    ~Geometry2D() = default;

    [[nodiscard]] std::size_t vertex_count() const noexcept {
        return vertex_coordinates_m_.size();
    }
    [[nodiscard]] std::size_t cell_count() const noexcept {
        return cell_centroids_m_.size();
    }
    [[nodiscard]] std::size_t face_count() const noexcept {
        return face_centroids_m_.size();
    }

    [[nodiscard]] std::span<const Coordinate2D> vertex_coordinates_m() const noexcept {
        return vertex_coordinates_m_;
    }
    [[nodiscard]] std::span<const Coordinate2D> cell_centroids_m() const noexcept {
        return cell_centroids_m_;
    }
    [[nodiscard]] std::span<const double> cell_areas_m2() const noexcept {
        return cell_areas_m2_;
    }
    [[nodiscard]] std::span<const Coordinate2D> face_centroids_m() const noexcept {
        return face_centroids_m_;
    }
    [[nodiscard]] std::span<const double> face_lengths_m() const noexcept {
        return face_lengths_m_;
    }
    [[nodiscard]] std::span<const LocalIndex> face_owners() const noexcept {
        return face_owners_;
    }
    [[nodiscard]] std::span<const UnitVector2D> face_owner_unit_normals() const noexcept {
        return face_owner_unit_normals_;
    }

    [[nodiscard]] Coordinate2D vertex_coordinate_m(LocalIndex vertex) const {
        return vertex_coordinates_m_.at(static_cast<std::size_t>(vertex.value()));
    }
    [[nodiscard]] Coordinate2D cell_centroid_m(LocalIndex cell) const {
        return cell_centroids_m_.at(static_cast<std::size_t>(cell.value()));
    }
    [[nodiscard]] double cell_area_m2(LocalIndex cell) const {
        return cell_areas_m2_.at(static_cast<std::size_t>(cell.value()));
    }
    [[nodiscard]] Coordinate2D face_centroid_m(LocalIndex face) const {
        return face_centroids_m_.at(static_cast<std::size_t>(face.value()));
    }
    [[nodiscard]] double face_length_m(LocalIndex face) const {
        return face_lengths_m_.at(static_cast<std::size_t>(face.value()));
    }
    [[nodiscard]] LocalIndex face_owner(LocalIndex face) const {
        return face_owners_.at(static_cast<std::size_t>(face.value()));
    }
    [[nodiscard]] UnitVector2D face_owner_unit_normal(LocalIndex face) const {
        return face_owner_unit_normals_.at(static_cast<std::size_t>(face.value()));
    }

private:
    static void require_finite_coordinate(Coordinate2D value, const char* message) {
        if (!std::isfinite(value.x_m) || !std::isfinite(value.y_m)) {
            throw std::invalid_argument(message);
        }
    }

    void validate() const {
        if (cell_centroids_m_.size() != cell_areas_m2_.size()) {
            throw std::invalid_argument(
                "mpmc::mesh::Geometry2D: cell geometry array sizes do not match");
        }
        if (face_centroids_m_.size() != face_lengths_m_.size() ||
            face_centroids_m_.size() != face_owners_.size() ||
            face_centroids_m_.size() != face_owner_unit_normals_.size()) {
            throw std::invalid_argument(
                "mpmc::mesh::Geometry2D: face geometry array sizes do not match");
        }

        for (const Coordinate2D coordinate : vertex_coordinates_m_) {
            require_finite_coordinate(
                coordinate, "mpmc::mesh::Geometry2D: non-finite vertex coordinate");
        }
        for (const Coordinate2D centroid : cell_centroids_m_) {
            require_finite_coordinate(
                centroid, "mpmc::mesh::Geometry2D: non-finite cell centroid");
        }
        for (const double area : cell_areas_m2_) {
            if (!std::isfinite(area) || area <= 0.0) {
                throw std::invalid_argument(
                    "mpmc::mesh::Geometry2D: cell area must be finite and positive");
            }
        }

        constexpr double normal_tolerance =
            64.0 * std::numeric_limits<double>::epsilon();
        for (std::size_t face = 0; face < face_centroids_m_.size(); ++face) {
            require_finite_coordinate(
                face_centroids_m_[face],
                "mpmc::mesh::Geometry2D: non-finite face centroid");
            const double length = face_lengths_m_[face];
            if (!std::isfinite(length) || length <= 0.0) {
                throw std::invalid_argument(
                    "mpmc::mesh::Geometry2D: face length must be finite and positive");
            }
            if (static_cast<std::size_t>(face_owners_[face].value()) >=
                cell_centroids_m_.size()) {
                throw std::out_of_range(
                    "mpmc::mesh::Geometry2D: face owner index out of range");
            }

            const UnitVector2D normal = face_owner_unit_normals_[face];
            if (!std::isfinite(normal.x) || !std::isfinite(normal.y)) {
                throw std::invalid_argument(
                    "mpmc::mesh::Geometry2D: non-finite face unit normal");
            }
            const double magnitude = std::hypot(normal.x, normal.y);
            if (!std::isfinite(magnitude) ||
                std::abs(magnitude - 1.0) > normal_tolerance) {
                throw std::invalid_argument(
                    "mpmc::mesh::Geometry2D: face normal is not unit length");
            }
        }
    }

    std::vector<Coordinate2D> vertex_coordinates_m_;
    std::vector<Coordinate2D> cell_centroids_m_;
    std::vector<double> cell_areas_m2_;
    std::vector<Coordinate2D> face_centroids_m_;
    std::vector<double> face_lengths_m_;
    std::vector<LocalIndex> face_owners_;
    std::vector<UnitVector2D> face_owner_unit_normals_;
};

namespace detail {

[[nodiscard]] inline std::size_t geometry_checked_add(std::size_t left,
                                                       std::size_t right,
                                                       const char* message) {
    if (left > std::numeric_limits<std::size_t>::max() - right) {
        throw std::length_error(message);
    }
    return left + right;
}

[[nodiscard]] inline std::size_t geometry_checked_multiply(std::size_t left,
                                                            std::size_t right,
                                                            const char* message) {
    if (left != 0U && right > std::numeric_limits<std::size_t>::max() / left) {
        throw std::length_error(message);
    }
    return left * right;
}

[[nodiscard]] inline LocalIndex geometry_local_index(std::size_t value) {
    if (value > static_cast<std::size_t>(
                    std::numeric_limits<LocalIndex::value_type>::max())) {
        throw std::length_error(
            "mpmc::mesh::make_cartesian_geometry_2d: local index overflow");
    }
    return LocalIndex{static_cast<LocalIndex::value_type>(value)};
}

inline void require_axis(std::span<const double> axis, const char* message) {
    if (axis.size() < 2U) {
        throw std::invalid_argument(message);
    }
    if (!std::isfinite(axis.front())) {
        throw std::invalid_argument(message);
    }
    for (std::size_t i = 1U; i < axis.size(); ++i) {
        if (!std::isfinite(axis[i])) {
            throw std::invalid_argument(message);
        }
        const double width = axis[i] - axis[i - 1U];
        if (!std::isfinite(width) || width <= 0.0) {
            throw std::invalid_argument(message);
        }
    }
}

inline void require_row(std::span<const LocalIndex> actual,
                        std::span<const std::size_t> expected,
                        const char* message) {
    if (actual.size() != expected.size()) {
        throw std::invalid_argument(message);
    }
    for (std::size_t i = 0; i < expected.size(); ++i) {
        if (static_cast<std::size_t>(actual[i].value()) != expected[i]) {
            throw std::invalid_argument(message);
        }
    }
}

inline void require_cartesian_topology_2d(const Topology& topology,
                                          std::size_t nx,
                                          std::size_t ny) {
    const std::size_t nx_vertices =
        geometry_checked_add(nx, 1U,
                             "mpmc::mesh::make_cartesian_geometry_2d: nx overflow");
    const std::size_t ny_vertices =
        geometry_checked_add(ny, 1U,
                             "mpmc::mesh::make_cartesian_geometry_2d: ny overflow");
    const std::size_t vertex_count =
        geometry_checked_multiply(
            nx_vertices, ny_vertices,
            "mpmc::mesh::make_cartesian_geometry_2d: vertex count overflow");
    const std::size_t cell_count =
        geometry_checked_multiply(
            nx, ny,
            "mpmc::mesh::make_cartesian_geometry_2d: cell count overflow");
    const std::size_t vertical_face_count =
        geometry_checked_multiply(
            nx_vertices, ny,
            "mpmc::mesh::make_cartesian_geometry_2d: vertical face count overflow");
    const std::size_t horizontal_face_count =
        geometry_checked_multiply(
            nx, ny_vertices,
            "mpmc::mesh::make_cartesian_geometry_2d: horizontal face count overflow");
    const std::size_t face_count =
        geometry_checked_add(
            vertical_face_count, horizontal_face_count,
            "mpmc::mesh::make_cartesian_geometry_2d: face count overflow");

    if (topology.entity_count(EntityKind::vertex) != vertex_count ||
        topology.entity_count(EntityKind::edge) != 0U ||
        topology.entity_count(EntityKind::face) != face_count ||
        topology.entity_count(EntityKind::cell) != cell_count) {
        throw std::invalid_argument(
            "mpmc::mesh::make_cartesian_geometry_2d: topology entity counts "
            "do not match axis coordinates");
    }
    if (!topology.has_relation(EntityKind::cell, EntityKind::vertex) ||
        !topology.has_relation(EntityKind::cell, EntityKind::face) ||
        !topology.has_relation(EntityKind::face, EntityKind::vertex) ||
        !topology.has_relation(EntityKind::face, EntityKind::cell)) {
        throw std::invalid_argument(
            "mpmc::mesh::make_cartesian_geometry_2d: required Cartesian relation missing");
    }

    const auto& cell_vertex = topology.relation(EntityKind::cell, EntityKind::vertex);
    const auto& cell_face = topology.relation(EntityKind::cell, EntityKind::face);
    const auto& face_vertex = topology.relation(EntityKind::face, EntityKind::vertex);
    const auto& face_cell = topology.relation(EntityKind::face, EntityKind::cell);

    const auto vertex = [nx_vertices](std::size_t i, std::size_t j) {
        return j * nx_vertices + i;
    };
    const auto cell = [nx](std::size_t i, std::size_t j) {
        return j * nx + i;
    };
    const auto vertical_face = [nx_vertices](std::size_t i, std::size_t j) {
        return j * nx_vertices + i;
    };
    const auto horizontal_face =
        [nx, vertical_face_count](std::size_t i, std::size_t j) {
            return vertical_face_count + j * nx + i;
        };

    for (std::size_t j = 0; j < ny; ++j) {
        for (std::size_t i = 0; i < nx; ++i) {
            const std::size_t cell_index = cell(i, j);
            const std::array<std::size_t, 4> vertices{
                vertex(i, j), vertex(i + 1U, j),
                vertex(i + 1U, j + 1U), vertex(i, j + 1U)};
            require_row(
                cell_vertex.adjacent(geometry_local_index(cell_index)),
                vertices,
                "mpmc::mesh::make_cartesian_geometry_2d: cell-to-vertex "
                "relation is not canonical Cartesian topology");
            const std::array<std::size_t, 4> faces{
                vertical_face(i, j), vertical_face(i + 1U, j),
                horizontal_face(i, j), horizontal_face(i, j + 1U)};
            require_row(
                cell_face.adjacent(geometry_local_index(cell_index)),
                faces,
                "mpmc::mesh::make_cartesian_geometry_2d: cell-to-face "
                "relation is not canonical Cartesian topology");
        }
    }

    for (std::size_t j = 0; j < ny; ++j) {
        for (std::size_t i = 0; i < nx_vertices; ++i) {
            const std::size_t face = vertical_face(i, j);
            const std::array<std::size_t, 2> vertices{
                vertex(i, j), vertex(i, j + 1U)};
            require_row(
                face_vertex.adjacent(geometry_local_index(face)),
                vertices,
                "mpmc::mesh::make_cartesian_geometry_2d: vertical "
                "face-to-vertex relation is not canonical");
            if (i == 0U) {
                const std::array<std::size_t, 1> cells{cell(0U, j)};
                require_row(face_cell.adjacent(geometry_local_index(face)), cells,
                            "mpmc::mesh::make_cartesian_geometry_2d: left boundary owner mismatch");
            } else if (i == nx) {
                const std::array<std::size_t, 1> cells{cell(nx - 1U, j)};
                require_row(face_cell.adjacent(geometry_local_index(face)), cells,
                            "mpmc::mesh::make_cartesian_geometry_2d: right boundary owner mismatch");
            } else {
                const std::array<std::size_t, 2> cells{
                    cell(i - 1U, j), cell(i, j)};
                require_row(face_cell.adjacent(geometry_local_index(face)), cells,
                            "mpmc::mesh::make_cartesian_geometry_2d: vertical internal owner mismatch");
            }
        }
    }

    for (std::size_t j = 0; j < ny_vertices; ++j) {
        for (std::size_t i = 0; i < nx; ++i) {
            const std::size_t face = horizontal_face(i, j);
            const std::array<std::size_t, 2> vertices{
                vertex(i, j), vertex(i + 1U, j)};
            require_row(
                face_vertex.adjacent(geometry_local_index(face)),
                vertices,
                "mpmc::mesh::make_cartesian_geometry_2d: horizontal "
                "face-to-vertex relation is not canonical");
            if (j == 0U) {
                const std::array<std::size_t, 1> cells{cell(i, 0U)};
                require_row(face_cell.adjacent(geometry_local_index(face)), cells,
                            "mpmc::mesh::make_cartesian_geometry_2d: bottom boundary owner mismatch");
            } else if (j == ny) {
                const std::array<std::size_t, 1> cells{cell(i, ny - 1U)};
                require_row(face_cell.adjacent(geometry_local_index(face)), cells,
                            "mpmc::mesh::make_cartesian_geometry_2d: top boundary owner mismatch");
            } else {
                const std::array<std::size_t, 2> cells{
                    cell(i, j - 1U), cell(i, j)};
                require_row(face_cell.adjacent(geometry_local_index(face)), cells,
                            "mpmc::mesh::make_cartesian_geometry_2d: horizontal internal owner mismatch");
            }
        }
    }
}

} // namespace detail

/// Build metric geometry for a canonical 2D Cartesian Topology.
///
/// x_coordinates_m and y_coordinates_m are strictly increasing axis coordinates
/// in SI metres. The source topology is not modified and is not retained.
[[nodiscard]] inline Geometry2D make_cartesian_geometry_2d(
    const Topology& topology,
    std::span<const double> x_coordinates_m,
    std::span<const double> y_coordinates_m) {
    detail::require_axis(
        x_coordinates_m,
        "mpmc::mesh::make_cartesian_geometry_2d: x coordinates must be finite, "
        "strictly increasing, and contain at least two values");
    detail::require_axis(
        y_coordinates_m,
        "mpmc::mesh::make_cartesian_geometry_2d: y coordinates must be finite, "
        "strictly increasing, and contain at least two values");

    const std::size_t nx = x_coordinates_m.size() - 1U;
    const std::size_t ny = y_coordinates_m.size() - 1U;
    detail::require_cartesian_topology_2d(topology, nx, ny);

    const std::size_t vertex_count = topology.entity_count(EntityKind::vertex);
    const std::size_t cell_count = topology.entity_count(EntityKind::cell);
    const std::size_t face_count = topology.entity_count(EntityKind::face);
    const std::size_t nx_vertices = x_coordinates_m.size();

    std::vector<double> dx;
    std::vector<double> dy;
    dx.reserve(nx);
    dy.reserve(ny);
    for (std::size_t i = 0; i < nx; ++i) {
        dx.push_back(x_coordinates_m[i + 1U] - x_coordinates_m[i]);
    }
    for (std::size_t j = 0; j < ny; ++j) {
        dy.push_back(y_coordinates_m[j + 1U] - y_coordinates_m[j]);
    }

    std::vector<Coordinate2D> vertex_coordinates;
    vertex_coordinates.reserve(vertex_count);
    for (std::size_t j = 0; j < y_coordinates_m.size(); ++j) {
        for (std::size_t i = 0; i < x_coordinates_m.size(); ++i) {
            vertex_coordinates.push_back(
                Coordinate2D{x_coordinates_m[i], y_coordinates_m[j]});
        }
    }

    std::vector<Coordinate2D> cell_centroids;
    std::vector<double> cell_areas;
    cell_centroids.reserve(cell_count);
    cell_areas.reserve(cell_count);
    for (std::size_t j = 0; j < ny; ++j) {
        for (std::size_t i = 0; i < nx; ++i) {
            const Coordinate2D centroid{
                std::midpoint(x_coordinates_m[i], x_coordinates_m[i + 1U]),
                std::midpoint(y_coordinates_m[j], y_coordinates_m[j + 1U])};
            const double area = dx[i] * dy[j];
            if (!std::isfinite(area) || area <= 0.0) {
                throw std::invalid_argument(
                    "mpmc::mesh::make_cartesian_geometry_2d: cell area is not finite and positive");
            }
            cell_centroids.push_back(centroid);
            cell_areas.push_back(area);
        }
    }

    const auto& face_cell = topology.relation(EntityKind::face, EntityKind::cell);
    std::vector<Coordinate2D> face_centroids;
    std::vector<double> face_lengths;
    std::vector<LocalIndex> face_owners;
    std::vector<UnitVector2D> face_normals;
    face_centroids.reserve(face_count);
    face_lengths.reserve(face_count);
    face_owners.reserve(face_count);
    face_normals.reserve(face_count);

    const auto append_face =
        [&](Coordinate2D centroid, double length, std::size_t face_index) {
            if (!std::isfinite(length) || length <= 0.0) {
                throw std::invalid_argument(
                    "mpmc::mesh::make_cartesian_geometry_2d: face length is not finite and positive");
            }
            const auto adjacent_cells =
                face_cell.adjacent(detail::geometry_local_index(face_index));
            if (adjacent_cells.empty()) {
                throw std::invalid_argument(
                    "mpmc::mesh::make_cartesian_geometry_2d: face has no owner cell");
            }
            const LocalIndex owner = adjacent_cells.front();
            const Coordinate2D owner_centroid =
                cell_centroids.at(static_cast<std::size_t>(owner.value()));
            const double delta_x = centroid.x_m - owner_centroid.x_m;
            const double delta_y = centroid.y_m - owner_centroid.y_m;
            const double distance = std::hypot(delta_x, delta_y);
            if (!std::isfinite(distance) || distance <= 0.0) {
                throw std::invalid_argument(
                    "mpmc::mesh::make_cartesian_geometry_2d: invalid owner-to-face distance");
            }
            face_centroids.push_back(centroid);
            face_lengths.push_back(length);
            face_owners.push_back(owner);
            face_normals.push_back(UnitVector2D{delta_x / distance, delta_y / distance});
        };

    std::size_t face_index = 0U;
    for (std::size_t j = 0; j < ny; ++j) {
        for (std::size_t i = 0; i < nx_vertices; ++i) {
            append_face(
                Coordinate2D{
                    x_coordinates_m[i],
                    std::midpoint(y_coordinates_m[j], y_coordinates_m[j + 1U])},
                dy[j],
                face_index);
            ++face_index;
        }
    }
    for (std::size_t j = 0; j < y_coordinates_m.size(); ++j) {
        for (std::size_t i = 0; i < nx; ++i) {
            append_face(
                Coordinate2D{
                    std::midpoint(x_coordinates_m[i], x_coordinates_m[i + 1U]),
                    y_coordinates_m[j]},
                dx[i],
                face_index);
            ++face_index;
        }
    }
    if (face_index != face_count) {
        throw std::logic_error(
            "mpmc::mesh::make_cartesian_geometry_2d: internal face enumeration mismatch");
    }

    return Geometry2D{
        std::move(vertex_coordinates),
        std::move(cell_centroids),
        std::move(cell_areas),
        std::move(face_centroids),
        std::move(face_lengths),
        std::move(face_owners),
        std::move(face_normals)};
}

} // namespace mpmc::mesh

#endif // MPMC_MESH_GEOMETRY_2D_HPP
