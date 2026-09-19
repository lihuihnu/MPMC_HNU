#ifndef MPMC_MESH_GEOMETRY_1D_HPP
#define MPMC_MESH_GEOMETRY_1D_HPP

#include <mpmc/mesh/entity.hpp>
#include <mpmc/mesh/topology.hpp>

#include <array>
#include <cmath>
#include <cstddef>
#include <limits>
#include <span>
#include <stdexcept>
#include <utility>
#include <vector>

namespace mpmc::mesh {

/// Immutable one-dimensional metric snapshot aligned to a Topology ordering.
///
/// Coordinates, cell centroids and cell lengths use SI metres. A 1D face is a
/// point, so this contract deliberately stores no fabricated face "area".
/// face_owner_unit_normals are scalar directions and must be -1 or +1 within
/// floating-point tolerance.
class Geometry1D {
public:
    Geometry1D(
        std::vector<double> vertex_coordinates_m,
        std::vector<double> cell_centroids_m,
        std::vector<double> cell_lengths_m,
        std::vector<double> face_coordinates_m,
        std::vector<LocalIndex> face_owners,
        std::vector<double>
            face_owner_unit_normals)
        : vertex_coordinates_m_(
              std::move(
                  vertex_coordinates_m)),
          cell_centroids_m_(
              std::move(
                  cell_centroids_m)),
          cell_lengths_m_(
              std::move(
                  cell_lengths_m)),
          face_coordinates_m_(
              std::move(
                  face_coordinates_m)),
          face_owners_(
              std::move(face_owners)),
          face_owner_unit_normals_(
              std::move(
                  face_owner_unit_normals)) {
        validate();
    }

    [[nodiscard]] std::size_t
    vertex_count() const noexcept {
        return vertex_coordinates_m_.size();
    }

    [[nodiscard]] std::size_t
    cell_count() const noexcept {
        return cell_centroids_m_.size();
    }

    [[nodiscard]] std::size_t
    face_count() const noexcept {
        return face_coordinates_m_.size();
    }

    [[nodiscard]] std::span<const double>
    vertex_coordinates_m() const noexcept {
        return vertex_coordinates_m_;
    }

    [[nodiscard]] std::span<const double>
    cell_centroids_m() const noexcept {
        return cell_centroids_m_;
    }

    [[nodiscard]] std::span<const double>
    cell_lengths_m() const noexcept {
        return cell_lengths_m_;
    }

    [[nodiscard]] std::span<const double>
    face_coordinates_m() const noexcept {
        return face_coordinates_m_;
    }

    [[nodiscard]] std::span<const LocalIndex>
    face_owners() const noexcept {
        return face_owners_;
    }

    [[nodiscard]] std::span<const double>
    face_owner_unit_normals()
        const noexcept {
        return face_owner_unit_normals_;
    }

    [[nodiscard]] double
    vertex_coordinate_m(
        LocalIndex vertex) const {
        return vertex_coordinates_m_.at(
            static_cast<std::size_t>(
                vertex.value()));
    }

    [[nodiscard]] double
    cell_centroid_m(
        LocalIndex cell) const {
        return cell_centroids_m_.at(
            static_cast<std::size_t>(
                cell.value()));
    }

    [[nodiscard]] double
    cell_length_m(
        LocalIndex cell) const {
        return cell_lengths_m_.at(
            static_cast<std::size_t>(
                cell.value()));
    }

    [[nodiscard]] double
    face_coordinate_m(
        LocalIndex face) const {
        return face_coordinates_m_.at(
            static_cast<std::size_t>(
                face.value()));
    }

    [[nodiscard]] LocalIndex
    face_owner(
        LocalIndex face) const {
        return face_owners_.at(
            static_cast<std::size_t>(
                face.value()));
    }

    [[nodiscard]] double
    face_owner_unit_normal(
        LocalIndex face) const {
        return face_owner_unit_normals_.at(
            static_cast<std::size_t>(
                face.value()));
    }

private:
    void validate() const {
        if (cell_centroids_m_.size() !=
                cell_lengths_m_.size()) {
            throw std::invalid_argument(
                "mpmc::mesh::Geometry1D: cell geometry array sizes do not match");
        }
        if (face_coordinates_m_.size() !=
                face_owners_.size() ||
            face_coordinates_m_.size() !=
                face_owner_unit_normals_.size()) {
            throw std::invalid_argument(
                "mpmc::mesh::Geometry1D: face geometry array sizes do not match");
        }

        for (const double coordinate :
             vertex_coordinates_m_) {
            if (!std::isfinite(coordinate)) {
                throw std::invalid_argument(
                    "mpmc::mesh::Geometry1D: non-finite vertex coordinate");
            }
        }
        for (const double centroid :
             cell_centroids_m_) {
            if (!std::isfinite(centroid)) {
                throw std::invalid_argument(
                    "mpmc::mesh::Geometry1D: non-finite cell centroid");
            }
        }
        for (const double length :
             cell_lengths_m_) {
            if (!std::isfinite(length) ||
                length <= 0.0) {
                throw std::invalid_argument(
                    "mpmc::mesh::Geometry1D: cell length must be finite and positive");
            }
        }

        constexpr double normal_tolerance =
            64.0 *
            std::numeric_limits<double>::
                epsilon();
        for (std::size_t face = 0U;
             face <
             face_coordinates_m_.size();
             ++face) {
            if (!std::isfinite(
                    face_coordinates_m_[face])) {
                throw std::invalid_argument(
                    "mpmc::mesh::Geometry1D: non-finite face coordinate");
            }
            if (static_cast<std::size_t>(
                    face_owners_[face].value()) >=
                cell_centroids_m_.size()) {
                throw std::out_of_range(
                    "mpmc::mesh::Geometry1D: face owner index out of range");
            }
            const double normal =
                face_owner_unit_normals_[face];
            if (!std::isfinite(normal) ||
                std::abs(
                    std::abs(normal) -
                    1.0) >
                    normal_tolerance) {
                throw std::invalid_argument(
                    "mpmc::mesh::Geometry1D: face normal is not unit magnitude");
            }
        }
    }

    std::vector<double>
        vertex_coordinates_m_;
    std::vector<double>
        cell_centroids_m_;
    std::vector<double>
        cell_lengths_m_;
    std::vector<double>
        face_coordinates_m_;
    std::vector<LocalIndex>
        face_owners_;
    std::vector<double>
        face_owner_unit_normals_;
};

namespace geometry_1d_detail {

[[nodiscard]] inline LocalIndex
local_index(std::size_t value) {
    if (value >
        static_cast<std::size_t>(
            std::numeric_limits<
                LocalIndex::value_type>::max())) {
        throw std::length_error(
            "mpmc::mesh::make_cartesian_geometry_1d: local index overflow");
    }
    return LocalIndex{
        static_cast<
            LocalIndex::value_type>(
            value)};
}

inline void require_axis(
    std::span<const double> axis) {
    if (axis.size() < 2U) {
        throw std::invalid_argument(
            "mpmc::mesh::make_cartesian_geometry_1d: x coordinates must contain at least two values");
    }
    if (!std::isfinite(axis.front())) {
        throw std::invalid_argument(
            "mpmc::mesh::make_cartesian_geometry_1d: x coordinates must be finite and strictly increasing");
    }
    for (std::size_t i = 1U;
         i < axis.size();
         ++i) {
        if (!std::isfinite(axis[i])) {
            throw std::invalid_argument(
                "mpmc::mesh::make_cartesian_geometry_1d: x coordinates must be finite and strictly increasing");
        }
        const double width =
            axis[i] - axis[i - 1U];
        if (!std::isfinite(width) ||
            width <= 0.0) {
            throw std::invalid_argument(
                "mpmc::mesh::make_cartesian_geometry_1d: x coordinates must be finite and strictly increasing");
        }
    }
}

inline void require_row(
    std::span<const LocalIndex> actual,
    std::span<const std::size_t> expected,
    const char* message) {
    if (actual.size() !=
        expected.size()) {
        throw std::invalid_argument(message);
    }
    for (std::size_t i = 0U;
         i < expected.size();
         ++i) {
        if (static_cast<std::size_t>(
                actual[i].value()) !=
            expected[i]) {
            throw std::invalid_argument(
                message);
        }
    }
}

inline void require_cartesian_topology(
    const Topology& topology,
    std::size_t nx) {
    const std::size_t vertex_count =
        nx + 1U;
    const std::size_t face_count =
        vertex_count;
    if (topology.entity_count(
            EntityKind::vertex) !=
            vertex_count ||
        topology.entity_count(
            EntityKind::edge) != 0U ||
        topology.entity_count(
            EntityKind::face) !=
            face_count ||
        topology.entity_count(
            EntityKind::cell) != nx) {
        throw std::invalid_argument(
            "mpmc::mesh::make_cartesian_geometry_1d: topology entity counts do not match axis coordinates");
    }

    for (const auto& relation :
         {std::pair{
              EntityKind::cell,
              EntityKind::vertex},
          std::pair{
              EntityKind::cell,
              EntityKind::face},
          std::pair{
              EntityKind::face,
              EntityKind::vertex},
          std::pair{
              EntityKind::face,
              EntityKind::cell}}) {
        if (!topology.has_relation(
                relation.first,
                relation.second)) {
            throw std::invalid_argument(
                "mpmc::mesh::make_cartesian_geometry_1d: required Cartesian relation missing");
        }
    }

    const auto& cell_vertex =
        topology.relation(
            EntityKind::cell,
            EntityKind::vertex);
    const auto& cell_face =
        topology.relation(
            EntityKind::cell,
            EntityKind::face);
    const auto& face_vertex =
        topology.relation(
            EntityKind::face,
            EntityKind::vertex);
    const auto& face_cell =
        topology.relation(
            EntityKind::face,
            EntityKind::cell);

    for (std::size_t cell = 0U;
         cell < nx;
         ++cell) {
        const std::array<std::size_t, 2>
            vertices{
                cell,
                cell + 1U};
        require_row(
            cell_vertex.adjacent(
                local_index(cell)),
            vertices,
            "mpmc::mesh::make_cartesian_geometry_1d: cell-to-vertex relation is not canonical");
        require_row(
            cell_face.adjacent(
                local_index(cell)),
            vertices,
            "mpmc::mesh::make_cartesian_geometry_1d: cell-to-face relation is not canonical");
    }

    for (std::size_t face = 0U;
         face < face_count;
         ++face) {
        const std::array<std::size_t, 1>
            vertex{face};
        require_row(
            face_vertex.adjacent(
                local_index(face)),
            vertex,
            "mpmc::mesh::make_cartesian_geometry_1d: face-to-vertex relation is not canonical");

        if (face == 0U) {
            const std::array<std::size_t, 1>
                cells{0U};
            require_row(
                face_cell.adjacent(
                    local_index(face)),
                cells,
                "mpmc::mesh::make_cartesian_geometry_1d: left boundary owner mismatch");
        } else if (face == nx) {
            const std::array<std::size_t, 1>
                cells{nx - 1U};
            require_row(
                face_cell.adjacent(
                    local_index(face)),
                cells,
                "mpmc::mesh::make_cartesian_geometry_1d: right boundary owner mismatch");
        } else {
            const std::array<std::size_t, 2>
                cells{
                    face - 1U,
                    face};
            require_row(
                face_cell.adjacent(
                    local_index(face)),
                cells,
                "mpmc::mesh::make_cartesian_geometry_1d: internal face owner/neighbour mismatch");
        }
    }
}

} // namespace geometry_1d_detail

/// Build metric geometry for a canonical 1D Cartesian topology.
[[nodiscard]] inline Geometry1D
make_cartesian_geometry_1d(
    const Topology& topology,
    std::span<const double>
        x_coordinates_m) {
    using namespace geometry_1d_detail;

    require_axis(x_coordinates_m);
    const std::size_t nx =
        x_coordinates_m.size() - 1U;
    require_cartesian_topology(
        topology,
        nx);

    std::vector<double>
        cell_centroids;
    std::vector<double>
        cell_lengths;
    cell_centroids.reserve(nx);
    cell_lengths.reserve(nx);
    for (std::size_t cell = 0U;
         cell < nx;
         ++cell) {
        const double width =
            x_coordinates_m[cell + 1U] -
            x_coordinates_m[cell];
        const double centroid =
            x_coordinates_m[cell] +
            0.5 * width;
        if (!std::isfinite(width) ||
            width <= 0.0 ||
            !std::isfinite(centroid)) {
            throw std::invalid_argument(
                "mpmc::mesh::make_cartesian_geometry_1d: invalid cell metric");
        }
        cell_lengths.push_back(width);
        cell_centroids.push_back(
            centroid);
    }

    std::vector<LocalIndex>
        face_owners;
    std::vector<double>
        face_normals;
    face_owners.reserve(nx + 1U);
    face_normals.reserve(nx + 1U);
    for (std::size_t face = 0U;
         face <= nx;
         ++face) {
        if (face == 0U) {
            face_owners.push_back(
                local_index(0U));
            face_normals.push_back(-1.0);
        } else {
            face_owners.push_back(
                local_index(face - 1U));
            face_normals.push_back(1.0);
        }
    }

    return Geometry1D{
        std::vector<double>{
            x_coordinates_m.begin(),
            x_coordinates_m.end()},
        std::move(cell_centroids),
        std::move(cell_lengths),
        std::vector<double>{
            x_coordinates_m.begin(),
            x_coordinates_m.end()},
        std::move(face_owners),
        std::move(face_normals)};
}

} // namespace mpmc::mesh

#endif // MPMC_MESH_GEOMETRY_1D_HPP
