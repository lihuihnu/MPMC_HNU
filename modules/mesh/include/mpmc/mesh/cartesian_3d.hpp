#ifndef MPMC_MESH_CARTESIAN_3D_HPP
#define MPMC_MESH_CARTESIAN_3D_HPP

#include <mpmc/mesh/csr_adjacency.hpp>
#include <mpmc/mesh/linear_cell_mesh_3d.hpp>

#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <span>
#include <stdexcept>
#include <utility>
#include <vector>

namespace mpmc::mesh {

namespace cartesian_3d_detail {

[[nodiscard]] inline std::size_t checked_add(
    std::size_t left,
    std::size_t right,
    const char* message) {
    if (left >
        std::numeric_limits<std::size_t>::max() -
            right) {
        throw std::length_error(message);
    }
    return left + right;
}

[[nodiscard]] inline std::size_t checked_multiply(
    std::size_t left,
    std::size_t right,
    const char* message) {
    if (left != 0U &&
        right >
            std::numeric_limits<std::size_t>::max() /
                left) {
        throw std::length_error(message);
    }
    return left * right;
}

[[nodiscard]] inline std::size_t checked_product3(
    std::size_t first,
    std::size_t second,
    std::size_t third,
    const char* message) {
    return checked_multiply(
        checked_multiply(
            first,
            second,
            message),
        third,
        message);
}

inline void require_local_capacity(
    std::size_t count) {
    constexpr auto local_max =
        std::numeric_limits<
            LocalIndex::value_type>::max();
    const std::uint64_t capacity =
        static_cast<std::uint64_t>(
            local_max) +
        1ULL;
    if (static_cast<std::uint64_t>(
            count) >
        capacity) {
        throw std::length_error(
            "mpmc::mesh::make_cartesian_topology_3d: entity count exceeds LocalIndex capacity");
    }
}

inline void require_csr_capacity(
    std::size_t entries) {
    if (entries >
        static_cast<std::size_t>(
            std::numeric_limits<
                CsrAdjacency::Offset>::max())) {
        throw std::length_error(
            "mpmc::mesh::make_cartesian_topology_3d: relation exceeds CSR offset capacity");
    }
}

[[nodiscard]] inline LocalIndex local_index(
    std::size_t value) {
    return LocalIndex{
        static_cast<
            LocalIndex::value_type>(
            value)};
}

[[nodiscard]] inline CsrAdjacency::Offset
csr_offset(std::size_t value) {
    return static_cast<
        CsrAdjacency::Offset>(value);
}

inline void require_axis(
    std::span<const double> axis,
    const char* message) {
    if (axis.size() < 2U ||
        !std::isfinite(axis.front())) {
        throw std::invalid_argument(message);
    }
    for (std::size_t index = 1U;
         index < axis.size();
         ++index) {
        if (!std::isfinite(axis[index])) {
            throw std::invalid_argument(
                message);
        }
        const double width =
            axis[index] -
            axis[index - 1U];
        if (!std::isfinite(width) ||
            width <= 0.0) {
            throw std::invalid_argument(
                message);
        }
    }
}

[[nodiscard]] inline double midpoint(
    double lower,
    double upper) {
    const double value =
        lower +
        0.5 * (upper - lower);
    if (!std::isfinite(value)) {
        throw std::invalid_argument(
            "mpmc::mesh::make_cartesian_mesh_3d: non-finite midpoint");
    }
    return value;
}

[[nodiscard]] inline double positive_product(
    double first,
    double second,
    const char* message) {
    const double value =
        first * second;
    if (!std::isfinite(value) ||
        value <= 0.0) {
        throw std::invalid_argument(message);
    }
    return value;
}

[[nodiscard]] inline double positive_product3(
    double first,
    double second,
    double third,
    const char* message) {
    const double value =
        first * second * third;
    if (!std::isfinite(value) ||
        value <= 0.0) {
        throw std::invalid_argument(message);
    }
    return value;
}

} // namespace cartesian_3d_detail

/// Build the topology of an nx-by-ny-by-nz Cartesian hexahedral grid.
///
/// Numbering is deterministic with I fastest, then J, then K. Vertex and cell
/// GlobalEntityId values are zero-based local logical numbers. Faces are grouped
/// by normal direction: x-normal faces first, then y-normal, then z-normal.
/// EntityKind::edge remains unmaterialized.
[[nodiscard]] inline Topology
make_cartesian_topology_3d(
    std::size_t nx,
    std::size_t ny,
    std::size_t nz) {
    using namespace cartesian_3d_detail;

    if (nx == 0U ||
        ny == 0U ||
        nz == 0U) {
        throw std::invalid_argument(
            "mpmc::mesh::make_cartesian_topology_3d: nx, ny and nz must be positive");
    }

    const std::size_t nxv =
        checked_add(
            nx,
            1U,
            "mpmc::mesh::make_cartesian_topology_3d: nx overflow");
    const std::size_t nyv =
        checked_add(
            ny,
            1U,
            "mpmc::mesh::make_cartesian_topology_3d: ny overflow");
    const std::size_t nzv =
        checked_add(
            nz,
            1U,
            "mpmc::mesh::make_cartesian_topology_3d: nz overflow");

    const std::size_t vertex_count =
        checked_product3(
            nxv,
            nyv,
            nzv,
            "mpmc::mesh::make_cartesian_topology_3d: vertex count overflow");
    const std::size_t cell_count =
        checked_product3(
            nx,
            ny,
            nz,
            "mpmc::mesh::make_cartesian_topology_3d: cell count overflow");
    const std::size_t x_face_count =
        checked_product3(
            nxv,
            ny,
            nz,
            "mpmc::mesh::make_cartesian_topology_3d: x-face count overflow");
    const std::size_t y_face_count =
        checked_product3(
            nx,
            nyv,
            nz,
            "mpmc::mesh::make_cartesian_topology_3d: y-face count overflow");
    const std::size_t z_face_count =
        checked_product3(
            nx,
            ny,
            nzv,
            "mpmc::mesh::make_cartesian_topology_3d: z-face count overflow");
    const std::size_t face_count =
        checked_add(
            checked_add(
                x_face_count,
                y_face_count,
                "mpmc::mesh::make_cartesian_topology_3d: face count overflow"),
            z_face_count,
            "mpmc::mesh::make_cartesian_topology_3d: face count overflow");

    require_local_capacity(vertex_count);
    require_local_capacity(face_count);
    require_local_capacity(cell_count);

    const std::size_t cell_vertex_entries =
        checked_multiply(
            cell_count,
            8U,
            "mpmc::mesh::make_cartesian_topology_3d: cell-to-vertex size overflow");
    const std::size_t cell_face_entries =
        checked_multiply(
            cell_count,
            6U,
            "mpmc::mesh::make_cartesian_topology_3d: cell-to-face size overflow");
    const std::size_t face_vertex_entries =
        checked_multiply(
            face_count,
            4U,
            "mpmc::mesh::make_cartesian_topology_3d: face-to-vertex size overflow");
    const std::size_t face_cell_capacity =
        checked_multiply(
            face_count,
            2U,
            "mpmc::mesh::make_cartesian_topology_3d: face-to-cell size overflow");

    require_csr_capacity(cell_vertex_entries);
    require_csr_capacity(cell_face_entries);
    require_csr_capacity(face_vertex_entries);
    require_csr_capacity(face_cell_capacity);

    const auto vertex =
        [nxv, nyv](
            std::size_t i,
            std::size_t j,
            std::size_t k) {
            return (k * nyv + j) *
                       nxv +
                   i;
        };
    const auto cell =
        [nx, ny](
            std::size_t i,
            std::size_t j,
            std::size_t k) {
            return (k * ny + j) *
                       nx +
                   i;
        };
    const auto x_face =
        [nxv, ny](
            std::size_t i,
            std::size_t j,
            std::size_t k) {
            return (k * ny + j) *
                       nxv +
                   i;
        };
    const auto y_face =
        [nx, nyv, x_face_count](
            std::size_t i,
            std::size_t j,
            std::size_t k) {
            return x_face_count +
                   (k * nyv + j) *
                       nx +
                   i;
        };
    const auto z_face =
        [nx, ny, x_face_count,
         y_face_count](
            std::size_t i,
            std::size_t j,
            std::size_t k) {
            return x_face_count +
                   y_face_count +
                   (k * ny + j) *
                       nx +
                   i;
        };

    Topology::EntityIds ids;
    ids.vertices.reserve(vertex_count);
    ids.faces.reserve(face_count);
    ids.cells.reserve(cell_count);
    for (std::size_t id = 0U;
         id < vertex_count;
         ++id) {
        ids.vertices.emplace_back(
            static_cast<
                GlobalEntityId::value_type>(
                id));
    }
    for (std::size_t id = 0U;
         id < face_count;
         ++id) {
        ids.faces.emplace_back(
            static_cast<
                GlobalEntityId::value_type>(
                id));
    }
    for (std::size_t id = 0U;
         id < cell_count;
         ++id) {
        ids.cells.emplace_back(
            static_cast<
                GlobalEntityId::value_type>(
                id));
    }

    std::vector<CsrAdjacency::Offset>
        cell_vertex_offsets{
            CsrAdjacency::Offset{0U}};
    std::vector<LocalIndex>
        cell_vertices;
    std::vector<CsrAdjacency::Offset>
        cell_face_offsets{
            CsrAdjacency::Offset{0U}};
    std::vector<LocalIndex>
        cell_faces;
    cell_vertex_offsets.reserve(
        cell_count + 1U);
    cell_face_offsets.reserve(
        cell_count + 1U);
    cell_vertices.reserve(
        cell_vertex_entries);
    cell_faces.reserve(
        cell_face_entries);

    for (std::size_t k = 0U;
         k < nz;
         ++k) {
        for (std::size_t j = 0U;
             j < ny;
             ++j) {
            for (std::size_t i = 0U;
                 i < nx;
                 ++i) {
                cell_vertices.push_back(
                    local_index(
                        vertex(i, j, k)));
                cell_vertices.push_back(
                    local_index(
                        vertex(i + 1U, j, k)));
                cell_vertices.push_back(
                    local_index(
                        vertex(
                            i + 1U,
                            j + 1U,
                            k)));
                cell_vertices.push_back(
                    local_index(
                        vertex(
                            i,
                            j + 1U,
                            k)));
                cell_vertices.push_back(
                    local_index(
                        vertex(
                            i,
                            j,
                            k + 1U)));
                cell_vertices.push_back(
                    local_index(
                        vertex(
                            i + 1U,
                            j,
                            k + 1U)));
                cell_vertices.push_back(
                    local_index(
                        vertex(
                            i + 1U,
                            j + 1U,
                            k + 1U)));
                cell_vertices.push_back(
                    local_index(
                        vertex(
                            i,
                            j + 1U,
                            k + 1U)));
                cell_vertex_offsets.push_back(
                    csr_offset(
                        cell_vertices.size()));

                cell_faces.push_back(
                    local_index(
                        x_face(i, j, k)));
                cell_faces.push_back(
                    local_index(
                        x_face(
                            i + 1U,
                            j,
                            k)));
                cell_faces.push_back(
                    local_index(
                        y_face(i, j, k)));
                cell_faces.push_back(
                    local_index(
                        y_face(
                            i,
                            j + 1U,
                            k)));
                cell_faces.push_back(
                    local_index(
                        z_face(i, j, k)));
                cell_faces.push_back(
                    local_index(
                        z_face(
                            i,
                            j,
                            k + 1U)));
                cell_face_offsets.push_back(
                    csr_offset(
                        cell_faces.size()));
            }
        }
    }

    std::vector<CsrAdjacency::Offset>
        face_vertex_offsets{
            CsrAdjacency::Offset{0U}};
    std::vector<LocalIndex>
        face_vertices;
    face_vertex_offsets.reserve(
        face_count + 1U);
    face_vertices.reserve(
        face_vertex_entries);

    for (std::size_t k = 0U;
         k < nz;
         ++k) {
        for (std::size_t j = 0U;
             j < ny;
             ++j) {
            for (std::size_t i = 0U;
                 i < nxv;
                 ++i) {
                face_vertices.push_back(
                    local_index(
                        vertex(i, j, k)));
                face_vertices.push_back(
                    local_index(
                        vertex(
                            i,
                            j + 1U,
                            k)));
                face_vertices.push_back(
                    local_index(
                        vertex(
                            i,
                            j + 1U,
                            k + 1U)));
                face_vertices.push_back(
                    local_index(
                        vertex(
                            i,
                            j,
                            k + 1U)));
                face_vertex_offsets.push_back(
                    csr_offset(
                        face_vertices.size()));
            }
        }
    }

    for (std::size_t k = 0U;
         k < nz;
         ++k) {
        for (std::size_t j = 0U;
             j < nyv;
             ++j) {
            for (std::size_t i = 0U;
                 i < nx;
                 ++i) {
                face_vertices.push_back(
                    local_index(
                        vertex(i, j, k)));
                face_vertices.push_back(
                    local_index(
                        vertex(
                            i,
                            j,
                            k + 1U)));
                face_vertices.push_back(
                    local_index(
                        vertex(
                            i + 1U,
                            j,
                            k + 1U)));
                face_vertices.push_back(
                    local_index(
                        vertex(
                            i + 1U,
                            j,
                            k)));
                face_vertex_offsets.push_back(
                    csr_offset(
                        face_vertices.size()));
            }
        }
    }

    for (std::size_t k = 0U;
         k < nzv;
         ++k) {
        for (std::size_t j = 0U;
             j < ny;
             ++j) {
            for (std::size_t i = 0U;
                 i < nx;
                 ++i) {
                face_vertices.push_back(
                    local_index(
                        vertex(i, j, k)));
                face_vertices.push_back(
                    local_index(
                        vertex(
                            i + 1U,
                            j,
                            k)));
                face_vertices.push_back(
                    local_index(
                        vertex(
                            i + 1U,
                            j + 1U,
                            k)));
                face_vertices.push_back(
                    local_index(
                        vertex(
                            i,
                            j + 1U,
                            k)));
                face_vertex_offsets.push_back(
                    csr_offset(
                        face_vertices.size()));
            }
        }
    }

    std::vector<CsrAdjacency::Offset>
        face_cell_offsets{
            CsrAdjacency::Offset{0U}};
    std::vector<LocalIndex>
        face_cells;
    face_cell_offsets.reserve(
        face_count + 1U);
    face_cells.reserve(
        face_cell_capacity);

    for (std::size_t k = 0U;
         k < nz;
         ++k) {
        for (std::size_t j = 0U;
             j < ny;
             ++j) {
            for (std::size_t i = 0U;
                 i < nxv;
                 ++i) {
                if (i > 0U) {
                    face_cells.push_back(
                        local_index(
                            cell(
                                i - 1U,
                                j,
                                k)));
                }
                if (i < nx) {
                    face_cells.push_back(
                        local_index(
                            cell(i, j, k)));
                }
                face_cell_offsets.push_back(
                    csr_offset(
                        face_cells.size()));
            }
        }
    }

    for (std::size_t k = 0U;
         k < nz;
         ++k) {
        for (std::size_t j = 0U;
             j < nyv;
             ++j) {
            for (std::size_t i = 0U;
                 i < nx;
                 ++i) {
                if (j > 0U) {
                    face_cells.push_back(
                        local_index(
                            cell(
                                i,
                                j - 1U,
                                k)));
                }
                if (j < ny) {
                    face_cells.push_back(
                        local_index(
                            cell(i, j, k)));
                }
                face_cell_offsets.push_back(
                    csr_offset(
                        face_cells.size()));
            }
        }
    }

    for (std::size_t k = 0U;
         k < nzv;
         ++k) {
        for (std::size_t j = 0U;
             j < ny;
             ++j) {
            for (std::size_t i = 0U;
                 i < nx;
                 ++i) {
                if (k > 0U) {
                    face_cells.push_back(
                        local_index(
                            cell(
                                i,
                                j,
                                k - 1U)));
                }
                if (k < nz) {
                    face_cells.push_back(
                        local_index(
                            cell(i, j, k)));
                }
                face_cell_offsets.push_back(
                    csr_offset(
                        face_cells.size()));
            }
        }
    }

    std::vector<CsrAdjacency> relations;
    relations.reserve(4U);
    relations.emplace_back(
        EntityKind::cell,
        EntityKind::vertex,
        vertex_count,
        std::move(
            cell_vertex_offsets),
        std::move(cell_vertices));
    relations.emplace_back(
        EntityKind::cell,
        EntityKind::face,
        face_count,
        std::move(
            cell_face_offsets),
        std::move(cell_faces));
    relations.emplace_back(
        EntityKind::face,
        EntityKind::vertex,
        vertex_count,
        std::move(
            face_vertex_offsets),
        std::move(face_vertices));
    relations.emplace_back(
        EntityKind::face,
        EntityKind::cell,
        cell_count,
        std::move(
            face_cell_offsets),
        std::move(face_cells));

    return Topology{
        std::move(ids),
        std::move(relations)};
}

/// Build a nonuniform Cartesian hexahedral mesh from strictly increasing axes.
[[nodiscard]] inline LinearMesh3D
make_cartesian_mesh_3d(
    std::span<const double>
        x_coordinates_m,
    std::span<const double>
        y_coordinates_m,
    std::span<const double>
        z_coordinates_m) {
    using namespace cartesian_3d_detail;

    require_axis(
        x_coordinates_m,
        "mpmc::mesh::make_cartesian_mesh_3d: x coordinates must be finite, strictly increasing, and contain at least two values");
    require_axis(
        y_coordinates_m,
        "mpmc::mesh::make_cartesian_mesh_3d: y coordinates must be finite, strictly increasing, and contain at least two values");
    require_axis(
        z_coordinates_m,
        "mpmc::mesh::make_cartesian_mesh_3d: z coordinates must be finite, strictly increasing, and contain at least two values");

    const std::size_t nx =
        x_coordinates_m.size() - 1U;
    const std::size_t ny =
        y_coordinates_m.size() - 1U;
    const std::size_t nz =
        z_coordinates_m.size() - 1U;
    Topology topology =
        make_cartesian_topology_3d(
            nx,
            ny,
            nz);

    const std::size_t nxv =
        x_coordinates_m.size();
    const std::size_t nyv =
        y_coordinates_m.size();
    const std::size_t x_face_count =
        checked_product3(
            nxv,
            ny,
            nz,
            "mpmc::mesh::make_cartesian_mesh_3d: x-face count overflow");
    const std::size_t y_face_count =
        checked_product3(
            nx,
            nyv,
            nz,
            "mpmc::mesh::make_cartesian_mesh_3d: y-face count overflow");

    const auto cell =
        [nx, ny](
            std::size_t i,
            std::size_t j,
            std::size_t k) {
            return (k * ny + j) *
                       nx +
                   i;
        };

    std::vector<Coordinate3D>
        coordinates;
    coordinates.reserve(
        topology.entity_count(
            EntityKind::vertex));
    for (std::size_t k = 0U;
         k <= nz;
         ++k) {
        for (std::size_t j = 0U;
             j <= ny;
             ++j) {
            for (std::size_t i = 0U;
                 i <= nx;
                 ++i) {
                coordinates.push_back(
                    Coordinate3D{
                        x_coordinates_m[i],
                        y_coordinates_m[j],
                        z_coordinates_m[k]});
            }
        }
    }

    std::vector<double> dx;
    std::vector<double> dy;
    std::vector<double> dz;
    dx.reserve(nx);
    dy.reserve(ny);
    dz.reserve(nz);
    for (std::size_t i = 0U;
         i < nx;
         ++i) {
        dx.push_back(
            x_coordinates_m[i + 1U] -
            x_coordinates_m[i]);
    }
    for (std::size_t j = 0U;
         j < ny;
         ++j) {
        dy.push_back(
            y_coordinates_m[j + 1U] -
            y_coordinates_m[j]);
    }
    for (std::size_t k = 0U;
         k < nz;
         ++k) {
        dz.push_back(
            z_coordinates_m[k + 1U] -
            z_coordinates_m[k]);
    }

    std::vector<double> cell_volumes;
    cell_volumes.reserve(
        topology.entity_count(
            EntityKind::cell));
    for (std::size_t k = 0U;
         k < nz;
         ++k) {
        for (std::size_t j = 0U;
             j < ny;
             ++j) {
            for (std::size_t i = 0U;
                 i < nx;
                 ++i) {
                cell_volumes.push_back(
                    positive_product3(
                        dx[i],
                        dy[j],
                        dz[k],
                        "mpmc::mesh::make_cartesian_mesh_3d: cell volume overflow"));
            }
        }
    }

    std::vector<Coordinate3D>
        face_centroids;
    std::vector<double>
        face_areas;
    std::vector<LocalIndex>
        face_owners;
    std::vector<UnitVector3D>
        face_normals;
    const std::size_t face_count =
        topology.entity_count(
            EntityKind::face);
    face_centroids.reserve(face_count);
    face_areas.reserve(face_count);
    face_owners.reserve(face_count);
    face_normals.reserve(face_count);

    for (std::size_t k = 0U;
         k < nz;
         ++k) {
        for (std::size_t j = 0U;
             j < ny;
             ++j) {
            for (std::size_t i = 0U;
                 i <= nx;
                 ++i) {
                face_centroids.push_back(
                    Coordinate3D{
                        x_coordinates_m[i],
                        midpoint(
                            y_coordinates_m[j],
                            y_coordinates_m[
                                j + 1U]),
                        midpoint(
                            z_coordinates_m[k],
                            z_coordinates_m[
                                k + 1U])});
                face_areas.push_back(
                    positive_product(
                        dy[j],
                        dz[k],
                        "mpmc::mesh::make_cartesian_mesh_3d: x-face area overflow"));
                if (i == 0U) {
                    face_owners.push_back(
                        local_index(
                            cell(0U, j, k)));
                    face_normals.push_back(
                        UnitVector3D{
                            -1.0, 0.0, 0.0});
                } else {
                    face_owners.push_back(
                        local_index(
                            cell(
                                i - 1U,
                                j,
                                k)));
                    face_normals.push_back(
                        UnitVector3D{
                            1.0, 0.0, 0.0});
                }
            }
        }
    }

    for (std::size_t k = 0U;
         k < nz;
         ++k) {
        for (std::size_t j = 0U;
             j <= ny;
             ++j) {
            for (std::size_t i = 0U;
                 i < nx;
                 ++i) {
                face_centroids.push_back(
                    Coordinate3D{
                        midpoint(
                            x_coordinates_m[i],
                            x_coordinates_m[
                                i + 1U]),
                        y_coordinates_m[j],
                        midpoint(
                            z_coordinates_m[k],
                            z_coordinates_m[
                                k + 1U])});
                face_areas.push_back(
                    positive_product(
                        dx[i],
                        dz[k],
                        "mpmc::mesh::make_cartesian_mesh_3d: y-face area overflow"));
                if (j == 0U) {
                    face_owners.push_back(
                        local_index(
                            cell(i, 0U, k)));
                    face_normals.push_back(
                        UnitVector3D{
                            0.0, -1.0, 0.0});
                } else {
                    face_owners.push_back(
                        local_index(
                            cell(
                                i,
                                j - 1U,
                                k)));
                    face_normals.push_back(
                        UnitVector3D{
                            0.0, 1.0, 0.0});
                }
            }
        }
    }

    for (std::size_t k = 0U;
         k <= nz;
         ++k) {
        for (std::size_t j = 0U;
             j < ny;
             ++j) {
            for (std::size_t i = 0U;
                 i < nx;
                 ++i) {
                face_centroids.push_back(
                    Coordinate3D{
                        midpoint(
                            x_coordinates_m[i],
                            x_coordinates_m[
                                i + 1U]),
                        midpoint(
                            y_coordinates_m[j],
                            y_coordinates_m[
                                j + 1U]),
                        z_coordinates_m[k]});
                face_areas.push_back(
                    positive_product(
                        dx[i],
                        dy[j],
                        "mpmc::mesh::make_cartesian_mesh_3d: z-face area overflow"));
                if (k == 0U) {
                    face_owners.push_back(
                        local_index(
                            cell(i, j, 0U)));
                    face_normals.push_back(
                        UnitVector3D{
                            0.0, 0.0, -1.0});
                } else {
                    face_owners.push_back(
                        local_index(
                            cell(
                                i,
                                j,
                                k - 1U)));
                    face_normals.push_back(
                        UnitVector3D{
                            0.0, 0.0, 1.0});
                }
            }
        }
    }

    const std::size_t z_face_count =
        checked_product3(
            nx,
            ny,
            z_coordinates_m.size(),
            "mpmc::mesh::make_cartesian_mesh_3d: z-face count overflow");
    if (face_centroids.size() !=
            face_count ||
        face_centroids.size() !=
            x_face_count +
                y_face_count +
                z_face_count) {
        throw std::logic_error(
            "mpmc::mesh::make_cartesian_mesh_3d: face metric count drift");
    }

    FaceGeometry3D face_geometry{
        topology.entity_count(
            EntityKind::cell),
        std::move(face_centroids),
        std::move(face_areas),
        std::move(face_owners),
        std::move(face_normals)};
    auto boundary =
        make_face_boundary_snapshot(
            topology);

    return LinearMesh3D{
        std::move(topology),
        std::move(coordinates),
        std::move(cell_volumes),
        std::move(face_geometry),
        std::move(boundary)};
}

} // namespace mpmc::mesh

#endif // MPMC_MESH_CARTESIAN_3D_HPP
