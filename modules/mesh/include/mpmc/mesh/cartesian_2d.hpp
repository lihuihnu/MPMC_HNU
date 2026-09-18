#ifndef MPMC_MESH_CARTESIAN_2D_HPP
#define MPMC_MESH_CARTESIAN_2D_HPP

#include <mpmc/mesh/csr_adjacency.hpp>
#include <mpmc/mesh/entity.hpp>
#include <mpmc/mesh/topology.hpp>

#include <cstddef>
#include <cstdint>
#include <limits>
#include <stdexcept>
#include <utility>
#include <vector>

namespace mpmc::mesh {

namespace detail {

[[nodiscard]] inline std::size_t cartesian_checked_add(std::size_t left,
                                                        std::size_t right,
                                                        const char* message) {
    if (left > std::numeric_limits<std::size_t>::max() - right) {
        throw std::length_error(message);
    }
    return left + right;
}

[[nodiscard]] inline std::size_t cartesian_checked_multiply(std::size_t left,
                                                             std::size_t right,
                                                             const char* message) {
    if (left != 0U && right > std::numeric_limits<std::size_t>::max() / left) {
        throw std::length_error(message);
    }
    return left * right;
}

inline void require_local_capacity(std::size_t count) {
    constexpr auto local_max = std::numeric_limits<LocalIndex::value_type>::max();
    const std::uint64_t local_capacity = static_cast<std::uint64_t>(local_max) + 1ULL;
    if (static_cast<std::uint64_t>(count) > local_capacity) {
        throw std::length_error(
            "mpmc::mesh::make_cartesian_topology_2d: entity count exceeds LocalIndex capacity");
    }
}

inline void require_csr_capacity(std::size_t entries) {
    if (entries > static_cast<std::size_t>(std::numeric_limits<CsrAdjacency::Offset>::max())) {
        throw std::length_error(
            "mpmc::mesh::make_cartesian_topology_2d: relation exceeds CSR offset capacity");
    }
}

[[nodiscard]] inline LocalIndex local_index(std::size_t value) {
    return LocalIndex{static_cast<LocalIndex::value_type>(value)};
}

[[nodiscard]] inline CsrAdjacency::Offset csr_offset(std::size_t value) {
    return static_cast<CsrAdjacency::Offset>(value);
}

} // namespace detail

/// Build the topology of an nx-by-ny logical Cartesian cell grid.
///
/// This function creates topology only: no coordinates, metric geometry, boundary
/// tags, fields, DoFs or partition metadata are attached.
///
/// Stable local/global numbering is deterministic and zero-based within each kind:
/// - vertex(i,j) = j * (nx + 1) + i;
/// - cell(i,j) = j * nx + i;
/// - vertical faces are numbered first, row-major in (j,i);
/// - horizontal faces follow, row-major in (j,i).
///
/// In two dimensions, EntityKind::face represents the codimension-one interfaces
/// consumed by cell-centered discretizations. EntityKind::edge is left empty.
[[nodiscard]] inline Topology make_cartesian_topology_2d(std::size_t nx, std::size_t ny) {
    if (nx == 0U || ny == 0U) {
        throw std::invalid_argument(
            "mpmc::mesh::make_cartesian_topology_2d: nx and ny must be positive");
    }

    const std::size_t nx_vertices =
        detail::cartesian_checked_add(nx, 1U,
                                      "mpmc::mesh::make_cartesian_topology_2d: nx overflow");
    const std::size_t ny_vertices =
        detail::cartesian_checked_add(ny, 1U,
                                      "mpmc::mesh::make_cartesian_topology_2d: ny overflow");

    const std::size_t vertex_count =
        detail::cartesian_checked_multiply(
            nx_vertices, ny_vertices,
            "mpmc::mesh::make_cartesian_topology_2d: vertex count overflow");
    const std::size_t cell_count =
        detail::cartesian_checked_multiply(
            nx, ny,
            "mpmc::mesh::make_cartesian_topology_2d: cell count overflow");
    const std::size_t vertical_face_count =
        detail::cartesian_checked_multiply(
            nx_vertices, ny,
            "mpmc::mesh::make_cartesian_topology_2d: vertical face count overflow");
    const std::size_t horizontal_face_count =
        detail::cartesian_checked_multiply(
            nx, ny_vertices,
            "mpmc::mesh::make_cartesian_topology_2d: horizontal face count overflow");
    const std::size_t face_count =
        detail::cartesian_checked_add(
            vertical_face_count, horizontal_face_count,
            "mpmc::mesh::make_cartesian_topology_2d: face count overflow");

    detail::require_local_capacity(vertex_count);
    detail::require_local_capacity(face_count);
    detail::require_local_capacity(cell_count);

    const std::size_t cell_relation_entries =
        detail::cartesian_checked_multiply(
            cell_count, 4U,
            "mpmc::mesh::make_cartesian_topology_2d: cell relation size overflow");
    const std::size_t face_vertex_entries =
        detail::cartesian_checked_multiply(
            face_count, 2U,
            "mpmc::mesh::make_cartesian_topology_2d: face-to-vertex size overflow");
    const std::size_t boundary_face_count =
        detail::cartesian_checked_multiply(
            detail::cartesian_checked_add(
                nx, ny,
                "mpmc::mesh::make_cartesian_topology_2d: boundary count overflow"),
            2U,
            "mpmc::mesh::make_cartesian_topology_2d: boundary count overflow");
    const std::size_t doubled_face_count =
        detail::cartesian_checked_multiply(
            face_count, 2U,
            "mpmc::mesh::make_cartesian_topology_2d: face-to-cell size overflow");
    if (boundary_face_count > doubled_face_count) {
        throw std::length_error(
            "mpmc::mesh::make_cartesian_topology_2d: invalid face-to-cell size");
    }
    const std::size_t face_cell_entries = doubled_face_count - boundary_face_count;

    detail::require_csr_capacity(cell_relation_entries);
    detail::require_csr_capacity(face_vertex_entries);
    detail::require_csr_capacity(face_cell_entries);

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

    Topology::EntityIds ids;
    ids.vertices.reserve(vertex_count);
    ids.faces.reserve(face_count);
    ids.cells.reserve(cell_count);
    for (std::size_t id = 0; id < vertex_count; ++id) {
        ids.vertices.emplace_back(static_cast<GlobalEntityId::value_type>(id));
    }
    for (std::size_t id = 0; id < face_count; ++id) {
        ids.faces.emplace_back(static_cast<GlobalEntityId::value_type>(id));
    }
    for (std::size_t id = 0; id < cell_count; ++id) {
        ids.cells.emplace_back(static_cast<GlobalEntityId::value_type>(id));
    }

    std::vector<CsrAdjacency::Offset> cell_offsets;
    std::vector<LocalIndex> cell_vertices;
    std::vector<LocalIndex> cell_faces;
    cell_offsets.reserve(cell_count + 1U);
    cell_vertices.reserve(cell_relation_entries);
    cell_faces.reserve(cell_relation_entries);
    cell_offsets.push_back(CsrAdjacency::Offset{0U});

    for (std::size_t j = 0; j < ny; ++j) {
        for (std::size_t i = 0; i < nx; ++i) {
            // Canonical logical order: lower-left, lower-right, upper-right, upper-left.
            cell_vertices.push_back(detail::local_index(vertex(i, j)));
            cell_vertices.push_back(detail::local_index(vertex(i + 1U, j)));
            cell_vertices.push_back(detail::local_index(vertex(i + 1U, j + 1U)));
            cell_vertices.push_back(detail::local_index(vertex(i, j + 1U)));

            // Canonical logical order: left, right, bottom, top.
            cell_faces.push_back(detail::local_index(vertical_face(i, j)));
            cell_faces.push_back(detail::local_index(vertical_face(i + 1U, j)));
            cell_faces.push_back(detail::local_index(horizontal_face(i, j)));
            cell_faces.push_back(detail::local_index(horizontal_face(i, j + 1U)));

            cell_offsets.push_back(detail::csr_offset(cell_vertices.size()));
        }
    }

    std::vector<CsrAdjacency::Offset> face_vertex_offsets;
    std::vector<LocalIndex> face_vertices;
    face_vertex_offsets.reserve(face_count + 1U);
    face_vertices.reserve(face_vertex_entries);
    face_vertex_offsets.push_back(CsrAdjacency::Offset{0U});

    for (std::size_t j = 0; j < ny; ++j) {
        for (std::size_t i = 0; i < nx_vertices; ++i) {
            face_vertices.push_back(detail::local_index(vertex(i, j)));
            face_vertices.push_back(detail::local_index(vertex(i, j + 1U)));
            face_vertex_offsets.push_back(detail::csr_offset(face_vertices.size()));
        }
    }
    for (std::size_t j = 0; j < ny_vertices; ++j) {
        for (std::size_t i = 0; i < nx; ++i) {
            face_vertices.push_back(detail::local_index(vertex(i, j)));
            face_vertices.push_back(detail::local_index(vertex(i + 1U, j)));
            face_vertex_offsets.push_back(detail::csr_offset(face_vertices.size()));
        }
    }

    std::vector<CsrAdjacency::Offset> face_cell_offsets;
    std::vector<LocalIndex> face_cells;
    face_cell_offsets.reserve(face_count + 1U);
    face_cells.reserve(face_cell_entries);
    face_cell_offsets.push_back(CsrAdjacency::Offset{0U});

    for (std::size_t j = 0; j < ny; ++j) {
        for (std::size_t i = 0; i < nx_vertices; ++i) {
            if (i > 0U) {
                face_cells.push_back(detail::local_index(cell(i - 1U, j)));
            }
            if (i < nx) {
                face_cells.push_back(detail::local_index(cell(i, j)));
            }
            face_cell_offsets.push_back(detail::csr_offset(face_cells.size()));
        }
    }
    for (std::size_t j = 0; j < ny_vertices; ++j) {
        for (std::size_t i = 0; i < nx; ++i) {
            if (j > 0U) {
                face_cells.push_back(detail::local_index(cell(i, j - 1U)));
            }
            if (j < ny) {
                face_cells.push_back(detail::local_index(cell(i, j)));
            }
            face_cell_offsets.push_back(detail::csr_offset(face_cells.size()));
        }
    }

    std::vector<CsrAdjacency> relations;
    relations.reserve(4U);
    relations.emplace_back(EntityKind::cell, EntityKind::vertex, vertex_count,
                           cell_offsets, std::move(cell_vertices));
    relations.emplace_back(EntityKind::cell, EntityKind::face, face_count,
                           std::move(cell_offsets), std::move(cell_faces));
    relations.emplace_back(EntityKind::face, EntityKind::vertex, vertex_count,
                           std::move(face_vertex_offsets), std::move(face_vertices));
    relations.emplace_back(EntityKind::face, EntityKind::cell, cell_count,
                           std::move(face_cell_offsets), std::move(face_cells));

    return Topology{std::move(ids), std::move(relations)};
}

} // namespace mpmc::mesh

#endif // MPMC_MESH_CARTESIAN_2D_HPP
