#ifndef MPMC_MESH_CARTESIAN_1D_HPP
#define MPMC_MESH_CARTESIAN_1D_HPP

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

namespace cartesian_1d_detail {

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
            "mpmc::mesh::make_cartesian_topology_1d: entity count exceeds LocalIndex capacity");
    }
}

inline void require_csr_capacity(
    std::size_t entries) {
    if (entries >
        static_cast<std::size_t>(
            std::numeric_limits<
                CsrAdjacency::Offset>::max())) {
        throw std::length_error(
            "mpmc::mesh::make_cartesian_topology_1d: relation exceeds CSR offset capacity");
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

} // namespace cartesian_1d_detail

/// Build the topology of an nx-cell one-dimensional Cartesian grid.
///
/// Stable numbering is deterministic and zero-based within each entity kind:
/// vertex(i)=i, face(i)=i, cell(i)=i. In one dimension a topological face is a
/// point and therefore each face references exactly one vertex. Face and vertex
/// remain distinct EntityKind values even when they occupy the same coordinate.
[[nodiscard]] inline Topology
make_cartesian_topology_1d(
    std::size_t nx) {
    using namespace cartesian_1d_detail;

    if (nx == 0U) {
        throw std::invalid_argument(
            "mpmc::mesh::make_cartesian_topology_1d: nx must be positive");
    }

    const std::size_t vertex_count =
        checked_add(
            nx,
            1U,
            "mpmc::mesh::make_cartesian_topology_1d: vertex count overflow");
    const std::size_t face_count =
        vertex_count;
    const std::size_t cell_count =
        nx;

    require_local_capacity(vertex_count);
    require_local_capacity(face_count);
    require_local_capacity(cell_count);

    const std::size_t cell_entries =
        checked_multiply(
            cell_count,
            2U,
            "mpmc::mesh::make_cartesian_topology_1d: cell relation size overflow");
    const std::size_t face_vertex_entries =
        face_count;
    const std::size_t face_cell_entries =
        checked_multiply(
            cell_count,
            2U,
            "mpmc::mesh::make_cartesian_topology_1d: face-to-cell size overflow");

    require_csr_capacity(cell_entries);
    require_csr_capacity(face_vertex_entries);
    require_csr_capacity(face_cell_entries);

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
        cell_offsets{
            CsrAdjacency::Offset{0U}};
    std::vector<LocalIndex>
        cell_vertices;
    std::vector<LocalIndex>
        cell_faces;
    cell_offsets.reserve(cell_count + 1U);
    cell_vertices.reserve(cell_entries);
    cell_faces.reserve(cell_entries);
    for (std::size_t cell = 0U;
         cell < cell_count;
         ++cell) {
        cell_vertices.push_back(
            local_index(cell));
        cell_vertices.push_back(
            local_index(cell + 1U));
        cell_faces.push_back(
            local_index(cell));
        cell_faces.push_back(
            local_index(cell + 1U));
        cell_offsets.push_back(
            csr_offset(
                cell_vertices.size()));
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
    for (std::size_t face = 0U;
         face < face_count;
         ++face) {
        face_vertices.push_back(
            local_index(face));
        face_vertex_offsets.push_back(
            csr_offset(
                face_vertices.size()));
    }

    std::vector<CsrAdjacency::Offset>
        face_cell_offsets{
            CsrAdjacency::Offset{0U}};
    std::vector<LocalIndex>
        face_cells;
    face_cell_offsets.reserve(
        face_count + 1U);
    face_cells.reserve(
        face_cell_entries);
    for (std::size_t face = 0U;
         face < face_count;
         ++face) {
        if (face > 0U) {
            face_cells.push_back(
                local_index(face - 1U));
        }
        if (face < cell_count) {
            face_cells.push_back(
                local_index(face));
        }
        face_cell_offsets.push_back(
            csr_offset(
                face_cells.size()));
    }

    std::vector<CsrAdjacency> relations;
    relations.reserve(4U);
    relations.emplace_back(
        EntityKind::cell,
        EntityKind::vertex,
        vertex_count,
        cell_offsets,
        std::move(cell_vertices));
    relations.emplace_back(
        EntityKind::cell,
        EntityKind::face,
        face_count,
        std::move(cell_offsets),
        std::move(cell_faces));
    relations.emplace_back(
        EntityKind::face,
        EntityKind::vertex,
        vertex_count,
        std::move(face_vertex_offsets),
        std::move(face_vertices));
    relations.emplace_back(
        EntityKind::face,
        EntityKind::cell,
        cell_count,
        std::move(face_cell_offsets),
        std::move(face_cells));

    return Topology{
        std::move(ids),
        std::move(relations)};
}

} // namespace mpmc::mesh

#endif // MPMC_MESH_CARTESIAN_1D_HPP
