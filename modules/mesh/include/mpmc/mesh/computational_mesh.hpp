#ifndef MPMC_MESH_COMPUTATIONAL_MESH_HPP
#define MPMC_MESH_COMPUTATIONAL_MESH_HPP

#include <mpmc/mesh/linear_cell_geometric_operator_3d.hpp>
#include <mpmc/mesh/linear_cell_mesh_2d.hpp>
#include <mpmc/mesh/linear_cell_mesh_3d.hpp>
#include <mpmc/mesh/mesh_exchange_document.hpp>

#include <algorithm>
#include <array>
#include <cstddef>
#include <cstdint>
#include <map>
#include <span>
#include <stdexcept>
#include <utility>
#include <vector>

namespace mpmc::mesh {
namespace computational_mesh_detail {

template <typename T>
[[nodiscard]] inline std::vector<T> copy(std::span<const T> values) {
    return {values.begin(), values.end()};
}

[[nodiscard]] inline LocalIndex local(std::size_t value) {
    // Topology has already checked LocalIndex capacity for each entity kind.
    return LocalIndex{static_cast<LocalIndex::value_type>(value)};
}

inline void require_domain(const MeshExchangeDocument& document, int dimension) {
    if (document.dimension() != dimension) {
        throw std::invalid_argument("mpmc::mesh::prepare_linear_mesh: dimension mismatch");
    }
    const auto& topology = document.topology();
    if (topology.entity_count(EntityKind::edge) != 0U || topology.relation_count() != 4U ||
        !topology.has_relation(EntityKind::cell, EntityKind::vertex) ||
        !topology.has_relation(EntityKind::cell, EntityKind::face) ||
        !topology.has_relation(EntityKind::face, EntityKind::vertex) ||
        !topology.has_relation(EntityKind::face, EntityKind::cell)) {
        throw std::invalid_argument(
            "mpmc::mesh::prepare_linear_mesh: requires the four cell/face/vertex relations and no explicit edges or additional relations");
    }
}

[[nodiscard]] inline FaceBoundarySnapshot source_boundary(
    const MeshExchangeDocument& document) {
    const auto count = document.topology().entity_count(EntityKind::face);
    const auto tags = document.face_boundary() ?
        copy(document.face_boundary()->physical_tags()) :
        std::vector<PhysicalTag>(count, PhysicalTag{0U});
    auto derived = make_face_boundary_snapshot(document.topology(), tags);
    if (document.face_boundary()) {
        const auto source = document.face_boundary()->classifications();
        if (!std::equal(source.begin(), source.end(), derived.classifications().begin())) {
            throw std::invalid_argument(
                "mpmc::mesh::prepare_linear_mesh: face classification contradicts source topology");
        }
    }
    return derived;
}

[[nodiscard]] inline std::vector<GlobalEntityId::value_type> adjacency_ids(
    const Topology& topology, EntityKind source, EntityKind target, LocalIndex entity) {
    std::vector<GlobalEntityId::value_type> result;
    for (const auto adjacent : topology.relation(source, target).adjacent(entity)) {
        result.push_back(topology.global_id(target, adjacent).value());
    }
    std::sort(result.begin(), result.end());
    if (std::adjacent_find(result.begin(), result.end()) != result.end()) {
        throw std::invalid_argument(
            "mpmc::mesh::prepare_linear_mesh: repeated entity in source/rebuilt adjacency");
    }
    return result;
}

[[nodiscard]] inline bool same_face_cycle(
    const Topology& source, LocalIndex face, const Topology& rebuilt, LocalIndex other) {
    const auto before = source.relation(EntityKind::face, EntityKind::vertex).adjacent(face);
    const auto after = rebuilt.relation(EntityKind::face, EntityKind::vertex).adjacent(other);
    if (before.size() != after.size() || before.empty()) { return false; }
    // A face ring may start at any vertex and have either winding. A different
    // permutation (e.g. a bow-tie quad) must not be silently canonicalized.
    for (std::size_t start = 0; start < after.size(); ++start) {
        bool forward = true;
        bool reverse = true;
        for (std::size_t i = 0; i < before.size(); ++i) {
            const auto id = source.global_id(EntityKind::vertex, before[i]);
            forward = forward && id == rebuilt.global_id(
                EntityKind::vertex, after[(start + i) % after.size()]);
            reverse = reverse && id == rebuilt.global_id(
                EntityKind::vertex, after[(start + after.size() - i) % after.size()]);
        }
        if (forward || reverse) { return true; }
    }
    return false;
}

// Return original-local face -> rebuilt-local face. Cells and vertices keep
// input order in both builders. No repaired/renumbered topology is substituted.
[[nodiscard]] inline std::vector<LocalIndex> verify_topology(
    const Topology& source, const Topology& rebuilt) {
    for (const auto kind : {EntityKind::vertex, EntityKind::cell}) {
        const auto before = source.global_ids(kind);
        const auto after = rebuilt.global_ids(kind);
        if (before.size() != after.size() ||
            !std::equal(before.begin(), before.end(), after.begin())) {
            throw std::invalid_argument(
                "mpmc::mesh::prepare_linear_mesh: rebuilt vertex/cell ordering changed");
        }
    }
    const auto face_ids = source.global_ids(EntityKind::face);
    const auto rebuilt_ids = rebuilt.global_ids(EntityKind::face);
    if (face_ids.size() != rebuilt_ids.size()) {
        throw std::invalid_argument(
            "mpmc::mesh::prepare_linear_mesh: source is missing or has additional cell faces");
    }
    std::map<GlobalEntityId::value_type, LocalIndex> rebuilt_face_by_id;
    for (std::size_t face = 0; face < rebuilt_ids.size(); ++face) {
        rebuilt_face_by_id.emplace(rebuilt_ids[face].value(), local(face));
    }
    std::vector<LocalIndex> mapping;
    mapping.reserve(face_ids.size());
    for (std::size_t face = 0; face < face_ids.size(); ++face) {
        const auto found = rebuilt_face_by_id.find(face_ids[face].value());
        if (found == rebuilt_face_by_id.end()) {
            throw std::invalid_argument(
                "mpmc::mesh::prepare_linear_mesh: rebuilt face identity changed");
        }
        mapping.push_back(found->second);
        if (!same_face_cycle(source, local(face), rebuilt, found->second)) {
            throw std::invalid_argument(
                "mpmc::mesh::prepare_linear_mesh: source face vertex cycle contradicts cell boundary");
        }
        for (const auto target : {EntityKind::vertex, EntityKind::cell}) {
            if (adjacency_ids(source, EntityKind::face, target, local(face)) !=
                adjacency_ids(rebuilt, EntityKind::face, target, found->second)) {
                throw std::invalid_argument(
                    "mpmc::mesh::prepare_linear_mesh: face incidence contradicts cell connectivity");
            }
        }
    }
    for (std::size_t cell = 0; cell < source.entity_count(EntityKind::cell); ++cell) {
        for (const auto target : {EntityKind::vertex, EntityKind::face}) {
            if (adjacency_ids(source, EntityKind::cell, target, local(cell)) !=
                adjacency_ids(rebuilt, EntityKind::cell, target, local(cell))) {
                throw std::invalid_argument(
                    "mpmc::mesh::prepare_linear_mesh: cell incidence contradicts rebuilt topology");
            }
        }
    }
    return mapping;
}

[[nodiscard]] inline LinearCellType3D cell_type_3d(std::size_t vertices) {
    switch (vertices) {
    case 4U: return LinearCellType3D::tetrahedron;
    case 5U: return LinearCellType3D::pyramid;
    case 6U: return LinearCellType3D::wedge;
    case 8U: return LinearCellType3D::hexahedron;
    default:
        throw std::invalid_argument(
            "mpmc::mesh::prepare_linear_mesh_3d: unsupported linear cell");
    }
}

} // namespace computational_mesh_detail

/// Validate an exchange document as a conforming, planar, linear triangle/quad
/// computational mesh. Rebuild with the format-independent builder, then reject
/// every disagreement with the source cell/face incidence and classification.
/// The returned topology, entity IDs and all local orders are exactly the
/// source's, including face owner order; geometry is remapped by stable face ID.
/// Thus fields/groups kept in document remain aligned without copying them or
/// inventing physics. This baseline requires precisely the four builder-owned
/// relations (cell->vertex/face and face->vertex/cell), with no explicit edges.
/// It checks neither global intersections between disconnected cells nor any
/// discretization-specific admissibility. Input/output coordinates use SI.
[[nodiscard]] inline LinearMesh2D prepare_linear_mesh_2d(
    const MeshExchangeDocument& document) {
    using namespace computational_mesh_detail;
    require_domain(document, 2);
    const auto& source = document.topology();
    auto boundary = source_boundary(document);
    std::vector<Coordinate2D> coordinates;
    for (const auto point : document.vertex_coordinates_m()) {
        coordinates.push_back({point.x_m, point.y_m});
    }
    std::vector<LinearCell2D> cells;
    for (std::size_t cell = 0; cell < source.entity_count(EntityKind::cell); ++cell) {
        const auto vertices = source.relation(EntityKind::cell, EntityKind::vertex).adjacent(local(cell));
        if (vertices.size() != 3U && vertices.size() != 4U) {
            throw std::invalid_argument("mpmc::mesh::prepare_linear_mesh_2d: triangle/quad required");
        }
        cells.push_back({source.global_id(EntityKind::cell, local(cell)),
                         vertices.size() == 3U ? LinearCellType2D::triangle : LinearCellType2D::quadrilateral,
                         copy(vertices)});
    }
    std::vector<LinearFaceAnnotation2D> annotations;
    for (std::size_t face = 0; face < source.entity_count(EntityKind::face); ++face) {
        const auto vertices = source.relation(EntityKind::face, EntityKind::vertex).adjacent(local(face));
        if (vertices.size() != 2U) {
            throw std::invalid_argument("mpmc::mesh::prepare_linear_mesh_2d: two-vertex face required");
        }
        annotations.push_back({{vertices[0], vertices[1]},
            source.global_id(EntityKind::face, local(face)), boundary.physical_tag(local(face))});
    }
    auto rebuilt = make_linear_mesh_2d(copy(source.global_ids(EntityKind::vertex)),
                                     coordinates, cells, annotations);
    const auto mapping = verify_topology(source, rebuilt.topology);
    std::vector<Coordinate2D> centroids;
    std::vector<double> lengths;
    std::vector<LocalIndex> owners;
    std::vector<UnitVector2D> normals;
    for (std::size_t face = 0; face < mapping.size(); ++face) {
        const auto other = mapping[face];
        const auto owner = source.relation(EntityKind::face, EntityKind::cell).adjacent(local(face)).front();
        auto normal = rebuilt.geometry.face_owner_unit_normal(other);
        if (owner != rebuilt.geometry.face_owner(other)) {
            normal = {-normal.x, -normal.y};
        }
        centroids.push_back(rebuilt.geometry.face_centroid_m(other));
        lengths.push_back(rebuilt.geometry.face_length_m(other));
        owners.push_back(owner);
        normals.push_back(normal);
    }
    Geometry2D geometry{std::move(coordinates), copy(rebuilt.geometry.cell_centroids_m()),
        copy(rebuilt.geometry.cell_areas_m2()), std::move(centroids), std::move(lengths),
        std::move(owners), std::move(normals)};
    return {source, std::move(geometry), std::move(boundary), std::move(rebuilt.cell_types)};
}

/// Three-dimensional counterpart of prepare_linear_mesh_2d for linear
/// tetrahedra/hexahedra/wedges/pyramids with conforming planar faces. It also
/// validates positive two-sided distances for the builder's vertex-mean cell
/// reference points through the common geometric operator. These reference
/// points are not generally volume centroids; no TPFA suitability is promised.
[[nodiscard]] inline LinearMesh3D prepare_linear_mesh_3d(
    const MeshExchangeDocument& document) {
    using namespace computational_mesh_detail;
    require_domain(document, 3);
    const auto& source = document.topology();
    auto boundary = source_boundary(document);
    std::vector<LinearCell3D> cells;
    for (std::size_t cell = 0; cell < source.entity_count(EntityKind::cell); ++cell) {
        const auto vertices = source.relation(EntityKind::cell, EntityKind::vertex).adjacent(local(cell));
        cells.push_back({source.global_id(EntityKind::cell, local(cell)),
                         cell_type_3d(vertices.size()), copy(vertices)});
    }
    std::vector<LinearFaceAnnotation3D> annotations;
    for (std::size_t face = 0; face < source.entity_count(EntityKind::face); ++face) {
        annotations.push_back({
            copy(source.relation(EntityKind::face, EntityKind::vertex).adjacent(local(face))),
            source.global_id(EntityKind::face, local(face)), boundary.physical_tag(local(face))});
    }
    auto rebuilt = make_linear_mesh_3d(copy(source.global_ids(EntityKind::vertex)),
                                     document.vertex_coordinates_m(), cells, annotations);
    const auto mapping = verify_topology(source, rebuilt.topology);
    std::vector<Coordinate3D> centroids;
    std::vector<double> areas;
    std::vector<LocalIndex> owners;
    std::vector<UnitVector3D> normals;
    for (std::size_t face = 0; face < mapping.size(); ++face) {
        const auto other = mapping[face];
        const auto owner = source.relation(EntityKind::face, EntityKind::cell).adjacent(local(face)).front();
        auto normal = rebuilt.face_geometry.face_owner_unit_normal(other);
        if (owner != rebuilt.face_geometry.face_owner(other)) {
            normal = {-normal.x, -normal.y, -normal.z};
        }
        centroids.push_back(rebuilt.face_geometry.face_centroid_m(other));
        areas.push_back(rebuilt.face_geometry.face_area_m2(other));
        owners.push_back(owner);
        normals.push_back(normal);
    }
    FaceGeometry3D geometry{cells.size(), std::move(centroids), std::move(areas),
                            std::move(owners), std::move(normals)};
    LinearMesh3D result{source, copy(document.vertex_coordinates_m()),
        std::move(rebuilt.cell_volumes_m3), std::move(geometry), std::move(boundary)};
    (void)make_cell_face_geometric_operator_3d(result);
    return result;
}

} // namespace mpmc::mesh
#endif // MPMC_MESH_COMPUTATIONAL_MESH_HPP
