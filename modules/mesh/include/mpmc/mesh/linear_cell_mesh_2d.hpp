#ifndef MPMC_MESH_LINEAR_CELL_MESH_2D_HPP
#define MPMC_MESH_LINEAR_CELL_MESH_2D_HPP

#include <mpmc/mesh/face_boundary.hpp>
#include <mpmc/mesh/geometry_2d.hpp>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <map>
#include <numeric>
#include <optional>
#include <set>
#include <span>
#include <stdexcept>
#include <utility>
#include <vector>

namespace mpmc::mesh {

enum class LinearCellType2D : std::uint8_t { triangle, quadrilateral };

struct LinearCell2D {
    GlobalEntityId global_id;
    LinearCellType2D type;
    // Cyclic boundary order, either winding; never a sorted set of vertices.
    std::vector<LocalIndex> vertices;
};

struct LinearFaceAnnotation2D {
    std::array<LocalIndex, 2> vertices;
    GlobalEntityId global_id;
    PhysicalTag physical_tag;
};

/// Owning XY-plane mesh in SI units. Vertex/cell local order and cyclic cell
/// connectivity are retained. Faces use sorted stable vertex-ID keys, then
/// ascending face ID order; owner is the first incident input cell.
/// Triangle/strictly convex quad only; no coordinate welding or global overlap
/// search. Face means a line segment; separate edge entities are not created.
struct LinearMesh2D {
    Topology topology;
    Geometry2D geometry;
    FaceBoundarySnapshot face_boundary;
    std::vector<LinearCellType2D> cell_types;
};

namespace linear_cell_mesh_2d_detail {

[[nodiscard]] inline LocalIndex local_index(std::size_t value) {
    if (value > std::numeric_limits<LocalIndex::value_type>::max()) {
        throw std::length_error("mpmc::mesh::make_linear_mesh_2d: local index overflow");
    }
    return LocalIndex{static_cast<LocalIndex::value_type>(value)};
}

[[nodiscard]] inline CsrAdjacency::Offset offset(std::size_t value) {
    if (value > std::numeric_limits<CsrAdjacency::Offset>::max()) {
        throw std::length_error("mpmc::mesh::make_linear_mesh_2d: CSR overflow");
    }
    return static_cast<CsrAdjacency::Offset>(value);
}

struct PolygonMetric {
    Coordinate2D centroid;
    double area;
};

// Signed shoelace area and first moments (Green's theorem), evaluated in a
// translated/scaled frame. The relative threshold is dimensionless; an offset
// of the coordinate origin must not turn a valid cell into a degenerate cell.
[[nodiscard]] inline PolygonMetric polygon_metric(std::span<const Coordinate2D> points) {
    if (points.size() != 3U && points.size() != 4U) {
        throw std::invalid_argument("mpmc::mesh::make_linear_mesh_2d: triangle/quad required");
    }
    std::array<Coordinate2D, 4> p{};
    double scale = 0.0;
    const auto origin = points.front();
    for (std::size_t i = 0; i < points.size(); ++i) {
        p[i] = {points[i].x_m - origin.x_m, points[i].y_m - origin.y_m};
        if (!std::isfinite(p[i].x_m) || !std::isfinite(p[i].y_m)) {
            throw std::invalid_argument("mpmc::mesh::make_linear_mesh_2d: non-finite coordinate difference");
        }
        scale = std::max({scale, std::abs(p[i].x_m), std::abs(p[i].y_m)});
    }
    if (scale == 0.0) {
        throw std::invalid_argument("mpmc::mesh::make_linear_mesh_2d: degenerate cell");
    }
    for (auto& point : p) {
        point.x_m /= scale;
        point.y_m /= scale;
    }
    constexpr double tolerance = 256.0 * std::numeric_limits<double>::epsilon();
    double sign = 0.0;
    for (std::size_t i = 0; i < points.size(); ++i) {
        const auto a = p[i];
        const auto b = p[(i + 1U) % points.size()];
        const auto c = p[(i + 2U) % points.size()];
        const double turn = (b.x_m - a.x_m) * (c.y_m - a.y_m) -
                            (b.y_m - a.y_m) * (c.x_m - a.x_m);
        if (std::abs(turn) <= tolerance || (i != 0U && (turn > 0.0) != (sign > 0.0))) {
            throw std::invalid_argument("mpmc::mesh::make_linear_mesh_2d: degenerate, concave or self-intersecting cell");
        }
        sign = turn;
    }
    double twice_area = 0.0;
    double cx = 0.0;
    double cy = 0.0;
    for (std::size_t i = 0; i < points.size(); ++i) {
        const auto a = p[i];
        const auto b = p[(i + 1U) % points.size()];
        const double term = a.x_m * b.y_m - b.x_m * a.y_m;
        twice_area += term;
        cx += (a.x_m + b.x_m) * term;
        cy += (a.y_m + b.y_m) * term;
    }
    if (std::abs(twice_area) <= tolerance) {
        throw std::invalid_argument("mpmc::mesh::make_linear_mesh_2d: degenerate cell area");
    }
    const Coordinate2D centroid{origin.x_m + scale * (cx / (3.0 * twice_area)),
                                origin.y_m + scale * (cy / (3.0 * twice_area))};
    const double area = (0.5 * std::abs(twice_area) * scale) * scale;
    if (!std::isfinite(area) || area <= 0.0 ||
        !std::isfinite(centroid.x_m) || !std::isfinite(centroid.y_m)) {
        throw std::invalid_argument("mpmc::mesh::make_linear_mesh_2d: unrepresentable polygon metric");
    }
    return {centroid, area};
}

using EdgeKey = std::pair<std::uint64_t, std::uint64_t>;
struct FaceBuild {
    EdgeKey key;
    std::vector<LocalIndex> cells;
    std::optional<GlobalEntityId> id;
    PhysicalTag tag{0U};
};

} // namespace linear_cell_mesh_2d_detail

/// Construct a conforming linear 2D mesh. IDs are unique per entity kind,
/// including zero. Explicit annotations preserve face IDs; missing face IDs
/// are allocated above max(generated_face_id_floor, explicit face IDs).
/// Interior annotations must be untagged (the existing boundary contract).
/// Invalid geometry/annotations throw invalid_argument; connectivity outside
/// the vertex array throws out_of_range; capacity/ID exhaustion throws length_error.
[[nodiscard]] inline LinearMesh2D make_linear_mesh_2d(
    std::vector<GlobalEntityId> vertex_ids,
    std::vector<Coordinate2D> coordinates,
    std::span<const LinearCell2D> cells,
    std::span<const LinearFaceAnnotation2D> annotations = {},
    std::uint64_t generated_face_id_floor = 0U) {
    using namespace linear_cell_mesh_2d_detail;
    if (vertex_ids.size() != coordinates.size() || vertex_ids.empty() || cells.empty()) {
        throw std::invalid_argument("mpmc::mesh::make_linear_mesh_2d: nonempty matching vertex/cell arrays required");
    }
    (void)local_index(vertex_ids.size() - 1U);
    (void)local_index(cells.size() - 1U);
    std::map<std::uint64_t, LocalIndex> vertex_by_id;
    for (std::size_t i = 0; i < vertex_ids.size(); ++i) {
        if (!vertex_by_id.emplace(vertex_ids[i].value(), local_index(i)).second) {
            throw std::invalid_argument("mpmc::mesh::make_linear_mesh_2d: duplicate vertex ID");
        }
        if (!std::isfinite(coordinates[i].x_m) || !std::isfinite(coordinates[i].y_m)) {
            throw std::invalid_argument("mpmc::mesh::make_linear_mesh_2d: non-finite vertex coordinate");
        }
    }
    const auto key_for = [&](LocalIndex a, LocalIndex b) -> EdgeKey {
        if (a.value() >= vertex_ids.size() || b.value() >= vertex_ids.size()) {
            throw std::out_of_range("mpmc::mesh::make_linear_mesh_2d: vertex index out of range");
        }
        if (a == b) {
            throw std::invalid_argument("mpmc::mesh::make_linear_mesh_2d: repeated edge vertex");
        }
        const auto left = vertex_ids[a.value()].value();
        const auto right = vertex_ids[b.value()].value();
        return {std::min(left, right), std::max(left, right)};
    };
    Topology::EntityIds ids;
    std::set<std::uint64_t> cell_ids;
    std::vector<LinearCellType2D> cell_types;
    std::vector<Coordinate2D> centroids;
    std::vector<double> areas;
    std::map<EdgeKey, FaceBuild> by_edge;
    for (std::size_t i = 0; i < cells.size(); ++i) {
        const auto& cell = cells[i];
        const std::size_t width = cell.type == LinearCellType2D::triangle ? 3U :
            cell.type == LinearCellType2D::quadrilateral ? 4U : 0U;
        if (width == 0U || cell.vertices.size() != width ||
            !cell_ids.insert(cell.global_id.value()).second) {
            throw std::invalid_argument("mpmc::mesh::make_linear_mesh_2d: invalid cell type, width or duplicate ID");
        }
        std::set<LocalIndex::value_type> unique;
        std::array<Coordinate2D, 4> polygon{};
        for (std::size_t j = 0; j < width; ++j) {
            const auto vertex = cell.vertices[j].value();
            if (vertex >= coordinates.size()) {
                throw std::out_of_range("mpmc::mesh::make_linear_mesh_2d: cell vertex out of range");
            }
            if (!unique.insert(vertex).second) {
                throw std::invalid_argument("mpmc::mesh::make_linear_mesh_2d: repeated cell vertex");
            }
            polygon[j] = coordinates[vertex];
        }
        const auto metric = polygon_metric(std::span{polygon}.first(width));
        centroids.push_back(metric.centroid);
        areas.push_back(metric.area);
        ids.cells.push_back(cell.global_id);
        cell_types.push_back(cell.type);
        for (std::size_t j = 0; j < width; ++j) {
            const auto key = key_for(cell.vertices[j], cell.vertices[(j + 1U) % width]);
            auto [it, inserted] = by_edge.try_emplace(key, FaceBuild{key, {}, std::nullopt, PhysicalTag{0U}});
            (void)inserted;
            it->second.cells.push_back(local_index(i));
            if (it->second.cells.size() > 2U) {
                throw std::invalid_argument("mpmc::mesh::make_linear_mesh_2d: non-manifold face");
            }
        }
    }
    (void)local_index(by_edge.size() - 1U);
    std::set<std::uint64_t> annotated_ids;
    for (const auto& annotation : annotations) {
        const auto key = key_for(annotation.vertices[0], annotation.vertices[1]);
        auto it = by_edge.find(key);
        if (it == by_edge.end() || it->second.id.has_value() ||
            !annotated_ids.insert(annotation.global_id.value()).second) {
            throw std::invalid_argument("mpmc::mesh::make_linear_mesh_2d: missing or duplicate annotated face/ID");
        }
        if (it->second.cells.size() == 2U && annotation.physical_tag.is_tagged()) {
            throw std::invalid_argument("mpmc::mesh::make_linear_mesh_2d: tagged interior face");
        }
        it->second.id = annotation.global_id;
        it->second.tag = annotation.physical_tag;
        generated_face_id_floor = std::max(generated_face_id_floor, annotation.global_id.value());
    }
    std::vector<FaceBuild> faces;
    faces.reserve(by_edge.size());
    for (auto& [key, face] : by_edge) {
        (void)key;
        if (!face.id.has_value()) {
            if (generated_face_id_floor == std::numeric_limits<std::uint64_t>::max()) {
                throw std::length_error("mpmc::mesh::make_linear_mesh_2d: face ID exhaustion");
            }
            face.id = GlobalEntityId{++generated_face_id_floor};
        }
        faces.push_back(std::move(face));
    }
    std::sort(faces.begin(), faces.end(), [](const auto& a, const auto& b) {
        return a.id->value() < b.id->value();
    });
    std::map<EdgeKey, LocalIndex> face_by_edge;
    std::vector<CsrAdjacency::Offset> fv_offsets{0U}, fc_offsets{0U};
    std::vector<LocalIndex> fv, fc, owners;
    std::vector<Coordinate2D> face_centroids;
    std::vector<double> lengths;
    std::vector<UnitVector2D> normals;
    std::vector<PhysicalTag> tags;
    for (std::size_t i = 0; i < faces.size(); ++i) {
        const auto& face = faces[i];
        face_by_edge.emplace(face.key, local_index(i));
        ids.faces.push_back(*face.id);
        const auto first = vertex_by_id.at(face.key.first);
        const auto second = vertex_by_id.at(face.key.second);
        fv.insert(fv.end(), {first, second});
        fc.insert(fc.end(), face.cells.begin(), face.cells.end());
        fv_offsets.push_back(offset(fv.size()));
        fc_offsets.push_back(offset(fc.size()));
        const auto a = coordinates[first.value()];
        const auto b = coordinates[second.value()];
        const double dx = b.x_m - a.x_m;
        const double dy = b.y_m - a.y_m;
        const double length = std::hypot(dx, dy);
        if (!std::isfinite(length) || length <= 0.0) {
            throw std::invalid_argument("mpmc::mesh::make_linear_mesh_2d: invalid face length");
        }
        const Coordinate2D midpoint{std::midpoint(a.x_m, b.x_m), std::midpoint(a.y_m, b.y_m)};
        UnitVector2D normal{-dy / length, dx / length};
        const auto distance = [&](LocalIndex cell) {
            const auto c = centroids[cell.value()];
            return normal.x * (midpoint.x_m - c.x_m) + normal.y * (midpoint.y_m - c.y_m);
        };
        const double tolerance = 256.0 * std::numeric_limits<double>::epsilon() * length;
        const double owner_distance = distance(face.cells.front());
        if (!std::isfinite(owner_distance) || std::abs(owner_distance) <= tolerance) {
            throw std::invalid_argument("mpmc::mesh::make_linear_mesh_2d: cannot orient face normal");
        }
        if (owner_distance < 0.0) { normal.x = -normal.x; normal.y = -normal.y; }
        if (face.cells.size() == 2U) {
            const double neighbour_distance = distance(face.cells[1]);
            if (!std::isfinite(neighbour_distance) || neighbour_distance >= -tolerance) {
                throw std::invalid_argument("mpmc::mesh::make_linear_mesh_2d: adjacent cells must lie on opposite sides of face");
            }
        }
        face_centroids.push_back(midpoint);
        lengths.push_back(length);
        owners.push_back(face.cells.front());
        normals.push_back(normal);
        tags.push_back(face.tag);
    }
    std::vector<CsrAdjacency::Offset> cv_offsets{0U}, cf_offsets{0U};
    std::vector<LocalIndex> cv, cf;
    for (const auto& cell : cells) {
        cv.insert(cv.end(), cell.vertices.begin(), cell.vertices.end());
        for (std::size_t j = 0; j < cell.vertices.size(); ++j) {
            cf.push_back(face_by_edge.at(key_for(cell.vertices[j], cell.vertices[(j + 1U) % cell.vertices.size()])));
        }
        cv_offsets.push_back(offset(cv.size()));
        cf_offsets.push_back(offset(cf.size()));
    }
    std::vector<CsrAdjacency> relations;
    relations.emplace_back(EntityKind::cell, EntityKind::vertex, vertex_ids.size(), std::move(cv_offsets), std::move(cv));
    relations.emplace_back(EntityKind::cell, EntityKind::face, faces.size(), std::move(cf_offsets), std::move(cf));
    relations.emplace_back(EntityKind::face, EntityKind::vertex, vertex_ids.size(), std::move(fv_offsets), std::move(fv));
    relations.emplace_back(EntityKind::face, EntityKind::cell, cells.size(), std::move(fc_offsets), std::move(fc));
    ids.vertices = std::move(vertex_ids);
    Topology topology{std::move(ids), std::move(relations)};
    Geometry2D geometry{std::move(coordinates), std::move(centroids), std::move(areas),
        std::move(face_centroids), std::move(lengths), std::move(owners), std::move(normals)};
    auto boundary = make_face_boundary_snapshot(topology, tags);
    return {std::move(topology), std::move(geometry), std::move(boundary), std::move(cell_types)};
}

} // namespace mpmc::mesh

#endif
