#ifndef MPMC_MESH_ACTIVE_CORNER_POINT_HPP
#define MPMC_MESH_ACTIVE_CORNER_POINT_HPP

#include <mpmc/mesh/face_geometry_3d.hpp>
#include <mpmc/mesh/grdecl.hpp>

#include <algorithm>
#include <array>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <iterator>
#include <map>
#include <set>
#include <span>
#include <stdexcept>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

namespace mpmc::mesh {

struct ActiveCornerPointGrid {
    Topology topology;
    std::vector<Coordinate3D> vertex_coordinates_m;
    std::vector<double> cell_volumes_m3;
    FaceGeometry3D face_geometry;
    std::vector<DenseFieldSnapshot> cell_fields;
    std::vector<GlobalEntityId> source_logical_cell_ids;

    [[nodiscard]] std::size_t cell_count() const noexcept {
        return source_logical_cell_ids.size();
    }

    [[nodiscard]] GlobalEntityId source_logical_cell_id(
        LocalIndex processed_cell) const {
        const std::size_t local =
            static_cast<std::size_t>(
                processed_cell.value());
        if (local >= source_logical_cell_ids.size()) {
            throw std::out_of_range(
                "mpmc::mesh::ActiveCornerPointGrid: processed cell index out of range");
        }
        return source_logical_cell_ids[local];
    }
};

namespace active_corner_point_detail {

struct CoordinateKey {
    double x;
    double y;
    double z;

    [[nodiscard]] friend bool operator<(
        const CoordinateKey& left,
        const CoordinateKey& right) noexcept {
        return std::tie(left.x, left.y, left.z) <
               std::tie(right.x, right.y, right.z);
    }

    [[nodiscard]] friend bool operator==(
        const CoordinateKey& left,
        const CoordinateKey& right) noexcept {
        return left.x == right.x &&
               left.y == right.y &&
               left.z == right.z;
    }
};

[[nodiscard]] inline CoordinateKey coordinate_key(
    Coordinate3D coordinate) {
    return CoordinateKey{
        coordinate.x_m,
        coordinate.y_m,
        coordinate.z_m};
}

struct FaceKey {
    std::array<std::size_t, 4> vertices;

    [[nodiscard]] friend bool operator<(
        const FaceKey& left,
        const FaceKey& right) noexcept {
        return left.vertices < right.vertices;
    }
};

struct FaceOccurrence {
    std::array<std::size_t, 4> ordered_vertices;
    std::size_t processed_cell;
    std::size_t source_logical_cell;
    std::size_t face_slot;
};

struct FaceBuild {
    FaceKey key;
    std::array<std::size_t, 4> ordered_vertices;
    std::vector<std::size_t> processed_cells;
    std::vector<std::size_t> source_logical_cells;
};

[[nodiscard]] inline LocalIndex checked_local(
    std::size_t value,
    const char* message) {
    if (value >
        static_cast<std::size_t>(
            std::numeric_limits<
                LocalIndex::value_type>::max())) {
        throw std::length_error(message);
    }
    return LocalIndex{
        static_cast<
            LocalIndex::value_type>(value)};
}

[[nodiscard]] inline CsrAdjacency::Offset checked_offset(
    std::size_t value,
    const char* message) {
    if (value >
        static_cast<std::size_t>(
            std::numeric_limits<
                CsrAdjacency::Offset>::max())) {
        throw std::length_error(message);
    }
    return static_cast<
        CsrAdjacency::Offset>(value);
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

struct Vector3D {
    double x;
    double y;
    double z;
};

[[nodiscard]] inline Vector3D subtract(
    Coordinate3D left,
    Coordinate3D right) {
    return Vector3D{
        left.x_m - right.x_m,
        left.y_m - right.y_m,
        left.z_m - right.z_m};
}

[[nodiscard]] inline Vector3D cross(
    Vector3D left,
    Vector3D right) {
    return Vector3D{
        left.y * right.z - left.z * right.y,
        left.z * right.x - left.x * right.z,
        left.x * right.y - left.y * right.x};
}

[[nodiscard]] inline double dot(
    Vector3D left,
    Vector3D right) {
    return left.x * right.x +
           left.y * right.y +
           left.z * right.z;
}

[[nodiscard]] inline double magnitude(
    Vector3D value) {
    return std::sqrt(dot(value, value));
}

[[nodiscard]] inline Coordinate3D triangle_centroid(
    Coordinate3D a,
    Coordinate3D b,
    Coordinate3D c) {
    return Coordinate3D{
        (a.x_m + b.x_m + c.x_m) / 3.0,
        (a.y_m + b.y_m + c.y_m) / 3.0,
        (a.z_m + b.z_m + c.z_m) / 3.0};
}

[[nodiscard]] inline Coordinate3D cell_vertex_mean(
    const std::array<std::size_t, 8>& vertices,
    std::span<const Coordinate3D> coordinates) {
    Coordinate3D centroid{0.0, 0.0, 0.0};
    for (const auto vertex : vertices) {
        if (vertex >= coordinates.size()) {
            throw std::out_of_range(
                "mpmc::mesh::process_active_corner_point_grid: cell vertex index out of range while computing centroid");
        }
        centroid.x_m += coordinates[vertex].x_m;
        centroid.y_m += coordinates[vertex].y_m;
        centroid.z_m += coordinates[vertex].z_m;
    }
    centroid.x_m /= 8.0;
    centroid.y_m /= 8.0;
    centroid.z_m /= 8.0;
    return centroid;
}

struct QuadMetric {
    Coordinate3D centroid;
    double area_m2;
    UnitVector3D unit_normal;
};

[[nodiscard]] inline QuadMetric quad_metric(
    const std::array<std::size_t, 4>& vertices,
    std::span<const Coordinate3D> coordinates,
    Coordinate3D owner_centroid) {
    std::array<Coordinate3D, 4> points{};
    for (std::size_t i = 0U; i < 4U; ++i) {
        if (vertices[i] >= coordinates.size()) {
            throw std::out_of_range(
                "mpmc::mesh::process_active_corner_point_grid: face vertex index out of range while computing metric");
        }
        points[i] = coordinates[vertices[i]];
    }

    double max_edge = 0.0;
    for (std::size_t i = 0U; i < 4U; ++i) {
        const auto edge =
            subtract(
                points[(i + 1U) % 4U],
                points[i]);
        const double length =
            magnitude(edge);
        if (!std::isfinite(length)) {
            throw std::invalid_argument(
                "mpmc::mesh::process_active_corner_point_grid: non-finite quad edge length");
        }
        max_edge =
            std::max(max_edge, length);
    }
    if (max_edge <= 0.0) {
        throw std::invalid_argument(
            "mpmc::mesh::process_active_corner_point_grid: zero-size quad face");
    }

    const Vector3D doubled_area0 =
        cross(
            subtract(points[1], points[0]),
            subtract(points[2], points[0]));
    const Vector3D doubled_area1 =
        cross(
            subtract(points[2], points[0]),
            subtract(points[3], points[0]));
    const double magnitude0 =
        magnitude(doubled_area0);
    const double magnitude1 =
        magnitude(doubled_area1);
    const double area_tolerance =
        4096.0 *
        std::numeric_limits<double>::epsilon() *
        max_edge * max_edge;
    if (!std::isfinite(magnitude0) ||
        !std::isfinite(magnitude1) ||
        magnitude0 <= area_tolerance ||
        magnitude1 <= area_tolerance) {
        throw std::invalid_argument(
            "mpmc::mesh::process_active_corner_point_grid: degenerate triangle in quad face");
    }

    const double fold_tolerance =
        4096.0 *
        std::numeric_limits<double>::epsilon() *
        magnitude0 * magnitude1;
    if (dot(doubled_area0, doubled_area1) <
        -fold_tolerance) {
        throw std::invalid_argument(
            "mpmc::mesh::process_active_corner_point_grid: folded quad face has opposing triangle normals");
    }

    const double area0 =
        0.5 * magnitude0;
    const double area1 =
        0.5 * magnitude1;
    const double area =
        area0 + area1;
    const auto centroid0 =
        triangle_centroid(
            points[0], points[1], points[2]);
    const auto centroid1 =
        triangle_centroid(
            points[0], points[2], points[3]);
    const Coordinate3D centroid{
        (area0 * centroid0.x_m +
         area1 * centroid1.x_m) / area,
        (area0 * centroid0.y_m +
         area1 * centroid1.y_m) / area,
        (area0 * centroid0.z_m +
         area1 * centroid1.z_m) / area};

    Vector3D area_vector{
        doubled_area0.x + doubled_area1.x,
        doubled_area0.y + doubled_area1.y,
        doubled_area0.z + doubled_area1.z};
    const double area_vector_magnitude =
        magnitude(area_vector);
    if (!std::isfinite(area_vector_magnitude) ||
        area_vector_magnitude <= area_tolerance) {
        throw std::invalid_argument(
            "mpmc::mesh::process_active_corner_point_grid: quad face has near-zero resultant area vector");
    }

    Vector3D normal{
        area_vector.x / area_vector_magnitude,
        area_vector.y / area_vector_magnitude,
        area_vector.z / area_vector_magnitude};
    const Vector3D owner_to_face{
        centroid.x_m - owner_centroid.x_m,
        centroid.y_m - owner_centroid.y_m,
        centroid.z_m - owner_centroid.z_m};
    const double orientation =
        dot(normal, owner_to_face);
    const double orientation_tolerance =
        4096.0 *
        std::numeric_limits<double>::epsilon() *
        max_edge;
    if (std::abs(orientation) <=
        orientation_tolerance) {
        throw std::invalid_argument(
            "mpmc::mesh::process_active_corner_point_grid: cannot orient quad normal away from owner cell");
    }
    if (orientation < 0.0) {
        normal.x = -normal.x;
        normal.y = -normal.y;
        normal.z = -normal.z;
    }

    return QuadMetric{
        centroid,
        area,
        UnitVector3D{
            normal.x,
            normal.y,
            normal.z}};
}

[[nodiscard]] inline std::array<std::size_t, 4>
face_corner_slots(std::size_t face_slot) {
    // Cell corner ordering is
    // {LLL,HLL,LHL,HHL,LLH,HLH,LHH,HHH}.
    switch (face_slot) {
    case 0U: return {0U, 4U, 6U, 2U}; // I-
    case 1U: return {1U, 3U, 7U, 5U}; // I+
    case 2U: return {0U, 1U, 5U, 4U}; // J-
    case 3U: return {2U, 6U, 7U, 3U}; // J+
    case 4U: return {0U, 2U, 3U, 1U}; // K-
    case 5U: return {4U, 5U, 7U, 6U}; // K+
    default:
        throw std::out_of_range(
            "mpmc::mesh::process_active_corner_point_grid: invalid face slot");
    }
}

[[nodiscard]] inline std::array<std::size_t, 4>
opposite_face_slots(std::size_t axis,
                    bool positive_side) {
    if (axis == 0U) {
        return face_corner_slots(
            positive_side ? 1U : 0U);
    }
    if (axis == 1U) {
        return face_corner_slots(
            positive_side ? 3U : 2U);
    }
    if (axis == 2U) {
        return face_corner_slots(
            positive_side ? 5U : 4U);
    }
    throw std::out_of_range(
        "mpmc::mesh::process_active_corner_point_grid: invalid axis");
}

[[nodiscard]] inline std::array<CoordinateKey, 4>
sorted_face_coordinates(
    const std::array<Coordinate3D, 8>& corners,
    const std::array<std::size_t, 4>& slots) {
    std::array<CoordinateKey, 4> keys{
        coordinate_key(corners[slots[0]]),
        coordinate_key(corners[slots[1]]),
        coordinate_key(corners[slots[2]]),
        coordinate_key(corners[slots[3]])};
    std::sort(keys.begin(), keys.end());
    return keys;
}

[[nodiscard]] inline std::array<std::size_t, 3>
logical_ijk(std::size_t logical_cell,
            const std::array<std::size_t, 3>& dimensions) {
    const std::size_t nx = dimensions[0];
    const std::size_t ny = dimensions[1];
    const std::size_t plane = nx * ny;
    const std::size_t k = logical_cell / plane;
    const std::size_t remainder =
        logical_cell % plane;
    const std::size_t j =
        remainder / nx;
    const std::size_t i =
        remainder % nx;
    return {i, j, k};
}

[[nodiscard]] inline bool are_axis_neighbors(
    std::size_t first,
    std::size_t second,
    const std::array<std::size_t, 3>& dimensions) {
    const auto a =
        logical_ijk(first, dimensions);
    const auto b =
        logical_ijk(second, dimensions);
    const std::size_t di =
        a[0] > b[0]
            ? a[0] - b[0]
            : b[0] - a[0];
    const std::size_t dj =
        a[1] > b[1]
            ? a[1] - b[1]
            : b[1] - a[1];
    const std::size_t dk =
        a[2] > b[2]
            ? a[2] - b[2]
            : b[2] - a[2];
    return di + dj + dk == 1U;
}

inline void require_no_active_fault_split(
    const GrdeclImportResult& raw) {
    const auto dims = raw.dimensions;
    const std::size_t nx = dims[0];
    const std::size_t ny = dims[1];
    const std::size_t nz = dims[2];

    const auto logical =
        [nx, ny](std::size_t i,
                 std::size_t j,
                 std::size_t k) {
            return i + nx * (j + ny * k);
        };

    for (std::size_t k = 0U; k < nz; ++k) {
        for (std::size_t j = 0U; j < ny; ++j) {
            for (std::size_t i = 0U; i < nx; ++i) {
                const std::size_t first =
                    logical(i, j, k);
                if (!raw.is_active(
                        checked_local(
                            first,
                            "mpmc::mesh::process_active_corner_point_grid: logical cell index overflow"))) {
                    continue;
                }
                const auto first_corners =
                    raw.geometry.cell_corners_m(
                        checked_local(
                            first,
                            "mpmc::mesh::process_active_corner_point_grid: logical cell index overflow"));

                const auto check_neighbor =
                    [&](std::size_t second,
                        std::size_t axis) {
                        if (!raw.is_active(
                                checked_local(
                                    second,
                                    "mpmc::mesh::process_active_corner_point_grid: neighbor cell index overflow"))) {
                            return;
                        }
                        const auto second_corners =
                            raw.geometry.cell_corners_m(
                                checked_local(
                                    second,
                                    "mpmc::mesh::process_active_corner_point_grid: neighbor cell index overflow"));
                        const auto first_face =
                            sorted_face_coordinates(
                                first_corners,
                                opposite_face_slots(
                                    axis, true));
                        const auto second_face =
                            sorted_face_coordinates(
                                second_corners,
                                opposite_face_slots(
                                    axis, false));
                        if (first_face != second_face) {
                            throw std::invalid_argument(
                                "mpmc::mesh::process_active_corner_point_grid: active logical neighbors have split/noncoincident face geometry; faults/NNC are unsupported in this baseline");
                        }
                    };

                if (i + 1U < nx) {
                    check_neighbor(
                        logical(i + 1U, j, k),
                        0U);
                }
                if (j + 1U < ny) {
                    check_neighbor(
                        logical(i, j + 1U, k),
                        1U);
                }
                if (k + 1U < nz) {
                    check_neighbor(
                        logical(i, j, k + 1U),
                        2U);
                }
            }
        }
    }
}

} // namespace active_corner_point_detail

[[nodiscard]] inline ActiveCornerPointGrid
process_active_corner_point_grid(
    const GrdeclImportResult& raw) {
    using namespace active_corner_point_detail;

    if (raw.topology.entity_count(
            EntityKind::cell) !=
            raw.cell_count() ||
        raw.geometry.cell_count() !=
            raw.cell_count()) {
        throw std::invalid_argument(
            "mpmc::mesh::process_active_corner_point_grid: raw GRDECL cell counts are inconsistent");
    }

    const std::size_t logical_cell_count =
        checked_multiply(
            checked_multiply(
                raw.dimensions[0],
                raw.dimensions[1],
                "mpmc::mesh::process_active_corner_point_grid: logical XY cell count overflow"),
            raw.dimensions[2],
            "mpmc::mesh::process_active_corner_point_grid: logical XYZ cell count overflow");
    if (logical_cell_count != raw.cell_count()) {
        throw std::invalid_argument(
            "mpmc::mesh::process_active_corner_point_grid: dimensions do not match raw logical cell count");
    }
    if (!raw.topology.has_relation(
            EntityKind::cell,
            EntityKind::vertex)) {
        throw std::invalid_argument(
            "mpmc::mesh::process_active_corner_point_grid: raw cell->vertex relation is missing");
    }

    const auto& raw_cell_vertices =
        raw.topology.relation(
            EntityKind::cell,
            EntityKind::vertex);
    for (std::size_t cell = 0U;
         cell < raw.cell_count();
         ++cell) {
        if (raw_cell_vertices
                .adjacent(
                    checked_local(
                        cell,
                        "mpmc::mesh::process_active_corner_point_grid: raw cell local index overflow"))
                .size() != 8U) {
            throw std::invalid_argument(
                "mpmc::mesh::process_active_corner_point_grid: every raw cell must have eight corner vertices");
        }
    }

    require_no_active_fault_split(raw);

    std::vector<std::size_t>
        active_logical_cells;
    active_logical_cells.reserve(
        raw.active_cell_count());
    for (std::size_t logical_cell = 0U;
         logical_cell < raw.cell_count();
         ++logical_cell) {
        if (raw.is_active(
                checked_local(
                    logical_cell,
                    "mpmc::mesh::process_active_corner_point_grid: logical cell index overflow"))) {
            active_logical_cells.push_back(
                logical_cell);
        }
    }
    if (active_logical_cells.empty()) {
        throw std::invalid_argument(
            "mpmc::mesh::process_active_corner_point_grid: no active cells to process");
    }

    std::map<CoordinateKey, Coordinate3D>
        unique_coordinates;
    for (const std::size_t logical_cell :
         active_logical_cells) {
        const auto corners =
            raw.geometry.cell_corners_m(
                checked_local(
                    logical_cell,
                    "mpmc::mesh::process_active_corner_point_grid: logical cell index overflow"));
        for (const auto corner : corners) {
            unique_coordinates.emplace(
                coordinate_key(corner),
                corner);
        }
    }

    if (unique_coordinates.size() >
        static_cast<std::size_t>(
            std::numeric_limits<
                LocalIndex::value_type>::max()) +
            std::size_t{1U}) {
        throw std::length_error(
            "mpmc::mesh::process_active_corner_point_grid: merged vertex count exceeds LocalIndex capacity");
    }

    std::vector<Coordinate3D>
        processed_coordinates;
    processed_coordinates.reserve(
        unique_coordinates.size());
    std::map<CoordinateKey, std::size_t>
        vertex_local_by_coordinate;
    std::size_t next_vertex = 0U;
    for (const auto& [key, coordinate] :
         unique_coordinates) {
        vertex_local_by_coordinate.emplace(
            key, next_vertex++);
        processed_coordinates.push_back(
            coordinate);
    }

    const std::size_t processed_cell_count =
        active_logical_cells.size();
    std::vector<std::array<std::size_t, 8>>
        processed_cell_vertices;
    processed_cell_vertices.reserve(
        processed_cell_count);
    std::vector<double> processed_volumes;
    processed_volumes.reserve(
        processed_cell_count);
    std::vector<GlobalEntityId>
        source_logical_ids;
    source_logical_ids.reserve(
        processed_cell_count);

    std::map<FaceKey, FaceBuild>
        faces_by_key;
    std::vector<std::array<FaceKey, 6>>
        cell_face_keys;
    cell_face_keys.reserve(
        processed_cell_count);

    for (std::size_t processed_cell = 0U;
         processed_cell < processed_cell_count;
         ++processed_cell) {
        const std::size_t logical_cell =
            active_logical_cells[
                processed_cell];
        const auto raw_local =
            checked_local(
                logical_cell,
                "mpmc::mesh::process_active_corner_point_grid: logical cell index overflow");
        const auto corners =
            raw.geometry.cell_corners_m(
                raw_local);

        std::array<std::size_t, 8>
            merged_vertices{};
        for (std::size_t corner = 0U;
             corner < 8U;
             ++corner) {
            const auto found =
                vertex_local_by_coordinate.find(
                    coordinate_key(
                        corners[corner]));
            if (found ==
                vertex_local_by_coordinate.end()) {
                throw std::logic_error(
                    "mpmc::mesh::process_active_corner_point_grid: merged vertex lookup failed");
            }
            merged_vertices[corner] =
                found->second;
        }
        std::set<std::size_t> unique_cell_vertices(
            merged_vertices.begin(),
            merged_vertices.end());
        if (unique_cell_vertices.size() != 8U) {
            throw std::invalid_argument(
                "mpmc::mesh::process_active_corner_point_grid: active cell collapses to fewer than eight unique merged vertices");
        }

        // Raw GRDECL corners use
        // {LLL,HLL,LHL,HHL,LLH,HLH,LHH,HHH}. The processed computational
        // topology uses the common linear-hexa ordering
        // {LLL,HLL,HHL,LHL,LLH,HLH,HHH,LHH}, matching VTK/Gmsh connectivity.
        // Face construction below still consumes the raw corner slots so the
        // existing I/J/K face semantics are unchanged.
        processed_cell_vertices.push_back(
            std::array<std::size_t, 8>{
                merged_vertices[0],
                merged_vertices[1],
                merged_vertices[3],
                merged_vertices[2],
                merged_vertices[4],
                merged_vertices[5],
                merged_vertices[7],
                merged_vertices[6]});
        processed_volumes.push_back(
            raw.geometry.cell_volume_m3(
                raw_local));
        const auto logical_id =
            raw.topology.global_id(
                EntityKind::cell,
                raw_local);
        source_logical_ids.push_back(
            logical_id);

        std::array<FaceKey, 6>
            keys{};
        for (std::size_t face_slot = 0U;
             face_slot < 6U;
             ++face_slot) {
            const auto slots =
                face_corner_slots(face_slot);
            std::array<std::size_t, 4>
                ordered{
                    merged_vertices[slots[0]],
                    merged_vertices[slots[1]],
                    merged_vertices[slots[2]],
                    merged_vertices[slots[3]]};
            std::array<std::size_t, 4>
                sorted = ordered;
            std::sort(
                sorted.begin(),
                sorted.end());
            if (std::adjacent_find(
                    sorted.begin(),
                    sorted.end()) !=
                sorted.end()) {
                throw std::invalid_argument(
                    "mpmc::mesh::process_active_corner_point_grid: active face has duplicate merged vertices");
            }
            const FaceKey key{sorted};
            keys[face_slot] = key;

            auto [found, inserted] =
                faces_by_key.emplace(
                    key,
                    FaceBuild{
                        key,
                        ordered,
                        {},
                        {}});
            found->second.processed_cells
                .push_back(processed_cell);
            found->second.source_logical_cells
                .push_back(logical_cell);
            if (found->second
                    .processed_cells.size() >
                2U) {
                throw std::invalid_argument(
                    "mpmc::mesh::process_active_corner_point_grid: merged face is incident to more than two active cells");
            }
            (void)inserted;
        }
        cell_face_keys.push_back(keys);
    }

    std::vector<FaceBuild> faces;
    faces.reserve(
        faces_by_key.size());
    for (auto& [key, face] :
         faces_by_key) {
        if (face.processed_cells.size() == 2U &&
            !are_axis_neighbors(
                face.source_logical_cells[0],
                face.source_logical_cells[1],
                raw.dimensions)) {
            throw std::invalid_argument(
                "mpmc::mesh::process_active_corner_point_grid: non-neighbor active cells geometrically share a full face");
        }
        faces.push_back(
            std::move(face));
        (void)key;
    }

    std::map<FaceKey, std::size_t>
        face_local_by_key;
    for (std::size_t face = 0U;
         face < faces.size();
         ++face) {
        face_local_by_key.emplace(
            faces[face].key,
            face);
    }

    Topology::EntityIds ids;
    ids.vertices.reserve(
        processed_coordinates.size());
    for (std::size_t vertex = 0U;
         vertex < processed_coordinates.size();
         ++vertex) {
        ids.vertices.emplace_back(
            static_cast<std::uint64_t>(
                vertex) +
            1U);
    }
    ids.faces.reserve(faces.size());
    for (std::size_t face = 0U;
         face < faces.size();
         ++face) {
        ids.faces.emplace_back(
            static_cast<std::uint64_t>(
                face) +
            1U);
    }
    ids.cells.reserve(
        processed_cell_count);
    for (const auto id :
         source_logical_ids) {
        ids.cells.push_back(id);
    }

    std::vector<CsrAdjacency::Offset>
        cell_vertex_offsets{0U};
    std::vector<LocalIndex>
        cell_vertices;
    std::vector<CsrAdjacency::Offset>
        cell_face_offsets{0U};
    std::vector<LocalIndex>
        cell_faces;
    cell_vertex_offsets.reserve(
        processed_cell_count + 1U);
    cell_face_offsets.reserve(
        processed_cell_count + 1U);
    cell_vertices.reserve(
        checked_multiply(
            processed_cell_count,
            8U,
            "mpmc::mesh::process_active_corner_point_grid: cell->vertex storage overflow"));
    cell_faces.reserve(
        checked_multiply(
            processed_cell_count,
            6U,
            "mpmc::mesh::process_active_corner_point_grid: cell->face storage overflow"));

    for (std::size_t cell = 0U;
         cell < processed_cell_count;
         ++cell) {
        for (const auto vertex :
             processed_cell_vertices[cell]) {
            cell_vertices.push_back(
                checked_local(
                    vertex,
                    "mpmc::mesh::process_active_corner_point_grid: merged vertex local index overflow"));
        }
        for (const auto& key :
             cell_face_keys[cell]) {
            const auto found =
                face_local_by_key.find(key);
            if (found ==
                face_local_by_key.end()) {
                throw std::logic_error(
                    "mpmc::mesh::process_active_corner_point_grid: merged face lookup failed");
            }
            cell_faces.push_back(
                checked_local(
                    found->second,
                    "mpmc::mesh::process_active_corner_point_grid: merged face local index overflow"));
        }
        cell_vertex_offsets.push_back(
            checked_offset(
                cell_vertices.size(),
                "mpmc::mesh::process_active_corner_point_grid: cell->vertex CSR overflow"));
        cell_face_offsets.push_back(
            checked_offset(
                cell_faces.size(),
                "mpmc::mesh::process_active_corner_point_grid: cell->face CSR overflow"));
    }

    std::vector<CsrAdjacency::Offset>
        face_vertex_offsets{0U};
    std::vector<LocalIndex>
        face_vertices;
    std::vector<CsrAdjacency::Offset>
        face_cell_offsets{0U};
    std::vector<LocalIndex>
        face_cells;
    face_vertex_offsets.reserve(
        faces.size() + 1U);
    face_cell_offsets.reserve(
        faces.size() + 1U);
    face_vertices.reserve(
        checked_multiply(
            faces.size(),
            4U,
            "mpmc::mesh::process_active_corner_point_grid: face->vertex storage overflow"));

    for (const auto& face : faces) {
        for (const auto vertex :
             face.ordered_vertices) {
            face_vertices.push_back(
                checked_local(
                    vertex,
                    "mpmc::mesh::process_active_corner_point_grid: face vertex local index overflow"));
        }
        for (const auto cell :
             face.processed_cells) {
            face_cells.push_back(
                checked_local(
                    cell,
                    "mpmc::mesh::process_active_corner_point_grid: face cell local index overflow"));
        }
        face_vertex_offsets.push_back(
            checked_offset(
                face_vertices.size(),
                "mpmc::mesh::process_active_corner_point_grid: face->vertex CSR overflow"));
        face_cell_offsets.push_back(
            checked_offset(
                face_cells.size(),
                "mpmc::mesh::process_active_corner_point_grid: face->cell CSR overflow"));
    }

    std::vector<CsrAdjacency> relations;
    relations.reserve(4U);
    relations.emplace_back(
        EntityKind::cell,
        EntityKind::vertex,
        processed_coordinates.size(),
        std::move(cell_vertex_offsets),
        std::move(cell_vertices));
    relations.emplace_back(
        EntityKind::cell,
        EntityKind::face,
        faces.size(),
        std::move(cell_face_offsets),
        std::move(cell_faces));
    relations.emplace_back(
        EntityKind::face,
        EntityKind::vertex,
        processed_coordinates.size(),
        std::move(face_vertex_offsets),
        std::move(face_vertices));
    relations.emplace_back(
        EntityKind::face,
        EntityKind::cell,
        processed_cell_count,
        std::move(face_cell_offsets),
        std::move(face_cells));

    Topology topology{
        std::move(ids),
        std::move(relations)};

    std::vector<Coordinate3D>
        processed_cell_centroids;
    processed_cell_centroids.reserve(
        processed_cell_count);
    for (const auto& vertices :
         processed_cell_vertices) {
        processed_cell_centroids.push_back(
            cell_vertex_mean(
                vertices,
                processed_coordinates));
    }

    std::vector<Coordinate3D>
        face_centroids;
    std::vector<double>
        face_areas;
    std::vector<LocalIndex>
        face_owners;
    std::vector<UnitVector3D>
        face_owner_normals;
    face_centroids.reserve(faces.size());
    face_areas.reserve(faces.size());
    face_owners.reserve(faces.size());
    face_owner_normals.reserve(
        faces.size());

    for (const auto& face : faces) {
        if (face.processed_cells.empty()) {
            throw std::logic_error(
                "mpmc::mesh::process_active_corner_point_grid: merged face has no owner cell");
        }
        const std::size_t owner =
            face.processed_cells.front();
        if (owner >=
            processed_cell_centroids.size()) {
            throw std::logic_error(
                "mpmc::mesh::process_active_corner_point_grid: merged face owner index is invalid");
        }
        const auto metric =
            quad_metric(
                face.ordered_vertices,
                processed_coordinates,
                processed_cell_centroids[
                    owner]);
        face_centroids.push_back(
            metric.centroid);
        face_areas.push_back(
            metric.area_m2);
        face_owners.push_back(
            checked_local(
                owner,
                "mpmc::mesh::process_active_corner_point_grid: face owner local index overflow"));
        face_owner_normals.push_back(
            metric.unit_normal);
    }

    FaceGeometry3D face_geometry{
        processed_cell_count,
        std::move(face_centroids),
        std::move(face_areas),
        std::move(face_owners),
        std::move(face_owner_normals)};

    std::map<std::uint64_t, std::size_t>
        raw_cell_local_by_id;
    const auto raw_cell_ids =
        raw.topology.global_ids(
            EntityKind::cell);
    for (std::size_t local = 0U;
         local < raw_cell_ids.size();
         ++local) {
        if (!raw_cell_local_by_id
                 .emplace(
                     raw_cell_ids[local].value(),
                     local)
                 .second) {
            throw std::logic_error(
                "mpmc::mesh::process_active_corner_point_grid: duplicate raw cell GlobalEntityId");
        }
    }
    std::vector<DenseFieldSnapshot>
        projected_fields;
    projected_fields.reserve(
        raw.cell_fields.size());
    std::set<std::string>
        projected_field_ids;
    for (const auto& field :
         raw.cell_fields) {
        if (field.location() !=
                EntityKind::cell ||
            field.entity_count() !=
                raw.cell_count()) {
            throw std::invalid_argument(
                "mpmc::mesh::process_active_corner_point_grid: raw GRDECL field is not aligned to logical cells");
        }
        if (!projected_field_ids.insert(
                field.metadata().id)
                 .second) {
            throw std::invalid_argument(
                "mpmc::mesh::process_active_corner_point_grid: duplicate raw GRDECL cell field ID");
        }

        const std::size_t components =
            field.component_count();
        std::vector<double> values;
        values.reserve(
            checked_multiply(
                processed_cell_count,
                components,
                "mpmc::mesh::process_active_corner_point_grid: projected field storage overflow"));

        for (const auto logical_id :
             source_logical_ids) {
            const auto local_found =
                raw_cell_local_by_id.find(
                    logical_id.value());
            if (local_found ==
                raw_cell_local_by_id.end()) {
                throw std::logic_error(
                    "mpmc::mesh::process_active_corner_point_grid: source logical cell ID is absent from raw topology");
            }
            const auto raw_values =
                field.entity_values(
                    checked_local(
                        local_found->second,
                        "mpmc::mesh::process_active_corner_point_grid: raw field cell local index overflow"));
            values.insert(
                values.end(),
                raw_values.begin(),
                raw_values.end());
        }

        projected_fields.push_back(
            DenseFieldSnapshot::create(
                topology,
                EntityKind::cell,
                components,
                std::move(values),
                field.metadata()));
    }


                cell,
                "mpmc::mesh::process_active_corner_point_grid: processed cell local index overflow");
        if (topology.global_id(
                EntityKind::cell,
                local) !=
            source_logical_ids[cell]) {
            throw std::logic_error(
                "mpmc::mesh::process_active_corner_point_grid: processed cell identity mapping drift");
        }
    }

    return ActiveCornerPointGrid{
        std::move(topology),
        std::move(processed_coordinates),
        std::move(processed_volumes),
        std::move(face_geometry),
        std::move(projected_fields),
        std::move(source_logical_ids)};
}

} // namespace mpmc::mesh

#endif // MPMC_MESH_ACTIVE_CORNER_POINT_HPP
