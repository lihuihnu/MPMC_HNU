#ifndef MPMC_MESH_LINEAR_CELL_MESH_3D_HPP
#define MPMC_MESH_LINEAR_CELL_MESH_3D_HPP

#include <mpmc/mesh/csr_adjacency.hpp>
#include <mpmc/mesh/face_boundary.hpp>
#include <mpmc/mesh/face_geometry_3d.hpp>
#include <mpmc/mesh/topology.hpp>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <map>
#include <set>
#include <span>
#include <stdexcept>
#include <utility>
#include <vector>

namespace mpmc::mesh {

enum class LinearCellType3D : std::uint8_t {
    tetrahedron = 0,
    hexahedron = 1,
    wedge = 2,
    pyramid = 3,
};

struct LinearCell3D {
    GlobalEntityId global_id;
    LinearCellType3D type;
    std::vector<LocalIndex> vertices;
};

struct LinearFaceAnnotation3D {
    std::vector<LocalIndex> vertices;
    GlobalEntityId global_id;
    PhysicalTag physical_tag;
};

/// Common immutable output for linear 3D mesh importers.
///
/// Topology owns stable vertex/face/cell identities and compact CSR relations.
/// Coordinates are SI metres, cell volumes are m^3, and FaceGeometry3D stores
/// face centroid/area plus an owner-relative outward unit normal. Edges are not
/// materialized by this baseline.
struct LinearMesh3D {
    Topology topology;
    std::vector<Coordinate3D> vertex_coordinates_m;
    std::vector<double> cell_volumes_m3;
    FaceGeometry3D face_geometry;
    FaceBoundarySnapshot face_boundary;
};

namespace linear_cell_mesh_3d_detail {

struct Vector3D {
    double x;
    double y;
    double z;
};

struct FaceKey {
    std::uint8_t size;
    std::array<LocalIndex::value_type, 4> vertices;

    [[nodiscard]] friend bool operator<(
        const FaceKey& left,
        const FaceKey& right) noexcept {
        if (left.size != right.size) {
            return left.size < right.size;
        }
        return left.vertices < right.vertices;
    }
};

struct FaceBuild {
    FaceKey key;
    std::vector<LocalIndex> ordered_vertices;
    std::vector<LocalIndex> adjacent_cells;
    GlobalEntityId global_id{0U};
    PhysicalTag physical_tag{0U};
    bool annotated{false};
};

[[nodiscard]] inline LocalIndex checked_local(
    std::size_t value,
    const char* message) {
    if (value >
        static_cast<std::size_t>(
            std::numeric_limits<LocalIndex::value_type>::max())) {
        throw std::length_error(message);
    }
    return LocalIndex{
        static_cast<LocalIndex::value_type>(value)};
}

[[nodiscard]] inline CsrAdjacency::Offset checked_offset(
    std::size_t value,
    const char* message) {
    if (value >
        static_cast<std::size_t>(
            std::numeric_limits<CsrAdjacency::Offset>::max())) {
        throw std::length_error(message);
    }
    return static_cast<CsrAdjacency::Offset>(value);
}

[[nodiscard]] inline Vector3D subtract(
    Coordinate3D left,
    Coordinate3D right) noexcept {
    return Vector3D{
        left.x_m - right.x_m,
        left.y_m - right.y_m,
        left.z_m - right.z_m};
}

[[nodiscard]] inline Vector3D cross(
    Vector3D left,
    Vector3D right) noexcept {
    return Vector3D{
        left.y * right.z - left.z * right.y,
        left.z * right.x - left.x * right.z,
        left.x * right.y - left.y * right.x};
}

[[nodiscard]] inline double dot(
    Vector3D left,
    Vector3D right) noexcept {
    return left.x * right.x +
           left.y * right.y +
           left.z * right.z;
}

[[nodiscard]] inline double magnitude(
    Vector3D value) noexcept {
    return std::sqrt(dot(value, value));
}

[[nodiscard]] inline double distance(
    Coordinate3D left,
    Coordinate3D right) noexcept {
    return magnitude(subtract(left, right));
}

[[nodiscard]] inline double signed_tetra_volume(
    Coordinate3D a,
    Coordinate3D b,
    Coordinate3D c,
    Coordinate3D d) noexcept {
    return dot(
               subtract(b, a),
               cross(
                   subtract(c, a),
                   subtract(d, a))) /
           6.0;
}

[[nodiscard]] inline Coordinate3D average(
    std::span<const LocalIndex> vertices,
    std::span<const Coordinate3D> coordinates) {
    Coordinate3D result{0.0, 0.0, 0.0};
    for (const auto vertex : vertices) {
        const auto local =
            static_cast<std::size_t>(vertex.value());
        if (local >= coordinates.size()) {
            throw std::out_of_range(
                "mpmc::mesh::make_linear_mesh_3d: cell/face vertex index out of range");
        }
        const auto point = coordinates[local];
        result.x_m += point.x_m;
        result.y_m += point.y_m;
        result.z_m += point.z_m;
    }
    const double denominator =
        static_cast<double>(vertices.size());
    result.x_m /= denominator;
    result.y_m /= denominator;
    result.z_m /= denominator;
    return result;
}

[[nodiscard]] inline double characteristic_length(
    std::span<const LocalIndex> vertices,
    std::span<const Coordinate3D> coordinates) {
    double scale = 0.0;
    for (std::size_t i = 0U; i < vertices.size(); ++i) {
        for (std::size_t j = i + 1U;
             j < vertices.size();
             ++j) {
            const auto first =
                coordinates[static_cast<std::size_t>(
                    vertices[i].value())];
            const auto second =
                coordinates[static_cast<std::size_t>(
                    vertices[j].value())];
            scale = std::max(
                scale,
                distance(first, second));
        }
    }
    if (!std::isfinite(scale) || scale <= 0.0) {
        throw std::invalid_argument(
            "mpmc::mesh::make_linear_mesh_3d: degenerate cell/face characteristic length");
    }
    return scale;
}

[[nodiscard]] inline FaceKey face_key(
    std::span<const LocalIndex> vertices) {
    if (vertices.size() != 3U &&
        vertices.size() != 4U) {
        throw std::invalid_argument(
            "mpmc::mesh::make_linear_mesh_3d: face must contain three or four vertices");
    }

    FaceKey key{};
    key.size =
        static_cast<std::uint8_t>(vertices.size());
    key.vertices.fill(
        std::numeric_limits<LocalIndex::value_type>::max());
    for (std::size_t i = 0U;
         i < vertices.size();
         ++i) {
        key.vertices[i] = vertices[i].value();
    }
    // The face width is only three or four. Keep the ordering logic
    // explicitly bounded instead of invoking a general introsort over a
    // runtime-sized prefix of the fixed array; GCC 13 otherwise emits a
    // false-positive -Warray-bounds under Release inlining.
    const std::size_t count = vertices.size();
    for (std::size_t i = 1U; i < count; ++i) {
        const auto value = key.vertices[i];
        std::size_t j = i;
        while (j > 0U &&
               value < key.vertices[j - 1U]) {
            key.vertices[j] =
                key.vertices[j - 1U];
            --j;
        }
        key.vertices[j] = value;
    }
    for (std::size_t i = 1U; i < count; ++i) {
        if (key.vertices[i - 1U] ==
            key.vertices[i]) {
            throw std::invalid_argument(
                "mpmc::mesh::make_linear_mesh_3d: face contains repeated vertices");
        }
    }
    return key;
}

[[nodiscard]] inline std::size_t
cell_vertex_count(LinearCellType3D type) {
    switch (type) {
    case LinearCellType3D::tetrahedron:
        return 4U;
    case LinearCellType3D::hexahedron:
        return 8U;
    case LinearCellType3D::wedge:
        return 6U;
    case LinearCellType3D::pyramid:
        return 5U;
    }
    throw std::invalid_argument(
        "mpmc::mesh::make_linear_mesh_3d: invalid linear 3D cell type");
}

[[nodiscard]] inline std::size_t
cell_face_count(LinearCellType3D type) {
    switch (type) {
    case LinearCellType3D::tetrahedron:
        return 4U;
    case LinearCellType3D::hexahedron:
        return 6U;
    case LinearCellType3D::wedge:
    case LinearCellType3D::pyramid:
        return 5U;
    }
    throw std::invalid_argument(
        "mpmc::mesh::make_linear_mesh_3d: invalid linear 3D cell type");
}

template <typename Callback>
inline void for_each_cell_face(
    const LinearCell3D& cell,
    Callback&& callback) {
    if (cell.type ==
        LinearCellType3D::tetrahedron) {
        if (cell.vertices.size() != 4U) {
            throw std::invalid_argument(
                "mpmc::mesh::make_linear_mesh_3d: tetrahedron requires four vertices");
        }
        constexpr std::array<
            std::array<std::size_t, 3>, 4>
            faces{{
                {{0U, 2U, 1U}},
                {{0U, 1U, 3U}},
                {{1U, 2U, 3U}},
                {{2U, 0U, 3U}},
            }};
        for (const auto& slots : faces) {
            const std::array<LocalIndex, 3>
                vertices{
                    cell.vertices[slots[0]],
                    cell.vertices[slots[1]],
                    cell.vertices[slots[2]]};
            callback(
                std::span<const LocalIndex>{
                    vertices.data(),
                    vertices.size()});
        }
        return;
    }

    if (cell.type ==
        LinearCellType3D::wedge) {
        if (cell.vertices.size() != 6U) {
            throw std::invalid_argument(
                "mpmc::mesh::make_linear_mesh_3d: wedge requires six vertices");
        }
        constexpr std::array<
            std::array<std::size_t, 3>, 2>
            triangle_faces{{
                {{0U, 2U, 1U}},
                {{3U, 4U, 5U}},
            }};
        for (const auto& slots :
             triangle_faces) {
            const std::array<LocalIndex, 3>
                vertices{
                    cell.vertices[slots[0]],
                    cell.vertices[slots[1]],
                    cell.vertices[slots[2]]};
            callback(
                std::span<const LocalIndex>{
                    vertices.data(),
                    vertices.size()});
        }
        constexpr std::array<
            std::array<std::size_t, 4>, 3>
            quad_faces{{
                {{0U, 1U, 4U, 3U}},
                {{1U, 2U, 5U, 4U}},
                {{2U, 0U, 3U, 5U}},
            }};
        for (const auto& slots :
             quad_faces) {
            const std::array<LocalIndex, 4>
                vertices{
                    cell.vertices[slots[0]],
                    cell.vertices[slots[1]],
                    cell.vertices[slots[2]],
                    cell.vertices[slots[3]]};
            callback(
                std::span<const LocalIndex>{
                    vertices.data(),
                    vertices.size()});
        }
        return;
    }

    if (cell.type ==
        LinearCellType3D::pyramid) {
        if (cell.vertices.size() != 5U) {
            throw std::invalid_argument(
                "mpmc::mesh::make_linear_mesh_3d: pyramid requires five vertices");
        }
        const std::array<LocalIndex, 4>
            base{
                cell.vertices[0],
                cell.vertices[3],
                cell.vertices[2],
                cell.vertices[1]};
        callback(
            std::span<const LocalIndex>{
                base.data(),
                base.size()});
        constexpr std::array<
            std::array<std::size_t, 3>, 4>
            side_faces{{
                {{0U, 1U, 4U}},
                {{1U, 2U, 4U}},
                {{2U, 3U, 4U}},
                {{3U, 0U, 4U}},
            }};
        for (const auto& slots :
             side_faces) {
            const std::array<LocalIndex, 3>
                vertices{
                    cell.vertices[slots[0]],
                    cell.vertices[slots[1]],
                    cell.vertices[slots[2]]};
            callback(
                std::span<const LocalIndex>{
                    vertices.data(),
                    vertices.size()});
        }
        return;
    }

    if (cell.type !=
            LinearCellType3D::hexahedron ||
        cell.vertices.size() != 8U) {
        throw std::invalid_argument(
            "mpmc::mesh::make_linear_mesh_3d: hexahedron requires eight vertices");
    }
    constexpr std::array<
        std::array<std::size_t, 4>, 6>
        faces{{
            {{0U, 3U, 2U, 1U}},
            {{4U, 5U, 6U, 7U}},
            {{0U, 1U, 5U, 4U}},
            {{1U, 2U, 6U, 5U}},
            {{2U, 3U, 7U, 6U}},
            {{3U, 0U, 4U, 7U}},
        }};
    for (const auto& slots : faces) {
        const std::array<LocalIndex, 4>
            vertices{
                cell.vertices[slots[0]],
                cell.vertices[slots[1]],
                cell.vertices[slots[2]],
                cell.vertices[slots[3]]};
        callback(
            std::span<const LocalIndex>{
                vertices.data(),
                vertices.size()});
    }
}

[[nodiscard]] inline double cell_volume(
    const LinearCell3D& cell,
    std::span<const Coordinate3D> coordinates) {
    const double length =
        characteristic_length(
            cell.vertices, coordinates);
    const double tolerance =
        4096.0 *
        std::numeric_limits<double>::epsilon() *
        length * length * length;

    if (cell.type ==
        LinearCellType3D::tetrahedron) {
        const auto a =
            coordinates[static_cast<std::size_t>(
                cell.vertices[0].value())];
        const auto b =
            coordinates[static_cast<std::size_t>(
                cell.vertices[1].value())];
        const auto c =
            coordinates[static_cast<std::size_t>(
                cell.vertices[2].value())];
        const auto d =
            coordinates[static_cast<std::size_t>(
                cell.vertices[3].value())];
        const double volume =
            signed_tetra_volume(a, b, c, d);
        if (!std::isfinite(volume) ||
            volume <= tolerance) {
            throw std::invalid_argument(
                "mpmc::mesh::make_linear_mesh_3d: tetrahedron is inverted or degenerate");
        }
        return volume;
    }

    if (cell.type ==
        LinearCellType3D::wedge) {
        constexpr std::array<
            std::array<std::size_t, 4>, 3>
            tetrahedra{{
                {{0U, 1U, 2U, 3U}},
                {{1U, 2U, 3U, 4U}},
                {{2U, 3U, 4U, 5U}},
            }};
        double total = 0.0;
        for (const auto& slots :
             tetrahedra) {
            const double volume =
                signed_tetra_volume(
                    coordinates[static_cast<std::size_t>(
                        cell.vertices[slots[0]].value())],
                    coordinates[static_cast<std::size_t>(
                        cell.vertices[slots[1]].value())],
                    coordinates[static_cast<std::size_t>(
                        cell.vertices[slots[2]].value())],
                    coordinates[static_cast<std::size_t>(
                        cell.vertices[slots[3]].value())]);
            if (!std::isfinite(volume) ||
                volume <= tolerance) {
                throw std::invalid_argument(
                    "mpmc::mesh::make_linear_mesh_3d: wedge is inverted, folded, or degenerate");
            }
            total += volume;
        }
        return total;
    }

    if (cell.type ==
        LinearCellType3D::pyramid) {
        constexpr std::array<
            std::array<std::size_t, 4>, 2>
            tetrahedra{{
                {{0U, 1U, 2U, 4U}},
                {{0U, 2U, 3U, 4U}},
            }};
        double total = 0.0;
        for (const auto& slots :
             tetrahedra) {
            const double volume =
                signed_tetra_volume(
                    coordinates[static_cast<std::size_t>(
                        cell.vertices[slots[0]].value())],
                    coordinates[static_cast<std::size_t>(
                        cell.vertices[slots[1]].value())],
                    coordinates[static_cast<std::size_t>(
                        cell.vertices[slots[2]].value())],
                    coordinates[static_cast<std::size_t>(
                        cell.vertices[slots[3]].value())]);
            if (!std::isfinite(volume) ||
                volume <= tolerance) {
                throw std::invalid_argument(
                    "mpmc::mesh::make_linear_mesh_3d: pyramid is inverted, folded, or degenerate");
            }
            total += volume;
        }
        return total;
    }

    if (cell.type !=
        LinearCellType3D::hexahedron) {
        throw std::invalid_argument(
            "mpmc::mesh::make_linear_mesh_3d: invalid linear 3D cell type");
    }

    constexpr std::array<
        std::array<std::size_t, 4>, 5>
        tetrahedra{{
            {{0U, 1U, 3U, 4U}},
            {{1U, 2U, 3U, 6U}},
            {{1U, 3U, 4U, 6U}},
            {{1U, 4U, 5U, 6U}},
            {{3U, 4U, 6U, 7U}},
        }};

    double total = 0.0;
    for (const auto& slots : tetrahedra) {
        const double volume =
            signed_tetra_volume(
                coordinates[static_cast<std::size_t>(
                    cell.vertices[slots[0]].value())],
                coordinates[static_cast<std::size_t>(
                    cell.vertices[slots[1]].value())],
                coordinates[static_cast<std::size_t>(
                    cell.vertices[slots[2]].value())],
                coordinates[static_cast<std::size_t>(
                    cell.vertices[slots[3]].value())]);
        if (!std::isfinite(volume) ||
            volume <= tolerance) {
            throw std::invalid_argument(
                "mpmc::mesh::make_linear_mesh_3d: hexahedron is inverted, folded, or degenerate");
        }
        total += volume;
    }
    if (!std::isfinite(total) ||
        total <= tolerance) {
        throw std::invalid_argument(
            "mpmc::mesh::make_linear_mesh_3d: invalid hexahedron volume");
    }
    return total;
}

struct FaceMetric {
    Coordinate3D centroid_m;
    double area_m2;
    UnitVector3D owner_unit_normal;
};

[[nodiscard]] inline FaceMetric face_metric(
    std::span<const LocalIndex> vertices,
    std::span<const Coordinate3D> coordinates,
    Coordinate3D owner_centroid) {
    const double length =
        characteristic_length(
            vertices, coordinates);
    const double area_tolerance =
        4096.0 *
        std::numeric_limits<double>::epsilon() *
        length * length;
    const double distance_tolerance =
        4096.0 *
        std::numeric_limits<double>::epsilon() *
        length;

    const auto p0 =
        coordinates[static_cast<std::size_t>(
            vertices[0].value())];
    const auto p1 =
        coordinates[static_cast<std::size_t>(
            vertices[1].value())];
    const auto p2 =
        coordinates[static_cast<std::size_t>(
            vertices[2].value())];
    const auto first_cross =
        cross(
            subtract(p1, p0),
            subtract(p2, p0));
    const double first_twice_area =
        magnitude(first_cross);
    if (!std::isfinite(first_twice_area) ||
        0.5 * first_twice_area <=
            area_tolerance) {
        throw std::invalid_argument(
            "mpmc::mesh::make_linear_mesh_3d: degenerate triangle on face");
    }

    Coordinate3D centroid{};
    double area = 0.5 * first_twice_area;
    Vector3D area_vector = first_cross;

    if (vertices.size() == 3U) {
        centroid = Coordinate3D{
            (p0.x_m + p1.x_m + p2.x_m) / 3.0,
            (p0.y_m + p1.y_m + p2.y_m) / 3.0,
            (p0.z_m + p1.z_m + p2.z_m) / 3.0};
    } else {
        const auto p3 =
            coordinates[static_cast<std::size_t>(
                vertices[3].value())];
        const auto second_cross =
            cross(
                subtract(p2, p0),
                subtract(p3, p0));
        const double second_twice_area =
            magnitude(second_cross);
        if (!std::isfinite(second_twice_area) ||
            0.5 * second_twice_area <=
                area_tolerance ||
            dot(first_cross, second_cross) <= 0.0) {
            throw std::invalid_argument(
                "mpmc::mesh::make_linear_mesh_3d: quad face is degenerate or folded");
        }

        const double second_area =
            0.5 * second_twice_area;
        const Coordinate3D first_centroid{
            (p0.x_m + p1.x_m + p2.x_m) / 3.0,
            (p0.y_m + p1.y_m + p2.y_m) / 3.0,
            (p0.z_m + p1.z_m + p2.z_m) / 3.0};
        const Coordinate3D second_centroid{
            (p0.x_m + p2.x_m + p3.x_m) / 3.0,
            (p0.y_m + p2.y_m + p3.y_m) / 3.0,
            (p0.z_m + p2.z_m + p3.z_m) / 3.0};
        area += second_area;
        centroid = Coordinate3D{
            (first_centroid.x_m *
                 (0.5 * first_twice_area) +
             second_centroid.x_m * second_area) /
                area,
            (first_centroid.y_m *
                 (0.5 * first_twice_area) +
             second_centroid.y_m * second_area) /
                area,
            (first_centroid.z_m *
                 (0.5 * first_twice_area) +
             second_centroid.z_m * second_area) /
                area};
        area_vector.x += second_cross.x;
        area_vector.y += second_cross.y;
        area_vector.z += second_cross.z;
    }

    const double normal_magnitude =
        magnitude(area_vector);
    if (!std::isfinite(area) ||
        area <= area_tolerance ||
        !std::isfinite(normal_magnitude) ||
        normal_magnitude <=
            2.0 * area_tolerance) {
        throw std::invalid_argument(
            "mpmc::mesh::make_linear_mesh_3d: invalid face area vector");
    }

    Vector3D normal{
        area_vector.x / normal_magnitude,
        area_vector.y / normal_magnitude,
        area_vector.z / normal_magnitude};
    const auto owner_to_face =
        subtract(centroid, owner_centroid);
    double orientation =
        dot(normal, owner_to_face);
    if (!std::isfinite(orientation) ||
        std::abs(orientation) <=
            distance_tolerance) {
        throw std::invalid_argument(
            "mpmc::mesh::make_linear_mesh_3d: face normal cannot be oriented from owner centroid");
    }
    if (orientation < 0.0) {
        normal.x = -normal.x;
        normal.y = -normal.y;
        normal.z = -normal.z;
        orientation = -orientation;
    }
    (void)orientation;

    return FaceMetric{
        centroid,
        area,
        UnitVector3D{
            normal.x,
            normal.y,
            normal.z}};
}

} // namespace linear_cell_mesh_3d_detail

[[nodiscard]] inline LinearMesh3D make_linear_mesh_3d(
    std::span<const GlobalEntityId> vertex_global_ids,
    std::span<const Coordinate3D> vertex_coordinates_m,
    std::span<const LinearCell3D> cells,
    std::span<const LinearFaceAnnotation3D> face_annotations = {},
    GlobalEntityId::value_type generated_face_id_floor = 0U) {
    using namespace linear_cell_mesh_3d_detail;

    if (vertex_global_ids.empty() ||
        vertex_coordinates_m.empty() ||
        cells.empty() ||
        vertex_global_ids.size() !=
            vertex_coordinates_m.size()) {
        throw std::invalid_argument(
            "mpmc::mesh::make_linear_mesh_3d: nonempty aligned vertices and cells are required");
    }

    for (const auto coordinate :
         vertex_coordinates_m) {
        if (!std::isfinite(coordinate.x_m) ||
            !std::isfinite(coordinate.y_m) ||
            !std::isfinite(coordinate.z_m)) {
            throw std::invalid_argument(
                "mpmc::mesh::make_linear_mesh_3d: vertex coordinates must be finite");
        }
    }

    std::vector<Coordinate3D> cell_centroids;
    std::vector<double> cell_volumes;
    cell_centroids.reserve(cells.size());
    cell_volumes.reserve(cells.size());

    std::map<FaceKey, FaceBuild> faces;
    std::vector<std::vector<FaceKey>>
        cell_face_keys;
    cell_face_keys.reserve(cells.size());

    for (std::size_t cell_index = 0U;
         cell_index < cells.size();
         ++cell_index) {
        const auto& cell = cells[cell_index];
        const std::size_t expected =
            cell_vertex_count(cell.type);
        if (cell.vertices.size() != expected) {
            throw std::invalid_argument(
                "mpmc::mesh::make_linear_mesh_3d: cell vertex count does not match linear cell type");
        }
        std::set<LocalIndex::value_type> unique;
        for (const auto vertex : cell.vertices) {
            if (static_cast<std::size_t>(
                    vertex.value()) >=
                vertex_coordinates_m.size()) {
                throw std::out_of_range(
                    "mpmc::mesh::make_linear_mesh_3d: cell vertex index out of range");
            }
            if (!unique.insert(
                    vertex.value()).second) {
                throw std::invalid_argument(
                    "mpmc::mesh::make_linear_mesh_3d: cell contains repeated vertices");
            }
        }

        cell_centroids.push_back(
            average(
                cell.vertices,
                vertex_coordinates_m));
        cell_volumes.push_back(
            cell_volume(
                cell,
                vertex_coordinates_m));

        std::vector<FaceKey> keys;
        keys.reserve(
            cell_face_count(cell.type));
        const auto cell_local =
            checked_local(
                cell_index,
                "mpmc::mesh::make_linear_mesh_3d: cell local index overflow");
        for_each_cell_face(
            cell,
            [&](std::span<const LocalIndex> vertices) {
                const auto key =
                    face_key(vertices);
                auto [found, inserted] =
                    faces.emplace(
                        key,
                        FaceBuild{
                            key,
                            {},
                            {},
                            GlobalEntityId{0U},
                            PhysicalTag{0U},
                            false});
                if (inserted) {
                    found->second.ordered_vertices.assign(
                        vertices.begin(),
                        vertices.end());
                }
                auto& adjacent =
                    found->second.adjacent_cells;
                if (adjacent.size() >= 2U) {
                    throw std::invalid_argument(
                        "mpmc::mesh::make_linear_mesh_3d: non-manifold face belongs to more than two cells");
                }
                adjacent.push_back(cell_local);
                keys.push_back(key);
            });
        cell_face_keys.push_back(
            std::move(keys));
    }

    GlobalEntityId::value_type maximum_face_id =
        generated_face_id_floor;
    std::set<GlobalEntityId::value_type>
        explicit_face_ids;
    for (const auto& annotation :
         face_annotations) {
        const auto key =
            face_key(annotation.vertices);
        const auto found = faces.find(key);
        if (found == faces.end()) {
            throw std::invalid_argument(
                "mpmc::mesh::make_linear_mesh_3d: annotated face is not a cell face");
        }
        if (found->second.annotated) {
            throw std::invalid_argument(
                "mpmc::mesh::make_linear_mesh_3d: duplicate annotation for one face");
        }
        if (!explicit_face_ids.insert(
                annotation.global_id.value())
                 .second) {
            throw std::invalid_argument(
                "mpmc::mesh::make_linear_mesh_3d: duplicate annotated face GlobalEntityId");
        }
        found->second.global_id =
            annotation.global_id;
        found->second.physical_tag =
            annotation.physical_tag;
        found->second.annotated = true;
        maximum_face_id =
            std::max(
                maximum_face_id,
                annotation.global_id.value());
    }

    for (auto& [key, face] : faces) {
        (void)key;
        if (!face.annotated) {
            if (maximum_face_id ==
                std::numeric_limits<
                    GlobalEntityId::value_type>::max()) {
                throw std::length_error(
                    "mpmc::mesh::make_linear_mesh_3d: cannot allocate generated face GlobalEntityId");
            }
            ++maximum_face_id;
            face.global_id =
                GlobalEntityId{
                    maximum_face_id};
        }
        if (face.adjacent_cells.size() == 2U &&
            face.physical_tag.is_tagged()) {
            throw std::invalid_argument(
                "mpmc::mesh::make_linear_mesh_3d: interior face cannot carry a physical tag");
        }
    }

    Topology::EntityIds ids;
    ids.vertices.assign(
        vertex_global_ids.begin(),
        vertex_global_ids.end());
    ids.cells.reserve(cells.size());
    for (const auto& cell : cells) {
        ids.cells.push_back(cell.global_id);
    }
    ids.faces.reserve(faces.size());

    std::map<FaceKey, LocalIndex>
        face_local_by_key;
    std::vector<CsrAdjacency::Offset>
        face_vertex_offsets{CsrAdjacency::Offset{0U}};
    std::vector<LocalIndex>
        face_vertex_indices;
    std::vector<CsrAdjacency::Offset>
        face_cell_offsets{CsrAdjacency::Offset{0U}};
    std::vector<LocalIndex>
        face_cell_indices;
    std::vector<PhysicalTag>
        face_physical_tags;
    std::vector<Coordinate3D>
        face_centroids;
    std::vector<double>
        face_areas;
    std::vector<LocalIndex>
        face_owners;
    std::vector<UnitVector3D>
        face_normals;

    face_vertex_offsets.reserve(
        faces.size() + 1U);
    face_cell_offsets.reserve(
        faces.size() + 1U);
    face_physical_tags.reserve(faces.size());
    face_centroids.reserve(faces.size());
    face_areas.reserve(faces.size());
    face_owners.reserve(faces.size());
    face_normals.reserve(faces.size());

    std::size_t face_position = 0U;
    for (const auto& [key, face] : faces) {
        const auto face_local =
            checked_local(
                face_position,
                "mpmc::mesh::make_linear_mesh_3d: face local index overflow");
        face_local_by_key.emplace(
            key, face_local);
        ids.faces.push_back(
            face.global_id);
        face_vertex_indices.insert(
            face_vertex_indices.end(),
            face.ordered_vertices.begin(),
            face.ordered_vertices.end());
        face_vertex_offsets.push_back(
            checked_offset(
                face_vertex_indices.size(),
                "mpmc::mesh::make_linear_mesh_3d: face->vertex CSR overflow"));
        face_cell_indices.insert(
            face_cell_indices.end(),
            face.adjacent_cells.begin(),
            face.adjacent_cells.end());
        face_cell_offsets.push_back(
            checked_offset(
                face_cell_indices.size(),
                "mpmc::mesh::make_linear_mesh_3d: face->cell CSR overflow"));
        face_physical_tags.push_back(
            face.physical_tag);

        const auto owner =
            face.adjacent_cells.front();
        const auto metric =
            face_metric(
                face.ordered_vertices,
                vertex_coordinates_m,
                cell_centroids[
                    static_cast<std::size_t>(
                        owner.value())]);
        face_centroids.push_back(
            metric.centroid_m);
        face_areas.push_back(
            metric.area_m2);
        face_owners.push_back(owner);
        face_normals.push_back(
            metric.owner_unit_normal);
        ++face_position;
    }

    std::vector<CsrAdjacency::Offset>
        cell_vertex_offsets{CsrAdjacency::Offset{0U}};
    std::vector<LocalIndex>
        cell_vertex_indices;
    std::vector<CsrAdjacency::Offset>
        cell_face_offsets{CsrAdjacency::Offset{0U}};
    std::vector<LocalIndex>
        cell_face_indices;
    cell_vertex_offsets.reserve(
        cells.size() + 1U);
    cell_face_offsets.reserve(
        cells.size() + 1U);

    for (std::size_t cell = 0U;
         cell < cells.size();
         ++cell) {
        cell_vertex_indices.insert(
            cell_vertex_indices.end(),
            cells[cell].vertices.begin(),
            cells[cell].vertices.end());
        cell_vertex_offsets.push_back(
            checked_offset(
                cell_vertex_indices.size(),
                "mpmc::mesh::make_linear_mesh_3d: cell->vertex CSR overflow"));

        for (const auto& key :
             cell_face_keys[cell]) {
            const auto found =
                face_local_by_key.find(key);
            if (found ==
                face_local_by_key.end()) {
                throw std::logic_error(
                    "mpmc::mesh::make_linear_mesh_3d: internal face lookup failed");
            }
            cell_face_indices.push_back(
                found->second);
        }
        cell_face_offsets.push_back(
            checked_offset(
                cell_face_indices.size(),
                "mpmc::mesh::make_linear_mesh_3d: cell->face CSR overflow"));
    }

    const std::size_t vertex_count =
        vertex_global_ids.size();
    const std::size_t face_count =
        faces.size();
    const std::size_t cell_count =
        cells.size();

    std::vector<CsrAdjacency> relations;
    relations.emplace_back(
        EntityKind::cell,
        EntityKind::vertex,
        vertex_count,
        std::move(cell_vertex_offsets),
        std::move(cell_vertex_indices));
    relations.emplace_back(
        EntityKind::cell,
        EntityKind::face,
        face_count,
        std::move(cell_face_offsets),
        std::move(cell_face_indices));
    relations.emplace_back(
        EntityKind::face,
        EntityKind::vertex,
        vertex_count,
        std::move(face_vertex_offsets),
        std::move(face_vertex_indices));
    relations.emplace_back(
        EntityKind::face,
        EntityKind::cell,
        cell_count,
        std::move(face_cell_offsets),
        std::move(face_cell_indices));

    Topology topology{
        std::move(ids),
        std::move(relations)};
    auto boundary =
        make_face_boundary_snapshot(
            topology,
            face_physical_tags);
    FaceGeometry3D geometry{
        cell_count,
        std::move(face_centroids),
        std::move(face_areas),
        std::move(face_owners),
        std::move(face_normals)};

    return LinearMesh3D{
        std::move(topology),
        std::vector<Coordinate3D>{
            vertex_coordinates_m.begin(),
            vertex_coordinates_m.end()},
        std::move(cell_volumes),
        std::move(geometry),
        std::move(boundary)};
}

} // namespace mpmc::mesh

#endif // MPMC_MESH_LINEAR_CELL_MESH_3D_HPP
