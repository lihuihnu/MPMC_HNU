#include <mpmc/mesh/active_corner_point.hpp>
#include <mpmc/mesh/cartesian_2d.hpp>
#include <mpmc/mesh/cell_face_geometric_operator_3d.hpp>
#include <mpmc/mesh/corner_point_geometry_3d.hpp>
#include <mpmc/mesh/csr_adjacency.hpp>
#include <mpmc/mesh/dense_field.hpp>
#include <mpmc/mesh/dof_layout.hpp>
#include <mpmc/mesh/dof_numbering.hpp>
#include <mpmc/mesh/entity.hpp>
#include <mpmc/mesh/face_boundary.hpp>
#include <mpmc/mesh/geometry_2d.hpp>
#include <mpmc/mesh/grdecl.hpp>
#include <mpmc/mesh/gmsh_4_1.hpp>
#include <mpmc/mesh/partition_snapshot.hpp>
#include <mpmc/mesh/shared_entity_plan.hpp>
#include <mpmc/mesh/topology.hpp>
#include <mpmc/mesh/vtu.hpp>

#include <array>
#include <cmath>
#include <cstdint>
#include <exception>
#include <iostream>
#include <limits>
#include <source_location>
#include <span>
#include <stdexcept>
#include <string>
#include <string_view>
#include <type_traits>
#include <utility>
#include <vector>

namespace {
namespace mesh = mpmc::mesh;

void require(bool condition,
             std::string_view message,
             std::source_location where = std::source_location::current()) {
    if (!condition) {
        throw std::runtime_error(std::string(where.file_name()) + ":" +
                                 std::to_string(where.line()) + ": " + std::string(message));
    }
}


void require_close(double actual, double expected, double tolerance,
                   std::string_view message,
                   std::source_location where = std::source_location::current()) {
    if (!std::isfinite(actual) || std::abs(actual - expected) > tolerance) {
        throw std::runtime_error(std::string(where.file_name()) + ":" +
                                 std::to_string(where.line()) + ": " +
                                 std::string(message));
    }
}

template <typename Exception, typename Function>
void expect_throw(Function&& function) {
    bool caught = false;
    try {
        std::forward<Function>(function)();
    } catch (const Exception&) {
        caught = true;
    }
    require(caught, "expected exception was not thrown");
}

mesh::CsrAdjacency two_cell_quad_cell_to_vertex() {
    //  v3 ---- v4 ---- v5
    //   |  c0  |  c1  |
    //  v0 ---- v1 ---- v2
    //
    // Connectivity is synthetic software-test data, not a physical mesh dataset.
    std::vector<mesh::CsrAdjacency::Offset> offsets{0U, 4U, 8U};
    std::vector<mesh::LocalIndex> vertices{
        mesh::LocalIndex{0U}, mesh::LocalIndex{1U}, mesh::LocalIndex{4U},
        mesh::LocalIndex{3U}, mesh::LocalIndex{1U}, mesh::LocalIndex{2U},
        mesh::LocalIndex{5U}, mesh::LocalIndex{4U},
    };
    return mesh::CsrAdjacency{mesh::EntityKind::cell, mesh::EntityKind::vertex, 6U,
                              std::move(offsets), std::move(vertices)};
}

void strong_indices() {
    static_assert(sizeof(mesh::LocalIndex) == sizeof(std::uint32_t));
    static_assert(sizeof(mesh::GlobalEntityId) == sizeof(std::uint64_t));
    static_assert(std::is_trivially_copyable_v<mesh::LocalIndex>);
    static_assert(std::is_trivially_copyable_v<mesh::GlobalEntityId>);
    static_assert(!std::is_default_constructible_v<mesh::LocalIndex>);
    static_assert(!std::is_default_constructible_v<mesh::GlobalEntityId>);
    static_assert(!std::is_convertible_v<mesh::LocalIndex, std::uint32_t>);
    static_assert(!std::is_convertible_v<mesh::GlobalEntityId, std::uint64_t>);
    static_assert(!std::is_same_v<mesh::LocalIndex, mesh::GlobalEntityId>);

    constexpr mesh::LocalIndex local{17U};
    constexpr mesh::GlobalEntityId global{9000000001ULL};
    static_assert(local.value() == 17U);
    static_assert(global.value() == 9000000001ULL);
    require(mesh::EntityKind::cell != mesh::EntityKind::face, "entity kinds collapsed");
}

void two_cell_quad() {
    const auto adjacency = two_cell_quad_cell_to_vertex();
    require(adjacency.source_kind() == mesh::EntityKind::cell, "source kind lost");
    require(adjacency.target_kind() == mesh::EntityKind::vertex, "target kind lost");
    require(adjacency.source_count() == 2U, "wrong source count");
    require(adjacency.target_count() == 6U, "wrong target count");
    require(adjacency.entry_count() == 8U, "wrong entry count");

    const std::array<std::uint32_t, 4> expected0{0U, 1U, 4U, 3U};
    const std::array<std::uint32_t, 4> expected1{1U, 2U, 5U, 4U};
    const auto row0 = adjacency.adjacent(mesh::LocalIndex{0U});
    const auto row1 = adjacency.adjacent(mesh::LocalIndex{1U});
    require(row0.size() == expected0.size() && row1.size() == expected1.size(), "row width");
    for (std::size_t i = 0; i < expected0.size(); ++i) {
        require(row0[i].value() == expected0[i], "cell 0 connectivity mismatch");
        require(row1[i].value() == expected1[i], "cell 1 connectivity mismatch");
    }

    const auto offsets = adjacency.offsets();
    const auto indices = adjacency.indices();
    require(offsets.size() == 3U && offsets[0] == 0U && offsets[1] == 4U &&
                offsets[2] == 8U,
            "CSR offsets mismatch");
    require(indices.size() == 8U, "CSR values mismatch");
    require(row0.data() == indices.data(),
            "row view must alias contiguous value storage");
    require(row1.data() == indices.data() + 4,
            "second row must be contiguous without per-row storage");

    expect_throw<std::out_of_range>(
        [&] { (void)adjacency.adjacent(mesh::LocalIndex{2U}); });
}

void invalid_csr() {
    using Offset = mesh::CsrAdjacency::Offset;
    const auto make = [](std::size_t target_count,
                         std::vector<Offset> offsets,
                         std::vector<mesh::LocalIndex> indices) {
        return mesh::CsrAdjacency{mesh::EntityKind::cell, mesh::EntityKind::vertex,
                                  target_count, std::move(offsets), std::move(indices)};
    };

    expect_throw<std::invalid_argument>([&] { (void)make(0U, {}, {}); });
    expect_throw<std::invalid_argument>([&] { (void)make(2U, {1U}, {}); });
    expect_throw<std::invalid_argument>([&] {
        (void)make(3U, {0U, 2U, 1U},
                   {mesh::LocalIndex{0U}, mesh::LocalIndex{1U}});
    });
    expect_throw<std::invalid_argument>([&] {
        (void)make(3U, {0U, 1U},
                   {mesh::LocalIndex{0U}, mesh::LocalIndex{1U}});
    });
    expect_throw<std::out_of_range>([&] {
        (void)make(3U, {0U, 1U}, {mesh::LocalIndex{3U}});
    });

    const auto empty = make(0U, {0U}, {});
    require(empty.source_count() == 0U && empty.entry_count() == 0U,
            "empty relation rejected");
}


mesh::Topology two_cell_quad_topology() {
    mesh::Topology::EntityIds ids;
    ids.vertices = {
        mesh::GlobalEntityId{100U}, mesh::GlobalEntityId{101U},
        mesh::GlobalEntityId{102U}, mesh::GlobalEntityId{103U},
        mesh::GlobalEntityId{104U}, mesh::GlobalEntityId{105U},
    };
    ids.cells = {mesh::GlobalEntityId{200U}, mesh::GlobalEntityId{201U}};
    std::vector<mesh::CsrAdjacency> relations;
    relations.push_back(two_cell_quad_cell_to_vertex());
    return mesh::Topology{std::move(ids), std::move(relations)};
}

void topology_snapshot() {
    const auto topology = two_cell_quad_topology();
    require(topology.entity_count(mesh::EntityKind::vertex) == 6U, "vertex count");
    require(topology.entity_count(mesh::EntityKind::edge) == 0U, "edge count");
    require(topology.entity_count(mesh::EntityKind::face) == 0U, "face count");
    require(topology.entity_count(mesh::EntityKind::cell) == 2U, "cell count");
    require(topology.relation_count() == 1U, "relation count");

    const auto vertex_ids = topology.global_ids(mesh::EntityKind::vertex);
    require(vertex_ids.size() == 6U && vertex_ids.front().value() == 100U &&
                vertex_ids.back().value() == 105U,
            "global ID ordering changed");
    require(topology.global_id(mesh::EntityKind::cell, mesh::LocalIndex{1U}).value() == 201U,
            "local-to-global lookup mismatch");

    require(topology.has_relation(mesh::EntityKind::cell, mesh::EntityKind::vertex),
            "cell-to-vertex relation missing");
    require(!topology.has_relation(mesh::EntityKind::vertex, mesh::EntityKind::cell),
            "unexpected reverse relation materialized");
    const auto& relation =
        topology.relation(mesh::EntityKind::cell, mesh::EntityKind::vertex);
    require(relation.adjacent(mesh::LocalIndex{0U})[2].value() == 4U,
            "snapshot relation payload changed");

    expect_throw<std::out_of_range>([&] {
        (void)topology.global_id(mesh::EntityKind::cell, mesh::LocalIndex{2U});
    });
    expect_throw<std::out_of_range>([&] {
        (void)topology.relation(mesh::EntityKind::vertex, mesh::EntityKind::cell);
    });
}

void topology_invalid() {
    {
        mesh::Topology::EntityIds ids;
        ids.vertices = {mesh::GlobalEntityId{10U}, mesh::GlobalEntityId{10U}};
        expect_throw<std::invalid_argument>(
            [&] { (void)mesh::Topology{std::move(ids), {}}; });
    }
    {
        mesh::Topology::EntityIds ids;
        ids.vertices = {
            mesh::GlobalEntityId{100U}, mesh::GlobalEntityId{101U},
            mesh::GlobalEntityId{102U}, mesh::GlobalEntityId{103U},
            mesh::GlobalEntityId{104U}, mesh::GlobalEntityId{105U},
        };
        ids.cells = {mesh::GlobalEntityId{200U}, mesh::GlobalEntityId{201U}};
        std::vector<mesh::CsrAdjacency> relations;
        relations.emplace_back(
            mesh::EntityKind::cell, mesh::EntityKind::vertex, 6U,
            std::vector<mesh::CsrAdjacency::Offset>{0U, 4U},
            std::vector<mesh::LocalIndex>{
                mesh::LocalIndex{0U}, mesh::LocalIndex{1U},
                mesh::LocalIndex{4U}, mesh::LocalIndex{3U}});
        expect_throw<std::invalid_argument>(
            [&] { (void)mesh::Topology{std::move(ids), std::move(relations)}; });
    }
    {
        mesh::Topology::EntityIds ids;
        ids.vertices = {
            mesh::GlobalEntityId{100U}, mesh::GlobalEntityId{101U},
            mesh::GlobalEntityId{102U}, mesh::GlobalEntityId{103U},
            mesh::GlobalEntityId{104U}, mesh::GlobalEntityId{105U},
        };
        ids.cells = {mesh::GlobalEntityId{200U}, mesh::GlobalEntityId{201U}};
        std::vector<mesh::CsrAdjacency> relations;
        relations.emplace_back(
            mesh::EntityKind::cell, mesh::EntityKind::vertex, 5U,
            std::vector<mesh::CsrAdjacency::Offset>{0U, 4U, 8U},
            std::vector<mesh::LocalIndex>{
                mesh::LocalIndex{0U}, mesh::LocalIndex{1U},
                mesh::LocalIndex{4U}, mesh::LocalIndex{3U},
                mesh::LocalIndex{1U}, mesh::LocalIndex{2U},
                mesh::LocalIndex{4U}, mesh::LocalIndex{3U}});
        expect_throw<std::invalid_argument>(
            [&] { (void)mesh::Topology{std::move(ids), std::move(relations)}; });
    }
    {
        mesh::Topology::EntityIds ids;
        ids.vertices = {
            mesh::GlobalEntityId{100U}, mesh::GlobalEntityId{101U},
            mesh::GlobalEntityId{102U}, mesh::GlobalEntityId{103U},
            mesh::GlobalEntityId{104U}, mesh::GlobalEntityId{105U},
        };
        ids.cells = {mesh::GlobalEntityId{200U}, mesh::GlobalEntityId{201U}};
        std::vector<mesh::CsrAdjacency> relations;
        relations.push_back(two_cell_quad_cell_to_vertex());
        relations.push_back(two_cell_quad_cell_to_vertex());
        expect_throw<std::invalid_argument>(
            [&] { (void)mesh::Topology{std::move(ids), std::move(relations)}; });
    }

    const auto topology = two_cell_quad_topology();
    const auto invalid_kind = static_cast<mesh::EntityKind>(255U);
    expect_throw<std::invalid_argument>([&] { (void)topology.entity_count(invalid_kind); });
    expect_throw<std::invalid_argument>(
        [&] { (void)topology.has_relation(invalid_kind, mesh::EntityKind::vertex); });
}


void cartesian_2d_topology() {
    const auto topology = mesh::make_cartesian_topology_2d(2U, 2U);
    require(topology.entity_count(mesh::EntityKind::vertex) == 9U, "Cartesian vertex count");
    require(topology.entity_count(mesh::EntityKind::edge) == 0U, "2D edge kind must remain empty");
    require(topology.entity_count(mesh::EntityKind::face) == 12U, "Cartesian face count");
    require(topology.entity_count(mesh::EntityKind::cell) == 4U, "Cartesian cell count");
    require(topology.relation_count() == 4U, "Cartesian relation count");

    for (mesh::EntityKind kind : {mesh::EntityKind::vertex, mesh::EntityKind::face,
                                  mesh::EntityKind::cell}) {
        const auto ids = topology.global_ids(kind);
        for (std::size_t i = 0; i < ids.size(); ++i) {
            require(ids[i].value() == static_cast<mesh::GlobalEntityId::value_type>(i),
                    "Cartesian global IDs must be deterministic and zero-based");
        }
    }

    const auto& cell_vertex =
        topology.relation(mesh::EntityKind::cell, mesh::EntityKind::vertex);
    const auto cell3_vertices = cell_vertex.adjacent(mesh::LocalIndex{3U});
    const std::array<std::uint32_t, 4> expected_cell3_vertices{4U, 5U, 8U, 7U};
    require(cell3_vertices.size() == expected_cell3_vertices.size(), "cell vertex arity");
    for (std::size_t i = 0; i < expected_cell3_vertices.size(); ++i) {
        require(cell3_vertices[i].value() == expected_cell3_vertices[i],
                "Cartesian cell-to-vertex numbering mismatch");
    }

    const auto& cell_face =
        topology.relation(mesh::EntityKind::cell, mesh::EntityKind::face);
    const auto cell3_faces = cell_face.adjacent(mesh::LocalIndex{3U});
    const std::array<std::uint32_t, 4> expected_cell3_faces{4U, 5U, 9U, 11U};
    require(cell3_faces.size() == expected_cell3_faces.size(), "cell face arity");
    for (std::size_t i = 0; i < expected_cell3_faces.size(); ++i) {
        require(cell3_faces[i].value() == expected_cell3_faces[i],
                "Cartesian cell-to-face numbering mismatch");
    }

    const auto& face_vertex =
        topology.relation(mesh::EntityKind::face, mesh::EntityKind::vertex);
    const auto horizontal_internal_vertices = face_vertex.adjacent(mesh::LocalIndex{8U});
    require(horizontal_internal_vertices.size() == 2U &&
                horizontal_internal_vertices[0].value() == 3U &&
                horizontal_internal_vertices[1].value() == 4U,
            "horizontal internal face vertices");

    const auto& face_cell =
        topology.relation(mesh::EntityKind::face, mesh::EntityKind::cell);
    const auto vertical_internal_cells = face_cell.adjacent(mesh::LocalIndex{1U});
    require(vertical_internal_cells.size() == 2U &&
                vertical_internal_cells[0].value() == 0U &&
                vertical_internal_cells[1].value() == 1U,
            "vertical internal face sharing");
    const auto horizontal_internal_cells = face_cell.adjacent(mesh::LocalIndex{8U});
    require(horizontal_internal_cells.size() == 2U &&
                horizontal_internal_cells[0].value() == 0U &&
                horizontal_internal_cells[1].value() == 2U,
            "horizontal internal face sharing");

    for (mesh::LocalIndex face : {mesh::LocalIndex{0U}, mesh::LocalIndex{2U},
                                  mesh::LocalIndex{6U}, mesh::LocalIndex{10U}}) {
        require(face_cell.adjacent(face).size() == 1U,
                "boundary face must have exactly one adjacent cell");
    }
}

void cartesian_2d_invalid() {
    expect_throw<std::invalid_argument>(
        [] { (void)mesh::make_cartesian_topology_2d(0U, 1U); });
    expect_throw<std::invalid_argument>(
        [] { (void)mesh::make_cartesian_topology_2d(1U, 0U); });
    expect_throw<std::length_error>([] {
        (void)mesh::make_cartesian_topology_2d(
            std::numeric_limits<std::size_t>::max(), 1U);
    });
    expect_throw<std::length_error>([] {
        (void)mesh::make_cartesian_topology_2d(
            static_cast<std::size_t>(
                std::numeric_limits<mesh::LocalIndex::value_type>::max()),
            1U);
    });
}


void cartesian_2d_geometry() {
    const auto geometry = [] {
        const auto topology = mesh::make_cartesian_topology_2d(2U, 2U);
        const std::array<double, 3> x{-1.0, 1.0, 4.0};
        const std::array<double, 3> y{0.0, 2.0, 5.0};
        return mesh::make_cartesian_geometry_2d(topology, x, y);
    }();

    require(geometry.vertex_count() == 9U, "geometry vertex count");
    require(geometry.cell_count() == 4U, "geometry cell count");
    require(geometry.face_count() == 12U, "geometry face count");

    const auto v8 = geometry.vertex_coordinate_m(mesh::LocalIndex{8U});
    require_close(v8.x_m, 4.0, 0.0, "vertex x coordinate");
    require_close(v8.y_m, 5.0, 0.0, "vertex y coordinate");

    const std::array<double, 4> expected_areas{4.0, 6.0, 6.0, 9.0};
    const std::array<mesh::Coordinate2D, 4> expected_centroids{{
        {0.0, 1.0}, {2.5, 1.0}, {0.0, 3.5}, {2.5, 3.5}}};
    for (std::size_t i = 0; i < expected_areas.size(); ++i) {
        const auto index = mesh::LocalIndex{static_cast<mesh::LocalIndex::value_type>(i)};
        require_close(geometry.cell_area_m2(index), expected_areas[i], 0.0,
                      "cell area mismatch");
        const auto centroid = geometry.cell_centroid_m(index);
        require_close(centroid.x_m, expected_centroids[i].x_m, 0.0,
                      "cell centroid x mismatch");
        require_close(centroid.y_m, expected_centroids[i].y_m, 0.0,
                      "cell centroid y mismatch");
    }

    const auto check_face =
        [&](std::uint32_t face_id, double x, double y, double length,
            std::uint32_t owner, double nx, double ny) {
            const auto face = mesh::LocalIndex{face_id};
            const auto centroid = geometry.face_centroid_m(face);
            const auto normal = geometry.face_owner_unit_normal(face);
            require_close(centroid.x_m, x, 0.0, "face centroid x mismatch");
            require_close(centroid.y_m, y, 0.0, "face centroid y mismatch");
            require_close(geometry.face_length_m(face), length, 0.0,
                          "face length mismatch");
            require(geometry.face_owner(face).value() == owner, "face owner mismatch");
            require_close(normal.x, nx, 0.0, "face normal x mismatch");
            require_close(normal.y, ny, 0.0, "face normal y mismatch");
        };

    check_face(0U, -1.0, 1.0, 2.0, 0U, -1.0, 0.0);
    check_face(1U, 1.0, 1.0, 2.0, 0U, 1.0, 0.0);
    check_face(2U, 4.0, 1.0, 2.0, 1U, 1.0, 0.0);
    check_face(6U, 0.0, 0.0, 2.0, 0U, 0.0, -1.0);
    check_face(8U, 0.0, 2.0, 2.0, 0U, 0.0, 1.0);
    check_face(11U, 2.5, 5.0, 3.0, 3U, 0.0, 1.0);

    expect_throw<std::out_of_range>(
        [&] { (void)geometry.face_length_m(mesh::LocalIndex{12U}); });
}

void cartesian_2d_geometry_invalid() {
    const auto topology = mesh::make_cartesian_topology_2d(1U, 1U);
    const std::array<double, 2> valid_axis{0.0, 1.0};
    const std::array<double, 1> short_axis{0.0};
    const std::array<double, 2> duplicate_axis{0.0, 0.0};
    const std::array<double, 2> descending_axis{1.0, 0.0};
    const std::array<double, 2> nonfinite_axis{
        0.0, std::numeric_limits<double>::infinity()};

    expect_throw<std::invalid_argument>([&] {
        (void)mesh::make_cartesian_geometry_2d(topology, short_axis, valid_axis);
    });
    expect_throw<std::invalid_argument>([&] {
        (void)mesh::make_cartesian_geometry_2d(topology, duplicate_axis, valid_axis);
    });
    expect_throw<std::invalid_argument>([&] {
        (void)mesh::make_cartesian_geometry_2d(topology, descending_axis, valid_axis);
    });
    expect_throw<std::invalid_argument>([&] {
        (void)mesh::make_cartesian_geometry_2d(topology, nonfinite_axis, valid_axis);
    });

    const auto wrong_size_topology = mesh::make_cartesian_topology_2d(2U, 1U);
    expect_throw<std::invalid_argument>([&] {
        (void)mesh::make_cartesian_geometry_2d(
            wrong_size_topology, valid_axis, valid_axis);
    });

    mesh::Topology::EntityIds ids;
    for (std::uint64_t id = 0U; id < 4U; ++id) {
        ids.vertices.emplace_back(id);
        ids.faces.emplace_back(id);
    }
    ids.cells.emplace_back(0U);
    std::vector<mesh::CsrAdjacency> relations;
    relations.emplace_back(
        mesh::EntityKind::cell, mesh::EntityKind::vertex, 4U,
        std::vector<mesh::CsrAdjacency::Offset>{0U, 4U},
        std::vector<mesh::LocalIndex>{
            mesh::LocalIndex{0U}, mesh::LocalIndex{2U},
            mesh::LocalIndex{3U}, mesh::LocalIndex{1U}});
    relations.emplace_back(
        mesh::EntityKind::cell, mesh::EntityKind::face, 4U,
        std::vector<mesh::CsrAdjacency::Offset>{0U, 4U},
        std::vector<mesh::LocalIndex>{
            mesh::LocalIndex{0U}, mesh::LocalIndex{1U},
            mesh::LocalIndex{2U}, mesh::LocalIndex{3U}});
    relations.emplace_back(
        mesh::EntityKind::face, mesh::EntityKind::vertex, 4U,
        std::vector<mesh::CsrAdjacency::Offset>{0U, 2U, 4U, 6U, 8U},
        std::vector<mesh::LocalIndex>{
            mesh::LocalIndex{0U}, mesh::LocalIndex{2U},
            mesh::LocalIndex{1U}, mesh::LocalIndex{3U},
            mesh::LocalIndex{0U}, mesh::LocalIndex{1U},
            mesh::LocalIndex{2U}, mesh::LocalIndex{3U}});
    relations.emplace_back(
        mesh::EntityKind::face, mesh::EntityKind::cell, 1U,
        std::vector<mesh::CsrAdjacency::Offset>{0U, 1U, 2U, 3U, 4U},
        std::vector<mesh::LocalIndex>{
            mesh::LocalIndex{0U}, mesh::LocalIndex{0U},
            mesh::LocalIndex{0U}, mesh::LocalIndex{0U}});
    const mesh::Topology noncanonical{std::move(ids), std::move(relations)};
    expect_throw<std::invalid_argument>([&] {
        (void)mesh::make_cartesian_geometry_2d(
            noncanonical, valid_axis, valid_axis);
    });

    expect_throw<std::invalid_argument>([] {
        (void)mesh::Geometry2D{
            {},
            {mesh::Coordinate2D{0.0, 0.0}},
            {},
            {},
            {},
            {},
            {}};
    });
}


void face_boundary_snapshot() {
    const auto topology = mesh::make_cartesian_topology_2d(2U, 2U);
    std::vector<mesh::PhysicalTag> tags(12U, mesh::PhysicalTag{0U});
    for (std::uint32_t face : {0U, 3U}) tags[face] = mesh::PhysicalTag{11U};
    for (std::uint32_t face : {2U, 5U}) tags[face] = mesh::PhysicalTag{12U};
    for (std::uint32_t face : {6U, 7U}) tags[face] = mesh::PhysicalTag{21U};
    for (std::uint32_t face : {10U, 11U}) tags[face] = mesh::PhysicalTag{22U};

    const auto snapshot = mesh::make_face_boundary_snapshot(topology, tags);
    require(snapshot.face_count() == 12U, "boundary snapshot face count");
    require(snapshot.boundary_face_count() == 8U, "boundary face count");
    require(snapshot.interior_face_count() == 4U, "interior face count");

    for (std::uint32_t face : {1U, 4U, 8U, 9U}) {
        const auto index = mesh::LocalIndex{face};
        require(!snapshot.is_boundary(index), "interior face misclassified");
        require(!snapshot.has_physical_tag(index), "interior face unexpectedly tagged");
    }
    require(snapshot.is_boundary(mesh::LocalIndex{0U}), "boundary face misclassified");
    require(snapshot.physical_tag(mesh::LocalIndex{0U}).value() == 11U, "left tag");
    require(snapshot.physical_tag(mesh::LocalIndex{2U}).value() == 12U, "right tag");
    require(snapshot.physical_tag(mesh::LocalIndex{6U}).value() == 21U, "bottom tag");
    require(snapshot.physical_tag(mesh::LocalIndex{10U}).value() == 22U, "top tag");

    const auto untagged = mesh::make_face_boundary_snapshot(topology);
    require(untagged.boundary_face_count() == 8U, "untagged classification changed");
    require(!untagged.has_physical_tag(mesh::LocalIndex{0U}), "default tag must be absent");

    expect_throw<std::out_of_range>(
        [&] { (void)snapshot.classification(mesh::LocalIndex{12U}); });
}

void face_boundary_invalid() {
    const auto topology = mesh::make_cartesian_topology_2d(1U, 1U);

    expect_throw<std::invalid_argument>([&] {
        const std::array<mesh::PhysicalTag, 3> tags{
            mesh::PhysicalTag{0U}, mesh::PhysicalTag{0U}, mesh::PhysicalTag{0U}};
        (void)mesh::make_face_boundary_snapshot(topology, tags);
    });

    expect_throw<std::invalid_argument>([&] {
        const auto two_by_two = mesh::make_cartesian_topology_2d(2U, 2U);
        std::vector<mesh::PhysicalTag> tags(12U, mesh::PhysicalTag{0U});
        tags[1] = mesh::PhysicalTag{99U};
        (void)mesh::make_face_boundary_snapshot(two_by_two, tags);
    });

    mesh::Topology::EntityIds missing_ids;
    missing_ids.faces.emplace_back(1U);
    const mesh::Topology missing_relation{std::move(missing_ids), {}};
    expect_throw<std::invalid_argument>(
        [&] { (void)mesh::make_face_boundary_snapshot(missing_relation); });

    mesh::Topology::EntityIds zero_ids;
    zero_ids.faces.emplace_back(1U);
    zero_ids.cells.emplace_back(2U);
    std::vector<mesh::CsrAdjacency> zero_relations;
    zero_relations.emplace_back(
        mesh::EntityKind::face, mesh::EntityKind::cell, 1U,
        std::vector<mesh::CsrAdjacency::Offset>{0U, 0U},
        std::vector<mesh::LocalIndex>{});
    const mesh::Topology zero_adjacent{std::move(zero_ids), std::move(zero_relations)};
    expect_throw<std::invalid_argument>(
        [&] { (void)mesh::make_face_boundary_snapshot(zero_adjacent); });

    mesh::Topology::EntityIds three_ids;
    three_ids.faces.emplace_back(1U);
    three_ids.cells = {mesh::GlobalEntityId{10U}, mesh::GlobalEntityId{11U},
                       mesh::GlobalEntityId{12U}};
    std::vector<mesh::CsrAdjacency> three_relations;
    three_relations.emplace_back(
        mesh::EntityKind::face, mesh::EntityKind::cell, 3U,
        std::vector<mesh::CsrAdjacency::Offset>{0U, 3U},
        std::vector<mesh::LocalIndex>{
            mesh::LocalIndex{0U}, mesh::LocalIndex{1U}, mesh::LocalIndex{2U}});
    const mesh::Topology three_adjacent{std::move(three_ids), std::move(three_relations)};
    expect_throw<std::invalid_argument>(
        [&] { (void)mesh::make_face_boundary_snapshot(three_adjacent); });

    expect_throw<std::invalid_argument>([] {
        (void)mesh::FaceBoundarySnapshot{
            {mesh::FaceClassification::interior},
            {mesh::PhysicalTag{7U}}};
    });
    expect_throw<std::invalid_argument>([] {
        (void)mesh::FaceBoundarySnapshot{
            {static_cast<mesh::FaceClassification>(255U)},
            {mesh::PhysicalTag{0U}}};
    });
}


mesh::DenseFieldMetadata synthetic_field_metadata(std::string id, std::string unit) {
    return {
        std::move(id),
        std::move(unit),
        {mesh::FieldSourceKind::synthetic_test,
         "test://mesh/dense-field",
         "fixture-v1",
         "core_test.cpp"}};
}

void dense_field_snapshot() {
    const auto topology = mesh::make_cartesian_topology_2d(2U, 2U);

    const auto cell_scalar = mesh::DenseFieldSnapshot::create(
        topology,
        mesh::EntityKind::cell,
        1U,
        {10.0, 20.0, 30.0, 40.0},
        synthetic_field_metadata("test.cell.scalar", "1"));
    require(cell_scalar.location() == mesh::EntityKind::cell, "cell field location");
    require(cell_scalar.entity_count() == 4U, "cell field entity count");
    require(cell_scalar.component_count() == 1U, "cell scalar component count");
    require(cell_scalar.value_count() == 4U, "cell scalar value count");
    require_close(cell_scalar.value(mesh::LocalIndex{2U}, 0U), 30.0, 0.0,
                  "cell scalar lookup");
    require(cell_scalar.metadata().id == "test.cell.scalar", "field id lost");
    require(cell_scalar.metadata().unit == "1", "field unit lost");
    require(cell_scalar.metadata().source.kind == mesh::FieldSourceKind::synthetic_test,
            "field source kind lost");
    require(cell_scalar.metadata().source.reference == "test://mesh/dense-field",
            "field source reference lost");
    require(cell_scalar.metadata().source.revision == "fixture-v1",
            "field source revision lost");

    std::vector<double> face_values;
    face_values.reserve(24U);
    for (std::size_t face = 0; face < 12U; ++face) {
        face_values.push_back(static_cast<double>(100U + face));
        face_values.push_back(static_cast<double>(200U + face));
    }
    const auto face_vector = mesh::DenseFieldSnapshot::create(
        topology,
        mesh::EntityKind::face,
        2U,
        std::move(face_values),
        synthetic_field_metadata("test.face.vector", "m/s"));
    const auto face3 = face_vector.entity_values(mesh::LocalIndex{3U});
    require(face3.size() == 2U, "face vector width");
    require_close(face3[0], 103.0, 0.0, "face vector component zero");
    require_close(face3[1], 203.0, 0.0, "face vector component one");
    require(face3.data() == face_vector.values().data() + 6,
            "entity-major field storage must be contiguous");

    std::vector<double> vertex_values;
    vertex_values.reserve(27U);
    for (std::size_t vertex = 0; vertex < 9U; ++vertex) {
        vertex_values.push_back(static_cast<double>(vertex));
        vertex_values.push_back(static_cast<double>(vertex + 10U));
        vertex_values.push_back(static_cast<double>(vertex + 20U));
    }
    const auto vertex_vector = mesh::DenseFieldSnapshot::create(
        topology,
        mesh::EntityKind::vertex,
        3U,
        std::move(vertex_values),
        synthetic_field_metadata("test.vertex.vector", "m"));
    require(vertex_vector.entity_count() == 9U, "vertex field entity count");
    require(vertex_vector.component_count() == 3U, "vertex field component count");
    require_close(vertex_vector.value(mesh::LocalIndex{8U}, 2U), 28.0, 0.0,
                  "vertex component lookup");

    expect_throw<std::out_of_range>(
        [&] { (void)cell_scalar.value(mesh::LocalIndex{4U}, 0U); });
    expect_throw<std::out_of_range>(
        [&] { (void)cell_scalar.value(mesh::LocalIndex{0U}, 1U); });
}

void dense_field_invalid() {
    const auto topology = mesh::make_cartesian_topology_2d(1U, 1U);

    expect_throw<std::invalid_argument>([&] {
        (void)mesh::DenseFieldSnapshot::create(
            topology, mesh::EntityKind::cell, 0U, {},
            synthetic_field_metadata("test.zero-components", "1"));
    });
    expect_throw<std::invalid_argument>([&] {
        (void)mesh::DenseFieldSnapshot::create(
            topology, mesh::EntityKind::edge, 1U, {},
            synthetic_field_metadata("test.edge", "1"));
    });
    expect_throw<std::invalid_argument>([&] {
        (void)mesh::DenseFieldSnapshot::create(
            topology, static_cast<mesh::EntityKind>(255U), 1U, {},
            synthetic_field_metadata("test.bad-location", "1"));
    });
    expect_throw<std::invalid_argument>([&] {
        (void)mesh::DenseFieldSnapshot::create(
            topology, mesh::EntityKind::cell, 2U, {1.0},
            synthetic_field_metadata("test.bad-size", "1"));
    });
    expect_throw<std::invalid_argument>([&] {
        (void)mesh::DenseFieldSnapshot::create(
            topology, mesh::EntityKind::cell, 1U,
            {std::numeric_limits<double>::quiet_NaN()},
            synthetic_field_metadata("test.nan", "1"));
    });
    expect_throw<std::invalid_argument>([&] {
        (void)mesh::DenseFieldSnapshot::create(
            topology, mesh::EntityKind::cell, 1U,
            {std::numeric_limits<double>::infinity()},
            synthetic_field_metadata("test.inf", "1"));
    });
    expect_throw<std::length_error>([&] {
        (void)mesh::DenseFieldSnapshot::create(
            topology, mesh::EntityKind::face,
            std::numeric_limits<std::size_t>::max(), {},
            synthetic_field_metadata("test.overflow", "1"));
    });

    expect_throw<std::invalid_argument>([&] {
        (void)mesh::DenseFieldSnapshot::create(
            topology, mesh::EntityKind::cell, 1U, {1.0},
            mesh::DenseFieldMetadata{});
    });

    for (int invalid = 0; invalid < 4; ++invalid) {
        auto metadata = synthetic_field_metadata("test.metadata", "Pa");
        if (invalid == 0) metadata.id = " \t";
        if (invalid == 1) metadata.unit.clear();
        if (invalid == 2) metadata.source.reference.clear();
        if (invalid == 3) metadata.source.revision.clear();
        expect_throw<std::invalid_argument>([&] {
            (void)mesh::DenseFieldSnapshot::create(
                topology, mesh::EntityKind::cell, 1U, {1.0},
                std::move(metadata));
        });
    }

    auto bad_kind = synthetic_field_metadata("test.bad-source-kind", "1");
    bad_kind.source.kind = static_cast<mesh::FieldSourceKind>(255U);
    expect_throw<std::invalid_argument>([&] {
        (void)mesh::DenseFieldSnapshot::create(
            topology, mesh::EntityKind::cell, 1U, {1.0},
            std::move(bad_kind));
    });

    auto bad_locator = synthetic_field_metadata("test.bad-locator", "1");
    bad_locator.source.locator = std::string{"row\0hidden", 10U};
    expect_throw<std::invalid_argument>([&] {
        (void)mesh::DenseFieldSnapshot::create(
            topology, mesh::EntityKind::cell, 1U, {1.0},
            std::move(bad_locator));
    });
}


std::string gmsh_mixed_triangle_quad_fixture() {
    return R"msh($MeshFormat
4.1 0 8
$EndMeshFormat
$PhysicalNames
5
1 11 "left"
1 12 "bottom"
1 13 "right"
1 14 "top"
2 21 "domain"
$EndPhysicalNames
$Entities
0 5 1 0
1 0 0 0 1 0 0 1 12 0
2 0 0 0 1 1 0 1 11 0
3 1 0 0 3 0 0 1 12 0
4 3 0 0 3 1 0 1 13 0
5 1 1 0 3 1 0 1 14 0
100 0 0 0 3 1 0 1 21 5 1 2 3 4 5
$EndEntities
$Nodes
1 5 10 50
2 100 0 5
50
10
40
20
30
3 1 0
0 0 0
3 0 0
1 0 0
1 1 0
$EndNodes
$Elements
7 7 101 202
1 1 1 1
101 10 20
1 2 1 1
103 30 10
1 3 1 1
104 20 40
1 4 1 1
105 40 50
1 5 1 1
106 50 30
2 100 2 1
201 10 20 30
2 100 3 1
202 20 40 50 30
$EndElements
)msh";
}

void gmsh_4_1_import() {
    const auto imported =
        mesh::import_gmsh_4_1_ascii(
            gmsh_mixed_triangle_quad_fixture(),
            2.0);

    const auto& topology = imported.topology;
    require(
        topology.entity_count(mesh::EntityKind::vertex) == 5U,
        "Gmsh vertex count");
    require(
        topology.entity_count(mesh::EntityKind::edge) == 0U,
        "Gmsh edge kind must remain empty in 2D");
    require(
        topology.entity_count(mesh::EntityKind::face) == 6U,
        "Gmsh face count");
    require(
        topology.entity_count(mesh::EntityKind::cell) == 2U,
        "Gmsh cell count");

    const std::array<std::uint64_t, 5> expected_vertices{
        10U, 20U, 30U, 40U, 50U};
    const auto vertex_ids =
        topology.global_ids(mesh::EntityKind::vertex);
    for (std::size_t i = 0U;
         i < expected_vertices.size();
         ++i) {
        require(
            vertex_ids[i].value() ==
                expected_vertices[i],
            "Gmsh sparse node tags must become sorted stable vertex IDs");
    }

    const std::array<std::uint64_t, 6> expected_faces{
        101U, 103U, 104U, 105U, 106U, 203U};
    const auto face_ids =
        topology.global_ids(mesh::EntityKind::face);
    for (std::size_t i = 0U;
         i < expected_faces.size();
         ++i) {
        require(
            face_ids[i].value() ==
                expected_faces[i],
            "Gmsh boundary/generated face stable IDs");
    }
    require(
        topology.global_id(
            mesh::EntityKind::cell,
            mesh::LocalIndex{0U}).value() == 201U &&
        topology.global_id(
            mesh::EntityKind::cell,
            mesh::LocalIndex{1U}).value() == 202U,
        "Gmsh cell element tags must be stable cell IDs");

    const auto& cell_vertices =
        topology.relation(
            mesh::EntityKind::cell,
            mesh::EntityKind::vertex);
    require(
        cell_vertices.adjacent(
            mesh::LocalIndex{0U}).size() == 3U,
        "Gmsh triangle connectivity width");
    require(
        cell_vertices.adjacent(
            mesh::LocalIndex{1U}).size() == 4U,
        "Gmsh quad connectivity width");

    const auto& cell_faces =
        topology.relation(
            mesh::EntityKind::cell,
            mesh::EntityKind::face);
    require(
        cell_faces.adjacent(
            mesh::LocalIndex{0U}).size() == 3U &&
        cell_faces.adjacent(
            mesh::LocalIndex{1U}).size() == 4U,
        "Gmsh triangle/quad cell-to-face widths");

    const auto& face_cells =
        topology.relation(
            mesh::EntityKind::face,
            mesh::EntityKind::cell);
    require(
        face_cells.adjacent(
            mesh::LocalIndex{5U}).size() == 2U,
        "generated internal face must connect both cells");

    const auto& boundary =
        imported.face_boundary;
    require(
        boundary.face_count() == 6U &&
        boundary.boundary_face_count() == 5U &&
        boundary.interior_face_count() == 1U,
        "Gmsh boundary classification counts");
    const std::array<std::uint32_t, 6> expected_tags{
        12U, 11U, 12U, 13U, 14U, 0U};
    for (std::size_t face = 0U;
         face < expected_tags.size();
         ++face) {
        const auto local =
            mesh::LocalIndex{
                static_cast<
                    mesh::LocalIndex::value_type>(
                        face)};
        require(
            boundary.physical_tag(local).value() ==
                expected_tags[face],
            "Gmsh curve Physical Group to PhysicalTag mapping");
    }
    require(
        boundary.classification(
            mesh::LocalIndex{5U}) ==
            mesh::FaceClassification::interior &&
        !boundary.has_physical_tag(
            mesh::LocalIndex{5U}),
        "generated internal face must remain untagged interior");

    const auto& geometry = imported.geometry;
    const auto v0 =
        geometry.vertex_coordinate_m(
            mesh::LocalIndex{0U});
    const auto v4 =
        geometry.vertex_coordinate_m(
            mesh::LocalIndex{4U});
    require_close(v0.x_m, 0.0, 0.0,
                  "scaled Gmsh vertex x0");
    require_close(v0.y_m, 0.0, 0.0,
                  "scaled Gmsh vertex y0");
    require_close(v4.x_m, 6.0, 0.0,
                  "scaled Gmsh vertex x4");
    require_close(v4.y_m, 2.0, 0.0,
                  "scaled Gmsh vertex y4");

    const auto triangle_centroid =
        geometry.cell_centroid_m(
            mesh::LocalIndex{0U});
    require_close(
        triangle_centroid.x_m,
        4.0 / 3.0,
        1.0e-14,
        "Gmsh triangle centroid x");
    require_close(
        triangle_centroid.y_m,
        2.0 / 3.0,
        1.0e-14,
        "Gmsh triangle centroid y");
    require_close(
        geometry.cell_area_m2(
            mesh::LocalIndex{0U}),
        2.0,
        1.0e-14,
        "Gmsh triangle area");

    const auto quad_centroid =
        geometry.cell_centroid_m(
            mesh::LocalIndex{1U});
    require_close(
        quad_centroid.x_m,
        4.0,
        1.0e-14,
        "Gmsh quad centroid x");
    require_close(
        quad_centroid.y_m,
        1.0,
        1.0e-14,
        "Gmsh quad centroid y");
    require_close(
        geometry.cell_area_m2(
            mesh::LocalIndex{1U}),
        8.0,
        1.0e-14,
        "Gmsh quad area");

    const std::array<double, 6> expected_lengths{
        2.0,
        std::sqrt(8.0),
        4.0,
        2.0,
        4.0,
        2.0};
    for (std::size_t face = 0U;
         face < expected_lengths.size();
         ++face) {
        require_close(
            geometry.face_length_m(
                mesh::LocalIndex{
                    static_cast<
                        mesh::LocalIndex::value_type>(
                            face)}),
            expected_lengths[face],
            1.0e-14,
            "Gmsh face length");
    }

    require(
        imported.physical_names.size() == 5U,
        "Gmsh PhysicalNames retention");
    require(
        imported.physical_names.front().dimension == 1 &&
        imported.physical_names.front().tag == 11U &&
        imported.physical_names.front().name == "left",
        "Gmsh first PhysicalName");
    require(
        imported.physical_names.back().dimension == 2 &&
        imported.physical_names.back().tag == 21U &&
        imported.physical_names.back().name == "domain",
        "Gmsh surface PhysicalName retention");

    require(
        imported.cell_physical_groups.size() == 2U,
        "Gmsh cell Physical Group count");
    for (std::size_t cell = 0U;
         cell < imported.cell_physical_groups.size();
         ++cell) {
        const auto& groups =
            imported.cell_physical_groups[cell];
        require(
            groups.cell_global_id.value() ==
                201U + cell &&
            groups.physical_tags.size() == 1U &&
            groups.physical_tags.front() == 21U,
            "Gmsh surface Physical Group retention per cell");
    }
}

std::string gmsh_roundtrip_surface_groups_fixture() {
    auto content =
        gmsh_mixed_triangle_quad_fixture();

    const auto physical_count =
        content.find("$PhysicalNames\n5\n");
    require(
        physical_count != std::string::npos,
        "round-trip PhysicalNames count marker");
    content.replace(
        physical_count,
        std::string{"$PhysicalNames\n5\n"}.size(),
        "$PhysicalNames\n6\n");

    const auto physical_end =
        content.find("$EndPhysicalNames");
    require(
        physical_end != std::string::npos,
        "round-trip PhysicalNames end marker");
    content.insert(
        physical_end,
        "2 22 \"material\"\n");

    const std::string surface_entity =
        "100 0 0 0 3 1 0 1 21 5 1 2 3 4 5";
    const auto surface =
        content.find(surface_entity);
    require(
        surface != std::string::npos,
        "round-trip surface entity marker");
    content.replace(
        surface,
        surface_entity.size(),
        "100 0 0 0 3 1 0 2 21 22 5 1 2 3 4 5");
    return content;
}

void require_topology_equal(
    const mesh::Topology& actual,
    const mesh::Topology& expected) {
    for (const auto kind :
         {mesh::EntityKind::vertex,
          mesh::EntityKind::edge,
          mesh::EntityKind::face,
          mesh::EntityKind::cell}) {
        require(
            actual.entity_count(kind) ==
                expected.entity_count(kind),
            "round-trip topology entity count");
        const auto actual_ids =
            actual.global_ids(kind);
        const auto expected_ids =
            expected.global_ids(kind);
        require(
            actual_ids.size() ==
                expected_ids.size(),
            "round-trip topology GlobalEntityId count");
        for (std::size_t i = 0U;
             i < actual_ids.size();
             ++i) {
            require(
                actual_ids[i] ==
                    expected_ids[i],
                "round-trip topology GlobalEntityId");
        }
    }

    for (const auto& relation :
         {std::pair{
              mesh::EntityKind::cell,
              mesh::EntityKind::vertex},
          std::pair{
              mesh::EntityKind::cell,
              mesh::EntityKind::face},
          std::pair{
              mesh::EntityKind::face,
              mesh::EntityKind::vertex},
          std::pair{
              mesh::EntityKind::face,
              mesh::EntityKind::cell}}) {
        require(
            actual.has_relation(
                relation.first,
                relation.second) ==
                expected.has_relation(
                    relation.first,
                    relation.second),
            "round-trip topology relation presence");
        const auto& actual_relation =
            actual.relation(
                relation.first,
                relation.second);
        const auto& expected_relation =
            expected.relation(
                relation.first,
                relation.second);
        require(
            actual_relation.offsets().size() ==
                expected_relation.offsets().size() &&
            actual_relation.indices().size() ==
                expected_relation.indices().size(),
            "round-trip topology relation storage size");
        for (std::size_t i = 0U;
             i < actual_relation.offsets().size();
             ++i) {
            require(
                actual_relation.offsets()[i] ==
                    expected_relation.offsets()[i],
                "round-trip topology relation offsets");
        }
        for (std::size_t i = 0U;
             i < actual_relation.indices().size();
             ++i) {
            require(
                actual_relation.indices()[i] ==
                    expected_relation.indices()[i],
                "round-trip topology relation indices");
        }
    }
}

void require_geometry_equal(
    const mesh::Geometry2D& actual,
    const mesh::Geometry2D& expected) {
    constexpr double tolerance = 1.0e-13;
    require(
        actual.vertex_count() ==
            expected.vertex_count() &&
        actual.face_count() ==
            expected.face_count() &&
        actual.cell_count() ==
            expected.cell_count(),
        "round-trip geometry entity counts");

    for (std::size_t vertex = 0U;
         vertex < actual.vertex_count();
         ++vertex) {
        const auto local =
            mesh::LocalIndex{
                static_cast<
                    mesh::LocalIndex::value_type>(
                        vertex)};
        const auto a =
            actual.vertex_coordinate_m(local);
        const auto e =
            expected.vertex_coordinate_m(local);
        require_close(
            a.x_m, e.x_m, tolerance,
            "round-trip vertex x");
        require_close(
            a.y_m, e.y_m, tolerance,
            "round-trip vertex y");
    }

    for (std::size_t cell = 0U;
         cell < actual.cell_count();
         ++cell) {
        const auto local =
            mesh::LocalIndex{
                static_cast<
                    mesh::LocalIndex::value_type>(
                        cell)};
        const auto a =
            actual.cell_centroid_m(local);
        const auto e =
            expected.cell_centroid_m(local);
        require_close(
            a.x_m, e.x_m, tolerance,
            "round-trip cell centroid x");
        require_close(
            a.y_m, e.y_m, tolerance,
            "round-trip cell centroid y");
        require_close(
            actual.cell_area_m2(local),
            expected.cell_area_m2(local),
            tolerance,
            "round-trip cell area");
    }

    for (std::size_t face = 0U;
         face < actual.face_count();
         ++face) {
        const auto local =
            mesh::LocalIndex{
                static_cast<
                    mesh::LocalIndex::value_type>(
                        face)};
        const auto a =
            actual.face_centroid_m(local);
        const auto e =
            expected.face_centroid_m(local);
        require_close(
            a.x_m, e.x_m, tolerance,
            "round-trip face centroid x");
        require_close(
            a.y_m, e.y_m, tolerance,
            "round-trip face centroid y");
        require_close(
            actual.face_length_m(local),
            expected.face_length_m(local),
            tolerance,
            "round-trip face length");
        require(
            actual.face_owner(local) ==
                expected.face_owner(local),
            "round-trip face owner");
        const auto an =
            actual.face_owner_unit_normal(local);
        const auto en =
            expected.face_owner_unit_normal(local);
        require_close(
            an.x, en.x, tolerance,
            "round-trip face normal x");
        require_close(
            an.y, en.y, tolerance,
            "round-trip face normal y");
    }
}

void require_boundary_equal(
    const mesh::FaceBoundarySnapshot& actual,
    const mesh::FaceBoundarySnapshot& expected) {
    require(
        actual.face_count() ==
            expected.face_count(),
        "round-trip boundary face count");
    for (std::size_t face = 0U;
         face < actual.face_count();
         ++face) {
        const auto local =
            mesh::LocalIndex{
                static_cast<
                    mesh::LocalIndex::value_type>(
                        face)};
        require(
            actual.classification(local) ==
                expected.classification(local),
            "round-trip face classification");
        require(
            actual.physical_tag(local) ==
                expected.physical_tag(local),
            "round-trip PhysicalTag");
    }
}

void require_physical_metadata_equal(
    const mesh::Gmsh41ImportResult& actual,
    const mesh::Gmsh41ImportResult& expected) {
    require(
        actual.physical_names.size() ==
            expected.physical_names.size(),
        "round-trip PhysicalNames count");
    for (std::size_t i = 0U;
         i < actual.physical_names.size();
         ++i) {
        const auto& a =
            actual.physical_names[i];
        const auto& e =
            expected.physical_names[i];
        require(
            a.dimension == e.dimension &&
                a.tag == e.tag &&
                a.name == e.name,
            "round-trip PhysicalName");
    }

    require(
        actual.cell_physical_groups.size() ==
            expected.cell_physical_groups.size(),
        "round-trip surface Physical Group record count");
    for (std::size_t i = 0U;
         i < actual.cell_physical_groups.size();
         ++i) {
        const auto& a =
            actual.cell_physical_groups[i];
        const auto& e =
            expected.cell_physical_groups[i];
        require(
            a.cell_global_id ==
                e.cell_global_id &&
                a.physical_tags ==
                    e.physical_tags,
            "round-trip surface Physical Groups");
    }
}

void gmsh_4_1_roundtrip() {
    const auto first =
        mesh::import_gmsh_4_1_ascii(
            gmsh_roundtrip_surface_groups_fixture(),
            2.0);

    require(
        first.cell_physical_groups.size() == 2U,
        "round-trip source surface group count");
    for (const auto& record :
         first.cell_physical_groups) {
        require(
            record.physical_tags ==
                std::vector<std::uint32_t>{
                    21U, 22U},
            "round-trip source multiple surface groups");
    }

    const std::string exported =
        mesh::export_gmsh_4_1_ascii(first);
    require(
        exported.find("$MeshFormat\n4.1 0 8\n") !=
            std::string::npos,
        "round-trip export MSH4.1 ASCII header");
    require(
        exported.find("2 22 \"material\"") !=
            std::string::npos,
        "round-trip export PhysicalName");
    require(
        exported.find("\n203 ") !=
            std::string::npos,
        "round-trip export generated internal face element tag");

    const auto second =
        mesh::import_gmsh_4_1_ascii(
            exported,
            1.0);

    require_topology_equal(
        second.topology,
        first.topology);
    require_geometry_equal(
        second.geometry,
        first.geometry);
    require_boundary_equal(
        second.face_boundary,
        first.face_boundary);
    require_physical_metadata_equal(
        second,
        first);
}

void gmsh_4_1_invalid() {
    const auto valid =
        gmsh_mixed_triangle_quad_fixture();

    const auto export_source =
        mesh::import_gmsh_4_1_ascii(
            valid, 1.0);

    {
        auto bad = export_source;
        bad.physical_names.push_back(
            bad.physical_names.front());
        expect_throw<std::invalid_argument>([&] {
            (void)mesh::export_gmsh_4_1_ascii(
                bad);
        });
    }

    {
        auto bad = export_source;
        bad.cell_physical_groups.pop_back();
        expect_throw<std::invalid_argument>([&] {
            (void)mesh::export_gmsh_4_1_ascii(
                bad);
        });
    }

    {
        auto bad = export_source;
        bad.physical_names.push_back(
            mesh::GmshPhysicalName{
                1, 999U, "unused"});
        expect_throw<std::invalid_argument>([&] {
            (void)mesh::export_gmsh_4_1_ascii(
                bad);
        });
    }

    expect_throw<std::invalid_argument>([&] {
        (void)mesh::import_gmsh_4_1_ascii(
            valid, 0.0);
    });

    {
        auto binary = valid;
        const auto where =
            binary.find("4.1 0 8");
        require(where != std::string::npos,
                "binary fixture marker");
        binary.replace(
            where, std::string{"4.1 0 8"}.size(),
            "4.1 1 8");
        expect_throw<std::invalid_argument>([&] {
            (void)mesh::import_gmsh_4_1_ascii(
                binary, 1.0);
        });
    }

    {
        auto high_order = valid;
        const auto where =
            high_order.find("2 100 2 1\n201");
        require(where != std::string::npos,
                "high-order fixture marker");
        high_order.replace(
            where,
            std::string{"2 100 2 1"}.size(),
            "2 100 9 1");
        expect_throw<std::invalid_argument>([&] {
            (void)mesh::import_gmsh_4_1_ascii(
                high_order, 1.0);
        });
    }

    {
        auto missing_node = valid;
        const auto where =
            missing_node.find("201 10 20 30");
        require(where != std::string::npos,
                "missing-node fixture marker");
        missing_node.replace(
            where,
            std::string{"201 10 20 30"}.size(),
            "201 10 20 999");
        expect_throw<std::invalid_argument>([&] {
            (void)mesh::import_gmsh_4_1_ascii(
                missing_node, 1.0);
        });
    }

    {
        auto multi_group = valid;
        const auto where =
            multi_group.find(
                "1 0 0 0 1 0 0 1 12 0");
        require(where != std::string::npos,
                "multi-group fixture marker");
        multi_group.replace(
            where,
            std::string{
                "1 0 0 0 1 0 0 1 12 0"}.size(),
            "1 0 0 0 1 0 0 2 12 99 0");
        expect_throw<std::invalid_argument>([&] {
            (void)mesh::import_gmsh_4_1_ascii(
                multi_group, 1.0);
        });
    }

    {
        auto nonplanar = valid;
        const auto where =
            nonplanar.find("3 1 0\n0 0 0");
        require(where != std::string::npos,
                "nonplanar fixture marker");
        nonplanar.replace(
            where,
            std::string{"3 1 0"}.size(),
            "3 1 0.25");
        expect_throw<std::invalid_argument>([&] {
            (void)mesh::import_gmsh_4_1_ascii(
                nonplanar, 1.0);
        });
    }
}

std::string vtu_mixed_triangle_quad_fixture() {
    return R"vtu(<?xml version="1.0"?>
<VTKFile type="UnstructuredGrid" version="1.0" byte_order="LittleEndian">
  <UnstructuredGrid>
    <Piece NumberOfPoints="5" NumberOfCells="2">
      <PointData>
        <DataArray type="Float32" Name="point.velocity" NumberOfComponents="2" format="ascii"
                   mpmc_unit="m/s"
                   mpmc_source_kind="synthetic_test"
                   mpmc_source_reference="vtu&amp;fixture"
                   mpmc_source_revision="r1"
                   mpmc_source_locator="PointData/velocity">
          1 2
          3 4
          5 6
          7 8
          9 10
        </DataArray>
      </PointData>
      <CellData>
        <DataArray type="UInt64" Name="mpmc_global_cell_id" format="ascii">
          9223372036854775813 7000000003
        </DataArray>
        <DataArray type="Float64" Name="cell.state" NumberOfComponents="2" format="ascii"
                   mpmc_unit="Pa"
                   mpmc_source_kind="user_supplied"
                   mpmc_source_reference="lab-run"
                   mpmc_source_revision="v2"
                   mpmc_source_locator="CellData/state">
          101.5 0.25
          202.5 0.75
        </DataArray>
      </CellData>
      <Points>
        <DataArray type="Float64" NumberOfComponents="3" format="ascii">
          0 0 0
          1 0 0
          1 1 0
          3 0 0
          3 1 0
        </DataArray>
      </Points>
      <Cells>
        <DataArray type="Int64" Name="connectivity" format="ascii">
          0 1 2
          1 3 4 2
        </DataArray>
        <DataArray type="Int64" Name="offsets" format="ascii">
          3 7
        </DataArray>
        <DataArray type="UInt8" Name="types" format="ascii">
          5 9
        </DataArray>
      </Cells>
    </Piece>
  </UnstructuredGrid>
</VTKFile>
)vtu";
}

void require_dense_field_equal(
    const mesh::DenseFieldSnapshot& actual,
    const mesh::DenseFieldSnapshot& expected) {
    require(
        actual.location() ==
            expected.location(),
        "VTU round-trip field location");
    require(
        actual.entity_count() ==
            expected.entity_count(),
        "VTU round-trip field entity count");
    require(
        actual.component_count() ==
            expected.component_count(),
        "VTU round-trip field component count");
    require(
        actual.values().size() ==
            expected.values().size(),
        "VTU round-trip field value count");
    for (std::size_t i = 0U;
         i < actual.values().size();
         ++i) {
        require_close(
            actual.values()[i],
            expected.values()[i],
            1.0e-14,
            "VTU round-trip field value");
    }

    const auto& a = actual.metadata();
    const auto& e = expected.metadata();
    require(
        a.id == e.id &&
            a.unit == e.unit &&
            a.source.kind ==
                e.source.kind &&
            a.source.reference ==
                e.source.reference &&
            a.source.revision ==
                e.source.revision &&
            a.source.locator ==
                e.source.locator,
        "VTU round-trip field metadata");
}

void require_vtu_fields_equal(
    const std::vector<mesh::DenseFieldSnapshot>& actual,
    const std::vector<mesh::DenseFieldSnapshot>& expected) {
    require(
        actual.size() == expected.size(),
        "VTU round-trip field array count");
    for (std::size_t i = 0U;
         i < actual.size();
         ++i) {
        require_dense_field_equal(
            actual[i], expected[i]);
    }
}

void vtu_ascii_import() {
    const auto imported =
        mesh::import_vtu_ascii(
            vtu_mixed_triangle_quad_fixture());

    const auto& topology = imported.topology;
    require(
        topology.entity_count(
            mesh::EntityKind::vertex) == 5U &&
        topology.entity_count(
            mesh::EntityKind::face) == 6U &&
        topology.entity_count(
            mesh::EntityKind::cell) == 2U &&
        topology.entity_count(
            mesh::EntityKind::edge) == 0U,
        "VTU imported entity counts");

    require(
        topology.global_id(
            mesh::EntityKind::cell,
            mesh::LocalIndex{0U}).value() ==
                9223372036854775813ULL &&
        topology.global_id(
            mesh::EntityKind::cell,
            mesh::LocalIndex{1U}).value() ==
                7000000003ULL,
        "VTU UInt64 stable cell IDs");

    const auto& cell_vertices =
        topology.relation(
            mesh::EntityKind::cell,
            mesh::EntityKind::vertex);
    require(
        cell_vertices.adjacent(
            mesh::LocalIndex{0U}).size() == 3U &&
        cell_vertices.adjacent(
            mesh::LocalIndex{1U}).size() == 4U,
        "VTU triangle/quad connectivity widths");

    const auto& face_cells =
        topology.relation(
            mesh::EntityKind::face,
            mesh::EntityKind::cell);
    std::size_t interior_faces = 0U;
    for (std::size_t face = 0U;
         face < topology.entity_count(
             mesh::EntityKind::face);
         ++face) {
        if (face_cells.adjacent(
                mesh::LocalIndex{
                    static_cast<
                        mesh::LocalIndex::value_type>(
                            face)})
                .size() == 2U) {
            ++interior_faces;
        }
    }
    require(
        interior_faces == 1U,
        "VTU reconstructed shared face");

    const auto& geometry = imported.geometry;
    require_close(
        geometry.cell_area_m2(
            mesh::LocalIndex{0U}),
        0.5,
        1.0e-14,
        "VTU triangle area");
    require_close(
        geometry.cell_area_m2(
            mesh::LocalIndex{1U}),
        2.0,
        1.0e-14,
        "VTU quad area");
    const auto triangle_centroid =
        geometry.cell_centroid_m(
            mesh::LocalIndex{0U});
    require_close(
        triangle_centroid.x_m,
        2.0 / 3.0,
        1.0e-14,
        "VTU triangle centroid x");
    require_close(
        triangle_centroid.y_m,
        1.0 / 3.0,
        1.0e-14,
        "VTU triangle centroid y");

    require(
        imported.point_fields.size() == 1U &&
        imported.cell_fields.size() == 1U,
        "VTU imported point/cell field counts");
    const auto& point_field =
        imported.point_fields.front();
    require(
        point_field.location() ==
            mesh::EntityKind::vertex &&
        point_field.component_count() == 2U &&
        point_field.entity_count() == 5U,
        "VTU point field layout");
    require_close(
        point_field.value(
            mesh::LocalIndex{4U}, 1U),
        10.0,
        0.0,
        "VTU point field value");
    require(
        point_field.metadata().unit == "m/s" &&
        point_field.metadata().source.kind ==
            mesh::FieldSourceKind::synthetic_test &&
        point_field.metadata().source.reference ==
            "vtu&fixture",
        "VTU point field metadata/unescape");

    const auto& cell_field =
        imported.cell_fields.front();
    require(
        cell_field.location() ==
            mesh::EntityKind::cell &&
        cell_field.component_count() == 2U &&
        cell_field.entity_count() == 2U,
        "VTU cell field layout");
    require_close(
        cell_field.value(
            mesh::LocalIndex{1U}, 1U),
        0.75,
        0.0,
        "VTU cell field value");
}

void vtu_ascii_roundtrip() {
    const auto first =
        mesh::import_vtu_ascii(
            vtu_mixed_triangle_quad_fixture());
    const std::string exported =
        mesh::export_vtu_ascii(first);

    require(
        exported.find(
            "<VTKFile type=\"UnstructuredGrid\"") !=
            std::string::npos,
        "VTU export root type");
    require(
        exported.find(
            "Name=\"mpmc_global_cell_id\"") !=
            std::string::npos,
        "VTU export stable cell identity array");
    require(
        exported.find(
            "9223372036854775813") !=
            std::string::npos,
        "VTU export full UInt64 cell ID");
    require(
        exported.find(
            "mpmc_source_reference=\"vtu&amp;fixture\"") !=
            std::string::npos,
        "VTU export XML-escaped provenance");

    const auto second =
        mesh::import_vtu_ascii(exported);

    require_topology_equal(
        second.topology,
        first.topology);
    require_geometry_equal(
        second.geometry,
        first.geometry);
    require_vtu_fields_equal(
        second.point_fields,
        first.point_fields);
    require_vtu_fields_equal(
        second.cell_fields,
        first.cell_fields);
}

void vtu_ascii_invalid() {
    const auto valid =
        vtu_mixed_triangle_quad_fixture();

    const auto expect_import_invalid =
        [](std::string_view label,
           const std::string& content) {
            bool caught = false;
            try {
                (void)mesh::import_vtu_ascii(content);
            } catch (const std::invalid_argument&) {
                caught = true;
            } catch (const std::out_of_range&) {
                caught = true;
            }
            require(
                caught,
                std::string{"expected invalid VTU rejection: "} +
                    std::string{label});
        };

    {
        auto appended = valid;
        const auto where =
            appended.find("format=\"ascii\"");
        require(
            where != std::string::npos,
            "VTU appended marker");
        appended.replace(
            where,
            std::string{"format=\"ascii\""}.size(),
            "format=\"appended\"");
        expect_import_invalid(
            "appended format", appended);
    }

    {
        auto compressed = valid;
        const auto where =
            compressed.find(
                "byte_order=\"LittleEndian\"");
        require(
            where != std::string::npos,
            "VTU compressor marker");
        compressed.replace(
            where,
            std::string{
                "byte_order=\"LittleEndian\""}.size(),
            "byte_order=\"LittleEndian\" compressor=\"vtkZLibDataCompressor\"");
        expect_import_invalid(
            "compressor attribute", compressed);
    }

    {
        auto unsupported = valid;
        const auto where =
            unsupported.find("5 9");
        require(
            where != std::string::npos,
            "VTU unsupported type marker");
        unsupported.replace(
            where,
            std::string{"5 9"}.size(),
            "5 10");
        expect_import_invalid(
            "unsupported cell type", unsupported);
    }

    {
        auto nonplanar = valid;
        const auto where =
            nonplanar.find("3 1 0");
        require(
            where != std::string::npos,
            "VTU nonplanar point marker");
        nonplanar.replace(
            where,
            std::string{"3 1 0"}.size(),
            "3 1 0.5");
        expect_import_invalid(
            "nonplanar coordinates", nonplanar);
    }

    {
        auto bad_offsets = valid;
        const auto offsets_name =
            bad_offsets.find(
                "Name=\"offsets\"");
        require(
            offsets_name != std::string::npos,
            "VTU offsets array marker");
        const auto where =
            bad_offsets.find(
                "3 7", offsets_name);
        require(
            where != std::string::npos,
            "VTU offsets values marker");
        bad_offsets.replace(
            where,
            std::string{"3 7"}.size(),
            "3 6");
        expect_import_invalid(
            "bad offsets", bad_offsets);
    }

    {
        auto duplicate_ids = valid;
        const auto where =
            duplicate_ids.find(
                "9223372036854775813 7000000003");
        require(
            where != std::string::npos,
            "VTU stable ID marker");
        duplicate_ids.replace(
            where,
            std::string{
                "9223372036854775813 7000000003"}.size(),
            "7000000003 7000000003");
        expect_import_invalid(
            "duplicate global cell IDs",
            duplicate_ids);
    }

    {
        auto imported =
            mesh::import_vtu_ascii(valid);
        imported.cell_fields.push_back(
            mesh::DenseFieldSnapshot::create(
                imported.topology,
                mesh::EntityKind::cell,
                1U,
                {1.0, 2.0},
                mesh::DenseFieldMetadata{
                    "mpmc_global_cell_id",
                    "1",
                    mesh::FieldSourceMetadata{
                        mesh::FieldSourceKind::synthetic_test,
                        "reserved-name-test",
                        "v1",
                        "CellData"}}));
        bool caught = false;
        try {
            (void)mesh::export_vtu_ascii(
                imported);
        } catch (const std::invalid_argument&) {
            caught = true;
        }
        require(
            caught,
            "expected invalid VTU rejection: reserved cell field name");
    }
}

std::string grdecl_two_cell_fixture() {
    return R"grdecl(-- minimal 2x1x1 Cartesian corner-point deck
SPECGRID
  2 1 1 1 F /
COORD
  0 0 0   0 0 1
  1 0 0   1 0 1
  2 0 0   2 0 1
  0 1 0   0 1 1
  1 1 0   1 1 1
  2 1 0   2 1 1 /
ZCORN
  8*0 8*1 /
ACTNUM
  1 0 /
PORO
  0.20 0.35 /
PERMX
  100 200 /
PERMY
  50 75 /
PERMZ
  1.0D+1 2.0D+1 /
)grdecl";
}

std::string grdecl_single_cell_fixture(
    std::string_view zcorn,
    int actnum) {
    std::ostringstream output;
    output
        << "SPECGRID\n"
        << " 1 1 1 1 F /\n"
        << "COORD\n"
        << " 0 0 0  0 0 1\n"
        << " 1 0 0  1 0 1\n"
        << " 0 1 0  0 1 1\n"
        << " 1 1 0  1 1 1 /\n"
        << "ZCORN\n "
        << zcorn << " /\n"
        << "ACTNUM\n "
        << actnum << " /\n"
        << "PORO\n 0.25 /\n"
        << "PERMX\n 100 /\n"
        << "PERMY\n 100 /\n"
        << "PERMZ\n 10 /\n";
    return output.str();
}

const mesh::DenseFieldSnapshot&
grdecl_field(const mesh::GrdeclImportResult& imported,
             std::string_view id) {
    const auto found =
        std::find_if(
            imported.cell_fields.begin(),
            imported.cell_fields.end(),
            [id](const auto& field) {
                return field.metadata().id == id;
            });
    require(
        found != imported.cell_fields.end(),
        "GRDECL field ID missing");
    return *found;
}

void grdecl_import() {
    constexpr double permeability_scale =
        1.0e-15;
    const auto imported =
        mesh::import_grdecl(
            grdecl_two_cell_fixture(),
            mesh::GrdeclImportOptions{
                2.0,
                permeability_scale});

    require(
        imported.dimensions ==
            std::array<std::size_t, 3>{
                2U, 1U, 1U},
        "GRDECL SPECGRID dimensions");
    require(
        imported.cell_count() == 2U &&
        imported.active_cell_count() == 1U,
        "GRDECL ACTNUM counts");
    require(
        imported.is_active(
            mesh::LocalIndex{0U}) &&
        !imported.is_active(
            mesh::LocalIndex{1U}),
        "GRDECL active/inactive alignment");

    const auto& topology =
        imported.topology;
    require(
        topology.entity_count(
            mesh::EntityKind::cell) == 2U &&
        topology.entity_count(
            mesh::EntityKind::vertex) == 16U &&
        topology.entity_count(
            mesh::EntityKind::face) == 0U &&
        topology.entity_count(
            mesh::EntityKind::edge) == 0U,
        "GRDECL minimal topology entity counts");
    require(
        topology.relation_count() == 1U &&
        topology.has_relation(
            mesh::EntityKind::cell,
            mesh::EntityKind::vertex),
        "GRDECL baseline only materializes cell-to-corner relation");
    require(
        topology.global_id(
            mesh::EntityKind::cell,
            mesh::LocalIndex{0U}).value() == 1U &&
        topology.global_id(
            mesh::EntityKind::cell,
            mesh::LocalIndex{1U}).value() == 2U,
        "GRDECL logical cell stable IDs");
    require(
        topology.global_id(
            mesh::EntityKind::vertex,
            mesh::LocalIndex{15U}).value() == 16U,
        "GRDECL cell-local corner stable IDs");

    const auto& cell_vertices =
        topology.relation(
            mesh::EntityKind::cell,
            mesh::EntityKind::vertex);
    for (std::size_t cell = 0U;
         cell < 2U;
         ++cell) {
        const auto row =
            cell_vertices.adjacent(
                mesh::LocalIndex{
                    static_cast<
                        mesh::LocalIndex::value_type>(
                            cell)});
        require(
            row.size() == 8U,
            "GRDECL cell must retain eight local corners");
        for (std::size_t corner = 0U;
             corner < 8U;
             ++corner) {
            require(
                row[corner].value() ==
                    static_cast<
                        mesh::LocalIndex::value_type>(
                            cell * 8U + corner),
                "GRDECL cell-local corner ordering");
        }
    }

    const auto& geometry =
        imported.geometry;
    require(
        geometry.cell_count() == 2U &&
        geometry.vertex_count() == 16U,
        "GRDECL corner geometry counts");
    require_close(
        geometry.cell_volume_m3(
            mesh::LocalIndex{0U}),
        8.0,
        1.0e-13,
        "GRDECL active cell volume");
    require_close(
        geometry.cell_volume_m3(
            mesh::LocalIndex{1U}),
        8.0,
        1.0e-13,
        "GRDECL inactive nondegenerate cell volume");

    const auto first =
        geometry.cell_corners_m(
            mesh::LocalIndex{0U});
    require_close(
        first[0].x_m, 0.0, 0.0,
        "GRDECL first corner x");
    require_close(
        first[0].y_m, 0.0, 0.0,
        "GRDECL first corner y");
    require_close(
        first[0].z_m, 0.0, 0.0,
        "GRDECL first corner z");
    require_close(
        first[7].x_m, 2.0, 0.0,
        "GRDECL last corner x");
    require_close(
        first[7].y_m, 2.0, 0.0,
        "GRDECL last corner y");
    require_close(
        first[7].z_m, 2.0, 0.0,
        "GRDECL last corner z");

    require(
        imported.cell_fields.size() == 4U,
        "GRDECL PORO/PERM field count");
    const auto& poro =
        grdecl_field(imported, "PORO");
    const auto& permx =
        grdecl_field(imported, "PERMX");
    const auto& permy =
        grdecl_field(imported, "PERMY");
    const auto& permz =
        grdecl_field(imported, "PERMZ");

    require(
        poro.location() ==
            mesh::EntityKind::cell &&
        poro.component_count() == 1U &&
        poro.metadata().unit == "1",
        "GRDECL PORO field contract");
    require_close(
        poro.value(
            mesh::LocalIndex{1U}, 0U),
        0.35,
        0.0,
        "GRDECL inactive PORO retained");

    for (const auto* field :
         {&permx, &permy, &permz}) {
        require(
            field->metadata().unit == "m2" &&
            field->metadata().source.kind ==
                mesh::FieldSourceKind::file_import &&
            field->entity_count() == 2U,
            "GRDECL permeability metadata/alignment");
    }
    require_close(
        permx.value(
            mesh::LocalIndex{0U}, 0U),
        100.0 * permeability_scale,
        1.0e-28,
        "GRDECL PERMX SI conversion");
    require_close(
        permy.value(
            mesh::LocalIndex{1U}, 0U),
        75.0 * permeability_scale,
        1.0e-28,
        "GRDECL inactive PERMY retained");
    require_close(
        permz.value(
            mesh::LocalIndex{1U}, 0U),
        20.0 * permeability_scale,
        1.0e-28,
        "GRDECL Fortran-D exponent parsing");
}

void grdecl_inactive_degenerate() {
    const auto inactive =
        mesh::import_grdecl(
            grdecl_single_cell_fixture(
                "8*0", 0),
            mesh::GrdeclImportOptions{
                1.0, 1.0e-15});
    require(
        inactive.cell_count() == 1U &&
        !inactive.is_active(
            mesh::LocalIndex{0U}),
        "GRDECL degenerate inactive identity");
    require(
        inactive.topology.global_id(
            mesh::EntityKind::cell,
            mesh::LocalIndex{0U}).value() == 1U,
        "GRDECL degenerate inactive stable cell ID");
    require_close(
        inactive.geometry.cell_volume_m3(
            mesh::LocalIndex{0U}),
        0.0,
        0.0,
        "GRDECL degenerate inactive zero volume");
    require_close(
        grdecl_field(
            inactive, "PORO")
            .value(
                mesh::LocalIndex{0U}, 0U),
        0.25,
        0.0,
        "GRDECL degenerate inactive property retained");
}

void grdecl_invalid() {
    const auto options =
        mesh::GrdeclImportOptions{
            1.0, 1.0e-15};

    {
        auto radial =
            grdecl_two_cell_fixture();
        const auto where =
            radial.find("2 1 1 1 F");
        require(
            where != std::string::npos,
            "GRDECL radial marker");
        radial.replace(
            where,
            std::string{"2 1 1 1 F"}.size(),
            "2 1 1 1 T");
        expect_throw<std::invalid_argument>(
            [&] {
                (void)mesh::import_grdecl(
                    radial, options);
            });
    }

    {
        auto bad_actnum =
            grdecl_two_cell_fixture();
        const auto where =
            bad_actnum.find(
                "ACTNUM\n  1 0");
        require(
            where != std::string::npos,
            "GRDECL ACTNUM marker");
        bad_actnum.replace(
            where,
            std::string{
                "ACTNUM\n  1 0"}.size(),
            "ACTNUM\n  1 2");
        expect_throw<std::invalid_argument>(
            [&] {
                (void)mesh::import_grdecl(
                    bad_actnum, options);
            });
    }

    {
        auto missing =
            grdecl_two_cell_fixture();
        const auto begin =
            missing.find("PERMZ");
        require(
            begin != std::string::npos,
            "GRDECL missing-keyword marker");
        missing.erase(begin);
        expect_throw<std::invalid_argument>(
            [&] {
                (void)mesh::import_grdecl(
                    missing, options);
            });
    }

    {
        auto bad_poro =
            grdecl_two_cell_fixture();
        const auto where =
            bad_poro.find(
                "0.20 0.35");
        require(
            where != std::string::npos,
            "GRDECL PORO marker");
        bad_poro.replace(
            where,
            std::string{"0.20 0.35"}.size(),
            "0.20 1.35");
        expect_throw<std::invalid_argument>(
            [&] {
                (void)mesh::import_grdecl(
                    bad_poro, options);
            });
    }

    {
        auto bad_perm =
            grdecl_two_cell_fixture();
        const auto where =
            bad_perm.find(
                "100 200");
        require(
            where != std::string::npos,
            "GRDECL PERMX marker");
        bad_perm.replace(
            where,
            std::string{"100 200"}.size(),
            "100 -1");
        expect_throw<std::invalid_argument>(
            [&] {
                (void)mesh::import_grdecl(
                    bad_perm, options);
            });
    }

    expect_throw<std::invalid_argument>(
        [&] {
            (void)mesh::import_grdecl(
                grdecl_single_cell_fixture(
                    "8*0", 1),
                options);
        });

    expect_throw<std::invalid_argument>(
        [&] {
            (void)mesh::import_grdecl(
                grdecl_single_cell_fixture(
                    "4*1 4*0", 0),
                options);
        });

    expect_throw<std::invalid_argument>(
        [&] {
            (void)mesh::import_grdecl(
                grdecl_two_cell_fixture() +
                    "\nINCLUDE\n 'other.inc' /\n",
                options);
        });

    expect_throw<std::invalid_argument>(
        [&] {
            (void)mesh::import_grdecl(
                grdecl_two_cell_fixture(),
                mesh::GrdeclImportOptions{
                    0.0, 1.0e-15});
        });
    expect_throw<std::invalid_argument>(
        [&] {
            (void)mesh::import_grdecl(
                grdecl_two_cell_fixture(),
                mesh::GrdeclImportOptions{
                    1.0, 0.0});
        });
}

std::string grdecl_two_cell_all_active_fixture() {
    auto content =
        grdecl_two_cell_fixture();
    const auto where =
        content.find("ACTNUM\n  1 0");
    require(
        where != std::string::npos,
        "active face-processor ACTNUM marker");
    content.replace(
        where,
        std::string{"ACTNUM\n  1 0"}.size(),
        "ACTNUM\n  1 1");
    return content;
}

std::string grdecl_second_cell_only_fixture() {
    auto content =
        grdecl_two_cell_fixture();
    const auto where =
        content.find("ACTNUM\n  1 0");
    require(
        where != std::string::npos,
        "filtered face-processor ACTNUM marker");
    content.replace(
        where,
        std::string{"ACTNUM\n  1 0"}.size(),
        "ACTNUM\n  0 1");
    return content;
}

std::string grdecl_vertical_two_cell_fixture() {
    return R"grdecl(SPECGRID
  1 1 2 1 F /
COORD
  0 0 0   0 0 2
  1 0 0   1 0 2
  0 1 0   0 1 2
  1 1 0   1 1 2 /
ZCORN
  4*0 8*1 4*2 /
ACTNUM
  2*1 /
PORO
  0.20 0.30 /
PERMX
  100 110 /
PERMY
  90 95 /
PERMZ
  10 12 /
)grdecl";
}

std::string grdecl_skewed_two_cell_fixture() {
    return R"grdecl(-- valid skewed 2x1x1 corner-point deck
SPECGRID
  2 1 1 1 F /
COORD
  0.0 0 0   0.0 0 1
  1.2 0 0   1.2 0 1
  2.0 0 0   2.0 0 1
  0.0 1 0   0.0 1 1
  0.8 1 0   0.8 1 1
  2.0 1 0   2.0 1 1 /
ZCORN
  8*0 8*1 /
ACTNUM
  2*1 /
PORO
  0.20 0.35 /
PERMX
  100 200 /
PERMY
  50 75 /
PERMZ
  10 20 /
)grdecl";
}

std::string grdecl_fault_split_two_cell_fixture() {
    return R"grdecl(SPECGRID
  2 1 1 1 F /
COORD
  0 0 0   0 0 2
  1 0 0   1 0 2
  2 0 0   2 0 2
  0 1 0   0 1 2
  1 1 0   1 1 2
  2 1 0   2 1 2 /
ZCORN
  0 0 0.1 0
  0 0 0.1 0
  1 1 1.1 1
  1 1 1.1 1 /
ACTNUM
  2*1 /
PORO
  0.20 0.25 /
PERMX
  100 120 /
PERMY
  80 90 /
PERMZ
  10 11 /
)grdecl";
}

const mesh::DenseFieldSnapshot&
active_corner_point_field(
    const mesh::ActiveCornerPointGrid& processed,
    std::string_view id) {
    const auto found =
        std::find_if(
            processed.cell_fields.begin(),
            processed.cell_fields.end(),
            [id](const auto& field) {
                return field.metadata().id == id;
            });
    require(
        found != processed.cell_fields.end(),
        "processed active field ID missing");
    return *found;
}

void require_metadata_equal(
    const mesh::DenseFieldMetadata& actual,
    const mesh::DenseFieldMetadata& expected,
    std::string_view message) {
    require(
        actual.id == expected.id &&
        actual.unit == expected.unit &&
        actual.source.kind ==
            expected.source.kind &&
        actual.source.reference ==
            expected.source.reference &&
        actual.source.revision ==
            expected.source.revision &&
        actual.source.locator ==
            expected.source.locator,
        std::string{message});
}

mesh::Coordinate3D processed_cell_vertex_mean(
    const mesh::ActiveCornerPointGrid& processed,
    mesh::LocalIndex cell) {
    const auto& cell_vertices =
        processed.topology.relation(
            mesh::EntityKind::cell,
            mesh::EntityKind::vertex);
    const auto vertices =
        cell_vertices.adjacent(cell);
    require(
        vertices.size() == 8U,
        "processed cell mean requires eight vertices");

    mesh::Coordinate3D mean{0.0, 0.0, 0.0};
    for (const auto vertex : vertices) {
        const auto& coordinate =
            processed.vertex_coordinates_m[
                static_cast<std::size_t>(
                    vertex.value())];
        mean.x_m += coordinate.x_m;
        mean.y_m += coordinate.y_m;
        mean.z_m += coordinate.z_m;
    }
    mean.x_m /= 8.0;
    mean.y_m /= 8.0;
    mean.z_m /= 8.0;
    return mean;
}

std::size_t shared_face_count(
    const mesh::Topology& topology) {
    const auto& face_cells =
        topology.relation(
            mesh::EntityKind::face,
            mesh::EntityKind::cell);
    std::size_t shared = 0U;
    for (std::size_t face = 0U;
         face < topology.entity_count(
             mesh::EntityKind::face);
         ++face) {
        if (face_cells.adjacent(
                mesh::LocalIndex{
                    static_cast<
                        mesh::LocalIndex::value_type>(
                            face)})
                .size() == 2U) {
            ++shared;
        }
    }
    return shared;
}

void require_active_processor_shape(
    const mesh::ActiveCornerPointGrid& processed) {
    const auto& topology =
        processed.topology;
    require(
        topology.entity_count(
            mesh::EntityKind::cell) == 2U &&
        topology.entity_count(
            mesh::EntityKind::vertex) == 12U &&
        topology.entity_count(
            mesh::EntityKind::face) == 11U &&
        topology.entity_count(
            mesh::EntityKind::edge) == 0U,
        "active corner-point merged entity counts");
    require(
        topology.relation_count() == 4U,
        "active corner-point four core relations");
    require(
        shared_face_count(topology) == 1U,
        "active corner-point exactly one shared face");
    require(
        processed.source_logical_cell_ids.size() == 2U &&
        processed.source_logical_cell_ids[0].value() == 1U &&
        processed.source_logical_cell_ids[1].value() == 2U,
        "active corner-point source logical mapping");
    require(
        topology.global_id(
            mesh::EntityKind::cell,
            mesh::LocalIndex{0U}).value() == 1U &&
        topology.global_id(
            mesh::EntityKind::cell,
            mesh::LocalIndex{1U}).value() == 2U,
        "active corner-point processed cell identities");

    require(
        processed.face_geometry.cell_count() == 2U &&
        processed.face_geometry.face_count() == 11U,
        "active corner-point face metric counts");
    require(
        processed.cell_fields.size() == 4U,
        "active corner-point projected field count");
    for (const auto& field :
         processed.cell_fields) {
        require(
            field.location() ==
                mesh::EntityKind::cell &&
            field.entity_count() == 2U,
            "active corner-point projected field alignment");
    }

    const auto& cell_vertices =
        topology.relation(
            mesh::EntityKind::cell,
            mesh::EntityKind::vertex);
    const auto& cell_faces =
        topology.relation(
            mesh::EntityKind::cell,
            mesh::EntityKind::face);
    const auto& face_vertices =
        topology.relation(
            mesh::EntityKind::face,
            mesh::EntityKind::vertex);
    const auto& face_cells =
        topology.relation(
            mesh::EntityKind::face,
            mesh::EntityKind::cell);
    for (std::size_t cell = 0U;
         cell < 2U;
         ++cell) {
        const auto local =
            mesh::LocalIndex{
                static_cast<
                    mesh::LocalIndex::value_type>(
                        cell)};
        require(
            cell_vertices.adjacent(local).size() == 8U,
            "processed hexa cell has eight shared vertices");
        require(
            cell_faces.adjacent(local).size() == 6U,
            "processed hexa cell has six faces");
    }
    for (std::size_t face = 0U;
         face < topology.entity_count(
             mesh::EntityKind::face);
         ++face) {
        const auto local =
            mesh::LocalIndex{
                static_cast<
                    mesh::LocalIndex::value_type>(
                        face)};
        require(
            face_vertices.adjacent(local).size() == 4U,
            "processed quad face has four vertices");
        const auto support =
            face_cells.adjacent(local).size();
        require(
            support == 1U || support == 2U,
            "processed face support cardinality");

        require_close(
            processed.face_geometry.face_area_m2(local),
            1.0,
            1.0e-14,
            "processed unit-cube face area");

        const auto owner =
            processed.face_geometry.face_owner(local);
        require(
            face_cells.adjacent(local).front() == owner,
            "face metric owner follows canonical face support owner");
        const auto centroid =
            processed.face_geometry.face_centroid_m(local);
        const auto owner_centroid =
            processed_cell_vertex_mean(
                processed, owner);
        const auto normal =
            processed.face_geometry.face_owner_unit_normal(
                local);
        const double outward =
            normal.x *
                (centroid.x_m -
                 owner_centroid.x_m) +
            normal.y *
                (centroid.y_m -
                 owner_centroid.y_m) +
            normal.z *
                (centroid.z_m -
                 owner_centroid.z_m);
        require(
            outward > 0.0,
            "processed face normal points away from owner cell");
        require_close(
            std::sqrt(
                normal.x * normal.x +
                normal.y * normal.y +
                normal.z * normal.z),
            1.0,
            1.0e-14,
            "processed face normal magnitude");
    }
}

void active_corner_point_i_neighbor() {
    const auto raw =
        mesh::import_grdecl(
            grdecl_two_cell_all_active_fixture(),
            mesh::GrdeclImportOptions{
                1.0, 1.0e-15});
    const auto processed =
        mesh::process_active_corner_point_grid(
            raw);
    require_active_processor_shape(
        processed);
    require_close(
        processed.cell_volumes_m3[0],
        1.0,
        1.0e-14,
        "processed I-neighbor first volume");
    require_close(
        processed.cell_volumes_m3[1],
        1.0,
        1.0e-14,
        "processed I-neighbor second volume");

    const auto& face_cells =
        processed.topology.relation(
            mesh::EntityKind::face,
            mesh::EntityKind::cell);
    const auto& face_vertices =
        processed.topology.relation(
            mesh::EntityKind::face,
            mesh::EntityKind::vertex);
    bool found_interface = false;
    for (std::size_t face = 0U;
         face < processed.topology.entity_count(
             mesh::EntityKind::face);
         ++face) {
        const auto local =
            mesh::LocalIndex{
                static_cast<
                    mesh::LocalIndex::value_type>(
                        face)};
        if (face_cells.adjacent(local).size() != 2U) {
            continue;
        }
        found_interface = true;
        const auto centroid =
            processed.face_geometry.face_centroid_m(local);
        const auto normal =
            processed.face_geometry.face_owner_unit_normal(
                local);
        require_close(
            centroid.x_m, 1.0, 0.0,
            "processed I-interface centroid x");
        require_close(
            centroid.y_m, 0.5, 1.0e-14,
            "processed I-interface centroid y");
        require_close(
            centroid.z_m, 0.5, 1.0e-14,
            "processed I-interface centroid z");
        require_close(
            processed.face_geometry.face_area_m2(local),
            1.0,
            1.0e-14,
            "processed I-interface area");
        require(
            processed.face_geometry.face_owner(local) ==
                mesh::LocalIndex{0U},
            "processed I-interface owner");
        require_close(
            normal.x, 1.0, 1.0e-14,
            "processed I-interface normal x");
        require_close(
            normal.y, 0.0, 1.0e-14,
            "processed I-interface normal y");
        require_close(
            normal.z, 0.0, 1.0e-14,
            "processed I-interface normal z");
        for (const auto vertex :
             face_vertices.adjacent(local)) {
            require_close(
                processed.vertex_coordinates_m[
                    static_cast<std::size_t>(
                        vertex.value())]
                    .x_m,
                1.0,
                0.0,
                "processed I-interface x coordinate");
        }
    }
    require(
        found_interface,
        "processed I-interface face exists");

    const auto& poro =
        active_corner_point_field(
            processed, "PORO");
    const auto& permx =
        active_corner_point_field(
            processed, "PERMX");
    require_close(
        poro.value(
            mesh::LocalIndex{0U}, 0U),
        0.20,
        0.0,
        "processed first PORO");
    require_close(
        poro.value(
            mesh::LocalIndex{1U}, 0U),
        0.35,
        0.0,
        "processed second PORO");
    require_close(
        permx.value(
            mesh::LocalIndex{0U}, 0U),
        100.0e-15,
        1.0e-28,
        "processed first PERMX");
    require_close(
        permx.value(
            mesh::LocalIndex{1U}, 0U),
        200.0e-15,
        1.0e-28,
        "processed second PERMX");
    require_metadata_equal(
        poro.metadata(),
        grdecl_field(raw, "PORO").metadata(),
        "processed PORO metadata preserved");
    require_metadata_equal(
        permx.metadata(),
        grdecl_field(raw, "PERMX").metadata(),
        "processed PERMX metadata preserved");

    const auto repeated =
        mesh::process_active_corner_point_grid(
            raw);
    require_topology_equal(
        repeated.topology,
        processed.topology);
    require(
        repeated.vertex_coordinates_m.size() ==
            processed.vertex_coordinates_m.size(),
        "active corner-point deterministic vertex count");
    for (std::size_t vertex = 0U;
         vertex < processed.vertex_coordinates_m.size();
         ++vertex) {
        const auto a =
            repeated.vertex_coordinates_m[vertex];
        const auto e =
            processed.vertex_coordinates_m[vertex];
        require(
            a.x_m == e.x_m &&
            a.y_m == e.y_m &&
            a.z_m == e.z_m,
            "active corner-point deterministic vertex ordering");
    }
}

void active_corner_point_k_neighbor() {
    const auto raw =
        mesh::import_grdecl(
            grdecl_vertical_two_cell_fixture(),
            mesh::GrdeclImportOptions{
                1.0, 1.0e-15});
    const auto processed =
        mesh::process_active_corner_point_grid(
            raw);
    require_active_processor_shape(
        processed);

    const auto& face_cells =
        processed.topology.relation(
            mesh::EntityKind::face,
            mesh::EntityKind::cell);
    const auto& face_vertices =
        processed.topology.relation(
            mesh::EntityKind::face,
            mesh::EntityKind::vertex);
    bool found_interface = false;
    for (std::size_t face = 0U;
         face < processed.topology.entity_count(
             mesh::EntityKind::face);
         ++face) {
        const auto local =
            mesh::LocalIndex{
                static_cast<
                    mesh::LocalIndex::value_type>(
                        face)};
        if (face_cells.adjacent(local).size() != 2U) {
            continue;
        }
        found_interface = true;
        const auto centroid =
            processed.face_geometry.face_centroid_m(local);
        const auto normal =
            processed.face_geometry.face_owner_unit_normal(
                local);
        require_close(
            centroid.x_m, 0.5, 1.0e-14,
            "processed K-interface centroid x");
        require_close(
            centroid.y_m, 0.5, 1.0e-14,
            "processed K-interface centroid y");
        require_close(
            centroid.z_m, 1.0, 0.0,
            "processed K-interface centroid z");
        require_close(
            processed.face_geometry.face_area_m2(local),
            1.0,
            1.0e-14,
            "processed K-interface area");
        require(
            processed.face_geometry.face_owner(local) ==
                mesh::LocalIndex{0U},
            "processed K-interface owner");
        require_close(
            normal.x, 0.0, 1.0e-14,
            "processed K-interface normal x");
        require_close(
            normal.y, 0.0, 1.0e-14,
            "processed K-interface normal y");
        require_close(
            normal.z, 1.0, 1.0e-14,
            "processed K-interface normal z");
        for (const auto vertex :
             face_vertices.adjacent(local)) {
            require_close(
                processed.vertex_coordinates_m[
                    static_cast<std::size_t>(
                        vertex.value())]
                    .z_m,
                1.0,
                0.0,
                "processed K-interface z coordinate");
        }
    }
    require(
        found_interface,
        "processed K-interface face exists");
}

void active_corner_point_activity_mapping() {
    const auto raw =
        mesh::import_grdecl(
            grdecl_second_cell_only_fixture(),
            mesh::GrdeclImportOptions{
                1.0, 1.0e-15});
    const auto processed =
        mesh::process_active_corner_point_grid(
            raw);

    require(
        processed.cell_count() == 1U &&
        processed.topology.entity_count(
            mesh::EntityKind::cell) == 1U &&
        processed.topology.entity_count(
            mesh::EntityKind::vertex) == 8U &&
        processed.topology.entity_count(
            mesh::EntityKind::face) == 6U,
        "active processor filters inactive logical cell");
    require(
        processed.source_logical_cell_id(
            mesh::LocalIndex{0U}).value() == 2U &&
        processed.topology.global_id(
            mesh::EntityKind::cell,
            mesh::LocalIndex{0U}).value() == 2U,
        "processed local cell maps to original logical GlobalEntityId");

    require(
        processed.face_geometry.cell_count() == 1U &&
        processed.face_geometry.face_count() == 6U &&
        processed.cell_fields.size() == 4U,
        "filtered active geometry/field counts");
    require_close(
        active_corner_point_field(
            processed, "PORO")
            .value(
                mesh::LocalIndex{0U}, 0U),
        0.35,
        0.0,
        "filtered logical cell 2 PORO");
    require_close(
        active_corner_point_field(
            processed, "PERMX")
            .value(
                mesh::LocalIndex{0U}, 0U),
        200.0e-15,
        1.0e-28,
        "filtered logical cell 2 PERMX");
    require_close(
        active_corner_point_field(
            processed, "PERMY")
            .value(
                mesh::LocalIndex{0U}, 0U),
        75.0e-15,
        1.0e-28,
        "filtered logical cell 2 PERMY");
    require_close(
        active_corner_point_field(
            processed, "PERMZ")
            .value(
                mesh::LocalIndex{0U}, 0U),
        20.0e-15,
        1.0e-28,
        "filtered logical cell 2 PERMZ");

    const auto& face_cells =
        processed.topology.relation(
            mesh::EntityKind::face,
            mesh::EntityKind::cell);
    for (std::size_t face = 0U;
         face < 6U;
         ++face) {
        require(
            face_cells.adjacent(
                mesh::LocalIndex{
                    static_cast<
                        mesh::LocalIndex::value_type>(
                            face)})
                .size() == 1U,
            "active-inactive interface becomes processed boundary");
    }
}

void active_corner_point_invalid() {
    const auto faulted_raw =
        mesh::import_grdecl(
            grdecl_fault_split_two_cell_fixture(),
            mesh::GrdeclImportOptions{
                1.0, 1.0e-15});
    expect_throw<std::invalid_argument>(
        [&] {
            (void)mesh::process_active_corner_point_grid(
                faulted_raw);
        });

    const auto no_active =
        mesh::import_grdecl(
            grdecl_single_cell_fixture(
                "4*0 4*1", 0),
            mesh::GrdeclImportOptions{
                1.0, 1.0e-15});
    expect_throw<std::invalid_argument>(
        [&] {
            (void)mesh::process_active_corner_point_grid(
                no_active);
        });
}


mesh::LocalIndex only_shared_face(
    const mesh::Topology& topology) {
    const auto& face_cells =
        topology.relation(
            mesh::EntityKind::face,
            mesh::EntityKind::cell);
    std::optional<mesh::LocalIndex> found;
    for (std::size_t face = 0U;
         face < topology.entity_count(
             mesh::EntityKind::face);
         ++face) {
        const auto local =
            mesh::LocalIndex{
                static_cast<
                    mesh::LocalIndex::value_type>(
                        face)};
        if (face_cells.adjacent(local).size() != 2U) {
            continue;
        }
        require(
            !found.has_value(),
            "expected exactly one shared face");
        found = local;
    }
    require(
        found.has_value(),
        "expected one shared face");
    return *found;
}

void require_shared_operator_axis(
    const mesh::ActiveCornerPointGrid& processed,
    std::array<double, 3> expected_owner_displacement,
    std::array<double, 3> expected_neighbour_displacement,
    std::array<double, 3> expected_normal) {
    const auto op =
        mesh::make_cell_face_geometric_operator_3d(
            processed.topology,
            processed.vertex_coordinates_m,
            processed.face_geometry);

    require(
        op.cell_count() == 2U &&
            op.face_count() == 11U,
        "3D cell-face geometric operator counts");
    require_close(
        op.cell_centroid_m(
            mesh::LocalIndex{0U}).x_m,
        processed_cell_vertex_mean(
            processed,
            mesh::LocalIndex{0U}).x_m,
        0.0,
        "3D operator owner cell centroid x");
    require_close(
        op.cell_centroid_m(
            mesh::LocalIndex{1U}).z_m,
        processed_cell_vertex_mean(
            processed,
            mesh::LocalIndex{1U}).z_m,
        0.0,
        "3D operator neighbour cell centroid z");

    const auto interface =
        only_shared_face(
            processed.topology);
    require(
        op.face_owner(interface) ==
            mesh::LocalIndex{0U},
        "3D operator shared face owner");
    require(
        op.face_neighbour(interface).has_value() &&
            *op.face_neighbour(interface) ==
                mesh::LocalIndex{1U},
        "3D operator shared face neighbour");

    const auto owner_displacement =
        op.owner_to_face_displacement_m(
            interface);
    const auto neighbour_displacement =
        op.neighbour_to_face_displacement_m(
            interface);
    require(
        neighbour_displacement.has_value(),
        "3D internal face neighbour displacement");
    require_close(
        owner_displacement.x_m,
        expected_owner_displacement[0],
        1.0e-14,
        "3D owner displacement x");
    require_close(
        owner_displacement.y_m,
        expected_owner_displacement[1],
        1.0e-14,
        "3D owner displacement y");
    require_close(
        owner_displacement.z_m,
        expected_owner_displacement[2],
        1.0e-14,
        "3D owner displacement z");
    require_close(
        neighbour_displacement->x_m,
        expected_neighbour_displacement[0],
        1.0e-14,
        "3D neighbour displacement x");
    require_close(
        neighbour_displacement->y_m,
        expected_neighbour_displacement[1],
        1.0e-14,
        "3D neighbour displacement y");
    require_close(
        neighbour_displacement->z_m,
        expected_neighbour_displacement[2],
        1.0e-14,
        "3D neighbour displacement z");

    require_close(
        op.owner_normal_distance_m(interface),
        0.5,
        1.0e-14,
        "3D owner normal distance");
    require(
        op.neighbour_normal_distance_m(
            interface).has_value(),
        "3D internal face neighbour distance");
    require_close(
        *op.neighbour_normal_distance_m(
            interface),
        0.5,
        1.0e-14,
        "3D neighbour normal distance");

    const auto inputs =
        op.transmissibility_geometry(
            interface);
    require_close(
        inputs.area_m2,
        1.0,
        1.0e-14,
        "3D transmissibility geometry area");
    require_close(
        inputs.owner_unit_normal.x,
        expected_normal[0],
        1.0e-14,
        "3D transmissibility geometry normal x");
    require_close(
        inputs.owner_unit_normal.y,
        expected_normal[1],
        1.0e-14,
        "3D transmissibility geometry normal y");
    require_close(
        inputs.owner_unit_normal.z,
        expected_normal[2],
        1.0e-14,
        "3D transmissibility geometry normal z");
    require_close(
        inputs.owner_normal_distance_m,
        0.5,
        1.0e-14,
        "3D transmissibility owner distance");
    require(
        inputs.neighbour_normal_distance_m.has_value(),
        "3D transmissibility neighbour distance");

    const auto non_orthogonality =
        op.internal_non_orthogonality(
            interface);
    require(
        non_orthogonality.has_value() &&
            inputs.internal_non_orthogonality.has_value(),
        "3D internal face exposes non-orthogonality geometry");
    require_close(
        non_orthogonality->center_distance_m,
        1.0,
        1.0e-14,
        "3D orthogonal center distance");
    require_close(
        non_orthogonality->normal_alignment_cosine,
        1.0,
        1.0e-14,
        "3D orthogonal normal alignment cosine");
    require_close(
        non_orthogonality->angle_rad,
        0.0,
        1.0e-14,
        "3D orthogonal non-orthogonality angle");

    const auto& face_cells =
        processed.topology.relation(
            mesh::EntityKind::face,
            mesh::EntityKind::cell);
    bool found_boundary = false;
    for (std::size_t face = 0U;
         face < processed.topology.entity_count(
             mesh::EntityKind::face);
         ++face) {
        const auto local =
            mesh::LocalIndex{
                static_cast<
                    mesh::LocalIndex::value_type>(
                        face)};
        if (face_cells.adjacent(local).size() != 1U) {
            continue;
        }
        found_boundary = true;
        require(
            !op.face_neighbour(local).has_value() &&
                !op.neighbour_to_face_displacement_m(
                    local).has_value() &&
                !op.neighbour_normal_distance_m(
                    local).has_value() &&
                !op.internal_non_orthogonality(
                    local).has_value(),
            "3D boundary face has no neighbour/non-orthogonality geometry");
        require(
            op.owner_normal_distance_m(local) > 0.0,
            "3D boundary owner normal distance positive");
        break;
    }
    require(
        found_boundary,
        "3D operator boundary face exists");
}

void cell_face_geometric_operator_3d() {
    const auto horizontal =
        mesh::process_active_corner_point_grid(
            mesh::import_grdecl(
                grdecl_two_cell_all_active_fixture(),
                mesh::GrdeclImportOptions{
                    1.0, 1.0e-15}));
    require_shared_operator_axis(
        horizontal,
        {0.5, 0.0, 0.0},
        {-0.5, 0.0, 0.0},
        {1.0, 0.0, 0.0});

    const auto vertical =
        mesh::process_active_corner_point_grid(
            mesh::import_grdecl(
                grdecl_vertical_two_cell_fixture(),
                mesh::GrdeclImportOptions{
                    1.0, 1.0e-15}));
    require_shared_operator_axis(
        vertical,
        {0.0, 0.0, 0.5},
        {0.0, 0.0, -0.5},
        {0.0, 0.0, 1.0});
}

void cell_face_geometric_operator_3d_skewed() {
    const auto processed =
        mesh::process_active_corner_point_grid(
            mesh::import_grdecl(
                grdecl_skewed_two_cell_fixture(),
                mesh::GrdeclImportOptions{
                    1.0, 1.0e-15}));

    require_active_processor_shape(
        processed);
    require_close(
        processed.cell_volumes_m3[0],
        1.0,
        1.0e-14,
        "skewed first cell volume");
    require_close(
        processed.cell_volumes_m3[1],
        1.0,
        1.0e-14,
        "skewed second cell volume");

    const auto op =
        mesh::make_cell_face_geometric_operator_3d(
            processed.topology,
            processed.vertex_coordinates_m,
            processed.face_geometry);
    const auto interface =
        only_shared_face(
            processed.topology);

    require_close(
        op.cell_centroid_m(
            mesh::LocalIndex{0U}).x_m,
        0.5,
        1.0e-14,
        "skewed owner 8-corner centroid x");
    require_close(
        op.cell_centroid_m(
            mesh::LocalIndex{0U}).y_m,
        0.5,
        1.0e-14,
        "skewed owner 8-corner centroid y");
    require_close(
        op.cell_centroid_m(
            mesh::LocalIndex{1U}).x_m,
        1.5,
        1.0e-14,
        "skewed neighbour 8-corner centroid x");
    require_close(
        op.cell_centroid_m(
            mesh::LocalIndex{1U}).y_m,
        0.5,
        1.0e-14,
        "skewed neighbour 8-corner centroid y");

    const auto face_centroid =
        processed.face_geometry
            .face_centroid_m(interface);
    require_close(
        face_centroid.x_m,
        1.0,
        1.0e-14,
        "skewed interface centroid x");
    require_close(
        face_centroid.y_m,
        0.5,
        1.0e-14,
        "skewed interface centroid y");
    require_close(
        face_centroid.z_m,
        0.5,
        1.0e-14,
        "skewed interface centroid z");

    const double sqrt29 =
        std::sqrt(29.0);
    const double expected_area =
        sqrt29 / 5.0;
    const double expected_normal_x =
        5.0 / sqrt29;
    const double expected_normal_y =
        2.0 / sqrt29;
    const double expected_distance =
        2.5 / sqrt29;
    const double expected_angle =
        std::acos(expected_normal_x);

    require_close(
        processed.face_geometry
            .face_area_m2(interface),
        expected_area,
        1.0e-14,
        "skewed interface area");
    const auto normal =
        processed.face_geometry
            .face_owner_unit_normal(
                interface);
    require_close(
        normal.x,
        expected_normal_x,
        1.0e-14,
        "skewed interface normal x");
    require_close(
        normal.y,
        expected_normal_y,
        1.0e-14,
        "skewed interface normal y");
    require_close(
        normal.z,
        0.0,
        1.0e-14,
        "skewed interface normal z");

    const auto owner_displacement =
        op.owner_to_face_displacement_m(
            interface);
    const auto neighbour_displacement =
        op.neighbour_to_face_displacement_m(
            interface);
    require(
        neighbour_displacement.has_value(),
        "skewed interface neighbour displacement");
    require_close(
        owner_displacement.x_m,
        0.5,
        1.0e-14,
        "skewed owner displacement x");
    require_close(
        owner_displacement.y_m,
        0.0,
        1.0e-14,
        "skewed owner displacement y");
    require_close(
        neighbour_displacement->x_m,
        -0.5,
        1.0e-14,
        "skewed neighbour displacement x");
    require_close(
        neighbour_displacement->y_m,
        0.0,
        1.0e-14,
        "skewed neighbour displacement y");

    require_close(
        op.owner_normal_distance_m(
            interface),
        expected_distance,
        1.0e-14,
        "skewed owner positive normal distance");
    require(
        op.neighbour_normal_distance_m(
            interface).has_value(),
        "skewed neighbour normal distance exists");
    require_close(
        *op.neighbour_normal_distance_m(
            interface),
        expected_distance,
        1.0e-14,
        "skewed neighbour positive normal distance");
    require(
        op.owner_normal_distance_m(interface) > 0.0 &&
            *op.neighbour_normal_distance_m(interface) > 0.0,
        "skewed two-sided normal distances remain positive");

    const auto non_orthogonality =
        op.internal_non_orthogonality(
            interface);
    require(
        non_orthogonality.has_value(),
        "skewed internal face non-orthogonality exists");
    require_close(
        non_orthogonality
            ->owner_to_neighbour_displacement_m.x_m,
        1.0,
        1.0e-14,
        "skewed centre vector x");
    require_close(
        non_orthogonality
            ->owner_to_neighbour_displacement_m.y_m,
        0.0,
        1.0e-14,
        "skewed centre vector y");
    require_close(
        non_orthogonality->center_distance_m,
        1.0,
        1.0e-14,
        "skewed centre distance");
    require_close(
        non_orthogonality
            ->normal_alignment_cosine,
        expected_normal_x,
        1.0e-14,
        "skewed normal alignment cosine");
    require_close(
        non_orthogonality->angle_rad,
        expected_angle,
        1.0e-14,
        "skewed non-orthogonality angle");
    require(
        non_orthogonality->angle_rad >
                20.0 * std::acos(-1.0) / 180.0 &&
            non_orthogonality->angle_rad <
                25.0 * std::acos(-1.0) / 180.0,
        "skewed fixture must be visibly non-orthogonal");

    const auto transmissibility_geometry =
        op.transmissibility_geometry(
            interface);
    require(
        transmissibility_geometry
            .internal_non_orthogonality
            .has_value(),
        "skewed transmissibility geometry exposes non-orthogonality");
    require_close(
        transmissibility_geometry
            .internal_non_orthogonality
            ->angle_rad,
        expected_angle,
        1.0e-14,
        "skewed transmissibility non-orthogonality angle");
}

mesh::FaceGeometry3D modified_face_geometry(
    const mesh::ActiveCornerPointGrid& processed,
    std::vector<mesh::Coordinate3D> centroids,
    std::vector<mesh::UnitVector3D> normals) {
    return mesh::FaceGeometry3D{
        processed.face_geometry.cell_count(),
        std::move(centroids),
        std::vector<double>{
            processed.face_geometry.face_areas_m2().begin(),
            processed.face_geometry.face_areas_m2().end()},
        std::vector<mesh::LocalIndex>{
            processed.face_geometry.face_owners().begin(),
            processed.face_geometry.face_owners().end()},
        std::move(normals)};
}

void cell_face_geometric_operator_3d_invalid() {
    const auto processed =
        mesh::process_active_corner_point_grid(
            mesh::import_grdecl(
                grdecl_two_cell_all_active_fixture(),
                mesh::GrdeclImportOptions{
                    1.0, 1.0e-15}));
    const auto baseline =
        mesh::make_cell_face_geometric_operator_3d(
            processed.topology,
            processed.vertex_coordinates_m,
            processed.face_geometry);
    const auto interface =
        only_shared_face(
            processed.topology);
    const std::size_t face =
        static_cast<std::size_t>(
            interface.value());
    const auto owner =
        baseline.face_owner(interface);
    const auto neighbour =
        *baseline.face_neighbour(interface);

    const std::vector<mesh::Coordinate3D>
        original_centroids{
            processed.face_geometry.face_centroids_m().begin(),
            processed.face_geometry.face_centroids_m().end()};
    const std::vector<mesh::UnitVector3D>
        original_normals{
            processed.face_geometry
                .face_owner_unit_normals()
                .begin(),
            processed.face_geometry
                .face_owner_unit_normals()
                .end()};

    {
        auto centroids =
            original_centroids;
        centroids[face] =
            baseline.cell_centroid_m(owner);
        const auto geometry =
            modified_face_geometry(
                processed,
                std::move(centroids),
                original_normals);
        expect_throw<std::invalid_argument>(
            [&] {
                (void)mesh::make_cell_face_geometric_operator_3d(
                    processed.topology,
                    processed.vertex_coordinates_m,
                    geometry);
            });
    }

    {
        auto normals =
            original_normals;
        normals[face].x =
            -normals[face].x;
        normals[face].y =
            -normals[face].y;
        normals[face].z =
            -normals[face].z;
        const auto geometry =
            modified_face_geometry(
                processed,
                original_centroids,
                std::move(normals));
        expect_throw<std::invalid_argument>(
            [&] {
                (void)mesh::make_cell_face_geometric_operator_3d(
                    processed.topology,
                    processed.vertex_coordinates_m,
                    geometry);
            });
    }

    {
        auto centroids =
            original_centroids;
        centroids[face] =
            baseline.cell_centroid_m(
                neighbour);
        const auto geometry =
            modified_face_geometry(
                processed,
                std::move(centroids),
                original_normals);
        expect_throw<std::invalid_argument>(
            [&] {
                (void)mesh::make_cell_face_geometric_operator_3d(
                    processed.topology,
                    processed.vertex_coordinates_m,
                    geometry);
            });
    }
}

void dof_layout_snapshot() {
    const auto topology = mesh::make_cartesian_topology_2d(2U, 1U);
    const auto layout = mesh::DofLayout::create(
        topology,
        {
            {"cell.primary", mesh::EntityKind::cell, 2U},
            {"vertex.aux", mesh::EntityKind::vertex, 1U},
            {"cell.secondary", mesh::EntityKind::cell, 1U},
            {"face.trace", mesh::EntityKind::face, 2U},
        });

    require(layout.variable_count() == 4U, "DoF variable count");
    require(layout.total_dof_count() == 26U, "total DoF count");

    require(layout.entity_count(mesh::EntityKind::cell) == 2U, "cell entity count");
    require(layout.entity_count(mesh::EntityKind::face) == 7U, "face entity count");
    require(layout.entity_count(mesh::EntityKind::vertex) == 6U, "vertex entity count");

    require(layout.dofs_per_entity(mesh::EntityKind::cell) == 3U,
            "cell DoFs per entity");
    require(layout.dofs_per_entity(mesh::EntityKind::face) == 2U,
            "face DoFs per entity");
    require(layout.dofs_per_entity(mesh::EntityKind::vertex) == 1U,
            "vertex DoFs per entity");

    require(layout.location_offset(mesh::EntityKind::cell) == 0U, "cell block offset");
    require(layout.location_scalar_count(mesh::EntityKind::cell) == 6U,
            "cell block size");
    require(layout.location_offset(mesh::EntityKind::face) == 6U, "face block offset");
    require(layout.location_scalar_count(mesh::EntityKind::face) == 14U,
            "face block size");
    require(layout.location_offset(mesh::EntityKind::vertex) == 20U,
            "vertex block offset");
    require(layout.location_scalar_count(mesh::EntityKind::vertex) == 6U,
            "vertex block size");

    require(layout.entity_offset(mesh::EntityKind::cell, mesh::LocalIndex{0U}) == 0U,
            "cell zero entity offset");
    require(layout.entity_offset(mesh::EntityKind::cell, mesh::LocalIndex{1U}) == 3U,
            "cell one entity offset");
    require(layout.entity_offset(mesh::EntityKind::face, mesh::LocalIndex{3U}) == 12U,
            "face entity offset");
    require(layout.entity_offset(mesh::EntityKind::vertex, mesh::LocalIndex{5U}) == 25U,
            "vertex entity offset");

    require(layout.variable(0U).id == "cell.primary", "variable order changed");
    require(layout.variable(1U).id == "vertex.aux", "variable order changed");
    require(layout.variable(2U).id == "cell.secondary", "variable order changed");
    require(layout.variable(3U).id == "face.trace", "variable order changed");
    require(layout.contains("cell.secondary"), "variable ID lookup");
    require(layout.variable_index("cell.secondary") == 2U, "variable ID index");

    require(layout.scalar_offset(0U, mesh::LocalIndex{0U}, 0U) == 0U,
            "cell primary c0 component0");
    require(layout.scalar_offset(0U, mesh::LocalIndex{0U}, 1U) == 1U,
            "cell primary c0 component1");
    require(layout.scalar_offset(2U, mesh::LocalIndex{0U}, 0U) == 2U,
            "cell secondary c0");
    require(layout.scalar_offset(0U, mesh::LocalIndex{1U}, 0U) == 3U,
            "cell primary c1 component0");
    require(layout.scalar_offset(2U, mesh::LocalIndex{1U}, 0U) == 5U,
            "cell secondary c1");

    require(layout.scalar_offset(3U, mesh::LocalIndex{0U}, 0U) == 6U,
            "face trace first component");
    require(layout.scalar_offset(3U, mesh::LocalIndex{6U}, 1U) == 19U,
            "face trace final component");

    require(layout.scalar_offset(1U, mesh::LocalIndex{0U}, 0U) == 20U,
            "vertex block first");
    require(layout.scalar_offset("vertex.aux", mesh::LocalIndex{5U}, 0U) == 25U,
            "vertex block final");

    // The snapshot owns its variable IDs and remains usable after input destruction.
    require(layout.variables().front().id == "cell.primary", "owned variable metadata");
}

void dof_layout_invalid() {
    const auto topology = mesh::make_cartesian_topology_2d(1U, 1U);

    expect_throw<std::invalid_argument>(
        [&] { (void)mesh::DofLayout::create(topology, {}); });
    expect_throw<std::invalid_argument>([&] {
        (void)mesh::DofLayout::create(
            topology, {{"", mesh::EntityKind::cell, 1U}});
    });
    expect_throw<std::invalid_argument>([&] {
        (void)mesh::DofLayout::create(
            topology, {{" \t", mesh::EntityKind::cell, 1U}});
    });
    expect_throw<std::invalid_argument>([&] {
        (void)mesh::DofLayout::create(
            topology,
            {{std::string{"bad\0id", 6U}, mesh::EntityKind::cell, 1U}});
    });
    expect_throw<std::invalid_argument>([&] {
        (void)mesh::DofLayout::create(
            topology, {{"zero", mesh::EntityKind::cell, 0U}});
    });
    expect_throw<std::invalid_argument>([&] {
        (void)mesh::DofLayout::create(
            topology,
            {{"same", mesh::EntityKind::cell, 1U},
             {"same", mesh::EntityKind::face, 1U}});
    });
    expect_throw<std::invalid_argument>([&] {
        (void)mesh::DofLayout::create(
            topology, {{"edge", mesh::EntityKind::edge, 1U}});
    });
    expect_throw<std::invalid_argument>([&] {
        (void)mesh::DofLayout::create(
            topology,
            {{"bad-kind", static_cast<mesh::EntityKind>(255U), 1U}});
    });
    expect_throw<std::length_error>([&] {
        (void)mesh::DofLayout::create(
            topology,
            {{"huge", mesh::EntityKind::cell,
              std::numeric_limits<std::size_t>::max()},
             {"overflow", mesh::EntityKind::cell, 1U}});
    });
    expect_throw<std::length_error>([&] {
        (void)mesh::DofLayout::create(
            topology,
            {{"face-huge", mesh::EntityKind::face,
              std::numeric_limits<std::size_t>::max()}});
    });

    const auto layout = mesh::DofLayout::create(
        topology, {{"cell.x", mesh::EntityKind::cell, 2U}});
    expect_throw<std::out_of_range>(
        [&] { (void)layout.variable(1U); });
    expect_throw<std::out_of_range>(
        [&] { (void)layout.variable_index("missing"); });
    expect_throw<std::out_of_range>(
        [&] { (void)layout.entity_offset(mesh::EntityKind::cell, mesh::LocalIndex{1U}); });
    expect_throw<std::out_of_range>(
        [&] { (void)layout.scalar_offset(0U, mesh::LocalIndex{0U}, 2U); });
    expect_throw<std::out_of_range>(
        [&] { (void)layout.scalar_offset(1U, mesh::LocalIndex{0U}, 0U); });
    expect_throw<std::invalid_argument>(
        [&] { (void)layout.dofs_per_entity(mesh::EntityKind::edge); });
    expect_throw<std::invalid_argument>([&] {
        (void)layout.location_offset(static_cast<mesh::EntityKind>(255U));
    });
}


void partition_serial_snapshot() {
    const auto snapshot = [] {
        const auto topology = mesh::make_cartesian_topology_2d(2U, 1U);
        return mesh::make_serial_partition_snapshot(topology);
    }();

    require(snapshot.is_serial(), "serial partition flag");
    require(snapshot.rank_count() == 1U, "serial rank count");
    require(snapshot.local_rank().value() == 0U, "serial local rank");

    for (const mesh::EntityKind kind :
         {mesh::EntityKind::vertex, mesh::EntityKind::edge,
          mesh::EntityKind::face, mesh::EntityKind::cell}) {
        require(snapshot.owned_count(kind) == snapshot.entity_count(kind),
                "serial snapshot must own every entity");
        require(snapshot.ghost_count(kind) == 0U,
                "serial snapshot cannot contain ghosts");
        const auto ids = snapshot.global_ids(kind);
        const auto owners = snapshot.owner_ranks(kind);
        require(ids.size() == owners.size(), "serial owner/global alignment");
        for (std::size_t local = 0; local < ids.size(); ++local) {
            const auto index =
                mesh::LocalIndex{static_cast<mesh::LocalIndex::value_type>(local)};
            require(owners[local].value() == 0U, "serial owner must be rank zero");
            require(snapshot.is_owned(kind, index), "serial entity must be owned");
            require(!snapshot.is_ghost(kind, index), "serial entity cannot be ghost");
            require(snapshot.ownership(kind, index) == mesh::EntityOwnership::owned,
                    "serial ownership enum");
            require(snapshot.local_index(kind, ids[local]) == index,
                    "serial global-to-local roundtrip");
            require(snapshot.global_id(kind, index) == ids[local],
                    "serial local-to-global roundtrip");
            require(snapshot.contains_global(kind, ids[local]),
                    "serial global presence");
        }
    }
}

void partition_local_snapshot() {
    mesh::Topology::EntityIds ids;
    ids.vertices = {
        mesh::GlobalEntityId{900U},
        mesh::GlobalEntityId{100U},
        mesh::GlobalEntityId{500U}};
    ids.faces = {
        mesh::GlobalEntityId{42U},
        mesh::GlobalEntityId{5U}};
    ids.cells = {
        mesh::GlobalEntityId{77U},
        mesh::GlobalEntityId{9U}};
    const mesh::Topology topology{std::move(ids), {}};

    mesh::EntityOwnerRanks owners;
    owners.vertices = {
        mesh::PartitionRank{1U}, mesh::PartitionRank{0U}, mesh::PartitionRank{2U}};
    owners.faces = {
        mesh::PartitionRank{1U}, mesh::PartitionRank{2U}};
    owners.cells = {
        mesh::PartitionRank{0U}, mesh::PartitionRank{1U}};

    const auto snapshot = mesh::PartitionSnapshot::create(
        topology, mesh::PartitionRank{1U}, 3U, std::move(owners));

    require(!snapshot.is_serial(), "local multi-rank snapshot misclassified");
    require(snapshot.local_rank().value() == 1U, "local rank");
    require(snapshot.rank_count() == 3U, "rank count");

    require(snapshot.owned_count(mesh::EntityKind::vertex) == 1U,
            "vertex owned count");
    require(snapshot.ghost_count(mesh::EntityKind::vertex) == 2U,
            "vertex ghost count");
    require(snapshot.owned_count(mesh::EntityKind::face) == 1U,
            "face owned count");
    require(snapshot.ghost_count(mesh::EntityKind::face) == 1U,
            "face ghost count");
    require(snapshot.owned_count(mesh::EntityKind::cell) == 1U,
            "cell owned count");
    require(snapshot.ghost_count(mesh::EntityKind::cell) == 1U,
            "cell ghost count");

    require(snapshot.global_id(mesh::EntityKind::vertex, mesh::LocalIndex{0U}).value() == 900U,
            "local-to-global must preserve topology ordering");
    require(snapshot.local_index(mesh::EntityKind::vertex, mesh::GlobalEntityId{100U}).value() == 1U,
            "global-to-local must not assume numeric identity");
    require(snapshot.local_index(mesh::EntityKind::vertex, mesh::GlobalEntityId{500U}).value() == 2U,
            "global-to-local sorted lookup");
    require(snapshot.local_index(mesh::EntityKind::face, mesh::GlobalEntityId{5U}).value() == 1U,
            "face global-to-local lookup");
    require(snapshot.local_index(mesh::EntityKind::cell, mesh::GlobalEntityId{77U}).value() == 0U,
            "cell global-to-local lookup");

    require(snapshot.owner_rank(mesh::EntityKind::vertex, mesh::LocalIndex{0U}).value() == 1U,
            "owned vertex owner rank");
    require(snapshot.owner_rank(mesh::EntityKind::vertex, mesh::LocalIndex{1U}).value() == 0U,
            "ghost vertex owner rank");
    require(snapshot.is_owned(mesh::EntityKind::vertex, mesh::LocalIndex{0U}),
            "owned vertex classification");
    require(snapshot.is_ghost(mesh::EntityKind::vertex, mesh::LocalIndex{1U}),
            "ghost vertex classification");
    require(snapshot.is_ghost(mesh::EntityKind::vertex, mesh::LocalIndex{2U}),
            "second ghost vertex classification");
    require(snapshot.is_owned(mesh::EntityKind::face, mesh::LocalIndex{0U}),
            "owned face classification");
    require(snapshot.is_ghost(mesh::EntityKind::face, mesh::LocalIndex{1U}),
            "ghost face classification");
    require(snapshot.is_ghost(mesh::EntityKind::cell, mesh::LocalIndex{0U}),
            "ghost cell classification");
    require(snapshot.is_owned(mesh::EntityKind::cell, mesh::LocalIndex{1U}),
            "owned cell classification");

    require(!snapshot.contains_global(
                mesh::EntityKind::vertex, mesh::GlobalEntityId{12345U}),
            "missing global ID reported present");
    expect_throw<std::out_of_range>([&] {
        (void)snapshot.local_index(
            mesh::EntityKind::vertex, mesh::GlobalEntityId{12345U});
    });
}

void partition_invalid() {
    const auto topology = mesh::make_cartesian_topology_2d(1U, 1U);

    mesh::EntityOwnerRanks serial_owners;
    serial_owners.vertices.assign(4U, mesh::PartitionRank{0U});
    serial_owners.faces.assign(4U, mesh::PartitionRank{0U});
    serial_owners.cells.assign(1U, mesh::PartitionRank{0U});

    expect_throw<std::invalid_argument>([&] {
        auto owners = serial_owners;
        (void)mesh::PartitionSnapshot::create(
            topology, mesh::PartitionRank{0U}, 0U, std::move(owners));
    });
    expect_throw<std::out_of_range>([&] {
        auto owners = serial_owners;
        (void)mesh::PartitionSnapshot::create(
            topology, mesh::PartitionRank{1U}, 1U, std::move(owners));
    });
    expect_throw<std::invalid_argument>([&] {
        auto owners = serial_owners;
        owners.faces.pop_back();
        (void)mesh::PartitionSnapshot::create(
            topology, mesh::PartitionRank{0U}, 1U, std::move(owners));
    });
    expect_throw<std::out_of_range>([&] {
        auto owners = serial_owners;
        owners.vertices[2] = mesh::PartitionRank{2U};
        (void)mesh::PartitionSnapshot::create(
            topology, mesh::PartitionRank{0U}, 2U, std::move(owners));
    });

    const auto snapshot = mesh::make_serial_partition_snapshot(topology);
    expect_throw<std::out_of_range>([&] {
        (void)snapshot.global_id(
            mesh::EntityKind::cell, mesh::LocalIndex{1U});
    });
    expect_throw<std::out_of_range>([&] {
        (void)snapshot.owner_rank(
            mesh::EntityKind::face, mesh::LocalIndex{4U});
    });
    expect_throw<std::out_of_range>([&] {
        (void)snapshot.local_index(
            mesh::EntityKind::cell, mesh::GlobalEntityId{999U});
    });
    expect_throw<std::invalid_argument>([&] {
        (void)snapshot.entity_count(static_cast<mesh::EntityKind>(255U));
    });
    expect_throw<std::invalid_argument>([&] {
        (void)snapshot.contains_global(
            static_cast<mesh::EntityKind>(255U), mesh::GlobalEntityId{0U});
    });
}


void dof_numbering_serial() {
    const auto topology = mesh::make_cartesian_topology_2d(2U, 1U);
    const auto partition = mesh::make_serial_partition_snapshot(topology);
    const auto layout = mesh::DofLayout::create(
        topology,
        {
            {"cell.primary", mesh::EntityKind::cell, 2U},
            {"face.trace", mesh::EntityKind::face, 1U},
            {"vertex.aux", mesh::EntityKind::vertex, 1U},
            {"cell.secondary", mesh::EntityKind::cell, 1U},
        });
    const auto numbering =
        mesh::DofNumberingSnapshot::create_serial(layout, partition);

    require(numbering.is_serial(), "serial DoF numbering flag");
    require(numbering.local_rank().value() == 0U, "serial numbering local rank");
    require(numbering.rank_count() == 1U, "serial numbering rank count");
    require(numbering.local_dof_count() == layout.total_dof_count(),
            "serial local DoF count");
    require(numbering.global_dof_count() ==
                static_cast<std::uint64_t>(layout.total_dof_count()),
            "serial global DoF count");
    require(numbering.owned_dof_count() == layout.total_dof_count(),
            "serial owned DoF count");
    require(numbering.ghost_dof_count() == 0U,
            "serial ghost DoF count");

    require(numbering.global_entity_count(mesh::EntityKind::cell) == 2U,
            "serial global cell count");
    require(numbering.global_entity_count(mesh::EntityKind::face) == 7U,
            "serial global face count");
    require(numbering.global_entity_count(mesh::EntityKind::vertex) == 6U,
            "serial global vertex count");

    require(numbering.global_location_offset(mesh::EntityKind::cell) == 0U,
            "serial cell global block");
    require(numbering.global_location_dof_count(mesh::EntityKind::cell) == 6U,
            "serial cell global block size");
    require(numbering.global_location_offset(mesh::EntityKind::face) == 6U,
            "serial face global block");
    require(numbering.global_location_offset(mesh::EntityKind::vertex) == 13U,
            "serial vertex global block");

    for (std::size_t local = 0; local < layout.total_dof_count(); ++local) {
        const auto global = numbering.global_index(local);
        require(global.value() == static_cast<std::uint64_t>(local),
                "serial global numbering must equal local scalar ordering");
        require(numbering.local_scalar(global) == local,
                "serial scalar roundtrip");
        require(numbering.contains_global(global),
                "serial global scalar presence");
        require(numbering.is_owned(local),
                "serial local scalar must be owned");
        require(!numbering.is_ghost(local),
                "serial local scalar cannot be ghost");
    }
}

void dof_numbering_local() {
    mesh::Topology::EntityIds ids;
    ids.vertices = {
        mesh::GlobalEntityId{900U},
        mesh::GlobalEntityId{100U},
        mesh::GlobalEntityId{500U}};
    ids.faces = {
        mesh::GlobalEntityId{42U},
        mesh::GlobalEntityId{5U}};
    ids.cells = {
        mesh::GlobalEntityId{77U},
        mesh::GlobalEntityId{9U}};
    const mesh::Topology topology{std::move(ids), {}};

    mesh::EntityOwnerRanks owners;
    owners.vertices = {
        mesh::PartitionRank{1U}, mesh::PartitionRank{0U}, mesh::PartitionRank{2U}};
    owners.faces = {
        mesh::PartitionRank{1U}, mesh::PartitionRank{2U}};
    owners.cells = {
        mesh::PartitionRank{0U}, mesh::PartitionRank{1U}};
    const auto partition = mesh::PartitionSnapshot::create(
        topology, mesh::PartitionRank{1U}, 3U, std::move(owners));

    const auto layout = mesh::DofLayout::create(
        topology,
        {
            {"cell.v", mesh::EntityKind::cell, 2U},
            {"face.v", mesh::EntityKind::face, 1U},
            {"vertex.v", mesh::EntityKind::vertex, 1U},
        });

    mesh::GlobalEntityNumberingInput input;
    input.global_cell_count = 5U;
    input.global_face_count = 7U;
    input.global_vertex_count = 8U;
    input.cells = {
        {mesh::GlobalEntityId{77U}, mesh::GlobalEntityOrdinal{4U}},
        {mesh::GlobalEntityId{9U}, mesh::GlobalEntityOrdinal{1U}}};
    input.faces = {
        {mesh::GlobalEntityId{42U}, mesh::GlobalEntityOrdinal{6U}},
        {mesh::GlobalEntityId{5U}, mesh::GlobalEntityOrdinal{2U}}};
    input.vertices = {
        {mesh::GlobalEntityId{900U}, mesh::GlobalEntityOrdinal{7U}},
        {mesh::GlobalEntityId{100U}, mesh::GlobalEntityOrdinal{0U}},
        {mesh::GlobalEntityId{500U}, mesh::GlobalEntityOrdinal{4U}}};

    const auto numbering = mesh::DofNumberingSnapshot::create_local(
        layout, partition, std::move(input));

    require(!numbering.is_serial(), "local numbering misclassified as serial");
    require(numbering.local_dof_count() == 9U, "local numbering local count");
    require(numbering.global_dof_count() == 25U, "local numbering global count");
    require(numbering.owned_dof_count() == 4U, "owned local DoF count");
    require(numbering.ghost_dof_count() == 5U, "ghost local DoF count");

    const std::array<std::uint64_t, 9> expected_global{
        8U, 9U, 2U, 3U, 16U, 12U, 24U, 17U, 21U};
    const std::array<mesh::EntityOwnership, 9> expected_ownership{
        mesh::EntityOwnership::ghost,
        mesh::EntityOwnership::ghost,
        mesh::EntityOwnership::owned,
        mesh::EntityOwnership::owned,
        mesh::EntityOwnership::owned,
        mesh::EntityOwnership::ghost,
        mesh::EntityOwnership::owned,
        mesh::EntityOwnership::ghost,
        mesh::EntityOwnership::ghost};

    for (std::size_t local = 0; local < expected_global.size(); ++local) {
        const auto global = numbering.global_index(local);
        require(global.value() == expected_global[local],
                "generic local-to-global DoF mapping");
        require(numbering.ownership(local) == expected_ownership[local],
                "generic DoF ownership");
        require(numbering.local_scalar(global) == local,
                "generic global-to-local DoF roundtrip");
    }

    require(numbering.global_location_offset(mesh::EntityKind::cell) == 0U,
            "generic cell global block");
    require(numbering.global_location_offset(mesh::EntityKind::face) == 10U,
            "generic face global block");
    require(numbering.global_location_offset(mesh::EntityKind::vertex) == 17U,
            "generic vertex global block");

    require(!numbering.contains_global(mesh::GlobalDofIndex{0U}),
            "nonlocal global cell DoF reported local");
    require(!numbering.contains_global(mesh::GlobalDofIndex{10U}),
            "nonlocal global face DoF reported local");
    require(!numbering.contains_global(mesh::GlobalDofIndex{18U}),
            "nonlocal global vertex DoF reported local");
    expect_throw<std::out_of_range>(
        [&] { (void)numbering.local_scalar(mesh::GlobalDofIndex{0U}); });
    expect_throw<std::out_of_range>(
        [&] { (void)numbering.local_scalar(mesh::GlobalDofIndex{25U}); });
}

void dof_numbering_invalid() {
    const auto topology = mesh::make_cartesian_topology_2d(1U, 1U);
    const auto serial_partition = mesh::make_serial_partition_snapshot(topology);
    const auto layout = mesh::DofLayout::create(
        topology, {{"cell.x", mesh::EntityKind::cell, 2U}});

    mesh::EntityOwnerRanks parallel_owners;
    parallel_owners.vertices.assign(4U, mesh::PartitionRank{1U});
    parallel_owners.faces.assign(4U, mesh::PartitionRank{1U});
    parallel_owners.cells.assign(1U, mesh::PartitionRank{0U});
    const auto parallel_partition = mesh::PartitionSnapshot::create(
        topology, mesh::PartitionRank{0U}, 2U, std::move(parallel_owners));
    expect_throw<std::invalid_argument>([&] {
        (void)mesh::DofNumberingSnapshot::create_serial(
            layout, parallel_partition);
    });

    const auto larger_topology = mesh::make_cartesian_topology_2d(2U, 1U);
    const auto larger_partition =
        mesh::make_serial_partition_snapshot(larger_topology);
    expect_throw<std::invalid_argument>([&] {
        mesh::GlobalEntityNumberingInput input;
        input.global_cell_count = 2U;
        input.global_face_count = 7U;
        input.global_vertex_count = 6U;
        for (std::size_t i = 0; i < 2U; ++i) {
            input.cells.push_back({
                larger_partition.global_id(
                    mesh::EntityKind::cell,
                    mesh::LocalIndex{
                        static_cast<mesh::LocalIndex::value_type>(i)}),
                mesh::GlobalEntityOrdinal{static_cast<std::uint64_t>(i)}});
        }
        for (std::size_t i = 0; i < 7U; ++i) {
            input.faces.push_back({
                larger_partition.global_id(
                    mesh::EntityKind::face,
                    mesh::LocalIndex{
                        static_cast<mesh::LocalIndex::value_type>(i)}),
                mesh::GlobalEntityOrdinal{static_cast<std::uint64_t>(i)}});
        }
        for (std::size_t i = 0; i < 6U; ++i) {
            input.vertices.push_back({
                larger_partition.global_id(
                    mesh::EntityKind::vertex,
                    mesh::LocalIndex{
                        static_cast<mesh::LocalIndex::value_type>(i)}),
                mesh::GlobalEntityOrdinal{static_cast<std::uint64_t>(i)}});
        }
        (void)mesh::DofNumberingSnapshot::create_local(
            layout, larger_partition, std::move(input));
    });

    expect_throw<std::invalid_argument>([&] {
        mesh::GlobalEntityNumberingInput input;
        input.global_cell_count = 1U;
        input.global_face_count = 4U;
        input.global_vertex_count = 4U;
        input.cells = {};
        input.faces = {
            {mesh::GlobalEntityId{0U}, mesh::GlobalEntityOrdinal{0U}},
            {mesh::GlobalEntityId{1U}, mesh::GlobalEntityOrdinal{1U}},
            {mesh::GlobalEntityId{2U}, mesh::GlobalEntityOrdinal{2U}},
            {mesh::GlobalEntityId{3U}, mesh::GlobalEntityOrdinal{3U}}};
        input.vertices = {
            {mesh::GlobalEntityId{0U}, mesh::GlobalEntityOrdinal{0U}},
            {mesh::GlobalEntityId{1U}, mesh::GlobalEntityOrdinal{1U}},
            {mesh::GlobalEntityId{2U}, mesh::GlobalEntityOrdinal{2U}},
            {mesh::GlobalEntityId{3U}, mesh::GlobalEntityOrdinal{3U}}};
        (void)mesh::DofNumberingSnapshot::create_local(
            layout, serial_partition, std::move(input));
    });

    expect_throw<std::invalid_argument>([&] {
        mesh::GlobalEntityNumberingInput input;
        input.global_cell_count = 1U;
        input.global_face_count = 4U;
        input.global_vertex_count = 4U;
        input.cells = {
            {mesh::GlobalEntityId{999U}, mesh::GlobalEntityOrdinal{0U}}};
        for (std::uint64_t i = 0U; i < 4U; ++i) {
            input.faces.push_back({
                mesh::GlobalEntityId{i}, mesh::GlobalEntityOrdinal{i}});
            input.vertices.push_back({
                mesh::GlobalEntityId{i}, mesh::GlobalEntityOrdinal{i}});
        }
        (void)mesh::DofNumberingSnapshot::create_local(
            layout, serial_partition, std::move(input));
    });

    expect_throw<std::out_of_range>([&] {
        mesh::GlobalEntityNumberingInput input;
        input.global_cell_count = 1U;
        input.global_face_count = 4U;
        input.global_vertex_count = 4U;
        input.cells = {
            {mesh::GlobalEntityId{0U}, mesh::GlobalEntityOrdinal{1U}}};
        for (std::uint64_t i = 0U; i < 4U; ++i) {
            input.faces.push_back({
                mesh::GlobalEntityId{i}, mesh::GlobalEntityOrdinal{i}});
            input.vertices.push_back({
                mesh::GlobalEntityId{i}, mesh::GlobalEntityOrdinal{i}});
        }
        (void)mesh::DofNumberingSnapshot::create_local(
            layout, serial_partition, std::move(input));
    });

    mesh::Topology::EntityIds duplicate_ids;
    duplicate_ids.cells = {
        mesh::GlobalEntityId{10U}, mesh::GlobalEntityId{11U}};
    const mesh::Topology duplicate_topology{std::move(duplicate_ids), {}};
    mesh::EntityOwnerRanks duplicate_owners;
    duplicate_owners.cells = {
        mesh::PartitionRank{0U}, mesh::PartitionRank{0U}};
    const auto duplicate_partition = mesh::PartitionSnapshot::create(
        duplicate_topology, mesh::PartitionRank{0U}, 1U,
        std::move(duplicate_owners));
    const auto duplicate_layout = mesh::DofLayout::create(
        duplicate_topology, {{"cell.x", mesh::EntityKind::cell, 1U}});
    expect_throw<std::invalid_argument>([&] {
        mesh::GlobalEntityNumberingInput input;
        input.global_cell_count = 2U;
        input.cells = {
            {mesh::GlobalEntityId{10U}, mesh::GlobalEntityOrdinal{0U}},
            {mesh::GlobalEntityId{11U}, mesh::GlobalEntityOrdinal{0U}}};
        (void)mesh::DofNumberingSnapshot::create_local(
            duplicate_layout, duplicate_partition, std::move(input));
    });

    expect_throw<std::length_error>([&] {
        mesh::GlobalEntityNumberingInput input;
        input.global_cell_count =
            std::numeric_limits<std::uint64_t>::max();
        input.global_face_count = 4U;
        input.global_vertex_count = 4U;
        input.cells = {
            {mesh::GlobalEntityId{0U}, mesh::GlobalEntityOrdinal{0U}}};
        for (std::uint64_t i = 0U; i < 4U; ++i) {
            input.faces.push_back({
                mesh::GlobalEntityId{i}, mesh::GlobalEntityOrdinal{i}});
            input.vertices.push_back({
                mesh::GlobalEntityId{i}, mesh::GlobalEntityOrdinal{i}});
        }
        (void)mesh::DofNumberingSnapshot::create_local(
            layout, serial_partition, std::move(input));
    });

    const auto numbering =
        mesh::DofNumberingSnapshot::create_serial(layout, serial_partition);
    expect_throw<std::out_of_range>(
        [&] { (void)numbering.global_index(numbering.local_dof_count()); });
    expect_throw<std::out_of_range>(
        [&] { (void)numbering.ownership(numbering.local_dof_count()); });
    expect_throw<std::invalid_argument>([&] {
        (void)numbering.global_entity_count(mesh::EntityKind::edge);
    });
    expect_throw<std::invalid_argument>([&] {
        (void)numbering.global_location_offset(
            static_cast<mesh::EntityKind>(255U));
    });
}


mesh::PartitionSnapshot shared_fixture_partition(std::uint32_t rank) {
    mesh::Topology::EntityIds ids;
    if (rank == 0U) {
        ids.vertices = {
            mesh::GlobalEntityId{10U}, mesh::GlobalEntityId{20U}};
        ids.faces = {
            mesh::GlobalEntityId{100U}, mesh::GlobalEntityId{200U}};
        ids.cells = {
            mesh::GlobalEntityId{1000U}, mesh::GlobalEntityId{2000U}};
    } else if (rank == 1U) {
        ids.vertices = {
            mesh::GlobalEntityId{20U}, mesh::GlobalEntityId{10U}};
        ids.faces = {
            mesh::GlobalEntityId{200U}, mesh::GlobalEntityId{100U}};
        ids.cells = {
            mesh::GlobalEntityId{2000U}, mesh::GlobalEntityId{1000U}};
    } else {
        throw std::invalid_argument("shared fixture rank must be 0 or 1");
    }

    const mesh::Topology topology{std::move(ids), {}};
    mesh::EntityOwnerRanks owners;
    if (rank == 0U) {
        owners.vertices = {mesh::PartitionRank{0U}, mesh::PartitionRank{1U}};
        owners.faces = {mesh::PartitionRank{0U}, mesh::PartitionRank{1U}};
        owners.cells = {mesh::PartitionRank{0U}, mesh::PartitionRank{1U}};
    } else {
        owners.vertices = {mesh::PartitionRank{1U}, mesh::PartitionRank{0U}};
        owners.faces = {mesh::PartitionRank{1U}, mesh::PartitionRank{0U}};
        owners.cells = {mesh::PartitionRank{1U}, mesh::PartitionRank{0U}};
    }
    return mesh::PartitionSnapshot::create(
        topology, mesh::PartitionRank{rank}, 2U, std::move(owners));
}

std::vector<mesh::SharedEntityLink> shared_fixture_links() {
    return {
        {mesh::EntityKind::vertex, mesh::GlobalEntityId{10U},
         mesh::PartitionRank{0U}, mesh::LocalIndex{0U},
         mesh::PartitionRank{1U}, mesh::LocalIndex{1U}},
        {mesh::EntityKind::vertex, mesh::GlobalEntityId{20U},
         mesh::PartitionRank{1U}, mesh::LocalIndex{0U},
         mesh::PartitionRank{0U}, mesh::LocalIndex{1U}},
        {mesh::EntityKind::face, mesh::GlobalEntityId{100U},
         mesh::PartitionRank{0U}, mesh::LocalIndex{0U},
         mesh::PartitionRank{1U}, mesh::LocalIndex{1U}},
        {mesh::EntityKind::face, mesh::GlobalEntityId{200U},
         mesh::PartitionRank{1U}, mesh::LocalIndex{0U},
         mesh::PartitionRank{0U}, mesh::LocalIndex{1U}},
        {mesh::EntityKind::cell, mesh::GlobalEntityId{1000U},
         mesh::PartitionRank{0U}, mesh::LocalIndex{0U},
         mesh::PartitionRank{1U}, mesh::LocalIndex{1U}},
        {mesh::EntityKind::cell, mesh::GlobalEntityId{2000U},
         mesh::PartitionRank{1U}, mesh::LocalIndex{0U},
         mesh::PartitionRank{0U}, mesh::LocalIndex{1U}},
    };
}

void shared_entity_plan_pair() {
    const auto rank0_partition = shared_fixture_partition(0U);
    const auto rank1_partition = shared_fixture_partition(1U);
    const auto links = shared_fixture_links();

    const auto rank0 = mesh::SharedEntityPlan::create(rank0_partition, links);
    const auto rank1 = mesh::SharedEntityPlan::create(rank1_partition, links);

    require(rank0.local_rank().value() == 0U, "rank0 plan local rank");
    require(rank1.local_rank().value() == 1U, "rank1 plan local rank");
    require(rank0.rank_count() == 2U && rank1.rank_count() == 2U,
            "shared plan rank count");
    require(!rank0.is_serial() && !rank1.is_serial(),
            "2-rank plans cannot be serial");

    require(rank0.neighbor_count() == 1U, "rank0 neighbor count");
    require(rank1.neighbor_count() == 1U, "rank1 neighbor count");
    require(rank0.contains_neighbor(mesh::PartitionRank{1U}),
            "rank0 missing rank1 neighbor");
    require(rank1.contains_neighbor(mesh::PartitionRank{0U}),
            "rank1 missing rank0 neighbor");
    require(!rank0.contains_neighbor(mesh::PartitionRank{0U}),
            "rank0 cannot be its own neighbor");
    require(!rank1.contains_neighbor(mesh::PartitionRank{1U}),
            "rank1 cannot be its own neighbor");

    const auto rank0_sends = rank0.send_entities_to(mesh::PartitionRank{1U});
    const auto rank0_receives =
        rank0.receive_entities_from(mesh::PartitionRank{1U});
    const auto rank1_sends = rank1.send_entities_to(mesh::PartitionRank{0U});
    const auto rank1_receives =
        rank1.receive_entities_from(mesh::PartitionRank{0U});

    require(rank0_sends.size() == 3U && rank0_receives.size() == 3U,
            "rank0 send/receive counts");
    require(rank1_sends.size() == 3U && rank1_receives.size() == 3U,
            "rank1 send/receive counts");
    require(rank0.send_count() == 3U && rank0.receive_count() == 3U,
            "rank0 flat counts");
    require(rank1.send_count() == 3U && rank1.receive_count() == 3U,
            "rank1 flat counts");

    const auto rank0_neighbor = rank0.neighbors().front();
    require(rank0_neighbor.rank.value() == 1U, "rank0 neighbor rank");
    require(rank0_neighbor.send_begin == 0U &&
                rank0_neighbor.send_count == 3U &&
                rank0_neighbor.receive_begin == 0U &&
                rank0_neighbor.receive_count == 3U,
            "rank0 compact neighbor ranges");

    for (std::size_t i = 0; i < rank0_sends.size(); ++i) {
        const auto& send = rank0_sends[i];
        const auto& receive = rank1_receives[i];
        require(send.kind == receive.kind, "rank0->rank1 kind symmetry");
        require(send.global_id == receive.global_id,
                "rank0->rank1 global ID symmetry");
        require(send.local == receive.remote_local,
                "rank0 owner-local must match rank1 remote owner-local");
        require(send.remote_local == receive.local,
                "rank0 remote ghost-local must match rank1 local ghost");
    }
    for (std::size_t i = 0; i < rank1_sends.size(); ++i) {
        const auto& send = rank1_sends[i];
        const auto& receive = rank0_receives[i];
        require(send.kind == receive.kind, "rank1->rank0 kind symmetry");
        require(send.global_id == receive.global_id,
                "rank1->rank0 global ID symmetry");
        require(send.local == receive.remote_local,
                "rank1 owner-local must match rank0 remote owner-local");
        require(send.remote_local == receive.local,
                "rank1 remote ghost-local must match rank0 local ghost");
    }

    require(rank0_receives[0].kind == mesh::EntityKind::vertex &&
                rank0_receives[0].global_id.value() == 20U &&
                rank0_receives[0].local.value() == 1U &&
                rank0_receives[0].remote_local.value() == 0U,
            "rank0 vertex receive owner/ghost contract");
    require(rank0_receives[1].kind == mesh::EntityKind::face &&
                rank0_receives[1].global_id.value() == 200U,
            "rank0 face receive ordering");
    require(rank0_receives[2].kind == mesh::EntityKind::cell &&
                rank0_receives[2].global_id.value() == 2000U,
            "rank0 cell receive ordering");

    const auto serial_topology = mesh::make_cartesian_topology_2d(1U, 1U);
    const auto serial_partition =
        mesh::make_serial_partition_snapshot(serial_topology);
    const std::vector<mesh::SharedEntityLink> no_links;
    const auto serial =
        mesh::SharedEntityPlan::create(serial_partition, no_links);
    require(serial.is_serial(), "empty serial halo plan");
    require(serial.neighbor_count() == 0U &&
                serial.send_count() == 0U &&
                serial.receive_count() == 0U,
            "serial halo plan must be empty");
}

void shared_entity_plan_invalid() {
    const auto rank0_partition = shared_fixture_partition(0U);
    const auto links = shared_fixture_links();

    {
        auto missing = links;
        missing.erase(missing.begin() + 1);
        expect_throw<std::invalid_argument>([&] {
            (void)mesh::SharedEntityPlan::create(rank0_partition, missing);
        });
    }
    {
        auto duplicate = links;
        duplicate.push_back(links.front());
        expect_throw<std::invalid_argument>([&] {
            (void)mesh::SharedEntityPlan::create(rank0_partition, duplicate);
        });
    }
    {
        auto duplicate_identity = links;
        auto repeated = links.front();
        repeated.ghost_local = mesh::LocalIndex{0U};
        duplicate_identity.push_back(repeated);
        expect_throw<std::invalid_argument>([&] {
            (void)mesh::SharedEntityPlan::create(
                rank0_partition, duplicate_identity);
        });
    }
    {
        auto wrong_global = links;
        wrong_global.front().global_id = mesh::GlobalEntityId{999U};
        expect_throw<std::invalid_argument>([&] {
            (void)mesh::SharedEntityPlan::create(rank0_partition, wrong_global);
        });
    }
    {
        auto wrong_owner_local = links;
        wrong_owner_local.front().owner_local = mesh::LocalIndex{1U};
        expect_throw<std::invalid_argument>([&] {
            (void)mesh::SharedEntityPlan::create(
                rank0_partition, wrong_owner_local);
        });
    }
    {
        auto wrong_ghost_owner = links;
        wrong_ghost_owner[1].owner_rank = mesh::PartitionRank{0U};
        expect_throw<std::invalid_argument>([&] {
            (void)mesh::SharedEntityPlan::create(
                rank0_partition, wrong_ghost_owner);
        });
    }
    {
        auto same_rank = links;
        same_rank.front().ghost_rank = mesh::PartitionRank{0U};
        expect_throw<std::invalid_argument>([&] {
            (void)mesh::SharedEntityPlan::create(rank0_partition, same_rank);
        });
    }
    {
        auto outside_rank = links;
        outside_rank.front().ghost_rank = mesh::PartitionRank{2U};
        expect_throw<std::out_of_range>([&] {
            (void)mesh::SharedEntityPlan::create(
                rank0_partition, outside_rank);
        });
    }
    {
        mesh::Topology::EntityIds ids;
        ids.cells = {mesh::GlobalEntityId{1U}};
        const mesh::Topology topology{std::move(ids), {}};
        mesh::EntityOwnerRanks owners;
        owners.cells = {mesh::PartitionRank{0U}};
        const auto partition = mesh::PartitionSnapshot::create(
            topology, mesh::PartitionRank{0U}, 3U, std::move(owners));
        const std::array<mesh::SharedEntityLink, 1> unrelated{{
            {mesh::EntityKind::cell, mesh::GlobalEntityId{20U},
             mesh::PartitionRank{1U}, mesh::LocalIndex{0U},
             mesh::PartitionRank{2U}, mesh::LocalIndex{0U}}}};
        expect_throw<std::invalid_argument>([&] {
            (void)mesh::SharedEntityPlan::create(partition, unrelated);
        });
    }

    const auto plan = mesh::SharedEntityPlan::create(rank0_partition, links);
    require(!plan.contains_neighbor(mesh::PartitionRank{0U}),
            "exact neighbor lookup regression");
    require(!plan.contains_neighbor(mesh::PartitionRank{2U}),
            "out-of-range rank must not alias a neighbor");
    expect_throw<std::out_of_range>([&] {
        (void)plan.send_entities_to(mesh::PartitionRank{0U});
    });
    expect_throw<std::out_of_range>([&] {
        (void)plan.receive_entities_from(mesh::PartitionRank{2U});
    });
}

void headers() {
    static_assert(
        std::is_same_v<decltype(std::declval<const mesh::CsrAdjacency&>().indices()),
                       std::span<const mesh::LocalIndex>>);
    static_assert(
        std::is_same_v<decltype(std::declval<const mesh::CsrAdjacency&>().offsets()),
                       std::span<const mesh::CsrAdjacency::Offset>>);
    require(true, "headers");
}

} // namespace

int main(int argc, char** argv) {
    try {
        require(argc == 2, "provide one named mesh core test");
        const std::string_view name{argv[1]};
        if (name == "strong_indices") { strong_indices(); }
        else if (name == "two_cell_quad") { two_cell_quad(); }
        else if (name == "invalid_csr") { invalid_csr(); }
        else if (name == "topology_snapshot") { topology_snapshot(); }
        else if (name == "topology_invalid") { topology_invalid(); }
        else if (name == "cartesian_2d_topology") { cartesian_2d_topology(); }
        else if (name == "cartesian_2d_invalid") { cartesian_2d_invalid(); }
        else if (name == "cartesian_2d_geometry") { cartesian_2d_geometry(); }
        else if (name == "cartesian_2d_geometry_invalid") { cartesian_2d_geometry_invalid(); }
        else if (name == "face_boundary_snapshot") { face_boundary_snapshot(); }
        else if (name == "face_boundary_invalid") { face_boundary_invalid(); }
        else if (name == "dense_field_snapshot") { dense_field_snapshot(); }
        else if (name == "dense_field_invalid") { dense_field_invalid(); }
        else if (name == "gmsh_4_1_import") { gmsh_4_1_import(); }
        else if (name == "gmsh_4_1_roundtrip") { gmsh_4_1_roundtrip(); }
        else if (name == "gmsh_4_1_invalid") { gmsh_4_1_invalid(); }
        else if (name == "vtu_ascii_import") { vtu_ascii_import(); }
        else if (name == "vtu_ascii_roundtrip") { vtu_ascii_roundtrip(); }
        else if (name == "vtu_ascii_invalid") { vtu_ascii_invalid(); }
        else if (name == "grdecl_import") { grdecl_import(); }
        else if (name == "grdecl_inactive_degenerate") { grdecl_inactive_degenerate(); }
        else if (name == "grdecl_invalid") { grdecl_invalid(); }
        else if (name == "active_corner_point_i_neighbor") { active_corner_point_i_neighbor(); }
        else if (name == "active_corner_point_k_neighbor") { active_corner_point_k_neighbor(); }
        else if (name == "active_corner_point_activity_mapping") { active_corner_point_activity_mapping(); }
        else if (name == "active_corner_point_invalid") { active_corner_point_invalid(); }
        else if (name == "cell_face_geometric_operator_3d") { cell_face_geometric_operator_3d(); }
        else if (name == "cell_face_geometric_operator_3d_skewed") { cell_face_geometric_operator_3d_skewed(); }
        else if (name == "cell_face_geometric_operator_3d_invalid") { cell_face_geometric_operator_3d_invalid(); }
        else if (name == "dof_layout_snapshot") { dof_layout_snapshot(); }
        else if (name == "dof_layout_invalid") { dof_layout_invalid(); }
        else if (name == "partition_serial_snapshot") { partition_serial_snapshot(); }
        else if (name == "partition_local_snapshot") { partition_local_snapshot(); }
        else if (name == "partition_invalid") { partition_invalid(); }
        else if (name == "dof_numbering_serial") { dof_numbering_serial(); }
        else if (name == "dof_numbering_local") { dof_numbering_local(); }
        else if (name == "dof_numbering_invalid") { dof_numbering_invalid(); }
        else if (name == "shared_entity_plan_pair") { shared_entity_plan_pair(); }
        else if (name == "shared_entity_plan_invalid") { shared_entity_plan_invalid(); }
        else if (name == "headers") { headers(); }
        else { throw std::invalid_argument("unknown mesh core test"); }
        std::cout << "[PASS] " << name << '\n';
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "[FAIL] " << error.what() << '\n';
        return 1;
    }
}
