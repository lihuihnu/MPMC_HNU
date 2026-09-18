#include <mpmc/mesh/cartesian_2d.hpp>
#include <mpmc/mesh/csr_adjacency.hpp>
#include <mpmc/mesh/entity.hpp>
#include <mpmc/mesh/geometry_2d.hpp>
#include <mpmc/mesh/topology.hpp>

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
        else if (name == "headers") { headers(); }
        else { throw std::invalid_argument("unknown mesh core test"); }
        std::cout << "[PASS] " << name << '\n';
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "[FAIL] " << error.what() << '\n';
        return 1;
    }
}
