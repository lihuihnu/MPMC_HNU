#include <mpmc/mesh/linear_cell_geometric_operator_3d.hpp>
#include <mpmc/mesh/gmsh_4_1_3d.hpp>
#include <mpmc/mesh/linear_cell_mesh_3d.hpp>
#include <mpmc/mesh/vtu_3d.hpp>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <iostream>
#include <limits>
#include <span>
#include <stdexcept>
#include <string>
#include <string_view>
#include <vector>

namespace mesh = mpmc::mesh;

namespace {

void require(bool condition, const char* message) {
    if (!condition) {
        throw std::runtime_error(message);
    }
}

void require_close(
    double actual,
    double expected,
    double tolerance,
    const char* message) {
    if (!std::isfinite(actual) ||
        std::abs(actual - expected) >
            tolerance) {
        throw std::runtime_error(message);
    }
}

bool same_ids(
    std::span<const mesh::GlobalEntityId> left,
    std::span<const mesh::GlobalEntityId> right) {
    return left.size() == right.size() &&
           std::equal(
               left.begin(),
               left.end(),
               right.begin(),
               right.end());
}

std::string gmsh_tetra_hexa_fixture() {
    return R"MSH($MeshFormat
4.1 0 8
$EndMeshFormat
$PhysicalNames
4
2 11 "tetra_wall"
2 12 "hex_wall"
3 21 "tetra_region"
3 22 "hex_region"
$EndPhysicalNames
$Entities
12 0 2 2
1 0 0 0 0
2 1 0 0 0
3 0 1 0 0
4 0 0 1 0
5 2 0 0 0
6 3 0 0 0
7 3 1 0 0
8 2 1 0 0
9 2 0 1 0
10 3 0 1 0
11 3 1 1 0
12 2 1 1 0
1 0 0 0 1 1 0 1 11 0
2 2 0 0 3 1 0 1 12 0
1 0 0 0 1 1 1 1 21 0
2 2 0 0 3 1 1 1 22 0
$EndEntities
$Nodes
12 12 1 12
0 1 0 1
1
0 0 0
0 2 0 1
2
1 0 0
0 3 0 1
3
0 1 0
0 4 0 1
4
0 0 1
0 5 0 1
5
2 0 0
0 6 0 1
6
3 0 0
0 7 0 1
7
3 1 0
0 8 0 1
8
2 1 0
0 9 0 1
9
2 0 1
0 10 0 1
10
3 0 1
0 11 0 1
11
3 1 1
0 12 0 1
12
2 1 1
$EndNodes
$Elements
4 4 1001 2002
2 1 2 1
1001 1 3 2
2 2 3 1
1002 5 8 7 6
3 1 4 1
2001 1 2 3 4
3 2 5 1
2002 5 6 7 8 9 10 11 12
$EndElements
)MSH";
}


std::string gmsh_reversed_hexa_fixture() {
    return R"MSH($MeshFormat
4.1 0 8
$EndMeshFormat
$Nodes
1 8 1 8
3 1 0 8
1
2
3
4
5
6
7
8
0 0 0
1 0 0
1 0 1
0 0 1
0 1 0
1 1 0
1 1 1
0 1 1
$EndNodes
$Elements
1 1 100 100
3 1 5 1
100 1 2 3 4 5 6 7 8
$EndElements
)MSH";
}

std::string vtu_tetra_hexa_fixture() {
    return R"VTU(<?xml version="1.0"?>
<VTKFile type="UnstructuredGrid" version="1.0" byte_order="LittleEndian">
  <UnstructuredGrid>
    <Piece NumberOfPoints="12" NumberOfCells="2">
      <PointData>
        <DataArray type="Float64" Name="point_marker" mpmc_unit="1" format="ascii">
          0 1 2 3 4 5 6 7 8 9 10 11
        </DataArray>
      </PointData>
      <CellData>
        <DataArray type="UInt64" Name="mpmc_global_cell_id" format="ascii">
          3001 3002
        </DataArray>
        <DataArray type="Float64" Name="PORO" mpmc_unit="1" format="ascii">
          0.1 0.2
        </DataArray>
      </CellData>
      <Points>
        <DataArray type="Float64" NumberOfComponents="3" format="ascii">
          0 0 0
          1 0 0
          0 1 0
          0 0 1
          2 0 0
          3 0 0
          3 1 0
          2 1 0
          2 0 1
          3 0 1
          3 1 1
          2 1 1
        </DataArray>
      </Points>
      <Cells>
        <DataArray type="Int64" Name="connectivity" format="ascii">
          0 1 2 3 4 5 6 7 8 9 10 11
        </DataArray>
        <DataArray type="Int64" Name="offsets" format="ascii">
          4 12
        </DataArray>
        <DataArray type="UInt8" Name="types" format="ascii">
          10 12
        </DataArray>
      </Cells>
    </Piece>
  </UnstructuredGrid>
</VTKFile>
)VTU";
}


std::string gmsh_wedge_pyramid_fixture() {
    return R"MSH($MeshFormat
4.1 0 8
$EndMeshFormat
$PhysicalNames
2
3 31 "wedge_region"
3 32 "pyramid_region"
$EndPhysicalNames
$Entities
11 0 0 2
1 0 0 0 0
2 1 0 0 0
3 0 1 0 0
4 0 0 1 0
5 1 0 1 0
6 0 1 1 0
7 2 0 0 0
8 3 0 0 0
9 3 1 0 0
10 2 1 0 0
11 2.5 0.5 1 0
1 0 0 0 1 1 1 1 31 0
2 2 0 0 3 1 1 1 32 0
$EndEntities
$Nodes
11 11 1 11
0 1 0 1
1
0 0 0
0 2 0 1
2
1 0 0
0 3 0 1
3
0 1 0
0 4 0 1
4
0 0 1
0 5 0 1
5
1 0 1
0 6 0 1
6
0 1 1
0 7 0 1
7
2 0 0
0 8 0 1
8
3 0 0
0 9 0 1
9
3 1 0
0 10 0 1
10
2 1 0
0 11 0 1
11
2.5 0.5 1
$EndNodes
$Elements
2 2 4001 4002
3 1 6 1
4001 1 2 3 4 5 6
3 2 7 1
4002 7 8 9 10 11
$EndElements
)MSH";
}

std::string vtu_wedge_pyramid_fixture() {
    return R"VTU(<?xml version="1.0"?>
<VTKFile type="UnstructuredGrid" version="1.0" byte_order="LittleEndian">
  <UnstructuredGrid>
    <Piece NumberOfPoints="11" NumberOfCells="2">
      <PointData>
        <DataArray type="Float64" Name="point_marker" mpmc_unit="1" format="ascii">
          0 1 2 3 4 5 6 7 8 9 10
        </DataArray>
      </PointData>
      <CellData>
        <DataArray type="UInt64" Name="mpmc_global_cell_id" format="ascii">
          5001 5002
        </DataArray>
        <DataArray type="Float64" Name="PORO" mpmc_unit="1" format="ascii">
          0.3 0.4
        </DataArray>
      </CellData>
      <Points>
        <DataArray type="Float64" NumberOfComponents="3" format="ascii">
          0 0 0
          1 0 0
          0 1 0
          0 0 1
          1 0 1
          0 1 1
          2 0 0
          3 0 0
          3 1 0
          2 1 0
          2.5 0.5 1
        </DataArray>
      </Points>
      <Cells>
        <DataArray type="Int64" Name="connectivity" format="ascii">
          0 1 2 3 4 5 6 7 8 9 10
        </DataArray>
        <DataArray type="Int64" Name="offsets" format="ascii">
          6 11
        </DataArray>
        <DataArray type="UInt8" Name="types" format="ascii">
          13 14
        </DataArray>
      </Cells>
    </Piece>
  </UnstructuredGrid>
</VTKFile>
)VTU";
}

void verify_shared_tetra_face() {
    const std::vector<mesh::GlobalEntityId>
        vertex_ids{
            mesh::GlobalEntityId{1U},
            mesh::GlobalEntityId{2U},
            mesh::GlobalEntityId{3U},
            mesh::GlobalEntityId{4U},
            mesh::GlobalEntityId{5U}};
    const std::vector<mesh::Coordinate3D>
        coordinates{
            {0.0, 0.0, 0.0},
            {1.0, 0.0, 0.0},
            {0.0, 1.0, 0.0},
            {0.0, 0.0, 1.0},
            {0.0, 0.0, -1.0}};
    const std::vector<mesh::LinearCell3D>
        cells{
            {mesh::GlobalEntityId{10U},
             mesh::LinearCellType3D::tetrahedron,
             {mesh::LocalIndex{0U},
              mesh::LocalIndex{1U},
              mesh::LocalIndex{2U},
              mesh::LocalIndex{3U}}},
            {mesh::GlobalEntityId{11U},
             mesh::LinearCellType3D::tetrahedron,
             {mesh::LocalIndex{0U},
              mesh::LocalIndex{2U},
              mesh::LocalIndex{1U},
              mesh::LocalIndex{4U}}}};

    const auto mesh3d =
        mesh::make_linear_mesh_3d(
            vertex_ids,
            coordinates,
            cells);
    require(
        mesh3d.topology.entity_count(
            mesh::EntityKind::face) == 7U,
        "two tetrahedra must share exactly one face");
    require(
        mesh3d.face_boundary.interior_face_count() == 1U &&
            mesh3d.face_boundary.boundary_face_count() == 6U,
        "tetra shared-face classification");
    require_close(
        mesh3d.cell_volumes_m3[0],
        1.0 / 6.0,
        1.0e-14,
        "first tetra volume");
    require_close(
        mesh3d.cell_volumes_m3[1],
        1.0 / 6.0,
        1.0e-14,
        "second tetra volume");
}

void verify_shared_hexa_face() {
    const std::vector<mesh::GlobalEntityId>
        vertex_ids{
            mesh::GlobalEntityId{1U},
            mesh::GlobalEntityId{2U},
            mesh::GlobalEntityId{3U},
            mesh::GlobalEntityId{4U},
            mesh::GlobalEntityId{5U},
            mesh::GlobalEntityId{6U},
            mesh::GlobalEntityId{7U},
            mesh::GlobalEntityId{8U},
            mesh::GlobalEntityId{9U},
            mesh::GlobalEntityId{10U},
            mesh::GlobalEntityId{11U},
            mesh::GlobalEntityId{12U}};
    const std::vector<mesh::Coordinate3D>
        coordinates{
            {0.0, 0.0, 0.0},
            {1.0, 0.0, 0.0},
            {1.0, 1.0, 0.0},
            {0.0, 1.0, 0.0},
            {0.0, 0.0, 1.0},
            {1.0, 0.0, 1.0},
            {1.0, 1.0, 1.0},
            {0.0, 1.0, 1.0},
            {2.0, 0.0, 0.0},
            {2.0, 1.0, 0.0},
            {2.0, 0.0, 1.0},
            {2.0, 1.0, 1.0}};
    const std::vector<mesh::LinearCell3D>
        cells{
            {mesh::GlobalEntityId{20U},
             mesh::LinearCellType3D::hexahedron,
             {mesh::LocalIndex{0U},
              mesh::LocalIndex{1U},
              mesh::LocalIndex{2U},
              mesh::LocalIndex{3U},
              mesh::LocalIndex{4U},
              mesh::LocalIndex{5U},
              mesh::LocalIndex{6U},
              mesh::LocalIndex{7U}}},
            {mesh::GlobalEntityId{21U},
             mesh::LinearCellType3D::hexahedron,
             {mesh::LocalIndex{1U},
              mesh::LocalIndex{8U},
              mesh::LocalIndex{9U},
              mesh::LocalIndex{2U},
              mesh::LocalIndex{5U},
              mesh::LocalIndex{10U},
              mesh::LocalIndex{11U},
              mesh::LocalIndex{6U}}}};

    const auto mesh3d =
        mesh::make_linear_mesh_3d(
            vertex_ids,
            coordinates,
            cells);
    require(
        mesh3d.topology.entity_count(
            mesh::EntityKind::face) == 11U,
        "two hexahedra must share exactly one face");
    require(
        mesh3d.face_boundary.interior_face_count() == 1U &&
            mesh3d.face_boundary.boundary_face_count() == 10U,
        "hexa shared-face classification");
    require_close(
        mesh3d.cell_volumes_m3[0],
        1.0,
        1.0e-14,
        "first hexa volume");
    require_close(
        mesh3d.cell_volumes_m3[1],
        1.0,
        1.0e-14,
        "second hexa volume");

    const auto& face_cells =
        mesh3d.topology.relation(
            mesh::EntityKind::face,
            mesh::EntityKind::cell);
    bool found_internal = false;
    for (std::size_t face = 0U;
         face <
         mesh3d.topology.entity_count(
             mesh::EntityKind::face);
         ++face) {
        const auto local =
            mesh::LocalIndex{
                static_cast<
                    mesh::LocalIndex::value_type>(
                    face)};
        if (face_cells.adjacent(local).size() ==
            2U) {
            found_internal = true;
            require_close(
                mesh3d.face_geometry
                    .face_area_m2(local),
                1.0,
                1.0e-14,
                "shared hexa face area");
            const auto normal =
                mesh3d.face_geometry
                    .face_owner_unit_normal(
                        local);
            require_close(
                normal.x,
                1.0,
                1.0e-14,
                "shared hexa owner normal");
            require_close(
                normal.y,
                0.0,
                1.0e-14,
                "shared hexa owner normal y");
            require_close(
                normal.z,
                0.0,
                1.0e-14,
                "shared hexa owner normal z");
        }
    }
    require(
        found_internal,
        "hexa shared face must be materialized");

    // Exact compatibility with the original eight-vertex reference convention.
    const auto legacy = mesh::make_cell_face_geometric_operator_3d(
        mesh3d.topology, mesh3d.vertex_coordinates_m, mesh3d.face_geometry);
    const auto generic = mesh::make_cell_face_geometric_operator_3d(mesh3d);
    for (std::size_t face = 0; face < legacy.face_count(); ++face) {
        const auto local = mesh::LocalIndex{
            static_cast<mesh::LocalIndex::value_type>(face)};
        require(legacy.owner_normal_distance_m(local) ==
                    generic.owner_normal_distance_m(local) &&
                legacy.neighbour_normal_distance_m(local) ==
                    generic.neighbour_normal_distance_m(local),
                "linear hexahedron bridge must preserve legacy distances exactly");
    }
}

// Synthetic, conforming analytic fixture: a unit cube capped by a square
// pyramid, and a right triangular prism capped by a tetrahedron. It exercises
// both mixed quadrilateral and triangular interfaces and all four cell types.
// Analytic volumes and reference points below are independent of the builder.
struct MixedGeometryFixture {
    std::vector<mesh::GlobalEntityId> ids;
    std::vector<mesh::Coordinate3D> points;
    std::vector<mesh::LinearCell3D> cells;
};

MixedGeometryFixture mixed_geometry_fixture() {
    using I = mesh::LocalIndex;
    using G = mesh::GlobalEntityId;
    using T = mesh::LinearCellType3D;
    MixedGeometryFixture fixture{
        {},
        {{0,0,0}, {1,0,0}, {1,1,0}, {0,1,0},
         {0,0,1}, {1,0,1}, {1,1,1}, {0,1,1}, {.5,.5,2},
         {3,0,0}, {4,0,0}, {3,1,0},
         {3,0,1}, {4,0,1}, {3,1,1}, {3,0,2}},
        {{G{100}, T::hexahedron, {I{0},I{1},I{2},I{3},I{4},I{5},I{6},I{7}}},
         {G{200}, T::pyramid, {I{4},I{5},I{6},I{7},I{8}}},
         {G{300}, T::wedge, {I{9},I{10},I{11},I{12},I{13},I{14}}},
         {G{400}, T::tetrahedron, {I{12},I{13},I{14},I{15}}}}};
    for (std::size_t i = 0; i < fixture.points.size(); ++i) {
        fixture.ids.push_back(G{101U + 17U * static_cast<std::uint64_t>(i)});
    }
    return fixture;
}

void verify_mixed_connection_geometry(const MixedGeometryFixture& fixture) {
    const auto grid = mesh::make_linear_mesh_3d(
        fixture.ids, fixture.points, fixture.cells);
    const auto op = mesh::make_cell_face_geometric_operator_3d(grid);
    require(grid.face_boundary.interior_face_count() == 2U,
            "mixed mesh has one quad and one triangular interface");

    const std::array<mesh::Coordinate3D, 4> reference_points{{
        {.5,.5,.5}, {.5,.5,1.2}, {10.0/3.0,1.0/3.0,.5}, {3.25,.25,1.25}}};
    constexpr std::array<double, 4> volumes{1.0, 1.0/3.0, .5, 1.0/6.0};
    for (std::size_t cell = 0; cell < fixture.cells.size(); ++cell) {
        const auto local = mesh::LocalIndex{
            static_cast<mesh::LocalIndex::value_type>(cell)};
        const auto oracle = static_cast<std::size_t>(
            fixture.cells[cell].global_id.value() / 100U - 1U);
        const auto point = op.cell_reference_point_m(local);
        require_close(point.x_m, reference_points[oracle].x_m, 1e-14,
                      "mixed cell reference x and stable identity");
        require_close(point.y_m, reference_points[oracle].y_m, 1e-14,
                      "mixed cell reference y and stable identity");
        require_close(point.z_m, reference_points[oracle].z_m, 1e-14,
                      "mixed cell reference z and stable identity");
        require_close(grid.cell_volumes_m3[cell], volumes[oracle], 1e-14,
                      "mixed cell analytic volume");
        require(op.cell_centroid_m(local).z_m == point.z_m,
                "legacy accessor remains a reference-point alias");

        std::array<double, 3> closure{};
        std::array<std::array<double, 3>, 3> moment{};
        for (const auto face : grid.topology.relation(
                 mesh::EntityKind::cell, mesh::EntityKind::face).adjacent(local)) {
            const auto connection = op.face_connection_geometry(face);
            const bool is_owner = connection.owner == local;
            const double sign = is_owner ? 1.0 : -1.0;
            const std::array<double, 3> normal{
                sign * connection.owner_unit_normal.x,
                sign * connection.owner_unit_normal.y,
                sign * connection.owner_unit_normal.z};
            const auto displacement = is_owner ?
                op.owner_to_face_displacement_m(face) :
                *op.neighbour_to_face_displacement_m(face);
            const std::array<double, 3> d{
                displacement.x_m, displacement.y_m, displacement.z_m};
            for (std::size_t i = 0; i < 3U; ++i) {
                closure[i] += connection.area_m2 * normal[i];
                for (std::size_t j = 0; j < 3U; ++j) {
                    moment[i][j] += connection.area_m2 * normal[i] * d[j];
                }
            }
            require(connection.owner_normal_distance_m > 0.0,
                    "all mixed owner normal distances positive");
            if (connection.neighbour.has_value()) {
                require(*connection.neighbour_normal_distance_m > 0.0,
                        "all mixed neighbour normal distances positive");
            }
        }
        // Divergence theorem: zero flux for a constant vector field; affine
        // surface moments integral n_i (x_j-reference_j) dA = V delta_ij.
        // Face-centroid quadrature is exact for these planar linear fields.
        for (std::size_t i = 0; i < 3U; ++i) {
            require_close(closure[i], 0.0, 1e-13, "constant-field closure");
            for (std::size_t j = 0; j < 3U; ++j) {
                require_close(moment[i][j], i == j ? volumes[oracle] : 0.0,
                              1e-13, "analytic affine surface moment");
            }
        }
    }
    for (std::size_t face = 0; face < op.face_count(); ++face) {
        const auto local = mesh::LocalIndex{
            static_cast<mesh::LocalIndex::value_type>(face)};
        const auto neighbour = op.face_neighbour(local);
        if (!neighbour) {
            require(!op.internal_non_orthogonality(local),
                    "boundary faces have no internal non-orthogonality");
            continue;
        }
        const auto degree = grid.topology.relation(
            mesh::EntityKind::face, mesh::EntityKind::vertex).adjacent(local).size();
        const auto measure = *op.internal_non_orthogonality(local);
        require_close(op.owner_normal_distance_m(local) +
                          *op.neighbour_normal_distance_m(local),
                      degree == 4U ? .7 : .75, 1e-14,
                      "mixed interface two-sided analytic normal distance");
        require_close(measure.center_distance_m,
                      degree == 4U ? .7 : std::sqrt(83.0)/12.0, 1e-14,
                      "mixed interface analytic centre separation");
        require_close(measure.normal_alignment_cosine,
                      degree == 4U ? 1.0 : 9.0/std::sqrt(83.0), 1e-14,
                      "mixed interface analytic non-orthogonality");
    }
}

void verify_linear_geometry_bridge() {
    auto fixture = mixed_geometry_fixture();
    verify_mixed_connection_geometry(fixture);
    // Reverse local vertex/cell ordering, while keeping element node ordering
    // and stable identities. This also reverses interface owner/neighbour.
    std::reverse(fixture.ids.begin(), fixture.ids.end());
    std::reverse(fixture.points.begin(), fixture.points.end());
    std::reverse(fixture.cells.begin(), fixture.cells.end());
    for (auto& cell : fixture.cells) {
        for (auto& vertex : cell.vertices) {
            vertex = mesh::LocalIndex{static_cast<mesh::LocalIndex::value_type>(
                fixture.points.size() - 1U - vertex.value())};
        }
    }
    verify_mixed_connection_geometry(fixture);
}

void verify_explicit_reference_points() {
    const auto fixture = mixed_geometry_fixture();
    const auto grid = mesh::make_linear_mesh_3d(
        fixture.ids, fixture.points, fixture.cells);
    std::vector<mesh::Coordinate3D> points{
        {.5,.5,.5}, {.5,.5,1.25}, {10.0/3.0,1.0/3.0,.5}, {3.25,.25,1.25}};
    // For a height-one pyramid, the uniform-volume centroid is 1/4 of the
    // height above the base, whereas its five-vertex mean is 1/5. The explicit
    // API must honor the selected point instead of silently averaging vertices.
    const auto op = mesh::make_cell_face_geometric_operator_3d_from_reference_points(
        grid.topology, points, grid.face_geometry);
    require_close(op.cell_reference_point_m(mesh::LocalIndex{1}).z_m,
                  1.25, 0.0, "explicit pyramid volume-centroid reference");
    bool found_quad = false;
    for (std::size_t face = 0; face < op.face_count(); ++face) {
        const auto local = mesh::LocalIndex{
            static_cast<mesh::LocalIndex::value_type>(face)};
        if (op.face_neighbour(local) == mesh::LocalIndex{1}) {
            found_quad = true;
            require_close(*op.neighbour_normal_distance_m(local), .25, 1e-14,
                          "explicit pyramid reference changes base distance");
        }
    }
    require(found_quad, "pyramid/cube explicit interface exists");

    const auto rejects = [&](const std::vector<mesh::Coordinate3D>& invalid) {
        bool rejected = false;
        try {
            (void)mesh::make_cell_face_geometric_operator_3d_from_reference_points(
                grid.topology, invalid, grid.face_geometry);
        } catch (const std::invalid_argument&) {
            rejected = true;
        }
        require(rejected, "invalid explicit reference points must be rejected");
    };
    auto invalid = points;
    invalid.pop_back();
    rejects(invalid);
    invalid = points;
    invalid[0].x_m = std::numeric_limits<double>::quiet_NaN();
    rejects(invalid);
    invalid = points;
    invalid[0].x_m = std::numeric_limits<double>::infinity();
    rejects(invalid);
    invalid = points;
    invalid[1].z_m = 1.0; // On the common face: zero neighbour distance.
    rejects(invalid);
    invalid[1].z_m = .9; // Wrong side of the common face.
    rejects(invalid);
    invalid = points;
    invalid[0].z_m = 1.1; // Wrong side for the owner.
    rejects(invalid);
    bool legacy_rejected = false;
    try {
        (void)mesh::make_cell_face_geometric_operator_3d(
            grid.topology, grid.vertex_coordinates_m, grid.face_geometry);
    } catch (const std::invalid_argument&) {
        legacy_rejected = true;
    }
    require(legacy_rejected, "legacy factory still requires eight-vertex cells");
}


void verify_wedge_pyramid_geometry() {
    const auto imported =
        mesh::import_vtu_ascii_3d(
            vtu_wedge_pyramid_fixture());
    require(
        imported.topology.entity_count(
            mesh::EntityKind::cell) == 2U &&
            imported.topology.entity_count(
                mesh::EntityKind::face) == 10U &&
            imported.topology.entity_count(
                mesh::EntityKind::vertex) == 11U,
        "wedge/pyramid entity counts");
    require_close(
        imported.cell_volumes_m3[0],
        0.5,
        1.0e-14,
        "wedge volume");
    require_close(
        imported.cell_volumes_m3[1],
        1.0 / 3.0,
        1.0e-14,
        "pyramid volume");

    const auto& face_vertices =
        imported.topology.relation(
            mesh::EntityKind::face,
            mesh::EntityKind::vertex);
    std::size_t triangle_faces = 0U;
    std::size_t quad_faces = 0U;
    for (std::size_t face = 0U;
         face <
         imported.topology.entity_count(
             mesh::EntityKind::face);
         ++face) {
        const auto degree =
            face_vertices.adjacent(
                mesh::LocalIndex{
                    static_cast<
                        mesh::LocalIndex::value_type>(
                        face)})
                .size();
        if (degree == 3U) {
            ++triangle_faces;
        } else if (degree == 4U) {
            ++quad_faces;
        }
    }
    require(
        triangle_faces == 6U &&
            quad_faces == 4U,
        "wedge/pyramid triangle/quad face families");
}

void verify_gmsh_wedge_pyramid_roundtrip() {
    const auto first =
        mesh::import_gmsh_4_1_ascii_3d(
            gmsh_wedge_pyramid_fixture(),
            1.0);
    require(
        first.topology.entity_count(
            mesh::EntityKind::cell) == 2U &&
            first.physical_names.size() == 2U &&
            first.cell_physical_groups.size() == 2U,
        "Gmsh wedge/pyramid physical metadata");
    require_close(
        first.cell_volumes_m3[0],
        0.5,
        1.0e-14,
        "Gmsh wedge volume");
    require_close(
        first.cell_volumes_m3[1],
        1.0 / 3.0,
        1.0e-14,
        "Gmsh pyramid volume");

    const auto text =
        mesh::export_gmsh_4_1_ascii_3d(
            first);
    const auto second =
        mesh::import_gmsh_4_1_ascii_3d(
            text,
            1.0);
    require(
        same_ids(
            second.topology.global_ids(
                mesh::EntityKind::face),
            first.topology.global_ids(
                mesh::EntityKind::face)) &&
            same_ids(
                second.topology.global_ids(
                    mesh::EntityKind::cell),
                first.topology.global_ids(
                    mesh::EntityKind::cell)),
        "Gmsh wedge/pyramid stable IDs roundtrip");
    require(
        second.physical_names.size() ==
                first.physical_names.size() &&
            second.cell_physical_groups.size() ==
                first.cell_physical_groups.size(),
        "Gmsh wedge/pyramid metadata roundtrip");
}

void verify_vtu_wedge_pyramid_roundtrip() {
    const auto first =
        mesh::import_vtu_ascii_3d(
            vtu_wedge_pyramid_fixture());
    require(
        first.point_fields.size() == 1U &&
            first.cell_fields.size() == 1U,
        "VTU wedge/pyramid fields imported");
    require_close(
        first.cell_fields[0].value(
            mesh::LocalIndex{0U},
            0U),
        0.3,
        1.0e-14,
        "VTU wedge PORO");
    require_close(
        first.cell_fields[0].value(
            mesh::LocalIndex{1U},
            0U),
        0.4,
        1.0e-14,
        "VTU pyramid PORO");

    const auto text =
        mesh::export_vtu_ascii_3d(
            first);
    const auto second =
        mesh::import_vtu_ascii_3d(
            text);
    require(
        same_ids(
            second.topology.global_ids(
                mesh::EntityKind::cell),
            first.topology.global_ids(
                mesh::EntityKind::cell)) &&
            same_ids(
                second.topology.global_ids(
                    mesh::EntityKind::face),
                first.topology.global_ids(
                    mesh::EntityKind::face)),
        "VTU wedge/pyramid stable IDs roundtrip");
    require_close(
        second.cell_volumes_m3[0],
        0.5,
        1.0e-14,
        "VTU wedge roundtrip volume");
    require_close(
        second.cell_volumes_m3[1],
        1.0 / 3.0,
        1.0e-14,
        "VTU pyramid roundtrip volume");
}


void verify_gmsh_reversed_orientation_canonicalization() {
    const auto first =
        mesh::import_gmsh_4_1_ascii_3d(
            gmsh_reversed_hexa_fixture(),
            1.0);
    require(
        first.topology.entity_count(
            mesh::EntityKind::cell) == 1U &&
            first.topology.entity_count(
                mesh::EntityKind::face) == 6U &&
            first.topology.entity_count(
                mesh::EntityKind::vertex) == 8U,
        "reversed Gmsh hexa entity counts");
    require_close(
        first.cell_volumes_m3[0],
        1.0,
        1.0e-14,
        "reversed Gmsh hexa canonical volume");

    const auto exported =
        mesh::export_gmsh_4_1_ascii_3d(first);
    const auto second =
        mesh::import_gmsh_4_1_ascii_3d(
            exported,
            1.0);
    require_close(
        second.cell_volumes_m3[0],
        1.0,
        1.0e-14,
        "canonicalized Gmsh hexa roundtrip volume");
    require(
        same_ids(
            first.topology.global_ids(
                mesh::EntityKind::cell),
            second.topology.global_ids(
                mesh::EntityKind::cell)),
        "canonicalized Gmsh hexa stable cell ID");
}

void verify_gmsh_roundtrip() {
    const auto first =
        mesh::import_gmsh_4_1_ascii_3d(
            gmsh_tetra_hexa_fixture(),
            1.0);
    require(
        first.topology.entity_count(
            mesh::EntityKind::cell) == 2U &&
            first.topology.entity_count(
                mesh::EntityKind::face) == 10U &&
            first.topology.entity_count(
                mesh::EntityKind::vertex) == 12U,
        "Gmsh tetra/hexa entity counts");
    require_close(
        first.cell_volumes_m3[0],
        1.0 / 6.0,
        1.0e-14,
        "Gmsh tetra volume");
    require_close(
        first.cell_volumes_m3[1],
        1.0,
        1.0e-14,
        "Gmsh hexa volume");
    require(
        first.face_boundary.boundary_face_count() == 10U &&
            first.physical_names.size() == 4U &&
            first.cell_physical_groups.size() == 2U,
        "Gmsh physical metadata");

    const auto text =
        mesh::export_gmsh_4_1_ascii_3d(
            first);
    const auto second =
        mesh::import_gmsh_4_1_ascii_3d(
            text,
            1.0);
    require(
        same_ids(
            second.topology.global_ids(
                mesh::EntityKind::vertex),
            first.topology.global_ids(
                mesh::EntityKind::vertex)),
        "Gmsh vertex stable IDs roundtrip");
    require(
        same_ids(
            second.topology.global_ids(
                mesh::EntityKind::face),
            first.topology.global_ids(
                mesh::EntityKind::face)),
        "Gmsh face stable IDs roundtrip");
    require(
        same_ids(
            second.topology.global_ids(
                mesh::EntityKind::cell),
            first.topology.global_ids(
                mesh::EntityKind::cell)),
        "Gmsh cell stable IDs roundtrip");
    require(
        second.physical_names.size() ==
            first.physical_names.size() &&
            second.cell_physical_groups.size() ==
                first.cell_physical_groups.size(),
        "Gmsh physical metadata roundtrip");
    require_close(
        second.cell_volumes_m3[0],
        first.cell_volumes_m3[0],
        1.0e-14,
        "Gmsh tetra roundtrip volume");
    require_close(
        second.cell_volumes_m3[1],
        first.cell_volumes_m3[1],
        1.0e-14,
        "Gmsh hexa roundtrip volume");
}

void verify_vtu_roundtrip() {
    const auto first =
        mesh::import_vtu_ascii_3d(
            vtu_tetra_hexa_fixture());
    require(
        first.topology.entity_count(
            mesh::EntityKind::cell) == 2U &&
            first.topology.entity_count(
                mesh::EntityKind::face) == 10U &&
            first.topology.entity_count(
                mesh::EntityKind::vertex) == 12U,
        "VTU tetra/hexa entity counts");
    require(
        first.point_fields.size() == 1U &&
            first.cell_fields.size() == 1U,
        "VTU 3D point/cell fields imported");
    require(
        first.cell_fields[0].metadata().id ==
                "PORO" &&
            first.cell_fields[0].metadata().unit ==
                "1",
        "VTU 3D cell field metadata");
    require_close(
        first.cell_fields[0].value(
            mesh::LocalIndex{0U},
            0U),
        0.1,
        1.0e-14,
        "VTU PORO first cell");

    const auto text =
        mesh::export_vtu_ascii_3d(
            first);
    const auto second =
        mesh::import_vtu_ascii_3d(
            text);
    require(
        same_ids(
            second.topology.global_ids(
                mesh::EntityKind::cell),
            first.topology.global_ids(
                mesh::EntityKind::cell)),
        "VTU stable cell IDs roundtrip");
    require(
        same_ids(
            second.topology.global_ids(
                mesh::EntityKind::face),
            first.topology.global_ids(
                mesh::EntityKind::face)),
        "VTU deterministic face IDs roundtrip");
    require(
        second.point_fields.size() == 1U &&
            second.cell_fields.size() == 1U,
        "VTU fields roundtrip");
    require_close(
        second.cell_volumes_m3[0],
        1.0 / 6.0,
        1.0e-14,
        "VTU tetra volume roundtrip");
    require_close(
        second.cell_volumes_m3[1],
        1.0,
        1.0e-14,
        "VTU hexa volume roundtrip");
    require_close(
        second.cell_fields[0].value(
            mesh::LocalIndex{1U},
            0U),
        0.2,
        1.0e-14,
        "VTU PORO second cell roundtrip");
}

void verify_invalid_cases() {
    bool inverted_rejected = false;
    try {
        const std::vector<mesh::GlobalEntityId>
            vertex_ids{
                mesh::GlobalEntityId{1U},
                mesh::GlobalEntityId{2U},
                mesh::GlobalEntityId{3U},
                mesh::GlobalEntityId{4U}};
        const std::vector<mesh::Coordinate3D>
            coordinates{
                {0.0, 0.0, 0.0},
                {1.0, 0.0, 0.0},
                {0.0, 1.0, 0.0},
                {0.0, 0.0, 1.0}};
        const std::vector<mesh::LinearCell3D>
            cells{
                {mesh::GlobalEntityId{1U},
                 mesh::LinearCellType3D::tetrahedron,
                 {mesh::LocalIndex{0U},
                  mesh::LocalIndex{2U},
                  mesh::LocalIndex{1U},
                  mesh::LocalIndex{3U}}}};
        (void)mesh::make_linear_mesh_3d(
            vertex_ids,
            coordinates,
            cells);
    } catch (const std::invalid_argument&) {
        inverted_rejected = true;
    }
    require(
        inverted_rejected,
        "inverted tetrahedron must be rejected");

    std::string unsupported =
        vtu_tetra_hexa_fixture();
    const auto position =
        unsupported.find(
            "10 12");
    require(
        position != std::string::npos,
        "VTU invalid fixture anchor");
    unsupported.replace(
        position,
        5U,
        "15 12");
    bool unsupported_cell_rejected = false;
    try {
        (void)mesh::import_vtu_ascii_3d(
            unsupported);
    } catch (const std::invalid_argument&) {
        unsupported_cell_rejected = true;
    }
    require(
        unsupported_cell_rejected,
        "unsupported VTK 3D cell type must be rejected");

    const std::vector<mesh::GlobalEntityId>
        wedge_ids{
            mesh::GlobalEntityId{1U},
            mesh::GlobalEntityId{2U},
            mesh::GlobalEntityId{3U},
            mesh::GlobalEntityId{4U},
            mesh::GlobalEntityId{5U},
            mesh::GlobalEntityId{6U}};
    const std::vector<mesh::Coordinate3D>
        wedge_coordinates{
            {0.0, 0.0, 0.0},
            {1.0, 0.0, 0.0},
            {0.0, 1.0, 0.0},
            {0.0, 0.0, 1.0},
            {1.0, 0.0, 1.0},
            {0.0, 1.0, 1.0}};
    const std::vector<mesh::LinearCell3D>
        inverted_wedge{
            {mesh::GlobalEntityId{1U},
             mesh::LinearCellType3D::wedge,
             {mesh::LocalIndex{0U},
              mesh::LocalIndex{2U},
              mesh::LocalIndex{1U},
              mesh::LocalIndex{3U},
              mesh::LocalIndex{5U},
              mesh::LocalIndex{4U}}}};
    bool inverted_wedge_rejected = false;
    try {
        (void)mesh::make_linear_mesh_3d(
            wedge_ids,
            wedge_coordinates,
            inverted_wedge);
    } catch (const std::invalid_argument&) {
        inverted_wedge_rejected = true;
    }
    require(
        inverted_wedge_rejected,
        "inverted wedge must be rejected");

    const std::vector<mesh::GlobalEntityId>
        pyramid_ids{
            mesh::GlobalEntityId{1U},
            mesh::GlobalEntityId{2U},
            mesh::GlobalEntityId{3U},
            mesh::GlobalEntityId{4U},
            mesh::GlobalEntityId{5U}};
    const std::vector<mesh::Coordinate3D>
        pyramid_coordinates{
            {0.0, 0.0, 0.0},
            {1.0, 0.0, 0.0},
            {1.0, 1.0, 0.0},
            {0.0, 1.0, 0.0},
            {0.5, 0.5, 1.0}};
    const std::vector<mesh::LinearCell3D>
        inverted_pyramid{
            {mesh::GlobalEntityId{1U},
             mesh::LinearCellType3D::pyramid,
             {mesh::LocalIndex{0U},
              mesh::LocalIndex{3U},
              mesh::LocalIndex{2U},
              mesh::LocalIndex{1U},
              mesh::LocalIndex{4U}}}};
    bool inverted_pyramid_rejected = false;
    try {
        (void)mesh::make_linear_mesh_3d(
            pyramid_ids,
            pyramid_coordinates,
            inverted_pyramid);
    } catch (const std::invalid_argument&) {
        inverted_pyramid_rejected = true;
    }
    require(
        inverted_pyramid_rejected,
        "inverted pyramid must be rejected");
}

} // namespace

int main() {
    try {
        verify_shared_tetra_face();
        verify_shared_hexa_face();
        verify_linear_geometry_bridge();
        verify_explicit_reference_points();
        verify_wedge_pyramid_geometry();
        verify_gmsh_reversed_orientation_canonicalization();
        verify_gmsh_roundtrip();
        verify_gmsh_wedge_pyramid_roundtrip();
        verify_vtu_roundtrip();
        verify_vtu_wedge_pyramid_roundtrip();
        verify_invalid_cases();
        std::cout
            << "[PASS] mesh.core.io_3d\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr
            << "[FAIL] mesh.core.io_3d: "
            << error.what()
            << '\n';
        return 1;
    }
}
