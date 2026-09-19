#include <mpmc/mesh/gmsh_4_1_3d.hpp>
#include <mpmc/mesh/linear_cell_mesh_3d.hpp>
#include <mpmc/mesh/vtu_3d.hpp>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <iostream>
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
        "13 12");
    bool wedge_rejected = false;
    try {
        (void)mesh::import_vtu_ascii_3d(
            unsupported);
    } catch (const std::invalid_argument&) {
        wedge_rejected = true;
    }
    require(
        wedge_rejected,
        "unsupported VTK wedge must be rejected");
}

} // namespace

int main() {
    try {
        verify_shared_tetra_face();
        verify_shared_hexa_face();
        verify_gmsh_roundtrip();
        verify_vtu_roundtrip();
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
