#include <mpmc/mesh/mesh_exchange_io.hpp>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <iostream>
#include <span>
#include <stdexcept>
#include <string_view>

namespace mesh = mpmc::mesh;

namespace {

void require(
    bool condition,
    const char* message) {
    if (!condition) {
        throw std::runtime_error(message);
    }
}

bool same_values(
    std::span<const double> left,
    std::span<const double> right,
    double tolerance) {
    if (left.size() != right.size()) {
        return false;
    }
    for (std::size_t index = 0U;
         index < left.size();
         ++index) {
        if (!std::isfinite(left[index]) ||
            !std::isfinite(right[index]) ||
            std::abs(
                left[index] -
                right[index]) >
                tolerance) {
            return false;
        }
    }
    return true;
}

const mesh::DenseFieldSnapshot&
field(
    const mesh::GrdeclImportResult& result,
    std::string_view id) {
    const auto found =
        std::find_if(
            result.cell_fields.begin(),
            result.cell_fields.end(),
            [&](const auto& candidate) {
                return candidate.metadata().id ==
                    id;
            });
    if (found == result.cell_fields.end()) {
        throw std::runtime_error(
            "field missing");
    }
    return *found;
}

std::string_view gmsh_group_fixture() {
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

std::string_view grdecl_fixture() {
    return R"GRDECL(SPECGRID
1 1 1 1 F /
COORD
0 0 0 0 0 1
1 0 0 1 0 1
0 1 0 0 1 1
1 1 0 1 1 1
/
ZCORN
0 0 0 0 1 1 1 1
/
ACTNUM
1
/
PORO
0.25
/
PERMX
10
/
PERMY
20
/
PERMZ
30
/
)GRDECL";
}

void verify_grdecl_canonical_roundtrip() {
    const mesh::GrdeclImportOptions options{
        2.0,
        1.0e-15};
    const auto first =
        mesh::import_grdecl(
            grdecl_fixture(),
            options);
    require(
        first.coord_m.size() == 24U &&
            first.zcorn_m.size() == 8U,
        "GRDECL source semantics retained");

    const auto document =
        mesh::make_mesh_exchange_document(
            first);
    require(
        document.source_format() ==
                mesh::MeshExchangeFormat::grdecl &&
            document.dimension() == 3 &&
            document.logical_corner_point()
                .has_value(),
        "GRDECL canonical document");

    const auto report =
        mesh::analyze_conversion(
            document,
            mesh::MeshExchangeFormat::grdecl);
    require(
        report.lossless() &&
            report.issues().empty(),
        "GRDECL canonical export must be lossless");

    const auto exported =
        mesh::export_grdecl_ascii(
            document);
    require(
        exported.exported() &&
            exported.report.lossless(),
        "native GRDECL export");

    const auto second =
        mesh::import_grdecl(
            *exported.content,
            options);
    require(
        first.dimensions ==
            second.dimensions &&
            first.active ==
                second.active,
        "GRDECL logical identity roundtrip");
    require(
        same_values(
            first.coord_m,
            second.coord_m,
            1.0e-14) &&
            same_values(
                first.zcorn_m,
                second.zcorn_m,
                1.0e-14),
        "GRDECL COORD/ZCORN roundtrip");

    for (const auto id :
         {"PORO", "PERMX", "PERMY", "PERMZ"}) {
        require(
            same_values(
                field(first, id).values(),
                field(second, id).values(),
                1.0e-28),
            "GRDECL property roundtrip");
    }

    const auto as_vtu =
        mesh::export_vtu_ascii(
            document);
    require(
        as_vtu.exported() &&
            as_vtu.report.disposition() ==
                mesh::ConversionDisposition::lossy,
        "GRDECL canonical VTU export");
    const auto vtu_second =
        mesh::import_vtu_ascii_3d(
            *as_vtu.content);
    require(
        vtu_second.topology.entity_count(
            mesh::EntityKind::cell) == 1U &&
            vtu_second.cell_fields.size() == 4U,
        "GRDECL canonical VTU re-import");

    const auto as_gmsh =
        mesh::export_gmsh_4_1_ascii(
            document);
    require(
        as_gmsh.exported() &&
            as_gmsh.report.disposition() ==
                mesh::ConversionDisposition::lossy,
        "GRDECL canonical Gmsh export");
    const auto gmsh_second =
        mesh::import_gmsh_4_1_ascii_3d(
            *as_gmsh.content,
            1.0);
    require(
        gmsh_second.topology.entity_count(
            mesh::EntityKind::cell) == 1U,
        "GRDECL canonical Gmsh re-import");
}

void verify_gmsh_group_bridge() {
    const auto first =
        mesh::import_gmsh_4_1_ascii(
            gmsh_group_fixture(),
            1.0);
    const auto document =
        mesh::make_mesh_exchange_document(
            first);
    const auto exported =
        mesh::export_gmsh_4_1_ascii(
            document);
    require(
        exported.exported() &&
            exported.report.lossless(),
        "Gmsh group bridge must be lossless");
    const auto second =
        mesh::import_gmsh_4_1_ascii(
            *exported.content,
            1.0);
    require(
        second.physical_names.size() ==
            first.physical_names.size() &&
            second.cell_physical_groups.size() ==
                first.cell_physical_groups.size(),
        "Gmsh physical group bridge");
    const auto first_tags =
        first.face_boundary.physical_tags();
    const auto second_tags =
        second.face_boundary.physical_tags();
    require(
        first_tags.size() ==
                second_tags.size() &&
            std::equal(
                first_tags.begin(),
                first_tags.end(),
                second_tags.begin()),
        "Gmsh boundary tags bridge");
}

void verify_non_corner_point_report() {
    const auto source =
        mesh::import_vtu_ascii_3d(
            R"VTU(<?xml version="1.0"?>
<VTKFile type="UnstructuredGrid" version="1.0" byte_order="LittleEndian">
  <UnstructuredGrid>
    <Piece NumberOfPoints="4" NumberOfCells="1">
      <PointData/>
      <CellData/>
      <Points>
        <DataArray type="Float64" NumberOfComponents="3" format="ascii">
          0 0 0 1 0 0 0 1 0 0 0 1
        </DataArray>
      </Points>
      <Cells>
        <DataArray type="Int64" Name="connectivity" format="ascii">0 1 2 3</DataArray>
        <DataArray type="Int64" Name="offsets" format="ascii">4</DataArray>
        <DataArray type="UInt8" Name="types" format="ascii">10</DataArray>
      </Cells>
    </Piece>
  </UnstructuredGrid>
</VTKFile>)VTU");
    const auto document =
        mesh::make_mesh_exchange_document(
            source);
    const auto vtu_export =
        mesh::export_vtu_ascii(
            document);
    require(
        vtu_export.exported() &&
            vtu_export.report.lossless(),
        "VTU canonical VTU bridge");
    const auto vtu_roundtrip =
        mesh::import_vtu_ascii_3d(
            *vtu_export.content);
    require(
        vtu_roundtrip.topology.entity_count(
            mesh::EntityKind::cell) == 1U,
        "VTU canonical VTU re-import");

    const auto gmsh_export =
        mesh::export_gmsh_4_1_ascii(
            document);
    require(
        gmsh_export.exported() &&
            gmsh_export.report.disposition() !=
                mesh::ConversionDisposition::
                    unsupported,
        "VTU canonical Gmsh bridge");
    const auto gmsh_roundtrip =
        mesh::import_gmsh_4_1_ascii_3d(
            *gmsh_export.content,
            1.0);
    require(
        gmsh_roundtrip.topology.entity_count(
            mesh::EntityKind::cell) == 1U,
        "VTU canonical Gmsh re-import");

    const auto report =
        mesh::analyze_conversion(
            document,
            mesh::MeshExchangeFormat::grdecl);
    require(
        report.disposition() ==
                mesh::ConversionDisposition::
                    unsupported &&
            !report.issues().empty(),
        "arbitrary VTU must not silently convert to GRDECL");

    const auto exported =
        mesh::export_grdecl_ascii(
            document,
            mesh::GrdeclExportOptions{
                1.0,
                1.0});
    require(
        !exported.exported() &&
            exported.report.disposition() ==
                mesh::ConversionDisposition::
                    unsupported,
        "unsupported GRDECL export must not produce content");
}

} // namespace

int main() {
    try {
        verify_grdecl_canonical_roundtrip();
        verify_gmsh_group_bridge();
        verify_non_corner_point_report();
        std::cout
            << "[PASS] mesh.core.exchange_io\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr
            << "[FAIL] mesh.core.exchange_io: "
            << error.what()
            << '\n';
        return 1;
    }
}
