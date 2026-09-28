#include <mpmc/mesh/mesh_exchange_io.hpp>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <iostream>
#include <limits>
#include <optional>
#include <span>
#include <stdexcept>
#include <string_view>
#include <utility>
#include <vector>

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

// Stable IDs need not be sorted, contiguous, or fit a local index. The writer
// must map group members to the original local order before canonical renaming.
void verify_sparse_group_lookup() {
    const std::array ids{
        mesh::GlobalEntityId{900U}, mesh::GlobalEntityId{0U},
        mesh::GlobalEntityId{std::numeric_limits<std::uint64_t>::max()},
        mesh::GlobalEntityId{3U}};
    mesh::mesh_exchange_io_detail::CanonicalEntityLookup lookup(ids);
    for (int repeat = 0; repeat < 2; ++repeat) {
        for (std::size_t local = 0U; local < ids.size(); ++local) {
            require(lookup.local(ids[local]) == local,
                    "sparse global ID must retain its original local index");
        }
    }
    const auto rejects_missing = [](auto& index) {
        try {
            (void)index.local(mesh::GlobalEntityId{42U});
        } catch (const std::invalid_argument& error) {
            return std::string_view(error.what()) ==
                "mpmc::mesh::canonical writer: group member is absent from topology";
        }
        return false;
    };
    require(rejects_missing(lookup), "missing group member must be rejected");
    mesh::mesh_exchange_io_detail::CanonicalEntityLookup empty({});
    require(rejects_missing(empty), "empty topology must reject a group member");

    const auto raw = mesh::import_gmsh_4_1_ascii(gmsh_group_fixture(), 1.0);
    const auto original = mesh::make_mesh_exchange_document(raw);
    mesh::Topology::EntityIds remapped;
    const auto copy_ids = [&](mesh::EntityKind kind) {
        const auto source = original.topology().global_ids(kind);
        return std::vector<mesh::GlobalEntityId>(source.begin(), source.end());
    };
    remapped.vertices = copy_ids(mesh::EntityKind::vertex);
    remapped.edges = copy_ids(mesh::EntityKind::edge);
    remapped.faces = copy_ids(mesh::EntityKind::face);
    remapped.cells = {ids[0], ids[3]};
    for (std::size_t local = 0U; local < remapped.faces.size(); ++local) {
        remapped.faces[local] = mesh::GlobalEntityId{
            1000U + static_cast<std::uint64_t>(remapped.faces.size() - local) * 7U};
    }
    const std::vector<mesh::MeshExchangeGroup> groups{
        {mesh::EntityKind::face, 11U, "boundary", {remapped.faces[0], remapped.faces[2]}},
        {mesh::EntityKind::cell, 21U, "second", {remapped.cells[1]}},
        {mesh::EntityKind::cell, 22U, "both", {remapped.cells[1], remapped.cells[0]}}};
    const auto coordinates = original.vertex_coordinates_m();
    const auto document = mesh::MeshExchangeDocument::create(
        original.source_format(), original.dimension(),
        mesh::Topology{std::move(remapped),
            mesh::mesh_exchange_io_detail::copy_relations(original.topology())},
        {coordinates.begin(), coordinates.end()}, std::nullopt, {}, groups, std::nullopt);
    const auto boundary = mesh::mesh_exchange_io_detail::canonical_face_boundary(
        document, original.topology());
    const auto tags = boundary.physical_tags();
    for (std::size_t local = 0U; local < tags.size(); ++local) {
        require(tags[local].value() == (local == 0U || local == 2U ? 11U : 0U),
                "face groups must follow local order after canonical renaming");
    }
    const auto cells = mesh::mesh_exchange_io_detail::canonical_cell_groups(
        document, original.topology());
    require(cells.size() == 2U &&
            cells[0].physical_tags == std::vector<std::uint32_t>{22U} &&
            cells[1].physical_tags == std::vector<std::uint32_t>{21U, 22U},
            "cell groups must retain memberships and sorted tags after renaming");
}

void verify_rectilinear_grdecl_reconstruction() {
    const auto source =
        mesh::import_vtu_ascii_3d(
            R"VTU(<?xml version="1.0"?>
<VTKFile type="UnstructuredGrid" version="1.0" byte_order="LittleEndian">
  <UnstructuredGrid>
    <Piece NumberOfPoints="12" NumberOfCells="2">
      <PointData/>
      <CellData>
        <DataArray type="Float64" Name="PORO" mpmc_unit="1" format="ascii">0.2 0.1</DataArray>
      </CellData>
      <Points>
        <DataArray type="Float64" NumberOfComponents="3" format="ascii">
          0 0 0  1 0 0  3 0 0
          0 2 0  1 2 0  3 2 0
          0 0 1  1 0 1  3 0 1
          0 2 1  1 2 1  3 2 1
        </DataArray>
      </Points>
      <Cells>
        <DataArray type="Int64" Name="connectivity" format="ascii">
          1 2 5 4 7 8 11 10
          0 1 4 3 6 7 10 9
        </DataArray>
        <DataArray type="Int64" Name="offsets" format="ascii">8 16</DataArray>
        <DataArray type="UInt8" Name="types" format="ascii">12 12</DataArray>
      </Cells>
    </Piece>
  </UnstructuredGrid>
</VTKFile>)VTU");
    const auto document =
        mesh::make_mesh_exchange_document(
            source);

    const auto reconstruction =
        mesh::reconstruct_structured_logical_grid_3d(
            document);
    require(
        reconstruction.report.representable() &&
            reconstruction.report.dimensions()
                .has_value() &&
            *reconstruction.report.dimensions() ==
                std::array<std::size_t, 3>{
                    2U, 1U, 1U} &&
            reconstruction.logical_grid
                .has_value() &&
            reconstruction.source_cell_by_logical
                .size() == 2U &&
            reconstruction.source_cell_by_logical[0]
                    .value() == 1U &&
            reconstruction.source_cell_by_logical[1]
                    .value() == 0U,
        "reversed source cells must reconstruct as I-fastest 2x1x1 GRDECL");
    require(
        reconstruction.logical_grid->
                porosity.size() == 2U &&
            reconstruction.logical_grid->
                porosity[0] == 0.1 &&
            reconstruction.logical_grid->
                porosity[1] == 0.2,
        "reconstructed PORO must follow logical cell ordering");

    const auto exported =
        mesh::export_grdecl_ascii(
            document);
    require(
        exported.exported() &&
            exported.report.lossless(),
        "rectilinear canonical mesh must export exact reconstructed GRDECL");

    const auto imported =
        mesh::import_grdecl(
            *exported.content,
            mesh::GrdeclImportOptions{
                1.0,
                1.0});
    require(
        imported.dimensions ==
                std::array<std::size_t, 3>{
                    2U, 1U, 1U} &&
            imported.cell_fields.size() == 1U &&
            imported.active_cell_count() == 2U,
        "reconstructed GRDECL import");
    require(
        field(imported, "PORO")
                .value(mesh::LocalIndex{0U}, 0U) ==
                0.1 &&
            field(imported, "PORO")
                .value(mesh::LocalIndex{1U}, 0U) ==
                0.2,
        "re-imported PORO logical ordering");
    const auto processed =
        mesh::process_active_corner_point_grid(
            imported);
    require(
        processed.cell_count() == 2U &&
            processed.cell_fields.size() == 1U,
        "reconstructed GRDECL active processing");
}

void verify_slanted_hexa_reconstruction_rejection() {
    const auto source =
        mesh::import_vtu_ascii_3d(
            R"VTU(<?xml version="1.0"?>
<VTKFile type="UnstructuredGrid" version="1.0" byte_order="LittleEndian">
  <UnstructuredGrid>
    <Piece NumberOfPoints="8" NumberOfCells="1">
      <PointData/>
      <CellData/>
      <Points>
        <DataArray type="Float64" NumberOfComponents="3" format="ascii">
          0 0 0 1 0 0 1 1 0 0 1 0
          0.1 0 1 1.1 0 1 1.1 1 1 0.1 1 1
        </DataArray>
      </Points>
      <Cells>
        <DataArray type="Int64" Name="connectivity" format="ascii">0 1 2 3 4 5 6 7</DataArray>
        <DataArray type="Int64" Name="offsets" format="ascii">8</DataArray>
        <DataArray type="UInt8" Name="types" format="ascii">12</DataArray>
      </Cells>
    </Piece>
  </UnstructuredGrid>
</VTKFile>)VTU");
    const auto document =
        mesh::make_mesh_exchange_document(
            source);
    const auto reconstruction =
        mesh::reconstruct_structured_logical_grid_3d(
            document);
    require(
        !reconstruction.report.representable() &&
            !reconstruction.report.issues().empty(),
        "slanted hexa is outside rectilinear GRDECL reconstruction baseline");

    const auto exported =
        mesh::export_grdecl_ascii(
            document);
    require(
        !exported.exported() &&
            exported.report.disposition() ==
                mesh::ConversionDisposition::
                    unsupported,
        "slanted hexa GRDECL export must remain unsupported");
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
        verify_sparse_group_lookup();
        verify_rectilinear_grdecl_reconstruction();
        verify_slanted_hexa_reconstruction_rejection();
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
