#include <mpmc/mesh/mesh_exchange_io.hpp>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <iostream>
#include <limits>
#include <map>
#include <optional>
#include <span>
#include <stdexcept>
#include <string>
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

// Validation must distinguish entity kinds even when stable IDs overlap, and
// retain the first diagnostic when an input violates more than one rule.
void verify_group_validation_lookup() {
    using Kind = mesh::EntityKind;
    using Id = mesh::GlobalEntityId;
    const auto create = [](std::vector<mesh::MeshExchangeGroup> groups) {
        mesh::Topology::EntityIds ids;
        ids.vertices = {Id{0U}, Id{std::numeric_limits<std::uint64_t>::max()}};
        ids.edges = {Id{900U}};
        ids.faces = {Id{900U}};
        ids.cells = {Id{3U}};
        return mesh::MeshExchangeDocument::create(
            mesh::MeshExchangeFormat::gmsh_4_1_ascii, 2,
            mesh::Topology{std::move(ids), {}}, {{0.0, 0.0, 0.0}, {1.0, 0.0, 0.0}},
            std::nullopt, {}, std::move(groups), std::nullopt);
    };
    const auto valid = create({
        {Kind::vertex, 1U, "vertices", {Id{std::numeric_limits<std::uint64_t>::max()}, Id{0U}, Id{0U}}},
        {Kind::vertex, 2U, "first", {Id{0U}}},
        {Kind::edge, 1U, "edge", {Id{900U}}},
        {Kind::face, 1U, "face", {Id{900U}}},
        {Kind::cell, 1U, "cell", {Id{3U}}}});
    require(valid.groups().size() == 5U && valid.groups()[0].members.size() == 3U,
            "validation must preserve group order and repeated members");
    require(create({}).groups().empty(), "empty groups must remain valid");
    const auto rejects = [&](std::vector<mesh::MeshExchangeGroup> groups,
                             std::string_view expected) {
        try {
            (void)create(std::move(groups));
        } catch (const std::invalid_argument& error) {
            return std::string_view(error.what()) == expected;
        }
        return false;
    };
    constexpr std::string_view missing =
        "mpmc::mesh::MeshExchangeDocument: group member is absent from topology";
    require(rejects({{Kind::vertex, 1U, "wrong kind", {Id{900U}}}}, missing),
            "an ID from another entity kind must be rejected");
    require(rejects({{Kind::cell, 1U, "valid", {Id{3U}}},
                     {Kind::cell, 2U, "missing", {Id{42U}}}}, missing),
            "reused index must still reject absent members");
    require(rejects({{Kind::cell, 0U, "zero", {Id{42U}}}},
            "mpmc::mesh::MeshExchangeDocument: group tag zero is reserved"),
            "zero tag must precede missing-member diagnostic");
    require(rejects({{Kind::cell, 1U, std::string("bad\0name", 8U), {Id{42U}}}},
            "mpmc::mesh::MeshExchangeDocument: group name cannot contain NUL"),
            "invalid name must precede missing-member diagnostic");
    require(rejects({{Kind::cell, 1U, "first", {}}, {Kind::cell, 1U, "duplicate", {Id{42U}}}},
            "mpmc::mesh::MeshExchangeDocument: duplicate group key"),
            "duplicate key must precede missing-member diagnostic");
    require(rejects({{static_cast<Kind>(99), 1U, "invalid kind", {}}},
            "mpmc::mesh::Topology: invalid entity kind"),
            "invalid kind must be rejected before indexing the lookup array");
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

// Identity regression oracle: compare the IDs attached to actual point/face
// entities after export/import, never reproduce the analyzer's numbering rule.
using FaceIdentities = std::map<std::vector<std::uint32_t>, std::uint64_t>;

FaceIdentities face_identities(const mesh::Topology& topology) {
    FaceIdentities result;
    const auto ids = topology.global_ids(mesh::EntityKind::face);
    const auto& relation = topology.relation(mesh::EntityKind::face, mesh::EntityKind::vertex);
    for (std::size_t i = 0; i < ids.size(); ++i) {
        std::vector<std::uint32_t> vertices;
        for (const auto vertex : relation.adjacent(mesh::LocalIndex{static_cast<std::uint32_t>(i)})) {
            vertices.push_back(vertex.value());
        }
        std::sort(vertices.begin(), vertices.end());
        require(result.emplace(std::move(vertices), ids[i].value()).second, "duplicate face in identity fixture");
    }
    return result;
}

bool same_ids(const mesh::Topology& a, const mesh::Topology& b, mesh::EntityKind kind) {
    const auto left = a.global_ids(kind);
    const auto right = b.global_ids(kind);
    return std::equal(left.begin(), left.end(), right.begin(), right.end());
}

void verify_identity_report(const mesh::MeshExchangeDocument& document,
                            bool vertices_changed, bool faces_changed) {
    const auto exported = mesh::export_vtu_ascii(document);
    require(exported.exported(), "identity fixture must export");
    const auto actual = document.dimension() == 2 ?
        mesh::make_mesh_exchange_document(mesh::import_vtu_ascii(*exported.content)) :
        mesh::make_mesh_exchange_document(mesh::import_vtu_ascii_3d(*exported.content));
    require(!same_ids(document.topology(), actual.topology(), mesh::EntityKind::vertex) == vertices_changed,
            "fixture vertex identity change differs from explicit expectation");
    require((face_identities(document.topology()) != face_identities(actual.topology())) == faces_changed,
            "fixture face identity change differs from explicit expectation");
    require(same_ids(document.topology(), actual.topology(), mesh::EntityKind::cell),
            "VTU must preserve exact UInt64 cell IDs");
    const auto before = document.vertex_coordinates_m();
    const auto after = actual.vertex_coordinates_m();
    require(before.size() == after.size(), "identity change must not drop points");
    for (std::size_t i = 0; i < before.size(); ++i) {
        require(before[i].x_m == after[i].x_m && before[i].y_m == after[i].y_m &&
                before[i].z_m == after[i].z_m, "identity change must not alter coordinates");
    }
    const auto& left = document.topology().relation(mesh::EntityKind::cell, mesh::EntityKind::vertex);
    const auto& right = actual.topology().relation(mesh::EntityKind::cell, mesh::EntityKind::vertex);
    for (std::size_t i = 0; i < document.topology().entity_count(mesh::EntityKind::cell); ++i) {
        const auto a = left.adjacent(mesh::LocalIndex{static_cast<std::uint32_t>(i)});
        const auto b = right.adjacent(mesh::LocalIndex{static_cast<std::uint32_t>(i)});
        require(std::equal(a.begin(), a.end(), b.begin(), b.end()), "cyclic cell connectivity changed");
    }
    for (const auto& report : {mesh::analyze_conversion(document, mesh::MeshExchangeFormat::vtu_ascii),
                               exported.report}) {
        bool vertex_issue = false;
        bool face_issue = false;
        for (const auto& issue : report.issues()) {
            if (issue.code == "vtu.vertex_ids_remapped") vertex_issue = true;
            else if (issue.code == "vtu.face_ids_remapped") face_issue = true;
            else throw std::runtime_error("unrelated loss masks identity regression: " + issue.code);
        }
        require(vertex_issue == vertices_changed, "VTU report must match actual vertex identity loss");
        require(face_issue == faces_changed, "VTU report must match actual face identity loss");
        require(report.issues().size() == static_cast<std::size_t>(vertices_changed) +
                                            static_cast<std::size_t>(faces_changed), "duplicate identity issue");
        require(report.lossless() == (!vertices_changed && !faces_changed),
                "VTU lossless promise differs from actual entity identities");
    }
}

mesh::MeshExchangeDocument identity_variant(const mesh::MeshExchangeDocument& original,
    std::vector<mesh::GlobalEntityId> vertices, std::vector<mesh::GlobalEntityId> faces,
    bool reverse_face_storage = false) {
    // Reorder storage and all affected relations together, preserving bindings.
    const auto& topology = original.topology();
    const auto face_count = faces.size();
    if (reverse_face_storage) std::reverse(faces.begin(), faces.end());
    mesh::Topology::EntityIds ids;
    ids.vertices = std::move(vertices);
    ids.faces = std::move(faces);
    const auto cell_count = topology.entity_count(mesh::EntityKind::cell);
    for (std::size_t i = 0; i < cell_count; ++i) {
        ids.cells.emplace_back(9007199254740993ULL + static_cast<std::uint64_t>(i));
    }
    constexpr std::array kinds{mesh::EntityKind::vertex, mesh::EntityKind::face, mesh::EntityKind::cell};
    std::vector<mesh::CsrAdjacency> relations;
    for (const auto from : kinds) for (const auto to : kinds) {
        if (!topology.has_relation(from, to)) continue;
        std::vector<mesh::CsrAdjacency::Offset> offsets{0};
        std::vector<mesh::LocalIndex> indices;
        for (std::size_t i = 0; i < topology.entity_count(from); ++i) {
            const auto source = reverse_face_storage && from == mesh::EntityKind::face ? face_count-1U-i : i;
            for (const auto index : topology.relation(from, to).adjacent(
                     mesh::LocalIndex{static_cast<std::uint32_t>(source)})) {
                const auto target = reverse_face_storage && to == mesh::EntityKind::face ?
                    face_count-1U-index.value() : index.value();
                indices.emplace_back(static_cast<std::uint32_t>(target));
            }
            offsets.push_back(static_cast<mesh::CsrAdjacency::Offset>(indices.size()));
        }
        relations.emplace_back(from, to, topology.entity_count(to), std::move(offsets), std::move(indices));
    }
    return mesh::MeshExchangeDocument::create(mesh::MeshExchangeFormat::gmsh_4_1_ascii,
        original.dimension(), mesh::Topology{std::move(ids), std::move(relations)},
        {original.vertex_coordinates_m().begin(), original.vertex_coordinates_m().end()},
        std::nullopt, {}, {}, std::nullopt);
}

// Same textual contract for both public importers, independent of writer internals.
void verify_vertex_id_contract(const mesh::MeshExchangeDocument& source) {
    const auto parse = [&](const std::string& text) {
        return source.dimension() == 2 ?
            mesh::make_mesh_exchange_document(mesh::import_vtu_ascii(text)) :
            mesh::make_mesh_exchange_document(mesh::import_vtu_ascii_3d(text));
    };
    const auto xml = *mesh::export_vtu_ascii(source).content;
    const auto name = xml.find("Name=\"mpmc_global_vertex_id\"");
    require(name != std::string::npos, "writer must emit vertex identity");
    const auto begin = xml.rfind("<DataArray", name);
    const auto end = xml.find("</DataArray>", name) + std::string_view{"</DataArray>"}.size();
    const auto make = [&](const std::string& attributes, const std::string& values) {
        return "<DataArray Name=\"mpmc_global_vertex_id\" " + attributes + ">" + values + "</DataArray>";
    };
    const auto replace = [&](const std::string& array) {
        auto text = xml;
        text.replace(begin, end-begin, array);
        return text;
    };
    const auto missing = parse(replace(""));
    for (std::size_t i = 0; i < 5U; ++i) {
        require(missing.topology().global_ids(mesh::EntityKind::vertex)[i].value() == i+1U,
                "legacy VTU missing identity must retain 1..N convention");
    }
    for (const auto& [type, payload] : std::array<std::pair<std::string, std::string>, 2>{{
            {"UInt64", "0 18446744073709551615 9007199254740993 12 7"},
            {"Int64", "0 9223372036854775807 9007199254740993 12 7"}}}) {
        const auto imported = parse(replace(make("type=\""+type+"\" NumberOfComponents=\"1\" NumberOfTuples=\"5\" format=\"ascii\"", payload)));
        const auto ids = imported.topology().global_ids(mesh::EntityKind::vertex);
        require(ids[0].value() == 0U && ids[2].value() == 9007199254740993ULL &&
                ids[3].value() == 12U && ids[4].value() == 7U, "exact vertex identity values");
        require(ids[1].value() == (type == "UInt64" ? std::numeric_limits<std::uint64_t>::max() :
                static_cast<std::uint64_t>(std::numeric_limits<std::int64_t>::max())), "full-width vertex identity");
        require(imported.fields().fields().empty(), "identity must not become a floating field");
        verify_identity_report(imported, false, false);
    }
    const std::string valid = make("type=\"UInt64\" format=\"ascii\"", "5 4 3 2 1");
    const std::vector<std::string> invalid{
        make("type=\"Float64\" format=\"ascii\"", "5 4 3 2 1"),
        make("type=\"UInt32\" format=\"ascii\"", "5 4 3 2 1"),
        make("type=\"UInt64\" format=\"ascii\"", "1 1 3 4 5"),
        make("type=\"UInt64\" format=\"ascii\"", "1 2 3 4"),
        make("type=\"UInt64\" format=\"ascii\"", "1 2 3 4 5 6"),
        make("type=\"UInt64\" format=\"ascii\"", "-1 2 3 4 5"),
        make("type=\"Int64\" format=\"ascii\"", "-1 2 3 4 5"),
        make("type=\"UInt64\" format=\"ascii\"", "18446744073709551616 2 3 4 5"),
        make("type=\"Int64\" format=\"ascii\"", "9223372036854775808 2 3 4 5"),
        make("type=\"UInt64\" format=\"ascii\"", "1.0 2 3 4 5"),
        make("type=\"UInt64\" NumberOfComponents=\"2\" format=\"ascii\"", "1 2 3 4 5"),
        make("type=\"UInt64\" NumberOfComponents=\"0\" format=\"ascii\"", "1 2 3 4 5"),
        make("type=\"UInt64\" NumberOfTuples=\"4\" format=\"ascii\"", "1 2 3 4 5"),
        make("type=\"UInt64\" format=\"binary\"", "1 2 3 4 5"), valid+valid};
    for (const auto& array : invalid) {
        bool rejected = false;
        try { (void)parse(replace(array)); } catch (const std::invalid_argument&) { rejected = true; }
        require(rejected, "invalid vertex identity array must be rejected");
    }
    auto misplaced = replace("");
    misplaced.insert(misplaced.find("</CellData>"), valid);
    bool rejected = false;
    try { (void)parse(misplaced); } catch (const std::invalid_argument&) { rejected = true; }
    require(rejected, "vertex identity is reserved to PointData");

    const auto collisions = [](auto imported, auto writer) {
        for (const auto location : {mesh::EntityKind::vertex, mesh::EntityKind::cell}) {
            for (const auto* reserved : {"mpmc_global_vertex_id", "mpmc_global_cell_id"}) {
                imported.point_fields.clear(); imported.cell_fields.clear();
                auto value = mesh::DenseFieldSnapshot::create(imported.topology, location, 1U,
                    std::vector<double>(imported.topology.entity_count(location), 1.0),
                    {reserved, "1", {mesh::FieldSourceKind::synthetic_test, "identity-contract", "v1", ""}});
                (location == mesh::EntityKind::vertex ? imported.point_fields : imported.cell_fields).push_back(value);
                bool failed = false;
                try { (void)writer(imported); } catch (const std::invalid_argument&) { failed = true; }
                require(failed, "scientific fields cannot shadow reserved identity names");
            }
        }
    };
    if (source.dimension() == 2) {
        collisions(mesh::import_vtu_ascii(xml), [](const auto& input) { return mesh::export_vtu_ascii(input); });
    } else {
        collisions(mesh::import_vtu_ascii_3d(xml), [](const auto& input) { return mesh::export_vtu_ascii_3d(input); });
    }
    std::cout << "[PASS] mesh.core.exchange_io.vtu_vertex_ids dimension=" << source.dimension()
              << " invalid_arrays=" << invalid.size()+1U << " reserved_collisions=4\n";
}

void verify_vtu_identity_losses() {
    // Handwritten geometry; mixed 2D has a shared edge, pyramid mixes triangle
    // and quad faces. Default IDs are checked against actual readers below.
    const auto two = mesh::make_mesh_exchange_document(mesh::import_vtu_ascii(R"VTU(
<VTKFile type="UnstructuredGrid"><UnstructuredGrid><Piece NumberOfPoints="5" NumberOfCells="2">
<Points><DataArray type="Float64" NumberOfComponents="3" format="ascii">0 0 0 3 0 0 2 2 0 0 2 0 4 0 0</DataArray></Points>
<Cells><DataArray type="Int64" Name="connectivity" format="ascii">1 4 2 3 2 1 0</DataArray>
<DataArray type="Int64" Name="offsets" format="ascii">3 7</DataArray>
<DataArray type="UInt8" Name="types" format="ascii">5 9</DataArray></Cells>
</Piece></UnstructuredGrid></VTKFile>)VTU"));
    const auto three = mesh::make_mesh_exchange_document(mesh::import_vtu_ascii_3d(R"VTU(
<VTKFile type="UnstructuredGrid"><UnstructuredGrid><Piece NumberOfPoints="5" NumberOfCells="1">
<Points><DataArray type="Float64" NumberOfComponents="3" format="ascii">0 0 0 1 0 0 1 1 0 0 1 0 0.5 0.5 1</DataArray></Points>
<Cells><DataArray type="Int64" Name="connectivity" format="ascii">0 1 2 3 4</DataArray>
<DataArray type="Int64" Name="offsets" format="ascii">5</DataArray>
<DataArray type="UInt8" Name="types" format="ascii">14</DataArray></Cells>
</Piece></UnstructuredGrid></VTKFile>)VTU"));
    for (const auto* source : {&two, &three}) {
        verify_vertex_id_contract(*source);
        const auto v = source->topology().global_ids(mesh::EntityKind::vertex);
        const auto f = source->topology().global_ids(mesh::EntityKind::face);
        const std::vector<mesh::GlobalEntityId> vertices{v.begin(), v.end()}, faces{f.begin(), f.end()};
        verify_identity_report(identity_variant(*source, vertices, faces), false, false);
        verify_identity_report(identity_variant(*source, vertices, faces, true), false, false);
        auto sparse_vertices = vertices;
        sparse_vertices[0] = mesh::GlobalEntityId{0};
        sparse_vertices[1] = mesh::GlobalEntityId{std::numeric_limits<std::uint64_t>::max()};
        verify_identity_report(identity_variant(*source, sparse_vertices, faces), false, source->dimension() == 2);
        auto permuted_vertices = vertices;
        std::swap(permuted_vertices[0], permuted_vertices[1]);
        verify_identity_report(identity_variant(*source, permuted_vertices, faces), false, source->dimension() == 2);
        auto sparse_faces = faces;
        sparse_faces[0] = mesh::GlobalEntityId{9007199254740999ULL};
        verify_identity_report(identity_variant(*source, vertices, sparse_faces), false, true);
        auto permuted_faces = faces;
        std::swap(permuted_faces[0], permuted_faces[1]);
        verify_identity_report(identity_variant(*source, vertices, permuted_faces), false, true);
        verify_identity_report(identity_variant(*source, sparse_vertices, sparse_faces), false, true);
        verify_identity_report(identity_variant(*source, sparse_vertices, sparse_faces, true), false, true);
    }
    std::cout << "[PASS] mesh.core.exchange_io.vtu_identity_losses cases=16\n";
}

} // namespace

int main() {
    try {
        verify_vtu_identity_losses();
        verify_grdecl_canonical_roundtrip();
        verify_gmsh_group_bridge();
        verify_sparse_group_lookup();
        verify_group_validation_lookup();
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
