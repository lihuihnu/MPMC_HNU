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
