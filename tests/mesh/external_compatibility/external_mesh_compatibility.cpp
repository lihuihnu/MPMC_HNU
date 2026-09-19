#include <mpmc/mesh/active_corner_point.hpp>
#include <mpmc/mesh/gmsh_4_1.hpp>
#include <mpmc/mesh/gmsh_4_1_3d.hpp>
#include <mpmc/mesh/grdecl.hpp>
#include <mpmc/mesh/mesh_exchange_io.hpp>
#include <mpmc/mesh/vtu.hpp>
#include <mpmc/mesh/vtu_3d.hpp>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <fstream>
#include <iostream>
#include <span>
#include <stdexcept>
#include <string>
#include <string_view>
#include <vector>

namespace mesh = mpmc::mesh;

namespace {

[[nodiscard]] std::string read_file(
    const char* path) {
    std::ifstream input(
        path,
        std::ios::binary);
    if (!input) {
        throw std::runtime_error(
            std::string{
                "cannot open external sample: "} +
            path);
    }
    return std::string{
        std::istreambuf_iterator<char>{input},
        std::istreambuf_iterator<char>{}};
}

void require(
    bool condition,
    std::string_view message) {
    if (!condition) {
        throw std::runtime_error(
            std::string{message});
    }
}

void require_close(
    double actual,
    double expected,
    double tolerance,
    std::string_view message) {
    if (!std::isfinite(actual) ||
        !std::isfinite(expected) ||
        std::abs(actual - expected) >
            tolerance) {
        throw std::runtime_error(
            std::string{message});
    }
}

[[nodiscard]] std::vector<std::string>
split_vtu_documents(
    std::string_view bundle) {
    std::vector<std::string> documents;
    std::size_t search_from = 0U;
    constexpr std::string_view begin_tag =
        "<VTKFile";
    constexpr std::string_view end_tag =
        "</VTKFile>";

    while (true) {
        const auto begin =
            bundle.find(
                begin_tag,
                search_from);
        if (begin ==
            std::string_view::npos) {
            break;
        }
        const auto end =
            bundle.find(
                end_tag,
                begin);
        if (end ==
            std::string_view::npos) {
            throw std::runtime_error(
                "unterminated VTKFile document in external bundle");
        }
        const std::size_t final =
            end + end_tag.size();
        documents.emplace_back(
            bundle.substr(
                begin,
                final - begin));
        search_from = final;
    }
    return documents;
}

void verify_gmsh_2d(
    std::string_view content) {
    const auto source =
        mesh::import_gmsh_4_1_ascii(
            content,
            1.0);
    require(
        source.topology.entity_count(
            mesh::EntityKind::vertex) == 4U &&
            source.topology.entity_count(
                mesh::EntityKind::face) == 4U &&
            source.topology.entity_count(
                mesh::EntityKind::cell) == 1U,
        "deal.II Gmsh 2D entity counts");

    const auto canonical =
        mesh::make_mesh_exchange_document(
            source);

    const auto as_gmsh =
        mesh::export_gmsh_4_1_ascii(
            canonical);
    require(
        as_gmsh.exported() &&
            as_gmsh.report.lossless(),
        "Gmsh 2D canonical->Gmsh lossless");
    const auto gmsh_second =
        mesh::import_gmsh_4_1_ascii(
            *as_gmsh.content,
            1.0);
    require(
        gmsh_second.topology.entity_count(
            mesh::EntityKind::cell) == 1U &&
            gmsh_second.physical_names.size() ==
                source.physical_names.size(),
        "Gmsh 2D canonical Gmsh re-import");

    const auto as_vtu =
        mesh::export_vtu_ascii(
            canonical);
    require(
        as_vtu.exported() &&
            as_vtu.report.disposition() ==
                mesh::ConversionDisposition::lossy,
        "Gmsh 2D canonical->VTU must report group/tag loss");
    const auto vtu_second =
        mesh::import_vtu_ascii(
            *as_vtu.content);
    require(
        vtu_second.topology.entity_count(
            mesh::EntityKind::cell) == 1U,
        "Gmsh 2D canonical VTU re-import");

    const auto as_grdecl =
        mesh::export_grdecl_ascii(
            canonical);
    require(
        !as_grdecl.exported() &&
            as_grdecl.report.disposition() ==
                mesh::ConversionDisposition::
                    unsupported,
        "Gmsh 2D canonical->GRDECL unsupported");

    std::cout
        << "[PASS] external.gmsh.2d.dealii\n"
        << "matrix.gmsh2d.gmsh=lossless\n"
        << "matrix.gmsh2d.vtu=lossy\n"
        << "matrix.gmsh2d.grdecl=unsupported\n";
}

void verify_gmsh_3d(
    std::string_view content) {
    const auto source =
        mesh::import_gmsh_4_1_ascii_3d(
            content,
            1.0);
    require(
        source.topology.entity_count(
            mesh::EntityKind::vertex) == 8U &&
            source.topology.entity_count(
                mesh::EntityKind::face) == 6U &&
            source.topology.entity_count(
                mesh::EntityKind::cell) == 1U,
        "deal.II Gmsh 3D entity counts");

    const auto canonical =
        mesh::make_mesh_exchange_document(
            source);

    const auto as_gmsh =
        mesh::export_gmsh_4_1_ascii(
            canonical);
    require(
        as_gmsh.exported() &&
            as_gmsh.report.lossless(),
        "Gmsh 3D canonical->Gmsh lossless");
    const auto gmsh_second =
        mesh::import_gmsh_4_1_ascii_3d(
            *as_gmsh.content,
            1.0);
    require_close(
        gmsh_second.cell_volumes_m3[0],
        1.0,
        1.0e-14,
        "Gmsh 3D canonical Gmsh volume");

    const auto as_vtu =
        mesh::export_vtu_ascii(
            canonical);
    require(
        as_vtu.exported() &&
            as_vtu.report.disposition() ==
                mesh::ConversionDisposition::lossy,
        "Gmsh 3D canonical->VTU must report group/tag loss");
    const auto vtu_second =
        mesh::import_vtu_ascii_3d(
            *as_vtu.content);
    require_close(
        vtu_second.cell_volumes_m3[0],
        1.0,
        1.0e-14,
        "Gmsh 3D canonical VTU volume");

    const auto as_grdecl =
        mesh::export_grdecl_ascii(
            canonical);
    require(
        !as_grdecl.exported() &&
            as_grdecl.report.disposition() ==
                mesh::ConversionDisposition::
                    unsupported,
        "Gmsh 3D canonical->GRDECL unsupported");

    std::cout
        << "[PASS] external.gmsh.3d.dealii\n"
        << "matrix.gmsh3d.gmsh=lossless\n"
        << "matrix.gmsh3d.vtu=lossy\n"
        << "matrix.gmsh3d.grdecl=unsupported\n";
}

void verify_vtu_bundle(
    std::string_view content) {
    const auto documents =
        split_vtu_documents(content);
    require(
        documents.size() == 2U,
        "deal.II VTU external bundle must contain exactly two VTKFile documents");

    const auto source_2d =
        mesh::import_vtu_ascii(
            documents[0]);
    const auto canonical_2d =
        mesh::make_mesh_exchange_document(
            source_2d);

    const auto vtu_2d =
        mesh::export_vtu_ascii(
            canonical_2d);
    require(
        vtu_2d.exported() &&
            vtu_2d.report.lossless(),
        "VTU 2D canonical->VTU lossless");
    const auto vtu_2d_second =
        mesh::import_vtu_ascii(
            *vtu_2d.content);
    require(
        vtu_2d_second.point_fields.size() ==
            source_2d.point_fields.size(),
        "VTU 2D canonical field preservation");

    const auto gmsh_2d =
        mesh::export_gmsh_4_1_ascii(
            canonical_2d);
    require(
        gmsh_2d.exported() &&
            gmsh_2d.report.disposition() ==
                mesh::ConversionDisposition::lossy,
        "VTU 2D canonical->Gmsh field loss report");
    const auto gmsh_2d_second =
        mesh::import_gmsh_4_1_ascii(
            *gmsh_2d.content,
            1.0);
    require(
        gmsh_2d_second.topology.entity_count(
            mesh::EntityKind::cell) == 1U,
        "VTU 2D canonical Gmsh re-import");

    const auto grdecl_2d =
        mesh::export_grdecl_ascii(
            canonical_2d);
    require(
        !grdecl_2d.exported() &&
            grdecl_2d.report.disposition() ==
                mesh::ConversionDisposition::
                    unsupported,
        "VTU 2D canonical->GRDECL unsupported");

    const auto source_3d =
        mesh::import_vtu_ascii_3d(
            documents[1]);
    const auto canonical_3d =
        mesh::make_mesh_exchange_document(
            source_3d);

    const auto vtu_3d =
        mesh::export_vtu_ascii(
            canonical_3d);
    require(
        vtu_3d.exported() &&
            vtu_3d.report.lossless(),
        "VTU 3D canonical->VTU lossless");
    const auto vtu_3d_second =
        mesh::import_vtu_ascii_3d(
            *vtu_3d.content);
    require(
        vtu_3d_second.point_fields.size() ==
            source_3d.point_fields.size(),
        "VTU 3D canonical field preservation");

    const auto gmsh_3d =
        mesh::export_gmsh_4_1_ascii(
            canonical_3d);
    require(
        gmsh_3d.exported() &&
            gmsh_3d.report.disposition() ==
                mesh::ConversionDisposition::lossy,
        "VTU 3D canonical->Gmsh field loss report");
    const auto gmsh_3d_second =
        mesh::import_gmsh_4_1_ascii_3d(
            *gmsh_3d.content,
            1.0);
    require_close(
        gmsh_3d_second.cell_volumes_m3[0],
        1.0,
        1.0e-14,
        "VTU 3D canonical Gmsh volume");

    const auto grdecl_3d =
        mesh::export_grdecl_ascii(
            canonical_3d);
    require(
        !grdecl_3d.exported() &&
            grdecl_3d.report.disposition() ==
                mesh::ConversionDisposition::
                    unsupported,
        "VTU 3D canonical->GRDECL unsupported");

    std::cout
        << "[PASS] external.vtu.2d3d.dealii\n"
        << "matrix.vtu2d.gmsh=lossy\n"
        << "matrix.vtu2d.vtu=lossless\n"
        << "matrix.vtu2d.grdecl=unsupported\n"
        << "matrix.vtu3d.gmsh=lossy\n"
        << "matrix.vtu3d.vtu=lossless\n"
        << "matrix.vtu3d.grdecl=unsupported\n";
}

void verify_grdecl_tube(
    std::string_view content) {
    const auto raw =
        mesh::import_grdecl(
            content,
            mesh::GrdeclImportOptions{
                1.0,
                1.0});
    require(
        raw.dimensions[0] == 3U &&
            raw.dimensions[1] == 1U &&
            raw.dimensions[2] == 2U &&
            raw.cell_count() == 6U,
        "OPM tube.grdecl dimensions");

    const auto canonical =
        mesh::make_mesh_exchange_document(
            raw);
    require(
        canonical.topology().entity_count(
            mesh::EntityKind::cell) == 6U &&
            canonical.topology().entity_count(
                mesh::EntityKind::face) != 0U &&
            canonical.logical_corner_point()
                .has_value(),
        "GRDECL canonical generic/logical dual representation");

    const auto native =
        mesh::export_grdecl_ascii(
            canonical);
    require(
        native.exported() &&
            native.report.lossless(),
        "GRDECL canonical->GRDECL lossless");
    const auto native_second =
        mesh::import_grdecl(
            *native.content,
            raw.source_options);
    require(
        native_second.dimensions ==
                raw.dimensions &&
            native_second.active ==
                raw.active,
        "GRDECL native structural roundtrip");
    for (std::size_t index = 0U;
         index < raw.coord_m.size();
         ++index) {
        require_close(
            native_second.coord_m[index],
            raw.coord_m[index],
            1.0e-12,
            "GRDECL native COORD");
    }
    for (std::size_t index = 0U;
         index < raw.zcorn_m.size();
         ++index) {
        require_close(
            native_second.zcorn_m[index],
            raw.zcorn_m[index],
            1.0e-12,
            "GRDECL native ZCORN");
    }

    const auto as_vtu =
        mesh::export_vtu_ascii(
            canonical);
    require(
        as_vtu.exported() &&
            as_vtu.report.disposition() ==
                mesh::ConversionDisposition::lossy,
        "GRDECL canonical->VTU logical semantic loss report");
    const auto vtu_second =
        mesh::import_vtu_ascii_3d(
            *as_vtu.content);
    require(
        vtu_second.topology.entity_count(
            mesh::EntityKind::cell) == 6U &&
            vtu_second.cell_fields.size() == 4U,
        "GRDECL canonical VTU re-import");

    const auto as_gmsh =
        mesh::export_gmsh_4_1_ascii(
            canonical);
    require(
        as_gmsh.exported() &&
            as_gmsh.report.disposition() ==
                mesh::ConversionDisposition::lossy,
        "GRDECL canonical->Gmsh field/logical loss report");
    const auto gmsh_second =
        mesh::import_gmsh_4_1_ascii_3d(
            *as_gmsh.content,
            1.0);
    require(
        gmsh_second.topology.entity_count(
            mesh::EntityKind::cell) == 6U,
        "GRDECL canonical Gmsh re-import");

    std::cout
        << "[PASS] external.grdecl.opm_tube_native_roundtrip\n"
        << "[PASS] external.grdecl.opm_tube_to_vtu\n"
        << "[PASS] external.grdecl.opm_tube_to_gmsh\n"
        << "matrix.grdecl.gmsh=lossy\n"
        << "matrix.grdecl.vtu=lossy\n"
        << "matrix.grdecl.grdecl=lossless\n"
        << "grdecl_native_export=implemented_minimal_baseline\n"
        << "grdecl_test_scale_assumption=coordinate:1,permeability:1\n";
}

void verify_grdecl_wrapper_rejection(
    std::string_view content) {
    bool rejected = false;
    std::string message;
    try {
        (void)mesh::import_grdecl(
            content,
            mesh::GrdeclImportOptions{
                1.0,
                1.0});
    } catch (
        const std::invalid_argument& error) {
        rejected = true;
        message = error.what();
    }
    require(
        rejected &&
            message.find(
                "unsupported keyword 'GRID'") !=
                std::string::npos,
        "OPM 27cellsAniso.grdecl must expose current GRID-wrapper limitation");

    std::cout
        << "[PASS] external.grdecl.opm_27cells_expected_grid_wrapper_rejection"
        << '\n';
}

} // namespace

int main(
    int argc,
    char** argv) {
    if (argc != 6) {
        std::cerr
            << "[FAIL] external.mesh.compatibility: "
            << "usage: external_mesh_compatibility <gmsh2d> <gmsh3d> <vtu_bundle> <grdecl_tube> <grdecl_27cells>"
            << '\n';
        return 1;
    }

    int failures = 0;
    const auto run =
        [&](std::string_view name,
            auto&& function) {
            try {
                function();
            } catch (const std::exception& error) {
                ++failures;
                std::cerr
                    << "[FAIL] "
                    << name
                    << ": "
                    << error.what()
                    << '\n';
            }
        };

    run(
        "external.gmsh.2d.dealii",
        [&] {
            verify_gmsh_2d(
                read_file(argv[1]));
        });
    run(
        "external.gmsh.3d.dealii",
        [&] {
            verify_gmsh_3d(
                read_file(argv[2]));
        });
    run(
        "external.vtu.2d3d.dealii",
        [&] {
            verify_vtu_bundle(
                read_file(argv[3]));
        });
    run(
        "external.grdecl.opm_tube_to_vtu",
        [&] {
            verify_grdecl_tube(
                read_file(argv[4]));
        });
    run(
        "external.grdecl.opm_27cells_expected_grid_wrapper_rejection",
        [&] {
            verify_grdecl_wrapper_rejection(
                read_file(argv[5]));
        });

    if (failures != 0) {
        std::cerr
            << "[FAIL] external.mesh.compatibility failures="
            << failures
            << '\n';
        return 1;
    }

    std::cout
        << "[PASS] external.mesh.compatibility"
        << '\n';
    return 0;
}
