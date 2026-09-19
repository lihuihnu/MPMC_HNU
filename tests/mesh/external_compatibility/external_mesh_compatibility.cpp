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

[[nodiscard]] bool same_ids(
    std::span<const mesh::GlobalEntityId> left,
    std::span<const mesh::GlobalEntityId> right) {
    return left.size() == right.size() &&
           std::equal(
               left.begin(),
               left.end(),
               right.begin());
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
    const auto first =
        mesh::import_gmsh_4_1_ascii(
            content,
            1.0);
    require(
        first.topology.entity_count(
            mesh::EntityKind::vertex) == 4U &&
            first.topology.entity_count(
                mesh::EntityKind::face) == 4U &&
            first.topology.entity_count(
                mesh::EntityKind::cell) == 1U,
        "deal.II Gmsh 2D entity counts");

    const auto canonical =
        mesh::make_mesh_exchange_document(
            first);
    require(
        canonical.source_format() ==
                mesh::MeshExchangeFormat::
                    gmsh_4_1_ascii &&
            canonical.dimension() == 2,
        "deal.II Gmsh 2D canonical adapter");
    require(
        mesh::analyze_conversion(
            canonical,
            mesh::MeshExchangeFormat::grdecl)
                .disposition() ==
            mesh::ConversionDisposition::
                unsupported,
        "generic Gmsh 2D must not silently become GRDECL");

    const auto exported =
        mesh::export_gmsh_4_1_ascii(
            first);
    const auto second =
        mesh::import_gmsh_4_1_ascii(
            exported,
            1.0);

    for (const auto kind :
         {mesh::EntityKind::vertex,
          mesh::EntityKind::face,
          mesh::EntityKind::cell}) {
        require(
            same_ids(
                first.topology.global_ids(kind),
                second.topology.global_ids(kind)),
            "deal.II Gmsh 2D stable IDs roundtrip");
    }
    require(
        second.face_boundary.face_count() ==
            first.face_boundary.face_count(),
        "deal.II Gmsh 2D boundary roundtrip");

    std::cout
        << "[PASS] external.gmsh.2d.dealii"
        << '\n';
}

void verify_gmsh_3d(
    std::string_view content) {
    const auto first =
        mesh::import_gmsh_4_1_ascii_3d(
            content,
            1.0);
    require(
        first.topology.entity_count(
            mesh::EntityKind::vertex) == 8U &&
            first.topology.entity_count(
                mesh::EntityKind::face) == 6U &&
            first.topology.entity_count(
                mesh::EntityKind::cell) == 1U,
        "deal.II Gmsh 3D entity counts");
    require_close(
        first.cell_volumes_m3[0],
        1.0,
        1.0e-14,
        "deal.II Gmsh 3D unit-cube volume");

    const auto canonical =
        mesh::make_mesh_exchange_document(
            first);
    require(
        canonical.source_format() ==
                mesh::MeshExchangeFormat::
                    gmsh_4_1_ascii &&
            canonical.dimension() == 3,
        "deal.II Gmsh 3D canonical adapter");
    require(
        mesh::analyze_conversion(
            canonical,
            mesh::MeshExchangeFormat::grdecl)
                .disposition() ==
            mesh::ConversionDisposition::
                unsupported,
        "generic Gmsh 3D must not silently become GRDECL");

    const auto exported =
        mesh::export_gmsh_4_1_ascii_3d(
            first);
    const auto second =
        mesh::import_gmsh_4_1_ascii_3d(
            exported,
            1.0);

    for (const auto kind :
         {mesh::EntityKind::vertex,
          mesh::EntityKind::face,
          mesh::EntityKind::cell}) {
        require(
            same_ids(
                first.topology.global_ids(kind),
                second.topology.global_ids(kind)),
            "deal.II Gmsh 3D stable IDs roundtrip");
    }
    require_close(
        second.cell_volumes_m3[0],
        first.cell_volumes_m3[0],
        1.0e-14,
        "deal.II Gmsh 3D volume roundtrip");

    std::cout
        << "[PASS] external.gmsh.3d.dealii"
        << '\n';
}

void verify_vtu_bundle(
    std::string_view content) {
    const auto documents =
        split_vtu_documents(content);
    require(
        documents.size() == 2U,
        "deal.II VTU external bundle must contain exactly two VTKFile documents");

    const auto first_2d =
        mesh::import_vtu_ascii(
            documents[0]);
    require(
        first_2d.topology.entity_count(
            mesh::EntityKind::vertex) == 4U &&
            first_2d.topology.entity_count(
                mesh::EntityKind::cell) == 1U &&
            first_2d.point_fields.size() == 2U,
        "deal.II VTU 2D quad import");
    const auto canonical_2d =
        mesh::make_mesh_exchange_document(
            first_2d);
    require(
        canonical_2d.source_format() ==
                mesh::MeshExchangeFormat::vtu_ascii &&
            canonical_2d.dimension() == 2 &&
            canonical_2d.fields().size() ==
                first_2d.point_fields.size() +
                    first_2d.cell_fields.size(),
        "deal.II VTU 2D canonical adapter");

    const auto exported_2d =
        mesh::export_vtu_ascii(
            first_2d);
    const auto second_2d =
        mesh::import_vtu_ascii(
            exported_2d);
    require(
        second_2d.topology.entity_count(
            mesh::EntityKind::cell) == 1U &&
            second_2d.point_fields.size() ==
                first_2d.point_fields.size(),
        "deal.II VTU 2D roundtrip");

    const auto first_3d =
        mesh::import_vtu_ascii_3d(
            documents[1]);
    require(
        first_3d.topology.entity_count(
            mesh::EntityKind::vertex) == 8U &&
            first_3d.topology.entity_count(
                mesh::EntityKind::face) == 6U &&
            first_3d.topology.entity_count(
                mesh::EntityKind::cell) == 1U &&
            first_3d.point_fields.size() == 2U,
        "deal.II VTU 3D hexa import");
    require_close(
        first_3d.cell_volumes_m3[0],
        1.0,
        1.0e-14,
        "deal.II VTU 3D unit-cube volume");

    const auto canonical_3d =
        mesh::make_mesh_exchange_document(
            first_3d);
    require(
        canonical_3d.source_format() ==
                mesh::MeshExchangeFormat::vtu_ascii &&
            canonical_3d.dimension() == 3 &&
            canonical_3d.fields().size() ==
                first_3d.point_fields.size() +
                    first_3d.cell_fields.size(),
        "deal.II VTU 3D canonical adapter");
    require(
        mesh::analyze_conversion(
            canonical_3d,
            mesh::MeshExchangeFormat::grdecl)
                .disposition() ==
            mesh::ConversionDisposition::
                unsupported,
        "generic VTU 3D must not silently become GRDECL");

    const auto exported_3d =
        mesh::export_vtu_ascii_3d(
            first_3d);
    const auto second_3d =
        mesh::import_vtu_ascii_3d(
            exported_3d);
    require(
        same_ids(
            first_3d.topology.global_ids(
                mesh::EntityKind::cell),
            second_3d.topology.global_ids(
                mesh::EntityKind::cell)) &&
            second_3d.point_fields.size() ==
                first_3d.point_fields.size(),
        "deal.II VTU 3D roundtrip");

    std::cout
        << "[PASS] external.vtu.2d3d.dealii"
        << '\n';
}

void verify_grdecl_tube(
    std::string_view content) {
    // This is a format-compatibility test, not a unit-system validation.
    // The standalone upstream GRDECL does not declare coordinate/permeability
    // units, so explicit no-conversion scales are used here and reported as
    // a test assumption rather than as physical interpretation.
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
    const auto native_report =
        mesh::analyze_conversion(
            canonical,
            mesh::MeshExchangeFormat::grdecl);
    require(
        native_report.lossless() &&
            native_report.issues().empty(),
        "OPM tube canonical GRDECL report");

    const auto native_export =
        mesh::export_grdecl_ascii(
            canonical);
    require(
        native_export.exported() &&
            native_export.report.lossless(),
        "OPM tube native GRDECL export");
    const auto native_second =
        mesh::import_grdecl(
            *native_export.content,
            raw.source_options);
    require(
        native_second.dimensions ==
                raw.dimensions &&
            native_second.active ==
                raw.active &&
            native_second.coord_m.size() ==
                raw.coord_m.size() &&
            native_second.zcorn_m.size() ==
                raw.zcorn_m.size(),
        "OPM native GRDECL structural roundtrip");
    for (std::size_t index = 0U;
         index < raw.coord_m.size();
         ++index) {
        require_close(
            native_second.coord_m[index],
            raw.coord_m[index],
            1.0e-12,
            "OPM native GRDECL COORD roundtrip");
    }
    for (std::size_t index = 0U;
         index < raw.zcorn_m.size();
         ++index) {
        require_close(
            native_second.zcorn_m[index],
            raw.zcorn_m[index],
            1.0e-12,
            "OPM native GRDECL ZCORN roundtrip");
    }
    for (const auto id :
         {"PORO", "PERMX", "PERMY", "PERMZ"}) {
        const auto find_field =
            [&](const mesh::GrdeclImportResult& value)
                -> const mesh::DenseFieldSnapshot& {
                const auto found =
                    std::find_if(
                        value.cell_fields.begin(),
                        value.cell_fields.end(),
                        [&](const auto& candidate) {
                            return candidate.metadata().id ==
                                id;
                        });
                if (found ==
                    value.cell_fields.end()) {
                    throw std::runtime_error(
                        "OPM native GRDECL field missing");
                }
                return *found;
            };
        const auto& expected =
            find_field(raw);
        const auto& actual =
            find_field(native_second);
        require(
            expected.values().size() ==
                actual.values().size(),
            "OPM native GRDECL field size");
        for (std::size_t index = 0U;
             index < expected.values().size();
             ++index) {
            require_close(
                actual.values()[index],
                expected.values()[index],
                1.0e-12,
                "OPM native GRDECL property roundtrip");
        }
    }

    const auto processed =
        mesh::process_active_corner_point_grid(
            raw);
    require(
        processed.cell_count() == 6U &&
            processed.cell_fields.size() == 4U,
        "OPM tube.grdecl active processing");

    const mesh::VtuImportResult3D as_vtu{
        processed.topology,
        processed.vertex_coordinates_m,
        processed.cell_volumes_m3,
        processed.face_geometry,
        mesh::make_face_boundary_snapshot(
            processed.topology),
        {},
        processed.cell_fields};

    const auto exported =
        mesh::export_vtu_ascii_3d(
            as_vtu);
    const auto second =
        mesh::import_vtu_ascii_3d(
            exported);

    require(
        second.topology.entity_count(
            mesh::EntityKind::cell) ==
            processed.cell_count(),
        "OPM GRDECL to VTU cell count");
    require(
        same_ids(
            processed.topology.global_ids(
                mesh::EntityKind::cell),
            second.topology.global_ids(
                mesh::EntityKind::cell)),
        "OPM GRDECL to VTU stable cell IDs");
    require(
        second.cell_fields.size() ==
            processed.cell_fields.size(),
        "OPM GRDECL to VTU property fields");

    for (std::size_t cell = 0U;
         cell < processed.cell_volumes_m3.size();
         ++cell) {
        require_close(
            second.cell_volumes_m3[cell],
            processed.cell_volumes_m3[cell],
            1.0e-10,
            "OPM GRDECL to VTU volume roundtrip");
    }

    std::cout
        << "[PASS] external.grdecl.opm_tube_native_roundtrip"
        << '\n';
    std::cout
        << "[PASS] external.grdecl.opm_tube_to_vtu"
        << '\n';
    std::cout
        << "grdecl_native_export=implemented_minimal_baseline"
        << '\n';
    std::cout
        << "grdecl_test_scale_assumption=coordinate:1,permeability:1"
        << '\n';
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
