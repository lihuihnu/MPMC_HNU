#include <mpmc/mesh/active_corner_point.hpp>
#include <mpmc/mesh/gmsh_4_1.hpp>
#include <mpmc/mesh/gmsh_4_1_3d.hpp>
#include <mpmc/mesh/grdecl.hpp>
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
        << "[PASS] external.grdecl.opm_tube_to_vtu"
        << '\n';
    std::cout
        << "grdecl_native_export=not_implemented"
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
    try {
        require(
            argc == 6,
            "usage: external_mesh_compatibility <gmsh2d> <gmsh3d> <vtu_bundle> <grdecl_tube> <grdecl_27cells>");

        verify_gmsh_2d(
            read_file(argv[1]));
        verify_gmsh_3d(
            read_file(argv[2]));
        verify_vtu_bundle(
            read_file(argv[3]));
        verify_grdecl_tube(
            read_file(argv[4]));
        verify_grdecl_wrapper_rejection(
            read_file(argv[5]));

        std::cout
            << "[PASS] external.mesh.compatibility"
            << '\n';
        return 0;
    } catch (const std::exception& error) {
        std::cerr
            << "[FAIL] external.mesh.compatibility: "
            << error.what()
            << '\n';
        return 1;
    }
}
