#include <mpmc/mesh/gmsh_4_1.hpp>
#include <mpmc/mesh/linear_cell_mesh_2d.hpp>
#include <mpmc/mesh/mesh_exchange_io.hpp>
#include <mpmc/mesh/vtu.hpp>

#include <array>
#include <filesystem>
#include <fstream>
#include <stdexcept>
#include <string>
#include <string_view>
#include <vector>

namespace {
namespace mesh = mpmc::mesh;
using mesh::GlobalEntityId;
using mesh::LocalIndex;
using mesh::LinearCellType2D;

void write(const std::filesystem::path& path, const std::string& text) {
    std::ofstream stream(path, std::ios::binary);
    stream.exceptions(std::ios::badbit | std::ios::failbit);
    stream << text;
}

void write_conversion(const std::filesystem::path& path, const mesh::MeshTextExportResult& result) {
    if (!result.exported()) throw std::runtime_error("2D test export unexpectedly unsupported");
    write(path, *result.content);
    std::string report = result.report.lossless() ? "lossless\n" : "lossy\n";
    for (const auto& issue : result.report.issues()) report += issue.code + "\n";
    write(path.string() + ".report", report);
}

mesh::DenseFieldMetadata metadata(const char* id, const char* unit) {
    return {id, unit, {mesh::FieldSourceKind::synthetic_test,
        "analytic://triangle-trapezoid", "v1", "emit_2d_exports.cpp"}};
}

void emit_case(const std::filesystem::path& directory, std::string_view name) {
    // Explicit synthetic inputs. The Python oracle has independent analytic
    // expectations; no expected topology/geometry is serialized by this driver.
    std::vector<GlobalEntityId> ids{GlobalEntityId{50}, GlobalEntityId{7}, GlobalEntityId{90}};
    std::vector<mesh::Coordinate2D> coordinates{{0,0}, {2,0}, {0,1}};
    std::vector<mesh::LinearCell2D> cells{{GlobalEntityId{901}, LinearCellType2D::triangle,
        {LocalIndex{0}, LocalIndex{1}, LocalIndex{2}}}};
    std::array<mesh::LinearFaceAnnotation2D, 2> annotations{{
        {{LocalIndex{0}, LocalIndex{1}}, GlobalEntityId{701}, mesh::PhysicalTag{11}},
        {{LocalIndex{1}, LocalIndex{2}}, GlobalEntityId{702}, mesh::PhysicalTag{12}}}};
    if (name != "triangle") {
        ids.emplace_back(12);
        coordinates = {{0,0}, {3,0}, {2,2}, {0,2}};
        cells = {{GlobalEntityId{901}, LinearCellType2D::quadrilateral,
            {LocalIndex{3}, LocalIndex{2}, LocalIndex{1}, LocalIndex{0}}}};
        if (name == "mixed") {
            ids.emplace_back(110);
            coordinates.push_back({4,0});
            cells.insert(cells.begin(), {GlobalEntityId{31}, LinearCellType2D::triangle,
                {LocalIndex{1}, LocalIndex{4}, LocalIndex{2}}});
            annotations[1].vertices = {LocalIndex{1}, LocalIndex{4}};
        }
    }
    auto grid = mesh::make_linear_mesh_2d(ids, coordinates, cells, annotations, 10000);
    std::vector<mesh::GmshCellPhysicalGroups> groups;
    for (std::size_t i = 0; i < cells.size(); ++i) {
        groups.push_back({cells[i].global_id, i == 0U ? std::vector<std::uint32_t>{21,22} :
                                                     std::vector<std::uint32_t>{21}});
    }
    mesh::Gmsh41ImportResult gmsh{grid.topology, grid.geometry, grid.face_boundary,
        {{1,11,"inlet & lower"}, {1,12,"outlet upper"},
         {2,21,"rock & sand"}, {2,22,"selected region"}}, std::move(groups)};
    std::vector<double> point_scalars, point_vectors, cell_scalars, cell_vectors;
    for (std::size_t i = 0; i < coordinates.size(); ++i) {
        point_scalars.push_back(300.0 + static_cast<double>(i));
        point_vectors.insert(point_vectors.end(),
            {coordinates[i].x_m, coordinates[i].y_m, coordinates[i].x_m - 2.0 * coordinates[i].y_m});
    }
    for (const auto& cell : cells) {
        const double marker = cell.global_id.value() == 31U ? 3.1 : 90.1;
        cell_scalars.push_back(marker);
        cell_vectors.insert(cell_vectors.end(), {marker, -marker});
    }
    mesh::VtuImportResult vtu{grid.topology, grid.geometry,
        {mesh::DenseFieldSnapshot::create(grid.topology, mesh::EntityKind::vertex, 1U,
             std::move(point_scalars), metadata("temperature", "K")),
         mesh::DenseFieldSnapshot::create(grid.topology, mesh::EntityKind::vertex, 3U,
             std::move(point_vectors), metadata("point_vector", "1"))},
        {mesh::DenseFieldSnapshot::create(grid.topology, mesh::EntityKind::cell, 1U,
             std::move(cell_scalars), metadata("marker", "1")),
         mesh::DenseFieldSnapshot::create(grid.topology, mesh::EntityKind::cell, 2U,
             std::move(cell_vectors), metadata("cell_pair", "1"))}};
    const std::string stem{name};
    write(directory / (stem + ".msh"), mesh::export_gmsh_4_1_ascii(gmsh));
    write(directory / (stem + ".vtu"), mesh::export_vtu_ascii(vtu));
    write_conversion(directory / (stem + "_from_gmsh.vtu"),
        mesh::export_vtu_ascii(mesh::make_mesh_exchange_document(gmsh)));
    write_conversion(directory / (stem + "_from_vtu.msh"),
        mesh::export_gmsh_4_1_ascii(mesh::make_mesh_exchange_document(vtu)));
}
} // namespace

void emit_2d_exports(const char* output_directory) {
    const std::filesystem::path directory{output_directory};
    std::filesystem::create_directories(directory);
    for (const auto name : {"triangle", "quad", "mixed"}) emit_case(directory, name);
}
