#include <mpmc/mesh/mesh_exchange_io.hpp>

#include <array>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <initializer_list>
#include <stdexcept>
#include <string>
#include <string_view>
#include <vector>

namespace {
namespace mesh = mpmc::mesh;

void write(const std::filesystem::path& path, const std::string& text) {
    std::ofstream stream(path, std::ios::binary);
    stream.exceptions(std::ios::badbit | std::ios::failbit);
    stream << text;
}

void write_conversion(const std::filesystem::path& path, const mesh::MeshTextExportResult& result) {
    if (!result.exported()) throw std::runtime_error("3D test export unexpectedly unsupported");
    write(path, *result.content);
    std::string report = result.report.lossless() ? "lossless\n" : "lossy\n";
    for (const auto& issue : result.report.issues()) report += issue.code + "\n";
    write(path.string() + ".report", report);
}

mesh::DenseFieldMetadata metadata(const char* id, const char* unit) {
    return {id, unit, {mesh::FieldSourceKind::synthetic_test,
        "analytic://linear-solids", "v1", "emit_3d_exports.cpp"}};
}

void write_grid(const std::filesystem::path& directory, std::string_view name,
                const mesh::LinearMesh3D& grid) {
    using mesh::GlobalEntityId;
    const auto& coordinates = grid.vertex_coordinates_m;
    mesh::Gmsh41ImportResult3D gmsh{grid.topology, grid.vertex_coordinates_m, grid.cell_volumes_m3,
        grid.face_geometry, grid.face_boundary,
        {{2,11,"base & inlet"}, {2,12,"translated base"}, {3,21,"rock & sand"}, {3,22,"selected region"}},
        {{GlobalEntityId{31}, {21,22}}, {GlobalEntityId{901}, {21}}}};
    std::vector<double> temperature, position;
    for (std::size_t i = 0; i < coordinates.size(); ++i) {
        temperature.push_back(300.0 + static_cast<double>(i));
        const auto p = coordinates[i];
        position.insert(position.end(), {p.x_m, p.y_m, p.z_m});
    }
    mesh::VtuImportResult3D vtu{grid.topology, grid.vertex_coordinates_m, grid.cell_volumes_m3,
        grid.face_geometry, mesh::make_face_boundary_snapshot(grid.topology),
        {mesh::DenseFieldSnapshot::create(grid.topology, mesh::EntityKind::vertex, 1U,
             std::move(temperature), metadata("temperature", "K")),
         mesh::DenseFieldSnapshot::create(grid.topology, mesh::EntityKind::vertex, 3U,
             std::move(position), metadata("position", "m"))},
        {mesh::DenseFieldSnapshot::create(grid.topology, mesh::EntityKind::cell, 1U,
             {3.1,90.1}, metadata("marker", "1")),
         mesh::DenseFieldSnapshot::create(grid.topology, mesh::EntityKind::cell, 2U,
             {3.1,-3.1,90.1,-90.1}, metadata("cell_pair", "1"))}};
    const std::string stem{name};
    write(directory / (stem + ".msh"), mesh::export_gmsh_4_1_ascii_3d(gmsh));
    write(directory / (stem + ".vtu"), mesh::export_vtu_ascii_3d(vtu));
    write_conversion(directory / (stem + "_from_gmsh.vtu"),
        mesh::export_vtu_ascii(mesh::make_mesh_exchange_document(gmsh)));
    write_conversion(directory / (stem + "_from_vtu.msh"),
        mesh::export_gmsh_4_1_ascii(mesh::make_mesh_exchange_document(vtu)));
}

void emit_case(const std::filesystem::path& directory, std::string_view name) {
    using mesh::GlobalEntityId;
    using mesh::LocalIndex;
    using mesh::LinearCellType3D;
    // Two disconnected solids per file: the second is scaled by 2 and
    // translated. Cell order deliberately differs from point block order.
    std::vector<mesh::Coordinate3D> base;
    LinearCellType3D type;
    if (name == "tetrahedron") {
        type = LinearCellType3D::tetrahedron;
        base = {{0,0,0}, {2,0,0}, {0,3,0}, {0,0,4}};
    } else if (name == "hexahedron") {
        type = LinearCellType3D::hexahedron;
        base = {{0,0,0}, {2,0,0}, {2,3,0}, {0,3,0},
                {0,0,4}, {2,0,4}, {2,3,4}, {0,3,4}};
    } else if (name == "wedge") {
        type = LinearCellType3D::wedge;
        base = {{0,0,0}, {2,0,0}, {0,3,0}, {0,0,4}, {2,0,4}, {0,3,4}};
    } else if (name == "pyramid") {
        type = LinearCellType3D::pyramid;
        base = {{0,0,0}, {2,0,0}, {2,3,0}, {0,3,0}, {1,1.5,4}};
    } else {
        throw std::invalid_argument("unknown 3D export fixture");
    }
    const std::array<std::uint64_t, 8> sparse{50,7,90,12,110,19,44,3};
    std::vector<GlobalEntityId> ids;
    std::vector<mesh::Coordinate3D> coordinates;
    std::array<std::vector<LocalIndex>, 2> rows;
    for (std::size_t block = 0; block < 2U; ++block) {
        for (std::size_t i = 0; i < base.size(); ++i) {
            ids.emplace_back(sparse[i] + 1000U * block);
            const auto p = base[i];
            coordinates.push_back(block == 0U ? p : mesh::Coordinate3D{10+2*p.x_m, -5+2*p.y_m, 2+2*p.z_m});
            rows[block].emplace_back(static_cast<LocalIndex::value_type>(coordinates.size()-1U));
        }
    }
    const std::vector<mesh::LinearCell3D> cells{
        {GlobalEntityId{31}, type, rows[1]}, {GlobalEntityId{901}, type, rows[0]}};
    const std::size_t base_width = (name == "tetrahedron" || name == "wedge") ? 3U : 4U;
    std::vector<mesh::LinearFaceAnnotation3D> annotations;
    for (std::size_t block = 0; block < 2U; ++block) {
        std::vector<LocalIndex> face;
        for (std::size_t i = 0; i < base_width; ++i) face.push_back(rows[block][i]);
        annotations.push_back({std::move(face), GlobalEntityId{701U+block},
            mesh::PhysicalTag{static_cast<std::uint32_t>(11U+block)}});
    }
    auto grid = mesh::make_linear_mesh_3d(ids, coordinates, cells, annotations, 10000);
    write_grid(directory, name, grid);
}

void emit_mixed(const std::filesystem::path& directory, bool triangle) {
    using mesh::GlobalEntityId;
    using mesh::LocalIndex;
    using mesh::LinearCellType3D;
    std::vector<mesh::Coordinate3D> coordinates = triangle
        ? std::vector<mesh::Coordinate3D>{{0,0,0}, {2,0,0}, {0,3,0},
              {0,0,4}, {2,0,4}, {0,3,4}, {0,0,8}}
        : std::vector<mesh::Coordinate3D>{{0,0,0}, {2,0,0}, {2,3,0}, {0,3,0},
              {0,0,4}, {2,0,4}, {2,3,4}, {0,3,4}, {1,1.5,8}};
    const std::array<std::uint64_t, 9> sparse{50,7,90,12,110,19,44,3,250};
    std::vector<GlobalEntityId> ids;
    for (std::size_t i = 0; i < coordinates.size(); ++i) ids.emplace_back(sparse[i]);
    const auto row = [](std::initializer_list<LocalIndex::value_type> indices) {
        std::vector<LocalIndex> result;
        for (const auto i : indices) result.emplace_back(i);
        return result;
    };
    const std::vector<mesh::LinearCell3D> cells = triangle
        ? std::vector<mesh::LinearCell3D>{
            {GlobalEntityId{31}, LinearCellType3D::tetrahedron, row({3,4,5,6})},
            {GlobalEntityId{901}, LinearCellType3D::wedge, row({0,1,2,3,4,5})}}
        : std::vector<mesh::LinearCell3D>{
            {GlobalEntityId{31}, LinearCellType3D::pyramid, row({4,5,6,7,8})},
            {GlobalEntityId{901}, LinearCellType3D::hexahedron, row({0,1,2,3,4,5,6,7})}};
    const std::vector<mesh::LinearFaceAnnotation3D> annotations{
        {triangle ? row({0,1,2}) : row({0,1,2,3}), GlobalEntityId{701}, mesh::PhysicalTag{11}},
        {triangle ? row({3,4,6}) : row({4,5,8}), GlobalEntityId{702}, mesh::PhysicalTag{12}},
        {triangle ? row({3,4,5}) : row({4,5,6,7}), GlobalEntityId{703}, mesh::PhysicalTag{0}}};
    const auto grid = mesh::make_linear_mesh_3d(ids, coordinates, cells, annotations, 10000);
    write_grid(directory, triangle ? "tetra_wedge" : "hexa_pyramid", grid);
}

} // namespace

void emit_3d_exports(const char* output_directory) {
    const std::filesystem::path directory{output_directory};
    std::filesystem::create_directories(directory);
    for (const auto name : {"tetrahedron", "hexahedron", "wedge", "pyramid"}) emit_case(directory, name);
    emit_mixed(directory, true);
    emit_mixed(directory, false);
}
