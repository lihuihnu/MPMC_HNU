#include <mpmc/mesh/mesh_exchange_io.hpp>

#include <fstream>
#include <iomanip>
#include <iterator>
#include <limits>
#include <stdexcept>
#include <string>
#include <string_view>

namespace {
namespace mesh = mpmc::mesh;

void write_export(const std::string& path, const mesh::MeshTextExportResult& result) {
    if (!result.exported()) throw std::runtime_error("2D chain export unsupported");
    std::ofstream file(path, std::ios::binary);
    file.exceptions(std::ios::badbit | std::ios::failbit);
    file << *result.content;
    std::ofstream report(path + ".report", std::ios::binary);
    report.exceptions(std::ios::badbit | std::ios::failbit);
    report << (result.report.lossless() ? "lossless\n" : "lossy\n");
    for (const auto& issue : result.report.issues()) report << issue.code << '\n';
}

// Raw imported state, before the canonical writer can reconstruct any geometry.
// This is evidence consumed by the Python oracle, never an expected-value file.
template<class Imported>
void convert(const Imported& imported, const std::string& stem) {
    using mesh::EntityKind;
    const auto& topology = imported.topology;
    std::ofstream audit(stem + ".import", std::ios::binary);
    audit.exceptions(std::ios::badbit | std::ios::failbit);
    audit << std::setprecision(std::numeric_limits<double>::max_digits10);
    for (const auto kind : {EntityKind::vertex, EntityKind::cell, EntityKind::face}) {
        const auto ids = topology.global_ids(kind);
        for (std::size_t i = 0; i < ids.size(); ++i) {
            const mesh::LocalIndex local{static_cast<mesh::LocalIndex::value_type>(i)};
            audit << static_cast<unsigned>(kind) << ' ' << ids[i].value();
            if (kind == EntityKind::vertex) {
                const auto p = imported.geometry.vertex_coordinates_m()[i];
                audit << ' ' << p.x_m << ' ' << p.y_m;
            } else {
                const auto vertices = topology.relation(kind, EntityKind::vertex).adjacent(local);
                audit << ' ' << vertices.size();
                for (const auto v : vertices)
                    audit << ' ' << topology.global_ids(EntityKind::vertex)[v.value()].value();
                if (kind == EntityKind::cell) {
                    const auto center = imported.geometry.cell_centroid_m(local);
                    audit << ' ' << imported.geometry.cell_area_m2(local)
                          << ' ' << center.x_m << ' ' << center.y_m;
                    const auto faces = topology.relation(kind, EntityKind::face).adjacent(local);
                    audit << ' ' << faces.size();
                    for (const auto face : faces)
                        audit << ' ' << topology.global_ids(EntityKind::face)[face.value()].value();
                } else {
                    const auto cells = topology.relation(kind, EntityKind::cell).adjacent(local);
                    audit << ' ' << cells.size();
                    for (const auto cell : cells)
                        audit << ' ' << topology.global_ids(EntityKind::cell)[cell.value()].value();
                    const auto owner = imported.geometry.face_owner(local);
                    const auto normal = imported.geometry.face_owner_unit_normal(local);
                    audit << ' ' << imported.face_boundary.is_boundary(local)
                          << ' ' << imported.face_boundary.physical_tag(local).value();
                    audit << ' ' << imported.geometry.face_length_m(local)
                          << ' ' << topology.global_ids(EntityKind::cell)[owner.value()].value()
                          << ' ' << normal.x << ' ' << normal.y;
                    const auto center = imported.geometry.face_centroid_m(local);
                    audit << ' ' << center.x_m << ' ' << center.y_m;
                }
            }
            audit << '\n';
        }
    }
    const auto document = mesh::make_mesh_exchange_document(imported);
    write_export(stem + ".msh", mesh::export_gmsh_4_1_ascii(document));
    write_export(stem + ".vtu", mesh::export_vtu_ascii(document));
}
} // namespace

void convert_2d_file(const char* format, const char* input, const char* output_stem) {
    std::ifstream file(input, std::ios::binary);
    if (!file) throw std::runtime_error("cannot open 2D chain input");
    const std::string content{std::istreambuf_iterator<char>{file}, std::istreambuf_iterator<char>{}};
    if (std::string_view{format} == "gmsh") {
        convert(mesh::import_gmsh_4_1_ascii(content, 1.0), output_stem);
    } else if (std::string_view{format} == "vtu") {
        convert(mesh::import_vtu_ascii(content), output_stem);
    } else {
        throw std::invalid_argument("2D chain input format must be gmsh or vtu");
    }
}
