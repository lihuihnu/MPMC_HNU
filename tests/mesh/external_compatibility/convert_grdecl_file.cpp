#include <mpmc/mesh/mesh_exchange_io.hpp>

#include <fstream>
#include <iomanip>
#include <iterator>
#include <limits>
#include <stdexcept>
#include <string>

namespace {
namespace mesh = mpmc::mesh;

template<class Values>
void write_values(std::ostream& output, const Values& values) {
    output << '[';
    bool first = true;
    for (const auto value : values) {
        if (!first) output << ',';
        first = false;
        output << +value;
    }
    output << ']';
}
} // namespace

// The audit is the imported state, not an expected-value fixture. The independent
// XTGeo oracle checks both this state and the exported GRDECL against its input.
void convert_grdecl_file(const char* input, const char* output_stem,
                         double coordinate_scale, double permeability_scale) {
    std::ifstream file(input, std::ios::binary);
    if (!file) throw std::runtime_error("cannot open GRDECL chain input");
    const std::string content{std::istreambuf_iterator<char>{file},
                              std::istreambuf_iterator<char>{}};
    const auto imported = mesh::import_grdecl(
        content, {coordinate_scale, permeability_scale});
    const auto processed = mesh::process_active_corner_point_grid(imported);
    const auto document = mesh::make_mesh_exchange_document(imported);
    const auto result = mesh::export_grdecl_ascii(document);
    if (!result.exported()) throw std::runtime_error("GRDECL chain export unsupported");
    const std::string stem{output_stem};
    std::ofstream output(stem + ".grdecl", std::ios::binary);
    output.exceptions(std::ios::badbit | std::ios::failbit);
    output << *result.content;
    std::ofstream report(stem + ".report", std::ios::binary);
    report.exceptions(std::ios::badbit | std::ios::failbit);
    report << (result.report.lossless() ? "lossless\n" : "lossy\n");
    for (const auto& issue : result.report.issues()) report << issue.code << '\n';

    std::ofstream audit(stem + ".import.json", std::ios::binary);
    audit.exceptions(std::ios::badbit | std::ios::failbit);
    audit << std::setprecision(std::numeric_limits<double>::max_digits10);
    audit << "{\"dimensions\":";
    write_values(audit, imported.dimensions);
    audit << ",\"coord_m\":";
    write_values(audit, imported.coord_m);
    audit << ",\"zcorn_m\":";
    write_values(audit, imported.zcorn_m);
    audit << ",\"actnum\":";
    write_values(audit, imported.active);
    audit << ",\"fields\":{";
    bool first = true;
    for (const auto& field : imported.cell_fields) {
        if (!first) audit << ',';
        first = false;
        audit << '"' << field.metadata().id << "\":";
        write_values(audit, field.values());
    }
    audit << "},\"active_logical_ids\":[";
    first = true;
    for (const auto id : processed.source_logical_cell_ids) {
        if (!first) audit << ',';
        first = false;
        audit << id.value();
    }
    audit << "],\"active_volumes_m3\":";
    write_values(audit, processed.cell_volumes_m3);
    audit << ",\"active_fields\":{";
    first = true;
    for (const auto& field : processed.cell_fields) {
        if (!first) audit << ',';
        first = false;
        audit << '"' << field.metadata().id << "\":";
        write_values(audit, field.values());
    }
    audit << "}}\n";
}
