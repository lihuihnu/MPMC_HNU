#include <mpmc/mesh/mesh_exchange_io.hpp>
#include <fstream>
#include <iomanip>
#include <iterator>
#include <stdexcept>
#include <string>

void convert_vtu_groups(const char* dimension, const char* input, const char* stem_arg) {
    std::ifstream file(input, std::ios::binary);
    if (!file) throw std::runtime_error("cannot open group-chain input");
    const std::string text{std::istreambuf_iterator<char>{file},std::istreambuf_iterator<char>{}};
    namespace mesh = mpmc::mesh;
    if (std::string{dimension} != "2" && std::string{dimension} != "3")
        throw std::invalid_argument("group chain dimension must be 2 or 3");
    const auto document = std::string{dimension} == "2" ?
        mesh::make_mesh_exchange_document(mesh::import_vtu_ascii(text)) :
        mesh::make_mesh_exchange_document(mesh::import_vtu_ascii_3d(text));
    const std::string stem{stem_arg};
    std::ofstream audit(stem+".groups.json");
    audit.exceptions(std::ios::badbit | std::ios::failbit);
    audit << '[';
    bool first=true;
    for (const auto& group:document.groups()) {
        if (!first) audit << ',';
        first=false;
        audit << "{\"location\":" << static_cast<unsigned>(group.location)
              << ",\"dimension\":" << mesh::vtu_detail::group_dimension(group.location,document.dimension())
              << ",\"tag\":" << group.tag << ",\"name\":" << std::quoted(group.name) << ",\"members\":[";
        bool first_member=true;
        for (const auto id:group.members) {
            if (!first_member) audit << ',';
            first_member=false;
            audit << id.value();
        }
        audit << "]}";
    }
    audit << "]\n";
    const auto result=mesh::export_vtu_ascii(document);
    if (!result.exported() || !result.report.lossless()) throw std::runtime_error("group chain must be lossless");
    std::ofstream output(stem+".vtu",std::ios::binary);
    output.exceptions(std::ios::badbit | std::ios::failbit);
    output << *result.content;
    std::ofstream report(stem+".report");
    report.exceptions(std::ios::badbit | std::ios::failbit);
    report << "lossless\n";
}
