#include <mpmc/mesh/mesh_exchange_io.hpp>

#include <cstdint>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <string_view>
#include <vector>

bool mesh_exchange_group_header_self_contained();

namespace {
namespace mesh = mpmc::mesh;
void require(bool value, const char* message) {
    if (!value) throw std::runtime_error(message);
}
bool equal_groups(std::span<const mesh::MeshExchangeGroup> a,
                  std::span<const mesh::MeshExchangeGroup> b) {
    if (a.size() != b.size()) return false;
    for (std::size_t i=0; i<a.size(); ++i) {
        if (a[i].location != b[i].location || a[i].tag != b[i].tag ||
            a[i].name != b[i].name || a[i].members != b[i].members) return false;
    }
    return true;
}
std::string change_array(std::string xml, std::string_view name, std::string body,
                         std::string attributes = "") {
    const auto marker=xml.find("Name=\""+std::string{name}+"\"");
    const auto start=xml.rfind("<DataArray",marker);
    const auto end=xml.find("</DataArray>",marker)+12U;
    require(marker!=std::string::npos && start!=std::string::npos,"group array fixture missing");
    if (attributes.empty()) attributes="type=\"UInt64\" format=\"ascii\"";
    xml.replace(start,end-start,"<DataArray Name=\""+std::string{name}+"\" "+attributes+">"+body+"</DataArray>");
    return xml;
}
} // namespace

void verify_vtu_group_contract(const mpmc::mesh::MeshExchangeDocument& source) {
    require(mesh_exchange_group_header_self_contained(),"group public header probe");
    using K=mesh::EntityKind;
    const auto& topology=source.topology();
    const auto vertex=topology.global_ids(K::vertex).front();
    const auto face=topology.global_ids(K::face).front();
    const auto cell=topology.global_ids(K::cell).front();
    // Preserve old repeated-member semantics; overlaps and tag reuse across
    // entity kinds are intentional. Numeric face labels remain independent.
    const std::vector<mesh::MeshExchangeGroup> groups{
        {K::face,17U,"inlet & <wall> \"A\"",{face}},
        {K::face,18U,"overlap",{face}},
        {K::cell,17U,"rock",{cell,cell}},
        {K::vertex,std::numeric_limits<std::uint32_t>::max(),"\xe7\xbd\x91\xe6\xa0\xbc",{vertex}},
        {K::edge,33U,"empty edge set",{}},
        {K::cell,44U,"",{}}};
    const auto document=mesh::MeshExchangeDocument::create(source.source_format(),source.dimension(),
        topology,{source.vertex_coordinates_m().begin(),source.vertex_coordinates_m().end()},
        source.face_boundary(),{},groups,std::nullopt);
    const auto exported=mesh::export_vtu_ascii(document);
    require(exported.exported() && exported.report.lossless(),"VTU group-only export must be lossless");
    const auto gmsh_report=mesh::analyze_conversion(document,mesh::MeshExchangeFormat::gmsh_4_1_ascii);
    const auto gmsh_export=mesh::export_gmsh_4_1_ascii(document);
    require(gmsh_report.disposition()==mesh::ConversionDisposition::unsupported && !gmsh_export.exported(),
            "unsupported overlapping Gmsh face groups must agree between preflight and export");
    const auto parse=[&](std::string_view text) {
        return source.dimension()==2 ? mesh::make_mesh_exchange_document(mesh::import_vtu_ascii(text)) :
            mesh::make_mesh_exchange_document(mesh::import_vtu_ascii_3d(text));
    };
    const auto read=parse(*exported.content);
    require(equal_groups(read.groups(),groups),"VTU group name/kind/tag/membership binding mismatch");
    const auto native=source.dimension()==2 ? mesh::export_vtu_ascii(mesh::import_vtu_ascii(*exported.content)) :
        mesh::export_vtu_ascii_3d(mesh::import_vtu_ascii_3d(*exported.content));
    require(equal_groups(parse(native).groups(),groups),"VTU native group roundtrip mismatch");
    const auto second=mesh::export_vtu_ascii(read);
    require(second.report.lossless() && equal_groups(parse(*second.content).groups(),groups),
            "VTU canonical group roundtrip mismatch");
    const auto& xml=*exported.content;
    std::vector<std::string> invalid;
    for (const auto name:mesh::vtu_detail::group_array_names) {
        auto missing=xml;
        const auto marker=missing.find("Name=\""+std::string{name}+"\"");
        const auto start=missing.rfind("<DataArray",marker);
        const auto end=missing.find("</DataArray>",marker)+12U;
        const auto original=missing.substr(start,end-start);
        missing.erase(start,end-start);
        invalid.push_back(missing);
        auto duplicate=xml;
        duplicate.insert(start,original);
        invalid.push_back(duplicate);
        for (const auto* association:{"</PointData>","</CellData>"}) {
            auto misplaced=missing;
            misplaced.insert(misplaced.find(association),original);
            invalid.push_back(misplaced);
        }
    }
    invalid.push_back(change_array(xml,"mpmc_group_schema_version","2","type=\"UInt32\" format=\"ascii\""));
    invalid.push_back(change_array(xml,"mpmc_group_schema_version","1","type=\"UInt32\" format=\"binary\""));
    invalid.push_back(change_array(xml,"mpmc_group_schema_version","1","type=\"UInt32\" format=\"ascii\" NumberOfComponents=\"2\""));
    invalid.push_back(change_array(xml,"mpmc_group_schema_version","1","type=\"UInt32\" format=\"ascii\" NumberOfTuples=\"2\""));
    invalid.push_back(change_array(xml,"mpmc_group_location","7 2 3 0 1 3","type=\"UInt8\" format=\"ascii\""));
    invalid.push_back(change_array(xml,"mpmc_group_dimension","9 9 9 9 9 9","type=\"UInt8\" format=\"ascii\""));
    invalid.push_back(change_array(xml,"mpmc_group_tag","0 18 17 4294967295 33 44","type=\"UInt32\" format=\"ascii\""));
    invalid.push_back(change_array(xml,"mpmc_group_tag","17 17 17 4294967295 33 44","type=\"UInt32\" format=\"ascii\""));
    invalid.push_back(change_array(xml,"mpmc_group_tag","4294967296 18 17 4294967295 33 44","type=\"UInt32\" format=\"ascii\""));
    invalid.push_back(change_array(xml,"mpmc_group_name_offsets","99999 99999 99999 99999 99999 99999"));
    invalid.push_back(change_array(xml,"mpmc_group_member_offsets","1 0 4 5 5 5"));
    invalid.push_back(change_array(xml,"mpmc_group_member_ids","-1 1 1 1 1"));
    invalid.push_back(change_array(xml,"mpmc_group_member_ids","1.5 1 1 1 1"));
    invalid.push_back(change_array(xml,"mpmc_group_member_ids","18446744073709551616 1 1 1 1"));
    auto bad_utf8=xml;
    const auto utf8_marker=bad_utf8.find("Name=\"mpmc_group_name_utf8\"");
    const auto begin=bad_utf8.find('>',utf8_marker)+1U;
    const auto number=bad_utf8.find_first_not_of(" \n\r\t",begin);
    const auto after=bad_utf8.find_first_of(" \n\r\t<",number);
    bad_utf8.replace(number,after-number,"255");
    invalid.push_back(bad_utf8);
    auto unknown=xml;
    unknown.insert(unknown.find("</FieldData>"),"<DataArray Name=\"mpmc_group_unknown\" type=\"UInt8\" format=\"ascii\">1</DataArray>");
    invalid.push_back(unknown);
    for (std::size_t i=0; i<invalid.size(); ++i) {
        bool rejected=false;
        try { (void)parse(invalid[i]); } catch (const std::invalid_argument&) { rejected=true; }
        require(rejected,("invalid VTU group table accepted: "+std::to_string(i)).c_str());
    }
    std::cout<<"[PASS] mesh.core.vtu_groups dimension="<<source.dimension()<<" invalid_tables="<<invalid.size()<<'\n';
}
