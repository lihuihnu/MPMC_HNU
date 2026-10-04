#include <mpmc/mesh/face_mesh_geometry.hpp>
#include <mpmc/mesh/face_mesh_vtu.hpp>
#include <mpmc/mesh/face_mesh_exchange.hpp>
#include <iostream>
#include <sstream>
#include <source_location>

namespace {
using namespace mpmc::mesh;
void check(bool v,const std::source_location where=std::source_location::current()) {
    if(!v) throw std::runtime_error("face mesh test failed at line "+std::to_string(where.line()));
}
template<class F> void invalid(F f) { bool rejected=false; try {f();} catch(const std::exception&) {rejected=true;} check(rejected); }
FaceMesh cube() {
    FaceMesh m;
    m.points={{0,0,0},{1,0,0},{1,1,0},{0,1,0},{0,0,1},{1,0,1},{1,1,1},{0,1,1}};
    m.node_ids={0,1,2,3,4,5,6,18446744073709551615ULL};
    m.cell_ids={9007199254740993ULL}; m.face_ids={10,30,50,70,90,110};
    m.face_offsets={0,4,8,12,16,20,24};
    m.face_nodes={0,3,2,1,4,5,6,7,0,1,5,4,1,2,6,5,2,3,7,6,3,0,4,7};
    m.cell_offsets={0,6}; m.cell_faces={0,1,2,3,4,5}; m.signs={1,1,1,1,1,1};
    m.arrays["mrst_index_map"]={"cell","1",1,std::vector<MeshIndex>{20}};
    m.arrays["rock_perm"]={"cell","m2",3,std::vector<double>{1e-12,2e-12,3e-12}};
    return m;
}
}
int main() {
    try {
        auto m=cube();
        auto g=prepare_face_mesh_geometry(m);
        check(std::abs(g.cell_volumes[0]-1)<1e-14);
        for(auto x:g.cell_centroids[0]) check(std::abs(x-0.5)<1e-14);
        check(face_mesh_topology(m).entity_count(EntityKind::cell)==1);
        check(face_mesh_cell_graph(m).entry_count()==1);
        auto split=m;
        split.points.push_back({0.5,0,0}); split.node_ids.push_back(1000);
        split.face_nodes.insert(split.face_nodes.begin()+4,8);
        for(std::size_t i=1;i<split.face_offsets.size();++i) ++split.face_offsets[i];
        check(std::abs(prepare_face_mesh_geometry(split).cell_volumes[0]-1)<1e-14);
        split.points.back()[2]=0.02;
        (void)validate_face_mesh(split); // Lossless transport still meaningful.
        invalid([&]{(void)prepare_face_mesh_geometry(split);});
        // NNC is a graph edge between geometrically disconnected cells. It
        // does not invent a shared face or geometrical normal.
        auto nnc=m; nnc.arrays.clear();
        std::iota(nnc.node_ids.begin(),nnc.node_ids.end(),0ULL);
        for(auto point:m.points) { point[0]+=3; nnc.points.push_back(point);nnc.node_ids.push_back(nnc.node_ids.size()); }
        for(auto id:m.face_ids) nnc.face_ids.push_back(id+1000);
        for(auto node:m.face_nodes) nnc.face_nodes.push_back(node+8);
        for(std::size_t i=1;i<m.face_offsets.size();++i) nnc.face_offsets.push_back(m.face_offsets[i]+24);
        nnc.cell_ids={10,20}; nnc.cell_offsets.push_back(12);
        for(auto f:m.cell_faces) nnc.cell_faces.push_back(f+6);
        nnc.signs.insert(nnc.signs.end(),6,1);
        check(face_mesh_cell_graph(nnc).entry_count()==2);
        nnc.arrays["nnc_cells"]={"nnc","1",2,std::vector<MeshIndex>{0,1}};
        nnc.arrays["nnc_transmissibility"]={"nnc","m3/(Pa*s)",1,std::vector<double>{1e-10}};
        check(face_mesh_cell_graph(nnc).entry_count()==4 && nnc.face_ids.size()==12);
        check(prepare_face_mesh_geometry(nnc).cell_volumes.size()==2);
        auto linear_source=m;
        std::iota(linear_source.node_ids.begin(),linear_source.node_ids.end(),0ULL);
        linear_source.cell_ids={123};
        const auto legacy=export_face_mesh_legacy(linear_source,MeshExchangeFormat::gmsh_4_1_ascii);
        check(legacy.exported() && !legacy.report.lossless());
        const auto legacy_back=make_face_mesh(make_mesh_exchange_document(import_gmsh_4_1_ascii_3d(*legacy.content,1.0)));
        check(legacy_back.cell_ids.size()==1 && legacy_back.face_ids.size()==6);
        const auto old_doc=make_mesh_exchange_document(import_gmsh_4_1_ascii_3d(*legacy.content,1.0));
        const auto old_vtu=export_vtu_ascii(old_doc);
        check(old_vtu.exported());
        const auto old_vtu_back=import_face_mesh_vtu(*old_vtu.content,3);
        check(old_vtu_back.face_ids==legacy_back.face_ids && old_vtu_back.face_nodes==legacy_back.face_nodes);
        const auto deck=export_face_mesh_legacy(linear_source,MeshExchangeFormat::grdecl);
        check(deck.exported() && !deck.report.lossless());
        std::ostringstream xml; write_face_mesh_vtu(xml,m);
        auto r=read_face_mesh_vtu(xml.str(),3);
        check(r.node_ids==m.node_ids && r.face_ids==m.face_ids && r.cell_ids==m.cell_ids);
        check(r.face_nodes==m.face_nodes && r.cell_faces==m.cell_faces && r.signs==m.signs);
        check(r.arrays.at("rock_perm").values==m.arrays.at("rock_perm").values);
        auto compressed=xml.str();
        const auto vtk_tag=compressed.find("<VTKFile");
        compressed.insert(vtk_tag+8," compressor=\"unsupported\"");
        invalid([&]{(void)read_face_mesh_vtu(compressed,3);});
        auto conflict=m;
        conflict.arrays["mpmc_global_cell_id"]={"cell","1",1,std::vector<MeshIndex>{10}};
        invalid([&]{std::ostringstream out;write_face_mesh_vtu(out,conflict);});
        auto bad=m; bad.signs[0]=-1; invalid([&]{(void)validate_face_mesh(bad,true);});
        bad=m; bad.face_nodes[1]=bad.face_nodes[0]; invalid([&]{(void)validate_face_mesh(bad,true);});
        bad=m; bad.cell_offsets.back()=999; invalid([&]{(void)validate_face_mesh(bad,true);});
        bad=m; bad.node_ids[1]=bad.node_ids[0]; invalid([&]{(void)validate_face_mesh(bad,true);});
        bad=m; bad.arrays["rock_perm"].components=2; invalid([&]{(void)validate_face_mesh(bad,true);});
        bad=m; for(auto& p:bad.points) p[0]=-p[0]; invalid([&]{(void)prepare_face_mesh_geometry(bad);});
        bad=m; bad.arrays["reference_cell_volumes"]={"cell","m3",1,std::vector<double>{2}};
        invalid([&]{(void)prepare_face_mesh_geometry(bad);});
        // Translation stability; warped shared surface has scalar area distinct
        // from vector area. Divergence integration still gives positive volume.
        auto warped=m; warped.points[6][2]=1.2;
        const auto wg=prepare_face_mesh_geometry(warped);
        check(wg.cell_volumes[0]>1 && wg.cell_volumes[0]<1.2);
        check(wg.face_areas[1]>face_mesh_detail::norm(wg.area_vectors[1]));
        for(auto& p:warped.points) for(auto& x:p) x+=1e6;
        check(std::abs(prepare_face_mesh_geometry(warped).cell_volumes[0]-wg.cell_volumes[0])<1e-9);
        // Genuine pentagon, including collinear-independent area oracle.
        FaceMesh p; p.dimension=2; p.points={{0,0,0},{2,0,0},{2,1,0},{1,2,0},{0,1,0}};
        p.node_ids={4,8,12,16,20}; p.face_ids={2,3,4,5,6}; p.cell_ids={99};
        p.face_offsets={0,2,4,6,8,10}; p.face_nodes={0,1,1,2,2,3,3,4,4,0};
        p.cell_offsets={0,5};p.cell_faces={0,1,2,3,4};p.signs={1,1,1,1,1};
        check(std::abs(prepare_face_mesh_geometry(p).cell_volumes[0]-3)<1e-14);
        check(!export_face_mesh_legacy(p,MeshExchangeFormat::gmsh_4_1_ascii).exported());
        std::ostringstream px;write_face_mesh_vtu(px,p);check(read_face_mesh_vtu(px.str(),2).face_ids==p.face_ids);
        std::cout<<"[PASS] face mesh topology, geometry, typed IDs, VTU, failures\n";
    } catch(const std::exception& e) {std::cerr<<e.what()<<'\n';return 1;}
}
