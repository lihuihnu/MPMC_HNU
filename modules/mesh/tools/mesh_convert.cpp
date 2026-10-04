#include <mpmc/mesh/face_mesh_hdf5.hpp>
#include <mpmc/mesh/face_mesh_geometry.hpp>
#include <mpmc/mesh/face_mesh_vtu.hpp>
#include <mpmc/mesh/face_mesh_exchange.hpp>
#include <fstream>
#include <iostream>
#include <iterator>
#include <chrono>

int main(int argc,char** argv) {
    try {
        if(argc<3) throw std::invalid_argument("usage: mpmc_mesh_convert input output [--2d] [--geometry] [--graph] [--grdecl-si]; formats: h5, vtu, msh, grdecl");
        bool geometry=false, graph=false, grdecl_si=false; int dim=3;
        for(int i=3;i<argc;++i) {
            const std::string arg=argv[i];
            if(arg=="--2d") dim=2; else if(arg=="--geometry") geometry=true; else if(arg=="--graph") graph=true;
            else if(arg=="--grdecl-si") grdecl_si=true;
            else throw std::invalid_argument("unknown argument: "+arg);
        }
        const auto start=std::chrono::steady_clock::now();
        const std::filesystem::path input=argv[1],output=argv[2];
        auto m=[&] {
            if(input.extension()==".h5") return mpmc::mesh::read_face_mesh_hdf5(input);
            std::ifstream in(input); if(!in) throw std::runtime_error("cannot open input");
            const std::string xml{std::istreambuf_iterator<char>(in),{}};
            if(input.extension()==".vtu") return mpmc::mesh::import_face_mesh_vtu(xml,dim);
            if(input.extension()==".msh") {
                if(dim==2) return mpmc::mesh::make_face_mesh(mpmc::mesh::make_mesh_exchange_document(mpmc::mesh::import_gmsh_4_1_ascii(xml,1.0)));
                return mpmc::mesh::make_face_mesh(mpmc::mesh::make_mesh_exchange_document(mpmc::mesh::import_gmsh_4_1_ascii_3d(xml,1.0)));
            }
            if(input.extension()==".grdecl" || input.extension()==".GRDECL") {
                if(!grdecl_si) throw std::invalid_argument("CLI GRDECL import requires explicit --grdecl-si (metres, permeability m2); use library scale options for other units");
                return mpmc::mesh::make_face_mesh(mpmc::mesh::make_mesh_exchange_document(mpmc::mesh::import_grdecl(xml,{1.0,1.0})));
            }
            throw std::invalid_argument("unsupported input extension");
        }();
        if(std::filesystem::exists(output)) throw std::invalid_argument("refusing to overwrite output");
        if(output.extension()==".h5") mpmc::mesh::write_face_mesh_hdf5(output,m);
        else if(output.extension()==".vtu") { std::ofstream out(output); mpmc::mesh::write_face_mesh_vtu(out,m); }
        else if(output.extension()==".msh" || output.extension()==".grdecl" || output.extension()==".GRDECL") {
            const auto format=output.extension()==".msh"?mpmc::mesh::MeshExchangeFormat::gmsh_4_1_ascii:mpmc::mesh::MeshExchangeFormat::grdecl;
            const auto converted=mpmc::mesh::export_face_mesh_legacy(m,format);
            for(const auto& issue:converted.report.issues()) std::cout<<"conversion_issue="<<issue.code<<" "<<issue.message<<'\n';
            if(!converted.exported()) throw std::invalid_argument("target cannot represent this face mesh");
            std::ofstream out(output);out<<*converted.content;if(!out) throw std::runtime_error("legacy export failed");
        } else throw std::invalid_argument("unsupported output extension");
        std::cout<<"transport=pass nodes="<<m.points.size()<<" faces="<<m.face_ids.size()<<" cells="<<m.cell_ids.size()<<'\n';
        int result=0;
        if(geometry) {
            try {
                const auto g=mpmc::mesh::prepare_face_mesh_geometry(m);
                const auto total=std::accumulate(g.cell_volumes.begin(),g.cell_volumes.end(),0.0);
                std::cout<<std::setprecision(17)<<"geometry=pass volume="<<total<<" closure="<<g.maximum_relative_closure<<" reference_error="<<g.maximum_reference_relative_error<<" reference_values="<<g.compared_reference_values<<'\n';
            } catch(const std::exception& e) {
                std::cerr<<"geometry=fail "<<e.what()<<'\n'; result=2;
            }
        }
        if(graph) {
            const auto adj=mpmc::mesh::face_mesh_cell_graph(m);
            std::cout<<"graph=pass rows="<<adj.source_count()<<" entries="<<adj.entry_count()<<'\n';
        }
        std::cout<<"seconds="<<std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count()<<'\n';
        return result;
    } catch(const std::exception& e) { std::cerr<<e.what()<<'\n'; return 1; }
}
