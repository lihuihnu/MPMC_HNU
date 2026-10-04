#include <mpmc/mesh/face_mesh_hdf5.hpp>
#include <chrono>
#include <iostream>

int main(int argc,char** argv) {
    try {
        if(argc!=2) throw std::runtime_error("test output directory required");
        auto filename=std::filesystem::path(u8"网格-");
        filename+=std::to_string(std::chrono::steady_clock::now().time_since_epoch().count())+".h5";
        const auto path=std::filesystem::path(argv[1])/filename;
        mpmc::mesh::FaceMesh m; m.dimension=2;
        m.points={{0,0,0},{1,0,0},{0,1,0}}; m.node_ids={1,2,UINT64_MAX};
        m.face_ids={100,300,500}; m.cell_ids={9007199254740993ULL};
        m.face_offsets={0,2,4,6}; m.face_nodes={0,1,1,2,2,0};
        m.cell_offsets={0,3}; m.cell_faces={0,1,2}; m.signs={1,1,1};
        m.arrays["region"]={"cell","1",1,std::vector<std::int64_t>{-10}};
        m.arrays["nnc_cells"]={"nnc","1",2,std::vector<std::uint64_t>{}};
        mpmc::mesh::write_face_mesh_hdf5(path,m);
        const auto r=mpmc::mesh::read_face_mesh_hdf5(path);
        if(r.node_ids!=m.node_ids || r.face_nodes!=m.face_nodes || r.points!=m.points ||
           r.cell_ids!=m.cell_ids || r.signs!=m.signs || r.arrays.at("region").values!=m.arrays.at("region").values)
            throw std::runtime_error("HDF5 roundtrip mismatch");
        bool refused=false;
        try { (void)mpmc::mesh::read_face_mesh_hdf5(path,8); } catch(const std::length_error&) { refused=true; }
        if(!refused) throw std::runtime_error("memory budget not enforced");
        refused=false;
        try { mpmc::mesh::write_face_mesh_hdf5(path,m); } catch(const std::exception&) { refused=true; }
        if(!refused) throw std::runtime_error("existing file overwritten");
        std::filesystem::remove(path);
        std::cout<<"[PASS] HDF5 typed roundtrip, empty arrays, budget, no overwrite\n";
    } catch(const std::exception& e) { std::cerr<<e.what()<<'\n'; return 1; }
}
