#include <mpmc/mesh/face_mesh_hdf5.hpp>
#include <hdf5.h>
#include <algorithm>
#include <cstring>
#include <type_traits>

namespace mpmc::mesh {
namespace {
struct Handle {
    hid_t id; herr_t (*close)(hid_t);
    Handle(hid_t value,herr_t (*closer)(hid_t)):id(value),close(closer) {
        if(id<0) throw std::runtime_error("HDF5 operation failed");
    }
    ~Handle() { (void)close(id); }
    Handle(const Handle&)=delete;
    Handle& operator=(const Handle&)=delete;
    operator hid_t() const { return id; }
};
void check(herr_t value) { if(value<0) throw std::runtime_error("HDF5 operation failed"); }
std::string filename_utf8(const std::filesystem::path& file) {
    // HDF5 1.14.6 uses Wopen_utf8 on Windows, independently of the active ANSI code page.
    const auto bytes=file.u8string();
    return {reinterpret_cast<const char*>(bytes.data()),bytes.size()};
}
template<class T> hid_t native() {
    if constexpr(std::is_same_v<T,double>) return H5T_NATIVE_DOUBLE;
    else if constexpr(std::is_same_v<T,std::uint64_t>) return H5T_NATIVE_UINT64;
    else if constexpr(std::is_same_v<T,std::int64_t>) return H5T_NATIVE_INT64;
    else return H5T_NATIVE_UCHAR;
}
template<class T> hid_t disk() {
    if constexpr(std::is_same_v<T,double>) return H5T_IEEE_F64LE;
    else if constexpr(std::is_same_v<T,std::uint64_t>) return H5T_STD_U64LE;
    else if constexpr(std::is_same_v<T,std::int64_t>) return H5T_STD_I64LE;
    else return H5T_STD_U8LE;
}
void hard_link(hid_t f,const std::string& path) {
    H5L_info_t info{};
    check(H5Lget_info(f,path.c_str(),&info,H5P_DEFAULT));
    if(info.type!=H5L_TYPE_HARD) throw std::invalid_argument("HDF5 external/soft links are not mesh data");
}
template<class T> void put(hid_t f,const std::string& path,const T* data,std::size_t size) {
    const hsize_t n=size;
    Handle space(H5Screate_simple(1,&n,nullptr),H5Sclose);
    Handle link(H5Pcreate(H5P_LINK_CREATE),H5Pclose);
    check(H5Pset_create_intermediate_group(link,1));
    Handle dataset(H5Dcreate2(f,path.c_str(),disk<T>(),space,link,H5P_DEFAULT,H5P_DEFAULT),H5Dclose);
    if(size) check(H5Dwrite(dataset,native<T>(),H5S_ALL,H5S_ALL,H5P_DEFAULT,data));
}
template<class T> void put(hid_t f,const std::string& p,const std::vector<T>& v) { put(f,p,v.data(),v.size()); }
void text(hid_t f,const std::string& p,const std::string& v) {
    put(f,p,reinterpret_cast<const unsigned char*>(v.data()),v.size());
}
void put_scalar(hid_t f,const std::string& p,std::uint64_t v) { put(f,p,&v,1); }
template<class T> std::vector<T> get(hid_t f,const std::string& path,std::uint64_t& budget) {
    hard_link(f,path);
    Handle d(H5Dopen2(f,path.c_str(),H5P_DEFAULT),H5Dclose);
    Handle type(H5Dget_type(d),H5Tclose);
    const bool floating=std::is_same_v<T,double>;
    if(H5Tget_size(type)!=sizeof(T) || H5Tget_class(type)!=(floating?H5T_FLOAT:H5T_INTEGER) ||
        (!floating && H5Tget_sign(type)!=(std::is_signed_v<T>?H5T_SGN_2:H5T_SGN_NONE)))
        throw std::invalid_argument("HDF5 mesh datatype mismatch: "+path);
    Handle s(H5Dget_space(d),H5Sclose);
    if(H5Sget_simple_extent_ndims(s)!=1) throw std::invalid_argument("HDF5 mesh arrays must be flat rank 1: "+path);
    hsize_t count{}; check(H5Sget_simple_extent_dims(s,&count,nullptr));
    if(count>budget/sizeof(T) || count>std::numeric_limits<std::size_t>::max()/sizeof(T))
        throw std::length_error("HDF5 mesh read exceeds memory budget");
    budget-=count*sizeof(T);
    std::vector<T> v(static_cast<std::size_t>(count));
    // Bounded I/O blocks, directly into final storage; no dataset-sized staging.
    for(hsize_t begin=0;begin<count;) {
        const hsize_t n=std::min<hsize_t>(count-begin,1024*1024);
        check(H5Sselect_hyperslab(s,H5S_SELECT_SET,&begin,nullptr,&n,nullptr));
        Handle memory(H5Screate_simple(1,&n,nullptr),H5Sclose);
        check(H5Dread(d,native<T>(),memory,s,H5P_DEFAULT,v.data()+begin)); begin+=n;
    }
    return v;
}
std::string text(hid_t f,const std::string& p,std::uint64_t& budget) {
    const auto v=get<unsigned char>(f,p,budget);
    return {v.begin(),v.end()};
}
std::uint64_t get_scalar(hid_t f,const std::string& p,std::uint64_t& budget) {
    const auto v=get<std::uint64_t>(f,p,budget);
    if(v.size()!=1) throw std::invalid_argument("HDF5 expected scalar: "+p);
    return v[0];
}
std::vector<std::string> children(hid_t f,const std::string& p) {
    Handle g(H5Gopen2(f,p.c_str(),H5P_DEFAULT),H5Gclose);
    H5G_info_t info{}; check(H5Gget_info(g,&info));
    std::vector<std::string> result;
    for(hsize_t i=0;i<info.nlinks;++i) {
        auto n=H5Lget_name_by_idx(g,".",H5_INDEX_NAME,H5_ITER_INC,i,nullptr,0,H5P_DEFAULT);
        if(n<0 || n>4096) throw std::invalid_argument("HDF5 invalid link name");
        std::string name(static_cast<std::size_t>(n)+1,'\0');
        if(H5Lget_name_by_idx(g,".",H5_INDEX_NAME,H5_ITER_INC,i,name.data(),name.size(),H5P_DEFAULT)<0)
            throw std::runtime_error("HDF5 link read failed");
        name.resize(static_cast<std::size_t>(n)); hard_link(g,name); result.push_back(std::move(name));
    }
    return result;
}
}
void write_face_mesh_hdf5(const std::filesystem::path& file,const FaceMesh& m) {
    (void)validate_face_mesh(m);
    if(std::filesystem::exists(file)) throw std::invalid_argument("refusing to overwrite HDF5 mesh");
    Handle f(H5Fcreate(filename_utf8(file).c_str(),H5F_ACC_EXCL,H5P_DEFAULT,H5P_DEFAULT),H5Fclose);
    put_scalar(f,"/schema_version",1); put_scalar(f,"/dimension",static_cast<std::uint64_t>(m.dimension));
    text(f,"/length_unit",m.length_unit); text(f,"/z_convention",m.z_convention); text(f,"/provenance",m.provenance);
    static_assert(sizeof(MeshPoint)==3*sizeof(double));
    put(f,"/points",m.points.front().data(),m.points.size()*3);
    put(f,"/node_ids",m.node_ids); put(f,"/face_ids",m.face_ids); put(f,"/cell_ids",m.cell_ids);
    put(f,"/face_offsets",m.face_offsets); put(f,"/face_nodes",m.face_nodes);
    put(f,"/cell_offsets",m.cell_offsets); put(f,"/cell_faces",m.cell_faces); put(f,"/signs",m.signs);
    Handle group(H5Gcreate2(f,"/arrays",H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT),H5Gclose);
    for(const auto& [name,a]:m.arrays) {
        const auto p="/arrays/"+name;
        std::visit([&](const auto& v){put(f,p+"/values",v);},a.values);
        put_scalar(f,p+"/components",a.components); text(f,p+"/location",a.location); text(f,p+"/unit",a.unit);
    }
    check(H5Fflush(f,H5F_SCOPE_GLOBAL));
}
FaceMesh read_face_mesh_hdf5(const std::filesystem::path& file,std::uint64_t budget) {
    Handle f(H5Fopen(filename_utf8(file).c_str(),H5F_ACC_RDONLY,H5P_DEFAULT),H5Fclose);
    const std::vector<std::string> expected{"arrays","cell_faces","cell_ids","cell_offsets","dimension","face_ids","face_nodes","face_offsets","length_unit","node_ids","points","provenance","schema_version","signs","z_convention"};
    if(children(f,"/")!=expected) throw std::invalid_argument("unknown/missing HDF5 mesh schema entries");
    if(get_scalar(f,"/schema_version",budget)!=1) throw std::invalid_argument("unsupported HDF5 mesh schema version");
    FaceMesh m; const auto dim=get_scalar(f,"/dimension",budget);
    if(dim!=2 && dim!=3) throw std::invalid_argument("HDF5 invalid mesh dimension");
    m.dimension=static_cast<int>(dim);
    m.length_unit=text(f,"/length_unit",budget); m.z_convention=text(f,"/z_convention",budget); m.provenance=text(f,"/provenance",budget);
    {
        const auto v=get<double>(f,"/points",budget);
        if(v.size()%3) throw std::invalid_argument("HDF5 xyz extent");
        m.points.resize(v.size()/3);
        for(std::size_t i=0;i<m.points.size();++i) std::copy_n(v.data()+i*3,3,m.points[i].begin());
    }
    m.node_ids=get<MeshIndex>(f,"/node_ids",budget); m.face_ids=get<MeshIndex>(f,"/face_ids",budget); m.cell_ids=get<MeshIndex>(f,"/cell_ids",budget);
    m.face_offsets=get<MeshIndex>(f,"/face_offsets",budget); m.face_nodes=get<MeshIndex>(f,"/face_nodes",budget);
    m.cell_offsets=get<MeshIndex>(f,"/cell_offsets",budget); m.cell_faces=get<MeshIndex>(f,"/cell_faces",budget); m.signs=get<std::int64_t>(f,"/signs",budget);
    for(const auto& name:children(f,"/arrays")) {
        const auto p="/arrays/"+name;
        if(children(f,p)!=std::vector<std::string>{"components","location","unit","values"})
            throw std::invalid_argument("unknown attribute schema entries");
        MeshArray a; a.location=text(f,p+"/location",budget); a.unit=text(f,p+"/unit",budget);
        const auto count=get_scalar(f,p+"/components",budget);
        if(count>std::numeric_limits<std::size_t>::max()) throw std::length_error("attribute components overflow");
        a.components=static_cast<std::size_t>(count);
        Handle d(H5Dopen2(f,(p+"/values").c_str(),H5P_DEFAULT),H5Dclose);
        Handle type(H5Dget_type(d),H5Tclose);
        if(H5Tget_class(type)==H5T_FLOAT) a.values=get<double>(f,p+"/values",budget);
        else if(H5Tget_sign(type)==H5T_SGN_NONE) a.values=get<std::uint64_t>(f,p+"/values",budget);
        else a.values=get<std::int64_t>(f,p+"/values",budget);
        m.arrays.emplace(name,std::move(a));
    }
    (void)validate_face_mesh(m); return m;
}
} // namespace mpmc::mesh
