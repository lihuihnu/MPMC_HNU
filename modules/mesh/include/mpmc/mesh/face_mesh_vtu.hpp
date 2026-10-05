#ifndef MPMC_MESH_FACE_MESH_VTU_HPP
#define MPMC_MESH_FACE_MESH_VTU_HPP
#include <mpmc/mesh/face_mesh.hpp>
#include <mpmc/mesh/vtu.hpp>
#include <mpmc/mesh/linear_cell_mesh_3d.hpp>
#include <ostream>

namespace mpmc::mesh {
namespace face_vtu_detail {
using Ring=std::vector<MeshIndex>;
inline Ring cycle(Ring r) {
    if(!r.empty()) std::rotate(r.begin(),std::min_element(r.begin(),r.end()),r.end());
    return r;
}
inline int orientation(std::span<const MeshIndex> a,std::span<const MeshIndex> b) {
    if(a.size()!=b.size()) return 0;
    if(a.size()==2) return a[0]==b[0] && a[1]==b[1]?1:a[0]==b[1] && a[1]==b[0]?-1:0;
    const auto c=cycle(Ring(a.begin(),a.end()));
    Ring d(b.begin(),b.end());
    if(c==cycle(d)) return 1;
    std::reverse(d.begin(),d.end()); return c==cycle(d)?-1:0;
}
template<class V> void array(std::ostream& o,const std::string& name,const V& values,std::size_t components=1) {
    using T=typename V::value_type;
    const char* type=std::is_floating_point_v<T>?"Float64":std::is_signed_v<T>?"Int64":"UInt64";
    o<<"<DataArray type=\""<<type<<"\" Name=\""<<vtu_detail::xml_escape(name)
     <<"\" NumberOfComponents=\""<<components<<"\" NumberOfTuples=\""<<values.size()/components<<"\" format=\"ascii\">\n";
    for(auto x:values) o<<x<<' ';
    o<<"\n</DataArray>\n";
}
inline void text(std::ostream& o,const std::string& name,const std::string& value) {
    std::vector<std::uint64_t> bytes; for(char c:value) bytes.push_back(static_cast<unsigned char>(c)); array(o,name,bytes);
}
inline std::vector<Ring> cell_rings(const FaceMesh& m,std::size_t c) {
    std::vector<Ring> r; const auto faces=m.cell(c);
    for(std::size_t j=0;j<faces.size();++j) {
        const auto face=m.face(static_cast<std::size_t>(faces[j])); r.emplace_back(face.begin(),face.end());
        if(m.signs[static_cast<std::size_t>(m.cell_offsets[c])+j]<0) std::reverse(r.back().begin(),r.back().end());
    }
    return r;
}
inline Ring polygon(const std::vector<Ring>& edges) {
    std::map<MeshIndex,MeshIndex> next;
    for(const auto& e:edges) face_mesh_detail::require(e.size()==2 && next.emplace(e[0],e[1]).second,"invalid polygon ring");
    Ring result; auto node=next.begin()->first;
    for(std::size_t i=0;i<edges.size();++i) { result.push_back(node); node=next.at(node); }
    face_mesh_detail::require(node==result.front(),"open polygon"); return result;
}
}

/// Standard VTK_POLYGON / VTK_POLYHEDRON plus dataset-level schema arrays.
/// Writes to a stream: full meshes need not form one giant XML string.
/// Geometry is transported as declared, not silently repaired/triangulated.
inline void write_face_mesh_vtu(std::ostream& out,const FaceMesh& m) {
    using namespace face_vtu_detail;
    (void)validate_face_mesh(m);
    for(const auto& [name,a]:m.arrays)
        face_mesh_detail::require(!((a.location=="node" && name=="mpmc_global_vertex_id") ||
            (a.location=="cell" && name=="mpmc_global_cell_id")),"attribute collides with VTU stable ID");
    out<<std::setprecision(std::numeric_limits<double>::max_digits10);
    out<<"<?xml version=\"1.0\"?>\n<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\"><UnstructuredGrid>\n<FieldData>\n";
    array(out,"mpmc_fm_version",std::vector<MeshIndex>{1});
    array(out,"mpmc_fm_dimension",std::vector<MeshIndex>{static_cast<MeshIndex>(m.dimension)});
    text(out,"mpmc_fm_z",m.z_convention); text(out,"mpmc_fm_provenance",m.provenance);
    array(out,"mpmc_fm_face_ids",m.face_ids); array(out,"mpmc_fm_face_offsets",m.face_offsets); array(out,"mpmc_fm_face_nodes",m.face_nodes);
    array(out,"mpmc_fm_cell_offsets",m.cell_offsets); array(out,"mpmc_fm_cell_faces",m.cell_faces); array(out,"mpmc_fm_signs",m.signs);
    for(const auto& [name,a]:m.arrays) {
        std::visit([&](const auto& v){array(out,"mpmc_fm_values_"+name,v,a.components);},a.values);
        text(out,"mpmc_fm_location_"+name,a.location); text(out,"mpmc_fm_unit_"+name,a.unit);
    }
    out<<"</FieldData><Piece NumberOfPoints=\""<<m.points.size()<<"\" NumberOfCells=\""<<m.cell_ids.size()<<"\">\n<Points><DataArray type=\"Float64\" NumberOfComponents=\"3\" format=\"ascii\">\n";
    for(auto p:m.points) out<<p[0]<<' '<<p[1]<<' '<<p[2]<<'\n';
    out<<"</DataArray></Points><PointData>\n"; array(out,"mpmc_global_vertex_id",m.node_ids);
    for(const auto& [name,a]:m.arrays) if(a.location=="node") std::visit([&](const auto& v){array(out,name,v,a.components);},a.values);
    out<<"</PointData><CellData>\n"; array(out,"mpmc_global_cell_id",m.cell_ids);
    for(const auto& [name,a]:m.arrays) if(a.location=="cell") std::visit([&](const auto& v){array(out,name,v,a.components);},a.values);
    out<<"</CellData><Cells>\n<DataArray type=\"Int64\" Name=\"connectivity\" format=\"ascii\">\n";
    std::vector<MeshIndex> offsets; offsets.reserve(m.cell_ids.size()); MeshIndex end=0;
    for(std::size_t c=0;c<m.cell_ids.size();++c) {
        const auto rings=cell_rings(m,c); Ring nodes;
        if(m.dimension==2) nodes=polygon(rings);
        else {
            for(const auto& r:rings) nodes.insert(nodes.end(),r.begin(),r.end());
            std::sort(nodes.begin(),nodes.end()); nodes.erase(std::unique(nodes.begin(),nodes.end()),nodes.end());
        }
        for(auto n:nodes) out<<n<<' ';
        out<<'\n'; end+=nodes.size(); offsets.push_back(end);
    }
    out<<"</DataArray>\n"; array(out,"offsets",offsets);
    out<<"<DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">\n";
    for(std::size_t c=0;c<m.cell_ids.size();++c) out<<(m.dimension==2?7:42)<<' ';
    out<<"\n</DataArray>\n";
    if(m.dimension==3) {
        out<<"<DataArray type=\"Int64\" Name=\"faces\" format=\"ascii\">\n"; end=0;
        for(std::size_t c=0;c<m.cell_ids.size();++c) {
            const auto rings=cell_rings(m,c); out<<rings.size()<<' '; ++end;
            for(const auto& r:rings) { out<<r.size()<<' '; for(auto n:r) out<<n<<' '; end+=r.size()+1; }
            offsets[c]=end; out<<'\n';
        }
        out<<"</DataArray>\n"; array(out,"faceoffsets",offsets);
    }
    out<<"</Cells></Piece></UnstructuredGrid></VTKFile>\n";
    if(!out) throw std::runtime_error("face mesh VTU write failed");
}

[[nodiscard]] inline FaceMesh read_face_mesh_vtu(std::string_view xml,int dimension) {
    using namespace vtu_detail; using namespace face_vtu_detail;
    using face_mesh_detail::require;
    require(dimension==2 || dimension==3,"VTU dimension must be explicit");
    const auto root=require_unique_element(xml,"VTKFile","missing VTKFile");
    require(require_attribute(root,"type","missing VTK type")=="UnstructuredGrid","VTU type unsupported");
    require(!root.attributes.contains("compressor") && !find_element(root.body,"AppendedData"),"compressed/appended VTU unsupported");
    const auto grid=require_unique_element(root.body,"UnstructuredGrid","missing UnstructuredGrid");
    const auto piece=require_unique_element(grid.body,"Piece","single Piece required");
    const auto np=parse_size(require_attribute(piece,"NumberOfPoints","point count required"),"invalid point count");
    const auto nc=parse_size(require_attribute(piece,"NumberOfCells","cell count required"),"invalid cell count");
    const auto points=data_arrays(require_unique_element(piece.body,"Points","missing Points").body);
    require(points.size()==1 && component_count(points[0])==3,"VTU xyz required"); require_ascii_data_array(points[0]);
    const auto coordinate_type=require_attribute(points[0],"type","coordinate type required");
    require(coordinate_type=="Float32" || coordinate_type=="Float64","VTU coordinates must be floating point");
    const auto coords=parse_double_values(points[0].body,"invalid VTU coordinates");
    require(coords.size()%3==0 && coords.size()/3==np,"VTU point extent mismatch");
    FaceMesh m; m.dimension=dimension; m.provenance="external VTU ASCII, explicit SI metres";
    for(std::size_t i=0;i<np;++i) m.points.push_back({coords[3*i],coords[3*i+1],coords[3*i+2]});
    m.node_ids.resize(np); m.cell_ids.resize(nc); std::iota(m.node_ids.begin(),m.node_ids.end(),0ULL); std::iota(m.cell_ids.begin(),m.cell_ids.end(),0ULL);
    std::map<std::string,XmlElement> cell_arrays,metadata;
    const auto collect=[&](const XmlElement& section,auto& result) {
        for(const auto& a:data_arrays(section.body)) { require_ascii_data_array(a); require(result.emplace(require_attribute(a,"Name","VTU array name required"),a).second,"duplicate VTU array"); }
    };
    collect(require_unique_element(piece.body,"Cells","missing Cells"),cell_arrays);
    if(find_element(grid.body,"FieldData")) collect(require_unique_element(grid.body,"FieldData","unique FieldData required"),metadata);
    require(!find_element(piece.body,"FieldData"),"FieldData must be dataset-level");
    const auto unsigned_values=[&](const XmlElement& a) {
        const auto type=require_attribute(a,"type","array type required");
        require(type=="Int64" || type=="UInt64" || type=="Int32" || type=="UInt32" || type=="UInt8","integer VTU array required");
        return parse_integer_values<MeshIndex>(a.body,"invalid/nonpositive VTU integer");
    };
    const auto typed=[&](const XmlElement& a)->MeshArrayValues {
        const auto type=require_attribute(a,"type","array type required");
        if(type=="Float64" || type=="Float32") return parse_double_values(a.body,"invalid field");
        if(type=="UInt64" || type=="UInt32" || type=="UInt8") return unsigned_values(a);
        require(type=="Int64" || type=="Int32","unsupported field type");
        return parse_integer_values<std::int64_t>(a.body,"invalid field integer");
    };
    for(const auto& [section,location,idname]:std::vector<std::array<std::string,3>>{{"PointData","node","mpmc_global_vertex_id"},{"CellData","cell","mpmc_global_cell_id"}}) {
        bool seen_id=false;
        if(find_element(piece.body,section)) for(const auto& a:data_arrays(require_unique_element(piece.body,section,"unique data section required").body)) {
            require_ascii_data_array(a); const auto name=require_attribute(a,"Name","field name required");
            if(name==idname) { require(!seen_id && component_count(a)==1,"duplicate/non-scalar stable ID array"); seen_id=true; if(location=="node") m.node_ids=unsigned_values(a); else m.cell_ids=unsigned_values(a); }
            else require(m.arrays.emplace(name,MeshArray{location,"unspecified",component_count(a),typed(a)}).second,"duplicate VTU field");
        }
    }
    require(m.node_ids.size()==np && m.cell_ids.size()==nc,"stable ID count mismatch");
    const auto conn=unsigned_values(cell_arrays.at("connectivity")), offsets=unsigned_values(cell_arrays.at("offsets")), types=unsigned_values(cell_arrays.at("types"));
    require(offsets.size()==nc && types.size()==nc,"VTU cell array extent");
    std::vector<MeshIndex> stream; std::vector<std::int64_t> faceends;
    if(cell_arrays.contains("faces")) {
        stream=unsigned_values(cell_arrays.at("faces"));
        faceends=parse_integer_values<std::int64_t>(cell_arrays.at("faceoffsets").body,"invalid face offsets");
        require(faceends.size()==nc,"faceoffsets extent");
    } else if(cell_arrays.contains("face_connectivity")) {
        // VTK 9.4+ also writes explicit face tables. Normalize their CSR
        // representation into the legacy stream with all bounds checked.
        const auto fn=unsigned_values(cell_arrays.at("face_connectivity"));
        const auto fo=unsigned_values(cell_arrays.at("face_offsets"));
        const auto cf=unsigned_values(cell_arrays.at("polyhedron_to_faces"));
        const auto co=unsigned_values(cell_arrays.at("polyhedron_offsets"));
        require(co.size()==nc && std::is_sorted(co.begin(),co.end()) &&
                (co.empty()?cf.empty():co.back()==cf.size()),"invalid polyhedron table offsets");
        require(std::is_sorted(fo.begin(),fo.end()) &&
                (fo.empty()?fn.empty():fo.back()==fn.size()),"invalid face table offsets");
        MeshIndex first=0;
        for(std::size_t c=0;c<nc;++c) {
            if(types[c]!=42) { require(co[c]==first,"nonpolyhedron has explicit faces"); faceends.push_back(-1); continue; }
            stream.push_back(co[c]-first);
            for(auto j=first;j<co[c];++j) {
                const auto face=cf[static_cast<std::size_t>(j)]; require(face<fo.size(),"invalid polyhedron face reference");
                const auto lo=face==0?0:fo[static_cast<std::size_t>(face)-1],hi=fo[static_cast<std::size_t>(face)];
                stream.push_back(hi-lo);
                stream.insert(stream.end(),fn.begin()+static_cast<std::ptrdiff_t>(lo),fn.begin()+static_cast<std::ptrdiff_t>(hi));
            }
            first=co[c]; faceends.push_back(static_cast<std::int64_t>(stream.size()));
        }
    }
    require(cell_arrays.size()==(cell_arrays.contains("faces")?5U:cell_arrays.contains("face_connectivity")?7U:3U),"uncontracted VTU Cells arrays");
    std::map<Ring,std::size_t> face_map;
    MeshIndex begin=0, facebegin=0;
    const bool own=metadata.contains("mpmc_fm_version");
    FaceMesh original;
    if(own) {
        require(unsigned_values(metadata.at("mpmc_fm_version"))==std::vector<MeshIndex>{1},"unsupported face VTU schema");
        require(unsigned_values(metadata.at("mpmc_fm_dimension"))==std::vector<MeshIndex>{static_cast<MeshIndex>(dimension)},"VTU dimension mismatch");
        original.dimension=dimension; original.points=m.points; original.node_ids=m.node_ids; original.cell_ids=m.cell_ids;
        original.face_ids=unsigned_values(metadata.at("mpmc_fm_face_ids")); original.face_offsets=unsigned_values(metadata.at("mpmc_fm_face_offsets")); original.face_nodes=unsigned_values(metadata.at("mpmc_fm_face_nodes"));
        original.cell_offsets=unsigned_values(metadata.at("mpmc_fm_cell_offsets")); original.cell_faces=unsigned_values(metadata.at("mpmc_fm_cell_faces"));
        original.signs=parse_integer_values<std::int64_t>(metadata.at("mpmc_fm_signs").body,"invalid signs");
        const auto string_value=[&](const std::string& name) { std::string value; for(auto x:unsigned_values(metadata.at(name))) { require(x<256,"invalid UTF8 byte"); value.push_back(static_cast<char>(x)); } return value; };
        original.z_convention=string_value("mpmc_fm_z"); original.provenance=string_value("mpmc_fm_provenance");
        for(const auto& [name,a]:metadata) if(name.starts_with("mpmc_fm_values_")) {
            const auto key=name.substr(15); MeshArray value{string_value("mpmc_fm_location_"+key),string_value("mpmc_fm_unit_"+key),component_count(a),typed(a)};
            if(value.location=="node" || value.location=="cell") {
                const auto it=m.arrays.find(key); require(it!=m.arrays.end() && it->second.location==value.location && it->second.components==value.components && it->second.values==value.values,"VTU visible field contradicts schema");
                m.arrays.erase(it);
            }
            original.arrays.emplace(key,std::move(value));
        }
        require(metadata.size()==10+3*original.arrays.size(),"uncontracted canonical VTU metadata");
        require(m.arrays.empty(),"uncontracted visible VTU fields");
        (void)validate_face_mesh(original);
    } else require(metadata.empty(),"uncontracted VTU FieldData; preserve through explicit schema");
    for(std::size_t c=0;c<nc;++c) {
        require(offsets[c]>=begin && offsets[c]<=conn.size(),"invalid cell offset");
        Ring nodes(conn.begin()+static_cast<std::ptrdiff_t>(begin),conn.begin()+static_cast<std::ptrdiff_t>(offsets[c])); begin=offsets[c];
        for(auto n:nodes) require(n<np,"VTU node out of range");
        std::vector<Ring> rings;
        if(dimension==2) {
            require(types[c]==5 || types[c]==7 || types[c]==9,"unsupported 2D VTU type");
            require(nodes.size()>=3 && (types[c]!=5 || nodes.size()==3) && (types[c]!=9 || nodes.size()==4),"invalid polygon arity");
            for(std::size_t j=0;j<nodes.size();++j) rings.push_back({nodes[j],nodes[(j+1)%nodes.size()]});
        } else if(types[c]==42) {
            require(c<faceends.size() && faceends[c]>=0,"missing polyhedron faceoffset");
            const auto stop=static_cast<MeshIndex>(faceends[c]); require(stop<=stream.size() && stop>facebegin,"invalid polyhedron stream");
            const auto count=stream[static_cast<std::size_t>(facebegin++)]; require(count<=stop-facebegin,"invalid face count");
            for(MeshIndex j=0;j<count;++j) {
                require(facebegin<stop,"truncated face stream"); const auto n=stream[static_cast<std::size_t>(facebegin++)]; require(n>=3 && n<=stop-facebegin,"invalid face length");
                rings.emplace_back(stream.begin()+static_cast<std::ptrdiff_t>(facebegin),stream.begin()+static_cast<std::ptrdiff_t>(facebegin+n)); facebegin+=n;
            }
            require(facebegin==stop,"trailing face stream");
            Ring used; for(const auto& r:rings) used.insert(used.end(),r.begin(),r.end());
            std::sort(used.begin(),used.end()); used.erase(std::unique(used.begin(),used.end()),used.end());
            auto listed=nodes; std::sort(listed.begin(),listed.end()); require(used==listed,"polyhedron connectivity contradicts faces");
        } else {
            LinearCellType3D type;
            if(types[c]==10) type=LinearCellType3D::tetrahedron;
            else if(types[c]==12) type=LinearCellType3D::hexahedron;
            else if(types[c]==13) type=LinearCellType3D::wedge;
            else { require(types[c]==14,"unsupported 3D VTU cell"); type=LinearCellType3D::pyramid; }
            LinearCell3D cell{GlobalEntityId{m.cell_ids[c]},type,{}};
            for(auto n:nodes) { require(n<=std::numeric_limits<std::uint32_t>::max(),"VTU linear index capacity"); cell.vertices.emplace_back(static_cast<std::uint32_t>(n)); }
            linear_cell_mesh_3d_detail::for_each_cell_face(cell,[&](auto r) { Ring v; for(auto n:r) v.push_back(n.value()); rings.push_back(std::move(v)); });
        }
        if(own) {
            auto expected=cell_rings(original,c);
            const auto normalize=[](auto& rr) { for(auto& r:rr) if(r.size()>2) r=cycle(std::move(r)); std::sort(rr.begin(),rr.end()); };
            normalize(expected); normalize(rings); require(expected==rings,"VTU cell faces contradict canonical schema");
        } else for(const auto& r:rings) {
            auto key=r; std::sort(key.begin(),key.end());
            const auto [it,inserted]=face_map.emplace(std::move(key),m.face_ids.size());
            if(inserted) { m.face_ids.push_back(it->second); m.face_nodes.insert(m.face_nodes.end(),r.begin(),r.end()); m.face_offsets.push_back(m.face_nodes.size()); }
            const auto sign=orientation(m.face(it->second),r); require(sign!=0,"shared face cycles disagree");
            m.cell_faces.push_back(it->second); m.signs.push_back(sign);
        }
        if(!own) m.cell_offsets.push_back(m.cell_faces.size());
    }
    require(begin==conn.size() && facebegin==stream.size(),"unused VTU connectivity");
    if(own) return original;
    (void)validate_face_mesh(m); return m;
}
} // namespace mpmc::mesh
#endif
