#ifndef MPMC_MESH_FACE_MESH_EXCHANGE_HPP
#define MPMC_MESH_FACE_MESH_EXCHANGE_HPP
#include <mpmc/mesh/face_mesh_vtu.hpp>
#include <mpmc/mesh/mesh_exchange_io.hpp>
#include <mpmc/mesh/computational_mesh.hpp>

namespace mpmc::mesh {
namespace face_exchange_detail {
using face_vtu_detail::Ring;
inline LinearCellType3D linear_type(std::size_t n) {
    if(n==4) return LinearCellType3D::tetrahedron;
    if(n==5) return LinearCellType3D::pyramid;
    if(n==6) return LinearCellType3D::wedge;
    if(n==8) return LinearCellType3D::hexahedron;
    throw std::invalid_argument("cell has no supported linear-element representation");
}
inline std::vector<Ring> linear_rings(const Ring& nodes,int dim) {
    std::vector<Ring> result;
    if(dim==2) {
        face_mesh_detail::require(nodes.size()==3 || nodes.size()==4,"legacy 2D export requires triangles/quadrilaterals");
        for(std::size_t i=0;i<nodes.size();++i) result.push_back({nodes[i],nodes[(i+1)%nodes.size()]});
    } else {
        LinearCell3D cell{GlobalEntityId{0},linear_type(nodes.size()),{}};
        for(auto n:nodes) {
            face_mesh_detail::require(n<=UINT32_MAX,"legacy topology index capacity");
            cell.vertices.emplace_back(static_cast<std::uint32_t>(n));
        }
        linear_cell_mesh_3d_detail::for_each_cell_face(cell,[&](auto r) {
            Ring face; for(auto n:r) face.push_back(n.value()); result.push_back(std::move(face));
        });
    }
    return result;
}
inline Ring ordered_cell(const FaceMesh& m,std::size_t c) {
    auto rings=face_vtu_detail::cell_rings(m,c);
    if(m.dimension==2) return face_vtu_detail::polygon(rings);
    Ring vertices; for(const auto& r:rings) vertices.insert(vertices.end(),r.begin(),r.end());
    std::sort(vertices.begin(),vertices.end());vertices.erase(std::unique(vertices.begin(),vertices.end()),vertices.end());
    const auto type=linear_type(vertices.size());
    const auto base_size=(type==LinearCellType3D::tetrahedron || type==LinearCellType3D::wedge)?3U:4U;
    const auto base=std::find_if(rings.begin(),rings.end(),[&](const auto& r){return r.size()==base_size;});
    face_mesh_detail::require(base!=rings.end(),"linear cell base not found");
    Ring order=*base; std::reverse(order.begin(),order.end());
    Ring outside;
    for(auto v:vertices) if(std::find(order.begin(),order.end(),v)==order.end()) outside.push_back(v);
    if(outside.size()==1) order.push_back(outside[0]);
    else {
        for(std::size_t i=0;i<base_size;++i) {
            Ring matches;
            for(const auto& r:rings) for(std::size_t j=0;j<r.size();++j) {
                const auto a=r[j],b=r[(j+1)%r.size()];
                if(a==order[i] && std::find(outside.begin(),outside.end(),b)!=outside.end()) matches.push_back(b);
                if(b==order[i] && std::find(outside.begin(),outside.end(),a)!=outside.end()) matches.push_back(a);
            }
            std::sort(matches.begin(),matches.end()); matches.erase(std::unique(matches.begin(),matches.end()),matches.end());
            face_mesh_detail::require(matches.size()==1,"cell cannot be represented by an unsplit linear element");order.push_back(matches[0]);
        }
    }
    const auto normalize=[](auto& faces) { for(auto& r:faces) r=face_vtu_detail::cycle(std::move(r)); std::sort(faces.begin(),faces.end()); };
    auto predicted=linear_rings(order,3);normalize(predicted);normalize(rings);
    face_mesh_detail::require(predicted==rings,"split/general polyhedron unsupported by legacy format subset");
    return order;
}
}

/// Import bridge for existing Gmsh/VTU/GRDECL exchange documents. Preserve
/// topology IDs and numeric fields. Original logical corner-point payload and
/// named groups are represented as typed metadata arrays, not inferred later.
[[nodiscard]] inline FaceMesh make_face_mesh(const MeshExchangeDocument& doc) {
    using namespace face_exchange_detail;
    FaceMesh m; m.dimension=doc.dimension(); m.provenance="MPMC legacy exchange document";
    const auto insert=[&](const std::string& name,MeshArray a) {
        face_mesh_detail::require(m.arrays.emplace(name,std::move(a)).second,"field name collides with exchange metadata: "+name);
    };
    const auto text=[&](const std::string& name,const std::string& value) {
        std::vector<MeshIndex> bytes;for(char ch:value) bytes.push_back(static_cast<unsigned char>(ch));
        insert(name,{"metadata","UTF8",1,std::move(bytes)});
    };
    const auto& t=doc.topology();
    face_mesh_detail::require(t.entity_count(EntityKind::edge)==0,"explicit edge entities need a separate contract");
    for(auto p:doc.vertex_coordinates_m()) m.points.push_back({p.x_m,p.y_m,p.z_m});
    for(auto id:t.global_ids(EntityKind::vertex)) m.node_ids.push_back(id.value());
    for(auto id:t.global_ids(EntityKind::face)) m.face_ids.push_back(id.value());
    for(auto id:t.global_ids(EntityKind::cell)) m.cell_ids.push_back(id.value());
    const auto& fv=t.relation(EntityKind::face,EntityKind::vertex);
    m.face_offsets.assign(fv.offsets().begin(),fv.offsets().end());
    for(auto n:fv.indices()) m.face_nodes.push_back(n.value());
    const auto& cf=t.relation(EntityKind::cell,EntityKind::face);
    m.cell_offsets.assign(cf.offsets().begin(),cf.offsets().end());
    for(auto f:cf.indices()) m.cell_faces.push_back(f.value());
    const auto& cv=t.relation(EntityKind::cell,EntityKind::vertex);
    for(std::size_t c=0;c<m.cell_ids.size();++c) {
        Ring nodes;for(auto n:cv.adjacent(LocalIndex{static_cast<std::uint32_t>(c)})) nodes.push_back(n.value());
        const auto rings=linear_rings(nodes,m.dimension);
        for(auto f:m.cell(c)) {
            int sign=0;for(const auto& r:rings) { const auto value=face_vtu_detail::orientation(m.face(static_cast<std::size_t>(f)),r);if(value) { sign=value;break; } }
            face_mesh_detail::require(sign!=0,"source document faces contradict cell order");m.signs.push_back(sign);
        }
    }
    for(const auto& field:doc.fields().fields()) {
        const auto location=field.location()==EntityKind::vertex?"node":field.location()==EntityKind::face?"face":"cell";
        face_mesh_detail::require(field.location()!=EntityKind::edge,"edge fields unsupported");
        insert(field.metadata().id,MeshArray{location,field.metadata().unit,field.component_count(),
            std::vector<double>(field.values().begin(),field.values().end())});
    }
    for(const auto& field:doc.fields().fields()) {
        const auto& meta=field.metadata();
        insert("field_source_kind_"+meta.id,{"metadata","1",1,std::vector<MeshIndex>{static_cast<MeshIndex>(meta.source.kind)}});
        text("field_source_reference_"+meta.id,meta.source.reference);
        text("field_source_revision_"+meta.id,meta.source.revision);
        text("field_source_locator_"+meta.id,meta.source.locator);
    }
    if(doc.face_boundary()) {
        std::vector<MeshIndex> tags;
        for(auto tag:doc.face_boundary()->physical_tags()) tags.push_back(tag.value());
        insert("physical_tag",{"face","1",1,std::move(tags)});
    }
    for(const auto& group:doc.groups()) {
        const auto suffix=std::to_string(static_cast<int>(group.location))+"_"+std::to_string(group.tag);
        std::vector<MeshIndex> members,bytes;
        for(auto id:group.members) members.push_back(id.value());
        for(char ch:group.name) bytes.push_back(static_cast<unsigned char>(ch));
        insert("group_members_"+suffix,{"metadata","1",1,std::move(members)});
        insert("group_name_"+suffix,{"metadata","UTF8",1,std::move(bytes)});
    }
    if(doc.logical_corner_point()) {
        const auto& q=*doc.logical_corner_point();
        insert("grdecl_dimensions",{"metadata","1",1,std::vector<MeshIndex>(q.dimensions.begin(),q.dimensions.end())});
        insert("grdecl_coord",{"metadata","m",1,q.coord_m});insert("grdecl_zcorn",{"metadata","m",1,q.zcorn_m});
        insert("grdecl_active",{"metadata","1",1,std::vector<MeshIndex>(q.active.begin(),q.active.end())});
        insert("grdecl_poro",{"metadata","1",1,q.porosity});
        insert("grdecl_permx",{"metadata","m2",1,q.permx_m2});insert("grdecl_permy",{"metadata","m2",1,q.permy_m2});insert("grdecl_permz",{"metadata","m2",1,q.permz_m2});
        insert("grdecl_scales",{"metadata","source_to_SI",2,std::vector<double>{q.source_coordinate_scale_to_m,q.source_permeability_scale_to_m2}});
    }
    (void)validate_face_mesh(m);return m;
}

/// Version-directed VTU entry: preserve existing linear MPMC face/group/field
/// semantics through their original parser; use the face schema for polygons
/// and polyhedra. Never retry a failed parse as a less strict format.
[[nodiscard]] inline FaceMesh import_face_mesh_vtu(std::string_view xml,int dimension) {
    using namespace vtu_detail;
    const auto root=require_unique_element(xml,"VTKFile","missing VTKFile");
    const auto grid=require_unique_element(root.body,"UnstructuredGrid","missing UnstructuredGrid");
    const auto piece=require_unique_element(grid.body,"Piece","single Piece required");
    bool legacy=false,face=false;
    if(const auto section=find_element(grid.body,"FieldData")) for(const auto& a:data_arrays(section->body)) {
        const auto name=require_attribute(a,"Name","FieldData name required");
        face=face || name.starts_with("mpmc_fm_");
        legacy=legacy || name.starts_with("mpmc_group_") || name.starts_with("mpmc_face_") || name=="mpmc_global_face_id";
    }
    for(const auto name:{"PointData","CellData"}) if(const auto section=find_element(piece.body,name))
        for(const auto& a:data_arrays(section->body)) legacy=legacy || a.attributes.contains("mpmc_unit");
    face_mesh_detail::require(!(face && legacy),"mixed incompatible VTU metadata schemas");
    if(!legacy) return read_face_mesh_vtu(xml,dimension);
    face_mesh_detail::require(dimension==2 || dimension==3,"VTU dimension must be explicit");
    return dimension==2?make_face_mesh(make_mesh_exchange_document(import_vtu_ascii(xml))):
        make_face_mesh(make_mesh_exchange_document(import_vtu_ascii_3d(xml)));
}

/// Bounded legacy export, with explicit refusal for general/split polyhedra.
/// HDF5/face-VTU remain the lossless routes; this operation reports metadata
/// omission even when coordinates and connectivity fit the legacy subset.
[[nodiscard]] inline MeshTextExportResult export_face_mesh_legacy(const FaceMesh& m,MeshExchangeFormat format) {
    ConversionReport rejected{format};
    try {
        face_mesh_detail::require(format==MeshExchangeFormat::gmsh_4_1_ascii || format==MeshExchangeFormat::grdecl,"legacy target must be Gmsh or GRDECL");
        const auto source=face_mesh_topology(m);
        Topology::EntityIds ids;
        for(auto x:m.node_ids) ids.vertices.emplace_back(x);
        for(auto x:m.face_ids) ids.faces.emplace_back(x);
        for(auto x:m.cell_ids) ids.cells.emplace_back(x);
        std::vector<CsrAdjacency> relations;
        relations.push_back(source.relation(EntityKind::cell,EntityKind::face));
        relations.push_back(source.relation(EntityKind::face,EntityKind::vertex));
        relations.push_back(source.relation(EntityKind::face,EntityKind::cell));
        std::vector<CsrAdjacency::Offset> offsets{0};std::vector<LocalIndex> vertices;
        for(std::size_t c=0;c<m.cell_ids.size();++c) {
            for(auto n:face_exchange_detail::ordered_cell(m,c)) {
                face_mesh_detail::require(n<=UINT32_MAX,"legacy index capacity");vertices.emplace_back(static_cast<std::uint32_t>(n));
            }
            face_mesh_detail::require(vertices.size()<=UINT32_MAX,"legacy CSR capacity");offsets.push_back(static_cast<std::uint32_t>(vertices.size()));
        }
        relations.emplace_back(EntityKind::cell,EntityKind::vertex,m.points.size(),std::move(offsets),std::move(vertices));
        Topology topology{std::move(ids),std::move(relations)};
        std::vector<Coordinate3D> coordinates;for(auto p:m.points) coordinates.push_back({p[0],p[1],p[2]});
        std::vector<DenseFieldSnapshot> fields;
        for(const auto& [name,a]:m.arrays) if(const auto* v=std::get_if<std::vector<double>>(&a.values)) {
            if(a.location=="node" || a.location=="face" || a.location=="cell") fields.push_back(DenseFieldSnapshot::create(
                topology,a.location=="node"?EntityKind::vertex:a.location=="face"?EntityKind::face:EntityKind::cell,
                a.components,*v,{name,a.unit,{FieldSourceKind::user_supplied,
                    m.provenance.empty()?"FaceMesh caller-supplied data":m.provenance,"1",name}}));
        }
        const auto doc=MeshExchangeDocument::create(MeshExchangeFormat::face_based,m.dimension,std::move(topology),std::move(coordinates),std::nullopt,std::move(fields),{},std::nullopt);
        if(m.dimension==2) (void)prepare_linear_mesh_2d(doc); else (void)prepare_linear_mesh_3d(doc);
        auto result=format==MeshExchangeFormat::gmsh_4_1_ascii?export_gmsh_4_1_ascii(doc):export_grdecl_ascii(doc);
        if(result.exported()) result.report.note_lossy("face_mesh_metadata","Legacy format subset does not preserve the complete face-mesh schema, typed attributes, geometry policy and provenance; use HDF5 or face-VTU for lossless exchange.");
        return result;
    } catch(const std::invalid_argument& e) {
        rejected.note_unsupported("linear_subset_incompatible",e.what());return {std::nullopt,std::move(rejected)};
    }
}
} // namespace mpmc::mesh
#endif
