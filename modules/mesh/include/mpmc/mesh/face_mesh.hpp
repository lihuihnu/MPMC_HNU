#ifndef MPMC_MESH_FACE_MESH_HPP
#define MPMC_MESH_FACE_MESH_HPP

#include <mpmc/mesh/topology.hpp>
#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <limits>
#include <map>
#include <numeric>
#include <span>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <utility>
#include <variant>
#include <vector>

namespace mpmc::mesh {
using MeshIndex = std::uint64_t;
using MeshPoint = std::array<double, 3>;
using MeshArrayValues = std::variant<std::vector<double>, std::vector<std::int64_t>,
                                     std::vector<std::uint64_t>>;
/// Row-major typed attribute, associated with node, face, cell, incidence, nnc,
/// or metadata. Integer values never pass through double. Units are explicit.
struct MeshArray {
    std::string location;
    std::string unit;
    std::size_t components{1};
    MeshArrayValues values;
    [[nodiscard]] std::size_t size() const {
        return std::visit([](const auto& v) { return v.size(); }, values);
    }
};

/// Source-preserving, face-based interchange storage. Zero-based references;
/// stable IDs are independent of array positions. A face ring's right-hand
/// normal (right normal of a 2D edge) is outward for incidence sign +1.
/// Ordinary connected closed manifold cells only; embedded surfaces excluded.
struct FaceMesh {
    int dimension{3};
    std::string length_unit{"m"};
    std::string z_convention{"depth"};
    std::string provenance;
    std::vector<MeshPoint> points;
    std::vector<MeshIndex> node_ids, face_ids, cell_ids;
    std::vector<MeshIndex> face_offsets{0}, face_nodes;
    std::vector<MeshIndex> cell_offsets{0}, cell_faces;
    std::vector<std::int64_t> signs;
    // Optional arrays include mrst_index_map (1-based external numbering),
    // cart_dims, cell_face_tags, reference geometry, material fields, and NNC.
    std::map<std::string, MeshArray> arrays;
    [[nodiscard]] std::span<const MeshIndex> face(std::size_t f) const {
        return std::span<const MeshIndex>(face_nodes).subspan(
            static_cast<std::size_t>(face_offsets.at(f)),
            static_cast<std::size_t>(face_offsets.at(f+1)-face_offsets.at(f)));
    }
    [[nodiscard]] std::span<const MeshIndex> cell(std::size_t c) const {
        return std::span<const MeshIndex>(cell_faces).subspan(
            static_cast<std::size_t>(cell_offsets.at(c)),
            static_cast<std::size_t>(cell_offsets.at(c+1)-cell_offsets.at(c)));
    }
};

namespace face_mesh_detail {
inline void require(bool valid, const std::string& message) {
    if (!valid) throw std::invalid_argument("face mesh: " + message);
}
inline void ids(const std::vector<MeshIndex>& v) {
    if (std::is_sorted(v.begin(),v.end()) && std::adjacent_find(v.begin(),v.end())==v.end()) return;
    auto sorted=v;
    std::sort(sorted.begin(),sorted.end());
    require(std::adjacent_find(sorted.begin(),sorted.end())==sorted.end(), "duplicate stable ID");
}
inline void csr(const std::vector<MeshIndex>& o, const std::vector<MeshIndex>& v,
                std::size_t rows, std::size_t targets) {
    require(o.size()==rows+1 && o.front()==0 && o.back()==v.size(), "invalid CSR extent");
    require(std::is_sorted(o.begin(),o.end()), "nonmonotone CSR offsets");
    for (auto x:v) require(x<targets,"CSR index out of range");
}
inline MeshPoint add(MeshPoint a, MeshPoint b) { for(std::size_t k=0;k<3;++k) a[k]+=b[k]; return a; }
inline MeshPoint sub(MeshPoint a, MeshPoint b) { for(std::size_t k=0;k<3;++k) a[k]-=b[k]; return a; }
inline MeshPoint mul(MeshPoint a, double s) { for(auto& x:a) x*=s; return a; }
inline double dot(MeshPoint a, MeshPoint b) { return a[0]*b[0]+a[1]*b[1]+a[2]*b[2]; }
inline MeshPoint cross(MeshPoint a, MeshPoint b) {
    return {a[1]*b[2]-a[2]*b[1],a[2]*b[0]-a[0]*b[2],a[0]*b[1]-a[1]*b[0]};
}
inline double norm(MeshPoint a) { return std::sqrt(dot(a,a)); }

inline void validate_reference_geometry_array_contract(
    const std::string& name, const MeshArray& a, int dimension) {
    const auto require_contract=[&](const char* location,std::size_t components,const char* unit) {
        require(a.location==location && a.components==components && a.unit==unit &&
                std::holds_alternative<std::vector<double>>(a.values),
                "invalid reference geometry array contract for "+name);
    };
    if(name=="reference_face_areas")
        require_contract("face",1,dimension==2 ? "m" : "m2");
    else if(name=="reference_face_normals")
        require_contract("face",3,dimension==2 ? "m" : "m2");
    else if(name=="reference_face_centroids")
        require_contract("face",3,"m");
    else if(name=="reference_cell_volumes")
        require_contract("cell",1,dimension==2 ? "m2" : "m3");
    else if(name=="reference_cell_centroids")
        require_contract("cell",3,"m");
}
}

/// Checks transport invariants. Optional closed-boundary certification is a
/// computational gate; preserving a source mesh need not certify its geometry.
/// Without it, no claim is made about positive volume or numerical
/// method admissibility. Returns signed face neighbors; -1 denotes exterior.
[[nodiscard]] inline std::vector<std::array<std::int64_t,2>> validate_face_mesh(const FaceMesh& m, bool require_closed_boundary=false) {
    using namespace face_mesh_detail;
    require(m.dimension==2 || m.dimension==3,"dimension must be 2 or 3");
    require(m.length_unit=="m", "coordinates must be SI metres");
    require(m.z_convention=="depth" || m.z_convention=="elevation","unknown Z convention");
    require(m.points.size()==m.node_ids.size(),"point/ID extent mismatch");
    require(!m.points.empty() && !m.cell_ids.empty(),"empty mesh");
    for(const auto p:m.points) for(std::size_t k=0;k<3;++k)
        require(std::isfinite(p[k]) && (m.dimension!=2 || k!=2 || p[k]==0),"invalid coordinate");
    ids(m.node_ids); ids(m.face_ids); ids(m.cell_ids);
    csr(m.face_offsets,m.face_nodes,m.face_ids.size(),m.points.size());
    csr(m.cell_offsets,m.cell_faces,m.cell_ids.size(),m.face_ids.size());
    require(m.signs.size()==m.cell_faces.size(),"incidence sign extent mismatch");
    std::vector<std::array<std::int64_t,2>> neighbors(m.face_ids.size(),{-1,-1});
    for(std::size_t f=0;f<m.face_ids.size();++f) {
        auto ring=m.face(f);
        require(m.dimension==2 ? ring.size()==2 : ring.size()>=3,"invalid face arity");
        auto v=std::vector<MeshIndex>(ring.begin(),ring.end());
        std::sort(v.begin(),v.end());
        require(std::adjacent_find(v.begin(),v.end())==v.end(),"repeated vertex in face");
    }
    for(std::size_t c=0;c<m.cell_ids.size();++c) {
        const auto faces=m.cell(c);
        require(faces.size()>=static_cast<std::size_t>(m.dimension+1),"cell has too few faces");
        std::vector<std::array<MeshIndex,4>> edges;
        std::vector<MeshIndex> sorted_faces(faces.begin(),faces.end());
        std::sort(sorted_faces.begin(),sorted_faces.end());
        require(std::adjacent_find(sorted_faces.begin(),sorted_faces.end())==sorted_faces.end(),"repeated face in cell");
        for(std::size_t j=0;j<faces.size();++j) {
            const auto f=static_cast<std::size_t>(faces[j]);
            const auto sign=m.signs[static_cast<std::size_t>(m.cell_offsets[c])+j];
            require(sign==1 || sign==-1,"incidence sign must be +/-1");
            auto& owner=neighbors[f][sign==1 ? 0U:1U];
            require(owner==-1,"nonmanifold face or equal shared-face orientation");
            owner=static_cast<std::int64_t>(c);
            const auto ring=m.face(f);
            if(!require_closed_boundary) continue;
            if(m.dimension==2) {
                const auto a=ring[sign==1?0U:1U], b=ring[sign==1?1U:0U];
                edges.push_back({a,a,0,j}); edges.push_back({b,b,1,j});
            } else {
                for(std::size_t k=0;k<ring.size();++k) {
                    auto a=ring[k], b=ring[(k+1)%ring.size()];
                    if(sign<0) std::swap(a,b);
                    edges.push_back({std::min(a,b),std::max(a,b),a<b?0U:1U,j});
                }
            }
        }
        if(!require_closed_boundary) continue;
        const auto closed=[&](auto& boundary) {
            std::sort(boundary.begin(),boundary.end());
            if(boundary.size()%2) return false;
            std::vector<std::size_t> parent(faces.size());
            std::iota(parent.begin(),parent.end(),0U);
            const auto root=[&parent](std::size_t x) { while(parent[x]!=x) x=parent[x]; return x; };
            for(std::size_t e=0;e<boundary.size();e+=2) {
                const auto& a=boundary[e]; const auto& b=boundary[e+1];
                if(!(a[0]==b[0] && a[1]==b[1] && a[2]!=b[2] &&
                    (e+2==boundary.size() || boundary[e+2][0]!=a[0] || boundary[e+2][1]!=a[1]))) return false;
                parent[root(static_cast<std::size_t>(a[3]))]=root(static_cast<std::size_t>(b[3]));
            }
            for(std::size_t j=1;j<faces.size();++j) if(root(j)!=root(0)) return false;
            return true;
        };
        if(!closed(edges)) {
            require(m.dimension==3,"open polygon boundary");
            // MRST fault faces may subdivide the same geometrical cell edge
            // differently. Construct a virtual common refinement for validation
            // only. Original face rings, IDs and metrics remain untouched.
            std::vector<MeshIndex> vertices;
            for(auto f:faces) { const auto r=m.face(static_cast<std::size_t>(f)); vertices.insert(vertices.end(),r.begin(),r.end()); }
            std::sort(vertices.begin(),vertices.end()); vertices.erase(std::unique(vertices.begin(),vertices.end()),vertices.end());
            double coordinate_scale=1, length_scale=0;
            const auto origin=m.points[static_cast<std::size_t>(vertices.front())];
            for(auto v:vertices) {
                const auto p=m.points[static_cast<std::size_t>(v)];
                length_scale=std::max(length_scale,norm(sub(p,origin)));
                for(auto x:p) coordinate_scale=std::max(coordinate_scale,std::abs(x));
            }
            const double tolerance=std::numeric_limits<double>::epsilon()*(64*coordinate_scale+4096*length_scale);
            std::vector<std::array<MeshIndex,4>> refined;
            for(const auto& e:edges) {
                const auto a=m.points[static_cast<std::size_t>(e[0])], b=m.points[static_cast<std::size_t>(e[1])];
                const auto d=sub(b,a); const double square=dot(d,d);
                require(square>tolerance*tolerance,"collapsed cell edge");
                std::vector<std::pair<double,MeshIndex>> cuts{{0,e[0]},{1,e[1]}};
                for(auto v:vertices) if(v!=e[0] && v!=e[1]) {
                    const auto q=sub(m.points[static_cast<std::size_t>(v)],a);
                    const double t=dot(q,d)/square;
                    if(t>0 && t<1 && norm(sub(q,mul(d,t)))<=tolerance) cuts.emplace_back(t,v);
                }
                std::sort(cuts.begin(),cuts.end());
                for(std::size_t k=1;k<cuts.size();++k) {
                    auto u=cuts[k-1].second,v=cuts[k].second;
                    if(e[2]) std::swap(u,v);
                    refined.push_back({std::min(u,v),std::max(u,v),u<v?0U:1U,e[3]});
                }
            }
            require(closed(refined),"cell boundary does not close after collinear edge refinement at cell "+std::to_string(c));
        }
    }
    for(auto n:neighbors) require(n[0]>=0 || n[1]>=0,"orphan face");
    std::size_t nnc_count=0;
    if(auto it=m.arrays.find("nnc_cells");it!=m.arrays.end()) {
        const auto& a=it->second;
        require(a.location=="nnc" && a.components==2 && a.size()%2==0 &&
            std::holds_alternative<std::vector<MeshIndex>>(a.values),"invalid NNC endpoint array");
        const auto& endpoints=std::get<std::vector<MeshIndex>>(a.values);
        nnc_count=endpoints.size()/2;
        for(std::size_t i=0;i<endpoints.size();i+=2)
            require(endpoints[i]<m.cell_ids.size() && endpoints[i+1]<m.cell_ids.size() &&
                    endpoints[i]!=endpoints[i+1],"invalid NNC endpoints");
    }
    for(const auto& [name,a]:m.arrays) {
        require(!name.empty() && name.find_first_not_of("abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789_")==std::string::npos,"invalid array name");
        require(a.components>0 && a.size()%a.components==0 && !a.unit.empty(),"invalid attribute shape/unit");
        const auto rows=a.size()/a.components;
        if(a.location=="node") require(rows==m.node_ids.size(),"node attribute extent");
        else if(a.location=="face") require(rows==m.face_ids.size(),"face attribute extent");
        else if(a.location=="cell") require(rows==m.cell_ids.size(),"cell attribute extent");
        else if(a.location=="incidence") require(rows==m.signs.size(),"incidence attribute extent");
        else if(a.location=="nnc") require(rows==nnc_count,"NNC attribute extent");
        else require(a.location=="metadata","unknown attribute location");
        validate_reference_geometry_array_contract(name,a,m.dimension);
        if(const auto* v=std::get_if<std::vector<double>>(&a.values))
            for(auto x:*v) require(std::isfinite(x),"nonfinite attribute");
    }
    if(auto it=m.arrays.find("mrst_index_map");it!=m.arrays.end()) {
        require(it->second.location=="cell" && it->second.components==1 &&
            std::holds_alternative<std::vector<MeshIndex>>(it->second.values),"invalid MRST indexMap");
        const auto& v=std::get<std::vector<MeshIndex>>(it->second.values); ids(v);
        for(auto x:v) require(x>0,"MRST external indexMap is 1-based");
    }
    if(auto it=m.arrays.find("cart_dims");it!=m.arrays.end()) {
        const auto& a=it->second;
        require(a.location=="metadata" && a.components==1 &&
            a.size()==static_cast<std::size_t>(m.dimension) &&
            std::holds_alternative<std::vector<MeshIndex>>(a.values),"invalid logical dimensions");
        MeshIndex count=1;
        for(auto n:std::get<std::vector<MeshIndex>>(a.values)) {
            require(n>0 && n<=std::numeric_limits<MeshIndex>::max()/count,"logical dimension overflow");count*=n;
        }
        if(auto map=m.arrays.find("mrst_index_map");map!=m.arrays.end())
            for(auto n:std::get<std::vector<MeshIndex>>(map->second.values)) require(n<=count,"indexMap exceeds logical grid");
    }
    return neighbors;
}

/// Materialize existing MPMC Topology only when a consumer needs it. The
/// interchange path avoids a second topology-sized copy. Checked narrowing
/// protects the existing 32-bit local-index and CSR contracts.
[[nodiscard]] inline Topology face_mesh_topology(const FaceMesh& m) {
    const auto neighbors=validate_face_mesh(m);
    const auto narrow=[](MeshIndex x) {
        if(x>std::numeric_limits<std::uint32_t>::max()) throw std::length_error("face mesh exceeds local topology capacity");
        return static_cast<std::uint32_t>(x);
    };
    Topology::EntityIds ids;
    for(auto x:m.node_ids) ids.vertices.emplace_back(x);
    for(auto x:m.face_ids) ids.faces.emplace_back(x);
    for(auto x:m.cell_ids) ids.cells.emplace_back(x);
    std::vector<CsrAdjacency> relations;
    const auto append=[&](EntityKind from,EntityKind to,std::size_t count,
                         const auto& offsets,const auto& values) {
        std::vector<CsrAdjacency::Offset> o; std::vector<LocalIndex> v;
        o.reserve(offsets.size()); v.reserve(values.size());
        for(auto x:offsets) o.push_back(narrow(x));
        for(auto x:values) v.emplace_back(narrow(x));
        relations.emplace_back(from,to,count,std::move(o),std::move(v));
    };
    append(EntityKind::face,EntityKind::vertex,m.points.size(),m.face_offsets,m.face_nodes);
    append(EntityKind::cell,EntityKind::face,m.face_ids.size(),m.cell_offsets,m.cell_faces);
    std::vector<MeshIndex> offsets{0},values;
    for(auto n:neighbors) {
        for(auto c:n) if(c>=0) values.push_back(static_cast<MeshIndex>(c));
        offsets.push_back(values.size());
    }
    append(EntityKind::face,EntityKind::cell,m.cell_ids.size(),offsets,values);
    offsets={0}; values.clear();
    for(std::size_t c=0;c<m.cell_ids.size();++c) {
        std::vector<MeshIndex> nodes;
        for(auto f:m.cell(c)) { auto r=m.face(static_cast<std::size_t>(f)); nodes.insert(nodes.end(),r.begin(),r.end()); }
        std::sort(nodes.begin(),nodes.end()); nodes.erase(std::unique(nodes.begin(),nodes.end()),nodes.end());
        values.insert(values.end(),nodes.begin(),nodes.end()); offsets.push_back(values.size());
    }
    append(EntityKind::cell,EntityKind::vertex,m.points.size(),offsets,values);
    return Topology{std::move(ids),std::move(relations)};
}

/// Scalar cell-centered matrix graph: diagonal, geometric neighbors and NNCs.
/// No transmissibility/TPFA admissibility is implied by this structural graph.
[[nodiscard]] inline CsrAdjacency face_mesh_cell_graph(const FaceMesh& m) {
    const auto neighbors=validate_face_mesh(m);
    face_mesh_detail::require(m.cell_ids.size()<=std::numeric_limits<std::uint32_t>::max(),"matrix graph capacity");
    std::vector<std::pair<std::uint32_t,std::uint32_t>> edges;
    edges.reserve(m.cell_ids.size()+2*m.face_ids.size());
    for(std::size_t c=0;c<m.cell_ids.size();++c) edges.emplace_back(static_cast<std::uint32_t>(c),static_cast<std::uint32_t>(c));
    const auto connect=[&](MeshIndex a,MeshIndex b) {
        edges.emplace_back(static_cast<std::uint32_t>(a),static_cast<std::uint32_t>(b));
        edges.emplace_back(static_cast<std::uint32_t>(b),static_cast<std::uint32_t>(a));
    };
    for(auto n:neighbors) if(n[0]>=0 && n[1]>=0) connect(static_cast<MeshIndex>(n[0]),static_cast<MeshIndex>(n[1]));
    if(auto it=m.arrays.find("nnc_cells");it!=m.arrays.end()) {
        const auto& v=std::get<std::vector<MeshIndex>>(it->second.values);
        for(std::size_t i=0;i<v.size();i+=2) connect(v[i],v[i+1]);
    }
    std::sort(edges.begin(),edges.end()); edges.erase(std::unique(edges.begin(),edges.end()),edges.end());
    face_mesh_detail::require(edges.size()<=std::numeric_limits<std::uint32_t>::max(),"matrix graph CSR capacity");
    std::vector<CsrAdjacency::Offset> offsets(m.cell_ids.size()+1,0);
    std::vector<LocalIndex> values; values.reserve(edges.size());
    for(auto [a,b]:edges) { ++offsets[static_cast<std::size_t>(a)+1]; values.emplace_back(b); }
    std::partial_sum(offsets.begin(),offsets.end(),offsets.begin());
    return {EntityKind::cell,EntityKind::cell,m.cell_ids.size(),std::move(offsets),std::move(values)};
}
} // namespace mpmc::mesh
#endif
