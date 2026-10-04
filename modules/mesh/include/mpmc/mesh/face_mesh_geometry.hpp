#ifndef MPMC_MESH_FACE_MESH_GEOMETRY_HPP
#define MPMC_MESH_FACE_MESH_GEOMETRY_HPP
#include <mpmc/mesh/face_mesh.hpp>

namespace mpmc::mesh {
struct FaceMeshGeometry {
    std::vector<double> face_areas, cell_volumes;
    std::vector<MeshPoint> face_centroids, area_vectors, cell_centroids;
    double maximum_relative_closure{};
    double maximum_reference_relative_error{};
    std::size_t compared_reference_values{};
};

/// Piecewise-planar surface: each 3D ring is triangulated about its arithmetic
/// node mean, once per shared face. Scalar area is sum(|triangle area vector|),
/// while flux uses the SUM of the oriented area vectors (not area*unit normal).
/// Volume and first moments use oriented tetrahedra and the divergence theorem.
/// Concave cells are allowed; negative sub-tetrahedra are NOT clipped/abs'd.
/// 2D uses signed boundary integration. SI m, m^2, m^3; not a TPFA certificate.
[[nodiscard]] inline FaceMeshGeometry prepare_face_mesh_geometry(
    const FaceMesh& m, bool compare_reference=true, double reference_tolerance=1e-8) {
    using namespace face_mesh_detail;
    (void)validate_face_mesh(m,true);
    require(std::isfinite(reference_tolerance) && reference_tolerance>0,"invalid geometry tolerance");
    FaceMeshGeometry g;
    const auto nf=m.face_ids.size(), nc=m.cell_ids.size();
    g.face_areas.resize(nf); g.face_centroids.resize(nf); g.area_vectors.resize(nf);
    g.cell_volumes.resize(nc); g.cell_centroids.resize(nc);
    std::vector<MeshPoint> centers(nf);
    std::vector<double> face_roundoff(nf),cell_roundoff(nc);
    constexpr double eps=4096*std::numeric_limits<double>::epsilon();
    for(std::size_t f=0;f<nf;++f) {
        const auto r=m.face(f);
        const auto base=m.points[static_cast<std::size_t>(r[0])];
        MeshPoint center{};
        for(auto n:r) center=add(center,sub(m.points[static_cast<std::size_t>(n)],base));
        center=add(base,mul(center,1/static_cast<double>(r.size()))); centers[f]=center;
        if(m.dimension==2) {
            const auto d=sub(m.points[static_cast<std::size_t>(r[1])],base);
            g.face_areas[f]=norm(d); g.area_vectors[f]={d[1],-d[0],0}; g.face_centroids[f]=center;
        } else {
            MeshPoint moment{}, av{}; double area=0;
            for(std::size_t j=0;j<r.size();++j) {
                const auto a=sub(m.points[static_cast<std::size_t>(r[j])],center);
                const auto b=sub(m.points[static_cast<std::size_t>(r[(j+1)%r.size()])],center);
                const auto n=mul(cross(a,b),0.5); const auto weight=norm(n);
                av=add(av,n); area+=weight;
                moment=add(moment,mul(add(a,b),weight/3));
            }
            require(std::isfinite(area) && area>0 && norm(av)>eps*area,"degenerate or folded face "+std::to_string(f));
            // Oppositely oriented fan triangles require a different surface
            // definition: retain transport but refuse silent geometry repair.
            for(std::size_t j=0;j<r.size();++j) {
                const auto a=sub(m.points[static_cast<std::size_t>(r[j])],center);
                const auto b=sub(m.points[static_cast<std::size_t>(r[(j+1)%r.size()])],center);
                const auto n=cross(a,b);
                require(dot(n,av)>=-eps*norm(n)*norm(av),"face fan folds at face "+std::to_string(f));
            }
            g.face_areas[f]=area; g.area_vectors[f]=av; g.face_centroids[f]=add(center,mul(moment,1/area));
        }
        double coordinate_magnitude=1, perimeter=0;
        for(std::size_t j=0;j<r.size();++j) {
            const auto p=m.points[static_cast<std::size_t>(r[j])];
            for(auto x:p) coordinate_magnitude=std::max(coordinate_magnitude,std::abs(x));
            perimeter+=norm(sub(p,m.points[static_cast<std::size_t>(r[(j+1)%r.size()])]));
        }
        face_roundoff[f]=64*std::numeric_limits<double>::epsilon()*coordinate_magnitude*(m.dimension==3?perimeter:1);
        require(std::isfinite(g.face_areas[f]) && g.face_areas[f]>0,"zero/nonfinite face measure");
    }
    for(std::size_t c=0;c<nc;++c) {
        const auto faces=m.cell(c);
        const auto base=g.face_centroids[static_cast<std::size_t>(faces[0])];
        MeshPoint origin{};
        for(auto f:faces) origin=add(origin,sub(g.face_centroids[static_cast<std::size_t>(f)],base));
        origin=add(base,mul(origin,1/static_cast<double>(faces.size())));
        MeshPoint closure{}, moment{}; double volume=0, surface=0, scale=0, edge_length=0, coordinate_scale=1;
        for(std::size_t j=0;j<faces.size();++j) {
            const auto f=static_cast<std::size_t>(faces[j]);
            const double sign=static_cast<double>(m.signs[static_cast<std::size_t>(m.cell_offsets[c])+j]);
            closure=add(closure,mul(g.area_vectors[f],sign)); surface+=g.face_areas[f];
            const auto ring=m.face(f);
            for(std::size_t k=0;k<ring.size();++k) {
                const auto p=m.points[static_cast<std::size_t>(ring[k])];
                scale=std::max(scale,norm(sub(p,origin)));
                for(auto x:p) coordinate_scale=std::max(coordinate_scale,std::abs(x));
                edge_length+=norm(sub(p,m.points[static_cast<std::size_t>(ring[(k+1)%ring.size()])]));
            }
            if(m.dimension==2) {
                const auto a=sub(m.points[static_cast<std::size_t>(ring[0])],origin);
                const auto b=sub(m.points[static_cast<std::size_t>(ring[1])],origin);
                const double v=sign*(a[0]*b[1]-a[1]*b[0])/2;
                volume+=v; moment=add(moment,mul(add(a,b),v/3));
            } else {
                const auto h=sub(centers[f],origin);
                for(std::size_t k=0;k<ring.size();++k) {
                    const auto a=sub(m.points[static_cast<std::size_t>(ring[k])],origin);
                    const auto b=sub(m.points[static_cast<std::size_t>(ring[(k+1)%ring.size()])],origin);
                    const double v=sign*dot(h,cross(a,b))/6;
                    volume+=v; moment=add(moment,mul(add(add(h,a),b),v/4));
                }
            }
        }
        const auto relative=norm(closure)/surface;
        g.maximum_relative_closure=std::max(g.maximum_relative_closure,relative);
        // Absolute-coordinate quantization perturbs cross products by O(u*X*L).
        // Include that forward-error bound for large map coordinates, in
        // addition to relative accumulation error; never scale by source residual.
        const double closure_bound=eps+64*std::numeric_limits<double>::epsilon()*coordinate_scale*
            (m.dimension==3?edge_length:static_cast<double>(faces.size()))/surface;
        require(std::isfinite(relative) && relative<=closure_bound,"cell area-vector closure failed "+std::to_string(c));
        const auto dimensional_scale=m.dimension==2?scale*scale:scale*scale*scale;
        require(std::isfinite(volume) && volume>eps*dimensional_scale,"inverted/degenerate cell "+std::to_string(c));
        cell_roundoff[c]=64*std::numeric_limits<double>::epsilon()*coordinate_scale*surface;
        g.cell_volumes[c]=volume; g.cell_centroids[c]=add(origin,mul(moment,1/volume));
        for(auto x:g.cell_centroids[c]) require(std::isfinite(x),"nonfinite cell centroid");
    }
    if(compare_reference) {
        const auto compare=[&](const char* name,const auto& values,std::size_t components,const char* unit) {
            const auto it=m.arrays.find(name); if(it==m.arrays.end()) return;
            const auto& a=it->second;
            require(a.components==components && a.unit==unit && std::holds_alternative<std::vector<double>>(a.values),"invalid reference geometry metadata");
            const auto& expected=std::get<std::vector<double>>(a.values);
            require(expected.size()==values.size()*components,"reference geometry extent");
            for(std::size_t i=0;i<expected.size();++i) {
                ++g.compared_reference_values;
                double actual;
                if constexpr(std::is_same_v<typename std::decay_t<decltype(values)>::value_type,double>) actual=values[i];
                else actual=values[i/components][i%components];
                const double error=std::abs(actual-expected[i])/std::max({1.0,std::abs(actual),std::abs(expected[i])});
                g.maximum_reference_relative_error=std::max(g.maximum_reference_relative_error,error);
                const std::string key=name;
                const double roundoff=key=="reference_face_areas" || key=="reference_face_normals" ? face_roundoff[i/components] :
                    key=="reference_cell_volumes" ? cell_roundoff[i/components] : 0;
                const double allowed=reference_tolerance*std::max({1.0,std::abs(actual),std::abs(expected[i])})+roundoff;
                require(std::abs(actual-expected[i])<=allowed,"reference geometry mismatch in "+std::string(name)+" at "+std::to_string(i));
            }
        };
        compare("reference_face_areas",g.face_areas,1,m.dimension==2?"m":"m2");
        compare("reference_face_normals",g.area_vectors,3,m.dimension==2?"m":"m2");
        compare("reference_face_centroids",g.face_centroids,3,"m");
        compare("reference_cell_volumes",g.cell_volumes,1,m.dimension==2?"m2":"m3");
        compare("reference_cell_centroids",g.cell_centroids,3,"m");
    }
    return g;
}
} // namespace mpmc::mesh
#endif
