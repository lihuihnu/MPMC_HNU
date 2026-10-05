#include <mpmc/mesh/computational_mesh.hpp>
#include <mpmc/mesh/dof_numbering.hpp>
#include <mpmc/mesh/mesh_exchange_io.hpp>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <iostream>
#include <optional>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

void computational_mesh_header_probe();
bool mesh_exchange_document_header_self_contained();
bool linear_cell_geometric_operator_3d_header_self_contained();
namespace mesh = mpmc::mesh;
using K = mesh::EntityKind;
using I = mesh::LocalIndex;
using G = mesh::GlobalEntityId;

namespace {
void require(bool condition, const char* message) {
    if (!condition) { throw std::runtime_error(message); }
}
void close(double actual, double expected, const char* message) {
    require(std::isfinite(actual) && std::abs(actual - expected) < 1e-12, message);
}
I local(std::size_t value) { return I{static_cast<I::value_type>(value)}; }
template <typename T> std::vector<T> copy(std::span<const T> data) {
    return {data.begin(), data.end()};
}
mesh::FaceBoundarySnapshot tagged_boundary(const mesh::Topology& topology) {
    std::vector<mesh::PhysicalTag> tags(topology.entity_count(K::face),mesh::PhysicalTag{0U});
    for (std::size_t f = 0; f < tags.size(); ++f) {
        if (topology.relation(K::face,K::cell).adjacent(local(f)).size() == 1U) {
            tags[f] = mesh::PhysicalTag{72U};
            break;
        }
    }
    return mesh::make_face_boundary_snapshot(topology,tags);
}

mesh::MeshExchangeDocument document(
    int dimension, mesh::Topology topology, std::vector<mesh::Coordinate3D> points,
    std::optional<mesh::FaceBoundarySnapshot> boundary = std::nullopt) {
    std::vector<double> values;
    for (const auto id : topology.global_ids(K::face)) {
        values.push_back(static_cast<double>(id.value()));
    }
    std::vector<mesh::DenseFieldSnapshot> fields;
    fields.push_back(mesh::DenseFieldSnapshot::create(
        topology, K::face, 1U, std::move(values),
        {"face_marker", "1", {mesh::FieldSourceKind::synthetic_test,
                              "analytic mesh fixture", "1", "face stable ID"}}));
    std::vector<mesh::MeshExchangeGroup> groups{
        {K::cell, 71U, "material", {topology.global_id(K::cell, I{0})}},
        {K::face, 72U, "face_set", {topology.global_id(K::face, I{0})}}};
    return mesh::MeshExchangeDocument::create(mesh::MeshExchangeFormat::vtu_ascii,
        dimension, std::move(topology), std::move(points), std::move(boundary),
        std::move(fields), std::move(groups), std::nullopt);
}

mesh::MeshExchangeDocument fixture_2d() {
    using T = mesh::LinearCellType2D;
    const std::vector<mesh::LinearCell2D> cells{
        {G{101}, T::quadrilateral, {I{0},I{1},I{2},I{3}}},
        {G{307}, T::triangle, {I{1},I{4},I{2}}}};
    const auto grid = mesh::make_linear_mesh_2d(
        {G{91},G{3},G{507},G{22},G{66}},
        {{0,0},{1,0},{1,1},{0,1},{2,0}}, cells);
    std::vector<mesh::Coordinate3D> points;
    for (const auto p : grid.geometry.vertex_coordinates_m()) { points.push_back({p.x_m,p.y_m,0}); }
    return document(2, grid.topology, std::move(points), tagged_boundary(grid.topology));
}

mesh::MeshExchangeDocument fixture_3d() {
    using T = mesh::LinearCellType3D;
    const std::vector<mesh::Coordinate3D> points{
        {0,0,0},{1,0,0},{1,1,0},{0,1,0},
        {0,0,1},{1,0,1},{1,1,1},{0,1,1},{.5,.5,2},
        {3,0,0},{4,0,0},{3,1,0},{3,0,1},{4,0,1},{3,1,1},{3,0,2}};
    std::vector<G> ids;
    for (std::size_t v = 0; v < points.size(); ++v) { ids.push_back(G{51U + 17U * v}); }
    const std::vector<mesh::LinearCell3D> cells{
        {G{101},T::hexahedron,{I{0},I{1},I{2},I{3},I{4},I{5},I{6},I{7}}},
        {G{307},T::pyramid,{I{4},I{5},I{6},I{7},I{8}}},
        {G{415},T::wedge,{I{9},I{10},I{11},I{12},I{13},I{14}}},
        {G{599},T::tetrahedron,{I{12},I{13},I{14},I{15}}}};
    const auto grid = mesh::make_linear_mesh_3d(ids, points, cells);
    return document(3, grid.topology, points, tagged_boundary(grid.topology));
}

enum class Mutation { reorder_faces, wrong_support, duplicate_support,
                      wrong_cell_faces, missing_relation, duplicate_annotation,
                      crossed_face, extra_relation, explicit_edges };

mesh::Topology changed_topology(const mesh::Topology& source, Mutation mutation) {
    const auto faces = source.entity_count(K::face);
    mesh::Topology::EntityIds ids{copy(source.global_ids(K::vertex)), {},
                                 copy(source.global_ids(K::face)), copy(source.global_ids(K::cell))};
    if (mutation == Mutation::explicit_edges) { ids.edges.push_back(G{909}); }
    if (mutation == Mutation::reorder_faces) { std::reverse(ids.faces.begin(), ids.faces.end()); }
    const std::array<std::pair<K,K>,4> kinds{{
        {K::cell,K::vertex}, {K::cell,K::face}, {K::face,K::vertex}, {K::face,K::cell}}};
    std::vector<mesh::CsrAdjacency> relations;
    for (const auto& [from,to] : kinds) {
        if (mutation == Mutation::missing_relation && from == K::cell && to == K::face) { continue; }
        std::vector<mesh::CsrAdjacency::Offset> offsets{0U};
        std::vector<I> entries;
        bool modified = false;
        for (std::size_t index = 0; index < source.entity_count(from); ++index) {
            const auto source_index = mutation == Mutation::reorder_faces && from == K::face ?
                                      faces - 1U - index : index;
            auto row = copy(source.relation(from,to).adjacent(local(source_index)));
            if (mutation == Mutation::reorder_faces) {
                if (to == K::face) {
                    for (auto& entity : row) { entity = local(faces - 1U - entity.value()); }
                }
                if (from == K::face && to == K::cell) { std::reverse(row.begin(),row.end()); }
            } else if (from == K::face && to == K::cell && !modified) {
                if (mutation == Mutation::wrong_support && row.size() == 1U) {
                    row[0] = row[0] == I{0} ? I{1} : I{0};
                    modified = true;
                } else if (mutation == Mutation::duplicate_support && row.size() == 2U) {
                    row[1] = row[0];
                    modified = true;
                }
            } else if (mutation == Mutation::wrong_cell_faces && from == K::cell && to == K::face && index == 0U) {
                row[1] = row[0];
            } else if (mutation == Mutation::duplicate_annotation && from == K::face && to == K::vertex && index == 1U) {
                row = copy(source.relation(from,to).adjacent(I{0}));
            } else if (mutation == Mutation::crossed_face && from == K::face && to == K::vertex &&
                       !modified && row.size() == 4U) {
                std::swap(row[1],row[2]);
                modified = true;
            }
            entries.insert(entries.end(),row.begin(),row.end());
            offsets.push_back(static_cast<mesh::CsrAdjacency::Offset>(entries.size()));
        }
        relations.emplace_back(from,to,source.entity_count(to),std::move(offsets),std::move(entries));
    }
    if (mutation == Mutation::extra_relation) {
        relations.emplace_back(K::vertex,K::cell,source.entity_count(K::cell),
            std::vector<mesh::CsrAdjacency::Offset>(source.entity_count(K::vertex)+1U,0U),std::vector<I>{});
    }
    return {std::move(ids),std::move(relations)};
}

void verify_binding(const mesh::MeshExchangeDocument& doc, const mesh::Topology& prepared) {
    const auto& source = doc.topology();
    for (const auto kind : {K::vertex,K::face,K::cell}) {
        require(std::equal(source.global_ids(kind).begin(),source.global_ids(kind).end(),
                           prepared.global_ids(kind).begin(),prepared.global_ids(kind).end()),
                "preparation must preserve all entity local identities");
    }
    const auto& marker = doc.fields().at(K::face, "face_marker");
    const std::vector<mesh::DofVariable> variables{{"cell_unknown",K::cell,2U},{"face_unknown",K::face,1U}};
    const auto before = mesh::DofLayout::create(source, variables);
    const auto after = mesh::DofLayout::create(prepared, variables);
    const auto numbering_before = mesh::DofNumberingSnapshot::create_serial(
        before, mesh::make_serial_partition_snapshot(source));
    const auto numbering_after = mesh::DofNumberingSnapshot::create_serial(
        after, mesh::make_serial_partition_snapshot(prepared));
    for (std::size_t f = 0; f < prepared.entity_count(K::face); ++f) {
        close(marker.value(local(f),0U), static_cast<double>(prepared.global_id(K::face,local(f)).value()),
              "source face field remains attached to stable face ID");
        const auto offset = before.scalar_offset("face_unknown",local(f),0U);
        require(offset == after.scalar_offset("face_unknown",local(f),0U) &&
                numbering_before.global_index(offset) == numbering_after.global_index(offset),
                "face DoF local/global index binding unchanged");
    }
    for (const auto& group : doc.groups()) {
        for (const auto id : group.members) {
            const auto ids = prepared.global_ids(group.location);
            require(std::find(ids.begin(),ids.end(),id) != ids.end(), "group entity ID remains present");
        }
    }
}

void verify_prepared(const mesh::MeshExchangeDocument& doc, bool metadata) {
    const auto verify_tags = [&](const mesh::FaceBoundarySnapshot& ready) {
        if (doc.face_boundary()) {
            for (std::size_t f = 0; f < ready.face_count(); ++f) {
                require(ready.physical_tag(local(f)) == doc.face_boundary()->physical_tag(local(f)),
                        "face physical tag remains aligned to source ID");
            }
        }
    };
    if (doc.dimension() == 2) {
        const auto ready = mesh::prepare_linear_mesh_2d(doc);
        verify_tags(ready.face_boundary);
        if (metadata) { verify_binding(doc,ready.topology); }
        close(ready.geometry.cell_area_m2(I{0}),1.0,"quad analytic area");
        close(ready.geometry.cell_area_m2(I{1}),.5,"triangle analytic area");
        require(ready.cell_types[0] == mesh::LinearCellType2D::quadrilateral &&
                ready.cell_types[1] == mesh::LinearCellType2D::triangle,"original cell types retained");
        for (std::size_t f = 0; f < ready.geometry.face_count(); ++f) {
            const auto face = local(f);
            const auto owner = ready.topology.relation(K::face,K::cell).adjacent(face).front();
            require(ready.geometry.face_owner(face) == owner, "source 2D owner retained after remap");
            const auto c = ready.geometry.cell_centroid_m(owner);
            const auto p = ready.geometry.face_centroid_m(face);
            const auto vertices = ready.topology.relation(K::face,K::vertex).adjacent(face);
            const auto a = doc.vertex_coordinates_m()[vertices[0].value()];
            const auto b = doc.vertex_coordinates_m()[vertices[1].value()];
            close(p.x_m,.5*(a.x_m+b.x_m),"remapped 2D face-centroid x");
            close(p.y_m,.5*(a.y_m+b.y_m),"remapped 2D face-centroid y");
            const auto n = ready.geometry.face_owner_unit_normal(face);
            require((p.x_m-c.x_m)*n.x + (p.y_m-c.y_m)*n.y > 0.0,"remapped 2D normal outward");
        }
    } else {
        const auto ready = mesh::prepare_linear_mesh_3d(doc);
        verify_tags(ready.face_boundary);
        if (metadata) { verify_binding(doc,ready.topology); }
        const auto op = mesh::make_cell_face_geometric_operator_3d(ready);
        const std::array<double,4> volumes{1.0,1.0/3.0,.5,1.0/6.0};
        for (std::size_t c = 0; c < volumes.size(); ++c) {
            close(ready.cell_volumes_m3[c],volumes[c],"prepared analytic volume");
        }
        for (std::size_t f = 0; f < op.face_count(); ++f) {
            const auto face = local(f);
            require(op.face_owner(face) == doc.topology().relation(K::face,K::cell).adjacent(face).front(),
                    "source 3D owner retained after remap");
            const auto vertices = ready.topology.relation(K::face,K::vertex).adjacent(face);
            mesh::Coordinate3D mean{0,0,0};
            for (const auto vertex : vertices) {
                const auto p = doc.vertex_coordinates_m()[vertex.value()];
                mean.x_m += p.x_m; mean.y_m += p.y_m; mean.z_m += p.z_m;
            }
            const double count = static_cast<double>(vertices.size());
            const auto p = ready.face_geometry.face_centroid_m(face);
            // All fixture quads are rectangles; triangles have exact vertex-mean centroids.
            close(p.x_m,mean.x_m/count,"remapped 3D face-centroid x");
            close(p.y_m,mean.y_m/count,"remapped 3D face-centroid y");
            close(p.z_m,mean.z_m/count,"remapped 3D face-centroid z");
        }
    }
}

template <typename F> void rejects(F&& run, const char* message) {
    bool rejected = false;
    try { run(); } catch (const std::invalid_argument&) { rejected = true; }
    require(rejected,message);
}

void contract_and_failures(const mesh::MeshExchangeDocument& base) {
    verify_prepared(base,true);
    auto reordered_class = copy(base.face_boundary()->classifications());
    auto reordered_tags = copy(base.face_boundary()->physical_tags());
    std::reverse(reordered_class.begin(),reordered_class.end());
    std::reverse(reordered_tags.begin(),reordered_tags.end());
    const auto reordered = document(base.dimension(),changed_topology(base.topology(),Mutation::reorder_faces),
                                    copy(base.vertex_coordinates_m()),
                                    mesh::FaceBoundarySnapshot{std::move(reordered_class),std::move(reordered_tags)});
    verify_prepared(reordered,true);
    const auto prepare = [&](const mesh::MeshExchangeDocument& doc) {
        if (doc.dimension() == 2) { (void)mesh::prepare_linear_mesh_2d(doc); }
        else { (void)mesh::prepare_linear_mesh_3d(doc); }
    };
    for (const auto mutation : {Mutation::wrong_support,Mutation::duplicate_support,
                               Mutation::wrong_cell_faces,Mutation::missing_relation,Mutation::duplicate_annotation,
                               Mutation::extra_relation,Mutation::explicit_edges}) {
        const auto bad = document(base.dimension(),changed_topology(base.topology(),mutation),copy(base.vertex_coordinates_m()));
        rejects([&] { prepare(bad); },"malformed topology must be rejected, not repaired");
    }
    if (base.dimension() == 3) {
        const auto crossed = document(3,changed_topology(base.topology(),Mutation::crossed_face),
                                      copy(base.vertex_coordinates_m()));
        rejects([&] { prepare(crossed); },"noncyclic quad face must not be silently repaired");
    }
    auto classifications = copy(base.face_boundary()->classifications());
    const auto found = std::find(classifications.begin(),classifications.end(),mesh::FaceClassification::boundary);
    require(found != classifications.end(),"fixture has boundary faces");
    *found = mesh::FaceClassification::interior;
    const auto wrong_class = document(base.dimension(),base.topology(),copy(base.vertex_coordinates_m()),
        mesh::FaceBoundarySnapshot{std::move(classifications),
            std::vector<mesh::PhysicalTag>(base.face_boundary()->face_count(),mesh::PhysicalTag{0U})});
    rejects([&] { prepare(wrong_class); },"wrong source boundary classification must be rejected");
    auto degenerate = copy(base.vertex_coordinates_m());
    degenerate[1] = degenerate[0];
    const auto bad_geometry = document(base.dimension(),base.topology(),std::move(degenerate));
    rejects([&] { prepare(bad_geometry); },"degenerate geometry rejected");
    rejects([&] {
        if (base.dimension() == 2) { (void)mesh::prepare_linear_mesh_3d(base); }
        else { (void)mesh::prepare_linear_mesh_2d(base); }
    },"dimension-specific API rejects mismatched input");
}

void supported_format_paths(const mesh::MeshExchangeDocument& source) {
    // This is an integration test of each native importer and preparation.
    // It is not claimed as an independent external-software reader oracle.
    const auto gmsh = mesh::export_gmsh_4_1_ascii(source);
    const auto vtu = mesh::export_vtu_ascii(source);
    require(gmsh.exported() && vtu.exported(),"fixture export available");
    if (source.dimension() == 2) {
        verify_prepared(mesh::make_mesh_exchange_document(mesh::import_gmsh_4_1_ascii(*gmsh.content, 1.0)),false);
        verify_prepared(mesh::make_mesh_exchange_document(mesh::import_vtu_ascii(*vtu.content)),false);
    } else {
        verify_prepared(mesh::make_mesh_exchange_document(mesh::import_gmsh_4_1_ascii_3d(*gmsh.content, 1.0)),false);
        verify_prepared(mesh::make_mesh_exchange_document(mesh::import_vtu_ascii_3d(*vtu.content)),false);
    }
}

void grdecl_path() {
    // Explicit unit cube, plus material data, with no exporter used to create it.
    constexpr auto input = R"(SPECGRID
1 1 1 1 F /
COORD
0 0 0 0 0 1 1 0 0 1 0 1 0 1 0 0 1 1 1 1 0 1 1 1 /
ZCORN
0 0 0 0 1 1 1 1 /
ACTNUM
1 /
PORO
0.25 /
PERMX
10 /
PERMY
20 /
PERMZ
30 /
)";
    const auto doc = mesh::make_mesh_exchange_document(mesh::import_grdecl(input,{1.0,1e-15}));
    const auto grid = mesh::prepare_linear_mesh_3d(doc);
    close(grid.cell_volumes_m3[0],1.0,"GRDECL cube ready volume");
    close(doc.fields().at(K::cell,"PORO").value(I{0},0),.25,"GRDECL material alignment");
    const auto geometry = mesh::make_cell_face_geometric_operator_3d(grid);
    for (std::size_t f = 0; f < geometry.face_count(); ++f) {
        close(geometry.owner_normal_distance_m(local(f)),.5,"GRDECL analytic half distance");
    }
}
} // namespace

int main() {
    try {
        computational_mesh_header_probe();
        require(linear_cell_geometric_operator_3d_header_self_contained(),"linear geometry public header probe");
        require(mesh_exchange_document_header_self_contained(),"exchange document public header probe");
        const auto two = fixture_2d();
        const auto three = fixture_3d();
        contract_and_failures(two);
        contract_and_failures(three);
        supported_format_paths(two);
        supported_format_paths(three);
        grdecl_path();
        std::cout << "[PASS] mesh.core.computational_mesh\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "[FAIL] mesh.core.computational_mesh: " << error.what() << '\n';
        return 1;
    }
}
