#include <mpmc/mesh/linear_cell_mesh_2d.hpp>
#include <mpmc/mesh/gmsh_4_1.hpp>
#include <mpmc/mesh/vtu.hpp>

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <string>

namespace {
namespace mesh = mpmc::mesh;
using mesh::GlobalEntityId;
using mesh::LocalIndex;
using mesh::LinearCell2D;
using mesh::LinearCellType2D;

void check(bool condition) {
    if (!condition) throw std::runtime_error("linear 2D contract regression");
}
void near(double value, double expected, double tolerance = 1e-12) {
    check(std::isfinite(value) && std::abs(value - expected) <= tolerance);
}
void replace_once(std::string& text, const std::string& from, const std::string& to) {
    const auto at = text.find(from);
    check(at != std::string::npos);
    text.replace(at, from.size(), to);
}
template <typename Exception = std::invalid_argument, typename F>
void rejects(F&& action) {
    bool caught = false;
    try { action(); } catch (const Exception&) { caught = true; }
    check(caught);
}
std::vector<GlobalEntityId> ids() {
    return {GlobalEntityId{1}, GlobalEntityId{2}, GlobalEntityId{3}, GlobalEntityId{4}, GlobalEntityId{5}};
}
std::vector<mesh::Coordinate2D> points() {
    // Rectangle [0,2]x[0,1] plus right triangle on its right edge.
    return {{0,0}, {2,0}, {2,1}, {0,1}, {3,0}};
}
std::vector<LinearCell2D> cells() {
    return {{GlobalEntityId{900}, LinearCellType2D::quadrilateral,
             {LocalIndex{0}, LocalIndex{1}, LocalIndex{2}, LocalIndex{3}}},
            {GlobalEntityId{30}, LinearCellType2D::triangle,
             {LocalIndex{1}, LocalIndex{4}, LocalIndex{2}}}};
}
void check_closure(const mesh::Topology& topology, const mesh::Geometry2D& geometry) {
    const auto& cf = topology.relation(mesh::EntityKind::cell, mesh::EntityKind::face);
    for (std::size_t i = 0; i < geometry.cell_count(); ++i) {
        const LocalIndex cell{static_cast<LocalIndex::value_type>(i)};
        double sx = 0, sy = 0;
        for (const auto face : cf.adjacent(cell)) {
            const auto n = geometry.face_owner_unit_normal(face);
            const double sign = geometry.face_owner(face) == cell ? 1.0 : -1.0;
            sx += sign * geometry.face_length_m(face) * n.x;
            sy += sign * geometry.face_length_m(face) * n.y;
        }
        near(sx, 0); near(sy, 0);
    }
}

// Independent ASCII fixtures, not produced by either writer under test.
std::string gmsh_fixture() {
    return R"($MeshFormat
4.1 0 8
$EndMeshFormat
$Nodes
1 5 1 5
2 1 0 5
1 2 3 4 5
0 0 0
2 0 0
2 1 0
0 1 0
3 0 0
$EndNodes
$Elements
2 2 30 900
2 1 3 1
900 1 2 3 4
2 1 2 1
30 2 5 3
$EndElements
)";
}
std::string vtu_fixture() {
    return R"(<VTKFile type="UnstructuredGrid" version="0.1" byte_order="LittleEndian">
<UnstructuredGrid><Piece NumberOfPoints="5" NumberOfCells="2">
<Points><DataArray type="Float64" NumberOfComponents="3" format="ascii">
0 0 0 2 0 0 2 1 0 0 1 0 3 0 0
</DataArray></Points>
<Cells>
<DataArray type="Int64" Name="connectivity" format="ascii">0 1 2 3 1 4 2</DataArray>
<DataArray type="Int64" Name="offsets" format="ascii">4 7</DataArray>
<DataArray type="UInt8" Name="types" format="ascii">9 5</DataArray>
</Cells>
<PointData><DataArray type="Float64" Name="p" format="ascii">10 20 30 40 50</DataArray></PointData>
<CellData>
<DataArray type="UInt64" Name="mpmc_global_cell_id" format="ascii">900 30</DataArray>
<DataArray type="Float64" Name="marker" format="ascii">90 3</DataArray>
</CellData>
</Piece></UnstructuredGrid></VTKFile>)";
}

void compare_semantics(const mesh::Topology& a, const mesh::Geometry2D& ag,
                       const mesh::Topology& b, const mesh::Geometry2D& bg) {
    for (const auto kind : {mesh::EntityKind::vertex, mesh::EntityKind::cell, mesh::EntityKind::face}) {
        check(a.entity_count(kind) == b.entity_count(kind));
    }
    const auto a_ids = a.global_ids(mesh::EntityKind::cell);
    const auto b_ids = b.global_ids(mesh::EntityKind::cell);
    for (std::size_t i = 0; i < a_ids.size(); ++i) {
        const auto found = std::find(b_ids.begin(), b_ids.end(), a_ids[i]);
        check(found != b_ids.end());
        const LocalIndex ai{static_cast<LocalIndex::value_type>(i)};
        const LocalIndex bi{static_cast<LocalIndex::value_type>(found - b_ids.begin())};
        near(ag.cell_area_m2(ai), bg.cell_area_m2(bi));
        near(ag.cell_centroid_m(ai).x_m, bg.cell_centroid_m(bi).x_m);
        near(ag.cell_centroid_m(ai).y_m, bg.cell_centroid_m(bi).y_m);
        const auto av = a.relation(mesh::EntityKind::cell, mesh::EntityKind::vertex).adjacent(ai);
        const auto bv = b.relation(mesh::EntityKind::cell, mesh::EntityKind::vertex).adjacent(bi);
        check(av.size() == bv.size());
        for (std::size_t j = 0; j < av.size(); ++j) {
            check(a.global_ids(mesh::EntityKind::vertex)[av[j].value()] ==
                  b.global_ids(mesh::EntityKind::vertex)[bv[j].value()]);
        }
    }
    // Face IDs need not survive conversion when the format has no face IDs.
    // Compare by stable endpoint identity and normals relative to the same cell.
    const auto& afv = a.relation(mesh::EntityKind::face, mesh::EntityKind::vertex);
    const auto& bfv = b.relation(mesh::EntityKind::face, mesh::EntityKind::vertex);
    const auto& afc = a.relation(mesh::EntityKind::face, mesh::EntityKind::cell);
    const auto& bfc = b.relation(mesh::EntityKind::face, mesh::EntityKind::cell);
    for (std::size_t i = 0; i < ag.face_count(); ++i) {
        const LocalIndex ai{static_cast<LocalIndex::value_type>(i)};
        const auto av = afv.adjacent(ai);
        bool found = false;
        for (std::size_t j = 0; j < bg.face_count(); ++j) {
            const LocalIndex bi{static_cast<LocalIndex::value_type>(j)};
            const auto bv = bfv.adjacent(bi);
            if (a.global_ids(mesh::EntityKind::vertex)[av[0].value()] != b.global_ids(mesh::EntityKind::vertex)[bv[0].value()] ||
                a.global_ids(mesh::EntityKind::vertex)[av[1].value()] != b.global_ids(mesh::EntityKind::vertex)[bv[1].value()]) continue;
            found = true;
            near(ag.face_length_m(ai), bg.face_length_m(bi));
            near(ag.face_centroid_m(ai).x_m, bg.face_centroid_m(bi).x_m);
            near(ag.face_centroid_m(ai).y_m, bg.face_centroid_m(bi).y_m);
            std::vector<std::uint64_t> ac, bc;
            for (auto c : afc.adjacent(ai)) ac.push_back(a_ids[c.value()].value());
            for (auto c : bfc.adjacent(bi)) bc.push_back(b_ids[c.value()].value());
            std::sort(ac.begin(), ac.end()); std::sort(bc.begin(), bc.end());
            check(ac == bc);
            const double sign = a_ids[ag.face_owner(ai).value()] == b_ids[bg.face_owner(bi).value()] ? 1.0 : -1.0;
            near(ag.face_owner_unit_normal(ai).x, sign * bg.face_owner_unit_normal(bi).x);
            near(ag.face_owner_unit_normal(ai).y, sign * bg.face_owner_unit_normal(bi).y);
        }
        check(found);
    }
}
} // namespace

void linear_cell_mesh_2d_contract() {
    const auto input = cells();
    const std::array<mesh::LinearFaceAnnotation2D, 1> annotations{{
        {{LocalIndex{0}, LocalIndex{1}}, GlobalEntityId{77}, mesh::PhysicalTag{8}}}};
    const auto grid = mesh::make_linear_mesh_2d(ids(), points(), input, annotations, 100);
    check(grid.cell_types == std::vector{LinearCellType2D::quadrilateral, LinearCellType2D::triangle});
    check(grid.topology.entity_count(mesh::EntityKind::face) == 6U);
    check(grid.topology.global_ids(mesh::EntityKind::face)[0].value() == 77U);
    check(grid.face_boundary.physical_tag(LocalIndex{0}).value() == 8U);
    near(grid.geometry.cell_area_m2(LocalIndex{0}), 2.0);
    near(grid.geometry.cell_area_m2(LocalIndex{1}), 0.5);
    near(grid.geometry.cell_centroid_m(LocalIndex{0}).x_m, 1.0);
    near(grid.geometry.cell_centroid_m(LocalIndex{0}).y_m, 0.5);
    near(grid.geometry.cell_centroid_m(LocalIndex{1}).x_m, 7.0 / 3.0);
    near(grid.geometry.cell_centroid_m(LocalIndex{1}).y_m, 1.0 / 3.0);
    check_closure(grid.topology, grid.geometry);
    // Trapezoid: integrating x in [0, 3-y/2], y in [0,2] gives
    // A=5, integral(x dA)=19/3, integral(y dA)=14/3. Not vertex mean.
    const std::array trapezoid_cells{input.front()};
    const auto trapezoid = mesh::make_linear_mesh_2d(
        {GlobalEntityId{0}, GlobalEntityId{40}, GlobalEntityId{5}, GlobalEntityId{80}},
        {{0,0}, {3,0}, {2,2}, {0,2}}, trapezoid_cells);
    near(trapezoid.geometry.cell_area_m2(LocalIndex{0}), 5.0);
    near(trapezoid.geometry.cell_centroid_m(LocalIndex{0}).x_m, 19.0 / 15.0);
    near(trapezoid.geometry.cell_centroid_m(LocalIndex{0}).y_m, 14.0 / 15.0);
    check_closure(trapezoid.topology, trapezoid.geometry);
    for (const double scale : {1e-6, 1.0, 1e6}) {
        auto transformed = points();
        // Rotation by 90 degrees, scale, and translation in the same units.
        for (auto& p : transformed) p = {scale * (10.0 - p.y_m), scale * (20.0 + p.x_m)};
        const auto changed = mesh::make_linear_mesh_2d(ids(), transformed, input);
        near(changed.geometry.cell_area_m2(LocalIndex{0}) / (scale * scale), 2.0);
        near(changed.geometry.cell_centroid_m(LocalIndex{0}).x_m / scale, 9.5);
        near(changed.geometry.cell_centroid_m(LocalIndex{0}).y_m / scale, 21.0);
    }
    auto translated = points();
    for (auto& p : translated) { p.x_m += 1e8; p.y_m -= 1e8; }
    const auto shifted = mesh::make_linear_mesh_2d(ids(), translated, input);
    near(shifted.geometry.cell_area_m2(LocalIndex{0}), 2.0);
    near(shifted.geometry.cell_centroid_m(LocalIndex{0}).x_m, 1e8 + 1.0);
    auto reversed = input;
    for (auto& cell : reversed) std::reverse(cell.vertices.begin(), cell.vertices.end());
    const auto clockwise = mesh::make_linear_mesh_2d(ids(), points(), reversed);
    near(clockwise.geometry.cell_area_m2(LocalIndex{1}), 0.5);
    check_closure(clockwise.topology, clockwise.geometry);
    auto vertex_ids = ids();
    auto reordered_points = points();
    std::swap(vertex_ids[0], vertex_ids[4]); std::swap(reordered_points[0], reordered_points[4]);
    auto reordered_cells = input;
    for (auto& cell : reordered_cells) for (auto& v : cell.vertices) {
        if (v.value() == 0) v = LocalIndex{4}; else if (v.value() == 4) v = LocalIndex{0};
    }
    const auto reordered = mesh::make_linear_mesh_2d(vertex_ids, reordered_points, reordered_cells);
    compare_semantics(grid.topology, grid.geometry, reordered.topology, reordered.geometry);
}

void linear_cell_mesh_2d_invalid() {
    rejects([&] { auto c = cells(); c[0].vertices[1] = c[0].vertices[0]; (void)mesh::make_linear_mesh_2d(ids(), points(), c); });
    rejects<std::out_of_range>([&] { auto c = cells(); c[0].vertices[1] = LocalIndex{99}; (void)mesh::make_linear_mesh_2d(ids(), points(), c); });
    rejects([&] { auto c = cells(); c[0].type = LinearCellType2D::triangle; (void)mesh::make_linear_mesh_2d(ids(), points(), c); });
    rejects([&] { auto c = cells(); c[0].type = static_cast<LinearCellType2D>(99); (void)mesh::make_linear_mesh_2d(ids(), points(), c); });
    rejects([&] { auto c = cells(); c[1].global_id = c[0].global_id; (void)mesh::make_linear_mesh_2d(ids(), points(), c); });
    rejects([&] { auto v = ids(); v[1] = v[0]; (void)mesh::make_linear_mesh_2d(v, points(), cells()); });
    rejects([&] { auto p = points(); p[0].x_m = std::numeric_limits<double>::infinity(); (void)mesh::make_linear_mesh_2d(ids(), p, cells()); });
    rejects([&] { auto p = points(); p[2] = {0.5, 0.25}; (void)mesh::make_linear_mesh_2d(ids(), p, cells()); });
    rejects([&] { auto c = cells(); std::swap(c[0].vertices[1], c[0].vertices[2]); (void)mesh::make_linear_mesh_2d(ids(), points(), c); });
    rejects([&] { auto p = points(); p[2] = {2,0}; (void)mesh::make_linear_mesh_2d(ids(), p, cells()); });
    rejects([&] { auto c = cells(); c[1].vertices = {LocalIndex{1}, LocalIndex{0}, LocalIndex{2}}; (void)mesh::make_linear_mesh_2d(ids(), points(), c); });
    rejects([&] { auto c = cells(); c.push_back(c[1]); c.back().global_id = GlobalEntityId{42}; (void)mesh::make_linear_mesh_2d(ids(), points(), c); });
    rejects<std::length_error>([&] { (void)mesh::make_linear_mesh_2d(ids(), points(), cells(), {}, std::numeric_limits<std::uint64_t>::max()); });
    const mesh::LinearFaceAnnotation2D boundary{{LocalIndex{0}, LocalIndex{1}}, GlobalEntityId{77}, mesh::PhysicalTag{8}};
    rejects([&] { std::array a{boundary, boundary}; (void)mesh::make_linear_mesh_2d(ids(), points(), cells(), a); });
    rejects([&] { auto b = boundary; b.vertices = {LocalIndex{2}, LocalIndex{3}}; std::array a{boundary, b}; (void)mesh::make_linear_mesh_2d(ids(), points(), cells(), a); });
    rejects([&] { auto b = boundary; b.vertices = {LocalIndex{1}, LocalIndex{2}}; std::array a{b}; (void)mesh::make_linear_mesh_2d(ids(), points(), cells(), a); });
    rejects([&] { auto b = boundary; b.vertices = {LocalIndex{0}, LocalIndex{4}}; std::array a{b}; (void)mesh::make_linear_mesh_2d(ids(), points(), cells(), a); });
}

void linear_cell_mesh_2d_import_parity() {
    const auto direct = mesh::make_linear_mesh_2d(ids(), points(), cells());
    const auto gmsh = mesh::import_gmsh_4_1_ascii(gmsh_fixture(), 1.0);
    const auto vtu = mesh::import_vtu_ascii(vtu_fixture());
    check(gmsh.topology.global_ids(mesh::EntityKind::cell)[0].value() == 30U);
    check(vtu.topology.global_ids(mesh::EntityKind::cell)[0].value() == 900U);
    compare_semantics(direct.topology, direct.geometry, gmsh.topology, gmsh.geometry);
    compare_semantics(direct.topology, direct.geometry, vtu.topology, vtu.geometry);
    check_closure(gmsh.topology, gmsh.geometry);
    check_closure(vtu.topology, vtu.geometry);
    check(vtu.cell_fields[0].values()[0] == 90.0 && vtu.cell_fields[0].values()[1] == 3.0);
    check(vtu.point_fields[0].values()[4] == 50.0);
    const auto gmsh_again = mesh::import_gmsh_4_1_ascii(mesh::export_gmsh_4_1_ascii(gmsh), 1.0);
    const auto vtu_again = mesh::import_vtu_ascii(mesh::export_vtu_ascii(vtu));
    compare_semantics(gmsh.topology, gmsh.geometry, gmsh_again.topology, gmsh_again.geometry);
    compare_semantics(vtu.topology, vtu.geometry, vtu_again.topology, vtu_again.geometry);
    // Both format adapters must reach the common geometry rejection, not just
    // agree for valid rectangles. The fifth point now lies inside the quad.
    auto bad_gmsh = gmsh_fixture();
    auto bad_vtu = vtu_fixture();
    replace_once(bad_gmsh, "3 0 0\n", "1 0 0\n");
    replace_once(bad_vtu, "0 1 0 3 0 0", "0 1 0 1 0 0");
    rejects([&] { (void)mesh::import_gmsh_4_1_ascii(bad_gmsh, 1.0); });
    rejects([&] { (void)mesh::import_vtu_ascii(bad_vtu); });
    auto shifted_gmsh = gmsh_fixture();
    auto shifted_vtu = vtu_fixture();
    const std::string shifted_points =
        "100000000 100000000 0 100000002 100000000 0 100000002 100000001 0 "
        "100000000 100000001 0 100000003 100000000 0";
    replace_once(shifted_gmsh, "0 0 0\n2 0 0\n2 1 0\n0 1 0\n3 0 0", shifted_points);
    replace_once(shifted_vtu, "0 0 0 2 0 0 2 1 0 0 1 0 3 0 0", shifted_points);
    const auto gs = mesh::import_gmsh_4_1_ascii(shifted_gmsh, 1.0);
    const auto vs = mesh::import_vtu_ascii(shifted_vtu);
    compare_semantics(gs.topology, gs.geometry, vs.topology, vs.geometry);
    near(gs.geometry.cell_area_m2(LocalIndex{0}), 0.5);
    near(vs.geometry.cell_area_m2(LocalIndex{0}), 2.0);
}
