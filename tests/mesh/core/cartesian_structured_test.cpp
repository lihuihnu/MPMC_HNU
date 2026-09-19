#include <mpmc/mesh/cartesian_1d.hpp>
#include <mpmc/mesh/cartesian_3d.hpp>
#include <mpmc/mesh/geometry_1d.hpp>

#include <array>
#include <cmath>
#include <cstddef>
#include <iostream>
#include <span>
#include <stdexcept>

namespace mesh = mpmc::mesh;

namespace {

void require(
    bool condition,
    const char* message) {
    if (!condition) {
        throw std::runtime_error(message);
    }
}

void require_close(
    double actual,
    double expected,
    double tolerance,
    const char* message) {
    if (!std::isfinite(actual) ||
        std::abs(actual - expected) >
            tolerance) {
        throw std::runtime_error(message);
    }
}

void require_row(
    std::span<const mesh::LocalIndex> actual,
    std::span<const std::size_t> expected,
    const char* message) {
    if (actual.size() != expected.size()) {
        throw std::runtime_error(message);
    }
    for (std::size_t i = 0U;
         i < expected.size();
         ++i) {
        if (static_cast<std::size_t>(
                actual[i].value()) !=
            expected[i]) {
            throw std::runtime_error(message);
        }
    }
}

void verify_cartesian_1d() {
    const auto topology =
        mesh::make_cartesian_topology_1d(
            3U);
    require(
        topology.entity_count(
            mesh::EntityKind::vertex) == 4U &&
            topology.entity_count(
                mesh::EntityKind::edge) == 0U &&
            topology.entity_count(
                mesh::EntityKind::face) == 4U &&
            topology.entity_count(
                mesh::EntityKind::cell) == 3U,
        "1D Cartesian entity counts");

    const auto& cell_vertices =
        topology.relation(
            mesh::EntityKind::cell,
            mesh::EntityKind::vertex);
    const auto& cell_faces =
        topology.relation(
            mesh::EntityKind::cell,
            mesh::EntityKind::face);
    const auto& face_vertices =
        topology.relation(
            mesh::EntityKind::face,
            mesh::EntityKind::vertex);
    const auto& face_cells =
        topology.relation(
            mesh::EntityKind::face,
            mesh::EntityKind::cell);

    const std::array<std::size_t, 2>
        middle_vertices{1U, 2U};
    require_row(
        cell_vertices.adjacent(
            mesh::LocalIndex{1U}),
        middle_vertices,
        "1D middle cell vertices");
    require_row(
        cell_faces.adjacent(
            mesh::LocalIndex{1U}),
        middle_vertices,
        "1D middle cell faces");
    const std::array<std::size_t, 1>
        face_vertex{2U};
    require_row(
        face_vertices.adjacent(
            mesh::LocalIndex{2U}),
        face_vertex,
        "1D face point identity");
    const std::array<std::size_t, 2>
        internal_cells{1U, 2U};
    require_row(
        face_cells.adjacent(
            mesh::LocalIndex{2U}),
        internal_cells,
        "1D internal face owner/neighbour");

    const std::array<double, 4>
        x{-2.0, 0.5, 3.5, 10.0};
    const auto geometry =
        mesh::make_cartesian_geometry_1d(
            topology,
            x);
    const std::array<double, 3>
        lengths{2.5, 3.0, 6.5};
    const std::array<double, 3>
        centroids{-0.75, 2.0, 6.75};
    for (std::size_t cell = 0U;
         cell < 3U;
         ++cell) {
        const auto local =
            mesh::LocalIndex{
                static_cast<
                    mesh::LocalIndex::value_type>(
                    cell)};
        require_close(
            geometry.cell_length_m(local),
            lengths[cell],
            1.0e-14,
            "1D cell length");
        require_close(
            geometry.cell_centroid_m(local),
            centroids[cell],
            1.0e-14,
            "1D cell centroid");
    }
    require_close(
        geometry.face_coordinate_m(
            mesh::LocalIndex{0U}),
        -2.0,
        1.0e-14,
        "1D left face coordinate");
    require_close(
        geometry.face_owner_unit_normal(
            mesh::LocalIndex{0U}),
        -1.0,
        0.0,
        "1D left outward normal");
    require_close(
        geometry.face_owner_unit_normal(
            mesh::LocalIndex{2U}),
        1.0,
        0.0,
        "1D internal owner normal");
    require(
        geometry.face_owner(
            mesh::LocalIndex{2U}) ==
            mesh::LocalIndex{1U},
        "1D internal face owner");
}

void verify_cartesian_3d() {
    const std::array<double, 3>
        x{0.0, 1.0, 3.0};
    const std::array<double, 2>
        y{-1.0, 2.0};
    const std::array<double, 3>
        z{0.0, 2.0, 5.0};
    const auto grid =
        mesh::make_cartesian_mesh_3d(
            x,
            y,
            z);

    require(
        grid.topology.entity_count(
            mesh::EntityKind::vertex) == 18U &&
            grid.topology.entity_count(
                mesh::EntityKind::edge) == 0U &&
            grid.topology.entity_count(
                mesh::EntityKind::face) == 20U &&
            grid.topology.entity_count(
                mesh::EntityKind::cell) == 4U,
        "3D Cartesian entity counts");
    require(
        grid.face_boundary
                .boundary_face_count() == 16U &&
            grid.face_boundary
                .interior_face_count() == 4U,
        "3D Cartesian boundary counts");

    const std::array<double, 4>
        volumes{6.0, 12.0, 9.0, 18.0};
    for (std::size_t cell = 0U;
         cell < volumes.size();
         ++cell) {
        require_close(
            grid.cell_volumes_m3[cell],
            volumes[cell],
            1.0e-14,
            "3D Cartesian cell volume");
    }

    const auto& cell_vertices =
        grid.topology.relation(
            mesh::EntityKind::cell,
            mesh::EntityKind::vertex);
    const auto& cell_faces =
        grid.topology.relation(
            mesh::EntityKind::cell,
            mesh::EntityKind::face);
    const std::array<std::size_t, 8>
        first_vertices{
            0U, 1U, 4U, 3U,
            6U, 7U, 10U, 9U};
    const std::array<std::size_t, 6>
        first_faces{
            0U, 1U, 6U, 8U, 14U, 16U};
    require_row(
        cell_vertices.adjacent(
            mesh::LocalIndex{0U}),
        first_vertices,
        "3D Cartesian first-cell vertices");
    require_row(
        cell_faces.adjacent(
            mesh::LocalIndex{0U}),
        first_faces,
        "3D Cartesian first-cell faces");

    const auto x_boundary =
        mesh::LocalIndex{0U};
    const auto x_internal =
        mesh::LocalIndex{1U};
    const auto y_boundary =
        mesh::LocalIndex{6U};
    const auto z_boundary =
        mesh::LocalIndex{14U};

    const auto x0 =
        grid.face_geometry
            .face_centroid_m(
                x_boundary);
    require_close(
        x0.x_m,
        0.0,
        1.0e-14,
        "3D x-face centroid x");
    require_close(
        x0.y_m,
        0.5,
        1.0e-14,
        "3D x-face centroid y");
    require_close(
        x0.z_m,
        1.0,
        1.0e-14,
        "3D x-face centroid z");
    require_close(
        grid.face_geometry
            .face_area_m2(
                x_boundary),
        6.0,
        1.0e-14,
        "3D x-face area");
    const auto x0_normal =
        grid.face_geometry
            .face_owner_unit_normal(
                x_boundary);
    require_close(
        x0_normal.x,
        -1.0,
        0.0,
        "3D left x normal");

    const auto xi_normal =
        grid.face_geometry
            .face_owner_unit_normal(
                x_internal);
    require_close(
        xi_normal.x,
        1.0,
        0.0,
        "3D internal x owner normal");
    require(
        grid.face_geometry.face_owner(
            x_internal) ==
            mesh::LocalIndex{0U},
        "3D internal x face owner");

    require_close(
        grid.face_geometry
            .face_area_m2(
                y_boundary),
        2.0,
        1.0e-14,
        "3D y-face area");
    const auto y0_normal =
        grid.face_geometry
            .face_owner_unit_normal(
                y_boundary);
    require_close(
        y0_normal.y,
        -1.0,
        0.0,
        "3D lower y normal");

    require_close(
        grid.face_geometry
            .face_area_m2(
                z_boundary),
        3.0,
        1.0e-14,
        "3D z-face area");
    const auto z0_normal =
        grid.face_geometry
            .face_owner_unit_normal(
                z_boundary);
    require_close(
        z0_normal.z,
        -1.0,
        0.0,
        "3D lower z normal");

    const auto last_vertex =
        grid.vertex_coordinates_m.back();
    require_close(
        last_vertex.x_m,
        3.0,
        0.0,
        "3D final vertex x");
    require_close(
        last_vertex.y_m,
        2.0,
        0.0,
        "3D final vertex y");
    require_close(
        last_vertex.z_m,
        5.0,
        0.0,
        "3D final vertex z");
}

void verify_invalid() {
    bool rejected = false;
    try {
        (void)mesh::make_cartesian_topology_1d(
            0U);
    } catch (const std::invalid_argument&) {
        rejected = true;
    }
    require(
        rejected,
        "1D zero-cell grid must be rejected");

    rejected = false;
    try {
        (void)mesh::make_cartesian_topology_3d(
            1U,
            0U,
            1U);
    } catch (const std::invalid_argument&) {
        rejected = true;
    }
    require(
        rejected,
        "3D zero dimension must be rejected");

    const auto one_d_topology =
        mesh::make_cartesian_topology_1d(
            2U);
    const std::array<double, 3>
        bad_x{0.0, 1.0, 1.0};
    rejected = false;
    try {
        (void)mesh::make_cartesian_geometry_1d(
            one_d_topology,
            bad_x);
    } catch (const std::invalid_argument&) {
        rejected = true;
    }
    require(
        rejected,
        "1D non-increasing axis must be rejected");

    const std::array<double, 2>
        x{0.0, 1.0};
    const std::array<double, 2>
        y{0.0, 1.0};
    const std::array<double, 2>
        bad_z{2.0, 2.0};
    rejected = false;
    try {
        (void)mesh::make_cartesian_mesh_3d(
            x,
            y,
            bad_z);
    } catch (const std::invalid_argument&) {
        rejected = true;
    }
    require(
        rejected,
        "3D non-increasing axis must be rejected");
}

} // namespace

int main() {
    try {
        verify_cartesian_1d();
        verify_cartesian_3d();
        verify_invalid();
        std::cout
            << "[PASS] mesh.core.cartesian_structured\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr
            << "[FAIL] mesh.core.cartesian_structured: "
            << error.what()
            << '\n';
        return 1;
    }
}
