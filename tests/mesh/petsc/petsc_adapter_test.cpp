#include <mpmc/mesh/active_corner_point.hpp>
#include <mpmc/mesh/cartesian_2d.hpp>
#include <mpmc/mesh/dense_field.hpp>
#include <mpmc/mesh/dof_layout.hpp>
#include <mpmc/mesh/face_boundary.hpp>
#include <mpmc/mesh/geometry_2d.hpp>
#include <mpmc/mesh/gmsh_4_1.hpp>
#include <mpmc/mesh/grdecl.hpp>
#include <mpmc/mesh/linear_cell_mesh_3d.hpp>
#include <mpmc/mesh/dof_numbering.hpp>
#include <mpmc/mesh/partition_snapshot.hpp>
#include <mpmc/mesh/shared_entity_plan.hpp>
#include <mpmc/mesh/topology.hpp>
#include <mpmc/mesh_petsc/adapter.hpp>
#include <mpmc/discretization/transmissibility_admissibility_3d.hpp>
#include <mpmc/discretization_petsc/adapter.hpp>

#include <petscsys.h>

#include <algorithm>
#include <array>
#include <cstdint>
#include <cmath>
#include <cstring>
#include <exception>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

namespace {
namespace mesh = mpmc::mesh;
namespace mesh_petsc = mpmc::mesh_petsc;
namespace discretization = mpmc::discretization;
namespace discretization_petsc = mpmc::discretization_petsc;

void require(bool condition, std::string_view message) {
    if (!condition) {
        throw std::runtime_error(std::string(message));
    }
}

void require_petsc(PetscErrorCode error, std::string_view message) {
    if (error != PETSC_SUCCESS) {
        throw std::runtime_error(
            std::string(message) + " PETSc error=" +
            std::to_string(static_cast<long long>(error)));
    }
}

mesh::Topology two_by_one_cartesian_with_stable_ids() {
    const auto base = mesh::make_cartesian_topology_2d(2U, 1U);

    mesh::Topology::EntityIds ids;
    for (std::size_t i = 0U;
         i < base.entity_count(mesh::EntityKind::vertex);
         ++i) {
        ids.vertices.emplace_back(
            5000000000ULL + static_cast<std::uint64_t>(i));
    }
    for (std::size_t i = 0U;
         i < base.entity_count(mesh::EntityKind::face);
         ++i) {
        ids.faces.emplace_back(
            6000000000ULL + static_cast<std::uint64_t>(i));
    }
    for (std::size_t i = 0U;
         i < base.entity_count(mesh::EntityKind::cell);
         ++i) {
        ids.cells.emplace_back(
            7000000000ULL + static_cast<std::uint64_t>(i));
    }

    std::vector<mesh::CsrAdjacency> relations;
    relations.emplace_back(
        base.relation(
            mesh::EntityKind::cell,
            mesh::EntityKind::vertex));
    relations.emplace_back(
        base.relation(
            mesh::EntityKind::cell,
            mesh::EntityKind::face));
    relations.emplace_back(
        base.relation(
            mesh::EntityKind::face,
            mesh::EntityKind::vertex));
    relations.emplace_back(
        base.relation(
            mesh::EntityKind::face,
            mesh::EntityKind::cell));

    return mesh::Topology{
        std::move(ids), std::move(relations)};
}


mesh::ActiveCornerPointGrid processed_grdecl_two_by_one_all_active() {
    constexpr std::string_view deck = R"grdecl(
SPECGRID
  2 1 1 1 F /
COORD
  0 0 0   0 0 1
  1 0 0   1 0 1
  2 0 0   2 0 1
  0 1 0   0 1 1
  1 1 0   1 1 1
  2 1 0   2 1 1 /
ZCORN
  8*0 8*1 /
ACTNUM
  2*1 /
PORO
  0.20 0.35 /
PERMX
  100 200 /
PERMY
  50 75 /
PERMZ
  10 20 /
)grdecl";

    const auto raw = mesh::import_grdecl(
        deck,
        mesh::GrdeclImportOptions{
            1.0,
            1.0e-15});
    return mesh::process_active_corner_point_grid(raw);
}


struct SerialCoordinateView3D {
    PetscSection section;
    Vec values;
    const PetscScalar* array;
};

SerialCoordinateView3D get_serial_coordinate_view_3d(DM dm) {
    PetscInt coordinate_dim = -1;
    require_petsc(
        DMGetCoordinateDim(dm, &coordinate_dim),
        "serial 3D DMGetCoordinateDim");
    require(
        coordinate_dim == 3,
        "serial processed GRDECL coordinate dimension");

    PetscSection section = nullptr;
    Vec values = nullptr;
    require_petsc(
        DMGetCoordinateSection(dm, &section),
        "serial 3D DMGetCoordinateSection");
    require(
        section != nullptr,
        "serial 3D coordinate section");
    require_petsc(
        DMGetCoordinatesLocal(dm, &values),
        "serial 3D DMGetCoordinatesLocal");
    require(
        values != nullptr,
        "serial 3D local coordinate vector");

    const PetscScalar* array = nullptr;
    require_petsc(
        VecGetArrayRead(values, &array),
        "serial 3D VecGetArrayRead coordinates");
    return SerialCoordinateView3D{
        section, values, array};
}

void restore_serial_coordinate_view_3d(
    SerialCoordinateView3D* view) {
    require(
        view != nullptr,
        "serial 3D coordinate view pointer");
    require_petsc(
        VecRestoreArrayRead(
            view->values, &view->array),
        "serial 3D VecRestoreArrayRead coordinates");
}

mesh::Coordinate3D serial_vertex_coordinate_3d(
    const SerialCoordinateView3D& view,
    PetscInt point) {
    PetscInt dof = -1;
    PetscInt offset = -1;
    require_petsc(
        PetscSectionGetDof(
            view.section, point, &dof),
        "serial 3D vertex coordinate dof");
    require_petsc(
        PetscSectionGetOffset(
            view.section, point, &offset),
        "serial 3D vertex coordinate offset");
    require(
        dof == 3 && offset >= 0,
        "serial DMPlex vertex must carry three coordinate DoFs");
    return mesh::Coordinate3D{
        static_cast<double>(
            PetscRealPart(
                view.array[
                    static_cast<std::size_t>(
                        offset)])),
        static_cast<double>(
            PetscRealPart(
                view.array[
                    static_cast<std::size_t>(
                        offset + 1)])),
        static_cast<double>(
            PetscRealPart(
                view.array[
                    static_cast<std::size_t>(
                        offset + 2)]))};
}

struct TestVector3D {
    double x;
    double y;
    double z;
};

TestVector3D subtract(
    mesh::Coordinate3D left,
    mesh::Coordinate3D right) {
    return TestVector3D{
        left.x_m - right.x_m,
        left.y_m - right.y_m,
        left.z_m - right.z_m};
}

TestVector3D cross(
    TestVector3D left,
    TestVector3D right) {
    return TestVector3D{
        left.y * right.z - left.z * right.y,
        left.z * right.x - left.x * right.z,
        left.x * right.y - left.y * right.x};
}

double dot(
    TestVector3D left,
    TestVector3D right) {
    return left.x * right.x +
           left.y * right.y +
           left.z * right.z;
}

double magnitude(TestVector3D value) {
    return std::sqrt(dot(value, value));
}

TestVector3D scaled(
    TestVector3D value,
    double factor) {
    return TestVector3D{
        value.x * factor,
        value.y * factor,
        value.z * factor};
}

mesh::Coordinate3D triangle_centroid(
    mesh::Coordinate3D a,
    mesh::Coordinate3D b,
    mesh::Coordinate3D c) {
    return mesh::Coordinate3D{
        (a.x_m + b.x_m + c.x_m) / 3.0,
        (a.y_m + b.y_m + c.y_m) / 3.0,
        (a.z_m + b.z_m + c.z_m) / 3.0};
}

TestVector3D triangle_area_vector(
    mesh::Coordinate3D a,
    mesh::Coordinate3D b,
    mesh::Coordinate3D c) {
    return scaled(
        cross(subtract(b, a), subtract(c, a)),
        0.5);
}

mesh::Coordinate3D serial_cell_vertex_mean(
    DM dm,
    const SerialCoordinateView3D& view,
    PetscInt cell) {
    PetscInt vertex_start = -1;
    PetscInt vertex_end = -1;
    require_petsc(
        DMPlexGetDepthStratum(
            dm, 0, &vertex_start, &vertex_end),
        "serial 3D vertex stratum for cell centroid");

    PetscInt closure_size = 0;
    PetscInt* closure = nullptr;
    require_petsc(
        DMPlexGetTransitiveClosure(
            dm, cell, PETSC_TRUE,
            &closure_size, &closure),
        "serial 3D cell transitive closure");

    std::vector<PetscInt> vertices;
    vertices.reserve(8U);
    for (PetscInt i = 0;
         i < closure_size;
         ++i) {
        const PetscInt point = closure[2 * i];
        if (point < vertex_start ||
            point >= vertex_end) {
            continue;
        }
        if (std::find(
                vertices.begin(),
                vertices.end(),
                point) == vertices.end()) {
            vertices.push_back(point);
        }
    }
    require_petsc(
        DMPlexRestoreTransitiveClosure(
            dm, cell, PETSC_TRUE,
            &closure_size, &closure),
        "serial 3D restore cell transitive closure");
    require(
        vertices.size() == 8U,
        "serial hexa closure must contain eight unique vertices");

    mesh::Coordinate3D mean{
        0.0, 0.0, 0.0};
    for (const PetscInt vertex : vertices) {
        const auto coordinate =
            serial_vertex_coordinate_3d(
                view, vertex);
        mean.x_m += coordinate.x_m;
        mean.y_m += coordinate.y_m;
        mean.z_m += coordinate.z_m;
    }
    mean.x_m /= 8.0;
    mean.y_m /= 8.0;
    mean.z_m /= 8.0;
    return mean;
}

struct SerialFaceMetric3D {
    mesh::Coordinate3D centroid;
    double area;
    mesh::UnitVector3D owner_unit_normal;
};

SerialFaceMetric3D serial_face_metric_3d(
    DM dm,
    const SerialCoordinateView3D& view,
    PetscInt face,
    PetscInt owner_cell) {
    PetscInt cone_size = -1;
    const PetscInt* cone = nullptr;
    require_petsc(
        DMPlexGetConeSize(
            dm, face, &cone_size),
        "serial 3D face cone size for metric");
    require_petsc(
        DMPlexGetCone(dm, face, &cone),
        "serial 3D face cone for metric");
    require(
        cone_size == 4 && cone != nullptr,
        "serial 3D metric face must be a quad");

    std::array<mesh::Coordinate3D, 4> p{};
    for (std::size_t i = 0U;
         i < p.size();
         ++i) {
        p[i] = serial_vertex_coordinate_3d(
            view, cone[i]);
    }

    const auto area_vector0 =
        triangle_area_vector(
            p[0], p[1], p[2]);
    const auto area_vector1 =
        triangle_area_vector(
            p[0], p[2], p[3]);
    const double area0 =
        magnitude(area_vector0);
    const double area1 =
        magnitude(area_vector1);
    require(
        area0 > 0.0 && area1 > 0.0,
        "serial 3D face triangles must be nondegenerate");

    const auto centroid0 =
        triangle_centroid(
            p[0], p[1], p[2]);
    const auto centroid1 =
        triangle_centroid(
            p[0], p[2], p[3]);
    const double area = area0 + area1;
    const mesh::Coordinate3D centroid{
        (area0 * centroid0.x_m +
         area1 * centroid1.x_m) / area,
        (area0 * centroid0.y_m +
         area1 * centroid1.y_m) / area,
        (area0 * centroid0.z_m +
         area1 * centroid1.z_m) / area};

    TestVector3D total{
        area_vector0.x + area_vector1.x,
        area_vector0.y + area_vector1.y,
        area_vector0.z + area_vector1.z};
    const double total_magnitude =
        magnitude(total);
    require(
        total_magnitude > 0.0,
        "serial 3D quad resultant area vector");
    total = scaled(
        total, 1.0 / total_magnitude);

    const auto owner_centroid =
        serial_cell_vertex_mean(
            dm, view, owner_cell);
    if (dot(
            total,
            subtract(
                centroid,
                owner_centroid)) < 0.0) {
        total = scaled(total, -1.0);
    }

    return SerialFaceMetric3D{
        centroid,
        area,
        mesh::UnitVector3D{
            total.x, total.y, total.z}};
}

double serial_cell_volume_3d(
    DM dm,
    const SerialCoordinateView3D& view,
    PetscInt cell) {
    const auto cell_centroid =
        serial_cell_vertex_mean(
            dm, view, cell);

    PetscInt cone_size = -1;
    const PetscInt* faces = nullptr;
    require_petsc(
        DMPlexGetConeSize(
            dm, cell, &cone_size),
        "serial 3D cell cone size for volume");
    require_petsc(
        DMPlexGetCone(dm, cell, &faces),
        "serial 3D cell cone for volume");
    require(
        cone_size == 6 && faces != nullptr,
        "serial hexa volume requires six quad faces");

    double signed_volume = 0.0;
    for (PetscInt f = 0;
         f < cone_size;
         ++f) {
        PetscInt face_cone_size = -1;
        const PetscInt* face_vertices = nullptr;
        require_petsc(
            DMPlexGetConeSize(
                dm, faces[f], &face_cone_size),
            "serial 3D volume face cone size");
        require_petsc(
            DMPlexGetCone(
                dm, faces[f], &face_vertices),
            "serial 3D volume face cone");
        require(
            face_cone_size == 4 &&
                face_vertices != nullptr,
            "serial hexa volume face must be quad");

        std::array<mesh::Coordinate3D, 4> p{};
        for (std::size_t i = 0U;
             i < p.size();
             ++i) {
            p[i] =
                serial_vertex_coordinate_3d(
                    view, face_vertices[i]);
        }

        for (const auto triangle :
             {std::array<std::size_t, 3>{0U, 1U, 2U},
              std::array<std::size_t, 3>{0U, 2U, 3U}}) {
            const auto a = p[triangle[0]];
            const auto b = p[triangle[1]];
            const auto c = p[triangle[2]];
            auto area_vector =
                triangle_area_vector(a, b, c);
            const auto centroid =
                triangle_centroid(a, b, c);
            if (dot(
                    area_vector,
                    subtract(
                        centroid,
                        cell_centroid)) < 0.0) {
                area_vector =
                    scaled(area_vector, -1.0);
            }
            signed_volume +=
                (centroid.x_m * area_vector.x +
                 centroid.y_m * area_vector.y +
                 centroid.z_m * area_vector.z) /
                3.0;
        }
    }
    return std::abs(signed_volume);
}

void verify_processed_grdecl_3d_identity_coordinates(
    DM dm,
    const std::vector<mesh_petsc::DMPlexPointIdentity>& identities,
    const std::array<double, 36>& expected_vertex_coordinates) {
    PetscInt dimension = -1;
    PetscInt depth = -1;
    require_petsc(
        DMGetDimension(dm, &dimension),
        "processed GRDECL distributed DM dimension");
    require_petsc(
        DMPlexGetDepth(dm, &depth),
        "processed GRDECL distributed DMPlex depth");
    require(
        dimension == 3 && depth == 2,
        "processed GRDECL distributed mesh must remain 3D partial interpolation");

    auto coordinate_view =
        get_serial_coordinate_view_3d(dm);
    for (const auto& identity : identities) {
        switch (identity.kind) {
        case mesh::EntityKind::cell:
            require(
                identity.global.value() >= 1U &&
                    identity.global.value() <= 2U,
                "processed GRDECL distributed cell stable ID range");
            break;
        case mesh::EntityKind::face:
            require(
                identity.global.value() >= 1U &&
                    identity.global.value() <= 11U,
                "processed GRDECL distributed face stable ID range");
            break;
        case mesh::EntityKind::vertex: {
            require(
                identity.global.value() >= 1U &&
                    identity.global.value() <= 12U,
                "processed GRDECL distributed vertex stable ID range");
            const auto actual =
                serial_vertex_coordinate_3d(
                    coordinate_view,
                    identity.point);
            const std::size_t slot =
                static_cast<std::size_t>(
                    identity.global.value() - 1U) *
                3U;
            constexpr double tolerance = 1.0e-12;
            require(
                std::abs(
                    actual.x_m -
                    expected_vertex_coordinates[slot]) <=
                        tolerance &&
                    std::abs(
                        actual.y_m -
                        expected_vertex_coordinates[slot + 1U]) <=
                        tolerance &&
                    std::abs(
                        actual.z_m -
                        expected_vertex_coordinates[slot + 2U]) <=
                        tolerance,
                "processed GRDECL distributed vertex coordinate by stable ID");
            break;
        }
        case mesh::EntityKind::edge:
            throw std::runtime_error(
                "processed GRDECL distributed gate must not invent edge identities");
        }
    }
    restore_serial_coordinate_view_3d(
        &coordinate_view);
}

void require_processed_grdecl_identity_owner_counts(
    const mesh::PartitionSnapshot& partition) {
    std::array<int, 2> cells{0, 0};
    std::array<int, 11> faces{
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
    std::array<int, 12> vertices{
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};

    const auto record =
        [&](mesh::EntityKind kind, auto& counts) {
            for (std::size_t local = 0U;
                 local < partition.entity_count(kind);
                 ++local) {
                const auto index =
                    mesh::LocalIndex{
                        static_cast<
                            mesh::LocalIndex::value_type>(
                                local)};
                const auto id =
                    partition.global_id(kind, index);
                require(
                    id.value() >= 1U &&
                        id.value() <= counts.size(),
                    "processed GRDECL partition stable ID range");
                if (partition.is_owned(kind, index)) {
                    ++counts[
                        static_cast<std::size_t>(
                            id.value() - 1U)];
                }
            }
        };

    record(mesh::EntityKind::cell, cells);
    record(mesh::EntityKind::face, faces);
    record(mesh::EntityKind::vertex, vertices);

    require(
        MPI_Allreduce(
            MPI_IN_PLACE,
            cells.data(),
            static_cast<int>(cells.size()),
            MPI_INT,
            MPI_SUM,
            PETSC_COMM_WORLD) == MPI_SUCCESS,
        "MPI_Allreduce processed GRDECL cell owners");
    require(
        MPI_Allreduce(
            MPI_IN_PLACE,
            faces.data(),
            static_cast<int>(faces.size()),
            MPI_INT,
            MPI_SUM,
            PETSC_COMM_WORLD) == MPI_SUCCESS,
        "MPI_Allreduce processed GRDECL face owners");
    require(
        MPI_Allreduce(
            MPI_IN_PLACE,
            vertices.data(),
            static_cast<int>(vertices.size()),
            MPI_INT,
            MPI_SUM,
            PETSC_COMM_WORLD) == MPI_SUCCESS,
        "MPI_Allreduce processed GRDECL vertex owners");

    for (const int count : cells) {
        require(
            count == 1,
            "each processed GRDECL stable cell ID must have exactly one owner");
    }
    for (const int count : faces) {
        require(
            count == 1,
            "each processed GRDECL stable face ID must have exactly one owner");
    }
    for (const int count : vertices) {
        require(
            count == 1,
            "each processed GRDECL stable vertex ID must have exactly one owner");
    }
}

void verify_serial_dmplex_processed_grdecl() {
    const auto processed =
        processed_grdecl_two_by_one_all_active();
    const auto& topology = processed.topology;

    require(
        topology.entity_count(mesh::EntityKind::cell) == 2U &&
        topology.entity_count(mesh::EntityKind::face) == 11U &&
        topology.entity_count(mesh::EntityKind::vertex) == 12U &&
        topology.entity_count(mesh::EntityKind::edge) == 0U,
        "processed GRDECL 2x1x1 entity counts");

    DM dm = nullptr;
    std::vector<mesh_petsc::DMPlexPointIdentity> identities;
    require_petsc(
        mesh_petsc::create_serial_dmplex_topology(
            topology, &dm, &identities),
        "create_serial_dmplex_topology processed GRDECL");
    require_petsc(
        mesh_petsc::attach_serial_vertex_coordinates_3d(
            dm,
            processed.vertex_coordinates_m,
            identities),
        "attach processed GRDECL serial 3D coordinates");

    PetscInt dimension = -1;
    PetscInt depth = -1;
    require_petsc(
        DMGetDimension(dm, &dimension),
        "processed GRDECL DM dimension");
    require_petsc(
        DMPlexGetDepth(dm, &depth),
        "processed GRDECL DMPlex depth");
    require(
        dimension == 3 && depth == 2,
        "processed GRDECL must form a 3D partially interpolated DMPlex");

    PetscInt chart_start = -1;
    PetscInt chart_end = -1;
    require_petsc(
        DMPlexGetChart(dm, &chart_start, &chart_end),
        "processed GRDECL DMPlex chart");
    require(
        chart_start == 0 && chart_end == 25,
        "processed GRDECL DMPlex chart size");

    PetscInt cell_start = -1;
    PetscInt cell_end = -1;
    PetscInt face_start = -1;
    PetscInt face_end = -1;
    PetscInt vertex_start = -1;
    PetscInt vertex_end = -1;
    require_petsc(
        DMPlexGetHeightStratum(
            dm, 0, &cell_start, &cell_end),
        "processed GRDECL cell stratum");
    require_petsc(
        DMPlexGetHeightStratum(
            dm, 1, &face_start, &face_end),
        "processed GRDECL face stratum");
    require_petsc(
        DMPlexGetDepthStratum(
            dm, 0, &vertex_start, &vertex_end),
        "processed GRDECL vertex stratum");
    require(
        cell_start == 0 && cell_end == 2,
        "processed GRDECL hexa cell points");
    require(
        face_start == 2 && face_end == 13,
        "processed GRDECL quad face points");
    require(
        vertex_start == 13 && vertex_end == 25,
        "processed GRDECL vertex points");

    require(
        identities.size() == 25U,
        "processed GRDECL DMPlex identity count");
    for (std::size_t point = 0U;
         point < identities.size();
         ++point) {
        const auto& identity = identities[point];
        require(
            identity.point == static_cast<PetscInt>(point),
            "processed GRDECL identity point alignment");

        mesh::EntityKind expected_kind =
            mesh::EntityKind::cell;
        std::size_t expected_local = point;
        if (point >= 13U) {
            expected_kind = mesh::EntityKind::vertex;
            expected_local = point - 13U;
        } else if (point >= 2U) {
            expected_kind = mesh::EntityKind::face;
            expected_local = point - 2U;
        }

        const auto local = mesh::LocalIndex{
            static_cast<mesh::LocalIndex::value_type>(
                expected_local)};
        require(
            identity.kind == expected_kind &&
            identity.local == local &&
            identity.global ==
                topology.global_id(expected_kind, local),
            "processed GRDECL DMPlex stable identity");
    }

    for (std::size_t cell = 0U; cell < 2U; ++cell) {
        const auto local = mesh::LocalIndex{
            static_cast<mesh::LocalIndex::value_type>(cell)};
        require(
            processed.source_logical_cell_id(local) ==
                topology.global_id(mesh::EntityKind::cell, local),
            "processed GRDECL cell identity must remain source logical identity");

        DMPolytopeType type = DM_POLYTOPE_UNKNOWN;
        require_petsc(
            DMPlexGetCellType(
                dm, static_cast<PetscInt>(cell), &type),
            "processed GRDECL hexa cell type");
        require(
            type == DM_POLYTOPE_HEXAHEDRON,
            "processed GRDECL cell must be hexahedron");
    }
    for (PetscInt face = face_start;
         face < face_end;
         ++face) {
        DMPolytopeType type = DM_POLYTOPE_UNKNOWN;
        require_petsc(
            DMPlexGetCellType(dm, face, &type),
            "processed GRDECL quad face type");
        require(
            type == DM_POLYTOPE_QUADRILATERAL,
            "processed GRDECL face must be quadrilateral");
    }
    for (PetscInt vertex = vertex_start;
         vertex < vertex_end;
         ++vertex) {
        DMPolytopeType type = DM_POLYTOPE_UNKNOWN;
        require_petsc(
            DMPlexGetCellType(dm, vertex, &type),
            "processed GRDECL point type");
        require(
            type == DM_POLYTOPE_POINT,
            "processed GRDECL vertex must be point");
    }

    const auto& cell_faces = topology.relation(
        mesh::EntityKind::cell,
        mesh::EntityKind::face);
    const auto& face_vertices = topology.relation(
        mesh::EntityKind::face,
        mesh::EntityKind::vertex);
    const auto& face_cells = topology.relation(
        mesh::EntityKind::face,
        mesh::EntityKind::cell);

    for (std::size_t cell = 0U; cell < 2U; ++cell) {
        const auto local = mesh::LocalIndex{
            static_cast<mesh::LocalIndex::value_type>(cell)};
        const auto expected_faces =
            cell_faces.adjacent(local);

        PetscInt cone_size = -1;
        const PetscInt* cone = nullptr;
        require_petsc(
            DMPlexGetConeSize(
                dm,
                static_cast<PetscInt>(cell),
                &cone_size),
            "processed GRDECL cell cone size");
        require_petsc(
            DMPlexGetCone(
                dm,
                static_cast<PetscInt>(cell),
                &cone),
            "processed GRDECL cell cone");
        require(
            cone_size == 6 && cone != nullptr,
            "processed GRDECL hexa must have six face cone points");
        for (std::size_t i = 0U;
             i < expected_faces.size();
             ++i) {
            require(
                cone[i] ==
                    2 + static_cast<PetscInt>(
                            expected_faces[i].value()),
                "processed GRDECL cell-to-face cone identity");
        }
    }

    std::size_t shared_face_count = 0U;
    for (std::size_t face = 0U; face < 11U; ++face) {
        const auto local = mesh::LocalIndex{
            static_cast<mesh::LocalIndex::value_type>(face)};
        const auto expected_vertices =
            face_vertices.adjacent(local);
        const auto expected_cells =
            face_cells.adjacent(local);

        PetscInt cone_size = -1;
        const PetscInt* cone = nullptr;
        require_petsc(
            DMPlexGetConeSize(
                dm,
                2 + static_cast<PetscInt>(face),
                &cone_size),
            "processed GRDECL face cone size");
        require_petsc(
            DMPlexGetCone(
                dm,
                2 + static_cast<PetscInt>(face),
                &cone),
            "processed GRDECL face cone");
        require(
            cone_size == 4 && cone != nullptr,
            "processed GRDECL quad must have four vertex cone points");
        for (std::size_t i = 0U; i < 4U; ++i) {
            require(
                cone[i] ==
                    13 + static_cast<PetscInt>(
                             expected_vertices[i].value()),
                "processed GRDECL face-to-vertex cone identity");
        }

        PetscInt support_size = -1;
        const PetscInt* support = nullptr;
        require_petsc(
            DMPlexGetSupportSize(
                dm,
                2 + static_cast<PetscInt>(face),
                &support_size),
            "processed GRDECL face support size");
        require_petsc(
            DMPlexGetSupport(
                dm,
                2 + static_cast<PetscInt>(face),
                &support),
            "processed GRDECL face support");
        require(
            support_size ==
                static_cast<PetscInt>(
                    expected_cells.size()),
            "processed GRDECL face-to-cell support width");

        std::vector<PetscInt> actual_support(
            support,
            support +
                static_cast<std::ptrdiff_t>(
                    support_size));
        std::vector<PetscInt> expected_support;
        for (const auto expected_cell : expected_cells) {
            expected_support.push_back(
                static_cast<PetscInt>(
                    expected_cell.value()));
        }
        std::sort(
            actual_support.begin(),
            actual_support.end());
        std::sort(
            expected_support.begin(),
            expected_support.end());
        require(
            actual_support == expected_support,
            "processed GRDECL face-to-cell support identity");

        if (support_size == 2) {
            ++shared_face_count;
        }
    }
    require(
        shared_face_count == 1U,
        "processed GRDECL DMPlex must preserve exactly one shared face");

    constexpr double geometry_tolerance =
        1.0e-12;
    auto coordinate_view =
        get_serial_coordinate_view_3d(dm);

    for (const auto& identity : identities) {
        PetscInt dof = -1;
        require_petsc(
            PetscSectionGetDof(
                coordinate_view.section,
                identity.point,
                &dof),
            "processed GRDECL coordinate point dof");

        if (identity.kind !=
            mesh::EntityKind::vertex) {
            require(
                dof == 0,
                "only processed DMPlex vertices may carry coordinates");
            continue;
        }

        require(
            dof == 3,
            "processed DMPlex vertex coordinate width");
        const auto actual =
            serial_vertex_coordinate_3d(
                coordinate_view,
                identity.point);
        const auto expected =
            processed.vertex_coordinates_m[
                static_cast<std::size_t>(
                    identity.local.value())];
        require(
            std::abs(actual.x_m - expected.x_m) <=
                    geometry_tolerance &&
                std::abs(actual.y_m - expected.y_m) <=
                    geometry_tolerance &&
                std::abs(actual.z_m - expected.z_m) <=
                    geometry_tolerance,
            "serial DMPlex coordinates must match processed GRDECL vertex_coordinates_m");
    }

    for (std::size_t face = 0U;
         face < 11U;
         ++face) {
        const auto local = mesh::LocalIndex{
            static_cast<
                mesh::LocalIndex::value_type>(
                    face)};
        const PetscInt point =
            2 + static_cast<PetscInt>(face);

        PetscInt support_size = -1;
        const PetscInt* support = nullptr;
        require_petsc(
            DMPlexGetSupportSize(
                dm, point, &support_size),
            "processed GRDECL metric face support size");
        require_petsc(
            DMPlexGetSupport(
                dm, point, &support),
            "processed GRDECL metric face support");
        require(
            support_size >= 1 &&
                support != nullptr,
            "processed GRDECL metric face owner support");

        const auto expected_owner =
            processed.face_geometry.face_owner(
                local);
        const auto owner_identity =
            std::find_if(
                identities.begin(),
                identities.end(),
                [expected_owner](const auto& identity) {
                    return identity.kind ==
                               mesh::EntityKind::cell &&
                           identity.local ==
                               expected_owner;
                });
        require(
            owner_identity != identities.end(),
            "processed GRDECL face owner must have DMPlex identity");
        require(
            std::find(
                support,
                support +
                    static_cast<std::ptrdiff_t>(
                        support_size),
                owner_identity->point) !=
                support +
                    static_cast<std::ptrdiff_t>(
                        support_size),
            "processed GRDECL face owner must be in DMPlex support");

        const auto actual =
            serial_face_metric_3d(
                dm,
                coordinate_view,
                point,
                owner_identity->point);
        const auto expected_centroid =
            processed.face_geometry
                .face_centroid_m(local);
        const double expected_area =
            processed.face_geometry
                .face_area_m2(local);
        const auto expected_normal =
            processed.face_geometry
                .face_owner_unit_normal(local);

        require(
            std::abs(
                actual.centroid.x_m -
                expected_centroid.x_m) <=
                    geometry_tolerance &&
                std::abs(
                    actual.centroid.y_m -
                    expected_centroid.y_m) <=
                    geometry_tolerance &&
                std::abs(
                    actual.centroid.z_m -
                    expected_centroid.z_m) <=
                    geometry_tolerance,
            "serial DMPlex face centroid must match FaceGeometry3D");
        require(
            std::abs(
                actual.area -
                expected_area) <=
                    geometry_tolerance,
            "serial DMPlex face area must match FaceGeometry3D");
        require(
            std::abs(
                actual.owner_unit_normal.x -
                expected_normal.x) <=
                    geometry_tolerance &&
                std::abs(
                    actual.owner_unit_normal.y -
                    expected_normal.y) <=
                    geometry_tolerance &&
                std::abs(
                    actual.owner_unit_normal.z -
                    expected_normal.z) <=
                    geometry_tolerance,
            "serial DMPlex owner-relative face normal must match FaceGeometry3D");
    }

    for (std::size_t cell = 0U;
         cell < 2U;
         ++cell) {
        const double actual_volume =
            serial_cell_volume_3d(
                dm,
                coordinate_view,
                static_cast<PetscInt>(cell));
        const double expected_volume =
            processed.cell_volumes_m3[cell];
        require(
            std::abs(
                actual_volume -
                expected_volume) <=
                    geometry_tolerance,
            "serial DMPlex cell volume must match processed cell_volumes_m3");
    }

    restore_serial_coordinate_view_3d(
        &coordinate_view);

    require_petsc(
        DMDestroy(&dm),
        "DMDestroy processed GRDECL serial DMPlex");
    require(
        dm == nullptr,
        "processed GRDECL DMPlex destroy must clear handle");
}

void verify_serial_dmplex_topology() {
    const auto topology =
        two_by_one_cartesian_with_stable_ids();

    DM dm = nullptr;
    std::vector<mesh_petsc::DMPlexPointIdentity> identities;
    require_petsc(
        mesh_petsc::create_serial_dmplex_topology(
            topology, &dm, &identities),
        "create_serial_dmplex_topology");

    PetscBool is_plex = PETSC_FALSE;
    require_petsc(
        PetscObjectTypeCompare(
            reinterpret_cast<PetscObject>(dm),
            DMPLEX,
            &is_plex),
        "DMPlex object type");
    require(is_plex == PETSC_TRUE,
            "serial topology adapter must create DMPLEX");

    MPI_Comm dm_comm = MPI_COMM_NULL;
    require_petsc(
        PetscObjectGetComm(
            reinterpret_cast<PetscObject>(dm),
            &dm_comm),
        "DMPlex communicator");
    int dm_size = -1;
    require(
        MPI_Comm_size(dm_comm, &dm_size) == MPI_SUCCESS,
        "DMPlex communicator size");
    require(dm_size == 1,
            "serial DMPlex must live on PETSC_COMM_SELF");

    PetscInt dimension = -1;
    PetscInt depth = -1;
    require_petsc(
        DMGetDimension(dm, &dimension),
        "DMGetDimension");
    require_petsc(
        DMPlexGetDepth(dm, &depth),
        "DMPlexGetDepth");
    require(dimension == 2 && depth == 2,
            "2D fully interpolated DMPlex depth");

    PetscInt chart_start = -1;
    PetscInt chart_end = -1;
    require_petsc(
        DMPlexGetChart(dm, &chart_start, &chart_end),
        "DMPlexGetChart");
    require(chart_start == 0 && chart_end == 15,
            "2x1 quad DMPlex chart");

    PetscInt cell_start = -1;
    PetscInt cell_end = -1;
    PetscInt face_start = -1;
    PetscInt face_end = -1;
    PetscInt vertex_start = -1;
    PetscInt vertex_end = -1;
    require_petsc(
        DMPlexGetHeightStratum(
            dm, 0, &cell_start, &cell_end),
        "DMPlex cell stratum");
    require_petsc(
        DMPlexGetHeightStratum(
            dm, 1, &face_start, &face_end),
        "DMPlex face stratum");
    require_petsc(
        DMPlexGetDepthStratum(
            dm, 0, &vertex_start, &vertex_end),
        "DMPlex vertex stratum");

    require(cell_start == 0 && cell_end == 2,
            "DMPlex cell points");
    require(face_start == 2 && face_end == 9,
            "DMPlex face points");
    require(vertex_start == 9 && vertex_end == 15,
            "DMPlex vertex points");

    require(identities.size() == 15U,
            "DMPlex identity map size");
    for (std::size_t point = 0U;
         point < identities.size();
         ++point) {
        const auto& identity = identities[point];
        require(
            identity.point ==
                static_cast<PetscInt>(point),
            "identity map point alignment");

        mesh::EntityKind expected_kind =
            mesh::EntityKind::cell;
        std::size_t expected_local = point;
        if (point >= 9U) {
            expected_kind = mesh::EntityKind::vertex;
            expected_local = point - 9U;
        } else if (point >= 2U) {
            expected_kind = mesh::EntityKind::face;
            expected_local = point - 2U;
        }

        require(identity.kind == expected_kind,
                "DMPlex stable identity kind");
        require(
            identity.local.value() ==
                static_cast<mesh::LocalIndex::value_type>(
                    expected_local),
            "DMPlex stable identity local index");
        require(
            identity.global ==
                topology.global_id(
                    expected_kind,
                    mesh::LocalIndex{
                        static_cast<
                            mesh::LocalIndex::value_type>(
                                expected_local)}),
            "DMPlex stable GlobalEntityId");
        require(identity.global.value() > 4000000000ULL,
                "stable ID regression must exercise 64-bit identity");
    }

    for (PetscInt cell = cell_start;
         cell < cell_end;
         ++cell) {
        DMPolytopeType type = DM_POLYTOPE_UNKNOWN;
        require_petsc(
            DMPlexGetCellType(dm, cell, &type),
            "DMPlex quadrilateral cell type");
        require(type == DM_POLYTOPE_QUADRILATERAL,
                "DMPlex cell must be quadrilateral");
    }
    for (PetscInt face = face_start;
         face < face_end;
         ++face) {
        DMPolytopeType type = DM_POLYTOPE_UNKNOWN;
        require_petsc(
            DMPlexGetCellType(dm, face, &type),
            "DMPlex segment face type");
        require(type == DM_POLYTOPE_SEGMENT,
                "DMPlex face must be segment");
    }
    for (PetscInt vertex = vertex_start;
         vertex < vertex_end;
         ++vertex) {
        DMPolytopeType type = DM_POLYTOPE_UNKNOWN;
        require_petsc(
            DMPlexGetCellType(dm, vertex, &type),
            "DMPlex point vertex type");
        require(type == DM_POLYTOPE_POINT,
                "DMPlex vertex must be point");
    }

    const auto& cell_faces = topology.relation(
        mesh::EntityKind::cell,
        mesh::EntityKind::face);
    const auto& face_vertices = topology.relation(
        mesh::EntityKind::face,
        mesh::EntityKind::vertex);
    const auto& face_cells = topology.relation(
        mesh::EntityKind::face,
        mesh::EntityKind::cell);

    for (std::size_t cell = 0U; cell < 2U; ++cell) {
        const auto local = mesh::LocalIndex{
            static_cast<mesh::LocalIndex::value_type>(
                cell)};
        const auto expected = cell_faces.adjacent(local);

        PetscInt cone_size = -1;
        const PetscInt* cone = nullptr;
        require_petsc(
            DMPlexGetConeSize(
                dm,
                static_cast<PetscInt>(cell),
                &cone_size),
            "DMPlex cell cone size");
        require_petsc(
            DMPlexGetCone(
                dm,
                static_cast<PetscInt>(cell),
                &cone),
            "DMPlex cell cone");
        require(cone_size == 4 && cone != nullptr,
                "DMPlex cell cone width");
        for (std::size_t i = 0U; i < 4U; ++i) {
            require(
                cone[i] ==
                    2 + static_cast<PetscInt>(
                        expected[i].value()),
                "DMPlex cell-to-face cone identity");
        }
    }

    for (std::size_t face = 0U; face < 7U; ++face) {
        const auto local = mesh::LocalIndex{
            static_cast<mesh::LocalIndex::value_type>(
                face)};
        const auto expected_vertices =
            face_vertices.adjacent(local);

        PetscInt cone_size = -1;
        const PetscInt* cone = nullptr;
        require_petsc(
            DMPlexGetConeSize(
                dm,
                2 + static_cast<PetscInt>(face),
                &cone_size),
            "DMPlex face cone size");
        require_petsc(
            DMPlexGetCone(
                dm,
                2 + static_cast<PetscInt>(face),
                &cone),
            "DMPlex face cone");
        require(cone_size == 2 && cone != nullptr,
                "DMPlex face cone width");
        require(
            cone[0] ==
                    9 + static_cast<PetscInt>(
                        expected_vertices[0].value()) &&
                cone[1] ==
                    9 + static_cast<PetscInt>(
                        expected_vertices[1].value()),
            "DMPlex face-to-vertex cone identity");

        PetscInt support_size = -1;
        const PetscInt* support = nullptr;
        require_petsc(
            DMPlexGetSupportSize(
                dm,
                2 + static_cast<PetscInt>(face),
                &support_size),
            "DMPlex face support size");
        require_petsc(
            DMPlexGetSupport(
                dm,
                2 + static_cast<PetscInt>(face),
                &support),
            "DMPlex face support");

        const auto expected_cells =
            face_cells.adjacent(local);
        require(
            support_size ==
                static_cast<PetscInt>(
                    expected_cells.size()),
            "DMPlex face-to-cell support width");

        std::vector<PetscInt> actual_support(
            support,
            support +
                static_cast<std::ptrdiff_t>(
                    support_size));
        std::vector<PetscInt> expected_support;
        for (const auto expected_cell :
             expected_cells) {
            expected_support.push_back(
                static_cast<PetscInt>(
                    expected_cell.value()));
        }
        std::sort(
            actual_support.begin(),
            actual_support.end());
        std::sort(
            expected_support.begin(),
            expected_support.end());
        require(
            actual_support == expected_support,
            "DMPlex face-to-cell support identity");
    }

    for (std::size_t vertex = 0U;
         vertex < 6U;
         ++vertex) {
        std::vector<PetscInt> expected_faces;
        for (std::size_t face = 0U;
             face < 7U;
             ++face) {
            const auto face_local = mesh::LocalIndex{
                static_cast<
                    mesh::LocalIndex::value_type>(
                        face)};
            for (const auto incident_vertex :
                 face_vertices.adjacent(face_local)) {
                if (incident_vertex.value() ==
                    static_cast<
                        mesh::LocalIndex::value_type>(
                            vertex)) {
                    expected_faces.push_back(
                        2 + static_cast<PetscInt>(
                            face));
                }
            }
        }

        PetscInt support_size = -1;
        const PetscInt* support = nullptr;
        require_petsc(
            DMPlexGetSupportSize(
                dm,
                9 + static_cast<PetscInt>(vertex),
                &support_size),
            "DMPlex vertex support size");
        require_petsc(
            DMPlexGetSupport(
                dm,
                9 + static_cast<PetscInt>(vertex),
                &support),
            "DMPlex vertex support");
        require(
            support_size ==
                static_cast<PetscInt>(
                    expected_faces.size()),
            "DMPlex vertex-to-face support width");

        std::vector<PetscInt> actual_support(
            support,
            support +
                static_cast<std::ptrdiff_t>(
                    support_size));
        std::sort(
            actual_support.begin(),
            actual_support.end());
        std::sort(
            expected_faces.begin(),
            expected_faces.end());
        require(
            actual_support == expected_faces,
            "DMPlex vertex-to-face support identity");
    }

    require_petsc(
        DMDestroy(&dm),
        "DMDestroy serial DMPlex");
    require(dm == nullptr,
            "DMDestroy must clear serial DMPlex handle");

    mesh::Topology::EntityIds invalid_ids;
    invalid_ids.vertices = {
        mesh::GlobalEntityId{1U},
        mesh::GlobalEntityId{2U},
        mesh::GlobalEntityId{3U},
        mesh::GlobalEntityId{4U}};
    invalid_ids.faces = {
        mesh::GlobalEntityId{5U}};
    invalid_ids.cells = {
        mesh::GlobalEntityId{6U}};
    const mesh::Topology missing_relations{
        std::move(invalid_ids), {}};

    identities.push_back(
        mesh_petsc::DMPlexPointIdentity{
            0,
            mesh::EntityKind::cell,
            mesh::LocalIndex{0U},
            mesh::GlobalEntityId{0U}});
    const PetscErrorCode invalid_error =
        mesh_petsc::create_serial_dmplex_topology(
            missing_relations, &dm, &identities);
    require(
        invalid_error == PETSC_ERR_ARG_INCOMP,
        "DMPlex adapter must reject missing core relations");
    require(dm == nullptr && identities.empty(),
            "failed DMPlex creation must leave clean outputs");
}

mesh::Topology two_rank_topology(std::uint32_t rank) {
    mesh::Topology::EntityIds ids;
    if (rank == 0U) {
        ids.vertices = {
            mesh::GlobalEntityId{10U}, mesh::GlobalEntityId{20U}};
        ids.faces = {
            mesh::GlobalEntityId{100U}, mesh::GlobalEntityId{200U}};
        ids.cells = {
            mesh::GlobalEntityId{1000U}, mesh::GlobalEntityId{2000U}};
    } else if (rank == 1U) {
        ids.vertices = {
            mesh::GlobalEntityId{20U}, mesh::GlobalEntityId{10U}};
        ids.faces = {
            mesh::GlobalEntityId{200U}, mesh::GlobalEntityId{100U}};
        ids.cells = {
            mesh::GlobalEntityId{2000U}, mesh::GlobalEntityId{1000U}};
    } else {
        throw std::invalid_argument("fixture rank must be 0 or 1");
    }
    return mesh::Topology{std::move(ids), {}};
}

mesh::PartitionSnapshot two_rank_partition(
    const mesh::Topology& topology,
    std::uint32_t rank) {
    mesh::EntityOwnerRanks owners;
    if (rank == 0U) {
        owners.vertices = {
            mesh::PartitionRank{0U}, mesh::PartitionRank{1U}};
        owners.faces = {
            mesh::PartitionRank{0U}, mesh::PartitionRank{1U}};
        owners.cells = {
            mesh::PartitionRank{0U}, mesh::PartitionRank{1U}};
    } else {
        owners.vertices = {
            mesh::PartitionRank{1U}, mesh::PartitionRank{0U}};
        owners.faces = {
            mesh::PartitionRank{1U}, mesh::PartitionRank{0U}};
        owners.cells = {
            mesh::PartitionRank{1U}, mesh::PartitionRank{0U}};
    }
    return mesh::PartitionSnapshot::create(
        topology, mesh::PartitionRank{rank}, 2U, std::move(owners));
}

std::vector<mesh::SharedEntityLink> shared_links() {
    return {
        {mesh::EntityKind::vertex, mesh::GlobalEntityId{10U},
         mesh::PartitionRank{0U}, mesh::LocalIndex{0U},
         mesh::PartitionRank{1U}, mesh::LocalIndex{1U}},
        {mesh::EntityKind::vertex, mesh::GlobalEntityId{20U},
         mesh::PartitionRank{1U}, mesh::LocalIndex{0U},
         mesh::PartitionRank{0U}, mesh::LocalIndex{1U}},
        {mesh::EntityKind::face, mesh::GlobalEntityId{100U},
         mesh::PartitionRank{0U}, mesh::LocalIndex{0U},
         mesh::PartitionRank{1U}, mesh::LocalIndex{1U}},
        {mesh::EntityKind::face, mesh::GlobalEntityId{200U},
         mesh::PartitionRank{1U}, mesh::LocalIndex{0U},
         mesh::PartitionRank{0U}, mesh::LocalIndex{1U}},
        {mesh::EntityKind::cell, mesh::GlobalEntityId{1000U},
         mesh::PartitionRank{0U}, mesh::LocalIndex{0U},
         mesh::PartitionRank{1U}, mesh::LocalIndex{1U}},
        {mesh::EntityKind::cell, mesh::GlobalEntityId{2000U},
         mesh::PartitionRank{1U}, mesh::LocalIndex{0U},
         mesh::PartitionRank{0U}, mesh::LocalIndex{1U}},
    };
}

mesh::GlobalEntityNumberingInput global_entity_numbering(
    const mesh::PartitionSnapshot& partition) {
    mesh::GlobalEntityNumberingInput input;
    input.global_cell_count = 2U;
    input.global_face_count = 2U;
    input.global_vertex_count = 2U;

    const auto append = [&](mesh::EntityKind kind,
                            std::vector<mesh::GlobalEntityOrdinalRecord>& records,
                            std::uint64_t first_id,
                            std::uint64_t second_id) {
        for (std::size_t local = 0U;
             local < partition.entity_count(kind);
             ++local) {
            const auto index = mesh::LocalIndex{
                static_cast<mesh::LocalIndex::value_type>(local)};
            const auto id = partition.global_id(kind, index);
            std::uint64_t ordinal = 0U;
            if (id.value() == first_id) {
                ordinal = 0U;
            } else if (id.value() == second_id) {
                ordinal = 1U;
            } else {
                throw std::runtime_error("unexpected fixture GlobalEntityId");
            }
            records.push_back(
                {id, mesh::GlobalEntityOrdinal{ordinal}});
        }
    };

    append(mesh::EntityKind::cell, input.cells, 1000U, 2000U);
    append(mesh::EntityKind::face, input.faces, 100U, 200U);
    append(mesh::EntityKind::vertex, input.vertices, 10U, 20U);
    return input;
}

void verify_section(
    const mesh::DofLayout& layout,
    const mesh::DofNumberingSnapshot& numbering,
    int mpi_rank) {
    PetscSection section = nullptr;
    std::vector<PetscInt> local_to_global;
    require_petsc(
        mesh_petsc::create_section_mapping(
            PETSC_COMM_WORLD, layout, numbering,
            &section, &local_to_global),
        "create_section_mapping");

    PetscInt chart_start = -1;
    PetscInt chart_end = -1;
    require_petsc(
        PetscSectionGetChart(section, &chart_start, &chart_end),
        "PetscSectionGetChart");
    require(chart_start == 0 && chart_end == 6,
            "PetscSection chart must be [0,6)");

    PetscInt field_count = 0;
    require_petsc(
        PetscSectionGetNumFields(section, &field_count),
        "PetscSectionGetNumFields");
    require(field_count == 4, "PetscSection field count");

    PetscBool point_major = PETSC_FALSE;
    require_petsc(
        PetscSectionGetPointMajor(section, &point_major),
        "PetscSectionGetPointMajor");
    require(point_major == PETSC_TRUE,
            "PetscSection must be explicitly point-major");

    const std::array<const char*, 4> expected_names{
        "cell.primary", "vertex.aux", "cell.secondary", "face.trace"};
    const std::array<PetscInt, 4> expected_components{2, 1, 1, 1};
    for (PetscInt field = 0; field < field_count; ++field) {
        const char* name = nullptr;
        PetscInt components = 0;
        require_petsc(
            PetscSectionGetFieldName(section, field, &name),
            "PetscSectionGetFieldName");
        require_petsc(
            PetscSectionGetFieldComponents(section, field, &components),
            "PetscSectionGetFieldComponents");
        require(name != nullptr &&
                    std::strcmp(name, expected_names[
                        static_cast<std::size_t>(field)]) == 0,
                "PetscSection field name");
        require(components == expected_components[
                    static_cast<std::size_t>(field)],
                "PetscSection field components");
    }

    const std::array<PetscInt, 6> expected_dofs{3, 3, 1, 1, 1, 1};
    const std::array<PetscInt, 6> expected_offsets{0, 3, 6, 7, 8, 9};
    for (PetscInt point = 0; point < 6; ++point) {
        PetscInt dof = 0;
        PetscInt offset = -1;
        require_petsc(
            PetscSectionGetDof(section, point, &dof),
            "PetscSectionGetDof");
        require_petsc(
            PetscSectionGetOffset(section, point, &offset),
            "PetscSectionGetOffset");
        require(dof == expected_dofs[
                    static_cast<std::size_t>(point)],
                "PetscSection point dof");
        require(offset == expected_offsets[
                    static_cast<std::size_t>(point)],
                "PetscSection point offset");
    }

    PetscInt field_dof = 0;
    PetscInt field_offset = -1;
    require_petsc(
        PetscSectionGetFieldDof(section, 0, 0, &field_dof),
        "cell primary field dof");
    require_petsc(
        PetscSectionGetFieldOffset(section, 0, 0, &field_offset),
        "cell primary field offset");
    require(field_dof == 2 && field_offset == 0,
            "cell primary field layout");
    require_petsc(
        PetscSectionGetFieldDof(section, 0, 2, &field_dof),
        "cell secondary field dof");
    require_petsc(
        PetscSectionGetFieldOffset(section, 0, 2, &field_offset),
        "cell secondary field offset");
    require(field_dof == 1 && field_offset == 2,
            "cell secondary field layout");
    require_petsc(
        PetscSectionGetFieldDof(section, 2, 3, &field_dof),
        "face trace field dof");
    require_petsc(
        PetscSectionGetFieldOffset(section, 2, 3, &field_offset),
        "face trace field offset");
    require(field_dof == 1 && field_offset == 6,
            "face trace field layout");
    require_petsc(
        PetscSectionGetFieldDof(section, 4, 1, &field_dof),
        "vertex aux field dof");
    require_petsc(
        PetscSectionGetFieldOffset(section, 4, 1, &field_offset),
        "vertex aux field offset");
    require(field_dof == 1 && field_offset == 8,
            "vertex aux field layout");

    require(local_to_global.size() == 10U,
            "PETSc local-to-global scalar map size");
    const std::array<PetscInt, 10> rank0_expected{
        0, 1, 2, 3, 4, 5, 6, 7, 8, 9};
    const std::array<PetscInt, 10> rank1_expected{
        3, 4, 5, 0, 1, 2, 7, 6, 9, 8};
    const auto& expected =
        mpi_rank == 0 ? rank0_expected : rank1_expected;
    for (std::size_t local = 0U;
         local < local_to_global.size();
         ++local) {
        require(local_to_global[local] == expected[local],
                "PETSc-compatible local-to-global scalar map");
    }

    require_petsc(PetscSectionDestroy(&section), "PetscSectionDestroy");
    require(section == nullptr, "PetscSectionDestroy must clear handle");
}

void verify_sf(
    const mesh::PartitionSnapshot& partition,
    const mesh::SharedEntityPlan& plan,
    int mpi_rank) {
    for (const auto kind :
         {mesh::EntityKind::cell,
          mesh::EntityKind::face,
          mesh::EntityKind::vertex}) {
        PetscSF sf = nullptr;
        require_petsc(
            mesh_petsc::create_entity_sf(
                PETSC_COMM_WORLD, partition, plan, kind, &sf),
            "create_entity_sf");

        PetscInt nroots = -1;
        PetscInt nleaves = -1;
        const PetscInt* ilocal = nullptr;
        const PetscSFNode* remote = nullptr;
        require_petsc(
            PetscSFGetGraph(
                sf, &nroots, &nleaves, &ilocal, &remote),
            "PetscSFGetGraph");

        require(nroots == 2, "PetscSF root count");
        require(nleaves == 1, "PetscSF leaf count");
        require(ilocal != nullptr && ilocal[0] == 1,
                "PetscSF local ghost leaf index");
        require(remote != nullptr,
                "PetscSF remote root array");
        require(remote[0].rank == (mpi_rank == 0 ? 1 : 0),
                "PetscSF remote owner rank");
        require(remote[0].index == 0,
                "PetscSF remote owner-local root index");

        if (kind == mesh::EntityKind::cell) {
            std::array<PetscInt, 2> root_values{
                mpi_rank == 0 ? 1000 : 2000,
                -77};
            std::array<PetscInt, 2> leaf_values{-1, -1};

            require_petsc(
                PetscSFBcastBegin(
                    sf, MPIU_INT,
                    root_values.data(),
                    leaf_values.data(),
                    MPI_REPLACE),
                "PetscSFBcastBegin");
            require_petsc(
                PetscSFBcastEnd(
                    sf, MPIU_INT,
                    root_values.data(),
                    leaf_values.data(),
                    MPI_REPLACE),
                "PetscSFBcastEnd");

            const PetscInt expected =
                mpi_rank == 0 ? 2000 : 1000;
            require(leaf_values[0] == -1,
                    "non-leaf local slot must remain untouched");
            require(leaf_values[1] == expected,
                    "PetscSF scalar broadcast owner-to-ghost");
        }

        require_petsc(PetscSFDestroy(&sf), "PetscSFDestroy");
        require(sf == nullptr, "PetscSFDestroy must clear handle");
    }
}


void verify_global_section_and_section_sf(
    const mesh::DofLayout& layout,
    const mesh::DofNumberingSnapshot& numbering,
    const mesh::PartitionSnapshot& partition,
    const mesh::SharedEntityPlan& plan) {
    PetscSection local_section = nullptr;
    std::vector<PetscInt> local_to_global;
    require_petsc(
        mesh_petsc::create_section_mapping(
            PETSC_COMM_WORLD, layout, numbering,
            &local_section, &local_to_global),
        "create_section_mapping for global section");

    PetscSF point_sf = nullptr;
    require_petsc(
        mesh_petsc::create_point_sf(
            PETSC_COMM_WORLD, layout, partition, plan, &point_sf),
        "create_point_sf");

    PetscInt point_roots = -1;
    PetscInt point_leaves = -1;
    require_petsc(
        PetscSFGetGraph(
            point_sf, &point_roots, &point_leaves, nullptr, nullptr),
        "PetscSFGetGraph point SF");
    require(point_roots == 6, "flattened point SF root count");
    require(point_leaves == 3, "flattened point SF ghost leaf count");

    PetscSection global_section = nullptr;
    require_petsc(
        PetscSectionCreateGlobalSection(
            local_section,
            point_sf,
            PETSC_FALSE,
            PETSC_FALSE,
            &global_section),
        "PetscSectionCreateGlobalSection");

    PetscInt local_storage = 0;
    PetscInt owned_storage = 0;
    require_petsc(
        PetscSectionGetStorageSize(local_section, &local_storage),
        "local PetscSection storage size");
    require_petsc(
        PetscSectionGetConstrainedStorageSize(
            global_section, &owned_storage),
        "global PetscSection owned storage size");
    require(local_storage == 10, "local section storage must contain all local DoFs");
    require(owned_storage == 5, "each rank must own five scalar DoFs");

    int mpi_rank = -1;
    require(
        MPI_Comm_rank(PETSC_COMM_WORLD, &mpi_rank) == MPI_SUCCESS,
        "MPI_Comm_rank for global section");

    PetscInt rank_global_begin = 0;
    require(
        MPI_Exscan(
            &owned_storage,
            &rank_global_begin,
            1,
            MPIU_INT,
            MPI_SUM,
            PETSC_COMM_WORLD) == MPI_SUCCESS,
        "MPI_Exscan global PETSc DoF range");
    if (mpi_rank == 0) rank_global_begin = 0;

    PetscInt global_storage = 0;
    require(
        MPI_Allreduce(
            &owned_storage,
            &global_storage,
            1,
            MPIU_INT,
            MPI_SUM,
            PETSC_COMM_WORLD) == MPI_SUCCESS,
        "MPI_Allreduce global PETSc DoF size");
    require(global_storage == 10, "global PETSc storage size");

    const std::size_t global_storage_size =
        static_cast<std::size_t>(global_storage);
    std::vector<std::uint64_t> petsc_global_to_core(
        global_storage_size,
        std::numeric_limits<std::uint64_t>::max());
    std::vector<int> owner_coverage(global_storage_size, 0);

    const std::size_t cell_count =
        layout.entity_count(mesh::EntityKind::cell);
    const std::size_t face_count =
        layout.entity_count(mesh::EntityKind::face);

    for (PetscInt point = 0; point < 6; ++point) {
        PetscInt local_dof = 0;
        PetscInt local_offset = -1;
        PetscInt global_dof = 0;
        PetscInt global_offset = 0;
        require_petsc(
            PetscSectionGetDof(local_section, point, &local_dof),
            "local point DoF");
        require_petsc(
            PetscSectionGetOffset(local_section, point, &local_offset),
            "local point offset");
        require_petsc(
            PetscSectionGetDof(global_section, point, &global_dof),
            "global point DoF");
        require_petsc(
            PetscSectionGetOffset(global_section, point, &global_offset),
            "global point offset");

        mesh::EntityKind kind = mesh::EntityKind::cell;
        std::size_t local_entity = 0U;
        const std::size_t point_index =
            static_cast<std::size_t>(point);
        if (point_index < cell_count) {
            kind = mesh::EntityKind::cell;
            local_entity = point_index;
        } else if (point_index < cell_count + face_count) {
            kind = mesh::EntityKind::face;
            local_entity = point_index - cell_count;
        } else {
            kind = mesh::EntityKind::vertex;
            local_entity = point_index - cell_count - face_count;
        }

        const auto local_index = mesh::LocalIndex{
            static_cast<mesh::LocalIndex::value_type>(local_entity)};
        const bool owned = partition.is_owned(kind, local_index);

        if (owned) {
            require(global_dof == local_dof,
                    "owned point global DoF sign/width");
            require(global_offset >= 0,
                    "owned point global offset must be nonnegative");
            for (PetscInt d = 0; d < local_dof; ++d) {
                const PetscInt petsc_global = global_offset + d;
                require(
                    petsc_global >= 0 &&
                        petsc_global < global_storage,
                    "owned PETSc global offset range");
                const PetscInt core_global =
                    local_to_global[
                        static_cast<std::size_t>(local_offset + d)];
                require(core_global >= 0,
                        "core global DoF must fit nonnegative PetscInt");
                petsc_global_to_core[
                    static_cast<std::size_t>(petsc_global)] =
                    static_cast<std::uint64_t>(core_global);
                owner_coverage[
                    static_cast<std::size_t>(petsc_global)] = 1;
            }
        } else {
            require(global_dof == -(local_dof + 1),
                    "ghost point global DoF must use PETSc negative encoding");
            require(global_offset < 0,
                    "ghost point global offset must be negative");
        }
    }

    require(
        global_storage <=
            static_cast<PetscInt>(std::numeric_limits<int>::max()),
        "fixture MPI collective count must fit int");
    const int collective_count =
        static_cast<int>(global_storage);
    require(
        MPI_Allreduce(
            MPI_IN_PLACE,
            petsc_global_to_core.data(),
            collective_count,
            MPI_UINT64_T,
            MPI_MIN,
            PETSC_COMM_WORLD) == MPI_SUCCESS,
        "MPI_Allreduce PETSc-global to core-global map");
    require(
        MPI_Allreduce(
            MPI_IN_PLACE,
            owner_coverage.data(),
            collective_count,
            MPI_INT,
            MPI_SUM,
            PETSC_COMM_WORLD) == MPI_SUCCESS,
        "MPI_Allreduce PETSc global owner coverage");

    for (std::size_t i = 0U;
         i < global_storage_size;
         ++i) {
        require(owner_coverage[i] == 1,
                "each PETSc global DoF must have exactly one owner");
        require(
            petsc_global_to_core[i] !=
                std::numeric_limits<std::uint64_t>::max(),
            "PETSc global DoF must resolve to one core GlobalDofIndex");
    }

    for (PetscInt point = 0; point < 6; ++point) {
        PetscInt local_dof = 0;
        PetscInt local_offset = -1;
        PetscInt global_dof = 0;
        PetscInt global_offset = 0;
        require_petsc(
            PetscSectionGetDof(local_section, point, &local_dof),
            "local DoF for core mapping");
        require_petsc(
            PetscSectionGetOffset(local_section, point, &local_offset),
            "local offset for core mapping");
        require_petsc(
            PetscSectionGetDof(global_section, point, &global_dof),
            "global DoF for core mapping");
        require_petsc(
            PetscSectionGetOffset(global_section, point, &global_offset),
            "global offset for core mapping");

        const PetscInt decoded_dof =
            global_dof < 0 ? -(global_dof + 1) : global_dof;
        const PetscInt owner_global_offset =
            global_offset < 0 ? -(global_offset + 1) : global_offset;
        require(decoded_dof == local_dof,
                "global section DoF width must decode to local width");

        for (PetscInt d = 0; d < local_dof; ++d) {
            const std::size_t petsc_global =
                static_cast<std::size_t>(
                    owner_global_offset + d);
            const PetscInt local_core =
                local_to_global[
                    static_cast<std::size_t>(local_offset + d)];
            require(
                petsc_global_to_core[petsc_global] ==
                    static_cast<std::uint64_t>(local_core),
                "owned/ghost PETSc offset must recover the same core GlobalDofIndex");
        }
    }

    PetscSF section_sf = nullptr;
    require_petsc(
        PetscSFCreate(PETSC_COMM_WORLD, &section_sf),
        "PetscSFCreate section SF");
    require_petsc(
        PetscSFSetGraphSection(
            section_sf, local_section, global_section),
        "PetscSFSetGraphSection");
    require_petsc(
        PetscSFSetUp(section_sf),
        "PetscSFSetUp section SF");

    PetscInt section_roots = -1;
    PetscInt section_leaves = -1;
    require_petsc(
        PetscSFGetGraph(
            section_sf,
            &section_roots,
            &section_leaves,
            nullptr,
            nullptr),
        "PetscSFGetGraph section SF");
    require(section_roots == owned_storage,
            "section SF roots must match local owned global storage");
    require(section_leaves == local_storage,
            "section SF leaves must cover the full local DoF vector");

    std::vector<PetscInt> root_values(
        static_cast<std::size_t>(section_roots), -1);
    for (PetscInt root = 0; root < section_roots; ++root) {
        const PetscInt petsc_global =
            rank_global_begin + root;
        const std::uint64_t core_global =
            petsc_global_to_core[
                static_cast<std::size_t>(petsc_global)];
        require(
            core_global <=
                static_cast<std::uint64_t>(
                    std::numeric_limits<PetscInt>::max() - 1000),
            "fixture encoded core GlobalDofIndex must fit PetscInt");
        root_values[static_cast<std::size_t>(root)] =
            static_cast<PetscInt>(1000U + core_global);
    }

    std::vector<PetscInt> local_values(
        static_cast<std::size_t>(local_storage), -777);
    require_petsc(
        PetscSFBcastBegin(
            section_sf,
            MPIU_INT,
            root_values.data(),
            local_values.data(),
            MPI_REPLACE),
        "PetscSFBcastBegin section SF");
    require_petsc(
        PetscSFBcastEnd(
            section_sf,
            MPIU_INT,
            root_values.data(),
            local_values.data(),
            MPI_REPLACE),
        "PetscSFBcastEnd section SF");

    for (PetscInt local = 0; local < local_storage; ++local) {
        const PetscInt expected =
            static_cast<PetscInt>(
                1000 + local_to_global[
                    static_cast<std::size_t>(local)]);
        require(
            local_values[static_cast<std::size_t>(local)] == expected,
            "section SF must broadcast global-layout values into every local DoF");
    }

    require(
        local_values[0] ==
            static_cast<PetscInt>(1000 + local_to_global[0]) &&
        local_values[1] ==
            static_cast<PetscInt>(1000 + local_to_global[1]) &&
        local_values[2] ==
            static_cast<PetscInt>(1000 + local_to_global[2]),
        "multi-DoF cell point broadcast");

    Vec global_vec = nullptr;
    Vec local_vec = nullptr;
    require_petsc(
        mesh_petsc::create_section_vecs(
            PETSC_COMM_WORLD,
            local_section,
            global_section,
            &global_vec,
            &local_vec),
        "create_section_vecs");

    PetscInt global_local_size = -1;
    PetscInt global_size = -1;
    PetscInt local_local_size = -1;
    PetscInt local_size = -1;
    require_petsc(
        VecGetLocalSize(global_vec, &global_local_size),
        "VecGetLocalSize global");
    require_petsc(
        VecGetSize(global_vec, &global_size),
        "VecGetSize global");
    require_petsc(
        VecGetLocalSize(local_vec, &local_local_size),
        "VecGetLocalSize local");
    require_petsc(
        VecGetSize(local_vec, &local_size),
        "VecGetSize local");
    require(global_local_size == owned_storage,
            "global Vec must store only locally owned DoFs");
    require(global_size == global_storage,
            "global Vec global size");
    require(local_local_size == local_storage &&
                local_size == local_storage,
            "local Vec must store owned plus ghost DoFs");

    PetscBool global_is_mpi = PETSC_FALSE;
    PetscBool local_is_seq = PETSC_FALSE;
    require_petsc(
        PetscObjectTypeCompare(
            reinterpret_cast<PetscObject>(global_vec),
            VECMPI,
            &global_is_mpi),
        "global Vec type");
    require_petsc(
        PetscObjectTypeCompare(
            reinterpret_cast<PetscObject>(local_vec),
            VECSEQ,
            &local_is_seq),
        "local Vec type");
    require(global_is_mpi == PETSC_TRUE,
            "global Vec must use VECMPI");
    require(local_is_seq == PETSC_TRUE,
            "local Vec must use VECSEQ");

    PetscInt ownership_begin = -1;
    PetscInt ownership_end = -1;
    require_petsc(
        VecGetOwnershipRange(
            global_vec, &ownership_begin, &ownership_end),
        "VecGetOwnershipRange global");
    require(ownership_begin == rank_global_begin,
            "global Vec ownership start must match global section");
    require(ownership_end - ownership_begin == owned_storage,
            "global Vec ownership width");

    PetscScalar* global_array = nullptr;
    require_petsc(
        VecGetArray(global_vec, &global_array),
        "VecGetArray global initialization");
    for (PetscInt local_root = 0;
         local_root < owned_storage;
         ++local_root) {
        const PetscInt petsc_global =
            ownership_begin + local_root;
        const std::uint64_t core_global =
            petsc_global_to_core[
                static_cast<std::size_t>(petsc_global)];
        global_array[
            static_cast<std::size_t>(local_root)] =
            static_cast<PetscScalar>(1000U + core_global);
    }
    require_petsc(
        VecRestoreArray(global_vec, &global_array),
        "VecRestoreArray global initialization");

    require_petsc(
        VecSet(local_vec, static_cast<PetscScalar>(-777.0)),
        "VecSet local sentinel");
    require_petsc(
        mesh_petsc::global_to_local(
            section_sf, global_vec, local_vec),
        "global_to_local Vec broadcast");

    const PetscScalar* local_array = nullptr;
    require_petsc(
        VecGetArrayRead(local_vec, &local_array),
        "VecGetArrayRead local broadcast");
    for (PetscInt local = 0;
         local < local_storage;
         ++local) {
        const PetscInt core_global =
            local_to_global[
                static_cast<std::size_t>(local)];
        const PetscReal actual =
            PetscRealPart(
                local_array[
                    static_cast<std::size_t>(local)]);
        const PetscReal expected =
            static_cast<PetscReal>(1000 + core_global);
        require(actual == expected,
                "global Vec to local Vec must preserve core GlobalDofIndex identity");
    }
    require_petsc(
        VecRestoreArrayRead(local_vec, &local_array),
        "VecRestoreArrayRead local broadcast");

    std::vector<int> local_copy_counts(
        static_cast<std::size_t>(global_storage), 0);
    for (const PetscInt core_global : local_to_global) {
        require(core_global >= 0 && core_global < global_storage,
                "fixture core GlobalDofIndex range");
        ++local_copy_counts[
            static_cast<std::size_t>(core_global)];
    }
    require(
        MPI_Allreduce(
            MPI_IN_PLACE,
            local_copy_counts.data(),
            collective_count,
            MPI_INT,
            MPI_SUM,
            PETSC_COMM_WORLD) == MPI_SUCCESS,
        "MPI_Allreduce local copy counts");
    for (const int copies : local_copy_counts) {
        require(copies == 2,
                "fixture must contain exactly one owner and one ghost copy per global DoF");
    }

    PetscScalar* local_contributions = nullptr;
    require_petsc(
        VecGetArray(local_vec, &local_contributions),
        "VecGetArray local contributions");
    for (PetscInt local = 0;
         local < local_storage;
         ++local) {
        const PetscInt core_global =
            local_to_global[
                static_cast<std::size_t>(local)];
        local_contributions[
            static_cast<std::size_t>(local)] =
            static_cast<PetscScalar>(
                (mpi_rank + 1) * 100 + core_global);
    }
    require_petsc(
        VecRestoreArray(
            local_vec, &local_contributions),
        "VecRestoreArray local contributions");

    require_petsc(
        VecSet(global_vec, static_cast<PetscScalar>(10.0)),
        "VecSet global ADD baseline");
    require_petsc(
        mesh_petsc::local_to_global_add(
            section_sf, local_vec, global_vec),
        "local_to_global_add Vec reduction");

    const PetscScalar* assembled_global = nullptr;
    require_petsc(
        VecGetArrayRead(global_vec, &assembled_global),
        "VecGetArrayRead assembled global");
    for (PetscInt local_root = 0;
         local_root < owned_storage;
         ++local_root) {
        const PetscInt petsc_global =
            ownership_begin + local_root;
        const std::uint64_t core_global =
            petsc_global_to_core[
                static_cast<std::size_t>(petsc_global)];
        const PetscReal actual =
            PetscRealPart(
                assembled_global[
                    static_cast<std::size_t>(local_root)]);
        const PetscReal expected =
            static_cast<PetscReal>(
                310U + 2U * core_global);
        require(actual == expected,
                "ADD_VALUES must sum owner and ghost contributions exactly once onto the unique owner");
    }
    require_petsc(
        VecRestoreArrayRead(
            global_vec, &assembled_global),
        "VecRestoreArrayRead assembled global");

    require_petsc(
        VecSet(local_vec, static_cast<PetscScalar>(-999.0)),
        "VecSet local post-assembly sentinel");
    require_petsc(
        mesh_petsc::global_to_local(
            section_sf, global_vec, local_vec),
        "global_to_local assembled Vec broadcast");
    require_petsc(
        VecGetArrayRead(local_vec, &local_array),
        "VecGetArrayRead post-assembly local");
    for (PetscInt local = 0;
         local < local_storage;
         ++local) {
        const PetscInt core_global =
            local_to_global[
                static_cast<std::size_t>(local)];
        const PetscReal actual =
            PetscRealPart(
                local_array[
                    static_cast<std::size_t>(local)]);
        const PetscReal expected =
            static_cast<PetscReal>(
                310 + 2 * core_global);
        require(actual == expected,
                "assembled global Vec must broadcast back to owner and ghost local slots");
    }
    require_petsc(
        VecRestoreArrayRead(local_vec, &local_array),
        "VecRestoreArrayRead post-assembly local");

    require_petsc(
        VecDestroy(&local_vec),
        "VecDestroy local Vec");
    require_petsc(
        VecDestroy(&global_vec),
        "VecDestroy global Vec");

    require_petsc(
        PetscSFDestroy(&section_sf),
        "PetscSFDestroy section SF");
    require_petsc(
        PetscSectionDestroy(&global_section),
        "PetscSectionDestroy global section");
    require_petsc(
        PetscSFDestroy(&point_sf),
        "PetscSFDestroy point SF");
    require_petsc(
        PetscSectionDestroy(&local_section),
        "PetscSectionDestroy local section");
}


struct PlexStrata {
    PetscInt cell_start;
    PetscInt cell_end;
    PetscInt face_start;
    PetscInt face_end;
    PetscInt vertex_start;
    PetscInt vertex_end;
};

PlexStrata plex_strata(DM dm) {
    PlexStrata strata{-1, -1, -1, -1, -1, -1};
    require_petsc(
        DMPlexGetHeightStratum(
            dm, 0, &strata.cell_start, &strata.cell_end),
        "DMPlex distributed cell stratum");
    require_petsc(
        DMPlexGetHeightStratum(
            dm, 1, &strata.face_start, &strata.face_end),
        "DMPlex distributed face stratum");
    require_petsc(
        DMPlexGetDepthStratum(
            dm, 0, &strata.vertex_start, &strata.vertex_end),
        "DMPlex distributed vertex stratum");
    return strata;
}

const mesh_petsc::DMPlexPointIdentity& identity_for_point(
    const std::vector<mesh_petsc::DMPlexPointIdentity>& identities,
    PetscInt point) {
    const auto found = std::find_if(
        identities.begin(),
        identities.end(),
        [point](const auto& identity) {
            return identity.point == point;
        });
    require(found != identities.end(),
            "distributed DMPlex point missing stable identity");
    return *found;
}

mesh::PartitionSnapshot partition_from_dm_point_sf(
    DM dm,
    const std::vector<mesh_petsc::DMPlexPointIdentity>& identities,
    int mpi_rank,
    int mpi_size) {
    const auto strata = plex_strata(dm);

    const auto cell_count =
        static_cast<std::size_t>(
            strata.cell_end - strata.cell_start);
    const auto face_count =
        static_cast<std::size_t>(
            strata.face_end - strata.face_start);
    const auto vertex_count =
        static_cast<std::size_t>(
            strata.vertex_end - strata.vertex_start);

    mesh::Topology::EntityIds ids;
    ids.cells.resize(
        cell_count, mesh::GlobalEntityId{0U});
    ids.faces.resize(
        face_count, mesh::GlobalEntityId{0U});
    ids.vertices.resize(
        vertex_count, mesh::GlobalEntityId{0U});

    mesh::EntityOwnerRanks owners;
    owners.cells.assign(
        cell_count,
        mesh::PartitionRank{
            static_cast<
                mesh::PartitionRank::value_type>(
                    mpi_rank)});
    owners.faces.assign(
        face_count,
        mesh::PartitionRank{
            static_cast<
                mesh::PartitionRank::value_type>(
                    mpi_rank)});
    owners.vertices.assign(
        vertex_count,
        mesh::PartitionRank{
            static_cast<
                mesh::PartitionRank::value_type>(
                    mpi_rank)});

    std::vector<std::uint8_t> cell_seen(
        cell_count, std::uint8_t{0U});
    std::vector<std::uint8_t> face_seen(
        face_count, std::uint8_t{0U});
    std::vector<std::uint8_t> vertex_seen(
        vertex_count, std::uint8_t{0U});

    for (const auto& identity : identities) {
        const std::size_t local =
            static_cast<std::size_t>(
                identity.local.value());
        switch (identity.kind) {
        case mesh::EntityKind::cell:
            require(local < ids.cells.size(),
                    "distributed cell identity local index");
            require(cell_seen[local] == 0U,
                    "duplicate distributed cell identity");
            cell_seen[local] = std::uint8_t{1U};
            ids.cells[local] = identity.global;
            break;
        case mesh::EntityKind::face:
            require(local < ids.faces.size(),
                    "distributed face identity local index");
            require(face_seen[local] == 0U,
                    "duplicate distributed face identity");
            face_seen[local] = std::uint8_t{1U};
            ids.faces[local] = identity.global;
            break;
        case mesh::EntityKind::vertex:
            require(local < ids.vertices.size(),
                    "distributed vertex identity local index");
            require(vertex_seen[local] == 0U,
                    "duplicate distributed vertex identity");
            vertex_seen[local] = std::uint8_t{1U};
            ids.vertices[local] = identity.global;
            break;
        case mesh::EntityKind::edge:
            throw std::runtime_error(
                "distributed DMPlex unexpectedly contains core edge identity");
        }
    }

    require(
        std::find(
            cell_seen.begin(), cell_seen.end(),
            std::uint8_t{0U}) == cell_seen.end(),
        "every distributed cell needs stable identity");
    require(
        std::find(
            face_seen.begin(), face_seen.end(),
            std::uint8_t{0U}) == face_seen.end(),
        "every distributed face needs stable identity");
    require(
        std::find(
            vertex_seen.begin(), vertex_seen.end(),
            std::uint8_t{0U}) == vertex_seen.end(),
        "every distributed vertex needs stable identity");

    PetscSF point_sf = nullptr;
    require_petsc(
        DMGetPointSF(dm, &point_sf),
        "DMGetPointSF distributed ownership");
    require(point_sf != nullptr,
            "distributed DMPlex must expose point SF");

    PetscInt nroots = -1;
    PetscInt nleaves = -1;
    const PetscInt* ilocal = nullptr;
    const PetscSFNode* remote = nullptr;
    require_petsc(
        PetscSFGetGraph(
            point_sf,
            &nroots,
            &nleaves,
            &ilocal,
            &remote),
        "PetscSFGetGraph distributed ownership");
    require(nroots >= 0 && nleaves >= 0,
            "distributed point SF graph must be set");

    for (PetscInt leaf = 0; leaf < nleaves; ++leaf) {
        const PetscInt point =
            ilocal != nullptr ? ilocal[leaf] : leaf;
        require(remote != nullptr,
                "distributed point SF remote roots");

        const auto& identity =
            identity_for_point(identities, point);
        require(
            remote[leaf].rank >= 0 &&
                remote[leaf].rank < mpi_size,
            "distributed point SF owner rank range");

        const auto owner =
            mesh::PartitionRank{
                static_cast<
                    mesh::PartitionRank::value_type>(
                        remote[leaf].rank)};
        const std::size_t local =
            static_cast<std::size_t>(
                identity.local.value());

        switch (identity.kind) {
        case mesh::EntityKind::cell:
            owners.cells[local] = owner;
            break;
        case mesh::EntityKind::face:
            owners.faces[local] = owner;
            break;
        case mesh::EntityKind::vertex:
            owners.vertices[local] = owner;
            break;
        case mesh::EntityKind::edge:
            throw std::runtime_error(
                "edge identity cannot be point-SF ghost");
        }
    }

    const mesh::Topology identity_topology{
        std::move(ids), {}};
    return mesh::PartitionSnapshot::create(
        identity_topology,
        mesh::PartitionRank{
            static_cast<
                mesh::PartitionRank::value_type>(
                    mpi_rank)},
        static_cast<std::uint32_t>(mpi_size),
        std::move(owners));
}

void require_expected_identity_owner_counts(
    const mesh::PartitionSnapshot& partition) {
    std::array<int, 2> local_cells{0, 0};
    std::array<int, 7> local_faces{0, 0, 0, 0, 0, 0, 0};
    std::array<int, 6> local_vertices{0, 0, 0, 0, 0, 0};

    const auto record = [&](mesh::EntityKind kind,
                            std::uint64_t base,
                            auto& counts) {
        for (std::size_t local = 0U;
             local < partition.entity_count(kind);
             ++local) {
            const auto index = mesh::LocalIndex{
                static_cast<
                    mesh::LocalIndex::value_type>(
                        local)};
            const auto id =
                partition.global_id(kind, index);
            require(id.value() >= base,
                    "distributed stable ID lower bound");
            const std::uint64_t ordinal =
                id.value() - base;
            require(
                ordinal <
                    static_cast<std::uint64_t>(
                        counts.size()),
                "distributed stable ID range");
            if (partition.is_owned(kind, index)) {
                ++counts[
                    static_cast<std::size_t>(
                        ordinal)];
            }
        }
    };

    record(
        mesh::EntityKind::cell,
        7000000000ULL,
        local_cells);
    record(
        mesh::EntityKind::face,
        6000000000ULL,
        local_faces);
    record(
        mesh::EntityKind::vertex,
        5000000000ULL,
        local_vertices);

    require(
        MPI_Allreduce(
            MPI_IN_PLACE,
            local_cells.data(),
            static_cast<int>(local_cells.size()),
            MPI_INT,
            MPI_SUM,
            PETSC_COMM_WORLD) == MPI_SUCCESS,
        "MPI_Allreduce cell stable owners");
    require(
        MPI_Allreduce(
            MPI_IN_PLACE,
            local_faces.data(),
            static_cast<int>(local_faces.size()),
            MPI_INT,
            MPI_SUM,
            PETSC_COMM_WORLD) == MPI_SUCCESS,
        "MPI_Allreduce face stable owners");
    require(
        MPI_Allreduce(
            MPI_IN_PLACE,
            local_vertices.data(),
            static_cast<int>(local_vertices.size()),
            MPI_INT,
            MPI_SUM,
            PETSC_COMM_WORLD) == MPI_SUCCESS,
        "MPI_Allreduce vertex stable owners");

    for (const int count : local_cells) {
        require(count == 1,
                "each stable cell ID must have exactly one owner");
    }
    for (const int count : local_faces) {
        require(count == 1,
                "each stable face ID must have exactly one owner");
    }
    for (const int count : local_vertices) {
        require(count == 1,
                "each stable vertex ID must have exactly one owner");
    }
}

std::vector<mesh::SharedEntityLink>
shared_links_from_dm_point_sf(
    DM dm,
    const std::vector<mesh_petsc::DMPlexPointIdentity>& identities,
    int mpi_rank,
    int mpi_size) {
    const auto strata = plex_strata(dm);
    const std::array<PetscInt, 6> local_ranges{
        strata.cell_start,
        strata.cell_end,
        strata.face_start,
        strata.face_end,
        strata.vertex_start,
        strata.vertex_end};

    std::vector<PetscInt> all_ranges(
        static_cast<std::size_t>(mpi_size) *
        local_ranges.size());
    require(
        MPI_Allgather(
            local_ranges.data(),
            static_cast<int>(local_ranges.size()),
            MPIU_INT,
            all_ranges.data(),
            static_cast<int>(local_ranges.size()),
            MPIU_INT,
            PETSC_COMM_WORLD) == MPI_SUCCESS,
        "MPI_Allgather DMPlex strata ranges");

    PetscSF point_sf = nullptr;
    require_petsc(
        DMGetPointSF(dm, &point_sf),
        "DMGetPointSF overlap links");

    PetscInt nroots = -1;
    PetscInt nleaves = -1;
    const PetscInt* ilocal = nullptr;
    const PetscSFNode* remote = nullptr;
    require_petsc(
        PetscSFGetGraph(
            point_sf,
            &nroots,
            &nleaves,
            &ilocal,
            &remote),
        "PetscSFGetGraph overlap links");
    require(nroots >= 0 && nleaves >= 0,
            "overlap point SF graph");
    if (nleaves > 0) {
        require(remote != nullptr,
                "overlap point SF remote roots");
    }

    constexpr std::size_t words_per_link = 6U;
    std::vector<std::uint64_t> local_words;
    local_words.reserve(
        static_cast<std::size_t>(nleaves) *
        words_per_link);

    for (PetscInt leaf = 0; leaf < nleaves; ++leaf) {
        const PetscInt point =
            ilocal != nullptr ? ilocal[leaf] : leaf;
        const auto& identity =
            identity_for_point(identities, point);
        const PetscMPIInt owner_rank =
            remote[leaf].rank;
        require(
            owner_rank >= 0 &&
                owner_rank < mpi_size &&
                owner_rank != mpi_rank,
            "overlap leaf must reference remote owner");

        const std::size_t remote_slot =
            static_cast<std::size_t>(owner_rank) *
            local_ranges.size();
        PetscInt remote_start = -1;
        PetscInt remote_end = -1;
        switch (identity.kind) {
        case mesh::EntityKind::cell:
            remote_start =
                all_ranges[remote_slot];
            remote_end =
                all_ranges[remote_slot + 1U];
            break;
        case mesh::EntityKind::face:
            remote_start =
                all_ranges[remote_slot + 2U];
            remote_end =
                all_ranges[remote_slot + 3U];
            break;
        case mesh::EntityKind::vertex:
            remote_start =
                all_ranges[remote_slot + 4U];
            remote_end =
                all_ranges[remote_slot + 5U];
            break;
        case mesh::EntityKind::edge:
            throw std::runtime_error(
                "overlap identity cannot be edge");
        }

        require(
            remote[leaf].index >= remote_start &&
                remote[leaf].index < remote_end,
            "remote root point must lie in matching kind stratum");
        const PetscInt owner_local =
            remote[leaf].index - remote_start;
        require(owner_local >= 0,
                "remote owner local index");

        local_words.push_back(
            static_cast<std::uint64_t>(
                identity.kind));
        local_words.push_back(
            identity.global.value());
        local_words.push_back(
            static_cast<std::uint64_t>(
                owner_rank));
        local_words.push_back(
            static_cast<std::uint64_t>(
                owner_local));
        local_words.push_back(
            static_cast<std::uint64_t>(
                mpi_rank));
        local_words.push_back(
            static_cast<std::uint64_t>(
                identity.local.value()));
    }

    require(
        local_words.size() <=
            static_cast<std::size_t>(
                std::numeric_limits<int>::max()),
        "local shared-link wire size");
    const int local_word_count =
        static_cast<int>(local_words.size());
    std::vector<int> word_counts(
        static_cast<std::size_t>(mpi_size), 0);
    require(
        MPI_Allgather(
            &local_word_count,
            1,
            MPI_INT,
            word_counts.data(),
            1,
            MPI_INT,
            PETSC_COMM_WORLD) == MPI_SUCCESS,
        "MPI_Allgather shared-link word counts");

    std::vector<int> displacements(
        static_cast<std::size_t>(mpi_size), 0);
    int total_words = 0;
    for (int rank = 0; rank < mpi_size; ++rank) {
        require(word_counts[
                    static_cast<std::size_t>(rank)] >= 0,
                "shared-link word count nonnegative");
        displacements[
            static_cast<std::size_t>(rank)] =
            total_words;
        require(
            word_counts[
                static_cast<std::size_t>(rank)] <=
                std::numeric_limits<int>::max() -
                    total_words,
            "shared-link gathered wire size overflow");
        total_words +=
            word_counts[
                static_cast<std::size_t>(rank)];
    }
    require(
        total_words %
            static_cast<int>(words_per_link) == 0,
        "shared-link gathered wire alignment");

    std::vector<std::uint64_t> all_words(
        static_cast<std::size_t>(total_words));
    require(
        MPI_Allgatherv(
            local_words.empty()
                ? nullptr
                : local_words.data(),
            local_word_count,
            MPI_UINT64_T,
            all_words.empty()
                ? nullptr
                : all_words.data(),
            word_counts.data(),
            displacements.data(),
            MPI_UINT64_T,
            PETSC_COMM_WORLD) == MPI_SUCCESS,
        "MPI_Allgatherv canonical shared links");

    std::vector<mesh::SharedEntityLink> links;
    links.reserve(
        all_words.size() / words_per_link);
    for (std::size_t offset = 0U;
         offset < all_words.size();
         offset += words_per_link) {
        const std::uint64_t raw_kind =
            all_words[offset];
        mesh::EntityKind kind;
        switch (raw_kind) {
        case static_cast<std::uint64_t>(
                 mesh::EntityKind::cell):
            kind = mesh::EntityKind::cell;
            break;
        case static_cast<std::uint64_t>(
                 mesh::EntityKind::face):
            kind = mesh::EntityKind::face;
            break;
        case static_cast<std::uint64_t>(
                 mesh::EntityKind::vertex):
            kind = mesh::EntityKind::vertex;
            break;
        default:
            throw std::runtime_error(
                "canonical shared-link kind");
        }

        require(
            all_words[offset + 2U] <
                static_cast<std::uint64_t>(
                    mpi_size) &&
                all_words[offset + 4U] <
                    static_cast<std::uint64_t>(
                        mpi_size),
            "canonical shared-link rank range");
        require(
            all_words[offset + 3U] <=
                static_cast<std::uint64_t>(
                    std::numeric_limits<
                        mesh::LocalIndex::value_type>::max()) &&
                all_words[offset + 5U] <=
                    static_cast<std::uint64_t>(
                        std::numeric_limits<
                            mesh::LocalIndex::value_type>::max()),
            "canonical shared-link local index range");

        links.push_back(
            mesh::SharedEntityLink{
                kind,
                mesh::GlobalEntityId{
                    all_words[offset + 1U]},
                mesh::PartitionRank{
                    static_cast<
                        mesh::PartitionRank::value_type>(
                            all_words[offset + 2U])},
                mesh::LocalIndex{
                    static_cast<
                        mesh::LocalIndex::value_type>(
                            all_words[offset + 3U])},
                mesh::PartitionRank{
                    static_cast<
                        mesh::PartitionRank::value_type>(
                            all_words[offset + 4U])},
                mesh::LocalIndex{
                    static_cast<
                        mesh::LocalIndex::value_type>(
                            all_words[offset + 5U])}});
    }
    return links;
}

mesh::FaceBoundarySnapshot two_by_one_reference_boundary(
    const mesh::Topology& topology) {
    const std::array<mesh::PhysicalTag, 7> tags{
        mesh::PhysicalTag{101U},
        mesh::PhysicalTag{0U},
        mesh::PhysicalTag{102U},
        mesh::PhysicalTag{103U},
        mesh::PhysicalTag{104U},
        mesh::PhysicalTag{105U},
        mesh::PhysicalTag{106U}};
    return mesh::make_face_boundary_snapshot(
        topology, tags);
}

mesh::DenseFieldMetadata migration_metadata(
    std::string id,
    std::string unit,
    std::string locator) {
    return mesh::DenseFieldMetadata{
        std::move(id),
        std::move(unit),
        mesh::FieldSourceMetadata{
            mesh::FieldSourceKind::synthetic_test,
            "mpmc.mesh_petsc.distributed_boundary_property_gate",
            "gate-v1",
            std::move(locator)}};
}

mesh::DenseFieldSnapshot reference_cell_scalar_field(
    const mesh::Topology& topology) {
    std::vector<double> values;
    values.reserve(
        topology.entity_count(
            mesh::EntityKind::cell));
    for (std::size_t i = 0U;
         i < topology.entity_count(
             mesh::EntityKind::cell);
         ++i) {
        values.push_back(
            10.0 +
            static_cast<double>(i) * 0.5);
    }
    return mesh::DenseFieldSnapshot::create(
        topology,
        mesh::EntityKind::cell,
        1U,
        std::move(values),
        migration_metadata(
            "cell.synthetic.scalar",
            "Pa",
            "cell/scalar"));
}

mesh::DenseFieldSnapshot reference_face_vector_field(
    const mesh::Topology& topology) {
    std::vector<double> values;
    values.reserve(
        topology.entity_count(
            mesh::EntityKind::face) *
        2U);
    for (std::size_t i = 0U;
         i < topology.entity_count(
             mesh::EntityKind::face);
         ++i) {
        values.push_back(
            100.0 +
            static_cast<double>(i));
        values.push_back(
            -200.0 -
            static_cast<double>(i) * 2.0);
    }
    return mesh::DenseFieldSnapshot::create(
        topology,
        mesh::EntityKind::face,
        2U,
        std::move(values),
        migration_metadata(
            "face.synthetic.vector",
            "m/s",
            "face/vector2"));
}

mesh::DenseFieldSnapshot reference_vertex_vector_field(
    const mesh::Topology& topology) {
    std::vector<double> values;
    values.reserve(
        topology.entity_count(
            mesh::EntityKind::vertex) *
        3U);
    for (std::size_t i = 0U;
         i < topology.entity_count(
             mesh::EntityKind::vertex);
         ++i) {
        values.push_back(
            1000.0 +
            static_cast<double>(i));
        values.push_back(
            2000.0 +
            static_cast<double>(i) * 3.0);
        values.push_back(
            -3000.0 -
            static_cast<double>(i) * 5.0);
    }
    return mesh::DenseFieldSnapshot::create(
        topology,
        mesh::EntityKind::vertex,
        3U,
        std::move(values),
        migration_metadata(
            "vertex.synthetic.vector",
            "kg/mol",
            "vertex/vector3"));
}

std::uint64_t stable_id_base(
    mesh::EntityKind kind) {
    switch (kind) {
    case mesh::EntityKind::vertex:
        return 5000000000ULL;
    case mesh::EntityKind::face:
        return 6000000000ULL;
    case mesh::EntityKind::cell:
        return 7000000000ULL;
    case mesh::EntityKind::edge:
        break;
    }
    throw std::runtime_error(
        "unsupported stable ID kind");
}

void require_metadata_equal(
    const mesh::DenseFieldMetadata& actual,
    const mesh::DenseFieldMetadata& expected) {
    require(actual.id == expected.id,
            "migrated field metadata id");
    require(actual.unit == expected.unit,
            "migrated field metadata unit");
    require(
        actual.source.kind ==
            expected.source.kind,
        "migrated field metadata source kind");
    require(
        actual.source.reference ==
            expected.source.reference,
        "migrated field metadata reference");
    require(
        actual.source.revision ==
            expected.source.revision,
        "migrated field metadata revision");
    require(
        actual.source.locator ==
            expected.source.locator,
        "migrated field metadata locator");
}

void verify_migrated_dense_field(
    const mesh::DenseFieldSnapshot& actual,
    const mesh::DenseFieldSnapshot& reference,
    const std::vector<
        mesh_petsc::DMPlexPointIdentity>& identities) {
    require(
        actual.location() ==
            reference.location(),
        "migrated field location");
    require(
        actual.component_count() ==
            reference.component_count(),
        "migrated field component count");
    require_metadata_equal(
        actual.metadata(),
        reference.metadata());

    std::size_t expected_entity_count = 0U;
    for (const auto& identity : identities) {
        if (identity.kind !=
            actual.location()) {
            continue;
        }
        ++expected_entity_count;

        const std::uint64_t base =
            stable_id_base(identity.kind);
        require(
            identity.global.value() >= base,
            "migrated field stable ID base");
        const std::uint64_t ordinal =
            identity.global.value() - base;
        require(
            ordinal <
                static_cast<std::uint64_t>(
                    reference.entity_count()),
            "migrated field stable ID range");

        const auto source_local =
            mesh::LocalIndex{
                static_cast<
                    mesh::LocalIndex::value_type>(
                        ordinal)};
        for (std::size_t component = 0U;
             component <
                 actual.component_count();
             ++component) {
            require(
                actual.value(
                    identity.local,
                    component) ==
                    reference.value(
                        source_local,
                        component),
                "migrated DenseFieldSnapshot value by stable GlobalEntityId");
        }
    }
    require(
        actual.entity_count() ==
            expected_entity_count,
        "migrated field target entity count");
}

void verify_migrated_boundary(
    const mesh::FaceBoundarySnapshot& actual,
    const mesh::FaceBoundarySnapshot& reference,
    const std::vector<
        mesh_petsc::DMPlexPointIdentity>& identities) {
    std::size_t face_count = 0U;
    for (const auto& identity : identities) {
        if (identity.kind !=
            mesh::EntityKind::face) {
            continue;
        }
        ++face_count;
        require(
            identity.global.value() >=
                6000000000ULL,
            "migrated boundary stable face ID base");
        const std::uint64_t ordinal =
            identity.global.value() -
            6000000000ULL;
        require(
            ordinal <
                static_cast<std::uint64_t>(
                    reference.face_count()),
            "migrated boundary stable face ID range");
        const auto source_face =
            mesh::LocalIndex{
                static_cast<
                    mesh::LocalIndex::value_type>(
                        ordinal)};
        require(
            actual.classification(
                identity.local) ==
                reference.classification(
                    source_face),
            "migrated face classification by stable GlobalEntityId");
        require(
            actual.physical_tag(
                identity.local) ==
                reference.physical_tag(
                    source_face),
            "migrated PhysicalTag by stable GlobalEntityId");
    }
    require(
        actual.face_count() == face_count,
        "migrated FaceBoundarySnapshot face count");
}

mesh::Geometry2D two_by_one_reference_geometry(
    const mesh::Topology& topology) {
    const std::array<double, 3> x{
        0.0, 1.25, 3.75};
    const std::array<double, 2> y{
        -2.0, 2.0};
    return mesh::make_cartesian_geometry_2d(
        topology, x, y);
}

struct DMPlexCoordinateView {
    PetscSection section;
    Vec values;
    const PetscScalar* array;
};

DMPlexCoordinateView get_dmplex_coordinate_view(DM dm) {
    PetscInt coordinate_dim = -1;
    require_petsc(
        DMGetCoordinateDim(dm, &coordinate_dim),
        "DMGetCoordinateDim distributed geometry");
    require(coordinate_dim == 2,
            "distributed DMPlex coordinate dimension");

    PetscSection section = nullptr;
    Vec values = nullptr;
    require_petsc(
        DMGetCoordinateSection(dm, &section),
        "DMGetCoordinateSection distributed geometry");
    require(section != nullptr,
            "distributed DMPlex coordinate section");
    require_petsc(
        DMGetCoordinatesLocal(dm, &values),
        "DMGetCoordinatesLocal distributed geometry");
    require(values != nullptr,
            "distributed DMPlex local coordinate vector");

    const PetscScalar* array = nullptr;
    require_petsc(
        VecGetArrayRead(values, &array),
        "VecGetArrayRead distributed coordinates");
    return DMPlexCoordinateView{
        section, values, array};
}

void restore_dmplex_coordinate_view(
    DMPlexCoordinateView* view) {
    require(view != nullptr,
            "coordinate view pointer");
    require_petsc(
        VecRestoreArrayRead(
            view->values, &view->array),
        "VecRestoreArrayRead distributed coordinates");
}

mesh::Coordinate2D coordinate_for_dmplex_vertex(
    const DMPlexCoordinateView& view,
    PetscInt point) {
    PetscInt dof = -1;
    PetscInt offset = -1;
    require_petsc(
        PetscSectionGetDof(
            view.section, point, &dof),
        "coordinate section vertex dof");
    require_petsc(
        PetscSectionGetOffset(
            view.section, point, &offset),
        "coordinate section vertex offset");
    require(dof == 2 && offset >= 0,
            "DMPlex vertex must carry two coordinate DoFs");
    return mesh::Coordinate2D{
        static_cast<double>(
            PetscRealPart(
                view.array[
                    static_cast<std::size_t>(
                        offset)])),
        static_cast<double>(
            PetscRealPart(
                view.array[
                    static_cast<std::size_t>(
                        offset + 1)]))};
}

void verify_dmplex_geometry_against_core(
    DM dm,
    const std::vector<
        mesh_petsc::DMPlexPointIdentity>& identities,
    const mesh::Geometry2D& reference) {
    constexpr double tolerance = 1.0e-12;

    auto view = get_dmplex_coordinate_view(dm);

    for (const auto& identity : identities) {
        if (identity.kind !=
            mesh::EntityKind::vertex) {
            PetscInt dof = -1;
            require_petsc(
                PetscSectionGetDof(
                    view.section,
                    identity.point,
                    &dof),
                "nonvertex coordinate dof");
            require(dof == 0,
                    "only DMPlex vertices may carry coordinates");
            continue;
        }

        const auto actual =
            coordinate_for_dmplex_vertex(
                view, identity.point);
        require(
            identity.global.value() >=
                5000000000ULL,
            "vertex stable ID geometry base");
        const std::uint64_t ordinal =
            identity.global.value() -
            5000000000ULL;
        require(
            ordinal <
                static_cast<std::uint64_t>(
                    reference.vertex_count()),
            "vertex stable ID geometry range");
        const auto expected =
            reference.vertex_coordinate_m(
                mesh::LocalIndex{
                    static_cast<
                        mesh::LocalIndex::value_type>(
                            ordinal)});
        require(
            std::abs(actual.x_m - expected.x_m) <=
                    tolerance &&
                std::abs(actual.y_m - expected.y_m) <=
                    tolerance,
            "distributed vertex coordinates must match Geometry2D by stable GlobalEntityId");
    }

    for (const auto& identity : identities) {
        if (identity.kind ==
            mesh::EntityKind::face) {
            PetscInt cone_size = -1;
            const PetscInt* cone = nullptr;
            require_petsc(
                DMPlexGetConeSize(
                    dm,
                    identity.point,
                    &cone_size),
                "geometry face cone size");
            require_petsc(
                DMPlexGetCone(
                    dm,
                    identity.point,
                    &cone),
                "geometry face cone");
            require(cone_size == 2 && cone != nullptr,
                    "geometry face must have two vertices");

            const auto a =
                coordinate_for_dmplex_vertex(
                    view, cone[0]);
            const auto b =
                coordinate_for_dmplex_vertex(
                    view, cone[1]);
            const double dx = b.x_m - a.x_m;
            const double dy = b.y_m - a.y_m;
            const double actual_length =
                std::hypot(dx, dy);

            require(
                identity.global.value() >=
                    6000000000ULL,
                "face stable ID geometry base");
            const std::uint64_t ordinal =
                identity.global.value() -
                6000000000ULL;
            require(
                ordinal <
                    static_cast<std::uint64_t>(
                        reference.face_count()),
                "face stable ID geometry range");
            const double expected_length =
                reference.face_length_m(
                    mesh::LocalIndex{
                        static_cast<
                            mesh::LocalIndex::value_type>(
                                ordinal)});
            require(
                std::abs(
                    actual_length -
                    expected_length) <=
                    tolerance,
                "distributed face length recomputation must match Geometry2D");
        }
    }

    for (const auto& identity : identities) {
        if (identity.kind !=
            mesh::EntityKind::cell) {
            continue;
        }

        PetscInt closure_size = 0;
        PetscInt* closure = nullptr;
        require_petsc(
            DMPlexGetTransitiveClosure(
                dm,
                identity.point,
                PETSC_TRUE,
                &closure_size,
                &closure),
            "DMPlex cell transitive closure geometry");

        std::array<PetscInt, 4> vertices{
            -1, -1, -1, -1};
        std::size_t vertex_count = 0U;
        for (PetscInt i = 0;
             i < closure_size;
             ++i) {
            const PetscInt point =
                closure[2 * i];
            const auto& closure_identity =
                identity_for_point(
                    identities, point);
            if (closure_identity.kind !=
                mesh::EntityKind::vertex) {
                continue;
            }

            const bool already_present =
                std::find(
                    vertices.begin(),
                    vertices.begin() +
                        static_cast<
                            std::ptrdiff_t>(
                                vertex_count),
                    point) !=
                vertices.begin() +
                    static_cast<
                        std::ptrdiff_t>(
                            vertex_count);
            if (!already_present) {
                require(
                    vertex_count <
                        vertices.size(),
                    "quad closure vertex overflow");
                vertices[
                    vertex_count++] = point;
            }
        }
        require_petsc(
            DMPlexRestoreTransitiveClosure(
                dm,
                identity.point,
                PETSC_TRUE,
                &closure_size,
                &closure),
            "DMPlexRestoreTransitiveClosure geometry");
        require(vertex_count == 4U,
                "quad cell closure must contain four unique vertices");

        std::array<mesh::Coordinate2D, 4>
            coordinates{};
        double centroid_x = 0.0;
        double centroid_y = 0.0;
        for (std::size_t i = 0U;
             i < coordinates.size();
             ++i) {
            coordinates[i] =
                coordinate_for_dmplex_vertex(
                    view, vertices[i]);
            centroid_x += coordinates[i].x_m;
            centroid_y += coordinates[i].y_m;
        }
        centroid_x /= 4.0;
        centroid_y /= 4.0;

        std::sort(
            coordinates.begin(),
            coordinates.end(),
            [centroid_x, centroid_y](
                const mesh::Coordinate2D& left,
                const mesh::Coordinate2D& right) {
                return std::atan2(
                           left.y_m - centroid_y,
                           left.x_m - centroid_x) <
                       std::atan2(
                           right.y_m - centroid_y,
                           right.x_m - centroid_x);
            });

        double twice_area = 0.0;
        for (std::size_t i = 0U;
             i < coordinates.size();
             ++i) {
            const auto& a = coordinates[i];
            const auto& b =
                coordinates[
                    (i + 1U) %
                    coordinates.size()];
            twice_area +=
                a.x_m * b.y_m -
                b.x_m * a.y_m;
        }
        const double actual_area =
            0.5 * std::abs(twice_area);

        require(
            identity.global.value() >=
                7000000000ULL,
            "cell stable ID geometry base");
        const std::uint64_t ordinal =
            identity.global.value() -
            7000000000ULL;
        require(
            ordinal <
                static_cast<std::uint64_t>(
                    reference.cell_count()),
            "cell stable ID geometry range");
        const auto cell_local =
            mesh::LocalIndex{
                static_cast<
                    mesh::LocalIndex::value_type>(
                        ordinal)};
        const auto expected_centroid =
            reference.cell_centroid_m(cell_local);
        const double expected_area =
            reference.cell_area_m2(cell_local);

        require(
            std::abs(
                centroid_x -
                expected_centroid.x_m) <=
                    tolerance &&
                std::abs(
                    centroid_y -
                    expected_centroid.y_m) <=
                    tolerance,
            "distributed cell centroid recomputation must match Geometry2D");
        require(
            std::abs(
                actual_area -
                expected_area) <=
                    tolerance,
            "distributed cell area recomputation must match Geometry2D");
    }

    restore_dmplex_coordinate_view(&view);
}

void verify_dmplex_distribute_overlap_identity() {
    int mpi_rank = -1;
    int mpi_size = -1;
    require(
        MPI_Comm_rank(
            PETSC_COMM_WORLD, &mpi_rank) == MPI_SUCCESS,
        "MPI_Comm_rank DMPlex distribute");
    require(
        MPI_Comm_size(
            PETSC_COMM_WORLD, &mpi_size) == MPI_SUCCESS,
        "MPI_Comm_size DMPlex distribute");
    require(mpi_size == 2,
            "DMPlex distribute gate requires exactly two ranks");

    const auto root_topology =
        two_by_one_cartesian_with_stable_ids();
    const auto reference_geometry =
        two_by_one_reference_geometry(
            root_topology);
    const auto reference_boundary =
        two_by_one_reference_boundary(
            root_topology);
    const auto reference_cell_field =
        reference_cell_scalar_field(
            root_topology);
    const auto reference_face_field =
        reference_face_vector_field(
            root_topology);
    const auto reference_vertex_field =
        reference_vertex_vector_field(
            root_topology);

    DM source_dm = nullptr;
    std::vector<mesh_petsc::DMPlexPointIdentity>
        source_identities;
    require_petsc(
        mesh_petsc::create_root_dmplex_topology(
            PETSC_COMM_WORLD,
            0,
            mpi_rank == 0 ? &root_topology : nullptr,
            &source_dm,
            &source_identities),
        "create_root_dmplex_topology");
    require(source_dm != nullptr,
            "rooted DMPlex source");

    require_petsc(
        mesh_petsc::attach_root_geometry2d_coordinates(
            source_dm,
            0,
            mpi_rank == 0
                ? &reference_geometry
                : nullptr,
            source_identities),
        "attach Geometry2D coordinates to rooted DMPlex");
    verify_dmplex_geometry_against_core(
        source_dm,
        source_identities,
        reference_geometry);

    PetscInt source_start = -1;
    PetscInt source_end = -1;
    require_petsc(
        DMPlexGetChart(
            source_dm,
            &source_start,
            &source_end),
        "rooted DMPlex source chart");
    if (mpi_rank == 0) {
        require(
            source_start == 0 &&
                source_end == 15 &&
                source_identities.size() == 15U,
            "rank0 must own complete serial source DAG");
    } else {
        require(
            source_start == 0 &&
                source_end == 0 &&
                source_identities.empty(),
            "non-root rank must start with empty source DAG");
    }

    PetscPartitioner partitioner = nullptr;
    require_petsc(
        DMPlexGetPartitioner(
            source_dm, &partitioner),
        "DMPlexGetPartitioner");
    require_petsc(
        PetscPartitionerSetType(
            partitioner,
            PETSCPARTITIONERSIMPLE),
        "PetscPartitionerSetType simple");

    PetscSF migration_sf = nullptr;
    DM distributed_dm = nullptr;
    require_petsc(
        DMPlexDistribute(
            source_dm,
            0,
            &migration_sf,
            &distributed_dm),
        "DMPlexDistribute overlap0");
    require(
        distributed_dm != nullptr &&
            migration_sf != nullptr,
        "DMPlexDistribute must produce two-rank mesh and migration SF");

    std::vector<mesh_petsc::DMPlexPointIdentity>
        distributed_identities;
    require_petsc(
        mesh_petsc::migrate_dmplex_identities(
            source_dm,
            migration_sf,
            source_identities,
            distributed_dm,
            &distributed_identities),
        "migrate DMPlex identities after distribute");
    verify_dmplex_geometry_against_core(
        distributed_dm,
        distributed_identities,
        reference_geometry);

    std::optional<mesh::FaceBoundarySnapshot>
        distributed_boundary;
    std::optional<mesh::DenseFieldSnapshot>
        distributed_cell_field;
    std::optional<mesh::DenseFieldSnapshot>
        distributed_face_field;
    std::optional<mesh::DenseFieldSnapshot>
        distributed_vertex_field;

    require_petsc(
        mesh_petsc::migrate_face_boundary_snapshot(
            source_dm,
            migration_sf,
            reference_boundary,
            source_identities,
            distributed_dm,
            distributed_identities,
            &distributed_boundary),
        "migrate face boundary after distribute");
    require_petsc(
        mesh_petsc::migrate_dense_field_snapshot(
            source_dm,
            migration_sf,
            reference_cell_field,
            source_identities,
            distributed_dm,
            distributed_identities,
            &distributed_cell_field),
        "migrate cell field after distribute");
    require_petsc(
        mesh_petsc::migrate_dense_field_snapshot(
            source_dm,
            migration_sf,
            reference_face_field,
            source_identities,
            distributed_dm,
            distributed_identities,
            &distributed_face_field),
        "migrate face field after distribute");
    require_petsc(
        mesh_petsc::migrate_dense_field_snapshot(
            source_dm,
            migration_sf,
            reference_vertex_field,
            source_identities,
            distributed_dm,
            distributed_identities,
            &distributed_vertex_field),
        "migrate vertex field after distribute");

    require(
        distributed_boundary.has_value() &&
            distributed_cell_field.has_value() &&
            distributed_face_field.has_value() &&
            distributed_vertex_field.has_value(),
        "distributed boundary and field snapshots must be reconstructed");
    verify_migrated_boundary(
        *distributed_boundary,
        reference_boundary,
        distributed_identities);
    verify_migrated_dense_field(
        *distributed_cell_field,
        reference_cell_field,
        distributed_identities);
    verify_migrated_dense_field(
        *distributed_face_field,
        reference_face_field,
        distributed_identities);
    verify_migrated_dense_field(
        *distributed_vertex_field,
        reference_vertex_field,
        distributed_identities);

    require_petsc(
        PetscSFDestroy(&migration_sf),
        "PetscSFDestroy distribution migration SF");
    require_petsc(
        DMDestroy(&source_dm),
        "DMDestroy rooted source DM");

    const auto distributed_strata =
        plex_strata(distributed_dm);
    require(
        distributed_strata.cell_end -
                distributed_strata.cell_start ==
            1,
        "simple partitioner must assign one cell per rank");

    PetscInt distributed_overlap = -1;
    require_petsc(
        DMPlexGetOverlap(
            distributed_dm,
            &distributed_overlap),
        "DMPlexGetOverlap distributed");
    require(distributed_overlap == 0,
            "first distributed mesh must have zero overlap");

    const auto distributed_partition =
        partition_from_dm_point_sf(
            distributed_dm,
            distributed_identities,
            mpi_rank,
            mpi_size);
    require(
        distributed_partition.owned_count(
            mesh::EntityKind::cell) == 1U &&
            distributed_partition.ghost_count(
                mesh::EntityKind::cell) == 0U,
        "overlap0 cell ownership must be unique");

    const std::size_t local_shared_closure_ghosts =
        distributed_partition.ghost_count(
            mesh::EntityKind::face) +
        distributed_partition.ghost_count(
            mesh::EntityKind::vertex);
    std::uint64_t shared_closure_ghosts =
        static_cast<std::uint64_t>(
            local_shared_closure_ghosts);
    require(
        MPI_Allreduce(
            MPI_IN_PLACE,
            &shared_closure_ghosts,
            1,
            MPI_UINT64_T,
            MPI_SUM,
            PETSC_COMM_WORLD) == MPI_SUCCESS,
        "MPI_Allreduce distributed closure ghosts");
    require(shared_closure_ghosts > 0U,
            "distributed mesh must expose shared face/vertex point-SF leaves");

    require_expected_identity_owner_counts(
        distributed_partition);

    PetscSF distributed_point_sf = nullptr;
    require_petsc(
        DMGetPointSF(
            distributed_dm,
            &distributed_point_sf),
        "DMGetPointSF distributed");
    PetscInt distributed_roots = -1;
    PetscInt distributed_leaves = -1;
    require_petsc(
        PetscSFGetGraph(
            distributed_point_sf,
            &distributed_roots,
            &distributed_leaves,
            nullptr,
            nullptr),
        "PetscSFGetGraph distributed point SF");
    PetscInt distributed_chart_start = -1;
    PetscInt distributed_chart_end = -1;
    require_petsc(
        DMPlexGetChart(
            distributed_dm,
            &distributed_chart_start,
            &distributed_chart_end),
        "DMPlexGetChart distributed point SF");
    require(
        distributed_roots == distributed_chart_end,
        "distributed point SF root space must use DMPlex point-index upper bound");
    require(distributed_leaves >= 0,
            "distributed point SF leaf count");

    const PetscInt* distributed_ilocal = nullptr;
    const PetscSFNode* distributed_remote = nullptr;
    require_petsc(
        PetscSFGetGraph(
            distributed_point_sf,
            &distributed_roots,
            &distributed_leaves,
            &distributed_ilocal,
            &distributed_remote),
        "PetscSFGetGraph distributed point SF leaves");
    for (PetscInt leaf = 0;
         leaf < distributed_leaves;
         ++leaf) {
        const PetscInt point =
            distributed_ilocal != nullptr
                ? distributed_ilocal[leaf]
                : leaf;
        require(
            point >= distributed_chart_start &&
                point < distributed_chart_end,
            "distributed point SF leaf must lie inside local DMPlex chart");
        require(distributed_remote != nullptr,
                "distributed point SF remote roots");
    }

    PetscSF overlap_migration_sf = nullptr;
    DM overlap_dm = nullptr;
    require_petsc(
        DMPlexDistributeOverlap(
            distributed_dm,
            1,
            &overlap_migration_sf,
            &overlap_dm),
        "DMPlexDistributeOverlap depth1");
    require(
        overlap_dm != nullptr &&
            overlap_migration_sf != nullptr,
        "DMPlexDistributeOverlap must produce overlap mesh and migration SF");

    std::vector<mesh_petsc::DMPlexPointIdentity>
        overlap_identities;
    require_petsc(
        mesh_petsc::migrate_dmplex_identities(
            distributed_dm,
            overlap_migration_sf,
            distributed_identities,
            overlap_dm,
            &overlap_identities),
        "migrate DMPlex identities into overlap");
    verify_dmplex_geometry_against_core(
        overlap_dm,
        overlap_identities,
        reference_geometry);

    std::optional<mesh::FaceBoundarySnapshot>
        overlap_boundary;
    std::optional<mesh::DenseFieldSnapshot>
        overlap_cell_field;
    std::optional<mesh::DenseFieldSnapshot>
        overlap_face_field;
    std::optional<mesh::DenseFieldSnapshot>
        overlap_vertex_field;

    require_petsc(
        mesh_petsc::migrate_face_boundary_snapshot(
            distributed_dm,
            overlap_migration_sf,
            *distributed_boundary,
            distributed_identities,
            overlap_dm,
            overlap_identities,
            &overlap_boundary),
        "migrate face boundary into overlap");
    require_petsc(
        mesh_petsc::migrate_dense_field_snapshot(
            distributed_dm,
            overlap_migration_sf,
            *distributed_cell_field,
            distributed_identities,
            overlap_dm,
            overlap_identities,
            &overlap_cell_field),
        "migrate cell field into overlap");
    require_petsc(
        mesh_petsc::migrate_dense_field_snapshot(
            distributed_dm,
            overlap_migration_sf,
            *distributed_face_field,
            distributed_identities,
            overlap_dm,
            overlap_identities,
            &overlap_face_field),
        "migrate face field into overlap");
    require_petsc(
        mesh_petsc::migrate_dense_field_snapshot(
            distributed_dm,
            overlap_migration_sf,
            *distributed_vertex_field,
            distributed_identities,
            overlap_dm,
            overlap_identities,
            &overlap_vertex_field),
        "migrate vertex field into overlap");

    require(
        overlap_boundary.has_value() &&
            overlap_cell_field.has_value() &&
            overlap_face_field.has_value() &&
            overlap_vertex_field.has_value(),
        "overlap boundary and field snapshots must be reconstructed");
    verify_migrated_boundary(
        *overlap_boundary,
        reference_boundary,
        overlap_identities);
    verify_migrated_dense_field(
        *overlap_cell_field,
        reference_cell_field,
        overlap_identities);
    verify_migrated_dense_field(
        *overlap_face_field,
        reference_face_field,
        overlap_identities);
    verify_migrated_dense_field(
        *overlap_vertex_field,
        reference_vertex_field,
        overlap_identities);

    require_petsc(
        PetscSFDestroy(&overlap_migration_sf),
        "PetscSFDestroy overlap migration SF");

    PetscInt overlap_depth = -1;
    require_petsc(
        DMPlexGetOverlap(
            overlap_dm,
            &overlap_depth),
        "DMPlexGetOverlap overlap mesh");
    require(overlap_depth == 1,
            "overlap mesh must record depth one");

    const auto overlap_strata =
        plex_strata(overlap_dm);
    require(
        overlap_strata.cell_end -
                overlap_strata.cell_start ==
            2,
        "depth-one overlap must expose both adjacent cells on each rank");
    require(
        overlap_identities.size() == 15U,
        "depth-one overlap of 2x1 mesh must expose full stable identity set");

    const auto overlap_partition =
        partition_from_dm_point_sf(
            overlap_dm,
            overlap_identities,
            mpi_rank,
            mpi_size);
    require(
        overlap_partition.owned_count(
            mesh::EntityKind::cell) == 1U &&
            overlap_partition.ghost_count(
                mesh::EntityKind::cell) == 1U,
        "overlap partition must contain one owned and one ghost cell");
    require_expected_identity_owner_counts(
        overlap_partition);

    const auto all_links =
        shared_links_from_dm_point_sf(
            overlap_dm,
            overlap_identities,
            mpi_rank,
            mpi_size);
    const auto shared_plan =
        mesh::SharedEntityPlan::create(
            overlap_partition,
            all_links);

    const std::size_t expected_receive_count =
        overlap_partition.ghost_count(
            mesh::EntityKind::cell) +
        overlap_partition.ghost_count(
            mesh::EntityKind::face) +
        overlap_partition.ghost_count(
            mesh::EntityKind::vertex);
    require(
        shared_plan.receive_count() ==
            expected_receive_count,
        "SharedEntityPlan receives must equal DMPlex point-SF ghosts");
    require(
        shared_plan.neighbor_count() == 1U,
        "two-rank overlap must have one halo neighbor");

    std::array<std::uint64_t, 2> local_exchange{
        static_cast<std::uint64_t>(
            shared_plan.send_count()),
        static_cast<std::uint64_t>(
            shared_plan.receive_count())};
    std::array<std::uint64_t, 4> all_exchange{
        0U, 0U, 0U, 0U};
    require(
        MPI_Allgather(
            local_exchange.data(),
            2,
            MPI_UINT64_T,
            all_exchange.data(),
            2,
            MPI_UINT64_T,
            PETSC_COMM_WORLD) == MPI_SUCCESS,
        "MPI_Allgather SharedEntityPlan send receive counts");
    require(
        all_exchange[0] == all_exchange[3] &&
            all_exchange[2] == all_exchange[1],
        "SharedEntityPlan send/receive symmetry across two ranks");

    PetscSF overlap_point_sf = nullptr;
    require_petsc(
        DMGetPointSF(
            overlap_dm,
            &overlap_point_sf),
        "DMGetPointSF overlap");
    PetscInt overlap_roots = -1;
    PetscInt overlap_leaves = -1;
    const PetscInt* overlap_ilocal = nullptr;
    const PetscSFNode* overlap_remote = nullptr;
    require_petsc(
        PetscSFGetGraph(
            overlap_point_sf,
            &overlap_roots,
            &overlap_leaves,
            &overlap_ilocal,
            &overlap_remote),
        "PetscSFGetGraph overlap point SF");
    require(
        overlap_leaves ==
            static_cast<PetscInt>(
                expected_receive_count),
        "overlap point SF leaves must match core ghost count");

    for (PetscInt leaf = 0;
         leaf < overlap_leaves;
         ++leaf) {
        const PetscInt point =
            overlap_ilocal != nullptr
                ? overlap_ilocal[leaf]
                : leaf;
        const auto& identity =
            identity_for_point(
                overlap_identities,
                point);
        require(
            overlap_partition.is_ghost(
                identity.kind,
                identity.local),
            "DMPlex point-SF leaf must be core ghost");
        require(
            overlap_partition.owner_rank(
                identity.kind,
                identity.local).value() ==
                static_cast<
                    mesh::PartitionRank::value_type>(
                        overlap_remote[leaf].rank),
            "DMPlex point-SF owner rank must match PartitionSnapshot");
    }

    require_petsc(
        DMDestroy(&overlap_dm),
        "DMDestroy overlap DMPlex");
    require_petsc(
        DMDestroy(&distributed_dm),
        "DMDestroy distributed DMPlex");
}

mesh::LocalIndex reference_local_by_global(
    const mesh::Topology& topology,
    mesh::EntityKind kind,
    mesh::GlobalEntityId global);

void require_field_metadata_same(
    const mesh::DenseFieldMetadata& left,
    const mesh::DenseFieldMetadata& right);

const mesh::DenseFieldSnapshot&
processed_grdecl_cell_field(
    const mesh::ActiveCornerPointGrid& processed,
    std::string_view id) {
    const auto found =
        std::find_if(
            processed.cell_fields.begin(),
            processed.cell_fields.end(),
            [id](const auto& field) {
                return field.metadata().id == id;
            });
    require(
        found != processed.cell_fields.end(),
        "processed GRDECL cell field ID missing");
    require(
        found->location() ==
                mesh::EntityKind::cell &&
            found->component_count() == 1U,
        "processed GRDECL property must remain scalar cell field");
    return *found;
}


mesh::CellCartesianDiagonalPermeability3D
processed_grdecl_diagonal_permeability(
    const mesh::ActiveCornerPointGrid& processed) {
    return mesh::make_cell_cartesian_diagonal_permeability_3d(
        processed.topology,
        processed_grdecl_cell_field(
            processed, "PERMX"),
        processed_grdecl_cell_field(
            processed, "PERMY"),
        processed_grdecl_cell_field(
            processed, "PERMZ"));
}

void verify_stable_gated_tpfa_transport_stage(
    const discretization_petsc::StableFaceGatedTpfaSnapshot3D& actual,
    const discretization::TpfaInternalFaceTransmissibilitySnapshot3D& reference,
    const mesh::Topology& reference_topology,
    const std::vector<
        mesh_petsc::DMPlexPointIdentity>& identities,
    const mesh::PartitionSnapshot& partition) {
    require(
        actual.geometry_policy
                .max_direct_normal_projection_angle_rad ==
            reference.geometry_policy()
                .max_direct_normal_projection_angle_rad &&
        actual.k_policy
                .max_half_face_co_normal_angle_rad ==
            reference.k_policy()
                .max_half_face_co_normal_angle_rad,
        "migrated gated TPFA policies must match source snapshot");

    require(
        actual.dispositions.size() ==
                actual.entry_count() &&
            actual
                    .materialized_face_transmissibilities_m3
                    .size() ==
                actual.entry_count(),
        "migrated gated TPFA transport array sizes");

    std::size_t expected_entries = 0U;
    std::uint64_t local_owned = 0U;
    std::uint64_t local_ghost = 0U;

    for (const auto& identity : identities) {
        if (identity.kind !=
            mesh::EntityKind::face) {
            continue;
        }

        const auto source_local =
            reference_local_by_global(
                reference_topology,
                mesh::EntityKind::face,
                identity.global);
        const bool expected_internal =
            reference.contains_internal_face(
                source_local);
        const auto found =
            std::find(
                actual.face_global_ids.begin(),
                actual.face_global_ids.end(),
                identity.global);

        if (!expected_internal) {
            require(
                found ==
                    actual.face_global_ids.end(),
                "boundary stable face ID must be absent from gated TPFA transport");
            continue;
        }

        ++expected_entries;
        require(
            found !=
                actual.face_global_ids.end(),
            "internal stable face ID must survive gated TPFA migration");
        const std::size_t entry =
            static_cast<std::size_t>(
                std::distance(
                    actual.face_global_ids.begin(),
                    found));
        const auto& expected =
            reference.entry(
                source_local);
        require(
            actual.dispositions[entry] ==
                expected.disposition,
            "gated TPFA disposition must match by stable face GlobalEntityId");

        if (expected.static_transmissibility.has_value()) {
            require(
                actual
                    .materialized_face_transmissibilities_m3[
                        entry]
                    .has_value(),
                "materialized gated TPFA face must retain numeric T_f after migration");
            require(
                *actual
                     .materialized_face_transmissibilities_m3[
                         entry] ==
                    expected
                        .static_transmissibility
                        ->face_transmissibility_m3,
                "materialized gated TPFA T_f must match by stable face GlobalEntityId");
        } else {
            require(
                !actual
                     .materialized_face_transmissibilities_m3[
                         entry]
                     .has_value(),
                "blocked gated TPFA face must remain value-free after migration");
        }

        if (partition.is_owned(
                mesh::EntityKind::face,
                identity.local)) {
            ++local_owned;
        } else {
            require(
                partition.is_ghost(
                    mesh::EntityKind::face,
                    identity.local),
                "non-owned gated TPFA internal face copy must be ghost");
            ++local_ghost;
        }
    }

    require(
        actual.entry_count() ==
            expected_entries,
        "gated TPFA transport target internal-face count");

    std::array<std::uint64_t, 3> global_counts{
        static_cast<std::uint64_t>(
            expected_entries),
        local_owned,
        local_ghost};
    require(
        MPI_Allreduce(
            MPI_IN_PLACE,
            global_counts.data(),
            static_cast<int>(
                global_counts.size()),
            MPI_UINT64_T,
            MPI_SUM,
            PETSC_COMM_WORLD) == MPI_SUCCESS,
        "MPI_Allreduce gated TPFA stable face copy counts");
    require(
        global_counts[0] == 2U &&
            global_counts[1] == 1U &&
            global_counts[2] == 1U,
        "one stable internal face must appear as one owner and one ghost copy");
}


void verify_target_local_gated_tpfa_view_stage(
    const discretization_petsc::TargetLocalGatedTpfaTransmissibilityView3D& view,
    const discretization_petsc::StableFaceGatedTpfaSnapshot3D& transport,
    const discretization::TpfaInternalFaceTransmissibilitySnapshot3D& reference,
    const mesh::Topology& reference_topology,
    const std::vector<
        mesh_petsc::DMPlexPointIdentity>& identities,
    const mesh::PartitionSnapshot& partition) {
    require(
        view.geometry_policy()
                .max_direct_normal_projection_angle_rad ==
            transport.geometry_policy
                .max_direct_normal_projection_angle_rad &&
        view.k_policy()
                .max_half_face_co_normal_angle_rad ==
            transport.k_policy
                .max_half_face_co_normal_angle_rad,
        "target-local gated TPFA view retains transport policies");

    std::size_t target_face_count = 0U;
    std::size_t materialized_count = 0U;
    std::uint64_t local_owned = 0U;
    std::uint64_t local_ghost = 0U;

    for (const auto& identity : identities) {
        if (identity.kind !=
            mesh::EntityKind::face) {
            continue;
        }
        ++target_face_count;

        const auto source_local =
            reference_local_by_global(
                reference_topology,
                mesh::EntityKind::face,
                identity.global);
        const bool expected_internal =
            reference.contains_internal_face(
                source_local);

        if (!expected_internal) {
            require(
                !view.contains_internal_face(
                    identity.local),
                "boundary target face must remain absent from assembly-facing TPFA view");
            bool boundary_entry_rejected = false;
            try {
                (void)view.entry(
                    identity.local);
            } catch (const std::invalid_argument&) {
                boundary_entry_rejected = true;
            }
            require(
                boundary_entry_rejected,
                "boundary target face entry accessor must reject");
            continue;
        }

        require(
            view.contains_internal_face(
                identity.local),
            "transported internal face must be addressable by target LocalIndex");
        const auto& local_entry =
            view.entry(identity.local);
        const auto& expected =
            reference.entry(source_local);

        require(
            local_entry.face ==
                    identity.local &&
                local_entry.global ==
                    identity.global &&
                view.global_id(
                    identity.local) ==
                    identity.global &&
                local_entry.disposition ==
                    expected.disposition &&
                view.disposition(
                    identity.local) ==
                    expected.disposition,
            "target-local gated TPFA entry identity/disposition");

        const auto optional_value =
            view.optional_transmissibility_m3(
                identity.local);
        if (expected.static_transmissibility.has_value()) {
            ++materialized_count;
            require(
                optional_value.has_value(),
                "materialized target-local TPFA view exposes optional T_f");
            const double expected_value =
                expected
                    .static_transmissibility
                    ->face_transmissibility_m3;
            require(
                *optional_value ==
                    expected_value &&
                view.transmissibility_m3(
                    identity.local) ==
                    expected_value,
                "materialized target-local TPFA numeric accessor");
        } else {
            require(
                !optional_value.has_value(),
                "blocked target-local TPFA view exposes no optional T_f");
            bool blocked_value_rejected = false;
            try {
                (void)view.transmissibility_m3(
                    identity.local);
            } catch (const std::invalid_argument&) {
                blocked_value_rejected = true;
            }
            require(
                blocked_value_rejected,
                "blocked target-local TPFA numeric accessor must reject");
        }

        if (partition.is_owned(
                mesh::EntityKind::face,
                identity.local)) {
            ++local_owned;
        } else {
            require(
                partition.is_ghost(
                    mesh::EntityKind::face,
                    identity.local),
                "target-local gated TPFA non-owned face copy must be ghost");
            ++local_ghost;
        }
    }

    require(
        view.target_face_count() ==
                target_face_count &&
            view.internal_face_count() ==
                transport.entry_count() &&
            view.materialized_face_count() ==
                materialized_count &&
            view.blocked_face_count() ==
                view.internal_face_count() -
                    materialized_count,
        "target-local gated TPFA view counts");

    std::array<std::uint64_t, 3> global_counts{
        static_cast<std::uint64_t>(
            view.internal_face_count()),
        local_owned,
        local_ghost};
    require(
        MPI_Allreduce(
            MPI_IN_PLACE,
            global_counts.data(),
            static_cast<int>(
                global_counts.size()),
            MPI_UINT64_T,
            MPI_SUM,
            PETSC_COMM_WORLD) == MPI_SUCCESS,
        "MPI_Allreduce target-local gated TPFA view copies");
    require(
        global_counts[0] == 2U &&
            global_counts[1] == 1U &&
            global_counts[2] == 1U,
        "target-local gated TPFA view must expose one owned and one ghost stable internal-face copy");
}


mesh::GlobalEntityId
processed_grdecl_face_owner_global_id(
    const mesh::ActiveCornerPointGrid& reference,
    mesh::LocalIndex face);

void verify_assembly_ready_internal_connection_table_stage(
    const discretization_petsc::AssemblyReadyInternalConnectionTable3D& table,
    const discretization_petsc::TargetLocalGatedTpfaTransmissibilityView3D& view,
    const mesh::ActiveCornerPointGrid& reference,
    const std::vector<
        mesh_petsc::DMPlexPointIdentity>& identities) {
    std::size_t target_face_count = 0U;
    std::size_t target_cell_count = 0U;
    for (const auto& identity : identities) {
        if (identity.kind == mesh::EntityKind::face) {
            ++target_face_count;
        } else if (
            identity.kind == mesh::EntityKind::cell) {
            ++target_cell_count;
        }
    }
    require(
        table.target_face_count() ==
                target_face_count &&
            table.target_cell_count() ==
                target_cell_count &&
            table.row_count() ==
                view.materialized_face_count(),
        "assembly-ready internal connection table counts");

    const auto& face_cells =
        reference.topology.relation(
            mesh::EntityKind::face,
            mesh::EntityKind::cell);

    for (const auto& gated :
         view.entries()) {
        const auto source_face =
            reference_local_by_global(
                reference.topology,
                mesh::EntityKind::face,
                gated.global);
        const bool materialized =
            gated.disposition ==
            discretization::TpfaInternalFaceTransmissibilityDisposition3D::
                materialized;

        require(
            table.contains_face(
                gated.face) ==
                materialized,
            "only materialized gated faces enter active connection table");

        if (!materialized) {
            bool blocked_row_rejected = false;
            try {
                (void)table.row(
                    gated.face);
            } catch (const std::invalid_argument&) {
                blocked_row_rejected = true;
            }
            require(
                blocked_row_rejected,
                "blocked gated face active-row accessor must reject");
            continue;
        }

        const auto& row =
            table.row(gated.face);
        require(
            row.face ==
                    gated.face &&
                row.face_global ==
                    gated.global &&
                row.transmissibility_m3 ==
                    view.transmissibility_m3(
                        gated.face) &&
                table.transmissibility_m3(
                    gated.face) ==
                    row.transmissibility_m3,
            "assembly-ready row face identity and T_f");

        const auto adjacent =
            face_cells.adjacent(
                source_face);
        require(
            adjacent.size() == 2U,
            "reference active connection must have two adjacent cells");

        const auto expected_owner_global =
            processed_grdecl_face_owner_global_id(
                reference,
                source_face);
        mesh::GlobalEntityId
            expected_neighbour_global{
                0U};
        bool found_neighbour = false;
        for (const auto cell :
             adjacent) {
            const auto global =
                reference.topology.global_id(
                    mesh::EntityKind::cell,
                    cell);
            if (global != expected_owner_global) {
                expected_neighbour_global =
                    global;
                found_neighbour = true;
            }
        }
        require(
            found_neighbour &&
                row.owner_cell_global ==
                    expected_owner_global &&
                row.neighbour_cell_global ==
                    expected_neighbour_global &&
                row.owner_cell !=
                    row.neighbour_cell,
            "assembly-ready canonical owner/neighbour stable cell identities");

        const auto owner_identity =
            std::find_if(
                identities.begin(),
                identities.end(),
                [&row](const auto& identity) {
                    return identity.kind ==
                               mesh::EntityKind::cell &&
                           identity.local ==
                               row.owner_cell;
                });
        const auto neighbour_identity =
            std::find_if(
                identities.begin(),
                identities.end(),
                [&row](const auto& identity) {
                    return identity.kind ==
                               mesh::EntityKind::cell &&
                           identity.local ==
                               row.neighbour_cell;
                });
        require(
            owner_identity != identities.end() &&
                neighbour_identity != identities.end() &&
                owner_identity->global ==
                    row.owner_cell_global &&
                neighbour_identity->global ==
                    row.neighbour_cell_global,
            "assembly-ready row local cell indices resolve stable cell IDs");
    }

    for (const auto& identity : identities) {
        if (identity.kind !=
            mesh::EntityKind::face ||
            view.contains_internal_face(
                identity.local)) {
            continue;
        }
        require(
            !table.contains_face(
                identity.local),
            "boundary/untransported face never enters active connection table");
    }
}


void verify_parallel_owned_connection_schedule_stage(
    const discretization_petsc::ParallelOwnedConnectionSchedule3D& schedule,
    const discretization_petsc::AssemblyReadyInternalConnectionTable3D& table,
    const mesh::PartitionSnapshot& partition,
    bool expect_one_global_materialized_connection) {
    require(
        schedule.local_rank() ==
                partition.local_rank() &&
            schedule.rank_count() ==
                partition.rank_count() &&
            schedule.local_copy_count() ==
                table.row_count(),
        "parallel owned connection schedule rank metadata and local copy count");

    const auto verify_row_against_table =
        [&](const discretization_petsc::AssemblyReadyInternalConnectionRow3D& row,
            bool should_be_owned) {
            require(
                table.contains_face(
                    row.face),
                "scheduled row must come from active connection table");
            const auto& reference =
                table.row(
                    row.face);
            require(
                reference.face_global ==
                        row.face_global &&
                    reference.owner_cell_global ==
                        row.owner_cell_global &&
                    reference.neighbour_cell_global ==
                        row.neighbour_cell_global &&
                    reference.transmissibility_m3 ==
                        row.transmissibility_m3,
                "scheduled row must preserve table stable IDs and T_f");
            require(
                partition.global_id(
                    mesh::EntityKind::face,
                    row.face) ==
                    row.face_global &&
                partition.is_owned(
                    mesh::EntityKind::face,
                    row.face) ==
                    should_be_owned &&
                partition.is_ghost(
                    mesh::EntityKind::face,
                    row.face) !=
                    should_be_owned,
                "scheduled row role must match PartitionSnapshot face ownership");
        };

    for (const auto& row :
         schedule.assembly_rows()) {
        verify_row_against_table(
            row, true);
    }
    for (const auto& row :
         schedule.ghost_rows()) {
        verify_row_against_table(
            row, false);
    }

    std::array<std::uint64_t, 2> global_counts{
        static_cast<std::uint64_t>(
            schedule.authoritative_row_count()),
        static_cast<std::uint64_t>(
            schedule.ghost_row_count())};
    require(
        MPI_Allreduce(
            MPI_IN_PLACE,
            global_counts.data(),
            static_cast<int>(
                global_counts.size()),
            MPI_UINT64_T,
            MPI_SUM,
            PETSC_COMM_WORLD) == MPI_SUCCESS,
        "MPI_Allreduce authoritative/ghost connection schedule counts");

    if (!expect_one_global_materialized_connection) {
        require(
            global_counts[0] == 0U &&
                global_counts[1] == 0U,
            "blocked-only schedule must have no authoritative or ghost active rows");
        return;
    }

    require(
        global_counts[0] == 1U &&
            global_counts[1] == 1U,
        "one stable materialized face must have exactly one authoritative and one ghost row globally");
    require(
        schedule.local_copy_count() == 1U,
        "two-rank fixture exposes exactly one local copy of materialized connection per rank");

    const bool authoritative =
        schedule.authoritative_row_count() == 1U;
    require(
        authoritative !=
            (schedule.ghost_row_count() == 1U),
        "local materialized connection copy must be exactly authoritative or ghost");

    const auto& local_row =
        authoritative
            ? schedule.assembly_rows().front()
            : schedule.ghost_rows().front();
    std::array<std::uint64_t, 4> local_words{
        authoritative ? 0U : 1U,
        local_row.face_global.value(),
        local_row.owner_cell_global.value(),
        local_row.neighbour_cell_global.value()};
    std::array<std::uint64_t, 8> all_words{};
    require(
        MPI_Allgather(
            local_words.data(),
            static_cast<int>(
                local_words.size()),
            MPI_UINT64_T,
            all_words.data(),
            static_cast<int>(
                local_words.size()),
            MPI_UINT64_T,
            PETSC_COMM_WORLD) == MPI_SUCCESS,
        "MPI_Allgather authoritative/ghost connection stable IDs");

    double local_t =
        local_row.transmissibility_m3;
    std::array<double, 2> all_t{};
    require(
        MPI_Allgather(
            &local_t,
            1,
            MPI_DOUBLE,
            all_t.data(),
            1,
            MPI_DOUBLE,
            PETSC_COMM_WORLD) == MPI_SUCCESS,
        "MPI_Allgather authoritative/ghost connection T_f");

    require(
        all_words[0] != all_words[4] &&
            all_words[1] == all_words[5] &&
            all_words[2] == all_words[6] &&
            all_words[3] == all_words[7] &&
            all_t[0] == all_t[1],
        "authoritative and ghost rows must be identical except for ownership role");
}


void verify_cell_pair_sparsity_stencil_stage(
    const discretization_petsc::CellPairSparsityStencilSnapshot3D& snapshot,
    const mesh::PartitionSnapshot& partition,
    bool expect_materialized_pair) {
    require(
        snapshot.local_rank() ==
                partition.local_rank() &&
            snapshot.rank_count() ==
                partition.rank_count() &&
            snapshot.local_cell_count() ==
                partition.entity_count(
                    mesh::EntityKind::cell) &&
            snapshot.owned_cell_count() ==
                partition.owned_count(
                    mesh::EntityKind::cell),
        "cell-pair sparsity snapshot partition metadata");

    require(
        partition.owned_count(
            mesh::EntityKind::cell) == 1U &&
            partition.ghost_count(
                mesh::EntityKind::cell) == 1U,
        "cell-pair sparsity fixture requires one owned and one ghost cell");

    mesh::LocalIndex owned_cell{0U};
    mesh::LocalIndex ghost_cell{0U};
    bool found_owned = false;
    bool found_ghost = false;
    for (std::size_t local = 0U;
         local < partition.entity_count(
             mesh::EntityKind::cell);
         ++local) {
        const auto cell =
            mesh::LocalIndex{
                static_cast<
                    mesh::LocalIndex::value_type>(
                        local)};
        if (partition.is_owned(
                mesh::EntityKind::cell,
                cell)) {
            owned_cell = cell;
            found_owned = true;
        } else {
            require(
                partition.is_ghost(
                    mesh::EntityKind::cell,
                    cell),
                "non-owned sparsity fixture cell must be ghost");
            ghost_cell = cell;
            found_ghost = true;
        }
    }
    require(
        found_owned && found_ghost,
        "cell-pair sparsity fixture resolves owned/ghost cells");

    require(
        snapshot.contains_owned_cell(
            owned_cell) &&
            !snapshot.contains_owned_cell(
                ghost_cell),
        "structural counts exist only for owned cells");

    bool ghost_counts_rejected = false;
    try {
        (void)snapshot.structural_counts(
            ghost_cell);
    } catch (const std::invalid_argument&) {
        ghost_counts_rejected = true;
    }
    require(
        ghost_counts_rejected,
        "ghost cell structural-count accessor must reject");

    const auto& counts =
        snapshot.structural_counts(
            owned_cell);
    require(
        counts.cell ==
                owned_cell &&
            counts.cell_global ==
                partition.global_id(
                    mesh::EntityKind::cell,
                    owned_cell),
        "owned-cell structural count identity");

    if (!expect_materialized_pair) {
        require(
            snapshot.coupling_count() == 0U &&
                counts.diagonal_block_nnz == 1U &&
                counts.off_diagonal_block_nnz == 0U,
            "blocked-only sparsity snapshot contains self diagonal structure only");
        return;
    }

    require(
        snapshot.coupling_count() == 1U &&
            counts.diagonal_block_nnz == 1U &&
            counts.off_diagonal_block_nnz == 1U,
        "materialized two-cell sparsity snapshot has one self diagonal and one remote off-diagonal");

    const auto& pair =
        snapshot.couplings().front();
    require(
        pair.first_cell_global <
                pair.second_cell_global &&
            partition.global_id(
                mesh::EntityKind::cell,
                pair.first_cell) ==
                pair.first_cell_global &&
            partition.global_id(
                mesh::EntityKind::cell,
                pair.second_cell) ==
                pair.second_cell_global,
        "cell-pair coupling uses canonical stable-ID order and target-local indices");

    const bool owned_is_first =
        pair.first_cell == owned_cell;
    const bool owned_is_second =
        pair.second_cell == owned_cell;
    require(
        owned_is_first != owned_is_second &&
            (pair.first_cell == ghost_cell ||
             pair.second_cell == ghost_cell),
        "local coupling connects exactly one owned and one ghost cell");

    std::array<std::uint64_t, 2> local_pair{
        pair.first_cell_global.value(),
        pair.second_cell_global.value()};
    std::array<std::uint64_t, 4> all_pairs{};
    require(
        MPI_Allgather(
            local_pair.data(),
            static_cast<int>(
                local_pair.size()),
            MPI_UINT64_T,
            all_pairs.data(),
            static_cast<int>(
                local_pair.size()),
            MPI_UINT64_T,
            PETSC_COMM_WORLD) == MPI_SUCCESS,
        "MPI_Allgather local cell-pair sparsity graph");
    require(
        all_pairs[0] == all_pairs[2] &&
            all_pairs[1] == all_pairs[3],
        "both cell-row owner ranks reconstruct the same stable cell pair from authoritative rows only");
}


void verify_petsc_mpiaij_symbolic_preallocation_stage(
    const discretization_petsc::PetscMpiAijSymbolicPreallocation3D& bridge,
    const discretization_petsc::CellPairSparsityStencilSnapshot3D& sparsity,
    const mesh::PartitionSnapshot& partition,
    bool expect_remote_coupling) {
    require(
        bridge.local_rank() ==
                partition.local_rank() &&
            bridge.rank_count() ==
                partition.rank_count() &&
            bridge.local_owned_row_count() ==
                static_cast<PetscInt>(
                    partition.owned_count(
                        mesh::EntityKind::cell)) &&
            bridge.global_row_count() == 2 &&
            bridge.global_row_end() -
                    bridge.global_row_start() ==
                1,
        "symbolic MPIAIJ bridge PETSc row ownership metadata");

    require(
        bridge.owned_cells_in_petsc_row_order().size() ==
                1U &&
            bridge.owned_cell_global_ids().size() ==
                1U &&
            bridge.owned_global_rows().size() ==
                1U &&
            bridge.diagonal_nnz().size() ==
                1U &&
            bridge.off_diagonal_nnz().size() ==
                1U,
        "symbolic MPIAIJ bridge array sizes");

    const auto owned_cell =
        bridge.owned_cells_in_petsc_row_order()
            .front();
    const auto owned_global =
        bridge.owned_cell_global_ids()
            .front();
    const PetscInt owned_row =
        bridge.owned_global_rows()
            .front();

    require(
        partition.is_owned(
            mesh::EntityKind::cell,
            owned_cell) &&
            partition.global_id(
                mesh::EntityKind::cell,
                owned_cell) ==
                owned_global &&
            owned_row ==
                bridge.global_row_start() &&
            bridge.global_row(
                owned_cell) ==
                owned_row,
        "symbolic MPIAIJ owned stable cell maps to PETSc local row range");

    const auto& structural =
        sparsity.structural_counts(
            owned_cell);
    require(
        bridge.diagonal_nnz().front() ==
                static_cast<PetscInt>(
                    structural.diagonal_block_nnz) &&
            bridge.off_diagonal_nnz().front() ==
                static_cast<PetscInt>(
                    structural.off_diagonal_block_nnz),
        "symbolic MPIAIJ arrays exactly mirror structural d_nnz/o_nnz counts");

    require(
        bridge.diagonal_nnz().front() == 1 &&
            bridge.off_diagonal_nnz().front() ==
                (expect_remote_coupling ? 1 : 0),
        "symbolic MPIAIJ fixture structural counts");

    mesh::LocalIndex ghost_cell{0U};
    bool found_ghost = false;
    for (std::size_t local = 0U;
         local < partition.entity_count(
             mesh::EntityKind::cell);
         ++local) {
        const auto cell =
            mesh::LocalIndex{
                static_cast<
                    mesh::LocalIndex::value_type>(
                        local)};
        if (partition.is_ghost(
                mesh::EntityKind::cell,
                cell)) {
            ghost_cell = cell;
            found_ghost = true;
            break;
        }
    }
    require(
        found_ghost,
        "symbolic MPIAIJ fixture resolves one ghost cell");

    const PetscInt ghost_row =
        bridge.global_row(
            ghost_cell);
    require(
        ghost_row >= 0 &&
            ghost_row <
                bridge.global_row_count() &&
            (ghost_row <
                 bridge.global_row_start() ||
             ghost_row >=
                 bridge.global_row_end()),
        "ghost stable cell row is outside local PETSc ownership range");

    std::array<std::uint64_t, 2>
        local_owned_record{
            owned_global.value(),
            static_cast<std::uint64_t>(
                owned_row)};
    std::array<std::uint64_t, 4>
        all_owned_records{};
    require(
        MPI_Allgather(
            local_owned_record.data(),
            static_cast<int>(
                local_owned_record.size()),
            MPI_UINT64_T,
            all_owned_records.data(),
            static_cast<int>(
                local_owned_record.size()),
            MPI_UINT64_T,
            PETSC_COMM_WORLD) == MPI_SUCCESS,
        "MPI_Allgather stable cell to PETSc global-row numbering");

    require(
        all_owned_records[0] !=
                all_owned_records[2] &&
            all_owned_records[1] !=
                all_owned_records[3] &&
            ((all_owned_records[1] == 0U &&
              all_owned_records[3] == 1U) ||
             (all_owned_records[1] == 1U &&
              all_owned_records[3] == 0U)),
        "owned stable cell IDs map one-to-one onto PETSc global rows [0,2)");

    const auto ghost_global =
        partition.global_id(
            mesh::EntityKind::cell,
            ghost_cell);
    bool ghost_matches_remote_owner = false;
    for (std::size_t record = 0U;
         record < 2U;
         ++record) {
        const std::size_t base =
            2U * record;
        if (all_owned_records[base] ==
            ghost_global.value()) {
            ghost_matches_remote_owner =
                all_owned_records[
                    base + 1U] ==
                static_cast<std::uint64_t>(
                    ghost_row);
        }
    }
    require(
        ghost_matches_remote_owner,
        "point-SF ghost row resolves to the remote owner's stable-cell PETSc row");
}

void verify_owned_cell_structural_column_pattern_stage(
    const discretization_petsc::OwnedCellStructuralColumnPatternSnapshot3D& snapshot,
    const discretization_petsc::PetscMpiAijSymbolicPreallocation3D& bridge,
    const discretization_petsc::CellPairSparsityStencilSnapshot3D& sparsity,
    const mesh::PartitionSnapshot& partition,
    bool expect_remote_coupling) {
    require(
        snapshot.local_rank() ==
                partition.local_rank() &&
            snapshot.rank_count() ==
                partition.rank_count() &&
            snapshot.local_cell_count() ==
                partition.entity_count(
                    mesh::EntityKind::cell) &&
            snapshot.global_row_start() ==
                bridge.global_row_start() &&
            snapshot.global_row_end() ==
                bridge.global_row_end() &&
            snapshot.global_row_count() ==
                bridge.global_row_count() &&
            snapshot.row_count() ==
                static_cast<std::size_t>(
                    bridge.local_owned_row_count()),
        "owned-cell structural column-pattern snapshot metadata");

    const auto owned_cells =
        bridge.owned_cells_in_petsc_row_order();
    const auto owned_ids =
        bridge.owned_cell_global_ids();
    const auto owned_rows =
        bridge.owned_global_rows();
    const auto d_nnz =
        bridge.diagonal_nnz();
    const auto o_nnz =
        bridge.off_diagonal_nnz();

    require(
        snapshot.rows().size() ==
            owned_cells.size(),
        "column-pattern row count follows symbolic bridge");

    std::size_t expected_diagonal_storage = 0U;
    std::size_t expected_off_diagonal_storage = 0U;

    for (std::size_t index = 0U;
         index < snapshot.rows().size();
         ++index) {
        const auto& pattern =
            snapshot.rows()[index];
        const auto cell =
            owned_cells[index];
        const auto diagonal_columns =
            snapshot.diagonal_global_columns(
                cell);
        const auto off_diagonal_columns =
            snapshot.off_diagonal_global_columns(
                cell);

        require(
            snapshot.contains_owned_cell(
                cell) &&
                &snapshot.row(cell) ==
                    &pattern &&
                pattern.cell ==
                    cell &&
                pattern.cell_global ==
                    owned_ids[index] &&
                pattern.global_row ==
                    owned_rows[index] &&
                pattern.diagonal_column_count ==
                    diagonal_columns.size() &&
                pattern.off_diagonal_column_count ==
                    off_diagonal_columns.size(),
            "column-pattern row preserves owned stable-cell PETSc row identity");

        const auto& structural =
            sparsity.structural_counts(
                cell);
        require(
            diagonal_columns.size() ==
                    structural.diagonal_block_nnz &&
                off_diagonal_columns.size() ==
                    structural.off_diagonal_block_nnz &&
                diagonal_columns.size() ==
                    static_cast<std::size_t>(
                        d_nnz[index]) &&
                off_diagonal_columns.size() ==
                    static_cast<std::size_t>(
                        o_nnz[index]) &&
                diagonal_columns.size() +
                        off_diagonal_columns.size() ==
                    static_cast<std::size_t>(
                        d_nnz[index] +
                        o_nnz[index]),
            "N_columns,row must equal d_nnz + o_nnz with exact block counts");

        require(
            std::is_sorted(
                diagonal_columns.begin(),
                diagonal_columns.end()) &&
                std::adjacent_find(
                    diagonal_columns.begin(),
                    diagonal_columns.end()) ==
                    diagonal_columns.end() &&
                std::binary_search(
                    diagonal_columns.begin(),
                    diagonal_columns.end(),
                    pattern.global_row),
            "diagonal structural columns are sorted, unique, and include self");

        for (const PetscInt column :
             diagonal_columns) {
            require(
                column >=
                        bridge.global_row_start() &&
                    column <
                        bridge.global_row_end(),
                "diagonal structural column is locally owned");
        }

        require(
            std::is_sorted(
                off_diagonal_columns.begin(),
                off_diagonal_columns.end()) &&
                std::adjacent_find(
                    off_diagonal_columns.begin(),
                    off_diagonal_columns.end()) ==
                    off_diagonal_columns.end(),
            "off-diagonal structural columns are sorted and unique");

        for (const PetscInt column :
             off_diagonal_columns) {
            require(
                column >= 0 &&
                    column <
                        bridge.global_row_count() &&
                    (column <
                         bridge.global_row_start() ||
                     column >=
                         bridge.global_row_end()),
                "off-diagonal structural column is remote-owned");
        }

        require(
            diagonal_columns.size() ==
                    1U &&
                diagonal_columns.front() ==
                    pattern.global_row &&
                off_diagonal_columns.size() ==
                    (expect_remote_coupling
                         ? 1U
                         : 0U),
            "two-cell fixture structural column pattern");

        expected_diagonal_storage +=
            diagonal_columns.size();
        expected_off_diagonal_storage +=
            off_diagonal_columns.size();
    }

    require(
        snapshot.diagonal_column_storage().size() ==
                expected_diagonal_storage &&
            snapshot.off_diagonal_column_storage().size() ==
                expected_off_diagonal_storage,
        "finalized column-pattern storage is compact and has no unreferenced entries");

    mesh::LocalIndex ghost_cell{0U};
    bool found_ghost = false;
    for (std::size_t local = 0U;
         local < partition.entity_count(
             mesh::EntityKind::cell);
         ++local) {
        const auto cell =
            mesh::LocalIndex{
                static_cast<
                    mesh::LocalIndex::value_type>(
                        local)};
        if (partition.is_ghost(
                mesh::EntityKind::cell,
                cell)) {
            ghost_cell = cell;
            found_ghost = true;
            break;
        }
    }
    require(
        found_ghost &&
            !snapshot.contains_owned_cell(
                ghost_cell),
        "column-pattern snapshot excludes ghost rows");

    bool ghost_row_rejected = false;
    try {
        (void)snapshot.row(
            ghost_cell);
    } catch (const std::invalid_argument&) {
        ghost_row_rejected = true;
    }
    require(
        ghost_row_rejected,
        "ghost cell column-pattern accessor must reject");
}


void verify_empty_petsc_mpiaij_symbolic_matrix_stage(
    Mat matrix,
    const discretization_petsc::PetscMpiAijSymbolicPreallocation3D& bridge) {
    require(
        matrix != nullptr,
        "empty symbolic MPIAIJ matrix exists");

    PetscBool is_mpiaij = PETSC_FALSE;
    require_petsc(
        PetscObjectTypeCompare(
            reinterpret_cast<PetscObject>(
                matrix),
            MATMPIAIJ,
            &is_mpiaij),
        "empty symbolic matrix type");
    require(
        is_mpiaij == PETSC_TRUE,
        "symbolic matrix must be MATMPIAIJ");

    PetscInt local_rows = -1;
    PetscInt local_columns = -1;
    PetscInt global_rows = -1;
    PetscInt global_columns = -1;
    PetscInt row_start = -1;
    PetscInt row_end = -1;
    PetscInt column_start = -1;
    PetscInt column_end = -1;

    require_petsc(
        MatGetLocalSize(
            matrix,
            &local_rows,
            &local_columns),
        "empty symbolic MatGetLocalSize");
    require_petsc(
        MatGetSize(
            matrix,
            &global_rows,
            &global_columns),
        "empty symbolic MatGetSize");
    require_petsc(
        MatGetOwnershipRange(
            matrix,
            &row_start,
            &row_end),
        "empty symbolic MatGetOwnershipRange");
    require_petsc(
        MatGetOwnershipRangeColumn(
            matrix,
            &column_start,
            &column_end),
        "empty symbolic MatGetOwnershipRangeColumn");

    require(
        local_rows ==
                bridge.local_owned_row_count() &&
            local_columns ==
                bridge.local_owned_row_count() &&
            global_rows ==
                bridge.global_row_count() &&
            global_columns ==
                bridge.global_row_count() &&
            row_start ==
                bridge.global_row_start() &&
            row_end ==
                bridge.global_row_end() &&
            column_start ==
                row_start &&
            column_end ==
                row_end,
        "empty MPIAIJ matrix sizes and ownership exactly match symbolic bridge");

    const auto owned_cells =
        bridge.owned_cells_in_petsc_row_order();
    const auto owned_ids =
        bridge.owned_cell_global_ids();
    const auto owned_rows =
        bridge.owned_global_rows();
    require(
        owned_cells.size() ==
                owned_rows.size() &&
            owned_cells.size() ==
                owned_ids.size() &&
            owned_cells.size() ==
                static_cast<std::size_t>(
                    local_rows),
        "empty symbolic matrix stable-cell row map sizes");

    for (std::size_t index = 0U;
         index < owned_cells.size();
         ++index) {
        const PetscInt expected_row =
            row_start +
            static_cast<PetscInt>(
                index);
        require(
            owned_rows[index] ==
                    expected_row &&
                bridge.global_row(
                    owned_cells[index]) ==
                    expected_row &&
                (index == 0U ||
                 owned_ids[index - 1U] <
                     owned_ids[index]),
            "empty symbolic matrix ownership range matches stable-cell row map");
    }

    MatInfo info{};
    require_petsc(
        MatGetInfo(
            matrix,
            MAT_LOCAL,
            &info),
        "empty symbolic MatGetInfo");
    require(
        info.nz_used == 0.0,
        "empty symbolic MPIAIJ matrix must remain value-free");
}


void verify_processed_grdecl_cell_field_stage(
    const mesh::DenseFieldSnapshot& actual,
    const mesh::DenseFieldSnapshot& reference,
    const mesh::Topology& reference_topology,
    const std::vector<
        mesh_petsc::DMPlexPointIdentity>& identities,
    const mesh::PartitionSnapshot& partition,
    std::size_t expected_owned,
    std::size_t expected_ghost) {
    require(
        actual.location() ==
                mesh::EntityKind::cell &&
            actual.component_count() ==
                reference.component_count(),
        "processed GRDECL migrated cell field layout");
    require_field_metadata_same(
        actual.metadata(),
        reference.metadata());

    std::size_t cell_count = 0U;
    std::size_t owned_count = 0U;
    std::size_t ghost_count = 0U;
    for (const auto& identity : identities) {
        if (identity.kind !=
            mesh::EntityKind::cell) {
            continue;
        }
        ++cell_count;

        const auto source_local =
            reference_local_by_global(
                reference_topology,
                mesh::EntityKind::cell,
                identity.global);
        require(
            actual.value(
                identity.local, 0U) ==
                reference.value(
                    source_local, 0U),
            "processed GRDECL migrated value by stable cell GlobalEntityId");

        if (partition.is_owned(
                mesh::EntityKind::cell,
                identity.local)) {
            ++owned_count;
        } else {
            require(
                partition.is_ghost(
                    mesh::EntityKind::cell,
                    identity.local),
                "processed GRDECL non-owned cell copy must be ghost");
            ++ghost_count;
        }
    }

    require(
        actual.entity_count() == cell_count,
        "processed GRDECL migrated cell field entity count");
    require(
        owned_count == expected_owned &&
            ghost_count == expected_ghost,
        "processed GRDECL migrated field owned/ghost coverage");
}

mesh::GlobalEntityId
processed_grdecl_face_owner_global_id(
    const mesh::ActiveCornerPointGrid& reference,
    mesh::LocalIndex face) {
    const auto owner_local =
        reference.face_geometry.face_owner(face);
    return reference.topology.global_id(
        mesh::EntityKind::cell,
        owner_local);
}

void verify_stable_owner_face_geometry_stage(
    const mesh_petsc::StableOwnerFaceGeometry3D& actual,
    const mesh::ActiveCornerPointGrid& reference,
    const std::vector<
        mesh_petsc::DMPlexPointIdentity>& identities) {
    constexpr double tolerance = 1.0e-12;

    std::size_t face_count = 0U;
    for (const auto& identity : identities) {
        if (identity.kind !=
            mesh::EntityKind::face) {
            continue;
        }
        ++face_count;

        const auto source_local =
            reference_local_by_global(
                reference.topology,
                mesh::EntityKind::face,
                identity.global);
        const std::size_t target_local =
            static_cast<std::size_t>(
                identity.local.value());
        require(
            target_local < actual.face_count(),
            "migrated stable-owner face geometry local range");

        const auto actual_centroid =
            actual.face_centroids_m[target_local];
        const auto expected_centroid =
            reference.face_geometry
                .face_centroid_m(source_local);
        require(
            std::abs(
                actual_centroid.x_m -
                expected_centroid.x_m) <=
                    tolerance &&
                std::abs(
                    actual_centroid.y_m -
                    expected_centroid.y_m) <=
                    tolerance &&
                std::abs(
                    actual_centroid.z_m -
                    expected_centroid.z_m) <=
                    tolerance,
            "migrated face centroid by stable face GlobalEntityId");

        require(
            std::abs(
                actual.face_areas_m2[target_local] -
                reference.face_geometry
                    .face_area_m2(source_local)) <=
                    tolerance,
            "migrated face area by stable face GlobalEntityId");

        require(
            actual.face_owner_global_ids[target_local] ==
                processed_grdecl_face_owner_global_id(
                    reference,
                    source_local),
            "migrated face owner must remain stable cell GlobalEntityId");

        const auto actual_normal =
            actual.face_owner_unit_normals[
                target_local];
        const auto expected_normal =
            reference.face_geometry
                .face_owner_unit_normal(
                    source_local);
        require(
            std::abs(
                actual_normal.x -
                expected_normal.x) <=
                    tolerance &&
                std::abs(
                    actual_normal.y -
                    expected_normal.y) <=
                    tolerance &&
                std::abs(
                    actual_normal.z -
                    expected_normal.z) <=
                    tolerance,
            "migrated owner-relative normal by stable face GlobalEntityId");
    }

    require(
        actual.face_count() == face_count &&
            actual.face_areas_m2.size() == face_count &&
            actual.face_owner_global_ids.size() == face_count &&
            actual.face_owner_unit_normals.size() == face_count,
        "migrated stable-owner face geometry target face count");
}

void verify_materialized_face_geometry_stage(
    const mesh::FaceGeometry3D& actual,
    const mesh_petsc::StableOwnerFaceGeometry3D& stable,
    const mesh::ActiveCornerPointGrid& reference,
    const std::vector<
        mesh_petsc::DMPlexPointIdentity>& identities) {
    constexpr double tolerance = 1.0e-12;
    require(
        actual.face_count() == stable.face_count(),
        "materialized face geometry face count");

    for (const auto& identity : identities) {
        if (identity.kind !=
            mesh::EntityKind::face) {
            continue;
        }
        const auto source_local =
            reference_local_by_global(
                reference.topology,
                mesh::EntityKind::face,
                identity.global);
        const auto target_face =
            identity.local;
        const std::size_t target_local =
            static_cast<std::size_t>(
                target_face.value());

        const auto owner_local =
            actual.face_owner(target_face);
        const auto owner_identity =
            std::find_if(
                identities.begin(),
                identities.end(),
                [owner_local](const auto& candidate) {
                    return candidate.kind ==
                               mesh::EntityKind::cell &&
                           candidate.local ==
                               owner_local;
                });
        require(
            owner_identity != identities.end(),
            "materialized face owner must resolve to target-local cell identity");
        require(
            owner_identity->global ==
                stable.face_owner_global_ids[
                    target_local] &&
                owner_identity->global ==
                    processed_grdecl_face_owner_global_id(
                        reference,
                        source_local),
            "materialized face owner local index must resolve stable owner cell ID");

        const auto actual_centroid =
            actual.face_centroid_m(
                target_face);
        const auto expected_centroid =
            reference.face_geometry
                .face_centroid_m(source_local);
        require(
            std::abs(
                actual_centroid.x_m -
                expected_centroid.x_m) <=
                    tolerance &&
                std::abs(
                    actual_centroid.y_m -
                    expected_centroid.y_m) <=
                    tolerance &&
                std::abs(
                    actual_centroid.z_m -
                    expected_centroid.z_m) <=
                    tolerance,
            "materialized face centroid");

        require(
            std::abs(
                actual.face_area_m2(
                    target_face) -
                reference.face_geometry
                    .face_area_m2(
                        source_local)) <=
                    tolerance,
            "materialized face area");

        const auto actual_normal =
            actual.face_owner_unit_normal(
                target_face);
        const auto expected_normal =
            reference.face_geometry
                .face_owner_unit_normal(
                    source_local);
        require(
            std::abs(
                actual_normal.x -
                expected_normal.x) <=
                    tolerance &&
                std::abs(
                    actual_normal.y -
                expected_normal.y) <=
                    tolerance &&
                std::abs(
                    actual_normal.z -
                    expected_normal.z) <=
                    tolerance,
            "materialized owner-relative face normal");
    }
}

void verify_processed_grdecl_3d_dmplex_distribute_overlap() {
    int mpi_rank = -1;
    int mpi_size = -1;
    require(
        MPI_Comm_rank(
            PETSC_COMM_WORLD,
            &mpi_rank) == MPI_SUCCESS,
        "MPI_Comm_rank processed GRDECL 3D distribute");
    require(
        MPI_Comm_size(
            PETSC_COMM_WORLD,
            &mpi_size) == MPI_SUCCESS,
        "MPI_Comm_size processed GRDECL 3D distribute");
    require(
        mpi_size == 2,
        "processed GRDECL 3D distribute gate requires exactly two ranks");

    const auto field_reference =
        processed_grdecl_two_by_one_all_active();

    std::optional<mesh::ActiveCornerPointGrid>
        root_processed;
    if (mpi_rank == 0) {
        root_processed.emplace(
            field_reference);
    }

    const auto reference_operator =
        mesh::make_cell_face_geometric_operator_3d(
            field_reference.topology,
            field_reference.vertex_coordinates_m,
            field_reference.face_geometry);
    const auto reference_permeability =
        processed_grdecl_diagonal_permeability(
            field_reference);
    const auto strict_geometry_policy =
        discretization::TransmissibilityGeometryAdmissibilityPolicy3D{
            0.0};
    const auto strict_k_policy =
        discretization::KOrthogonalityAdmissibilityPolicy3D{
            0.0};
    const auto reference_gated_materialized =
        discretization::make_admissibility_gated_internal_face_transmissibility_snapshot_3d(
            reference_operator,
            reference_permeability,
            strict_geometry_policy,
            strict_k_policy);
    const auto zero_permeability =
        mesh::CellCartesianDiagonalPermeability3D{
            std::vector<
                mesh::CartesianDiagonalPermeabilityTensor3D>(
                field_reference.cell_count(),
                mesh::CartesianDiagonalPermeabilityTensor3D{
                    0.0, 0.0, 0.0})};
    const auto reference_gated_blocked =
        discretization::make_admissibility_gated_internal_face_transmissibility_snapshot_3d(
            reference_operator,
            zero_permeability,
            strict_geometry_policy,
            strict_k_policy);

    const std::array<std::string_view, 4>
        property_ids{
            "PORO", "PERMX", "PERMY", "PERMZ"};

    std::array<double, 36>
        expected_vertex_coordinates{};
    if (mpi_rank == 0) {
        require(
            root_processed->topology.entity_count(
                mesh::EntityKind::cell) == 2U &&
                root_processed->topology.entity_count(
                    mesh::EntityKind::face) == 11U &&
                root_processed->topology.entity_count(
                    mesh::EntityKind::vertex) == 12U,
            "root processed GRDECL entity counts");

        std::array<std::uint8_t, 12> seen{};
        for (std::size_t local = 0U;
             local < 12U;
             ++local) {
            const auto local_index =
                mesh::LocalIndex{
                    static_cast<
                        mesh::LocalIndex::value_type>(
                            local)};
            const auto id =
                root_processed->topology.global_id(
                    mesh::EntityKind::vertex,
                    local_index);
            require(
                id.value() >= 1U &&
                    id.value() <= 12U,
                "root processed GRDECL vertex ID range");
            const std::size_t slot =
                static_cast<std::size_t>(
                    id.value() - 1U);
            require(
                seen[slot] == std::uint8_t{0U},
                "root processed GRDECL duplicate vertex ID");
            seen[slot] = std::uint8_t{1U};

            const auto coordinate =
                root_processed
                    ->vertex_coordinates_m[local];
            expected_vertex_coordinates[
                3U * slot] = coordinate.x_m;
            expected_vertex_coordinates[
                3U * slot + 1U] = coordinate.y_m;
            expected_vertex_coordinates[
                3U * slot + 2U] = coordinate.z_m;
        }
    }
    require(
        MPI_Bcast(
            expected_vertex_coordinates.data(),
            static_cast<int>(
                expected_vertex_coordinates.size()),
            MPI_DOUBLE,
            0,
            PETSC_COMM_WORLD) == MPI_SUCCESS,
        "MPI_Bcast processed GRDECL coordinate reference");

    DM source_dm = nullptr;
    std::vector<mesh_petsc::DMPlexPointIdentity>
        source_identities;
    require_petsc(
        mesh_petsc::create_root_dmplex_topology(
            PETSC_COMM_WORLD,
            0,
            mpi_rank == 0
                ? &root_processed->topology
                : nullptr,
            &source_dm,
            &source_identities),
        "create rooted processed GRDECL 3D DMPlex");

    std::span<const mesh::Coordinate3D>
        root_coordinates;
    if (mpi_rank == 0) {
        root_coordinates =
            root_processed->vertex_coordinates_m;
    }
    require_petsc(
        mesh_petsc::attach_root_vertex_coordinates_3d(
            source_dm,
            0,
            root_coordinates,
            source_identities),
        "attach rooted processed GRDECL 3D coordinates");
    verify_processed_grdecl_3d_identity_coordinates(
        source_dm,
        source_identities,
        expected_vertex_coordinates);

    PetscInt source_start = -1;
    PetscInt source_end = -1;
    require_petsc(
        DMPlexGetChart(
            source_dm,
            &source_start,
            &source_end),
        "rooted processed GRDECL source chart");
    if (mpi_rank == 0) {
        require(
            source_start == 0 &&
                source_end == 25 &&
                source_identities.size() == 25U,
            "root rank must hold complete processed GRDECL 3D DAG");
    } else {
        require(
            source_start == 0 &&
                source_end == 0 &&
                source_identities.empty(),
            "non-root rank must begin with empty processed GRDECL DAG");
    }

    mesh_petsc::StableOwnerFaceGeometry3D
        source_face_geometry;
    if (mpi_rank == 0) {
        require_petsc(
            mesh_petsc::make_stable_owner_face_geometry_3d(
                root_processed->face_geometry,
                source_identities,
                &source_face_geometry),
            "freeze processed GRDECL face owners as stable cell IDs");
        verify_stable_owner_face_geometry_stage(
            source_face_geometry,
            field_reference,
            source_identities);
    }

    discretization_petsc::StableFaceGatedTpfaSnapshot3D
        source_gated_materialized;
    discretization_petsc::StableFaceGatedTpfaSnapshot3D
        source_gated_blocked;
    require_petsc(
        discretization_petsc::make_root_stable_face_gated_tpfa_snapshot_3d(
            PETSC_COMM_WORLD,
            0,
            mpi_rank == 0
                ? &reference_gated_materialized
                : nullptr,
            source_identities,
            &source_gated_materialized),
        "freeze materialized gated TPFA snapshot by stable face ID");
    require_petsc(
        discretization_petsc::make_root_stable_face_gated_tpfa_snapshot_3d(
            PETSC_COMM_WORLD,
            0,
            mpi_rank == 0
                ? &reference_gated_blocked
                : nullptr,
            source_identities,
            &source_gated_blocked),
        "freeze blocked gated TPFA snapshot by stable face ID");
    require(
        source_gated_materialized.geometry_policy
                .max_direct_normal_projection_angle_rad ==
                0.0 &&
            source_gated_materialized.k_policy
                .max_half_face_co_normal_angle_rad ==
                0.0 &&
            source_gated_blocked.geometry_policy
                .max_direct_normal_projection_angle_rad ==
                0.0 &&
            source_gated_blocked.k_policy
                .max_half_face_co_normal_angle_rad ==
                0.0,
        "root gated TPFA transport retains strict policies");
    if (mpi_rank == 0) {
        require(
            source_gated_materialized.entry_count() == 1U &&
                source_gated_materialized
                    .materialized_face_transmissibilities_m3[
                        0U]
                    .has_value(),
            "root materialized gated TPFA transport");
        require(
            source_gated_blocked.entry_count() == 1U &&
                source_gated_blocked.dispositions[0U] ==
                    discretization::TpfaInternalFaceTransmissibilityDisposition3D::
                        blocked_degenerate_permeability_direction &&
                !source_gated_blocked
                     .materialized_face_transmissibilities_m3[
                         0U]
                     .has_value(),
            "root blocked gated TPFA transport contains no T_f");
    } else {
        require(
            source_gated_materialized.entry_count() == 0U &&
                source_gated_blocked.entry_count() == 0U,
            "non-root rooted gated TPFA transports contain policy only");
    }

    PetscPartitioner partitioner = nullptr;
    require_petsc(
        DMPlexGetPartitioner(
            source_dm,
            &partitioner),
        "processed GRDECL DMPlexGetPartitioner");
    require_petsc(
        PetscPartitionerSetType(
            partitioner,
            PETSCPARTITIONERSIMPLE),
        "processed GRDECL simple partitioner");

    PetscSF migration_sf = nullptr;
    DM distributed_dm = nullptr;
    require_petsc(
        DMPlexDistribute(
            source_dm,
            0,
            &migration_sf,
            &distributed_dm),
        "distribute processed GRDECL 3D overlap0");
    require(
        distributed_dm != nullptr &&
            migration_sf != nullptr,
        "processed GRDECL distribution outputs");

    std::vector<mesh_petsc::DMPlexPointIdentity>
        distributed_identities;
    require_petsc(
        mesh_petsc::migrate_dmplex_identities(
            source_dm,
            migration_sf,
            source_identities,
            distributed_dm,
            &distributed_identities),
        "migrate processed GRDECL identities after distribute");
    verify_processed_grdecl_3d_identity_coordinates(
        distributed_dm,
        distributed_identities,
        expected_vertex_coordinates);

    const auto distributed_partition =
        partition_from_dm_point_sf(
            distributed_dm,
            distributed_identities,
            mpi_rank,
            mpi_size);
    require(
        distributed_partition.owned_count(
            mesh::EntityKind::cell) == 1U &&
            distributed_partition.ghost_count(
                mesh::EntityKind::cell) == 0U,
        "processed GRDECL overlap0 cell ownership");
    require_processed_grdecl_identity_owner_counts(
        distributed_partition);

    discretization_petsc::StableFaceGatedTpfaSnapshot3D
        distributed_gated_materialized;
    discretization_petsc::StableFaceGatedTpfaSnapshot3D
        distributed_gated_blocked;
    require_petsc(
        discretization_petsc::migrate_stable_face_gated_tpfa_snapshot_3d(
            source_dm,
            migration_sf,
            source_gated_materialized,
            source_identities,
            distributed_dm,
            distributed_identities,
            &distributed_gated_materialized),
        "migrate materialized gated TPFA snapshot after distribute");
    require_petsc(
        discretization_petsc::migrate_stable_face_gated_tpfa_snapshot_3d(
            source_dm,
            migration_sf,
            source_gated_blocked,
            source_identities,
            distributed_dm,
            distributed_identities,
            &distributed_gated_blocked),
        "migrate blocked gated TPFA snapshot after distribute");
    verify_stable_gated_tpfa_transport_stage(
        distributed_gated_materialized,
        reference_gated_materialized,
        field_reference.topology,
        distributed_identities,
        distributed_partition);
    verify_stable_gated_tpfa_transport_stage(
        distributed_gated_blocked,
        reference_gated_blocked,
        field_reference.topology,
        distributed_identities,
        distributed_partition);

    std::optional<
        discretization_petsc::TargetLocalGatedTpfaTransmissibilityView3D>
        distributed_materialized_view;
    std::optional<
        discretization_petsc::TargetLocalGatedTpfaTransmissibilityView3D>
        distributed_blocked_view;
    require_petsc(
        discretization_petsc::make_target_local_gated_tpfa_transmissibility_view_3d(
            distributed_dm,
            distributed_gated_materialized,
            distributed_identities,
            &distributed_materialized_view),
        "build distributed target-local materialized TPFA view");
    require_petsc(
        discretization_petsc::make_target_local_gated_tpfa_transmissibility_view_3d(
            distributed_dm,
            distributed_gated_blocked,
            distributed_identities,
            &distributed_blocked_view),
        "build distributed target-local blocked TPFA view");
    require(
        distributed_materialized_view.has_value() &&
            distributed_blocked_view.has_value(),
        "distributed target-local TPFA views materialized");
    verify_target_local_gated_tpfa_view_stage(
        *distributed_materialized_view,
        distributed_gated_materialized,
        reference_gated_materialized,
        field_reference.topology,
        distributed_identities,
        distributed_partition);
    verify_target_local_gated_tpfa_view_stage(
        *distributed_blocked_view,
        distributed_gated_blocked,
        reference_gated_blocked,
        field_reference.topology,
        distributed_identities,
        distributed_partition);

    std::array<
        std::optional<mesh::DenseFieldSnapshot>,
        4>
        distributed_fields;
    for (std::size_t field_index = 0U;
         field_index < property_ids.size();
         ++field_index) {
        const auto& source_field =
            processed_grdecl_cell_field(
                field_reference,
                property_ids[field_index]);
        require_petsc(
            mesh_petsc::migrate_dense_field_snapshot(
                source_dm,
                migration_sf,
                source_field,
                source_identities,
                distributed_dm,
                distributed_identities,
                &distributed_fields[field_index]),
            "migrate processed GRDECL cell field after distribute");
        require(
            distributed_fields[field_index].has_value(),
            "distributed processed GRDECL field must be reconstructed");
        verify_processed_grdecl_cell_field_stage(
            *distributed_fields[field_index],
            source_field,
            field_reference.topology,
            distributed_identities,
            distributed_partition,
            1U,
            0U);
    }

    mesh_petsc::StableOwnerFaceGeometry3D
        distributed_face_geometry;
    require_petsc(
        mesh_petsc::migrate_stable_owner_face_geometry_3d(
            source_dm,
            migration_sf,
            source_face_geometry,
            source_identities,
            distributed_dm,
            distributed_identities,
            &distributed_face_geometry),
        "migrate processed GRDECL face geometry after distribute");
    verify_stable_owner_face_geometry_stage(
        distributed_face_geometry,
        field_reference,
        distributed_identities);

    std::optional<
        discretization_petsc::AssemblyReadyInternalConnectionTable3D>
        distributed_materialized_table;
    const PetscErrorCode
        distributed_materialized_table_error =
            discretization_petsc::make_assembly_ready_internal_connection_table_3d(
                distributed_dm,
                *distributed_materialized_view,
                distributed_face_geometry,
                distributed_identities,
                &distributed_materialized_table);
    require(
        distributed_materialized_table_error ==
                PETSC_ERR_ARG_WRONGSTATE &&
            !distributed_materialized_table.has_value(),
        "overlap0 materialized connection table must reject missing remote cell LocalIndex");

    int local_overlap0_wrongstate =
        distributed_materialized_table_error ==
                PETSC_ERR_ARG_WRONGSTATE
            ? 1
            : 0;
    int overlap0_wrongstate_ranks = 0;
    require(
        MPI_Allreduce(
            &local_overlap0_wrongstate,
            &overlap0_wrongstate_ranks,
            1,
            MPI_INT,
            MPI_SUM,
            PETSC_COMM_WORLD) == MPI_SUCCESS,
        "MPI_Allreduce overlap0 active connection incompleteness");
    require(
        overlap0_wrongstate_ranks == 2,
        "both overlap0 ranks must refuse a cross-rank materialized connection without both local cells");

    std::optional<
        discretization_petsc::AssemblyReadyInternalConnectionTable3D>
        distributed_blocked_table;
    require_petsc(
        discretization_petsc::make_assembly_ready_internal_connection_table_3d(
            distributed_dm,
            *distributed_blocked_view,
            distributed_face_geometry,
            distributed_identities,
            &distributed_blocked_table),
        "build overlap0 blocked-only active connection table");
    require(
        distributed_blocked_table.has_value() &&
            distributed_blocked_table
                    ->row_count() ==
                0U,
        "blocked overlap0 view produces empty active connection table");
    verify_assembly_ready_internal_connection_table_stage(
        *distributed_blocked_table,
        *distributed_blocked_view,
        field_reference,
        distributed_identities);

    std::optional<mesh::FaceGeometry3D>
        distributed_materialized_face_geometry;
    const PetscErrorCode distributed_materialize_error =
        mesh_petsc::materialize_face_geometry_3d(
            distributed_face_geometry,
            distributed_identities,
            &distributed_materialized_face_geometry);
    require(
        distributed_materialize_error == PETSC_SUCCESS ||
            distributed_materialize_error ==
                PETSC_ERR_ARG_WRONGSTATE,
        "overlap0 face geometry materialization may fail only when canonical owner cell is remote");
    const int local_unresolved_owner =
        distributed_materialize_error ==
                PETSC_ERR_ARG_WRONGSTATE
            ? 1
            : 0;
    int unresolved_owner_ranks = 0;
    require(
        MPI_Allreduce(
            &local_unresolved_owner,
            &unresolved_owner_ranks,
            1,
            MPI_INT,
            MPI_SUM,
            PETSC_COMM_WORLD) == MPI_SUCCESS,
        "MPI_Allreduce unresolved stable face owners");
    require(
        unresolved_owner_ranks >= 1,
        "overlap0 must exercise at least one face whose stable canonical owner is remote");

    require_petsc(
        PetscSFDestroy(&migration_sf),
        "destroy processed GRDECL distribution migration SF");
    require_petsc(
        DMDestroy(&source_dm),
        "destroy rooted processed GRDECL source DM");

    PetscInt distributed_overlap = -1;
    require_petsc(
        DMPlexGetOverlap(
            distributed_dm,
            &distributed_overlap),
        "processed GRDECL overlap0 query");
    require(
        distributed_overlap == 0,
        "processed GRDECL first distributed mesh must have overlap zero");

    const auto distributed_strata =
        plex_strata(distributed_dm);
    require(
        distributed_strata.cell_end -
                distributed_strata.cell_start ==
            1,
        "processed GRDECL simple partitioner must assign one hexa per rank");

    PetscSF overlap_migration_sf = nullptr;
    DM overlap_dm = nullptr;
    require_petsc(
        DMPlexDistributeOverlap(
            distributed_dm,
            1,
            &overlap_migration_sf,
            &overlap_dm),
        "processed GRDECL depth-one overlap");
    require(
        overlap_dm != nullptr &&
            overlap_migration_sf != nullptr,
        "processed GRDECL depth-one overlap outputs");

    std::vector<mesh_petsc::DMPlexPointIdentity>
        overlap_identities;
    require_petsc(
        mesh_petsc::migrate_dmplex_identities(
            distributed_dm,
            overlap_migration_sf,
            distributed_identities,
            overlap_dm,
            &overlap_identities),
        "migrate processed GRDECL identities into overlap");
    verify_processed_grdecl_3d_identity_coordinates(
        overlap_dm,
        overlap_identities,
        expected_vertex_coordinates);

    PetscInt overlap_depth = -1;
    require_petsc(
        DMPlexGetOverlap(
            overlap_dm,
            &overlap_depth),
        "processed GRDECL overlap depth query");
    require(
        overlap_depth == 1,
        "processed GRDECL overlap mesh must record depth one");

    const auto overlap_strata =
        plex_strata(overlap_dm);
    require(
        overlap_strata.cell_end -
                overlap_strata.cell_start ==
            2,
        "processed GRDECL depth-one overlap must expose both hexa cells");
    require(
        overlap_identities.size() == 25U,
        "processed GRDECL depth-one overlap must expose full stable identity set");

    const auto overlap_partition =
        partition_from_dm_point_sf(
            overlap_dm,
            overlap_identities,
            mpi_rank,
            mpi_size);
    require(
        overlap_partition.owned_count(
            mesh::EntityKind::cell) == 1U &&
            overlap_partition.ghost_count(
                mesh::EntityKind::cell) == 1U,
        "processed GRDECL overlap must contain one owned and one ghost cell");
    require_processed_grdecl_identity_owner_counts(
        overlap_partition);

    discretization_petsc::StableFaceGatedTpfaSnapshot3D
        overlap_gated_materialized;
    discretization_petsc::StableFaceGatedTpfaSnapshot3D
        overlap_gated_blocked;
    require_petsc(
        discretization_petsc::migrate_stable_face_gated_tpfa_snapshot_3d(
            distributed_dm,
            overlap_migration_sf,
            distributed_gated_materialized,
            distributed_identities,
            overlap_dm,
            overlap_identities,
            &overlap_gated_materialized),
        "migrate materialized gated TPFA snapshot into overlap");
    require_petsc(
        discretization_petsc::migrate_stable_face_gated_tpfa_snapshot_3d(
            distributed_dm,
            overlap_migration_sf,
            distributed_gated_blocked,
            distributed_identities,
            overlap_dm,
            overlap_identities,
            &overlap_gated_blocked),
        "migrate blocked gated TPFA snapshot into overlap");
    verify_stable_gated_tpfa_transport_stage(
        overlap_gated_materialized,
        reference_gated_materialized,
        field_reference.topology,
        overlap_identities,
        overlap_partition);
    verify_stable_gated_tpfa_transport_stage(
        overlap_gated_blocked,
        reference_gated_blocked,
        field_reference.topology,
        overlap_identities,
        overlap_partition);

    std::optional<
        discretization_petsc::TargetLocalGatedTpfaTransmissibilityView3D>
        overlap_materialized_view;
    std::optional<
        discretization_petsc::TargetLocalGatedTpfaTransmissibilityView3D>
        overlap_blocked_view;
    require_petsc(
        discretization_petsc::make_target_local_gated_tpfa_transmissibility_view_3d(
            overlap_dm,
            overlap_gated_materialized,
            overlap_identities,
            &overlap_materialized_view),
        "build overlap target-local materialized TPFA view");
    require_petsc(
        discretization_petsc::make_target_local_gated_tpfa_transmissibility_view_3d(
            overlap_dm,
            overlap_gated_blocked,
            overlap_identities,
            &overlap_blocked_view),
        "build overlap target-local blocked TPFA view");
    require(
        overlap_materialized_view.has_value() &&
            overlap_blocked_view.has_value(),
        "overlap target-local TPFA views materialized");
    verify_target_local_gated_tpfa_view_stage(
        *overlap_materialized_view,
        overlap_gated_materialized,
        reference_gated_materialized,
        field_reference.topology,
        overlap_identities,
        overlap_partition);
    verify_target_local_gated_tpfa_view_stage(
        *overlap_blocked_view,
        overlap_gated_blocked,
        reference_gated_blocked,
        field_reference.topology,
        overlap_identities,
        overlap_partition);

    std::array<
        std::optional<mesh::DenseFieldSnapshot>,
        4>
        overlap_fields;
    for (std::size_t field_index = 0U;
         field_index < property_ids.size();
         ++field_index) {
        require(
            distributed_fields[field_index].has_value(),
            "distributed processed GRDECL source field must exist");
        const auto& reference_field =
            processed_grdecl_cell_field(
                field_reference,
                property_ids[field_index]);
        require_petsc(
            mesh_petsc::migrate_dense_field_snapshot(
                distributed_dm,
                overlap_migration_sf,
                *distributed_fields[field_index],
                distributed_identities,
                overlap_dm,
                overlap_identities,
                &overlap_fields[field_index]),
            "migrate processed GRDECL cell field into overlap");
        require(
            overlap_fields[field_index].has_value(),
            "overlap processed GRDECL field must be reconstructed");
        verify_processed_grdecl_cell_field_stage(
            *overlap_fields[field_index],
            reference_field,
            field_reference.topology,
            overlap_identities,
            overlap_partition,
            1U,
            1U);
    }

    mesh_petsc::StableOwnerFaceGeometry3D
        overlap_face_geometry;
    require_petsc(
        mesh_petsc::migrate_stable_owner_face_geometry_3d(
            distributed_dm,
            overlap_migration_sf,
            distributed_face_geometry,
            distributed_identities,
            overlap_dm,
            overlap_identities,
            &overlap_face_geometry),
        "migrate processed GRDECL face geometry into overlap");
    verify_stable_owner_face_geometry_stage(
        overlap_face_geometry,
        field_reference,
        overlap_identities);

    std::optional<
        discretization_petsc::AssemblyReadyInternalConnectionTable3D>
        overlap_materialized_table;
    std::optional<
        discretization_petsc::AssemblyReadyInternalConnectionTable3D>
        overlap_blocked_table;
    require_petsc(
        discretization_petsc::make_assembly_ready_internal_connection_table_3d(
            overlap_dm,
            *overlap_materialized_view,
            overlap_face_geometry,
            overlap_identities,
            &overlap_materialized_table),
        "build depth-one assembly-ready materialized connection table");
    require_petsc(
        discretization_petsc::make_assembly_ready_internal_connection_table_3d(
            overlap_dm,
            *overlap_blocked_view,
            overlap_face_geometry,
            overlap_identities,
            &overlap_blocked_table),
        "build depth-one blocked-only active connection table");
    require(
        overlap_materialized_table.has_value() &&
            overlap_materialized_table
                    ->row_count() ==
                1U &&
            overlap_blocked_table.has_value() &&
            overlap_blocked_table
                    ->row_count() ==
                0U,
        "depth-one active connection row counts");
    verify_assembly_ready_internal_connection_table_stage(
        *overlap_materialized_table,
        *overlap_materialized_view,
        field_reference,
        overlap_identities);
    verify_assembly_ready_internal_connection_table_stage(
        *overlap_blocked_table,
        *overlap_blocked_view,
        field_reference,
        overlap_identities);

    std::optional<
        discretization_petsc::ParallelOwnedConnectionSchedule3D>
        overlap_materialized_schedule;
    std::optional<
        discretization_petsc::ParallelOwnedConnectionSchedule3D>
        overlap_blocked_schedule;
    require_petsc(
        discretization_petsc::make_parallel_owned_connection_schedule_3d(
            *overlap_materialized_table,
            overlap_partition,
            &overlap_materialized_schedule),
        "build parallel owned materialized connection schedule");
    require_petsc(
        discretization_petsc::make_parallel_owned_connection_schedule_3d(
            *overlap_blocked_table,
            overlap_partition,
            &overlap_blocked_schedule),
        "build parallel owned blocked-only connection schedule");
    require(
        overlap_materialized_schedule.has_value() &&
            overlap_blocked_schedule.has_value(),
        "parallel owned connection schedules constructed");
    verify_parallel_owned_connection_schedule_stage(
        *overlap_materialized_schedule,
        *overlap_materialized_table,
        overlap_partition,
        true);
    verify_parallel_owned_connection_schedule_stage(
        *overlap_blocked_schedule,
        *overlap_blocked_table,
        overlap_partition,
        false);

    std::optional<
        discretization_petsc::CellPairSparsityStencilSnapshot3D>
        materialized_sparsity;
    std::optional<
        discretization_petsc::CellPairSparsityStencilSnapshot3D>
        blocked_sparsity;
    require_petsc(
        discretization_petsc::make_cell_pair_sparsity_stencil_snapshot_3d(
            PETSC_COMM_WORLD,
            *overlap_materialized_schedule,
            overlap_partition,
            &materialized_sparsity),
        "build materialized cell-pair sparsity/stencil snapshot");
    require_petsc(
        discretization_petsc::make_cell_pair_sparsity_stencil_snapshot_3d(
            PETSC_COMM_WORLD,
            *overlap_blocked_schedule,
            overlap_partition,
            &blocked_sparsity),
        "build blocked-only cell-pair sparsity/stencil snapshot");
    require(
        materialized_sparsity.has_value() &&
            blocked_sparsity.has_value(),
        "cell-pair sparsity/stencil snapshots constructed");
    verify_cell_pair_sparsity_stencil_stage(
        *materialized_sparsity,
        overlap_partition,
        true);
    verify_cell_pair_sparsity_stencil_stage(
        *blocked_sparsity,
        overlap_partition,
        false);

    std::optional<
        discretization_petsc::PetscMpiAijSymbolicPreallocation3D>
        materialized_preallocation;
    std::optional<
        discretization_petsc::PetscMpiAijSymbolicPreallocation3D>
        blocked_preallocation;
    require_petsc(
        discretization_petsc::make_petsc_mpiaij_symbolic_preallocation_3d(
            overlap_dm,
            *materialized_sparsity,
            overlap_partition,
            overlap_identities,
            &materialized_preallocation),
        "build materialized symbolic MPIAIJ preallocation bridge");
    require_petsc(
        discretization_petsc::make_petsc_mpiaij_symbolic_preallocation_3d(
            overlap_dm,
            *blocked_sparsity,
            overlap_partition,
            overlap_identities,
            &blocked_preallocation),
        "build blocked-only symbolic MPIAIJ preallocation bridge");
    require(
        materialized_preallocation.has_value() &&
            blocked_preallocation.has_value(),
        "symbolic MPIAIJ preallocation bridges constructed");
    verify_petsc_mpiaij_symbolic_preallocation_stage(
        *materialized_preallocation,
        *materialized_sparsity,
        overlap_partition,
        true);
    verify_petsc_mpiaij_symbolic_preallocation_stage(
        *blocked_preallocation,
        *blocked_sparsity,
        overlap_partition,
        false);

    std::optional<
        discretization_petsc::OwnedCellStructuralColumnPatternSnapshot3D>
        materialized_column_pattern;
    std::optional<
        discretization_petsc::OwnedCellStructuralColumnPatternSnapshot3D>
        blocked_column_pattern;
    require_petsc(
        discretization_petsc::make_owned_cell_structural_column_pattern_snapshot_3d(
                *materialized_sparsity,
                *materialized_preallocation,
                overlap_partition,
                &materialized_column_pattern),
        "build materialized owned-cell structural column pattern");
    require_petsc(
        discretization_petsc::make_owned_cell_structural_column_pattern_snapshot_3d(
                *blocked_sparsity,
                *blocked_preallocation,
                overlap_partition,
                &blocked_column_pattern),
        "build blocked-only owned-cell structural column pattern");
    require(
        materialized_column_pattern
                .has_value() &&
            blocked_column_pattern
                .has_value(),
        "owned-cell structural column patterns constructed");
    verify_owned_cell_structural_column_pattern_stage(
        *materialized_column_pattern,
        *materialized_preallocation,
        *materialized_sparsity,
        overlap_partition,
        true);
    verify_owned_cell_structural_column_pattern_stage(
        *blocked_column_pattern,
        *blocked_preallocation,
        *blocked_sparsity,
        overlap_partition,
        false);

    Mat materialized_symbolic_matrix = nullptr;
    require_petsc(
        discretization_petsc::create_empty_petsc_mpiaij_symbolic_matrix_3d(
                PETSC_COMM_WORLD,
                *materialized_preallocation,
                &materialized_symbolic_matrix),
        "create materialized empty symbolic MPIAIJ matrix");
    verify_empty_petsc_mpiaij_symbolic_matrix_stage(
        materialized_symbolic_matrix,
        *materialized_preallocation);
    require_petsc(
        MatDestroy(
            &materialized_symbolic_matrix),
        "destroy materialized empty symbolic MPIAIJ matrix");

    Mat blocked_symbolic_matrix = nullptr;
    require_petsc(
        discretization_petsc::create_empty_petsc_mpiaij_symbolic_matrix_3d(
                PETSC_COMM_WORLD,
                *blocked_preallocation,
                &blocked_symbolic_matrix),
        "create blocked-only empty symbolic MPIAIJ matrix");
    verify_empty_petsc_mpiaij_symbolic_matrix_stage(
        blocked_symbolic_matrix,
        *blocked_preallocation);
    require_petsc(
        MatDestroy(
            &blocked_symbolic_matrix),
        "destroy blocked-only empty symbolic MPIAIJ matrix");

    std::vector<
        discretization_petsc::AssemblyReadyInternalConnectionRow3D>
        duplicate_pair_authoritative_rows;
    if (!overlap_materialized_schedule
             ->assembly_rows()
             .empty()) {
        const auto source_row =
            overlap_materialized_schedule
                ->assembly_rows()
                .front();
        auto duplicate_row =
            source_row;
        duplicate_row.face_global =
            mesh::GlobalEntityId{
                source_row.face_global.value() +
                1000000ULL};
        duplicate_pair_authoritative_rows.push_back(
            source_row);
        duplicate_pair_authoritative_rows.push_back(
            duplicate_row);
    }
    const auto duplicate_pair_schedule =
        discretization_petsc::ParallelOwnedConnectionSchedule3D{
            overlap_partition.local_rank(),
            overlap_partition.rank_count(),
            std::move(
                duplicate_pair_authoritative_rows),
            {}};
    std::optional<
        discretization_petsc::CellPairSparsityStencilSnapshot3D>
        deduplicated_sparsity;
    require_petsc(
        discretization_petsc::make_cell_pair_sparsity_stencil_snapshot_3d(
            PETSC_COMM_WORLD,
            duplicate_pair_schedule,
            overlap_partition,
            &deduplicated_sparsity),
        "deduplicate multiple authoritative face rows sharing one stable cell pair");
    require(
        deduplicated_sparsity.has_value(),
        "deduplicated cell-pair sparsity snapshot constructed");
    verify_cell_pair_sparsity_stencil_stage(
        *deduplicated_sparsity,
        overlap_partition,
        true);

    std::optional<mesh::FaceGeometry3D>
        overlap_materialized_face_geometry;
    require_petsc(
        mesh_petsc::materialize_face_geometry_3d(
            overlap_face_geometry,
            overlap_identities,
            &overlap_materialized_face_geometry),
        "materialize processed GRDECL overlap FaceGeometry3D");
    require(
        overlap_materialized_face_geometry.has_value(),
        "depth-one overlap must resolve every stable face owner to target-local cell");
    verify_materialized_face_geometry_stage(
        *overlap_materialized_face_geometry,
        overlap_face_geometry,
        field_reference,
        overlap_identities);

    require_petsc(
        PetscSFDestroy(&overlap_migration_sf),
        "destroy processed GRDECL overlap migration SF");

    require_petsc(
        DMDestroy(&overlap_dm),
        "destroy processed GRDECL overlap DM");
    require_petsc(
        DMDestroy(&distributed_dm),
        "destroy processed GRDECL distributed DM");
}

std::string gmsh_petsc_chain_fixture() {
    return R"msh($MeshFormat
4.1 0 8
$EndMeshFormat
$PhysicalNames
5
1 11 "left"
1 12 "bottom"
1 13 "right"
1 14 "top"
2 21 "domain"
$EndPhysicalNames
$Entities
0 5 1 0
1 0 0 0 1 0 0 1 12 0
2 0 0 0 1 1 0 1 11 0
3 1 0 0 3 0 0 1 12 0
4 3 0 0 3 1 0 1 13 0
5 1 1 0 3 1 0 1 14 0
100 0 0 0 3 1 0 1 21 5 1 2 3 4 5
$EndEntities
$Nodes
1 5 10 50
2 100 0 5
50
10
40
20
30
3 1 0
0 0 0
3 0 0
1 0 0
1 1 0
$EndNodes
$Elements
7 7 101 202
1 1 1 1
101 10 20
1 2 1 1
103 30 10
1 3 1 1
104 20 40
1 4 1 1
105 40 50
1 5 1 1
106 50 30
2 100 2 1
201 10 20 30
2 100 3 1
202 20 40 50 30
$EndElements
)msh";
}

mesh::LocalIndex reference_local_by_global(
    const mesh::Topology& topology,
    mesh::EntityKind kind,
    mesh::GlobalEntityId global) {
    const auto ids = topology.global_ids(kind);
    for (std::size_t local = 0U;
         local < ids.size();
         ++local) {
        if (ids[local] == global) {
            return mesh::LocalIndex{
                static_cast<
                    mesh::LocalIndex::value_type>(
                        local)};
        }
    }
    throw std::runtime_error(
        "stable GlobalEntityId missing from imported reference topology");
}

void verify_imported_boundary_by_global(
    const mesh::FaceBoundarySnapshot& actual,
    const std::vector<
        mesh_petsc::DMPlexPointIdentity>& identities,
    const mesh::Topology& reference_topology,
    const mesh::FaceBoundarySnapshot& reference) {
    std::size_t count = 0U;
    for (const auto& identity : identities) {
        if (identity.kind !=
            mesh::EntityKind::face) {
            continue;
        }
        ++count;
        const auto source_local =
            reference_local_by_global(
                reference_topology,
                mesh::EntityKind::face,
                identity.global);
        require(
            actual.classification(
                identity.local) ==
                reference.classification(
                    source_local),
            "imported DMPlex face classification by stable ID");
        require(
            actual.physical_tag(
                identity.local) ==
                reference.physical_tag(
                    source_local),
            "imported DMPlex PhysicalTag by stable ID");
    }
    require(
        actual.face_count() == count,
        "imported DMPlex boundary target face count");
}

void require_field_metadata_same(
    const mesh::DenseFieldMetadata& left,
    const mesh::DenseFieldMetadata& right) {
    require(
        left.id == right.id &&
            left.unit == right.unit &&
            left.source.kind ==
                right.source.kind &&
            left.source.reference ==
                right.source.reference &&
            left.source.revision ==
                right.source.revision &&
            left.source.locator ==
                right.source.locator,
        "imported DMPlex property metadata");
}

void verify_imported_cell_field_by_global(
    const mesh::DenseFieldSnapshot& actual,
    const std::vector<
        mesh_petsc::DMPlexPointIdentity>& identities,
    const mesh::Topology& reference_topology,
    const mesh::DenseFieldSnapshot& reference) {
    require(
        actual.location() ==
            mesh::EntityKind::cell &&
        actual.component_count() ==
            reference.component_count(),
        "imported DMPlex cell field layout");
    require_field_metadata_same(
        actual.metadata(),
        reference.metadata());

    std::size_t count = 0U;
    for (const auto& identity : identities) {
        if (identity.kind !=
            mesh::EntityKind::cell) {
            continue;
        }
        ++count;
        const auto source_local =
            reference_local_by_global(
                reference_topology,
                mesh::EntityKind::cell,
                identity.global);
        for (std::size_t component = 0U;
             component <
                 actual.component_count();
             ++component) {
            require(
                actual.value(
                    identity.local,
                    component) ==
                    reference.value(
                        source_local,
                        component),
                "imported DMPlex cell property by stable ID");
        }
    }
    require(
        actual.entity_count() == count,
        "imported DMPlex cell field entity count");
}

void verify_imported_geometry_by_global(
    DM dm,
    const std::vector<
        mesh_petsc::DMPlexPointIdentity>& identities,
    const mesh::Topology& reference_topology,
    const mesh::Geometry2D& reference_geometry) {
    constexpr double tolerance = 1.0e-12;
    auto view = get_dmplex_coordinate_view(dm);

    for (const auto& identity : identities) {
        if (identity.kind !=
            mesh::EntityKind::vertex) {
            continue;
        }
        const auto source_local =
            reference_local_by_global(
                reference_topology,
                mesh::EntityKind::vertex,
                identity.global);
        const auto expected =
            reference_geometry.vertex_coordinate_m(
                source_local);
        const auto actual =
            coordinate_for_dmplex_vertex(
                view, identity.point);
        require(
            std::abs(actual.x_m - expected.x_m) <=
                    tolerance &&
                std::abs(actual.y_m - expected.y_m) <=
                    tolerance,
            "imported distributed vertex coordinate by stable ID");
    }

    for (const auto& identity : identities) {
        if (identity.kind !=
            mesh::EntityKind::face) {
            continue;
        }
        PetscInt cone_size = -1;
        const PetscInt* cone = nullptr;
        require_petsc(
            DMPlexGetConeSize(
                dm,
                identity.point,
                &cone_size),
            "imported face cone size");
        require_petsc(
            DMPlexGetCone(
                dm,
                identity.point,
                &cone),
            "imported face cone");
        require(
            cone_size == 2 &&
                cone != nullptr,
            "imported face has two vertices");
        const auto a =
            coordinate_for_dmplex_vertex(
                view, cone[0]);
        const auto b =
            coordinate_for_dmplex_vertex(
                view, cone[1]);
        const double length =
            std::hypot(
                b.x_m - a.x_m,
                b.y_m - a.y_m);
        const auto source_local =
            reference_local_by_global(
                reference_topology,
                mesh::EntityKind::face,
                identity.global);
        require(
            std::abs(
                length -
                reference_geometry.face_length_m(
                    source_local)) <=
                tolerance,
            "imported distributed face length by stable ID");
    }

    for (const auto& identity : identities) {
        if (identity.kind !=
            mesh::EntityKind::cell) {
            continue;
        }
        const auto source_local =
            reference_local_by_global(
                reference_topology,
                mesh::EntityKind::cell,
                identity.global);
        const std::size_t expected_vertex_count =
            reference_topology
                .relation(
                    mesh::EntityKind::cell,
                    mesh::EntityKind::vertex)
                .adjacent(source_local)
                .size();
        require(
            expected_vertex_count == 3U ||
                expected_vertex_count == 4U,
            "imported reference cell degree");

        DMPolytopeType type =
            DM_POLYTOPE_UNKNOWN;
        require_petsc(
            DMPlexGetCellType(
                dm,
                identity.point,
                &type),
            "imported distributed cell type");
        require(
            type ==
                (expected_vertex_count == 3U
                     ? DM_POLYTOPE_TRIANGLE
                     : DM_POLYTOPE_QUADRILATERAL),
            "imported triangle/quad cell type after DMPlex distribution");

        PetscInt closure_size = 0;
        PetscInt* closure = nullptr;
        require_petsc(
            DMPlexGetTransitiveClosure(
                dm,
                identity.point,
                PETSC_TRUE,
                &closure_size,
                &closure),
            "imported cell closure");

        std::vector<PetscInt> vertices;
        vertices.reserve(expected_vertex_count);
        for (PetscInt i = 0;
             i < closure_size;
             ++i) {
            const PetscInt point =
                closure[2 * i];
            const auto& closure_identity =
                identity_for_point(
                    identities, point);
            if (closure_identity.kind ==
                    mesh::EntityKind::vertex &&
                std::find(
                    vertices.begin(),
                    vertices.end(),
                    point) ==
                    vertices.end()) {
                vertices.push_back(point);
            }
        }
        require_petsc(
            DMPlexRestoreTransitiveClosure(
                dm,
                identity.point,
                PETSC_TRUE,
                &closure_size,
                &closure),
            "restore imported cell closure");
        require(
            vertices.size() ==
                expected_vertex_count,
            "imported cell closure vertex count");

        std::vector<mesh::Coordinate2D>
            coordinates;
        coordinates.reserve(vertices.size());
        double mean_x = 0.0;
        double mean_y = 0.0;
        for (const PetscInt vertex :
             vertices) {
            const auto coordinate =
                coordinate_for_dmplex_vertex(
                    view, vertex);
            coordinates.push_back(coordinate);
            mean_x += coordinate.x_m;
            mean_y += coordinate.y_m;
        }
        mean_x /=
            static_cast<double>(
                coordinates.size());
        mean_y /=
            static_cast<double>(
                coordinates.size());

        std::sort(
            coordinates.begin(),
            coordinates.end(),
            [mean_x, mean_y](
                const mesh::Coordinate2D& left,
                const mesh::Coordinate2D& right) {
                return std::atan2(
                           left.y_m - mean_y,
                           left.x_m - mean_x) <
                       std::atan2(
                           right.y_m - mean_y,
                           right.x_m - mean_x);
            });

        double twice_area = 0.0;
        double centroid_numerator_x = 0.0;
        double centroid_numerator_y = 0.0;
        for (std::size_t i = 0U;
             i < coordinates.size();
             ++i) {
            const auto& a = coordinates[i];
            const auto& b =
                coordinates[
                    (i + 1U) %
                    coordinates.size()];
            const double cross =
                a.x_m * b.y_m -
                b.x_m * a.y_m;
            twice_area += cross;
            centroid_numerator_x +=
                (a.x_m + b.x_m) * cross;
            centroid_numerator_y +=
                (a.y_m + b.y_m) * cross;
        }
        require(
            std::abs(twice_area) >
                tolerance,
            "imported distributed cell nonzero area");
        const double centroid_x =
            centroid_numerator_x /
            (3.0 * twice_area);
        const double centroid_y =
            centroid_numerator_y /
            (3.0 * twice_area);
        const double area =
            0.5 * std::abs(twice_area);

        const auto expected_centroid =
            reference_geometry.cell_centroid_m(
                source_local);
        require(
            std::abs(
                centroid_x -
                expected_centroid.x_m) <=
                    tolerance &&
                std::abs(
                    centroid_y -
                expected_centroid.y_m) <=
                    tolerance,
            "imported distributed cell centroid by stable ID");
        require(
            std::abs(
                area -
                reference_geometry.cell_area_m2(
                    source_local)) <=
                tolerance,
            "imported distributed cell area by stable ID");
    }

    restore_dmplex_coordinate_view(&view);
}

mesh::DenseFieldMetadata gmsh_chain_field_metadata() {
    return mesh::DenseFieldMetadata{
        "gmsh.import.cell.scalar",
        "kg/m3",
        mesh::FieldSourceMetadata{
            mesh::FieldSourceKind::synthetic_test,
            "gmsh_mixed_triangle_quad_fixture",
            "msh4.1-v1",
            "$Elements/2D"}};
}

mesh::DenseFieldSnapshot gmsh_chain_source_field(
    const mesh::Topology& topology) {
    std::vector<double> values;
    values.reserve(
        topology.entity_count(
            mesh::EntityKind::cell));
    for (const auto id :
         topology.global_ids(
             mesh::EntityKind::cell)) {
        values.push_back(
            0.01 *
            static_cast<double>(
                id.value()));
    }
    return mesh::DenseFieldSnapshot::create(
        topology,
        mesh::EntityKind::cell,
        1U,
        std::move(values),
        gmsh_chain_field_metadata());
}

void verify_gmsh_import_through_dmplex_chain() {
    int mpi_rank = -1;
    int mpi_size = -1;
    require(
        MPI_Comm_rank(
            PETSC_COMM_WORLD,
            &mpi_rank) == MPI_SUCCESS,
        "MPI rank for Gmsh PETSc chain");
    require(
        MPI_Comm_size(
            PETSC_COMM_WORLD,
            &mpi_size) == MPI_SUCCESS &&
            mpi_size == 2,
        "Gmsh PETSc chain requires two ranks");

    const auto reference_import =
        mesh::import_gmsh_4_1_ascii(
            gmsh_petsc_chain_fixture(),
            2.0);

    mesh::Topology empty_topology{
        mesh::Topology::EntityIds{}, {}};
    const auto source_boundary =
        mpi_rank == 0
            ? reference_import.face_boundary
            : mesh::FaceBoundarySnapshot{
                  {}, {}};
    const auto source_field =
        mpi_rank == 0
            ? gmsh_chain_source_field(
                  reference_import.topology)
            : mesh::DenseFieldSnapshot::create(
                  empty_topology,
                  mesh::EntityKind::cell,
                  1U,
                  {},
                  gmsh_chain_field_metadata());

    DM source_dm = nullptr;
    std::vector<mesh_petsc::DMPlexPointIdentity>
        source_identities;
    require_petsc(
        mesh_petsc::create_root_dmplex_topology(
            PETSC_COMM_WORLD,
            0,
            mpi_rank == 0
                ? &reference_import.topology
                : nullptr,
            &source_dm,
            &source_identities),
        "create rooted Gmsh DMPlex");

    require_petsc(
        mesh_petsc::attach_root_geometry2d_coordinates(
            source_dm,
            0,
            mpi_rank == 0
                ? &reference_import.geometry
                : nullptr,
            source_identities),
        "attach imported Gmsh coordinates");

    PetscPartitioner partitioner = nullptr;
    require_petsc(
        DMPlexGetPartitioner(
            source_dm,
            &partitioner),
        "get Gmsh DMPlex partitioner");
    require_petsc(
        PetscPartitionerSetType(
            partitioner,
            PETSCPARTITIONERSIMPLE),
        "set simple Gmsh partitioner");

    PetscSF migration_sf = nullptr;
    DM distributed_dm = nullptr;
    require_petsc(
        DMPlexDistribute(
            source_dm,
            0,
            &migration_sf,
            &distributed_dm),
        "distribute imported Gmsh DMPlex");
    require(
        distributed_dm != nullptr &&
            migration_sf != nullptr,
        "Gmsh distribution outputs");

    std::vector<mesh_petsc::DMPlexPointIdentity>
        distributed_identities;
    require_petsc(
        mesh_petsc::migrate_dmplex_identities(
            source_dm,
            migration_sf,
            source_identities,
            distributed_dm,
            &distributed_identities),
        "migrate imported Gmsh identities");

    std::optional<mesh::FaceBoundarySnapshot>
        distributed_boundary;
    std::optional<mesh::DenseFieldSnapshot>
        distributed_field;
    require_petsc(
        mesh_petsc::migrate_face_boundary_snapshot(
            source_dm,
            migration_sf,
            source_boundary,
            source_identities,
            distributed_dm,
            distributed_identities,
            &distributed_boundary),
        "migrate imported Gmsh boundary");
    require_petsc(
        mesh_petsc::migrate_dense_field_snapshot(
            source_dm,
            migration_sf,
            source_field,
            source_identities,
            distributed_dm,
            distributed_identities,
            &distributed_field),
        "migrate imported Gmsh property");

    require(
        distributed_boundary.has_value() &&
            distributed_field.has_value(),
        "distributed imported snapshots");

    const auto& reference_topology =
        reference_import.topology;
    const auto& reference_geometry =
        reference_import.geometry;
    const auto& reference_boundary =
        reference_import.face_boundary;
    const auto reference_field =
        gmsh_chain_source_field(
            reference_topology);

    verify_imported_geometry_by_global(
        distributed_dm,
        distributed_identities,
        reference_topology,
        reference_geometry);
    verify_imported_boundary_by_global(
        *distributed_boundary,
        distributed_identities,
        reference_topology,
        reference_boundary);
    verify_imported_cell_field_by_global(
        *distributed_field,
        distributed_identities,
        reference_topology,
        reference_field);

    PetscInt local_triangles = 0;
    PetscInt local_quads = 0;
    for (const auto& identity :
         distributed_identities) {
        if (identity.kind !=
            mesh::EntityKind::cell) {
            continue;
        }
        DMPolytopeType type =
            DM_POLYTOPE_UNKNOWN;
        require_petsc(
            DMPlexGetCellType(
                distributed_dm,
                identity.point,
                &type),
            "imported distributed polytope");
        if (type == DM_POLYTOPE_TRIANGLE) {
            ++local_triangles;
        } else if (
            type ==
            DM_POLYTOPE_QUADRILATERAL) {
            ++local_quads;
        } else {
            throw std::runtime_error(
                "imported distributed cell is neither triangle nor quad");
        }
    }
    PetscInt global_triangles =
        local_triangles;
    PetscInt global_quads =
        local_quads;
    require(
        MPI_Allreduce(
            MPI_IN_PLACE,
            &global_triangles,
            1,
            MPIU_INT,
            MPI_SUM,
            PETSC_COMM_WORLD) == MPI_SUCCESS,
        "count distributed triangles");
    require(
        MPI_Allreduce(
            MPI_IN_PLACE,
            &global_quads,
            1,
            MPIU_INT,
            MPI_SUM,
            PETSC_COMM_WORLD) == MPI_SUCCESS,
        "count distributed quads");
    require(
        global_triangles == 1 &&
            global_quads == 1,
        "Gmsh mixed triangle/quad types survive distribution");

    PetscSF overlap_sf = nullptr;
    DM overlap_dm = nullptr;
    require_petsc(
        DMPlexDistributeOverlap(
            distributed_dm,
            1,
            &overlap_sf,
            &overlap_dm),
        "overlap imported Gmsh DMPlex");
    require(
        overlap_sf != nullptr &&
            overlap_dm != nullptr,
        "Gmsh overlap outputs");

    std::vector<mesh_petsc::DMPlexPointIdentity>
        overlap_identities;
    require_petsc(
        mesh_petsc::migrate_dmplex_identities(
            distributed_dm,
            overlap_sf,
            distributed_identities,
            overlap_dm,
            &overlap_identities),
        "migrate Gmsh identities into overlap");

    std::optional<mesh::FaceBoundarySnapshot>
        overlap_boundary;
    std::optional<mesh::DenseFieldSnapshot>
        overlap_field;
    require_petsc(
        mesh_petsc::migrate_face_boundary_snapshot(
            distributed_dm,
            overlap_sf,
            *distributed_boundary,
            distributed_identities,
            overlap_dm,
            overlap_identities,
            &overlap_boundary),
        "migrate Gmsh boundary into overlap");
    require_petsc(
        mesh_petsc::migrate_dense_field_snapshot(
            distributed_dm,
            overlap_sf,
            *distributed_field,
            distributed_identities,
            overlap_dm,
            overlap_identities,
            &overlap_field),
        "migrate Gmsh property into overlap");

    require(
        overlap_boundary.has_value() &&
            overlap_field.has_value(),
        "overlap imported snapshots");
    verify_imported_geometry_by_global(
        overlap_dm,
        overlap_identities,
        reference_topology,
        reference_geometry);
    verify_imported_boundary_by_global(
        *overlap_boundary,
        overlap_identities,
        reference_topology,
        reference_boundary);
    verify_imported_cell_field_by_global(
        *overlap_field,
        overlap_identities,
        reference_topology,
        reference_field);

    const auto generated_internal =
        reference_local_by_global(
            reference_topology,
            mesh::EntityKind::face,
            mesh::GlobalEntityId{203U});
    require(
        reference_boundary.classification(
            generated_internal) ==
            mesh::FaceClassification::interior &&
        !reference_boundary.has_physical_tag(
            generated_internal),
        "generated Gmsh internal face remains interior through PETSc chain");

    require_petsc(
        PetscSFDestroy(&overlap_sf),
        "destroy Gmsh overlap migration SF");
    require_petsc(
        DMDestroy(&overlap_dm),
        "destroy Gmsh overlap DM");
    require_petsc(
        PetscSFDestroy(&migration_sf),
        "destroy Gmsh distribution migration SF");
    require_petsc(
        DMDestroy(&distributed_dm),
        "destroy Gmsh distributed DM");
    require_petsc(
        DMDestroy(&source_dm),
        "destroy Gmsh source DM");
}

void run_three_rank_sparsity_test() {
    int mpi_rank = -1;
    int mpi_size = -1;
    require(
        MPI_Comm_rank(
            PETSC_COMM_WORLD,
            &mpi_rank) == MPI_SUCCESS,
        "3-rank sparsity MPI_Comm_rank");
    require(
        MPI_Comm_size(
            PETSC_COMM_WORLD,
            &mpi_size) == MPI_SUCCESS,
        "3-rank sparsity MPI_Comm_size");
    require(
        mpi_size == 3,
        "owner-targeted sparsity test requires exactly three MPI ranks");
    require(
        mpi_rank >= 0 && mpi_rank < 3,
        "unexpected 3-rank sparsity MPI rank");

    const auto rank =
        static_cast<std::uint32_t>(
            mpi_rank);

    mesh::Topology::EntityIds ids;
    ids.cells = {
        mesh::GlobalEntityId{1000U},
        mesh::GlobalEntityId{2000U}};
    const mesh::Topology topology{
        std::move(ids),
        {}};

    mesh::EntityOwnerRanks owners;
    owners.cells = {
        mesh::PartitionRank{0U},
        mesh::PartitionRank{2U}};
    const auto partition =
        mesh::PartitionSnapshot::create(
            topology,
            mesh::PartitionRank{rank},
            3U,
            std::move(owners));

    const auto make_row =
        [](mesh::LocalIndex face,
           std::uint64_t face_global,
           double transmissibility_m3) {
            return discretization_petsc::AssemblyReadyInternalConnectionRow3D{
                    face,
                    mesh::GlobalEntityId{
                        face_global},
                    mesh::LocalIndex{0U},
                    mesh::GlobalEntityId{1000U},
                    mesh::LocalIndex{1U},
                    mesh::GlobalEntityId{2000U},
                    transmissibility_m3};
        };

    std::vector<
        discretization_petsc::AssemblyReadyInternalConnectionRow3D>
        authoritative_rows;
    std::vector<
        discretization_petsc::AssemblyReadyInternalConnectionRow3D>
        ghost_rows;

    if (rank == 1U) {
        authoritative_rows.push_back(
            make_row(
                mesh::LocalIndex{0U},
                9000U,
                1.0));
        authoritative_rows.push_back(
            make_row(
                mesh::LocalIndex{1U},
                9001U,
                2.0));
    } else {
        ghost_rows.push_back(
            make_row(
                mesh::LocalIndex{0U},
                9000U,
                11.0));
        ghost_rows.push_back(
            make_row(
                mesh::LocalIndex{1U},
                9001U,
                22.0));
    }

    const discretization_petsc::ParallelOwnedConnectionSchedule3D
        schedule{
            mesh::PartitionRank{rank},
            3U,
            std::move(authoritative_rows),
            std::move(ghost_rows)};

    std::optional<
        discretization_petsc::CellPairSparsityStencilSnapshot3D>
        snapshot;
    require_petsc(
        discretization_petsc::make_cell_pair_sparsity_stencil_snapshot_3d(
                PETSC_COMM_WORLD,
                schedule,
                partition,
                &snapshot),
        "owner-targeted 3-rank cell-pair sparsity");
    require(
        snapshot.has_value(),
        "owner-targeted 3-rank sparsity snapshot");

    if (rank == 1U) {
        require(
            partition.owned_count(
                mesh::EntityKind::cell) == 0U &&
                snapshot->coupling_count() == 0U &&
                snapshot->owned_cell_count() == 0U,
            "authoritative face rank that owns neither endpoint stores no cell-pair graph");
    } else {
        require(
            partition.owned_count(
                mesh::EntityKind::cell) == 1U &&
                partition.ghost_count(
                    mesh::EntityKind::cell) == 1U &&
                snapshot->coupling_count() == 1U &&
                snapshot->owned_cell_count() == 1U,
            "endpoint owner receives exactly one deduplicated coupling");

        const auto& pair =
            snapshot->couplings().front();
        require(
            pair.first_cell ==
                    mesh::LocalIndex{0U} &&
                pair.second_cell ==
                    mesh::LocalIndex{1U} &&
                pair.first_cell_global ==
                    mesh::GlobalEntityId{1000U} &&
                pair.second_cell_global ==
                    mesh::GlobalEntityId{2000U},
            "owner-targeted coupling preserves canonical stable cell pair");

        const auto owned_cell =
            rank == 0U
                ? mesh::LocalIndex{0U}
                : mesh::LocalIndex{1U};
        const auto ghost_cell =
            rank == 0U
                ? mesh::LocalIndex{1U}
                : mesh::LocalIndex{0U};
        require(
            partition.is_owned(
                mesh::EntityKind::cell,
                owned_cell) &&
                partition.is_ghost(
                    mesh::EntityKind::cell,
                    ghost_cell),
            "3-rank fixture endpoint ownership");

        const auto& counts =
            snapshot->structural_counts(
                owned_cell);
        require(
            counts.cell ==
                    owned_cell &&
                counts.cell_global ==
                    partition.global_id(
                        mesh::EntityKind::cell,
                        owned_cell) &&
                counts.diagonal_block_nnz == 1U &&
                counts.off_diagonal_block_nnz == 1U,
            "endpoint owner structural counts remain one self plus one remote neighbour");
    }

    const std::uint64_t local_couplings =
        static_cast<std::uint64_t>(
            snapshot->coupling_count());
    std::array<std::uint64_t, 3>
        all_couplings{};
    require(
        MPI_Allgather(
            &local_couplings,
            1,
            MPI_UINT64_T,
            all_couplings.data(),
            1,
            MPI_UINT64_T,
            PETSC_COMM_WORLD) == MPI_SUCCESS,
        "gather 3-rank owner-targeted local coupling counts");
    require(
        all_couplings[0] == 1U &&
            all_couplings[1] == 0U &&
            all_couplings[2] == 1U,
        "stable pair payload is present only on endpoint owner ranks");

    require(
        MPI_Barrier(
            PETSC_COMM_WORLD) == MPI_SUCCESS,
        "3-rank sparsity MPI_Barrier");
}

void verify_structural_column_pattern_sorting_fixture(
    int mpi_rank) {
    require(
        mpi_rank == 0 ||
            mpi_rank == 1,
        "column-pattern sorting fixture requires rank 0 or 1");

    const auto rank =
        static_cast<std::uint32_t>(
            mpi_rank);

    mesh::Topology::EntityIds ids;
    ids.cells = {
        mesh::GlobalEntityId{100U},
        mesh::GlobalEntityId{200U},
        mesh::GlobalEntityId{300U},
        mesh::GlobalEntityId{400U}};
    const mesh::Topology topology{
        std::move(ids),
        {}};

    mesh::EntityOwnerRanks owners;
    owners.cells = {
        mesh::PartitionRank{0U},
        mesh::PartitionRank{1U},
        mesh::PartitionRank{0U},
        mesh::PartitionRank{1U}};
    const auto partition =
        mesh::PartitionSnapshot::create(
            topology,
            mesh::PartitionRank{rank},
            2U,
            std::move(owners));

    const auto coupling =
        [](std::uint32_t first_local,
           std::uint64_t first_global,
           std::uint32_t second_local,
           std::uint64_t second_global) {
            return discretization_petsc::CellPairCoupling3D{
                    mesh::LocalIndex{
                        first_local},
                    mesh::GlobalEntityId{
                        first_global},
                    mesh::LocalIndex{
                        second_local},
                    mesh::GlobalEntityId{
                        second_global}};
        };

    std::vector<
        discretization_petsc::CellPairCoupling3D>
        couplings;
    std::vector<
        discretization_petsc::OwnedCellStructuralCounts3D>
        counts;
    std::vector<mesh::LocalIndex>
        owned_cells;
    std::vector<mesh::GlobalEntityId>
        owned_ids;
    std::vector<PetscInt>
        owned_rows;
    std::vector<PetscInt>
        d_nnz;
    std::vector<PetscInt>
        o_nnz;

    if (rank == 0U) {
        couplings = {
            coupling(2U, 300U, 3U, 400U),
            coupling(0U, 100U, 3U, 400U),
            coupling(0U, 100U, 2U, 300U),
            coupling(0U, 100U, 1U, 200U)};
        counts = {
            {mesh::LocalIndex{0U},
             mesh::GlobalEntityId{100U},
             2U,
             2U},
            {mesh::LocalIndex{2U},
             mesh::GlobalEntityId{300U},
             2U,
             1U}};
        owned_cells = {
            mesh::LocalIndex{0U},
            mesh::LocalIndex{2U}};
        owned_ids = {
            mesh::GlobalEntityId{100U},
            mesh::GlobalEntityId{300U}};
        owned_rows = {0, 1};
        d_nnz = {2, 2};
        o_nnz = {2, 1};
    } else {
        couplings = {
            coupling(0U, 100U, 3U, 400U),
            coupling(1U, 200U, 3U, 400U),
            coupling(0U, 100U, 1U, 200U),
            coupling(2U, 300U, 3U, 400U)};
        counts = {
            {mesh::LocalIndex{1U},
             mesh::GlobalEntityId{200U},
             2U,
             1U},
            {mesh::LocalIndex{3U},
             mesh::GlobalEntityId{400U},
             2U,
             2U}};
        owned_cells = {
            mesh::LocalIndex{1U},
            mesh::LocalIndex{3U}};
        owned_ids = {
            mesh::GlobalEntityId{200U},
            mesh::GlobalEntityId{400U}};
        owned_rows = {2, 3};
        d_nnz = {2, 2};
        o_nnz = {1, 2};
    }

    const discretization_petsc::CellPairSparsityStencilSnapshot3D
        sparsity{
            mesh::PartitionRank{rank},
            2U,
            4U,
            std::move(couplings),
            std::move(counts)};

    const discretization_petsc::PetscMpiAijSymbolicPreallocation3D
        bridge{
            mesh::PartitionRank{rank},
            2U,
            rank == 0U ? 0 : 2,
            rank == 0U ? 2 : 4,
            4,
            std::move(owned_cells),
            std::move(owned_ids),
            std::move(owned_rows),
            std::move(d_nnz),
            std::move(o_nnz),
            {0, 2, 1, 3}};

    std::optional<
        discretization_petsc::OwnedCellStructuralColumnPatternSnapshot3D>
        snapshot;
    require_petsc(
        discretization_petsc::make_owned_cell_structural_column_pattern_snapshot_3d(
                sparsity,
                bridge,
                partition,
                &snapshot),
        "build multi-neighbour sorted structural column pattern");
    require(
        snapshot.has_value() &&
            snapshot->row_count() == 2U,
        "multi-neighbour structural column-pattern row count");

    const auto matches =
        [](std::span<const PetscInt> actual,
           std::span<const PetscInt> expected) {
            return actual.size() ==
                       expected.size() &&
                   std::equal(
                       actual.begin(),
                       actual.end(),
                       expected.begin(),
                       expected.end());
        };

    if (rank == 0U) {
        const std::array<PetscInt, 2>
            first_diagonal{0, 1};
        const std::array<PetscInt, 2>
            first_off_diagonal{2, 3};
        const std::array<PetscInt, 2>
            second_diagonal{0, 1};
        const std::array<PetscInt, 1>
            second_off_diagonal{3};
        require(
            matches(
                snapshot->diagonal_global_columns(
                    snapshot->rows()[0].cell),
                first_diagonal) &&
            matches(
                snapshot->off_diagonal_global_columns(
                    snapshot->rows()[0].cell),
                first_off_diagonal) &&
            matches(
                snapshot->diagonal_global_columns(
                    snapshot->rows()[1].cell),
                second_diagonal) &&
            matches(
                snapshot->off_diagonal_global_columns(
                    snapshot->rows()[1].cell),
                second_off_diagonal),
            "rank 0 structural columns are globally sorted and block-classified");
    } else {
        const std::array<PetscInt, 2>
            first_diagonal{2, 3};
        const std::array<PetscInt, 1>
            first_off_diagonal{0};
        const std::array<PetscInt, 2>
            second_diagonal{2, 3};
        const std::array<PetscInt, 2>
            second_off_diagonal{0, 1};
        require(
            matches(
                snapshot->diagonal_global_columns(
                    snapshot->rows()[0].cell),
                first_diagonal) &&
            matches(
                snapshot->off_diagonal_global_columns(
                    snapshot->rows()[0].cell),
                first_off_diagonal) &&
            matches(
                snapshot->diagonal_global_columns(
                    snapshot->rows()[1].cell),
                second_diagonal) &&
            matches(
                snapshot->off_diagonal_global_columns(
                    snapshot->rows()[1].cell),
                second_off_diagonal),
            "rank 1 structural columns are globally sorted and block-classified");
    }

    for (std::size_t index = 0U;
         index < snapshot->rows().size();
         ++index) {
        const auto cell =
            snapshot->rows()[index].cell;
        const auto diagonal_columns =
            snapshot->diagonal_global_columns(
                cell);
        const auto off_diagonal_columns =
            snapshot->off_diagonal_global_columns(
                cell);
        require(
            diagonal_columns.size() +
                    off_diagonal_columns.size() ==
                static_cast<std::size_t>(
                    bridge.diagonal_nnz()[index] +
                    bridge.off_diagonal_nnz()[index]),
            "multi-neighbour N_columns,row equals d_nnz + o_nnz");
    }

    require(
        snapshot->diagonal_column_storage().size() ==
                4U &&
            snapshot->off_diagonal_column_storage().size() ==
                3U,
        "multi-neighbour finalized column storage is compact across rows");
}

void verify_linear_3d_cell_families_dmplex(
    int mpi_rank,
    int mpi_size) {
    require(
        mpi_size == 2 &&
            (mpi_rank == 0 || mpi_rank == 1),
        "linear 3D PETSc fixture requires two MPI ranks");

    const std::vector<mesh::GlobalEntityId>
        vertex_ids{
            mesh::GlobalEntityId{1U},
            mesh::GlobalEntityId{2U},
            mesh::GlobalEntityId{3U},
            mesh::GlobalEntityId{4U},
            mesh::GlobalEntityId{5U},
            mesh::GlobalEntityId{6U},
            mesh::GlobalEntityId{7U},
            mesh::GlobalEntityId{8U},
            mesh::GlobalEntityId{9U},
            mesh::GlobalEntityId{10U},
            mesh::GlobalEntityId{11U},
            mesh::GlobalEntityId{12U},
            mesh::GlobalEntityId{13U},
            mesh::GlobalEntityId{14U},
            mesh::GlobalEntityId{15U},
            mesh::GlobalEntityId{16U},
            mesh::GlobalEntityId{17U},
            mesh::GlobalEntityId{18U},
            mesh::GlobalEntityId{19U},
            mesh::GlobalEntityId{20U},
            mesh::GlobalEntityId{21U},
            mesh::GlobalEntityId{22U},
            mesh::GlobalEntityId{23U}};
    const std::vector<mesh::Coordinate3D>
        coordinates{
            {0.0, 0.0, 0.0},
            {1.0, 0.0, 0.0},
            {0.0, 1.0, 0.0},
            {0.0, 0.0, 1.0},
            {2.0, 0.0, 0.0},
            {3.0, 0.0, 0.0},
            {3.0, 1.0, 0.0},
            {2.0, 1.0, 0.0},
            {2.0, 0.0, 1.0},
            {3.0, 0.0, 1.0},
            {3.0, 1.0, 1.0},
            {2.0, 1.0, 1.0},
            {4.0, 0.0, 0.0},
            {5.0, 0.0, 0.0},
            {4.0, 1.0, 0.0},
            {4.0, 0.0, 1.0},
            {5.0, 0.0, 1.0},
            {4.0, 1.0, 1.0},
            {6.0, 0.0, 0.0},
            {7.0, 0.0, 0.0},
            {7.0, 1.0, 0.0},
            {6.0, 1.0, 0.0},
            {6.5, 0.5, 1.0}};
    const std::vector<mesh::LinearCell3D>
        cells{
            {mesh::GlobalEntityId{101U},
             mesh::LinearCellType3D::tetrahedron,
             {mesh::LocalIndex{0U},
              mesh::LocalIndex{1U},
              mesh::LocalIndex{2U},
              mesh::LocalIndex{3U}}},
            {mesh::GlobalEntityId{102U},
             mesh::LinearCellType3D::hexahedron,
             {mesh::LocalIndex{4U},
              mesh::LocalIndex{5U},
              mesh::LocalIndex{6U},
              mesh::LocalIndex{7U},
              mesh::LocalIndex{8U},
              mesh::LocalIndex{9U},
              mesh::LocalIndex{10U},
              mesh::LocalIndex{11U}}},
            {mesh::GlobalEntityId{103U},
             mesh::LinearCellType3D::wedge,
             {mesh::LocalIndex{12U},
              mesh::LocalIndex{13U},
              mesh::LocalIndex{14U},
              mesh::LocalIndex{15U},
              mesh::LocalIndex{16U},
              mesh::LocalIndex{17U}}},
            {mesh::GlobalEntityId{104U},
             mesh::LinearCellType3D::pyramid,
             {mesh::LocalIndex{18U},
              mesh::LocalIndex{19U},
              mesh::LocalIndex{20U},
              mesh::LocalIndex{21U},
              mesh::LocalIndex{22U}}}};
    const auto fixture =
        mesh::make_linear_mesh_3d(
            vertex_ids,
            coordinates,
            cells);

    DM source_dm = nullptr;
    std::vector<mesh_petsc::DMPlexPointIdentity>
        source_identities;
    require_petsc(
        mesh_petsc::create_root_dmplex_topology(
            PETSC_COMM_WORLD,
            0,
            mpi_rank == 0
                ? &fixture.topology
                : nullptr,
            &source_dm,
            &source_identities),
        "create rooted mixed linear 3D cell families DMPlex");

    PetscInt source_dimension = -1;
    require_petsc(
        DMGetDimension(
            source_dm,
            &source_dimension),
        "mixed linear 3D cell families source dimension");
    require(
        source_dimension == 3,
        "mixed linear 3D cell families source DMPlex is 3D");

    PetscPartitioner partitioner = nullptr;
    require_petsc(
        DMPlexGetPartitioner(
            source_dm,
            &partitioner),
        "mixed linear 3D cell families DMPlexGetPartitioner");
    require_petsc(
        PetscPartitionerSetType(
            partitioner,
            PETSCPARTITIONERSIMPLE),
        "mixed linear 3D cell families simple partitioner");

    PetscSF migration_sf = nullptr;
    DM distributed_dm = nullptr;
    require_petsc(
        DMPlexDistribute(
            source_dm,
            0,
            &migration_sf,
            &distributed_dm),
        "distribute mixed linear 3D cell families DMPlex");
    require(
        migration_sf != nullptr &&
            distributed_dm != nullptr,
        "mixed linear 3D cell families distribution outputs");

    PetscInt distributed_dimension = -1;
    require_petsc(
        DMGetDimension(
            distributed_dm,
            &distributed_dimension),
        "mixed linear 3D cell families distributed dimension");
    require(
        distributed_dimension == 3,
        "distributed mixed linear 3D cell families DMPlex is 3D");

    PetscInt cell_start = 0;
    PetscInt cell_end = 0;
    require_petsc(
        DMPlexGetHeightStratum(
            distributed_dm,
            0,
            &cell_start,
            &cell_end),
        "mixed linear 3D cell families distributed cell stratum");
    const int local_cells =
        static_cast<int>(
            cell_end - cell_start);
    int global_cells = 0;
    require(
        MPI_Allreduce(
            &local_cells,
            &global_cells,
            1,
            MPI_INT,
            MPI_SUM,
            PETSC_COMM_WORLD) ==
            MPI_SUCCESS,
        "mixed linear 3D cell families global cell count");
    require(
        global_cells == 4,
        "mixed linear 3D cell families distribution preserves all cells");

    PetscSF overlap_sf = nullptr;
    DM overlap_dm = nullptr;
    require_petsc(
        DMPlexDistributeOverlap(
            distributed_dm,
            1,
            &overlap_sf,
            &overlap_dm),
        "mixed linear 3D cell families depth-1 overlap");
    require(
        overlap_sf != nullptr &&
            overlap_dm != nullptr,
        "mixed linear 3D cell families overlap outputs");

    require_petsc(
        PetscSFDestroy(&overlap_sf),
        "destroy mixed linear 3D cell families overlap SF");
    require_petsc(
        DMDestroy(&overlap_dm),
        "destroy mixed linear 3D cell families overlap DM");
    require_petsc(
        PetscSFDestroy(&migration_sf),
        "destroy mixed linear 3D cell families migration SF");
    require_petsc(
        DMDestroy(&distributed_dm),
        "destroy mixed linear 3D cell families distributed DM");
    require_petsc(
        DMDestroy(&source_dm),
        "destroy mixed linear 3D cell families source DM");
}

void run_two_rank_test() {
    int mpi_rank = -1;
    int mpi_size = -1;
    require(
        MPI_Comm_rank(PETSC_COMM_WORLD, &mpi_rank) == MPI_SUCCESS,
        "MPI_Comm_rank");
    require(
        MPI_Comm_size(PETSC_COMM_WORLD, &mpi_size) == MPI_SUCCESS,
        "MPI_Comm_size");
    require(mpi_size == 2, "integration test requires exactly two MPI ranks");
    require(mpi_rank == 0 || mpi_rank == 1, "unexpected MPI rank");

    const auto rank = static_cast<std::uint32_t>(mpi_rank);
    verify_linear_3d_cell_families_dmplex(
        mpi_rank,
        mpi_size);
    const auto topology = two_rank_topology(rank);
    const auto partition = two_rank_partition(topology, rank);
    const auto plan = mesh::SharedEntityPlan::create(
        partition, shared_links());

    const auto layout = mesh::DofLayout::create(
        topology,
        {
            {"cell.primary", mesh::EntityKind::cell, 2U},
            {"vertex.aux", mesh::EntityKind::vertex, 1U},
            {"cell.secondary", mesh::EntityKind::cell, 1U},
            {"face.trace", mesh::EntityKind::face, 1U},
        });
    const auto numbering =
        mesh::DofNumberingSnapshot::create_local(
            layout, partition,
            global_entity_numbering(partition));

    verify_structural_column_pattern_sorting_fixture(
        mpi_rank);
    verify_serial_dmplex_topology();
    verify_serial_dmplex_processed_grdecl();
    verify_dmplex_distribute_overlap_identity();
    verify_processed_grdecl_3d_dmplex_distribute_overlap();
    verify_gmsh_import_through_dmplex_chain();
    verify_section(layout, numbering, mpi_rank);
    verify_sf(partition, plan, mpi_rank);
    verify_global_section_and_section_sf(
        layout, numbering, partition, plan);

    require(
        MPI_Barrier(PETSC_COMM_WORLD) == MPI_SUCCESS,
        "MPI_Barrier");
}

} // namespace

int main(int argc, char** argv) {
    const PetscErrorCode initialize_error =
        PetscInitialize(&argc, &argv, nullptr, nullptr);
    if (initialize_error != PETSC_SUCCESS) {
        return static_cast<int>(initialize_error);
    }

    int result = 0;
    try {
        int mpi_size = -1;
        require(
            MPI_Comm_size(
                PETSC_COMM_WORLD,
                &mpi_size) == MPI_SUCCESS,
            "main MPI_Comm_size");

        std::string_view pass_label;
        if (mpi_size == 2) {
            run_two_rank_test();
            pass_label =
                "mesh.petsc.synthetic_2rank";
        } else if (mpi_size == 3) {
            run_three_rank_sparsity_test();
            pass_label =
                "mesh.petsc.sparsity_owner_targeted_3rank";
        } else {
            throw std::runtime_error(
                "PETSc integration test requires exactly two or three MPI ranks");
        }

        int rank = -1;
        MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
        if (rank == 0) {
            std::cout << "[PASS] "
                      << pass_label
                      << '\n';
        }
    } catch (const std::exception& error) {
        int rank = -1;
        MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
        std::cerr << "[FAIL][rank " << rank << "] "
                  << error.what() << '\n';
        result = 1;
    }

    const PetscErrorCode finalize_error = PetscFinalize();
    if (finalize_error != PETSC_SUCCESS && result == 0) {
        result = static_cast<int>(finalize_error);
    }
    return result;
}
