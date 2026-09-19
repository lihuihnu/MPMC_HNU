#include <mpmc/mesh/cartesian_2d.hpp>
#include <mpmc/mesh/dof_layout.hpp>
#include <mpmc/mesh/geometry_2d.hpp>
#include <mpmc/mesh/dof_numbering.hpp>
#include <mpmc/mesh/partition_snapshot.hpp>
#include <mpmc/mesh/shared_entity_plan.hpp>
#include <mpmc/mesh/topology.hpp>
#include <mpmc/mesh_petsc/adapter.hpp>

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

    verify_serial_dmplex_topology();
    verify_dmplex_distribute_overlap_identity();
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
        run_two_rank_test();
        int rank = -1;
        MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
        if (rank == 0) {
            std::cout << "[PASS] mesh.petsc.synthetic_2rank\n";
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
