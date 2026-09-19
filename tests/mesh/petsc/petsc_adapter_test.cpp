#include <mpmc/mesh/cartesian_2d.hpp>
#include <mpmc/mesh/dof_layout.hpp>
#include <mpmc/mesh/dof_numbering.hpp>
#include <mpmc/mesh/partition_snapshot.hpp>
#include <mpmc/mesh/shared_entity_plan.hpp>
#include <mpmc/mesh/topology.hpp>
#include <mpmc/mesh_petsc/adapter.hpp>

#include <petscsys.h>

#include <algorithm>
#include <array>
#include <cstdint>
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
