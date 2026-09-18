#include <mpmc/mesh/dof_layout.hpp>
#include <mpmc/mesh/dof_numbering.hpp>
#include <mpmc/mesh/partition_snapshot.hpp>
#include <mpmc/mesh/shared_entity_plan.hpp>
#include <mpmc/mesh/topology.hpp>
#include <mpmc/mesh_petsc/adapter.hpp>

#include <petscsys.h>

#include <array>
#include <cstdint>
#include <cstring>
#include <exception>
#include <iostream>
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

    verify_section(layout, numbering, mpi_rank);
    verify_sf(partition, plan, mpi_rank);

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
