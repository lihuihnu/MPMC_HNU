#include <mpmc/flow_discretization_petsc/variable_cardinality_natural_variable_numbering.hpp>

#include <petscsection.h>
#include <petscvec.h>

#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <optional>
#include <stdexcept>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

namespace {

namespace fdp =
    mpmc::flow_discretization_petsc;
namespace mesh = mpmc::mesh;

void require_variable_collective(
    bool condition,
    std::string_view message) {
    int local = condition ? 1 : 0;
    int global = 0;
    if (MPI_Allreduce(
            &local,
            &global,
            1,
            MPI_INT,
            MPI_MIN,
            PETSC_COMM_WORLD) !=
        MPI_SUCCESS) {
        throw std::runtime_error(
            "MPI_Allreduce failed in variable-cardinality assertion");
    }
    if (global == 0) {
        throw std::runtime_error(
            std::string{message});
    }
}

mesh::Topology variable_topology(
    int rank) {
    mesh::Topology::EntityIds ids;
    ids.cells =
        rank == 0
            ? std::vector<mesh::GlobalEntityId>{
                  mesh::GlobalEntityId{UINT64_C(10)},
                  mesh::GlobalEntityId{UINT64_C(20)},
                  mesh::GlobalEntityId{UINT64_C(30)}}
            : std::vector<mesh::GlobalEntityId>{
                  mesh::GlobalEntityId{UINT64_C(20)},
                  mesh::GlobalEntityId{UINT64_C(30)},
                  mesh::GlobalEntityId{UINT64_C(10)}};
    return {
        std::move(ids),
        {}};
}

mesh::PartitionSnapshot
variable_partition(
    int rank) {
    auto topology =
        variable_topology(rank);
    mesh::EntityOwnerRanks owners;
    owners.cells =
        rank == 0
            ? std::vector<mesh::PartitionRank>{
                  mesh::PartitionRank{0U},
                  mesh::PartitionRank{1U},
                  mesh::PartitionRank{1U}}
            : std::vector<mesh::PartitionRank>{
                  mesh::PartitionRank{1U},
                  mesh::PartitionRank{1U},
                  mesh::PartitionRank{0U}};
    return mesh::PartitionSnapshot::create(
        topology,
        mesh::PartitionRank{
            static_cast<
                mesh::PartitionRank::value_type>(
                    rank)},
        2U,
        std::move(owners));
}

std::vector<std::size_t>
variable_phase_counts(
    int rank) {
    return rank == 0
        ? std::vector<std::size_t>{
              1U, 2U, 3U}
        : std::vector<std::size_t>{
              2U, 3U, 1U};
}

PetscInt expected_global_start(
    std::uint64_t stable_cell) {
    switch (stable_cell) {
    case UINT64_C(10):
        return 0;
    case UINT64_C(20):
        return 4;
    case UINT64_C(30):
        return 11;
    default:
        throw std::invalid_argument(
            "unknown variable-cardinality stable cell");
    }
}

std::size_t expected_phase_count(
    std::uint64_t stable_cell) {
    switch (stable_cell) {
    case UINT64_C(10):
        return 1U;
    case UINT64_C(20):
        return 2U;
    case UINT64_C(30):
        return 3U;
    default:
        throw std::invalid_argument(
            "unknown variable-cardinality stable cell");
    }
}

std::size_t expected_width(
    std::uint64_t stable_cell) {
    return expected_phase_count(
               stable_cell) *
            3U +
        1U;
}

} // namespace

void variable_cardinality_natural_variable_numbering_test() {
    int rank = -1;
    int size = -1;
    if (MPI_Comm_rank(
            PETSC_COMM_WORLD,
            &rank) != MPI_SUCCESS ||
        MPI_Comm_size(
            PETSC_COMM_WORLD,
            &size) != MPI_SUCCESS) {
        throw std::runtime_error(
            "MPI rank/size query failed in variable-cardinality test");
    }
    require_variable_collective(
        size == 2,
        "variable-cardinality test requires exactly two ranks");

    const auto partition =
        variable_partition(rank);
    const auto phase_counts =
        variable_phase_counts(rank);

    std::optional<
        fdp::VariableCardinalityNaturalVariableNumbering3D>
        numbering;
    PetscErrorCode error =
        fdp::
            make_variable_cardinality_natural_variable_numbering_3d(
                PETSC_COMM_WORLD,
                partition,
                3U,
                phase_counts,
                &numbering);
    require_variable_collective(
        error == PETSC_SUCCESS &&
            numbering.has_value(),
        "failed to build variable-cardinality natural-variable numbering");

    require_variable_collective(
        numbering->component_count() == 3U &&
            numbering->local_cell_count() == 3U &&
            numbering->local_packed_scalar_count() ==
                21U &&
            numbering->petsc_global_scalar_count() ==
                21,
        "variable-cardinality global size does not equal 4+7+10");

    const PetscInt expected_start =
        rank == 0 ? 0 : 4;
    const PetscInt expected_end =
        rank == 0 ? 4 : 21;
    require_variable_collective(
        numbering->petsc_owned_scalar_start() ==
                expected_start &&
            numbering->petsc_owned_scalar_end() ==
                expected_end &&
            numbering->petsc_local_owned_scalar_count() ==
                expected_end -
                    expected_start,
        "variable-cardinality PETSc ownership range mismatch");

    const std::array<std::uint64_t, 3>
        local_stable =
            rank == 0
                ? std::array<std::uint64_t, 3>{
                      UINT64_C(10),
                      UINT64_C(20),
                      UINT64_C(30)}
                : std::array<std::uint64_t, 3>{
                      UINT64_C(20),
                      UINT64_C(30),
                      UINT64_C(10)};
    const std::array<std::size_t, 3>
        expected_offsets =
            rank == 0
                ? std::array<std::size_t, 3>{
                      0U, 4U, 11U}
                : std::array<std::size_t, 3>{
                      0U, 7U, 17U};

    bool local_cell_metadata_ok = true;
    bool local_scalar_mapping_ok = true;
    for (std::size_t local = 0U;
         local < local_stable.size();
         ++local) {
        const auto cell =
            mesh::LocalIndex{
                static_cast<
                    mesh::LocalIndex::value_type>(
                        local)};
        const auto& record =
            numbering->cell(cell);
        const std::uint64_t stable =
            local_stable[local];

        local_cell_metadata_ok =
            local_cell_metadata_ok &&
            record.cell_global.value() ==
                stable &&
            record.phase_count ==
                expected_phase_count(stable) &&
            record.scalar_count ==
                expected_width(stable) &&
            record.local_scalar_offset ==
                expected_offsets[local] &&
            record.petsc_global_scalar_start ==
                expected_global_start(stable);

        for (std::size_t slot = 0U;
             slot < record.scalar_count;
             ++slot) {
            const PetscInt expected =
                expected_global_start(
                    stable) +
                static_cast<PetscInt>(slot);
            local_scalar_mapping_ok =
                local_scalar_mapping_ok &&
                numbering->petsc_global_scalar(
                    cell,
                    slot) ==
                    expected &&
                numbering
                        ->local_packed_scalar_global_index(
                            record.local_scalar_offset +
                            slot) ==
                    expected;
        }
    }
    require_variable_collective(
        local_cell_metadata_ok,
        "variable-cardinality cell DoF metadata mismatch");
    require_variable_collective(
        local_scalar_mapping_ok,
        "owner/ghost scalar numbering disagrees");

    require_variable_collective(
        numbering->is_owned_cell(
            partition.local_index(
                mesh::EntityKind::cell,
                mesh::GlobalEntityId{
                    rank == 0
                        ? UINT64_C(10)
                        : UINT64_C(20)})),
        "owned variable-cardinality cell was not marked owned");
    require_variable_collective(
        numbering->is_ghost_cell(
            partition.local_index(
                mesh::EntityKind::cell,
                mesh::GlobalEntityId{
                    rank == 0
                        ? UINT64_C(20)
                        : UINT64_C(10)})),
        "ghost variable-cardinality cell was not marked ghost");

    PetscSection section = nullptr;
    error =
        fdp::
            create_variable_cardinality_natural_variable_section_3d(
                PETSC_COMM_WORLD,
                *numbering,
                &section);
    require_variable_collective(
        error == PETSC_SUCCESS &&
            section != nullptr,
        "failed to create variable-cardinality PetscSection");

    PetscInt section_storage = -1;
    require_variable_collective(
        PetscSectionGetStorageSize(
            section,
            &section_storage) ==
                PETSC_SUCCESS &&
            section_storage == 21,
        "ragged PetscSection storage does not equal local 4+7+10");

    for (PetscInt point = 0;
         point < 3;
         ++point) {
        PetscInt dof = -1;
        PetscInt offset = -1;
        require_variable_collective(
            PetscSectionGetDof(
                section,
                point,
                &dof) ==
                    PETSC_SUCCESS &&
                PetscSectionGetOffset(
                    section,
                    point,
                    &offset) ==
                    PETSC_SUCCESS,
            "failed to inspect variable-cardinality PetscSection");
        const auto& record =
            numbering->cell(
                mesh::LocalIndex{
                    static_cast<
                        mesh::LocalIndex::value_type>(
                            point)});
        require_variable_collective(
            dof ==
                    static_cast<PetscInt>(
                        record.scalar_count) &&
                offset ==
                    static_cast<PetscInt>(
                        record.local_scalar_offset),
            "PetscSection ragged dof/offset mismatch");
    }

    Vec vector = nullptr;
    error =
        fdp::
            create_variable_cardinality_natural_variable_vec_3d(
                PETSC_COMM_WORLD,
                *numbering,
                &vector);
    require_variable_collective(
        error == PETSC_SUCCESS &&
            vector != nullptr,
        "failed to create mixed-cardinality PETSc Vec");

    bool local_insert_ok = true;
    for (const auto& record :
         numbering->cells()) {
        if (record.owner_rank !=
            numbering->local_rank()) {
            continue;
        }
        for (std::size_t slot = 0U;
             slot < record.scalar_count;
             ++slot) {
            const PetscInt index =
                record.petsc_global_scalar_start +
                static_cast<PetscInt>(
                    slot);
            const PetscScalar value =
                static_cast<PetscScalar>(
                    index + 1);
            local_insert_ok =
                local_insert_ok &&
                VecSetValues(
                    vector,
                    1,
                    &index,
                    &value,
                    INSERT_VALUES) ==
                    PETSC_SUCCESS;
        }
    }
    require_variable_collective(
        local_insert_ok,
        "failed to insert mixed-cardinality owned scalar");
    require_variable_collective(
        VecAssemblyBegin(vector) ==
                PETSC_SUCCESS &&
            VecAssemblyEnd(vector) ==
                PETSC_SUCCESS,
        "failed to assemble mixed-cardinality PETSc Vec");

    const PetscScalar* local_values =
        nullptr;
    require_variable_collective(
        VecGetArrayRead(
            vector,
            &local_values) ==
            PETSC_SUCCESS,
        "failed to access mixed-cardinality PETSc Vec");
    bool local_values_ok = true;
    for (PetscInt local = 0;
         local <
             numbering
                 ->petsc_local_owned_scalar_count();
         ++local) {
        const PetscInt global =
            numbering
                ->petsc_owned_scalar_start() +
            local;
        const double value =
            static_cast<double>(
                PetscRealPart(
                    local_values[local]));
        local_values_ok =
            local_values_ok &&
            value ==
                static_cast<double>(
                    global + 1);
    }
    require_variable_collective(
        local_values_ok,
        "mixed-cardinality PETSc Vec ownership/value mismatch");
    require_variable_collective(
        VecRestoreArrayRead(
            vector,
            &local_values) ==
            PETSC_SUCCESS,
        "failed to restore mixed-cardinality PETSc Vec");

    require_variable_collective(
        VecDestroy(&vector) ==
                PETSC_SUCCESS &&
            PetscSectionDestroy(&section) ==
                PETSC_SUCCESS,
        "variable-cardinality PETSc fixture cleanup failed");

    auto inconsistent_counts =
        phase_counts;
    if (rank == 1) {
        inconsistent_counts[2] =
            2U;
    }
    std::optional<
        fdp::VariableCardinalityNaturalVariableNumbering3D>
        rejected;
    error =
        fdp::
            make_variable_cardinality_natural_variable_numbering_3d(
                PETSC_COMM_WORLD,
                partition,
                3U,
                inconsistent_counts,
                &rejected);
    require_variable_collective(
        error == PETSC_ERR_ARG_INCOMP &&
            !rejected.has_value(),
        "owner/ghost phase-cardinality mismatch was not rejected collectively");
}
