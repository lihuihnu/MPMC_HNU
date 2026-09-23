#include <mpmc/flow_discretization_petsc/variable_cardinality_snes_assembly.hpp>

#include <petscmat.h>
#include <petscvec.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <optional>
#include <span>
#include <stdexcept>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

namespace {

namespace dp = mpmc::discretization_petsc;
namespace fdp = mpmc::flow_discretization_petsc;
namespace mesh = mpmc::mesh;

void require_mixed_collective(
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
            "MPI_Allreduce failed in mixed-cardinality assertion");
    }
    if (global == 0) {
        throw std::runtime_error(
            std::string{message});
    }
}

void near_mixed(
    double actual,
    double expected,
    double relative = 2.0e-9,
    double absolute = 2.0e-10) {
    if (!std::isfinite(actual) ||
        !std::isfinite(expected)) {
        throw std::runtime_error(
            "non-finite mixed-cardinality comparison");
    }
    const double scale =
        std::max(
            std::abs(actual),
            std::abs(expected));
    if (std::abs(actual - expected) >
        absolute + relative * scale) {
        throw std::runtime_error(
            "mixed-cardinality numeric mismatch");
    }
}

mesh::Topology mixed_topology(
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
mixed_partition(
    int rank) {
    auto topology =
        mixed_topology(rank);
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
mixed_phase_counts(
    int rank) {
    return rank == 0
        ? std::vector<std::size_t>{
              1U, 2U, 3U}
        : std::vector<std::size_t>{
              2U, 3U, 1U};
}

dp::PetscMpiAijSymbolicPreallocation3D
mixed_cell_bridge(
    int rank) {
    if (rank == 0) {
        return {
            mesh::PartitionRank{0U},
            2U,
            0,
            1,
            3,
            {mesh::LocalIndex{0U}},
            {mesh::GlobalEntityId{UINT64_C(10)}},
            {0},
            {1},
            {1},
            {0, 1, 2}};
    }

    return {
        mesh::PartitionRank{1U},
        2U,
        1,
        3,
        3,
        {
            mesh::LocalIndex{0U},
            mesh::LocalIndex{1U}},
        {
            mesh::GlobalEntityId{UINT64_C(20)},
            mesh::GlobalEntityId{UINT64_C(30)}},
        {1, 2},
        {2, 2},
        {1, 0},
        {1, 2, 0}};
}

dp::OwnedCellStructuralColumnPatternSnapshot3D
mixed_cell_pattern(
    int rank) {
    if (rank == 0) {
        return {
            mesh::PartitionRank{0U},
            2U,
            3U,
            0,
            1,
            3,
            {
                {
                    mesh::LocalIndex{0U},
                    mesh::GlobalEntityId{
                        UINT64_C(10)},
                    0,
                    0U,
                    1U,
                    0U,
                    1U}},
            {0},
            {1}};
    }

    return {
        mesh::PartitionRank{1U},
        2U,
        3U,
        1,
        3,
        3,
        {
            {
                mesh::LocalIndex{0U},
                mesh::GlobalEntityId{
                    UINT64_C(20)},
                1,
                0U,
                2U,
                0U,
                1U},
            {
                mesh::LocalIndex{1U},
                mesh::GlobalEntityId{
                    UINT64_C(30)},
                2,
                2U,
                2U,
                1U,
                0U}},
        {1, 2, 1, 2},
        {0}};
}

std::vector<double>
target_for(
    std::uint64_t stable_cell) {
    switch (stable_cell) {
    case UINT64_C(10):
        return {
            10.0,
            8.0,
            0.20,
            0.30};
    case UINT64_C(20):
        return {
            11.0,
            8.5,
            0.30,
            0.20,
            0.25,
            0.25,
            0.30};
    case UINT64_C(30):
        return {
            12.0,
            9.0,
            0.25,
            0.30,
            0.20,
            0.25,
            0.20,
            0.30,
            0.15,
            0.25};
    default:
        throw std::invalid_argument(
            "unknown mixed-cardinality stable cell");
    }
}

std::vector<std::uint64_t>
neighbours_for(
    std::uint64_t stable_cell) {
    switch (stable_cell) {
    case UINT64_C(10):
        return {UINT64_C(20)};
    case UINT64_C(20):
        return {
            UINT64_C(10),
            UINT64_C(30)};
    case UINT64_C(30):
        return {UINT64_C(20)};
    default:
        throw std::invalid_argument(
            "unknown mixed-cardinality stable cell");
    }
}

mesh::LocalIndex
local_cell_for_stable(
    const fdp::VariableCardinalityNaturalVariableNumbering3D&
        numbering,
    std::uint64_t stable_cell) {
    for (const auto& cell :
         numbering.cells()) {
        if (cell.cell_global.value() ==
            stable_cell) {
            return cell.cell;
        }
    }
    throw std::invalid_argument(
        "stable cell missing from local mixed-cardinality overlap");
}

struct ControlledMixedClosure {
    static constexpr double coupling =
        1.0e-3;
    PetscInt calls{};
};

PetscErrorCode
controlled_mixed_evaluator(
    mesh::LocalIndex row_cell,
    std::size_t phase_count,
    std::span<const double>
        local_packed_state,
    const fdp::
        VariableCardinalityNaturalVariableNumbering3D&
            numbering,
    void* raw_context,
    fdp::
        VariableCardinalityNaturalVariableCellAssembly3D*
            output,
    fdp::NaturalVariableSnesEvaluationStatus3D*
        status) {
    if (raw_context == nullptr ||
        output == nullptr ||
        status == nullptr) {
        return PETSC_ERR_ARG_NULL;
    }
    auto* context =
        static_cast<ControlledMixedClosure*>(
            raw_context);
    ++context->calls;

    const auto& row =
        numbering.cell(
            row_cell);
    if (phase_count !=
            row.phase_count ||
        local_packed_state.size() !=
            numbering
                .local_packed_scalar_count()) {
        return PETSC_ERR_ARG_INCOMP;
    }

    const auto target =
        target_for(
            row.cell_global.value());
    if (target.size() !=
        row.scalar_count) {
        return PETSC_ERR_ARG_INCOMP;
    }

    const auto state =
        local_packed_state.subspan(
            row.local_scalar_offset,
            row.scalar_count);
    for (const double value :
         state) {
        if (!std::isfinite(value) ||
            !(value > 0.0)) {
            *status =
                fdp::
                    NaturalVariableSnesEvaluationStatus3D::
                        domain_error;
            output->residual.clear();
            output->jacobian_blocks.clear();
            return PETSC_SUCCESS;
        }
    }

    output->residual.resize(
        row.scalar_count);
    for (std::size_t slot = 0U;
         slot < row.scalar_count;
         ++slot) {
        output->residual[slot] =
            std::log(
                state[slot] /
                target[slot]);
    }

    auto self =
        fdp::
            VariableCardinalityNaturalVariableJacobianBlock3D{};
    self.column_cell =
        row_cell;
    self.values_row_major.assign(
        row.scalar_count *
            row.scalar_count,
        0.0);
    for (std::size_t slot = 0U;
         slot < row.scalar_count;
         ++slot) {
        self.values_row_major[
            slot * row.scalar_count +
            slot] =
            1.0 / state[slot];
    }
    output->jacobian_blocks.clear();
    output->jacobian_blocks.push_back(
        std::move(self));

    for (const std::uint64_t
             neighbour_stable :
         neighbours_for(
             row.cell_global.value())) {
        const auto neighbour_cell =
            local_cell_for_stable(
                numbering,
                neighbour_stable);
        const auto& neighbour =
            numbering.cell(
                neighbour_cell);
        const auto neighbour_target =
            target_for(
                neighbour_stable);
        const double neighbour_value =
            local_packed_state[
                neighbour.local_scalar_offset];

        for (double& value :
             output->residual) {
            value +=
                ControlledMixedClosure::coupling *
                (neighbour_value -
                 neighbour_target[0]);
        }

        fdp::
            VariableCardinalityNaturalVariableJacobianBlock3D
            block;
        block.column_cell =
            neighbour_cell;
        block.values_row_major.assign(
            row.scalar_count *
                neighbour.scalar_count,
            0.0);
        for (std::size_t row_slot = 0U;
             row_slot <
                row.scalar_count;
             ++row_slot) {
            block.values_row_major[
                row_slot *
                    neighbour.scalar_count] =
                ControlledMixedClosure::
                    coupling;
        }
        output->jacobian_blocks.push_back(
            std::move(block));
    }

    *status =
        fdp::
            NaturalVariableSnesEvaluationStatus3D::
                success;
    return PETSC_SUCCESS;
}

void insert_initial_state(
    Vec state,
    const fdp::
        VariableCardinalityNaturalVariableNumbering3D&
            numbering) {
    bool local_ok = true;
    for (const auto& cell :
         numbering.cells()) {
        if (cell.owner_rank !=
            numbering.local_rank()) {
            continue;
        }
        const auto target =
            target_for(
                cell.cell_global.value());
        for (std::size_t slot = 0U;
             slot <
                cell.scalar_count;
             ++slot) {
            const PetscInt index =
                cell.petsc_global_scalar_start +
                static_cast<PetscInt>(
                    slot);
            const PetscScalar value =
                static_cast<PetscScalar>(
                    1.2 *
                    target[slot]);
            local_ok =
                local_ok &&
                VecSetValues(
                    state,
                    1,
                    &index,
                    &value,
                    INSERT_VALUES) ==
                    PETSC_SUCCESS;
        }
    }
    require_mixed_collective(
        local_ok,
        "failed to insert mixed-cardinality initial state");
    require_mixed_collective(
        VecAssemblyBegin(state) ==
                PETSC_SUCCESS &&
            VecAssemblyEnd(state) ==
                PETSC_SUCCESS,
        "failed to assemble mixed-cardinality initial state");
}

} // namespace

void variable_cardinality_snes_assembly_test() {
    int rank = -1;
    int size = -1;
    if (MPI_Comm_rank(
            PETSC_COMM_WORLD,
            &rank) != MPI_SUCCESS ||
        MPI_Comm_size(
            PETSC_COMM_WORLD,
            &size) != MPI_SUCCESS) {
        throw std::runtime_error(
            "failed to query MPI rank/size for mixed-cardinality SNES test");
    }
    require_mixed_collective(
        size == 2,
        "mixed-cardinality SNES test requires two ranks");

    const auto partition =
        mixed_partition(
            rank);
    const auto phase_counts =
        mixed_phase_counts(
            rank);
    std::optional<
        fdp::
            VariableCardinalityNaturalVariableNumbering3D>
        numbering;
    PetscErrorCode error =
        fdp::
            make_variable_cardinality_natural_variable_numbering_3d(
                PETSC_COMM_WORLD,
                partition,
                3U,
                phase_counts,
                &numbering);
    require_mixed_collective(
        error == PETSC_SUCCESS &&
            numbering.has_value(),
        "failed to build mixed-cardinality numbering");

    const auto bridge =
        mixed_cell_bridge(
            rank);
    const auto pattern =
        mixed_cell_pattern(
            rank);

    ControlledMixedClosure closure;
    std::optional<
        fdp::
            VariableCardinalityNaturalVariableSnesAssemblyContext3D>
        context;
    error =
        fdp::
            VariableCardinalityNaturalVariableSnesAssemblyContext3D::
                create(
                    PETSC_COMM_WORLD,
                    *numbering,
                    bridge,
                    pattern,
                    {
                        &controlled_mixed_evaluator,
                        &closure},
                    &context);
    require_mixed_collective(
        error == PETSC_SUCCESS &&
            context.has_value(),
        "failed to create mixed-cardinality SNES assembly context");

    Mat structure = nullptr;
    error =
        context->create_jacobian_structure(
            &structure);
    require_mixed_collective(
        error == PETSC_SUCCESS &&
            structure != nullptr,
        "failed to create q-ragged MPIAIJ structure");

    bool local_structure_ok = true;
    if (rank == 0) {
        PetscInt count = -1;
        const PetscInt row = 0;
        error =
            MatGetRow(
                structure,
                row,
                &count,
                nullptr,
                nullptr);
        local_structure_ok =
            error == PETSC_SUCCESS &&
            count == 11;
        if (error == PETSC_SUCCESS) {
            local_structure_ok =
                local_structure_ok &&
                MatRestoreRow(
                    structure,
                    row,
                    &count,
                    nullptr,
                    nullptr) ==
                    PETSC_SUCCESS;
        }
    } else {
        PetscInt count20 = -1;
        PetscInt count30 = -1;
        const PetscInt row20 = 4;
        const PetscInt row30 = 11;
        error =
            MatGetRow(
                structure,
                row20,
                &count20,
                nullptr,
                nullptr);
        if (error == PETSC_SUCCESS) {
            local_structure_ok =
                count20 == 21 &&
                MatRestoreRow(
                    structure,
                    row20,
                    &count20,
                    nullptr,
                    nullptr) ==
                    PETSC_SUCCESS;
        } else {
            local_structure_ok = false;
        }
        error =
            MatGetRow(
                structure,
                row30,
                &count30,
                nullptr,
                nullptr);
        if (error == PETSC_SUCCESS) {
            local_structure_ok =
                local_structure_ok &&
                count30 == 17 &&
                MatRestoreRow(
                    structure,
                    row30,
                    &count30,
                    nullptr,
                    nullptr) ==
                    PETSC_SUCCESS;
        } else {
            local_structure_ok = false;
        }
    }
    require_mixed_collective(
        local_structure_ok,
        "ragged MPIAIJ scalar nnz does not match 4/7/10 block widths");

    Vec initial = nullptr;
    error =
        fdp::
            create_variable_cardinality_natural_variable_vec_3d(
                PETSC_COMM_WORLD,
                *numbering,
                &initial);
    require_mixed_collective(
        error == PETSC_SUCCESS &&
            initial != nullptr,
        "failed to create mixed-cardinality SNES state");
    insert_initial_state(
        initial,
        *numbering);

    Vec row_scaling = nullptr;
    error =
        fdp::
            make_variable_cardinality_initial_row_equilibration_3d(
                PETSC_COMM_WORLD,
                *numbering,
                initial,
                structure,
                context->snes_evaluator(),
                &row_scaling);
    require_mixed_collective(
        error == PETSC_SUCCESS &&
            row_scaling != nullptr,
        "failed to construct q-ragged frozen row equilibration");

    Vec solution = nullptr;
    std::optional<
        fdp::
            VariableCardinalityNaturalVariableSnesSolveReport3D>
        report;
    error =
        fdp::
            solve_variable_cardinality_natural_variable_snes_3d(
                PETSC_COMM_WORLD,
                *numbering,
                initial,
                structure,
                context->snes_evaluator(),
                &solution,
                &report,
                row_scaling);
    require_mixed_collective(
        error == PETSC_SUCCESS &&
            solution != nullptr &&
            report.has_value(),
        "mixed-cardinality PETSc SNES solve failed");

    require_mixed_collective(
        static_cast<int>(
            report->converged_reason) >
                0 &&
            report->snes_type ==
                SNESNEWTONLS &&
            report->line_search_type ==
                SNESLINESEARCHBT &&
            report->ksp_type ==
                KSPGMRES &&
            report->pc_type ==
                PCASM &&
            report->asm_overlap == 1 &&
            report->function_evaluations >
                0 &&
            report->jacobian_evaluations >
                0 &&
            report->final_function_l2_norm <
                1.0e-8,
        "mixed-cardinality SNES solver policy/convergence mismatch");

    bool local_solution_ok = true;
    for (const auto& entry :
         report->locally_owned_solution) {
        const auto target =
            target_for(
                entry.cell_global.value());
        local_solution_ok =
            local_solution_ok &&
            entry.phase_count >= 1U &&
            entry.phase_count <= 3U &&
            entry.natural_variable_slot <
                target.size() &&
            std::abs(
                entry.value -
                target[
                    entry.natural_variable_slot]) <=
                2.0e-8 *
                    std::max(
                        1.0,
                        std::abs(
                            target[
                                entry.natural_variable_slot]));
    }
    require_mixed_collective(
        local_solution_ok,
        "mixed-cardinality SNES solution did not recover q-ragged target");
    require_mixed_collective(
        VecDestroy(
            &row_scaling) ==
            PETSC_SUCCESS,
        "mixed-cardinality row scaling cleanup failed");

    std::vector<double>
        local_packed;
    error =
        context->copy_local_packed_state(
            solution,
            &local_packed);
    require_mixed_collective(
        error == PETSC_SUCCESS &&
            local_packed.size() ==
                numbering
                    ->local_packed_scalar_count(),
        "mixed-cardinality PetscSF state copy failed");

    bool local_ghost_ok = true;
    for (const auto& cell :
         numbering->cells()) {
        const auto target =
            target_for(
                cell.cell_global.value());
        for (std::size_t slot = 0U;
             slot <
                cell.scalar_count;
             ++slot) {
            const double actual =
                local_packed[
                    cell.local_scalar_offset +
                    slot];
            local_ghost_ok =
                local_ghost_ok &&
                std::abs(
                    actual -
                    target[slot]) <=
                    2.0e-8 *
                        std::max(
                            1.0,
                            std::abs(
                                target[slot]));
        }
    }
    require_mixed_collective(
        local_ghost_ok,
        "q-ragged PetscSF owner/ghost values disagree after solve");

    Mat evaluated = nullptr;
    error =
        MatDuplicate(
            structure,
            MAT_COPY_VALUES,
            &evaluated);
    if (error == PETSC_SUCCESS) {
        error =
            MatZeroEntries(
                evaluated);
    }

    auto evaluator =
        context->snes_evaluator();
    fdp::NaturalVariableSnesEvaluationStatus3D
        status =
            fdp::
                NaturalVariableSnesEvaluationStatus3D::
                    success;
    if (error == PETSC_SUCCESS) {
        error =
            evaluator.jacobian(
                solution,
                evaluated,
                evaluator.user_context,
                &status);
    }
    require_mixed_collective(
        error == PETSC_SUCCESS &&
            status ==
                fdp::
                    NaturalVariableSnesEvaluationStatus3D::
                        success,
        "failed to evaluate mixed-cardinality rectangular Jacobian");
    require_mixed_collective(
        MatAssemblyBegin(
            evaluated,
            MAT_FINAL_ASSEMBLY) ==
                PETSC_SUCCESS &&
            MatAssemblyEnd(
                evaluated,
                MAT_FINAL_ASSEMBLY) ==
                PETSC_SUCCESS,
        "failed to assemble mixed-cardinality rectangular Jacobian");

    bool local_rectangular_ok = true;
    if (rank == 0) {
        const PetscInt row = 0;
        const PetscInt column = 4;
        PetscScalar value = 0.0;
        local_rectangular_ok =
            MatGetValues(
                evaluated,
                1,
                &row,
                1,
                &column,
                &value) ==
                PETSC_SUCCESS;
        if (local_rectangular_ok) {
            near_mixed(
                static_cast<double>(
                    PetscRealPart(value)),
                ControlledMixedClosure::
                    coupling);
        }
    } else {
        const std::array<
            std::pair<PetscInt, PetscInt>,
            3>
            probes{{
                {4, 0},
                {4, 11},
                {11, 4}}};
        for (const auto& [row, column] :
             probes) {
            PetscScalar value = 0.0;
            local_rectangular_ok =
                local_rectangular_ok &&
                MatGetValues(
                    evaluated,
                    1,
                    &row,
                    1,
                    &column,
                    &value) ==
                    PETSC_SUCCESS;
            if (local_rectangular_ok) {
                near_mixed(
                    static_cast<double>(
                        PetscRealPart(value)),
                    ControlledMixedClosure::
                        coupling);
            }
        }
    }
    require_mixed_collective(
        local_rectangular_ok,
        "cross-cardinality rectangular Jacobian coupling is missing");

    require_mixed_collective(
        MatDestroy(&evaluated) ==
                PETSC_SUCCESS &&
            VecDestroy(&solution) ==
                PETSC_SUCCESS &&
            VecDestroy(&initial) ==
                PETSC_SUCCESS &&
            MatDestroy(&structure) ==
                PETSC_SUCCESS,
        "mixed-cardinality SNES fixture cleanup failed");
}
