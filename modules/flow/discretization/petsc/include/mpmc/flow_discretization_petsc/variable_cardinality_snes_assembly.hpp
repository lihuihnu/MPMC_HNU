#ifndef MPMC_FLOW_DISCRETIZATION_PETSC_VARIABLE_CARDINALITY_SNES_ASSEMBLY_HPP
#define MPMC_FLOW_DISCRETIZATION_PETSC_VARIABLE_CARDINALITY_SNES_ASSEMBLY_HPP

#include <mpmc/discretization_petsc/adapter.hpp>
#include <mpmc/flow_discretization_petsc/natural_variable_snes_solver.hpp>
#include <mpmc/flow_discretization_petsc/variable_cardinality_natural_variable_numbering.hpp>

#include <petscmat.h>
#include <petscsf.h>
#include <petscsnes.h>
#include <petscvec.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <optional>
#include <span>
#include <stdexcept>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

namespace mpmc::flow_discretization_petsc {

inline constexpr std::string_view
    variable_cardinality_snes_assembly_convention =
        "flow_discretization_petsc/variable-cardinality-snes-assembly/v1";

struct VariableCardinalityNaturalVariableJacobianBlock3D {
    mpmc::mesh::LocalIndex column_cell{
        mpmc::mesh::LocalIndex::value_type{0}};
    std::vector<double> values_row_major;
};

struct VariableCardinalityNaturalVariableCellAssembly3D {
    std::vector<double> residual;
    std::vector<
        VariableCardinalityNaturalVariableJacobianBlock3D>
        jacobian_blocks;
};

using VariableCardinalityNaturalVariableCellEvaluator3D =
    PetscErrorCode (*)(
        mpmc::mesh::LocalIndex row_cell,
        std::size_t phase_count,
        std::span<const double> local_packed_state,
        const VariableCardinalityNaturalVariableNumbering3D&
            numbering,
        void* user_context,
        VariableCardinalityNaturalVariableCellAssembly3D*
            output,
        NaturalVariableSnesEvaluationStatus3D* status);

struct VariableCardinalityNaturalVariableCellEvaluatorBinding3D {
    VariableCardinalityNaturalVariableCellEvaluator3D evaluator{};
    void* user_context{};
};

struct VariableCardinalityNaturalVariableSnesSolutionEntry3D {
    PetscInt petsc_global_scalar{-1};
    mpmc::mesh::GlobalEntityId cell_global{
        mpmc::mesh::GlobalEntityId::value_type{0}};
    std::size_t phase_count{};
    std::size_t natural_variable_slot{};
    double value{};
};

struct VariableCardinalityNaturalVariableSnesSolveReport3D {
    SNESConvergedReason converged_reason{
        SNES_CONVERGED_ITERATING};
    PetscInt nonlinear_iterations{};
    PetscInt function_evaluations{};
    PetscInt jacobian_evaluations{};
    PetscInt line_search_prechecks{};
    PetscInt line_search_direction_changes{};
    double final_function_l2_norm{};
    std::string snes_type;
    std::string line_search_type;
    std::string ksp_type;
    std::string pc_type;
    PetscInt asm_overlap{};
    std::vector<
        VariableCardinalityNaturalVariableSnesSolutionEntry3D>
        locally_owned_solution;
};

namespace variable_cardinality_snes_assembly_detail {

[[nodiscard]] inline PetscErrorCode
collective_error(
    MPI_Comm comm,
    PetscErrorCode local_error) {
    int local =
        static_cast<int>(local_error);
    int global = 0;
    if (MPI_Allreduce(
            &local,
            &global,
            1,
            MPI_INT,
            MPI_MAX,
            comm) != MPI_SUCCESS) {
        return PETSC_ERR_MPI;
    }
    return static_cast<PetscErrorCode>(
        global);
}

[[nodiscard]] inline bool
contains_cell_row(
    std::span<const PetscInt> rows,
    PetscInt row) {
    return std::binary_search(
        rows.begin(),
        rows.end(),
        row);
}

[[nodiscard]] inline bool
communicators_are_compatible(
    MPI_Comm expected,
    MPI_Comm actual) {
    int comparison = MPI_UNEQUAL;
    return MPI_Comm_compare(
               expected,
               actual,
               &comparison) ==
               MPI_SUCCESS &&
        (comparison == MPI_IDENT ||
         comparison == MPI_CONGRUENT);
}

} // namespace variable_cardinality_snes_assembly_detail

/// Distributed mixed-cardinality PETSc assembly context.
///
/// The context is numerical/assembly infrastructure only. It does not select
/// phase sets, evaluate EOS/flash, or invent cross-phase identity. A caller
/// supplies one evaluator that receives the q-ragged owned+ghost state and
/// returns an owned row-cell residual plus rectangular q_row x q_col Jacobian
/// blocks. Structural block membership is checked against the existing
/// cell-pair sparsity snapshot.
///
/// State exchange uses PetscSF. State values are never allgathered.
class VariableCardinalityNaturalVariableSnesAssemblyContext3D {
public:
    static constexpr std::string_view convention =
        variable_cardinality_snes_assembly_convention;

    VariableCardinalityNaturalVariableSnesAssemblyContext3D(
        const VariableCardinalityNaturalVariableSnesAssemblyContext3D&) =
        delete;
    VariableCardinalityNaturalVariableSnesAssemblyContext3D& operator=(
        const VariableCardinalityNaturalVariableSnesAssemblyContext3D&) =
        delete;
    VariableCardinalityNaturalVariableSnesAssemblyContext3D& operator=(
        VariableCardinalityNaturalVariableSnesAssemblyContext3D&&) =
        delete;

    VariableCardinalityNaturalVariableSnesAssemblyContext3D(
        VariableCardinalityNaturalVariableSnesAssemblyContext3D&&
            other) noexcept
        : comm_(other.comm_),
          numbering_(other.numbering_),
          cell_bridge_(other.cell_bridge_),
          cell_pattern_(other.cell_pattern_),
          evaluator_(other.evaluator_),
          cell_row_to_local_(
              std::move(
                  other.cell_row_to_local_)),
          state_sf_(other.state_sf_),
          local_state_(other.local_state_) {
        other.state_sf_ = nullptr;
        other.local_state_ = nullptr;
    }

    ~VariableCardinalityNaturalVariableSnesAssemblyContext3D() {
        if (local_state_ != nullptr) {
            (void)VecDestroy(
                &local_state_);
        }
        if (state_sf_ != nullptr) {
            (void)PetscSFDestroy(
                &state_sf_);
        }
    }

    [[nodiscard]] static PetscErrorCode
    create(
        MPI_Comm comm,
        const VariableCardinalityNaturalVariableNumbering3D&
            numbering,
        const mpmc::discretization_petsc::
            PetscMpiAijSymbolicPreallocation3D&
                cell_bridge,
        const mpmc::discretization_petsc::
            OwnedCellStructuralColumnPatternSnapshot3D&
                cell_pattern,
        VariableCardinalityNaturalVariableCellEvaluatorBinding3D
            evaluator,
        std::optional<
            VariableCardinalityNaturalVariableSnesAssemblyContext3D>*
            output) {
        using namespace
            variable_cardinality_snes_assembly_detail;

        if (output == nullptr) {
            return PETSC_ERR_ARG_NULL;
        }
        output->reset();

        int mpi_rank = -1;
        int mpi_size = -1;
        if (MPI_Comm_rank(
                comm,
                &mpi_rank) != MPI_SUCCESS ||
            MPI_Comm_size(
                comm,
                &mpi_size) != MPI_SUCCESS ||
            mpi_rank < 0 ||
            mpi_size <= 0) {
            return PETSC_ERR_MPI;
        }

        PetscErrorCode local_error =
            PETSC_SUCCESS;
        std::vector<
            std::pair<
                PetscInt,
                mpmc::mesh::LocalIndex>>
            row_to_local;

        try {
            if (evaluator.evaluator == nullptr ||
                numbering.local_rank().value() !=
                    static_cast<std::uint32_t>(
                        mpi_rank) ||
                numbering.rank_count() !=
                    static_cast<std::uint32_t>(
                        mpi_size) ||
                cell_bridge.local_rank() !=
                    numbering.local_rank() ||
                cell_bridge.rank_count() !=
                    numbering.rank_count() ||
                cell_pattern.local_rank() !=
                    numbering.local_rank() ||
                cell_pattern.rank_count() !=
                    numbering.rank_count() ||
                cell_pattern.local_cell_count() !=
                    numbering.local_cell_count() ||
                cell_pattern.global_row_start() !=
                    cell_bridge.global_row_start() ||
                cell_pattern.global_row_end() !=
                    cell_bridge.global_row_end() ||
                cell_pattern.global_row_count() !=
                    cell_bridge.global_row_count()) {
                throw std::invalid_argument(
                    "mixed-cardinality SNES metadata mismatch");
            }

            row_to_local.reserve(
                numbering.local_cell_count());
            for (std::size_t local = 0U;
                 local < numbering.local_cell_count();
                 ++local) {
                const auto cell =
                    mpmc::mesh::LocalIndex{
                        static_cast<
                            mpmc::mesh::LocalIndex::value_type>(
                                local)};
                const auto& record =
                    numbering.cell(cell);
                const PetscInt cell_row =
                    cell_bridge.global_row(
                        cell);
                row_to_local.emplace_back(
                    cell_row,
                    cell);

                if (record.owner_rank ==
                    numbering.local_rank()) {
                    if (!cell_pattern.contains_owned_cell(
                            cell)) {
                        throw std::invalid_argument(
                            "owned mixed-cardinality cell has no structural row");
                    }
                    const auto& pattern =
                        cell_pattern.row(
                            cell);
                    if (pattern.cell_global !=
                            record.cell_global ||
                        pattern.global_row !=
                            cell_row) {
                        throw std::invalid_argument(
                            "mixed-cardinality owned cell identity disagrees with structural row");
                    }
                }
            }

            std::sort(
                row_to_local.begin(),
                row_to_local.end(),
                [](const auto& left,
                   const auto& right) {
                    return left.first <
                        right.first;
                });
            if (std::adjacent_find(
                    row_to_local.begin(),
                    row_to_local.end(),
                    [](const auto& left,
                       const auto& right) {
                        return left.first ==
                            right.first;
                    }) !=
                row_to_local.end()) {
                throw std::invalid_argument(
                    "local overlap maps two cells to one structural global row");
            }
        } catch (...) {
            local_error =
                PETSC_ERR_ARG_INCOMP;
        }

        PetscErrorCode error =
            collective_error(
                comm,
                local_error);
        if (error != PETSC_SUCCESS) {
            return error;
        }

        VariableCardinalityNaturalVariableSnesAssemblyContext3D
            context{
                comm,
                numbering,
                cell_bridge,
                cell_pattern,
                evaluator,
                std::move(row_to_local)};

        error =
            context.initialize_state_exchange();
        error =
            collective_error(
                comm,
                error);
        if (error != PETSC_SUCCESS) {
            return error;
        }

        output->emplace(
            std::move(context));
        return PETSC_SUCCESS;
    }

    [[nodiscard]] NaturalVariableSnesEvaluator3D
    snes_evaluator() noexcept {
        return {
            &snes_function,
            &snes_jacobian,
            &snes_precheck,
            this};
    }

    [[nodiscard]] PetscErrorCode
    copy_local_packed_state(
        Vec global_state,
        std::vector<double>* output) {
        if (output == nullptr) {
            return PETSC_ERR_ARG_NULL;
        }
        output->clear();

        PetscErrorCode error =
            broadcast_state(
                global_state);
        error =
            variable_cardinality_snes_assembly_detail::
                collective_error(
                    comm_,
                    error);
        if (error != PETSC_SUCCESS) {
            return error;
        }

        const PetscScalar* values =
            nullptr;
        error =
            VecGetArrayRead(
                local_state_,
                &values);
        if (error != PETSC_SUCCESS) {
            return error;
        }

        output->resize(
            numbering_
                ->local_packed_scalar_count());
        for (std::size_t index = 0U;
             index < output->size();
             ++index) {
            (*output)[index] =
                static_cast<double>(
                    PetscRealPart(
                        values[index]));
        }

        return VecRestoreArrayRead(
            local_state_,
            &values);
    }

    [[nodiscard]] PetscErrorCode
    create_jacobian_structure(
        Mat* output) const {
        using namespace
            variable_cardinality_snes_assembly_detail;

        if (output == nullptr) {
            return PETSC_ERR_ARG_NULL;
        }
        if (*output != nullptr) {
            return PETSC_ERR_ARG_WRONGSTATE;
        }

        const PetscInt local_scalar_count =
            numbering_
                ->petsc_local_owned_scalar_count();
        const PetscInt global_scalar_count =
            numbering_
                ->petsc_global_scalar_count();
        if (local_scalar_count < 0 ||
            global_scalar_count <
                local_scalar_count) {
            return PETSC_ERR_ARG_INCOMP;
        }

        std::vector<PetscInt>
            diagonal_nnz;
        std::vector<PetscInt>
            off_diagonal_nnz;
        diagonal_nnz.reserve(
            static_cast<std::size_t>(
                local_scalar_count));
        off_diagonal_nnz.reserve(
            static_cast<std::size_t>(
                local_scalar_count));

        try {
            for (const auto& row :
                 cell_pattern_->rows()) {
                const auto& row_record =
                    numbering_->cell(
                        row.cell);
                std::size_t diagonal_width = 0U;
                std::size_t off_diagonal_width = 0U;

                for (const PetscInt column_row :
                     cell_pattern_
                         ->diagonal_global_columns(
                             row.cell)) {
                    diagonal_width +=
                        scalar_width_for_cell_row(
                            column_row);
                }
                for (const PetscInt column_row :
                     cell_pattern_
                         ->off_diagonal_global_columns(
                             row.cell)) {
                    off_diagonal_width +=
                        scalar_width_for_cell_row(
                            column_row);
                }

                if (diagonal_width == 0U ||
                    diagonal_width >
                        static_cast<std::size_t>(
                            std::numeric_limits<
                                PetscInt>::max()) ||
                    off_diagonal_width >
                        static_cast<std::size_t>(
                            std::numeric_limits<
                                PetscInt>::max())) {
                    throw std::length_error(
                        "mixed-cardinality scalar preallocation overflow");
                }

                const PetscInt diagonal =
                    static_cast<PetscInt>(
                        diagonal_width);
                const PetscInt off_diagonal =
                    static_cast<PetscInt>(
                        off_diagonal_width);
                for (std::size_t slot = 0U;
                     slot <
                        row_record.scalar_count;
                     ++slot) {
                    diagonal_nnz.push_back(
                        diagonal);
                    off_diagonal_nnz.push_back(
                        off_diagonal);
                }
            }

            if (diagonal_nnz.size() !=
                    static_cast<std::size_t>(
                        local_scalar_count) ||
                off_diagonal_nnz.size() !=
                    diagonal_nnz.size()) {
                throw std::runtime_error(
                    "mixed-cardinality scalar preallocation cardinality mismatch");
            }
        } catch (...) {
            return PETSC_ERR_ARG_INCOMP;
        }

        Mat matrix = nullptr;
        PetscErrorCode error =
            MatCreate(
                comm_,
                &matrix);
        if (error == PETSC_SUCCESS) {
            error =
                MatSetSizes(
                    matrix,
                    local_scalar_count,
                    local_scalar_count,
                    global_scalar_count,
                    global_scalar_count);
        }
        if (error == PETSC_SUCCESS) {
            error =
                MatSetType(
                    matrix,
                    MATMPIAIJ);
        }
        if (error == PETSC_SUCCESS) {
            error =
                MatMPIAIJSetPreallocation(
                    matrix,
                    0,
                    diagonal_nnz.empty()
                        ? nullptr
                        : diagonal_nnz.data(),
                    0,
                    off_diagonal_nnz.empty()
                        ? nullptr
                        : off_diagonal_nnz.data());
        }
        if (error == PETSC_SUCCESS) {
            error =
                MatSetOption(
                    matrix,
                    MAT_IGNORE_ZERO_ENTRIES,
                    PETSC_FALSE);
        }
        if (error != PETSC_SUCCESS) {
            (void)MatDestroy(
                &matrix);
            return error;
        }

        PetscErrorCode local_error =
            PETSC_SUCCESS;
        try {
            for (const auto& row :
                 cell_pattern_->rows()) {
                const auto& row_record =
                    numbering_->cell(
                        row.cell);
                std::vector<PetscInt>
                    scalar_rows(
                        row_record.scalar_count);
                for (std::size_t slot = 0U;
                     slot <
                        row_record.scalar_count;
                     ++slot) {
                    scalar_rows[slot] =
                        row_record
                            .petsc_global_scalar_start +
                        static_cast<PetscInt>(
                            slot);
                }

                const auto seed_columns =
                    [&](std::span<const PetscInt>
                            cell_columns) {
                        for (const PetscInt cell_row :
                             cell_columns) {
                            const auto column_cell =
                                local_cell_for_cell_row(
                                    cell_row);
                            const auto& column_record =
                                numbering_->cell(
                                    column_cell);
                            std::vector<PetscInt>
                                scalar_columns(
                                    column_record
                                        .scalar_count);
                            for (std::size_t slot = 0U;
                                 slot <
                                    column_record
                                        .scalar_count;
                                 ++slot) {
                                scalar_columns[slot] =
                                    column_record
                                        .petsc_global_scalar_start +
                                    static_cast<PetscInt>(
                                        slot);
                            }
                            std::vector<PetscScalar>
                                zeros(
                                    row_record.scalar_count *
                                        column_record
                                            .scalar_count,
                                    PetscScalar{0.0});
                            const PetscErrorCode seed_error =
                                MatSetValues(
                                    matrix,
                                    static_cast<PetscInt>(
                                        scalar_rows.size()),
                                    scalar_rows.data(),
                                    static_cast<PetscInt>(
                                        scalar_columns.size()),
                                    scalar_columns.data(),
                                    zeros.data(),
                                    INSERT_VALUES);
                            if (seed_error !=
                                PETSC_SUCCESS) {
                                return seed_error;
                            }
                        }
                        return PETSC_SUCCESS;
                    };

                local_error =
                    seed_columns(
                        cell_pattern_
                            ->diagonal_global_columns(
                                row.cell));
                if (local_error ==
                    PETSC_SUCCESS) {
                    local_error =
                        seed_columns(
                            cell_pattern_
                                ->off_diagonal_global_columns(
                                    row.cell));
                }
                if (local_error !=
                    PETSC_SUCCESS) {
                    break;
                }
            }
        } catch (...) {
            local_error =
                PETSC_ERR_ARG_INCOMP;
        }

        error =
            collective_error(
                comm_,
                local_error);
        if (error == PETSC_SUCCESS) {
            error =
                MatAssemblyBegin(
                    matrix,
                    MAT_FINAL_ASSEMBLY);
        }
        if (error == PETSC_SUCCESS) {
            error =
                MatAssemblyEnd(
                    matrix,
                    MAT_FINAL_ASSEMBLY);
        }
        if (error == PETSC_SUCCESS) {
            error =
                MatSetOption(
                    matrix,
                    MAT_NEW_NONZERO_LOCATION_ERR,
                    PETSC_TRUE);
        }
        if (error == PETSC_SUCCESS) {
            error =
                MatSetOption(
                    matrix,
                    MAT_NEW_NONZERO_ALLOCATION_ERR,
                    PETSC_TRUE);
        }
        if (error != PETSC_SUCCESS) {
            (void)MatDestroy(
                &matrix);
            return error;
        }

        *output = matrix;
        return PETSC_SUCCESS;
    }

private:
    VariableCardinalityNaturalVariableSnesAssemblyContext3D(
        MPI_Comm comm,
        const VariableCardinalityNaturalVariableNumbering3D&
            numbering,
        const mpmc::discretization_petsc::
            PetscMpiAijSymbolicPreallocation3D&
                cell_bridge,
        const mpmc::discretization_petsc::
            OwnedCellStructuralColumnPatternSnapshot3D&
                cell_pattern,
        VariableCardinalityNaturalVariableCellEvaluatorBinding3D
            evaluator,
        std::vector<
            std::pair<
                PetscInt,
                mpmc::mesh::LocalIndex>>
            cell_row_to_local)
        : comm_(comm),
          numbering_(&numbering),
          cell_bridge_(&cell_bridge),
          cell_pattern_(&cell_pattern),
          evaluator_(evaluator),
          cell_row_to_local_(
              std::move(
                  cell_row_to_local)) {}

    [[nodiscard]] mpmc::mesh::LocalIndex
    local_cell_for_cell_row(
        PetscInt cell_row) const {
        const auto found =
            std::lower_bound(
                cell_row_to_local_.begin(),
                cell_row_to_local_.end(),
                cell_row,
                [](const auto& entry,
                   PetscInt value) {
                    return entry.first <
                        value;
                });
        if (found ==
                cell_row_to_local_.end() ||
            found->first !=
                cell_row) {
            throw std::invalid_argument(
                "structural cell row is not present in local overlap");
        }
        return found->second;
    }

    [[nodiscard]] std::size_t
    scalar_width_for_cell_row(
        PetscInt cell_row) const {
        return numbering_
            ->cell(
                local_cell_for_cell_row(
                    cell_row))
            .scalar_count;
    }

    [[nodiscard]] PetscErrorCode
    initialize_state_exchange() {
        const PetscInt nroots =
            numbering_
                ->petsc_local_owned_scalar_count();
        if (nroots < 0 ||
            numbering_
                    ->local_packed_scalar_count() >
                static_cast<std::size_t>(
                    std::numeric_limits<
                        PetscInt>::max())) {
            return PETSC_ERR_ARG_OUTOFRANGE;
        }
        const PetscInt nleaves =
            static_cast<PetscInt>(
                numbering_
                    ->local_packed_scalar_count());

        std::array<PetscInt, 2>
            local_range{
                numbering_
                    ->petsc_owned_scalar_start(),
                numbering_
                    ->petsc_owned_scalar_end()};
        std::vector<PetscInt>
            all_ranges(
                static_cast<std::size_t>(
                    numbering_->rank_count()) *
                2U);
        if (MPI_Allgather(
                local_range.data(),
                2,
                MPIU_INT,
                all_ranges.data(),
                2,
                MPIU_INT,
                comm_) != MPI_SUCCESS) {
            return PETSC_ERR_MPI;
        }

        std::vector<PetscSFNode>
            remote(
                static_cast<std::size_t>(
                    nleaves));
        for (const auto& record :
             numbering_->cells()) {
            const std::size_t owner =
                static_cast<std::size_t>(
                    record.owner_rank.value());
            if (owner >=
                numbering_->rank_count()) {
                return PETSC_ERR_ARG_INCOMP;
            }
            const PetscInt owner_start =
                all_ranges[
                    owner * 2U];
            const PetscInt owner_end =
                all_ranges[
                    owner * 2U + 1U];
            if (record.petsc_global_scalar_start <
                    owner_start ||
                record.petsc_global_scalar_start +
                        static_cast<PetscInt>(
                            record.scalar_count) >
                    owner_end) {
                return PETSC_ERR_ARG_INCOMP;
            }

            for (std::size_t slot = 0U;
                 slot <
                    record.scalar_count;
                 ++slot) {
                const std::size_t leaf =
                    record.local_scalar_offset +
                    slot;
                remote[leaf] =
                    PetscSFNode{
                        static_cast<PetscInt>(
                            record.owner_rank
                                .value()),
                        record
                                .petsc_global_scalar_start -
                            owner_start +
                            static_cast<PetscInt>(
                                slot)};
            }
        }

        PetscErrorCode error =
            PetscSFCreate(
                comm_,
                &state_sf_);
        if (error == PETSC_SUCCESS) {
            error =
                PetscSFSetGraph(
                    state_sf_,
                    nroots,
                    nleaves,
                    nullptr,
                    PETSC_COPY_VALUES,
                    remote.data(),
                    PETSC_COPY_VALUES);
        }
        if (error == PETSC_SUCCESS) {
            error =
                PetscSFSetUp(
                    state_sf_);
        }
        if (error == PETSC_SUCCESS) {
            error =
                VecCreateSeq(
                    PETSC_COMM_SELF,
                    nleaves,
                    &local_state_);
        }
        return error;
    }

    [[nodiscard]] PetscErrorCode
    broadcast_state(
        Vec global_state) {
        if (global_state == nullptr ||
            state_sf_ == nullptr ||
            local_state_ == nullptr) {
            return PETSC_ERR_ARG_NULL;
        }

        PetscInt local_size = -1;
        PetscErrorCode error =
            VecGetLocalSize(
                global_state,
                &local_size);
        if (error != PETSC_SUCCESS) {
            return error;
        }
        if (local_size !=
            numbering_
                ->petsc_local_owned_scalar_count()) {
            return PETSC_ERR_ARG_SIZ;
        }

        const PetscScalar* roots =
            nullptr;
        PetscScalar* leaves =
            nullptr;
        error =
            VecGetArrayRead(
                global_state,
                &roots);
        if (error == PETSC_SUCCESS) {
            error =
                VecGetArray(
                    local_state_,
                    &leaves);
        }
        if (error != PETSC_SUCCESS) {
            if (roots != nullptr) {
                (void)VecRestoreArrayRead(
                    global_state,
                    &roots);
            }
            return error;
        }

        PetscErrorCode communication =
            PetscSFBcastBegin(
                state_sf_,
                MPIU_SCALAR,
                roots,
                leaves,
                MPI_REPLACE);
        if (communication ==
            PETSC_SUCCESS) {
            communication =
                PetscSFBcastEnd(
                    state_sf_,
                    MPIU_SCALAR,
                    roots,
                    leaves,
                    MPI_REPLACE);
        }

        const PetscErrorCode leaf_restore =
            VecRestoreArray(
                local_state_,
                &leaves);
        const PetscErrorCode root_restore =
            VecRestoreArrayRead(
                global_state,
                &roots);
        if (communication !=
            PETSC_SUCCESS) {
            return communication;
        }
        if (leaf_restore !=
            PETSC_SUCCESS) {
            return leaf_restore;
        }
        return root_restore;
    }

    [[nodiscard]] PetscErrorCode
    evaluate_and_insert(
        Vec global_state,
        Vec residual,
        Mat jacobian,
        bool insert_residual,
        bool insert_jacobian,
        NaturalVariableSnesEvaluationStatus3D*
            status) {
        using namespace
            variable_cardinality_snes_assembly_detail;

        if (global_state == nullptr ||
            status == nullptr ||
            (insert_residual &&
             residual == nullptr) ||
            (insert_jacobian &&
             jacobian == nullptr)) {
            return PETSC_ERR_ARG_NULL;
        }
        *status =
            NaturalVariableSnesEvaluationStatus3D::
                success;

        PetscErrorCode error =
            broadcast_state(
                global_state);
        error =
            collective_error(
                comm_,
                error);
        if (error != PETSC_SUCCESS) {
            return error;
        }

        const PetscScalar* local_values =
            nullptr;
        error =
            VecGetArrayRead(
                local_state_,
                &local_values);
        if (error != PETSC_SUCCESS) {
            return error;
        }

        std::vector<double>
            local_state(
                numbering_
                    ->local_packed_scalar_count());
        for (std::size_t index = 0U;
             index < local_state.size();
             ++index) {
            local_state[index] =
                static_cast<double>(
                    PetscRealPart(
                        local_values[index]));
        }
        const PetscErrorCode restore_error =
            VecRestoreArrayRead(
                local_state_,
                &local_values);
        if (restore_error !=
            PETSC_SUCCESS) {
            return restore_error;
        }

        PetscErrorCode local_error =
            PETSC_SUCCESS;
        bool local_domain_error = false;

        try {
            for (const auto& row_record :
                 numbering_->cells()) {
                if (row_record.owner_rank !=
                    numbering_->local_rank()) {
                    continue;
                }

                VariableCardinalityNaturalVariableCellAssembly3D
                    assembly;
                NaturalVariableSnesEvaluationStatus3D
                    cell_status =
                        NaturalVariableSnesEvaluationStatus3D::
                            success;
                const PetscErrorCode cell_error =
                    evaluator_.evaluator(
                        row_record.cell,
                        row_record.phase_count,
                        local_state,
                        *numbering_,
                        evaluator_.user_context,
                        &assembly,
                        &cell_status);
                if (cell_error !=
                    PETSC_SUCCESS) {
                    local_error =
                        cell_error;
                    break;
                }
                if (cell_status ==
                    NaturalVariableSnesEvaluationStatus3D::
                        domain_error) {
                    local_domain_error = true;
                    continue;
                }
                if (cell_status !=
                    NaturalVariableSnesEvaluationStatus3D::
                        success) {
                    local_error =
                        PETSC_ERR_ARG_INCOMP;
                    break;
                }

                if (assembly.residual.size() !=
                    row_record.scalar_count) {
                    local_error =
                        PETSC_ERR_ARG_SIZ;
                    break;
                }
                for (const double value :
                     assembly.residual) {
                    if (!std::isfinite(value)) {
                        local_error =
                            PETSC_ERR_FP;
                        break;
                    }
                }
                if (local_error !=
                    PETSC_SUCCESS) {
                    break;
                }

                const auto diagonal_columns =
                    cell_pattern_
                        ->diagonal_global_columns(
                            row_record.cell);
                const auto off_diagonal_columns =
                    cell_pattern_
                        ->off_diagonal_global_columns(
                            row_record.cell);
                std::vector<PetscInt>
                    seen_columns;
                seen_columns.reserve(
                    assembly
                        .jacobian_blocks
                        .size());
                bool has_self = false;

                for (const auto& block :
                     assembly.jacobian_blocks) {
                    const auto& column_record =
                        numbering_->cell(
                            block.column_cell);
                    const PetscInt column_cell_row =
                        cell_bridge_->global_row(
                            block.column_cell);
                    if (!contains_cell_row(
                            diagonal_columns,
                            column_cell_row) &&
                        !contains_cell_row(
                            off_diagonal_columns,
                            column_cell_row)) {
                        local_error =
                            PETSC_ERR_ARG_INCOMP;
                        break;
                    }

                    const std::size_t expected =
                        row_record.scalar_count *
                        column_record.scalar_count;
                    if (block.values_row_major
                            .size() !=
                        expected) {
                        local_error =
                            PETSC_ERR_ARG_SIZ;
                        break;
                    }
                    for (const double value :
                         block.values_row_major) {
                        if (!std::isfinite(value)) {
                            local_error =
                                PETSC_ERR_FP;
                            break;
                        }
                    }
                    if (local_error !=
                        PETSC_SUCCESS) {
                        break;
                    }

                    if (std::find(
                            seen_columns.begin(),
                            seen_columns.end(),
                            column_cell_row) !=
                        seen_columns.end()) {
                        local_error =
                            PETSC_ERR_ARG_INCOMP;
                        break;
                    }
                    seen_columns.push_back(
                        column_cell_row);
                    has_self =
                        has_self ||
                        block.column_cell ==
                            row_record.cell;
                }
                if (local_error !=
                    PETSC_SUCCESS) {
                    break;
                }
                if (!has_self) {
                    local_error =
                        PETSC_ERR_ARG_INCOMP;
                    break;
                }

                if (insert_residual) {
                    std::vector<PetscInt>
                        rows(
                            row_record.scalar_count);
                    std::vector<PetscScalar>
                        values(
                            row_record.scalar_count);
                    for (std::size_t slot = 0U;
                         slot <
                            row_record.scalar_count;
                         ++slot) {
                        rows[slot] =
                            row_record
                                .petsc_global_scalar_start +
                            static_cast<PetscInt>(
                                slot);
                        values[slot] =
                            static_cast<PetscScalar>(
                                assembly
                                    .residual[slot]);
                    }
                    local_error =
                        VecSetValues(
                            residual,
                            static_cast<PetscInt>(
                                rows.size()),
                            rows.data(),
                            values.data(),
                            INSERT_VALUES);
                    if (local_error !=
                        PETSC_SUCCESS) {
                        break;
                    }
                }

                if (insert_jacobian) {
                    std::vector<PetscInt>
                        rows(
                            row_record.scalar_count);
                    for (std::size_t slot = 0U;
                         slot <
                            row_record.scalar_count;
                         ++slot) {
                        rows[slot] =
                            row_record
                                .petsc_global_scalar_start +
                            static_cast<PetscInt>(
                                slot);
                    }

                    for (const auto& block :
                         assembly
                             .jacobian_blocks) {
                        const auto& column_record =
                            numbering_->cell(
                                block.column_cell);
                        std::vector<PetscInt>
                            columns(
                                column_record
                                    .scalar_count);
                        for (std::size_t slot = 0U;
                             slot <
                                column_record
                                    .scalar_count;
                             ++slot) {
                            columns[slot] =
                                column_record
                                    .petsc_global_scalar_start +
                                static_cast<PetscInt>(
                                    slot);
                        }
                        std::vector<PetscScalar>
                            values;
                        values.reserve(
                            block
                                .values_row_major
                                .size());
                        for (const double value :
                             block
                                 .values_row_major) {
                            values.push_back(
                                static_cast<
                                    PetscScalar>(
                                        value));
                        }
                        local_error =
                            MatSetValues(
                                jacobian,
                                static_cast<PetscInt>(
                                    rows.size()),
                                rows.data(),
                                static_cast<PetscInt>(
                                    columns.size()),
                                columns.data(),
                                values.data(),
                                INSERT_VALUES);
                        if (local_error !=
                            PETSC_SUCCESS) {
                            break;
                        }
                    }
                    if (local_error !=
                        PETSC_SUCCESS) {
                        break;
                    }
                }
            }
        } catch (...) {
            local_error =
                PETSC_ERR_ARG_INCOMP;
        }

        error =
            collective_error(
                comm_,
                local_error);
        if (error != PETSC_SUCCESS) {
            return error;
        }

        int local_domain =
            local_domain_error ? 1 : 0;
        int global_domain = 0;
        if (MPI_Allreduce(
                &local_domain,
                &global_domain,
                1,
                MPI_INT,
                MPI_MAX,
                comm_) != MPI_SUCCESS) {
            return PETSC_ERR_MPI;
        }
        if (global_domain != 0) {
            *status =
                NaturalVariableSnesEvaluationStatus3D::
                    domain_error;
        }
        return PETSC_SUCCESS;
    }

    [[nodiscard]] PetscErrorCode
    raw_positive_support_precheck(
        Vec state,
        Vec search_direction,
        PetscBool* changed_direction) {
        if (state == nullptr ||
            search_direction == nullptr ||
            changed_direction == nullptr) {
            return PETSC_ERR_ARG_NULL;
        }
        *changed_direction =
            PETSC_FALSE;

        PetscInt state_local = -1;
        PetscInt direction_local = -1;
        PetscErrorCode error =
            VecGetLocalSize(
                state,
                &state_local);
        if (error == PETSC_SUCCESS) {
            error =
                VecGetLocalSize(
                    search_direction,
                    &direction_local);
        }
        if (error != PETSC_SUCCESS) {
            return error;
        }
        if (state_local !=
                direction_local ||
            state_local !=
                numbering_
                    ->petsc_local_owned_scalar_count()) {
            return PETSC_ERR_ARG_SIZ;
        }

        const PetscScalar* state_values =
            nullptr;
        const PetscScalar* direction_values =
            nullptr;
        error =
            VecGetArrayRead(
                state,
                &state_values);
        if (error == PETSC_SUCCESS) {
            error =
                VecGetArrayRead(
                    search_direction,
                    &direction_values);
        }
        if (error != PETSC_SUCCESS) {
            if (state_values != nullptr) {
                (void)VecRestoreArrayRead(
                    state,
                    &state_values);
            }
            return error;
        }

        double local_scale = 1.0;
        bool invalid = false;
        for (PetscInt index = 0;
             index < state_local;
             ++index) {
            const double x =
                static_cast<double>(
                    PetscRealPart(
                        state_values[index]));
            const double y =
                static_cast<double>(
                    PetscRealPart(
                        direction_values[index]));
            if (!std::isfinite(x) ||
                !std::isfinite(y) ||
                !(x > 0.0)) {
                invalid = true;
                break;
            }
            if (y > 0.0 &&
                x - y <= 0.0) {
                local_scale =
                    std::min(
                        local_scale,
                        0.8 * x / y);
            }
        }

        const PetscErrorCode direction_restore =
            VecRestoreArrayRead(
                search_direction,
                &direction_values);
        const PetscErrorCode state_restore =
            VecRestoreArrayRead(
                state,
                &state_values);
        if (direction_restore !=
            PETSC_SUCCESS) {
            return direction_restore;
        }
        if (state_restore !=
            PETSC_SUCCESS) {
            return state_restore;
        }
        if (invalid) {
            return PETSC_ERR_ARG_OUTOFRANGE;
        }

        double global_scale = 1.0;
        if (MPI_Allreduce(
                &local_scale,
                &global_scale,
                1,
                MPI_DOUBLE,
                MPI_MIN,
                comm_) != MPI_SUCCESS) {
            return PETSC_ERR_MPI;
        }
        if (!std::isfinite(
                global_scale) ||
            !(global_scale > 0.0)) {
            return PETSC_ERR_FP;
        }
        if (global_scale < 1.0) {
            error =
                VecScale(
                    search_direction,
                    static_cast<PetscScalar>(
                        global_scale));
            if (error != PETSC_SUCCESS) {
                return error;
            }
            *changed_direction =
                PETSC_TRUE;
        }
        return PETSC_SUCCESS;
    }

    static PetscErrorCode
    snes_function(
        Vec state,
        Vec residual,
        void* raw_context,
        NaturalVariableSnesEvaluationStatus3D*
            status) {
        if (raw_context == nullptr) {
            return PETSC_ERR_ARG_NULL;
        }
        return static_cast<
            VariableCardinalityNaturalVariableSnesAssemblyContext3D*>(
                raw_context)
            ->evaluate_and_insert(
                state,
                residual,
                nullptr,
                true,
                false,
                status);
    }

    static PetscErrorCode
    snes_jacobian(
        Vec state,
        Mat jacobian,
        void* raw_context,
        NaturalVariableSnesEvaluationStatus3D*
            status) {
        if (raw_context == nullptr) {
            return PETSC_ERR_ARG_NULL;
        }
        return static_cast<
            VariableCardinalityNaturalVariableSnesAssemblyContext3D*>(
                raw_context)
            ->evaluate_and_insert(
                state,
                nullptr,
                jacobian,
                false,
                true,
                status);
    }

    static PetscErrorCode
    snes_precheck(
        Vec state,
        Vec search_direction,
        void* raw_context,
        PetscBool* changed_direction) {
        if (raw_context == nullptr) {
            return PETSC_ERR_ARG_NULL;
        }
        return static_cast<
            VariableCardinalityNaturalVariableSnesAssemblyContext3D*>(
                raw_context)
            ->raw_positive_support_precheck(
                state,
                search_direction,
                changed_direction);
    }

    MPI_Comm comm_;
    const VariableCardinalityNaturalVariableNumbering3D*
        numbering_;
    const mpmc::discretization_petsc::
        PetscMpiAijSymbolicPreallocation3D*
            cell_bridge_;
    const mpmc::discretization_petsc::
        OwnedCellStructuralColumnPatternSnapshot3D*
            cell_pattern_;
    VariableCardinalityNaturalVariableCellEvaluatorBinding3D
        evaluator_;
    std::vector<
        std::pair<
            PetscInt,
            mpmc::mesh::LocalIndex>>
        cell_row_to_local_;
    PetscSF state_sf_{};
    Vec local_state_{};
};

namespace variable_cardinality_snes_solver_detail {

[[nodiscard]] inline PetscErrorCode
validate_layout(
    MPI_Comm comm,
    const VariableCardinalityNaturalVariableNumbering3D&
        numbering,
    Vec state,
    Mat jacobian) {
    using namespace
        variable_cardinality_snes_assembly_detail;

    if (state == nullptr ||
        jacobian == nullptr) {
        return PETSC_ERR_ARG_NULL;
    }
    if (!communicators_are_compatible(
            comm,
            PetscObjectComm(
                reinterpret_cast<PetscObject>(
                    state))) ||
        !communicators_are_compatible(
            comm,
            PetscObjectComm(
                reinterpret_cast<PetscObject>(
                    jacobian)))) {
        return PETSC_ERR_ARG_INCOMP;
    }

    PetscBool assembled =
        PETSC_FALSE;
    PetscErrorCode error =
        MatAssembled(
            jacobian,
            &assembled);
    if (error != PETSC_SUCCESS ||
        assembled != PETSC_TRUE) {
        return error != PETSC_SUCCESS
            ? error
            : PETSC_ERR_ARG_WRONGSTATE;
    }

    PetscInt vec_local = -1;
    PetscInt vec_global = -1;
    PetscInt vec_start = -1;
    PetscInt vec_end = -1;
    PetscInt mat_local_rows = -1;
    PetscInt mat_local_columns = -1;
    PetscInt mat_global_rows = -1;
    PetscInt mat_global_columns = -1;
    PetscInt mat_start = -1;
    PetscInt mat_end = -1;
    PetscInt mat_column_start = -1;
    PetscInt mat_column_end = -1;

    error =
        VecGetLocalSize(
            state,
            &vec_local);
    if (error == PETSC_SUCCESS) {
        error =
            VecGetSize(
                state,
                &vec_global);
    }
    if (error == PETSC_SUCCESS) {
        error =
            VecGetOwnershipRange(
                state,
                &vec_start,
                &vec_end);
    }
    if (error == PETSC_SUCCESS) {
        error =
            MatGetLocalSize(
                jacobian,
                &mat_local_rows,
                &mat_local_columns);
    }
    if (error == PETSC_SUCCESS) {
        error =
            MatGetSize(
                jacobian,
                &mat_global_rows,
                &mat_global_columns);
    }
    if (error == PETSC_SUCCESS) {
        error =
            MatGetOwnershipRange(
                jacobian,
                &mat_start,
                &mat_end);
    }
    if (error == PETSC_SUCCESS) {
        error =
            MatGetOwnershipRangeColumn(
                jacobian,
                &mat_column_start,
                &mat_column_end);
    }
    if (error != PETSC_SUCCESS) {
        return error;
    }

    const PetscInt expected_local =
        numbering
            .petsc_local_owned_scalar_count();
    const PetscInt expected_global =
        numbering
            .petsc_global_scalar_count();
    const PetscInt expected_start =
        numbering
            .petsc_owned_scalar_start();
    const PetscInt expected_end =
        numbering
            .petsc_owned_scalar_end();

    if (vec_local !=
            expected_local ||
        vec_global !=
            expected_global ||
        vec_start !=
            expected_start ||
        vec_end !=
            expected_end ||
        mat_local_rows !=
            expected_local ||
        mat_local_columns !=
            expected_local ||
        mat_global_rows !=
            expected_global ||
        mat_global_columns !=
            expected_global ||
        mat_start !=
            expected_start ||
        mat_end !=
            expected_end ||
        mat_column_start !=
            expected_start ||
        mat_column_end !=
            expected_end) {
        return PETSC_ERR_ARG_INCOMP;
    }
    return PETSC_SUCCESS;
}

[[nodiscard]] inline PetscErrorCode
validate_row_scaling(
    MPI_Comm comm,
    const VariableCardinalityNaturalVariableNumbering3D&
        numbering,
    Vec row_scaling) {
    using namespace
        natural_variable_snes_detail;

    if (row_scaling == nullptr) {
        return PETSC_SUCCESS;
    }
    if (!communicators_are_compatible(
            comm,
            PetscObjectComm(
                reinterpret_cast<PetscObject>(
                    row_scaling)))) {
        return PETSC_ERR_ARG_INCOMP;
    }

    PetscInt local = -1;
    PetscInt global = -1;
    PetscInt start = -1;
    PetscInt end = -1;
    PetscErrorCode error =
        VecGetLocalSize(
            row_scaling,
            &local);
    if (error == PETSC_SUCCESS) {
        error =
            VecGetSize(
                row_scaling,
                &global);
    }
    if (error == PETSC_SUCCESS) {
        error =
            VecGetOwnershipRange(
                row_scaling,
                &start,
                &end);
    }
    if (error != PETSC_SUCCESS) {
        return error;
    }

    if (local !=
            numbering
                .petsc_local_owned_scalar_count() ||
        global !=
            numbering
                .petsc_global_scalar_count() ||
        start !=
            numbering
                .petsc_owned_scalar_start() ||
        end !=
            numbering
                .petsc_owned_scalar_end()) {
        return PETSC_ERR_ARG_INCOMP;
    }

    const PetscScalar* values = nullptr;
    error =
        VecGetArrayRead(
            row_scaling,
            &values);
    if (error != PETSC_SUCCESS) {
        return error;
    }
    bool valid = true;
    for (PetscInt index = 0;
         index < local;
         ++index) {
        const double value =
            static_cast<double>(
                PetscRealPart(
                    values[index]));
        if (!std::isfinite(value) ||
            !(value > 0.0)) {
            valid = false;
            break;
        }
    }
    const PetscErrorCode restore =
        VecRestoreArrayRead(
            row_scaling,
            &values);
    if (restore != PETSC_SUCCESS) {
        return restore;
    }
    return valid
        ? PETSC_SUCCESS
        : PETSC_ERR_ARG_OUTOFRANGE;
}

} // namespace variable_cardinality_snes_solver_detail

/// Construct a frozen positive left row-equilibration vector from the analytic
/// q-ragged Jacobian at one initial state.
///
/// For each owned scalar equation row i:
///   D_i = 1 / max_j |J_ij(q0)|.
///
/// The evaluator's physical residual/Jacobian is not mutated. D is caller-owned
/// and is intended to remain frozen for the subsequent nonlinear solve.
[[nodiscard]] inline PetscErrorCode
make_variable_cardinality_initial_row_equilibration_3d(
    MPI_Comm comm,
    const VariableCardinalityNaturalVariableNumbering3D&
        numbering,
    Vec initial_state,
    Mat jacobian_structure_template,
    NaturalVariableSnesEvaluator3D evaluator,
    Vec* row_scaling) {
    using namespace
        variable_cardinality_snes_solver_detail;
    using variable_cardinality_snes_assembly_detail::
        collective_error;

    if (row_scaling == nullptr) {
        return PETSC_ERR_ARG_NULL;
    }
    if (*row_scaling != nullptr) {
        return PETSC_ERR_ARG_WRONGSTATE;
    }
    if (evaluator.jacobian == nullptr) {
        return PETSC_ERR_ARG_NULL;
    }

    PetscErrorCode error =
        validate_layout(
            comm,
            numbering,
            initial_state,
            jacobian_structure_template);
    error =
        collective_error(
            comm,
            error);
    if (error != PETSC_SUCCESS) {
        return error;
    }

    Mat jacobian = nullptr;
    Vec scaling = nullptr;
    error =
        MatDuplicate(
            jacobian_structure_template,
            MAT_DO_NOT_COPY_VALUES,
            &jacobian);
    if (error == PETSC_SUCCESS) {
        error =
            MatZeroEntries(
                jacobian);
    }

    NaturalVariableSnesEvaluationStatus3D
        status =
            NaturalVariableSnesEvaluationStatus3D::
                success;
    if (error == PETSC_SUCCESS) {
        error =
            evaluator.jacobian(
                initial_state,
                jacobian,
                evaluator.user_context,
                &status);
    }
    if (error == PETSC_SUCCESS) {
        int local_status =
            status ==
                    NaturalVariableSnesEvaluationStatus3D::
                        success
                ? 0
                : (status ==
                           NaturalVariableSnesEvaluationStatus3D::
                               domain_error
                       ? 1
                       : 2);
        int global_status = 0;
        if (MPI_Allreduce(
                &local_status,
                &global_status,
                1,
                MPI_INT,
                MPI_MAX,
                comm) != MPI_SUCCESS) {
            error =
                PETSC_ERR_MPI;
        } else if (global_status != 0) {
            error =
                global_status == 1
                ? PETSC_ERR_ARG_OUTOFRANGE
                : PETSC_ERR_ARG_INCOMP;
        }
    }
    if (error == PETSC_SUCCESS) {
        error =
            MatAssemblyBegin(
                jacobian,
                MAT_FINAL_ASSEMBLY);
    }
    if (error == PETSC_SUCCESS) {
        error =
            MatAssemblyEnd(
                jacobian,
                MAT_FINAL_ASSEMBLY);
    }
    if (error == PETSC_SUCCESS) {
        error =
            VecDuplicate(
                initial_state,
                &scaling);
    }

    PetscInt start = -1;
    PetscInt end = -1;
    if (error == PETSC_SUCCESS) {
        error =
            MatGetOwnershipRange(
                jacobian,
                &start,
                &end);
    }
    if (error == PETSC_SUCCESS &&
        (start !=
             numbering
                 .petsc_owned_scalar_start() ||
         end !=
             numbering
                 .petsc_owned_scalar_end())) {
        error =
            PETSC_ERR_ARG_INCOMP;
    }

    PetscScalar* values = nullptr;
    if (error == PETSC_SUCCESS) {
        error =
            VecGetArray(
                scaling,
                &values);
    }
    if (error == PETSC_SUCCESS) {
        for (PetscInt row = start;
             row < end;
             ++row) {
            PetscInt count = 0;
            const PetscInt* columns =
                nullptr;
            const PetscScalar* row_values =
                nullptr;
            error =
                MatGetRow(
                    jacobian,
                    row,
                    &count,
                    &columns,
                    &row_values);
            if (error != PETSC_SUCCESS) {
                break;
            }
            double maximum = 0.0;
            for (PetscInt entry = 0;
                 entry < count;
                 ++entry) {
                const double magnitude =
                    std::abs(
                        static_cast<double>(
                            PetscRealPart(
                                row_values[
                                    entry])));
                if (!std::isfinite(
                        magnitude)) {
                    error =
                        PETSC_ERR_FP;
                    break;
                }
                maximum =
                    std::max(
                        maximum,
                        magnitude);
            }
            const PetscErrorCode restore =
                MatRestoreRow(
                    jacobian,
                    row,
                    &count,
                    &columns,
                    &row_values);
            if (error == PETSC_SUCCESS &&
                restore != PETSC_SUCCESS) {
                error =
                    restore;
            }
            if (error != PETSC_SUCCESS) {
                break;
            }
            if (!(maximum > 0.0) ||
                !std::isfinite(maximum)) {
                error =
                    PETSC_ERR_ARG_WRONGSTATE;
                break;
            }
            values[
                row - start] =
                static_cast<PetscScalar>(
                    1.0 / maximum);
        }
    }

    if (values != nullptr) {
        const PetscErrorCode restore =
            VecRestoreArray(
                scaling,
                &values);
        if (error == PETSC_SUCCESS &&
            restore != PETSC_SUCCESS) {
            error =
                restore;
        }
    }
    const PetscErrorCode destroy =
        MatDestroy(
            &jacobian);
    if (error == PETSC_SUCCESS &&
        destroy != PETSC_SUCCESS) {
        error =
            destroy;
    }
    if (error != PETSC_SUCCESS) {
        if (scaling != nullptr) {
            (void)VecDestroy(
                &scaling);
        }
        return error;
    }

    error =
        validate_row_scaling(
            comm,
            numbering,
            scaling);
    if (error != PETSC_SUCCESS) {
        (void)VecDestroy(
            &scaling);
        return error;
    }

    *row_scaling =
        scaling;
    return PETSC_SUCCESS;
}

/// Solve one distributed q-ragged natural-variable system with PETSc-owned
/// nonlinear orchestration.
///
/// Solver policy intentionally matches the existing fixed-cardinality path:
/// SNESNEWTONLS -> SNESLINESEARCHBT -> KSPGMRES -> PCASM(restrict, overlap=1).
/// Runtime SNESSetFromOptions() remains disabled so an unvalidated option
/// cannot replace the analytic/AD Jacobian path with finite differences.
///
/// An optional frozen positive left row-scaling vector D applies the same
/// audited contract as the fixed-cardinality solver: PETSc solves D*R=0 with
/// D*J while the physical evaluator continues to publish native R and J.
/// With D enabled, ASM sub-block LU uses the same weak-diagonal reordering and
/// MAT_SHIFT_NONZERO stabilization as the fixed-cardinality correctness path.
/// Optional failure diagnostics preserve SNES/KSP/PC reasons without changing
/// non-convergence semantics.
inline PetscErrorCode
solve_variable_cardinality_natural_variable_snes_3d(
    MPI_Comm comm,
    const VariableCardinalityNaturalVariableNumbering3D&
        numbering,
    Vec initial_state,
    Mat jacobian_structure_template,
    NaturalVariableSnesEvaluator3D evaluator,
    Vec* solution,
    std::optional<
        VariableCardinalityNaturalVariableSnesSolveReport3D>*
        report,
    Vec row_scaling = nullptr,
    std::optional<
        NaturalVariableSnesFailureDiagnostics3D>*
            failure_diagnostics = nullptr) {
    using namespace
        variable_cardinality_snes_solver_detail;
    using namespace natural_variable_snes_detail;
    using variable_cardinality_snes_assembly_detail::
        collective_error;

    int mpi_rank = -1;
    int mpi_size = -1;
    if (MPI_Comm_rank(
            comm,
            &mpi_rank) != MPI_SUCCESS ||
        MPI_Comm_size(
            comm,
            &mpi_size) != MPI_SUCCESS ||
        mpi_rank < 0 ||
        mpi_size <= 0) {
        return PETSC_ERR_MPI;
    }

    PetscErrorCode local_error =
        solution == nullptr ||
                report == nullptr
            ? PETSC_ERR_ARG_NULL
            : PETSC_SUCCESS;
    if (local_error == PETSC_SUCCESS &&
        *solution != nullptr) {
        local_error =
            PETSC_ERR_ARG_WRONGSTATE;
    }
    if (local_error == PETSC_SUCCESS &&
        (evaluator.function == nullptr ||
         evaluator.jacobian == nullptr ||
         numbering.local_rank().value() !=
             static_cast<std::uint32_t>(
                 mpi_rank) ||
         numbering.rank_count() !=
             static_cast<std::uint32_t>(
                 mpi_size))) {
        local_error =
            PETSC_ERR_ARG_INCOMP;
    }

    PetscErrorCode error =
        collective_error(
            comm,
            local_error);
    if (error != PETSC_SUCCESS) {
        return error;
    }

    *solution = nullptr;
    report->reset();
    if (failure_diagnostics != nullptr) {
        failure_diagnostics->reset();
    }

    local_error =
        validate_layout(
            comm,
            numbering,
            initial_state,
            jacobian_structure_template);
    error =
        collective_error(
            comm,
            local_error);
    if (error == PETSC_SUCCESS) {
        local_error =
            validate_row_scaling(
                comm,
                numbering,
                row_scaling);
        error =
            collective_error(
                comm,
                local_error);
    }
    if (error != PETSC_SUCCESS) {
        return error;
    }

    Vec solved_state = nullptr;
    Vec residual = nullptr;
    Mat jacobian = nullptr;
    SNES snes = nullptr;

    const auto cleanup =
        [&]() {
            if (snes != nullptr) {
                (void)SNESDestroy(
                    &snes);
            }
            if (jacobian != nullptr) {
                (void)MatDestroy(
                    &jacobian);
            }
            if (residual != nullptr) {
                (void)VecDestroy(
                    &residual);
            }
            if (solved_state != nullptr) {
                (void)VecDestroy(
                    &solved_state);
            }
        };

    error =
        VecDuplicate(
            initial_state,
            &solved_state);
    if (error == PETSC_SUCCESS) {
        error =
            VecCopy(
                initial_state,
                solved_state);
    }
    if (error == PETSC_SUCCESS) {
        error =
            VecDuplicate(
                initial_state,
                &residual);
    }
    if (error == PETSC_SUCCESS) {
        error =
            MatDuplicate(
                jacobian_structure_template,
                MAT_COPY_VALUES,
                &jacobian);
    }
    if (error == PETSC_SUCCESS) {
        error =
            SNESCreate(
                comm,
                &snes);
    }
    if (error != PETSC_SUCCESS) {
        cleanup();
        return error;
    }

    CallbackContext callback_context{
        evaluator,
        row_scaling};

    error =
        SNESSetFunction(
            snes,
            residual,
            form_function,
            &callback_context);
    if (error == PETSC_SUCCESS) {
        error =
            SNESSetJacobian(
                snes,
                jacobian,
                jacobian,
                form_jacobian,
                &callback_context);
    }
    if (error == PETSC_SUCCESS) {
        error =
            SNESSetType(
                snes,
                SNESNEWTONLS);
    }

    SNESLineSearch line_search = nullptr;
    if (error == PETSC_SUCCESS) {
        error =
            SNESGetLineSearch(
                snes,
                &line_search);
    }
    if (error == PETSC_SUCCESS) {
        error =
            SNESLineSearchSetType(
                line_search,
                SNESLINESEARCHBT);
    }
    if (error == PETSC_SUCCESS &&
        evaluator.step_precheck !=
            nullptr) {
        error =
            SNESLineSearchSetPreCheck(
                line_search,
                line_search_precheck,
                &callback_context);
    }

    KSP ksp = nullptr;
    if (error == PETSC_SUCCESS) {
        error =
            SNESGetKSP(
                snes,
                &ksp);
    }
    if (error == PETSC_SUCCESS) {
        error =
            KSPSetType(
                ksp,
                KSPGMRES);
    }

    PC pc = nullptr;
    constexpr PetscInt asm_overlap = 1;
    if (error == PETSC_SUCCESS) {
        error =
            KSPGetPC(
                ksp,
                &pc);
    }
    if (error == PETSC_SUCCESS) {
        error =
            PCSetType(
                pc,
                PCASM);
    }
    if (error == PETSC_SUCCESS) {
        error =
            PCASMSetType(
                pc,
                PC_ASM_RESTRICT);
    }
    if (error == PETSC_SUCCESS) {
        error =
            PCASMSetOverlap(
                pc,
                asm_overlap);
    }

    constexpr const char* asm_prefix =
        "mpmc_variable_cardinality_";
    constexpr const char* sub_pc_option =
        "-mpmc_variable_cardinality_sub_pc_type";
    constexpr const char* sub_pc_reorder_option =
        "-mpmc_variable_cardinality_sub_pc_factor_nonzeros_along_diagonal";
    constexpr const char* sub_pc_shift_type_option =
        "-mpmc_variable_cardinality_sub_pc_factor_shift_type";
    bool sub_pc_option_installed = false;
    bool sub_pc_reorder_installed = false;
    bool sub_pc_shift_type_installed = false;

    if (error == PETSC_SUCCESS) {
        error =
            PCSetOptionsPrefix(
                pc,
                asm_prefix);
    }
    if (error == PETSC_SUCCESS) {
        error =
            PetscOptionsSetValue(
                nullptr,
                sub_pc_option,
                "lu");
        sub_pc_option_installed =
            error == PETSC_SUCCESS;
    }
    if (error == PETSC_SUCCESS) {
        const char* reorder_tolerance =
            row_scaling != nullptr
                ? "1.0e-10"
                : "0.0";
        error =
            PetscOptionsSetValue(
                nullptr,
                sub_pc_reorder_option,
                reorder_tolerance);
        sub_pc_reorder_installed =
            error == PETSC_SUCCESS;
    }
    if (error == PETSC_SUCCESS &&
        row_scaling != nullptr) {
        error =
            PetscOptionsSetValue(
                nullptr,
                sub_pc_shift_type_option,
                "nonzero");
        sub_pc_shift_type_installed =
            error == PETSC_SUCCESS;
    }
    if (error == PETSC_SUCCESS) {
        error =
            KSPSetTolerances(
                ksp,
                PetscReal{1.0e-12},
                PetscReal{1.0e-14},
                PETSC_DEFAULT,
                PetscInt{200});
    }
    if (error == PETSC_SUCCESS) {
        error =
            SNESSetTolerances(
                snes,
                PetscReal{1.0e-10},
                PetscReal{1.0e-10},
                PetscReal{1.0e-12},
                PetscInt{30},
                PetscInt{2000});
    }
    if (error == PETSC_SUCCESS) {
        error =
            SNESSetErrorIfNotConverged(
                snes,
                PETSC_FALSE);
    }
    if (error == PETSC_SUCCESS) {
        error =
            SNESSetCheckJacobianDomainError(
                snes,
                PETSC_TRUE);
    }

    const auto clear_options =
        [&]() {
            PetscErrorCode first =
                PETSC_SUCCESS;
            if (sub_pc_shift_type_installed) {
                first =
                    PetscOptionsClearValue(
                        nullptr,
                        sub_pc_shift_type_option);
                sub_pc_shift_type_installed =
                    false;
            }
            if (sub_pc_reorder_installed) {
                const PetscErrorCode current =
                    PetscOptionsClearValue(
                        nullptr,
                        sub_pc_reorder_option);
                sub_pc_reorder_installed =
                    false;
                if (first == PETSC_SUCCESS &&
                    current != PETSC_SUCCESS) {
                    first = current;
                }
            }
            if (sub_pc_option_installed) {
                const PetscErrorCode current =
                    PetscOptionsClearValue(
                        nullptr,
                        sub_pc_option);
                sub_pc_option_installed =
                    false;
                if (first == PETSC_SUCCESS &&
                    current != PETSC_SUCCESS) {
                    first = current;
                }
            }
            return first;
        };

    if (error != PETSC_SUCCESS) {
        (void)clear_options();
        cleanup();
        return error;
    }

    error =
        SNESSolve(
            snes,
            nullptr,
            solved_state);
    const PetscErrorCode clear_error =
        clear_options();
    if (error == PETSC_SUCCESS &&
        clear_error != PETSC_SUCCESS) {
        error =
            clear_error;
    }
    if (error != PETSC_SUCCESS) {
        cleanup();
        return error;
    }

    SNESConvergedReason reason{
        SNES_CONVERGED_ITERATING};
    PetscInt nonlinear_iterations = -1;
    error =
        SNESGetConvergedReason(
            snes,
            &reason);
    if (error == PETSC_SUCCESS) {
        error =
            SNESGetIterationNumber(
                snes,
                &nonlinear_iterations);
    }
    if (error != PETSC_SUCCESS) {
        cleanup();
        return error;
    }
    if (static_cast<int>(reason) <= 0) {
        if (failure_diagnostics != nullptr) {
            // Diagnostic-only replay of the actual frozen linear problem. Never
            // replace the failed Newton correction or alter the SNES outcome.
            Mat frozen_matrix = nullptr;
            Vec frozen_rhs = nullptr;
            Vec failed_solution = nullptr;
            PetscErrorCode replay_error = KSPGetOperators(ksp, &frozen_matrix, nullptr);
            if (replay_error == PETSC_SUCCESS) replay_error = KSPGetRhs(ksp, &frozen_rhs);
            if (replay_error == PETSC_SUCCESS) replay_error = KSPGetSolution(ksp, &failed_solution);
            const auto print_linear_diagnostic = [&](KSP solver, const char* label) {
                Vec x = nullptr;
                Vec residual = nullptr;
                PetscInt iterations = -1;
                PetscReal rhs_norm = 0.0;
                PetscReal true_norm = std::numeric_limits<PetscReal>::quiet_NaN();
                KSPConvergedReason linear_reason = KSP_CONVERGED_ITERATING;
                PetscErrorCode diagnostic_error = KSPGetSolution(solver, &x);
                if (diagnostic_error == PETSC_SUCCESS) diagnostic_error = KSPGetIterationNumber(solver, &iterations);
                if (diagnostic_error == PETSC_SUCCESS) diagnostic_error = KSPGetConvergedReason(solver, &linear_reason);
                if (diagnostic_error == PETSC_SUCCESS) diagnostic_error = VecNorm(frozen_rhs, NORM_2, &rhs_norm);
                if (diagnostic_error == PETSC_SUCCESS) diagnostic_error = VecDuplicate(frozen_rhs, &residual);
                if (diagnostic_error == PETSC_SUCCESS) diagnostic_error = MatMult(frozen_matrix, x, residual);
                if (diagnostic_error == PETSC_SUCCESS) diagnostic_error = VecAXPY(residual, -1.0, frozen_rhs);
                if (diagnostic_error == PETSC_SUCCESS) diagnostic_error = VecNorm(residual, NORM_2, &true_norm);
                (void)PetscPrintf(comm,
                    "[frozen linear diagnostic] path=%s error=%d reason=%d iterations=%d "
                    "rhs_l2=%.17g true_residual_l2=%.17g relative_true_residual=%.17g\n",
                    label, static_cast<int>(diagnostic_error), static_cast<int>(linear_reason),
                    static_cast<int>(iterations), static_cast<double>(rhs_norm),
                    static_cast<double>(true_norm), static_cast<double>(
                        rhs_norm > 0.0 ? true_norm / rhs_norm : true_norm));
                (void)KSPView(solver, PETSC_VIEWER_STDOUT_(comm));
                if (residual != nullptr) (void)VecDestroy(&residual);
            };
            if (replay_error == PETSC_SUCCESS && frozen_matrix != nullptr &&
                frozen_rhs != nullptr && failed_solution != nullptr) {
                print_linear_diagnostic(ksp, "actual-gmres-asm");
                PetscMPIInt comm_size = 0;
                if (MPI_Comm_size(comm, &comm_size) == MPI_SUCCESS && comm_size == 1) {
                    KSP direct = nullptr;
                    Vec direct_solution = nullptr;
                    PC direct_pc = nullptr;
                    // Return diagnostic LU errors instead of aborting the original failure report.
                    replay_error = PetscPushErrorHandler(PetscReturnErrorHandler, nullptr);
                    const bool handler_installed = replay_error == PETSC_SUCCESS;
                    if (replay_error == PETSC_SUCCESS) replay_error = KSPCreate(comm, &direct);
                    if (replay_error == PETSC_SUCCESS) replay_error = KSPSetOperators(direct, frozen_matrix, frozen_matrix);
                    if (replay_error == PETSC_SUCCESS) replay_error = KSPSetType(direct, KSPPREONLY);
                    if (replay_error == PETSC_SUCCESS) replay_error = KSPGetPC(direct, &direct_pc);
                    if (replay_error == PETSC_SUCCESS) replay_error = PCSetType(direct_pc, PCLU);
                    if (replay_error == PETSC_SUCCESS) replay_error = PCFactorReorderForNonzeroDiagonal(direct_pc, 1.0e-10);
                    if (replay_error == PETSC_SUCCESS) replay_error = PCFactorSetShiftType(direct_pc, MAT_SHIFT_NONE);
                    if (replay_error == PETSC_SUCCESS) replay_error = VecDuplicate(failed_solution, &direct_solution);
                    if (replay_error == PETSC_SUCCESS) replay_error = VecSet(direct_solution, 0.0);
                    if (replay_error == PETSC_SUCCESS) replay_error = KSPSolve(direct, frozen_rhs, direct_solution);
                    (void)PetscPrintf(comm, "[frozen linear diagnostic] path=direct-lu solve_error=%d\n",
                        static_cast<int>(replay_error));
                    if (replay_error == PETSC_SUCCESS) {
                        print_linear_diagnostic(direct, "direct-lu-no-shift");
                    } else if (direct != nullptr) {
                        (void)KSPView(direct, PETSC_VIEWER_STDOUT_(comm));
                    }
                    if (direct != nullptr) (void)KSPDestroy(&direct);
                    if (direct_solution != nullptr) (void)VecDestroy(&direct_solution);
                    if (handler_installed) (void)PetscPopErrorHandler();
                }
            } else {
                (void)PetscPrintf(comm, "[frozen linear diagnostic] unavailable error=%d\n",
                    static_cast<int>(replay_error));
            }
            KSPConvergedReason ksp_reason{
                KSP_CONVERGED_ITERATING};
            PCFailedReason pc_reason{};
            KSPConvergedReason sub_ksp_reason{
                KSP_CONVERGED_ITERATING};
            PCFailedReason sub_pc_reason{};
            PetscReal function_norm = 0.0;

            const PetscErrorCode ksp_reason_error =
                KSPGetConvergedReason(
                    ksp,
                    &ksp_reason);
            const PetscErrorCode pc_reason_error =
                PCGetFailedReason(
                    pc,
                    &pc_reason);
            const PetscErrorCode norm_error =
                SNESGetFunctionNorm(
                    snes,
                    &function_norm);

            PetscErrorCode sub_reason_error =
                PETSC_SUCCESS;
            PetscInt sub_count = 0;
            KSP* sub_ksp = nullptr;
            if (pc_reason_error == PETSC_SUCCESS) {
                sub_reason_error =
                    PCASMGetSubKSP(
                        pc,
                        &sub_count,
                        nullptr,
                        &sub_ksp);
            }
            if (sub_reason_error == PETSC_SUCCESS &&
                sub_count > 0 &&
                sub_ksp != nullptr) {
                for (PetscInt index = 0;
                     index < sub_count;
                     ++index) {
                    KSPConvergedReason candidate_ksp{
                        KSP_CONVERGED_ITERATING};
                    PC sub_pc = nullptr;
                    PCFailedReason candidate_pc{};
                    if (KSPGetConvergedReason(
                            sub_ksp[index],
                            &candidate_ksp) !=
                            PETSC_SUCCESS ||
                        KSPGetPC(
                            sub_ksp[index],
                            &sub_pc) !=
                            PETSC_SUCCESS ||
                        sub_pc == nullptr ||
                        PCGetFailedReason(
                            sub_pc,
                            &candidate_pc) !=
                            PETSC_SUCCESS) {
                        sub_reason_error =
                            PETSC_ERR_LIB;
                        break;
                    }
                    if (static_cast<int>(
                            candidate_ksp) < 0) {
                        sub_ksp_reason =
                            candidate_ksp;
                    }
                    if (static_cast<int>(
                            candidate_pc) != 0) {
                        sub_pc_reason =
                            candidate_pc;
                    }
                }
            }

            if (ksp_reason_error == PETSC_SUCCESS &&
                pc_reason_error == PETSC_SUCCESS &&
                norm_error == PETSC_SUCCESS &&
                sub_reason_error == PETSC_SUCCESS &&
                std::isfinite(
                    static_cast<double>(
                        function_norm)) &&
                function_norm >= 0.0) {
                failure_diagnostics->emplace(
                    NaturalVariableSnesFailureDiagnostics3D{
                        reason,
                        ksp_reason,
                        static_cast<int>(
                            pc_reason),
                        sub_ksp_reason,
                        static_cast<int>(
                            sub_pc_reason),
                        nonlinear_iterations,
                        callback_context.function_evaluations,
                        callback_context.jacobian_evaluations,
                        callback_context.function_domain_errors,
                        callback_context.jacobian_domain_errors,
                        callback_context.line_search_prechecks,
                        callback_context.line_search_direction_changes,
                        static_cast<double>(
                            function_norm)});
            }
        }
        cleanup();
        return PETSC_ERR_NOT_CONVERGED;
    }

    error =
        SNESComputeFunction(
            snes,
            solved_state,
            residual);
    PetscReal final_norm = 0.0;
    if (error == PETSC_SUCCESS) {
        error =
            VecNorm(
                residual,
                NORM_2,
                &final_norm);
    }

    const char* snes_name = nullptr;
    const char* line_search_name = nullptr;
    const char* ksp_name = nullptr;
    const char* pc_name = nullptr;
    if (error == PETSC_SUCCESS) {
        error =
            SNESGetType(
                snes,
                &snes_name);
    }
    if (error == PETSC_SUCCESS) {
        error =
            SNESLineSearchGetType(
                line_search,
                &line_search_name);
    }
    if (error == PETSC_SUCCESS) {
        error =
            KSPGetType(
                ksp,
                &ksp_name);
    }
    if (error == PETSC_SUCCESS) {
        error =
            PCGetType(
                pc,
                &pc_name);
    }
    if (error != PETSC_SUCCESS ||
        snes_name == nullptr ||
        line_search_name == nullptr ||
        ksp_name == nullptr ||
        pc_name == nullptr) {
        cleanup();
        return error != PETSC_SUCCESS
            ? error
            : PETSC_ERR_PLIB;
    }

    const PetscScalar* solved_values =
        nullptr;
    error =
        VecGetArrayRead(
            solved_state,
            &solved_values);
    if (error != PETSC_SUCCESS) {
        cleanup();
        return error;
    }

    std::vector<
        VariableCardinalityNaturalVariableSnesSolutionEntry3D>
        entries(
            static_cast<std::size_t>(
                numbering
                    .petsc_local_owned_scalar_count()));
    bool entries_ok = true;
    for (const auto& cell :
         numbering.cells()) {
        if (cell.owner_rank !=
            numbering.local_rank()) {
            continue;
        }
        for (std::size_t slot = 0U;
             slot < cell.scalar_count;
             ++slot) {
            const PetscInt global =
                cell.petsc_global_scalar_start +
                static_cast<PetscInt>(
                    slot);
            const PetscInt local =
                global -
                numbering
                    .petsc_owned_scalar_start();
            if (local < 0 ||
                local >=
                    numbering
                        .petsc_local_owned_scalar_count()) {
                entries_ok = false;
                continue;
            }
            entries[
                static_cast<std::size_t>(
                    local)] =
                {
                    global,
                    cell.cell_global,
                    cell.phase_count,
                    slot,
                    static_cast<double>(
                        PetscRealPart(
                            solved_values[local]))};
        }
    }
    const PetscErrorCode restore_error =
        VecRestoreArrayRead(
            solved_state,
            &solved_values);
    if (restore_error !=
        PETSC_SUCCESS) {
        cleanup();
        return restore_error;
    }
    if (!entries_ok) {
        cleanup();
        return PETSC_ERR_PLIB;
    }
    for (std::size_t index = 0U;
         index < entries.size();
         ++index) {
        if (entries[index]
                .petsc_global_scalar !=
                numbering
                    .petsc_owned_scalar_start() +
                    static_cast<PetscInt>(
                        index) ||
            !std::isfinite(
                entries[index].value)) {
            cleanup();
            return PETSC_ERR_PLIB;
        }
    }

    VariableCardinalityNaturalVariableSnesSolveReport3D
        completed;
    completed.converged_reason =
        reason;
    completed.nonlinear_iterations =
        nonlinear_iterations;
    completed.function_evaluations =
        callback_context
            .function_evaluations;
    completed.jacobian_evaluations =
        callback_context
            .jacobian_evaluations;
    completed.line_search_prechecks =
        callback_context
            .line_search_prechecks;
    completed.line_search_direction_changes =
        callback_context
            .line_search_direction_changes;
    completed.final_function_l2_norm =
        static_cast<double>(
            final_norm);
    completed.snes_type =
        snes_name;
    completed.line_search_type =
        line_search_name;
    completed.ksp_type =
        ksp_name;
    completed.pc_type =
        pc_name;
    completed.asm_overlap =
        asm_overlap;
    completed.locally_owned_solution =
        std::move(entries);

    error =
        SNESDestroy(
            &snes);
    if (error == PETSC_SUCCESS) {
        error =
            MatDestroy(
                &jacobian);
    }
    if (error == PETSC_SUCCESS) {
        error =
            VecDestroy(
                &residual);
    }
    if (error != PETSC_SUCCESS) {
        if (solved_state != nullptr) {
            (void)VecDestroy(
                &solved_state);
        }
        return error;
    }

    *solution = solved_state;
    solved_state = nullptr;
    report->emplace(
        std::move(completed));
    return PETSC_SUCCESS;
}

} // namespace mpmc::flow_discretization_petsc

#endif // MPMC_FLOW_DISCRETIZATION_PETSC_VARIABLE_CARDINALITY_SNES_ASSEMBLY_HPP
