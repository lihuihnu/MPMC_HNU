#ifndef MPMC_FLOW_DISCRETIZATION_PETSC_COMPLETE_NATURAL_VARIABLE_PETSC_MATERIALIZATION_HPP
#define MPMC_FLOW_DISCRETIZATION_PETSC_COMPLETE_NATURAL_VARIABLE_PETSC_MATERIALIZATION_HPP

#include <mpmc/flow_discretization_petsc/complete_natural_variable_assembly_snapshot.hpp>

#include <petscmat.h>
#include <petscvec.h>

#include <cstddef>
#include <cstdint>
#include <limits>
#include <span>
#include <stdexcept>
#include <string_view>
#include <vector>

namespace mpmc::flow_discretization_petsc {

inline constexpr std::string_view
    complete_natural_variable_petsc_materialization_convention =
        "flow_discretization_petsc/complete-natural-variable-petsc-materialization/v1";

namespace complete_petsc_materialization_detail {

[[nodiscard]] inline bool
checked_nnz_product(
    PetscInt cell_nnz,
    std::size_t block_width,
    PetscInt* output) {
    if (output == nullptr ||
        cell_nnz < 0 ||
        block_width == 0U ||
        block_width >
            static_cast<std::size_t>(
                std::numeric_limits<PetscInt>::max())) {
        return false;
    }
    const PetscInt width =
        static_cast<PetscInt>(
            block_width);
    if (cell_nnz >
        std::numeric_limits<PetscInt>::max() /
            width) {
        return false;
    }
    *output =
        cell_nnz * width;
    return true;
}

[[nodiscard]] inline PetscErrorCode
insert_residual_values(
    Vec residual,
    std::span<
        const CompleteNaturalVariableResidualEntry3D>
        entries) {
    if (residual == nullptr) {
        return PETSC_ERR_ARG_NULL;
    }
    if (entries.empty()) {
        return PETSC_SUCCESS;
    }
    if (entries.size() >
        static_cast<std::size_t>(
            std::numeric_limits<PetscInt>::max())) {
        return PETSC_ERR_ARG_OUTOFRANGE;
    }

    std::vector<PetscInt> rows;
    std::vector<PetscScalar> values;
    rows.reserve(entries.size());
    values.reserve(entries.size());
    for (const auto& entry : entries) {
        rows.push_back(
            entry.petsc_global_row);
        values.push_back(
            static_cast<PetscScalar>(
                entry.native_value));
    }

    return VecSetValues(
        residual,
        static_cast<PetscInt>(
            rows.size()),
        rows.data(),
        values.data(),
        INSERT_VALUES);
}

[[nodiscard]] inline PetscErrorCode
insert_jacobian_values(
    Mat jacobian,
    std::span<
        const CompleteNaturalVariableJacobianEntry3D>
        entries,
    bool structural_zero_seed) {
    if (jacobian == nullptr) {
        return PETSC_ERR_ARG_NULL;
    }

    for (const auto& entry : entries) {
        const PetscInt row =
            entry.petsc_global_row;
        const PetscInt column =
            entry.petsc_global_column;
        const PetscScalar value =
            structural_zero_seed
                ? PetscScalar{0.0}
                : static_cast<PetscScalar>(
                      entry.value);
        const PetscErrorCode error =
            MatSetValues(
                jacobian,
                1,
                &row,
                1,
                &column,
                &value,
                INSERT_VALUES);
        if (error != PETSC_SUCCESS) {
            return error;
        }
    }
    return PETSC_SUCCESS;
}

} // namespace complete_petsc_materialization_detail

/// Materialize and assemble the complete natural-variable residual/Jacobian.
///
/// Ownership:
/// - on success, *residual and *jacobian are owned by the caller and must be
///   destroyed with VecDestroy()/MatDestroy().
/// - on failure, both output handles remain null.
///
/// The residual vector intentionally contains heterogeneous native equation
/// units. No scaling, norm definition, linear solve or nonlinear solve is
/// performed here.
inline PetscErrorCode
materialize_complete_natural_variable_petsc_system_3d(
    MPI_Comm comm,
    const CompleteNaturalVariableAssemblySnapshot3D&
        snapshot,
    const mpmc::discretization_petsc::
        PetscMpiAijSymbolicPreallocation3D&
            cell_bridge,
    Vec* residual,
    Mat* jacobian) {
    using namespace complete_petsc_materialization_detail;

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
        residual == nullptr ||
                jacobian == nullptr
            ? PETSC_ERR_ARG_NULL
            : PETSC_SUCCESS;
    if (local_error == PETSC_SUCCESS &&
        (*residual != nullptr ||
         *jacobian != nullptr)) {
        local_error =
            PETSC_ERR_ARG_WRONGSTATE;
    }

    PetscErrorCode error =
        detail::collective_error(
            comm,
            local_error);
    if (error != PETSC_SUCCESS) {
        return error;
    }

    *residual = nullptr;
    *jacobian = nullptr;

    const std::size_t q =
        snapshot.natural_variable_count();
    PetscInt expected_scalar_start = -1;
    PetscInt expected_scalar_end = -1;
    PetscInt expected_scalar_count = -1;
    PetscInt expected_local_scalar_rows = -1;

    std::vector<PetscInt>
        scalar_diagonal_nnz;
    std::vector<PetscInt>
        scalar_off_diagonal_nnz;

    try {
        if (snapshot.local_rank() !=
                cell_bridge.local_rank() ||
            snapshot.rank_count() !=
                cell_bridge.rank_count() ||
            snapshot.local_rank().value() !=
                static_cast<std::uint32_t>(
                    mpi_rank) ||
            snapshot.rank_count() !=
                static_cast<std::uint32_t>(
                    mpi_size) ||
            q == 0U) {
            throw std::invalid_argument(
                "complete PETSc materialization rank/block metadata mismatch");
        }

        if (!global_component_assembly_mapping_detail::
                checked_scalar_count(
                    cell_bridge.global_row_start(),
                    q,
                    &expected_scalar_start) ||
            !global_component_assembly_mapping_detail::
                checked_scalar_count(
                    cell_bridge.global_row_end(),
                    q,
                    &expected_scalar_end) ||
            !global_component_assembly_mapping_detail::
                checked_scalar_count(
                    cell_bridge.global_row_count(),
                    q,
                    &expected_scalar_count) ||
            !global_component_assembly_mapping_detail::
                checked_scalar_count(
                    cell_bridge.local_owned_row_count(),
                    q,
                    &expected_local_scalar_rows) ||
            expected_scalar_start !=
                snapshot.petsc_scalar_row_start() ||
            expected_scalar_end !=
                snapshot.petsc_scalar_row_end() ||
            expected_scalar_count !=
                snapshot.petsc_scalar_row_count() ||
            expected_local_scalar_rows !=
                expected_scalar_end -
                    expected_scalar_start) {
            throw std::invalid_argument(
                "complete PETSc materialization scalar ownership mismatch");
        }

        if (expected_local_scalar_rows < 0) {
            throw std::invalid_argument(
                "complete PETSc materialization local scalar row count is negative");
        }

        const auto diagonal =
            cell_bridge.diagonal_nnz();
        const auto off_diagonal =
            cell_bridge.off_diagonal_nnz();
        if (diagonal.size() !=
                off_diagonal.size() ||
            diagonal.size() !=
                cell_bridge
                    .owned_cells_in_petsc_row_order()
                    .size()) {
            throw std::invalid_argument(
                "complete PETSc materialization cell preallocation arrays disagree");
        }

        const std::size_t local_scalar_count =
            static_cast<std::size_t>(
                expected_local_scalar_rows);
        scalar_diagonal_nnz.reserve(
            local_scalar_count);
        scalar_off_diagonal_nnz.reserve(
            local_scalar_count);

        for (std::size_t cell = 0U;
             cell < diagonal.size();
             ++cell) {
            PetscInt diagonal_scalar_nnz = 0;
            PetscInt off_diagonal_scalar_nnz = 0;
            if (!checked_nnz_product(
                    diagonal[cell],
                    q,
                    &diagonal_scalar_nnz) ||
                !checked_nnz_product(
                    off_diagonal[cell],
                    q,
                    &off_diagonal_scalar_nnz) ||
                diagonal_scalar_nnz <= 0 ||
                off_diagonal_scalar_nnz < 0) {
                throw std::length_error(
                    "complete PETSc materialization scalar preallocation overflow");
            }

            for (std::size_t slot = 0U;
                 slot < q;
                 ++slot) {
                scalar_diagonal_nnz.push_back(
                    diagonal_scalar_nnz);
                scalar_off_diagonal_nnz.push_back(
                    off_diagonal_scalar_nnz);
            }
        }

        if (scalar_diagonal_nnz.size() !=
                local_scalar_count ||
            scalar_off_diagonal_nnz.size() !=
                local_scalar_count) {
            throw std::runtime_error(
                "complete PETSc materialization scalar preallocation cardinality mismatch");
        }
    } catch (...) {
        local_error =
            PETSC_ERR_ARG_INCOMP;
    }

    error =
        detail::collective_error(
            comm,
            local_error);
    if (error != PETSC_SUCCESS) {
        return error;
    }

    Vec created_residual = nullptr;
    Mat created_jacobian = nullptr;

    const auto cleanup =
        [&]() {
            PetscErrorCode first =
                PETSC_SUCCESS;
            if (created_residual !=
                nullptr) {
                const PetscErrorCode destroy =
                    VecDestroy(
                        &created_residual);
                if (first ==
                        PETSC_SUCCESS &&
                    destroy !=
                        PETSC_SUCCESS) {
                    first = destroy;
                }
            }
            if (created_jacobian !=
                nullptr) {
                const PetscErrorCode destroy =
                    MatDestroy(
                        &created_jacobian);
                if (first ==
                        PETSC_SUCCESS &&
                    destroy !=
                        PETSC_SUCCESS) {
                    first = destroy;
                }
            }
            return first;
        };

    error =
        VecCreateMPI(
            comm,
            expected_local_scalar_rows,
            expected_scalar_count,
            &created_residual);
    if (error != PETSC_SUCCESS) {
        (void)cleanup();
        return error;
    }

    error =
        MatCreate(
            comm,
            &created_jacobian);
    if (error == PETSC_SUCCESS) {
        error =
            MatSetSizes(
                created_jacobian,
                expected_local_scalar_rows,
                expected_local_scalar_rows,
                expected_scalar_count,
                expected_scalar_count);
    }
    if (error == PETSC_SUCCESS) {
        error =
            MatSetType(
                created_jacobian,
                MATMPIAIJ);
    }
    if (error == PETSC_SUCCESS) {
        error =
            MatMPIAIJSetPreallocation(
                created_jacobian,
                0,
                scalar_diagonal_nnz.empty()
                    ? nullptr
                    : scalar_diagonal_nnz.data(),
                0,
                scalar_off_diagonal_nnz.empty()
                    ? nullptr
                    : scalar_off_diagonal_nnz.data());
    }
    if (error == PETSC_SUCCESS) {
        error =
            MatSetOption(
                created_jacobian,
                MAT_IGNORE_ZERO_ENTRIES,
                PETSC_FALSE);
    }
    if (error == PETSC_SUCCESS) {
        error =
            MatSetOption(
                created_jacobian,
                MAT_NEW_NONZERO_ALLOCATION_ERR,
                PETSC_TRUE);
    }
    if (error != PETSC_SUCCESS) {
        (void)cleanup();
        return error;
    }

    // Local insertions are checked collectively before any assembly collective.
    local_error =
        insert_residual_values(
            created_residual,
            snapshot.residual_entries());
    if (local_error ==
        PETSC_SUCCESS) {
        local_error =
            insert_jacobian_values(
                created_jacobian,
                snapshot.jacobian_entries(),
                true);
    }

    error =
        detail::collective_error(
            comm,
            local_error);
    if (error != PETSC_SUCCESS) {
        (void)cleanup();
        return error;
    }

    error =
        VecAssemblyBegin(
            created_residual);
    if (error == PETSC_SUCCESS) {
        error =
            VecAssemblyEnd(
                created_residual);
    }
    if (error == PETSC_SUCCESS) {
        error =
            MatAssemblyBegin(
                created_jacobian,
                MAT_FINAL_ASSEMBLY);
    }
    if (error == PETSC_SUCCESS) {
        error =
            MatAssemblyEnd(
                created_jacobian,
                MAT_FINAL_ASSEMBLY);
    }
    if (error != PETSC_SUCCESS) {
        (void)cleanup();
        return error;
    }

    // Freeze the actual equation-aware locations seeded from the complete
    // snapshot, not the wider cell-level capacity pattern.
    error =
        MatSetOption(
            created_jacobian,
            MAT_NEW_NONZERO_LOCATION_ERR,
            PETSC_TRUE);
    if (error != PETSC_SUCCESS) {
        (void)cleanup();
        return error;
    }

    local_error =
        insert_jacobian_values(
            created_jacobian,
            snapshot.jacobian_entries(),
            false);
    error =
        detail::collective_error(
            comm,
            local_error);
    if (error != PETSC_SUCCESS) {
        (void)cleanup();
        return error;
    }

    error =
        MatAssemblyBegin(
            created_jacobian,
            MAT_FINAL_ASSEMBLY);
    if (error == PETSC_SUCCESS) {
        error =
            MatAssemblyEnd(
                created_jacobian,
                MAT_FINAL_ASSEMBLY);
    }
    if (error != PETSC_SUCCESS) {
        (void)cleanup();
        return error;
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
    PetscBool is_mpiaij =
        PETSC_FALSE;

    error =
        VecGetLocalSize(
            created_residual,
            &vec_local);
    if (error == PETSC_SUCCESS) {
        error =
            VecGetSize(
                created_residual,
                &vec_global);
    }
    if (error == PETSC_SUCCESS) {
        error =
            VecGetOwnershipRange(
                created_residual,
                &vec_start,
                &vec_end);
    }
    if (error == PETSC_SUCCESS) {
        error =
            MatGetLocalSize(
                created_jacobian,
                &mat_local_rows,
                &mat_local_columns);
    }
    if (error == PETSC_SUCCESS) {
        error =
            MatGetSize(
                created_jacobian,
                &mat_global_rows,
                &mat_global_columns);
    }
    if (error == PETSC_SUCCESS) {
        error =
            MatGetOwnershipRange(
                created_jacobian,
                &mat_start,
                &mat_end);
    }
    if (error == PETSC_SUCCESS) {
        error =
            MatGetOwnershipRangeColumn(
                created_jacobian,
                &mat_column_start,
                &mat_column_end);
    }
    if (error == PETSC_SUCCESS) {
        error =
            PetscObjectTypeCompare(
                reinterpret_cast<PetscObject>(
                    created_jacobian),
                MATMPIAIJ,
                &is_mpiaij);
    }
    if (error != PETSC_SUCCESS) {
        (void)cleanup();
        return error;
    }

    if (vec_local !=
            expected_local_scalar_rows ||
        vec_global !=
            expected_scalar_count ||
        vec_start !=
            expected_scalar_start ||
        vec_end !=
            expected_scalar_end ||
        mat_local_rows !=
            expected_local_scalar_rows ||
        mat_local_columns !=
            expected_local_scalar_rows ||
        mat_global_rows !=
            expected_scalar_count ||
        mat_global_columns !=
            expected_scalar_count ||
        mat_start !=
            expected_scalar_start ||
        mat_end !=
            expected_scalar_end ||
        mat_column_start !=
            expected_scalar_start ||
        mat_column_end !=
            expected_scalar_end ||
        is_mpiaij !=
            PETSC_TRUE) {
        (void)cleanup();
        return PETSC_ERR_ARG_INCOMP;
    }

    *residual =
        created_residual;
    *jacobian =
        created_jacobian;
    return PETSC_SUCCESS;
}

} // namespace mpmc::flow_discretization_petsc

#endif // MPMC_FLOW_DISCRETIZATION_PETSC_COMPLETE_NATURAL_VARIABLE_PETSC_MATERIALIZATION_HPP
