#ifndef MPMC_FLOW_DISCRETIZATION_PETSC_NATURAL_VARIABLE_NEWTON_LINEAR_SYSTEM_HPP
#define MPMC_FLOW_DISCRETIZATION_PETSC_NATURAL_VARIABLE_NEWTON_LINEAR_SYSTEM_HPP

#include <mpmc/flow_discretization_petsc/complete_natural_variable_assembly_snapshot.hpp>

#include <petscksp.h>

#include <algorithm>
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
    frozen_natural_variable_newton_linear_system_convention =
        "flow_discretization_petsc/frozen-natural-variable-newton-linear-system/v1";

struct NaturalVariableLinearCorrectionEntry3D {
    PetscInt petsc_global_scalar{-1};
    mpmc::mesh::GlobalDofIndex mesh_global_dof{
        mpmc::mesh::GlobalDofIndex::value_type{0}};
    mpmc::mesh::GlobalEntityId cell_global{
        mpmc::mesh::GlobalEntityId::value_type{0}};
    std::size_t natural_variable_slot{};
    double value{};
};

class FrozenNaturalVariableNewtonLinearCorrection3D {
public:
    static constexpr std::string_view convention =
        frozen_natural_variable_newton_linear_system_convention;

    FrozenNaturalVariableNewtonLinearCorrection3D(
        mpmc::mesh::PartitionRank local_rank,
        std::uint32_t rank_count,
        std::string natural_variable_id,
        std::size_t component_count,
        std::size_t natural_variable_count,
        PetscInt petsc_scalar_row_start,
        PetscInt petsc_scalar_row_end,
        PetscInt petsc_scalar_row_count,
        std::string ksp_type,
        std::string pc_type,
        PetscInt iteration_count,
        KSPConvergedReason converged_reason,
        double residual_l2_norm,
        double rhs_l2_norm,
        double linear_residual_l2_norm,
        std::vector<NaturalVariableLinearCorrectionEntry3D>
            locally_owned_correction)
        : local_rank_(local_rank),
          rank_count_(rank_count),
          natural_variable_id_(
              std::move(natural_variable_id)),
          component_count_(component_count),
          natural_variable_count_(
              natural_variable_count),
          petsc_scalar_row_start_(
              petsc_scalar_row_start),
          petsc_scalar_row_end_(
              petsc_scalar_row_end),
          petsc_scalar_row_count_(
              petsc_scalar_row_count),
          ksp_type_(std::move(ksp_type)),
          pc_type_(std::move(pc_type)),
          iteration_count_(iteration_count),
          converged_reason_(converged_reason),
          residual_l2_norm_(
              residual_l2_norm),
          rhs_l2_norm_(
              rhs_l2_norm),
          linear_residual_l2_norm_(
              linear_residual_l2_norm),
          locally_owned_correction_(
              std::move(
                  locally_owned_correction)) {
        validate();
    }

    [[nodiscard]] mpmc::mesh::PartitionRank
    local_rank() const noexcept {
        return local_rank_;
    }

    [[nodiscard]] std::uint32_t
    rank_count() const noexcept {
        return rank_count_;
    }

    [[nodiscard]] std::string_view
    natural_variable_id() const noexcept {
        return natural_variable_id_;
    }

    [[nodiscard]] std::size_t
    component_count() const noexcept {
        return component_count_;
    }

    [[nodiscard]] std::size_t
    natural_variable_count() const noexcept {
        return natural_variable_count_;
    }

    [[nodiscard]] PetscInt
    petsc_scalar_row_start() const noexcept {
        return petsc_scalar_row_start_;
    }

    [[nodiscard]] PetscInt
    petsc_scalar_row_end() const noexcept {
        return petsc_scalar_row_end_;
    }

    [[nodiscard]] PetscInt
    petsc_scalar_row_count() const noexcept {
        return petsc_scalar_row_count_;
    }

    [[nodiscard]] std::string_view
    ksp_type() const noexcept {
        return ksp_type_;
    }

    [[nodiscard]] std::string_view
    pc_type() const noexcept {
        return pc_type_;
    }

    [[nodiscard]] PetscInt
    iteration_count() const noexcept {
        return iteration_count_;
    }

    [[nodiscard]] KSPConvergedReason
    converged_reason() const noexcept {
        return converged_reason_;
    }

    [[nodiscard]] double
    residual_l2_norm() const noexcept {
        return residual_l2_norm_;
    }

    [[nodiscard]] double
    rhs_l2_norm() const noexcept {
        return rhs_l2_norm_;
    }

    /// Unscaled algebraic norm ||J*delta_q + R||_2.
    ///
    /// The natural-variable residual contains rows with different native
    /// physical dimensions. This norm is therefore an algebraic linear-solve
    /// verification quantity only; it is not a physical Newton convergence
    /// metric.
    [[nodiscard]] double
    linear_residual_l2_norm() const noexcept {
        return linear_residual_l2_norm_;
    }

    [[nodiscard]] std::span<
        const NaturalVariableLinearCorrectionEntry3D>
    locally_owned_correction() const noexcept {
        return locally_owned_correction_;
    }

private:
    void validate() const {
        if (rank_count_ == 0U ||
            local_rank_.value() >=
                rank_count_ ||
            natural_variable_id_.empty() ||
            component_count_ < 2U ||
            natural_variable_count_ !=
                mpmc::flow::fixed_three_phase_count *
                    component_count_ +
                    1U ||
            petsc_scalar_row_start_ < 0 ||
            petsc_scalar_row_end_ <
                petsc_scalar_row_start_ ||
            petsc_scalar_row_count_ <
                petsc_scalar_row_end_ ||
            ksp_type_.empty() ||
            pc_type_.empty() ||
            iteration_count_ < 0 ||
            static_cast<int>(
                converged_reason_) <= 0 ||
            !std::isfinite(
                residual_l2_norm_) ||
            residual_l2_norm_ < 0.0 ||
            !std::isfinite(
                rhs_l2_norm_) ||
            rhs_l2_norm_ < 0.0 ||
            !std::isfinite(
                linear_residual_l2_norm_) ||
            linear_residual_l2_norm_ < 0.0 ||
            locally_owned_correction_.size() !=
                static_cast<std::size_t>(
                    petsc_scalar_row_end_ -
                    petsc_scalar_row_start_)) {
            throw std::invalid_argument(
                "mpmc::flow_discretization_petsc: malformed frozen Newton linear-correction report");
        }

        for (std::size_t index = 0U;
             index <
             locally_owned_correction_.size();
             ++index) {
            const auto& entry =
                locally_owned_correction_[
                    index];
            const PetscInt expected =
                petsc_scalar_row_start_ +
                static_cast<PetscInt>(
                    index);
            if (entry.petsc_global_scalar !=
                    expected ||
                entry.natural_variable_slot !=
                    static_cast<std::size_t>(
                        expected %
                        static_cast<PetscInt>(
                            natural_variable_count_)) ||
                !std::isfinite(entry.value)) {
                throw std::invalid_argument(
                    "mpmc::flow_discretization_petsc: frozen Newton correction numbering/value mismatch");
            }
        }
    }

    mpmc::mesh::PartitionRank local_rank_;
    std::uint32_t rank_count_{};
    std::string natural_variable_id_;
    std::size_t component_count_{};
    std::size_t natural_variable_count_{};
    PetscInt petsc_scalar_row_start_{};
    PetscInt petsc_scalar_row_end_{};
    PetscInt petsc_scalar_row_count_{};
    std::string ksp_type_;
    std::string pc_type_;
    PetscInt iteration_count_{};
    KSPConvergedReason converged_reason_{
        KSP_CONVERGED_ITERATING};
    double residual_l2_norm_{};
    double rhs_l2_norm_{};
    double linear_residual_l2_norm_{};
    std::vector<NaturalVariableLinearCorrectionEntry3D>
        locally_owned_correction_;
};

namespace frozen_newton_linear_detail {

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

[[nodiscard]] inline PetscErrorCode
validate_petsc_linear_objects(
    MPI_Comm comm,
    const CompleteNaturalVariableAssemblySnapshot3D&
        numbering,
    Mat jacobian,
    Vec residual) {
    if (jacobian == nullptr ||
        residual == nullptr) {
        return PETSC_ERR_ARG_NULL;
    }

    if (!communicators_are_compatible(
            comm,
            PetscObjectComm(
                reinterpret_cast<PetscObject>(
                    jacobian))) ||
        !communicators_are_compatible(
            comm,
            PetscObjectComm(
                reinterpret_cast<PetscObject>(
                    residual)))) {
        return PETSC_ERR_ARG_INCOMP;
    }

    PetscBool assembled =
        PETSC_FALSE;
    PetscErrorCode error =
        MatAssembled(
            jacobian,
            &assembled);
    if (error != PETSC_SUCCESS) {
        return error;
    }
    if (assembled != PETSC_TRUE) {
        return PETSC_ERR_ARG_WRONGSTATE;
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
            residual,
            &vec_local);
    if (error == PETSC_SUCCESS) {
        error =
            VecGetSize(
                residual,
                &vec_global);
    }
    if (error == PETSC_SUCCESS) {
        error =
            VecGetOwnershipRange(
                residual,
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
        numbering.petsc_scalar_row_end() -
        numbering.petsc_scalar_row_start();
    if (expected_local < 0 ||
        vec_local != expected_local ||
        vec_global !=
            numbering.petsc_scalar_row_count() ||
        vec_start !=
            numbering.petsc_scalar_row_start() ||
        vec_end !=
            numbering.petsc_scalar_row_end() ||
        mat_local_rows != expected_local ||
        mat_local_columns != expected_local ||
        mat_global_rows !=
            numbering.petsc_scalar_row_count() ||
        mat_global_columns !=
            numbering.petsc_scalar_row_count() ||
        mat_start !=
            numbering.petsc_scalar_row_start() ||
        mat_end !=
            numbering.petsc_scalar_row_end() ||
        mat_column_start !=
            numbering.petsc_scalar_row_start() ||
        mat_column_end !=
            numbering.petsc_scalar_row_end()) {
        return PETSC_ERR_ARG_INCOMP;
    }

    return PETSC_SUCCESS;
}

} // namespace frozen_newton_linear_detail

/// Solve one frozen natural-variable Newton linearization:
///
///     J(q^k) * delta_q = -R(q^k)
///
/// The input Mat and residual Vec are treated as already assembled and
/// authoritative. This function does not evaluate or update any physical
/// state, EOS, flash, flux, accumulation, phase set or natural-variable chart.
///
/// The baseline KSP is deliberately package-independent: full-restart GMRES
/// with PCNONE and strict tolerances. It exists to verify algebraic correctness
/// of the assembled small system, not to select a production reservoir solver.
///
/// On success, *correction is owned by the caller and must be destroyed with
/// VecDestroy(). On failure, it remains null.
inline PetscErrorCode
solve_frozen_natural_variable_newton_linear_system_3d(
    MPI_Comm comm,
    const CompleteNaturalVariableAssemblySnapshot3D&
        numbering,
    Mat assembled_jacobian,
    Vec assembled_residual,
    Vec* correction,
    std::optional<
        FrozenNaturalVariableNewtonLinearCorrection3D>*
            report) {
    using namespace frozen_newton_linear_detail;

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
        correction == nullptr ||
                report == nullptr
            ? PETSC_ERR_ARG_NULL
            : PETSC_SUCCESS;
    if (local_error == PETSC_SUCCESS &&
        *correction != nullptr) {
        local_error =
            PETSC_ERR_ARG_WRONGSTATE;
    }
    if (local_error == PETSC_SUCCESS &&
        (numbering.local_rank().value() !=
             static_cast<std::uint32_t>(
                 mpi_rank) ||
         numbering.rank_count() !=
             static_cast<std::uint32_t>(
                 mpi_size))) {
        local_error =
            PETSC_ERR_ARG_INCOMP;
    }

    PetscErrorCode error =
        detail::collective_error(
            comm,
            local_error);
    if (error != PETSC_SUCCESS) {
        return error;
    }
    *correction = nullptr;
    report->reset();

    local_error =
        validate_petsc_linear_objects(
            comm,
            numbering,
            assembled_jacobian,
            assembled_residual);
    error =
        detail::collective_error(
            comm,
            local_error);
    if (error != PETSC_SUCCESS) {
        return error;
    }

    Vec rhs = nullptr;
    Vec delta = nullptr;
    Vec linear_residual = nullptr;
    KSP ksp = nullptr;

    const auto cleanup =
        [&]() {
            PetscErrorCode first =
                PETSC_SUCCESS;
            if (ksp != nullptr) {
                const PetscErrorCode destroy =
                    KSPDestroy(&ksp);
                if (first ==
                        PETSC_SUCCESS &&
                    destroy !=
                        PETSC_SUCCESS) {
                    first = destroy;
                }
            }
            if (rhs != nullptr) {
                const PetscErrorCode destroy =
                    VecDestroy(&rhs);
                if (first ==
                        PETSC_SUCCESS &&
                    destroy !=
                        PETSC_SUCCESS) {
                    first = destroy;
                }
            }
            if (linear_residual != nullptr) {
                const PetscErrorCode destroy =
                    VecDestroy(
                        &linear_residual);
                if (first ==
                        PETSC_SUCCESS &&
                    destroy !=
                        PETSC_SUCCESS) {
                    first = destroy;
                }
            }
            if (delta != nullptr) {
                const PetscErrorCode destroy =
                    VecDestroy(&delta);
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
        VecDuplicate(
            assembled_residual,
            &rhs);
    if (error == PETSC_SUCCESS) {
        error =
            VecCopy(
                assembled_residual,
                rhs);
    }
    if (error == PETSC_SUCCESS) {
        error =
            VecScale(
                rhs,
                PetscScalar{-1.0});
    }
    if (error == PETSC_SUCCESS) {
        error =
            VecDuplicate(
                assembled_residual,
                &delta);
    }
    if (error == PETSC_SUCCESS) {
        error =
            VecSet(
                delta,
                PetscScalar{0.0});
    }
    if (error == PETSC_SUCCESS) {
        error =
            KSPCreate(
                comm,
                &ksp);
    }
    if (error == PETSC_SUCCESS) {
        error =
            KSPSetOperators(
                ksp,
                assembled_jacobian,
                assembled_jacobian);
    }
    if (error == PETSC_SUCCESS) {
        error =
            KSPSetType(
                ksp,
                KSPGMRES);
    }

    PC pc = nullptr;
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
                PCNONE);
    }

    PetscInt restart = 0;
    const PetscInt global_count =
        numbering.petsc_scalar_row_count();
    if (error == PETSC_SUCCESS) {
        if (global_count <= 0) {
            error =
                PETSC_ERR_ARG_INCOMP;
        } else if (
            global_count >
                PetscInt{200}) {
            // This contract is intentionally the small-system algebraic gate,
            // not a production large-scale solver policy.
            error =
                PETSC_ERR_SUP;
        } else {
            restart =
                global_count;
            error =
                KSPGMRESSetRestart(
                    ksp,
                    restart);
        }
    }

    PetscInt max_iterations = 0;
    if (error == PETSC_SUCCESS) {
        if (global_count >
            std::numeric_limits<PetscInt>::max() /
                2) {
            error =
                PETSC_ERR_ARG_OUTOFRANGE;
        } else {
            max_iterations =
                std::max(
                    PetscInt{50},
                    std::min(
                        PetscInt{400},
                        PetscInt{2} *
                            global_count));
            error =
                KSPSetTolerances(
                    ksp,
                    PetscReal{1.0e-12},
                    PetscReal{1.0e-14},
                    PETSC_DEFAULT,
                    max_iterations);
        }
    }
    if (error == PETSC_SUCCESS) {
        error =
            KSPSetInitialGuessNonzero(
                ksp,
                PETSC_FALSE);
    }
    if (error == PETSC_SUCCESS) {
        error =
            KSPSetErrorIfNotConverged(
                ksp,
                PETSC_FALSE);
    }
    if (error != PETSC_SUCCESS) {
        (void)cleanup();
        return error;
    }

    double residual_norm = 0.0;
    double rhs_norm = 0.0;
    PetscReal residual_norm_petsc = 0.0;
    PetscReal rhs_norm_petsc = 0.0;
    error =
        VecNorm(
            assembled_residual,
            NORM_2,
            &residual_norm_petsc);
    if (error == PETSC_SUCCESS) {
        error =
            VecNorm(
                rhs,
                NORM_2,
                &rhs_norm_petsc);
    }
    if (error != PETSC_SUCCESS) {
        (void)cleanup();
        return error;
    }
    residual_norm =
        static_cast<double>(
            residual_norm_petsc);
    rhs_norm =
        static_cast<double>(
            rhs_norm_petsc);

    error =
        KSPSolve(
            ksp,
            rhs,
            delta);
    if (error != PETSC_SUCCESS) {
        (void)cleanup();
        return error;
    }

    KSPConvergedReason reason{
        KSP_CONVERGED_ITERATING};
    PetscInt iterations = -1;
    error =
        KSPGetConvergedReason(
            ksp,
            &reason);
    if (error == PETSC_SUCCESS) {
        error =
            KSPGetIterationNumber(
                ksp,
                &iterations);
    }
    if (error != PETSC_SUCCESS) {
        (void)cleanup();
        return error;
    }
    if (static_cast<int>(reason) <= 0) {
        (void)cleanup();
        return PETSC_ERR_NOT_CONVERGED;
    }

    error =
        VecDuplicate(
            assembled_residual,
            &linear_residual);
    if (error == PETSC_SUCCESS) {
        error =
            MatMult(
                assembled_jacobian,
                delta,
                linear_residual);
    }
    if (error == PETSC_SUCCESS) {
        error =
            VecAXPY(
                linear_residual,
                PetscScalar{1.0},
                assembled_residual);
    }

    PetscReal linear_norm_petsc = 0.0;
    if (error == PETSC_SUCCESS) {
        error =
            VecNorm(
                linear_residual,
                NORM_2,
                &linear_norm_petsc);
    }
    if (error != PETSC_SUCCESS) {
        (void)cleanup();
        return error;
    }

    const double linear_norm =
        static_cast<double>(
            linear_norm_petsc);
    const double accepted =
        32.0 *
        std::max(
            1.0e-14,
            1.0e-12 *
                std::max(
                    1.0,
                    rhs_norm));
    if (!std::isfinite(linear_norm) ||
        linear_norm > accepted) {
        (void)cleanup();
        return PETSC_ERR_NOT_CONVERGED;
    }

    const PetscScalar* local_values =
        nullptr;
    error =
        VecGetArrayRead(
            delta,
            &local_values);
    if (error != PETSC_SUCCESS) {
        (void)cleanup();
        return error;
    }

    std::vector<NaturalVariableLinearCorrectionEntry3D>
        entries;
    try {
        const auto residual_entries =
            numbering.residual_entries();
        const std::size_t local_count =
            residual_entries.size();
        entries.reserve(local_count);

        for (std::size_t index = 0U;
             index < local_count;
             ++index) {
            const auto& provenance =
                residual_entries[index];
            const PetscInt expected =
                numbering
                    .petsc_scalar_row_start() +
                static_cast<PetscInt>(
                    index);
            if (provenance.petsc_global_row !=
                expected) {
                throw std::invalid_argument(
                    "frozen Newton correction provenance is not in local PETSc scalar order");
            }

            entries.push_back(
                NaturalVariableLinearCorrectionEntry3D{
                    expected,
                    provenance.mesh_global_row_dof,
                    provenance.row_cell_global,
                    static_cast<std::size_t>(
                        expected %
                        static_cast<PetscInt>(
                            numbering
                                .natural_variable_count())),
                    static_cast<double>(
                        PetscRealPart(
                            local_values[index]))});
        }
    } catch (...) {
        local_error =
            PETSC_ERR_ARG_INCOMP;
    }

    const PetscErrorCode restore_error =
        VecRestoreArrayRead(
            delta,
            &local_values);
    if (local_error ==
            PETSC_SUCCESS &&
        restore_error !=
            PETSC_SUCCESS) {
        local_error =
            restore_error;
    }

    error =
        detail::collective_error(
            comm,
            local_error);
    if (error != PETSC_SUCCESS) {
        (void)cleanup();
        return error;
    }

    const char* ksp_name = nullptr;
    const char* pc_name = nullptr;
    error =
        KSPGetType(
            ksp,
            &ksp_name);
    if (error == PETSC_SUCCESS) {
        error =
            PCGetType(
                pc,
                &pc_name);
    }
    if (error != PETSC_SUCCESS ||
        ksp_name == nullptr ||
        pc_name == nullptr) {
        (void)cleanup();
        return error != PETSC_SUCCESS
            ? error
            : PETSC_ERR_PLIB;
    }

    try {
        report->emplace(
            numbering.local_rank(),
            numbering.rank_count(),
            std::string{
                numbering.natural_variable_id()},
            numbering.component_count(),
            numbering.natural_variable_count(),
            numbering.petsc_scalar_row_start(),
            numbering.petsc_scalar_row_end(),
            numbering.petsc_scalar_row_count(),
            std::string{ksp_name},
            std::string{pc_name},
            iterations,
            reason,
            residual_norm,
            rhs_norm,
            linear_norm,
            std::move(entries));
    } catch (...) {
        local_error =
            PETSC_ERR_ARG_INCOMP;
    }

    error =
        detail::collective_error(
            comm,
            local_error);
    if (error != PETSC_SUCCESS) {
        (void)cleanup();
        report->reset();
        return error;
    }

    const PetscErrorCode ksp_destroy_error =
        KSPDestroy(&ksp);
    const PetscErrorCode rhs_destroy_error =
        VecDestroy(&rhs);
    const PetscErrorCode linear_destroy_error =
        VecDestroy(
            &linear_residual);
    local_error =
        ksp_destroy_error !=
                PETSC_SUCCESS
            ? ksp_destroy_error
            : (rhs_destroy_error !=
                       PETSC_SUCCESS
                   ? rhs_destroy_error
                   : linear_destroy_error);
    error =
        detail::collective_error(
            comm,
            local_error);
    if (error != PETSC_SUCCESS) {
        (void)VecDestroy(&delta);
        report->reset();
        return error;
    }

    *correction =
        delta;
    return PETSC_SUCCESS;
}

} // namespace mpmc::flow_discretization_petsc

#endif // MPMC_FLOW_DISCRETIZATION_PETSC_NATURAL_VARIABLE_NEWTON_LINEAR_SYSTEM_HPP
