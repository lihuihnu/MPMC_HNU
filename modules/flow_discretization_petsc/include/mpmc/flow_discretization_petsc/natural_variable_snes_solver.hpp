#ifndef MPMC_FLOW_DISCRETIZATION_PETSC_NATURAL_VARIABLE_SNES_SOLVER_HPP
#define MPMC_FLOW_DISCRETIZATION_PETSC_NATURAL_VARIABLE_SNES_SOLVER_HPP

#include <mpmc/flow_discretization_petsc/complete_natural_variable_assembly_snapshot.hpp>

#include <petscsnes.h>

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

namespace mpmc::flow_discretization_petsc {

inline constexpr std::string_view
    natural_variable_snes_solver_convention =
        "flow_discretization_petsc/natural-variable-snes-solver/v1";

enum class NaturalVariableSnesEvaluationStatus3D :
    std::uint8_t {
    success,
    domain_error
};

using NaturalVariableSnesFunctionEvaluator3D =
    PetscErrorCode (*)(
        Vec state,
        Vec residual,
        void* user_context,
        NaturalVariableSnesEvaluationStatus3D* status);

using NaturalVariableSnesJacobianEvaluator3D =
    PetscErrorCode (*)(
        Vec state,
        Mat jacobian,
        void* user_context,
        NaturalVariableSnesEvaluationStatus3D* status);

/// PETSc SNES line-search precheck hook.
///
/// PETSc supplies the current state X and Newton/search direction Y and forms
/// trial states as X - lambda*Y. The hook may modify Y and set changed_direction
/// to PETSC_TRUE. This is the correct PETSc 3.19 path for recoverable
/// trial-step domain protection; SNESSetFunctionDomainError() is reserved for
/// a genuinely non-evaluable state and terminates SNESSolve().
using NaturalVariableSnesStepPrecheck3D =
    PetscErrorCode (*)(
        Vec state,
        Vec search_direction,
        void* user_context,
        PetscBool* changed_direction);

struct NaturalVariableSnesEvaluator3D {
    NaturalVariableSnesFunctionEvaluator3D
        function{};
    NaturalVariableSnesJacobianEvaluator3D
        jacobian{};
    NaturalVariableSnesStepPrecheck3D
        step_precheck{};
    void* user_context{};
};

struct NaturalVariableSnesFailureDiagnostics3D {
    SNESConvergedReason snes_reason{
        SNES_CONVERGED_ITERATING};
    KSPConvergedReason ksp_reason{
        KSP_CONVERGED_ITERATING};
    int pc_failed_reason{};
    KSPConvergedReason asm_sub_ksp_reason{
        KSP_CONVERGED_ITERATING};
    int asm_sub_pc_failed_reason{};
    PetscInt nonlinear_iterations{};
    PetscInt function_evaluations{};
    PetscInt jacobian_evaluations{};
    PetscInt function_domain_errors{};
    PetscInt jacobian_domain_errors{};
    PetscInt line_search_prechecks{};
    PetscInt line_search_direction_changes{};
    double function_l2_norm{};
};

struct NaturalVariableSnesSolutionEntry3D {
    PetscInt petsc_global_scalar{-1};
    mpmc::mesh::GlobalDofIndex mesh_global_dof{
        mpmc::mesh::GlobalDofIndex::value_type{0}};
    mpmc::mesh::GlobalEntityId cell_global{
        mpmc::mesh::GlobalEntityId::value_type{0}};
    std::size_t natural_variable_slot{};
    double value{};
};

class NaturalVariableSnesSolveReport3D {
public:
    static constexpr std::string_view convention =
        natural_variable_snes_solver_convention;

    NaturalVariableSnesSolveReport3D(
        mpmc::mesh::PartitionRank local_rank,
        std::uint32_t rank_count,
        std::string natural_variable_id,
        std::size_t component_count,
        std::size_t natural_variable_count,
        PetscInt petsc_scalar_row_start,
        PetscInt petsc_scalar_row_end,
        PetscInt petsc_scalar_row_count,
        std::string snes_type,
        std::string line_search_type,
        std::string ksp_type,
        std::string pc_type,
        PetscInt default_asm_overlap,
        SNESConvergedReason converged_reason,
        PetscInt nonlinear_iterations,
        PetscInt function_evaluations,
        PetscInt jacobian_evaluations,
        PetscInt function_domain_errors,
        PetscInt jacobian_domain_errors,
        PetscInt line_search_prechecks,
        PetscInt line_search_direction_changes,
        double final_function_l2_norm,
        std::vector<NaturalVariableSnesSolutionEntry3D>
            locally_owned_solution)
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
          snes_type_(std::move(snes_type)),
          line_search_type_(
              std::move(line_search_type)),
          ksp_type_(std::move(ksp_type)),
          pc_type_(std::move(pc_type)),
          default_asm_overlap_(
              default_asm_overlap),
          converged_reason_(converged_reason),
          nonlinear_iterations_(
              nonlinear_iterations),
          function_evaluations_(
              function_evaluations),
          jacobian_evaluations_(
              jacobian_evaluations),
          function_domain_errors_(
              function_domain_errors),
          jacobian_domain_errors_(
              jacobian_domain_errors),
          line_search_prechecks_(
              line_search_prechecks),
          line_search_direction_changes_(
              line_search_direction_changes),
          final_function_l2_norm_(
              final_function_l2_norm),
          locally_owned_solution_(
              std::move(
                  locally_owned_solution)) {
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
    snes_type() const noexcept {
        return snes_type_;
    }

    [[nodiscard]] std::string_view
    line_search_type() const noexcept {
        return line_search_type_;
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
    default_asm_overlap() const noexcept {
        return default_asm_overlap_;
    }

    [[nodiscard]] SNESConvergedReason
    converged_reason() const noexcept {
        return converged_reason_;
    }

    [[nodiscard]] PetscInt
    nonlinear_iterations() const noexcept {
        return nonlinear_iterations_;
    }

    [[nodiscard]] PetscInt
    function_evaluations() const noexcept {
        return function_evaluations_;
    }

    [[nodiscard]] PetscInt
    jacobian_evaluations() const noexcept {
        return jacobian_evaluations_;
    }

    [[nodiscard]] PetscInt
    function_domain_errors() const noexcept {
        return function_domain_errors_;
    }

    [[nodiscard]] PetscInt
    jacobian_domain_errors() const noexcept {
        return jacobian_domain_errors_;
    }

    [[nodiscard]] PetscInt
    line_search_prechecks() const noexcept {
        return line_search_prechecks_;
    }

    [[nodiscard]] PetscInt
    line_search_direction_changes() const noexcept {
        return line_search_direction_changes_;
    }

    /// PETSc function norm used by the nonlinear solve.
    ///
    /// Without explicit row scaling this is the native mixed-equation norm and
    /// component/energy/fugacity rows retain different physical units. When an
    /// explicit frozen left row-scaling vector is supplied to the solver, this
    /// is the norm of D*R. Physical acceptance must still be checked against an
    /// independently reassembled unscaled residual.
    [[nodiscard]] double
    final_function_l2_norm() const noexcept {
        return final_function_l2_norm_;
    }

    [[nodiscard]] std::span<
        const NaturalVariableSnesSolutionEntry3D>
    locally_owned_solution() const noexcept {
        return locally_owned_solution_;
    }

private:
    void validate() const {
        if (rank_count_ == 0U ||
            local_rank_.value() >= rank_count_ ||
            natural_variable_id_.empty() ||
            component_count_ < 2U ||
            natural_variable_count_ <= 1U ||
            (natural_variable_count_ - 1U) %
                    component_count_ !=
                0U ||
            (natural_variable_count_ - 1U) /
                    component_count_ ==
                0U ||
            (natural_variable_count_ - 1U) /
                    component_count_ >
                mpmc::flow::fixed_three_phase_count ||
            petsc_scalar_row_start_ < 0 ||
            petsc_scalar_row_end_ <
                petsc_scalar_row_start_ ||
            petsc_scalar_row_count_ <
                petsc_scalar_row_end_ ||
            snes_type_.empty() ||
            line_search_type_.empty() ||
            ksp_type_.empty() ||
            pc_type_.empty() ||
            default_asm_overlap_ < 0 ||
            static_cast<int>(
                converged_reason_) <= 0 ||
            nonlinear_iterations_ < 0 ||
            function_evaluations_ <= 0 ||
            jacobian_evaluations_ <= 0 ||
            function_domain_errors_ < 0 ||
            jacobian_domain_errors_ < 0 ||
            line_search_prechecks_ < 0 ||
            line_search_direction_changes_ < 0 ||
            line_search_direction_changes_ >
                line_search_prechecks_ ||
            !std::isfinite(
                final_function_l2_norm_) ||
            final_function_l2_norm_ < 0.0 ||
            locally_owned_solution_.size() !=
                static_cast<std::size_t>(
                    petsc_scalar_row_end_ -
                    petsc_scalar_row_start_)) {
            throw std::invalid_argument(
                "mpmc::flow_discretization_petsc: malformed natural-variable SNES solve report");
        }

        for (std::size_t index = 0U;
             index < locally_owned_solution_.size();
             ++index) {
            const auto& entry =
                locally_owned_solution_[index];
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
                    "mpmc::flow_discretization_petsc: SNES solution numbering/value mismatch");
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
    std::string snes_type_;
    std::string line_search_type_;
    std::string ksp_type_;
    std::string pc_type_;
    PetscInt default_asm_overlap_{};
    SNESConvergedReason converged_reason_{
        SNES_CONVERGED_ITERATING};
    PetscInt nonlinear_iterations_{};
    PetscInt function_evaluations_{};
    PetscInt jacobian_evaluations_{};
    PetscInt function_domain_errors_{};
    PetscInt jacobian_domain_errors_{};
    PetscInt line_search_prechecks_{};
    PetscInt line_search_direction_changes_{};
    double final_function_l2_norm_{};
    std::vector<NaturalVariableSnesSolutionEntry3D>
        locally_owned_solution_;
};

namespace natural_variable_snes_detail {

struct CallbackContext {
    NaturalVariableSnesEvaluator3D evaluator;
    Vec row_scaling{};
    PetscInt function_evaluations{};
    PetscInt jacobian_evaluations{};
    PetscInt function_domain_errors{};
    PetscInt jacobian_domain_errors{};
    PetscInt line_search_prechecks{};
    PetscInt line_search_direction_changes{};
};

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
validate_linear_layout(
    MPI_Comm comm,
    const CompleteNaturalVariableAssemblySnapshot3D&
        numbering,
    Vec state,
    Mat jacobian_template) {
    if (state == nullptr ||
        jacobian_template == nullptr) {
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
                    jacobian_template)))) {
        return PETSC_ERR_ARG_INCOMP;
    }

    PetscBool assembled =
        PETSC_FALSE;
    PetscErrorCode error =
        MatAssembled(
            jacobian_template,
            &assembled);
    if (error != PETSC_SUCCESS) {
        return error;
    }
    if (assembled != PETSC_TRUE) {
        return PETSC_ERR_ARG_WRONGSTATE;
    }

    PetscBool is_mpiaij =
        PETSC_FALSE;
    error =
        PetscObjectTypeCompare(
            reinterpret_cast<PetscObject>(
                jacobian_template),
            MATMPIAIJ,
            &is_mpiaij);
    if (error != PETSC_SUCCESS) {
        return error;
    }
    if (is_mpiaij != PETSC_TRUE) {
        return PETSC_ERR_ARG_INCOMP;
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
                jacobian_template,
                &mat_local_rows,
                &mat_local_columns);
    }
    if (error == PETSC_SUCCESS) {
        error =
            MatGetSize(
                jacobian_template,
                &mat_global_rows,
                &mat_global_columns);
    }
    if (error == PETSC_SUCCESS) {
        error =
            MatGetOwnershipRange(
                jacobian_template,
                &mat_start,
                &mat_end);
    }
    if (error == PETSC_SUCCESS) {
        error =
            MatGetOwnershipRangeColumn(
                jacobian_template,
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

[[nodiscard]] inline PetscErrorCode
validate_row_scaling(
    MPI_Comm comm,
    const CompleteNaturalVariableAssemblySnapshot3D& numbering,
    Vec row_scaling) {
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
        VecGetLocalSize(row_scaling, &local);
    if (error == PETSC_SUCCESS) {
        error =
            VecGetSize(row_scaling, &global);
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

    const PetscInt expected_local =
        numbering.petsc_scalar_row_end() -
        numbering.petsc_scalar_row_start();
    if (local != expected_local ||
        global != numbering.petsc_scalar_row_count() ||
        start != numbering.petsc_scalar_row_start() ||
        end != numbering.petsc_scalar_row_end()) {
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
                PetscRealPart(values[index]));
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

inline PetscErrorCode
apply_row_scaling_to_residual(
    Vec residual,
    Vec row_scaling) {
    if (row_scaling == nullptr) {
        return PETSC_SUCCESS;
    }
    PetscInt residual_local = -1;
    PetscInt scaling_local = -1;
    PetscErrorCode error =
        VecGetLocalSize(
            residual,
            &residual_local);
    if (error == PETSC_SUCCESS) {
        error =
            VecGetLocalSize(
                row_scaling,
                &scaling_local);
    }
    if (error != PETSC_SUCCESS) {
        return error;
    }
    if (residual_local != scaling_local) {
        return PETSC_ERR_ARG_SIZ;
    }

    PetscScalar* residual_values = nullptr;
    const PetscScalar* scaling_values = nullptr;
    error =
        VecGetArray(
            residual,
            &residual_values);
    if (error == PETSC_SUCCESS) {
        error =
            VecGetArrayRead(
                row_scaling,
                &scaling_values);
    }
    if (error != PETSC_SUCCESS) {
        if (residual_values != nullptr) {
            (void)VecRestoreArray(
                residual,
                &residual_values);
        }
        return error;
    }

    for (PetscInt index = 0;
         index < residual_local;
         ++index) {
        residual_values[index] *=
            scaling_values[index];
    }

    const PetscErrorCode restore_scaling =
        VecRestoreArrayRead(
            row_scaling,
            &scaling_values);
    const PetscErrorCode restore_residual =
        VecRestoreArray(
            residual,
            &residual_values);
    return restore_scaling != PETSC_SUCCESS
        ? restore_scaling
        : restore_residual;
}

inline PetscErrorCode
form_function(
    SNES snes,
    Vec state,
    Vec residual,
    void* raw_context) {
    if (raw_context == nullptr) {
        return PETSC_ERR_ARG_NULL;
    }
    auto* context =
        static_cast<CallbackContext*>(
            raw_context);
    ++context->function_evaluations;

    PetscErrorCode error =
        VecSet(
            residual,
            PetscScalar{0.0});
    if (error != PETSC_SUCCESS) {
        return error;
    }

    NaturalVariableSnesEvaluationStatus3D status =
        NaturalVariableSnesEvaluationStatus3D::
            success;
    error =
        context->evaluator.function(
            state,
            residual,
            context->evaluator.user_context,
            &status);
    if (error != PETSC_SUCCESS) {
        return error;
    }

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
            PetscObjectComm(
                reinterpret_cast<PetscObject>(
                    state))) !=
        MPI_SUCCESS) {
        return PETSC_ERR_MPI;
    }

    if (global_status == 2) {
        return PETSC_ERR_ARG_INCOMP;
    }
    if (global_status == 1) {
        ++context->function_domain_errors;
        return SNESSetFunctionDomainError(
            snes);
    }

    error =
        VecAssemblyBegin(
            residual);
    if (error == PETSC_SUCCESS) {
        error =
            VecAssemblyEnd(
                residual);
    }
    if (error == PETSC_SUCCESS) {
        error =
            apply_row_scaling_to_residual(
                residual,
                context->row_scaling);
    }
    return error;
}

inline PetscErrorCode
line_search_precheck(
    SNESLineSearch,
    Vec state,
    Vec search_direction,
    PetscBool* changed_direction,
    void* raw_context) {
    if (raw_context == nullptr ||
        changed_direction == nullptr) {
        return PETSC_ERR_ARG_NULL;
    }

    auto* context =
        static_cast<CallbackContext*>(
            raw_context);
    ++context->line_search_prechecks;
    *changed_direction =
        PETSC_FALSE;

    if (context->evaluator.step_precheck ==
        nullptr) {
        return PETSC_SUCCESS;
    }

    PetscBool changed =
        PETSC_FALSE;
    const PetscErrorCode error =
        context->evaluator.step_precheck(
            state,
            search_direction,
            context->evaluator.user_context,
            &changed);
    if (error != PETSC_SUCCESS) {
        return error;
    }
    if (changed != PETSC_FALSE &&
        changed != PETSC_TRUE) {
        return PETSC_ERR_ARG_INCOMP;
    }
    if (changed == PETSC_TRUE) {
        ++context
              ->line_search_direction_changes;
    }
    *changed_direction =
        changed;
    return PETSC_SUCCESS;
}

inline PetscErrorCode
form_jacobian(
    SNES snes,
    Vec state,
    Mat jacobian,
    Mat preconditioner,
    void* raw_context) {
    if (raw_context == nullptr ||
        jacobian == nullptr ||
        preconditioner == nullptr) {
        return PETSC_ERR_ARG_NULL;
    }
    if (jacobian != preconditioner) {
        return PETSC_ERR_SUP;
    }

    auto* context =
        static_cast<CallbackContext*>(
            raw_context);
    ++context->jacobian_evaluations;

    PetscErrorCode error =
        MatZeroEntries(
            jacobian);
    if (error != PETSC_SUCCESS) {
        return error;
    }

    NaturalVariableSnesEvaluationStatus3D status =
        NaturalVariableSnesEvaluationStatus3D::
            success;
    error =
        context->evaluator.jacobian(
            state,
            jacobian,
            context->evaluator.user_context,
            &status);
    if (error != PETSC_SUCCESS) {
        return error;
    }

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
            PetscObjectComm(
                reinterpret_cast<PetscObject>(
                    state))) !=
        MPI_SUCCESS) {
        return PETSC_ERR_MPI;
    }

    if (global_status == 2) {
        return PETSC_ERR_ARG_INCOMP;
    }
    if (global_status == 1) {
        ++context->jacobian_domain_errors;
        return SNESSetJacobianDomainError(
            snes);
    }

    error =
        MatAssemblyBegin(
            jacobian,
            MAT_FINAL_ASSEMBLY);
    if (error == PETSC_SUCCESS) {
        error =
            MatAssemblyEnd(
                jacobian,
                MAT_FINAL_ASSEMBLY);
    }
    if (error == PETSC_SUCCESS &&
        context->row_scaling != nullptr) {
        error =
            MatDiagonalScale(
                jacobian,
                context->row_scaling,
                nullptr);
    }
    return error;
}

} // namespace natural_variable_snes_detail

/// Construct a frozen positive left row-equilibration vector from one
/// independently assembled analytic Jacobian.
///
/// For each owned equation row i:
///   D_i = 1 / max_j |J_ij(q0)|.
///
/// This is algebraic solver scaling only. It does not mutate the physical
/// residual/Jacobian snapshot and it is frozen for the subsequent nonlinear
/// solve. Rows without a finite nonzero analytic coefficient are rejected.
inline PetscErrorCode
make_natural_variable_initial_row_equilibration_3d(
    MPI_Comm comm,
    const CompleteNaturalVariableAssemblySnapshot3D& snapshot,
    Vec* row_scaling) {
    if (row_scaling == nullptr) {
        return PETSC_ERR_ARG_NULL;
    }
    if (*row_scaling != nullptr) {
        return PETSC_ERR_ARG_WRONGSTATE;
    }

    const PetscInt start =
        snapshot.petsc_scalar_row_start();
    const PetscInt end =
        snapshot.petsc_scalar_row_end();
    const PetscInt local_count =
        end - start;
    if (local_count <= 0 ||
        snapshot.petsc_scalar_row_count() <= 0) {
        return PETSC_ERR_ARG_SIZ;
    }

    std::vector<double> row_max(
        static_cast<std::size_t>(
            local_count),
        0.0);
    for (const auto& entry :
         snapshot.jacobian_entries()) {
        if (entry.petsc_global_row < start ||
            entry.petsc_global_row >= end) {
            return PETSC_ERR_ARG_INCOMP;
        }
        const std::size_t local =
            static_cast<std::size_t>(
                entry.petsc_global_row - start);
        const double magnitude =
            std::abs(entry.value);
        if (!std::isfinite(magnitude)) {
            return PETSC_ERR_ARG_OUTOFRANGE;
        }
        row_max[local] =
            std::max(
                row_max[local],
                magnitude);
    }
    for (const double value : row_max) {
        if (!(value > 0.0) ||
            !std::isfinite(value)) {
            return PETSC_ERR_ARG_WRONGSTATE;
        }
    }

    Vec scaling = nullptr;
    PetscErrorCode error =
        VecCreateMPI(
            comm,
            local_count,
            snapshot.petsc_scalar_row_count(),
            &scaling);
    if (error != PETSC_SUCCESS) {
        return error;
    }

    PetscInt actual_start = -1;
    PetscInt actual_end = -1;
    error =
        VecGetOwnershipRange(
            scaling,
            &actual_start,
            &actual_end);
    if (error != PETSC_SUCCESS ||
        actual_start != start ||
        actual_end != end) {
        (void)VecDestroy(&scaling);
        return error != PETSC_SUCCESS
            ? error
            : PETSC_ERR_ARG_INCOMP;
    }

    PetscScalar* values = nullptr;
    error =
        VecGetArray(
            scaling,
            &values);
    if (error != PETSC_SUCCESS) {
        (void)VecDestroy(&scaling);
        return error;
    }
    for (PetscInt index = 0;
         index < local_count;
         ++index) {
        values[index] =
            PetscScalar{
                1.0 /
                row_max[
                    static_cast<std::size_t>(
                        index)]};
    }
    const PetscErrorCode restore =
        VecRestoreArray(
            scaling,
            &values);
    if (restore != PETSC_SUCCESS) {
        (void)VecDestroy(&scaling);
        return restore;
    }

    *row_scaling = scaling;
    return PETSC_SUCCESS;
}

/// Solve the nonlinear natural-variable system with PETSc-owned orchestration.
///
/// Project defaults:
///   SNESNEWTONLS + backtracking line search
///   KSPGMRES + PCASM(overlap=1, restricted)
///   ASM local solve: KSPPREONLY + PCLU, with exact-zero diagonal reordering
///
/// An optional frozen positive left row-scaling vector D may be supplied.
/// PETSc then solves D*R=0 with D*J while the caller-owned physical evaluator
/// continues to publish the unscaled R and J. No scaling is inferred by default.
///
/// This first production contract deliberately does not call
/// SNESSetFromOptions(): unrestricted PETSc options could replace the audited
/// analytic Jacobian path with matrix-free/finite-difference behavior. Solver
/// configurability must later enter through an explicit validated settings
/// contract.
///
/// The supplied evaluators are the only path from PETSc state x to F(x)
/// and J(x). An optional PETSc line-search precheck may modify the search
/// direction before BT evaluates a trial state. Hard function/Jacobian domain
/// errors use PETSc domain-error semantics and terminate the solve in 3.19.
/// This adapter performs no EOS/flash/flux/accumulation evaluation itself and
/// contains no custom Newton or line-search loop.
///
/// On success, *solution is caller-owned and must be destroyed with
/// VecDestroy(). On failure, *solution remains null and report is empty.
inline PetscErrorCode
solve_natural_variable_snes_3d(
    MPI_Comm comm,
    const CompleteNaturalVariableAssemblySnapshot3D&
        numbering,
    Vec initial_state,
    Mat jacobian_structure_template,
    NaturalVariableSnesEvaluator3D evaluator,
    Vec* solution,
    std::optional<NaturalVariableSnesSolveReport3D>*
        report,
    Vec row_scaling = nullptr,
    std::optional<
        NaturalVariableSnesFailureDiagnostics3D>*
            failure_diagnostics = nullptr) {
    using namespace natural_variable_snes_detail;

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
        detail::collective_error(
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
        validate_linear_layout(
            comm,
            numbering,
            initial_state,
            jacobian_structure_template);
    if (local_error == PETSC_SUCCESS) {
        local_error =
            validate_row_scaling(
                comm,
                numbering,
                row_scaling);
    }
    error =
        detail::collective_error(
            comm,
            local_error);
    if (error != PETSC_SUCCESS) {
        return error;
    }

    Vec solved_state = nullptr;
    Vec residual = nullptr;
    Mat jacobian = nullptr;
    SNES snes = nullptr;

    const auto cleanup =
        [&]() {
            PetscErrorCode first =
                PETSC_SUCCESS;
            if (snes != nullptr) {
                const PetscErrorCode destroy =
                    SNESDestroy(&snes);
                if (first == PETSC_SUCCESS &&
                    destroy != PETSC_SUCCESS) {
                    first = destroy;
                }
            }
            if (residual != nullptr) {
                const PetscErrorCode destroy =
                    VecDestroy(&residual);
                if (first == PETSC_SUCCESS &&
                    destroy != PETSC_SUCCESS) {
                    first = destroy;
                }
            }
            if (jacobian != nullptr) {
                const PetscErrorCode destroy =
                    MatDestroy(&jacobian);
                if (first == PETSC_SUCCESS &&
                    destroy != PETSC_SUCCESS) {
                    first = destroy;
                }
            }
            if (solved_state != nullptr) {
                const PetscErrorCode destroy =
                    VecDestroy(&solved_state);
                if (first == PETSC_SUCCESS &&
                    destroy != PETSC_SUCCESS) {
                    first = destroy;
                }
            }
            return first;
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
                MAT_DO_NOT_COPY_VALUES,
                &jacobian);
    }
    if (error == PETSC_SUCCESS) {
        error =
            MatSetOption(
                jacobian,
                MAT_NEW_NONZERO_LOCATION_ERR,
                PETSC_TRUE);
    }
    if (error == PETSC_SUCCESS) {
        error =
            MatSetOption(
                jacobian,
                MAT_NEW_NONZERO_ALLOCATION_ERR,
                PETSC_TRUE);
    }
    if (error == PETSC_SUCCESS) {
        error =
            SNESCreate(
                comm,
                &snes);
    }
    if (error != PETSC_SUCCESS) {
        (void)cleanup();
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

    // Natural-variable equation ordering does not guarantee a nonzero scalar
    // diagonal even when J is nonsingular. PETSc ASM defaults to PREONLY+ILU(0)
    // on each local block, so a valid coupled Jacobian can fail during ILU
    // setup before GMRES performs an iteration. Keep the audited top-level
    // GMRES+ASM contract, but use an exact PETSc LU solve inside each ASM block
    // for this correctness baseline. Before local factorization, PETSc reorders
    // exact-zero diagonal entries onto usable pivots; no Jacobian value is
    // shifted or regularized. The private prefix prevents these internal
    // sub-solver options from replacing the top-level KSP/PC or analytic J.
    constexpr const char* asm_options_prefix =
        "mpmc_natural_variable_";
    constexpr const char* asm_sub_pc_option =
        "-mpmc_natural_variable_sub_pc_type";
    constexpr const char* asm_sub_pc_reorder_option =
        "-mpmc_natural_variable_sub_pc_factor_nonzeros_along_diagonal";
    bool asm_sub_pc_option_installed = false;
    bool asm_sub_pc_reorder_option_installed = false;
    if (error == PETSC_SUCCESS) {
        error =
            PCSetOptionsPrefix(
                pc,
                asm_options_prefix);
    }
    if (error == PETSC_SUCCESS) {
        error =
            PetscOptionsSetValue(
                nullptr,
                asm_sub_pc_option,
                "lu");
        asm_sub_pc_option_installed =
            error == PETSC_SUCCESS;
    }
    if (error == PETSC_SUCCESS) {
        // With explicit frozen row equilibration, real-SI natural variables can
        // leave a numerically tiny but nonzero diagonal after row scaling
        // (notably pressure derivatives against composition/energy rows).
        // PETSc's sparse LU does not perform numerical pivoting; allow its
        // documented nonzero-diagonal reordering to treat values below the
        // normalized 1e-10 threshold as weak pivots. Preserve the historical
        // exact-zero-only behavior for unscaled callers.
        const char* reorder_tolerance =
            row_scaling != nullptr
                ? "1.0e-10"
                : "0.0";
        error =
            PetscOptionsSetValue(
                nullptr,
                asm_sub_pc_reorder_option,
                reorder_tolerance);
        asm_sub_pc_reorder_option_installed =
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
    if (error != PETSC_SUCCESS) {
        if (asm_sub_pc_reorder_option_installed) {
            (void)PetscOptionsClearValue(
                nullptr,
                asm_sub_pc_reorder_option);
        }
        if (asm_sub_pc_option_installed) {
            (void)PetscOptionsClearValue(
                nullptr,
                asm_sub_pc_option);
        }
        (void)cleanup();
        return error;
    }

    error =
        SNESSolve(
            snes,
            nullptr,
            solved_state);

    PetscErrorCode options_clear_error =
        PETSC_SUCCESS;
    if (asm_sub_pc_reorder_option_installed) {
        options_clear_error =
            PetscOptionsClearValue(
                nullptr,
                asm_sub_pc_reorder_option);
        asm_sub_pc_reorder_option_installed =
            false;
    }
    if (asm_sub_pc_option_installed) {
        const PetscErrorCode clear_pc_type_error =
            PetscOptionsClearValue(
                nullptr,
                asm_sub_pc_option);
        asm_sub_pc_option_installed =
            false;
        if (options_clear_error == PETSC_SUCCESS &&
            clear_pc_type_error != PETSC_SUCCESS) {
            options_clear_error =
                clear_pc_type_error;
        }
    }
    if (error == PETSC_SUCCESS &&
        options_clear_error !=
            PETSC_SUCCESS) {
        error =
            options_clear_error;
    }
    if (error != PETSC_SUCCESS) {
        (void)cleanup();
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
        (void)cleanup();
        return error;
    }
    if (static_cast<int>(reason) <= 0) {
        if (failure_diagnostics != nullptr) {
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
        (void)cleanup();
        return PETSC_ERR_NOT_CONVERGED;
    }

    error =
        SNESComputeFunction(
            snes,
            solved_state,
            residual);
    PetscReal final_function_norm = 0.0;
    if (error == PETSC_SUCCESS) {
        error =
            VecNorm(
                residual,
                NORM_2,
                &final_function_norm);
    }
    if (error != PETSC_SUCCESS) {
        (void)cleanup();
        return error;
    }

    const char* snes_name = nullptr;
    const char* line_search_name = nullptr;
    const char* ksp_name = nullptr;
    const char* pc_name = nullptr;
    error =
        SNESGetType(
            snes,
            &snes_name);
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
        (void)cleanup();
        return error != PETSC_SUCCESS
            ? error
            : PETSC_ERR_PLIB;
    }

    const PetscScalar* local_values = nullptr;
    error =
        VecGetArrayRead(
            solved_state,
            &local_values);
    if (error != PETSC_SUCCESS) {
        (void)cleanup();
        return error;
    }

    std::vector<NaturalVariableSnesSolutionEntry3D>
        solution_entries;
    try {
        const auto provenance =
            numbering.residual_entries();
        solution_entries.reserve(
            provenance.size());
        for (std::size_t index = 0U;
             index < provenance.size();
             ++index) {
            const auto& row =
                provenance[index];
            const PetscInt expected =
                numbering
                    .petsc_scalar_row_start() +
                static_cast<PetscInt>(
                    index);
            if (row.petsc_global_row !=
                expected) {
                throw std::invalid_argument(
                    "SNES solution provenance is not in local PETSc scalar order");
            }
            solution_entries.push_back(
                NaturalVariableSnesSolutionEntry3D{
                    expected,
                    row.mesh_global_row_dof,
                    row.row_cell_global,
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
            solved_state,
            &local_values);
    if (local_error == PETSC_SUCCESS &&
        restore_error != PETSC_SUCCESS) {
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
            std::string{snes_name},
            std::string{line_search_name},
            std::string{ksp_name},
            std::string{pc_name},
            asm_overlap,
            reason,
            nonlinear_iterations,
            callback_context.function_evaluations,
            callback_context.jacobian_evaluations,
            callback_context.function_domain_errors,
            callback_context.jacobian_domain_errors,
            callback_context.line_search_prechecks,
            callback_context.line_search_direction_changes,
            static_cast<double>(
                final_function_norm),
            std::move(solution_entries));
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

    const PetscErrorCode snes_destroy_error =
        SNESDestroy(&snes);
    const PetscErrorCode residual_destroy_error =
        VecDestroy(&residual);
    const PetscErrorCode jacobian_destroy_error =
        MatDestroy(&jacobian);

    local_error =
        snes_destroy_error != PETSC_SUCCESS
            ? snes_destroy_error
            : (residual_destroy_error !=
                       PETSC_SUCCESS
                   ? residual_destroy_error
                   : jacobian_destroy_error);
    error =
        detail::collective_error(
            comm,
            local_error);
    if (error != PETSC_SUCCESS) {
        (void)VecDestroy(&solved_state);
        report->reset();
        return error;
    }

    *solution =
        solved_state;
    return PETSC_SUCCESS;
}

} // namespace mpmc::flow_discretization_petsc

#endif // MPMC_FLOW_DISCRETIZATION_PETSC_NATURAL_VARIABLE_SNES_SOLVER_HPP
