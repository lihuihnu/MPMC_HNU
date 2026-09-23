#ifndef MPMC_WELL_DISCRETIZATION_PETSC_FIXED_TOTAL_MOLAR_RATE_TIMESTEP_DRIVER_HPP
#define MPMC_WELL_DISCRETIZATION_PETSC_FIXED_TOTAL_MOLAR_RATE_TIMESTEP_DRIVER_HPP

#include <mpmc/flow_discretization_petsc/accepted_physical_time.hpp>
#include <mpmc/well_discretization_petsc/fixed_total_molar_rate_control.hpp>

#include <petscvec.h>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <memory>
#include <optional>
#include <string_view>
#include <utility>

namespace mpmc::well_discretization_petsc {

inline constexpr std::string_view
    fixed_total_molar_rate_timestep_driver_convention =
        "well-discretization-petsc/single-well/fixed-total-molar-rate-physical-timestep-driver/v1";

struct FixedTotalMolarRatePhysicalTimestepDriverOptions3D {
    mpmc::flow_discretization_petsc::
        AdaptiveTimestepControllerOptions3D
            adaptive;
};

struct FixedTotalMolarRatePhysicalTimestepDriverReport3D {
    mpmc::flow_discretization_petsc::
        AdaptiveTimestepControllerReport3D
            adaptive;
    std::optional<
        mpmc::flow_discretization_petsc::
            AcceptedPhysicalTimestepRecord3D>
        accepted_record;
    std::optional<
        FixedTotalMolarRateWellControlSolveReport3D>
        accepted_solve;
    double entry_bottom_hole_pressure_pa{};
    double accepted_bottom_hole_pressure_pa{};

    [[nodiscard]] bool accepted() const noexcept {
        return adaptive.accepted() &&
            accepted_record.has_value() &&
            accepted_solve.has_value();
    }
};

namespace fixed_total_molar_rate_timestep_driver_detail {

[[nodiscard]] inline bool same_timestep(
    double first,
    double second) noexcept {
    if (!std::isfinite(first) ||
        !std::isfinite(second)) {
        return false;
    }
    const double scale =
        std::max(
            {1.0,
             std::abs(first),
             std::abs(second)});
    return std::abs(first - second) <=
        64.0 *
            std::numeric_limits<double>::epsilon() *
            scale;
}

struct DriverContext3D {
    MPI_Comm comm{};
    mpmc::flow_discretization_petsc::
        PhaseTransitionRebuiltNaturalVariableSystem3D*
            reservoir_system{};
    FixedTotalMolarRateWellSourceEvaluatorContext3D*
        source_context{};
    double target_total_molar_rate_mol_per_s{};
    double entry_bottom_hole_pressure_pa{};
    const mpmc::flow_discretization_petsc::
        AdaptiveTimestepControllerOptions3D*
            adaptive_options{};
    mpmc::flow_discretization_petsc::
        AcceptedPhysicalTimeClock3D*
            clock{};
    double* accepted_bottom_hole_pressure_pa{};

    std::unique_ptr<
        FixedTotalMolarRateWellControlSystem3D>
        pending_control_system;
    Vec pending_augmented_state{};
    std::optional<
        FixedTotalMolarRateWellControlSolveReport3D>
        pending_solve_report;

    std::optional<
        mpmc::flow_discretization_petsc::
            AcceptedPhysicalTimestepRecord3D>
        accepted_record;
    std::optional<
        FixedTotalMolarRateWellControlSolveReport3D>
        accepted_solve;
};

[[nodiscard]] inline PetscErrorCode
clear_pending(
    DriverContext3D* context) noexcept {
    if (context == nullptr) {
        return PETSC_ERR_ARG_NULL;
    }
    PetscErrorCode error =
        PETSC_SUCCESS;
    if (context->pending_augmented_state !=
        nullptr) {
        error =
            VecDestroy(
                &context
                     ->pending_augmented_state);
    }
    context->pending_control_system.reset();
    context->pending_solve_report.reset();
    return error;
}

[[nodiscard]] inline PetscErrorCode
restore_entry(
    DriverContext3D* context,
    double entry_timestep_seconds,
    PetscErrorCode primary_error) noexcept {
    if (context == nullptr ||
        context->reservoir_system == nullptr ||
        context->source_context == nullptr) {
        return primary_error != PETSC_SUCCESS
            ? primary_error
            : PETSC_ERR_ARG_NULL;
    }

    const PetscErrorCode dt_error =
        context->reservoir_system
            ->set_trial_timestep_seconds(
                entry_timestep_seconds);

    PetscErrorCode bhp_error =
        PETSC_SUCCESS;
    try {
        context->source_context
            ->begin_evaluation(
                context
                    ->entry_bottom_hole_pressure_pa);
    } catch (...) {
        bhp_error =
            PETSC_ERR_ARG_INCOMP;
    }

    if (primary_error != PETSC_SUCCESS) {
        return primary_error;
    }
    if (dt_error != PETSC_SUCCESS) {
        return dt_error;
    }
    return bhp_error;
}

inline PetscErrorCode
attempt(
    const mpmc::flow_discretization_petsc::
        AdaptiveTimestepAttemptRequest3D&
            request,
    void* raw_context,
    mpmc::flow_discretization_petsc::
        AdaptiveTimestepAttemptResult3D*
            result) {
    namespace fdp =
        mpmc::flow_discretization_petsc;

    if (raw_context == nullptr ||
        result == nullptr) {
        return PETSC_ERR_ARG_NULL;
    }
    auto* context =
        static_cast<DriverContext3D*>(
            raw_context);
    if (context->reservoir_system == nullptr ||
        context->source_context == nullptr) {
        return PETSC_ERR_ARG_NULL;
    }

    PetscErrorCode error =
        clear_pending(context);
    if (error != PETSC_SUCCESS) {
        return error;
    }

    error =
        context->reservoir_system
            ->set_trial_timestep_seconds(
                request.timestep_seconds);
    if (error != PETSC_SUCCESS) {
        return error;
    }

    std::unique_ptr<
        FixedTotalMolarRateWellControlSystem3D>
        control_system;
    error =
        FixedTotalMolarRateWellControlSystem3D::
            create(
                context->comm,
                context->reservoir_system,
                context->source_context,
                context
                    ->entry_bottom_hole_pressure_pa,
                context
                    ->target_total_molar_rate_mol_per_s,
                &control_system);
    if (error != PETSC_SUCCESS ||
        control_system == nullptr) {
        return error != PETSC_SUCCESS
            ? error
            : PETSC_ERR_PLIB;
    }

    Vec solution = nullptr;
    std::optional<
        FixedTotalMolarRateWellControlSolveReport3D>
        solve_report;
    error =
        control_system->solve(
            &solution,
            &solve_report);

    if (error == PETSC_ERR_NOT_CONVERGED) {
        if (solution != nullptr) {
            (void)VecDestroy(&solution);
        }
        *result = {
            fdp::AdaptiveTimestepAttemptOutcome3D::
                nonlinear_solve_diverged,
            0,
            0,
            0,
            0,
            0U};
        return PETSC_SUCCESS;
    }
    if (error != PETSC_SUCCESS) {
        if (solution != nullptr) {
            (void)VecDestroy(&solution);
        }
        return error;
    }
    if (solution == nullptr ||
        !solve_report.has_value() ||
        !solve_report->converged()) {
        if (solution != nullptr) {
            (void)VecDestroy(&solution);
        }
        return PETSC_ERR_PLIB;
    }

    *result = {
        fdp::AdaptiveTimestepAttemptOutcome3D::
            stable_phase_set,
        solve_report->nonlinear_iterations,
        solve_report->function_domain_errors,
        solve_report->jacobian_domain_errors,
        solve_report
            ->line_search_direction_changes,
        0U};

    context->pending_control_system =
        std::move(control_system);
    context->pending_augmented_state =
        solution;
    context->pending_solve_report =
        std::move(solve_report);
    return PETSC_SUCCESS;
}

inline PetscErrorCode
commit(
    const mpmc::flow_discretization_petsc::
        AdaptiveTimestepAttemptRequest3D&
            request,
    const mpmc::flow_discretization_petsc::
        AdaptiveTimestepAttemptResult3D&
            result,
    void* raw_context) {
    namespace fdp =
        mpmc::flow_discretization_petsc;

    if (raw_context == nullptr) {
        return PETSC_ERR_ARG_NULL;
    }
    auto* context =
        static_cast<DriverContext3D*>(
            raw_context);
    if (context->reservoir_system == nullptr ||
        context->source_context == nullptr ||
        context->adaptive_options == nullptr ||
        context->clock == nullptr ||
        context
                ->accepted_bottom_hole_pressure_pa ==
            nullptr ||
        context->pending_control_system ==
            nullptr ||
        context->pending_augmented_state ==
            nullptr ||
        !context
             ->pending_solve_report
             .has_value()) {
        return PETSC_ERR_ARG_WRONGSTATE;
    }

    fdp::AdaptiveTimestepDecisionResult3D
        decision;
    try {
        decision =
            fdp::decide_adaptive_timestep_3d(
                request.timestep_seconds,
                request.retry_index,
                result,
                *context
                     ->adaptive_options);
    } catch (...) {
        return PETSC_ERR_ARG_INCOMP;
    }
    if (decision.decision !=
            fdp::AdaptiveTimestepDecision3D::
                accept_and_grow &&
        decision.decision !=
            fdp::AdaptiveTimestepDecision3D::
                accept_and_hold) {
        return PETSC_ERR_ARG_INCOMP;
    }

    Vec reservoir_state = nullptr;
    PetscErrorCode error =
        context->pending_control_system
            ->copy_reservoir_state(
                context
                    ->pending_augmented_state,
                &reservoir_state);
    if (error != PETSC_SUCCESS ||
        reservoir_state == nullptr) {
        if (reservoir_state != nullptr) {
            (void)VecDestroy(
                &reservoir_state);
        }
        return error != PETSC_SUCCESS
            ? error
            : PETSC_ERR_PLIB;
    }

    std::optional<
        fdp::AcceptedPhysicalTimestepRecord3D>
        accepted_record;
    error =
        fdp::commit_accepted_physical_timestep_3d(
            *context->reservoir_system,
            reservoir_state,
            request.timestep_seconds,
            result,
            decision,
            context->clock,
            &accepted_record);

    const PetscErrorCode destroy_error =
        VecDestroy(
            &reservoir_state);
    if (error == PETSC_SUCCESS &&
        destroy_error != PETSC_SUCCESS) {
        error =
            destroy_error;
    }
    if (error != PETSC_SUCCESS) {
        return error;
    }
    if (!accepted_record.has_value()) {
        return PETSC_ERR_PLIB;
    }

    const double accepted_bhp =
        context
            ->pending_solve_report
            ->bottom_hole_pressure_pa;
    if (!std::isfinite(accepted_bhp) ||
        !(accepted_bhp > 0.0)) {
        return PETSC_ERR_PLIB;
    }

    try {
        context->source_context
            ->begin_evaluation(
                accepted_bhp);
    } catch (...) {
        return PETSC_ERR_ARG_INCOMP;
    }

    *context
         ->accepted_bottom_hole_pressure_pa =
        accepted_bhp;
    context->accepted_record =
        std::move(accepted_record);
    context->accepted_solve =
        context->pending_solve_report;

    return clear_pending(context);
}

} // namespace fixed_total_molar_rate_timestep_driver_detail

/// Advance one frozen-phase physical timestep under a fixed whole-well total
/// molar production-rate target.
///
/// The reservoir state/history and physical clock retain their existing
/// ownership. BHP is a nonlinear unknown only inside each augmented trial
/// system and is committed separately as the next timestep's initial guess.
/// Every rejected retry is rebuilt from the same accepted reservoir baseline
/// and entry BHP. No BHP value is inserted into backward-Euler accumulation
/// history.
///
/// A terminal adaptive rejection is reported with PETSC_SUCCESS and
/// report->accepted()==false, matching the model-neutral adaptive controller.
/// The entry trial timestep and BHP evaluator state are restored in that case.
inline PetscErrorCode
advance_fixed_total_molar_rate_controlled_physical_timestep_3d(
    MPI_Comm comm,
    mpmc::flow_discretization_petsc::
        PhaseTransitionRebuiltNaturalVariableSystem3D*
            reservoir_system,
    FixedTotalMolarRateWellSourceEvaluatorContext3D*
        source_context,
    double target_total_molar_rate_mol_per_s,
    const FixedTotalMolarRatePhysicalTimestepDriverOptions3D&
        options,
    mpmc::flow_discretization_petsc::
        AcceptedPhysicalTimeClock3D*
            clock,
    double* accepted_bottom_hole_pressure_pa,
    std::optional<
        FixedTotalMolarRatePhysicalTimestepDriverReport3D>*
            report) {
    namespace fdp =
        mpmc::flow_discretization_petsc;
    using namespace
        fixed_total_molar_rate_timestep_driver_detail;

    if (reservoir_system == nullptr ||
        source_context == nullptr ||
        clock == nullptr ||
        accepted_bottom_hole_pressure_pa ==
            nullptr ||
        report == nullptr) {
        return PETSC_ERR_ARG_NULL;
    }
    report->reset();

    const double entry_bhp =
        *accepted_bottom_hole_pressure_pa;
    const double entry_dt =
        clock->next_timestep_seconds();
    if (!std::isfinite(entry_bhp) ||
        !(entry_bhp > 0.0) ||
        !std::isfinite(
            target_total_molar_rate_mol_per_s) ||
        !(target_total_molar_rate_mol_per_s >
          0.0)) {
        return PETSC_ERR_ARG_OUTOFRANGE;
    }
    if (!same_timestep(
            reservoir_system
                ->time_step_seconds(),
            entry_dt)) {
        return PETSC_ERR_ARG_INCOMP;
    }

    DriverContext3D context{
        comm,
        reservoir_system,
        source_context,
        target_total_molar_rate_mol_per_s,
        entry_bhp,
        &options.adaptive,
        clock,
        accepted_bottom_hole_pressure_pa};

    std::optional<
        fdp::AdaptiveTimestepControllerReport3D>
        adaptive_report;
    PetscErrorCode error =
        fdp::solve_adaptive_timestep_3d(
            entry_dt,
            options.adaptive,
            {
                &attempt,
                &context,
                &commit,
                &context},
            &adaptive_report);

    const PetscErrorCode clear_error =
        clear_pending(&context);
    if (error == PETSC_SUCCESS &&
        clear_error != PETSC_SUCCESS) {
        error =
            clear_error;
    }

    if (error != PETSC_SUCCESS) {
        return restore_entry(
            &context,
            entry_dt,
            error);
    }
    if (!adaptive_report.has_value()) {
        return restore_entry(
            &context,
            entry_dt,
            PETSC_ERR_PLIB);
    }

    if (!adaptive_report->accepted()) {
        const PetscErrorCode restore_error =
            restore_entry(
                &context,
                entry_dt,
                PETSC_SUCCESS);
        if (restore_error != PETSC_SUCCESS) {
            return restore_error;
        }
        report->emplace(
            FixedTotalMolarRatePhysicalTimestepDriverReport3D{
                std::move(
                    *adaptive_report),
                std::nullopt,
                std::nullopt,
                entry_bhp,
                entry_bhp});
        return PETSC_SUCCESS;
    }

    if (!context.accepted_record.has_value() ||
        !context.accepted_solve.has_value() ||
        !std::isfinite(
            *accepted_bottom_hole_pressure_pa) ||
        !(*accepted_bottom_hole_pressure_pa >
          0.0)) {
        return PETSC_ERR_PLIB;
    }

    report->emplace(
        FixedTotalMolarRatePhysicalTimestepDriverReport3D{
            std::move(
                *adaptive_report),
            std::move(
                context.accepted_record),
            std::move(
                context.accepted_solve),
            entry_bhp,
            *accepted_bottom_hole_pressure_pa});
    return PETSC_SUCCESS;
}

} // namespace mpmc::well_discretization_petsc

#endif // MPMC_WELL_DISCRETIZATION_PETSC_FIXED_TOTAL_MOLAR_RATE_TIMESTEP_DRIVER_HPP
