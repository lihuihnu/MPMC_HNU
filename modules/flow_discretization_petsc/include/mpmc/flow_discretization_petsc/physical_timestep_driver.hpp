#ifndef MPMC_FLOW_DISCRETIZATION_PETSC_PHYSICAL_TIMESTEP_DRIVER_HPP
#define MPMC_FLOW_DISCRETIZATION_PETSC_PHYSICAL_TIMESTEP_DRIVER_HPP

#include <mpmc/flow_discretization_petsc/accepted_physical_time.hpp>
#include <mpmc/flow_discretization_petsc/post_snes_phase_transition_controller.hpp>

#include <petscvec.h>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <memory>
#include <optional>
#include <string_view>
#include <utility>
#include <vector>

namespace mpmc::flow_discretization_petsc {

inline constexpr std::string_view
    physical_timestep_driver_convention =
        "flow_discretization_petsc/physical-timestep-driver/v1";

struct PhysicalTimestepDriverOptions3D {
    AdaptiveTimestepControllerOptions3D adaptive;
    PostSnesPhaseTransitionControllerOptions3D
        phase_transition;

    /// Optional cap for the first nonlinear attempt of this physical timestep.
    /// This is intended for an outer time loop that must land on an exact
    /// schedule/end-time boundary without mutating the accepted clock first.
    std::optional<double>
        initial_timestep_cap_seconds;
};

struct PhysicalTimestepDriverReport3D {
    AdaptiveTimestepControllerReport3D adaptive;
    std::optional<AcceptedPhysicalTimestepRecord3D>
        accepted_record;

    [[nodiscard]] bool accepted() const noexcept {
        return adaptive.accepted();
    }
};

namespace physical_timestep_driver_detail {

[[nodiscard]] inline PetscErrorCode
destroy_state(
    Vec* state) noexcept {
    if (state == nullptr) {
        return PETSC_ERR_ARG_NULL;
    }
    if (*state == nullptr) {
        return PETSC_SUCCESS;
    }
    return VecDestroy(state);
}

[[nodiscard]] inline PetscErrorCode
restore_trial_timestep(
    std::unique_ptr<
        PhaseTransitionRebuiltNaturalVariableSystem3D>*
            accepted_system,
    double timestep_seconds,
    PetscErrorCode primary_error) noexcept {
    if (accepted_system == nullptr ||
        *accepted_system == nullptr) {
        return primary_error != PETSC_SUCCESS
            ? primary_error
            : PETSC_ERR_ARG_NULL;
    }
    const PetscErrorCode restore_error =
        (*accepted_system)
            ->set_trial_timestep_seconds(
                timestep_seconds);
    return primary_error != PETSC_SUCCESS
        ? primary_error
        : restore_error;
}

[[nodiscard]] inline bool
same_timestep(
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

[[nodiscard]] inline
AdaptiveTimestepAttemptResult3D
nonconverged_attempt(
    const VariableCardinalityNaturalVariableSnesSolveReport3D&
        report,
    const std::optional<
        NaturalVariableSnesFailureDiagnostics3D>&
        diagnostics,
    std::size_t phase_transition_restarts) {
    AdaptiveTimestepAttemptResult3D result;
    if (diagnostics.has_value()) {
        result =
            make_adaptive_timestep_attempt_result(
                *diagnostics);
    } else {
        result.outcome =
            AdaptiveTimestepAttemptOutcome3D::
                nonlinear_solve_diverged;
        result.nonlinear_iterations =
            std::max<PetscInt>(
                report.nonlinear_iterations,
                0);
        result.line_search_direction_changes =
            std::max<PetscInt>(
                report
                    .line_search_direction_changes,
                0);
    }
    result.phase_transition_restarts =
        phase_transition_restarts;
    return result;
}

} // namespace physical_timestep_driver_detail

/// Advance exactly one physical backward-Euler timestep from the supplied
/// accepted natural-variable system.
///
/// The accepted system remains the rollback anchor for every rejected attempt.
/// Trial dt changes never rewrite accepted state/history. If a converged trial
/// requests a phase-set change, the first target system is rebuilt separately
/// and all restarted topology generations are solved on that disposable trial
/// object. A failed transition trial is discarded before cutback, so the
/// accepted topology is preserved. Only a stable candidate reaches
/// commit_accepted_physical_timestep_3d(), which advances history and physical
/// time exactly once.
///
/// initial_phase_transition_restarts records transitions already materialized
/// before this entry but still belonging to this physical timestep. It is
/// included in adaptive effort/growth and in the total transition-restart
/// budget. Normal calls from an already accepted system pass zero.
inline PetscErrorCode
advance_one_physical_timestep_3d(
    MPI_Comm comm,
    std::unique_ptr<
        PhaseTransitionRebuiltNaturalVariableSystem3D>*
            accepted_system,
    std::size_t initial_phase_transition_restarts,
    PostSnesPhaseTransitionControllerBindings3D
        transition_bindings,
    const PhysicalTimestepDriverOptions3D& options,
    AcceptedPhysicalTimeClock3D* clock,
    std::optional<
        PhysicalTimestepDriverReport3D>*
            report) {
    using namespace
        physical_timestep_driver_detail;
    using post_snes_transition_detail::
        collective_error;
    using post_snes_transition_detail::
        gather_accepted_batch;
    using post_snes_transition_detail::
        gather_phase_signature;
    using post_snes_transition_detail::
        validate_local_proposals;
    using post_snes_transition_detail::
        validate_rebuilt_signature_against_batch;

    if (accepted_system == nullptr ||
        clock == nullptr ||
        report == nullptr ||
        transition_bindings.scanner == nullptr ||
        transition_bindings.rebuild_factory ==
            nullptr) {
        return PETSC_ERR_ARG_NULL;
    }
    report->reset();
    if (*accepted_system == nullptr ||
        initial_phase_transition_restarts >
            options.phase_transition
                .max_transition_restarts) {
        return PETSC_ERR_ARG_INCOMP;
    }

    const double entry_timestep_seconds =
        clock->next_timestep_seconds();
    if (!same_timestep(
            (*accepted_system)
                ->time_step_seconds(),
            entry_timestep_seconds)) {
        return PETSC_ERR_ARG_INCOMP;
    }

    double initial_timestep_seconds =
        entry_timestep_seconds;
    if (options
            .initial_timestep_cap_seconds
            .has_value()) {
        const double cap =
            *options
                 .initial_timestep_cap_seconds;
        if (!std::isfinite(cap) ||
            !(cap > 0.0)) {
            return PETSC_ERR_ARG_OUTOFRANGE;
        }
        initial_timestep_seconds =
            std::min(
                initial_timestep_seconds,
                cap);
    }

    try {
        adaptive_timestep_detail::
            validate_options(
                initial_timestep_seconds,
                options.adaptive);
    } catch (const std::exception&) {
        return PETSC_ERR_ARG_OUTOFRANGE;
    }

    AdaptiveTimestepControllerReport3D
        adaptive;
    adaptive.initial_timestep_seconds =
        initial_timestep_seconds;

    double timestep_seconds =
        initial_timestep_seconds;
    std::size_t retries = 0U;

    for (std::size_t attempt_index = 0U;;
         ++attempt_index) {
        const AdaptiveTimestepAttemptRequest3D
            request{
                attempt_index,
                retries,
                timestep_seconds};

        PetscErrorCode error =
            (*accepted_system)
                ->set_trial_timestep_seconds(
                    timestep_seconds);
        if (error != PETSC_SUCCESS) {
            return restore_trial_timestep(
                accepted_system,
                entry_timestep_seconds,
                error);
        }

        Vec solved = nullptr;
        std::optional<
            VariableCardinalityNaturalVariableSnesSolveReport3D>
            solve_report;
        std::optional<
            NaturalVariableSnesFailureDiagnostics3D>
            failure_diagnostics;
        error =
            (*accepted_system)
                ->solve(
                    &solved,
                    &solve_report,
                    &failure_diagnostics);
        if (error != PETSC_SUCCESS) {
            (void)destroy_state(
                &solved);
            return restore_trial_timestep(
                accepted_system,
                entry_timestep_seconds,
                error);
        }
        if (solved == nullptr ||
            !solve_report.has_value()) {
            (void)destroy_state(
                &solved);
            return restore_trial_timestep(
                accepted_system,
                entry_timestep_seconds,
                PETSC_ERR_PLIB);
        }

        AdaptiveTimestepAttemptResult3D
            attempt_result;
        std::unique_ptr<
            PhaseTransitionRebuiltNaturalVariableSystem3D>
            transitioned_candidate;
        Vec candidate_state = nullptr;

        if (static_cast<int>(
                solve_report
                    ->converged_reason) <= 0) {
            attempt_result =
                nonconverged_attempt(
                    *solve_report,
                    failure_diagnostics,
                    initial_phase_transition_restarts);
            error =
                destroy_state(
                    &solved);
            if (error != PETSC_SUCCESS) {
                return restore_trial_timestep(
                    accepted_system,
                    initial_timestep_seconds,
                    error);
            }
        } else {
            attempt_result =
                make_adaptive_timestep_attempt_result(
                    *solve_report,
                    initial_phase_transition_restarts);

            PostSnesPhaseTransitionScanStatus3D
                local_scan_status =
                    PostSnesPhaseTransitionScanStatus3D::
                        complete;
            std::vector<
                PostSnesPhaseTransitionProposal3D>
                local_proposals;
            PetscErrorCode local_error =
                transition_bindings.scanner(
                    **accepted_system,
                    solved,
                    *solve_report,
                    transition_bindings
                        .scanner_context,
                    &local_scan_status,
                    &local_proposals);
            error =
                collective_error(
                    comm,
                    local_error);
            if (error != PETSC_SUCCESS) {
                (void)destroy_state(
                    &solved);
                return restore_trial_timestep(
                    accepted_system,
                    initial_timestep_seconds,
                    error);
            }

            const int local_indeterminate =
                local_scan_status ==
                        PostSnesPhaseTransitionScanStatus3D::
                            indeterminate
                    ? 1
                    : 0;
            int global_indeterminate = 0;
            if (MPI_Allreduce(
                    &local_indeterminate,
                    &global_indeterminate,
                    1,
                    MPI_INT,
                    MPI_MAX,
                    comm) != MPI_SUCCESS) {
                (void)destroy_state(
                    &solved);
                return restore_trial_timestep(
                    accepted_system,
                    initial_timestep_seconds,
                    PETSC_ERR_MPI);
            }

            if (global_indeterminate != 0) {
                attempt_result.outcome =
                    AdaptiveTimestepAttemptOutcome3D::
                        phase_set_scan_indeterminate;
                error =
                    destroy_state(
                        &solved);
                if (error != PETSC_SUCCESS) {
                    return restore_trial_timestep(
                        accepted_system,
                        initial_timestep_seconds,
                        error);
                }
            } else {
                local_error =
                    validate_local_proposals(
                        (*accepted_system)
                            ->numbering(),
                        &local_proposals);
                error =
                    collective_error(
                        comm,
                        local_error);
                if (error != PETSC_SUCCESS) {
                    (void)destroy_state(
                        &solved);
                    return restore_trial_timestep(
                        accepted_system,
                        initial_timestep_seconds,
                        error);
                }

                std::vector<
                    AcceptedPhaseTransitionSummary3D>
                    accepted_batch;
                error =
                    gather_accepted_batch(
                        comm,
                        (*accepted_system)
                            ->numbering(),
                        local_proposals,
                        &accepted_batch);
                if (error != PETSC_SUCCESS) {
                    (void)destroy_state(
                        &solved);
                    return restore_trial_timestep(
                        accepted_system,
                        initial_timestep_seconds,
                        error);
                }

                if (accepted_batch.empty()) {
                    candidate_state =
                        solved;
                    solved = nullptr;
                } else if (
                    initial_phase_transition_restarts >=
                    options.phase_transition
                        .max_transition_restarts) {
                    attempt_result.outcome =
                        AdaptiveTimestepAttemptOutcome3D::
                            transition_restart_budget_exhausted;
                    error =
                        destroy_state(
                            &solved);
                    if (error != PETSC_SUCCESS) {
                        return restore_trial_timestep(
                            accepted_system,
                            initial_timestep_seconds,
                            error);
                    }
                } else {
                    std::unique_ptr<
                        PhaseTransitionRebuiltNaturalVariableSystem3D>
                        first_rebuilt;
                    local_error =
                        transition_bindings
                            .rebuild_factory(
                                **accepted_system,
                                solved,
                                *solve_report,
                                local_proposals,
                                accepted_batch,
                                transition_bindings
                                    .rebuild_context,
                                &first_rebuilt);
                    error =
                        collective_error(
                            comm,
                            local_error);
                    const PetscErrorCode
                        solved_destroy =
                            destroy_state(
                                &solved);
                    if (error == PETSC_SUCCESS &&
                        solved_destroy !=
                            PETSC_SUCCESS) {
                        error =
                            solved_destroy;
                    }
                    if (error != PETSC_SUCCESS) {
                        return restore_trial_timestep(
                            accepted_system,
                            initial_timestep_seconds,
                            error);
                    }
                    if (first_rebuilt == nullptr) {
                        return restore_trial_timestep(
                            accepted_system,
                            initial_timestep_seconds,
                            PETSC_ERR_PLIB);
                    }

                    std::vector<
                        GlobalPhaseSetSignatureEntry3D>
                        current_signature;
                    std::vector<
                        GlobalPhaseSetSignatureEntry3D>
                        rebuilt_signature;
                    error =
                        gather_phase_signature(
                            comm,
                            **accepted_system,
                            &current_signature);
                    if (error == PETSC_SUCCESS) {
                        error =
                            gather_phase_signature(
                                comm,
                                *first_rebuilt,
                                &rebuilt_signature);
                    }
                    if (error == PETSC_SUCCESS) {
                        error =
                            validate_rebuilt_signature_against_batch(
                                current_signature,
                                accepted_batch,
                                rebuilt_signature);
                    }
                    if (error != PETSC_SUCCESS) {
                        return restore_trial_timestep(
                            accepted_system,
                            initial_timestep_seconds,
                            error);
                    }

                    PostSnesPhaseTransitionControllerOptions3D
                        remaining_transition_options =
                            options.phase_transition;
                    remaining_transition_options
                        .max_transition_restarts =
                        options.phase_transition
                            .max_transition_restarts -
                        initial_phase_transition_restarts -
                        1U;

                    Vec final_state = nullptr;
                    std::optional<
                        PostSnesPhaseTransitionControllerReport3D>
                        transition_report;
                    error =
                        solve_nonlinear_timestep_with_phase_transitions_3d(
                            comm,
                            std::move(
                                first_rebuilt),
                            transition_bindings,
                            remaining_transition_options,
                            &transitioned_candidate,
                            &final_state,
                            &transition_report);
                    if (error != PETSC_SUCCESS) {
                        (void)destroy_state(
                            &final_state);
                        return restore_trial_timestep(
                            accepted_system,
                            initial_timestep_seconds,
                            error);
                    }
                    if (!transition_report.has_value() ||
                        transitioned_candidate ==
                            nullptr ||
                        final_state == nullptr) {
                        (void)destroy_state(
                            &final_state);
                        return restore_trial_timestep(
                            accepted_system,
                            initial_timestep_seconds,
                            PETSC_ERR_PLIB);
                    }

                    attempt_result =
                        make_adaptive_timestep_attempt_result(
                            *transition_report);
                    const std::size_t
                        completed_before_outer =
                            initial_phase_transition_restarts +
                            1U;
                    if (attempt_result
                            .phase_transition_restarts >
                        options.phase_transition
                                .max_transition_restarts -
                            completed_before_outer) {
                        (void)destroy_state(
                            &final_state);
                        return restore_trial_timestep(
                            accepted_system,
                            initial_timestep_seconds,
                            PETSC_ERR_PLIB);
                    }
                    attempt_result
                        .phase_transition_restarts +=
                        completed_before_outer;

                    if (transition_report
                            ->timestep_accepted()) {
                        candidate_state =
                            final_state;
                        final_state =
                            nullptr;
                    } else {
                        error =
                            destroy_state(
                                &final_state);
                        transitioned_candidate
                            .reset();
                        if (error != PETSC_SUCCESS) {
                            return restore_trial_timestep(
                                accepted_system,
                                initial_timestep_seconds,
                                error);
                        }
                    }
                }
            }
        }

        AdaptiveTimestepDecisionResult3D
            decision;
        try {
            decision =
                decide_adaptive_timestep_3d(
                    timestep_seconds,
                    retries,
                    attempt_result,
                    options.adaptive);
        } catch (const std::exception&) {
            (void)destroy_state(
                &candidate_state);
            transitioned_candidate.reset();
            return restore_trial_timestep(
                accepted_system,
                entry_timestep_seconds,
                PETSC_ERR_ARG_INCOMP);
        }

        adaptive.attempts.push_back(
            {
                request,
                attempt_result,
                decision.decision,
                decision.next_timestep_seconds});

        if (decision.decision ==
                AdaptiveTimestepDecision3D::
                    accept_and_grow ||
            decision.decision ==
                AdaptiveTimestepDecision3D::
                    accept_and_hold) {
            if (candidate_state == nullptr) {
                return restore_trial_timestep(
                    accepted_system,
                    initial_timestep_seconds,
                    PETSC_ERR_PLIB);
            }

            std::optional<
                AcceptedPhysicalTimestepRecord3D>
                accepted_record;
            PhaseTransitionRebuiltNaturalVariableSystem3D*
                system_to_commit =
                    transitioned_candidate != nullptr
                    ? transitioned_candidate.get()
                    : accepted_system->get();
            error =
                commit_accepted_physical_timestep_3d(
                    *system_to_commit,
                    candidate_state,
                    timestep_seconds,
                    attempt_result,
                    decision,
                    clock,
                    &accepted_record);
            if (error == PETSC_SUCCESS &&
                transitioned_candidate != nullptr) {
                *accepted_system =
                    std::move(
                        transitioned_candidate);
            }
            const PetscErrorCode destroy_error =
                destroy_state(
                    &candidate_state);
            if (error == PETSC_SUCCESS &&
                destroy_error != PETSC_SUCCESS) {
                error = destroy_error;
            }
            if (error != PETSC_SUCCESS) {
                return error;
            }
            if (!accepted_record.has_value()) {
                return PETSC_ERR_PLIB;
            }

            adaptive.outcome =
                AdaptiveTimestepControllerOutcome3D::
                    timestep_accepted;
            adaptive.accepted_timestep_seconds =
                timestep_seconds;
            adaptive.next_timestep_seconds =
                decision.next_timestep_seconds;
            adaptive.retries =
                retries;

            PhysicalTimestepDriverReport3D
                completed;
            completed.adaptive =
                std::move(adaptive);
            completed.accepted_record =
                std::move(
                    accepted_record);
            report->emplace(
                std::move(completed));
            return PETSC_SUCCESS;
        }

        (void)destroy_state(
            &candidate_state);
        transitioned_candidate.reset();

        if (decision.decision ==
            AdaptiveTimestepDecision3D::
                reject_retry_budget_exhausted) {
            adaptive.outcome =
                AdaptiveTimestepControllerOutcome3D::
                    retry_budget_exhausted;
            adaptive.retries =
                retries;
            const PetscErrorCode restore_error =
                restore_trial_timestep(
                    accepted_system,
                    initial_timestep_seconds,
                    PETSC_SUCCESS);
            if (restore_error != PETSC_SUCCESS) {
                return restore_error;
            }
            PhysicalTimestepDriverReport3D
                completed;
            completed.adaptive =
                std::move(adaptive);
            report->emplace(
                std::move(completed));
            return PETSC_SUCCESS;
        }

        if (decision.decision ==
            AdaptiveTimestepDecision3D::
                reject_minimum_timestep_reached) {
            adaptive.outcome =
                AdaptiveTimestepControllerOutcome3D::
                    minimum_timestep_reached;
            adaptive.retries =
                retries;
            const PetscErrorCode restore_error =
                restore_trial_timestep(
                    accepted_system,
                    initial_timestep_seconds,
                    PETSC_SUCCESS);
            if (restore_error != PETSC_SUCCESS) {
                return restore_error;
            }
            PhysicalTimestepDriverReport3D
                completed;
            completed.adaptive =
                std::move(adaptive);
            report->emplace(
                std::move(completed));
            return PETSC_SUCCESS;
        }

        if (decision.decision !=
                AdaptiveTimestepDecision3D::
                    reject_and_cutback ||
            !decision.next_timestep_seconds
                 .has_value()) {
            return restore_trial_timestep(
                accepted_system,
                entry_timestep_seconds,
                PETSC_ERR_PLIB);
        }

        timestep_seconds =
            *decision
                 .next_timestep_seconds;
        ++retries;
    }
}

} // namespace mpmc::flow_discretization_petsc

#endif // MPMC_FLOW_DISCRETIZATION_PETSC_PHYSICAL_TIMESTEP_DRIVER_HPP
