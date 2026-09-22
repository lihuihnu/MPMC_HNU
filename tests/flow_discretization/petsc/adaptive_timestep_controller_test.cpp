#include <mpmc/flow_discretization_petsc/adaptive_timestep_controller.hpp>

#include <petscsys.h>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <optional>
#include <stdexcept>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

namespace {

namespace fdp = mpmc::flow_discretization_petsc;

void require_adaptive_collective(
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
            PETSC_COMM_WORLD) != MPI_SUCCESS) {
        throw std::runtime_error(
            "MPI_Allreduce failed in adaptive timestep regression");
    }
    if (global == 0) {
        throw std::runtime_error(
            std::string{message});
    }
}

void near_adaptive_collective(
    double actual,
    double expected) {
    const double scale =
        std::max(
            {1.0,
             std::abs(actual),
             std::abs(expected)});
    require_adaptive_collective(
        std::isfinite(actual) &&
            std::isfinite(expected) &&
            std::abs(actual - expected) <=
                32.0 *
                    std::numeric_limits<double>::
                        epsilon() *
                    scale,
        "adaptive timestep numeric mismatch");
}

struct ScriptedAdaptiveContext {
    std::vector<
        fdp::AdaptiveTimestepAttemptResult3D>
        scripted;
    std::vector<double> observed_dt;
    std::size_t next{};
    std::size_t commits{};
    double committed_dt{};
    bool fail_attempt{};
    bool fail_commit{};
};

PetscErrorCode scripted_attempt(
    const fdp::AdaptiveTimestepAttemptRequest3D&
        request,
    void* raw_context,
    fdp::AdaptiveTimestepAttemptResult3D*
        result) {
    if (raw_context == nullptr ||
        result == nullptr) {
        return PETSC_ERR_ARG_NULL;
    }
    auto* context =
        static_cast<
            ScriptedAdaptiveContext*>(
                raw_context);
    if (context->fail_attempt) {
        return PETSC_ERR_LIB;
    }
    if (context->next >=
        context->scripted.size()) {
        return PETSC_ERR_ARG_OUTOFRANGE;
    }
    context->observed_dt.push_back(
        request.timestep_seconds);
    *result =
        context->scripted[
            context->next++];
    return PETSC_SUCCESS;
}

PetscErrorCode scripted_commit(
    const fdp::AdaptiveTimestepAttemptRequest3D&
        request,
    const fdp::AdaptiveTimestepAttemptResult3D&,
    void* raw_context) {
    if (raw_context == nullptr) {
        return PETSC_ERR_ARG_NULL;
    }
    auto* context =
        static_cast<
            ScriptedAdaptiveContext*>(
                raw_context);
    if (context->fail_commit) {
        return PETSC_ERR_LIB;
    }
    ++context->commits;
    context->committed_dt =
        request.timestep_seconds;
    return PETSC_SUCCESS;
}

fdp::AdaptiveTimestepControllerOptions3D
default_options() {
    fdp::AdaptiveTimestepControllerOptions3D
        options;
    options.minimum_timestep_seconds = 1.0;
    options.maximum_timestep_seconds = 16.0;
    options.cutback_factor = 0.5;
    options.growth_factor = 2.0;
    options.maximum_retries = 4U;
    options.growth_nonlinear_iteration_limit = 4;
    options
        .growth_line_search_direction_change_limit =
        0;
    options.growth_transition_restart_limit = 0U;
    return options;
}

fdp::AdaptiveTimestepAttemptResult3D
stable_easy() {
    return {
        fdp::AdaptiveTimestepAttemptOutcome3D::
            stable_phase_set,
        2,
        0,
        0,
        0,
        0U};
}

fdp::AdaptiveTimestepAttemptResult3D
nonlinear_domain() {
    return {
        fdp::AdaptiveTimestepAttemptOutcome3D::
            nonlinear_domain_error,
        1,
        1,
        0,
        0,
        0U};
}

fdp::AdaptiveTimestepAttemptResult3D
scan_indeterminate() {
    return {
        fdp::AdaptiveTimestepAttemptOutcome3D::
            phase_set_scan_indeterminate,
        3,
        0,
        0,
        0,
        0U};
}

fdp::AdaptiveTimestepAttemptResult3D
phase_cycle() {
    return {
        fdp::AdaptiveTimestepAttemptOutcome3D::
            phase_set_cycle_detected,
        4,
        0,
        0,
        1,
        1U};
}

void retry_cutback_growth_contract() {
    auto options =
        default_options();
    ScriptedAdaptiveContext context;
    context.scripted = {
        nonlinear_domain(),
        scan_indeterminate(),
        stable_easy()};

    std::optional<
        fdp::AdaptiveTimestepControllerReport3D>
        report;
    const PetscErrorCode error =
        fdp::solve_adaptive_timestep_3d(
            8.0,
            options,
            {
                &scripted_attempt,
                &context,
                &scripted_commit,
                &context},
            &report);

    require_adaptive_collective(
        error == PETSC_SUCCESS &&
            report.has_value() &&
            report->accepted() &&
            report->retries == 2U &&
            report->attempts.size() == 3U &&
            context.commits == 1U &&
            context.observed_dt.size() == 3U,
        "adaptive retry/cutback/growth lifecycle failed");
    near_adaptive_collective(
        context.observed_dt[0],
        8.0);
    near_adaptive_collective(
        context.observed_dt[1],
        4.0);
    near_adaptive_collective(
        context.observed_dt[2],
        2.0);
    near_adaptive_collective(
        *report->accepted_timestep_seconds,
        2.0);
    near_adaptive_collective(
        *report->next_timestep_seconds,
        4.0);
    near_adaptive_collective(
        context.committed_dt,
        2.0);

    require_adaptive_collective(
        report->attempts[0].decision ==
                fdp::AdaptiveTimestepDecision3D::
                    reject_and_cutback &&
            report->attempts[1].decision ==
                fdp::AdaptiveTimestepDecision3D::
                    reject_and_cutback &&
            report->attempts[2].decision ==
                fdp::AdaptiveTimestepDecision3D::
                    accept_and_grow,
        "adaptive timestep decisions changed");
}

void accepted_hard_step_holds_contract() {
    auto options =
        default_options();
    auto hard =
        stable_easy();
    hard.nonlinear_iterations = 7;
    hard.line_search_direction_changes = 1;
    hard.phase_transition_restarts = 1U;

    ScriptedAdaptiveContext context;
    context.scripted = {hard};

    std::optional<
        fdp::AdaptiveTimestepControllerReport3D>
        report;
    const PetscErrorCode error =
        fdp::solve_adaptive_timestep_3d(
            4.0,
            options,
            {
                &scripted_attempt,
                &context,
                &scripted_commit,
                &context},
            &report);

    require_adaptive_collective(
        error == PETSC_SUCCESS &&
            report.has_value() &&
            report->accepted() &&
            report->retries == 0U &&
            report->attempts.size() == 1U &&
            report->attempts.front().decision ==
                fdp::AdaptiveTimestepDecision3D::
                    accept_and_hold &&
            context.commits == 1U,
        "hard accepted timestep did not hold dt");
    near_adaptive_collective(
        *report->accepted_timestep_seconds,
        4.0);
    near_adaptive_collective(
        *report->next_timestep_seconds,
        4.0);
}

void retry_budget_contract() {
    auto options =
        default_options();
    options.maximum_retries = 1U;

    ScriptedAdaptiveContext context;
    context.scripted = {
        phase_cycle(),
        phase_cycle()};

    std::optional<
        fdp::AdaptiveTimestepControllerReport3D>
        report;
    const PetscErrorCode error =
        fdp::solve_adaptive_timestep_3d(
            8.0,
            options,
            {
                &scripted_attempt,
                &context,
                &scripted_commit,
                &context},
            &report);

    require_adaptive_collective(
        error == PETSC_SUCCESS &&
            report.has_value() &&
            !report->accepted() &&
            report->outcome ==
                fdp::AdaptiveTimestepControllerOutcome3D::
                    retry_budget_exhausted &&
            report->retries == 1U &&
            report->attempts.size() == 2U &&
            context.commits == 0U,
        "adaptive retry budget semantics mismatch");
    near_adaptive_collective(
        context.observed_dt[0],
        8.0);
    near_adaptive_collective(
        context.observed_dt[1],
        4.0);
}

void minimum_timestep_contract() {
    auto options =
        default_options();
    options.cutback_factor = 0.25;
    options.maximum_retries = 8U;

    ScriptedAdaptiveContext context;
    context.scripted = {
        phase_cycle(),
        phase_cycle()};

    std::optional<
        fdp::AdaptiveTimestepControllerReport3D>
        report;
    const PetscErrorCode error =
        fdp::solve_adaptive_timestep_3d(
            2.0,
            options,
            {
                &scripted_attempt,
                &context,
                &scripted_commit,
                &context},
            &report);

    require_adaptive_collective(
        error == PETSC_SUCCESS &&
            report.has_value() &&
            !report->accepted() &&
            report->outcome ==
                fdp::AdaptiveTimestepControllerOutcome3D::
                    minimum_timestep_reached &&
            report->retries == 1U &&
            report->attempts.size() == 2U &&
            context.commits == 0U,
        "adaptive minimum timestep semantics mismatch");
    near_adaptive_collective(
        context.observed_dt[0],
        2.0);
    near_adaptive_collective(
        context.observed_dt[1],
        1.0);
}

void fatal_callback_errors_propagate() {
    auto options =
        default_options();

    ScriptedAdaptiveContext context;
    context.scripted = {stable_easy()};
    context.fail_attempt = true;

    std::optional<
        fdp::AdaptiveTimestepControllerReport3D>
        report;
    require_adaptive_collective(
        fdp::solve_adaptive_timestep_3d(
            2.0,
            options,
            {
                &scripted_attempt,
                &context,
                &scripted_commit,
                &context},
            &report) ==
                PETSC_ERR_LIB &&
            !report.has_value() &&
            context.commits == 0U,
        "fatal adaptive attempt error was incorrectly converted to cutback");

    context.fail_attempt = false;
    context.fail_commit = true;
    context.next = 0U;
    context.observed_dt.clear();
    report.reset();
    require_adaptive_collective(
        fdp::solve_adaptive_timestep_3d(
            2.0,
            options,
            {
                &scripted_attempt,
                &context,
                &scripted_commit,
                &context},
            &report) ==
                PETSC_ERR_LIB &&
            !report.has_value(),
        "fatal adaptive commit error was swallowed");
}

void outcome_adapter_contract() {
    fdp::NaturalVariableSnesFailureDiagnostics3D
        domain_diagnostics;
    domain_diagnostics.snes_reason =
        SNES_DIVERGED_FUNCTION_DOMAIN;
    domain_diagnostics.function_domain_errors = 2;
    domain_diagnostics.nonlinear_iterations = 1;
    domain_diagnostics.line_search_direction_changes = 1;

    const auto domain =
        fdp::make_adaptive_timestep_attempt_result(
            domain_diagnostics);
    require_adaptive_collective(
        domain.outcome ==
                fdp::AdaptiveTimestepAttemptOutcome3D::
                    nonlinear_domain_error &&
            domain.function_domain_errors == 2,
        "SNES domain-error adaptive mapping failed");

    fdp::NaturalVariableSnesFailureDiagnostics3D
        divergence_diagnostics;
    divergence_diagnostics.snes_reason =
        SNES_DIVERGED_LINEAR_SOLVE;
    divergence_diagnostics.nonlinear_iterations = 0;
    const auto divergence =
        fdp::make_adaptive_timestep_attempt_result(
            divergence_diagnostics);
    require_adaptive_collective(
        divergence.outcome ==
            fdp::AdaptiveTimestepAttemptOutcome3D::
                nonlinear_solve_diverged,
        "SNES divergence adaptive mapping failed");

    fdp::PostSnesPhaseTransitionControllerReport3D
        phase_report;
    phase_report.outcome =
        fdp::PostSnesPhaseTransitionOutcome3D::
            phase_set_scan_indeterminate;
    phase_report.transition_restarts = 2U;
    fdp::PostSnesPhaseTransitionGenerationReport3D
        generation;
    generation.nonlinear_solve
        .nonlinear_iterations = 5;
    generation.nonlinear_solve
        .line_search_direction_changes = 1;
    phase_report.generations.push_back(
        std::move(generation));

    const auto phase =
        fdp::make_adaptive_timestep_attempt_result(
            phase_report);
    require_adaptive_collective(
        phase.outcome ==
                fdp::AdaptiveTimestepAttemptOutcome3D::
                    phase_set_scan_indeterminate &&
            phase.phase_transition_restarts == 2U &&
            phase.nonlinear_iterations == 5 &&
            phase.line_search_direction_changes == 1,
        "phase-scan adaptive mapping failed");

    phase_report.outcome =
        fdp::PostSnesPhaseTransitionOutcome3D::
            transition_restart_budget_exhausted;
    require_adaptive_collective(
        fdp::make_adaptive_timestep_attempt_result(
            phase_report)
                .outcome ==
            fdp::AdaptiveTimestepAttemptOutcome3D::
                transition_restart_budget_exhausted,
        "phase restart-budget adaptive mapping failed");

    phase_report.outcome =
        fdp::PostSnesPhaseTransitionOutcome3D::
            phase_set_cycle_detected;
    require_adaptive_collective(
        fdp::make_adaptive_timestep_attempt_result(
            phase_report)
                .outcome ==
            fdp::AdaptiveTimestepAttemptOutcome3D::
                phase_set_cycle_detected,
        "phase cycle adaptive mapping failed");
}

} // namespace

void adaptive_timestep_controller_test() {
    retry_cutback_growth_contract();
    accepted_hard_step_holds_contract();
    retry_budget_contract();
    minimum_timestep_contract();
    fatal_callback_errors_propagate();
    outcome_adapter_contract();
}
