#ifndef MPMC_FLOW_DISCRETIZATION_PETSC_ADAPTIVE_TIMESTEP_CONTROLLER_HPP
#define MPMC_FLOW_DISCRETIZATION_PETSC_ADAPTIVE_TIMESTEP_CONTROLLER_HPP

#include <mpmc/flow_discretization_petsc/natural_variable_snes_solver.hpp>
#include <mpmc/flow_discretization_petsc/post_snes_phase_transition_controller.hpp>

#include <petscsys.h>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <optional>
#include <stdexcept>
#include <string_view>
#include <utility>
#include <vector>

namespace mpmc::flow_discretization_petsc {

inline constexpr std::string_view
    adaptive_timestep_controller_convention =
        "flow_discretization_petsc/adaptive-timestep-controller/v1";

enum class AdaptiveTimestepAttemptOutcome3D {
    stable_phase_set,
    phase_transition_proposed,
    nonlinear_solve_diverged,
    nonlinear_domain_error,
    phase_set_scan_indeterminate,
    transition_restart_budget_exhausted,
    phase_set_cycle_detected
};

enum class AdaptiveTimestepControllerOutcome3D {
    timestep_accepted,
    phase_transition_handoff_required,
    retry_budget_exhausted,
    minimum_timestep_reached
};

enum class AdaptiveTimestepDecision3D {
    accept_and_grow,
    accept_and_hold,
    handoff_phase_transition,
    reject_and_cutback,
    reject_retry_budget_exhausted,
    reject_minimum_timestep_reached
};

struct AdaptiveTimestepControllerOptions3D {
    double minimum_timestep_seconds{1.0e-6};
    double maximum_timestep_seconds{1.0e30};
    double cutback_factor{0.5};
    double growth_factor{2.0};
    std::size_t maximum_retries{8U};

    /// Growth is deliberately conservative. A successful step grows only when
    /// all three effort gates are satisfied.
    PetscInt growth_nonlinear_iteration_limit{4};
    PetscInt growth_line_search_direction_change_limit{0};
    std::size_t growth_transition_restart_limit{0U};
};

struct AdaptiveTimestepAttemptRequest3D {
    std::size_t attempt_index{};
    std::size_t retry_index{};
    double timestep_seconds{};
};

struct AdaptiveTimestepAttemptResult3D {
    AdaptiveTimestepAttemptOutcome3D outcome{
        AdaptiveTimestepAttemptOutcome3D::
            nonlinear_solve_diverged};
    PetscInt nonlinear_iterations{};
    PetscInt function_domain_errors{};
    PetscInt jacobian_domain_errors{};
    PetscInt line_search_direction_changes{};
    std::size_t phase_transition_restarts{};
};

struct AdaptiveTimestepAttemptReport3D {
    AdaptiveTimestepAttemptRequest3D request;
    AdaptiveTimestepAttemptResult3D result;
    AdaptiveTimestepDecision3D decision{
        AdaptiveTimestepDecision3D::
            reject_and_cutback};
    std::optional<double>
        next_attempt_timestep_seconds;
};

struct AdaptiveTimestepControllerReport3D {
    AdaptiveTimestepControllerOutcome3D outcome{
        AdaptiveTimestepControllerOutcome3D::
            retry_budget_exhausted};
    double initial_timestep_seconds{};
    std::optional<double>
        accepted_timestep_seconds;
    std::optional<double>
        next_timestep_seconds;
    std::size_t retries{};
    std::vector<AdaptiveTimestepAttemptReport3D>
        attempts;

    [[nodiscard]] bool accepted() const noexcept {
        return outcome ==
            AdaptiveTimestepControllerOutcome3D::
                timestep_accepted;
    }
};

using AdaptiveTimestepAttemptEvaluator3D =
    PetscErrorCode (*)(
        const AdaptiveTimestepAttemptRequest3D& request,
        void* user_context,
        AdaptiveTimestepAttemptResult3D* result);

/// Called exactly once and only after a stable accepted attempt. The attempt
/// context may retain the pending converged state/system; commit is the only
/// point at which accepted-history ownership may advance.
using AdaptiveTimestepAcceptedCommit3D =
    PetscErrorCode (*)(
        const AdaptiveTimestepAttemptRequest3D& request,
        const AdaptiveTimestepAttemptResult3D& result,
        void* user_context);

struct AdaptiveTimestepControllerBindings3D {
    AdaptiveTimestepAttemptEvaluator3D attempt{};
    void* attempt_context{};
    AdaptiveTimestepAcceptedCommit3D commit{};
    void* commit_context{};
};

namespace adaptive_timestep_detail {

inline void validate_options(
    double initial_timestep_seconds,
    const AdaptiveTimestepControllerOptions3D&
        options) {
    if (!std::isfinite(initial_timestep_seconds) ||
        !(initial_timestep_seconds > 0.0) ||
        !std::isfinite(
            options.minimum_timestep_seconds) ||
        !(options.minimum_timestep_seconds > 0.0) ||
        !std::isfinite(
            options.maximum_timestep_seconds) ||
        !(options.maximum_timestep_seconds >=
          options.minimum_timestep_seconds) ||
        initial_timestep_seconds <
            options.minimum_timestep_seconds ||
        initial_timestep_seconds >
            options.maximum_timestep_seconds ||
        !std::isfinite(options.cutback_factor) ||
        !(options.cutback_factor > 0.0) ||
        !(options.cutback_factor < 1.0) ||
        !std::isfinite(options.growth_factor) ||
        !(options.growth_factor > 1.0) ||
        options.growth_nonlinear_iteration_limit < 0 ||
        options.growth_line_search_direction_change_limit <
            0) {
        throw std::invalid_argument(
            "mpmc::flow_discretization_petsc: invalid adaptive timestep options");
    }
}

inline void validate_attempt_result(
    const AdaptiveTimestepAttemptResult3D&
        result) {
    if (result.nonlinear_iterations < 0 ||
        result.function_domain_errors < 0 ||
        result.jacobian_domain_errors < 0 ||
        result.line_search_direction_changes < 0) {
        throw std::invalid_argument(
            "mpmc::flow_discretization_petsc: adaptive timestep attempt metrics must be nonnegative");
    }

    if ((result.outcome ==
             AdaptiveTimestepAttemptOutcome3D::
                 stable_phase_set ||
         result.outcome ==
             AdaptiveTimestepAttemptOutcome3D::
                 phase_transition_proposed) &&
        (result.function_domain_errors != 0 ||
         result.jacobian_domain_errors != 0)) {
        throw std::invalid_argument(
            "mpmc::flow_discretization_petsc: converged timestep/transition handoff cannot contain nonlinear domain errors");
    }

    if (result.outcome ==
            AdaptiveTimestepAttemptOutcome3D::
                nonlinear_domain_error &&
        result.function_domain_errors == 0 &&
        result.jacobian_domain_errors == 0) {
        throw std::invalid_argument(
            "mpmc::flow_discretization_petsc: nonlinear-domain timestep outcome requires domain-error evidence");
    }
}

[[nodiscard]] inline bool easy_acceptance(
    const AdaptiveTimestepAttemptResult3D& result,
    const AdaptiveTimestepControllerOptions3D&
        options) noexcept {
    return result.outcome ==
               AdaptiveTimestepAttemptOutcome3D::
                   stable_phase_set &&
        result.nonlinear_iterations <=
            options.growth_nonlinear_iteration_limit &&
        result.line_search_direction_changes <=
            options
                .growth_line_search_direction_change_limit &&
        result.phase_transition_restarts <=
            options.growth_transition_restart_limit &&
        result.function_domain_errors == 0 &&
        result.jacobian_domain_errors == 0;
}

[[nodiscard]] inline double grown_timestep(
    double accepted_timestep_seconds,
    const AdaptiveTimestepControllerOptions3D&
        options) {
    if (accepted_timestep_seconds >
        options.maximum_timestep_seconds /
            options.growth_factor) {
        return options.maximum_timestep_seconds;
    }
    return std::min(
        options.maximum_timestep_seconds,
        accepted_timestep_seconds *
            options.growth_factor);
}

[[nodiscard]] inline std::optional<double>
cutback_timestep(
    double current_timestep_seconds,
    const AdaptiveTimestepControllerOptions3D&
        options) {
    if (current_timestep_seconds <=
        options.minimum_timestep_seconds) {
        return std::nullopt;
    }

    double candidate =
        current_timestep_seconds *
        options.cutback_factor;
    if (!std::isfinite(candidate) ||
        !(candidate > 0.0) ||
        !(candidate < current_timestep_seconds)) {
        throw std::range_error(
            "mpmc::flow_discretization_petsc: adaptive timestep cutback is not representable");
    }
    candidate =
        std::max(
            candidate,
            options.minimum_timestep_seconds);
    if (!(candidate < current_timestep_seconds)) {
        return std::nullopt;
    }
    return candidate;
}

} // namespace adaptive_timestep_detail

struct AdaptiveTimestepDecisionResult3D {
    AdaptiveTimestepDecision3D decision{
        AdaptiveTimestepDecision3D::
            reject_and_cutback};
    std::optional<double>
        next_timestep_seconds;
};

[[nodiscard]] inline
AdaptiveTimestepDecisionResult3D
decide_adaptive_timestep_3d(
    double current_timestep_seconds,
    std::size_t retries_already_used,
    const AdaptiveTimestepAttemptResult3D& result,
    const AdaptiveTimestepControllerOptions3D&
        options) {
    adaptive_timestep_detail::validate_options(
        current_timestep_seconds,
        options);
    adaptive_timestep_detail::
        validate_attempt_result(result);

    if (result.outcome ==
        AdaptiveTimestepAttemptOutcome3D::
            phase_transition_proposed) {
        return {
            AdaptiveTimestepDecision3D::
                handoff_phase_transition,
            std::nullopt};
    }

    if (result.outcome ==
        AdaptiveTimestepAttemptOutcome3D::
            stable_phase_set) {
        if (adaptive_timestep_detail::
                easy_acceptance(
                    result,
                    options)) {
            return {
                AdaptiveTimestepDecision3D::
                    accept_and_grow,
                adaptive_timestep_detail::
                    grown_timestep(
                        current_timestep_seconds,
                        options)};
        }
        return {
            AdaptiveTimestepDecision3D::
                accept_and_hold,
            current_timestep_seconds};
    }

    if (retries_already_used >=
        options.maximum_retries) {
        return {
            AdaptiveTimestepDecision3D::
                reject_retry_budget_exhausted,
            std::nullopt};
    }

    const auto cutback =
        adaptive_timestep_detail::
            cutback_timestep(
                current_timestep_seconds,
                options);
    if (!cutback) {
        return {
            AdaptiveTimestepDecision3D::
                reject_minimum_timestep_reached,
            std::nullopt};
    }
    return {
        AdaptiveTimestepDecision3D::
            reject_and_cutback,
        *cutback};
}

[[nodiscard]] inline
AdaptiveTimestepAttemptResult3D
make_adaptive_timestep_attempt_result(
    const NaturalVariableSnesSolveReport3D&
        report) {
    if (static_cast<int>(
            report.converged_reason()) <= 0) {
        throw std::invalid_argument(
            "mpmc::flow_discretization_petsc: successful fixed-cardinality SNES report is not converged");
    }
    return {
        AdaptiveTimestepAttemptOutcome3D::
            stable_phase_set,
        report.nonlinear_iterations(),
        report.function_domain_errors(),
        report.jacobian_domain_errors(),
        report.line_search_direction_changes(),
        0U};
}

[[nodiscard]] inline
AdaptiveTimestepAttemptResult3D
make_adaptive_timestep_attempt_result(
    const NaturalVariableSnesFailureDiagnostics3D&
        diagnostics) {
    const bool domain =
        diagnostics.function_domain_errors > 0 ||
        diagnostics.jacobian_domain_errors > 0;
    return {
        domain
            ? AdaptiveTimestepAttemptOutcome3D::
                  nonlinear_domain_error
            : AdaptiveTimestepAttemptOutcome3D::
                  nonlinear_solve_diverged,
        std::max<PetscInt>(
            diagnostics.nonlinear_iterations,
            0),
        std::max<PetscInt>(
            diagnostics.function_domain_errors,
            0),
        std::max<PetscInt>(
            diagnostics.jacobian_domain_errors,
            0),
        std::max<PetscInt>(
            diagnostics.line_search_direction_changes,
            0),
        0U};
}

[[nodiscard]] inline
AdaptiveTimestepAttemptResult3D
make_adaptive_timestep_attempt_result(
    const VariableCardinalityNaturalVariableSnesSolveReport3D&
        report,
    std::size_t phase_transition_restarts = 0U) {
    if (static_cast<int>(
            report.converged_reason) <= 0 ||
        report.nonlinear_iterations < 0 ||
        report.line_search_direction_changes < 0 ||
        !std::isfinite(
            report.final_function_l2_norm)) {
        throw std::invalid_argument(
            "mpmc::flow_discretization_petsc: successful variable-cardinality SNES report is invalid");
    }
    return {
        AdaptiveTimestepAttemptOutcome3D::
            stable_phase_set,
        report.nonlinear_iterations,
        0,
        0,
        report.line_search_direction_changes,
        phase_transition_restarts};
}

[[nodiscard]] inline
AdaptiveTimestepAttemptResult3D
make_adaptive_timestep_attempt_result(
    const PostSnesPhaseTransitionControllerReport3D&
        report) {
    AdaptiveTimestepAttemptResult3D result;
    result.phase_transition_restarts =
        report.transition_restarts;

    for (const auto& generation :
         report.generations) {
        result.nonlinear_iterations =
            std::max(
                result.nonlinear_iterations,
                generation
                    .nonlinear_solve
                    .nonlinear_iterations);
        result.line_search_direction_changes +=
            generation
                .nonlinear_solve
                .line_search_direction_changes;
    }

    switch (report.outcome) {
    case PostSnesPhaseTransitionOutcome3D::
        stable_phase_set:
        result.outcome =
            AdaptiveTimestepAttemptOutcome3D::
                stable_phase_set;
        break;
    case PostSnesPhaseTransitionOutcome3D::
        nonlinear_solve_diverged:
        result.outcome =
            AdaptiveTimestepAttemptOutcome3D::
                nonlinear_solve_diverged;
        break;
    case PostSnesPhaseTransitionOutcome3D::
        phase_set_scan_indeterminate:
        result.outcome =
            AdaptiveTimestepAttemptOutcome3D::
                phase_set_scan_indeterminate;
        break;
    case PostSnesPhaseTransitionOutcome3D::
        transition_restart_budget_exhausted:
        result.outcome =
            AdaptiveTimestepAttemptOutcome3D::
                transition_restart_budget_exhausted;
        break;
    case PostSnesPhaseTransitionOutcome3D::
        phase_set_cycle_detected:
        result.outcome =
            AdaptiveTimestepAttemptOutcome3D::
                phase_set_cycle_detected;
        break;
    }
    return result;
}

/// Repeatedly attempt one physical timestep. Every retry uses the same accepted
/// history owned by the caller's attempt context; only the stable accepted
/// attempt reaches the commit callback. The controller never mutates dt inside
/// SNES callbacks and never commits a rejected trial state.
inline PetscErrorCode
solve_adaptive_timestep_3d(
    double initial_timestep_seconds,
    const AdaptiveTimestepControllerOptions3D&
        options,
    AdaptiveTimestepControllerBindings3D
        bindings,
    std::optional<
        AdaptiveTimestepControllerReport3D>*
            report) {
    if (report == nullptr ||
        bindings.attempt == nullptr ||
        bindings.commit == nullptr) {
        return PETSC_ERR_ARG_NULL;
    }
    report->reset();

    try {
        adaptive_timestep_detail::
            validate_options(
                initial_timestep_seconds,
                options);
    } catch (const std::exception&) {
        return PETSC_ERR_ARG_OUTOFRANGE;
    }

    AdaptiveTimestepControllerReport3D
        completed;
    completed.initial_timestep_seconds =
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
        AdaptiveTimestepAttemptResult3D result;
        const PetscErrorCode attempt_error =
            bindings.attempt(
                request,
                bindings.attempt_context,
                &result);
        if (attempt_error != PETSC_SUCCESS) {
            return attempt_error;
        }

        AdaptiveTimestepDecisionResult3D
            decision;
        try {
            decision =
                decide_adaptive_timestep_3d(
                    timestep_seconds,
                    retries,
                    result,
                    options);
        } catch (const std::exception&) {
            return PETSC_ERR_ARG_INCOMP;
        }

        completed.attempts.push_back(
            {
                request,
                result,
                decision.decision,
                decision.next_timestep_seconds});

        if (decision.decision ==
            AdaptiveTimestepDecision3D::
                handoff_phase_transition) {
            completed.outcome =
                AdaptiveTimestepControllerOutcome3D::
                    phase_transition_handoff_required;
            completed.retries =
                retries;
            report->emplace(
                std::move(completed));
            return PETSC_SUCCESS;
        }

        if (decision.decision ==
                AdaptiveTimestepDecision3D::
                    accept_and_grow ||
            decision.decision ==
                AdaptiveTimestepDecision3D::
                    accept_and_hold) {
            const PetscErrorCode commit_error =
                bindings.commit(
                    request,
                    result,
                    bindings.commit_context);
            if (commit_error != PETSC_SUCCESS) {
                return commit_error;
            }

            completed.outcome =
                AdaptiveTimestepControllerOutcome3D::
                    timestep_accepted;
            completed.accepted_timestep_seconds =
                timestep_seconds;
            completed.next_timestep_seconds =
                decision.next_timestep_seconds;
            completed.retries =
                retries;
            report->emplace(
                std::move(completed));
            return PETSC_SUCCESS;
        }

        if (decision.decision ==
            AdaptiveTimestepDecision3D::
                reject_retry_budget_exhausted) {
            completed.outcome =
                AdaptiveTimestepControllerOutcome3D::
                    retry_budget_exhausted;
            completed.retries =
                retries;
            report->emplace(
                std::move(completed));
            return PETSC_SUCCESS;
        }

        if (decision.decision ==
            AdaptiveTimestepDecision3D::
                reject_minimum_timestep_reached) {
            completed.outcome =
                AdaptiveTimestepControllerOutcome3D::
                    minimum_timestep_reached;
            completed.retries =
                retries;
            report->emplace(
                std::move(completed));
            return PETSC_SUCCESS;
        }

        if (!decision.next_timestep_seconds) {
            return PETSC_ERR_PLIB;
        }
        timestep_seconds =
            *decision.next_timestep_seconds;
        ++retries;
    }
}

} // namespace mpmc::flow_discretization_petsc

#endif // MPMC_FLOW_DISCRETIZATION_PETSC_ADAPTIVE_TIMESTEP_CONTROLLER_HPP
