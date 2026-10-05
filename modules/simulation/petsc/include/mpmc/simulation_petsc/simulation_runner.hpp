#ifndef MPMC_SIMULATION_PETSC_SIMULATION_RUNNER_HPP
#define MPMC_SIMULATION_PETSC_SIMULATION_RUNNER_HPP

#include <mpmc/flow_discretization_petsc/physical_timestep_driver.hpp>
#include <mpmc/simulation/simulation_cursor.hpp>
#include <mpmc/simulation/simulation_report.hpp>

#include <petscsys.h>

#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <memory>
#include <optional>
#include <stdexcept>
#include <string_view>
#include <utility>

namespace mpmc::simulation_petsc {

inline constexpr std::string_view simulation_runner_convention =
    "simulation_petsc/multi-timestep-runner/v1";

namespace simulation_runner_detail {

namespace fdp = mpmc::flow_discretization_petsc;
namespace sim = mpmc::simulation;

[[nodiscard]] inline bool same_time(double first, double second) noexcept {
    return sim::simulation_time_detail::same_time(first, second);
}

[[nodiscard]] inline bool same_positive_value(
    double first,
    double second) noexcept {
    return first > 0.0 && second > 0.0 && same_time(first, second);
}

[[nodiscard]] inline bool accepted_report_matches_clock(
    const fdp::PhysicalTimestepDriverReport3D& step_report,
    const fdp::AcceptedPhysicalTimeClock3D& clock,
    double time_before_seconds,
    std::size_t step_count_before,
    double effective_initial_cap_seconds) noexcept {
    if (!step_report.accepted() ||
        !step_report.accepted_record.has_value() ||
        !step_report.adaptive.accepted_timestep_seconds.has_value() ||
        !step_report.adaptive.next_timestep_seconds.has_value() ||
        step_report.adaptive.attempts.empty() ||
        step_report.adaptive.retries ==
            std::numeric_limits<std::size_t>::max() ||
        step_report.adaptive.attempts.size() !=
            step_report.adaptive.retries + 1U) {
        return false;
    }

    const auto& record = *step_report.accepted_record;
    const auto& last_attempt = step_report.adaptive.attempts.back();
    const bool accepted_decision =
        record.decision == fdp::AdaptiveTimestepDecision3D::accept_and_grow ||
        record.decision == fdp::AdaptiveTimestepDecision3D::accept_and_hold;

    return accepted_decision &&
        step_count_before != std::numeric_limits<std::size_t>::max() &&
        record.accepted_step_index == step_count_before &&
        clock.accepted_step_count() == step_count_before + 1U &&
        same_time(record.time_n_seconds, time_before_seconds) &&
        same_time(record.time_np1_seconds, clock.accepted_time_seconds()) &&
        same_positive_value(
            record.accepted_timestep_seconds,
            *step_report.adaptive.accepted_timestep_seconds) &&
        same_positive_value(
            record.next_timestep_seconds,
            *step_report.adaptive.next_timestep_seconds) &&
        same_positive_value(
            record.next_timestep_seconds,
            clock.next_timestep_seconds()) &&
        record.phase_transition_restarts ==
            last_attempt.result.phase_transition_restarts &&
        record.decision == last_attempt.decision &&
        std::isfinite(step_report.adaptive.initial_timestep_seconds) &&
        step_report.adaptive.initial_timestep_seconds > 0.0 &&
        (step_report.adaptive.initial_timestep_seconds <
             effective_initial_cap_seconds ||
         same_time(
             step_report.adaptive.initial_timestep_seconds,
             effective_initial_cap_seconds)) &&
        (record.accepted_timestep_seconds <
             effective_initial_cap_seconds ||
         same_time(
             record.accepted_timestep_seconds,
             effective_initial_cap_seconds));
}

[[nodiscard]] inline sim::AcceptedSimulationStepSummary project_accepted_step(
    const fdp::PhysicalTimestepDriverReport3D& step_report) {
    if (!step_report.accepted_record.has_value()) {
        throw std::logic_error(
            "mpmc::simulation_petsc: accepted physical timestep lacks accepted record");
    }
    const auto& record = *step_report.accepted_record;
    return {
        record.accepted_step_index,
        record.time_n_seconds,
        record.time_np1_seconds,
        record.accepted_timestep_seconds,
        record.next_timestep_seconds,
        step_report.adaptive.retries,
        step_report.adaptive.attempts.size(),
        record.phase_transition_restarts};
}

[[nodiscard]] inline bool terminal_rejection(
    const fdp::PhysicalTimestepDriverReport3D& step_report) noexcept {
    return !step_report.accepted() &&
        !step_report.accepted_record.has_value() &&
        (step_report.adaptive.outcome ==
             fdp::AdaptiveTimestepControllerOutcome3D::
                 retry_budget_exhausted ||
         step_report.adaptive.outcome ==
             fdp::AdaptiveTimestepControllerOutcome3D::
                 minimum_timestep_reached);
}

[[nodiscard]] inline bool clock_unchanged(
    const fdp::AcceptedPhysicalTimeClock3D& clock,
    double time_before_seconds,
    std::size_t step_count_before,
    double next_timestep_before_seconds) noexcept {
    return clock.accepted_step_count() == step_count_before &&
        same_time(clock.accepted_time_seconds(), time_before_seconds) &&
        same_positive_value(
            clock.next_timestep_seconds(),
            next_timestep_before_seconds);
}

[[nodiscard]] inline PetscErrorCode publish_terminal(
    sim::SimulationReport* lifecycle,
    sim::SimulationRunStatus status,
    std::optional<std::int64_t> lower_layer_error_code,
    std::optional<sim::SimulationReport>* report,
    PetscErrorCode return_code) {
    if (lifecycle == nullptr || report == nullptr) {
        return PETSC_ERR_ARG_NULL;
    }
    try {
        lifecycle->terminate(status, lower_layer_error_code);
        report->emplace(std::move(*lifecycle));
        return return_code;
    } catch (const std::exception&) {
        report->reset();
        return PETSC_ERR_PLIB;
    }
}

} // namespace simulation_runner_detail

struct SimulationRunnerOptions3D {
    mpmc::flow_discretization_petsc::PhysicalTimestepDriverOptions3D
        physical_timestep;
    std::size_t first_step_initial_phase_transition_restarts{};
};

/// Advance an already materialized accepted reservoir system until the supplied
/// SimulationCursor consumes its terminal hard boundary.
///
/// This function owns only lifecycle orchestration. Every physical step is
/// delegated to advance_one_physical_timestep_3d(), which remains the sole
/// owner of nonlinear solve, adaptive retry/cutback, phase-transition restart,
/// accepted-history rebase and AcceptedPhysicalTimeClock3D commit.
///
/// Lower-level timestep reports live only for one call, are validated against
/// the authoritative accepted clock, projected into the bounded simulation
/// report, and then discarded.
[[nodiscard]] inline PetscErrorCode advance_simulation_timeline_3d(
    MPI_Comm comm,
    std::unique_ptr<
        mpmc::flow_discretization_petsc::
            PhaseTransitionRebuiltNaturalVariableSystem3D>*
        accepted_system,
    mpmc::flow_discretization_petsc::
        PostSnesPhaseTransitionControllerBindings3D transition_bindings,
    const SimulationRunnerOptions3D& options,
    mpmc::flow_discretization_petsc::AcceptedPhysicalTimeClock3D* clock,
    mpmc::simulation::SimulationCursor* cursor,
    std::optional<mpmc::simulation::SimulationReport>* report) {
    namespace fdp = mpmc::flow_discretization_petsc;
    namespace sim = mpmc::simulation;
    using namespace simulation_runner_detail;

    if (accepted_system == nullptr || *accepted_system == nullptr ||
        clock == nullptr || cursor == nullptr || report == nullptr) {
        return PETSC_ERR_ARG_NULL;
    }
    report->reset();

    sim::SimulationReport lifecycle(
        {clock->accepted_time_seconds(),
         clock->accepted_step_count(),
         cursor->consumed_boundary_count()});

    const auto contract_failure = [&](PetscErrorCode code) {
        return publish_terminal(
            &lifecycle,
            sim::SimulationRunStatus::contract_violation,
            std::nullopt,
            report,
            code);
    };

    try {
        if (cursor->complete()) {
            if (!same_time(
                    clock->accepted_time_seconds(),
                    cursor->timeline().terminal_time_seconds())) {
                return contract_failure(PETSC_ERR_ARG_INCOMP);
            }
            return publish_terminal(
                &lifecycle,
                sim::SimulationRunStatus::completed,
                std::nullopt,
                report,
                PETSC_SUCCESS);
        }

        bool first_physical_step = true;
        while (!cursor->complete()) {
            const double time_before_seconds =
                clock->accepted_time_seconds();
            const std::size_t step_count_before =
                clock->accepted_step_count();
            const double next_timestep_before_seconds =
                clock->next_timestep_seconds();
            const auto next_boundary =
                cursor->next_hard_boundary();
            if (!next_boundary.has_value()) {
                return contract_failure(PETSC_ERR_PLIB);
            }

            const double effective_cap_seconds =
                cursor->cap_initial_timestep_seconds(
                    time_before_seconds,
                    options.physical_timestep
                        .initial_timestep_cap_seconds);

            auto step_options = options.physical_timestep;
            step_options.initial_timestep_cap_seconds =
                effective_cap_seconds;
            if (effective_cap_seconds <
                step_options.adaptive.minimum_timestep_seconds) {
                // A hard schedule clip is allowed to be smaller than the
                // normal adaptive floor. It becomes the retry floor for this
                // physical call, matching the existing Flow target-time loop.
                step_options.adaptive.minimum_timestep_seconds =
                    effective_cap_seconds;
            }

            lifecycle.record_physical_timestep_call();

            std::optional<fdp::PhysicalTimestepDriverReport3D>
                step_report;
            const PetscErrorCode error =
                fdp::advance_one_physical_timestep_3d(
                    comm,
                    accepted_system,
                    first_physical_step
                        ? options
                              .first_step_initial_phase_transition_restarts
                        : 0U,
                    transition_bindings,
                    step_options,
                    clock,
                    &step_report);
            first_physical_step = false;

            if (error != PETSC_SUCCESS) {
                if (!clock_unchanged(
                        *clock,
                        time_before_seconds,
                        step_count_before,
                        next_timestep_before_seconds)) {
                    return contract_failure(PETSC_ERR_PLIB);
                }
                return publish_terminal(
                    &lifecycle,
                    sim::SimulationRunStatus::execution_error,
                    static_cast<std::int64_t>(error),
                    report,
                    error);
            }

            if (!step_report.has_value()) {
                return contract_failure(PETSC_ERR_PLIB);
            }

            if (!step_report->accepted()) {
                if (!terminal_rejection(*step_report) ||
                    !clock_unchanged(
                        *clock,
                        time_before_seconds,
                        step_count_before,
                        next_timestep_before_seconds)) {
                    return contract_failure(PETSC_ERR_PLIB);
                }
                return publish_terminal(
                    &lifecycle,
                    sim::SimulationRunStatus::
                        physical_timestep_rejected,
                    std::nullopt,
                    report,
                    PETSC_SUCCESS);
            }

            if (!accepted_report_matches_clock(
                    *step_report,
                    *clock,
                    time_before_seconds,
                    step_count_before,
                    effective_cap_seconds)) {
                return contract_failure(PETSC_ERR_PLIB);
            }

            lifecycle.record_accepted_step(
                project_accepted_step(*step_report));

            const double time_after_seconds =
                clock->accepted_time_seconds();
            if (time_after_seconds > next_boundary->time_seconds &&
                !same_time(
                    time_after_seconds,
                    next_boundary->time_seconds)) {
                return contract_failure(PETSC_ERR_PLIB);
            }

            const auto reached =
                cursor->accept_time(time_after_seconds);
            if (reached.has_value()) {
                lifecycle.record_reached_boundary(
                    {*reached,
                     step_report->accepted_record
                         ->accepted_step_index});
            }
        }

        if (!same_time(
                clock->accepted_time_seconds(),
                cursor->timeline().terminal_time_seconds())) {
            return contract_failure(PETSC_ERR_PLIB);
        }

        return publish_terminal(
            &lifecycle,
            sim::SimulationRunStatus::completed,
            std::nullopt,
            report,
            PETSC_SUCCESS);
    } catch (const std::exception&) {
        if (lifecycle.terminal()) {
            report->reset();
            return PETSC_ERR_PLIB;
        }
        return contract_failure(PETSC_ERR_ARG_INCOMP);
    }
}

} // namespace mpmc::simulation_petsc

#endif // MPMC_SIMULATION_PETSC_SIMULATION_RUNNER_HPP
