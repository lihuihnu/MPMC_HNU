#ifndef MPMC_SIMULATION_SIMULATION_REPORT_HPP
#define MPMC_SIMULATION_SIMULATION_REPORT_HPP

#include <mpmc/simulation/simulation_status.hpp>
#include <mpmc/simulation/simulation_timeline.hpp>

#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <optional>
#include <stdexcept>
#include <string_view>

namespace mpmc::simulation {

inline constexpr std::string_view simulation_report_convention =
    "simulation/bounded-run-report/v1";

struct AcceptedSimulationStepSummary {
    std::size_t accepted_step_index{};
    double time_n_seconds{};
    double time_np1_seconds{};
    double accepted_timestep_seconds{};
    double next_timestep_seconds{};
    std::size_t retry_count{};
    std::size_t attempt_count{1U};
    std::size_t phase_transition_restarts{};
};

struct ReachedSimulationBoundary {
    SimulationHardBoundary boundary;
    std::size_t accepted_step_index{};
};

struct SimulationReportAnchor {
    double accepted_time_seconds{};
    std::size_t accepted_step_count{};
    std::size_t consumed_boundary_count{};
};

struct SimulationLifecycleDiagnostics {
    std::size_t physical_timestep_calls{};
    std::size_t accepted_steps{};
    std::size_t reached_boundaries{};
    std::size_t accepted_retries{};
    std::size_t accepted_phase_transition_restarts{};
};

namespace simulation_report_detail {

inline void require_running(SimulationRunStatus status) {
    if (status != SimulationRunStatus::running) {
        throw std::logic_error(
            "mpmc::simulation: completed report cannot accept more lifecycle records");
    }
}

[[nodiscard]] inline std::size_t checked_add(
    std::size_t first,
    std::size_t second,
    const char* message) {
    if (second > std::numeric_limits<std::size_t>::max() - first) {
        throw std::overflow_error(message);
    }
    return first + second;
}

inline void validate_accepted_step(
    const AcceptedSimulationStepSummary& step) {
    simulation_time_detail::require_accepted_time(step.time_n_seconds);
    simulation_time_detail::require_accepted_time(step.time_np1_seconds);
    if (!(step.time_np1_seconds > step.time_n_seconds) ||
        simulation_time_detail::same_time(
            step.time_np1_seconds, step.time_n_seconds) ||
        !std::isfinite(step.accepted_timestep_seconds) ||
        !(step.accepted_timestep_seconds > 0.0) ||
        !std::isfinite(step.next_timestep_seconds) ||
        !(step.next_timestep_seconds > 0.0)) {
        throw std::invalid_argument(
            "mpmc::simulation: invalid accepted physical-step summary");
    }
    const double reconstructed =
        step.time_n_seconds + step.accepted_timestep_seconds;
    if (!std::isfinite(reconstructed) ||
        !simulation_time_detail::same_time(
            reconstructed, step.time_np1_seconds)) {
        throw std::invalid_argument(
            "mpmc::simulation: accepted timestep does not match accepted time interval");
    }
    if (step.attempt_count == 0U ||
        step.retry_count == std::numeric_limits<std::size_t>::max() ||
        step.attempt_count != step.retry_count + 1U) {
        throw std::invalid_argument(
            "mpmc::simulation: accepted step attempt/retry summary is inconsistent");
    }
}

} // namespace simulation_report_detail

class SimulationReport {
public:
    explicit SimulationReport(SimulationReportAnchor anchor)
        : anchor_(anchor) {
        simulation_time_detail::require_accepted_time(
            anchor_.accepted_time_seconds);
    }

    [[nodiscard]] const SimulationReportAnchor& anchor() const noexcept {
        return anchor_;
    }

    [[nodiscard]] SimulationRunStatus status() const noexcept {
        return status_;
    }

    [[nodiscard]] bool terminal() const noexcept {
        return simulation_status_is_terminal(status_);
    }

    [[nodiscard]] const SimulationLifecycleDiagnostics& diagnostics() const noexcept {
        return diagnostics_;
    }

    [[nodiscard]] const std::optional<AcceptedSimulationStepSummary>&
    last_accepted_step() const noexcept {
        return last_accepted_step_;
    }

    [[nodiscard]] const std::optional<ReachedSimulationBoundary>&
    last_reached_boundary() const noexcept {
        return last_reached_boundary_;
    }

    [[nodiscard]] const std::optional<std::int64_t>&
    lower_layer_error_code() const noexcept {
        return lower_layer_error_code_;
    }

    void record_physical_timestep_call() {
        simulation_report_detail::require_running(status_);
        diagnostics_.physical_timestep_calls =
            simulation_report_detail::checked_add(
                diagnostics_.physical_timestep_calls,
                1U,
                "mpmc::simulation: physical timestep call counter overflow");
    }

    void record_accepted_step(const AcceptedSimulationStepSummary& step) {
        simulation_report_detail::require_running(status_);
        simulation_report_detail::validate_accepted_step(step);
        const std::size_t expected_call_count =
            simulation_report_detail::checked_add(
                diagnostics_.accepted_steps,
                1U,
                "mpmc::simulation: physical timestep call expectation overflow");
        if (diagnostics_.physical_timestep_calls != expected_call_count) {
            throw std::logic_error(
                "mpmc::simulation: accepted step must correspond to exactly one new physical-timestep call");
        }
        const std::size_t expected_index =
            simulation_report_detail::checked_add(
                anchor_.accepted_step_count,
                diagnostics_.accepted_steps,
                "mpmc::simulation: accepted step index overflow");
        if (step.accepted_step_index != expected_index) {
            throw std::invalid_argument(
                "mpmc::simulation: accepted step index is not sequential");
        }
        const double expected_time = last_accepted_step_.has_value()
            ? last_accepted_step_->time_np1_seconds
            : anchor_.accepted_time_seconds;
        if (!simulation_time_detail::same_time(
                step.time_n_seconds, expected_time)) {
            throw std::invalid_argument(
                "mpmc::simulation: accepted step does not continue the accepted time history");
        }

        diagnostics_.accepted_steps =
            simulation_report_detail::checked_add(
                diagnostics_.accepted_steps,
                1U,
                "mpmc::simulation: accepted step counter overflow");
        diagnostics_.accepted_retries =
            simulation_report_detail::checked_add(
                diagnostics_.accepted_retries,
                step.retry_count,
                "mpmc::simulation: accepted retry counter overflow");
        diagnostics_.accepted_phase_transition_restarts =
            simulation_report_detail::checked_add(
                diagnostics_.accepted_phase_transition_restarts,
                step.phase_transition_restarts,
                "mpmc::simulation: phase-transition restart counter overflow");
        last_accepted_step_ = step;
    }

    void record_reached_boundary(const ReachedSimulationBoundary& reached) {
        simulation_report_detail::require_running(status_);
        if (!last_accepted_step_.has_value()) {
            throw std::logic_error(
                "mpmc::simulation: hard boundary cannot be reached without an accepted step");
        }
        simulation_time_detail::require_accepted_time(
            reached.boundary.time_seconds);
        const std::size_t expected_boundary_index =
            simulation_report_detail::checked_add(
                anchor_.consumed_boundary_count,
                diagnostics_.reached_boundaries,
                "mpmc::simulation: hard-boundary index overflow");
        if (reached.boundary.index != expected_boundary_index ||
            reached.accepted_step_index !=
                last_accepted_step_->accepted_step_index ||
            !simulation_time_detail::same_time(
                reached.boundary.time_seconds,
                last_accepted_step_->time_np1_seconds)) {
            throw std::invalid_argument(
                "mpmc::simulation: reached hard boundary is inconsistent with accepted lifecycle state");
        }
        if (last_reached_boundary_.has_value() &&
            (!(reached.boundary.time_seconds >
               last_reached_boundary_->boundary.time_seconds) ||
             simulation_time_detail::same_time(
                 reached.boundary.time_seconds,
                 last_reached_boundary_->boundary.time_seconds))) {
            throw std::invalid_argument(
                "mpmc::simulation: reached hard boundaries are not strictly increasing");
        }
        diagnostics_.reached_boundaries =
            simulation_report_detail::checked_add(
                diagnostics_.reached_boundaries,
                1U,
                "mpmc::simulation: reached-boundary counter overflow");
        last_reached_boundary_ = reached;
    }

    void terminate(
        SimulationRunStatus terminal_status,
        std::optional<std::int64_t> lower_layer_error_code = std::nullopt) {
        simulation_report_detail::require_running(status_);
        if (!simulation_status_is_terminal(terminal_status)) {
            throw std::invalid_argument(
                "mpmc::simulation: terminal status cannot be running");
        }
        if (terminal_status == SimulationRunStatus::execution_error) {
            if (!lower_layer_error_code.has_value() ||
                *lower_layer_error_code == 0) {
                throw std::invalid_argument(
                    "mpmc::simulation: execution error requires a nonzero lower-layer code");
            }
        } else if (lower_layer_error_code.has_value()) {
            throw std::invalid_argument(
                "mpmc::simulation: lower-layer error code is valid only for execution_error");
        }
        status_ = terminal_status;
        lower_layer_error_code_ = lower_layer_error_code;
    }

private:
    SimulationReportAnchor anchor_;
    SimulationRunStatus status_{SimulationRunStatus::running};
    SimulationLifecycleDiagnostics diagnostics_;
    std::optional<AcceptedSimulationStepSummary> last_accepted_step_;
    std::optional<ReachedSimulationBoundary> last_reached_boundary_;
    std::optional<std::int64_t> lower_layer_error_code_;
};

} // namespace mpmc::simulation

#endif // MPMC_SIMULATION_SIMULATION_REPORT_HPP
