#ifndef MPMC_SIMULATION_SIMULATION_CURSOR_HPP
#define MPMC_SIMULATION_SIMULATION_CURSOR_HPP

#include <mpmc/simulation/simulation_timeline.hpp>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <optional>
#include <stdexcept>
#include <string_view>
#include <utility>

namespace mpmc::simulation {

inline constexpr std::string_view simulation_cursor_convention =
    "simulation/accepted-hard-boundary-cursor/v1";

class SimulationCursor {
public:
    SimulationCursor(SimulationTimeline timeline, double accepted_time_seconds)
        : timeline_(std::move(timeline)) {
        simulation_time_detail::require_accepted_time(accepted_time_seconds);
        const double terminal = timeline_.terminal_time_seconds();
        if (accepted_time_seconds > terminal &&
            !simulation_time_detail::same_time(accepted_time_seconds, terminal)) {
            throw std::invalid_argument(
                "mpmc::simulation: accepted time exceeds terminal hard boundary");
        }
        while (next_boundary_index_ < timeline_.size()) {
            const double boundary =
                timeline_.hard_boundary_seconds(next_boundary_index_);
            if (accepted_time_seconds < boundary &&
                !simulation_time_detail::same_time(accepted_time_seconds, boundary)) {
                break;
            }
            ++next_boundary_index_;
        }
    }

    [[nodiscard]] const SimulationTimeline& timeline() const noexcept {
        return timeline_;
    }

    [[nodiscard]] bool complete() const noexcept {
        return next_boundary_index_ == timeline_.size();
    }

    [[nodiscard]] std::size_t consumed_boundary_count() const noexcept {
        return next_boundary_index_;
    }

    [[nodiscard]] std::optional<SimulationHardBoundary> next_hard_boundary() const {
        if (complete()) {
            return std::nullopt;
        }
        return SimulationHardBoundary{
            next_boundary_index_,
            timeline_.hard_boundary_seconds(next_boundary_index_)};
    }

    [[nodiscard]] double remaining_to_next_boundary_seconds(
        double accepted_time_seconds) const {
        simulation_time_detail::require_accepted_time(accepted_time_seconds);
        const auto next = next_hard_boundary();
        if (!next.has_value()) {
            throw std::logic_error(
                "mpmc::simulation: completed cursor has no remaining hard-boundary time");
        }
        if (accepted_time_seconds > next->time_seconds &&
            !simulation_time_detail::same_time(
                accepted_time_seconds, next->time_seconds)) {
            throw std::invalid_argument(
                "mpmc::simulation: accepted time skipped an unconsumed hard boundary");
        }
        if (simulation_time_detail::same_time(
                accepted_time_seconds, next->time_seconds)) {
            return 0.0;
        }
        return next->time_seconds - accepted_time_seconds;
    }

    [[nodiscard]] double cap_initial_timestep_seconds(
        double accepted_time_seconds,
        std::optional<double> caller_cap_seconds = std::nullopt) const {
        const double remaining =
            remaining_to_next_boundary_seconds(accepted_time_seconds);
        if (!(remaining > 0.0)) {
            throw std::logic_error(
                "mpmc::simulation: consume the accepted hard boundary before requesting another timestep cap");
        }
        if (!caller_cap_seconds.has_value()) {
            return remaining;
        }
        if (!std::isfinite(*caller_cap_seconds) ||
            !(*caller_cap_seconds > 0.0)) {
            throw std::invalid_argument(
                "mpmc::simulation: caller timestep cap must be finite and positive");
        }
        return std::min(remaining, *caller_cap_seconds);
    }

    [[nodiscard]] std::optional<SimulationHardBoundary>
    accept_time(double accepted_time_seconds) {
        simulation_time_detail::require_accepted_time(accepted_time_seconds);
        if (complete()) {
            const double terminal = timeline_.terminal_time_seconds();
            if (accepted_time_seconds > terminal &&
                !simulation_time_detail::same_time(
                    accepted_time_seconds, terminal)) {
                throw std::invalid_argument(
                    "mpmc::simulation: accepted time exceeds completed timeline");
            }
            return std::nullopt;
        }

        const SimulationHardBoundary next{
            next_boundary_index_,
            timeline_.hard_boundary_seconds(next_boundary_index_)};
        if (simulation_time_detail::same_time(
                accepted_time_seconds, next.time_seconds)) {
            ++next_boundary_index_;
            return next;
        }
        if (accepted_time_seconds > next.time_seconds) {
            throw std::invalid_argument(
                "mpmc::simulation: accepted time skipped an unconsumed hard boundary");
        }
        return std::nullopt;
    }

private:
    SimulationTimeline timeline_;
    std::size_t next_boundary_index_{};
};

} // namespace mpmc::simulation

#endif // MPMC_SIMULATION_SIMULATION_CURSOR_HPP
