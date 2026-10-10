#ifndef MPMC_SIMULATION_SIMULATION_TIMELINE_HPP
#define MPMC_SIMULATION_SIMULATION_TIMELINE_HPP

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <stdexcept>
#include <string_view>
#include <utility>
#include <vector>

namespace mpmc::simulation {

inline constexpr std::string_view simulation_timeline_convention =
    "simulation/hard-timeline/v1";

namespace simulation_time_detail {

[[nodiscard]] inline bool same_time(double first, double second) noexcept {
    if (!std::isfinite(first) || !std::isfinite(second)) {
        return false;
    }
    const double scale = std::max({1.0, std::abs(first), std::abs(second)});
    return std::abs(first - second) <=
           64.0 * std::numeric_limits<double>::epsilon() * scale;
}

inline void require_accepted_time(double time_seconds) {
    if (!std::isfinite(time_seconds) || time_seconds < 0.0) {
        throw std::invalid_argument(
            "mpmc::simulation: accepted time must be finite and non-negative");
    }
}

} // namespace simulation_time_detail

struct SimulationHardBoundary {
    std::size_t index{};
    double time_seconds{};
};

class SimulationTimeline {
public:
    explicit SimulationTimeline(std::vector<double> hard_boundaries_seconds)
        : hard_boundaries_seconds_(std::move(hard_boundaries_seconds)) {
        if (hard_boundaries_seconds_.empty()) {
            throw std::invalid_argument(
                "mpmc::simulation: timeline requires at least one hard boundary");
        }
        for (std::size_t i = 0; i < hard_boundaries_seconds_.size(); ++i) {
            const double current = hard_boundaries_seconds_[i];
            simulation_time_detail::require_accepted_time(current);
            if (i == 0) {
                continue;
            }
            const double previous = hard_boundaries_seconds_[i - 1];
            if (!(current > previous) ||
                simulation_time_detail::same_time(current, previous)) {
                throw std::invalid_argument(
                    "mpmc::simulation: hard boundaries must be numerically distinct and strictly increasing");
            }
        }
    }

    [[nodiscard]] std::size_t size() const noexcept {
        return hard_boundaries_seconds_.size();
    }

    [[nodiscard]] double hard_boundary_seconds(std::size_t index) const {
        return hard_boundaries_seconds_.at(index);
    }

    [[nodiscard]] double terminal_time_seconds() const noexcept {
        return hard_boundaries_seconds_.back();
    }

private:
    std::vector<double> hard_boundaries_seconds_;
};

} // namespace mpmc::simulation

#endif // MPMC_SIMULATION_SIMULATION_TIMELINE_HPP
