#ifndef MPMC_SIMULATION_SIMULATION_STATUS_HPP
#define MPMC_SIMULATION_SIMULATION_STATUS_HPP

#include <string_view>

namespace mpmc::simulation {

inline constexpr std::string_view simulation_status_convention =
    "simulation/run-status/v1";

enum class SimulationRunStatus {
    running,
    completed,
    physical_timestep_rejected,
    execution_error,
    contract_violation
};

[[nodiscard]] constexpr bool simulation_status_is_terminal(
    SimulationRunStatus status) noexcept {
    return status != SimulationRunStatus::running;
}

[[nodiscard]] constexpr std::string_view simulation_status_name(
    SimulationRunStatus status) noexcept {
    switch (status) {
    case SimulationRunStatus::running:
        return "running";
    case SimulationRunStatus::completed:
        return "completed";
    case SimulationRunStatus::physical_timestep_rejected:
        return "physical_timestep_rejected";
    case SimulationRunStatus::execution_error:
        return "execution_error";
    case SimulationRunStatus::contract_violation:
        return "contract_violation";
    }
    return "unknown";
}

} // namespace mpmc::simulation

#endif // MPMC_SIMULATION_SIMULATION_STATUS_HPP
