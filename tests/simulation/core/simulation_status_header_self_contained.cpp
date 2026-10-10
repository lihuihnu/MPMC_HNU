#include <mpmc/simulation/simulation_status.hpp>

int main() {
    return mpmc::simulation::simulation_status_is_terminal(
               mpmc::simulation::SimulationRunStatus::completed)
        ? 0
        : 1;
}
