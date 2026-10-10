#include <mpmc/simulation/simulation_report.hpp>

int main() {
    mpmc::simulation::SimulationReport report({0.0, 0U, 0U});
    return report.terminal() ? 1 : 0;
}
