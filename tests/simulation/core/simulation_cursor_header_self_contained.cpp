#include <mpmc/simulation/simulation_cursor.hpp>

int main() {
    mpmc::simulation::SimulationCursor cursor(
        mpmc::simulation::SimulationTimeline({1.0}), 0.0);
    return cursor.complete() ? 1 : 0;
}
