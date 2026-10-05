#include <mpmc/simulation/simulation_timeline.hpp>

int main() {
    mpmc::simulation::SimulationTimeline timeline({1.0});
    return timeline.size() == 1U ? 0 : 1;
}
