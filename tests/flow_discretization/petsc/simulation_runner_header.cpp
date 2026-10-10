#include <mpmc/simulation_petsc/simulation_runner.hpp>

#include <string_view>

bool simulation_runner_header() {
    return
        mpmc::simulation_petsc::simulation_runner_convention ==
        std::string_view{
            "simulation_petsc/multi-timestep-runner/v1"};
}
