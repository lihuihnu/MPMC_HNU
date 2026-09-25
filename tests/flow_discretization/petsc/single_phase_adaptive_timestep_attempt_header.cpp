#include <mpmc/flow_discretization_petsc/single_phase_adaptive_timestep_attempt.hpp>

#include <string_view>

bool single_phase_adaptive_timestep_attempt_header() {
    return
        mpmc::flow_discretization_petsc::
            single_phase_adaptive_timestep_attempt_convention ==
        std::string_view{
            "flow_discretization_petsc/single-phase-adaptive-timestep-attempt/v1"};
}
