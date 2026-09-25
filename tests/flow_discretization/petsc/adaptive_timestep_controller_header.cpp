#include <mpmc/flow_discretization_petsc/adaptive_timestep_controller.hpp>

#include <string_view>

bool adaptive_timestep_controller_header() {
    return
        mpmc::flow_discretization_petsc::
            adaptive_timestep_controller_convention ==
        std::string_view{
            "flow_discretization_petsc/adaptive-timestep-controller/v1"};
}
