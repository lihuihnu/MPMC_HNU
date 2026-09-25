#include <mpmc/flow_discretization_petsc/physical_timestep_driver.hpp>

#include <string_view>

bool physical_timestep_driver_header() {
    return
        mpmc::flow_discretization_petsc::
            physical_timestep_driver_convention ==
        std::string_view{
            "flow_discretization_petsc/physical-timestep-driver/v1"};
}
