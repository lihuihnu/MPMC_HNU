#include <mpmc/flow_discretization_petsc/physical_time_loop.hpp>

#include <string_view>

bool physical_time_loop_header() {
    return
        mpmc::flow_discretization_petsc::
            physical_time_loop_convention ==
        std::string_view{
            "flow_discretization_petsc/physical-time-loop/v1"};
}
