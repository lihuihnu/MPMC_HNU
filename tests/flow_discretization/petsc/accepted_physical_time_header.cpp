#include <mpmc/flow_discretization_petsc/accepted_physical_time.hpp>

#include <string_view>

bool accepted_physical_time_header() {
    return
        mpmc::flow_discretization_petsc::
            accepted_physical_time_convention ==
        std::string_view{
            "flow_discretization_petsc/accepted-physical-time/v1"};
}
