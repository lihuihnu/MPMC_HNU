#include <mpmc/flow_discretization_petsc/single_phase_handoff_initial_system.hpp>

#include <string_view>

bool single_phase_handoff_initial_system_header() {
    return
        mpmc::flow_discretization_petsc::
            single_phase_handoff_initial_system_convention ==
        std::string_view{
            "flow_discretization_petsc/single-phase-handoff-initial-system/v1"};
}
