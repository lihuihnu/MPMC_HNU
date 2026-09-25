#include <mpmc/flow_discretization_petsc/pr76_production_cell_evaluator.hpp>

#include <string_view>

bool pr76_production_cell_evaluator_header() {
    return
        mpmc::flow_discretization_petsc::
            pr76_production_cell_evaluator_convention ==
        std::string_view{
            "flow_discretization_petsc/pr76-production-cell-evaluator/v1"};
}
