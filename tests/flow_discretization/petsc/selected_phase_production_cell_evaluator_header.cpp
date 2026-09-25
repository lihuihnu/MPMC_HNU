#include <mpmc/flow_discretization_petsc/selected_phase_production_cell_evaluator.hpp>

bool selected_phase_production_cell_evaluator_header() {
    return mpmc::flow_discretization_petsc::
               selected_phase_production_cell_evaluator_convention ==
           "flow_discretization_petsc/selected-phase-production-cell-evaluator/v1";
}
