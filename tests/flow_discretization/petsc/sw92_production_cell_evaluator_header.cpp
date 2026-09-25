#include <mpmc/flow_discretization_petsc/sw92_production_cell_evaluator.hpp>

bool sw92_production_cell_evaluator_header() {
    return mpmc::flow_discretization_petsc::
               sw92_production_cell_evaluator_convention ==
           "flow_discretization_petsc/sw92-production-cell-evaluator/zero-salinity-co2-water/v1";
}
