#include <mpmc/flash/sw92_dual_model.hpp>

bool sw92_dual_model_header() {
    return mpmc::flash::sw92_whitson_dual_model_algorithm ==
           "SW92-equilibrium/whitson-dual-model-observables/v1";
}
