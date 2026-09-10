#include <mpmc/flash/sw92_asymmetric_orchestration.hpp>

bool sw92_asymmetric_orchestration_header() {
    mpmc::flash::Sw92AsymmetricPairSelectionOptions options;
    return options.max_pair_attempts > 0;
}
