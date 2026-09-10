#include <mpmc/flash/sw92_asymmetric_max2.hpp>

bool sw92_asymmetric_max2_header() {
    mpmc::flash::Sw92AsymmetricMax2Options options;
    return options.selection.max_pair_attempts > 0;
}
