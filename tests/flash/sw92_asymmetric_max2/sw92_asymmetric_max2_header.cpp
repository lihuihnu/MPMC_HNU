#include <mpmc/flash/sw92_asymmetric_max2.hpp>

bool sw92_asymmetric_max2_header() {
    mpmc::flash::Sw92AsymmetricMax2Options options;
    return options.selection.max_pair_attempts > 0 &&
           mpmc::flash::detail::sw92_max2_effective_tpd_tolerance(0.0, 0.0) == 0.0 &&
           mpmc::flash::detail::sw92_max2_effective_tpd_tolerance(0.0, 0.25) == 0.25;
}
