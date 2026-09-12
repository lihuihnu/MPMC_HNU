#include <mpmc/flash/cpa_split.hpp>
#include <mpmc/thermodynamics/cpa_parameters.hpp>

int main() {
    static_assert(mpmc::thermodynamics::cpa_profile[0] != '\0');
    static_assert(mpmc::flash::cpa_pt_vle_convention[0] != '\0');
    return 0;
}
