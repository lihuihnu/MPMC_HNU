#include <mpmc/flash/cpa_three_phase.hpp>
#include <mpmc/flash/cpa_max3_phase_set.hpp>
#include <mpmc/flash/cpa_pt_flash_backend.hpp>

int main() {
    static_assert(mpmc::flash::cpa_pt_max3_convention[0] != '\0');
    return 0;
}
