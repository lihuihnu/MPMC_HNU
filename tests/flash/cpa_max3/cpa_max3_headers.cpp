#include <mpmc/flash/cpa_three_phase.hpp>
#include <mpmc/flash/cpa_max3_phase_set.hpp>
#include <mpmc/flash/cpa_pt_flash_backend.hpp>

int main() {
    static_assert(mpmc::flash::cpa_pt_max3_convention[0] != '\0');
    static_assert(!mpmc::flash::PtThreePhaseResult::final_stability_checked);
    static_assert(!mpmc::flash::CpaPtMax3Result::global_stability_proven);
    return 0;
}
