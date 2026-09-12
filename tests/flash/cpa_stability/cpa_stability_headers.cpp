#include <mpmc/flash/cpa_stability.hpp>

bool cpa_stability_public_header_self_contained() {
    return mpmc::flash::cpa_pt_stability_convention[0] != '\0';
}
