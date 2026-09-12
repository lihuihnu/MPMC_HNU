#include <mpmc/flash/cpa_split.hpp>

bool cpa_split_public_header_self_contained() {
    return mpmc::flash::cpa_pt_vle_convention[0] != '\0';
}
