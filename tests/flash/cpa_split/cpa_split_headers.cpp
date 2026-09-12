#include <mpmc/flash/cpa_split.hpp>

int main() {
    static_assert(mpmc::flash::cpa_pt_vle_convention.size() > 0U);
    return 0;
}
