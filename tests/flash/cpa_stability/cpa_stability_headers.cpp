#include <mpmc/flash/cpa_stability.hpp>

int main() {
    return mpmc::flash::cpa_pt_stability_convention[0] != '\0' ? 0 : 1;
}
