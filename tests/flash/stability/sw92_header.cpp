#include <mpmc/flash/sw92_stability.hpp>
#include <mpmc/flash/sw92_stability.hpp>

bool stability_sw92_header() {
    return sizeof(mpmc::flash::Sw92FamilyStabilityResult) > 0;
}
