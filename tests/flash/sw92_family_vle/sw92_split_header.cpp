#include <mpmc/flash/sw92_split.hpp>
#include <mpmc/flash/sw92_split.hpp>

bool sw92_split_header() {
    return sizeof(mpmc::flash::Sw92FamilyPtSplitResult) > 0;
}
