// Include the public header first, twice, and in a separate translation unit.
// This checks self-containment, include guards, and linkage of inline operators.
#include <mpmc/ad/dual.hpp>
#include <mpmc/ad/dual.hpp>

mpmc::ad::Dual<double, 2> evaluate_in_other_translation_unit() {
    const auto x = mpmc::ad::Dual<double, 2>::variable(3.0, 0);
    return 2.0 * x * x + 1.0;
}
