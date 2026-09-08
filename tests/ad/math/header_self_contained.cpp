#include <mpmc/ad/math.hpp>
#include <mpmc/ad/math.hpp> // Intentional: verify the include guard in a second translation unit.

mpmc::ad::Dual<double, 2> evaluate_math_in_other_translation_unit() {
    using Number = mpmc::ad::Dual<double, 2>;
    const auto x = Number::variable(4.0, 0);
    const auto y = Number::variable(9.0, 1);
    return mpmc::ad::exp(mpmc::ad::log(x)) + mpmc::ad::sqrt(y);
}
