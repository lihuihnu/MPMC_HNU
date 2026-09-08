#include <mpmc/ad/differentiate.hpp>
#include <mpmc/ad/differentiate.hpp>

// The public header is first and repeated. Link the same template specialization
// from two translation units without relying on math.hpp or test helpers.
mpmc::ad::ValueAndJacobian<double, 1, 1> jacobian_from_separate_translation_unit() {
    return mpmc::ad::value_and_jacobian(
        [](const auto& x) { return x; }, std::array<double, 1>{3.0});
}
