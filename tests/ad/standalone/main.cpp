#include <mpmc/ad/dual.hpp>

// No test helper or platform header: exercise exactly what a library user needs.
int main() {
    using Number = mpmc::ad::Dual<double, 2>;
    const auto x = Number::variable(2.0, 0);
    const auto y = Number::variable(3.0, 1);
    const auto result = x * x + x * y;
    return result.value() == 10.0 && result.derivative(0) == 7.0 && result.derivative(1) == 2.0
               ? 0
               : 1;
}
