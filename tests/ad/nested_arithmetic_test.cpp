#include <mpmc/ad/dual.hpp>

#include <algorithm>
#include <array>
#include <cmath>
#include <stdexcept>
#include <string>
#include <string_view>
#include <type_traits>

namespace ad = mpmc::ad;

namespace {
void require_nested(bool condition, std::string_view message) {
    if (!condition) {
        throw std::runtime_error(std::string{message});
    }
}

void near_nested(double actual, double expected) {
    const double scale = std::max({1.0, std::abs(actual), std::abs(expected)});
    require_nested(
        std::isfinite(actual) &&
            std::abs(actual - expected) <= 2.0e-13 * scale,
        "nested arithmetic value/derivative mismatch");
}
} // namespace

void test_nested_arithmetic() {
    using Inner = ad::Dual<double, 2>;
    using Outer = ad::Dual<Inner, 2>;

    static_assert(std::same_as<Outer::Scalar, Inner>);
    static_assert(std::same_as<Outer::BaseScalar, double>);
    static_assert(std::is_trivially_copyable_v<Outer>);

    const Outer x =
        Outer::variable(
            Inner::variable(2.0, 0),
            0);
    const Outer y =
        Outer::variable(
            Inner::variable(3.0, 1),
            1);

    const Outer f =
        x * x * x + x * y + 2.0 / y;

    near_nested(f.value().value(), 44.0 / 3.0);
    near_nested(f.value().derivative(0), 15.0);
    near_nested(f.value().derivative(1), 16.0 / 9.0);

    near_nested(f.derivative(0).value(), 15.0);
    near_nested(f.derivative(1).value(), 16.0 / 9.0);
    near_nested(f.derivative(0).derivative(0), 12.0);
    near_nested(f.derivative(0).derivative(1), 1.0);
    near_nested(f.derivative(1).derivative(0), 1.0);
    near_nested(f.derivative(1).derivative(1), 4.0 / 27.0);

    const Inner shift =
        Inner::variable(0.5, 0);
    const Outer shifted = x + shift;
    near_nested(shifted.value().value(), 2.5);
    near_nested(shifted.value().derivative(0), 2.0);
    near_nested(shifted.derivative(0).value(), 1.0);

    Outer state = x;
    const Outer original = state;
    const Inner zero_with_seed{
        0.0,
        {1.0, 0.0}};
    bool caught = false;
    try {
        state /= zero_with_seed;
    } catch (const std::domain_error&) {
        caught = true;
    }
    require_nested(
        caught,
        "nested denominator with zero deepest primal was accepted");
    near_nested(
        state.value().value(),
        original.value().value());
    near_nested(
        state.derivative(0).value(),
        original.derivative(0).value());
}
