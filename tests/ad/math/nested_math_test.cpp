#include <mpmc/ad/math.hpp>

#include <algorithm>
#include <array>
#include <cmath>
#include <numbers>
#include <stdexcept>
#include <string>
#include <string_view>

namespace ad = mpmc::ad;

namespace {
void require_nested_math(
    bool condition,
    std::string_view message) {
    if (!condition) {
        throw std::runtime_error(
            std::string{message});
    }
}

void near_nested_math(
    double actual,
    double expected) {
    const double scale =
        std::max(
            {1.0,
             std::abs(actual),
             std::abs(expected)});
    require_nested_math(
        std::isfinite(actual) &&
            std::abs(actual - expected) <=
                3.0e-12 * scale,
        "nested math value/derivative mismatch");
}
} // namespace

void test_nested_math() {
    using Inner = ad::Dual<double, 2>;
    using Outer = ad::Dual<Inner, 2>;

    const double xv = 1.2;
    const double yv = 0.7;
    const Outer x =
        Outer::variable(
            Inner::variable(xv, 0),
            0);
    const Outer y =
        Outer::variable(
            Inner::variable(yv, 1),
            1);

    using std::exp;
    using std::log;
    using std::sin;
    using std::sqrt;
    using std::pow;

    const Outer f =
        exp(x * y) +
        log(x) +
        sqrt(y) +
        sin(x) +
        pow(x, 1.5) +
        pow(2.0, y);

    const double e = std::exp(xv * yv);
    const double sqrt_x = std::sqrt(xv);
    const double sqrt_y = std::sqrt(yv);
    const double pow2y = std::pow(2.0, yv);
    const double ln2 = std::numbers::ln2_v<double>;

    const double expected_value =
        e + std::log(xv) + sqrt_y +
        std::sin(xv) +
        std::pow(xv, 1.5) + pow2y;
    const double gx =
        yv * e + 1.0 / xv +
        std::cos(xv) +
        1.5 * sqrt_x;
    const double gy =
        xv * e +
        1.0 / (2.0 * sqrt_y) +
        ln2 * pow2y;

    const double hxx =
        yv * yv * e -
        1.0 / (xv * xv) -
        std::sin(xv) +
        0.75 / sqrt_x;
    const double hxy =
        e * (1.0 + xv * yv);
    const double hyy =
        xv * xv * e -
        1.0 / (4.0 * yv * sqrt_y) +
        ln2 * ln2 * pow2y;

    near_nested_math(
        f.value().value(),
        expected_value);
    near_nested_math(
        f.value().derivative(0),
        gx);
    near_nested_math(
        f.value().derivative(1),
        gy);
    near_nested_math(
        f.derivative(0).derivative(0),
        hxx);
    near_nested_math(
        f.derivative(0).derivative(1),
        hxy);
    near_nested_math(
        f.derivative(1).derivative(0),
        hxy);
    near_nested_math(
        f.derivative(1).derivative(1),
        hyy);

    // Compile and evaluate every existing unary family on a nested scalar.
    const Outer positive =
        x + 1.0;
    const Outer unit =
        x / 4.0;
    const Outer above_one =
        x + 1.1;
    const Outer nonzero =
        x - 2.0;

    const std::array<Outer, 22> values{
        ad::exp(x),
        ad::exp2(x),
        ad::expm1(x),
        ad::log(positive),
        ad::log2(positive),
        ad::log10(positive),
        ad::log1p(x),
        ad::sqrt(positive),
        ad::cbrt(nonzero),
        ad::sin(x),
        ad::cos(x),
        ad::tan(unit),
        ad::asin(unit),
        ad::acos(unit),
        ad::atan(x),
        ad::sinh(x),
        ad::cosh(x),
        ad::tanh(x),
        ad::asinh(x),
        ad::acosh(above_one),
        ad::atanh(unit),
        ad::abs(nonzero)};

    for (const auto& value : values) {
        require_nested_math(
            std::isfinite(
                value.value().value()) &&
                std::isfinite(
                    value.value().derivative(0)) &&
                std::isfinite(
                    value.derivative(0).derivative(0)),
            "nested unary function produced non-finite test result");
    }

    const auto integer_power =
        ad::pow(x, 3);
    const auto variable_power =
        ad::pow(x, y);
    require_nested_math(
        std::isfinite(
            integer_power
                .derivative(0)
                .derivative(0)) &&
            std::isfinite(
                variable_power
                    .derivative(0)
                    .derivative(1)),
        "nested pow overload did not preserve second derivatives");

    bool domain_caught = false;
    try {
        const Outer negative{
            Inner::variable(-1.0, 0)};
        (void)ad::log(negative);
    } catch (const std::domain_error&) {
        domain_caught = true;
    }
    require_nested_math(
        domain_caught,
        "nested log did not apply domain check to deepest primal");
}
