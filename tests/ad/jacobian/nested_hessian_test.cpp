#include <mpmc/ad/differentiate.hpp>
#include <mpmc/ad/math.hpp>

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <string>
#include <string_view>

namespace ad = mpmc::ad;

namespace {
void require_hessian(
    bool condition,
    std::string_view message) {
    if (!condition) {
        throw std::runtime_error(
            std::string{message});
    }
}

void near_hessian(
    double actual,
    double expected) {
    const double scale =
        std::max(
            {1.0,
             std::abs(actual),
             std::abs(expected)});
    require_hessian(
        std::isfinite(actual) &&
            std::abs(actual - expected) <=
                3.0e-12 * scale,
        "nested Hessian differs from analytic reference");
}
} // namespace

void test_nested_hessian_driver() {
    const std::array<double, 2> input{
        1.2, 0.7};
    int calls = 0;

    const auto result =
        ad::value_gradient_hessian(
            [&](const auto& p) {
                ++calls;
                using std::exp;
                using std::log;
                const auto& x = p[0];
                const auto& y = p[1];
                return x * x * y +
                    exp(x * y) +
                    log(x);
            },
            input);

    require_hessian(
        calls == 1,
        "Hessian callback was not evaluated exactly once");

    const double x = input[0];
    const double y = input[1];
    const double e =
        std::exp(x * y);

    const double value =
        x * x * y +
        e +
        std::log(x);
    const std::array<double, 2> gradient{
        2.0 * x * y +
            y * e +
            1.0 / x,
        x * x + x * e};
    const std::array<
        std::array<double, 2>,
        2> hessian{{
            {
                2.0 * y +
                    y * y * e -
                    1.0 / (x * x),
                2.0 * x +
                    e * (1.0 + x * y)},
            {
                2.0 * x +
                    e * (1.0 + x * y),
                x * x * e}}};

    near_hessian(result.value, value);
    for (std::size_t row = 0;
         row < 2U;
         ++row) {
        near_hessian(
            result.gradient[row],
            gradient[row]);
        for (std::size_t column = 0;
             column < 2U;
             ++column) {
            near_hessian(
                result.hessian[row][column],
                hessian[row][column]);
        }
    }
    near_hessian(
        result.hessian[0][1],
        result.hessian[1][0]);

    int invalid_calls = 0;
    bool caught = false;
    try {
        (void)ad::value_gradient_hessian(
            [&](const auto& p) {
                ++invalid_calls;
                return p[0];
            },
            std::array<double, 1>{
                std::numeric_limits<double>::
                    quiet_NaN()});
    } catch (const std::domain_error&) {
        caught = true;
    }
    require_hessian(
        caught && invalid_calls == 0,
        "non-finite Hessian input reached callback");
}
