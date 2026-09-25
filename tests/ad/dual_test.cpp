#include <mpmc/ad/dual.hpp>

#include <array>
#include <cmath>
#include <concepts>
#include <cstddef>
#include <exception>
#include <iomanip>
#include <initializer_list>
#include <iostream>
#include <limits>
#include <source_location>
#include <stdexcept>
#include <string>
#include <string_view>
#include <type_traits>
#include <utility>

mpmc::ad::Dual<double, 2> evaluate_in_other_translation_unit();
void test_nested_arithmetic();

namespace {

using D = mpmc::ad::Dual<double, 2>;

template <typename A, typename B>
concept Addable = requires(A a, B b) { a + b; };

template <typename T, std::size_t N>
concept ValidDual = requires { typename mpmc::ad::Dual<T, N>; };

static_assert(ValidDual<float, 1> && ValidDual<double, 2> && ValidDual<long double, 4>);
static_assert(!ValidDual<int, 1> && !ValidDual<double, 0>);
static_assert(!ValidDual<const double, 1> && !ValidDual<volatile double, 1>);
static_assert(std::is_trivially_copyable_v<D>);
static_assert(!std::is_convertible_v<D, double> && !std::is_convertible_v<double, D>);
static_assert(!Addable<D, mpmc::ad::Dual<double, 3>>);
static_assert(!Addable<D, mpmc::ad::Dual<float, 2>>);
static_assert(std::is_same_v<decltype(std::declval<const D&>().derivatives()), const D::Gradient&>);
static_assert(std::is_same_v<decltype(std::declval<D>().derivatives()), D::Gradient>);
static_assert(noexcept(std::declval<D>() * std::declval<D>()));
static_assert(!noexcept(std::declval<D>() / std::declval<D>()));
constexpr auto compile_time_x = D::variable(2.0, 0);
constexpr auto compile_time_y =
    (compile_time_x * compile_time_x + 3.0 * compile_time_x - 1.0) / 2.0;
static_assert(compile_time_y.value() == 4.5);
static_assert(compile_time_y.derivative(0) == 3.5 && compile_time_y.derivative(1) == 0.0);

// These checks deliberately do not use assert: Release/NDEBUG must still test.
void require(bool condition, std::string_view message,
             const std::source_location location = std::source_location::current()) {
    if (!condition) {
        throw std::runtime_error(std::string(location.file_name()) + ":" +
                                 std::to_string(location.line()) + ": " + std::string(message));
    }
}

template <std::floating_point T>
void near(T actual, T expected, T absolute_tolerance = T{64} * std::numeric_limits<T>::epsilon(),
          const std::source_location location = std::source_location::current()) {
    // Dimensionless test functions: 64 eps covers a short chain's rounding error.
    // Extreme-scale regressions pass absolute_tolerance=0 to avoid accepting zero.
    const T relative_tolerance = T{64} * std::numeric_limits<T>::epsilon();
    const bool ok = std::isfinite(actual) && std::isfinite(expected) &&
                    std::abs(actual - expected) <=
                        absolute_tolerance + relative_tolerance * std::abs(expected);
    if (!ok) {
        std::cerr << std::setprecision(std::numeric_limits<T>::max_digits10)
                  << "actual=" << actual << ", expected=" << expected << '\n';
    }
    require(ok, "value/derivative outside tolerance", location);
}

template <typename Exception, typename Function>
void throws(Function function,
            const std::source_location location = std::source_location::current()) {
    try {
        function();
    } catch (const Exception&) {
        return;
    }
    require(false, "expected exception was not thrown", location);
}

template <std::floating_point T, std::size_t N>
void check(const mpmc::ad::Dual<T, N>& actual, T value, const std::array<T, N>& gradient,
           const std::source_location location = std::source_location::current()) {
    near(actual.value(), value, T{64} * std::numeric_limits<T>::epsilon(), location);
    for (std::size_t i = 0; i < N; ++i) {
        near(actual.derivative(i), gradient[i],
             T{64} * std::numeric_limits<T>::epsilon(), location);
    }
}

void construction_and_seeds() {
    check(D{}, 0.0, {0.0, 0.0});
    check(D{7.0}, 7.0, {0.0, 0.0});
    check(D::variable(2.0, 0), 2.0, {1.0, 0.0});
    check(D::variable(3.0, 1), 3.0, {0.0, 1.0});
    const D seeded{2.0, {4.0, -3.0}};
    auto copy = seeded;
    copy += D{1.0, {2.0, 5.0}};
    check(seeded, 2.0, {4.0, -3.0}); // Independent storage, not a shared tape/view.
    check(copy, 3.0, {6.0, 2.0});
    copy = D{9.0}; // Replacing with a constant must reset the old gradient.
    check(copy, 9.0, {0.0, 0.0});
    using Wide = mpmc::ad::Dual<double, 4>;
    check(Wide::variable(2.0, 3) * 3, 6.0, {0.0, 0.0, 0.0, 3.0});
    const auto owned_gradient = D::variable(5.0, 1).derivatives();
    require(owned_gradient[1] == 1.0, "temporary derivative storage must remain valid");
}

void arithmetic_and_scalar_paths() {
    const D a{6.0, {2.0, -1.0}};
    const D b{2.0, {3.0, 4.0}};
    check(+a, 6.0, {2.0, -1.0});
    check(-a, -6.0, {-2.0, 1.0});
    check(a + b, 8.0, {5.0, 3.0});
    check(a - b, 4.0, {-1.0, -5.0});
    check(a * b, 12.0, {22.0, 22.0});
    check(a / b, 3.0, {-3.5, -6.5});
    check(a / (-b), -3.0, {3.5, 6.5});
    check(a + 2, 8.0, {2.0, -1.0}); // Integer constants convert to the underlying T.
    check(2 + a, 8.0, {2.0, -1.0});
    check(a - 2.0, 4.0, {2.0, -1.0});
    check(2.0 - a, -4.0, {-2.0, 1.0});
    check(a * 2.0, 12.0, {4.0, -2.0});
    check(2.0 * a, 12.0, {4.0, -2.0});
    check(a / 2.0, 3.0, {1.0, -0.5});
    check(a / -2.0, -3.0, {-1.0, 0.5});
    check(12.0 / a, 2.0, {-2.0 / 3.0, 1.0 / 3.0});
    check(0.0 / a, 0.0, {0.0, 0.0});
    check(a * 0.0, 0.0, {0.0, 0.0});
    auto c = a;
    require(&(c += 2.0) == &c, "compound assignment must return *this");
    c -= 3.0;
    c *= 4.0;
    c /= 2.0;
    check(c, 10.0, {4.0, -2.0});
    check(a, 6.0, {2.0, -1.0});
}

void compound_aliasing() {
    const D initial{3.0, {2.0, -1.0}};
    auto a = initial;
    a += a;
    check(a, 6.0, {4.0, -2.0});
    a = initial;
    a -= a;
    check(a, 0.0, {0.0, 0.0});
    a = initial;
    const D& alias = a;
    a *= alias;
    check(a, 9.0, {12.0, -6.0});
    a = initial;
    a /= alias;
    check(a, 1.0, {0.0, 0.0});
    a = initial;
    a *= a + 1.0; // An owning temporary must not reference the modified lhs.
    check(a, 12.0, {14.0, -7.0});
}

void checked_indices() {
    throws<std::out_of_range>([] { (void)D::variable(1.0, D::derivative_count); });
    throws<std::out_of_range>([] {
        (void)D::variable(1.0, std::numeric_limits<std::size_t>::max());
    });
    const D a{};
    throws<std::out_of_range>([&] { (void)a.derivative(2); });
}

void zero_division_preserves_state() {
    for (const double zero : {0.0, -0.0}) {
        D a{3.0, {2.0, -1.0}};
        const D denominator{zero, {1.0, 4.0}};
        throws<std::domain_error>([&] { a /= denominator; });
        check(a, 3.0, {2.0, -1.0});
        throws<std::domain_error>([&] { a /= zero; });
        check(a, 3.0, {2.0, -1.0});
        throws<std::domain_error>([&] { (void)(a / denominator); });
        throws<std::domain_error>([&] { (void)(a / zero); });
        throws<std::domain_error>([&] { (void)(1.0 / denominator); });
        D self{zero, {2.0, -1.0}};
        throws<std::domain_error>([&] { self /= self; });
        check(self, zero, {2.0, -1.0});
    }
}

// One expression serves both ordinary floating point and AD numbers.
// f=(x*x+3*x*y-y)/(x+2), g=x*x*x-y*y.
template <typename Number>
std::array<Number, 2> residual(const Number& x, const Number& y) {
    return {(x * x + 3 * x * y - y) / (x + 2), x * x * x - y * y};
}

template <std::floating_point T>
void analytic_jacobian() {
    using Number = mpmc::ad::Dual<T, 2>;
    for (const T x : {T{-0.5}, T{0}, T{0.5}, T{2}}) {
        for (const T y : {T{-2}, T{0}, T{0.25}, T{3}}) {
            const auto result = residual(Number::variable(x, 0), Number::variable(y, 1));
            const auto ordinary = residual(x, y);
            // Closed-form derivatives, independently simplified by hand.
            const T denominator = x + T{2};
            check(result[0], ordinary[0],
                  {(x * x + T{4} * x + T{7} * y) / (denominator * denominator),
                   (T{3} * x - T{1}) / denominator});
            check(result[1], ordinary[1], {T{3} * x * x, -T{2} * y});
        }
    }
}

void directional_derivative() {
    using Number = mpmc::ad::Dual<double, 1>;
    const auto result = residual(Number{2.0, {2.0}}, Number{3.0, {-1.0}});
    // J*[2,-1] = [23/8,30], not merely one coordinate derivative.
    check(result[0], 19.0 / 4.0, {23.0 / 8.0});
    check(result[1], -1.0, {30.0});
}

void extreme_scale_division() {
    using Number = mpmc::ad::Dual<double, 1>;
    for (const double scale : {1e-200, 1e200}) {
        const Number a{2.0 * scale, {scale}};
        const Number b{scale, {2.0 * scale}};
        check(a / b, 2.0, {-3.0});
        check(a / scale, 2.0, {1.0});
        check(scale / b, 1.0, {-2.0});
        const auto small_derivative = Number{1.0, {1.0}} / scale;
        near(small_derivative.derivative(0), 1.0 / scale, 0.0);
    }
    if constexpr (std::numeric_limits<double>::is_iec559) {
        // A reciprocal would overflow, although each requested quotient is 1.
        const double tiny = 16.0 * std::numeric_limits<double>::denorm_min();
        check(Number{tiny, {tiny}} / tiny, 1.0, {1.0});
        check(Number{tiny, {tiny}} / Number{tiny}, 1.0, {1.0});
    }
}

void nonfinite_values_are_not_hidden() {
    if constexpr (std::numeric_limits<double>::is_iec559) {
        const double nan = std::numeric_limits<double>::quiet_NaN();
        const double infinity = std::numeric_limits<double>::infinity();
        require(std::isnan((D{nan} + 1.0).value()), "NaN primal was hidden");
        require(std::isinf((D{infinity} * 2.0).value()), "infinite primal was hidden");
        const auto result = D{2.0, {nan, 1.0}} + 1.0;
        require(result.value() == 3.0 && std::isnan(result.derivative(0)),
                "nonfinite derivative was hidden");
        const double largest = std::numeric_limits<double>::max();
        const auto overflow = D{largest, {largest, 0.0}} * 2.0;
        require(std::isinf(overflow.value()) && std::isinf(overflow.derivative(0)),
                "floating-point overflow was silently clipped");
    }
}

void separate_translation_unit() {
    check(evaluate_in_other_translation_unit(), 19.0, {12.0, 0.0});
}

struct TestCase {
    std::string_view name;
    void (*run)();
};

} // namespace

int main() {
    const std::array tests{
        TestCase{"construction_and_seeds", construction_and_seeds},
        TestCase{"arithmetic_and_scalar_paths", arithmetic_and_scalar_paths},
        TestCase{"compound_aliasing", compound_aliasing},
        TestCase{"checked_indices", checked_indices},
        TestCase{"zero_division_preserves_state", zero_division_preserves_state},
        TestCase{"analytic_jacobian_float", analytic_jacobian<float>},
        TestCase{"analytic_jacobian_double", analytic_jacobian<double>},
        TestCase{"analytic_jacobian_long_double", analytic_jacobian<long double>},
        TestCase{"directional_derivative", directional_derivative},
        TestCase{"extreme_scale_division", extreme_scale_division},
        TestCase{"nonfinite_values_are_not_hidden", nonfinite_values_are_not_hidden},
        TestCase{"nested_arithmetic", test_nested_arithmetic},
        TestCase{"separate_translation_unit", separate_translation_unit},
    };
    std::size_t failures = 0;
    for (const auto& test : tests) {
        try {
            test.run();
            std::cout << "[PASS] " << test.name << '\n';
        } catch (const std::exception& error) {
            ++failures;
            std::cerr << "[FAIL] " << test.name << ": " << error.what() << '\n';
        } catch (...) {
            ++failures;
            std::cerr << "[FAIL] " << test.name << ": unexpected non-standard exception\n";
        }
    }
    std::cout << tests.size() - failures << '/' << tests.size() << " cases passed\n";
    return failures == 0 ? 0 : 1;
}
