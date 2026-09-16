#include <mpmc/thermodynamics/pr76_pure.hpp>
#include <mpmc/ad/differentiate.hpp>
#include <mpmc/ad/math.hpp>
#include <mpmc/ad/runtime_differentiate.hpp>

#include <algorithm>
#include <array>
#include <cmath>
#include <concepts>
#include <cstddef>
#include <exception>
#include <iostream>
#include <limits>
#include <optional>
#include <source_location>
#include <span>
#include <stdexcept>
#include <string>
#include <string_view>
#include <type_traits>
#include <utility>
#include <vector>

long double pr76_header_gas_constant();
double pr76_plain_evaluation(const mpmc::thermodynamics::Pr76Pure<double>& pure, double t);

namespace {
namespace th = mpmc::thermodynamics;
namespace ad = mpmc::ad;

void require(bool condition, std::string_view message,
             std::source_location where = std::source_location::current()) {
    if (!condition) {
        throw std::runtime_error(std::string(where.file_name()) + ":" +
                                 std::to_string(where.line()) + ": " + std::string(message));
    }
}

template <typename Error, typename Function>
void expect_error(Function&& function) {
    bool caught = false;
    try {
        function();
    } catch (const Error&) {
        caught = true;
    }
    require(caught, "expected error was not thrown");
}

template <std::floating_point T>
void near(T actual, long double expected,
          std::source_location where = std::source_location::current()) {
    // Relative tolerance only; exact expected zeros must stay zero. Long-double
    // comparison avoids underflowing the tolerance for small float/double data.
    const long double tolerance = 128.0L * std::numeric_limits<T>::epsilon();
    require(std::isfinite(actual) && std::isfinite(expected) &&
                std::abs(static_cast<long double>(actual) - expected) <=
                    tolerance * std::abs(expected),
            "PR76 value/derivative differs from independent reference", where);
}

// Artificial inputs for software/formula verification, NOT real fluid parameters.
th::Provenance source(std::string key) {
    return {th::SourceKind::synthetic_test, "test://thermodynamics/pr76", "fixture-v1",
            std::move(key), "Manufactured inputs; no physical accuracy claim",
            "Constructed in this test", "No third-party property dataset copied"};
}
th::SourcedScalar datum(double value, th::Unit unit) {
    return {value, unit, source("parameter"), "already in specified SI unit", "identity"};
}
struct Spec {
    std::string id;
    double tc = 400.0;
    double pc = 4e6;
    double omega = 0.2;
};

th::PrParameterSet parameters(const std::vector<Spec>& specs,
                             std::optional<th::ClosedInterval> interval = std::nullopt) {
    std::vector<th::Component> catalog;
    std::vector<std::string> order;
    th::PrParameterInput input;
    input.model_id = th::pr76_profile;
    input.dataset_id = "manufactured-pr76";
    input.revision = "v1";
    input.applicability = {interval, std::nullopt, source("bounds-unknown-unless-explicit")};
    for (const auto& spec : specs) {
        catalog.push_back({spec.id, "Artificial", th::ComponentKind::pure,
                           source(spec.id), std::nullopt});
        order.push_back(spec.id);
        input.pure.push_back({spec.id, datum(spec.tc, th::Unit::kelvin),
                              datum(spec.pc, th::Unit::pascal),
                              datum(spec.omega, th::Unit::dimensionless)});
    }
    for (std::size_t i = 0; i < order.size(); ++i) {
        for (std::size_t j = i + 1; j < order.size(); ++j) {
            input.binary.push_back({order[i], order[j], datum(0, th::Unit::dimensionless)});
        }
    }
    std::reverse(input.pure.begin(), input.pure.end()); // Exercise ID-based consumption.
    return th::PrParameterSet::create(catalog, order, input, th::DataPolicy::allow_synthetic_tests);
}

template <typename T>
th::Pr76Pure<T> kernel(double tc = 400, double pc = 4e6, double omega = 0.2,
                     std::optional<th::ClosedInterval> interval = std::nullopt) {
    const auto data = parameters({{"test:a", tc, pc, omega}}, interval);
    return th::Pr76Pure<T>::from_parameters(data, 0); // Factory must not retain data views.
}

struct Reference { long double a, b, da; };
Reference reference(long double temperature, const Spec& spec = {"test:a"}) {
    // Independently expanded alpha and its analytic temperature derivative.
    // Uses original input parameters and paper literals, never production helpers.
    const long double tc = spec.tc;
    const long double pc = spec.pc;
    const long double w = spec.omega;
    const long double r = 8.31446261815324L;
    const long double k = 0.37464L + 1.54226L * w - 0.26992L * w * w;
    const long double ac = 0.45724L * r * r * tc * tc / pc;
    const long double tr = temperature / tc;
    const long double a = ac * ((1+k)*(1+k) - 2*k*(1+k)*std::sqrt(tr) + k*k*tr);
    const long double da = ac * (k*k/tc - k*(1+k)/std::sqrt(temperature*tc));
    return {a, 0.07780L * r * tc / pc, da};
}

template <typename Number>
concept CanEvaluate = requires(const th::Pr76Pure<double>& pure, const Number& number) {
    pure.evaluate(number);
};
static_assert(CanEvaluate<double> && CanEvaluate<ad::Dual<double, 3>>);
static_assert(!CanEvaluate<int> && !CanEvaluate<float> && !CanEvaluate<ad::Dual<float, 3>>);
static_assert(!std::is_default_constructible_v<th::Pr76Pure<double>>);
static_assert(!std::is_copy_assignable_v<th::Pr76Pure<double>>);

template <typename T>
void printed_reference() {
    const auto pure = kernel<T>();
    // Independent Decimal(80-digit) evaluation of the paper equations. Tc=400,
    // Pc=4000000 and R are exact decimals; omega is the stored binary64 0.2,
    // imported with Decimal.from_float(0.2), NOT the ideal decimal 1/5.
    // Promoting a double parameter to long double cannot recover its lost bits.
    // References are not fluid data or outputs of the tested C++ implementation.
    near(pure.critical_attraction(), 1.26436532690287949346109988523210496L);
    near(pure.covolume(), 0.0000646865191692322072L);
    near(pure.kappa(), 0.672295200000000015923840024356650214L);
    near(pure.gas_constant(), 8.31446261815324L);
    const std::array<T, 3> temperatures{T{100}, T{400}, T{900}};
    const std::array<long double, 3> values{
        2.257259291573855898805178948320131479L,
        1.26436532690287949346109988523210496L,
        0.557205810927382359273418873306236173L};
    const std::array<long double, 3> slopes{
        -0.005678805945093580204611390443345527L,
        -0.002125066850808091924414700093767369L,
        -0.000940487152712929164349136643907983L};
    for (std::size_t i = 0; i < temperatures.size(); ++i) {
        const auto value = pure.evaluate(temperatures[i]);
        const auto differentiated = pure.evaluate(ad::Dual<T>::variable(temperatures[i], 0));
        near(value.a, values[i]);
        near(differentiated.a.value(), values[i]);
        near(differentiated.a.derivative(0), slopes[i]);
        require(value.b == pure.covolume(), "b varies with temperature");
        require(differentiated.b.derivative(0) == T{0}, "db/dT must be exactly zero");
    }
    // At the critical temperature alpha=1; the attraction slope is NOT zero.
    near(pure.evaluate(pure.critical_temperature_k()).a, pure.critical_attraction());
}

template <typename T>
void analytic_temperature() {
    for (double w : std::array{-0.1, 0.0, 0.2, 0.6}) {
        const auto pure = kernel<T>(400, 4e6, w);
        for (T temperature : std::array<T, 5>{T{100}, T{256}, T{400}, T{576}, T{900}}) {
            const auto expected = reference(temperature, {"test:a", 400, 4e6, w});
            const auto result = pure.evaluate(ad::Dual<T>::variable(temperature, 0));
            near(result.a.value(), expected.a);
            near(result.a.derivative(0), expected.da);
            near(result.b.value(), expected.b);
            require(result.b.derivative(0) == T{0}, "temperature-independent b lost");
        }
    }
    // omega>0.49 still uses PR76, NOT PR78. Reference uses Decimal.from_float(0.6).
    near(kernel<T>(400, 4e6, 0.6).kappa(), 1.202824799999999972947062332195855428L);
}

template <typename T>
void chain_rule() {
    const auto pure = kernel<T>();
    using Number = ad::Dual<T, 3>;
    const Number independent{T{10}, {T{1}, T{-2}, T{0}}};
    const Number temperature = T{100} + independent * independent;
    const auto result = pure.evaluate(temperature);
    const auto expected = reference(200);
    near(result.a.value(), expected.a);
    near(result.a.derivative(0), expected.da * 20);
    near(result.a.derivative(1), expected.da * -40);
    require(result.a.derivative(2) == T{0}, "inactive direction became active");
    require(result.b.derivatives() == typename Number::Gradient{}, "b is not an AD constant");
    const auto constant = pure.evaluate(Number{T{200}});
    require(constant.a.derivatives() == typename Number::Gradient{}, "zero seeds changed");
}

template <typename T>
void fixed_jacobian() {
    const auto pure = kernel<T>();
    const auto result = ad::value_and_jacobian([&](const auto& input) {
        const auto coefficients = pure.evaluate(input[0]);
        return std::array{coefficients.a, coefficients.b};
    }, std::array<T, 1>{T{400}});
    const auto expected = reference(400);
    near(result.values[0], expected.a);
    near(result.values[1], expected.b);
    near(result.jacobian[0][0], expected.da);
    require(result.jacobian[1][0] == T{0}, "fixed Jacobian b row");
}

template <typename T>
void runtime_jacobian() {
    ad::RuntimeJacobianWorkspace<T, 2> workspace;
    // Runtime size changes consume NEW ordered parameter snapshots, not padded components.
    for (std::size_t n : std::array<std::size_t, 4>{1, 5, 2, 4}) {
        std::vector<Spec> specs;
        std::vector<T> temperatures(n, T{300});
        for (std::size_t i = 0; i < n; ++i) {
            specs.push_back({"test:" + std::to_string(i), 350.0 + 20.0 * static_cast<double>(i),
                             4e6, 0.2});
        }
        std::reverse(specs.begin(), specs.end());
        const auto data = parameters(specs);
        std::vector<th::Pr76Pure<T>> pure;
        pure.reserve(n);
        for (std::size_t i = 0; i < n; ++i) {
            pure.push_back(th::Pr76Pure<T>::from_parameters(data, i));
            require(pure.back().source_record().component_id == specs[i].id, "ID binding changed");
        }
        using Number = ad::Dual<T, 2>;
        std::size_t calls = 0;
        const auto result = ad::value_and_jacobian_runtime<2>(
            [&](std::span<const Number> input, std::span<Number> output) {
                ++calls;
                for (std::size_t i = 0; i < n; ++i) {
                    const auto coefficients = pure[i].evaluate(input[i]);
                    output[2*i] = coefficients.a;
                    output[2*i+1] = coefficients.b;
                }
            }, temperatures, 2*n, workspace);
        require(calls == (n+1)/2, "incorrect block count");
        require(result.input_count == n && result.output_count == 2*n, "incorrect shape");
        for (std::size_t i = 0; i < n; ++i) {
            const auto expected = reference(300, specs[i]);
            near(result.values[2*i], expected.a);
            near(result.values[2*i+1], expected.b);
            for (std::size_t j = 0; j < n; ++j) {
                near(result.jacobian[2*i*n+j], i == j ? expected.da : 0.0L);
                require(result.jacobian[(2*i+1)*n+j] == T{0}, "runtime b row not zero");
            }
        }
    }
}

template <typename T>
void snapshot_lifetime() {
    const auto pure = kernel<T>(); // Source snapshot is already destroyed.
    const auto copy = pure;
    auto to_move = pure;
    const auto moved = std::move(to_move);
    near(copy.evaluate(T{400}).a, reference(400).a);
    near(moved.evaluate(T{400}).a, reference(400).a);
    require(moved.source_record().component_id == "test:a", "source identity lost");
    require(moved.source_record().critical_temperature->source.kind ==
                th::SourceKind::synthetic_test, "synthetic provenance erased");
    require(moved.dataset_id() == "manufactured-pr76" && moved.revision() == "v1", "revision lost");
    const auto data = parameters({{"test:a"}});
    expect_error<std::out_of_range>([&] { (void)th::Pr76Pure<T>::from_parameters(data, 1); });
    expect_error<std::out_of_range>([&] {
        (void)th::Pr76Pure<T>::from_parameters(data, std::numeric_limits<std::size_t>::max());
    });
}

template <typename T>
void temperature_domain() {
    const auto pure = kernel<T>();
    const T inf = std::numeric_limits<T>::infinity();
    const T nan = std::numeric_limits<T>::quiet_NaN();
    for (T invalid : std::array<T, 6>{T{0}, -T{0}, T{-1}, inf, -inf, nan}) {
        expect_error<std::domain_error>([&] { (void)pure.evaluate(invalid); });
        expect_error<std::domain_error>([&] { (void)pure.evaluate(ad::Dual<T>{invalid}); });
    }
    for (T seed : std::array<T, 2>{inf, nan}) {
        const ad::Dual<T> input{T{300}, {seed}};
        expect_error<std::domain_error>([&] { (void)pure.evaluate(input); });
        require(input.value() == T{300}, "invalid input changed on error");
    }
    near(pure.evaluate(T{400}).a, reference(400).a); // Errors do not poison the prepared object.
}

template <typename T>
void declared_bounds() {
    const auto bounded = kernel<T>(400, 4e6, 0.2, th::ClosedInterval{150, 700});
    for (T temperature : std::array<T, 2>{T{150}, T{700}}) {
        near(bounded.evaluate(temperature).a, reference(temperature).a);
    }
    const T below = std::nextafter(T{150}, T{0});
    const T above = std::nextafter(T{700}, std::numeric_limits<T>::infinity());
    near(bounded.evaluate(below).a, reference(below).a);
    const auto extrapolated = bounded.evaluate(ad::Dual<T>::variable(above, 0));
    near(extrapolated.a.value(), reference(above).a);
    near(extrapolated.a.derivative(0), reference(above).da);
    require(bounded.applicability().assess(static_cast<double>(below), 1e6) ==
                th::RangeAssessment::outside_declared_bounds,
            "lower extrapolation must stay visibly outside declared applicability");
    require(bounded.applicability().assess(static_cast<double>(above), 1e6) ==
                th::RangeAssessment::outside_declared_bounds,
            "upper extrapolation must stay visibly outside declared applicability");
    require(bounded.applicability().assess(400, 1e6) == th::RangeAssessment::unknown,
            "missing pressure applicability must remain unknown, not proven inside");
    const auto unknown = kernel<T>();
    near(unknown.evaluate(T{900}).a, reference(900).a);
    require(unknown.applicability().assess(900, 1e6) == th::RangeAssessment::unknown,
            "absent empirical bounds were promoted to scientific validity");
}

template <typename T>
void algebraic_extension() {
    // Tr=100 makes the bracket negative. Squaring it is explicitly an algebraic
    // continuation, NOT physical validation of this extreme temperature.
    const auto pure = kernel<T>();
    const auto result = pure.evaluate(ad::Dual<T>::variable(T{40000}, 0));
    const auto expected = reference(40000);
    require(result.a.derivative(0) > T{0}, "bracket was silently clipped/modified");
    near(result.a.value(), expected.a);
    near(result.a.derivative(0), expected.da);
}

template <typename T>
void numerical_range() {
    const auto large = kernel<T>(1, 1, 0.2);
    expect_error<std::range_error>([&] {
        (void)large.evaluate(std::numeric_limits<T>::max());
    });
    expect_error<std::range_error>([&] {
        const ad::Dual<T> input{static_cast<T>(0.001L), {std::numeric_limits<T>::max()}};
        (void)large.evaluate(input);
    });
    if constexpr (std::numeric_limits<T>::max_exponent <
                  std::numeric_limits<double>::max_exponent) {
        expect_error<std::range_error>([] { (void)kernel<T>(1e100, 1e6, 0.2); });
        expect_error<std::range_error>([] { (void)kernel<T>(1e-100, 1e6, 0.2); });
    }
    if constexpr (std::numeric_limits<T>::max_exponent <=
                  std::numeric_limits<double>::max_exponent) {
        expect_error<std::range_error>([] { (void)kernel<T>(1e200, 1, 0.2); });
        expect_error<std::range_error>([] { (void)kernel<T>(1e-200, 1e200, 0.2); });
    }
    // When pc==tc, a_c = 0.45724*R^2*tc and b = 0.07780*R. This independent
    // identity checks large/small scales without an overflowing reference square.
    constexpr bool narrow = std::numeric_limits<T>::max_exponent < 1024;
    const std::array<double, 2> scales = narrow ? std::array{1e20, 1e-30}
                                               : std::array{1e160, 1e-200};
    const long double factor = 0.45724L * 8.31446261815324L * 8.31446261815324L;
    for (double scale : scales) {
        const auto pure = kernel<T>(scale, scale, 0.2);
        const auto result = pure.evaluate(ad::Dual<T>::variable(static_cast<T>(scale), 0));
        near(result.a.value(), factor * scale);
        near(result.a.derivative(0), -factor * 0.672295200000000015923840024356650214L);
        near(result.b.value(), 0.07780L * 8.31446261815324L);
    }
}

template <typename T>
void run_case(std::string_view name) {
    if (name == "printed_reference") printed_reference<T>();
    else if (name == "analytic_temperature") analytic_temperature<T>();
    else if (name == "chain_rule") chain_rule<T>();
    else if (name == "fixed_jacobian") fixed_jacobian<T>();
    else if (name == "runtime_jacobian") runtime_jacobian<T>();
    else if (name == "snapshot_lifetime") snapshot_lifetime<T>();
    else if (name == "temperature_domain") temperature_domain<T>();
    else if (name == "declared_bounds") declared_bounds<T>();
    else if (name == "algebraic_extension") algebraic_extension<T>();
    else if (name == "numerical_range") numerical_range<T>();
    else throw std::invalid_argument("unknown PR76 test case");
}
} // namespace

int main(int argc, char** argv) {
    try {
        require(argc == 2, "provide one named PR76 test case");
        const std::string_view name{argv[1]};
        if (name == "headers") {
            require(pr76_header_gas_constant() == th::Pr76Pure<long double>::gas_constant(),
                    "new header is not self-contained/ODR safe");
            near(pr76_plain_evaluation(kernel<double>(), 400.0), reference(400).a);
        } else {
            run_case<float>(name);
            run_case<double>(name);
            run_case<long double>(name);
        }
        std::cout << "[PASS] " << name << '\n';
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "[FAIL] " << error.what() << '\n';
        return 1;
    }
}
