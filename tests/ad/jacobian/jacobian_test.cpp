#include <mpmc/ad/differentiate.hpp>
#include <mpmc/ad/math.hpp>

#include <algorithm>
#include <array>
#include <cmath>
#include <concepts>
#include <cstddef>
#include <exception>
#include <functional>
#include <iostream>
#include <limits>
#include <source_location>
#include <stdexcept>
#include <string>
#include <string_view>
#include <type_traits>
#include <utility>

mpmc::ad::ValueAndJacobian<double, 1, 1> jacobian_from_separate_translation_unit();
void test_nested_hessian_driver();

namespace {
namespace ad = mpmc::ad;

void require(bool condition, std::string_view message,
             std::source_location where = std::source_location::current()) {
    if (!condition) {
        throw std::runtime_error(std::string(where.file_name()) + ":" +
                                 std::to_string(where.line()) + ": " + std::string(message));
    }
}

template <std::floating_point T>
void near(T actual, T expected,
          std::source_location where = std::source_location::current()) {
    // Dimensionless short expressions: absolute AND relative scale, not a
    // tolerance selected from the AD result. Expected derivatives are analytic.
    const T tolerance = T{64} * std::numeric_limits<T>::epsilon();
    require(std::isfinite(actual) && std::isfinite(expected) &&
                std::abs(actual - expected) <= tolerance * std::max(T{1}, std::abs(expected)),
            "value/Jacobian differs from the independent analytic reference", where);
}

template <typename Function, typename Inputs>
concept CanEvaluate = requires(Function&& function, const Inputs& inputs) {
    ad::value_and_jacobian(std::forward<Function>(function), inputs);
};

// Auto return deduction instantiates the wrapper body during callability checks.
// Define probes even though no test should execute them.
using Number2 = ad::Dual<double, 2>;
using Seed2 = std::array<Number2, 2>;
template <typename Output>
struct Returns {
    Output operator()(const Seed2&) const {
        throw std::logic_error("signature probe must not be executed");
    }
};
struct MutableInput {
    Seed2 operator()(Seed2& input) const { return input; }
};
struct Identity {
    template <typename Number, std::size_t N>
    auto operator()(const std::array<Number, N>& input) const { return input; }
};

static_assert(CanEvaluate<Returns<std::array<Number2, 3>>, std::array<double, 2>>);
static_assert(!CanEvaluate<Returns<std::array<double, 3>>, std::array<double, 2>>);
static_assert(!CanEvaluate<Returns<Number2>, std::array<double, 2>>);
static_assert(!CanEvaluate<Returns<std::array<ad::Dual<float, 2>, 3>>, std::array<double, 2>>);
static_assert(!CanEvaluate<Returns<std::array<ad::Dual<double, 3>, 3>>, std::array<double, 2>>);
static_assert(!CanEvaluate<Returns<std::array<Number2, 0>>, std::array<double, 2>>);
static_assert(!CanEvaluate<Returns<const Seed2&>, std::array<double, 2>>);
static_assert(!CanEvaluate<Returns<Seed2&&>, std::array<double, 2>>);
static_assert(!CanEvaluate<MutableInput, std::array<double, 2>>);
static_assert(!CanEvaluate<Identity, std::array<double, 0>>);
static_assert(!CanEvaluate<Identity, std::array<int, 2>>);
static_assert(!CanEvaluate<Identity, std::array<Number2, 2>>);
static_assert(!CanEvaluate<int, std::array<double, 2>>);

template <typename T>
void test_identity_seeds() {
    using Number = ad::Dual<T, 4>;
    const std::array<T, 4> input{T{-2}, T{0}, T{3}, T{5}};
    int calls = 0;
    const auto result = ad::value_and_jacobian([&](const auto& x) {
        static_assert(std::same_as<decltype(x), const std::array<Number, 4>&>);
        ++calls;
        for (std::size_t i = 0; i < 4; ++i) {
            require(x[i].value() == input[i], "input order changed");
            for (std::size_t j = 0; j < 4; ++j) {
                require(x[i].derivative(j) == (i == j ? T{1} : T{0}), "incorrect identity seed");
            }
        }
        return x;
    }, input);
    static_assert(std::same_as<decltype(result), const ad::ValueAndJacobian<T, 4, 4>>);
    require(calls == 1 && result.values == input, "identity must be evaluated exactly once");
    for (std::size_t i = 0; i < 4; ++i) {
        for (std::size_t j = 0; j < 4; ++j) {
            require(result.jacobian[i][j] == (i == j ? T{1} : T{0}), "identity extraction");
        }
    }
}

// The same user expression supports plain floating point and AD via ADL.
// These 3 outputs and 2 inputs deliberately form a NON-square Jacobian.
template <typename Number>
auto vector_function(const std::array<Number, 2>& p) {
    using std::exp;
    using std::log;
    const auto& x = p[0];
    const auto& y = p[1];
    return std::array{x * x + y, x * y, exp(x) + log(y)};
}

template <typename T>
void test_rectangular_analytic() {
    for (T x : std::array<T, 3>{T{-1}, T{0}, T{2}}) {
        for (T y : std::array<T, 3>{T{1}, T{2}, T{3}}) {
            const std::array<T, 2> input{x, y};
            const auto result = ad::value_and_jacobian(
                [](const auto& p) { return vector_function(p); }, input);
            static_assert(decltype(result)::input_count == 2);
            static_assert(decltype(result)::output_count == 3);
            const auto ordinary_values = vector_function(input);
            // Independently differentiated rows: [2x,1], [y,x], [exp(x),1/y].
            const std::array<std::array<T, 2>, 3> expected{{
                {T{2} * x, T{1}}, {y, x}, {std::exp(x), T{1} / y}}};
            for (std::size_t i = 0; i < 3; ++i) {
                near(result.values[i], ordinary_values[i]);
                for (std::size_t j = 0; j < 2; ++j) {
                    near(result.jacobian[i][j], expected[i][j]);
                }
            }
        }
    }
}

template <typename T>
void test_constants_unused_inputs() {
    using Number = ad::Dual<T, 3>;
    const auto result = ad::value_and_jacobian([](const auto& p) {
        return std::array{Number{T{7}}, p[0] * T{2}};
    }, std::array<T, 3>{T{2}, T{3}, T{4}});
    require(result.values == std::array<T, 2>{T{7}, T{4}}, "constant output values");
    require(result.jacobian[0] == std::array<T, 3>{}, "constant output must have a zero row");
    require(result.jacobian[1] == std::array<T, 3>{T{2}, T{0}, T{0}}, "unused input columns");
}

template <typename T>
void test_single_input_output() {
    const auto result = ad::value_and_jacobian([](const auto& p) {
        return std::array{p[0] * p[0] * p[0]};
    }, std::array<T, 1>{T{-2}});
    require(result.values[0] == T{-8} && result.jacobian[0][0] == T{12}, "1x1 Jacobian");
}

template <typename T>
struct NoncopyableFunction {
    explicit NoncopyableFunction(int& calls) : calls_(calls) {}
    NoncopyableFunction(const NoncopyableFunction&) = delete;
    NoncopyableFunction& operator=(const NoncopyableFunction&) = delete;

    auto operator()(const std::array<ad::Dual<T>, 1>& p) & {
        ++calls_;
        return std::array{p[0] * T{2}};
    }
    auto operator()(const std::array<ad::Dual<T>, 1>& p) && {
        ++calls_;
        return std::array{p[0] * T{3}};
    }
private:
    int& calls_;
};

template <typename T>
auto free_function(const std::array<ad::Dual<T>, 1>& p) { return p; }

template <typename T>
void test_callable_forwarding() {
    const std::array<T, 1> input{T{2}};
    int calls = 0;
    NoncopyableFunction<T> function{calls};
    const auto left = ad::value_and_jacobian(function, input);
    const auto right = ad::value_and_jacobian(NoncopyableFunction<T>{calls}, input);
    const auto reference = ad::value_and_jacobian(std::ref(function), input);
    const auto pointer = ad::value_and_jacobian(&free_function<T>, input);
    require(calls == 3, "callback copied, retained, or evaluated more than once");
    require(left.values[0] == T{4} && left.jacobian[0][0] == T{2}, "lvalue callable");
    require(right.values[0] == T{6} && right.jacobian[0][0] == T{3}, "rvalue callable");
    require(reference.jacobian[0][0] == T{2}, "reference_wrapper callable");
    require(pointer.jacobian[0][0] == T{1}, "function pointer callable");
}

template <typename T>
void test_repeated_calls_ownership() {
    std::array<T, 2> input{T{2}, T{3}};
    const auto original = input;
    const auto function = [](const auto& p) { return std::array{p[0] * p[1], p[0] + p[1]}; };
    const auto first = ad::value_and_jacobian(function, input);
    require(input == original, "evaluation modified caller inputs");
    input = {T{5}, T{7}};
    auto second = ad::value_and_jacobian(function, input);
    require(first.values[0] == T{6}, "result aliases input or later workspace");
    require(first.jacobian[0] == std::array<T, 2>{T{3}, T{2}}, "first result overwritten");
    require(second.values[0] == T{35}, "repeated primal evaluation");
    require(second.jacobian[0] == std::array<T, 2>{T{7}, T{5}}, "stale derivative seeds");
    second.jacobian[0][0] = T{-99};
    const auto third = ad::value_and_jacobian(function, original);
    require(third.values == first.values && third.jacobian == first.jacobian,
            "calls or returned arrays share mutable storage");
}

struct CallbackFailure : std::runtime_error {
    explicit CallbackFailure(int error_code)
        : std::runtime_error("callback failure"), code(error_code) {}
    int code;
};

template <typename T>
void test_exception_propagation() {
    const std::array<T, 1> input{T{-2}};
    int calls = 0;
    bool caught = false;
    try {
        (void)ad::value_and_jacobian([&](const auto&) -> std::array<ad::Dual<T>, 1> {
            ++calls;
            throw CallbackFailure{42};
        }, input);
    } catch (const CallbackFailure& error) {
        caught = error.code == 42;
    }
    require(caught && calls == 1, "callback exception type/payload must be preserved");
    caught = false;
    try {
        (void)ad::value_and_jacobian([](const auto& p) {
            using std::log;
            return std::array{log(p[0])};
        }, input);
    } catch (const std::domain_error&) {
        caught = true;
    }
    require(caught && input[0] == T{-2}, "math domain errors propagate without mutation");
}

template <typename T>
void test_nonfinite_inputs() {
    const T inf = std::numeric_limits<T>::infinity();
    for (T invalid : std::array<T, 3>{inf, -inf, std::numeric_limits<T>::quiet_NaN()}) {
        for (std::size_t index = 0; index < 2; ++index) {
            std::array<T, 2> input{T{2}, T{3}};
            input[index] = invalid;
            int calls = 0;
            bool caught = false;
            try {
                (void)ad::value_and_jacobian([&](const auto& p) { ++calls; return p; }, input);
            } catch (const std::domain_error&) {
                caught = true;
            }
            require(caught && calls == 0, "nonfinite inputs must fail before callback evaluation");
            require(std::isnan(invalid) ? std::isnan(input[index]) : input[index] == invalid,
                    "invalid input was modified");
        }
    }
}

template <typename T>
void test_output_passthrough() {
    using Number = ad::Dual<T>;
    const T inf = std::numeric_limits<T>::infinity();
    const T nan = std::numeric_limits<T>::quiet_NaN();
    const auto result = ad::value_and_jacobian([&](const auto&) {
        // Deliberate invalid outputs test extraction, not scientific validity.
        return std::array{Number{inf}, Number{nan}, Number{T{1}, {inf}}};
    }, std::array<T, 1>{T{2}});
    require(std::isinf(result.values[0]) && std::isnan(result.values[1]) &&
                std::isinf(result.jacobian[2][0]), "nonfinite outputs must not be repaired");
}

void test_header_odr() {
    const auto separate = jacobian_from_separate_translation_unit();
    const auto local = ad::value_and_jacobian(Identity{}, std::array<double, 1>{3.0});
    require(local.values == separate.values && local.jacobian == separate.jacobian,
            "header self-containment/multiple-translation-unit integration");
}

template <typename T>
void run_typed_case(std::string_view name) {
    if (name == "identity_seeds") {
        test_identity_seeds<T>();
    } else if (name == "rectangular_analytic") {
        test_rectangular_analytic<T>();
    } else if (name == "constants_unused_inputs") {
        test_constants_unused_inputs<T>();
    } else if (name == "single_input_output") {
        test_single_input_output<T>();
    } else if (name == "callable_forwarding") {
        test_callable_forwarding<T>();
    } else if (name == "repeated_calls_ownership") {
        test_repeated_calls_ownership<T>();
    } else if (name == "exception_propagation") {
        test_exception_propagation<T>();
    } else if (name == "nonfinite_inputs") {
        test_nonfinite_inputs<T>();
    } else if (name == "output_passthrough") {
        test_output_passthrough<T>();
    } else {
        throw std::invalid_argument("unknown Jacobian test case");
    }
}
} // namespace

int main(int argc, char** argv) {
    try {
        require(argc == 2, "provide one named Jacobian test case");
        const std::string_view name{argv[1]};
        if (name == "header_odr") {
            test_header_odr();
        } else if (name == "nested_hessian") {
            test_nested_hessian_driver();
        } else {
            run_typed_case<float>(name);
            run_typed_case<double>(name);
            run_typed_case<long double>(name);
        }
        std::cout << "[PASS] " << name << '\n';
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "[FAIL] " << error.what() << '\n';
        return 1;
    }
}
