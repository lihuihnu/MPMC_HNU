#include <mpmc/ad/runtime_differentiate.hpp>
#include <mpmc/ad/differentiate.hpp>
#include <mpmc/ad/math.hpp>

#include <algorithm>
#include <array>
#include <cmath>
#include <concepts>
#include <cstddef>
#include <functional>
#include <iostream>
#include <limits>
#include <source_location>
#include <span>
#include <stdexcept>
#include <string>
#include <string_view>
#include <type_traits>
#include <utility>
#include <vector>

mpmc::ad::RuntimeValueAndJacobian<double> runtime_from_separate_translation_unit();

namespace {
namespace ad = mpmc::ad;

void require(bool condition, std::string_view message,
             std::source_location where = std::source_location::current()) {
    if (!condition) {
        throw std::runtime_error(std::string(where.file_name()) + ":" +
                                 std::to_string(where.line()) + ": " + std::string(message));
    }
}

template <typename T>
void near(T actual, T expected, T absolute = T{64} * std::numeric_limits<T>::epsilon()) {
    // Short, dimensionless expressions; tiny nonzero slopes use absolute=0.
    const T relative = T{64} * std::numeric_limits<T>::epsilon();
    require(std::isfinite(actual) && std::isfinite(expected) &&
                std::abs(actual - expected) <= absolute + relative * std::abs(expected),
            "runtime result differs from analytic reference");
}

template <typename Error, typename Function>
void expect_throw(Function&& function) {
    bool caught = false;
    try {
        std::invoke(std::forward<Function>(function));
    } catch (const Error&) {
        caught = true;
    }
    require(caught, "expected exception was not thrown");
}

template <typename Number>
void identity(std::span<const Number> input, std::span<Number> output) {
    require(input.size() == output.size(), "identity shape");
    std::copy(input.begin(), input.end(), output.begin());
}

using ProbeNumber = ad::Dual<double, 4>;
struct ValidCallback {
    void operator()(std::span<const ProbeNumber> x, std::span<ProbeNumber> y) & {
        y[0] = x[0];
    }
};
struct WrongReturn {
    int operator()(std::span<const ProbeNumber>, std::span<ProbeNumber>) const { return 0; }
};
struct MutableInputs {
    void operator()(std::span<ProbeNumber>, std::span<ProbeNumber>) const {}
};
struct OrdinaryOutputs {
    void operator()(std::span<const ProbeNumber>, std::span<double>) const {}
};
struct RvalueOnly {
    void operator()(std::span<const ProbeNumber>, std::span<ProbeNumber>) && {}
};

template <typename Function, typename Inputs = std::vector<double>>
concept CanRun = requires(Function&& f, const Inputs& x,
                          ad::RuntimeJacobianWorkspace<double, 4>& work) {
    ad::value_and_jacobian_runtime<4>(std::forward<Function>(f), x, 1, work);
};
template <std::size_t K>
concept CanWorkspace = requires { typename ad::RuntimeJacobianWorkspace<double, K>; };

static_assert(CanRun<ValidCallback> && CanRun<ValidCallback&>);
static_assert(CanRun<ValidCallback, std::array<double, 3>>);
static_assert(CanRun<ValidCallback, std::span<const double>>);
static_assert(!CanRun<WrongReturn> && !CanRun<MutableInputs> && !CanRun<OrdinaryOutputs>);
static_assert(!CanRun<RvalueOnly> && !CanRun<int>);
static_assert(!CanRun<ValidCallback, std::vector<float>>);
static_assert(!CanRun<ValidCallback, std::vector<ProbeNumber>>);
static_assert(!CanWorkspace<0> && CanWorkspace<1>);
static_assert(!std::is_copy_constructible_v<ad::RuntimeJacobianWorkspace<double, 4>>);
static_assert(!std::is_move_constructible_v<ad::RuntimeJacobianWorkspace<double, 4>>);
static_assert(std::same_as<ad::RuntimeJacobianWorkspace<double, 4>::Number, ProbeNumber>);

template <typename T, std::size_t K>
void check_block_seeding() {
    using Number = ad::Dual<T, K>;
    ad::RuntimeJacobianWorkspace<T, K> work;
    for (std::size_t n : std::array<std::size_t, 8>{1, 2, 3, 4, 5, 8, 9, 10}) {
        std::vector<T> x(n, T{2});
        std::size_t calls = 0;
        const auto result = ad::value_and_jacobian_runtime<K>(
            [&](std::span<const Number> input, std::span<Number> output) {
                require(input.size() == n && output.size() == n, "full shape on every pass");
                const std::size_t first = calls * K;
                for (std::size_t j = 0; j < n; ++j) {
                    require(input[j].value() == x[j], "primal order changed");
                    for (std::size_t lane = 0; lane < K; ++lane) {
                        const T expected = j == first + lane ? T{1} : T{0};
                        require(input[j].derivative(lane) == expected,
                                "stale seed or nonzero inactive tail lane");
                    }
                    output[j] = input[j];
                }
                ++calls;
            }, x, n, work);
        require(calls == n / K + (n % K != 0 ? 1U : 0U), "incorrect number of passes");
        require(result.values == x && result.input_count == n && result.output_count == n,
                "identity value/shape");
        for (std::size_t i = 0; i < n; ++i) {
            for (std::size_t j = 0; j < n; ++j) {
                require(result.jacobian[i * n + j] == (i == j ? T{1} : T{0}),
                        "identity row-major layout");
            }
        }
    }
}

template <typename T>
void test_block_seeding() {
    check_block_seeding<T, 1>();
    check_block_seeding<T, 3>();
    check_block_seeding<T, 4>();
    check_block_seeding<T, 8>();
}

// Analytic variable-size family: F_i = (i+1)*sum(x_j^2) + x_a*x_b + exp(x_a),
// a=i mod n, b=(i+1) mod n. The n=1 case exercises coincident a/b contributions.
template <typename T>
void test_rectangular_analytic() {
    using Number = ad::Dual<T, 4>;
    ad::RuntimeJacobianWorkspace<T, 4> work;
    for (std::size_t n : std::array<std::size_t, 5>{1, 2, 3, 5, 9}) {
        std::vector<T> x(n);
        T sum{};
        for (std::size_t j = 0; j < n; ++j) {
            x[j] = static_cast<T>(j + 1) / T{4};
            sum += x[j] * x[j];
        }
        for (std::size_t m : std::array<std::size_t, 3>{1, 3, n + 2}) {
            const auto result = ad::value_and_jacobian_runtime<4>(
                [](std::span<const Number> p, std::span<Number> output) {
                    Number squares{T{0}};
                    for (const auto& value : p) {
                        squares += value * value;
                    }
                    using std::exp;
                    for (std::size_t i = 0; i < output.size(); ++i) {
                        const auto a = i % p.size();
                        const auto b = (i + 1) % p.size();
                        output[i] = squares * static_cast<T>(i + 1) + p[a] * p[b] + exp(p[a]);
                    }
                }, x, m, work);
            require(result.jacobian.size() == n * m, "runtime rectangular shape");
            for (std::size_t i = 0; i < m; ++i) {
                const auto a = i % n;
                const auto b = (i + 1) % n;
                near(result.values[i], static_cast<T>(i + 1) * sum + x[a] * x[b] + std::exp(x[a]));
                for (std::size_t j = 0; j < n; ++j) {
                    T slope = T{2} * static_cast<T>(i + 1) * x[j];
                    if (j == a) {
                        slope += x[b] + std::exp(x[a]);
                    }
                    if (j == b) {
                        slope += x[a];
                    }
                    near(result.jacobian[i * n + j], slope);
                }
            }
        }
    }
}

template <typename Number>
auto small_function(std::span<const Number> p) {
    using std::exp;
    using std::log;
    return std::array{p[0] * p[4] + p[1], exp(p[2]) + log(p[3]), p[4] * p[4]};
}

template <typename T>
void test_fixed_crosscheck() {
    const std::array<T, 5> x{T{1}, T{2}, T{0}, T{4}, T{3}};
    const auto fixed = ad::value_and_jacobian([](const auto& p) {
        using Number = typename std::remove_cvref_t<decltype(p)>::value_type;
        return small_function(std::span<const Number>{p});
    }, x);
    using Number = ad::Dual<T, 3>;
    ad::RuntimeJacobianWorkspace<T, 3> work;
    const auto result = ad::value_and_jacobian_runtime<3>(
        [](std::span<const Number> p, std::span<Number> output) {
            const auto values = small_function(p);
            std::copy(values.begin(), values.end(), output.begin());
        }, x, 3, work);
    // The fixed driver is an additional integration oracle, not the only oracle.
    const std::array<std::array<T, 5>, 3> analytic{{
        {T{3}, T{1}, T{0}, T{0}, T{1}}, {T{0}, T{0}, T{1}, T{0.25}, T{0}},
        {T{0}, T{0}, T{0}, T{0}, T{6}}}};
    for (std::size_t i = 0; i < 3; ++i) {
        near(result.values[i], fixed.values[i]);
        for (std::size_t j = 0; j < 5; ++j) {
            near(result.jacobian[i * 5 + j], analytic[i][j]);
            near(result.jacobian[i * 5 + j], fixed.jacobian[i][j]);
        }
    }
}

template <typename T>
void test_constants_unused() {
    using Number = ad::Dual<T, 4>;
    ad::RuntimeJacobianWorkspace<T, 4> work;
    const auto result = ad::value_and_jacobian_runtime<4>(
        [](std::span<const Number> p, std::span<Number> output) {
            output[0] = Number{T{7}};
            output[1] = p[0] * T{2};
        }, std::vector<T>(9, T{3}), 2, work);
    require(result.values == std::vector<T>{T{7}, T{6}}, "constant outputs");
    for (std::size_t j = 0; j < 9; ++j) {
        require(result.jacobian[j] == T{0}, "constant row");
        require(result.jacobian[9 + j] == (j == 0 ? T{2} : T{0}), "unused columns");
    }
}

template <typename T>
void test_resize_ownership() {
    using Number = ad::Dual<T, 4>;
    ad::RuntimeJacobianWorkspace<T, 4> work;
    std::vector<T> x(12, T{2});
    const auto first = ad::value_and_jacobian_runtime<4>(identity<Number>, x, x.size(), work);
    const auto in_capacity = work.input_capacity();
    const auto out_capacity = work.output_capacity();
    for (std::size_t n : std::array<std::size_t, 4>{2, 7, 1, 12}) {
        x.assign(n, static_cast<T>(n));
        const auto result = ad::value_and_jacobian_runtime<4>(identity<Number>, x, n, work);
        require(result.values == x, "resize used stale values");
        require(work.input_capacity() == in_capacity && work.output_capacity() == out_capacity,
                "smaller problems should reuse AD buffer capacity");
        require(first.values == std::vector<T>(12, T{2}), "result aliases workspace");
    }
    const auto owned = [] {
        ad::RuntimeJacobianWorkspace<T, 4> temporary;
        return ad::value_and_jacobian_runtime<4>(
            identity<Number>, std::array<T, 1>{T{5}}, 1, temporary);
    }();
    require(owned.values[0] == T{5} && owned.jacobian[0] == T{1}, "result lifetime");
}

template <typename T>
struct NoncopyableCallback {
    explicit NoncopyableCallback(std::size_t& calls) : calls_(calls) {}
    NoncopyableCallback(const NoncopyableCallback&) = delete;
    NoncopyableCallback& operator=(const NoncopyableCallback&) = delete;
    void operator()(std::span<const ad::Dual<T, 4>> x, std::span<ad::Dual<T, 4>> y) & {
        ++calls_;
        y[0] = x.front() + x.back();
    }
    void operator()(std::span<const ad::Dual<T, 4>>, std::span<ad::Dual<T, 4>>) && = delete;
    std::size_t& calls_;
};

template <typename T>
void test_callback_contract() {
    using Number = ad::Dual<T, 4>;
    ad::RuntimeJacobianWorkspace<T, 4> work;
    const std::vector<T> x(5, T{2});
    std::size_t calls = 0;
    NoncopyableCallback<T> callback{calls};
    const auto a = ad::value_and_jacobian_runtime<4>(callback, x, 1, work);
    const auto b = ad::value_and_jacobian_runtime<4>(NoncopyableCallback<T>{calls}, x, 1, work);
    const auto c = ad::value_and_jacobian_runtime<4>(std::ref(callback), x, 1, work);
    require(calls == 6 && a.jacobian == b.jacobian && b.jacobian == c.jacobian,
            "callback must be a repeated lvalue, without copying or moving it");
    const auto pointer = ad::value_and_jacobian_runtime<4>(&identity<Number>, x, 5, work);
    require(pointer.values == x, "function pointer");
    using Pointer = void (*)(std::span<const Number>, std::span<Number>);
    expect_throw<std::invalid_argument>([&] {
        (void)ad::value_and_jacobian_runtime<4>(Pointer{}, x, 1, work);
    });
}

struct CallbackFailure : std::runtime_error {
    explicit CallbackFailure(int value)
        : std::runtime_error("intentional callback failure"), code(value) {}
    int code;
};

template <typename T>
void test_failure_recovery() {
    using Number = ad::Dual<T, 4>;
    ad::RuntimeJacobianWorkspace<T, 4> work;
    const std::vector<T> x(9, T{2});
    auto published = ad::value_and_jacobian_runtime<4>(identity<Number>, x, 9, work);
    const auto original = published;
    std::size_t calls = 0;
    bool caught = false;
    try {
        published = ad::value_and_jacobian_runtime<4>(
            [&](std::span<const Number> p, std::span<Number> y) {
                if (++calls == 2) {
                    throw CallbackFailure{42};
                }
                identity(p, y);
            }, x, 9, work);
    } catch (const CallbackFailure& error) {
        caught = error.code == 42;
    }
    require(caught && calls == 2 && published.values == original.values &&
                published.jacobian == original.jacobian,
            "partial matrix published or exception changed");
    const auto recovered = ad::value_and_jacobian_runtime<4>(identity<Number>, x, 9, work);
    require(recovered.jacobian == original.jacobian, "failed call contaminated workspace");
    expect_throw<std::domain_error>([&] {
        (void)ad::value_and_jacobian_runtime<4>(
            [](std::span<const Number> p, std::span<Number> y) { y[0] = ad::log(p[0]); },
            std::array<T, 1>{T{-1}}, 1, work);
    });
}

template <typename T>
void test_input_validation() {
    using Number = ad::Dual<T, 4>;
    ad::RuntimeJacobianWorkspace<T, 4> work;
    std::size_t calls = 0;
    for (T bad : std::array<T, 3>{std::numeric_limits<T>::infinity(),
                                 -std::numeric_limits<T>::infinity(),
                                 std::numeric_limits<T>::quiet_NaN()}) {
        for (std::size_t j = 0; j < 5; ++j) {
            std::vector<T> input(5, T{2});
            input[j] = bad;
            expect_throw<std::domain_error>([&] {
                (void)ad::value_and_jacobian_runtime<4>(
                    [&](std::span<const Number> p, std::span<Number> y) {
                        ++calls;
                        identity(p, y);
                    }, input, 5, work);
            });
            require(std::isnan(bad) ? std::isnan(input[j]) : input[j] == bad, "input mutation");
        }
    }
    require(calls == 0 && work.input_capacity() == 0,
            "validation must precede allocation/callback");
}

template <typename T>
void test_size_limits() {
    using Number = ad::Dual<T, 4>;
    ad::RuntimeJacobianWorkspace<T, 4> work;
    std::size_t calls = 0;
    const auto callback = [&](std::span<const Number>, std::span<Number> y) {
        ++calls;
        std::fill(y.begin(), y.end(), Number{T{1}});
    };
    const std::vector<T> x(4, T{2});
    expect_throw<std::invalid_argument>([&] {
        (void)ad::value_and_jacobian_runtime<4>(callback, std::span<const T>{}, 1, work);
    });
    expect_throw<std::invalid_argument>([&] {
        (void)ad::value_and_jacobian_runtime<4>(callback, x, 0, work);
    });
    const auto max = std::numeric_limits<std::size_t>::max();
    for (auto limits : std::array<ad::RuntimeJacobianLimits, 4>{{
             {3, max, max}, {max, 2, max}, {max, max, 11}, {max, max, 0}}}) {
        expect_throw<std::length_error>([&] {
            (void)ad::value_and_jacobian_runtime<4>(callback, x, 3, work, limits);
        });
    }
    expect_throw<std::length_error>([&] {
        (void)ad::value_and_jacobian_runtime<4>(callback, x, max, work);
    });
    expect_throw<std::length_error>([&] {
        (void)ad::value_and_jacobian_runtime<4>(callback, std::array<T, 1>{T{2}}, max, work);
    });
    require(calls == 0 && work.input_capacity() == 0,
            "limits must fail before allocation/callback");
    const auto result = ad::value_and_jacobian_runtime<4>(callback, x, 3, work, {4, 3, 12});
    require(calls == 1 && result.jacobian.size() == 12, "inclusive size limits");
}

template <typename T>
void test_output_validation() {
    using Number = ad::Dual<T, 4>;
    ad::RuntimeJacobianWorkspace<T, 4> work;
    const std::vector<T> x(5, T{2});
    std::size_t calls = 0;
    // Deliberately omit a previously assigned output on the second pass.
    expect_throw<std::domain_error>([&] {
        (void)ad::value_and_jacobian_runtime<4>(
            [&](std::span<const Number> p, std::span<Number> y) {
                y[0] = p[0];
                if (++calls == 1) {
                    y[1] = Number{T{7}};
                }
            }, x, 2, work);
    });
    require(calls == 2, "missing output must not reuse a previous block's value");
    for (T bad : std::array<T, 3>{std::numeric_limits<T>::infinity(),
                                 -std::numeric_limits<T>::infinity(),
                                 std::numeric_limits<T>::quiet_NaN()}) {
        expect_throw<std::domain_error>([&] {
            (void)ad::value_and_jacobian_runtime<4>(
                [&](std::span<const Number>, std::span<Number> y) { y[0] = Number{bad}; },
                x, 1, work);
        });
        expect_throw<std::domain_error>([&] {
            (void)ad::value_and_jacobian_runtime<4>(
                [&](std::span<const Number>, std::span<Number> y) {
                    typename Number::Gradient gradient{};
                    gradient[0] = bad;
                    y[0] = Number{T{1}, gradient};
                }, x, 1, work);
        });
    }
    const auto recovered = ad::value_and_jacobian_runtime<4>(identity<Number>, x, 5, work);
    require(recovered.values == x, "output failure poisoned next call");
}

template <typename T>
void test_primal_consistency() {
    using Number = ad::Dual<T, 4>;
    ad::RuntimeJacobianWorkspace<T, 4> work;
    const std::vector<T> x(5, T{2});
    for (bool signed_zero : std::array<bool, 2>{false, true}) {
        std::size_t calls = 0;
        expect_throw<std::runtime_error>([&] {
            (void)ad::value_and_jacobian_runtime<4>(
                [&](std::span<const Number>, std::span<Number> y) {
                    const T value = signed_zero ? (calls == 0 ? T{0} : -T{0})
                                                : static_cast<T>(calls);
                    y[0] = Number{value};
                    ++calls;
                }, x, 1, work);
        });
        require(calls == 2, "primal mismatch not detected on second pass");
    }
}

template <typename T>
void test_workspace_isolation() {
    using Number = ad::Dual<T, 4>;
    ad::RuntimeJacobianWorkspace<T, 4> outer;
    ad::RuntimeJacobianWorkspace<T, 4> inner;
    const std::vector<T> x(5, T{3});
    expect_throw<std::logic_error>([&] {
        (void)ad::value_and_jacobian_runtime<4>(
            [&](std::span<const Number>, std::span<Number>) {
                (void)ad::value_and_jacobian_runtime<4>(identity<Number>, x, 5, outer);
            }, x, 1, outer);
    });
    const auto result = ad::value_and_jacobian_runtime<4>(
        [&](std::span<const Number> p, std::span<Number> y) {
            const auto separate = ad::value_and_jacobian_runtime<4>(
                identity<Number>, std::array<T, 1>{T{2}}, 1, inner);
            require(separate.jacobian[0] == T{1}, "independent nested workspace");
            y[0] = p.front() + p.back();
        }, x, 1, outer);
    require(result.jacobian == std::vector<T>{T{1}, T{0}, T{0}, T{0}, T{1}},
            "nested evaluation changed outer seeds");
}

template <typename T>
void test_scaled_derivatives() {
    using Number = ad::Dual<T, 4>;
    ad::RuntimeJacobianWorkspace<T, 4> work;
    const T slope = std::sqrt(std::numeric_limits<T>::min());
    const auto result = ad::value_and_jacobian_runtime<4>(
        [&](std::span<const Number> p, std::span<Number> y) {
            y[0] = p.front() * slope + p.back() * slope;
        }, std::vector<T>(5, T{1}), 1, work);
    near(result.jacobian[0], slope, T{0});
    near(result.jacobian[4], slope, T{0});
    near(result.values[0], T{2} * slope, T{0});
}

void test_header_odr() {
    const auto result = runtime_from_separate_translation_unit();
    require(result.input_count == 5 && result.output_count == 5, "separate TU shape");
    require(result.values == std::vector<double>(5, 2.0), "separate TU values");
    for (std::size_t i = 0; i < 5; ++i) {
        for (std::size_t j = 0; j < 5; ++j) {
            require(result.jacobian[i * 5 + j] == (i == j ? 1.0 : 0.0), "separate TU matrix");
        }
    }
}

template <typename T>
void run_case(std::string_view name) {
    if (name == "block_seeding") {
        test_block_seeding<T>();
    } else if (name == "rectangular_analytic") {
        test_rectangular_analytic<T>();
    } else if (name == "fixed_crosscheck") {
        test_fixed_crosscheck<T>();
    } else if (name == "constants_unused") {
        test_constants_unused<T>();
    } else if (name == "resize_ownership") {
        test_resize_ownership<T>();
    } else if (name == "callback_contract") {
        test_callback_contract<T>();
    } else if (name == "failure_recovery") {
        test_failure_recovery<T>();
    } else if (name == "input_validation") {
        test_input_validation<T>();
    } else if (name == "size_limits") {
        test_size_limits<T>();
    } else if (name == "output_validation") {
        test_output_validation<T>();
    } else if (name == "primal_consistency") {
        test_primal_consistency<T>();
    } else if (name == "workspace_isolation") {
        test_workspace_isolation<T>();
    } else if (name == "scaled_derivatives") {
        test_scaled_derivatives<T>();
    } else {
        throw std::invalid_argument("unknown runtime test case");
    }
}
} // namespace

int main(int argc, char** argv) {
    try {
        require(argc == 2, "provide one runtime test case");
        const std::string_view name{argv[1]};
        if (name == "header_odr") {
            test_header_odr();
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
