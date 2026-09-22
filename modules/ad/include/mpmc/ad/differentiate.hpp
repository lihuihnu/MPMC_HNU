#ifndef MPMC_AD_DIFFERENTIATE_HPP
#define MPMC_AD_DIFFERENTIATE_HPP

#include <mpmc/ad/dual.hpp>

#include <array>
#include <cmath>
#include <concepts>
#include <cstddef>
#include <functional>
#include <stdexcept>
#include <type_traits>
#include <utility>

namespace mpmc::ad {

template <std::floating_point T, std::size_t N, std::size_t M>
    requires(N > 0 && M > 0 &&
             std::same_as<T, std::remove_cv_t<T>>)
struct ValueAndJacobian {
    using Scalar = T;
    using Values = std::array<T, M>;
    using Jacobian = std::array<std::array<T, N>, M>;
    static constexpr std::size_t input_count = N;
    static constexpr std::size_t output_count = M;

    Values values{};
    Jacobian jacobian{};
};

template <std::floating_point T, std::size_t N>
    requires(N > 0 &&
             std::same_as<T, std::remove_cv_t<T>>)
struct ValueGradientHessian {
    using Scalar = T;
    using Gradient = std::array<T, N>;
    using Hessian = std::array<std::array<T, N>, N>;
    static constexpr std::size_t input_count = N;

    T value{};
    Gradient gradient{};
    Hessian hessian{};
};

namespace detail {

template <typename Output, typename Number>
struct JacobianOutput {
    static constexpr bool valid = false;
};

template <typename Number, std::size_t M>
struct JacobianOutput<std::array<Number, M>, Number> {
    static constexpr bool valid = M > 0;
    static constexpr std::size_t size = M;
};

} // namespace detail

template <typename Function, std::floating_point T, std::size_t N>
    requires(N > 0 &&
             std::same_as<T, std::remove_cv_t<T>> &&
             std::invocable<
                 Function,
                 const std::array<Dual<T, N>, N>&> &&
             detail::JacobianOutput<
                 std::invoke_result_t<
                     Function,
                     const std::array<Dual<T, N>, N>&>,
                 Dual<T, N>>::valid)
[[nodiscard]] auto value_and_jacobian(
    Function&& function,
    const std::array<T, N>& inputs) {
    using Number = Dual<T, N>;
    using SeedArray = std::array<Number, N>;
    using Output =
        std::invoke_result_t<
            Function,
            const SeedArray&>;
    constexpr std::size_t output_count =
        detail::JacobianOutput<
            Output,
            Number>::size;

    const SeedArray seeded = [&inputs] {
        SeedArray variables{};
        for (std::size_t j = 0; j < N; ++j) {
            if (!std::isfinite(inputs[j])) {
                throw std::domain_error(
                    "mpmc::ad::value_and_jacobian: nonfinite input");
            }
            variables[j] =
                Number::variable(
                    inputs[j],
                    j);
        }
        return variables;
    }();

    const auto outputs =
        std::invoke(
            std::forward<Function>(function),
            seeded);

    ValueAndJacobian<T, N, output_count>
        result{};
    for (std::size_t i = 0;
         i < output_count;
         ++i) {
        result.values[i] =
            outputs[i].value();
        result.jacobian[i] =
            outputs[i].derivatives();
    }
    return result;
}

/// Evaluate one fixed-size scalar function, its gradient and full Hessian using
/// nested forward AD. The callback is invoked exactly once.
template <typename Function, std::floating_point T, std::size_t N>
    requires(N > 0 &&
             std::same_as<T, std::remove_cv_t<T>>)
[[nodiscard]] ValueGradientHessian<T, N>
value_gradient_hessian(
    Function&& function,
    const std::array<T, N>& inputs) {
    using Inner = Dual<T, N>;
    using Outer = Dual<Inner, N>;
    using SeedArray = std::array<Outer, N>;

    static_assert(
        std::invocable<
            Function,
            const SeedArray&>,
        "value_gradient_hessian callback must accept the nested seed array");
    static_assert(
        std::same_as<
            std::invoke_result_t<
                Function,
                const SeedArray&>,
            Outer>,
        "value_gradient_hessian callback must return one nested Dual scalar");

    SeedArray seeded{};
    for (std::size_t j = 0;
         j < N;
         ++j) {
        if (!std::isfinite(inputs[j])) {
            throw std::domain_error(
                "mpmc::ad::value_gradient_hessian: nonfinite input");
        }
        seeded[j] =
            Outer::variable(
                Inner::variable(
                    inputs[j],
                    j),
                j);
    }

    const Outer output =
        std::invoke(
            std::forward<Function>(function),
            seeded);

    ValueGradientHessian<T, N> result{};
    result.value =
        output.value().value();
    result.gradient =
        output.value().derivatives();

    for (std::size_t row = 0;
         row < N;
         ++row) {
        const Inner first =
            output.derivative(row);
        for (std::size_t column = 0;
             column < N;
             ++column) {
            result.hessian[row][column] =
                first.derivative(column);
        }
    }
    return result;
}

} // namespace mpmc::ad

#endif // MPMC_AD_DIFFERENTIATE_HPP
