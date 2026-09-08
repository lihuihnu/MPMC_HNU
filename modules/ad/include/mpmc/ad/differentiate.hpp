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

/// Owned result for F: R^N -> R^M. Rows are outputs; columns are inputs.
/// jacobian[i][j] = d F_i / d x_j, in the caller's original input order.
/// Nested arrays are indexed by row; no flat-buffer or binary ABI is promised.
template <std::floating_point T, std::size_t N, std::size_t M>
    requires(N > 0 && M > 0 && std::same_as<T, std::remove_cv_t<T>>)
struct ValueAndJacobian {
    using Scalar = T;
    using Values = std::array<T, M>;
    using Jacobian = std::array<std::array<T, N>, M>;
    static constexpr std::size_t input_count = N;
    static constexpr std::size_t output_count = M;

    Values values{};
    Jacobian jacobian{};
};

namespace detail {

// Accept only an owning, nonempty array of the EXACT seeded scalar type.
// In particular, ordinary floating-point outputs must not silently acquire
// zero derivatives; reference/view outputs would introduce lifetime contracts.
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

/// Evaluate a fixed-size vector function and its full first-order Jacobian.
///
/// The callback receives const std::array<Dual<T, N>, N>& (or its own copy)
/// and returns std::array<Dual<T, N>, M> by value, M > 0. N and M are inferred.
/// Constants in the output must be explicitly constructed as Dual constants.
/// All N identity directions propagate in ONE callback invocation, not N calls.
///
/// Nonfinite inputs throw std::domain_error BEFORE invoking the callback.
/// Callback exceptions propagate unchanged. Output values/derivatives are copied
/// without clipping or finiteness repair; their validity remains the caller's
/// responsibility. An AD result type cannot detect deliberate value() stripping.
///
/// Inputs and results own separate storage. No callback is copied or retained;
/// forwarding preserves its value category. Callback side effects (including
/// writes through captures) cannot be rolled back. Null callable pointers and
/// callbacks retaining references into the temporary seed array are invalid use.
///
/// Only dual.hpp and the standard library are needed; include math.hpp in the
/// caller when needed. Local fixed-size arrays avoid wrapper heap allocation,
/// but occupy O(N*N + M*N) storage; intended for small/local dense problems.
/// No domain-specific normalization, constraints, or variable scaling is added.
template <typename Function, std::floating_point T, std::size_t N>
    requires(N > 0 && std::same_as<T, std::remove_cv_t<T>> &&
             std::invocable<Function, const std::array<Dual<T, N>, N>&> &&
             detail::JacobianOutput<
                 std::invoke_result_t<Function, const std::array<Dual<T, N>, N>&>,
                 Dual<T, N>>::valid)
[[nodiscard]] auto value_and_jacobian(Function&& function, const std::array<T, N>& inputs) {
    using Number = Dual<T, N>;
    using SeedArray = std::array<Number, N>;
    using Output = std::invoke_result_t<Function, const SeedArray&>;
    constexpr std::size_t output_count = detail::JacobianOutput<Output, Number>::size;

    const SeedArray seeded = [&inputs] {
        SeedArray variables{};
        for (std::size_t j = 0; j < N; ++j) {
            if (!std::isfinite(inputs[j])) {
                throw std::domain_error("mpmc::ad::value_and_jacobian: nonfinite input");
            }
            variables[j] = Number::variable(inputs[j], j);
        }
        return variables;
    }();

    // Result type inspection above is unevaluated; this is the only invocation.
    const auto outputs = std::invoke(std::forward<Function>(function), seeded);
    ValueAndJacobian<T, N, output_count> result{};
    for (std::size_t i = 0; i < output_count; ++i) {
        result.values[i] = outputs[i].value();
        result.jacobian[i] = outputs[i].derivatives();
    }
    return result;
}

} // namespace mpmc::ad

#endif // MPMC_AD_DIFFERENTIATE_HPP
