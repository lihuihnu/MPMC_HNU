#ifndef MPMC_AD_RUNTIME_DIFFERENTIATE_HPP
#define MPMC_AD_RUNTIME_DIFFERENTIATE_HPP

#include <mpmc/ad/dual.hpp>

#include <algorithm>
#include <cmath>
#include <concepts>
#include <cstddef>
#include <functional>
#include <limits>
#include <span>
#include <stdexcept>
#include <type_traits>
#include <utility>
#include <vector>

namespace mpmc::ad {

/// Logical problem limits, checked before allocation or callback invocation.
/// Defaults impose only representability/container limits, NOT a service quota.
/// The caller must set appropriate quotas for untrusted or interactive requests.
struct RuntimeJacobianLimits {
    std::size_t max_inputs = std::numeric_limits<std::size_t>::max();
    std::size_t max_outputs = std::numeric_limits<std::size_t>::max();
    std::size_t max_jacobian_entries = std::numeric_limits<std::size_t>::max();
};

/// Owned, complete result for F: R^n -> R^m, independent of the workspace.
/// jacobian[i * input_count + j] = d F_i / d x_j (contiguous row-major order).
/// Shape and vectors describe the returned result; callers modifying these
/// public fields are responsible for preserving their consistency.
template <std::floating_point T>
    requires std::same_as<T, std::remove_cv_t<T>>
struct RuntimeValueAndJacobian {
    std::size_t input_count{};
    std::size_t output_count{};
    std::vector<T> values;
    std::vector<T> jacobian;
};

namespace detail {

template <typename Function, typename T, std::size_t K>
concept RuntimeJacobianFunction =
    std::invocable<Function&, std::span<const Dual<T, K>>, std::span<Dual<T, K>>> &&
    std::same_as<std::invoke_result_t<Function&, std::span<const Dual<T, K>>,
                                    std::span<Dual<T, K>>>, void>;

} // namespace detail

/// Reusable input/output AD storage; K is a direction width, not a variable limit.
/// One workspace belongs to one synchronous evaluation at a time. Separate
/// workspaces may be used independently. Concurrent access to the SAME workspace
/// is forbidden; the reentry guard below is not a cross-thread synchronization.
/// Noncopyable/nonmovable so a callback cannot invalidate active storage by moving it.
template <std::floating_point T, std::size_t K>
    requires(K > 0 && std::same_as<T, std::remove_cv_t<T>>)
class RuntimeJacobianWorkspace {
public:
    using Scalar = T;
    using Number = Dual<T, K>;
    static constexpr std::size_t direction_count = K;

    // NaN sentinels detect missing output assignments, not just stale values.
    static_assert(std::numeric_limits<T>::has_quiet_NaN,
                  "Runtime Jacobian output validation requires quiet NaN support");

    RuntimeJacobianWorkspace() = default;
    RuntimeJacobianWorkspace(const RuntimeJacobianWorkspace&) = delete;
    RuntimeJacobianWorkspace& operator=(const RuntimeJacobianWorkspace&) = delete;
    RuntimeJacobianWorkspace(RuntimeJacobianWorkspace&&) = delete;
    RuntimeJacobianWorkspace& operator=(RuntimeJacobianWorkspace&&) = delete;

    [[nodiscard]] std::size_t input_capacity() const noexcept { return variables_.capacity(); }
    [[nodiscard]] std::size_t output_capacity() const noexcept { return outputs_.capacity(); }

    /// Implementation entry used by value_and_jacobian_runtime; same contract.
    /// The callback is invoked as an LVALUE on every pass, even if supplied as a
    /// temporary. It must be repeatable, preserve AD dependence and assign EVERY
    /// output (including constants) from scratch; spans may not be retained.
    ///
    /// This checked runtime driver requires finite values and active derivatives.
    /// It does not change the more permissive output policy of the fixed driver.
    /// Exact primal agreement (including signed zero) is required across passes;
    /// this detects some stateful callbacks, but is not a proof of callback purity.
    template <typename Function>
        requires detail::RuntimeJacobianFunction<Function, T, K>
    [[nodiscard]] RuntimeValueAndJacobian<T>
    evaluate(Function&& function, std::span<const T> inputs, std::size_t output_count,
             RuntimeJacobianLimits limits = {}) {
        ActiveCall active{active_};
        const std::size_t n = inputs.size();
        const std::size_t m = output_count;
        if (n == 0 || m == 0) {
            throw std::invalid_argument("mpmc::ad::runtime: input/output counts must be positive");
        }
        // Divide BEFORE multiplying: no n*m overflow and no artificial huge span.
        const std::size_t scalar_limit = std::vector<T>{}.max_size();
        if (n > limits.max_inputs || m > limits.max_outputs ||
            m > limits.max_jacobian_entries / n || n > variables_.max_size() ||
            m > outputs_.max_size() || m > scalar_limit / n) {
            throw std::length_error("mpmc::ad::runtime: problem exceeds size/resource limits");
        }
        for (T value : inputs) {
            if (!std::isfinite(value)) {
                throw std::domain_error("mpmc::ad::runtime: nonfinite input");
            }
        }
        if constexpr (std::is_pointer_v<std::remove_reference_t<Function>>) {
            if (function == nullptr) {
                throw std::invalid_argument("mpmc::ad::runtime: null callback");
            }
        }

        // Snapshot all primals once. Subsequent passes never reread the input view.
        // resize reuses capacity; a failed call leaves reusable, but unspecified,
        // scratch contents. The next call always overwrites them before use.
        variables_.resize(n);
        outputs_.resize(m);
        for (std::size_t j = 0; j < n; ++j) {
            variables_[j] = Number{inputs[j]};
        }
        RuntimeValueAndJacobian<T> result{n, m, std::vector<T>(m), std::vector<T>(m * n)};
        typename Number::Gradient unset_gradient{};
        const T nan = std::numeric_limits<T>::quiet_NaN();
        unset_gradient.fill(nan);
        const Number unset{nan, unset_gradient};

        for (std::size_t first = 0; first < n;) {
            const std::size_t width = std::min(K, n - first);
            for (std::size_t lane = 0; lane < width; ++lane) {
                const std::size_t j = first + lane;
                variables_[j] = Number::variable(variables_[j].value(), lane);
            }
            std::fill(outputs_.begin(), outputs_.end(), unset);
            std::invoke(function, std::span<const Number>{variables_}, std::span<Number>{outputs_});

            for (std::size_t i = 0; i < m; ++i) {
                const T value = outputs_[i].value();
                if (!std::isfinite(value)) {
                    throw std::domain_error("mpmc::ad::runtime: nonfinite or unwritten output");
                }
                if (first == 0) {
                    result.values[i] = value;
                } else if (value != result.values[i] ||
                           std::signbit(value) != std::signbit(result.values[i])) {
                    throw std::runtime_error("mpmc::ad::runtime: primals changed between passes");
                }
                for (std::size_t lane = 0; lane < width; ++lane) {
                    const T derivative = outputs_[i].derivatives()[lane];
                    if (!std::isfinite(derivative)) {
                        throw std::domain_error("mpmc::ad::runtime: nonfinite active derivative");
                    }
                    result.jacobian[i * n + first + lane] = derivative;
                }
            }
            // Only this block was active. Reset it before seeding the next block;
            // this avoids clearing all n*K lanes on every pass. Inactive lanes in
            // the final partial block remain zero; no fictitious variables exist.
            for (std::size_t lane = 0; lane < width; ++lane) {
                auto& variable = variables_[first + lane];
                variable = Number{variable.value()};
            }
            first += width; // width <= n-first, so the final increment cannot overflow.
        }
        return result; // Exceptions before here never publish a partial matrix.
    }

private:
    struct ActiveCall {
        explicit ActiveCall(bool& active) : active_(active) {
            if (active_) {
                throw std::logic_error("mpmc::ad::runtime: workspace already in use");
            }
            active_ = true;
        }
        ~ActiveCall() { active_ = false; }
        ActiveCall(const ActiveCall&) = delete;
        ActiveCall& operator=(const ActiveCall&) = delete;
        bool& active_;
    };

    std::vector<Number> variables_;
    std::vector<Number> outputs_;
    bool active_ = false;
};

/// Runtime n/m, fixed direction width K: ceil(n/K) callback invocations.
/// Callback signature: void(span<const Dual<T,K>>, span<Dual<T,K>>).
/// T is inferred from the explicit workspace; vector/array/span inputs can
/// convert to the read-only span without encoding n in the function type.
///
/// invalid_argument: zero dimensions or a null function pointer.
/// length_error: shape overflow, container limit or configured quota.
/// domain_error: nonfinite inputs, unwritten/nonfinite outputs or active slopes.
/// runtime_error: primal disagreement across blocks; logic_error: workspace reentry.
/// Allocation and callback exceptions propagate; external side effects cannot be
/// rolled back. Workspaces are reusable after an exception, results do not alias
/// them, and the wrapper never repairs, clips or normalizes mathematical data.
///
/// O(K*(n+m)) AD workspace, O(m*n) dense result; the callback's own costs are extra.
/// Results allocate per call; workspace growth may allocate, not per scalar or
/// per block. No realtime latency or relative performance guarantee is implied.
template <std::size_t K, typename Function, std::floating_point T>
    requires(K > 0 && std::same_as<T, std::remove_cv_t<T>> &&
             detail::RuntimeJacobianFunction<Function, T, K>)
[[nodiscard]] RuntimeValueAndJacobian<T> value_and_jacobian_runtime(
    Function&& function, std::type_identity_t<std::span<const T>> inputs,
    std::size_t output_count, RuntimeJacobianWorkspace<T, K>& workspace,
    RuntimeJacobianLimits limits = {}) {
    return workspace.evaluate(std::forward<Function>(function), inputs, output_count, limits);
}

} // namespace mpmc::ad

#endif // MPMC_AD_RUNTIME_DIFFERENTIATE_HPP
