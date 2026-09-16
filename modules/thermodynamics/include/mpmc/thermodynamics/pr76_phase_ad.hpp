#ifndef MPMC_THERMODYNAMICS_PR76_PHASE_AD_HPP
#define MPMC_THERMODYNAMICS_PR76_PHASE_AD_HPP

#include <mpmc/ad/math.hpp>
#include <mpmc/ad/runtime_differentiate.hpp>
#include <mpmc/thermodynamics/pr76_phase.hpp>

#include <concepts>
#include <cstddef>
#include <limits>
#include <span>
#include <stdexcept>
#include <type_traits>
#include <vector>

namespace mpmc::thermodynamics {

/// Optional direct first-order AD driver for Pr76Phase.
///
/// This header is intentionally separate from pr76_phase.hpp so plain scalar
/// thermodynamics users do not acquire an AD dependency. K is the block width,
/// not a component-count limit; runtime component counts are propagated in
/// ceil(input_count/K) Dual passes by the existing AD runtime driver.
///
/// Output order for both coordinate systems:
///   [Z, ln(phi_0), ..., ln(phi_{n-1})]
///
/// Full input order:
///   [p, T, x_0, ..., x_{n-1}]
/// and the x columns are the formal full-composition derivatives already
/// supported by Pr76Phase::evaluate_full (the primal composition is normalized,
/// but derivative seeds are independent).
///
/// Reduced input order:
///   [p, T, x_0, ..., x_{n-2}], x_{n-1}=1-sum(x_0...x_{n-2})
/// so every reduced composition column includes the dependent-last chain rule.
///
/// The selected root index follows Pr76Phase semantics and is re-resolved at the
/// same primal state on every AD block. Near-multiple/ill-conditioned root errors
/// are preserved; this adapter never finite-differences, clips, or changes a root.
template <std::floating_point T, std::size_t K>
    requires(K > 0 && std::same_as<T, std::remove_cv_t<T>>)
class Pr76PhaseAdWorkspace {
public:
    using Scalar = T;
    using Number = ad::Dual<T, K>;
    static constexpr std::size_t direction_count = K;

    Pr76PhaseAdWorkspace() = default;
    Pr76PhaseAdWorkspace(const Pr76PhaseAdWorkspace&) = delete;
    Pr76PhaseAdWorkspace& operator=(const Pr76PhaseAdWorkspace&) = delete;
    Pr76PhaseAdWorkspace(Pr76PhaseAdWorkspace&&) = delete;
    Pr76PhaseAdWorkspace& operator=(Pr76PhaseAdWorkspace&&) = delete;

    [[nodiscard]] std::size_t input_capacity() const noexcept { return inputs_.capacity(); }

    [[nodiscard]] ad::RuntimeValueAndJacobian<T> value_and_jacobian_full(
        const Pr76Phase<T>& model, T pressure_pa, T temperature_k,
        std::span<const T> fractions, std::size_t root_index,
        Pr76RootOptions options = {}, ad::RuntimeJacobianLimits limits = {}) {
        return evaluate_impl<false>(model, pressure_pa, temperature_k, fractions,
                                    root_index, options, limits);
    }

    [[nodiscard]] ad::RuntimeValueAndJacobian<T> value_and_jacobian_reduced(
        const Pr76Phase<T>& model, T pressure_pa, T temperature_k,
        std::span<const T> independent_fractions, std::size_t root_index,
        Pr76RootOptions options = {}, ad::RuntimeJacobianLimits limits = {}) {
        return evaluate_impl<true>(model, pressure_pa, temperature_k, independent_fractions,
                                   root_index, options, limits);
    }

private:
    template <bool Reduced>
    [[nodiscard]] ad::RuntimeValueAndJacobian<T> evaluate_impl(
        const Pr76Phase<T>& model, T pressure_pa, T temperature_k,
        std::span<const T> fractions, std::size_t root_index,
        Pr76RootOptions options, ad::RuntimeJacobianLimits limits) {
        const std::size_t n = model.size();
        if (n == 0) {
            throw std::invalid_argument("Pr76PhaseAdWorkspace: moved-from/empty model");
        }
        const std::size_t expected = Reduced ? n - 1 : n;
        if (fractions.size() != expected) {
            throw std::invalid_argument("Pr76PhaseAdWorkspace: composition size does not match coordinates");
        }
        if (expected > inputs_.max_size() - std::size_t{2} ||
            n == std::numeric_limits<std::size_t>::max()) {
            throw std::length_error("Pr76PhaseAdWorkspace: problem shape is not representable");
        }
        const std::size_t input_count = expected + 2;
        const std::size_t output_count = n + 1;
        if (input_count > limits.max_inputs || output_count > limits.max_outputs ||
            output_count > limits.max_jacobian_entries / input_count) {
            throw std::length_error("Pr76PhaseAdWorkspace: problem exceeds AD resource limits");
        }

        inputs_.resize(input_count);
        inputs_[0] = pressure_pa;
        inputs_[1] = temperature_k;
        for (std::size_t i = 0; i < expected; ++i) {
            inputs_[i + 2] = fractions[i];
        }

        const auto callback = [&](std::span<const Number> in, std::span<Number> out) {
            Pr76PhaseValues<Number, T> phase = [&] {
                if constexpr (Reduced) {
                    return model.evaluate_reduced(in[0], in[1], in.subspan(2), root_index,
                                                  phase_workspace_, options);
                } else {
                    return model.evaluate_full(in[0], in[1], in.subspan(2), root_index,
                                               phase_workspace_, options);
                }
            }();
            out[0] = phase.z;
            for (std::size_t i = 0; i < n; ++i) {
                out[i + 1] = phase.ln_phi[i];
            }
        };

        return jacobian_workspace_.evaluate(callback, inputs_, output_count, limits);
    }

    ad::RuntimeJacobianWorkspace<T, K> jacobian_workspace_;
    Pr76PhaseWorkspace<Number> phase_workspace_;
    std::vector<T> inputs_;
};

} // namespace mpmc::thermodynamics

#endif // MPMC_THERMODYNAMICS_PR76_PHASE_AD_HPP
