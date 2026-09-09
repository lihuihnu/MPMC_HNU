#ifndef MPMC_THERMODYNAMICS_PR76_MIXTURE_HPP
#define MPMC_THERMODYNAMICS_PR76_MIXTURE_HPP

#include <mpmc/thermodynamics/pr76_pure.hpp>

#include <cmath>
#include <concepts>
#include <cstddef>
#include <limits>
#include <span>
#include <stdexcept>
#include <type_traits>
#include <utility>
#include <vector>

namespace mpmc::thermodynamics {

template <std::floating_point T>
    requires std::same_as<T, std::remove_cv_t<T>>
class Pr76Mixture;

/// Private scratch storage, independent of the prepared model and returned values.
/// Reuse one workspace sequentially; concurrent evaluations need different workspaces.
/// No cached thermodynamic values survive a call: every root is recomputed, so AD
/// seeds from a previous Jacobian block or a different composition cannot leak.
template <typename Number>
class Pr76MixtureWorkspace {
public:
    Pr76MixtureWorkspace() = default;
    Pr76MixtureWorkspace(const Pr76MixtureWorkspace&) = delete;
    Pr76MixtureWorkspace& operator=(const Pr76MixtureWorkspace&) = delete;
    Pr76MixtureWorkspace(Pr76MixtureWorkspace&&) = delete;
    Pr76MixtureWorkspace& operator=(Pr76MixtureWorkspace&&) = delete;
    [[nodiscard]] std::size_t capacity() const noexcept { return roots_.capacity(); }
private:
    template <std::floating_point T>
        requires std::same_as<T, std::remove_cv_t<T>>
    friend class Pr76Mixture;
    std::vector<Number> roots_;
};

/// Molar SI coefficients, NOT the dimensionless EOS A and B. Ownership is by value.
template <typename Number>
struct Pr76MixtureValues {
    Number a; // Pa m^6 mol^-2.
    Number b; // m^3 mol^-1.
};

/// Classical PR76 mixing, journal p.60 Eqs.(20)-(22); paper delta_ij = project kij.
/// Immutable, owning preparation from ONE validated parameter snapshot. Includes
/// all selected components, including those with zero mole fraction at evaluation.
/// No EOS roots, phase decisions, normalization or implicit/default binary data.
/// See pr76_mixture.md for coordinate semantics, domains and independent derivatives.
template <std::floating_point T = double>
    requires std::same_as<T, std::remove_cv_t<T>>
class Pr76Mixture {
public:
    Pr76Mixture(const Pr76Mixture&) = default;
    Pr76Mixture(Pr76Mixture&&) noexcept = default;
    Pr76Mixture& operator=(const Pr76Mixture&) = delete;
    Pr76Mixture& operator=(Pr76Mixture&&) = delete;

    [[nodiscard]] static Pr76Mixture from_parameters(const PrParameterSet& parameters,
                                                    ContractLimits limits = {}) {
        const std::size_t n = parameters.components().size();
        if (n == 0 || n > limits.max_components ||
            parameters.binary_records().size() > limits.max_pair_records ||
            n > limits.max_matrix_entries / n || n > std::vector<T>{}.max_size() / n ||
            n > std::vector<Pr76Pure<T>>{}.max_size()) {
            throw std::length_error("Pr76Mixture: invalid size or preparation limit exceeded");
        }
        return Pr76Mixture(parameters);
    }

    [[nodiscard]] std::size_t size() const noexcept { return pure_.size(); }
    [[nodiscard]] const PrParameterSet& parameters() const & noexcept { return parameters_; }
    const PrParameterSet& parameters() const && = delete;

    /// Dimensionless ABSOLUTE acceptance threshold for sum(x)-1. Checking this
    /// tolerance does not modify x, normalize it or differentiate a correction.
    [[nodiscard]] static constexpr T sum_tolerance() noexcept {
        return T{64} * std::numeric_limits<T>::epsilon();
    }

    /// Full n-component vector, in parameters().components() order.
    /// Finite primals in [0,1], abs(sum(x)-1)<=sum_tolerance(), finite AD seeds.
    /// Each entry's seed is preserved INDEPENDENTLY: derivatives are those of the
    /// unconstrained mixing polynomial evaluated at the accepted composition.
    /// Seed sums need not vanish. For simplex-constrained derivatives use reduced
    /// coordinates below, or explicitly supply a tangent seed with zero sum.
    template <typename Number>
        requires detail::Pr76Number<Number, T>
    [[nodiscard]] Pr76MixtureValues<Number> evaluate_full(
        const Number& temperature_k, std::type_identity_t<std::span<const Number>> fractions,
        Pr76MixtureWorkspace<Number>& workspace) const {
        return evaluate_impl<false>(temperature_k, fractions, workspace);
    }

    /// n-1 independent mole fractions in the first n-1 positions of the snapshot.
    /// Construct x_last=1-sum(independent) in Number arithmetic: the dependent
    /// component's derivatives are NOT discarded. For n=1 the input span is empty.
    /// Reject an out-of-[0,1] remainder, even a small negative one; never clip it.
    /// Changing the last component changes the coordinate chart and its Jacobian.
    template <typename Number>
        requires detail::Pr76Number<Number, T>
    [[nodiscard]] Pr76MixtureValues<Number> evaluate_reduced(
        const Number& temperature_k, std::type_identity_t<std::span<const Number>> independent,
        Pr76MixtureWorkspace<Number>& workspace) const {
        return evaluate_impl<true>(temperature_k, independent, workspace);
    }

private:
    template <bool Reduced, typename Number>
    [[nodiscard]] Pr76MixtureValues<Number> evaluate_impl(
        const Number& temperature_k, std::span<const Number> fractions,
        Pr76MixtureWorkspace<Number>& workspace) const {
        const std::size_t n = size();
        if (n == 0 || fractions.size() != n - (Reduced ? 1U : 0U)) {
            throw std::invalid_argument("Pr76Mixture: composition dimension mismatch");
        }
        Number last{T{1}};
        if constexpr (Reduced) {
            Number independent_sum{T{0}};
            for (const auto& fraction : fractions) {
                check_fraction(fraction);
                independent_sum += fraction;
            }
            last -= independent_sum;
            check_fraction(last);
        } else {
            long double sum = 0, compensation = 0;
            for (const auto& fraction : fractions) {
                check_fraction(fraction);
                // Compensated primal-only validation; use ORIGINAL x in the expression.
                const long double increment = detail::pr76_value(fraction) - compensation;
                const long double next = sum + increment;
                compensation = (next - sum) - increment;
                sum = next;
            }
            if (std::abs(sum - 1.0L) > static_cast<long double>(sum_tolerance())) {
                throw std::domain_error("Pr76Mixture: mole fractions are not normalized");
            }
        }
        const auto x = [&](std::size_t i) -> const Number& {
            if constexpr (Reduced) {
                if (i == n - 1) {
                    return last;
                }
            }
            return fractions[i];
        };
        // Only capacity growth allocates; roots are never reused as a value cache.
        workspace.roots_.resize(n);
        Number attraction{T{0}}, covolume{T{0}};
        using std::sqrt;
        for (std::size_t i = 0; i < n; ++i) {
            const auto value = pure_[i].evaluate(temperature_k);
            attraction += (x(i) * value.a) * x(i); // Diagonal: use ai, not sqrt(ai)^2.
            covolume += x(i) * value.b;
            if (n > 1) {
                // At ai=0, sqrt(ai(T)) may have a cusp. This minimal differentiable
                // multi-component kernel rejects it uniformly, even when xi=0 or kij=1.
                if (!(detail::pr76_value(value.a) > T{0})) {
                    throw std::domain_error("Pr76Mixture: multi-component mixing requires ai > 0");
                }
                workspace.roots_[i] = sqrt(value.a); // ADL; principal NONNEGATIVE root.
                if (!detail::pr76_finite(workspace.roots_[i])) {
                    throw std::range_error("Pr76Mixture: nonrepresentable attraction root");
                }
            }
        }
        for (std::size_t i = 0; i < n; ++i) {
            for (std::size_t j = i + 1; j < n; ++j) {
                // Symmetry halves pair work. Never form sqrt(ai*aj), and never skip
                // a zero fraction: its derivative can be nonzero at the boundary.
                const Number cross = ((x(i) * workspace.roots_[i]) *
                                      (x(j) * workspace.roots_[j])) * factors_[i * n + j];
                attraction += cross;
                attraction += cross;
            }
        }
        if (!detail::pr76_finite(attraction) || !detail::pr76_finite(covolume) ||
            !(detail::pr76_value(covolume) > T{0})) {
            throw std::range_error("Pr76Mixture: nonrepresentable coefficient or derivative");
        }
        // Arbitrary finite kij can yield negative/zero a_mix. Preserve the algebraic
        // result, not an unrequested stability judgement or a clipped positive value.
        return {attraction, covolume};
    }

    template <typename Number>
    static void check_fraction(const Number& fraction) {
        const T value = detail::pr76_value(fraction);
        if (!detail::pr76_finite(fraction) || value < T{0} || value > T{1}) {
            throw std::domain_error("Pr76Mixture: fractions/seeds must be finite, x in [0,1]");
        }
    }

    explicit Pr76Mixture(const PrParameterSet& parameters) : parameters_(parameters) {
        const std::size_t n = parameters_.components().size();
        pure_.reserve(n);
        for (std::size_t i = 0; i < n; ++i) {
            pure_.push_back(Pr76Pure<T>::from_parameters(parameters_, i));
        }
        factors_.reserve(n * n); // Size multiplication was checked in the factory.
        for (const double kij : parameters_.kij_matrix()) {
            const long double limit = std::numeric_limits<T>::max();
            if (static_cast<long double>(kij) > limit || static_cast<long double>(kij) < -limit) {
                throw std::range_error("Pr76Mixture: kij outside selected scalar range");
            }
            const T converted = static_cast<T>(kij);
            if (converted == T{0} && kij != 0.0) {
                throw std::range_error("Pr76Mixture: nonzero kij underflows selected scalar");
            }
            const T factor = T{1} - converted;
            if (!std::isfinite(factor)) {
                throw std::range_error("Pr76Mixture: nonrepresentable binary mixing factor");
            }
            factors_.push_back(factor);
        }
    }

    PrParameterSet parameters_;
    std::vector<Pr76Pure<T>> pure_;
    std::vector<T> factors_; // Constant 1-kij, in the SAME snapshot order; no hot-loop ID lookup.
};

} // namespace mpmc::thermodynamics
#endif // MPMC_THERMODYNAMICS_PR76_MIXTURE_HPP
