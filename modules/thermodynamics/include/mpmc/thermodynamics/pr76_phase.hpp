#ifndef MPMC_THERMODYNAMICS_PR76_PHASE_HPP
#define MPMC_THERMODYNAMICS_PR76_PHASE_HPP

#include <mpmc/thermodynamics/pr76_mixture.hpp>
#include <mpmc/thermodynamics/pr76_roots.hpp>

#include <cmath>
#include <concepts>
#include <cstddef>
#include <limits>
#include <span>
#include <stdexcept>
#include <string_view>
#include <type_traits>
#include <utility>
#include <vector>

namespace mpmc::thermodynamics {

/// Eq.(4)-consistent radical factors, not the rounded 2.414/0.414 printed in Eq.(19).
inline constexpr std::string_view pr76_pt_convention =
    "PR76/printed-coefficients/R-SI-2019/exact-sqrt2/PT-v1";

enum class Pr76PhaseErrorCode {
    near_multiple, iteration_limit, unrepresentable_root, ill_conditioned_derivative
};
class Pr76PhaseError : public std::runtime_error {
public:
    Pr76PhaseError(Pr76PhaseErrorCode code, const char* message)
        : std::runtime_error(message), code_(code) {}
    [[nodiscard]] Pr76PhaseErrorCode code() const noexcept { return code_; }
private:
    Pr76PhaseErrorCode code_;
};

/// Scratch only. No cached values/seeds; one workspace per simultaneous call.
template <typename Number>
class Pr76PhaseWorkspace {
public:
    Pr76PhaseWorkspace() = default;
    Pr76PhaseWorkspace(const Pr76PhaseWorkspace&) = delete;
    Pr76PhaseWorkspace& operator=(const Pr76PhaseWorkspace&) = delete;
    Pr76PhaseWorkspace(Pr76PhaseWorkspace&&) = delete;
    Pr76PhaseWorkspace& operator=(Pr76PhaseWorkspace&&) = delete;
    [[nodiscard]] std::size_t capacity() const noexcept { return sums_.capacity(); }
private:
    template <std::floating_point U>
        requires std::same_as<U, std::remove_cv_t<U>>
    friend class Pr76Phase;
    Pr76MixtureWorkspace<Number> mixture_;
    std::vector<Number> sums_;
};

template <typename Number, std::floating_point T>
struct Pr76PhaseValues {
    Number z;
    std::vector<Number> ln_phi; // Dimensionless, in snapshot order; owns its storage.
    Pr76RootSet<T> root_set;
    std::size_t root_index;
};

namespace detail {
// Analytic continuation of log1p(u)/u. At |u|<=1/8, degree N=digits/3+2
// bounds the omitted value by |u|^(N+1)/[(N+2)(1-|u|)] and the omitted slope
// by |u|^N/(1-|u|)^2, both below epsilon. The outer quotient is then safely
// separated from zero. Coefficients are mathematical series terms, not fitted data.
template <std::floating_point T, typename Number>
[[nodiscard]] Number pr76_log1p_over_x(const Number& u) {
    if (std::abs(pr76_value(u)) <= T{1}/T{8}) {
        constexpr int degree = std::numeric_limits<T>::digits / 3 + 2;
        Number value{(degree % 2 == 0 ? T{1} : T{-1}) / static_cast<T>(degree + 1)};
        for (int k = degree - 1; k >= 0; --k) {
            value = (k % 2 == 0 ? T{1} : T{-1}) / static_cast<T>(k + 1) + u * value;
        }
        return value;
    }
    using std::log1p;
    return log1p(u) / u;
}

inline void pr76_require_resolved(Pr76RootStatus status) {
    switch (status) {
    case Pr76RootStatus::success: return;
    case Pr76RootStatus::near_multiple:
        throw Pr76PhaseError(Pr76PhaseErrorCode::near_multiple,
                            "Pr76Phase: unresolved near-multiple root topology");
    case Pr76RootStatus::iteration_limit:
        throw Pr76PhaseError(Pr76PhaseErrorCode::iteration_limit,
                            "Pr76Phase: root iteration limit exceeded");
    case Pr76RootStatus::unrepresentable:
        throw Pr76PhaseError(Pr76PhaseErrorCode::unrepresentable_root,
                            "Pr76Phase: root/domain/residual not representable");
    }
    throw std::invalid_argument("Pr76Phase: unknown root status");
}
} // namespace detail

/// PT candidate-phase properties for PR76; no flash, phase labels or TPD test.
/// Enumerate roots first, then explicitly select an index in INCREASING Z order.
/// Even the middle (mechanically unstable) candidate can be inspected; none is
/// silently discarded or called a stable equilibrium phase. Composition and AD
/// coordinate semantics are exactly those of Pr76Mixture. Declared T/p
/// applicability is advisory: finite positive states may be evaluated outside it,
/// while parameters().applicability().assess(T,p) carries inside/outside/unknown
/// information for a caller to publish as a model-validity warning.
template <std::floating_point T = double>
    requires std::same_as<T, std::remove_cv_t<T>>
class Pr76Phase {
public:
    Pr76Phase(const Pr76Phase&) = default;
    Pr76Phase(Pr76Phase&&) noexcept = default;
    Pr76Phase& operator=(const Pr76Phase&) = delete;
    Pr76Phase& operator=(Pr76Phase&&) = delete;

    [[nodiscard]] static Pr76Phase from_parameters(const PrParameterSet& parameters,
                                                   ContractLimits limits = {}) {
        return Pr76Phase(parameters, limits);
    }
    [[nodiscard]] std::size_t size() const noexcept { return mixture_.size(); }
    [[nodiscard]] const PrParameterSet& parameters() const & noexcept {
        return mixture_.parameters();
    }
    const PrParameterSet& parameters() const && = delete;

    [[nodiscard]] Pr76RootSet<T> roots_full(
        T pressure_pa, T temperature_k, std::span<const T> fractions,
        Pr76PhaseWorkspace<T>& workspace, Pr76RootOptions options = {}) const {
        check_pressure(pressure_pa);
        const auto mixed = mixture_.evaluate_full(temperature_k, fractions, workspace.mixture_);
        const auto ab = coefficients(pressure_pa, temperature_k, mixed);
        return pr76_roots(ab.a, ab.b, options);
    }
    [[nodiscard]] Pr76RootSet<T> roots_reduced(
        T pressure_pa, T temperature_k, std::span<const T> fractions,
        Pr76PhaseWorkspace<T>& workspace, Pr76RootOptions options = {}) const {
        check_pressure(pressure_pa);
        const auto mixed = mixture_.evaluate_reduced(temperature_k, fractions, workspace.mixture_);
        const auto ab = coefficients(pressure_pa, temperature_k, mixed);
        return pr76_roots(ab.a, ab.b, options);
    }

    template <typename Number>
        requires detail::Pr76PhaseNumber<Number, T>
    [[nodiscard]] Pr76PhaseValues<Number, T> evaluate_full(
        const Number& pressure_pa, const Number& temperature_k,
        std::type_identity_t<std::span<const Number>> fractions, std::size_t root_index,
        Pr76PhaseWorkspace<Number>& workspace, Pr76RootOptions options = {}) const {
        return evaluate_impl<false>(pressure_pa, temperature_k, fractions, root_index, workspace, options);
    }
    template <typename Number>
        requires detail::Pr76PhaseNumber<Number, T>
    [[nodiscard]] Pr76PhaseValues<Number, T> evaluate_reduced(
        const Number& pressure_pa, const Number& temperature_k,
        std::type_identity_t<std::span<const Number>> fractions, std::size_t root_index,
        Pr76PhaseWorkspace<Number>& workspace, Pr76RootOptions options = {}) const {
        return evaluate_impl<true>(pressure_pa, temperature_k, fractions, root_index, workspace, options);
    }

private:
    Pr76Phase(const PrParameterSet& parameters, ContractLimits limits)
        : mixture_(Pr76Mixture<T>::from_parameters(parameters, limits)) {}

    template <typename Number>
    void check_pressure(const Number& pressure) const {
        if (size() == 0) { throw std::invalid_argument("Pr76Phase: moved-from model"); }
        const T p = detail::pr76_value(pressure);
        if (!detail::pr76_finite(pressure) || !(p > T{0})) {
            throw std::domain_error("Pr76Phase: finite pressure/seeds and p>0 Pa required");
        }
        // Dataset/model applicability bounds are advisory metadata. Do not turn a
        // finite positive pressure into a numerical-domain error merely because it
        // is outside a declared validation interval; callers can query assess(T,p).
    }

    template <typename Number>
    [[nodiscard]] static Pr76MixtureValues<Number> coefficients(
        const Number& pressure, const Number& temperature, const Pr76MixtureValues<Number>& mixed) {
        const Number rt = Pr76Pure<T>::gas_constant() * temperature;
        if (!detail::pr76_finite(rt) || !(detail::pr76_value(rt) > T{0})) {
            throw std::range_error("Pr76Phase: nonrepresentable RT");
        }
        const Number p_over_rt = pressure / rt;
        const Number a = (mixed.a / rt) * p_over_rt;
        const Number b = mixed.b * p_over_rt;
        if (!detail::pr76_finite(a) || !detail::pr76_finite(b) ||
            !(detail::pr76_value(b) > T{0}) ||
            (detail::pr76_value(a) == T{0} && detail::pr76_value(mixed.a) != T{0})) {
            throw std::range_error("Pr76Phase: nonrepresentable A/B or derivative");
        }
        return {a, b};
    }

    template <bool Reduced, typename Number>
    [[nodiscard]] Pr76PhaseValues<Number, T> evaluate_impl(
        const Number& pressure, const Number& temperature, std::span<const Number> fractions,
        std::size_t root_index, Pr76PhaseWorkspace<Number>& workspace, Pr76RootOptions options) const {
        check_pressure(pressure);
        workspace.sums_.resize(size());
        // Reuse validated coordinates, pure kernels, principal roots and binary
        // factors in ONE mixing traversal. Output scratch cannot alias input.
        const auto mixed = mixture_.template evaluate_impl<Reduced, true>(
            temperature, fractions, workspace.mixture_, std::span<Number>{workspace.sums_});
        const auto ab = coefficients(pressure, temperature, mixed);
        const auto roots = pr76_roots(detail::pr76_value(ab.a), detail::pr76_value(ab.b), options);
        detail::pr76_require_resolved(roots.status);
        if (root_index >= roots.count) { throw std::out_of_range("Pr76Phase: invalid root index"); }
        const auto& root = roots.roots[root_index];
        Number free_volume{root.free_volume_z};
        if constexpr (!std::floating_point<Number>) {
            if (!root.derivative_valid) {
                throw Pr76PhaseError(Pr76PhaseErrorCode::ill_conditioned_derivative,
                                    "Pr76Phase: selected root has no reliable local derivative");
            }
            const T y0 = root.free_volume_z / roots.scale;
            const auto c = detail::pr76_cubic_coefficients(ab.a, ab.b, roots.scale);
            const Number residual = ((Number{y0} + c.c2) * y0 + c.c1) * y0 + c.c0;
            const T dh = (T{3} * y0 + T{2} * detail::pr76_value(c.c2)) * y0 +
                         detail::pr76_value(c.c1);
            // First-order implicit function theorem, NOT differentiation of bisection.
            // Remove only the residual's primal: returned root values stay unchanged.
            free_volume -= ((residual - Number{detail::pr76_value(residual)}) / dh) * roots.scale;
        }
        const Number z = free_volume + ab.b;
        const T sqrt2 = std::sqrt(T{2});
        const Number denominator = free_volume + (T{2} - sqrt2) * ab.b;
        const Number u = (T{2} * sqrt2) * (ab.b / denominator);
        const Number log_factor = detail::pr76_log1p_over_x<T>(u) / denominator;
        using std::log;
        using std::log1p;
        Number z_minus_one = z - T{1};
        Number log_free = log(free_volume);
        if (std::abs(detail::pr76_value(z_minus_one)) < T{1}/T{4} &&
            detail::pr76_value(ab.b) < T{1}/T{8}) {
            // EOS identity Z-1 = B/(Z-B) - A*Z/(Z^2+2BZ-B^2).
            // Preserve the O(p) terms even when Z rounds to one at low pressure.
            const Number bz = ab.b / z;
            z_minus_one = ab.b / free_volume - (ab.a / z) /
                          (T{1} + T{2} * bz - bz * bz);
            log_free = log1p(z_minus_one - ab.b);
        }
        const Number rt = Pr76Pure<T>::gas_constant() * temperature;
        const Number p_over_rt = pressure / rt;
        std::vector<Number> ln_phi(size());
        for (std::size_t i = 0; i < size(); ++i) {
            const Number ratio = Number{mixture_.pure_[i].covolume()} / mixed.b;
            // Cancel a and B symbolically from Eq.(19). a=0 needs no special
            // division or regularization, and w_i=0 still has finite ln(phi_i).
            const Number coefficient = ((T{2} * workspace.sums_[i] - mixed.a * ratio) / rt) * p_over_rt;
            ln_phi[i] = ratio * z_minus_one - log_free - coefficient * log_factor;
            if (!detail::pr76_finite(ln_phi[i])) {
                throw std::range_error("Pr76Phase: nonrepresentable ln(phi) or derivative");
            }
        }
        if (!detail::pr76_finite(z)) {
            throw std::range_error("Pr76Phase: nonrepresentable root derivative");
        }
        return {z, std::move(ln_phi), roots, root_index};
    }

    Pr76Mixture<T> mixture_;
};

} // namespace mpmc::thermodynamics
#endif // MPMC_THERMODYNAMICS_PR76_PHASE_HPP
