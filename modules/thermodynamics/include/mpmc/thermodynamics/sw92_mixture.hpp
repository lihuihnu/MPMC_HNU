#ifndef MPMC_THERMODYNAMICS_SW92_MIXTURE_HPP
#define MPMC_THERMODYNAMICS_SW92_MIXTURE_HPP

#include <mpmc/thermodynamics/sw92_pure.hpp>

#include <cmath>
#include <concepts>
#include <cstddef>
#include <limits>
#include <span>
#include <stdexcept>
#include <type_traits>
#include <vector>

namespace mpmc::thermodynamics {

template <std::floating_point T>
    requires std::same_as<T, std::remove_cv_t<T>>
class Sw92Phase;

template <typename Number>
struct Sw92MixtureValues { Number a, b; };

template <typename Number>
class Sw92MixtureWorkspace {
public:
    Sw92MixtureWorkspace() = default;
    Sw92MixtureWorkspace(const Sw92MixtureWorkspace&) = delete;
    Sw92MixtureWorkspace& operator=(const Sw92MixtureWorkspace&) = delete;
    [[nodiscard]] std::size_t capacity() const noexcept { return sqrt_ai_.capacity(); }
private:
    template <std::floating_point U>
        requires std::same_as<U, std::remove_cv_t<U>>
    friend class Sw92Mixture;
    template <std::floating_point U>
        requires std::same_as<U, std::remove_cv_t<U>>
    friend class Sw92Phase;
    std::vector<Number> sqrt_ai_;
};

/// Classical PR Eq.(10) with explicit AQ/NA water-pair behavior. No composition
/// normalization, stability search, root selection, or phase labeling. The
/// scalar model is immutable while local Number arithmetic may carry AD seeds.
template <std::floating_point T = double>
    requires std::same_as<T, std::remove_cv_t<T>>
class Sw92Mixture {
public:
    Sw92Mixture(const Sw92Mixture&) = default;
    Sw92Mixture(Sw92Mixture&&) noexcept = default;
    Sw92Mixture& operator=(const Sw92Mixture&) = delete;
    Sw92Mixture& operator=(Sw92Mixture&&) = delete;

    [[nodiscard]] static Sw92Mixture from_parameters(
        const Sw92ParameterSet& p, ContractLimits limits = {}) {
        const std::size_t n = p.components().size();
        if (n == 0 || n > limits.max_components || n > limits.max_matrix_entries/n ||
            n > std::vector<Sw92Pure<T>>{}.max_size())
            throw std::length_error("Sw92Mixture: size/preparation limit exceeded");
        return Sw92Mixture(p);
    }
    [[nodiscard]] std::size_t size() const noexcept { return pure_.size(); }
    [[nodiscard]] const Sw92ParameterSet& parameters() const & noexcept {
        return parameters_;
    }
    const Sw92ParameterSet& parameters() const && = delete;
    [[nodiscard]] static constexpr T sum_tolerance() noexcept {
        return T{64}*std::numeric_limits<T>::epsilon();
    }

    [[nodiscard]] T kij(T t, T molality, SwPhaseFamily family,
                        std::size_t i, std::size_t j) const {
        check_state(t, molality);
        check_family(family);
        if (i >= size() || j >= size())
            throw std::out_of_range("Sw92Mixture::kij: index");
        return detail::sw92_value(
            kij_unchecked<T>(t, molality, family, i, j));
    }

    // Backward-compatible full-composition floating-point entry point.
    [[nodiscard]] Sw92MixtureValues<T> evaluate(
        T t, std::span<const T> x, T molality, SwPhaseFamily family,
        Sw92MixtureWorkspace<T>& workspace) const {
        return evaluate_impl<false, false>(
            t, x, molality, family, workspace, {});
    }

    template <typename Number>
        requires detail::Sw92Number<Number, T>
    [[nodiscard]] Sw92MixtureValues<Number> evaluate_full(
        const Number& t, std::type_identity_t<std::span<const Number>> x,
        T molality, SwPhaseFamily family,
        Sw92MixtureWorkspace<Number>& workspace) const {
        return evaluate_impl<false, false>(
            t, x, molality, family, workspace, {});
    }

    template <typename Number>
        requires detail::Sw92Number<Number, T>
    [[nodiscard]] Sw92MixtureValues<Number> evaluate_reduced(
        const Number& t,
        std::type_identity_t<std::span<const Number>> independent,
        T molality, SwPhaseFamily family,
        Sw92MixtureWorkspace<Number>& workspace) const {
        return evaluate_impl<true, false>(
            t, independent, molality, family, workspace, {});
    }

private:
    template <std::floating_point U>
        requires std::same_as<U, std::remove_cv_t<U>>
    friend class Sw92Phase;

    template <bool Reduced, bool WithRows, typename Number>
        requires detail::Sw92Number<Number, T>
    [[nodiscard]] Sw92MixtureValues<Number> evaluate_impl(
        const Number& t, std::span<const Number> fractions,
        T molality, SwPhaseFamily family,
        Sw92MixtureWorkspace<Number>& workspace,
        std::span<Number> rows) const {
        check_state(t, molality);
        check_family(family);
        const std::size_t n = size();
        if (n == 0 || fractions.size() != n - (Reduced ? 1U : 0U))
            throw std::invalid_argument("Sw92Mixture: composition dimension");
        if constexpr (WithRows) {
            if (rows.size() != n)
                throw std::invalid_argument("Sw92Mixture: internal row dimension");
        }

        Number last{T{1}};
        if constexpr (Reduced) {
            Number sum{T{0}};
            for (const auto& xi : fractions) {
                check_fraction(xi);
                sum += xi;
            }
            last -= sum;
            check_fraction(last);
        } else {
            long double sum = 0.0L;
            long double correction = 0.0L;
            for (const auto& xi : fractions) {
                check_fraction(xi);
                const long double increment =
                    static_cast<long double>(detail::sw92_value(xi)) - correction;
                const long double next = sum + increment;
                correction = (next - sum) - increment;
                sum = next;
            }
            if (std::abs(sum - 1.0L) >
                static_cast<long double>(sum_tolerance()))
                throw std::domain_error("Sw92Mixture: mole fractions are not normalized");
        }
        const auto x = [&](std::size_t i) -> const Number& {
            if constexpr (Reduced) {
                if (i == n - 1U) return last;
            }
            return fractions[i];
        };

        workspace.sqrt_ai_.resize(n);
        Number attraction{T{0}};
        Number covolume{T{0}};
        using std::sqrt;
        for (std::size_t i = 0; i < n; ++i) {
            const auto value = pure_[i].evaluate(t, molality);
            if (!(detail::sw92_value(value.a) > T{0}))
                throw std::domain_error("Sw92Mixture: positive pure a required");
            workspace.sqrt_ai_[i] = sqrt(value.a);
            attraction += (x(i) * value.a) * x(i);
            covolume += x(i) * value.b;
            if constexpr (WithRows) rows[i] = x(i) * value.a;
        }
        for (std::size_t i = 0; i < n; ++i) {
            for (std::size_t j = i + 1; j < n; ++j) {
                const Number aij =
                    (workspace.sqrt_ai_[i] *
                     (T{1} - kij_unchecked<Number>(t, molality, family, i, j))) *
                    workspace.sqrt_ai_[j];
                const Number cross = (x(i) * x(j)) * aij;
                attraction += cross;
                attraction += cross;
                if constexpr (WithRows) {
                    rows[i] += x(j) * aij;
                    rows[j] += x(i) * aij;
                }
            }
        }
        if (!detail::sw92_finite(attraction) ||
            !detail::sw92_finite(covolume) ||
            !(detail::sw92_value(covolume) > T{0}))
            throw std::range_error(
                "Sw92Mixture: nonrepresentable mixed coefficient or derivative");
        if constexpr (WithRows) {
            for (const auto& row : rows) {
                if (!detail::sw92_finite(row))
                    throw std::range_error(
                        "Sw92Mixture: nonrepresentable row sum or derivative");
            }
        }
        return {attraction, covolume};
    }

    explicit Sw92Mixture(const Sw92ParameterSet& p) : parameters_(p) {
        pure_.reserve(p.components().size());
        for (std::size_t i = 0; i < p.components().size(); ++i)
            pure_.push_back(Sw92Pure<T>::from_parameters(p, i));
    }

    template <typename Number>
        requires detail::Sw92Number<Number, T>
    void check_state(const Number& t, T molality) const {
        detail::sw92_check_state(t, molality, parameters_.applicability());
    }

    template <typename Number>
        requires detail::Sw92Number<Number, T>
    static void check_fraction(const Number& xi) {
        const T value = detail::sw92_value(xi);
        if (!detail::sw92_finite(xi) || value < T{0} || value > T{1})
            throw std::domain_error(
                "Sw92Mixture: fractions/seeds must be finite with primal in [0,1]");
    }

    static void check_family(SwPhaseFamily family) {
        switch (family) {
        case SwPhaseFamily::aqueous:
        case SwPhaseFamily::nonaqueous:
            return;
        }
        throw std::invalid_argument("Sw92Mixture: unknown phase family");
    }

    template <typename Number>
        requires detail::Sw92Number<Number, T>
    [[nodiscard]] Number kij_unchecked(
        const Number& t, T molality, SwPhaseFamily family,
        std::size_t i, std::size_t j) const {
        if (i == j) return Number{T{0}};
        const std::size_t water = parameters_.water_index();
        if (i != water && j != water)
            return Number{detail::sw92_cast<T>(
                parameters_.nonwater_kij(family, i, j))};
        const std::size_t other = i == water ? j : i;
        const T tc = detail::sw92_cast<T>(
            parameters_.critical_temperatures_k()[other]);
        const Number tr = t / tc;
        if (family == SwPhaseFamily::aqueous)
            return sw92_aqueous_water_kij<Number, T>(
                parameters_.species(other), tr,
                detail::sw92_cast<T>(parameters_.acentric_factors()[other]),
                molality);
        if (parameters_.species(other) == Sw92Species::hydrogen_sulfide) {
            if (parameters_.water_nonaqueous_rule(other) !=
                Sw92NonAqueousWaterRule::hydrogen_sulfide_eq17)
                throw std::logic_error("SW92: H2S Eq.(17) routing lost");
            return sw92_h2s_nonaqueous_water_kij<Number, T>(tr);
        }
        return Number{detail::sw92_cast<T>(
            parameters_.water_nonaqueous_constant_kij(other))};
    }

    Sw92ParameterSet parameters_;
    std::vector<Sw92Pure<T>> pure_;
};

} // namespace mpmc::thermodynamics

#endif // MPMC_THERMODYNAMICS_SW92_MIXTURE_HPP
