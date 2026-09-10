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

template <std::floating_point T> requires std::same_as<T, std::remove_cv_t<T>> class Sw92Phase;
template <std::floating_point T> struct Sw92MixtureValues { T a, b; };

template <std::floating_point T>
class Sw92MixtureWorkspace {
public:
    Sw92MixtureWorkspace() = default;
    Sw92MixtureWorkspace(const Sw92MixtureWorkspace&) = delete;
    Sw92MixtureWorkspace& operator=(const Sw92MixtureWorkspace&) = delete;
    [[nodiscard]] std::size_t capacity() const noexcept { return sqrt_ai_.capacity(); }
private:
    template <std::floating_point U> requires std::same_as<U, std::remove_cv_t<U>> friend class Sw92Mixture;
    template <std::floating_point U> requires std::same_as<U, std::remove_cv_t<U>> friend class Sw92Phase;
    std::vector<T> sqrt_ai_;
};

/// Classical PR Eq.(10) with explicit AQ/NA water-pair behavior. No composition
/// normalization, stability search, root selection, or phase labeling.
template <std::floating_point T = double>
    requires std::same_as<T, std::remove_cv_t<T>>
class Sw92Mixture {
public:
    Sw92Mixture(const Sw92Mixture&) = default;
    Sw92Mixture(Sw92Mixture&&) noexcept = default;
    Sw92Mixture& operator=(const Sw92Mixture&) = delete;
    Sw92Mixture& operator=(Sw92Mixture&&) = delete;

    [[nodiscard]] static Sw92Mixture from_parameters(const Sw92ParameterSet& p,
                                                      ContractLimits limits = {}) {
        const std::size_t n = p.components().size();
        if (n == 0 || n > limits.max_components || n > limits.max_matrix_entries/n ||
            n > std::vector<Sw92Pure<T>>{}.max_size())
            throw std::length_error("Sw92Mixture: size/preparation limit exceeded");
        return Sw92Mixture(p);
    }
    [[nodiscard]] std::size_t size() const noexcept { return pure_.size(); }
    [[nodiscard]] const Sw92ParameterSet& parameters() const & noexcept { return parameters_; }
    const Sw92ParameterSet& parameters() const && = delete;
    [[nodiscard]] static constexpr T sum_tolerance() noexcept {
        return T{64}*std::numeric_limits<T>::epsilon();
    }

    [[nodiscard]] T kij(T t, T molality, SwPhaseFamily family,
                        std::size_t i, std::size_t j) const {
        check_state(t, molality); check_family(family);
        if (i >= size() || j >= size()) throw std::out_of_range("Sw92Mixture::kij: index");
        return kij_unchecked(t, molality, family, i, j);
    }
    [[nodiscard]] Sw92MixtureValues<T> evaluate(T t, std::span<const T> x,
                                                T molality, SwPhaseFamily family,
                                                Sw92MixtureWorkspace<T>& w) const {
        return evaluate_impl<false>(t, x, molality, family, w, {});
    }

private:
    template <std::floating_point U> requires std::same_as<U, std::remove_cv_t<U>> friend class Sw92Phase;

    template <bool WithRows>
    [[nodiscard]] Sw92MixtureValues<T> evaluate_impl(T t, std::span<const T> x,
        T molality, SwPhaseFamily family, Sw92MixtureWorkspace<T>& w,
        std::span<T> rows) const {
        check_state(t, molality); check_family(family);
        if (x.size() != size()) throw std::invalid_argument("Sw92Mixture: composition dimension");
        if constexpr (WithRows) if (rows.size() != size())
            throw std::invalid_argument("Sw92Mixture: internal row dimension");
        long double sum = 0, corr = 0;
        for (T xi : x) {
            if (!std::isfinite(xi) || xi < T{0} || xi > T{1})
                throw std::domain_error("Sw92Mixture: fractions must be finite in [0,1]");
            const long double y = static_cast<long double>(xi)-corr;
            const long double next = sum+y; corr=(next-sum)-y; sum=next;
        }
        if (std::abs(sum-1.0L) > static_cast<long double>(sum_tolerance()))
            throw std::domain_error("Sw92Mixture: mole fractions are not normalized");

        const std::size_t n = size(); w.sqrt_ai_.resize(n);
        T a{}, b{};
        for (std::size_t i=0;i<n;++i) {
            const auto v = pure_[i].evaluate(t, molality);
            if (!(v.a > T{0})) throw std::domain_error("Sw92Mixture: positive pure a required");
            w.sqrt_ai_[i]=std::sqrt(v.a);
            a += (x[i]*v.a)*x[i]; b += x[i]*v.b;
            if constexpr (WithRows) rows[i]=x[i]*v.a;
        }
        for (std::size_t i=0;i<n;++i) for (std::size_t j=i+1;j<n;++j) {
            const T aij=(w.sqrt_ai_[i]*(T{1}-kij_unchecked(t,molality,family,i,j)))*w.sqrt_ai_[j];
            const T cross=(x[i]*x[j])*aij; a += cross; a += cross;
            if constexpr (WithRows) { rows[i]+=x[j]*aij; rows[j]+=x[i]*aij; }
        }
        if (!std::isfinite(a) || !std::isfinite(b) || !(b > T{0}))
            throw std::range_error("Sw92Mixture: nonrepresentable mixed coefficient");
        if constexpr (WithRows) for (T r:rows) if (!std::isfinite(r))
            throw std::range_error("Sw92Mixture: nonrepresentable row sum");
        return {a,b};
    }

    explicit Sw92Mixture(const Sw92ParameterSet& p):parameters_(p) {
        pure_.reserve(p.components().size());
        for (std::size_t i=0;i<p.components().size();++i) pure_.push_back(Sw92Pure<T>::from_parameters(p,i));
    }
    void check_state(T t, T m) const { detail::sw92_check_state(t,m,parameters_.applicability()); }
    static void check_family(SwPhaseFamily f) {
        switch(f){case SwPhaseFamily::aqueous:case SwPhaseFamily::nonaqueous:return;}
        throw std::invalid_argument("Sw92Mixture: unknown phase family");
    }
    [[nodiscard]] T kij_unchecked(T t,T m,SwPhaseFamily f,std::size_t i,std::size_t j) const {
        if(i==j) return T{0};
        const std::size_t water=parameters_.water_index();
        if(i!=water && j!=water) return detail::sw92_cast<T>(parameters_.nonwater_kij(f,i,j));
        const std::size_t other=i==water?j:i;
        const T tc=detail::sw92_cast<T>(parameters_.critical_temperatures_k()[other]);
        const T tr=t/tc;
        if(f==SwPhaseFamily::aqueous)
            return sw92_aqueous_water_kij(parameters_.species(other),tr,
                detail::sw92_cast<T>(parameters_.acentric_factors()[other]),m);
        if(parameters_.species(other)==Sw92Species::hydrogen_sulfide) {
            if(parameters_.water_nonaqueous_rule(other)!=Sw92NonAqueousWaterRule::hydrogen_sulfide_eq17)
                throw std::logic_error("SW92: H2S Eq.(17) routing lost");
            return sw92_h2s_nonaqueous_water_kij(tr);
        }
        return detail::sw92_cast<T>(parameters_.water_nonaqueous_constant_kij(other));
    }

    Sw92ParameterSet parameters_;
    std::vector<Sw92Pure<T>> pure_;
};

} // namespace mpmc::thermodynamics

#endif // MPMC_THERMODYNAMICS_SW92_MIXTURE_HPP
