#ifndef MPMC_THERMODYNAMICS_SW92_PHASE_HPP
#define MPMC_THERMODYNAMICS_SW92_PHASE_HPP

#include <mpmc/thermodynamics/pr76_roots.hpp>
#include <mpmc/thermodynamics/sw92_mixture.hpp>

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

/// SW92 retains the Peng-Robinson cubic after replacing water alpha and
/// phase-family BIPs. Family selection and algebraic root selection are separate.
inline constexpr std::string_view sw92_pt_convention =
    "SW92/corrected-original/PR76-base/R-SI-2019/exact-sqrt2/PT-v1";

using Sw92RootStatus = Pr76RootStatus;
using Sw92RootOptions = Pr76RootOptions;
template <std::floating_point T> using Sw92Root = Pr76Root<T>;
template <std::floating_point T> using Sw92RootSet = Pr76RootSet<T>;

enum class Sw92PhaseErrorCode { near_multiple, iteration_limit, unrepresentable_root };
class Sw92PhaseError : public std::runtime_error {
public:
    Sw92PhaseError(Sw92PhaseErrorCode code, const char* message)
        : std::runtime_error(message), code_(code) {}
    [[nodiscard]] Sw92PhaseErrorCode code() const noexcept { return code_; }
private:
    Sw92PhaseErrorCode code_;
};

template <typename Number>
class Sw92PhaseWorkspace {
public:
    Sw92PhaseWorkspace() = default;
    Sw92PhaseWorkspace(const Sw92PhaseWorkspace&) = delete;
    Sw92PhaseWorkspace& operator=(const Sw92PhaseWorkspace&) = delete;
    [[nodiscard]] std::size_t capacity() const noexcept { return sums_.capacity(); }
private:
    template <std::floating_point U>
        requires std::same_as<U, std::remove_cv_t<U>>
    friend class Sw92Phase;
    Sw92MixtureWorkspace<Number> mixture_;
    std::vector<Number> sums_;
};

namespace detail {
template <typename Number, typename = void>
struct Sw92PhaseScalar { using type = Number; };
template <typename Number>
struct Sw92PhaseScalar<Number, std::void_t<typename Number::Scalar>> {
    using type = typename Number::Scalar;
};
template <typename Number>
using Sw92PhaseScalarT = typename Sw92PhaseScalar<Number>::type;
} // namespace detail

// The second template argument defaults from Number so the historical public
// spelling Sw92PhaseValues<double> remains source-compatible. AD callers may
// use Sw92PhaseValues<Dual<...>> while the stored algebraic-root set continues
// to use the underlying floating scalar.
template <typename Number,
          std::floating_point T = detail::Sw92PhaseScalarT<Number>>
struct Sw92PhaseValues {
    Number z{};
    std::vector<Number> ln_phi;
    Sw92RootSet<T> root_set;
    std::size_t root_index{};
    SwPhaseFamily family{SwPhaseFamily::aqueous};
};

namespace detail {
template <std::floating_point T, typename Number>
[[nodiscard]] Number sw92_log1p_over_x(const Number& u) {
    if (!(sw92_value(u) > T{-1}) || !sw92_finite(u))
        throw std::range_error("Sw92Phase: invalid log1p ratio");
    if (std::abs(sw92_value(u)) <= T{1}/T{8}) {
        constexpr int degree = std::numeric_limits<T>::digits/3 + 2;
        Number value{(degree%2==0?T{1}:T{-1})/static_cast<T>(degree+1)};
        for (int k=degree-1;k>=0;--k)
            value=(k%2==0?T{1}:T{-1})/static_cast<T>(k+1)+u*value;
        return value;
    }
    using std::log1p;
    return log1p(u)/u;
}

inline void sw92_require_resolved(Sw92RootStatus status) {
    switch(status) {
    case Sw92RootStatus::success: return;
    case Sw92RootStatus::near_multiple:
        throw Sw92PhaseError(
            Sw92PhaseErrorCode::near_multiple,
            "Sw92Phase: near-multiple PR roots");
    case Sw92RootStatus::iteration_limit:
        throw Sw92PhaseError(
            Sw92PhaseErrorCode::iteration_limit,
            "Sw92Phase: PR root iteration limit");
    case Sw92RootStatus::unrepresentable:
        throw Sw92PhaseError(
            Sw92PhaseErrorCode::unrepresentable_root,
            "Sw92Phase: unrepresentable PR root/domain");
    }
    throw std::invalid_argument("Sw92Phase: unknown root status");
}
} // namespace detail

/// Explicit-family PT candidate kernel. `family` fixes AQ/NA parameterization;
/// `root_index` then selects one increasing-Z algebraic branch. The selected
/// simple root has the same local implicit-derivative treatment as PR76; root
/// iteration itself is never differentiated.
template <std::floating_point T = double>
    requires std::same_as<T, std::remove_cv_t<T>>
class Sw92Phase {
public:
    Sw92Phase(const Sw92Phase&) = default;
    Sw92Phase(Sw92Phase&&) noexcept = default;
    Sw92Phase& operator=(const Sw92Phase&) = delete;
    Sw92Phase& operator=(Sw92Phase&&) = delete;

    [[nodiscard]] static Sw92Phase from_parameters(
        const Sw92ParameterSet& p, ContractLimits limits = {}) {
        return Sw92Phase(p,limits);
    }
    [[nodiscard]] std::size_t size() const noexcept { return mixture_.size(); }
    [[nodiscard]] const Sw92ParameterSet& parameters() const & noexcept {
        return mixture_.parameters();
    }
    const Sw92ParameterSet& parameters() const && = delete;

    [[nodiscard]] Sw92RootSet<T> roots(
        T pressure_pa, T temperature_k, std::span<const T> x,
        T molality, SwPhaseFamily family,
        Sw92PhaseWorkspace<T>& workspace, Sw92RootOptions options={}) const {
        check_pressure(pressure_pa);
        const auto mixed = mixture_.evaluate_full(
            temperature_k, x, molality, family, workspace.mixture_);
        const auto ab = coefficients(pressure_pa, temperature_k, mixed);
        return pr76_roots(ab.a, ab.b, options);
    }

    // Backward-compatible full-composition floating-point API.
    [[nodiscard]] Sw92PhaseValues<T> evaluate(
        T pressure_pa, T temperature_k, std::span<const T> x,
        T molality, SwPhaseFamily family, std::size_t root_index,
        Sw92PhaseWorkspace<T>& workspace, Sw92RootOptions options={}) const {
        return evaluate_impl<false>(
            pressure_pa, temperature_k, x, molality, family, root_index,
            workspace, options);
    }

    template <typename Number>
        requires detail::Sw92Number<Number, T>
    [[nodiscard]] Sw92PhaseValues<Number> evaluate_full(
        const Number& pressure_pa, const Number& temperature_k,
        std::type_identity_t<std::span<const Number>> x,
        T molality, SwPhaseFamily family, std::size_t root_index,
        Sw92PhaseWorkspace<Number>& workspace,
        Sw92RootOptions options={}) const {
        return evaluate_impl<false>(
            pressure_pa, temperature_k, x, molality, family, root_index,
            workspace, options);
    }

    template <typename Number>
        requires detail::Sw92Number<Number, T>
    [[nodiscard]] Sw92PhaseValues<Number> evaluate_reduced(
        const Number& pressure_pa, const Number& temperature_k,
        std::type_identity_t<std::span<const Number>> independent,
        T molality, SwPhaseFamily family, std::size_t root_index,
        Sw92PhaseWorkspace<Number>& workspace,
        Sw92RootOptions options={}) const {
        return evaluate_impl<true>(
            pressure_pa, temperature_k, independent, molality, family,
            root_index, workspace, options);
    }

private:
    Sw92Phase(const Sw92ParameterSet& p, ContractLimits limits)
        : mixture_(Sw92Mixture<T>::from_parameters(p,limits)) {
        covolumes_.reserve(p.components().size());
        for(std::size_t i=0;i<p.components().size();++i)
            covolumes_.push_back(Sw92Pure<T>::from_parameters(p,i).covolume());
    }

    template <typename Number>
        requires detail::Sw92Number<Number, T>
    void check_pressure(const Number& p) const {
        const T pressure = detail::sw92_value(p);
        if(!detail::sw92_finite(p) || !(pressure>T{0}))
            throw std::domain_error(
                "Sw92Phase: finite pressure/seeds and pressure>0 Pa required");
        if(const auto& b=parameters().applicability().state.pressure_pa; b &&
           (static_cast<long double>(pressure)<b->lower ||
            static_cast<long double>(pressure)>b->upper))
            throw std::domain_error(
                "Sw92Phase: pressure outside declared dataset interval");
    }

    template <typename Number>
        requires detail::Sw92Number<Number, T>
    [[nodiscard]] static Sw92MixtureValues<Number> coefficients(
        const Number& p, const Number& t,
        const Sw92MixtureValues<Number>& mixed) {
        if(!detail::sw92_finite(t) || !(detail::sw92_value(t)>T{0}))
            throw std::domain_error(
                "Sw92Phase: finite temperature/seeds and T>0 K required");
        const Number rt=Sw92Pure<T>::gas_constant()*t;
        const Number q=p/rt;
        const Number a=(mixed.a/rt)*q;
        const Number b=mixed.b*q;
        if(!detail::sw92_finite(a)||!detail::sw92_finite(b)||
           !(detail::sw92_value(b)>T{0})||
           (detail::sw92_value(a)==T{0} && detail::sw92_value(mixed.a)!=T{0}))
            throw std::range_error(
                "Sw92Phase: nonrepresentable PR A/B or derivative");
        return {a,b};
    }

    template <bool Reduced, typename Number>
        requires detail::Sw92Number<Number, T>
    [[nodiscard]] Sw92PhaseValues<Number> evaluate_impl(
        const Number& pressure_pa, const Number& temperature_k,
        std::span<const Number> x, T molality, SwPhaseFamily family,
        std::size_t root_index, Sw92PhaseWorkspace<Number>& workspace,
        Sw92RootOptions options) const {
        check_pressure(pressure_pa);
        workspace.sums_.resize(size());
        const auto mixed = mixture_.template evaluate_impl<Reduced, true>(
            temperature_k, x, molality, family, workspace.mixture_,
            std::span<Number>{workspace.sums_});
        const auto ab=coefficients(pressure_pa,temperature_k,mixed);
        const auto roots=pr76_roots(
            detail::sw92_value(ab.a), detail::sw92_value(ab.b), options);
        detail::sw92_require_resolved(roots.status);
        if(root_index>=roots.count)
            throw std::out_of_range("Sw92Phase: invalid root index");
        const auto& root=roots.roots[root_index];

        Number free_volume{root.free_volume_z};
        if constexpr (!std::floating_point<Number>) {
            if (!root.derivative_valid) {
                throw std::range_error(
                    "Sw92Phase: selected root has no reliable local derivative");
            }
            const T y0 = root.free_volume_z / roots.scale;
            const auto c = detail::pr76_cubic_coefficients(
                ab.a, ab.b, roots.scale);
            const Number residual =
                ((Number{y0}+c.c2)*y0+c.c1)*y0+c.c0;
            const T dh =
                (T{3}*y0+T{2}*detail::sw92_value(c.c2))*y0+
                detail::sw92_value(c.c1);
            free_volume -=
                ((residual-Number{detail::sw92_value(residual)})/dh)*roots.scale;
        }

        const Number z=free_volume+ab.b;
        const T sqrt2=std::sqrt(T{2});
        const Number denominator=free_volume+(T{2}-sqrt2)*ab.b;
        if(!detail::sw92_finite(denominator)||
           !(detail::sw92_value(denominator)>T{0}))
            throw std::range_error("Sw92Phase: PR log denominator");
        const Number u=(T{2}*sqrt2)*(ab.b/denominator);
        const Number log_factor=
            detail::sw92_log1p_over_x<T>(u)/denominator;
        if(!detail::sw92_finite(log_factor))
            throw std::range_error("Sw92Phase: PR log factor");

        using std::log;
        using std::log1p;
        Number zm1=z-T{1};
        Number log_free=log(free_volume);
        if(std::abs(detail::sw92_value(zm1))<T{1}/T{4} &&
           detail::sw92_value(ab.b)<T{1}/T{8}) {
            const Number bz=ab.b/z;
            zm1=ab.b/free_volume-(ab.a/z)/
                (T{1}+T{2}*bz-bz*bz);
            log_free=log1p(zm1-ab.b);
        }
        if(!detail::sw92_finite(zm1)||!detail::sw92_finite(log_free))
            throw std::range_error("Sw92Phase: PR free-volume log");

        const Number rt=Sw92Pure<T>::gas_constant()*temperature_k;
        const Number p_over_rt=pressure_pa/rt;
        std::vector<Number> lnphi(size());
        for(std::size_t i=0;i<size();++i) {
            const Number ratio=Number{covolumes_[i]}/mixed.b;
            const Number coefficient=
                ((T{2}*workspace.sums_[i]-mixed.a*ratio)/rt)*p_over_rt;
            lnphi[i]=ratio*zm1-log_free-coefficient*log_factor;
            if(!detail::sw92_finite(lnphi[i]))
                throw std::range_error(
                    "Sw92Phase: nonrepresentable ln(phi) or derivative");
        }
        if(!detail::sw92_finite(z))
            throw std::range_error(
                "Sw92Phase: nonrepresentable root derivative");
        return {z,std::move(lnphi),roots,root_index,family};
    }

    Sw92Mixture<T> mixture_;
    std::vector<T> covolumes_;
};

} // namespace mpmc::thermodynamics

#endif // MPMC_THERMODYNAMICS_SW92_PHASE_HPP
