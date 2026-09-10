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
private: Sw92PhaseErrorCode code_;
};

template <std::floating_point T>
class Sw92PhaseWorkspace {
public:
    Sw92PhaseWorkspace() = default;
    Sw92PhaseWorkspace(const Sw92PhaseWorkspace&) = delete;
    Sw92PhaseWorkspace& operator=(const Sw92PhaseWorkspace&) = delete;
    [[nodiscard]] std::size_t capacity() const noexcept { return sums_.capacity(); }
private:
    template <std::floating_point U> requires std::same_as<U, std::remove_cv_t<U>> friend class Sw92Phase;
    Sw92MixtureWorkspace<T> mixture_;
    std::vector<T> sums_;
};

template <std::floating_point T>
struct Sw92PhaseValues {
    T z{};
    std::vector<T> ln_phi;
    Sw92RootSet<T> root_set;
    std::size_t root_index{};
    SwPhaseFamily family{SwPhaseFamily::aqueous};
};

namespace detail {
template <std::floating_point T>
[[nodiscard]] T sw92_log1p_over_x(T u) {
    if (!std::isfinite(u) || !(u > T{-1}))
        throw std::range_error("Sw92Phase: invalid log1p ratio");
    if (std::abs(u) <= T{1}/T{8}) {
        constexpr int degree = std::numeric_limits<T>::digits/3 + 2;
        T value = (degree%2==0?T{1}:T{-1})/static_cast<T>(degree+1);
        for (int k=degree-1;k>=0;--k)
            value=(k%2==0?T{1}:T{-1})/static_cast<T>(k+1)+u*value;
        return value;
    }
    return std::log1p(u)/u;
}
inline void sw92_require_resolved(Sw92RootStatus status) {
    switch(status) {
    case Sw92RootStatus::success: return;
    case Sw92RootStatus::near_multiple:
        throw Sw92PhaseError(Sw92PhaseErrorCode::near_multiple,"Sw92Phase: near-multiple PR roots");
    case Sw92RootStatus::iteration_limit:
        throw Sw92PhaseError(Sw92PhaseErrorCode::iteration_limit,"Sw92Phase: PR root iteration limit");
    case Sw92RootStatus::unrepresentable:
        throw Sw92PhaseError(Sw92PhaseErrorCode::unrepresentable_root,"Sw92Phase: unrepresentable PR root/domain");
    }
    throw std::invalid_argument("Sw92Phase: unknown root status");
}
} // namespace detail

/// Explicit-family PT candidate kernel. `family` fixes AQ/NA parameterization;
/// `root_index` then selects one increasing-Z algebraic branch. This is not a
/// stability test, phase-count decision, flash, or physical phase label.
template <std::floating_point T = double>
    requires std::same_as<T, std::remove_cv_t<T>>
class Sw92Phase {
public:
    Sw92Phase(const Sw92Phase&) = default;
    Sw92Phase(Sw92Phase&&) noexcept = default;
    Sw92Phase& operator=(const Sw92Phase&) = delete;
    Sw92Phase& operator=(Sw92Phase&&) = delete;

    [[nodiscard]] static Sw92Phase from_parameters(const Sw92ParameterSet& p,
                                                    ContractLimits limits = {}) {
        return Sw92Phase(p,limits);
    }
    [[nodiscard]] std::size_t size() const noexcept { return mixture_.size(); }
    [[nodiscard]] const Sw92ParameterSet& parameters() const & noexcept { return mixture_.parameters(); }
    const Sw92ParameterSet& parameters() const && = delete;

    [[nodiscard]] Sw92RootSet<T> roots(T pressure_pa,T temperature_k,
        std::span<const T> x,T molality,SwPhaseFamily family,
        Sw92PhaseWorkspace<T>& workspace,Sw92RootOptions options={}) const {
        check_pressure(pressure_pa);
        const auto mixed=mixture_.evaluate(temperature_k,x,molality,family,workspace.mixture_);
        const auto ab=coefficients(pressure_pa,temperature_k,mixed);
        return pr76_roots(ab.a,ab.b,options);
    }

    [[nodiscard]] Sw92PhaseValues<T> evaluate(T pressure_pa,T temperature_k,
        std::span<const T> x,T molality,SwPhaseFamily family,std::size_t root_index,
        Sw92PhaseWorkspace<T>& workspace,Sw92RootOptions options={}) const {
        check_pressure(pressure_pa); workspace.sums_.resize(size());
        const auto mixed=mixture_.template evaluate_impl<true>(temperature_k,x,molality,family,
            workspace.mixture_,std::span<T>{workspace.sums_});
        const auto ab=coefficients(pressure_pa,temperature_k,mixed);
        const auto roots=pr76_roots(ab.a,ab.b,options); detail::sw92_require_resolved(roots.status);
        if(root_index>=roots.count) throw std::out_of_range("Sw92Phase: invalid root index");
        const auto& root=roots.roots[root_index]; const T z=root.z, free=root.free_volume_z;
        const T sqrt2=std::sqrt(T{2});
        const T denominator=free+(T{2}-sqrt2)*ab.b;
        if(!std::isfinite(denominator)||!(denominator>T{0}))
            throw std::range_error("Sw92Phase: PR log denominator");
        const T u=(T{2}*sqrt2)*(ab.b/denominator);
        const T log_factor=detail::sw92_log1p_over_x(u)/denominator;
        if(!std::isfinite(log_factor)) throw std::range_error("Sw92Phase: PR log factor");

        T zm1=z-T{1}; T log_free=std::log(free);
        if(std::abs(zm1)<T{1}/T{4} && ab.b<T{1}/T{8}) {
            const T bz=ab.b/z;
            zm1=ab.b/free-(ab.a/z)/(T{1}+T{2}*bz-bz*bz);
            log_free=std::log1p(zm1-ab.b);
        }
        if(!std::isfinite(zm1)||!std::isfinite(log_free))
            throw std::range_error("Sw92Phase: PR free-volume log");

        const T rt=Sw92Pure<T>::gas_constant()*temperature_k;
        const T p_over_rt=pressure_pa/rt;
        std::vector<T> lnphi(size());
        for(std::size_t i=0;i<size();++i) {
            const T ratio=covolumes_[i]/mixed.b;
            const T coefficient=((T{2}*workspace.sums_[i]-mixed.a*ratio)/rt)*p_over_rt;
            lnphi[i]=ratio*zm1-log_free-coefficient*log_factor;
            if(!std::isfinite(lnphi[i])) throw std::range_error("Sw92Phase: nonrepresentable ln(phi)");
        }
        return {z,std::move(lnphi),roots,root_index,family};
    }

private:
    Sw92Phase(const Sw92ParameterSet& p,ContractLimits limits)
        : mixture_(Sw92Mixture<T>::from_parameters(p,limits)) {
        covolumes_.reserve(p.components().size());
        for(std::size_t i=0;i<p.components().size();++i)
            covolumes_.push_back(Sw92Pure<T>::from_parameters(p,i).covolume());
    }
    void check_pressure(T p) const {
        if(!std::isfinite(p)||!(p>T{0})) throw std::domain_error("Sw92Phase: finite pressure >0 Pa required");
        if(const auto& b=parameters().applicability().state.pressure_pa; b &&
           (static_cast<long double>(p)<b->lower||static_cast<long double>(p)>b->upper))
            throw std::domain_error("Sw92Phase: pressure outside declared dataset interval");
    }
    [[nodiscard]] static Sw92MixtureValues<T> coefficients(T p,T t,const Sw92MixtureValues<T>& m) {
        if(!std::isfinite(t)||!(t>T{0})) throw std::domain_error("Sw92Phase: finite temperature >0 K required");
        const T rt=Sw92Pure<T>::gas_constant()*t, q=p/rt;
        const T A=(m.a/rt)*q, B=m.b*q;
        if(!std::isfinite(A)||!std::isfinite(B)||!(B>T{0})||(A==T{0}&&m.a!=T{0}))
            throw std::range_error("Sw92Phase: nonrepresentable PR A/B");
        return {A,B};
    }

    Sw92Mixture<T> mixture_;
    std::vector<T> covolumes_;
};

} // namespace mpmc::thermodynamics

#endif // MPMC_THERMODYNAMICS_SW92_PHASE_HPP
