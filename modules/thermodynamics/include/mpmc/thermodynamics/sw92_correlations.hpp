#ifndef MPMC_THERMODYNAMICS_SW92_CORRELATIONS_HPP
#define MPMC_THERMODYNAMICS_SW92_CORRELATIONS_HPP

#include <mpmc/thermodynamics/sw92_types.hpp>

#include <cmath>
#include <concepts>
#include <stdexcept>

namespace mpmc::thermodynamics {

namespace detail {
template <std::floating_point T>
[[nodiscard]] T sw92_pow_nonnegative(T x, T p) {
    if (!std::isfinite(x) || x < T{0} || !std::isfinite(p))
        throw std::domain_error("SW92 correlation: invalid real power input");
    if (x == T{0}) return T{0};
    const T y = std::pow(x, p);
    if (!std::isfinite(y)) throw std::range_error("SW92 correlation: nonrepresentable power");
    return y;
}
} // namespace detail

/// SW92 Eq.(9). c_sw is mol NaCl/kg H2O, not a composition fraction.
template <std::floating_point T>
[[nodiscard]] T sw92_water_alpha(T temperature_k, T water_tc_k,
                                 T nacl_molality_mol_per_kg_water) {
    if (!std::isfinite(temperature_k) || !(temperature_k > T{0}) ||
        !std::isfinite(water_tc_k) || !(water_tc_k > T{0}) ||
        !std::isfinite(nacl_molality_mol_per_kg_water) ||
        nacl_molality_mol_per_kg_water < T{0})
        throw std::domain_error("sw92_water_alpha: finite T,Tc>0 and molality>=0 required");
    const T tr = temperature_k / water_tc_k;
    const T c11 = detail::sw92_pow_nonnegative(nacl_molality_mol_per_kg_water,
                                               static_cast<T>(1.1L));
    const T tr3 = tr * tr * tr;
    if (!std::isfinite(tr3) || tr3 == T{0})
        throw std::range_error("sw92_water_alpha: nonrepresentable reduced state");
    const T q = T{1} + static_cast<T>(0.4530L) *
        (T{1} - tr * (T{1} - static_cast<T>(0.0103L) * c11)) +
        static_cast<T>(0.0034L) * (T{1}/tr3 - T{1});
    const T alpha = q * q;
    if (!std::isfinite(alpha)) throw std::range_error("sw92_water_alpha: nonrepresentable alpha");
    return alpha;
}

/// Corrected SW92 AQ water-pair rules: corrected Eq.(12)/Table 2 and Eq.(13),
/// plus Eqs.(14),(15). `tr` is T/Tc of the non-water component.
template <std::floating_point T>
[[nodiscard]] T sw92_aqueous_water_kij(Sw92Species species, T tr, T omega,
                                       T nacl_molality_mol_per_kg_water) {
    if (!std::isfinite(tr) || !(tr > T{0}) ||
        !std::isfinite(nacl_molality_mol_per_kg_water) ||
        nacl_molality_mol_per_kg_water < T{0})
        throw std::domain_error("SW92 AQ BIP: invalid reduced state/molality");
    const T c = nacl_molality_mol_per_kg_water;
    T k{};
    switch (species) {
    case Sw92Species::hydrocarbon: {
        if (!std::isfinite(omega) || !(omega > T{0}))
            throw std::domain_error("SW92 Eq.(12): hydrocarbon omega must be >0");
        const T w = std::pow(omega, static_cast<T>(-0.1L));
        if (!std::isfinite(w)) throw std::range_error("SW92 Eq.(12): omega power");
        const T a0 = static_cast<T>(1.1120L) - static_cast<T>(1.7369L)*w;
        const T a1 = static_cast<T>(1.1001L) + static_cast<T>(0.8360L)*omega;
        const T a2 = -static_cast<T>(0.15742L) - static_cast<T>(1.0988L)*omega;
        k = a0*(T{1}+static_cast<T>(0.017407L)*c) +
            a1*tr*(T{1}+static_cast<T>(0.033516L)*c) +
            a2*tr*tr*(T{1}+static_cast<T>(0.011478L)*c);
        break;
    }
    case Sw92Species::nitrogen: {
        const T cp = detail::sw92_pow_nonnegative(c, static_cast<T>(0.75L));
        k = -static_cast<T>(1.70235L)*(T{1}+static_cast<T>(0.025587L)*cp) +
             static_cast<T>(0.44338L)*(T{1}+static_cast<T>(0.08126L)*cp)*tr;
        break;
    }
    case Sw92Species::carbon_dioxide: {
        const T c1 = detail::sw92_pow_nonnegative(c, static_cast<T>(0.7505L));
        const T c2 = detail::sw92_pow_nonnegative(c, static_cast<T>(0.979L));
        const T e = std::exp(-static_cast<T>(6.7222L)*tr - c);
        k = -static_cast<T>(0.31092L)*(T{1}+static_cast<T>(0.15587L)*c1) +
             static_cast<T>(0.23580L)*(T{1}+static_cast<T>(0.17837L)*c2)*tr -
             static_cast<T>(21.2566L)*e;
        break;
    }
    case Sw92Species::hydrogen_sulfide:
        k = -static_cast<T>(0.20441L) + static_cast<T>(0.23426L)*tr;
        break;
    case Sw92Species::water:
    case Sw92Species::unspecified:
        throw std::invalid_argument("SW92 AQ BIP: invalid water-pair species");
    }
    if (!std::isfinite(k)) throw std::range_error("SW92 AQ BIP: nonrepresentable value");
    return k;
}

/// SW92 Eq.(17), used only for H2S/water in the non-aqueous family.
template <std::floating_point T>
[[nodiscard]] T sw92_h2s_nonaqueous_water_kij(T tr) {
    if (!std::isfinite(tr) || !(tr > T{0}))
        throw std::domain_error("SW92 Eq.(17): Tr must be finite and >0");
    return static_cast<T>(0.19031L) - static_cast<T>(0.05965L)*tr;
}

} // namespace mpmc::thermodynamics

#endif // MPMC_THERMODYNAMICS_SW92_CORRELATIONS_HPP
