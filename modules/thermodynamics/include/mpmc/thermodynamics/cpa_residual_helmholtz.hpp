#ifndef MPMC_THERMODYNAMICS_CPA_RESIDUAL_HELMHOLTZ_HPP
#define MPMC_THERMODYNAMICS_CPA_RESIDUAL_HELMHOLTZ_HPP

#include <mpmc/thermodynamics/cpa_association.hpp>

#include <cmath>
#include <cstddef>
#include <span>
#include <stdexcept>

namespace mpmc::thermodynamics {
namespace cpa_detail {

// Keep the thermodynamic kernel independent of the AD module. Unqualified
// dependent calls use std math for ordinary floating-point Scalars and allow
// ADL to find the corresponding operation for an external scalar type such as
// mpmc::ad::Dual when that scalar's math overloads are visible to the caller.
template <class Scalar>
[[nodiscard]] Scalar cpa_scalar_sqrt(const Scalar& value) {
    using std::sqrt;
    return sqrt(value);
}

template <class Scalar>
[[nodiscard]] Scalar cpa_scalar_log(const Scalar& value) {
    using std::log;
    return log(value);
}

template <class Scalar>
[[nodiscard]] Scalar cpa_scalar_log1p(const Scalar& value) {
    using std::log1p;
    return log1p(value);
}

template <class Scalar>
[[nodiscard]] Scalar cpa_scalar_expm1(const Scalar& value) {
    using std::expm1;
    return expm1(value);
}

template <class Scalar>
[[nodiscard]] Scalar cpa_ai_generic(
    const Scalar& temperature_k,
    const CpaPureParameters& pure) {
    const Scalar tr = temperature_k / Scalar{pure.critical_temperature_k};
    const Scalar alpha_base =
        Scalar{1.0} + Scalar{pure.c1_dimensionless} *
            (Scalar{1.0} - cpa_scalar_sqrt(tr));
    return Scalar{pure.a0_pa_m6_per_mol2} * alpha_base * alpha_base;
}

template <class Scalar>
[[nodiscard]] Scalar cpa_total_moles(
    std::span<const Scalar> mole_numbers) {
    Scalar total{};
    for (const auto& value : mole_numbers) { total += value; }
    return total;
}

template <class Scalar>
[[nodiscard]] Scalar cpa_extensive_b(
    std::span<const Scalar> mole_numbers,
    const CpaParameterSet& parameters) {
    if (mole_numbers.size() != parameters.size()) {
        throw std::invalid_argument(
            "CPA Helmholtz: mole-number dimension mismatch");
    }
    Scalar value{};
    for (std::size_t i = 0U; i < parameters.size(); ++i) {
        value += mole_numbers[i] * Scalar{parameters.pure(i).b_m3_per_mol};
    }
    return value;
}

template <class Scalar>
[[nodiscard]] Scalar cpa_extensive_a(
    const Scalar& temperature_k,
    std::span<const Scalar> mole_numbers,
    const CpaParameterSet& parameters) {
    if (mole_numbers.size() != parameters.size()) {
        throw std::invalid_argument(
            "CPA Helmholtz: mole-number dimension mismatch");
    }
    Scalar value{};
    for (std::size_t i = 0U; i < parameters.size(); ++i) {
        const Scalar ai = cpa_ai_generic(temperature_k, parameters.pure(i));
        for (std::size_t j = 0U; j < parameters.size(); ++j) {
            const Scalar aj = cpa_ai_generic(temperature_k, parameters.pure(j));
            const Scalar aij = cpa_scalar_sqrt(ai * aj) *
                Scalar{1.0 - parameters.kij(i, j)};
            value += mole_numbers[i] * mole_numbers[j] * aij;
        }
    }
    return value;
}

inline void validate_primal_association_for_helmholtz(
    const CpaAssociationResult& association,
    const CpaParameterSet& parameters) {
    if (!association.converged()) {
        throw std::invalid_argument(
            "CPA Helmholtz: converged primal association state required");
    }
    for (const auto& site : association.sites) {
        if (site.component_index >= parameters.size() ||
            site.multiplicity == 0U ||
            !std::isfinite(site.unbonded_fraction) ||
            !(site.unbonded_fraction > 0.0) ||
            site.unbonded_fraction > 1.0) {
            throw std::invalid_argument(
                "CPA Helmholtz: invalid primal association site state");
        }
    }
}

} // namespace cpa_detail

// Dimensionless EXTENSIVE SRK residual Helmholtz energy F=A^res/(R T).
// Independent variables are (T,V,n_i), not normalized composition. This is the
// canonical cubic term for the current CPA/SRK profile:
//
// F_cubic = -n ln(1-B/V) - A/(B R T) ln(1+B/V),
// B=sum_i n_i b_i, A=sum_ij n_i n_j a_ij(T).
template <class Scalar>
[[nodiscard]] Scalar cpa_cubic_residual_helmholtz_reduced(
    const Scalar& temperature_k,
    const Scalar& volume_m3,
    std::span<const Scalar> mole_numbers,
    const CpaParameterSet& parameters) {
    if (mole_numbers.size() != parameters.size() || mole_numbers.empty()) {
        throw std::invalid_argument(
            "CPA Helmholtz: nonempty mole-number vector must match parameter set");
    }

    const Scalar total_moles = cpa_detail::cpa_total_moles(mole_numbers);
    const Scalar extensive_b =
        cpa_detail::cpa_extensive_b(mole_numbers, parameters);
    const Scalar extensive_a = cpa_detail::cpa_extensive_a(
        temperature_k, mole_numbers, parameters);
    const Scalar reduced_covolume = extensive_b / volume_m3;
    const Scalar rt = Scalar{cpa_gas_constant_j_per_mol_k} * temperature_k;

    return -total_moles *
               cpa_detail::cpa_scalar_log1p(-reduced_covolume) -
           (extensive_a / (extensive_b * rt)) *
               cpa_detail::cpa_scalar_log1p(reduced_covolume);
}

// Stationary Michelsen/CPA association Q function, also dimensionless and
// extensive. The supplied association state is the converged PRIMAL X* at the
// same (T,V,n_i,parameters) state. X* is intentionally held constant while T,
// V and n_i may be active scalar variables. At the stationary point,
// dF_assoc/dz = partial Q/partial z for every first-order z in {T,V,n_i}; no
// differentiation through the fixed-point iteration is required.
template <class Scalar>
[[nodiscard]] Scalar cpa_association_q_reduced(
    const Scalar& temperature_k,
    const Scalar& volume_m3,
    std::span<const Scalar> mole_numbers,
    const CpaParameterSet& parameters,
    const CpaAssociationResult& primal_association) {
    if (mole_numbers.size() != parameters.size() || mole_numbers.empty()) {
        throw std::invalid_argument(
            "CPA Helmholtz: nonempty mole-number vector must match parameter set");
    }
    cpa_detail::validate_primal_association_for_helmholtz(
        primal_association, parameters);
    if (primal_association.sites.empty()) { return Scalar{}; }

    Scalar first_term{};
    for (const auto& site : primal_association.sites) {
        const Scalar x_site{site.unbonded_fraction};
        const Scalar site_count{static_cast<double>(site.multiplicity)};
        first_term += mole_numbers[site.component_index] * site_count *
            (cpa_detail::cpa_scalar_log(x_site) - x_site + Scalar{1.0});
    }

    const Scalar extensive_b =
        cpa_detail::cpa_extensive_b(mole_numbers, parameters);
    const Scalar g = Scalar{1.0} /
        (Scalar{1.0} - Scalar{0.475} * extensive_b / volume_m3);
    const Scalar rt = Scalar{cpa_gas_constant_j_per_mol_k} * temperature_k;

    Scalar pair_sum{};
    for (const auto& first : primal_association.sites) {
        for (const auto& second : primal_association.sites) {
            const auto* pair = parameters.association_pair(
                first.component_index, first.site_id,
                second.component_index, second.site_id);
            if (pair == nullptr || pair->beta_dimensionless == 0.0) { continue; }

            const double b_ij = 0.5 *
                (parameters.pure(first.component_index).b_m3_per_mol +
                 parameters.pure(second.component_index).b_m3_per_mol);
            const Scalar delta = g * cpa_detail::cpa_scalar_expm1(
                Scalar{pair->epsilon_j_per_mol} / rt) *
                Scalar{b_ij * pair->beta_dimensionless};
            const Scalar first_site_count{
                static_cast<double>(first.multiplicity)};
            const Scalar second_site_count{
                static_cast<double>(second.multiplicity)};
            pair_sum +=
                mole_numbers[first.component_index] *
                mole_numbers[second.component_index] *
                first_site_count * second_site_count *
                Scalar{first.unbonded_fraction} *
                Scalar{second.unbonded_fraction} * delta;
        }
    }

    return first_term - pair_sum / (Scalar{2.0} * volume_m3);
}

// Canonical first-order CPA residual potential for the current SRK+sCPA
// profile. The association state must be solved once at the primal state before
// this algebraic kernel is evaluated.
template <class Scalar>
[[nodiscard]] Scalar cpa_residual_helmholtz_reduced(
    const Scalar& temperature_k,
    const Scalar& volume_m3,
    std::span<const Scalar> mole_numbers,
    const CpaParameterSet& parameters,
    const CpaAssociationResult& primal_association) {
    return cpa_cubic_residual_helmholtz_reduced(
               temperature_k, volume_m3, mole_numbers, parameters) +
           cpa_association_q_reduced(
               temperature_k, volume_m3, mole_numbers,
               parameters, primal_association);
}

} // namespace mpmc::thermodynamics

#endif // MPMC_THERMODYNAMICS_CPA_RESIDUAL_HELMHOLTZ_HPP
