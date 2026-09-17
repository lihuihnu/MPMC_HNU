#ifndef MPMC_TEST_CPA_ASSOCIATION_OBSERVABLE_CERTIFICATE_ROUNDOFF_AUTO_HPP
#define MPMC_TEST_CPA_ASSOCIATION_OBSERVABLE_CERTIFICATE_ROUNDOFF_AUTO_HPP

#include "certificate_roundoff.hpp"

#include <mpmc/thermodynamics/cpa_phase.hpp>
#include <mpmc/thermodynamics/cpa_pt_phase.hpp>

#include <algorithm>
#include <cmath>
#include <span>
#include <vector>

namespace cpa_association_observable_certificate {

inline Result certify_with_pressure_roundoff_auto(
    double temperature_k,
    double molar_density_mol_per_m3,
    std::span<const double> composition,
    const th::CpaParameterSet& parameters,
    const th::CpaAssociationResult& association) {
    const double b_mix = th::cpa_detail::cpa_b_mix(composition, parameters);
    const double a_mix = th::cpa_detail::cpa_a_mix(
        temperature_k, composition, parameters);
    const double rt = th::cpa_gas_constant_j_per_mol_k * temperature_k;
    const double b_rho = b_mix * molar_density_mol_per_m3;
    if (!std::isfinite(b_rho) || !(b_rho < 1.0)) {
        return Result{};
    }
    const double pressure_physical_pa =
        rt * molar_density_mol_per_m3 / (1.0 - b_rho) -
        a_mix * molar_density_mol_per_m3 * molar_density_mol_per_m3 /
            (1.0 + b_rho);
    if (!std::isfinite(pressure_physical_pa)) {
        return Result{};
    }
    return certify_with_pressure_roundoff(
        temperature_k,
        molar_density_mol_per_m3,
        composition,
        parameters,
        association,
        pressure_physical_pa);
}

inline Result certify_with_full_roundoff_auto(
    double target_pressure_pa,
    double temperature_k,
    double molar_density_mol_per_m3,
    std::span<const double> composition,
    const th::CpaParameterSet& parameters,
    const th::CpaAssociationResult& association) {
    Result result = certify_with_pressure_roundoff_auto(
        temperature_k,
        molar_density_mol_per_m3,
        composition,
        parameters,
        association);
    if (!result.enclosure_certified ||
        association.status == th::CpaAssociationStatus::no_associating_sites) {
        return result;
    }

    const double rho = molar_density_mol_per_m3;
    const double rt = th::cpa_gas_constant_j_per_mol_k * temperature_k;
    const double b = th::cpa_detail::cpa_b_mix(composition, parameters);
    const double a = th::cpa_detail::cpa_a_mix(
        temperature_k, composition, parameters);
    const double b_rho = b * rho;
    const double z = target_pressure_pa / (rho * rt);
    const double pressure_physical_pa =
        rt * rho / (1.0 - b_rho) - a * rho * rho / (1.0 + b_rho);
    const double z_physical = pressure_physical_pa / (rho * rt);
    const double a_over_brt = a / (b * rt);
    const double log_free_volume = std::log1p(-b_rho);
    const double log_attraction_volume = std::log1p(b_rho);
    const double log_z = std::log(z);
    if (!std::isfinite(z) || !(z > 0.0) ||
        !std::isfinite(z_physical) || !std::isfinite(a_over_brt) ||
        !std::isfinite(log_free_volume) ||
        !std::isfinite(log_attraction_volume) ||
        !std::isfinite(log_z)) {
        return Result{};
    }

    const auto pure_a = th::cpa_detail::cpa_pure_a(temperature_k, parameters);
    const auto sums = th::cpa_detail::cpa_a_sums(
        composition, pure_a, parameters);

    std::vector<double> association_log_x(parameters.size(), 0.0);
    double association_sum = 0.0;
    for (const auto& site : association.sites) {
        const double multiplicity = static_cast<double>(site.multiplicity);
        association_log_x[site.component_index] +=
            multiplicity * std::log(site.unbonded_fraction);
        association_sum += composition[site.component_index] * multiplicity *
            (1.0 - site.unbonded_fraction);
    }
    const double g = association.radial_distribution;
    if (!std::isfinite(g) || !(g > 0.0) ||
        !std::isfinite(association_sum) || association_sum < 0.0) {
        return Result{};
    }

    // The X-dependent analytic and mu_assoc arithmetic terms were already
    // included by certify_with_pressure_roundoff_auto().  What remains is the
    // final production rounding of
    //     mu_cubic + mu_association - log(z).
    // All non-association operands are identical for warm/cold candidates at
    // fixed (p,T,rho,x).  Compare production double evaluation with the exact
    // long-double sum of those same stored double operands; this is fully
    // candidate-local and needs no cold reference.
    Real max_augmented_ln_phi = 0.0L;
    for (std::size_t component = 0U;
         component < parameters.size(); ++component) {
        const double b_i = parameters.pure(component).b_m3_per_mol;
        const double b_ratio = b_i / b;
        const double attraction_ratio =
            2.0 * sums[component] / a - b_ratio;
        const double mu_cubic =
            b_ratio * (z_physical - 1.0) - log_free_volume -
            a_over_brt * attraction_ratio * log_attraction_volume;
        const double mu_association = association_log_x[component] -
            (1.9 / 8.0) * rho * b_i * g * association_sum;
        const double value = mu_cubic + mu_association - log_z;
        if (!std::isfinite(mu_cubic) || !std::isfinite(mu_association) ||
            !std::isfinite(value)) {
            return Result{};
        }

        const Real exact_sum_of_double_operands =
            static_cast<Real>(mu_cubic) +
            static_cast<Real>(mu_association) -
            static_cast<Real>(log_z);
        const Real final_addition_gap = std::abs(
            static_cast<Real>(value) - exact_sum_of_double_operands);
        const Real diagnostic_pad =
            64.0L * std::numeric_limits<Real>::epsilon() *
            std::max(1.0L, std::abs(exact_sum_of_double_operands));
        const Real augmented =
            static_cast<Real>(result.ln_phi_error_bounds[component]) +
            final_addition_gap + diagnostic_pad;
        if (!finite(augmented) || augmented < 0.0L) {
            return Result{};
        }
        result.ln_phi_error_bounds[component] = static_cast<double>(augmented);
        max_augmented_ln_phi = std::max(max_augmented_ln_phi, augmented);
    }
    result.max_ln_phi_error_bound =
        static_cast<double>(max_augmented_ln_phi);
    result.ln_phi_scale_certified =
        max_augmented_ln_phi <= static_cast<Real>(contract::ln_phi_error_guard);
    return result;
}

} // namespace cpa_association_observable_certificate

#endif // MPMC_TEST_CPA_ASSOCIATION_OBSERVABLE_CERTIFICATE_ROUNDOFF_AUTO_HPP
