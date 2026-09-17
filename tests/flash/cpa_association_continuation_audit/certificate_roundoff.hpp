#ifndef MPMC_TEST_CPA_ASSOCIATION_OBSERVABLE_CERTIFICATE_ROUNDOFF_HPP
#define MPMC_TEST_CPA_ASSOCIATION_OBSERVABLE_CERTIFICATE_ROUNDOFF_HPP

#include "certificate.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <span>

namespace cpa_association_observable_certificate {

inline Result certify_with_pressure_roundoff(
    double temperature_k,
    double molar_density_mol_per_m3,
    std::span<const double> composition,
    const th::CpaParameterSet& parameters,
    const th::CpaAssociationResult& association,
    double pressure_physical_pa) {
    Result result = certify(
        temperature_k,
        molar_density_mol_per_m3,
        composition,
        parameters,
        association);
    if (!result.enclosure_certified ||
        association.status == th::CpaAssociationStatus::no_associating_sites) {
        return result;
    }

    // Reproduce the current production association arithmetic in double while
    // evaluating the same expressions from the same stored input doubles in
    // long double. These candidate-local gaps require no cold reference and
    // account for ULP-scale loss caused by dense-state association terms.
    double association_sum_double = 0.0;
    Real association_sum_real = 0.0L;
    for (const auto& site : association.sites) {
        const double weight_double =
            composition[site.component_index] *
            static_cast<double>(site.multiplicity);
        const double one_minus_x_double = 1.0 - site.unbonded_fraction;
        association_sum_double += weight_double * one_minus_x_double;

        const Real weight_real =
            static_cast<Real>(composition[site.component_index]) *
            static_cast<Real>(site.multiplicity);
        const Real one_minus_x_real =
            1.0L - static_cast<Real>(site.unbonded_fraction);
        association_sum_real += weight_real * one_minus_x_real;
    }

    const Real unit_roundoff =
        0.5L * std::numeric_limits<Real>::epsilon();
    constexpr Real operation_budget = 64.0L;
    const Real gamma =
        (operation_budget * unit_roundoff) /
        (1.0L - operation_budget * unit_roundoff);

    const double rt_double = th::cpa_gas_constant_j_per_mol_k * temperature_k;
    const double association_pressure_double =
        -0.5 * rt_double * molar_density_mol_per_m3 *
        (1.0 + association.rho_dln_g_drho) * association_sum_double;
    const double total_pressure_double =
        pressure_physical_pa + association_pressure_double;

    const Real rt_real =
        static_cast<Real>(th::cpa_gas_constant_j_per_mol_k) *
        static_cast<Real>(temperature_k);
    const Real association_pressure_real =
        -0.5L * rt_real * static_cast<Real>(molar_density_mol_per_m3) *
        (1.0L + static_cast<Real>(association.rho_dln_g_drho)) *
        association_sum_real;
    const Real total_pressure_real =
        static_cast<Real>(pressure_physical_pa) + association_pressure_real;
    if (!finite(association_pressure_real) || !finite(total_pressure_real) ||
        !std::isfinite(total_pressure_double)) {
        return Result{};
    }

    const Real observed_pressure_gap = std::abs(
        static_cast<Real>(total_pressure_double) - total_pressure_real);
    const Real pressure_diagnostic_pad = gamma *
        std::max(1.0L, std::abs(association_pressure_real));
    const Real pressure_floating_allowance =
        observed_pressure_gap + pressure_diagnostic_pad;

    const Real augmented_pressure_bound =
        static_cast<Real>(result.pressure_error_bound_pa) +
        pressure_floating_allowance;
    if (!finite(augmented_pressure_bound) || augmented_pressure_bound < 0.0L) {
        return Result{};
    }
    result.pressure_error_bound_pa =
        static_cast<double>(augmented_pressure_bound);
    result.pressure_scale_certified =
        augmented_pressure_bound <=
        static_cast<Real>(contract::pressure_error_guard_pa);

    // The association chemical-potential term is the only X-dependent part of
    // ln(phi) at fixed T, rho and composition. Add a candidate-local arithmetic
    // allowance to the analytic site-box propagation bound for each component.
    if (result.ln_phi_error_bounds.size() != parameters.size()) {
        return Result{};
    }
    Real max_augmented_ln_phi = 0.0L;
    for (std::size_t component = 0U;
         component < parameters.size(); ++component) {
        double mu_assoc_double = 0.0;
        Real mu_assoc_real = 0.0L;
        for (const auto& site : association.sites) {
            if (site.component_index != component) { continue; }
            const double multiplicity_double =
                static_cast<double>(site.multiplicity);
            mu_assoc_double +=
                multiplicity_double * std::log(site.unbonded_fraction);
            mu_assoc_real +=
                static_cast<Real>(site.multiplicity) *
                std::log(static_cast<Real>(site.unbonded_fraction));
        }

        const double sensitivity_double =
            (1.9 / 8.0) * molar_density_mol_per_m3 *
            parameters.pure(component).b_m3_per_mol *
            association.radial_distribution;
        mu_assoc_double -= sensitivity_double * association_sum_double;

        const Real sensitivity_real =
            (1.9L / 8.0L) * static_cast<Real>(molar_density_mol_per_m3) *
            static_cast<Real>(parameters.pure(component).b_m3_per_mol) *
            static_cast<Real>(association.radial_distribution);
        mu_assoc_real -= sensitivity_real * association_sum_real;
        if (!std::isfinite(mu_assoc_double) || !finite(mu_assoc_real)) {
            return Result{};
        }

        const Real observed_mu_gap = std::abs(
            static_cast<Real>(mu_assoc_double) - mu_assoc_real);
        const Real mu_diagnostic_pad = gamma *
            std::max(1.0L, std::abs(mu_assoc_real));
        const Real augmented =
            static_cast<Real>(result.ln_phi_error_bounds[component]) +
            observed_mu_gap + mu_diagnostic_pad;
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

#endif // MPMC_TEST_CPA_ASSOCIATION_OBSERVABLE_CERTIFICATE_ROUNDOFF_HPP
