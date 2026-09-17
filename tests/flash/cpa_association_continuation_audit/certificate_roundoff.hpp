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

    // Reproduce the current production association-pressure arithmetic in
    // double, while evaluating the same expression from the same stored input
    // doubles in long double.  The observed gap is candidate-local and requires
    // no cold reference.  It captures the ULP-scale loss caused by forming a
    // very large association pressure before cancellation with the physical
    // contribution at dense states.
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

    const Real observed_double_gap = std::abs(
        static_cast<Real>(total_pressure_double) - total_pressure_real);

    // Standard floating-point model pad for the long-double diagnostic path.
    // The pressure expression contains far fewer than 64 elementary arithmetic
    // operations for the current bounded site count.  gamma_64 is applied to
    // the large association-pressure magnitude, so cancellation in the final
    // total cannot make the absolute pad artificially small.
    const Real unit_roundoff =
        0.5L * std::numeric_limits<Real>::epsilon();
    constexpr Real operation_budget = 64.0L;
    const Real gamma =
        (operation_budget * unit_roundoff) /
        (1.0L - operation_budget * unit_roundoff);
    const Real diagnostic_pad = gamma *
        std::max(1.0L, std::abs(association_pressure_real));
    const Real floating_allowance = observed_double_gap + diagnostic_pad;

    const Real augmented_pressure_bound =
        static_cast<Real>(result.pressure_error_bound_pa) + floating_allowance;
    if (!finite(augmented_pressure_bound) || augmented_pressure_bound < 0.0L) {
        return Result{};
    }
    result.pressure_error_bound_pa =
        static_cast<double>(augmented_pressure_bound);
    result.pressure_scale_certified =
        augmented_pressure_bound <=
        static_cast<Real>(contract::pressure_error_guard_pa);
    return result;
}

} // namespace cpa_association_observable_certificate

#endif // MPMC_TEST_CPA_ASSOCIATION_OBSERVABLE_CERTIFICATE_ROUNDOFF_HPP
