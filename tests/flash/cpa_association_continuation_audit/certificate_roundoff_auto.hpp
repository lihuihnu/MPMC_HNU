#ifndef MPMC_TEST_CPA_ASSOCIATION_OBSERVABLE_CERTIFICATE_ROUNDOFF_AUTO_HPP
#define MPMC_TEST_CPA_ASSOCIATION_OBSERVABLE_CERTIFICATE_ROUNDOFF_AUTO_HPP

#include "certificate_roundoff.hpp"

#include <mpmc/thermodynamics/cpa_phase.hpp>

#include <cmath>
#include <span>

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

} // namespace cpa_association_observable_certificate

#endif // MPMC_TEST_CPA_ASSOCIATION_OBSERVABLE_CERTIFICATE_ROUNDOFF_AUTO_HPP
