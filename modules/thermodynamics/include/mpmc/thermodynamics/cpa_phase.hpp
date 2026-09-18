#ifndef MPMC_THERMODYNAMICS_CPA_PHASE_HPP
#define MPMC_THERMODYNAMICS_CPA_PHASE_HPP

#include <mpmc/thermodynamics/cpa_association.hpp>
#include <mpmc/thermodynamics/cpa_helmholtz_derivatives.hpp>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <span>
#include <stdexcept>
#include <string>
#include <vector>

namespace mpmc::thermodynamics {

struct CpaPhaseOptions {
    CpaAssociationOptions association;
};

struct CpaPhaseState {
    double temperature_k{};
    double molar_density_mol_per_m3{};
    double a_mix_pa_m6_per_mol2{};
    double b_mix_m3_per_mol{};
    double pressure_physical_pa{};
    double pressure_association_pa{};
    double pressure_pa{};
    CpaAssociationResult association;
};

namespace cpa_detail {

inline double cpa_ai(double temperature_k, const CpaPureParameters& pure) {
    if (!std::isfinite(temperature_k) || !(temperature_k > 0.0)) {
        throw std::domain_error("CPA phase: positive finite temperature required");
    }
    const double tr = temperature_k / pure.critical_temperature_k;
    if (!std::isfinite(tr) || !(tr > 0.0)) {
        throw std::range_error("CPA phase: nonrepresentable reduced temperature");
    }
    const double alpha_base =
        1.0 + pure.c1_dimensionless * (1.0 - std::sqrt(tr));
    const double value = pure.a0_pa_m6_per_mol2 * alpha_base * alpha_base;
    if (!std::isfinite(value) || !(value > 0.0)) {
        throw std::range_error("CPA phase: nonrepresentable pure attraction parameter");
    }
    return value;
}

inline double cpa_a_mix(double temperature_k,
                        std::span<const double> composition,
                        const CpaParameterSet& parameters) {
    std::vector<double> pure_a(parameters.size(), 0.0);
    for (std::size_t i = 0; i < parameters.size(); ++i) {
        pure_a[i] = cpa_ai(temperature_k, parameters.pure(i));
    }
    double value = 0.0;
    for (std::size_t i = 0; i < parameters.size(); ++i) {
        for (std::size_t j = 0; j < parameters.size(); ++j) {
            const double aij = std::sqrt(pure_a[i] * pure_a[j]) *
                               (1.0 - parameters.kij(i, j));
            value += composition[i] * composition[j] * aij;
        }
    }
    if (!std::isfinite(value) || !(value > 0.0)) {
        throw std::range_error("CPA phase: nonrepresentable mixture attraction parameter");
    }
    return value;
}

} // namespace cpa_detail

// Density-state CPA kernel. This is intentionally below the future p,T,x root
// solver/fugacity layer: callers provide molar density explicitly. The kernel
// supplies a validated pressure decomposition and converged association state.
[[nodiscard]] inline CpaPhaseState evaluate_cpa_phase_at_density(
    double temperature_k,
    double molar_density_mol_per_m3,
    std::span<const double> composition,
    const CpaParameterSet& parameters,
    CpaPhaseOptions options = {}) {
    (void)cpa_detail::validate_cpa_composition(composition, parameters.size());
    if (!std::isfinite(temperature_k) || !(temperature_k > 0.0) ||
        !std::isfinite(molar_density_mol_per_m3) ||
        !(molar_density_mol_per_m3 > 0.0)) {
        throw std::domain_error("CPA phase: positive finite T and molar density required");
    }

    CpaPhaseState result;
    result.temperature_k = temperature_k;
    result.molar_density_mol_per_m3 = molar_density_mol_per_m3;
    result.b_mix_m3_per_mol = cpa_detail::cpa_b_mix(composition, parameters);
    result.a_mix_pa_m6_per_mol2 =
        cpa_detail::cpa_a_mix(temperature_k, composition, parameters);

    const double b_rho = result.b_mix_m3_per_mol * molar_density_mol_per_m3;
    if (!std::isfinite(b_rho) || !(b_rho < 1.0)) {
        throw std::domain_error("CPA phase: density lies at or beyond SRK covolume singularity");
    }

    result.association = solve_cpa_association(
        temperature_k, molar_density_mol_per_m3,
        composition, parameters, options.association);
    if (!result.association.converged()) {
        throw std::runtime_error(
            "CPA phase: association state unavailable: " + result.association.diagnostic);
    }

    const auto pressure = cpa_helmholtz_pressure(
        temperature_k,
        molar_density_mol_per_m3,
        composition,
        parameters,
        result.association);
    result.pressure_physical_pa = pressure.pressure_physical_pa;
    result.pressure_association_pa = pressure.pressure_association_pa;
    result.pressure_pa = pressure.pressure_pa;
    return result;
}

} // namespace mpmc::thermodynamics

#endif // MPMC_THERMODYNAMICS_CPA_PHASE_HPP
