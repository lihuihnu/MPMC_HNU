#ifndef MPMC_THERMODYNAMICS_CPA_ASSOCIATION_HPP
#define MPMC_THERMODYNAMICS_CPA_ASSOCIATION_HPP

#include <mpmc/thermodynamics/cpa_parameters.hpp>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <span>
#include <stdexcept>
#include <string>
#include <vector>

namespace mpmc::thermodynamics {

inline constexpr double cpa_gas_constant_j_per_mol_k = 8.31446261815324;

struct CpaAssociationOptions {
    double site_fraction_tolerance{1e-12};
    double damping{0.5};
    int max_iterations{512};
    std::size_t max_site_classes{256U};
};

enum class CpaAssociationStatus {
    success,
    no_associating_sites,
    invalid_state,
    radial_distribution_singularity,
    nonrepresentable_strength,
    iteration_limit
};

struct CpaAssociationSiteState {
    std::size_t component_index{};
    std::string site_id;
    std::size_t multiplicity{};
    double unbonded_fraction{1.0};
};

struct CpaAssociationResult {
    CpaAssociationStatus status{CpaAssociationStatus::invalid_state};
    double temperature_k{};
    double molar_density_mol_per_m3{};
    double b_mix_m3_per_mol{};
    double eta{};
    double radial_distribution{};
    double rho_dln_g_drho{};
    std::vector<CpaAssociationSiteState> sites;
    double max_fixed_point_residual{};
    int iterations{};
    std::string diagnostic;

    [[nodiscard]] bool converged() const noexcept {
        return status == CpaAssociationStatus::success ||
               status == CpaAssociationStatus::no_associating_sites;
    }
};

namespace cpa_detail {

inline double validate_cpa_composition(std::span<const double> composition,
                                       std::size_t expected_size) {
    if (composition.size() != expected_size || composition.empty()) {
        throw std::invalid_argument("CPA: composition dimension mismatch");
    }
    double sum = 0.0;
    double correction = 0.0;
    for (const double value : composition) {
        if (!std::isfinite(value) || value < 0.0 || value > 1.0) {
            throw std::domain_error("CPA: composition entries must be finite in [0,1]");
        }
        const double y = value - correction;
        const double t = sum + y;
        correction = (t - sum) - y;
        sum = t;
    }
    const double tolerance = 64.0 * std::numeric_limits<double>::epsilon();
    if (!std::isfinite(sum) || std::abs(sum - 1.0) > tolerance) {
        throw std::domain_error("CPA: composition must sum to one within roundoff tolerance");
    }
    return sum;
}

inline void validate_association_options(const CpaAssociationOptions& options) {
    if (!std::isfinite(options.site_fraction_tolerance) ||
        !(options.site_fraction_tolerance > 0.0) ||
        !std::isfinite(options.damping) || !(options.damping > 0.0) ||
        options.damping > 1.0 || options.max_iterations < 0 ||
        options.max_site_classes == 0U) {
        throw std::invalid_argument("CPA association: invalid solver options");
    }
}

inline double cpa_b_mix(std::span<const double> composition,
                        const CpaParameterSet& parameters) {
    double value = 0.0;
    for (std::size_t i = 0; i < composition.size(); ++i) {
        value += composition[i] * parameters.pure(i).b_m3_per_mol;
    }
    if (!std::isfinite(value) || !(value > 0.0)) {
        throw std::range_error("CPA: nonrepresentable mixture covolume");
    }
    return value;
}

inline double cpa_association_strength(
    double temperature_k, double radial_distribution,
    double b_ij_m3_per_mol,
    const CpaAssociationPairParameters& pair) {
    const double exponent =
        pair.epsilon_j_per_mol / (cpa_gas_constant_j_per_mol_k * temperature_k);
    if (!std::isfinite(exponent) || exponent > 700.0) {
        throw std::range_error("CPA association: nonrepresentable association exponent");
    }
    const double strength = radial_distribution * std::expm1(exponent) *
                            b_ij_m3_per_mol * pair.beta_dimensionless;
    if (!std::isfinite(strength) || strength < 0.0) {
        throw std::range_error("CPA association: nonrepresentable association strength");
    }
    return strength;
}

} // namespace cpa_detail

[[nodiscard]] inline CpaAssociationResult solve_cpa_association(
    double temperature_k,
    double molar_density_mol_per_m3,
    std::span<const double> composition,
    const CpaParameterSet& parameters,
    CpaAssociationOptions options = {}) {
    cpa_detail::validate_association_options(options);
    (void)cpa_detail::validate_cpa_composition(composition, parameters.size());
    if (!std::isfinite(temperature_k) || !(temperature_k > 0.0) ||
        !std::isfinite(molar_density_mol_per_m3) ||
        !(molar_density_mol_per_m3 > 0.0)) {
        throw std::domain_error("CPA association: positive finite T and molar density required");
    }

    CpaAssociationResult result;
    result.temperature_k = temperature_k;
    result.molar_density_mol_per_m3 = molar_density_mol_per_m3;
    result.b_mix_m3_per_mol = cpa_detail::cpa_b_mix(composition, parameters);
    result.eta = 0.25 * result.b_mix_m3_per_mol * molar_density_mol_per_m3;
    const double denominator = 1.0 - 1.9 * result.eta;
    if (!std::isfinite(denominator) || !(denominator > 0.0)) {
        result.status = CpaAssociationStatus::radial_distribution_singularity;
        result.diagnostic = "CPA association: simplified radial-distribution denominator is nonpositive";
        return result;
    }
    result.radial_distribution = 1.0 / denominator;
    result.rho_dln_g_drho = (1.9 * result.eta) / denominator;

    for (std::size_t component = 0; component < parameters.size(); ++component) {
        for (const auto& site : parameters.pure(component).sites) {
            result.sites.push_back({component, site.id, site.multiplicity, 1.0});
        }
    }
    if (result.sites.empty()) {
        result.status = CpaAssociationStatus::no_associating_sites;
        result.diagnostic = "CPA association: configured component set has no association sites";
        return result;
    }
    if (result.sites.size() > options.max_site_classes) {
        throw std::length_error("CPA association: site-class quota exceeded");
    }

    std::vector<double> current(result.sites.size(), 1.0);
    std::vector<double> target(result.sites.size(), 1.0);
    for (int iteration = 0; iteration <= options.max_iterations; ++iteration) {
        double residual = 0.0;
        try {
            for (std::size_t a = 0; a < result.sites.size(); ++a) {
                const auto& first = result.sites[a];
                double association_sum = 0.0;
                for (std::size_t b = 0; b < result.sites.size(); ++b) {
                    const auto& second = result.sites[b];
                    const auto* pair = parameters.association_pair(
                        first.component_index, first.site_id,
                        second.component_index, second.site_id);
                    if (pair == nullptr || pair->beta_dimensionless == 0.0 ||
                        composition[second.component_index] == 0.0) {
                        continue;
                    }
                    const double b_ij = 0.5 *
                        (parameters.pure(first.component_index).b_m3_per_mol +
                         parameters.pure(second.component_index).b_m3_per_mol);
                    const double delta = cpa_detail::cpa_association_strength(
                        temperature_k, result.radial_distribution, b_ij, *pair);
                    association_sum += composition[second.component_index] *
                        static_cast<double>(second.multiplicity) * current[b] * delta;
                }
                const double value = 1.0 /
                    (1.0 + molar_density_mol_per_m3 * association_sum);
                if (!std::isfinite(value) || !(value > 0.0) || value > 1.0) {
                    result.status = CpaAssociationStatus::nonrepresentable_strength;
                    result.diagnostic = "CPA association: site fraction left (0,1]";
                    return result;
                }
                target[a] = value;
                residual = std::max(residual, std::abs(target[a] - current[a]));
            }
        } catch (const std::range_error& error) {
            result.status = CpaAssociationStatus::nonrepresentable_strength;
            result.diagnostic = error.what();
            return result;
        }

        result.iterations = iteration;
        result.max_fixed_point_residual = residual;
        if (residual <= options.site_fraction_tolerance) {
            for (std::size_t i = 0; i < result.sites.size(); ++i) {
                result.sites[i].unbonded_fraction = target[i];
            }
            result.status = CpaAssociationStatus::success;
            result.diagnostic = "CPA association site fractions converged";
            return result;
        }
        if (iteration == options.max_iterations) { break; }
        for (std::size_t i = 0; i < current.size(); ++i) {
            current[i] += options.damping * (target[i] - current[i]);
        }
    }

    result.status = CpaAssociationStatus::iteration_limit;
    result.diagnostic = "CPA association: fixed-point iteration limit";
    return result;
}

} // namespace mpmc::thermodynamics

#endif // MPMC_THERMODYNAMICS_CPA_ASSOCIATION_HPP
