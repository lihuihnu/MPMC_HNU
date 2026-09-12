#ifndef MPMC_THERMODYNAMICS_CPA_PT_PHASE_HPP
#define MPMC_THERMODYNAMICS_CPA_PT_PHASE_HPP

#include <mpmc/thermodynamics/cpa_phase.hpp>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <numbers>
#include <span>
#include <stdexcept>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

namespace mpmc::thermodynamics {

inline constexpr std::string_view cpa_pt_convention =
    "CPA/SRK-association/density-root-finite-scan/lnphi-helmholtz/v1";

struct CpaPtOptions {
    CpaPhaseOptions phase;
    std::size_t scan_intervals{512U};
    int max_bisection_iterations{128};
    std::size_t max_evaluations{4096U};
    std::size_t max_roots{16U};
    double pressure_absolute_tolerance_pa{1e-4};
    double pressure_relative_tolerance{1e-10};
};

enum class CpaPtRootStatus {
    success,
    near_multiple,
    no_root,
    property_failure,
    iteration_limit,
    evaluation_limit,
    root_limit
};

struct CpaPtRoot {
    double molar_density_mol_per_m3{};
    double compressibility_factor{};
    double pressure_residual_pa{};
    // +1: resolved negative-to-positive P(rho)-P crossing (mechanically stable
    // simple root); -1: positive-to-negative crossing; 0: tangent/unresolved.
    int pressure_slope_sign{};
    std::vector<double> ln_phi;
    double reduced_gibbs_offset{}; // sum_i x_i ln(phi_i), comparable at fixed p,T,x.
};

struct CpaPtRootSet {
    static constexpr std::string_view convention = cpa_pt_convention;

    CpaPtRootStatus status{CpaPtRootStatus::no_root};
    CpaPtOptions options;
    double pressure_pa{};
    double temperature_k{};
    std::vector<double> composition;
    std::vector<CpaPtRoot> roots;
    std::size_t evaluations{};
    std::string diagnostic;
};

namespace cpa_detail {

inline void validate_cpa_pt_options(const CpaPtOptions& options) {
    if (options.scan_intervals < 16U ||
        options.scan_intervals > options.max_evaluations ||
        options.max_bisection_iterations <= 0 ||
        options.max_evaluations == 0U || options.max_roots == 0U ||
        !std::isfinite(options.pressure_absolute_tolerance_pa) ||
        !(options.pressure_absolute_tolerance_pa > 0.0) ||
        !std::isfinite(options.pressure_relative_tolerance) ||
        !(options.pressure_relative_tolerance > 0.0)) {
        throw std::invalid_argument("CPA PT: invalid root-search options");
    }
    validate_association_options(options.phase.association);
}

inline int cpa_pressure_sign(double residual, double tolerance) noexcept {
    if (residual > tolerance) { return 1; }
    if (residual < -tolerance) { return -1; }
    return 0;
}

struct CpaPtSample {
    double reduced_density{}; // u=b_mix*rho
    double residual_pa{};
};

struct CpaPtRootCandidate {
    double reduced_density{};
    double residual_pa{};
    int slope_sign{};
};

inline std::vector<double> cpa_pure_a(
    double temperature_k, const CpaParameterSet& parameters) {
    std::vector<double> values(parameters.size(), 0.0);
    for (std::size_t i = 0; i < parameters.size(); ++i) {
        values[i] = cpa_ai(temperature_k, parameters.pure(i));
    }
    return values;
}

inline std::vector<double> cpa_a_sums(
    std::span<const double> composition,
    std::span<const double> pure_a,
    const CpaParameterSet& parameters) {
    std::vector<double> sums(parameters.size(), 0.0);
    for (std::size_t i = 0; i < parameters.size(); ++i) {
        for (std::size_t j = 0; j < parameters.size(); ++j) {
            const double aij = std::sqrt(pure_a[i] * pure_a[j]) *
                               (1.0 - parameters.kij(i, j));
            sums[i] += composition[j] * aij;
        }
    }
    return sums;
}

inline void cpa_fill_ln_phi(
    double target_pressure_pa,
    double temperature_k,
    std::span<const double> composition,
    const CpaParameterSet& parameters,
    const CpaPhaseState& state,
    CpaPtRoot& root) {
    const double rho = state.molar_density_mol_per_m3;
    const double rt = cpa_gas_constant_j_per_mol_k * temperature_k;
    const double z = target_pressure_pa / (rho * rt);
    const double b = state.b_mix_m3_per_mol;
    const double a = state.a_mix_pa_m6_per_mol2;
    const double b_rho = b * rho;
    if (!std::isfinite(z) || !(z > 0.0) ||
        !std::isfinite(b_rho) || !(b_rho > 0.0) || !(b_rho < 1.0)) {
        throw std::range_error("CPA PT: nonrepresentable Z or reduced density");
    }

    const auto pure_a = cpa_pure_a(temperature_k, parameters);
    const auto sums = cpa_a_sums(composition, pure_a, parameters);
    const double z_physical = state.pressure_physical_pa / (rho * rt);
    const double a_over_brt = a / (b * rt);
    const double log_free_volume = std::log1p(-b_rho);
    const double log_attraction_volume = std::log1p(b_rho);
    if (!std::isfinite(z_physical) || !std::isfinite(a_over_brt) ||
        !std::isfinite(log_free_volume) ||
        !std::isfinite(log_attraction_volume)) {
        throw std::range_error("CPA PT: nonrepresentable SRK chemical potential state");
    }

    std::vector<double> association_log_x(parameters.size(), 0.0);
    double association_sum = 0.0;
    for (const auto& site : state.association.sites) {
        if (!(site.unbonded_fraction > 0.0) ||
            !std::isfinite(site.unbonded_fraction)) {
            throw std::range_error("CPA PT: invalid converged association site fraction");
        }
        const double multiplicity = static_cast<double>(site.multiplicity);
        association_log_x[site.component_index] +=
            multiplicity * std::log(site.unbonded_fraction);
        association_sum += composition[site.component_index] * multiplicity *
            (1.0 - site.unbonded_fraction);
    }

    const double g = state.association.radial_distribution;
    if (!std::isfinite(g) || !(g > 0.0) ||
        !std::isfinite(association_sum) || association_sum < 0.0) {
        throw std::range_error("CPA PT: nonrepresentable association chemical potential state");
    }

    root.compressibility_factor = z;
    root.ln_phi.assign(parameters.size(), 0.0);
    root.reduced_gibbs_offset = 0.0;
    for (std::size_t i = 0; i < parameters.size(); ++i) {
        const double b_i = parameters.pure(i).b_m3_per_mol;
        const double b_ratio = b_i / b;
        const double attraction_ratio = 2.0 * sums[i] / a - b_ratio;
        const double mu_cubic =
            b_ratio * (z_physical - 1.0) - log_free_volume -
            a_over_brt * attraction_ratio * log_attraction_volume;
        // For g=1/(1-1.9*eta), (1/g) dg/deta = 1.9*g.
        const double mu_association = association_log_x[i] -
            (1.9 / 8.0) * rho * b_i * g * association_sum;
        const double value = mu_cubic + mu_association - std::log(z);
        if (!std::isfinite(mu_cubic) || !std::isfinite(mu_association) ||
            !std::isfinite(value)) {
            throw std::range_error("CPA PT: nonrepresentable ln(phi)");
        }
        root.ln_phi[i] = value;
        root.reduced_gibbs_offset += composition[i] * value;
    }
    if (!std::isfinite(root.reduced_gibbs_offset)) {
        throw std::range_error("CPA PT: nonrepresentable reduced Gibbs root offset");
    }
}

} // namespace cpa_detail

class CpaPtPhase {
public:
    CpaPtPhase(const CpaPtPhase&) = default;
    CpaPtPhase(CpaPtPhase&&) noexcept = default;
    CpaPtPhase& operator=(const CpaPtPhase&) = delete;
    CpaPtPhase& operator=(CpaPtPhase&&) = delete;

    [[nodiscard]] static CpaPtPhase from_parameters(
        const CpaParameterSet& parameters) {
        return CpaPtPhase(parameters);
    }

    [[nodiscard]] std::size_t size() const noexcept { return parameters_.size(); }
    [[nodiscard]] const CpaParameterSet& parameters() const & noexcept {
        return parameters_;
    }
    const CpaParameterSet& parameters() const && = delete;

    [[nodiscard]] CpaPtRootSet roots(
        double pressure_pa,
        double temperature_k,
        std::span<const double> composition,
        CpaPtOptions options = {}) const {
        cpa_detail::validate_cpa_pt_options(options);
        (void)cpa_detail::validate_cpa_composition(composition, size());
        if (!std::isfinite(pressure_pa) || !(pressure_pa > 0.0) ||
            !std::isfinite(temperature_k) || !(temperature_k > 0.0)) {
            throw std::domain_error("CPA PT: positive finite p and T required");
        }
        if (parameters_.applicability().assess(temperature_k, pressure_pa) ==
            RangeAssessment::outside_declared_bounds) {
            throw std::domain_error("CPA PT: state outside declared parameter applicability");
        }

        CpaPtRootSet result;
        result.options = options;
        result.pressure_pa = pressure_pa;
        result.temperature_k = temperature_k;
        result.composition.assign(composition.begin(), composition.end());

        const double b = cpa_detail::cpa_b_mix(composition, parameters_);
        const double pressure_tolerance = std::max(
            options.pressure_absolute_tolerance_pa,
            options.pressure_relative_tolerance * pressure_pa);
        constexpr double endpoint_guard =
            1024.0 * std::numeric_limits<double>::epsilon();
        const double u_max = 1.0 - endpoint_guard;

        auto evaluate_reduced_density = [&](double u, double& residual) -> bool {
            if (result.evaluations >= options.max_evaluations) {
                result.status = CpaPtRootStatus::evaluation_limit;
                result.diagnostic = "CPA PT: phase-property evaluation budget exhausted";
                return false;
            }
            ++result.evaluations;
            try {
                const double rho = u / b;
                const auto state = evaluate_cpa_phase_at_density(
                    temperature_k, rho, composition, parameters_, options.phase);
                residual = state.pressure_pa - pressure_pa;
                if (!std::isfinite(residual)) {
                    throw std::range_error("nonfinite pressure residual");
                }
                return true;
            } catch (const std::exception& error) {
                result.status = CpaPtRootStatus::property_failure;
                result.diagnostic = std::string("CPA PT: density-state property failure: ") +
                                    error.what();
                return false;
            }
        };

        std::vector<cpa_detail::CpaPtSample> samples;
        samples.reserve(options.scan_intervals + 1U);
        samples.push_back({0.0, -pressure_pa});
        for (std::size_t k = 1U; k <= options.scan_intervals; ++k) {
            const double angle = std::numbers::pi_v<double> *
                static_cast<double>(k) /
                (2.0 * static_cast<double>(options.scan_intervals));
            const double sine = std::sin(angle);
            const double u = u_max * sine * sine;
            double residual = 0.0;
            if (!evaluate_reduced_density(u, residual)) { return result; }
            samples.push_back({u, residual});
        }

        std::vector<cpa_detail::CpaPtRootCandidate> candidates;
        const auto append_candidate = [&](double u, double residual, int slope_sign) {
            if (candidates.size() >= options.max_roots) {
                result.status = CpaPtRootStatus::root_limit;
                result.diagnostic = "CPA PT: root-count quota exceeded";
                return false;
            }
            candidates.push_back({u, residual, slope_sign});
            return true;
        };

        for (std::size_t i = 0U; i + 1U < samples.size(); ++i) {
            const int left_sign = cpa_detail::cpa_pressure_sign(
                samples[i].residual_pa, pressure_tolerance);
            const int right_sign = cpa_detail::cpa_pressure_sign(
                samples[i + 1U].residual_pa, pressure_tolerance);
            if (left_sign == 0 || right_sign == 0 || left_sign == right_sign) {
                continue;
            }

            double left_u = samples[i].reduced_density;
            double right_u = samples[i + 1U].reduced_density;
            const double initial_left_r = samples[i].residual_pa;
            const double initial_right_r = samples[i + 1U].residual_pa;
            double best_u = std::abs(initial_left_r) < std::abs(initial_right_r)
                ? left_u : right_u;
            double best_r = std::abs(initial_left_r) < std::abs(initial_right_r)
                ? initial_left_r : initial_right_r;
            bool converged = false;
            for (int iteration = 0; iteration < options.max_bisection_iterations;
                 ++iteration) {
                const double mid_u = 0.5 * (left_u + right_u);
                if (mid_u == left_u || mid_u == right_u) {
                    converged = std::abs(best_r) <= pressure_tolerance;
                    break;
                }
                double mid_r = 0.0;
                if (!evaluate_reduced_density(mid_u, mid_r)) { return result; }
                if (std::abs(mid_r) < std::abs(best_r)) {
                    best_u = mid_u;
                    best_r = mid_r;
                }
                const int mid_sign =
                    cpa_detail::cpa_pressure_sign(mid_r, pressure_tolerance);
                if (mid_sign == 0) {
                    best_u = mid_u;
                    best_r = mid_r;
                    converged = true;
                    break;
                }
                if (mid_sign == left_sign) {
                    left_u = mid_u;
                } else {
                    right_u = mid_u;
                }
            }
            if (!converged && std::abs(best_r) > pressure_tolerance) {
                result.status = CpaPtRootStatus::iteration_limit;
                result.diagnostic = "CPA PT: bracketed density root did not reach pressure tolerance";
                return result;
            }
            const int slope_sign = left_sign < right_sign ? 1 : -1;
            if (!append_candidate(best_u, best_r, slope_sign)) { return result; }
        }

        // A scan point already within pressure tolerance is handled separately.
        // Opposite signs around the zero group indicate a resolved simple root;
        // equal signs indicate a tangent/near-multiple topology that is retained
        // explicitly instead of being published as a smooth root.
        for (std::size_t begin = 1U; begin < samples.size();) {
            if (cpa_detail::cpa_pressure_sign(
                    samples[begin].residual_pa, pressure_tolerance) != 0) {
                ++begin;
                continue;
            }
            std::size_t end = begin;
            std::size_t best = begin;
            while (end + 1U < samples.size() &&
                   cpa_detail::cpa_pressure_sign(
                       samples[end + 1U].residual_pa, pressure_tolerance) == 0) {
                ++end;
                if (std::abs(samples[end].residual_pa) <
                    std::abs(samples[best].residual_pa)) {
                    best = end;
                }
            }
            const int left_sign = begin > 0U
                ? cpa_detail::cpa_pressure_sign(
                      samples[begin - 1U].residual_pa, pressure_tolerance)
                : 0;
            const int right_sign = end + 1U < samples.size()
                ? cpa_detail::cpa_pressure_sign(
                      samples[end + 1U].residual_pa, pressure_tolerance)
                : 0;
            int slope_sign = 0;
            if (left_sign != 0 && right_sign != 0 && left_sign != right_sign) {
                slope_sign = left_sign < right_sign ? 1 : -1;
            }
            if (!append_candidate(
                    samples[best].reduced_density,
                    samples[best].residual_pa, slope_sign)) {
                return result;
            }
            begin = end + 1U;
        }

        if (candidates.empty()) {
            result.status = CpaPtRootStatus::no_root;
            result.diagnostic =
                "CPA PT: finite density scan found no representable pressure root";
            return result;
        }

        std::sort(candidates.begin(), candidates.end(),
                  [](const auto& first, const auto& second) {
                      return first.reduced_density < second.reduced_density;
                  });
        std::vector<cpa_detail::CpaPtRootCandidate> unique;
        unique.reserve(candidates.size());
        constexpr double root_merge_factor = 4096.0;
        for (const auto& candidate : candidates) {
            if (!unique.empty()) {
                const double scale = std::max(
                    {1.0, std::abs(unique.back().reduced_density),
                     std::abs(candidate.reduced_density)});
                if (std::abs(candidate.reduced_density -
                             unique.back().reduced_density) <=
                    root_merge_factor * std::numeric_limits<double>::epsilon() * scale) {
                    if (std::abs(candidate.residual_pa) <
                        std::abs(unique.back().residual_pa)) {
                        unique.back() = candidate;
                    }
                    continue;
                }
            }
            unique.push_back(candidate);
        }

        bool near_multiple = false;
        result.roots.reserve(unique.size());
        for (const auto& candidate : unique) {
            if (candidate.slope_sign == 0) { near_multiple = true; }
            if (result.evaluations >= options.max_evaluations) {
                result.status = CpaPtRootStatus::evaluation_limit;
                result.diagnostic = "CPA PT: phase-property evaluation budget exhausted";
                result.roots.clear();
                return result;
            }
            ++result.evaluations;
            try {
                const double rho = candidate.reduced_density / b;
                const auto state = evaluate_cpa_phase_at_density(
                    temperature_k, rho, composition, parameters_, options.phase);
                CpaPtRoot root;
                root.molar_density_mol_per_m3 = rho;
                root.pressure_residual_pa = state.pressure_pa - pressure_pa;
                if (!std::isfinite(root.pressure_residual_pa) ||
                    std::abs(root.pressure_residual_pa) > pressure_tolerance) {
                    result.status = CpaPtRootStatus::iteration_limit;
                    result.diagnostic =
                        "CPA PT: rebuilt density root no longer satisfies pressure tolerance";
                    result.roots.clear();
                    return result;
                }
                root.pressure_slope_sign = candidate.slope_sign;
                cpa_detail::cpa_fill_ln_phi(
                    pressure_pa, temperature_k, composition,
                    parameters_, state, root);
                result.roots.push_back(std::move(root));
            } catch (const std::exception& error) {
                result.status = CpaPtRootStatus::property_failure;
                result.diagnostic = std::string("CPA PT: accepted-root property failure: ") +
                                    error.what();
                result.roots.clear();
                return result;
            }
        }

        result.status = near_multiple
            ? CpaPtRootStatus::near_multiple
            : CpaPtRootStatus::success;
        result.diagnostic = near_multiple
            ? "CPA PT: finite density scan retained an unresolved tangent/near-multiple root"
            : "CPA PT: finite density scan resolved all detected simple pressure roots";
        return result;
    }

private:
    explicit CpaPtPhase(const CpaParameterSet& parameters)
        : parameters_(parameters) {}

    CpaParameterSet parameters_;
};

} // namespace mpmc::thermodynamics

#endif // MPMC_THERMODYNAMICS_CPA_PT_PHASE_HPP
