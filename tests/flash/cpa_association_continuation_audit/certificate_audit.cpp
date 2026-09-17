#include <mpmc/flash/cpa_split.hpp>

#include "certificate.hpp"
#include "contract.hpp"
#include "test_support.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <numbers>
#include <sstream>
#include <span>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace {
namespace fl = mpmc::flash;
namespace th = mpmc::thermodynamics;
namespace certificate = cpa_association_observable_certificate;
namespace contract = cpa_association_continuation_contract;

void require(bool condition, const std::string& message) {
    if (!condition) { throw std::runtime_error(message); }
}

std::size_t association_sweeps(const th::CpaAssociationResult& result) {
    if (result.status == th::CpaAssociationStatus::no_associating_sites) { return 0U; }
    if (!result.converged()) { return 0U; }
    return static_cast<std::size_t>(result.iterations) + 1U;
}

std::vector<double> site_fractions(const th::CpaAssociationResult& result) {
    std::vector<double> values;
    values.reserve(result.sites.size());
    for (const auto& site : result.sites) { values.push_back(site.unbonded_fraction); }
    return values;
}

struct ShadowAssociationSolve {
    th::CpaAssociationResult result;
    std::size_t sweeps{};
};

ShadowAssociationSolve solve_association_shadow(
    double temperature_k,
    double molar_density_mol_per_m3,
    std::span<const double> composition,
    const th::CpaParameterSet& parameters,
    const th::CpaAssociationOptions& options,
    std::span<const double> initial_site_fractions) {
    th::cpa_detail::validate_association_options(options);
    (void)th::cpa_detail::validate_cpa_composition(composition, parameters.size());

    ShadowAssociationSolve shadow;
    auto& result = shadow.result;
    result.temperature_k = temperature_k;
    result.molar_density_mol_per_m3 = molar_density_mol_per_m3;
    result.b_mix_m3_per_mol = th::cpa_detail::cpa_b_mix(composition, parameters);
    result.eta = 0.25 * result.b_mix_m3_per_mol * molar_density_mol_per_m3;
    const double denominator = 1.0 - 1.9 * result.eta;
    if (!std::isfinite(denominator) || !(denominator > 0.0)) {
        result.status = th::CpaAssociationStatus::radial_distribution_singularity;
        return shadow;
    }
    result.radial_distribution = 1.0 / denominator;
    result.rho_dln_g_drho = (1.9 * result.eta) / denominator;

    for (std::size_t component = 0U; component < parameters.size(); ++component) {
        for (const auto& site : parameters.pure(component).sites) {
            result.sites.push_back({component, site.id, site.multiplicity, 1.0});
        }
    }
    if (result.sites.empty()) {
        result.status = th::CpaAssociationStatus::no_associating_sites;
        return shadow;
    }

    std::vector<double> current(result.sites.size(), 1.0);
    require(initial_site_fractions.empty() ||
                initial_site_fractions.size() == current.size(),
            "certificate audit warm seed dimension mismatch");
    if (!initial_site_fractions.empty()) {
        for (std::size_t i = 0U; i < current.size(); ++i) {
            require(std::isfinite(initial_site_fractions[i]) &&
                        initial_site_fractions[i] > 0.0 &&
                        initial_site_fractions[i] <= 1.0,
                    "certificate audit warm seed left (0,1]");
            current[i] = initial_site_fractions[i];
        }
    }

    std::vector<double> target(result.sites.size(), 1.0);
    for (int iteration = 0; iteration <= options.max_iterations; ++iteration) {
        ++shadow.sweeps;
        double residual = 0.0;
        for (std::size_t a = 0U; a < result.sites.size(); ++a) {
            const auto& first = result.sites[a];
            double association_sum = 0.0;
            for (std::size_t b = 0U; b < result.sites.size(); ++b) {
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
                const double delta = th::cpa_detail::cpa_association_strength(
                    temperature_k, result.radial_distribution, b_ij, *pair);
                association_sum += composition[second.component_index] *
                    static_cast<double>(second.multiplicity) * current[b] * delta;
            }
            const double value = 1.0 /
                (1.0 + molar_density_mol_per_m3 * association_sum);
            if (!std::isfinite(value) || !(value > 0.0) || value > 1.0) {
                result.status = th::CpaAssociationStatus::nonrepresentable_strength;
                return shadow;
            }
            target[a] = value;
            residual = std::max(residual, std::abs(target[a] - current[a]));
        }
        result.iterations = iteration;
        result.max_fixed_point_residual = residual;
        if (residual <= options.site_fraction_tolerance) {
            for (std::size_t i = 0U; i < result.sites.size(); ++i) {
                result.sites[i].unbonded_fraction = target[i];
            }
            result.status = th::CpaAssociationStatus::success;
            return shadow;
        }
        if (iteration == options.max_iterations) { break; }
        for (std::size_t i = 0U; i < current.size(); ++i) {
            current[i] += options.damping * (target[i] - current[i]);
        }
    }
    result.status = th::CpaAssociationStatus::iteration_limit;
    return shadow;
}

th::CpaPhaseState phase_with_association(
    const th::CpaPhaseState& cold,
    std::span<const double> composition,
    const th::CpaAssociationResult& association) {
    th::CpaPhaseState result = cold;
    result.association = association;
    const double rho = result.molar_density_mol_per_m3;
    const double rt = th::cpa_gas_constant_j_per_mol_k * result.temperature_k;
    double association_sum = 0.0;
    for (const auto& site : association.sites) {
        association_sum += composition[site.component_index] *
            static_cast<double>(site.multiplicity) *
            (1.0 - site.unbonded_fraction);
    }
    result.pressure_association_pa = -0.5 * rt * rho *
        (1.0 + association.rho_dln_g_drho) * association_sum;
    result.pressure_pa = result.pressure_physical_pa + result.pressure_association_pa;
    return result;
}

std::vector<double> ln_phi_at_density(
    double target_pressure_pa,
    double temperature_k,
    std::span<const double> composition,
    const th::CpaParameterSet& parameters,
    const th::CpaPhaseState& phase) {
    th::CpaPtRoot root;
    root.molar_density_mol_per_m3 = phase.molar_density_mol_per_m3;
    th::cpa_detail::cpa_fill_ln_phi(
        target_pressure_pa, temperature_k, composition,
        parameters, phase, root);
    return root.ln_phi;
}

struct Request {
    double pressure_pa{};
    double temperature_k{};
    std::vector<double> composition;
};

class RecordingDelegate {
public:
    RecordingDelegate(
        const th::CpaPtPhase& model,
        th::CpaPtOptions options,
        std::vector<Request>& requests)
        : delegate_(model, std::move(options)), requests_(requests) {}

    fl::StabilityPhase operator()(
        double pressure_pa, double temperature_k,
        std::span<const double> composition) {
        record(pressure_pa, temperature_k, composition);
        return delegate_(pressure_pa, temperature_k, composition);
    }

    fl::PtSplitPhase operator()(
        double pressure_pa, double temperature_k,
        std::span<const double> composition,
        fl::PtPhaseRole role) {
        record(pressure_pa, temperature_k, composition);
        return delegate_(pressure_pa, temperature_k, composition, role);
    }

private:
    void record(
        double pressure_pa, double temperature_k,
        std::span<const double> composition) {
        requests_.push_back({
            pressure_pa, temperature_k,
            std::vector<double>(composition.begin(), composition.end())});
    }

    fl::CpaVleEvaluator delegate_;
    std::vector<Request>& requests_;
};

struct AuditStats {
    std::size_t requests{};
    std::size_t density_states{};
    std::size_t seeded_states{};
    std::size_t warm_converged{};
    std::size_t enclosure_certified{};
    std::size_t pressure_certified{};
    std::size_t ln_phi_certified{};
    std::size_t observable_certified{};
    std::size_t sign_certified{};
    std::size_t root_safe_certified{};
    std::size_t enclosure_failures{};
    std::size_t pressure_bound_failures{};
    std::size_t ln_phi_bound_failures{};
    std::size_t sign_ambiguities{};
    std::size_t cold_enclosure_certified{};
    std::size_t empirical_site_bound_checks{};
    std::size_t empirical_pressure_bound_checks{};
    std::size_t empirical_ln_phi_bound_checks{};
    std::size_t effective_sweeps{};
    std::size_t cold_sweeps{};
    std::size_t warm_sweeps{};
    double max_fresh_residual{};
    double max_site_error_bound{};
    double max_pressure_error_bound_pa{};
    double max_ln_phi_error_bound{};
    double max_actual_site_diff{};
    double max_actual_pressure_diff_pa{};
    double max_actual_ln_phi_diff{};
};

bool pressure_sign_certified(
    double warm_residual_pa,
    double pressure_error_bound_pa,
    double pressure_tolerance_pa) {
    const double lower = warm_residual_pa - pressure_error_bound_pa;
    const double upper = warm_residual_pa + pressure_error_bound_pa;
    if (lower > pressure_tolerance_pa) { return true; }
    if (upper < -pressure_tolerance_pa) { return true; }
    return lower >= -pressure_tolerance_pa &&
           upper <= pressure_tolerance_pa;
}

void validate_empirical_bounds(
    const th::CpaAssociationResult& warm,
    const th::CpaAssociationResult& cold,
    const certificate::Result& warm_certificate,
    const certificate::Result& cold_certificate,
    double actual_pressure_diff,
    double actual_ln_phi_diff,
    AuditStats& stats) {
    if (!warm_certificate.enclosure_certified ||
        !cold_certificate.enclosure_certified) {
        return;
    }
    require(warm.sites.size() == cold.sites.size(),
            "certificate validation site dimension mismatch");
    require(warm_certificate.site_error_bounds.size() == warm.sites.size() &&
                cold_certificate.site_error_bounds.size() == cold.sites.size(),
            "certificate validation radius dimension mismatch");
    for (std::size_t i = 0U; i < warm.sites.size(); ++i) {
        const double actual = std::abs(
            warm.sites[i].unbonded_fraction - cold.sites[i].unbonded_fraction);
        const double bound = warm_certificate.site_error_bounds[i] +
            cold_certificate.site_error_bounds[i];
        if (!(actual <= bound + 1.0e-15)) {
            std::ostringstream message;
            message << std::setprecision(17)
                    << "warm/cold site difference exceeded combined certificate bounds"
                    << " site=" << i
                    << " actual=" << actual
                    << " bound=" << bound
                    << " warm_radius=" << warm_certificate.site_error_bounds[i]
                    << " cold_radius=" << cold_certificate.site_error_bounds[i]
                    << " warm_residual=" << warm_certificate.fresh_fixed_point_residual_inf
                    << " cold_residual=" << cold_certificate.fresh_fixed_point_residual_inf
                    << " warm_ratio=" << warm_certificate.krawczyk_max_inclusion_ratio
                    << " cold_ratio=" << cold_certificate.krawczyk_max_inclusion_ratio;
            throw std::runtime_error(message.str());
        }
    }
    ++stats.empirical_site_bound_checks;

    const double pressure_bound = warm_certificate.pressure_error_bound_pa +
        cold_certificate.pressure_error_bound_pa;
    if (!(actual_pressure_diff <= pressure_bound + 1.0e-9)) {
        std::ostringstream message;
        message << std::setprecision(17)
                << "warm/cold pressure difference exceeded combined certificate bounds"
                << " actual=" << actual_pressure_diff
                << " bound=" << pressure_bound
                << " warm_bound=" << warm_certificate.pressure_error_bound_pa
                << " cold_bound=" << cold_certificate.pressure_error_bound_pa
                << " warm_max_site_radius=" << warm_certificate.max_site_error_bound
                << " cold_max_site_radius=" << cold_certificate.max_site_error_bound
                << " warm_residual=" << warm_certificate.fresh_fixed_point_residual_inf
                << " cold_residual=" << cold_certificate.fresh_fixed_point_residual_inf
                << " warm_ratio=" << warm_certificate.krawczyk_max_inclusion_ratio
                << " cold_ratio=" << cold_certificate.krawczyk_max_inclusion_ratio;
        throw std::runtime_error(message.str());
    }
    ++stats.empirical_pressure_bound_checks;

    const double ln_phi_bound = warm_certificate.max_ln_phi_error_bound +
        cold_certificate.max_ln_phi_error_bound;
    if (!(actual_ln_phi_diff <= ln_phi_bound + 1.0e-13)) {
        std::ostringstream message;
        message << std::setprecision(17)
                << "warm/cold ln(phi) difference exceeded combined certificate bounds"
                << " actual=" << actual_ln_phi_diff
                << " bound=" << ln_phi_bound
                << " warm_bound=" << warm_certificate.max_ln_phi_error_bound
                << " cold_bound=" << cold_certificate.max_ln_phi_error_bound;
        throw std::runtime_error(message.str());
    }
    ++stats.empirical_ln_phi_bound_checks;
}

fl::PtSplitOptions audit_split_options() {
    fl::PtSplitOptions options;
    options.initial_stability.automatic_starts = false;
    options.final_stability.automatic_starts = false;
    return options;
}

} // namespace

int main() {
    try {
        const auto parameters = cpa_physical_test::parameters(false);
        const auto model = th::CpaPtPhase::from_parameters(parameters);
        const auto phase_options = fl::cpa_pt_vle_default_phase_options();
        std::vector<Request> requests;

        for (const auto& state : cpa_physical_test::points()) {
            const auto feed = cpa_physical_test::feed(state, false);
            const auto starts = cpa_physical_test::starts(state, false);
            const std::vector<std::vector<double>> final_starts{feed};
            RecordingDelegate evaluator(model, phase_options, requests);
            const auto result = fl::solve_pt_vle(
                state.pressure_pa, cpa_physical_test::temperature_k,
                feed, evaluator, evaluator,
                audit_split_options(), starts, final_starts);
            require(result.status == fl::PtSplitStatus::two_phase_no_instability_found,
                    "certificate audit physical flash did not close two phase");
        }
        require(requests.size() == 497U,
                "certificate audit no longer sees the frozen 497 root requests");

        AuditStats stats;
        stats.requests = requests.size();
        constexpr double endpoint_guard =
            1024.0 * std::numeric_limits<double>::epsilon();
        const double u_max = 1.0 - endpoint_guard;

        for (const auto& request : requests) {
            const auto composition = std::span<const double>(request.composition);
            const double b_mix = th::cpa_detail::cpa_b_mix(composition, parameters);
            const double pressure_tolerance = std::max(
                phase_options.pressure_absolute_tolerance_pa,
                phase_options.pressure_relative_tolerance * request.pressure_pa);
            std::map<double, std::vector<double>> cache;

            for (std::size_t k = 1U; k <= phase_options.scan_intervals; ++k) {
                ++stats.density_states;
                const double angle = std::numbers::pi_v<double> *
                    static_cast<double>(k) /
                    (2.0 * static_cast<double>(phase_options.scan_intervals));
                const double sine = std::sin(angle);
                const double u = u_max * sine * sine;
                const double rho = u / b_mix;

                const auto cold = th::evaluate_cpa_phase_at_density(
                    request.temperature_k, rho, composition,
                    parameters, phase_options.phase);
                const std::size_t cold_sweep_count = association_sweeps(cold.association);
                stats.cold_sweeps += cold_sweep_count;

                if (cache.empty()) {
                    cache[u] = site_fractions(cold.association);
                    stats.effective_sweeps += cold_sweep_count;
                    continue;
                }
                ++stats.seeded_states;

                auto upper = cache.lower_bound(u);
                auto best = upper;
                if (upper == cache.end()) {
                    best = std::prev(cache.end());
                } else if (upper != cache.begin()) {
                    const auto lower = std::prev(upper);
                    if (std::abs(lower->first - u) <= std::abs(upper->first - u)) {
                        best = lower;
                    }
                }
                const auto warm = solve_association_shadow(
                    request.temperature_k, rho, composition, parameters,
                    phase_options.phase.association,
                    std::span<const double>(best->second));
                stats.warm_sweeps += warm.sweeps;
                require(warm.result.converged(),
                        "certificate audit warm solve failed to converge");
                ++stats.warm_converged;

                const auto warm_certificate = certificate::certify(
                    request.temperature_k, rho, composition, parameters, warm.result);
                const auto cold_certificate = certificate::certify(
                    request.temperature_k, rho, composition, parameters, cold.association);
                if (cold_certificate.enclosure_certified) {
                    ++stats.cold_enclosure_certified;
                }

                stats.max_fresh_residual = std::max(
                    stats.max_fresh_residual,
                    warm_certificate.fresh_fixed_point_residual_inf);
                stats.max_site_error_bound = std::max(
                    stats.max_site_error_bound,
                    warm_certificate.max_site_error_bound);
                stats.max_pressure_error_bound_pa = std::max(
                    stats.max_pressure_error_bound_pa,
                    warm_certificate.pressure_error_bound_pa);
                stats.max_ln_phi_error_bound = std::max(
                    stats.max_ln_phi_error_bound,
                    warm_certificate.max_ln_phi_error_bound);

                if (warm_certificate.enclosure_certified) {
                    ++stats.enclosure_certified;
                } else {
                    ++stats.enclosure_failures;
                }
                if (warm_certificate.pressure_scale_certified) {
                    ++stats.pressure_certified;
                } else if (warm_certificate.enclosure_certified) {
                    ++stats.pressure_bound_failures;
                }
                if (warm_certificate.ln_phi_scale_certified) {
                    ++stats.ln_phi_certified;
                } else if (warm_certificate.enclosure_certified) {
                    ++stats.ln_phi_bound_failures;
                }
                if (warm_certificate.observable_certified()) {
                    ++stats.observable_certified;
                }

                const auto warm_phase = phase_with_association(
                    cold, composition, warm.result);
                const double warm_residual = warm_phase.pressure_pa - request.pressure_pa;
                const bool sign_ok = warm_certificate.enclosure_certified &&
                    pressure_sign_certified(
                        warm_residual,
                        warm_certificate.pressure_error_bound_pa,
                        pressure_tolerance);
                if (sign_ok) {
                    ++stats.sign_certified;
                } else if (warm_certificate.observable_certified()) {
                    ++stats.sign_ambiguities;
                }
                const bool root_safe = warm_certificate.observable_certified() && sign_ok;
                if (root_safe) { ++stats.root_safe_certified; }

                const double actual_pressure_diff =
                    std::abs(warm_phase.pressure_pa - cold.pressure_pa);
                stats.max_actual_pressure_diff_pa = std::max(
                    stats.max_actual_pressure_diff_pa, actual_pressure_diff);
                double actual_site_diff = 0.0;
                for (std::size_t i = 0U; i < warm.result.sites.size(); ++i) {
                    actual_site_diff = std::max(
                        actual_site_diff,
                        std::abs(warm.result.sites[i].unbonded_fraction -
                                 cold.association.sites[i].unbonded_fraction));
                }
                stats.max_actual_site_diff = std::max(
                    stats.max_actual_site_diff, actual_site_diff);

                const auto warm_ln_phi = ln_phi_at_density(
                    request.pressure_pa, request.temperature_k,
                    composition, parameters, warm_phase);
                const auto cold_ln_phi = ln_phi_at_density(
                    request.pressure_pa, request.temperature_k,
                    composition, parameters, cold);
                double actual_ln_phi_diff = 0.0;
                for (std::size_t i = 0U; i < warm_ln_phi.size(); ++i) {
                    actual_ln_phi_diff = std::max(
                        actual_ln_phi_diff,
                        std::abs(warm_ln_phi[i] - cold_ln_phi[i]));
                }
                stats.max_actual_ln_phi_diff = std::max(
                    stats.max_actual_ln_phi_diff, actual_ln_phi_diff);

                validate_empirical_bounds(
                    warm.result, cold.association,
                    warm_certificate, cold_certificate,
                    actual_pressure_diff, actual_ln_phi_diff, stats);

                if (root_safe) {
                    cache[u] = site_fractions(warm.result);
                    stats.effective_sweeps += warm.sweeps;
                } else {
                    cache[u] = site_fractions(cold.association);
                    stats.effective_sweeps += warm.sweeps + cold_sweep_count;
                }
            }
        }

        const double seeded = static_cast<double>(stats.seeded_states);
        const double observable_rate = seeded > 0.0
            ? static_cast<double>(stats.observable_certified) / seeded
            : 0.0;
        const double root_safe_rate = seeded > 0.0
            ? static_cast<double>(stats.root_safe_certified) / seeded
            : 0.0;
        const double sweep_reduction = stats.cold_sweeps > 0U
            ? 1.0 - static_cast<double>(stats.effective_sweeps) /
                        static_cast<double>(stats.cold_sweeps)
            : 0.0;

        std::cout << std::setprecision(17)
            << "CPA_ASSOC_NO_COLD_CERTIFICATE_SUMMARY"
            << " requests=" << stats.requests
            << " density_states=" << stats.density_states
            << " seeded_states=" << stats.seeded_states
            << " warm_converged=" << stats.warm_converged
            << " enclosure_certified=" << stats.enclosure_certified
            << " pressure_certified=" << stats.pressure_certified
            << " lnphi_certified=" << stats.ln_phi_certified
            << " observable_certified=" << stats.observable_certified
            << " observable_certified_rate=" << observable_rate
            << " sign_certified=" << stats.sign_certified
            << " sign_ambiguities=" << stats.sign_ambiguities
            << " root_safe_certified=" << stats.root_safe_certified
            << " root_safe_certified_rate=" << root_safe_rate
            << " enclosure_failures=" << stats.enclosure_failures
            << " pressure_bound_failures=" << stats.pressure_bound_failures
            << " lnphi_bound_failures=" << stats.ln_phi_bound_failures
            << " cold_enclosure_certified=" << stats.cold_enclosure_certified
            << " empirical_site_bound_checks=" << stats.empirical_site_bound_checks
            << " empirical_pressure_bound_checks=" << stats.empirical_pressure_bound_checks
            << " empirical_lnphi_bound_checks=" << stats.empirical_ln_phi_bound_checks
            << " cold_sweeps=" << stats.cold_sweeps
            << " warm_sweeps=" << stats.warm_sweeps
            << " effective_sweeps=" << stats.effective_sweeps
            << " sweep_reduction_fraction=" << sweep_reduction
            << " max_fresh_residual=" << stats.max_fresh_residual
            << " max_site_error_bound=" << stats.max_site_error_bound
            << " max_pressure_error_bound_pa=" << stats.max_pressure_error_bound_pa
            << " max_lnphi_error_bound=" << stats.max_ln_phi_error_bound
            << " max_actual_dX=" << stats.max_actual_site_diff
            << " max_actual_dP_pa=" << stats.max_actual_pressure_diff_pa
            << " max_actual_dlnphi=" << stats.max_actual_ln_phi_diff
            << " pressure_guard_pa=" << contract::pressure_error_guard_pa
            << " lnphi_guard=" << contract::ln_phi_error_guard
            << " cold_reference_used_for_acceptance=false"
            << " production_solver_changes=none"
            << '\n';
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "[FAIL] " << error.what() << '\n';
        return 1;
    }
}
