#include <mpmc/flash/cpa_split.hpp>

#include "cpa_thermopack_parity_thresholds.hpp"
#include "test_support.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <optional>
#include <span>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace {
namespace fl = mpmc::flash;
namespace th = mpmc::thermodynamics;
namespace parity_gate = cpa_thermopack_parity_thresholds;

constexpr double site_fraction_audit_guard = 1.0e-10;
constexpr double replay_relative_density_guard = 2.0e-12;
constexpr double replay_abs_ln_phi_guard = 2.0e-12;

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
    if (!std::isfinite(temperature_k) || !(temperature_k > 0.0) ||
        !std::isfinite(molar_density_mol_per_m3) ||
        !(molar_density_mol_per_m3 > 0.0)) {
        throw std::domain_error("CPA continuation audit: positive finite T and density required");
    }

    ShadowAssociationSolve shadow;
    auto& result = shadow.result;
    result.temperature_k = temperature_k;
    result.molar_density_mol_per_m3 = molar_density_mol_per_m3;
    result.b_mix_m3_per_mol = th::cpa_detail::cpa_b_mix(composition, parameters);
    result.eta = 0.25 * result.b_mix_m3_per_mol * molar_density_mol_per_m3;
    const double denominator = 1.0 - 1.9 * result.eta;
    if (!std::isfinite(denominator) || !(denominator > 0.0)) {
        result.status = th::CpaAssociationStatus::radial_distribution_singularity;
        result.diagnostic = "shadow association: simplified radial-distribution denominator is nonpositive";
        return shadow;
    }
    result.radial_distribution = 1.0 / denominator;
    result.rho_dln_g_drho = (1.9 * result.eta) / denominator;

    for (std::size_t component = 0U; component < parameters.size(); ++component) {
        for (const auto& site : parameters.pure(component).sites) {
            if (result.sites.size() >= options.max_site_classes) {
                throw std::length_error("CPA continuation audit: site-class quota exceeded");
            }
            result.sites.push_back({component, site.id, site.multiplicity, 1.0});
        }
    }
    if (result.sites.empty()) {
        result.status = th::CpaAssociationStatus::no_associating_sites;
        result.diagnostic = "shadow association: no associating sites";
        return shadow;
    }

    std::vector<double> current(result.sites.size(), 1.0);
    if (!initial_site_fractions.empty()) {
        if (initial_site_fractions.size() != current.size()) {
            throw std::invalid_argument("CPA continuation audit: warm seed dimension mismatch");
        }
        for (std::size_t i = 0U; i < current.size(); ++i) {
            const double value = initial_site_fractions[i];
            if (!std::isfinite(value) || !(value > 0.0) || value > 1.0) {
                throw std::domain_error("CPA continuation audit: warm seed must lie in (0,1]");
            }
            current[i] = value;
        }
    }

    std::vector<double> target(result.sites.size(), 1.0);
    for (int iteration = 0; iteration <= options.max_iterations; ++iteration) {
        ++shadow.sweeps;
        double residual = 0.0;
        try {
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
                    result.diagnostic = "shadow association: site fraction left (0,1]";
                    return shadow;
                }
                target[a] = value;
                residual = std::max(residual, std::abs(target[a] - current[a]));
            }
        } catch (const std::range_error& error) {
            result.status = th::CpaAssociationStatus::nonrepresentable_strength;
            result.diagnostic = error.what();
            return shadow;
        }

        result.iterations = iteration;
        result.max_fixed_point_residual = residual;
        if (residual <= options.site_fraction_tolerance) {
            for (std::size_t i = 0U; i < result.sites.size(); ++i) {
                result.sites[i].unbonded_fraction = target[i];
            }
            result.status = th::CpaAssociationStatus::success;
            result.diagnostic = "shadow association site fractions converged";
            return shadow;
        }
        if (iteration == options.max_iterations) { break; }
        for (std::size_t i = 0U; i < current.size(); ++i) {
            current[i] += options.damping * (target[i] - current[i]);
        }
    }

    result.status = th::CpaAssociationStatus::iteration_limit;
    result.diagnostic = "shadow association: fixed-point iteration limit";
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
    for (const auto& site : result.association.sites) {
        association_sum += composition[site.component_index] *
            static_cast<double>(site.multiplicity) *
            (1.0 - site.unbonded_fraction);
    }
    result.pressure_association_pa = -0.5 * rt * rho *
        (1.0 + result.association.rho_dln_g_drho) * association_sum;
    result.pressure_pa = result.pressure_physical_pa + result.pressure_association_pa;
    return result;
}

struct SeedSelection {
    std::vector<double> values;
    bool exact{};
};

std::optional<SeedSelection> nearest_seed(
    const std::map<double, std::vector<double>>& cache, double reduced_density) {
    if (cache.empty()) { return std::nullopt; }
    auto upper = cache.lower_bound(reduced_density);
    auto best = upper;
    if (upper == cache.end()) {
        best = std::prev(cache.end());
    } else if (upper != cache.begin()) {
        const auto lower = std::prev(upper);
        if (std::abs(lower->first - reduced_density) <=
            std::abs(upper->first - reduced_density)) {
            best = lower;
        }
    }
    return SeedSelection{best->second, best->first == reduced_density};
}

struct ContinuationStats {
    std::size_t root_searches{};
    std::size_t density_evaluations{};
    std::size_t unseeded_evaluations{};
    std::size_t seeded_evaluations{};
    std::size_t exact_seed_hits{};
    std::size_t accepted_continuations{};
    std::size_t fallback_evaluations{};
    std::size_t nonconverged_fallbacks{};
    std::size_t site_guard_fallbacks{};
    std::size_t pressure_guard_fallbacks{};
    std::size_t root_lnphi_fallbacks{};
    std::size_t warm_better{};
    std::size_t warm_equal{};
    std::size_t warm_worse{};
    std::size_t cold_sweeps{};
    std::size_t seeded_warm_attempt_sweeps{};
    std::size_t effective_sweeps_with_fallback{};
    std::size_t max_cold_sweeps{};
    std::size_t max_warm_sweeps{};
    double max_raw_site_fraction_diff{};
    double max_raw_pressure_diff_pa{};
    double max_raw_root_ln_phi_diff{};
};

struct AuditEvaluation {
    double pressure_residual_pa{};
    th::CpaPhaseState cold;
    th::CpaPhaseState effective;
    bool had_seed{};
    bool accepted_continuation{};
    std::size_t cold_sweeps{};
};

AuditEvaluation audit_density_state(
    double target_pressure_pa,
    double temperature_k,
    double reduced_density,
    double b_mix,
    std::span<const double> composition,
    const th::CpaParameterSet& parameters,
    const th::CpaPhaseOptions& options,
    std::map<double, std::vector<double>>& cache,
    ContinuationStats& stats) {
    ++stats.density_evaluations;
    const double rho = reduced_density / b_mix;
    const auto cold = th::evaluate_cpa_phase_at_density(
        temperature_k, rho, composition, parameters, options);
    const std::size_t cold_sweep_count = association_sweeps(cold.association);
    stats.cold_sweeps += cold_sweep_count;
    stats.max_cold_sweeps = std::max(stats.max_cold_sweeps, cold_sweep_count);

    const auto seed = nearest_seed(cache, reduced_density);
    if (seed) {
        ++stats.seeded_evaluations;
        if (seed->exact) { ++stats.exact_seed_hits; }
    } else {
        ++stats.unseeded_evaluations;
    }

    const std::span<const double> seed_span = seed
        ? std::span<const double>(seed->values)
        : std::span<const double>();
    const auto warm = solve_association_shadow(
        temperature_k, rho, composition, parameters,
        options.association, seed_span);
    stats.max_warm_sweeps = std::max(stats.max_warm_sweeps, warm.sweeps);
    if (seed) { stats.seeded_warm_attempt_sweeps += warm.sweeps; }

    bool site_layout_match = warm.result.sites.size() == cold.association.sites.size();
    double max_site_diff = 0.0;
    if (site_layout_match) {
        for (std::size_t i = 0U; i < warm.result.sites.size(); ++i) {
            const auto& first = warm.result.sites[i];
            const auto& second = cold.association.sites[i];
            if (first.component_index != second.component_index ||
                first.site_id != second.site_id ||
                first.multiplicity != second.multiplicity) {
                site_layout_match = false;
                break;
            }
            max_site_diff = std::max(
                max_site_diff,
                std::abs(first.unbonded_fraction - second.unbonded_fraction));
        }
    }
    stats.max_raw_site_fraction_diff =
        std::max(stats.max_raw_site_fraction_diff, max_site_diff);

    const bool warm_converged = warm.result.converged() && site_layout_match;
    th::CpaPhaseState warm_phase = cold;
    double pressure_diff = std::numeric_limits<double>::infinity();
    if (warm_converged) {
        warm_phase = phase_with_association(cold, composition, warm.result);
        pressure_diff = std::abs(warm_phase.pressure_pa - cold.pressure_pa);
        stats.max_raw_pressure_diff_pa =
            std::max(stats.max_raw_pressure_diff_pa, pressure_diff);
    }

    const bool site_ok = warm_converged && max_site_diff <= site_fraction_audit_guard;
    const bool pressure_ok = warm_converged &&
        pressure_diff <= parity_gate::phase_max_abs_pressure_total_pa;
    const bool accepted = !seed || (site_ok && pressure_ok);

    th::CpaPhaseState effective = cold;
    if (!seed) {
        require(warm_converged, "cold-equivalent shadow association failed on first density state");
        require(site_ok && pressure_ok,
                "cold-equivalent shadow association disagreed with production cold solve");
        effective = std::move(warm_phase);
        stats.effective_sweeps_with_fallback += cold_sweep_count;
    } else if (accepted) {
        ++stats.accepted_continuations;
        effective = std::move(warm_phase);
        stats.effective_sweeps_with_fallback += warm.sweeps;
        if (warm.sweeps < cold_sweep_count) {
            ++stats.warm_better;
        } else if (warm.sweeps == cold_sweep_count) {
            ++stats.warm_equal;
        } else {
            ++stats.warm_worse;
        }
    } else {
        ++stats.fallback_evaluations;
        if (!warm_converged) { ++stats.nonconverged_fallbacks; }
        if (warm_converged && !site_ok) { ++stats.site_guard_fallbacks; }
        if (warm_converged && !pressure_ok) { ++stats.pressure_guard_fallbacks; }
        stats.effective_sweeps_with_fallback += warm.sweeps + cold_sweep_count;
    }

    cache[reduced_density] = site_fractions(effective.association);
    return {
        effective.pressure_pa - target_pressure_pa,
        std::move(const_cast<th::CpaPhaseState&>(cold)),
        std::move(effective),
        seed.has_value(),
        seed.has_value() && accepted,
        cold_sweep_count};
}

struct RootRecord {
    double pressure_pa{};
    double temperature_k{};
    std::vector<double> composition;
    th::CpaPtRootSet roots;
};

class RecordingEvaluator {
public:
    RecordingEvaluator(
        const th::CpaPtPhase& model,
        th::CpaPtOptions options,
        std::vector<RootRecord>& records)
        : model_(model), options_(std::move(options)), records_(records) {}

    fl::StabilityPhase operator()(
        double pressure_pa, double temperature_k,
        std::span<const double> composition) {
        const auto& roots = record_roots(pressure_pa, temperature_k, composition);
        fl::detail::cpa_flash_require_root_success(roots, "CPA continuation recording stability");
        const auto admissible = fl::detail::cpa_flash_admissible_roots(roots);

        const auto gibbs_offset = [&](std::size_t root_index) {
            const auto& root = roots.roots[root_index];
            double value = 0.0;
            double correction = 0.0;
            for (std::size_t i = 0U; i < composition.size(); ++i) {
                fl::detail::stability_add(
                    composition[i] * root.ln_phi[i], value, correction);
            }
            return value;
        };

        std::size_t best_position = 0U;
        double best_gibbs = gibbs_offset(admissible[0]);
        for (std::size_t position = 1U; position < admissible.size(); ++position) {
            const double value = gibbs_offset(admissible[position]);
            if (value < best_gibbs) {
                best_position = position;
                best_gibbs = value;
            }
        }
        const std::size_t best_root = admissible[best_position];
        fl::StabilityPhase phase{roots.roots[best_root].ln_phi, best_root, true};
        for (const std::size_t candidate : admissible) {
            if (candidate == best_root) { continue; }
            const double candidate_gibbs = gibbs_offset(candidate);
            double magnitude = 1.0;
            double correction = 0.0;
            for (std::size_t i = 0U; i < composition.size(); ++i) {
                fl::detail::stability_add(
                    composition[i] *
                        (std::abs(roots.roots[candidate].ln_phi[i]) +
                         std::abs(roots.roots[best_root].ln_phi[i])),
                    magnitude, correction);
            }
            if (std::abs(candidate_gibbs - best_gibbs) <=
                256.0 * fl::detail::stability_eps * magnitude) {
                phase.smooth = false;
            }
        }
        return phase;
    }

    fl::PtSplitPhase operator()(
        double pressure_pa, double temperature_k,
        std::span<const double> composition,
        fl::PtPhaseRole role) {
        const auto& roots = record_roots(pressure_pa, temperature_k, composition);
        fl::detail::cpa_flash_require_root_success(roots, "CPA continuation recording split");
        const auto admissible = fl::detail::cpa_flash_admissible_roots(roots);
        const auto side = role == fl::PtPhaseRole::liquid_candidate
            ? fl::CpaRootSide::upper_density_admissible
            : fl::CpaRootSide::lower_density_admissible;
        const std::size_t selected = fl::detail::cpa_flash_select_root_side(
            roots, admissible, side);
        const auto& root = roots.roots[selected];
        return {{root.ln_phi, selected, true}, root.compressibility_factor};
    }

private:
    const th::CpaPtRootSet& record_roots(
        double pressure_pa, double temperature_k,
        std::span<const double> composition) {
        RootRecord record;
        record.pressure_pa = pressure_pa;
        record.temperature_k = temperature_k;
        record.composition.assign(composition.begin(), composition.end());
        record.roots = model_.roots(
            pressure_pa, temperature_k, composition, options_);
        records_.push_back(std::move(record));
        return records_.back().roots;
    }

    const th::CpaPtPhase& model_;
    th::CpaPtOptions options_;
    std::vector<RootRecord>& records_;
};

void compare_outer_solution(
    const fl::PtSplitResult& production,
    const fl::PtSplitResult& recording) {
    require(production.status == recording.status,
            "recording provider changed PT split status");
    require(production.initial_stability.evaluations ==
                recording.initial_stability.evaluations,
            "recording provider changed initial stability evaluation count");
    require(production.split_evaluations == recording.split_evaluations,
            "recording provider changed split evaluation count");
    require(production.final_stability.has_value() ==
                recording.final_stability.has_value(),
            "recording provider changed final stability availability");
    if (production.final_stability && recording.final_stability) {
        require(production.final_stability->evaluations ==
                    recording.final_stability->evaluations,
                "recording provider changed final stability evaluation count");
    }
    const auto* first = production.candidate();
    const auto* second = recording.candidate();
    require((first != nullptr) == (second != nullptr),
            "recording provider changed candidate availability");
    if (first == nullptr || second == nullptr) { return; }
    require(std::abs(first->fractions.vapor_fraction -
                     second->fractions.vapor_fraction) <= 1.0e-14,
            "recording provider changed vapor fraction");
    for (std::size_t i = 0U; i < first->fractions.liquid.size(); ++i) {
        require(std::abs(first->fractions.liquid[i] -
                         second->fractions.liquid[i]) <= 1.0e-14,
                "recording provider changed liquid composition");
        require(std::abs(first->fractions.vapor[i] -
                         second->fractions.vapor[i]) <= 1.0e-14,
                "recording provider changed vapor composition");
    }
}

th::CpaPtRootSet replay_root_search(
    const RootRecord& record,
    const th::CpaParameterSet& parameters,
    ContinuationStats& stats) {
    const auto& options = record.roots.options;
    const auto composition = std::span<const double>(record.composition);
    th::CpaPtRootSet result;
    result.options = options;
    result.pressure_pa = record.pressure_pa;
    result.temperature_k = record.temperature_k;
    result.composition = record.composition;
    ++stats.root_searches;

    const double b = th::cpa_detail::cpa_b_mix(composition, parameters);
    const double pressure_tolerance = std::max(
        options.pressure_absolute_tolerance_pa,
        options.pressure_relative_tolerance * record.pressure_pa);
    constexpr double endpoint_guard =
        1024.0 * std::numeric_limits<double>::epsilon();
    const double u_max = 1.0 - endpoint_guard;
    std::map<double, std::vector<double>> cache;

    auto evaluate_reduced_density = [&](double u, double& residual) -> bool {
        if (result.evaluations >= options.max_evaluations) {
            result.status = th::CpaPtRootStatus::evaluation_limit;
            return false;
        }
        ++result.evaluations;
        const auto evaluated = audit_density_state(
            record.pressure_pa, record.temperature_k, u, b,
            composition, parameters, options.phase, cache, stats);
        residual = evaluated.pressure_residual_pa;
        return true;
    };

    std::vector<th::cpa_detail::CpaPtSample> samples;
    samples.reserve(options.scan_intervals + 1U);
    samples.push_back({0.0, -record.pressure_pa});
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

    std::vector<th::cpa_detail::CpaPtRootCandidate> candidates;
    const auto append_candidate = [&](double u, double residual, int slope_sign) {
        if (candidates.size() >= options.max_roots) {
            result.status = th::CpaPtRootStatus::root_limit;
            return false;
        }
        candidates.push_back({u, residual, slope_sign});
        return true;
    };

    for (std::size_t i = 0U; i + 1U < samples.size(); ++i) {
        const int left_sign = th::cpa_detail::cpa_pressure_sign(
            samples[i].residual_pa, pressure_tolerance);
        const int right_sign = th::cpa_detail::cpa_pressure_sign(
            samples[i + 1U].residual_pa, pressure_tolerance);
        if (left_sign == 0 || right_sign == 0 || left_sign == right_sign) { continue; }
        double left_u = samples[i].reduced_density;
        double right_u = samples[i + 1U].reduced_density;
        const double initial_left_r = samples[i].residual_pa;
        const double initial_right_r = samples[i + 1U].residual_pa;
        double best_u = std::abs(initial_left_r) < std::abs(initial_right_r)
            ? left_u : right_u;
        double best_r = std::abs(initial_left_r) < std::abs(initial_right_r)
            ? initial_left_r : initial_right_r;
        bool converged = false;
        for (int iteration = 0; iteration < options.max_bisection_iterations; ++iteration) {
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
            const int mid_sign = th::cpa_detail::cpa_pressure_sign(mid_r, pressure_tolerance);
            if (mid_sign == 0) {
                best_u = mid_u;
                best_r = mid_r;
                converged = true;
                break;
            }
            if (mid_sign == left_sign) { left_u = mid_u; }
            else { right_u = mid_u; }
        }
        if (!converged && std::abs(best_r) > pressure_tolerance) {
            result.status = th::CpaPtRootStatus::iteration_limit;
            return result;
        }
        const int slope_sign = left_sign < right_sign ? 1 : -1;
        if (!append_candidate(best_u, best_r, slope_sign)) { return result; }
    }

    for (std::size_t begin = 1U; begin < samples.size();) {
        if (th::cpa_detail::cpa_pressure_sign(
                samples[begin].residual_pa, pressure_tolerance) != 0) {
            ++begin;
            continue;
        }
        std::size_t end = begin;
        std::size_t best = begin;
        while (end + 1U < samples.size() &&
               th::cpa_detail::cpa_pressure_sign(
                   samples[end + 1U].residual_pa, pressure_tolerance) == 0) {
            ++end;
            if (std::abs(samples[end].residual_pa) <
                std::abs(samples[best].residual_pa)) {
                best = end;
            }
        }
        const int left_sign = begin > 0U
            ? th::cpa_detail::cpa_pressure_sign(
                  samples[begin - 1U].residual_pa, pressure_tolerance)
            : 0;
        const int right_sign = end + 1U < samples.size()
            ? th::cpa_detail::cpa_pressure_sign(
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
        result.status = th::CpaPtRootStatus::no_root;
        return result;
    }
    std::sort(candidates.begin(), candidates.end(),
              [](const auto& first, const auto& second) {
                  return first.reduced_density < second.reduced_density;
              });
    std::vector<th::cpa_detail::CpaPtRootCandidate> unique;
    unique.reserve(candidates.size());
    constexpr double root_merge_factor = 4096.0;
    for (const auto& candidate : candidates) {
        if (!unique.empty()) {
            const double scale = std::max(
                {1.0, std::abs(unique.back().reduced_density),
                 std::abs(candidate.reduced_density)});
            if (std::abs(candidate.reduced_density - unique.back().reduced_density) <=
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
            result.status = th::CpaPtRootStatus::evaluation_limit;
            result.roots.clear();
            return result;
        }
        ++result.evaluations;
        const auto evaluated = audit_density_state(
            record.pressure_pa, record.temperature_k,
            candidate.reduced_density, b, composition,
            parameters, options.phase, cache, stats);
        th::CpaPtRoot cold_root;
        cold_root.molar_density_mol_per_m3 =
            candidate.reduced_density / b;
        cold_root.pressure_residual_pa =
            evaluated.cold.pressure_pa - record.pressure_pa;
        cold_root.pressure_slope_sign = candidate.slope_sign;
        th::cpa_detail::cpa_fill_ln_phi(
            record.pressure_pa, record.temperature_k, composition,
            parameters, evaluated.cold, cold_root);

        th::CpaPtRoot warm_root;
        warm_root.molar_density_mol_per_m3 = cold_root.molar_density_mol_per_m3;
        warm_root.pressure_residual_pa =
            evaluated.effective.pressure_pa - record.pressure_pa;
        warm_root.pressure_slope_sign = candidate.slope_sign;
        th::cpa_detail::cpa_fill_ln_phi(
            record.pressure_pa, record.temperature_k, composition,
            parameters, evaluated.effective, warm_root);
        double ln_phi_diff = 0.0;
        for (std::size_t i = 0U; i < cold_root.ln_phi.size(); ++i) {
            ln_phi_diff = std::max(
                ln_phi_diff,
                std::abs(cold_root.ln_phi[i] - warm_root.ln_phi[i]));
        }
        stats.max_raw_root_ln_phi_diff =
            std::max(stats.max_raw_root_ln_phi_diff, ln_phi_diff);
        if (evaluated.accepted_continuation &&
            ln_phi_diff > parity_gate::phase_max_abs_ln_phi) {
            ++stats.root_lnphi_fallbacks;
            stats.effective_sweeps_with_fallback += evaluated.cold_sweeps;
            warm_root = cold_root;
        }
        result.roots.push_back(std::move(cold_root));
    }

    result.status = near_multiple
        ? th::CpaPtRootStatus::near_multiple
        : th::CpaPtRootStatus::success;
    return result;
}

void compare_replay(const th::CpaPtRootSet& expected, const th::CpaPtRootSet& actual) {
    require(expected.status == actual.status, "root replay changed root-search status");
    require(expected.evaluations == actual.evaluations,
            "root replay changed production density-evaluation count");
    require(expected.roots.size() == actual.roots.size(),
            "root replay changed root count");
    for (std::size_t k = 0U; k < expected.roots.size(); ++k) {
        const auto& first = expected.roots[k];
        const auto& second = actual.roots[k];
        const double density_scale = std::max(1.0, std::abs(first.molar_density_mol_per_m3));
        require(std::abs(first.molar_density_mol_per_m3 -
                         second.molar_density_mol_per_m3) <=
                    replay_relative_density_guard * density_scale,
                "root replay changed root density");
        require(first.pressure_slope_sign == second.pressure_slope_sign,
                "root replay changed root slope sign");
        require(first.ln_phi.size() == second.ln_phi.size(),
                "root replay changed ln(phi) dimension");
        for (std::size_t i = 0U; i < first.ln_phi.size(); ++i) {
            require(std::abs(first.ln_phi[i] - second.ln_phi[i]) <=
                        replay_abs_ln_phi_guard,
                    "root replay changed cold ln(phi)");
        }
    }
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
        std::vector<RootRecord> records;

        std::cout << std::setprecision(17);
        for (const auto& state : cpa_physical_test::points()) {
            const auto feed = cpa_physical_test::feed(state, false);
            const auto starts = cpa_physical_test::starts(state, false);
            const std::vector<std::vector<double>> final_starts{feed};
            const auto options = audit_split_options();

            fl::CpaVleEvaluator production_evaluator(model);
            const auto production = fl::solve_cpa_pt_vle(
                state.pressure_pa, cpa_physical_test::temperature_k,
                feed, production_evaluator, options, starts, final_starts);

            const std::size_t before = records.size();
            RecordingEvaluator recording_evaluator(model, phase_options, records);
            const auto recording = fl::solve_pt_vle(
                state.pressure_pa, cpa_physical_test::temperature_k,
                feed, recording_evaluator, recording_evaluator,
                options, starts, final_starts);
            compare_outer_solution(production.solution, recording);

            const std::size_t expected_calls =
                recording.initial_stability.evaluations +
                recording.split_evaluations +
                (recording.final_stability
                    ? recording.final_stability->evaluations
                    : 0U);
            require(records.size() - before == expected_calls,
                    "recorded root-search count does not match solver counters");
            std::cout
                << "CPA_ASSOC_CONTINUATION_CAPTURE"
                << " PPa=" << state.pressure_pa
                << " root_calls=" << expected_calls
                << " betaV=" << recording.candidate()->fractions.vapor_fraction
                << '\n';
        }

        ContinuationStats stats;
        std::size_t expected_density_evaluations = 0U;
        for (const auto& record : records) {
            expected_density_evaluations += record.roots.evaluations;
            const auto replay = replay_root_search(record, parameters, stats);
            compare_replay(record.roots, replay);
        }
        require(stats.root_searches == records.size(),
                "not every recorded root search was replayed");
        require(stats.density_evaluations == expected_density_evaluations,
                "shadow replay did not cover every production density evaluation");
        require(stats.unseeded_evaluations == stats.root_searches,
                "each root search should have exactly one unseeded first density evaluation");
        require(stats.seeded_evaluations + stats.unseeded_evaluations ==
                    stats.density_evaluations,
                "seeded/unseeded accounting mismatch");
        require(stats.max_raw_root_ln_phi_diff <= parity_gate::phase_max_abs_ln_phi ||
                    stats.root_lnphi_fallbacks > 0U,
                "root ln(phi) guard exceeded without fallback accounting");

        const double seeded = static_cast<double>(stats.seeded_evaluations);
        const double fallback_rate = seeded > 0.0
            ? static_cast<double>(stats.fallback_evaluations + stats.root_lnphi_fallbacks) / seeded
            : 0.0;
        const double sweep_reduction = stats.cold_sweeps > 0U
            ? 1.0 - static_cast<double>(stats.effective_sweeps_with_fallback) /
                        static_cast<double>(stats.cold_sweeps)
            : 0.0;
        const double mean_cold_sweeps = stats.density_evaluations > 0U
            ? static_cast<double>(stats.cold_sweeps) /
                static_cast<double>(stats.density_evaluations)
            : 0.0;
        const double mean_seeded_warm_sweeps = stats.seeded_evaluations > 0U
            ? static_cast<double>(stats.seeded_warm_attempt_sweeps) /
                static_cast<double>(stats.seeded_evaluations)
            : 0.0;

        std::cout
            << "CPA_ASSOC_CONTINUATION_SUMMARY"
            << " root_searches=" << stats.root_searches
            << " density_evaluations=" << stats.density_evaluations
            << " seeded_evaluations=" << stats.seeded_evaluations
            << " exact_seed_hits=" << stats.exact_seed_hits
            << " accepted_continuations=" << stats.accepted_continuations
            << " fallback_evaluations=" << stats.fallback_evaluations
            << " nonconverged_fallbacks=" << stats.nonconverged_fallbacks
            << " site_guard_fallbacks=" << stats.site_guard_fallbacks
            << " pressure_guard_fallbacks=" << stats.pressure_guard_fallbacks
            << " root_lnphi_fallbacks=" << stats.root_lnphi_fallbacks
            << " fallback_rate=" << fallback_rate
            << " cold_sweeps=" << stats.cold_sweeps
            << " effective_sweeps_with_fallback=" << stats.effective_sweeps_with_fallback
            << " sweep_reduction_fraction=" << sweep_reduction
            << " mean_cold_sweeps=" << mean_cold_sweeps
            << " mean_seeded_warm_sweeps=" << mean_seeded_warm_sweeps
            << " warm_better=" << stats.warm_better
            << " warm_equal=" << stats.warm_equal
            << " warm_worse=" << stats.warm_worse
            << " max_cold_sweeps=" << stats.max_cold_sweeps
            << " max_warm_sweeps=" << stats.max_warm_sweeps
            << " max_raw_dX=" << stats.max_raw_site_fraction_diff
            << " max_raw_dP_pa=" << stats.max_raw_pressure_diff_pa
            << " max_raw_dlnphi=" << stats.max_raw_root_ln_phi_diff
            << " site_guard=" << site_fraction_audit_guard
            << " pressure_guard_pa=" << parity_gate::phase_max_abs_pressure_total_pa
            << " lnphi_guard=" << parity_gate::phase_max_abs_ln_phi
            << " production_solver_changes=none"
            << " optimization=none_audit_only"
            << '\n';
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "[FAIL] " << error.what() << '\n';
        return 1;
    }
}
