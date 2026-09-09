#ifndef MPMC_FLASH_PT_SPLIT_HPP
#define MPMC_FLASH_PT_SPLIT_HPP

#include <mpmc/flash/rachford_rice.hpp>
#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <numeric>
#include <optional>
#include <span>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace mpmc::flash {

// Requested candidate branches, not a universal physical phase-identification rule.
enum class PtPhaseRole { liquid_candidate, vapor_candidate };
struct PtSplitPhase {
    StabilityPhase activity;
    double z{}; // p*v/(R*T), positive; used only for candidate density ordering.
};
enum class PtSplitAttemptStatus {
    converged, no_interior_rr_root, rr_failure, indistinguishable_phases,
    phase_disappearance, property_failure, iteration_limit, evaluation_limit,
    line_search_failed, unrepresentable_seed, balance_failure
};
struct PtSplitIterationOptions {
    double fugacity_tolerance{1e-11}; // Absolute max |ln(fL/fV)|, no mole-fraction weight.
    double mass_absolute_tolerance{1e-12};
    double mass_relative_tolerance{1e-10}; // Per POSITIVE feed component, with no floor.
    double minimum_phase_fraction{1e-10}; // Below this: unresolved disappearance, not clipping.
    double log_k_separation{1e-7};
    double relative_z_separation{1e-8};
    double max_log_step{2};
    double residual_decrease{1e-4}; // Unresolved-Gibbs step progress, not a stopping tolerance.
    int max_iterations{512};
    int max_backtracks{24};
    std::size_t max_evaluations{20000}; // Single-phase property calls, including failures.
    RachfordRiceOptions rr;
    double gibbs_armijo{1e-4};
};
struct PtSplitOptions {
    PtSplitIterationOptions iteration;
    StabilityOptions initial_stability;
    StabilityOptions final_stability;
    std::size_t max_split_attempts{16}; // Both role assignments of negative TPD witnesses.
};
struct PtSplitState {
    RachfordRiceResult fractions;
    PtSplitPhase liquid, vapor;
    std::vector<double> log_k, fugacity_residual, common_log_activity;
    double fugacity_norm{};
    double reduced_gibbs{}; // sum phase fractions*x*ln(x*phi); common reference terms omitted.
    double gibbs_roundoff_guard{};
};
struct PtSplitAttempt {
    PtSplitAttemptStatus status{PtSplitAttemptStatus::unrepresentable_seed};
    std::vector<double> initial_log_k;
    std::optional<PtSplitState> point; // Last accepted iterate; may NOT be converged.
    std::optional<std::size_t> witness_trial;
    bool witness_as_vapor{};
    int iterations{};
    std::size_t evaluations{}, backtracks{}, rejected_evaluations{};
    std::optional<StabilityPropertyIssue> property_issue;
    std::string diagnostic;
    std::size_t gibbs_descent_steps{}, residual_increase_steps{};
};
enum class PtSplitStatus {
    single_phase_no_instability_found, two_phase_no_instability_found,
    phase_set_unstable, indeterminate
};
struct PtSplitResult {
    PtSplitStatus status{PtSplitStatus::indeterminate};
    static constexpr bool global_stability_proven = false;
    static constexpr const char* convention = "PT-VLE/logK-SSI/RR/common-tangent-v1";
    PtSplitOptions options;
    StabilityResult initial_stability;
    std::vector<PtSplitAttempt> attempts;
    std::optional<std::size_t> selected_attempt;
    std::optional<StabilityResult> final_stability;
    std::size_t split_evaluations{};
    bool attempt_limit_reached{false};
    double gibbs_change{std::numeric_limits<double>::quiet_NaN()};
    // Half of the two computed phase log-activity discrepancy. This is a numerical
    // reference disagreement allowance, not a rigorous error bound to exact equilibrium.
    double common_reference_allowance{};
    std::string diagnostic;
    [[nodiscard]] bool equations_converged() const noexcept {
        return selected_attempt && *selected_attempt < attempts.size() &&
            attempts[*selected_attempt].status == PtSplitAttemptStatus::converged &&
            attempts[*selected_attempt].point.has_value();
    }
    [[nodiscard]] const PtSplitState* candidate() const & noexcept {
        return equations_converged() ? &*attempts[*selected_attempt].point : nullptr;
    }
    const PtSplitState* candidate() const && = delete;
};

namespace detail {
inline void split_check_options(const PtSplitIterationOptions& o) {
    const auto positive = [](double a) { return std::isfinite(a) && a > 0; };
    if (!positive(o.fugacity_tolerance) || !positive(o.mass_absolute_tolerance) ||
        !positive(o.mass_relative_tolerance) || !positive(o.minimum_phase_fraction) ||
        o.minimum_phase_fraction >= 0.5 || !positive(o.log_k_separation) ||
        !positive(o.relative_z_separation) || o.relative_z_separation >= 1 ||
        !positive(o.max_log_step) || !positive(o.residual_decrease) || o.residual_decrease >= 1 ||
        !positive(o.gibbs_armijo) || o.gibbs_armijo >= 1 ||
        o.max_iterations < 0 || o.max_backtracks <= 0) {
        throw std::invalid_argument("PT split: invalid tolerance or iteration option");
    }
    rr_check_options(o.rr);
}
inline void split_check_pt(double p, double t) {
    if (!std::isfinite(p) || p <= 0 || !std::isfinite(t) || t <= 0) {
        throw std::domain_error("PT split: finite p>0 Pa and T>0 K required");
    }
}
class SplitStepError : public std::runtime_error {
public:
    SplitStepError(PtSplitAttemptStatus status, const char* message)
        : std::runtime_error(message), status_(status) {}
    [[nodiscard]] PtSplitAttemptStatus status() const noexcept { return status_; }
private:
    PtSplitAttemptStatus status_;
};
inline void split_check_phase(const PtSplitPhase& phase, std::size_t n) {
    stability_check_phase(phase.activity, n);
    if (!std::isfinite(phase.z) || phase.z <= 0 || !phase.activity.smooth) {
        throw StabilityPropertyError(StabilityPropertyIssue::ill_conditioned_root,
                                     "PT split: invalid or nonsmooth candidate branch");
    }
}
inline bool split_balance_ok(const PtSplitState& state, const PtSplitIterationOptions& o) {
    return state.fractions.mass_absolute <= o.mass_absolute_tolerance &&
           state.fractions.mass_relative <= o.mass_relative_tolerance;
}
// For material-balanced RR states, c_i=x_i*y_i/z_i and H=sum((y_i-x_i)^2/z_i).
// RR gives d(beta)=sum(c_i*d(logK_i))/H. Writing vapor amounts v_i=beta*y_i,
// d(v_i)=c_i*[d(beta)+beta*(1-beta)*d(logK_i)]. Gibbs-Duhem then yields
// dG=-sum(r_i*dv_i). Along d(logK)=r/scale this is a nonpositive quadratic form.
// No EOS derivatives or numerical differentiation are needed for this slope.
inline double split_gibbs_slope(const PtSplitState& point, std::span<const double> feed, double scale) {
    double h = 0, hc = 0, first = 0, fc = 0, second = 0, sc = 0;
    for (std::size_t i = 0; i < feed.size(); ++i) {
        if (feed[i] == 0) { continue; }
        const double x = point.fractions.liquid[i], y = point.fractions.vapor[i];
        const double c = (x/feed[i])*y;
        const double difference = y-x;
        const double r = point.fugacity_residual[i];
        stability_add(difference*(difference/feed[i]), h, hc);
        stability_add(c*r, first, fc);
        stability_add((c*r)*r, second, sc);
    }
    if (!(h > 0) || !std::isfinite(h) || !std::isfinite(first) ||
        !std::isfinite(second) || !std::isfinite(scale) || !(scale >= 1)) {
        return std::numeric_limits<double>::quiet_NaN();
    }
    const double beta = point.fractions.vapor_fraction;
    return -(beta*(1-beta)*second + (first/h)*first)/scale;
}
inline bool split_accept_step(const PtSplitState& point, const PtSplitState& next,
                               double predicted, double relative_step, const PtSplitIterationOptions& o) {
    const double guard = point.gibbs_roundoff_guard + next.gibbs_roundoff_guard;
    if (predicted > guard) {
        return next.reduced_gibbs <= point.reduced_gibbs-o.gibbs_armijo*predicted;
    }
    // An unresolved objective reduction cannot be accepted merely because the
    // rounded Gibbs values compare equal. Require actual residual progress.
    return next.reduced_gibbs <= point.reduced_gibbs+guard &&
        (next.fugacity_norm <= o.fugacity_tolerance ||
         point.fugacity_norm-next.fugacity_norm >
             o.residual_decrease*relative_step*point.fugacity_norm);
}
inline std::optional<std::vector<double>> split_seed(
    std::span<const double> z, std::span<const double> w, bool as_vapor) {
    double amount = 0.1;
    for (std::size_t i = 0; i < z.size(); ++i) {
        if (z[i] > 0) { amount = std::min(amount, 0.5*(z[i]/w[i])); }
    }
    if (!(amount > 0 && amount < 1)) { return std::nullopt; }
    std::vector<double> other(z.size()), log_k(z.size());
    for (std::size_t i = 0; i < z.size(); ++i) {
        if (z[i] == 0) { continue; }
        other[i] = std::fma(-amount, w[i], z[i])/(1-amount);
        if (!(other[i] > 0)) { return std::nullopt; }
    }
    const double sum = stability_sum(other);
    if (!std::isfinite(sum) || std::abs(sum-1) > 64*stability_eps) { return std::nullopt; }
    for (std::size_t i = 0; i < z.size(); ++i) {
        if (z[i] == 0) { continue; }
        other[i] /= sum;
        log_k[i] = (as_vapor ? 1.0 : -1.0)*(std::log(w[i])-std::log(other[i]));
        if (!std::isfinite(log_k[i])) { return std::nullopt; }
    }
    return log_k;
}
} // namespace detail

// Low-level fixed-seed iteration, not a full flash/stability decision. Provider is
// provider(p_Pa,T_K,normalized_composition,PtPhaseRole)->PtSplitPhase. It must use
// common component reference states and mechanically admissible candidate branches.
// Input logK is not assumed physically valid and endpoint RR results are NOT single-phase evidence.
template <typename Provider>
[[nodiscard]] PtSplitAttempt iterate_pt_split(
    double pressure_pa, double temperature_k, std::span<const double> feed,
    std::span<const double> initial_log_k, Provider&& provider, PtSplitIterationOptions options = {}) {
    detail::split_check_options(options);
    detail::split_check_pt(pressure_pa, temperature_k);
    if (feed.empty() || feed.size() > options.rr.max_components) {
        throw std::length_error("PT split: component quota exceeded or empty feed");
    }
    if (feed.size() != initial_log_k.size()) { throw std::invalid_argument("PT split: logK dimension mismatch"); }
    (void)detail::stability_check_composition(feed);
    for (double value : initial_log_k) {
        if (!std::isfinite(value)) { throw std::domain_error("PT split: finite initial log K required"); }
    }
    const std::size_t n = feed.size();
    PtSplitAttempt result;
    result.initial_log_k.assign(initial_log_k.begin(), initial_log_k.end());
    const auto evaluate_phase = [&](std::span<const double> w, PtPhaseRole role) {
        if (result.evaluations >= options.max_evaluations) {
            throw detail::SplitStepError(PtSplitAttemptStatus::evaluation_limit, "PT split: property budget exhausted");
        }
        ++result.evaluations;
        auto phase = provider(pressure_pa, temperature_k, w, role);
        detail::split_check_phase(phase, n);
        return phase;
    };
    const auto evaluate = [&](std::span<const double> log_k) {
        PtSplitState state;
        state.fractions = solve_rachford_rice(feed, log_k, options.rr);
        if (state.fractions.status != RachfordRiceStatus::interior) {
            const auto reason = state.fractions.status == RachfordRiceStatus::degenerate
                ? PtSplitAttemptStatus::indistinguishable_phases
                : (state.fractions.status == RachfordRiceStatus::no_resolved_interior_root
                    ? PtSplitAttemptStatus::no_interior_rr_root : PtSplitAttemptStatus::rr_failure);
            throw detail::SplitStepError(reason, "PT split: RR has no usable interior solution");
        }
        state.liquid = evaluate_phase(state.fractions.liquid, PtPhaseRole::liquid_candidate);
        state.vapor = evaluate_phase(state.fractions.vapor, PtPhaseRole::vapor_candidate);
        state.log_k.resize(n);
        state.fugacity_residual.resize(n);
        state.common_log_activity.resize(n);
        double correction = 0, magnitude = 1, mc = 0;
        const double beta = state.fractions.vapor_fraction;
        for (std::size_t i = 0; i < n; ++i) {
            if (feed[i] == 0) { continue; }
            const double lx = std::log(state.fractions.liquid[i]);
            const double ly = std::log(state.fractions.vapor[i]);
            const double mu_l = lx + state.liquid.activity.ln_phi[i];
            const double mu_v = ly + state.vapor.activity.ln_phi[i];
            state.log_k[i] = ly-lx;
            state.fugacity_residual[i] = (lx-ly) +
                (state.liquid.activity.ln_phi[i]-state.vapor.activity.ln_phi[i]);
            state.common_log_activity[i] = std::midpoint(mu_l, mu_v);
            state.fugacity_norm = std::max(state.fugacity_norm, std::abs(state.fugacity_residual[i]));
            detail::stability_add((1-beta)*state.fractions.liquid[i]*mu_l +
                beta*state.fractions.vapor[i]*mu_v, state.reduced_gibbs, correction);
            detail::stability_add((1-beta)*state.fractions.liquid[i]*(std::abs(lx)+std::abs(state.liquid.activity.ln_phi[i])) +
                beta*state.fractions.vapor[i]*(std::abs(ly)+std::abs(state.vapor.activity.ln_phi[i])), magnitude, mc);
            if (!std::isfinite(state.fugacity_residual[i]) || !std::isfinite(state.common_log_activity[i])) {
                throw StabilityPropertyError(StabilityPropertyIssue::nonfinite_properties, "PT split: nonfinite activity arithmetic");
            }
        }
        state.gibbs_roundoff_guard = 256*detail::stability_eps*magnitude;
        if (!std::isfinite(state.reduced_gibbs) || !std::isfinite(state.gibbs_roundoff_guard)) {
            throw StabilityPropertyError(StabilityPropertyIssue::nonfinite_properties, "PT split: nonfinite Gibbs arithmetic");
        }
        return state;
    };
    try { result.point = evaluate(initial_log_k); }
    catch (const detail::SplitStepError& error) {
        result.status = error.status(); result.diagnostic = error.what(); return result;
    } catch (const StabilityPropertyError& error) {
        result.status = PtSplitAttemptStatus::property_failure;
        result.property_issue = error.issue(); result.diagnostic = error.what(); return result;
    }
    for (;;) {
        const auto& point = *result.point;
        if (point.fugacity_norm <= options.fugacity_tolerance) {
            const double beta = point.fractions.vapor_fraction;
            double contrast = 0;
            for (double value : point.log_k) { contrast = std::max(contrast, std::abs(value)); }
            if (!detail::split_balance_ok(point, options)) { result.status = PtSplitAttemptStatus::balance_failure; }
            else if (beta <= options.minimum_phase_fraction || 1-beta <= options.minimum_phase_fraction) {
                result.status = PtSplitAttemptStatus::phase_disappearance;
            } else if (contrast <= options.log_k_separation ||
                point.vapor.z-point.liquid.z <= options.relative_z_separation*std::max(point.vapor.z, point.liquid.z)) {
                result.status = PtSplitAttemptStatus::indistinguishable_phases;
            } else { result.status = PtSplitAttemptStatus::converged; }
            return result;
        }
        if (result.iterations >= options.max_iterations) {
            result.status = PtSplitAttemptStatus::iteration_limit; return result;
        }
        const double scale = std::max(1.0, point.fugacity_norm/options.max_log_step);
        const double slope = detail::split_gibbs_slope(point, feed, scale);
        if (!std::isfinite(slope) || !(slope < 0)) {
            result.status = PtSplitAttemptStatus::line_search_failed;
            result.diagnostic = "PT split: degenerate or unrepresentable Gibbs descent slope";
            return result;
        }
        bool accepted = false;
        double alpha = 1;
        result.status = PtSplitAttemptStatus::line_search_failed;
        for (int bt = 0; bt < options.max_backtracks; ++bt, alpha *= 0.5) {
            if (bt > 0) { ++result.backtracks; }
            std::vector<double> next_log_k(n);
            for (std::size_t i = 0; i < n; ++i) {
                next_log_k[i] = point.log_k[i] + (alpha/scale)*point.fugacity_residual[i];
            }
            try {
                PtSplitState next = evaluate(next_log_k);
                const double predicted = -alpha*slope;
                if (detail::split_accept_step(point, next, predicted, alpha/scale, options)) {
                    if (predicted > point.gibbs_roundoff_guard+next.gibbs_roundoff_guard) {
                        ++result.gibbs_descent_steps;
                    }
                    if (next.fugacity_norm > point.fugacity_norm) { ++result.residual_increase_steps; }
                    result.point = std::move(next); accepted = true; break;
                }
                ++result.rejected_evaluations;
            } catch (const detail::SplitStepError& error) {
                ++result.rejected_evaluations;
                result.diagnostic = error.what();
                if (error.status() == PtSplitAttemptStatus::evaluation_limit) {
                    result.status = error.status(); return result;
                }
            } catch (const StabilityPropertyError& error) {
                ++result.rejected_evaluations;
                result.property_issue = error.issue(); result.diagnostic = error.what();
            }
        }
        if (!accepted) { return result; }
        ++result.iterations;
    }
}

// Full baseline: feed stability -> witnessed, material-balanced seeds -> phase
// split -> common-tangent phase-set review. This remains a finite local search.
// Both providers must describe the SAME model/reference states. The stability
// provider must include every allowed trial branch, not just the two requested roles.
template <typename StabilityProvider, typename PhaseProvider>
[[nodiscard]] PtSplitResult solve_pt_vle(
    double pressure_pa, double temperature_k, std::span<const double> feed,
    StabilityProvider&& stability_provider, PhaseProvider&& phase_provider,
    PtSplitOptions options = {}, std::span<const std::vector<double>> initial_starts = {},
    std::span<const std::vector<double>> final_starts = {}) {
    detail::split_check_options(options.iteration);
    detail::stability_check_options(options.initial_stability);
    detail::stability_check_options(options.final_stability);
    detail::split_check_pt(pressure_pa, temperature_k);
    const std::size_t n = feed.size();
    if (n == 0 || n > options.iteration.rr.max_components || n > options.final_stability.max_components) {
        throw std::length_error("PT split: component quota exceeded or empty feed");
    }
    (void)detail::stability_check_composition(feed);
    // Validate final-search resources and starts BEFORE initial property calls.
    std::size_t active = 0;
    for (double value : feed) { if (value > 0) { ++active; } }
    if (active > std::numeric_limits<std::size_t>::max()-4) { throw std::length_error("PT split: start count overflow"); }
    const std::size_t generated = options.final_stability.automatic_starts ? (active == 1 ? 1 : active+2) : 0;
    const std::size_t required = generated+2; // Always review both computed phases as extra starts.
    if (required > options.final_stability.max_starts ||
        final_starts.size() > options.final_stability.max_starts-required ||
        required+final_starts.size() > options.final_stability.max_start_entries/n) {
        throw std::length_error("PT split: final stability start quota exceeded");
    }
    for (const auto& w : final_starts) {
        if (w.size() != n) { throw std::invalid_argument("PT split: final start dimension mismatch"); }
        (void)detail::stability_check_composition(w);
        for (std::size_t i = 0; i < n; ++i) {
            if ((feed[i] > 0) != (w[i] > 0)) { throw std::domain_error("PT split: final start support mismatch"); }
        }
    }
    PtSplitResult result;
    result.options = options;
    result.initial_stability = test_pt_stability(pressure_pa, temperature_k, feed,
        stability_provider, options.initial_stability, initial_starts);
    if (result.initial_stability.status == StabilityStatus::no_instability_found) {
        result.status = PtSplitStatus::single_phase_no_instability_found; return result;
    }
    if (result.initial_stability.status != StabilityStatus::unstable) {
        result.diagnostic = "initial feed stability is indeterminate"; return result;
    }
    const auto& z = result.initial_stability.feed;
    for (std::size_t k = 0; k < result.initial_stability.trials.size(); ++k) {
        const auto& witness = result.initial_stability.trials[k];
        if (witness.status != StabilityTrialStatus::negative_tpd || !witness.point ||
            !detail::stability_negative(*witness.point, options.initial_stability)) { continue; }
        for (bool as_vapor : {true, false}) {
            if (result.attempts.size() >= options.max_split_attempts ||
                result.split_evaluations >= options.iteration.max_evaluations) {
                result.attempt_limit_reached = true; break;
            }
            auto seed = detail::split_seed(z, witness.point->composition, as_vapor);
            PtSplitAttempt attempt;
            if (seed) {
                auto iteration = options.iteration;
                iteration.max_evaluations -= result.split_evaluations;
                attempt = iterate_pt_split(pressure_pa, temperature_k, z, *seed, phase_provider, iteration);
            } else { attempt.diagnostic = "finite, positive material-balanced witness seed is not representable"; }
            attempt.witness_trial = k;
            attempt.witness_as_vapor = as_vapor;
            result.split_evaluations += attempt.evaluations;
            const std::size_t index = result.attempts.size();
            result.attempts.push_back(std::move(attempt));
            if (result.attempts.back().status == PtSplitAttemptStatus::converged &&
                (!result.selected_attempt || result.attempts.back().point->reduced_gibbs < result.candidate()->reduced_gibbs)) {
                result.selected_attempt = index;
            }
        }
        if (result.attempt_limit_reached) { break; }
    }
    if (!result.selected_attempt) {
        result.diagnostic = "no witnessed seed produced a distinct, balanced, fugacity-converged pair"; return result;
    }
    const auto& candidate = *result.candidate();
    double feed_g = 0, correction = 0, magnitude = 1;
    for (std::size_t i = 0; i < n; ++i) {
        if (z[i] == 0) { continue; }
        const double lz = std::log(z[i]);
        const double phi = result.initial_stability.reference->ln_phi[i];
        detail::stability_add(z[i]*(lz+phi), feed_g, correction);
        magnitude += z[i]*(std::abs(lz)+std::abs(phi));
    }
    result.gibbs_change = candidate.reduced_gibbs-feed_g;
    const double gibbs_guard = candidate.gibbs_roundoff_guard+256*detail::stability_eps*magnitude;
    result.common_reference_allowance = 0.5*candidate.fugacity_norm;
    auto final_options = options.final_stability;
    final_options.tpd_tolerance += result.common_reference_allowance;
    if (!std::isfinite(result.gibbs_change) || !std::isfinite(gibbs_guard) ||
        !std::isfinite(final_options.tpd_tolerance)) {
        result.diagnostic = "nonrepresentable phase-set review scale"; return result;
    }
    std::vector<std::vector<double>> starts(final_starts.begin(), final_starts.end());
    starts.push_back(candidate.fractions.liquid);
    starts.push_back(candidate.fractions.vapor);
    result.final_stability = test_pt_stability_against(pressure_pa, temperature_k, z,
        candidate.common_log_activity, stability_provider, final_options, starts);
    if (result.final_stability->status == StabilityStatus::unstable) {
        result.status = PtSplitStatus::phase_set_unstable;
        result.diagnostic = "common-tangent review found a lower-Gibbs trial; not an accepted two-phase state";
    } else if (result.gibbs_change > gibbs_guard) {
        result.diagnostic = "converged pair has higher Gibbs energy than the original feed";
    } else if (result.final_stability->status == StabilityStatus::no_instability_found) {
        result.status = PtSplitStatus::two_phase_no_instability_found;
    } else { result.diagnostic = "split equations converged; final phase-set stability is indeterminate"; }
    return result;
}

} // namespace mpmc::flash
#endif // MPMC_FLASH_PT_SPLIT_HPP
