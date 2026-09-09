#ifndef MPMC_FLASH_PT_STABILITY_HPP
#define MPMC_FLASH_PT_STABILITY_HPP

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <optional>
#include <span>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace mpmc::flash {

// Finite local searches NEVER certify global stability or a phase count.
enum class StabilityStatus { unstable, no_instability_found, indeterminate };
enum class StabilityTrialStatus {
    negative_tpd, stationary, iteration_limit, evaluation_limit,
    line_search_failed, property_failure, nonsmooth, unrepresentable_composition
};
enum class StabilityPropertyIssue {
    root_topology, root_iteration_limit, root_range, ill_conditioned_root,
    no_admissible_branch, nonfinite_properties
};
class StabilityPropertyError : public std::runtime_error {
public:
    StabilityPropertyError(StabilityPropertyIssue issue, const std::string& message)
        : std::runtime_error(message), issue_(issue) {}
    [[nodiscard]] StabilityPropertyIssue issue() const noexcept { return issue_; }
private:
    StabilityPropertyIssue issue_;
};

// Provider must select the least Gibbs-energy admissible fluid branch at THIS
// composition, with common chemical-potential reference states and Gibbs-Duhem
// consistency. branch is a diagnostic, not a permanent phase identity.
struct StabilityPhase {
    std::vector<double> ln_phi;
    std::size_t branch{};
    bool smooth{true};
};
struct StabilityOptions {
    double tpd_tolerance{1e-10};         // Absolute, dimensionless D/(RT).
    double stationarity_tolerance{1e-8}; // Max unweighted chemical-potential residual/(RT).
    double armijo{1e-4};
    double log_step_limit{4.0};
    int max_iterations{512};
    int max_backtracks{32};
    std::size_t max_evaluations{100000}; // Includes reference and rejected evaluations.
    std::size_t max_components{256};
    std::size_t max_starts{1024};
    std::size_t max_start_entries{262144};
    bool automatic_starts{true}; // Feed, uniform, 0.1*feed+0.9*each active vertex.
};
struct TpdPoint {
    std::vector<double> composition;
    std::vector<double> residual; // q_i - weighted mean(q); 0 on absent feed components.
    double value{};
    double roundoff_guard{}; // Heuristic arithmetic guard, NOT an interval/error enclosure.
    double stationarity{};
    std::size_t branch{};
    bool smooth{true};
};
struct StabilityTrial {
    std::vector<double> initial_composition;
    std::optional<TpdPoint> point;
    StabilityTrialStatus status{StabilityTrialStatus::evaluation_limit};
    std::optional<StabilityPropertyIssue> property_issue;
    std::string diagnostic;
    int iterations{};
    std::size_t evaluations{};
    std::size_t backtracks{};
    std::size_t rejected_property_evaluations{};
};
struct StabilityResult {
    StabilityStatus status{StabilityStatus::indeterminate};
    static constexpr bool global_stability_proven = false;
    static constexpr const char* search_convention = "normalized-TPD/log-descent/multistart-v1";
    StabilityOptions options;
    double pressure_pa{}, temperature_k{};
    double input_feed_sum{};
    std::vector<double> feed; // Explicit roundoff-only normalization, not a cutoff or clipping.
    std::optional<StabilityPhase> reference;
    std::optional<StabilityPropertyIssue> reference_issue;
    std::string diagnostic;
    std::vector<StabilityTrial> trials;
    std::optional<TpdPoint> lowest_sampled; // NOT a global or necessarily stationary minimum.
    std::size_t evaluations{};
};

namespace detail {
inline constexpr double stability_eps = std::numeric_limits<double>::epsilon();
inline void stability_add(double term, double& sum, double& correction) {
    const double increment = term - correction;
    const double next = sum + increment;
    correction = (next - sum) - increment;
    sum = next;
}
inline double stability_sum(std::span<const double> x) {
    double sum = 0, correction = 0;
    for (double value : x) { stability_add(value, sum, correction); }
    return sum;
}
inline double stability_check_composition(std::span<const double> x) {
    if (x.empty()) { throw std::invalid_argument("stability: empty composition"); }
    for (double value : x) {
        if (!std::isfinite(value) || value < 0 || value > 1) {
            throw std::domain_error("stability: finite mole fractions in [0,1] required");
        }
    }
    const double sum = stability_sum(x);
    if (std::abs(sum - 1.0) > 64.0 * stability_eps) {
        throw std::domain_error("stability: composition is not normalized");
    }
    return sum;
}
inline void stability_check_options(const StabilityOptions& o) {
    if (!std::isfinite(o.tpd_tolerance) || o.tpd_tolerance < 0 ||
        !std::isfinite(o.stationarity_tolerance) || o.stationarity_tolerance <= 0 ||
        !std::isfinite(o.armijo) || o.armijo <= 0 || o.armijo >= 1 ||
        !std::isfinite(o.log_step_limit) || o.log_step_limit <= 0 ||
        o.max_iterations < 0 || o.max_backtracks <= 0 ||
        o.max_evaluations == 0 || o.max_components == 0 || o.max_starts == 0 ||
        o.max_start_entries == 0) {
        throw std::invalid_argument("stability: invalid tolerance, step or resource limit");
    }
}
inline void stability_check_phase(const StabilityPhase& phase, std::size_t n) {
    if (phase.ln_phi.size() != n) {
        throw std::logic_error("stability: provider output dimension mismatch");
    }
    for (double value : phase.ln_phi) {
        if (!std::isfinite(value)) {
            throw StabilityPropertyError(StabilityPropertyIssue::nonfinite_properties,
                                         "stability: provider returned nonfinite ln(phi)");
        }
    }
}
inline std::vector<double> stability_normalize(std::span<const double> x, double sum) {
    std::vector<double> result(x.begin(), x.end());
    for (double& value : result) { value /= sum; }
    return result;
}
inline bool stability_negative(const TpdPoint& point, const StabilityOptions& o) {
    return point.value < -o.tpd_tolerance - point.roundoff_guard;
}
} // namespace detail

// Normalized tangent-plane distance. Boundary w_i=0 uses lim w_i*log(w_i)=0;
// stationarity there is infinite on active feed support. w_i>0 with z_i=0 is
// infeasible for a closed, nonreactive feed and is rejected, not regularized.
// Inputs must already be normalized; this evaluator does not modify them.
[[nodiscard]] inline TpdPoint tangent_plane_distance(
    std::span<const double> w, std::span<const double> z,
    const StabilityPhase& trial, const StabilityPhase& reference) {
    if (w.size() != z.size()) { throw std::invalid_argument("TPD: dimension mismatch"); }
    const double sum = detail::stability_check_composition(w);
    (void)detail::stability_check_composition(z);
    detail::stability_check_phase(trial, w.size());
    detail::stability_check_phase(reference, w.size());
    TpdPoint result{std::vector<double>(w.begin(), w.end()), std::vector<double>(w.size()),
                    0, 0, 0, trial.branch, trial.smooth};
    double correction = 0, magnitude = 1, magnitude_correction = 0;
    bool boundary = false;
    for (std::size_t i = 0; i < w.size(); ++i) {
        if (z[i] == 0) {
            if (w[i] != 0) { throw std::domain_error("TPD: trial adds an absent feed component"); }
            continue;
        }
        if (w[i] == 0) { boundary = true; continue; }
        const double lw = std::log(w[i]), lz = std::log(z[i]);
        const double log_ratio = std::abs(w[i] - z[i]) < 0.5 * z[i]
            ? std::log1p((w[i] - z[i]) / z[i]) : lw - lz;
        const double q = log_ratio + (trial.ln_phi[i] - reference.ln_phi[i]);
        result.residual[i] = q;
        detail::stability_add(w[i] * q, result.value, correction);
        detail::stability_add(w[i] * (std::abs(lw) + std::abs(lz) +
            std::abs(trial.ln_phi[i]) + std::abs(reference.ln_phi[i])),
            magnitude, magnitude_correction);
    }
    const double mean = result.value / sum;
    for (std::size_t i = 0; i < w.size(); ++i) {
        if (z[i] == 0) { continue; }
        result.residual[i] = w[i] == 0 ? -std::numeric_limits<double>::infinity()
                                      : result.residual[i] - mean;
        result.stationarity = std::max(result.stationarity, std::abs(result.residual[i]));
    }
    result.roundoff_guard = 256.0 * detail::stability_eps * magnitude;
    if (!std::isfinite(result.value) || !std::isfinite(result.roundoff_guard) ||
        (!boundary && !std::isfinite(result.stationarity))) {
        throw StabilityPropertyError(StabilityPropertyIssue::nonfinite_properties,
                                     "TPD: nonrepresentable arithmetic");
    }
    return result;
}

// Provider is called synchronously as provider(p_Pa,T_K,normalized_w).
// It must be repeatable, must not retain spans, and must follow StabilityPhase's
// thermodynamic contract. Only StabilityPropertyError is converted to diagnostic
// failure; invalid inputs, programming errors and allocation exceptions propagate.
// double-only first increment; no AD through the search or solver sensitivities.
template <typename Provider>
[[nodiscard]] StabilityResult test_pt_stability(
    double pressure_pa, double temperature_k, std::span<const double> feed,
    Provider&& provider, StabilityOptions options = {},
    std::span<const std::vector<double>> extra_starts = {}) {
    detail::stability_check_options(options);
    if (!std::isfinite(pressure_pa) || pressure_pa <= 0 ||
        !std::isfinite(temperature_k) || temperature_k <= 0) {
        throw std::domain_error("stability: finite p>0 Pa and T>0 K required");
    }
    const std::size_t n = feed.size();
    if (n == 0 || n > options.max_components) {
        throw std::length_error("stability: invalid component count or quota exceeded");
    }
    const double input_sum = detail::stability_check_composition(feed);
    std::size_t active = 0;
    for (double value : feed) { if (value > 0) { ++active; } }
    if (active > std::numeric_limits<std::size_t>::max() - 2) {
        throw std::length_error("stability: start count overflow");
    }
    const std::size_t generated = options.automatic_starts ? (active == 1 ? 1 : active + 2) : 0;
    if (generated > options.max_starts || extra_starts.size() > options.max_starts - generated) {
        throw std::length_error("stability: start quota exceeded");
    }
    const std::size_t count = generated + extra_starts.size();
    if (count == 0 || count > options.max_start_entries / n ||
        count > std::vector<StabilityTrial>{}.max_size()) {
        throw std::length_error("stability: empty search or start storage quota exceeded");
    }
    // Validate ALL supplied starts before any provider call or search execution.
    for (const auto& start : extra_starts) {
        if (start.size() != n) { throw std::invalid_argument("stability: start dimension mismatch"); }
        (void)detail::stability_check_composition(start);
        for (std::size_t i = 0; i < n; ++i) {
            if ((feed[i] == 0 && start[i] != 0) || (feed[i] > 0 && start[i] == 0)) {
                throw std::domain_error("stability: starts must have exactly the active feed support");
            }
        }
    }
    StabilityResult result;
    result.options = options;
    result.pressure_pa = pressure_pa;
    result.temperature_k = temperature_k;
    result.input_feed_sum = input_sum;
    result.feed = detail::stability_normalize(feed, input_sum);
    for (std::size_t i = 0; i < n; ++i) {
        if (feed[i] > 0 && result.feed[i] == 0) {
            result.diagnostic = "roundoff normalization lost an active feed component";
            return result;
        }
    }
    result.trials.reserve(count);
    const auto add_start = [&](std::vector<double> start) {
        const double sum = detail::stability_sum(start);
        for (double& value : start) { value /= sum; }
        StabilityTrial trial;
        trial.initial_composition = std::move(start);
        result.trials.push_back(std::move(trial));
    };
    if (options.automatic_starts) {
        add_start(result.feed);
        if (active > 1) {
            std::vector<double> uniform(n);
            for (std::size_t i = 0; i < n; ++i) {
                if (result.feed[i] > 0) { uniform[i] = 1.0 / static_cast<double>(active); }
            }
            add_start(std::move(uniform));
            for (std::size_t i = 0; i < n; ++i) {
                if (result.feed[i] == 0) { continue; }
                std::vector<double> start = result.feed;
                for (double& value : start) { value *= 0.1; }
                start[i] += 0.9;
                add_start(std::move(start));
            }
        }
    }
    for (const auto& start : extra_starts) { add_start(start); }
    try {
        ++result.evaluations;
        result.reference = provider(pressure_pa, temperature_k, std::span<const double>{result.feed});
        detail::stability_check_phase(*result.reference, n);
    } catch (const StabilityPropertyError& error) {
        result.reference.reset();
        result.reference_issue = error.issue();
        result.diagnostic = error.what();
        for (auto& trial : result.trials) {
            trial.status = StabilityTrialStatus::property_failure;
            trial.property_issue = error.issue();
            trial.diagnostic = "reference phase unavailable";
        }
        return result;
    }
    const auto remember = [&](const TpdPoint& point) {
        if (!result.lowest_sampled || point.value < result.lowest_sampled->value) {
            result.lowest_sampled = point;
        }
    };
    bool negative_found = false, all_stationary = true;
    for (auto& trial : result.trials) {
        const auto evaluate = [&](std::span<const double> w) {
            ++result.evaluations;
            ++trial.evaluations;
            const StabilityPhase phase = provider(pressure_pa, temperature_k, w);
            return tangent_plane_distance(w, result.feed, phase, *result.reference);
        };
        bool support_lost = false;
        for (std::size_t i = 0; i < n; ++i) {
            support_lost = support_lost || (result.feed[i] > 0 && trial.initial_composition[i] == 0);
        }
        if (support_lost) {
            trial.status = StabilityTrialStatus::unrepresentable_composition;
            trial.diagnostic = "generated start lost an active component; no floor applied";
            all_stationary = false;
            continue;
        }
        if (result.evaluations >= options.max_evaluations) {
            trial.status = StabilityTrialStatus::evaluation_limit;
            all_stationary = false;
            continue;
        }
        try { trial.point = evaluate(trial.initial_composition); }
        catch (const StabilityPropertyError& error) {
            trial.status = StabilityTrialStatus::property_failure;
            trial.property_issue = error.issue();
            trial.diagnostic = error.what();
            all_stationary = false;
            continue;
        }
        for (;;) {
            remember(*trial.point);
            if (detail::stability_negative(*trial.point, options)) {
                trial.status = StabilityTrialStatus::negative_tpd;
                negative_found = true;
                break;
            }
            if (!trial.point->smooth) {
                trial.status = StabilityTrialStatus::nonsmooth;
                break;
            }
            if (trial.point->stationarity <= options.stationarity_tolerance) {
                trial.status = StabilityTrialStatus::stationary; // May be a saddle, NOT certified minimum.
                break;
            }
            if (trial.iterations >= options.max_iterations) {
                trial.status = StabilityTrialStatus::iteration_limit;
                break;
            }
            const auto& point = *trial.point;
            const double scale = std::max(1.0, point.stationarity / options.log_step_limit);
            std::vector<double> direction(n), candidate(n);
            double slope = 0, correction = 0;
            for (std::size_t i = 0; i < n; ++i) {
                if (result.feed[i] == 0) { continue; }
                direction[i] = -point.residual[i] / scale;
                detail::stability_add((point.composition[i] * point.residual[i]) * direction[i],
                                      slope, correction);
            }
            if (!std::isfinite(slope) || !(slope < 0)) {
                trial.status = StabilityTrialStatus::line_search_failed;
                trial.diagnostic = "nonrepresentable or zero descent slope with nonzero stationarity residual";
                break;
            }
            bool accepted = false;
            double alpha = 1;
            trial.status = StabilityTrialStatus::line_search_failed;
            for (int backtrack = 0; backtrack < options.max_backtracks; ++backtrack, alpha *= 0.5) {
                if (backtrack > 0) { ++trial.backtracks; }
                double largest = -std::numeric_limits<double>::infinity();
                for (std::size_t i = 0; i < n; ++i) {
                    if (result.feed[i] > 0) {
                        candidate[i] = std::log(point.composition[i]) + alpha * direction[i];
                        largest = std::max(largest, candidate[i]);
                    }
                }
                for (std::size_t i = 0; i < n; ++i) {
                    candidate[i] = result.feed[i] > 0 ? std::exp(candidate[i] - largest) : 0;
                }
                const double sum = detail::stability_sum(candidate);
                bool positive = std::isfinite(sum) && sum > 0;
                for (std::size_t i = 0; i < n; ++i) {
                    candidate[i] /= sum;
                    positive = positive && (result.feed[i] == 0 || candidate[i] > 0);
                }
                if (!positive) { continue; } // Reject, never replace by epsilon.
                if (result.evaluations >= options.max_evaluations) {
                    trial.status = StabilityTrialStatus::evaluation_limit;
                    break;
                }
                try {
                    TpdPoint next = evaluate(candidate);
                    remember(next); // Any feasible negative point is evidence, even before convergence.
                    if (detail::stability_negative(next, options)) {
                        trial.point = std::move(next);
                        accepted = true;
                        break;
                    }
                    const double arithmetic_guard = point.roundoff_guard + next.roundoff_guard;
                    const double predicted = -alpha * slope;
                    const bool roundoff_progress = predicted <= arithmetic_guard &&
                        next.value <= point.value + arithmetic_guard &&
                        next.stationarity < 0.9 * point.stationarity;
                    if (next.value <= point.value + options.armijo * alpha * slope || roundoff_progress) {
                        trial.point = std::move(next);
                        accepted = true;
                        break;
                    }
                } catch (const StabilityPropertyError& error) {
                    ++trial.rejected_property_evaluations;
                    trial.property_issue = error.issue();
                    trial.diagnostic = error.what();
                }
            }
            if (!accepted) { break; }
            ++trial.iterations;
        }
        all_stationary = all_stationary && trial.status == StabilityTrialStatus::stationary;
    }
    result.status = negative_found ? StabilityStatus::unstable :
        (all_stationary ? StabilityStatus::no_instability_found : StabilityStatus::indeterminate);
    return result;
}

} // namespace mpmc::flash
#endif // MPMC_FLASH_PT_STABILITY_HPP
