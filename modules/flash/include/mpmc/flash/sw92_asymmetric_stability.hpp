#ifndef MPMC_FLASH_SW92_ASYMMETRIC_STABILITY_HPP
#define MPMC_FLASH_SW92_ASYMMETRIC_STABILITY_HPP

#include <mpmc/flash/sw92_stability.hpp>

#include <cmath>
#include <cstddef>
#include <limits>
#include <optional>
#include <span>
#include <stdexcept>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

namespace mpmc::flash {

/// Reserved parent equilibrium profile. Gate 3A implements only its asymmetric
/// stability foundation; it does not publish or solve a joint AQ/NA phase split.
inline constexpr std::string_view sw92_xu_asymmetric_gibbs_profile =
    "SW92-equilibrium/xu-asymmetric-gibbs/v1";

/// Numerical identity of this finite two-family common-tangent search.
inline constexpr std::string_view sw92_xu_asymmetric_stability_convention =
    "SW92-equilibrium/xu-asymmetric-gibbs/stability-foundation/v1";

struct Sw92AsymmetricFamilyOptions {
    thermodynamics::Sw92RootOptions root_options;
    StabilityOptions stability;
};

struct Sw92AsymmetricStabilityOptions {
    Sw92AsymmetricFamilyOptions aqueous;
    Sw92AsymmetricFamilyOptions nonaqueous;
};

/// Non-owning extra-start views used only during one stability call.
struct Sw92AsymmetricStabilityStarts {
    std::span<const std::vector<double>> aqueous{};
    std::span<const std::vector<double>> nonaqueous{};
};

enum class Sw92AsymmetricFeedReferenceStatus {
    selected,
    family_gibbs_tie,
    family_reference_failure,
    family_reference_nonsmooth,
    normalization_failure
};

/// One required family leg of the asymmetric search. Feed-reference evaluation
/// and finite TPD-search accounting remain separate so neither family consumes
/// the other family's resource budget.
struct Sw92AsymmetricFamilySearchResult {
    thermodynamics::SwPhaseFamily family{thermodynamics::SwPhaseFamily::aqueous};
    std::size_t feed_reference_evaluations{};
    std::optional<StabilityPhase> feed_reference;
    std::optional<StabilityPropertyIssue> feed_reference_issue;
    double feed_reduced_gibbs{std::numeric_limits<double>::quiet_NaN()};
    std::optional<StabilityResult> search;

    [[nodiscard]] std::optional<std::size_t> total_property_evaluations() const noexcept {
        const std::size_t search_count = search ? search->evaluations : 0;
        if (search_count > std::numeric_limits<std::size_t>::max() -
                               feed_reference_evaluations) {
            return std::nullopt;
        }
        return feed_reference_evaluations + search_count;
    }
};

struct Sw92AsymmetricNegativeWitness {
    thermodynamics::SwPhaseFamily family{thermodynamics::SwPhaseFamily::aqueous};
    std::size_t trial_index{};
    TpdPoint point;
};

/// Two required finite family searches against one externally supplied common
/// reduced tangent. This is reusable by Gate 3A and later phase-set review; it
/// never decides or publishes a phase split.
struct Sw92AsymmetricCommonTangentSearchResult {
    StabilityStatus status{StabilityStatus::indeterminate};
    static constexpr bool global_stability_proven = false;

    double pressure_pa{};
    double temperature_k{};
    double input_feed_sum{};
    std::vector<double> feed;
    double nacl_molality_mol_per_kg_water{};

    std::string dataset_id;
    std::string revision;
    std::vector<std::string> component_ids;
    std::string model_profile{thermodynamics::sw92_corrected_profile};
    std::string phase_convention{thermodynamics::sw92_pt_convention};
    std::string equilibrium_profile{sw92_xu_asymmetric_gibbs_profile};
    std::string search_convention{sw92_xu_asymmetric_stability_convention};

    Sw92AsymmetricStabilityOptions options;
    std::vector<double> common_log_activity;
    StabilityResult aqueous;
    StabilityResult nonaqueous;
    std::vector<Sw92AsymmetricNegativeWitness> negative_witnesses;
    std::string diagnostic;
};

struct Sw92AsymmetricStabilityResult {
    StabilityStatus status{StabilityStatus::indeterminate};
    static constexpr bool global_stability_proven = false;

    double pressure_pa{};
    double temperature_k{};
    double input_feed_sum{};
    std::vector<double> feed;
    double nacl_molality_mol_per_kg_water{};

    std::string dataset_id;
    std::string revision;
    std::vector<std::string> component_ids;
    std::string model_profile{thermodynamics::sw92_corrected_profile};
    std::string phase_convention{thermodynamics::sw92_pt_convention};
    std::string equilibrium_profile{sw92_xu_asymmetric_gibbs_profile};
    std::string search_convention{sw92_xu_asymmetric_stability_convention};

    Sw92AsymmetricStabilityOptions options;
    Sw92AsymmetricFeedReferenceStatus feed_reference_status{
        Sw92AsymmetricFeedReferenceStatus::family_reference_failure};
    std::optional<thermodynamics::SwPhaseFamily> reference_family;
    double aqueous_minus_nonaqueous_feed_gibbs{
        std::numeric_limits<double>::quiet_NaN()};
    double feed_gibbs_roundoff_guard{
        std::numeric_limits<double>::quiet_NaN()};
    std::vector<double> common_log_activity;

    Sw92AsymmetricFamilySearchResult aqueous;
    Sw92AsymmetricFamilySearchResult nonaqueous;
    std::vector<Sw92AsymmetricNegativeWitness> negative_witnesses;
    std::string diagnostic;
};

namespace detail {

inline const char* sw92_family_name(thermodynamics::SwPhaseFamily family) noexcept {
    return family == thermodynamics::SwPhaseFamily::aqueous ? "AQ" : "NA";
}

inline void sw92_asymmetric_preflight_search(
    std::span<const double> feed, const StabilityOptions& options,
    std::span<const std::vector<double>> extra_starts) {
    stability_check_options(options);
    const std::size_t n = feed.size();
    if (n == 0 || n > options.max_components) {
        throw std::length_error("SW92 asymmetric stability: component quota exceeded or empty feed");
    }
    (void)stability_check_composition(feed);
    std::size_t active = 0;
    for (double value : feed) {
        if (value > 0.0) { ++active; }
    }
    if (active > std::numeric_limits<std::size_t>::max() - 2) {
        throw std::length_error("SW92 asymmetric stability: start count overflow");
    }
    const std::size_t generated = options.automatic_starts
        ? (active == 1 ? 1 : active + 2) : 0;
    if (generated > options.max_starts ||
        extra_starts.size() > options.max_starts - generated) {
        throw std::length_error("SW92 asymmetric stability: start quota exceeded");
    }
    const std::size_t count = generated + extra_starts.size();
    if (count == 0 || count > options.max_start_entries / n ||
        count > std::vector<StabilityTrial>{}.max_size()) {
        throw std::length_error("SW92 asymmetric stability: empty search or start storage quota exceeded");
    }
    for (const auto& start : extra_starts) {
        if (start.size() != n) {
            throw std::invalid_argument("SW92 asymmetric stability: start dimension mismatch");
        }
        (void)stability_check_composition(start);
        for (std::size_t i = 0; i < n; ++i) {
            if ((feed[i] == 0.0 && start[i] != 0.0) ||
                (feed[i] > 0.0 && start[i] == 0.0)) {
                throw std::domain_error(
                    "SW92 asymmetric stability: starts must match active feed support");
            }
        }
    }
}

inline double sw92_feed_reduced_gibbs(
    std::span<const double> feed, const StabilityPhase& phase) {
    stability_check_phase(phase, feed.size());
    double value = 0.0;
    double correction = 0.0;
    for (std::size_t i = 0; i < feed.size(); ++i) {
        if (feed[i] == 0.0) { continue; }
        const double term = feed[i] * (std::log(feed[i]) + phase.ln_phi[i]);
        stability_add(term, value, correction);
    }
    if (!std::isfinite(value)) {
        throw StabilityPropertyError(StabilityPropertyIssue::nonfinite_properties,
            "SW92 asymmetric stability: nonrepresentable feed Gibbs arithmetic");
    }
    return value;
}

inline std::pair<double, double> sw92_feed_family_difference(
    std::span<const double> feed, const StabilityPhase& aqueous,
    const StabilityPhase& nonaqueous) {
    stability_check_phase(aqueous, feed.size());
    stability_check_phase(nonaqueous, feed.size());
    double difference = 0.0;
    double correction = 0.0;
    double magnitude = 1.0;
    double magnitude_correction = 0.0;
    for (std::size_t i = 0; i < feed.size(); ++i) {
        if (feed[i] == 0.0) { continue; }
        stability_add(feed[i] * (aqueous.ln_phi[i] - nonaqueous.ln_phi[i]),
                      difference, correction);
        stability_add(feed[i] * (std::abs(aqueous.ln_phi[i]) +
                                 std::abs(nonaqueous.ln_phi[i])),
                      magnitude, magnitude_correction);
    }
    const double guard = 256.0 * stability_eps * magnitude;
    if (!std::isfinite(difference) || !std::isfinite(guard)) {
        throw StabilityPropertyError(StabilityPropertyIssue::nonfinite_properties,
            "SW92 asymmetric stability: nonrepresentable feed-family Gibbs comparison");
    }
    return {difference, guard};
}

inline void sw92_collect_negative_witnesses(
    thermodynamics::SwPhaseFamily family, const StabilityResult& search,
    std::vector<Sw92AsymmetricNegativeWitness>& output) {
    for (std::size_t trial_index = 0; trial_index < search.trials.size(); ++trial_index) {
        const auto& trial = search.trials[trial_index];
        if (trial.status != StabilityTrialStatus::negative_tpd || !trial.point ||
            !stability_negative(*trial.point, search.options)) {
            continue;
        }
        output.push_back({family, trial_index, *trial.point});
    }
}

inline StabilityStatus sw92_combine_asymmetric_searches(
    const StabilityResult& aqueous, const StabilityResult& nonaqueous) noexcept {
    if (aqueous.status == StabilityStatus::unstable ||
        nonaqueous.status == StabilityStatus::unstable) {
        return StabilityStatus::unstable;
    }
    if (aqueous.status == StabilityStatus::no_instability_found &&
        nonaqueous.status == StabilityStatus::no_instability_found) {
        return StabilityStatus::no_instability_found;
    }
    return StabilityStatus::indeterminate;
}

inline std::string sw92_asymmetric_search_diagnostic(StabilityStatus status) {
    switch (status) {
    case StabilityStatus::unstable:
        return "finite common-tangent search found a robust negative TPD witness in at least one SW92 family";
    case StabilityStatus::no_instability_found:
        return "both finite family searches found no instability against the common tangent; not a global proof";
    case StabilityStatus::indeterminate:
        return "no robust negative witness was found, but at least one required family search is indeterminate";
    }
    return "unknown SW92 asymmetric stability status";
}

} // namespace detail

/// Run both SW92 family searches against exactly one caller-supplied reduced
/// tangent. The caller owns the thermodynamic validity/uncertainty of that
/// tangent; this helper owns family symmetry, resource isolation and witnesses.
[[nodiscard]] inline Sw92AsymmetricCommonTangentSearchResult
    test_sw92_pt_asymmetric_stability_against(
        double pressure_pa, double temperature_k, std::span<const double> support_feed,
        std::span<const double> common_log_activity,
        const thermodynamics::Sw92Phase<double>& model,
        double nacl_molality_mol_per_kg_water,
        Sw92AsymmetricStabilityOptions options = {},
        Sw92AsymmetricStabilityStarts starts = {}) {
    if (!std::isfinite(pressure_pa) || pressure_pa <= 0.0 ||
        !std::isfinite(temperature_k) || temperature_k <= 0.0) {
        throw std::domain_error(
            "SW92 asymmetric stability: finite p>0 Pa and T>0 K required");
    }
    if (support_feed.size() != model.size()) {
        throw std::invalid_argument(
            "SW92 asymmetric stability: feed does not match ordered model snapshot");
    }
    if (common_log_activity.size() != support_feed.size()) {
        throw std::invalid_argument(
            "SW92 asymmetric stability: common tangent dimension mismatch");
    }
    for (double value : common_log_activity) {
        if (!std::isfinite(value)) {
            throw std::domain_error(
                "SW92 asymmetric stability: finite common tangent entries required");
        }
    }

    detail::sw92_asymmetric_preflight_search(
        support_feed, options.aqueous.stability, starts.aqueous);
    detail::sw92_asymmetric_preflight_search(
        support_feed, options.nonaqueous.stability, starts.nonaqueous);

    Sw92AsymmetricCommonTangentSearchResult result;
    result.pressure_pa = pressure_pa;
    result.temperature_k = temperature_k;
    result.input_feed_sum = detail::stability_check_composition(support_feed);
    result.feed = detail::stability_normalize(support_feed, result.input_feed_sum);
    result.nacl_molality_mol_per_kg_water = nacl_molality_mol_per_kg_water;
    result.options = options;
    result.common_log_activity.assign(common_log_activity.begin(), common_log_activity.end());

    const auto& parameters = model.parameters();
    result.dataset_id = parameters.dataset_id();
    result.revision = parameters.revision();
    for (const auto& component : parameters.components().items()) {
        result.component_ids.push_back(component.id);
    }

    Sw92FamilyStabilityEvaluator aqueous_evaluator(
        model, nacl_molality_mol_per_kg_water,
        thermodynamics::SwPhaseFamily::aqueous, options.aqueous.root_options);
    Sw92FamilyStabilityEvaluator nonaqueous_evaluator(
        model, nacl_molality_mol_per_kg_water,
        thermodynamics::SwPhaseFamily::nonaqueous, options.nonaqueous.root_options);

    result.aqueous = test_pt_stability_against(
        pressure_pa, temperature_k, support_feed, common_log_activity,
        aqueous_evaluator, options.aqueous.stability, starts.aqueous);
    result.nonaqueous = test_pt_stability_against(
        pressure_pa, temperature_k, support_feed, common_log_activity,
        nonaqueous_evaluator, options.nonaqueous.stability, starts.nonaqueous);

    detail::sw92_collect_negative_witnesses(
        thermodynamics::SwPhaseFamily::aqueous, result.aqueous,
        result.negative_witnesses);
    detail::sw92_collect_negative_witnesses(
        thermodynamics::SwPhaseFamily::nonaqueous, result.nonaqueous,
        result.negative_witnesses);
    result.status = detail::sw92_combine_asymmetric_searches(
        result.aqueous, result.nonaqueous);
    result.diagnostic = detail::sw92_asymmetric_search_diagnostic(result.status);
    return result;
}

/// Gate 3A: finite Xu-style lower-envelope stability foundation for SW92 AQ/NA.
///
/// The function evaluates both same-family minimum-Gibbs feed references, chooses
/// a resolved lower family only when the AQ/NA Gibbs gap exceeds the arithmetic
/// guard, constructs one common reduced tangent, and runs BOTH family TPD searches
/// against that exact tangent. It never solves or publishes a joint phase split.
[[nodiscard]] inline Sw92AsymmetricStabilityResult test_sw92_pt_asymmetric_stability(
    double pressure_pa, double temperature_k, std::span<const double> feed,
    const thermodynamics::Sw92Phase<double>& model,
    double nacl_molality_mol_per_kg_water,
    Sw92AsymmetricStabilityOptions options = {},
    Sw92AsymmetricStabilityStarts starts = {}) {
    if (!std::isfinite(pressure_pa) || pressure_pa <= 0.0 ||
        !std::isfinite(temperature_k) || temperature_k <= 0.0) {
        throw std::domain_error(
            "SW92 asymmetric stability: finite p>0 Pa and T>0 K required");
    }
    if (feed.size() != model.size()) {
        throw std::invalid_argument(
            "SW92 asymmetric stability: feed does not match ordered model snapshot");
    }

    // Validate both family search contracts before any thermodynamic evaluation,
    // so a bad NA budget/start cannot be hidden behind earlier AQ work or vice versa.
    detail::sw92_asymmetric_preflight_search(feed, options.aqueous.stability,
                                              starts.aqueous);
    detail::sw92_asymmetric_preflight_search(feed, options.nonaqueous.stability,
                                              starts.nonaqueous);

    Sw92AsymmetricStabilityResult result;
    result.pressure_pa = pressure_pa;
    result.temperature_k = temperature_k;
    result.input_feed_sum = detail::stability_check_composition(feed);
    result.feed = detail::stability_normalize(feed, result.input_feed_sum);
    result.nacl_molality_mol_per_kg_water = nacl_molality_mol_per_kg_water;
    result.options = options;
    result.aqueous.family = thermodynamics::SwPhaseFamily::aqueous;
    result.nonaqueous.family = thermodynamics::SwPhaseFamily::nonaqueous;

    const auto& parameters = model.parameters();
    result.dataset_id = parameters.dataset_id();
    result.revision = parameters.revision();
    for (const auto& component : parameters.components().items()) {
        result.component_ids.push_back(component.id);
    }

    for (std::size_t i = 0; i < feed.size(); ++i) {
        if (feed[i] > 0.0 && result.feed[i] == 0.0) {
            result.feed_reference_status =
                Sw92AsymmetricFeedReferenceStatus::normalization_failure;
            result.diagnostic =
                "roundoff normalization lost an active feed component";
            return result;
        }
    }

    Sw92FamilyStabilityEvaluator aqueous_evaluator(
        model, nacl_molality_mol_per_kg_water,
        thermodynamics::SwPhaseFamily::aqueous, options.aqueous.root_options);
    Sw92FamilyStabilityEvaluator nonaqueous_evaluator(
        model, nacl_molality_mol_per_kg_water,
        thermodynamics::SwPhaseFamily::nonaqueous, options.nonaqueous.root_options);

    const auto evaluate_feed = [&](auto& evaluator,
                                   Sw92AsymmetricFamilySearchResult& family_result) {
        ++family_result.feed_reference_evaluations;
        try {
            family_result.feed_reference = evaluator(
                pressure_pa, temperature_k, std::span<const double>{result.feed});
            detail::stability_check_phase(*family_result.feed_reference, result.feed.size());
            family_result.feed_reduced_gibbs = detail::sw92_feed_reduced_gibbs(
                result.feed, *family_result.feed_reference);
        } catch (const StabilityPropertyError& error) {
            family_result.feed_reference.reset();
            family_result.feed_reference_issue = error.issue();
        }
    };

    // Both required feed-family evaluations are attempted independently.
    evaluate_feed(aqueous_evaluator, result.aqueous);
    evaluate_feed(nonaqueous_evaluator, result.nonaqueous);

    if (!result.aqueous.feed_reference || !result.nonaqueous.feed_reference) {
        result.feed_reference_status =
            Sw92AsymmetricFeedReferenceStatus::family_reference_failure;
        result.diagnostic =
            "AQ/NA common feed tangent unavailable because at least one family reference failed";
        return result;
    }
    if (!result.aqueous.feed_reference->smooth ||
        !result.nonaqueous.feed_reference->smooth) {
        result.feed_reference_status =
            Sw92AsymmetricFeedReferenceStatus::family_reference_nonsmooth;
        result.diagnostic =
            "AQ/NA common feed tangent unavailable because a same-family minimum envelope is nonsmooth";
        return result;
    }

    try {
        const auto [difference, guard] = detail::sw92_feed_family_difference(
            result.feed, *result.aqueous.feed_reference,
            *result.nonaqueous.feed_reference);
        result.aqueous_minus_nonaqueous_feed_gibbs = difference;
        result.feed_gibbs_roundoff_guard = guard;
    } catch (const StabilityPropertyError& error) {
        result.feed_reference_status =
            Sw92AsymmetricFeedReferenceStatus::family_reference_failure;
        result.diagnostic = error.what();
        return result;
    }

    if (std::abs(result.aqueous_minus_nonaqueous_feed_gibbs) <=
        result.feed_gibbs_roundoff_guard) {
        result.feed_reference_status =
            Sw92AsymmetricFeedReferenceStatus::family_gibbs_tie;
        result.diagnostic =
            "AQ/NA feed Gibbs surfaces are tied at the arithmetic resolution; lower-envelope tangent is nonsmooth";
        return result;
    }

    const StabilityPhase* selected = nullptr;
    if (result.aqueous_minus_nonaqueous_feed_gibbs < 0.0) {
        result.reference_family = thermodynamics::SwPhaseFamily::aqueous;
        selected = &*result.aqueous.feed_reference;
    } else {
        result.reference_family = thermodynamics::SwPhaseFamily::nonaqueous;
        selected = &*result.nonaqueous.feed_reference;
    }
    result.feed_reference_status = Sw92AsymmetricFeedReferenceStatus::selected;
    result.common_log_activity.resize(result.feed.size());
    for (std::size_t i = 0; i < result.feed.size(); ++i) {
        result.common_log_activity[i] = result.feed[i] > 0.0
            ? std::log(result.feed[i]) + selected->ln_phi[i]
            : 0.0;
        if (!std::isfinite(result.common_log_activity[i])) {
            result.status = StabilityStatus::indeterminate;
            result.diagnostic = "common reduced feed tangent is nonrepresentable";
            result.common_log_activity.clear();
            return result;
        }
    }

    auto common_search = test_sw92_pt_asymmetric_stability_against(
        pressure_pa, temperature_k, result.feed, result.common_log_activity,
        model, nacl_molality_mol_per_kg_water, options, starts);
    result.aqueous.search = std::move(common_search.aqueous);
    result.nonaqueous.search = std::move(common_search.nonaqueous);
    result.negative_witnesses = std::move(common_search.negative_witnesses);
    result.status = common_search.status;
    result.diagnostic = std::move(common_search.diagnostic);
    return result;
}

} // namespace mpmc::flash

#endif // MPMC_FLASH_SW92_ASYMMETRIC_STABILITY_HPP
