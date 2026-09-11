#ifndef MPMC_FLASH_SW92_PHASE_ASSIGNED_NO_W_HPP
#define MPMC_FLASH_SW92_PHASE_ASSIGNED_NO_W_HPP

#include <mpmc/flash/sw92_phase_assigned_joint.hpp>
#include <mpmc/flash/sw92_split.hpp>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <numeric>
#include <optional>
#include <span>
#include <stdexcept>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

namespace mpmc::flash {

inline constexpr std::string_view sw92_phase_assigned_no_w_convention =
    "SW92-equilibrium/phase-assigned-aq-na-joint/no-w-na-flash-"
    "targeted-aq-appearance/v1";

struct Sw92PhaseAssignedNoWOptions {
    PtSplitOptions nonaqueous_flash;
    StabilityOptions aqueous_appearance;
    double log_composition_separation{1e-7};
    thermodynamics::Sw92RootOptions aqueous_root_options;
    thermodynamics::Sw92RootOptions nonaqueous_root_options;

    Sw92PhaseAssignedNoWOptions() {
        // Profile-C water appearance is a role-targeted search. Do not silently
        // turn it into whole-simplex AQ/NA model competition by default.
        aqueous_appearance.automatic_starts = false;
    }
};

enum class Sw92PhaseAssignedNoWStatus {
    no_w_single_h_locally_closed,
    no_w_two_h_locally_closed,
    aqueous_phase_witness_found,
    higher_h_multiplicity_or_wrong_candidate,
    phase_disappearance_unresolved,
    single_phase_role_unresolved,
    indeterminate
};

struct Sw92PhaseAssignedWaterAppearanceWitness {
    std::size_t trial_index{};
    TpdPoint point;
    std::vector<double> log_distance_from_retained_na_candidates;
    double trial_minus_min_na_water_fraction{};
    double water_role_roundoff_guard{};
    bool distinct_from_every_retained_na_candidate{};
    bool water_role_candidate_admissible{};

    /// Candidate-generation evidence only. A later joint W+H... solve must still
    /// require W to be water-richer than every FINAL retained H phase.
    [[nodiscard]] bool usable_w_seed() const noexcept {
        return distinct_from_every_retained_na_candidate &&
               water_role_candidate_admissible;
    }
};

struct Sw92PhaseAssignedNoWResult {
    static constexpr bool global_stability_proven = false;
    static constexpr bool accepted_phase_set_published = false;
    static constexpr bool morphology_resolved = false;

    Sw92PhaseAssignedNoWStatus status{Sw92PhaseAssignedNoWStatus::indeterminate};
    Sw92PhaseAssignedNoWOptions options;
    Sw92FamilyPtSplitResult hydrocarbon_flash;
    std::optional<StabilityResult> aqueous_appearance_search;
    std::vector<Sw92PhaseAssignedWaterAppearanceWitness> aqueous_witnesses;

    double pressure_pa{};
    double temperature_k{};
    double input_feed_sum{};
    std::vector<double> feed;
    double nacl_molality_mol_per_kg_water{};
    std::size_t water_index{};
    std::string dataset_id;
    std::string revision;
    std::vector<std::string> component_ids;
    std::string model_profile{thermodynamics::sw92_corrected_profile};
    std::string phase_convention{thermodynamics::sw92_pt_convention};
    std::string equilibrium_profile{sw92_phase_assigned_aq_na_joint_profile};
    std::string adapter_convention{sw92_phase_assigned_no_w_convention};

    // These are retained fixed-NA numerical candidates. Before a rival W
    // topology is solved they MUST NOT be interpreted as authoritative physical
    // hydrocarbon roles merely because they use the NA thermodynamic family.
    std::vector<std::vector<double>> retained_na_candidate_compositions;
    std::vector<double> retained_na_candidate_fractions;
    std::vector<double> common_log_activity;
    double min_retained_na_candidate_water_fraction{
        std::numeric_limits<double>::quiet_NaN()};
    double max_retained_na_candidate_water_fraction{
        std::numeric_limits<double>::quiet_NaN()};
    std::optional<std::vector<double>> targeted_water_start;
    std::optional<std::size_t> selected_water_witness_index;
    std::string diagnostic;

    [[nodiscard]] std::size_t retained_na_candidate_count() const noexcept {
        return retained_na_candidate_compositions.size();
    }
    [[nodiscard]] const Sw92PhaseAssignedWaterAppearanceWitness*
    selected_water_witness() const & noexcept {
        if (!selected_water_witness_index ||
            *selected_water_witness_index >= aqueous_witnesses.size()) {
            return nullptr;
        }
        return &aqueous_witnesses[*selected_water_witness_index];
    }
    const Sw92PhaseAssignedWaterAppearanceWitness* selected_water_witness() const && = delete;
    [[nodiscard]] bool no_w_locally_closed() const noexcept {
        return status == Sw92PhaseAssignedNoWStatus::no_w_single_h_locally_closed ||
               status == Sw92PhaseAssignedNoWStatus::no_w_two_h_locally_closed;
    }
};

namespace detail {

inline void sw92_phase_assigned_no_w_check_options(
    const Sw92PhaseAssignedNoWOptions& options) {
    split_check_options(options.nonaqueous_flash.iteration);
    stability_check_options(options.nonaqueous_flash.initial_stability);
    stability_check_options(options.nonaqueous_flash.final_stability);
    stability_check_options(options.aqueous_appearance);
    if (!std::isfinite(options.log_composition_separation) ||
        !(options.log_composition_separation > 0.0) ||
        options.aqueous_root_options.max_iterations <= 0 ||
        options.nonaqueous_root_options.max_iterations <= 0) {
        throw std::invalid_argument(
            "SW92 Profile-C no-W: invalid separation/root option");
    }
}

inline double sw92_phase_assigned_no_w_log_distance(
    std::span<const double> first, std::span<const double> second,
    std::span<const double> feed) {
    if (first.size() != second.size() || first.size() != feed.size()) {
        return std::numeric_limits<double>::infinity();
    }
    double distance = 0.0;
    for (std::size_t i = 0; i < feed.size(); ++i) {
        if (feed[i] == 0.0) {
            if (first[i] != 0.0 || second[i] != 0.0) {
                return std::numeric_limits<double>::infinity();
            }
            continue;
        }
        if (!(first[i] > 0.0) || !(second[i] > 0.0)) {
            return std::numeric_limits<double>::infinity();
        }
        distance = std::max(
            distance,
            std::abs(std::log(first[i]) - std::log(second[i])));
    }
    return distance;
}

inline std::optional<std::vector<double>> sw92_phase_assigned_targeted_water_start(
    std::span<const double> feed, std::size_t water_index,
    double max_na_candidate_water_fraction) {
    if (water_index >= feed.size() || !(feed[water_index] > 0.0) ||
        !std::isfinite(max_na_candidate_water_fraction) ||
        !(max_na_candidate_water_fraction >= 0.0) ||
        !(max_na_candidate_water_fraction < 1.0)) {
        return std::nullopt;
    }
    double nonwater_feed = 0.0;
    double correction = 0.0;
    for (std::size_t i = 0; i < feed.size(); ++i) {
        if (i != water_index) {
            stability_add(feed[i], nonwater_feed, correction);
        }
    }
    if (!(nonwater_feed > 0.0) || !std::isfinite(nonwater_feed)) {
        return std::nullopt;
    }

    const double target_water = std::midpoint(max_na_candidate_water_fraction, 1.0);
    const double guard = 256.0 * stability_eps *
        (1.0 + std::abs(max_na_candidate_water_fraction) +
         std::abs(target_water));
    if (!std::isfinite(target_water) || !(target_water < 1.0) ||
        !(target_water - max_na_candidate_water_fraction > guard)) {
        return std::nullopt;
    }

    std::vector<double> start(feed.size(), 0.0);
    start[water_index] = target_water;
    const double remainder = 1.0 - target_water;
    for (std::size_t i = 0; i < feed.size(); ++i) {
        if (i == water_index || feed[i] == 0.0) { continue; }
        start[i] = remainder * feed[i] / nonwater_feed;
        if (!(start[i] > 0.0) || !std::isfinite(start[i])) {
            return std::nullopt;
        }
    }
    const double sum = stability_sum(start);
    if (!std::isfinite(sum) || std::abs(sum - 1.0) > 64.0 * stability_eps) {
        return std::nullopt;
    }
    for (double& value : start) { value /= sum; }
    try {
        (void)stability_check_composition(start);
    } catch (const std::exception&) {
        return std::nullopt;
    }
    return start;
}

inline void sw92_phase_assigned_no_w_preflight_starts(
    std::span<const double> feed, const StabilityOptions& options,
    std::size_t mandatory, std::span<const std::vector<double>> extra_starts) {
    const std::size_t n = feed.size();
    if (n == 0U || n > options.max_components) {
        throw std::length_error(
            "SW92 Profile-C no-W: component quota exceeded or empty feed");
    }
    std::size_t active = 0U;
    for (double value : feed) { if (value > 0.0) { ++active; } }
    const std::size_t generated = options.automatic_starts
        ? (active == 1U ? 1U : active + 2U) : 0U;
    if (mandatory > options.max_starts || generated > options.max_starts - mandatory ||
        extra_starts.size() > options.max_starts - mandatory - generated) {
        throw std::length_error("SW92 Profile-C no-W: start quota exceeded");
    }
    const std::size_t total = mandatory + generated + extra_starts.size();
    if (total > options.max_start_entries / n ||
        total > std::vector<StabilityTrial>{}.max_size()) {
        throw std::length_error("SW92 Profile-C no-W: start storage quota exceeded");
    }
    for (const auto& start : extra_starts) {
        if (start.size() != n) {
            throw std::invalid_argument(
                "SW92 Profile-C no-W: AQ start dimension mismatch");
        }
        (void)stability_check_composition(start);
        for (std::size_t i = 0; i < n; ++i) {
            if ((feed[i] > 0.0) != (start[i] > 0.0)) {
                throw std::domain_error(
                    "SW92 Profile-C no-W: AQ starts must match active feed support");
            }
        }
    }
}

inline Sw92PhaseAssignedWaterAppearanceWitness
sw92_phase_assigned_classify_water_witness(
    std::size_t trial_index, const TpdPoint& point,
    std::span<const std::vector<double>> retained_na_candidates,
    std::span<const double> feed, std::size_t water_index,
    double min_na_candidate_water_fraction,
    double log_composition_separation) {
    Sw92PhaseAssignedWaterAppearanceWitness witness;
    witness.trial_index = trial_index;
    witness.point = point;
    witness.distinct_from_every_retained_na_candidate = true;
    for (const auto& candidate : retained_na_candidates) {
        const double distance = sw92_phase_assigned_no_w_log_distance(
            point.composition, candidate, feed);
        witness.log_distance_from_retained_na_candidates.push_back(distance);
        if (!std::isfinite(distance) || distance <= log_composition_separation) {
            witness.distinct_from_every_retained_na_candidate = false;
        }
    }
    if (water_index >= point.composition.size()) {
        witness.distinct_from_every_retained_na_candidate = false;
        return witness;
    }
    const double trial_water = point.composition[water_index];
    witness.trial_minus_min_na_water_fraction =
        trial_water - min_na_candidate_water_fraction;
    witness.water_role_roundoff_guard = 256.0 * stability_eps *
        (1.0 + std::abs(trial_water) +
         std::abs(min_na_candidate_water_fraction));
    witness.water_role_candidate_admissible =
        std::isfinite(witness.trial_minus_min_na_water_fraction) &&
        std::isfinite(witness.water_role_roundoff_guard) &&
        witness.trial_minus_min_na_water_fraction >
            witness.water_role_roundoff_guard;
    return witness;
}

} // namespace detail

[[nodiscard]] inline Sw92PhaseAssignedNoWResult solve_sw92_phase_assigned_no_w(
    double pressure_pa, double temperature_k, std::span<const double> feed,
    const thermodynamics::Sw92Phase<double>& model,
    double nacl_molality_mol_per_kg_water,
    Sw92PhaseAssignedNoWOptions options = {},
    std::span<const std::vector<double>> nonaqueous_initial_starts = {},
    std::span<const std::vector<double>> nonaqueous_final_starts = {},
    std::span<const std::vector<double>> aqueous_extra_starts = {}) {
    detail::sw92_phase_assigned_no_w_check_options(options);
    if (feed.size() != model.size()) {
        throw std::invalid_argument(
            "SW92 Profile-C no-W: feed does not match ordered model snapshot");
    }

    Sw92PhaseAssignedNoWResult result;
    result.options = options;
    result.pressure_pa = pressure_pa;
    result.temperature_k = temperature_k;
    result.input_feed_sum = detail::stability_check_composition(feed);
    result.nacl_molality_mol_per_kg_water = nacl_molality_mol_per_kg_water;
    result.water_index = model.parameters().water_index();
    const auto& parameters = model.parameters();
    result.dataset_id = parameters.dataset_id();
    result.revision = parameters.revision();
    for (const auto& component : parameters.components().items()) {
        result.component_ids.push_back(component.id);
    }

    Sw92FamilyVleEvaluator evaluator(
        model, nacl_molality_mol_per_kg_water,
        thermodynamics::SwPhaseFamily::nonaqueous,
        options.nonaqueous_root_options);
    result.hydrocarbon_flash = solve_sw92_pt_family_vle(
        pressure_pa, temperature_k, feed, evaluator,
        options.nonaqueous_flash, nonaqueous_initial_starts,
        nonaqueous_final_starts);
    result.feed = result.hydrocarbon_flash.solution.initial_stability.feed;

    const auto& solution = result.hydrocarbon_flash.solution;
    if (solution.status == PtSplitStatus::phase_set_unstable) {
        result.status = Sw92PhaseAssignedNoWStatus::
            higher_h_multiplicity_or_wrong_candidate;
        result.diagnostic =
            "fixed-NA maximum-two-phase state has a further same-family instability; no-W topology is not closed";
        return result;
    }
    if (solution.status == PtSplitStatus::indeterminate) {
        bool saw_disappearance = false;
        for (const auto& attempt : solution.attempts) {
            saw_disappearance = saw_disappearance ||
                attempt.status == PtSplitAttemptStatus::phase_disappearance;
        }
        result.status = saw_disappearance
            ? Sw92PhaseAssignedNoWStatus::phase_disappearance_unresolved
            : Sw92PhaseAssignedNoWStatus::indeterminate;
        result.diagnostic = solution.diagnostic.empty()
            ? "fixed-NA flash is indeterminate"
            : solution.diagnostic;
        return result;
    }

    if (solution.status == PtSplitStatus::single_phase_no_instability_found) {
        if (!solution.initial_stability.reference || result.feed.empty()) {
            result.status = Sw92PhaseAssignedNoWStatus::indeterminate;
            result.diagnostic = "fixed-NA single-phase result lacks reference activity";
            return result;
        }
        result.retained_na_candidate_compositions.push_back(result.feed);
        result.retained_na_candidate_fractions.push_back(1.0);
        result.common_log_activity.assign(result.feed.size(), 0.0);
        for (std::size_t i = 0; i < result.feed.size(); ++i) {
            if (result.feed[i] > 0.0) {
                result.common_log_activity[i] = std::log(result.feed[i]) +
                    solution.initial_stability.reference->ln_phi[i];
            }
        }
    } else if (solution.status == PtSplitStatus::two_phase_no_instability_found) {
        const auto* candidate = solution.candidate();
        if (candidate == nullptr || !solution.final_stability ||
            solution.final_stability->status != StabilityStatus::no_instability_found) {
            result.status = Sw92PhaseAssignedNoWStatus::indeterminate;
            result.diagnostic = "fixed-NA accepted pair lost its final-stability evidence";
            return result;
        }
        result.retained_na_candidate_compositions.push_back(candidate->fractions.liquid);
        result.retained_na_candidate_compositions.push_back(candidate->fractions.vapor);
        result.retained_na_candidate_fractions.push_back(
            1.0 - candidate->fractions.vapor_fraction);
        result.retained_na_candidate_fractions.push_back(
            candidate->fractions.vapor_fraction);
        result.common_log_activity = candidate->common_log_activity;
    } else {
        result.status = Sw92PhaseAssignedNoWStatus::indeterminate;
        result.diagnostic = "unexpected fixed-NA flash state";
        return result;
    }

    if (result.water_index >= result.feed.size() ||
        result.retained_na_candidate_compositions.empty()) {
        result.status = Sw92PhaseAssignedNoWStatus::indeterminate;
        result.diagnostic = "water index/candidate inventory is invalid";
        return result;
    }
    result.min_retained_na_candidate_water_fraction =
        std::numeric_limits<double>::infinity();
    result.max_retained_na_candidate_water_fraction = 0.0;
    for (const auto& candidate : result.retained_na_candidate_compositions) {
        if (candidate.size() != result.feed.size()) {
            result.status = Sw92PhaseAssignedNoWStatus::indeterminate;
            result.diagnostic = "retained fixed-NA candidate dimension mismatch";
            return result;
        }
        const double water = candidate[result.water_index];
        result.min_retained_na_candidate_water_fraction = std::min(
            result.min_retained_na_candidate_water_fraction, water);
        result.max_retained_na_candidate_water_fraction = std::max(
            result.max_retained_na_candidate_water_fraction, water);
    }
    if (!std::isfinite(result.min_retained_na_candidate_water_fraction) ||
        !std::isfinite(result.max_retained_na_candidate_water_fraction)) {
        result.status = Sw92PhaseAssignedNoWStatus::indeterminate;
        result.diagnostic = "retained fixed-NA water range is nonrepresentable";
        return result;
    }

    if (result.feed[result.water_index] == 0.0) {
        result.status = result.retained_na_candidate_count() == 1U
            ? Sw92PhaseAssignedNoWStatus::no_w_single_h_locally_closed
            : Sw92PhaseAssignedNoWStatus::no_w_two_h_locally_closed;
        result.diagnostic =
            "no AQ phase can appear because total water inventory is zero; retained NA morphology remains unresolved";
        return result;
    }

    result.targeted_water_start = detail::sw92_phase_assigned_targeted_water_start(
        result.feed, result.water_index,
        result.max_retained_na_candidate_water_fraction);
    if (!result.targeted_water_start) {
        result.status = Sw92PhaseAssignedNoWStatus::single_phase_role_unresolved;
        result.diagnostic =
            "a distinct water-richer AQ trial cannot be represented on active support; physical role/topology remains unresolved";
        return result;
    }

    const std::size_t mandatory = 1U + result.retained_na_candidate_count();
    detail::sw92_phase_assigned_no_w_preflight_starts(
        result.feed, options.aqueous_appearance, mandatory,
        aqueous_extra_starts);
    std::vector<std::vector<double>> starts;
    starts.reserve(mandatory + aqueous_extra_starts.size());
    starts.push_back(*result.targeted_water_start); // sole candidate-generating W start
    for (const auto& candidate : result.retained_na_candidate_compositions) {
        starts.push_back(candidate); // AQ-at-NA diagnostic only
    }
    for (const auto& extra : aqueous_extra_starts) {
        starts.push_back(extra); // diagnostic in v1; not candidate-generating
    }

    Sw92FamilyStabilityEvaluator aqueous(
        model, nacl_molality_mol_per_kg_water,
        thermodynamics::SwPhaseFamily::aqueous,
        options.aqueous_root_options);
    result.aqueous_appearance_search = test_pt_stability_against(
        pressure_pa, temperature_k, result.feed, result.common_log_activity,
        aqueous, options.aqueous_appearance, starts);
    const auto& search = *result.aqueous_appearance_search;
    for (std::size_t trial_index = 0; trial_index < search.trials.size(); ++trial_index) {
        const auto& trial = search.trials[trial_index];
        if (trial.status != StabilityTrialStatus::negative_tpd || !trial.point ||
            !detail::stability_negative(*trial.point, search.options)) {
            continue;
        }
        // automatic_starts is false by default and the targeted water start is
        // inserted first. All H-origin and caller starts remain diagnostic-only
        // so AQ model descent from an H-like seed cannot create a physical W.
        if (trial_index != 0U) { continue; }
        result.aqueous_witnesses.push_back(
            detail::sw92_phase_assigned_classify_water_witness(
                trial_index, *trial.point,
                result.retained_na_candidate_compositions, result.feed,
                result.water_index,
                result.min_retained_na_candidate_water_fraction,
                options.log_composition_separation));
    }

    for (std::size_t i = 0; i < result.aqueous_witnesses.size(); ++i) {
        if (result.aqueous_witnesses[i].usable_w_seed()) {
            result.selected_water_witness_index = i;
            result.status = Sw92PhaseAssignedNoWStatus::aqueous_phase_witness_found;
            result.diagnostic =
                "targeted AQ finite search found a robust distinct W/H-contrast phase-addition witness; joint Profile-C equations must decide the final W-vs-H role set";
            return result;
        }
    }

    // Negative diagnostics from retained NA candidates do not make the targeted
    // W search indeterminate; they are exactly the forbidden model-reassignment
    // evidence. The targeted start itself, however, must have a resolved result.
    if (!search.trials.empty()) {
        const auto targeted_status = search.trials.front().status;
        if (targeted_status == StabilityTrialStatus::stationary ||
            targeted_status == StabilityTrialStatus::negative_tpd) {
            result.status = result.retained_na_candidate_count() == 1U
                ? Sw92PhaseAssignedNoWStatus::no_w_single_h_locally_closed
                : Sw92PhaseAssignedNoWStatus::no_w_two_h_locally_closed;
            result.diagnostic =
                "fixed-NA multiplicity and targeted AQ appearance are locally closed; H-origin AQ negatives remain diagnostic-only and morphology unresolved";
            return result;
        }
    }

    result.status = Sw92PhaseAssignedNoWStatus::indeterminate;
    result.diagnostic =
        "targeted AQ water-phase search is numerically/property indeterminate";
    return result;
}

} // namespace mpmc::flash

#endif // MPMC_FLASH_SW92_PHASE_ASSIGNED_NO_W_HPP
