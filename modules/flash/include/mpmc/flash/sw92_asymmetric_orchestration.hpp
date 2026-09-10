#ifndef MPMC_FLASH_SW92_ASYMMETRIC_ORCHESTRATION_HPP
#define MPMC_FLASH_SW92_ASYMMETRIC_ORCHESTRATION_HPP

#include <mpmc/flash/sw92_asymmetric_fixed_pair.hpp>

#include <algorithm>
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

/// Gate 3B.2 numerical orchestration identity. This selects a family-aware
/// single/pair candidate but deliberately stops before Gate 3B.3 final
/// two-family phase-set stability and accepted-phase-set publication.
inline constexpr std::string_view sw92_xu_asymmetric_orchestration_convention =
    "SW92-equilibrium/xu-asymmetric-gibbs/max2-witness-orchestration-gibbs-selection/v1";

enum class Sw92AsymmetricPairPlanKind {
    reference_with_witness,
    witness_same_family_alternative
};

struct Sw92AsymmetricPairPlanEntry {
    thermodynamics::SwPhaseFamily witness_family{
        thermodynamics::SwPhaseFamily::aqueous};
    std::size_t witness_trial_index{};
    double witness_tpd{};
    double witness_roundoff_guard{};
    std::vector<double> witness_composition;
    Sw92AsymmetricPairPlanKind kind{
        Sw92AsymmetricPairPlanKind::reference_with_witness};
    Sw92AsymmetricFamilyPair family_pair;
    bool protected_family_representative{};
    std::optional<std::vector<double>> initial_log_k;
    std::string diagnostic;
};

struct Sw92AsymmetricPairAttemptRecord {
    std::size_t plan_index{};
    std::optional<Sw92AsymmetricFixedPairResult> fixed_pair;
    std::string diagnostic;
};

struct Sw92AsymmetricCandidateClass {
    std::size_t representative_attempt{};
    std::vector<std::size_t> equivalent_attempts;
    double reduced_gibbs{std::numeric_limits<double>::quiet_NaN()};
    double gibbs_roundoff_guard{std::numeric_limits<double>::quiet_NaN()};
};

struct Sw92AsymmetricSinglePhaseCandidate {
    thermodynamics::SwPhaseFamily family{
        thermodynamics::SwPhaseFamily::aqueous};
    std::vector<double> composition;
    StabilityPhase activity;
    double reduced_gibbs{std::numeric_limits<double>::quiet_NaN()};
};

struct Sw92AsymmetricPairSelectionOptions {
    Sw92AsymmetricStabilityOptions initial_stability;
    Sw92AsymmetricFixedPairOptions fixed_pair;
    std::size_t max_pair_attempts{32};
};

enum class Sw92AsymmetricPairSelectionStatus {
    single_phase_candidate_no_instability_found,
    pair_candidate_selected_pending_final_stability,
    initial_stability_indeterminate,
    pair_search_indeterminate,
    attempt_limit_reached,
    candidate_gibbs_tie
};

/// Gate 3B.2 result. There is intentionally no accepted_phase_set(): a selected
/// pair still requires Gate 3B.3 pair-vs-feed Gibbs and final AQ/NA stability.
struct Sw92AsymmetricPairSelectionResult {
    static constexpr std::size_t maximum_phase_count = 2;
    static constexpr bool global_stability_proven = false;

    Sw92AsymmetricPairSelectionStatus status{
        Sw92AsymmetricPairSelectionStatus::initial_stability_indeterminate};
    Sw92AsymmetricPairSelectionOptions options;

    double pressure_pa{};
    double temperature_k{};
    std::vector<double> feed;
    double nacl_molality_mol_per_kg_water{};

    std::string dataset_id;
    std::string revision;
    std::vector<std::string> component_ids;
    std::string model_profile{thermodynamics::sw92_corrected_profile};
    std::string phase_convention{thermodynamics::sw92_pt_convention};
    std::string equilibrium_profile{sw92_xu_asymmetric_gibbs_profile};
    std::string orchestration_convention{
        sw92_xu_asymmetric_orchestration_convention};

    Sw92AsymmetricStabilityResult initial_stability;
    std::optional<Sw92AsymmetricSinglePhaseCandidate> single_phase_candidate;

    std::vector<Sw92AsymmetricPairPlanEntry> plan;
    std::vector<Sw92AsymmetricPairAttemptRecord> attempts;
    bool attempt_limit_reached{};
    std::size_t total_pair_property_evaluations{};

    double candidate_fraction_equivalence_tolerance{};
    double candidate_log_composition_equivalence_tolerance{};
    std::vector<Sw92AsymmetricCandidateClass> candidate_classes;
    std::optional<std::size_t> selected_candidate_class;
    std::string diagnostic;

    [[nodiscard]] const Sw92AsymmetricFixedPairState*
        selected_pair_candidate_pending_final_stability() const & noexcept {
        if (status != Sw92AsymmetricPairSelectionStatus::
                          pair_candidate_selected_pending_final_stability ||
            !selected_candidate_class ||
            *selected_candidate_class >= candidate_classes.size()) {
            return nullptr;
        }
        const auto attempt_index =
            candidate_classes[*selected_candidate_class].representative_attempt;
        if (attempt_index >= attempts.size() || !attempts[attempt_index].fixed_pair ||
            !attempts[attempt_index].fixed_pair->candidate_admissible()) {
            return nullptr;
        }
        return &*attempts[attempt_index].fixed_pair->point;
    }
    const Sw92AsymmetricFixedPairState*
        selected_pair_candidate_pending_final_stability() const && = delete;
};

namespace detail {

inline std::optional<std::vector<double>> sw92_asymmetric_witness_seed(
    std::span<const double> feed, std::span<const double> witness) {
    if (feed.size() != witness.size() || feed.empty()) {
        return std::nullopt;
    }
    try {
        (void)stability_check_composition(feed);
        (void)stability_check_composition(witness);
    } catch (const std::exception&) {
        return std::nullopt;
    }

    double incipient_fraction = 0.1;
    for (std::size_t i = 0; i < feed.size(); ++i) {
        if (feed[i] == 0.0) {
            if (witness[i] != 0.0) { return std::nullopt; }
            continue;
        }
        if (!(witness[i] > 0.0)) { return std::nullopt; }
        incipient_fraction = std::min(
            incipient_fraction, 0.5 * (feed[i] / witness[i]));
    }
    if (!(incipient_fraction > 0.0 && incipient_fraction < 1.0)) {
        return std::nullopt;
    }

    std::vector<double> phase0(feed.size());
    std::vector<double> log_k(feed.size());
    for (std::size_t i = 0; i < feed.size(); ++i) {
        if (feed[i] == 0.0) { continue; }
        phase0[i] = std::fma(-incipient_fraction, witness[i], feed[i]) /
                    (1.0 - incipient_fraction);
        if (!(phase0[i] > 0.0) || !std::isfinite(phase0[i])) {
            return std::nullopt;
        }
    }
    const double sum = stability_sum(phase0);
    if (!std::isfinite(sum) || std::abs(sum - 1.0) > 64.0 * stability_eps) {
        return std::nullopt;
    }
    for (std::size_t i = 0; i < feed.size(); ++i) {
        if (feed[i] == 0.0) { continue; }
        log_k[i] = std::log(witness[i]) - std::log(phase0[i]);
        if (!std::isfinite(log_k[i])) { return std::nullopt; }
    }
    return log_k;
}

inline bool sw92_asymmetric_witness_less(
    const Sw92AsymmetricNegativeWitness* first,
    const Sw92AsymmetricNegativeWitness* second) {
    if (first->point.value != second->point.value) {
        return first->point.value < second->point.value;
    }
    if (first->point.roundoff_guard != second->point.roundoff_guard) {
        return first->point.roundoff_guard < second->point.roundoff_guard;
    }
    if (first->point.composition != second->point.composition) {
        return std::lexicographical_compare(
            first->point.composition.begin(), first->point.composition.end(),
            second->point.composition.begin(), second->point.composition.end());
    }
    if (first->trial_index != second->trial_index) {
        return first->trial_index < second->trial_index;
    }
    return static_cast<int>(first->family) < static_cast<int>(second->family);
}

inline Sw92AsymmetricPairPlanEntry sw92_make_pair_plan_entry(
    const Sw92AsymmetricNegativeWitness& witness,
    thermodynamics::SwPhaseFamily reference_family,
    Sw92AsymmetricPairPlanKind kind,
    std::span<const double> feed,
    bool protected_family_representative) {
    Sw92AsymmetricPairPlanEntry entry;
    entry.witness_family = witness.family;
    entry.witness_trial_index = witness.trial_index;
    entry.witness_tpd = witness.point.value;
    entry.witness_roundoff_guard = witness.point.roundoff_guard;
    entry.witness_composition = witness.point.composition;
    entry.kind = kind;
    entry.protected_family_representative = protected_family_representative;
    if (kind == Sw92AsymmetricPairPlanKind::reference_with_witness) {
        entry.family_pair = {reference_family, witness.family};
    } else {
        entry.family_pair = {witness.family, witness.family};
    }
    entry.initial_log_k = sw92_asymmetric_witness_seed(
        feed, witness.point.composition);
    if (!entry.initial_log_k) {
        entry.diagnostic =
            "family-tagged negative witness cannot form a finite positive material-balanced seed";
    }
    return entry;
}

inline double sw92_candidate_fraction_tolerance(
    const Sw92AsymmetricFixedPairOptions& options) noexcept {
    return std::max(512.0 * stability_eps,
                    8.0 * std::max(options.mass_absolute_tolerance,
                                   options.mass_relative_tolerance));
}

inline double sw92_candidate_log_composition_tolerance(
    const Sw92AsymmetricFixedPairOptions& options) noexcept {
    return std::max(512.0 * stability_eps, options.log_k_separation);
}

inline bool sw92_candidate_phase_equivalent(
    const Sw92AsymmetricFixedPairPhase& first,
    const Sw92AsymmetricFixedPairPhase& second,
    std::span<const double> feed,
    double fraction_tolerance,
    double log_composition_tolerance) {
    if (first.family != second.family ||
        std::abs(first.mole_phase_fraction - second.mole_phase_fraction) >
            fraction_tolerance ||
        first.composition.size() != feed.size() ||
        second.composition.size() != feed.size()) {
        return false;
    }
    for (std::size_t i = 0; i < feed.size(); ++i) {
        if (feed[i] == 0.0) { continue; }
        if (!(first.composition[i] > 0.0) || !(second.composition[i] > 0.0)) {
            return false;
        }
        const double difference =
            std::abs(std::log(first.composition[i]) -
                     std::log(second.composition[i]));
        if (!std::isfinite(difference) ||
            difference > log_composition_tolerance) {
            return false;
        }
    }
    return true;
}

inline bool sw92_candidate_state_equivalent(
    const Sw92AsymmetricFixedPairState& first,
    const Sw92AsymmetricFixedPairState& second,
    std::span<const double> feed,
    double fraction_tolerance,
    double log_composition_tolerance) {
    const bool direct = sw92_candidate_phase_equivalent(
                            first.phase0, second.phase0, feed,
                            fraction_tolerance, log_composition_tolerance) &&
                        sw92_candidate_phase_equivalent(
                            first.phase1, second.phase1, feed,
                            fraction_tolerance, log_composition_tolerance);
    if (direct) { return true; }
    return sw92_candidate_phase_equivalent(
               first.phase0, second.phase1, feed,
               fraction_tolerance, log_composition_tolerance) &&
           sw92_candidate_phase_equivalent(
               first.phase1, second.phase0, feed,
               fraction_tolerance, log_composition_tolerance);
}

inline bool sw92_candidate_representative_better(
    const Sw92AsymmetricFixedPairResult& candidate,
    const Sw92AsymmetricFixedPairResult& current) noexcept {
    const auto& a = *candidate.point;
    const auto& b = *current.point;
    if (a.chemical_potential_norm != b.chemical_potential_norm) {
        return a.chemical_potential_norm < b.chemical_potential_norm;
    }
    if (a.mass_relative != b.mass_relative) {
        return a.mass_relative < b.mass_relative;
    }
    if (a.mass_absolute != b.mass_absolute) {
        return a.mass_absolute < b.mass_absolute;
    }
    return a.gibbs_roundoff_guard < b.gibbs_roundoff_guard;
}

inline std::vector<Sw92AsymmetricCandidateClass> sw92_build_candidate_classes(
    const std::vector<Sw92AsymmetricPairAttemptRecord>& attempts,
    std::span<const double> feed,
    const Sw92AsymmetricFixedPairOptions& options) {
    const double fraction_tolerance = sw92_candidate_fraction_tolerance(options);
    const double log_composition_tolerance =
        sw92_candidate_log_composition_tolerance(options);
    std::vector<Sw92AsymmetricCandidateClass> classes;

    for (std::size_t attempt_index = 0; attempt_index < attempts.size(); ++attempt_index) {
        const auto& record = attempts[attempt_index];
        if (!record.fixed_pair || !record.fixed_pair->candidate_admissible()) {
            continue;
        }
        const auto& candidate = *record.fixed_pair;
        const auto& state = *candidate.point;
        std::optional<std::size_t> matched;
        for (std::size_t class_index = 0; class_index < classes.size(); ++class_index) {
            const auto representative_index = classes[class_index].representative_attempt;
            const auto& representative = *attempts[representative_index].fixed_pair;
            if (sw92_candidate_state_equivalent(
                    state, *representative.point, feed,
                    fraction_tolerance, log_composition_tolerance)) {
                matched = class_index;
                break;
            }
        }
        if (!matched) {
            Sw92AsymmetricCandidateClass candidate_class;
            candidate_class.representative_attempt = attempt_index;
            candidate_class.equivalent_attempts.push_back(attempt_index);
            candidate_class.reduced_gibbs = state.reduced_gibbs;
            candidate_class.gibbs_roundoff_guard = state.gibbs_roundoff_guard;
            classes.push_back(std::move(candidate_class));
            continue;
        }
        auto& candidate_class = classes[*matched];
        candidate_class.equivalent_attempts.push_back(attempt_index);
        const auto current_index = candidate_class.representative_attempt;
        const auto& current = *attempts[current_index].fixed_pair;
        if (sw92_candidate_representative_better(candidate, current)) {
            candidate_class.representative_attempt = attempt_index;
            candidate_class.reduced_gibbs = state.reduced_gibbs;
            candidate_class.gibbs_roundoff_guard = state.gibbs_roundoff_guard;
        }
    }
    return classes;
}

struct Sw92CandidateChoice {
    std::optional<std::size_t> selected_class;
    bool distinct_gibbs_tie{};
};

inline Sw92CandidateChoice sw92_choose_candidate_class(
    const std::vector<Sw92AsymmetricCandidateClass>& classes) {
    Sw92CandidateChoice choice;
    if (classes.empty()) { return choice; }

    std::size_t best = 0;
    for (std::size_t i = 1; i < classes.size(); ++i) {
        if (classes[i].reduced_gibbs < classes[best].reduced_gibbs) {
            best = i;
        }
    }
    for (std::size_t i = 0; i < classes.size(); ++i) {
        if (i == best) { continue; }
        const double guard = classes[best].gibbs_roundoff_guard +
                             classes[i].gibbs_roundoff_guard;
        if (std::abs(classes[i].reduced_gibbs - classes[best].reduced_gibbs) <=
            guard) {
            choice.distinct_gibbs_tie = true;
            return choice;
        }
    }
    choice.selected_class = best;
    return choice;
}

} // namespace detail

/// Build a deterministic family-neutral attempt plan from the complete Gate-3A
/// negative-witness set. For each witness rank, the primary entry from every
/// family is placed before same-family alternatives, preventing one family from
/// consuming all early slots merely because Gate 3A stores AQ diagnostics first.
[[nodiscard]] inline std::vector<Sw92AsymmetricPairPlanEntry>
    build_sw92_asymmetric_pair_plan(
        const Sw92AsymmetricStabilityResult& initial_stability) {
    if (initial_stability.status != StabilityStatus::unstable ||
        initial_stability.feed_reference_status !=
            Sw92AsymmetricFeedReferenceStatus::selected ||
        !initial_stability.reference_family) {
        throw std::invalid_argument(
            "SW92 asymmetric orchestration: unstable Gate-3A result with resolved reference family required");
    }

    using Witness = Sw92AsymmetricNegativeWitness;
    std::vector<const Witness*> aqueous;
    std::vector<const Witness*> nonaqueous;
    for (const auto& witness : initial_stability.negative_witnesses) {
        switch (witness.family) {
        case thermodynamics::SwPhaseFamily::aqueous:
            aqueous.push_back(&witness);
            break;
        case thermodynamics::SwPhaseFamily::nonaqueous:
            nonaqueous.push_back(&witness);
            break;
        }
    }
    std::sort(aqueous.begin(), aqueous.end(), detail::sw92_asymmetric_witness_less);
    std::sort(nonaqueous.begin(), nonaqueous.end(), detail::sw92_asymmetric_witness_less);

    std::vector<Sw92AsymmetricPairPlanEntry> plan;
    const std::size_t maximum_rank = std::max(aqueous.size(), nonaqueous.size());
    const std::size_t extra_alternatives =
        (*initial_stability.reference_family == thermodynamics::SwPhaseFamily::aqueous
             ? nonaqueous.size() : aqueous.size());
    if (initial_stability.negative_witnesses.size() >
        std::numeric_limits<std::size_t>::max() - extra_alternatives) {
        throw std::length_error("SW92 asymmetric orchestration: pair-plan size overflow");
    }
    plan.reserve(initial_stability.negative_witnesses.size() + extra_alternatives);

    for (std::size_t rank = 0; rank < maximum_rank; ++rank) {
        std::vector<const Witness*> round;
        if (rank < aqueous.size()) { round.push_back(aqueous[rank]); }
        if (rank < nonaqueous.size()) { round.push_back(nonaqueous[rank]); }
        std::sort(round.begin(), round.end(), detail::sw92_asymmetric_witness_less);

        for (const Witness* witness : round) {
            plan.push_back(detail::sw92_make_pair_plan_entry(
                *witness, *initial_stability.reference_family,
                Sw92AsymmetricPairPlanKind::reference_with_witness,
                initial_stability.feed, rank == 0));
        }
        for (const Witness* witness : round) {
            if (witness->family == *initial_stability.reference_family) { continue; }
            plan.push_back(detail::sw92_make_pair_plan_entry(
                *witness, *initial_stability.reference_family,
                Sw92AsymmetricPairPlanKind::witness_same_family_alternative,
                initial_stability.feed, false));
        }
    }
    return plan;
}

/// Gate 3B.2: initial asymmetric stability -> family-neutral witness plan ->
/// fixed-pair attempts -> slot-swap-aware deduplication -> guarded Gibbs choice.
/// A selected pair is explicitly pending Gate 3B.3 final stability.
[[nodiscard]] inline Sw92AsymmetricPairSelectionResult
    orchestrate_sw92_asymmetric_pair_candidates(
        double pressure_pa, double temperature_k, std::span<const double> feed,
        const thermodynamics::Sw92Phase<double>& model,
        double nacl_molality_mol_per_kg_water,
        Sw92AsymmetricPairSelectionOptions options = {},
        Sw92AsymmetricStabilityStarts initial_starts = {}) {
    if (options.max_pair_attempts == 0) {
        throw std::invalid_argument(
            "SW92 asymmetric orchestration: max_pair_attempts must be positive");
    }
    detail::sw92_fixed_pair_check_options(options.fixed_pair);
    if (feed.empty() || feed.size() > options.fixed_pair.rr.max_components) {
        throw std::length_error(
            "SW92 asymmetric orchestration: fixed-pair component quota exceeded or empty feed");
    }

    Sw92AsymmetricPairSelectionResult result;
    result.options = options;
    result.pressure_pa = pressure_pa;
    result.temperature_k = temperature_k;
    result.nacl_molality_mol_per_kg_water = nacl_molality_mol_per_kg_water;
    result.candidate_fraction_equivalence_tolerance =
        detail::sw92_candidate_fraction_tolerance(options.fixed_pair);
    result.candidate_log_composition_equivalence_tolerance =
        detail::sw92_candidate_log_composition_tolerance(options.fixed_pair);

    const auto& parameters = model.parameters();
    result.dataset_id = parameters.dataset_id();
    result.revision = parameters.revision();
    for (const auto& component : parameters.components().items()) {
        result.component_ids.push_back(component.id);
    }

    result.initial_stability = test_sw92_pt_asymmetric_stability(
        pressure_pa, temperature_k, feed, model,
        nacl_molality_mol_per_kg_water, options.initial_stability,
        initial_starts);
    result.feed = result.initial_stability.feed;

    if (result.initial_stability.status == StabilityStatus::no_instability_found) {
        if (result.initial_stability.feed_reference_status !=
                Sw92AsymmetricFeedReferenceStatus::selected ||
            !result.initial_stability.reference_family) {
            result.status =
                Sw92AsymmetricPairSelectionStatus::initial_stability_indeterminate;
            result.diagnostic =
                "Gate 3A reported no instability without a resolved lower-feed family";
            return result;
        }
        const auto family = *result.initial_stability.reference_family;
        const auto& leg = family == thermodynamics::SwPhaseFamily::aqueous
            ? result.initial_stability.aqueous
            : result.initial_stability.nonaqueous;
        if (!leg.feed_reference || !std::isfinite(leg.feed_reduced_gibbs)) {
            result.status =
                Sw92AsymmetricPairSelectionStatus::initial_stability_indeterminate;
            result.diagnostic =
                "resolved lower-feed family lacks a representable retained feed reference";
            return result;
        }
        result.single_phase_candidate = Sw92AsymmetricSinglePhaseCandidate{
            family, result.feed, *leg.feed_reference, leg.feed_reduced_gibbs};
        result.status = Sw92AsymmetricPairSelectionStatus::
            single_phase_candidate_no_instability_found;
        result.diagnostic =
            "Gate 3A finite search found no instability; retained family-aware single-phase candidate is not a global proof";
        return result;
    }

    if (result.initial_stability.status != StabilityStatus::unstable ||
        result.initial_stability.feed_reference_status !=
            Sw92AsymmetricFeedReferenceStatus::selected ||
        !result.initial_stability.reference_family) {
        result.status =
            Sw92AsymmetricPairSelectionStatus::initial_stability_indeterminate;
        result.diagnostic =
            "Gate 3A initial asymmetric stability is indeterminate; no pair orchestration started";
        return result;
    }

    result.plan = build_sw92_asymmetric_pair_plan(result.initial_stability);
    if (result.plan.empty()) {
        result.status = Sw92AsymmetricPairSelectionStatus::pair_search_indeterminate;
        result.diagnostic =
            "Gate 3A reported instability but retained no robust family-tagged witness for pair planning";
        return result;
    }

    const std::size_t attempt_count =
        std::min(options.max_pair_attempts, result.plan.size());
    result.attempts.reserve(attempt_count);
    for (std::size_t plan_index = 0; plan_index < attempt_count; ++plan_index) {
        const auto& entry = result.plan[plan_index];
        Sw92AsymmetricPairAttemptRecord record;
        record.plan_index = plan_index;
        if (!entry.initial_log_k) {
            record.diagnostic = entry.diagnostic;
            result.attempts.push_back(std::move(record));
            continue;
        }
        record.fixed_pair = iterate_sw92_asymmetric_fixed_pair(
            pressure_pa, temperature_k, result.feed, *entry.initial_log_k,
            entry.family_pair, model, nacl_molality_mol_per_kg_water,
            options.fixed_pair);
        if (record.fixed_pair->evaluations >
            std::numeric_limits<std::size_t>::max() -
                result.total_pair_property_evaluations) {
            throw std::length_error(
                "SW92 asymmetric orchestration: property-evaluation count overflow");
        }
        result.total_pair_property_evaluations += record.fixed_pair->evaluations;
        record.diagnostic = record.fixed_pair->diagnostic;
        result.attempts.push_back(std::move(record));
    }

    result.candidate_classes = detail::sw92_build_candidate_classes(
        result.attempts, result.feed, options.fixed_pair);
    result.attempt_limit_reached = result.attempts.size() < result.plan.size();
    if (result.attempt_limit_reached) {
        result.status = Sw92AsymmetricPairSelectionStatus::attempt_limit_reached;
        result.diagnostic =
            "pair-attempt budget truncated the family-neutral plan; candidate selection intentionally withheld";
        return result;
    }

    if (result.candidate_classes.empty()) {
        result.status = Sw92AsymmetricPairSelectionStatus::pair_search_indeterminate;
        result.diagnostic =
            "all planned fixed-family attempts completed, but none produced an admissible lower-envelope pair candidate";
        return result;
    }

    const auto choice = detail::sw92_choose_candidate_class(result.candidate_classes);
    if (choice.distinct_gibbs_tie) {
        result.status = Sw92AsymmetricPairSelectionStatus::candidate_gibbs_tie;
        result.diagnostic =
            "distinct deduplicated pair candidates are tied within combined Gibbs arithmetic guards";
        return result;
    }
    if (!choice.selected_class) {
        result.status = Sw92AsymmetricPairSelectionStatus::pair_search_indeterminate;
        result.diagnostic =
            "candidate classification completed without a selectable Gibbs minimum";
        return result;
    }

    result.selected_candidate_class = choice.selected_class;
    result.status = Sw92AsymmetricPairSelectionStatus::
        pair_candidate_selected_pending_final_stability;
    result.diagnostic =
        "lowest resolved deduplicated fixed-pair Gibbs candidate selected; Gate 3B.3 final two-family phase-set stability is still required";
    return result;
}

} // namespace mpmc::flash

#endif // MPMC_FLASH_SW92_ASYMMETRIC_ORCHESTRATION_HPP
