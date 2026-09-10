#ifndef MPMC_FLASH_SW92_ASYMMETRIC_MAX2_HPP
#define MPMC_FLASH_SW92_ASYMMETRIC_MAX2_HPP

#include <mpmc/flash/sw92_asymmetric_orchestration.hpp>

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

/// Gate 3B.3 maximum-two-phase algorithm identity. `max2` is an entry-point
/// capability, not a proof that the thermodynamic model can never admit >2 phases.
inline constexpr std::string_view sw92_xu_asymmetric_max2_convention =
    "SW92-equilibrium/xu-asymmetric-gibbs/max2-logK-SSI-RR/common-tangent/v1";

/// Rejected/indeterminate final-review paths retain the same phase data, so the
/// type deliberately says candidate. Publication is controlled only by
/// Sw92AsymmetricMax2Result::accepted_phase_set().
struct Sw92AsymmetricCandidatePhase {
    thermodynamics::SwPhaseFamily family{
        thermodynamics::SwPhaseFamily::aqueous};
    double mole_phase_fraction{};
    std::vector<double> composition;
    StabilityPhase activity;
    std::optional<double> compressibility_factor;
};

struct Sw92AsymmetricCandidatePhaseSet {
    std::vector<Sw92AsymmetricCandidatePhase> phases;
};

struct Sw92AsymmetricMax2Options {
    Sw92AsymmetricPairSelectionOptions selection;
    Sw92AsymmetricStabilityOptions final_stability;
};

enum class Sw92AsymmetricMax2Status {
    single_phase_no_instability_found,
    two_phase_no_instability_found,
    phase_set_unstable,
    pair_gibbs_above_feed,
    indeterminate
};

struct Sw92AsymmetricMax2Result {
    static constexpr std::size_t maximum_phase_count = 2;
    static constexpr bool global_stability_proven = false;

    Sw92AsymmetricMax2Status status{Sw92AsymmetricMax2Status::indeterminate};
    Sw92AsymmetricMax2Options options;
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
    std::string max2_convention{sw92_xu_asymmetric_max2_convention};

    Sw92AsymmetricPairSelectionResult selection;
    std::optional<Sw92AsymmetricCandidatePhaseSet> candidate_phase_set;

    double lower_feed_reduced_gibbs{std::numeric_limits<double>::quiet_NaN()};
    double lower_feed_gibbs_roundoff_guard{
        std::numeric_limits<double>::quiet_NaN()};
    double selected_pair_reduced_gibbs{
        std::numeric_limits<double>::quiet_NaN()};
    double selected_pair_gibbs_roundoff_guard{
        std::numeric_limits<double>::quiet_NaN()};
    double pair_minus_lower_feed_reduced_gibbs{
        std::numeric_limits<double>::quiet_NaN()};
    double pair_feed_gibbs_combined_guard{
        std::numeric_limits<double>::quiet_NaN()};

    double common_reference_allowance{
        std::numeric_limits<double>::quiet_NaN()};
    double final_aqueous_base_tpd_tolerance{
        std::numeric_limits<double>::quiet_NaN()};
    double final_nonaqueous_base_tpd_tolerance{
        std::numeric_limits<double>::quiet_NaN()};
    double final_aqueous_effective_tpd_tolerance{
        std::numeric_limits<double>::quiet_NaN()};
    double final_nonaqueous_effective_tpd_tolerance{
        std::numeric_limits<double>::quiet_NaN()};
    std::optional<Sw92AsymmetricCommonTangentSearchResult> final_stability;

    std::string diagnostic;

    [[nodiscard]] const Sw92AsymmetricCandidatePhaseSet*
        accepted_phase_set() const & noexcept {
        const bool accepted_single =
            status == Sw92AsymmetricMax2Status::single_phase_no_instability_found;
        const bool accepted_pair =
            status == Sw92AsymmetricMax2Status::two_phase_no_instability_found;
        if ((!accepted_single && !accepted_pair) || !candidate_phase_set) {
            return nullptr;
        }
        const std::size_t expected = accepted_single ? 1U : 2U;
        if (candidate_phase_set->phases.size() != expected ||
            expected > maximum_phase_count) {
            return nullptr;
        }
        double fraction_sum = 0.0;
        double fraction_correction = 0.0;
        for (const auto& phase : candidate_phase_set->phases) {
            if (!std::isfinite(phase.mole_phase_fraction) ||
                !(phase.mole_phase_fraction > 0.0) ||
                phase.composition.size() != feed.size()) {
                return nullptr;
            }
            detail::stability_add(
                phase.mole_phase_fraction, fraction_sum, fraction_correction);
            switch (phase.family) {
            case thermodynamics::SwPhaseFamily::aqueous:
            case thermodynamics::SwPhaseFamily::nonaqueous:
                break;
            default:
                return nullptr;
            }
            double composition_sum = 0.0;
            double composition_correction = 0.0;
            for (std::size_t i = 0; i < phase.composition.size(); ++i) {
                const double value = phase.composition[i];
                if (!std::isfinite(value) || value < 0.0 || value > 1.0 ||
                    (feed[i] == 0.0 && value != 0.0) ||
                    (feed[i] > 0.0 && !(value > 0.0))) {
                    return nullptr;
                }
                detail::stability_add(
                    value, composition_sum, composition_correction);
            }
            if (std::abs(composition_sum - 1.0) >
                64.0 * detail::stability_eps) {
                return nullptr;
            }
        }
        if (!std::isfinite(fraction_sum) ||
            std::abs(fraction_sum - 1.0) > 256.0 * detail::stability_eps) {
            return nullptr;
        }
        return &*candidate_phase_set;
    }
    const Sw92AsymmetricCandidatePhaseSet* accepted_phase_set() const && = delete;

    [[nodiscard]] std::size_t accepted_phase_count() const noexcept {
        const auto* accepted = accepted_phase_set();
        return accepted == nullptr ? 0U : accepted->phases.size();
    }
};

namespace detail {

inline double sw92_max2_phase_gibbs_guard(
    std::span<const double> composition, const StabilityPhase& phase) {
    stability_check_phase(phase, composition.size());
    double magnitude = 1.0;
    double correction = 0.0;
    for (std::size_t i = 0; i < composition.size(); ++i) {
        if (composition[i] == 0.0) { continue; }
        if (!(composition[i] > 0.0)) {
            throw std::domain_error(
                "SW92 asymmetric max2: active composition must be positive");
        }
        stability_add(
            composition[i] *
                (std::abs(std::log(composition[i])) + std::abs(phase.ln_phi[i])),
            magnitude, correction);
    }
    const double guard = 256.0 * stability_eps * magnitude;
    if (!std::isfinite(guard)) {
        throw StabilityPropertyError(
            StabilityPropertyIssue::nonfinite_properties,
            "SW92 asymmetric max2: nonrepresentable Gibbs roundoff guard");
    }
    return guard;
}

inline bool sw92_max2_model_matches_selection(
    const Sw92AsymmetricPairSelectionResult& selection,
    const thermodynamics::Sw92Phase<double>& model) {
    if (selection.feed.size() != model.size()) { return false; }
    const auto& parameters = model.parameters();
    if (selection.dataset_id != parameters.dataset_id() ||
        selection.revision != parameters.revision() ||
        selection.component_ids.size() != model.size()) {
        return false;
    }
    std::size_t index = 0;
    for (const auto& component : parameters.components().items()) {
        if (selection.component_ids[index++] != component.id) { return false; }
    }
    return true;
}

inline bool sw92_max2_stability_matches_selection(
    const Sw92AsymmetricPairSelectionResult& selection) {
    const auto& initial = selection.initial_stability;
    return initial.pressure_pa == selection.pressure_pa &&
           initial.temperature_k == selection.temperature_k &&
           initial.feed == selection.feed &&
           initial.nacl_molality_mol_per_kg_water ==
               selection.nacl_molality_mol_per_kg_water &&
           initial.dataset_id == selection.dataset_id &&
           initial.revision == selection.revision &&
           initial.component_ids == selection.component_ids &&
           std::string_view{initial.model_profile} ==
               thermodynamics::sw92_corrected_profile &&
           std::string_view{initial.phase_convention} ==
               thermodynamics::sw92_pt_convention &&
           std::string_view{initial.equilibrium_profile} ==
               sw92_xu_asymmetric_gibbs_profile;
}

inline Sw92AsymmetricCandidatePhase sw92_max2_from_fixed_pair_phase(
    const Sw92AsymmetricFixedPairPhase& phase) {
    return {phase.family, phase.mole_phase_fraction, phase.composition,
            phase.activity, phase.compressibility_factor};
}

inline std::vector<std::vector<double>> sw92_max2_final_starts(
    std::span<const std::vector<double>> caller_starts,
    const Sw92AsymmetricFixedPairState& pair) {
    if (caller_starts.size() > std::numeric_limits<std::size_t>::max() - 2U) {
        throw std::length_error("SW92 asymmetric max2: final-start count overflow");
    }
    std::vector<std::vector<double>> starts;
    starts.reserve(caller_starts.size() + 2U);
    for (const auto& start : caller_starts) { starts.push_back(start); }
    starts.push_back(pair.phase0.composition);
    starts.push_back(pair.phase1.composition);
    return starts;
}

inline double sw92_max2_effective_tpd_tolerance(double base, double allowance) {
    // Preserve the generic StabilityOptions contract: tpd_tolerance==0 is legal.
    if (!std::isfinite(base) || base < 0.0 ||
        !std::isfinite(allowance) || allowance < 0.0 ||
        allowance > std::numeric_limits<double>::max() - base) {
        throw std::domain_error(
            "SW92 asymmetric max2: nonrepresentable final TPD tolerance/allowance");
    }
    const double effective = base + allowance;
    if (!std::isfinite(effective) || effective < 0.0) {
        throw std::domain_error(
            "SW92 asymmetric max2: invalid effective final TPD tolerance");
    }
    return effective;
}

inline const Sw92AsymmetricFixedPairState* sw92_max2_revalidate_pair_choice(
    const Sw92AsymmetricPairSelectionResult& selection) {
    if (selection.status !=
            Sw92AsymmetricPairSelectionStatus::
                pair_candidate_selected_pending_final_stability ||
        selection.attempt_limit_reached ||
        selection.attempts.size() != selection.plan.size()) {
        return nullptr;
    }
    const auto classes = sw92_build_candidate_classes(
        selection.attempts, selection.feed, selection.options.fixed_pair);
    const auto choice = sw92_choose_candidate_class(classes);
    if (choice.distinct_gibbs_tie || !choice.selected_class ||
        *choice.selected_class >= classes.size()) {
        return nullptr;
    }
    const std::size_t attempt_index =
        classes[*choice.selected_class].representative_attempt;
    if (attempt_index >= selection.attempts.size() ||
        !selection.attempts[attempt_index].fixed_pair ||
        !selection.attempts[attempt_index].fixed_pair->candidate_admissible()) {
        return nullptr;
    }
    return &*selection.attempts[attempt_index].fixed_pair->point;
}

inline bool sw92_max2_single_phase_evidence_consistent(
    const Sw92AsymmetricPairSelectionResult& selection) {
    if (selection.status != Sw92AsymmetricPairSelectionStatus::
                                single_phase_candidate_no_instability_found ||
        !selection.single_phase_candidate ||
        selection.initial_stability.status != StabilityStatus::no_instability_found ||
        selection.initial_stability.feed_reference_status !=
            Sw92AsymmetricFeedReferenceStatus::selected ||
        !selection.initial_stability.reference_family ||
        *selection.initial_stability.reference_family !=
            selection.single_phase_candidate->family ||
        selection.single_phase_candidate->composition != selection.feed) {
        return false;
    }
    const auto& leg = *selection.initial_stability.reference_family ==
            thermodynamics::SwPhaseFamily::aqueous
        ? selection.initial_stability.aqueous
        : selection.initial_stability.nonaqueous;
    if (!leg.feed_reference || !std::isfinite(leg.feed_reduced_gibbs)) {
        return false;
    }
    const auto& candidate = *selection.single_phase_candidate;
    return candidate.activity.branch == leg.feed_reference->branch &&
           candidate.activity.smooth == leg.feed_reference->smooth &&
           candidate.activity.ln_phi == leg.feed_reference->ln_phi &&
           candidate.reduced_gibbs == leg.feed_reduced_gibbs;
}

struct Sw92Max2PairEvidence {
    double reduced_gibbs{};
    double gibbs_roundoff_guard{};
    double chemical_potential_norm{};
    double mass_absolute{};
    double mass_relative{};
    double log_k_contrast{};
    std::vector<double> common_log_activity;
};

/// Recompute the acceptance-critical continuous evidence from retained phase
/// fractions/compositions/activities. This avoids trusting stale cached Gibbs,
/// residual, mass-balance or common-tangent fields during authoritative review.
inline std::optional<Sw92Max2PairEvidence> sw92_max2_recompute_pair_evidence(
    const Sw92AsymmetricFixedPairState& pair,
    std::span<const double> feed,
    const Sw92AsymmetricFixedPairOptions& options) {
    if (pair.phase0.composition.size() != feed.size() ||
        pair.phase1.composition.size() != feed.size() ||
        pair.phase0.activity.ln_phi.size() != feed.size() ||
        pair.phase1.activity.ln_phi.size() != feed.size() ||
        !pair.phase0.activity.smooth || !pair.phase1.activity.smooth ||
        !std::isfinite(pair.phase0.mole_phase_fraction) ||
        !std::isfinite(pair.phase1.mole_phase_fraction) ||
        !(pair.phase0.mole_phase_fraction > options.minimum_phase_fraction) ||
        !(pair.phase1.mole_phase_fraction > options.minimum_phase_fraction)) {
        return std::nullopt;
    }
    const double fraction_sum =
        pair.phase0.mole_phase_fraction + pair.phase1.mole_phase_fraction;
    if (!std::isfinite(fraction_sum) ||
        std::abs(fraction_sum - 1.0) > 256.0 * stability_eps) {
        return std::nullopt;
    }
    try {
        (void)stability_check_composition(feed);
        (void)stability_check_composition(pair.phase0.composition);
        (void)stability_check_composition(pair.phase1.composition);
        stability_check_phase(pair.phase0.activity, feed.size());
        stability_check_phase(pair.phase1.activity, feed.size());
    } catch (const std::exception&) {
        return std::nullopt;
    }

    Sw92Max2PairEvidence evidence;
    evidence.common_log_activity.resize(feed.size());
    double gibbs_correction = 0.0;
    double magnitude = 1.0;
    double magnitude_correction = 0.0;
    for (std::size_t i = 0; i < feed.size(); ++i) {
        if (feed[i] == 0.0) {
            if (pair.phase0.composition[i] != 0.0 ||
                pair.phase1.composition[i] != 0.0) {
                return std::nullopt;
            }
            continue;
        }
        const double x0 = pair.phase0.composition[i];
        const double x1 = pair.phase1.composition[i];
        if (!(x0 > 0.0) || !(x1 > 0.0)) { return std::nullopt; }
        const double l0 = std::log(x0);
        const double l1 = std::log(x1);
        const double m0 = l0 + pair.phase0.activity.ln_phi[i];
        const double m1 = l1 + pair.phase1.activity.ln_phi[i];
        if (!std::isfinite(m0) || !std::isfinite(m1)) { return std::nullopt; }

        evidence.chemical_potential_norm = std::max(
            evidence.chemical_potential_norm, std::abs(m0 - m1));
        evidence.common_log_activity[i] = std::midpoint(m0, m1);
        evidence.log_k_contrast = std::max(
            evidence.log_k_contrast, std::abs(l1 - l0));

        const double recovered = std::fma(
            pair.phase1.mole_phase_fraction, x1,
            pair.phase0.mole_phase_fraction * x0);
        const double error = std::abs(recovered - feed[i]);
        evidence.mass_absolute = std::max(evidence.mass_absolute, error);
        evidence.mass_relative = std::max(
            evidence.mass_relative, error / feed[i]);

        stability_add(
            pair.phase0.mole_phase_fraction * x0 * m0 +
                pair.phase1.mole_phase_fraction * x1 * m1,
            evidence.reduced_gibbs, gibbs_correction);
        stability_add(
            pair.phase0.mole_phase_fraction * x0 *
                    (std::abs(l0) + std::abs(pair.phase0.activity.ln_phi[i])) +
                pair.phase1.mole_phase_fraction * x1 *
                    (std::abs(l1) + std::abs(pair.phase1.activity.ln_phi[i])),
            magnitude, magnitude_correction);
    }
    evidence.gibbs_roundoff_guard = 256.0 * stability_eps * magnitude;
    if (!std::isfinite(evidence.reduced_gibbs) ||
        !std::isfinite(evidence.gibbs_roundoff_guard) ||
        evidence.chemical_potential_norm > options.chemical_potential_tolerance ||
        evidence.mass_absolute > options.mass_absolute_tolerance ||
        evidence.mass_relative > options.mass_relative_tolerance ||
        evidence.log_k_contrast <= options.log_k_separation) {
        return std::nullopt;
    }
    return evidence;
}

} // namespace detail

[[nodiscard]] inline Sw92AsymmetricMax2Result finalize_sw92_asymmetric_max2_selection(
    Sw92AsymmetricPairSelectionResult selection,
    const thermodynamics::Sw92Phase<double>& model,
    Sw92AsymmetricStabilityOptions final_stability_options = {},
    Sw92AsymmetricStabilityStarts caller_final_starts = {}) {
    if (!detail::sw92_max2_model_matches_selection(selection, model) ||
        !detail::sw92_max2_stability_matches_selection(selection)) {
        throw std::invalid_argument(
            "SW92 asymmetric max2: selection/model or retained-stability identity mismatch");
    }
    if (std::string_view{selection.model_profile} !=
            thermodynamics::sw92_corrected_profile ||
        std::string_view{selection.phase_convention} !=
            thermodynamics::sw92_pt_convention ||
        std::string_view{selection.equilibrium_profile} !=
            sw92_xu_asymmetric_gibbs_profile ||
        std::string_view{selection.orchestration_convention} !=
            sw92_xu_asymmetric_orchestration_convention) {
        throw std::invalid_argument(
            "SW92 asymmetric max2: selection algorithm/model identity mismatch");
    }

    Sw92AsymmetricMax2Result result;
    result.options.selection = selection.options;
    result.options.final_stability = final_stability_options;
    result.pressure_pa = selection.pressure_pa;
    result.temperature_k = selection.temperature_k;
    result.feed = selection.feed;
    result.nacl_molality_mol_per_kg_water =
        selection.nacl_molality_mol_per_kg_water;
    result.dataset_id = selection.dataset_id;
    result.revision = selection.revision;
    result.component_ids = selection.component_ids;
    result.model_profile = selection.model_profile;
    result.phase_convention = selection.phase_convention;
    result.equilibrium_profile = selection.equilibrium_profile;
    result.selection = std::move(selection);

    if (result.selection.status ==
            Sw92AsymmetricPairSelectionStatus::single_phase_candidate_no_instability_found) {
        if (!detail::sw92_max2_single_phase_evidence_consistent(result.selection)) {
            result.status = Sw92AsymmetricMax2Status::indeterminate;
            result.diagnostic =
                "Gate 3B.2 single-phase status is inconsistent with retained Gate-3A evidence";
            return result;
        }
        const auto& source = *result.selection.single_phase_candidate;
        Sw92AsymmetricCandidatePhase phase;
        phase.family = source.family;
        phase.mole_phase_fraction = 1.0;
        phase.composition = source.composition;
        phase.activity = source.activity;
        result.candidate_phase_set = Sw92AsymmetricCandidatePhaseSet{{std::move(phase)}};
        result.status = Sw92AsymmetricMax2Status::single_phase_no_instability_found;
        result.diagnostic =
            "Gate 3A initial two-family finite stability found no instability; family-aware single phase accepted without forcing a split; not a global stability proof";
        return result;
    }

    const auto* pair = detail::sw92_max2_revalidate_pair_choice(result.selection);
    if (pair == nullptr) {
        result.status = Sw92AsymmetricMax2Status::indeterminate;
        result.diagnostic =
            "Gate 3B.2 did not retain one complete, uniquely revalidated pair candidate for final acceptance";
        return result;
    }
    const auto continuous = detail::sw92_max2_recompute_pair_evidence(
        *pair, result.feed, result.selection.options.fixed_pair);
    if (!continuous) {
        result.status = Sw92AsymmetricMax2Status::indeterminate;
        result.diagnostic =
            "selected pair retained stale or inconsistent equation/material-balance evidence";
        return result;
    }

    result.candidate_phase_set = Sw92AsymmetricCandidatePhaseSet{{
        detail::sw92_max2_from_fixed_pair_phase(pair->phase0),
        detail::sw92_max2_from_fixed_pair_phase(pair->phase1)}};
    result.selected_pair_reduced_gibbs = continuous->reduced_gibbs;
    result.selected_pair_gibbs_roundoff_guard =
        continuous->gibbs_roundoff_guard;

    const auto reference_family = result.selection.initial_stability.reference_family;
    if (!reference_family ||
        result.selection.initial_stability.feed_reference_status !=
            Sw92AsymmetricFeedReferenceStatus::selected) {
        result.status = Sw92AsymmetricMax2Status::indeterminate;
        result.diagnostic =
            "selected pair lacks the resolved lower-feed reference required for pair-vs-feed Gibbs comparison";
        return result;
    }
    const auto& feed_leg = *reference_family ==
            thermodynamics::SwPhaseFamily::aqueous
        ? result.selection.initial_stability.aqueous
        : result.selection.initial_stability.nonaqueous;
    if (!feed_leg.feed_reference || !std::isfinite(feed_leg.feed_reduced_gibbs)) {
        result.status = Sw92AsymmetricMax2Status::indeterminate;
        result.diagnostic =
            "pair-vs-feed Gibbs comparison lacks a representable lower-feed reference";
        return result;
    }

    try {
        result.lower_feed_reduced_gibbs = detail::sw92_feed_reduced_gibbs(
            result.feed, *feed_leg.feed_reference);
        result.lower_feed_gibbs_roundoff_guard =
            detail::sw92_max2_phase_gibbs_guard(
                result.feed, *feed_leg.feed_reference);
    } catch (const std::exception& error) {
        result.status = Sw92AsymmetricMax2Status::indeterminate;
        result.diagnostic = error.what();
        return result;
    }
    result.pair_minus_lower_feed_reduced_gibbs =
        result.selected_pair_reduced_gibbs - result.lower_feed_reduced_gibbs;
    result.pair_feed_gibbs_combined_guard =
        result.selected_pair_gibbs_roundoff_guard +
        result.lower_feed_gibbs_roundoff_guard;
    if (!std::isfinite(result.pair_minus_lower_feed_reduced_gibbs) ||
        !std::isfinite(result.pair_feed_gibbs_combined_guard)) {
        result.status = Sw92AsymmetricMax2Status::indeterminate;
        result.diagnostic = "pair-vs-feed Gibbs arithmetic is nonrepresentable";
        return result;
    }
    if (result.pair_minus_lower_feed_reduced_gibbs >
        result.pair_feed_gibbs_combined_guard) {
        result.status = Sw92AsymmetricMax2Status::pair_gibbs_above_feed;
        result.diagnostic =
            "selected pair reduced Gibbs is resolved above the feasible lower-envelope feed state; candidate retained but not accepted";
        return result;
    }

    result.common_reference_allowance =
        0.5 * continuous->chemical_potential_norm;
    if (!std::isfinite(result.common_reference_allowance) ||
        result.common_reference_allowance < 0.0) {
        result.status = Sw92AsymmetricMax2Status::indeterminate;
        result.diagnostic = "selected-pair common-reference allowance is nonrepresentable";
        return result;
    }

    detail::stability_check_options(final_stability_options.aqueous.stability);
    detail::stability_check_options(final_stability_options.nonaqueous.stability);
    result.final_aqueous_base_tpd_tolerance =
        final_stability_options.aqueous.stability.tpd_tolerance;
    result.final_nonaqueous_base_tpd_tolerance =
        final_stability_options.nonaqueous.stability.tpd_tolerance;

    Sw92AsymmetricStabilityOptions effective = final_stability_options;
    result.final_aqueous_effective_tpd_tolerance =
        detail::sw92_max2_effective_tpd_tolerance(
            result.final_aqueous_base_tpd_tolerance,
            result.common_reference_allowance);
    result.final_nonaqueous_effective_tpd_tolerance =
        detail::sw92_max2_effective_tpd_tolerance(
            result.final_nonaqueous_base_tpd_tolerance,
            result.common_reference_allowance);
    effective.aqueous.stability.tpd_tolerance =
        result.final_aqueous_effective_tpd_tolerance;
    effective.nonaqueous.stability.tpd_tolerance =
        result.final_nonaqueous_effective_tpd_tolerance;

    auto aqueous_starts = detail::sw92_max2_final_starts(
        caller_final_starts.aqueous, *pair);
    auto nonaqueous_starts = detail::sw92_max2_final_starts(
        caller_final_starts.nonaqueous, *pair);
    const Sw92AsymmetricStabilityStarts required_starts{
        aqueous_starts, nonaqueous_starts};

    result.final_stability = test_sw92_pt_asymmetric_stability_against(
        result.pressure_pa, result.temperature_k, result.feed,
        continuous->common_log_activity, model,
        result.nacl_molality_mol_per_kg_water, effective, required_starts);

    switch (result.final_stability->status) {
    case StabilityStatus::unstable:
        result.status = Sw92AsymmetricMax2Status::phase_set_unstable;
        result.diagnostic =
            "selected pair satisfies Gate 3B.1/3B.2 and pair-vs-feed Gibbs, but final AQ/NA common-tangent search found a robust negative TPD witness";
        return result;
    case StabilityStatus::indeterminate:
        result.status = Sw92AsymmetricMax2Status::indeterminate;
        result.diagnostic =
            "selected pair satisfies equations and Gibbs gates, but at least one required final family stability search is indeterminate";
        return result;
    case StabilityStatus::no_instability_found:
        result.status = Sw92AsymmetricMax2Status::two_phase_no_instability_found;
        result.diagnostic =
            "maximum-two-phase family-aware SW92/Xu candidate passed material balance, common chemical potentials, lower-envelope family checks, guarded Gibbs selection, pair-vs-feed Gibbs and both final finite family stability searches; not a global stability proof";
        return result;
    }

    result.status = Sw92AsymmetricMax2Status::indeterminate;
    result.diagnostic = "unknown final asymmetric stability status";
    return result;
}

[[nodiscard]] inline Sw92AsymmetricMax2Result solve_sw92_xu_asymmetric_max2(
    double pressure_pa, double temperature_k, std::span<const double> feed,
    const thermodynamics::Sw92Phase<double>& model,
    double nacl_molality_mol_per_kg_water,
    Sw92AsymmetricMax2Options options = {},
    Sw92AsymmetricStabilityStarts initial_starts = {},
    Sw92AsymmetricStabilityStarts final_starts = {}) {
    auto selection = orchestrate_sw92_asymmetric_pair_candidates(
        pressure_pa, temperature_k, feed, model,
        nacl_molality_mol_per_kg_water,
        options.selection, initial_starts);
    return finalize_sw92_asymmetric_max2_selection(
        std::move(selection), model, options.final_stability, final_starts);
}

} // namespace mpmc::flash

#endif // MPMC_FLASH_SW92_ASYMMETRIC_MAX2_HPP
