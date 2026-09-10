#ifndef MPMC_FLASH_SW92_PHASE_ASSIGNED_H_SIDE_WITNESS_HPP
#define MPMC_FLASH_SW92_PHASE_ASSIGNED_H_SIDE_WITNESS_HPP

#include <mpmc/flash/sw92_phase_assigned_joint.hpp>

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
#include <vector>

namespace mpmc::flash {

/// Gate C2a1 algorithm identity. This is an NA-only additional-phase witness
/// search from an already converged Profile-C W(AQ)+H(NA) C1 candidate. It is
/// not a phase-count decision and never publishes an accepted phase set.
inline constexpr std::string_view sw92_phase_assigned_h_side_na_witness_convention =
    "SW92-equilibrium/phase-assigned-aq-na-joint/h-side-na-phase-addition-witness/v1";

struct Sw92PhaseAssignedHSideWitnessOptions {
    thermodynamics::Sw92RootOptions nonaqueous_root_options;
    StabilityOptions stability;
    double log_composition_separation{1e-7};
};

enum class Sw92PhaseAssignedHSideWitnessStatus {
    source_candidate_unavailable,
    source_candidate_inconsistent,
    additional_nonaqueous_phase_witness_found,
    no_additional_nonaqueous_witness_found,
    indeterminate
};

/// One robustly negative NA trial from the finite search. The point is always
/// retained as diagnostic evidence. It is a usable H-split seed only when both
/// topology guards below pass; even then it is not an accepted new phase.
struct Sw92PhaseAssignedNaNegativeWitness {
    std::size_t trial_index{};
    TpdPoint point;
    double log_distance_from_retained_h{};
    double retained_w_minus_trial_water_fraction{};
    double water_role_roundoff_guard{};
    bool compositionally_distinct_from_retained_h{};
    bool water_role_admissible{};

    [[nodiscard]] bool usable_h_split_seed() const noexcept {
        return compositionally_distinct_from_retained_h && water_role_admissible;
    }
};

struct Sw92PhaseAssignedHSideWitnessResult {
    static constexpr bool global_stability_proven = false;
    static constexpr bool accepted_phase_set_published = false;

    Sw92PhaseAssignedHSideWitnessStatus status{
        Sw92PhaseAssignedHSideWitnessStatus::source_candidate_unavailable};
    Sw92PhaseAssignedHSideWitnessOptions options;

    double pressure_pa{};
    double temperature_k{};
    std::vector<double> feed;
    double nacl_molality_mol_per_kg_water{};
    std::string dataset_id;
    std::string revision;
    std::vector<std::string> component_ids;
    std::string model_profile{thermodynamics::sw92_corrected_profile};
    std::string phase_convention{thermodynamics::sw92_pt_convention};
    std::string equilibrium_profile{sw92_phase_assigned_aq_na_joint_profile};
    std::string witness_convention{sw92_phase_assigned_h_side_na_witness_convention};

    std::size_t water_index{};
    double retained_w_water_fraction{std::numeric_limits<double>::quiet_NaN()};
    double retained_h_water_fraction{std::numeric_limits<double>::quiet_NaN()};
    std::vector<double> common_log_activity;
    double common_reference_allowance{std::numeric_limits<double>::quiet_NaN()};
    double base_tpd_tolerance{std::numeric_limits<double>::quiet_NaN()};
    double effective_tpd_tolerance{std::numeric_limits<double>::quiet_NaN()};
    std::optional<TpdPoint> retained_h_trivial_point;

    std::optional<StabilityResult> nonaqueous_search;
    std::vector<Sw92PhaseAssignedNaNegativeWitness> negative_witnesses;
    std::string diagnostic;

    [[nodiscard]] std::size_t usable_witness_count() const noexcept {
        return static_cast<std::size_t>(std::count_if(
            negative_witnesses.begin(), negative_witnesses.end(),
            [](const Sw92PhaseAssignedNaNegativeWitness& witness) {
                return witness.usable_h_split_seed();
            }));
    }

    [[nodiscard]] bool has_usable_h_split_witness() const noexcept {
        return usable_witness_count() != 0U;
    }
};

namespace detail {

struct Sw92PhaseAssignedC1WitnessEvidence {
    std::vector<double> common_log_activity;
    double common_reference_allowance{};
    double retained_w_water_fraction{};
    double retained_h_water_fraction{};
    TpdPoint retained_h_trivial_point;
};

inline bool sw92_phase_assigned_h_side_model_matches(
    const Sw92PhaseAssignedJointResult& source,
    const thermodynamics::Sw92Phase<double>& model) {
    if (source.feed.size() != model.size() ||
        source.dataset_id != model.parameters().dataset_id() ||
        source.revision != model.parameters().revision() ||
        source.component_ids.size() != model.size()) {
        return false;
    }
    std::size_t index = 0;
    for (const auto& component : model.parameters().components().items()) {
        if (source.component_ids[index++] != component.id) { return false; }
    }
    return std::string_view{source.model_profile} == thermodynamics::sw92_corrected_profile &&
           std::string_view{source.phase_convention} == thermodynamics::sw92_pt_convention &&
           std::string_view{source.equilibrium_profile} ==
               sw92_phase_assigned_aq_na_joint_profile &&
           std::string_view{source.primitive_convention} ==
               sw92_phase_assigned_aq_na_joint_primitive;
}

inline std::optional<Sw92PhaseAssignedC1WitnessEvidence>
sw92_phase_assigned_revalidate_c1_for_h_side_witness(
    const Sw92PhaseAssignedJointResult& source,
    const thermodynamics::Sw92Phase<double>& model) {
    const auto* point = source.candidate();
    if (point == nullptr) { return std::nullopt; }
    const std::size_t n = source.feed.size();
    if (n == 0U || point->aqueous_phase.composition.size() != n ||
        point->nonaqueous_phase.composition.size() != n ||
        point->aqueous_phase.activity.ln_phi.size() != n ||
        point->nonaqueous_phase.activity.ln_phi.size() != n ||
        point->aqueous_phase.physical_role != Sw92PhysicalPhaseRole::aqueous ||
        point->nonaqueous_phase.physical_role != Sw92PhysicalPhaseRole::nonaqueous ||
        point->aqueous_phase.thermodynamic_family !=
            thermodynamics::SwPhaseFamily::aqueous ||
        point->nonaqueous_phase.thermodynamic_family !=
            thermodynamics::SwPhaseFamily::nonaqueous ||
        !point->aqueous_phase.activity.smooth ||
        !point->nonaqueous_phase.activity.smooth ||
        !std::isfinite(point->aqueous_phase.mole_phase_fraction) ||
        !std::isfinite(point->nonaqueous_phase.mole_phase_fraction) ||
        !(point->aqueous_phase.mole_phase_fraction > source.options.minimum_phase_fraction) ||
        !(point->nonaqueous_phase.mole_phase_fraction > source.options.minimum_phase_fraction)) {
        return std::nullopt;
    }

    try {
        (void)stability_check_composition(source.feed);
        (void)stability_check_composition(point->aqueous_phase.composition);
        (void)stability_check_composition(point->nonaqueous_phase.composition);
        stability_check_phase(point->aqueous_phase.activity, n);
        stability_check_phase(point->nonaqueous_phase.activity, n);
    } catch (const std::exception&) {
        return std::nullopt;
    }

    const double fraction_sum = point->aqueous_phase.mole_phase_fraction +
                                point->nonaqueous_phase.mole_phase_fraction;
    if (!std::isfinite(fraction_sum) ||
        std::abs(fraction_sum - 1.0) > 256.0 * stability_eps) {
        return std::nullopt;
    }

    Sw92PhaseAssignedC1WitnessEvidence evidence;
    evidence.common_log_activity.assign(n, 0.0);
    double residual_norm = 0.0;
    double mass_absolute = 0.0;
    double mass_relative = 0.0;
    double log_contrast = 0.0;
    for (std::size_t i = 0; i < n; ++i) {
        const double feed_i = source.feed[i];
        const double x = point->aqueous_phase.composition[i];
        const double y = point->nonaqueous_phase.composition[i];
        if (feed_i == 0.0) {
            if (x != 0.0 || y != 0.0) { return std::nullopt; }
            continue;
        }
        if (!(x > 0.0) || !(y > 0.0)) { return std::nullopt; }
        const double m_w = std::log(x) + point->aqueous_phase.activity.ln_phi[i];
        const double m_h = std::log(y) + point->nonaqueous_phase.activity.ln_phi[i];
        if (!std::isfinite(m_w) || !std::isfinite(m_h)) { return std::nullopt; }
        evidence.common_log_activity[i] = std::midpoint(m_w, m_h);
        residual_norm = std::max(residual_norm, std::abs(m_w - m_h));
        log_contrast = std::max(log_contrast, std::abs(std::log(y) - std::log(x)));

        const double reconstructed =
            point->aqueous_phase.mole_phase_fraction * x +
            point->nonaqueous_phase.mole_phase_fraction * y;
        const double error = std::abs(reconstructed - feed_i);
        mass_absolute = std::max(mass_absolute, error);
        mass_relative = std::max(mass_relative, error / feed_i);
    }
    if (residual_norm > source.options.chemical_potential_tolerance ||
        mass_absolute > source.options.mass_absolute_tolerance ||
        mass_relative > source.options.mass_relative_tolerance ||
        log_contrast <= source.options.log_k_separation) {
        return std::nullopt;
    }

    const std::size_t water_index = model.parameters().water_index();
    evidence.retained_w_water_fraction = point->aqueous_phase.composition[water_index];
    evidence.retained_h_water_fraction = point->nonaqueous_phase.composition[water_index];
    const double role_guard = 256.0 * stability_eps *
        (1.0 + std::abs(evidence.retained_w_water_fraction) +
         std::abs(evidence.retained_h_water_fraction));
    if (!std::isfinite(role_guard) ||
        !(evidence.retained_w_water_fraction - evidence.retained_h_water_fraction >
          role_guard)) {
        return std::nullopt;
    }

    evidence.common_reference_allowance = 0.5 * residual_norm;
    StabilityPhase equivalent_reference;
    equivalent_reference.ln_phi.assign(n, 0.0);
    for (std::size_t i = 0; i < n; ++i) {
        if (source.feed[i] > 0.0) {
            equivalent_reference.ln_phi[i] =
                evidence.common_log_activity[i] - std::log(source.feed[i]);
        }
    }
    try {
        evidence.retained_h_trivial_point = tangent_plane_distance(
            point->nonaqueous_phase.composition, source.feed,
            point->nonaqueous_phase.activity, equivalent_reference);
    } catch (const std::exception&) {
        return std::nullopt;
    }
    const double trivial_allowance = evidence.common_reference_allowance +
                                     evidence.retained_h_trivial_point.roundoff_guard;
    if (!std::isfinite(trivial_allowance) ||
        std::abs(evidence.retained_h_trivial_point.value) > trivial_allowance) {
        return std::nullopt;
    }
    return evidence;
}

inline double sw92_phase_assigned_h_side_effective_tpd_tolerance(
    double base, double allowance) {
    if (!std::isfinite(base) || base < 0.0 ||
        !std::isfinite(allowance) || allowance < 0.0 ||
        allowance > std::numeric_limits<double>::max() - base) {
        throw std::domain_error(
            "SW92 phase-assigned H-side witness: invalid TPD tolerance/reference allowance");
    }
    const double effective = base + allowance;
    if (!std::isfinite(effective)) {
        throw std::domain_error(
            "SW92 phase-assigned H-side witness: nonrepresentable effective TPD tolerance");
    }
    return effective;
}

inline Sw92PhaseAssignedNaNegativeWitness sw92_phase_assigned_classify_na_witness(
    std::size_t trial_index, const TpdPoint& point,
    std::span<const double> retained_h,
    double retained_w_water, std::size_t water_index,
    double log_composition_separation) {
    Sw92PhaseAssignedNaNegativeWitness witness;
    witness.trial_index = trial_index;
    witness.point = point;
    if (point.composition.size() != retained_h.size() ||
        water_index >= point.composition.size()) {
        return witness;
    }
    double distance = 0.0;
    for (std::size_t i = 0; i < retained_h.size(); ++i) {
        if (retained_h[i] == 0.0 && point.composition[i] == 0.0) { continue; }
        if (!(retained_h[i] > 0.0) || !(point.composition[i] > 0.0)) {
            distance = std::numeric_limits<double>::infinity();
            break;
        }
        distance = std::max(
            distance,
            std::abs(std::log(point.composition[i]) - std::log(retained_h[i])));
    }
    witness.log_distance_from_retained_h = distance;
    witness.compositionally_distinct_from_retained_h =
        std::isfinite(distance) && distance > log_composition_separation;

    const double trial_water = point.composition[water_index];
    witness.retained_w_minus_trial_water_fraction = retained_w_water - trial_water;
    witness.water_role_roundoff_guard = 256.0 * stability_eps *
        (1.0 + std::abs(retained_w_water) + std::abs(trial_water));
    witness.water_role_admissible =
        std::isfinite(witness.retained_w_minus_trial_water_fraction) &&
        std::isfinite(witness.water_role_roundoff_guard) &&
        witness.retained_w_minus_trial_water_fraction > witness.water_role_roundoff_guard;
    return witness;
}

} // namespace detail

/// Search for evidence that the retained Profile-C hydrocarbon phase H can add a
/// second non-aqueous phase while W remains. All trial properties use SW92 NA.
///
/// A negative trial is only a usable H-split seed when it is compositionally
/// distinct from retained H and remains water-poorer than retained W. This Gate
/// C2a1 adapter publishes no phase set and does not classify H as liquid/vapor.
[[nodiscard]] inline Sw92PhaseAssignedHSideWitnessResult
test_sw92_phase_assigned_h_side_na_witness(
    const Sw92PhaseAssignedJointResult& source,
    const thermodynamics::Sw92Phase<double>& model,
    Sw92PhaseAssignedHSideWitnessOptions options = {},
    std::span<const std::vector<double>> extra_starts = {}) {
    detail::stability_check_options(options.stability);
    if (!std::isfinite(options.log_composition_separation) ||
        !(options.log_composition_separation > 0.0)) {
        throw std::invalid_argument(
            "SW92 phase-assigned H-side witness: positive finite composition separation required");
    }
    if (!detail::sw92_phase_assigned_h_side_model_matches(source, model)) {
        throw std::invalid_argument(
            "SW92 phase-assigned H-side witness: C1 result/model snapshot mismatch");
    }

    Sw92PhaseAssignedHSideWitnessResult result;
    result.options = options;
    result.pressure_pa = source.pressure_pa;
    result.temperature_k = source.temperature_k;
    result.feed = source.feed;
    result.nacl_molality_mol_per_kg_water = source.nacl_molality_mol_per_kg_water;
    result.dataset_id = source.dataset_id;
    result.revision = source.revision;
    result.component_ids = source.component_ids;
    result.water_index = model.parameters().water_index();
    result.base_tpd_tolerance = options.stability.tpd_tolerance;

    if (source.candidate() == nullptr) {
        result.status = Sw92PhaseAssignedHSideWitnessStatus::source_candidate_unavailable;
        result.diagnostic =
            "Profile-C C2a1 requires an admissible C1 W(AQ)+H(NA) candidate";
        return result;
    }

    const auto evidence = detail::sw92_phase_assigned_revalidate_c1_for_h_side_witness(
        source, model);
    if (!evidence) {
        result.status = Sw92PhaseAssignedHSideWitnessStatus::source_candidate_inconsistent;
        result.diagnostic =
            "retained C1 state failed independent role/equilibrium/material-balance revalidation";
        return result;
    }
    result.common_log_activity = evidence->common_log_activity;
    result.common_reference_allowance = evidence->common_reference_allowance;
    result.retained_w_water_fraction = evidence->retained_w_water_fraction;
    result.retained_h_water_fraction = evidence->retained_h_water_fraction;
    result.retained_h_trivial_point = evidence->retained_h_trivial_point;
    result.effective_tpd_tolerance =
        detail::sw92_phase_assigned_h_side_effective_tpd_tolerance(
            result.base_tpd_tolerance, result.common_reference_allowance);

    if (extra_starts.size() > std::numeric_limits<std::size_t>::max() - 2U) {
        throw std::length_error(
            "SW92 phase-assigned H-side witness: start-count overflow");
    }
    std::vector<std::vector<double>> starts;
    starts.reserve(extra_starts.size() + 2U);
    const auto& c1 = *source.point;
    starts.push_back(c1.aqueous_phase.composition);   // diagnostic W composition under NA
    starts.push_back(c1.nonaqueous_phase.composition); // trivial retained H phase
    for (const auto& start : extra_starts) { starts.push_back(start); }

    StabilityOptions search_options = options.stability;
    search_options.tpd_tolerance = result.effective_tpd_tolerance;
    Sw92FamilyStabilityEvaluator evaluator(
        model, source.nacl_molality_mol_per_kg_water,
        thermodynamics::SwPhaseFamily::nonaqueous,
        options.nonaqueous_root_options);
    result.nonaqueous_search = test_pt_stability_against(
        source.pressure_pa, source.temperature_k, source.feed,
        result.common_log_activity, evaluator, search_options, starts);

    const auto& search = *result.nonaqueous_search;
    for (std::size_t trial_index = 0; trial_index < search.trials.size(); ++trial_index) {
        const auto& trial = search.trials[trial_index];
        if (trial.status != StabilityTrialStatus::negative_tpd || !trial.point ||
            !detail::stability_negative(*trial.point, search.options)) {
            continue;
        }
        result.negative_witnesses.push_back(
            detail::sw92_phase_assigned_classify_na_witness(
                trial_index, *trial.point, c1.nonaqueous_phase.composition,
                result.retained_w_water_fraction, result.water_index,
                options.log_composition_separation));
    }

    if (result.has_usable_h_split_witness()) {
        result.status = Sw92PhaseAssignedHSideWitnessStatus::
            additional_nonaqueous_phase_witness_found;
        result.diagnostic =
            "finite NA-only search found a robust, compositionally distinct and role-admissible additional nonaqueous witness; this is an H-split seed, not an accepted phase";
        return result;
    }
    if (search.status == StabilityStatus::no_instability_found) {
        result.status = Sw92PhaseAssignedHSideWitnessStatus::
            no_additional_nonaqueous_witness_found;
        result.diagnostic =
            "finite NA-only search found no additional nonaqueous instability; not a global stability or phase-count proof";
        return result;
    }

    result.status = Sw92PhaseAssignedHSideWitnessStatus::indeterminate;
    result.diagnostic = search.status == StabilityStatus::unstable
        ? "NA-only search found negative mathematical trials, but none qualifies as a role-admissible distinct H-split seed"
        : "NA-only search is numerically indeterminate; no authoritative topology conclusion is permitted";
    return result;
}

} // namespace mpmc::flash

#endif // MPMC_FLASH_SW92_PHASE_ASSIGNED_H_SIDE_WITNESS_HPP
