#ifndef MPMC_FLASH_SW92_PHASE_ASSIGNED_THREE_PHASE_CLOSURE_HPP
#define MPMC_FLASH_SW92_PHASE_ASSIGNED_THREE_PHASE_CLOSURE_HPP

#include <mpmc/flash/sw92_phase_assigned_three_phase.hpp>

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

inline constexpr std::string_view sw92_phase_assigned_c2b2_closure_convention =
    "SW92-equilibrium/phase-assigned-aq-na-joint/w-present-c2b1-"
    "nested-topology-na-review/v1";

struct Sw92PhaseAssignedC2b2Options {
    StabilityOptions stability;
    double log_composition_separation{1e-7};
};

enum class Sw92PhaseAssignedC2b2Status {
    source_chain_unavailable,
    source_chain_inconsistent,
    route_to_w_h,
    route_to_no_w_h0_h1,
    single_phase_endpoint_unresolved,
    nested_gibbs_inconsistent,
    higher_phase_count_witness_found,
    w_present_h_multiplicity_locally_closed,
    indeterminate
};

struct Sw92PhaseAssignedC2b2NaNegativeWitness {
    std::size_t trial_index{};
    TpdPoint point;
    double log_distance_from_h0{};
    double log_distance_from_h1{};
    double retained_w_minus_trial_water_fraction{};
    double water_role_roundoff_guard{};
    bool compositionally_distinct_from_h0{};
    bool compositionally_distinct_from_h1{};
    bool water_role_admissible{};

    [[nodiscard]] bool usable_additional_h_witness() const noexcept {
        return compositionally_distinct_from_h0 &&
               compositionally_distinct_from_h1 && water_role_admissible;
    }
};

struct Sw92PhaseAssignedC2b2Result {
    static constexpr bool global_stability_proven = false;
    static constexpr bool accepted_phase_set_published = false;
    static constexpr bool morphology_resolved = false;

    Sw92PhaseAssignedC2b2Status status{
        Sw92PhaseAssignedC2b2Status::source_chain_unavailable};
    Sw92PhaseAssignedC2b2Options options;

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
    std::string review_convention{sw92_phase_assigned_c2b2_closure_convention};

    std::optional<std::size_t> selected_c2a1_witness_index;
    std::size_t water_index{};
    std::size_t revalidation_property_evaluations{};

    double source_c1_reduced_gibbs{std::numeric_limits<double>::quiet_NaN()};
    double source_c1_gibbs_roundoff_guard{std::numeric_limits<double>::quiet_NaN()};
    double candidate_c2b1_reduced_gibbs{std::numeric_limits<double>::quiet_NaN()};
    double candidate_c2b1_gibbs_roundoff_guard{std::numeric_limits<double>::quiet_NaN()};
    double nested_gibbs_guard{std::numeric_limits<double>::quiet_NaN()};

    std::vector<double> common_log_activity;
    double common_reference_allowance{std::numeric_limits<double>::quiet_NaN()};
    double base_tpd_tolerance{std::numeric_limits<double>::quiet_NaN()};
    double effective_tpd_tolerance{std::numeric_limits<double>::quiet_NaN()};
    std::optional<TpdPoint> retained_h0_trivial_point;
    std::optional<TpdPoint> retained_h1_trivial_point;

    std::optional<StabilityResult> nonaqueous_search;
    std::vector<Sw92PhaseAssignedC2b2NaNegativeWitness> negative_witnesses;
    std::optional<StabilityPropertyIssue> property_issue;
    std::string diagnostic;

    [[nodiscard]] std::size_t usable_additional_h_witness_count() const noexcept {
        return static_cast<std::size_t>(std::count_if(
            negative_witnesses.begin(), negative_witnesses.end(),
            [](const Sw92PhaseAssignedC2b2NaNegativeWitness& witness) {
                return witness.usable_additional_h_witness();
            }));
    }

    [[nodiscard]] bool w_present_locally_closed() const noexcept {
        return status == Sw92PhaseAssignedC2b2Status::
                             w_present_h_multiplicity_locally_closed;
    }
};

namespace detail {

struct Sw92PhaseAssignedC2b2C1Evidence {
    double reduced_gibbs{};
    double gibbs_roundoff_guard{};
    double retained_w_water_fraction{};
};

struct Sw92PhaseAssignedC2b2ThreePhaseEvidence {
    std::vector<double> common_log_activity;
    double common_reference_allowance{};
    double reduced_gibbs{};
    double gibbs_roundoff_guard{};
    TpdPoint retained_h0_trivial_point;
    TpdPoint retained_h1_trivial_point;
    bool w_disappeared{};
    bool h0_disappeared{};
    bool h1_disappeared{};
};

inline void sw92_phase_assigned_c2b2_check_options(
    const Sw92PhaseAssignedC2b2Options& options) {
    stability_check_options(options.stability);
    if (!std::isfinite(options.log_composition_separation) ||
        !(options.log_composition_separation > 0.0)) {
        throw std::invalid_argument(
            "SW92 phase-assigned C2b.2: positive finite composition separation required");
    }
}

inline bool sw92_phase_assigned_c2b2_model_matches(
    const Sw92PhaseAssignedJointResult& c1,
    const thermodynamics::Sw92Phase<double>& model) {
    return sw92_phase_assigned_h_side_model_matches(c1, model);
}

inline bool sw92_phase_assigned_c2b2_three_phase_matches(
    const Sw92PhaseAssignedJointResult& c1,
    const Sw92PhaseAssignedThreePhaseResult& c2b1) {
    return c1.pressure_pa == c2b1.pressure_pa &&
           c1.temperature_k == c2b1.temperature_k &&
           c1.nacl_molality_mol_per_kg_water ==
               c2b1.nacl_molality_mol_per_kg_water &&
           c1.feed == c2b1.feed &&
           c1.dataset_id == c2b1.dataset_id &&
           c1.revision == c2b1.revision &&
           c1.component_ids == c2b1.component_ids &&
           std::string_view{c2b1.model_profile} == thermodynamics::sw92_corrected_profile &&
           std::string_view{c2b1.phase_convention} == thermodynamics::sw92_pt_convention &&
           std::string_view{c2b1.equilibrium_profile} ==
               sw92_phase_assigned_aq_na_joint_profile &&
           std::string_view{c2b1.primitive_convention} ==
               sw92_phase_assigned_c2b1_three_phase_convention;
}

inline bool sw92_phase_assigned_c2b2_close(
    double actual, double expected) noexcept {
    const double guard = 2048.0 * stability_eps *
        (1.0 + std::abs(actual) + std::abs(expected));
    return std::isfinite(actual) && std::isfinite(expected) &&
           std::abs(actual - expected) <= guard;
}

inline bool sw92_phase_assigned_c2b2_source_seed_matches(
    const Sw92PhaseAssignedJointResult& c1,
    const Sw92PhaseAssignedHSideWitnessResult& c2a1,
    const Sw92PhaseAssignedThreePhaseResult& c2b1) {
    if (!c2b1.source_witness_index ||
        *c2b1.source_witness_index >= c2a1.negative_witnesses.size()) {
        return false;
    }
    const auto* c1_point = c1.candidate();
    if (c1_point == nullptr ||
        c2b1.initial_log_k_h0.size() != c1.feed.size() ||
        c2b1.initial_log_k_h1.size() != c1.feed.size()) {
        return false;
    }
    const auto& witness = c2a1.negative_witnesses[*c2b1.source_witness_index];
    if (witness.point.composition.size() != c1.feed.size()) { return false; }
    for (std::size_t i = 0; i < c1.feed.size(); ++i) {
        if (c1.feed[i] == 0.0) { continue; }
        const double w = c1_point->aqueous_phase.composition[i];
        const double h = c1_point->nonaqueous_phase.composition[i];
        const double added = witness.point.composition[i];
        if (!(w > 0.0) || !(h > 0.0) || !(added > 0.0)) { return false; }
        const double expected_h0 = std::log(h) - std::log(w);
        const double expected_h1 = std::log(added) - std::log(w);
        if (!sw92_phase_assigned_c2b2_close(
                c2b1.initial_log_k_h0[i], expected_h0) ||
            !sw92_phase_assigned_c2b2_close(
                c2b1.initial_log_k_h1[i], expected_h1)) {
            return false;
        }
    }
    const double retained_h = c1_point->nonaqueous_phase.mole_phase_fraction;
    const double expected0 = retained_h *
        (1.0 - c2b1.options.new_hydrocarbon_seed_share);
    const double expected1 = retained_h * c2b1.options.new_hydrocarbon_seed_share;
    return sw92_phase_assigned_c2b2_close(
               c2b1.initial_hydrocarbon_fractions[0], expected0) &&
           sw92_phase_assigned_c2b2_close(
               c2b1.initial_hydrocarbon_fractions[1], expected1);
}

inline std::optional<Sw92PhaseAssignedC2b2C1Evidence>
sw92_phase_assigned_c2b2_revalidate_c1(
    const Sw92PhaseAssignedJointResult& c1,
    const thermodynamics::Sw92Phase<double>& model,
    std::size_t& evaluations) {
    const auto basic = sw92_phase_assigned_revalidate_c1_for_h_side_witness(c1, model);
    const auto* point = c1.candidate();
    if (!basic || point == nullptr) { return std::nullopt; }

    Sw92FamilyStabilityEvaluator aqueous(
        model, c1.nacl_molality_mol_per_kg_water,
        thermodynamics::SwPhaseFamily::aqueous,
        c1.options.aqueous_root_options);
    Sw92FamilyStabilityEvaluator nonaqueous(
        model, c1.nacl_molality_mol_per_kg_water,
        thermodynamics::SwPhaseFamily::nonaqueous,
        c1.options.nonaqueous_root_options);
    auto selected_w = aqueous.evaluate_selected(
        c1.pressure_pa, c1.temperature_k, point->aqueous_phase.composition);
    ++evaluations;
    auto selected_h = nonaqueous.evaluate_selected(
        c1.pressure_pa, c1.temperature_k, point->nonaqueous_phase.composition);
    ++evaluations;
    if (!selected_w.activity.smooth || !selected_h.activity.smooth) {
        throw StabilityPropertyError(
            StabilityPropertyIssue::root_topology,
            "SW92 phase-assigned C2b.2: C1 same-family minimum-root envelope is nonsmooth");
    }

    Sw92PhaseAssignedC2b2C1Evidence evidence;
    evidence.retained_w_water_fraction = basic->retained_w_water_fraction;
    double residual_norm = 0.0;
    double magnitude = 1.0;
    double magnitude_correction = 0.0;
    double gibbs_correction = 0.0;
    for (std::size_t i = 0; i < c1.feed.size(); ++i) {
        if (c1.feed[i] == 0.0) { continue; }
        const double xw = point->aqueous_phase.composition[i];
        const double xh = point->nonaqueous_phase.composition[i];
        if (!(xw > 0.0) || !(xh > 0.0)) { return std::nullopt; }
        const double lw = std::log(xw);
        const double lh = std::log(xh);
        const double mw = lw + selected_w.activity.ln_phi[i];
        const double mh = lh + selected_h.activity.ln_phi[i];
        if (!std::isfinite(mw) || !std::isfinite(mh)) { return std::nullopt; }
        residual_norm = std::max(residual_norm, std::abs(mw - mh));
        stability_add(
            point->aqueous_phase.mole_phase_fraction * xw * mw +
            point->nonaqueous_phase.mole_phase_fraction * xh * mh,
            evidence.reduced_gibbs, gibbs_correction);
        stability_add(
            point->aqueous_phase.mole_phase_fraction * xw *
                (std::abs(lw) + std::abs(selected_w.activity.ln_phi[i])) +
            point->nonaqueous_phase.mole_phase_fraction * xh *
                (std::abs(lh) + std::abs(selected_h.activity.ln_phi[i])),
            magnitude, magnitude_correction);
    }
    evidence.gibbs_roundoff_guard = 256.0 * stability_eps * magnitude;
    if (residual_norm > c1.options.chemical_potential_tolerance ||
        !std::isfinite(evidence.reduced_gibbs) ||
        !std::isfinite(evidence.gibbs_roundoff_guard)) {
        return std::nullopt;
    }
    return evidence;
}

inline std::optional<Sw92PhaseAssignedC2b2ThreePhaseEvidence>
sw92_phase_assigned_c2b2_revalidate_three_phase(
    const Sw92PhaseAssignedThreePhaseResult& c2b1,
    const thermodynamics::Sw92Phase<double>& model,
    std::size_t& evaluations) {
    if (!c2b1.point || !c2b1.equations_converged()) { return std::nullopt; }
    if (c2b1.status != Sw92PhaseAssignedThreePhaseStatus::converged_candidate &&
        c2b1.status != Sw92PhaseAssignedThreePhaseStatus::
                           hydrocarbon_phase_disappearance &&
        c2b1.status != Sw92PhaseAssignedThreePhaseStatus::
                           aqueous_phase_disappearance) {
        return std::nullopt;
    }
    const auto& point = *c2b1.point;
    const std::size_t n = c2b1.feed.size();
    if (n == 0U || point.aqueous_phase.composition.size() != n ||
        point.hydrocarbon0_phase.composition.size() != n ||
        point.hydrocarbon1_phase.composition.size() != n ||
        point.aqueous_phase.physical_role != Sw92PhysicalPhaseRole::aqueous ||
        point.hydrocarbon0_phase.physical_role != Sw92PhysicalPhaseRole::nonaqueous ||
        point.hydrocarbon1_phase.physical_role != Sw92PhysicalPhaseRole::nonaqueous ||
        point.aqueous_phase.thermodynamic_family !=
            thermodynamics::SwPhaseFamily::aqueous ||
        point.hydrocarbon0_phase.thermodynamic_family !=
            thermodynamics::SwPhaseFamily::nonaqueous ||
        point.hydrocarbon1_phase.thermodynamic_family !=
            thermodynamics::SwPhaseFamily::nonaqueous) {
        return std::nullopt;
    }
    try {
        (void)stability_check_composition(c2b1.feed);
        (void)stability_check_composition(point.aqueous_phase.composition);
        (void)stability_check_composition(point.hydrocarbon0_phase.composition);
        (void)stability_check_composition(point.hydrocarbon1_phase.composition);
    } catch (const std::exception&) {
        return std::nullopt;
    }

    const double fw = point.aqueous_phase.mole_phase_fraction;
    const double f0 = point.hydrocarbon0_phase.mole_phase_fraction;
    const double f1 = point.hydrocarbon1_phase.mole_phase_fraction;
    if (!std::isfinite(fw) || !std::isfinite(f0) || !std::isfinite(f1) ||
        fw < 0.0 || f0 < 0.0 || f1 < 0.0 ||
        std::abs((fw + f0 + f1) - 1.0) > 256.0 * stability_eps) {
        return std::nullopt;
    }

    Sw92FamilyStabilityEvaluator aqueous(
        model, c2b1.nacl_molality_mol_per_kg_water,
        thermodynamics::SwPhaseFamily::aqueous,
        c2b1.options.aqueous_root_options);
    Sw92FamilyStabilityEvaluator nonaqueous(
        model, c2b1.nacl_molality_mol_per_kg_water,
        thermodynamics::SwPhaseFamily::nonaqueous,
        c2b1.options.nonaqueous_root_options);
    auto selected_w = aqueous.evaluate_selected(
        c2b1.pressure_pa, c2b1.temperature_k, point.aqueous_phase.composition);
    ++evaluations;
    auto selected_h0 = nonaqueous.evaluate_selected(
        c2b1.pressure_pa, c2b1.temperature_k, point.hydrocarbon0_phase.composition);
    ++evaluations;
    auto selected_h1 = nonaqueous.evaluate_selected(
        c2b1.pressure_pa, c2b1.temperature_k, point.hydrocarbon1_phase.composition);
    ++evaluations;
    if (!selected_w.activity.smooth || !selected_h0.activity.smooth ||
        !selected_h1.activity.smooth) {
        throw StabilityPropertyError(
            StabilityPropertyIssue::root_topology,
            "SW92 phase-assigned C2b.2: C2b.1 same-family minimum-root envelope is nonsmooth");
    }

    Sw92PhaseAssignedC2b2ThreePhaseEvidence evidence;
    evidence.common_log_activity.assign(n, 0.0);
    double residual_norm = 0.0;
    double mass_absolute = 0.0;
    double mass_relative = 0.0;
    double magnitude = 1.0;
    double magnitude_correction = 0.0;
    double gibbs_correction = 0.0;
    for (std::size_t i = 0; i < n; ++i) {
        const double feed_i = c2b1.feed[i];
        const double xw = point.aqueous_phase.composition[i];
        const double x0 = point.hydrocarbon0_phase.composition[i];
        const double x1 = point.hydrocarbon1_phase.composition[i];
        if (feed_i == 0.0) {
            if (xw != 0.0 || x0 != 0.0 || x1 != 0.0) { return std::nullopt; }
            continue;
        }
        if (!(xw > 0.0) || !(x0 > 0.0) || !(x1 > 0.0)) {
            return std::nullopt;
        }
        const double lw = std::log(xw);
        const double l0 = std::log(x0);
        const double l1 = std::log(x1);
        const double mw = lw + selected_w.activity.ln_phi[i];
        const double m0 = l0 + selected_h0.activity.ln_phi[i];
        const double m1 = l1 + selected_h1.activity.ln_phi[i];
        if (!std::isfinite(mw) || !std::isfinite(m0) || !std::isfinite(m1)) {
            return std::nullopt;
        }
        const double common = (mw + m0 + m1) / 3.0;
        evidence.common_log_activity[i] = common;
        evidence.common_reference_allowance = std::max({
            evidence.common_reference_allowance,
            std::abs(mw - common), std::abs(m0 - common),
            std::abs(m1 - common)});
        residual_norm = std::max({residual_norm, std::abs(mw - m0),
                                  std::abs(mw - m1)});

        const double recovered = fw * xw + f0 * x0 + f1 * x1;
        const double error = std::abs(recovered - feed_i);
        mass_absolute = std::max(mass_absolute, error);
        mass_relative = std::max(mass_relative, error / feed_i);

        stability_add(fw * xw * mw + f0 * x0 * m0 + f1 * x1 * m1,
                      evidence.reduced_gibbs, gibbs_correction);
        stability_add(
            fw * xw * (std::abs(lw) + std::abs(selected_w.activity.ln_phi[i])) +
            f0 * x0 * (std::abs(l0) + std::abs(selected_h0.activity.ln_phi[i])) +
            f1 * x1 * (std::abs(l1) + std::abs(selected_h1.activity.ln_phi[i])),
            magnitude, magnitude_correction);
    }
    evidence.gibbs_roundoff_guard = 256.0 * stability_eps * magnitude;
    if (residual_norm > c2b1.options.chemical_potential_tolerance ||
        mass_absolute > c2b1.options.mass_absolute_tolerance ||
        mass_relative > c2b1.options.mass_relative_tolerance ||
        !std::isfinite(evidence.common_reference_allowance) ||
        !std::isfinite(evidence.reduced_gibbs) ||
        !std::isfinite(evidence.gibbs_roundoff_guard)) {
        return std::nullopt;
    }

    const std::size_t water_index = model.parameters().water_index();
    const double w_water = point.aqueous_phase.composition[water_index];
    const double h0_water = point.hydrocarbon0_phase.composition[water_index];
    const double h1_water = point.hydrocarbon1_phase.composition[water_index];
    const double role_guard = 256.0 * stability_eps *
        (1.0 + std::abs(w_water) + std::abs(h0_water) + std::abs(h1_water));
    if (!std::isfinite(role_guard) ||
        !(w_water - h0_water > role_guard) ||
        !(w_water - h1_water > role_guard)) {
        return std::nullopt;
    }

    evidence.w_disappeared = fw <= c2b1.options.minimum_phase_fraction;
    evidence.h0_disappeared = f0 <= c2b1.options.minimum_phase_fraction;
    evidence.h1_disappeared = f1 <= c2b1.options.minimum_phase_fraction;
    if (!evidence.w_disappeared && !evidence.h0_disappeared &&
        !evidence.h1_disappeared) {
        const double d_w0 = sw92_phase_assigned_log_distance(
            point.aqueous_phase.composition,
            point.hydrocarbon0_phase.composition, c2b1.feed);
        const double d_w1 = sw92_phase_assigned_log_distance(
            point.aqueous_phase.composition,
            point.hydrocarbon1_phase.composition, c2b1.feed);
        const double d_01 = sw92_phase_assigned_log_distance(
            point.hydrocarbon0_phase.composition,
            point.hydrocarbon1_phase.composition, c2b1.feed);
        if (!(d_w0 > c2b1.options.log_composition_separation) ||
            !(d_w1 > c2b1.options.log_composition_separation) ||
            !(d_01 > c2b1.options.log_composition_separation)) {
            return std::nullopt;
        }
    }

    if (c2b1.status == Sw92PhaseAssignedThreePhaseStatus::converged_candidate) {
        if (evidence.w_disappeared || evidence.h0_disappeared || evidence.h1_disappeared) {
            return std::nullopt;
        }
    } else if (c2b1.status ==
               Sw92PhaseAssignedThreePhaseStatus::hydrocarbon_phase_disappearance) {
        if (evidence.w_disappeared ||
            (!evidence.h0_disappeared && !evidence.h1_disappeared)) {
            return std::nullopt;
        }
    } else if (c2b1.status ==
               Sw92PhaseAssignedThreePhaseStatus::aqueous_phase_disappearance) {
        if (!evidence.w_disappeared) { return std::nullopt; }
    }

    StabilityPhase equivalent_reference;
    equivalent_reference.ln_phi.assign(n, 0.0);
    for (std::size_t i = 0; i < n; ++i) {
        if (c2b1.feed[i] > 0.0) {
            equivalent_reference.ln_phi[i] =
                evidence.common_log_activity[i] - std::log(c2b1.feed[i]);
        }
    }
    evidence.retained_h0_trivial_point = tangent_plane_distance(
        point.hydrocarbon0_phase.composition, c2b1.feed,
        selected_h0.activity, equivalent_reference);
    evidence.retained_h1_trivial_point = tangent_plane_distance(
        point.hydrocarbon1_phase.composition, c2b1.feed,
        selected_h1.activity, equivalent_reference);
    const double allowance0 = evidence.common_reference_allowance +
                              evidence.retained_h0_trivial_point.roundoff_guard;
    const double allowance1 = evidence.common_reference_allowance +
                              evidence.retained_h1_trivial_point.roundoff_guard;
    if (std::abs(evidence.retained_h0_trivial_point.value) > allowance0 ||
        std::abs(evidence.retained_h1_trivial_point.value) > allowance1) {
        return std::nullopt;
    }
    return evidence;
}

inline double sw92_phase_assigned_c2b2_effective_tpd_tolerance(
    double base, double allowance) {
    if (!std::isfinite(base) || base < 0.0 ||
        !std::isfinite(allowance) || allowance < 0.0 ||
        allowance > std::numeric_limits<double>::max() - base) {
        throw std::domain_error(
            "SW92 phase-assigned C2b.2: invalid TPD tolerance/reference allowance");
    }
    const double effective = base + allowance;
    if (!std::isfinite(effective)) {
        throw std::domain_error(
            "SW92 phase-assigned C2b.2: nonrepresentable effective TPD tolerance");
    }
    return effective;
}

inline void sw92_phase_assigned_c2b2_preflight_starts(
    std::span<const double> feed, const StabilityOptions& options,
    std::span<const std::vector<double>> extra_starts) {
    const std::size_t n = feed.size();
    if (n == 0U || n > options.max_components) {
        throw std::length_error(
            "SW92 phase-assigned C2b.2: component quota exceeded or empty feed");
    }
    std::size_t active = 0U;
    for (double value : feed) { if (value > 0.0) { ++active; } }
    constexpr std::size_t mandatory = 4U; // W diagnostic, H0, H1, C2a1 witness.
    const std::size_t generated = options.automatic_starts
        ? (active == 1U ? 1U : active + 2U) : 0U;
    if (extra_starts.size() > std::numeric_limits<std::size_t>::max() - mandatory) {
        throw std::length_error("SW92 phase-assigned C2b.2: start-count overflow");
    }
    const std::size_t supplied = mandatory + extra_starts.size();
    if (generated > options.max_starts || supplied > options.max_starts - generated) {
        throw std::length_error("SW92 phase-assigned C2b.2: start quota exceeded");
    }
    const std::size_t count = generated + supplied;
    if (count > options.max_start_entries / n ||
        count > std::vector<StabilityTrial>{}.max_size()) {
        throw std::length_error(
            "SW92 phase-assigned C2b.2: start storage quota exceeded");
    }
    for (const auto& start : extra_starts) {
        if (start.size() != n) {
            throw std::invalid_argument(
                "SW92 phase-assigned C2b.2: start dimension mismatch");
        }
        (void)stability_check_composition(start);
        for (std::size_t i = 0; i < n; ++i) {
            if ((feed[i] == 0.0 && start[i] != 0.0) ||
                (feed[i] > 0.0 && start[i] == 0.0)) {
                throw std::domain_error(
                    "SW92 phase-assigned C2b.2: starts must match active feed support");
            }
        }
    }
}

inline Sw92PhaseAssignedC2b2NaNegativeWitness
sw92_phase_assigned_c2b2_classify_na_witness(
    std::size_t trial_index, const TpdPoint& point,
    std::span<const double> h0, std::span<const double> h1,
    std::span<const double> feed, double retained_w_water,
    std::size_t water_index, double log_composition_separation) {
    Sw92PhaseAssignedC2b2NaNegativeWitness witness;
    witness.trial_index = trial_index;
    witness.point = point;
    if (point.composition.size() != h0.size() ||
        point.composition.size() != h1.size() ||
        point.composition.size() != feed.size() ||
        water_index >= point.composition.size()) {
        return witness;
    }
    witness.log_distance_from_h0 = sw92_phase_assigned_log_distance(
        point.composition, h0, feed);
    witness.log_distance_from_h1 = sw92_phase_assigned_log_distance(
        point.composition, h1, feed);
    witness.compositionally_distinct_from_h0 =
        std::isfinite(witness.log_distance_from_h0) &&
        witness.log_distance_from_h0 > log_composition_separation;
    witness.compositionally_distinct_from_h1 =
        std::isfinite(witness.log_distance_from_h1) &&
        witness.log_distance_from_h1 > log_composition_separation;
    const double trial_water = point.composition[water_index];
    witness.retained_w_minus_trial_water_fraction =
        retained_w_water - trial_water;
    witness.water_role_roundoff_guard = 256.0 * stability_eps *
        (1.0 + std::abs(retained_w_water) + std::abs(trial_water));
    witness.water_role_admissible =
        std::isfinite(witness.retained_w_minus_trial_water_fraction) &&
        std::isfinite(witness.water_role_roundoff_guard) &&
        witness.retained_w_minus_trial_water_fraction >
            witness.water_role_roundoff_guard;
    return witness;
}

inline bool sw92_phase_assigned_c2b2_search_complete(
    const StabilityResult& search) noexcept {
    return std::all_of(search.trials.begin(), search.trials.end(),
        [](const StabilityTrial& trial) {
            return trial.status == StabilityTrialStatus::stationary ||
                   trial.status == StabilityTrialStatus::negative_tpd;
        });
}

} // namespace detail

[[nodiscard]] inline Sw92PhaseAssignedC2b2Result
review_sw92_phase_assigned_c2b2(
    const Sw92PhaseAssignedJointResult& c1,
    const Sw92PhaseAssignedHSideWitnessResult& c2a1,
    const Sw92PhaseAssignedThreePhaseResult& c2b1,
    const thermodynamics::Sw92Phase<double>& model,
    Sw92PhaseAssignedC2b2Options options = {},
    std::span<const std::vector<double>> extra_starts = {}) {
    detail::sw92_phase_assigned_c2b2_check_options(options);
    if (!detail::sw92_phase_assigned_c2b2_model_matches(c1, model)) {
        throw std::invalid_argument(
            "SW92 phase-assigned C2b.2: C1 result/model snapshot mismatch");
    }

    Sw92PhaseAssignedC2b2Result result;
    result.options = options;
    result.pressure_pa = c1.pressure_pa;
    result.temperature_k = c1.temperature_k;
    result.feed = c1.feed;
    result.nacl_molality_mol_per_kg_water = c1.nacl_molality_mol_per_kg_water;
    result.dataset_id = c1.dataset_id;
    result.revision = c1.revision;
    result.component_ids = c1.component_ids;
    result.water_index = model.parameters().water_index();
    result.base_tpd_tolerance = options.stability.tpd_tolerance;

    if (c1.candidate() == nullptr ||
        c2a1.status != Sw92PhaseAssignedHSideWitnessStatus::
                           additional_nonaqueous_phase_witness_found ||
        !c2a1.nonaqueous_search || !c2b1.point ||
        !c2b1.source_witness_index) {
        result.status = Sw92PhaseAssignedC2b2Status::source_chain_unavailable;
        result.diagnostic =
            "Profile-C C2b.2 requires retained C1, C2a1 additional-NA and C2b.1 evidence";
        return result;
    }
    if (!detail::sw92_phase_assigned_c2b1_source_matches(c1, c2a1) ||
        !detail::sw92_phase_assigned_c2b2_three_phase_matches(c1, c2b1) ||
        !detail::sw92_phase_assigned_c2b2_source_seed_matches(c1, c2a1, c2b1)) {
        result.status = Sw92PhaseAssignedC2b2Status::source_chain_inconsistent;
        result.diagnostic =
            "Profile-C C2b.2 source metadata or C2a1->C2b.1 seed provenance is inconsistent";
        return result;
    }
    result.selected_c2a1_witness_index = c2b1.source_witness_index;
    const auto& source_witness =
        c2a1.negative_witnesses[*c2b1.source_witness_index];

    try {
        const auto c1_evidence = detail::sw92_phase_assigned_c2b2_revalidate_c1(
            c1, model, result.revalidation_property_evaluations);
        if (!c1_evidence) {
            result.status = Sw92PhaseAssignedC2b2Status::source_chain_inconsistent;
            result.diagnostic =
                "Profile-C C2b.2 C1 source failed fresh equilibrium/root/Gibbs revalidation";
            return result;
        }
        const auto* c1_point = c1.candidate();
        const auto reclassified_source = detail::sw92_phase_assigned_classify_na_witness(
            source_witness.trial_index, source_witness.point,
            c1_point->nonaqueous_phase.composition,
            c1_evidence->retained_w_water_fraction, result.water_index,
            c2a1.options.log_composition_separation);
        if (!detail::stability_negative(
                source_witness.point, c2a1.nonaqueous_search->options) ||
            !reclassified_source.usable_h_split_seed()) {
            result.status = Sw92PhaseAssignedC2b2Status::source_chain_inconsistent;
            result.diagnostic =
                "Profile-C C2b.2 selected C2a1 witness failed retained robustness/role guards";
            return result;
        }

        const auto c2b1_evidence =
            detail::sw92_phase_assigned_c2b2_revalidate_three_phase(
                c2b1, model, result.revalidation_property_evaluations);
        if (!c2b1_evidence) {
            result.status = Sw92PhaseAssignedC2b2Status::source_chain_inconsistent;
            result.diagnostic =
                "Profile-C C2b.2 C2b.1 source failed fresh three-phase equilibrium/topology revalidation";
            return result;
        }

        result.source_c1_reduced_gibbs = c1_evidence->reduced_gibbs;
        result.source_c1_gibbs_roundoff_guard = c1_evidence->gibbs_roundoff_guard;
        result.candidate_c2b1_reduced_gibbs = c2b1_evidence->reduced_gibbs;
        result.candidate_c2b1_gibbs_roundoff_guard =
            c2b1_evidence->gibbs_roundoff_guard;
        result.nested_gibbs_guard = c1_evidence->gibbs_roundoff_guard +
                                    c2b1_evidence->gibbs_roundoff_guard;
        result.common_log_activity = c2b1_evidence->common_log_activity;
        result.common_reference_allowance =
            c2b1_evidence->common_reference_allowance;
        result.retained_h0_trivial_point =
            c2b1_evidence->retained_h0_trivial_point;
        result.retained_h1_trivial_point =
            c2b1_evidence->retained_h1_trivial_point;

        if (c2b1_evidence->w_disappeared &&
            (c2b1_evidence->h0_disappeared || c2b1_evidence->h1_disappeared)) {
            result.status =
                Sw92PhaseAssignedC2b2Status::single_phase_endpoint_unresolved;
            result.diagnostic =
                "C2b.1 reached a one-phase endpoint; autonomous W-only versus H-only role selection remains unresolved";
            return result;
        }
        if (c2b1_evidence->w_disappeared) {
            result.status = Sw92PhaseAssignedC2b2Status::route_to_no_w_h0_h1;
            result.diagnostic =
                "C2b.1 W(AQ) disappeared; route to the dedicated no-W H0+H1 topology path without accepting it here";
            return result;
        }
        if (c2b1_evidence->h0_disappeared || c2b1_evidence->h1_disappeared) {
            result.status = Sw92PhaseAssignedC2b2Status::route_to_w_h;
            result.diagnostic =
                "C2b.1 one-H disappearance reached the neighboring W+H topology; route back to C1 review without auto-acceptance";
            return result;
        }

        if (result.candidate_c2b1_reduced_gibbs >
            result.source_c1_reduced_gibbs + result.nested_gibbs_guard) {
            result.status = Sw92PhaseAssignedC2b2Status::nested_gibbs_inconsistent;
            result.diagnostic =
                "C2b.1 phase addition has reduced Gibbs above its revalidated nested C1 source beyond arithmetic guard";
            return result;
        }

        detail::sw92_phase_assigned_c2b2_preflight_starts(
            c1.feed, options.stability, extra_starts);
        result.effective_tpd_tolerance =
            detail::sw92_phase_assigned_c2b2_effective_tpd_tolerance(
                result.base_tpd_tolerance,
                result.common_reference_allowance);

        const auto& state = *c2b1.point;
        std::vector<std::vector<double>> starts;
        starts.reserve(extra_starts.size() + 4U);
        starts.push_back(state.aqueous_phase.composition);
        starts.push_back(state.hydrocarbon0_phase.composition);
        starts.push_back(state.hydrocarbon1_phase.composition);
        starts.push_back(source_witness.point.composition);
        for (const auto& start : extra_starts) { starts.push_back(start); }

        StabilityOptions search_options = options.stability;
        search_options.tpd_tolerance = result.effective_tpd_tolerance;
        Sw92FamilyStabilityEvaluator evaluator(
            model, c1.nacl_molality_mol_per_kg_water,
            thermodynamics::SwPhaseFamily::nonaqueous,
            c2b1.options.nonaqueous_root_options);
        result.nonaqueous_search = test_pt_stability_against(
            c1.pressure_pa, c1.temperature_k, c1.feed,
            result.common_log_activity, evaluator, search_options, starts);

        const auto& search = *result.nonaqueous_search;
        for (std::size_t trial_index = 0; trial_index < search.trials.size();
             ++trial_index) {
            const auto& trial = search.trials[trial_index];
            if (trial.status != StabilityTrialStatus::negative_tpd || !trial.point ||
                !detail::stability_negative(*trial.point, search.options)) {
                continue;
            }
            result.negative_witnesses.push_back(
                detail::sw92_phase_assigned_c2b2_classify_na_witness(
                    trial_index, *trial.point,
                    state.hydrocarbon0_phase.composition,
                    state.hydrocarbon1_phase.composition,
                    c1.feed, state.aqueous_phase.composition[result.water_index],
                    result.water_index, options.log_composition_separation));
        }

        if (result.usable_additional_h_witness_count() != 0U) {
            result.status =
                Sw92PhaseAssignedC2b2Status::higher_phase_count_witness_found;
            result.diagnostic =
                "post-C2b.1 NA-only finite review found a robust additional-H witness distinct from both retained H phases; maximum-three-phase scope cannot accept this candidate";
            return result;
        }
        if (!detail::sw92_phase_assigned_c2b2_search_complete(search)) {
            result.status = Sw92PhaseAssignedC2b2Status::indeterminate;
            result.diagnostic =
                "post-C2b.1 NA-only finite review contains unresolved numerical/property trials";
            return result;
        }

        result.status =
            Sw92PhaseAssignedC2b2Status::w_present_h_multiplicity_locally_closed;
        result.diagnostic =
            "W-present H multiplicity is locally closed under the declared finite NA-only review; this is not global phase-number proof, morphology resolution, or accepted phase-set publication";
        return result;
    } catch (const StabilityPropertyError& error) {
        result.status = Sw92PhaseAssignedC2b2Status::indeterminate;
        result.property_issue = error.issue();
        result.diagnostic = error.what();
        return result;
    }
}

} // namespace mpmc::flash

#endif // MPMC_FLASH_SW92_PHASE_ASSIGNED_THREE_PHASE_CLOSURE_HPP
