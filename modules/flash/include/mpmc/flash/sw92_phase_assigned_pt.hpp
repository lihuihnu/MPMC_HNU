#ifndef MPMC_FLASH_SW92_PHASE_ASSIGNED_PT_HPP
#define MPMC_FLASH_SW92_PHASE_ASSIGNED_PT_HPP

#include <mpmc/flash/sw92_phase_assigned_no_w.hpp>
#include <mpmc/flash/sw92_phase_assigned_three_phase_closure.hpp>

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

inline constexpr std::string_view sw92_phase_assigned_pt_convention =
    "SW92-equilibrium/phase-assigned-aq-na-joint/topology-orchestration-pt/v1";

struct Sw92PhaseAssignedPtOptions {
    Sw92PhaseAssignedNoWOptions no_w;
    Sw92PhaseAssignedJointOptions c1;
    Sw92PhaseAssignedHSideWitnessOptions c2a1;
    Sw92PhaseAssignedThreePhaseOptions c2b1;
    Sw92PhaseAssignedC2b2Options c2b2;
    std::size_t max_c1_attempts{8};
    std::size_t max_c2b1_attempts{8};
};

enum class Sw92PhaseAssignedPtStatus {
    no_w_single_h_locally_closed,
    no_w_two_h_locally_closed,
    w_h_locally_closed,
    w_h0_h1_locally_closed,
    higher_phase_count_or_wrong_candidate,
    topology_unresolved,
    numerical_indeterminate
};

enum class Sw92PhaseAssignedPtPhysicalRole {
    aqueous,
    nonaqueous_unclassified
};

struct Sw92PhaseAssignedPtPhase {
    Sw92PhaseAssignedPtPhysicalRole physical_role{
        Sw92PhaseAssignedPtPhysicalRole::nonaqueous_unclassified};
    thermodynamics::SwPhaseFamily thermodynamic_family{
        thermodynamics::SwPhaseFamily::nonaqueous};
    double mole_phase_fraction{};
    std::vector<double> composition;
    std::optional<double> compressibility_factor;
};

struct Sw92PhaseAssignedPtResult {
    static constexpr std::size_t maximum_phase_count = 3;
    static constexpr bool global_stability_proven = false;
    static constexpr bool accepted_phase_set_published = false;
    static constexpr bool morphology_resolved = false;

    Sw92PhaseAssignedPtStatus status{
        Sw92PhaseAssignedPtStatus::numerical_indeterminate};
    Sw92PhaseAssignedPtOptions options;
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
    std::string equilibrium_profile{sw92_phase_assigned_aq_na_joint_profile};
    std::string orchestration_convention{sw92_phase_assigned_pt_convention};

    Sw92PhaseAssignedNoWResult no_w;
    std::optional<Sw92PhaseAssignedJointResult> c1;
    std::optional<Sw92PhaseAssignedHSideWitnessResult> c2a1;
    std::optional<Sw92PhaseAssignedThreePhaseResult> c2b1;
    std::optional<Sw92PhaseAssignedC2b2Result> c2b2;

    std::vector<Sw92PhaseAssignedPtPhase> phases;
    std::size_t c1_attempts{};
    std::size_t c2b1_attempts{};
    bool c1_attempt_limit_reached{false};
    bool c2b1_attempt_limit_reached{false};
    bool unstable_no_w_w_search_attempted{false};
    std::string diagnostic;

    [[nodiscard]] bool locally_closed_phase_candidate() const noexcept {
        switch (status) {
        case Sw92PhaseAssignedPtStatus::no_w_single_h_locally_closed:
        case Sw92PhaseAssignedPtStatus::no_w_two_h_locally_closed:
        case Sw92PhaseAssignedPtStatus::w_h_locally_closed:
        case Sw92PhaseAssignedPtStatus::w_h0_h1_locally_closed:
            return !phases.empty();
        case Sw92PhaseAssignedPtStatus::higher_phase_count_or_wrong_candidate:
        case Sw92PhaseAssignedPtStatus::topology_unresolved:
        case Sw92PhaseAssignedPtStatus::numerical_indeterminate:
            return false;
        }
        return false;
    }
};

namespace detail {

inline void sw92_phase_assigned_pt_check_options(
    const Sw92PhaseAssignedPtOptions& options) {
    if (options.max_c1_attempts == 0U || options.max_c2b1_attempts == 0U) {
        throw std::invalid_argument(
            "SW92 Profile-C PT: positive candidate-attempt quotas required");
    }
}

inline double sw92_phase_assigned_pt_log_distance(
    std::span<const double> first, std::span<const double> second,
    std::span<const double> feed) {
    return sw92_phase_assigned_no_w_log_distance(first, second, feed);
}

inline bool sw92_phase_assigned_pt_water_role_seed_ok(
    std::span<const double> w, std::span<const double> h,
    std::size_t water_index) {
    if (water_index >= w.size() || w.size() != h.size()) { return false; }
    const double guard = 256.0 * stability_eps *
        (1.0 + std::abs(w[water_index]) + std::abs(h[water_index]));
    return std::isfinite(guard) && w[water_index] - h[water_index] > guard;
}

inline std::optional<std::vector<double>> sw92_phase_assigned_pt_log_k(
    std::span<const double> w, std::span<const double> h,
    std::span<const double> feed) {
    if (w.size() != h.size() || w.size() != feed.size()) { return std::nullopt; }
    std::vector<double> log_k(feed.size(), 0.0);
    for (std::size_t i = 0; i < feed.size(); ++i) {
        if (feed[i] == 0.0) {
            if (w[i] != 0.0 || h[i] != 0.0) { return std::nullopt; }
            continue;
        }
        if (!(w[i] > 0.0) || !(h[i] > 0.0)) { return std::nullopt; }
        log_k[i] = std::log(h[i]) - std::log(w[i]);
        if (!std::isfinite(log_k[i])) { return std::nullopt; }
    }
    return log_k;
}

inline std::vector<std::vector<double>> sw92_phase_assigned_pt_h_side_starts(
    const Sw92PhaseAssignedNoWResult& no_w,
    std::span<const double> retained_h,
    double separation) {
    std::vector<std::vector<double>> starts;
    for (const auto& candidate : no_w.retained_na_candidate_compositions) {
        const double distance = sw92_phase_assigned_pt_log_distance(
            candidate, retained_h, no_w.feed);
        if (std::isfinite(distance) && distance > separation) {
            starts.push_back(candidate);
        }
    }
    return starts;
}

/// A fixed-NA two-phase candidate may itself be unstable to another NA phase.
/// That is NOT sufficient to terminate Profile-C: the physical solution may be
/// W+H0+H1. Reconstruct the selected fixed-NA pair and perform the same targeted
/// AQ appearance search used by the no-W adapter. If a W seed is found, the
/// top-level graph continues into joint Profile-C equations.
inline bool sw92_phase_assigned_pt_try_w_from_unstable_no_w(
    Sw92PhaseAssignedNoWResult& no_w,
    const thermodynamics::Sw92Phase<double>& model) {
    auto& solution = no_w.hydrocarbon_flash.solution;
    if (no_w.status !=
            Sw92PhaseAssignedNoWStatus::higher_h_multiplicity_or_wrong_candidate ||
        solution.status != PtSplitStatus::phase_set_unstable) {
        return false;
    }
    const auto* candidate = solution.candidate();
    if (candidate == nullptr || no_w.feed.empty() ||
        no_w.water_index >= no_w.feed.size()) {
        return false;
    }

    no_w.retained_na_candidate_compositions.clear();
    no_w.retained_na_candidate_fractions.clear();
    no_w.aqueous_witnesses.clear();
    no_w.selected_water_witness_index.reset();
    no_w.retained_na_candidate_compositions.push_back(candidate->fractions.liquid);
    no_w.retained_na_candidate_compositions.push_back(candidate->fractions.vapor);
    no_w.retained_na_candidate_fractions.push_back(
        1.0 - candidate->fractions.vapor_fraction);
    no_w.retained_na_candidate_fractions.push_back(candidate->fractions.vapor_fraction);
    no_w.common_log_activity = candidate->common_log_activity;

    no_w.min_retained_na_candidate_water_fraction =
        std::numeric_limits<double>::infinity();
    no_w.max_retained_na_candidate_water_fraction = 0.0;
    for (const auto& phase : no_w.retained_na_candidate_compositions) {
        if (phase.size() != no_w.feed.size()) { return false; }
        const double water = phase[no_w.water_index];
        no_w.min_retained_na_candidate_water_fraction = std::min(
            no_w.min_retained_na_candidate_water_fraction, water);
        no_w.max_retained_na_candidate_water_fraction = std::max(
            no_w.max_retained_na_candidate_water_fraction, water);
    }
    if (!std::isfinite(no_w.min_retained_na_candidate_water_fraction) ||
        !std::isfinite(no_w.max_retained_na_candidate_water_fraction) ||
        no_w.feed[no_w.water_index] == 0.0) {
        return false;
    }

    no_w.targeted_water_start = sw92_phase_assigned_targeted_water_start(
        no_w.feed, no_w.water_index,
        no_w.max_retained_na_candidate_water_fraction);
    if (!no_w.targeted_water_start) { return false; }

    const std::size_t mandatory = 1U + no_w.retained_na_candidate_count();
    sw92_phase_assigned_no_w_preflight_starts(
        no_w.feed, no_w.options.aqueous_appearance, mandatory, {});
    std::vector<std::vector<double>> starts;
    starts.reserve(mandatory);
    starts.push_back(*no_w.targeted_water_start);
    for (const auto& phase : no_w.retained_na_candidate_compositions) {
        starts.push_back(phase);
    }

    Sw92FamilyStabilityEvaluator aqueous(
        model, no_w.nacl_molality_mol_per_kg_water,
        thermodynamics::SwPhaseFamily::aqueous,
        no_w.options.aqueous_root_options);
    no_w.aqueous_appearance_search = test_pt_stability_against(
        no_w.pressure_pa, no_w.temperature_k, no_w.feed,
        no_w.common_log_activity, aqueous,
        no_w.options.aqueous_appearance, starts);
    const auto& search = *no_w.aqueous_appearance_search;
    if (search.trials.empty()) {
        no_w.status = Sw92PhaseAssignedNoWStatus::indeterminate;
        no_w.diagnostic =
            "unstable fixed-NA pair required W-rival review but the targeted AQ search produced no trial";
        return false;
    }

    const auto& targeted = search.trials.front();
    if (targeted.status == StabilityTrialStatus::negative_tpd && targeted.point &&
        stability_negative(*targeted.point, search.options)) {
        no_w.aqueous_witnesses.push_back(sw92_phase_assigned_classify_water_witness(
            0U, *targeted.point, no_w.retained_na_candidate_compositions,
            no_w.feed, no_w.water_index,
            no_w.min_retained_na_candidate_water_fraction,
            no_w.options.log_composition_separation));
        if (no_w.aqueous_witnesses.front().usable_w_seed()) {
            no_w.selected_water_witness_index = 0U;
            no_w.status = Sw92PhaseAssignedNoWStatus::aqueous_phase_witness_found;
            no_w.diagnostic =
                "fixed-NA pair has further NA instability and a targeted physical W witness; continue into joint Profile-C topology rather than terminating at all-NA multiplicity";
            return true;
        }
    }

    if (targeted.status != StabilityTrialStatus::stationary &&
        targeted.status != StabilityTrialStatus::negative_tpd) {
        no_w.status = Sw92PhaseAssignedNoWStatus::indeterminate;
        no_w.diagnostic =
            "unstable fixed-NA pair required W-rival review but targeted AQ search was numerically/property indeterminate";
    }
    return false;
}

inline void sw92_phase_assigned_pt_publish_no_w(
    Sw92PhaseAssignedPtResult& result) {
    const auto& no_w = result.no_w;
    const auto& solution = no_w.hydrocarbon_flash.solution;
    if (no_w.retained_na_candidate_compositions.empty() ||
        no_w.retained_na_candidate_compositions.size() !=
            no_w.retained_na_candidate_fractions.size()) {
        return;
    }

    const PtSplitState* pair = solution.candidate();
    for (std::size_t i = 0; i < no_w.retained_na_candidate_compositions.size(); ++i) {
        Sw92PhaseAssignedPtPhase phase;
        phase.physical_role = Sw92PhaseAssignedPtPhysicalRole::nonaqueous_unclassified;
        phase.thermodynamic_family = thermodynamics::SwPhaseFamily::nonaqueous;
        phase.mole_phase_fraction = no_w.retained_na_candidate_fractions[i];
        phase.composition = no_w.retained_na_candidate_compositions[i];
        if (pair != nullptr && no_w.retained_na_candidate_compositions.size() == 2U) {
            phase.compressibility_factor = i == 0U ? pair->liquid.z : pair->vapor.z;
        }
        result.phases.push_back(std::move(phase));
    }
}

inline void sw92_phase_assigned_pt_publish_c1(
    Sw92PhaseAssignedPtResult& result,
    const Sw92PhaseAssignedJointResult& c1) {
    const auto* point = c1.candidate();
    if (point == nullptr) { return; }
    Sw92PhaseAssignedPtPhase w;
    w.physical_role = Sw92PhaseAssignedPtPhysicalRole::aqueous;
    w.thermodynamic_family = thermodynamics::SwPhaseFamily::aqueous;
    w.mole_phase_fraction = point->aqueous_phase.mole_phase_fraction;
    w.composition = point->aqueous_phase.composition;
    w.compressibility_factor = point->aqueous_phase.compressibility_factor;
    result.phases.push_back(std::move(w));

    Sw92PhaseAssignedPtPhase h;
    h.physical_role = Sw92PhaseAssignedPtPhysicalRole::nonaqueous_unclassified;
    h.thermodynamic_family = thermodynamics::SwPhaseFamily::nonaqueous;
    h.mole_phase_fraction = point->nonaqueous_phase.mole_phase_fraction;
    h.composition = point->nonaqueous_phase.composition;
    h.compressibility_factor = point->nonaqueous_phase.compressibility_factor;
    result.phases.push_back(std::move(h));
}

inline void sw92_phase_assigned_pt_publish_c2b1(
    Sw92PhaseAssignedPtResult& result,
    const Sw92PhaseAssignedThreePhaseResult& c2b1) {
    const auto* point = c2b1.candidate();
    if (point == nullptr) { return; }
    const auto add = [&](const Sw92PhaseAssignedThreePhasePhase& source,
                         Sw92PhaseAssignedPtPhysicalRole role) {
        Sw92PhaseAssignedPtPhase phase;
        phase.physical_role = role;
        phase.thermodynamic_family = source.thermodynamic_family;
        phase.mole_phase_fraction = source.mole_phase_fraction;
        phase.composition = source.composition;
        phase.compressibility_factor = source.compressibility_factor;
        result.phases.push_back(std::move(phase));
    };
    add(point->aqueous_phase, Sw92PhaseAssignedPtPhysicalRole::aqueous);
    add(point->hydrocarbon0_phase,
        Sw92PhaseAssignedPtPhysicalRole::nonaqueous_unclassified);
    add(point->hydrocarbon1_phase,
        Sw92PhaseAssignedPtPhysicalRole::nonaqueous_unclassified);
}

} // namespace detail

[[nodiscard]] inline Sw92PhaseAssignedPtResult solve_sw92_phase_assigned_pt(
    double pressure_pa, double temperature_k, std::span<const double> feed,
    const thermodynamics::Sw92Phase<double>& model,
    double nacl_molality_mol_per_kg_water,
    Sw92PhaseAssignedPtOptions options = {}) {
    detail::sw92_phase_assigned_pt_check_options(options);

    Sw92PhaseAssignedPtResult result;
    result.options = options;
    result.pressure_pa = pressure_pa;
    result.temperature_k = temperature_k;
    result.input_feed_sum = detail::stability_check_composition(feed);
    result.feed = detail::stability_normalize(feed, result.input_feed_sum);
    result.nacl_molality_mol_per_kg_water = nacl_molality_mol_per_kg_water;
    const auto& parameters = model.parameters();
    result.dataset_id = parameters.dataset_id();
    result.revision = parameters.revision();
    for (const auto& component : parameters.components().items()) {
        result.component_ids.push_back(component.id);
    }

    result.no_w = solve_sw92_phase_assigned_no_w(
        pressure_pa, temperature_k, feed, model,
        nacl_molality_mol_per_kg_water, options.no_w);

    if (result.no_w.status ==
        Sw92PhaseAssignedNoWStatus::higher_h_multiplicity_or_wrong_candidate) {
        result.unstable_no_w_w_search_attempted = true;
        (void)detail::sw92_phase_assigned_pt_try_w_from_unstable_no_w(
            result.no_w, model);
    }

    switch (result.no_w.status) {
    case Sw92PhaseAssignedNoWStatus::no_w_single_h_locally_closed:
        result.status = Sw92PhaseAssignedPtStatus::no_w_single_h_locally_closed;
        detail::sw92_phase_assigned_pt_publish_no_w(result);
        result.diagnostic =
            "Profile-C finite topology graph locally closed at one no-W NA/H candidate; morphology remains unresolved";
        return result;
    case Sw92PhaseAssignedNoWStatus::no_w_two_h_locally_closed:
        result.status = Sw92PhaseAssignedPtStatus::no_w_two_h_locally_closed;
        detail::sw92_phase_assigned_pt_publish_no_w(result);
        result.diagnostic =
            "Profile-C finite topology graph locally closed at two no-W NA candidates; morphology remains unresolved";
        return result;
    case Sw92PhaseAssignedNoWStatus::higher_h_multiplicity_or_wrong_candidate:
        result.status = Sw92PhaseAssignedPtStatus::higher_phase_count_or_wrong_candidate;
        result.diagnostic = result.no_w.diagnostic;
        return result;
    case Sw92PhaseAssignedNoWStatus::phase_disappearance_unresolved:
    case Sw92PhaseAssignedNoWStatus::single_phase_role_unresolved:
        result.status = Sw92PhaseAssignedPtStatus::topology_unresolved;
        result.diagnostic = result.no_w.diagnostic;
        return result;
    case Sw92PhaseAssignedNoWStatus::indeterminate:
        result.status = Sw92PhaseAssignedPtStatus::numerical_indeterminate;
        result.diagnostic = result.no_w.diagnostic;
        return result;
    case Sw92PhaseAssignedNoWStatus::aqueous_phase_witness_found:
        break;
    }

    const auto* w_witness = result.no_w.selected_water_witness();
    if (w_witness == nullptr) {
        result.status = Sw92PhaseAssignedPtStatus::topology_unresolved;
        result.diagnostic =
            "no-W stage reported W appearance without a usable retained witness";
        return result;
    }
    const auto& w_seed = w_witness->point.composition;

    // Do NOT compare all-NA no-W Gibbs against W(AQ)+H(NA) Gibbs. Those are
    // different physical-role/model assignments in an asymmetric SW92 profile;
    // such a cross-topology ranking would recreate Profile-B model competition.
    // A targeted W instability generates the rival physical topology, and the
    // joint C1 equations plus role guard decide whether that topology exists.
    std::optional<Sw92PhaseAssignedJointResult> best_c1;
    double best_c1_gibbs = std::numeric_limits<double>::infinity();
    for (const auto& h_seed : result.no_w.retained_na_candidate_compositions) {
        if (result.c1_attempts >= options.max_c1_attempts) {
            result.c1_attempt_limit_reached = true;
            break;
        }
        if (!detail::sw92_phase_assigned_pt_water_role_seed_ok(
                w_seed, h_seed, parameters.water_index())) {
            continue;
        }
        const auto log_k = detail::sw92_phase_assigned_pt_log_k(
            w_seed, h_seed, result.feed);
        if (!log_k) { continue; }
        ++result.c1_attempts;
        auto candidate = iterate_sw92_phase_assigned_aq_na_joint(
            pressure_pa, temperature_k, result.feed, *log_k, model,
            nacl_molality_mol_per_kg_water, options.c1);
        const auto* point = candidate.candidate();
        if (point == nullptr) { continue; }
        if (point->reduced_gibbs < best_c1_gibbs) {
            best_c1_gibbs = point->reduced_gibbs;
            best_c1 = std::move(candidate);
        }
    }
    if (!best_c1) {
        result.status = Sw92PhaseAssignedPtStatus::topology_unresolved;
        result.diagnostic = result.c1_attempt_limit_reached
            ? "W witness exists but the C1 attempt quota was exhausted before an admissible W+H candidate was retained"
            : "W witness exists but no admissible W+H joint candidate was found";
        return result;
    }
    result.c1 = std::move(best_c1);

    const auto& retained_h = result.c1->candidate()->nonaqueous_phase.composition;
    const auto h_starts = detail::sw92_phase_assigned_pt_h_side_starts(
        result.no_w, retained_h, options.c2a1.log_composition_separation);
    result.c2a1 = test_sw92_phase_assigned_h_side_na_witness(
        *result.c1, model, options.c2a1, h_starts);

    switch (result.c2a1->status) {
    case Sw92PhaseAssignedHSideWitnessStatus::no_additional_nonaqueous_witness_found:
        result.status = Sw92PhaseAssignedPtStatus::w_h_locally_closed;
        detail::sw92_phase_assigned_pt_publish_c1(result, *result.c1);
        result.diagnostic =
            "Profile-C finite topology graph locally closed at W(AQ)+H(NA); H morphology remains unresolved";
        return result;
    case Sw92PhaseAssignedHSideWitnessStatus::source_candidate_unavailable:
    case Sw92PhaseAssignedHSideWitnessStatus::source_candidate_inconsistent:
        result.status = Sw92PhaseAssignedPtStatus::topology_unresolved;
        result.diagnostic = result.c2a1->diagnostic;
        return result;
    case Sw92PhaseAssignedHSideWitnessStatus::indeterminate:
        result.status = Sw92PhaseAssignedPtStatus::numerical_indeterminate;
        result.diagnostic = result.c2a1->diagnostic;
        return result;
    case Sw92PhaseAssignedHSideWitnessStatus::additional_nonaqueous_phase_witness_found:
        break;
    }

    std::optional<Sw92PhaseAssignedThreePhaseResult> best_c2b1;
    double best_c2b1_gibbs = std::numeric_limits<double>::infinity();
    for (std::size_t witness_index = 0;
         witness_index < result.c2a1->negative_witnesses.size(); ++witness_index) {
        const auto& witness = result.c2a1->negative_witnesses[witness_index];
        if (!witness.usable_h_split_seed()) { continue; }
        if (result.c2b1_attempts >= options.max_c2b1_attempts) {
            result.c2b1_attempt_limit_reached = true;
            break;
        }
        ++result.c2b1_attempts;
        auto candidate = solve_sw92_phase_assigned_c2b1_candidate(
            *result.c1, *result.c2a1, witness_index, model, options.c2b1);
        const auto* point = candidate.candidate();
        if (point == nullptr) { continue; }
        const double c1_guard = result.c1->candidate()->gibbs_roundoff_guard;
        if (point->reduced_gibbs > best_c1_gibbs +
                c1_guard + point->gibbs_roundoff_guard) {
            continue;
        }
        if (point->reduced_gibbs < best_c2b1_gibbs) {
            best_c2b1_gibbs = point->reduced_gibbs;
            best_c2b1 = std::move(candidate);
        }
    }
    if (!best_c2b1) {
        result.status = Sw92PhaseAssignedPtStatus::topology_unresolved;
        result.diagnostic = result.c2b1_attempt_limit_reached
            ? "additional-H evidence exists but the C2b.1 attempt quota was exhausted before a Gibbs-consistent three-phase candidate was retained"
            : "additional-H evidence exists but no Gibbs-consistent admissible W+H0+H1 candidate was found";
        return result;
    }
    result.c2b1 = std::move(best_c2b1);
    result.c2b2 = review_sw92_phase_assigned_c2b2(
        *result.c1, *result.c2a1, *result.c2b1, model, options.c2b2);

    switch (result.c2b2->status) {
    case Sw92PhaseAssignedC2b2Status::w_present_h_multiplicity_locally_closed:
        result.status = Sw92PhaseAssignedPtStatus::w_h0_h1_locally_closed;
        detail::sw92_phase_assigned_pt_publish_c2b1(result, *result.c2b1);
        result.diagnostic =
            "Profile-C finite topology graph locally closed at W(AQ)+H0(NA)+H1(NA); H morphology remains unresolved";
        return result;
    case Sw92PhaseAssignedC2b2Status::higher_phase_count_witness_found:
        result.status = Sw92PhaseAssignedPtStatus::higher_phase_count_or_wrong_candidate;
        result.diagnostic = result.c2b2->diagnostic;
        return result;
    case Sw92PhaseAssignedC2b2Status::indeterminate:
        result.status = Sw92PhaseAssignedPtStatus::numerical_indeterminate;
        result.diagnostic = result.c2b2->diagnostic;
        return result;
    case Sw92PhaseAssignedC2b2Status::source_chain_unavailable:
    case Sw92PhaseAssignedC2b2Status::source_chain_inconsistent:
    case Sw92PhaseAssignedC2b2Status::route_to_w_h:
    case Sw92PhaseAssignedC2b2Status::route_to_no_w_h0_h1:
    case Sw92PhaseAssignedC2b2Status::single_phase_endpoint_unresolved:
    case Sw92PhaseAssignedC2b2Status::nested_gibbs_inconsistent:
        result.status = Sw92PhaseAssignedPtStatus::topology_unresolved;
        result.diagnostic = result.c2b2->diagnostic;
        return result;
    }
    result.status = Sw92PhaseAssignedPtStatus::numerical_indeterminate;
    result.diagnostic = "unhandled Profile-C topology state";
    return result;
}

} // namespace mpmc::flash

#endif // MPMC_FLASH_SW92_PHASE_ASSIGNED_PT_HPP
