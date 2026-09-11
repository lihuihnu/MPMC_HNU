#ifndef MPMC_FLASH_SW92_PHASE_ASSIGNED_BOUNDARY_HPP
#define MPMC_FLASH_SW92_PHASE_ASSIGNED_BOUNDARY_HPP

#include <mpmc/flash/sw92_phase_assigned_pt.hpp>

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

inline constexpr std::string_view sw92_phase_assigned_boundary_convention =
    "SW92-equilibrium/phase-assigned-aq-na-joint/"
    "three-phase-disappearance-neighbor-resolve/v1";

enum class Sw92PhaseAssignedBoundaryStatus {
    no_boundary_route,
    resolved_to_w_h,
    resolved_to_no_w_single_h,
    resolved_to_no_w_two_h,
    neighbor_topology_not_closed,
    higher_phase_count_or_wrong_candidate,
    source_chain_inconsistent,
    numerical_indeterminate
};

struct Sw92PhaseAssignedBoundaryResult {
    static constexpr bool global_stability_proven = false;
    static constexpr bool accepted_phase_set_published = false;
    static constexpr bool morphology_resolved = false;

    Sw92PhaseAssignedBoundaryStatus status{
        Sw92PhaseAssignedBoundaryStatus::no_boundary_route};
    bool re_solve_attempted{false};
    std::string convention{sw92_phase_assigned_boundary_convention};

    std::optional<Sw92PhaseAssignedNoWResult> neighbor_no_w;
    std::optional<Sw92PhaseAssignedJointResult> neighbor_c1;
    std::optional<Sw92PhaseAssignedHSideWitnessResult> neighbor_c2a1;
    std::vector<Sw92PhaseAssignedPtPhase> phases;
    std::string diagnostic;

    [[nodiscard]] bool neighbor_locally_closed() const noexcept {
        return status == Sw92PhaseAssignedBoundaryStatus::resolved_to_w_h ||
               status == Sw92PhaseAssignedBoundaryStatus::resolved_to_no_w_single_h ||
               status == Sw92PhaseAssignedBoundaryStatus::resolved_to_no_w_two_h;
    }
};

struct Sw92PhaseAssignedBoundaryAwareResult {
    static constexpr bool global_stability_proven = false;
    static constexpr bool accepted_phase_set_published = false;
    static constexpr bool morphology_resolved = false;

    Sw92PhaseAssignedPtResult base;
    std::optional<Sw92PhaseAssignedThreePhaseResult> boundary_c2b1;
    std::optional<Sw92PhaseAssignedC2b2Result> boundary_c2b2;
    std::optional<Sw92PhaseAssignedBoundaryResult> boundary;
    Sw92PhaseAssignedPtStatus status{
        Sw92PhaseAssignedPtStatus::numerical_indeterminate};
    std::vector<Sw92PhaseAssignedPtPhase> phases;
    std::size_t disappearance_attempts{};
    bool disappearance_attempt_limit_reached{false};
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

inline void sw92_phase_assigned_boundary_publish_c1(
    Sw92PhaseAssignedBoundaryResult& result,
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

inline void sw92_phase_assigned_boundary_publish_no_w(
    Sw92PhaseAssignedBoundaryResult& result,
    const Sw92PhaseAssignedNoWResult& no_w) {
    if (no_w.retained_na_candidate_compositions.empty() ||
        no_w.retained_na_candidate_compositions.size() !=
            no_w.retained_na_candidate_fractions.size()) {
        return;
    }
    const auto* pair = no_w.hydrocarbon_flash.solution.candidate();
    for (std::size_t i = 0; i < no_w.retained_na_candidate_compositions.size(); ++i) {
        Sw92PhaseAssignedPtPhase phase;
        phase.physical_role =
            Sw92PhaseAssignedPtPhysicalRole::nonaqueous_unclassified;
        phase.thermodynamic_family = thermodynamics::SwPhaseFamily::nonaqueous;
        phase.mole_phase_fraction = no_w.retained_na_candidate_fractions[i];
        phase.composition = no_w.retained_na_candidate_compositions[i];
        if (pair != nullptr && no_w.retained_na_candidate_compositions.size() == 2U) {
            phase.compressibility_factor = i == 0U ? pair->liquid.z : pair->vapor.z;
        }
        result.phases.push_back(std::move(phase));
    }
}

inline std::vector<std::vector<double>> sw92_phase_assigned_boundary_unique_starts(
    std::span<const std::vector<double>> source,
    std::span<const double> retained,
    std::span<const double> feed,
    double separation) {
    std::vector<std::vector<double>> starts;
    for (const auto& candidate : source) {
        const double distance = sw92_phase_assigned_pt_log_distance(
            candidate, retained, feed);
        if (std::isfinite(distance) && distance > separation) {
            starts.push_back(candidate);
        }
    }
    return starts;
}

inline Sw92PhaseAssignedBoundaryResult sw92_phase_assigned_resolve_w_h_neighbor(
    double pressure_pa, double temperature_k, std::span<const double> feed,
    std::span<const double> w_seed, std::span<const double> h_seed,
    const thermodynamics::Sw92Phase<double>& model,
    double nacl_molality_mol_per_kg_water,
    const Sw92PhaseAssignedPtOptions& options,
    std::span<const std::vector<double>> h_side_extra_starts = {}) {
    Sw92PhaseAssignedBoundaryResult result;
    result.re_solve_attempted = true;

    const auto log_k = sw92_phase_assigned_pt_log_k(w_seed, h_seed, feed);
    if (!log_k) {
        result.status = Sw92PhaseAssignedBoundaryStatus::source_chain_inconsistent;
        result.diagnostic =
            "W+H boundary re-solve could not construct a finite positive logK seed";
        return result;
    }

    result.neighbor_c1 = iterate_sw92_phase_assigned_aq_na_joint(
        pressure_pa, temperature_k, feed, *log_k, model,
        nacl_molality_mol_per_kg_water, options.c1);
    if (!result.neighbor_c1->candidate_admissible()) {
        result.status = Sw92PhaseAssignedBoundaryStatus::neighbor_topology_not_closed;
        result.diagnostic =
            "W+H boundary re-solve did not produce an admissible joint C1 state; the disappearing H cannot be silently dropped";
        return result;
    }

    result.neighbor_c2a1 = test_sw92_phase_assigned_h_side_na_witness(
        *result.neighbor_c1, model, options.c2a1, h_side_extra_starts);
    switch (result.neighbor_c2a1->status) {
    case Sw92PhaseAssignedHSideWitnessStatus::no_additional_nonaqueous_witness_found:
        result.status = Sw92PhaseAssignedBoundaryStatus::resolved_to_w_h;
        sw92_phase_assigned_boundary_publish_c1(result, *result.neighbor_c1);
        if (result.phases.size() != 2U) {
            result.status = Sw92PhaseAssignedBoundaryStatus::source_chain_inconsistent;
            result.phases.clear();
            result.diagnostic =
                "W+H boundary re-solve converged but phase publication data are incomplete";
            return result;
        }
        result.diagnostic =
            "three-phase one-H disappearance was re-solved as W(AQ)+H(NA), and the surviving two-phase state has no robust additional-H witness under the declared finite review";
        return result;
    case Sw92PhaseAssignedHSideWitnessStatus::additional_nonaqueous_phase_witness_found:
        result.status = Sw92PhaseAssignedBoundaryStatus::neighbor_topology_not_closed;
        result.diagnostic =
            "re-solved W+H neighbor still has a robust additional-H witness; the nominally disappearing H remains required and no phase is dropped";
        return result;
    case Sw92PhaseAssignedHSideWitnessStatus::source_candidate_unavailable:
    case Sw92PhaseAssignedHSideWitnessStatus::source_candidate_inconsistent:
        result.status = Sw92PhaseAssignedBoundaryStatus::source_chain_inconsistent;
        result.diagnostic = result.neighbor_c2a1->diagnostic;
        return result;
    case Sw92PhaseAssignedHSideWitnessStatus::indeterminate:
        result.status = Sw92PhaseAssignedBoundaryStatus::numerical_indeterminate;
        result.diagnostic = result.neighbor_c2a1->diagnostic;
        return result;
    }
    result.status = Sw92PhaseAssignedBoundaryStatus::numerical_indeterminate;
    result.diagnostic = "unhandled W+H boundary stability state";
    return result;
}

inline Sw92PhaseAssignedBoundaryResult sw92_phase_assigned_resolve_no_w_neighbor(
    double pressure_pa, double temperature_k, std::span<const double> feed,
    std::span<const std::vector<double>> retained_h_starts,
    const thermodynamics::Sw92Phase<double>& model,
    double nacl_molality_mol_per_kg_water,
    const Sw92PhaseAssignedPtOptions& options) {
    Sw92PhaseAssignedBoundaryResult result;
    result.re_solve_attempted = true;
    if (retained_h_starts.empty()) {
        result.status = Sw92PhaseAssignedBoundaryStatus::source_chain_inconsistent;
        result.diagnostic =
            "no-W boundary re-solve requires at least one surviving NA phase composition";
        return result;
    }

    result.neighbor_no_w = solve_sw92_phase_assigned_no_w(
        pressure_pa, temperature_k, feed, model,
        nacl_molality_mol_per_kg_water, options.no_w,
        retained_h_starts, retained_h_starts);

    switch (result.neighbor_no_w->status) {
    case Sw92PhaseAssignedNoWStatus::no_w_single_h_locally_closed:
        result.status = Sw92PhaseAssignedBoundaryStatus::resolved_to_no_w_single_h;
        sw92_phase_assigned_boundary_publish_no_w(result, *result.neighbor_no_w);
        if (result.phases.size() != 1U) {
            result.status = Sw92PhaseAssignedBoundaryStatus::source_chain_inconsistent;
            result.phases.clear();
            result.diagnostic =
                "no-W single-H neighbor closed but phase publication data are incomplete";
            return result;
        }
        result.diagnostic =
            "three-phase disappearance was re-solved as a no-W single-NA/H state with targeted W appearance locally closed";
        return result;
    case Sw92PhaseAssignedNoWStatus::no_w_two_h_locally_closed:
        result.status = Sw92PhaseAssignedBoundaryStatus::resolved_to_no_w_two_h;
        sw92_phase_assigned_boundary_publish_no_w(result, *result.neighbor_no_w);
        if (result.phases.size() != 2U) {
            result.status = Sw92PhaseAssignedBoundaryStatus::source_chain_inconsistent;
            result.phases.clear();
            result.diagnostic =
                "no-W two-H neighbor closed but phase publication data are incomplete";
            return result;
        }
        result.diagnostic =
            "W disappearance was re-solved as no-W H0+H1, including fixed-NA split/final review and targeted AQ-W appearance review";
        return result;
    case Sw92PhaseAssignedNoWStatus::aqueous_phase_witness_found:
    case Sw92PhaseAssignedNoWStatus::phase_disappearance_unresolved:
    case Sw92PhaseAssignedNoWStatus::single_phase_role_unresolved:
        result.status = Sw92PhaseAssignedBoundaryStatus::neighbor_topology_not_closed;
        result.diagnostic =
            "re-solved no-W neighbor is not locally closed; the disappearing W cannot be silently dropped";
        return result;
    case Sw92PhaseAssignedNoWStatus::higher_h_multiplicity_or_wrong_candidate:
        result.status =
            Sw92PhaseAssignedBoundaryStatus::higher_phase_count_or_wrong_candidate;
        result.diagnostic = result.neighbor_no_w->diagnostic;
        return result;
    case Sw92PhaseAssignedNoWStatus::indeterminate:
        result.status = Sw92PhaseAssignedBoundaryStatus::numerical_indeterminate;
        result.diagnostic = result.neighbor_no_w->diagnostic;
        return result;
    }
    result.status = Sw92PhaseAssignedBoundaryStatus::numerical_indeterminate;
    result.diagnostic = "unhandled no-W boundary state";
    return result;
}

inline bool sw92_phase_assigned_boundary_disappearance_status(
    Sw92PhaseAssignedThreePhaseStatus status) noexcept {
    return status == Sw92PhaseAssignedThreePhaseStatus::hydrocarbon_phase_disappearance ||
           status == Sw92PhaseAssignedThreePhaseStatus::aqueous_phase_disappearance;
}

} // namespace detail

/// Re-solve a C2b.2 disappearance route in its neighboring topology.
///
/// This function never accepts a lower phase count merely because a C2b.1
/// fraction is small. The neighbor equations are solved again from the retained
/// boundary compositions and the neighbor-specific finite stability/topology
/// review must close before a phase is removed.
[[nodiscard]] inline Sw92PhaseAssignedBoundaryResult
resolve_sw92_phase_assigned_c2b2_boundary(
    const Sw92PhaseAssignedJointResult& c1,
    const Sw92PhaseAssignedHSideWitnessResult& c2a1,
    const Sw92PhaseAssignedThreePhaseResult& c2b1,
    const Sw92PhaseAssignedC2b2Result& c2b2,
    const thermodynamics::Sw92Phase<double>& model,
    const Sw92PhaseAssignedPtOptions& options = {}) {
    Sw92PhaseAssignedBoundaryResult result;
    if (!c2b1.point || !c2b1.equations_converged() ||
        !detail::sw92_phase_assigned_boundary_disappearance_status(c2b1.status)) {
        result.status = Sw92PhaseAssignedBoundaryStatus::source_chain_inconsistent;
        result.diagnostic =
            "boundary re-solve requires an equation-converged C2b.1 disappearance state";
        return result;
    }
    if (!detail::sw92_phase_assigned_c2b1_source_matches(c1, c2a1) ||
        !detail::sw92_phase_assigned_c2b2_three_phase_matches(c1, c2b1) ||
        !detail::sw92_phase_assigned_c2b2_source_seed_matches(c1, c2a1, c2b1) ||
        c1.pressure_pa != c2b2.pressure_pa ||
        c1.temperature_k != c2b2.temperature_k ||
        c1.nacl_molality_mol_per_kg_water != c2b2.nacl_molality_mol_per_kg_water ||
        c1.feed != c2b2.feed || c1.dataset_id != c2b2.dataset_id ||
        c1.revision != c2b2.revision || c1.component_ids != c2b2.component_ids ||
        c2b2.selected_c2a1_witness_index != c2b1.source_witness_index) {
        result.status = Sw92PhaseAssignedBoundaryStatus::source_chain_inconsistent;
        result.diagnostic =
            "boundary re-solve source metadata or C2a1->C2b.1 witness provenance is inconsistent";
        return result;
    }

    const auto& state = *c2b1.point;
    const double minimum = c2b1.options.minimum_phase_fraction;
    const bool w_gone = state.aqueous_phase.mole_phase_fraction <= minimum;
    const bool h0_gone = state.hydrocarbon0_phase.mole_phase_fraction <= minimum;
    const bool h1_gone = state.hydrocarbon1_phase.mole_phase_fraction <= minimum;

    if (c2b2.status == Sw92PhaseAssignedC2b2Status::route_to_w_h) {
        if (w_gone || h0_gone == h1_gone) {
            result.status = Sw92PhaseAssignedBoundaryStatus::source_chain_inconsistent;
            result.diagnostic =
                "C2b.2 W+H route does not identify exactly one surviving H phase";
            return result;
        }
        const auto& survivor = h0_gone
            ? state.hydrocarbon1_phase.composition
            : state.hydrocarbon0_phase.composition;
        const auto& disappeared = h0_gone
            ? state.hydrocarbon0_phase.composition
            : state.hydrocarbon1_phase.composition;
        std::vector<std::vector<double>> starts;
        starts.push_back(disappeared);
        auto inherited = detail::sw92_phase_assigned_boundary_unique_starts(
            std::span<const std::vector<double>>{}, survivor, c1.feed,
            options.c2a1.log_composition_separation);
        starts.insert(starts.end(), inherited.begin(), inherited.end());
        return detail::sw92_phase_assigned_resolve_w_h_neighbor(
            c1.pressure_pa, c1.temperature_k, c1.feed,
            state.aqueous_phase.composition, survivor, model,
            c1.nacl_molality_mol_per_kg_water, options, starts);
    }

    if (c2b2.status == Sw92PhaseAssignedC2b2Status::route_to_no_w_h0_h1) {
        if (!w_gone || h0_gone || h1_gone) {
            result.status = Sw92PhaseAssignedBoundaryStatus::source_chain_inconsistent;
            result.diagnostic =
                "C2b.2 no-W route does not retain exactly two positive H phases";
            return result;
        }
        const std::vector<std::vector<double>> starts{
            state.hydrocarbon0_phase.composition,
            state.hydrocarbon1_phase.composition};
        return detail::sw92_phase_assigned_resolve_no_w_neighbor(
            c1.pressure_pa, c1.temperature_k, c1.feed, starts, model,
            c1.nacl_molality_mol_per_kg_water, options);
    }

    if (c2b2.status ==
        Sw92PhaseAssignedC2b2Status::single_phase_endpoint_unresolved) {
        if (!w_gone) {
            result.status = Sw92PhaseAssignedBoundaryStatus::source_chain_inconsistent;
            result.diagnostic =
                "single-phase endpoint route is inconsistent with a retained W fraction";
            return result;
        }
        std::vector<std::vector<double>> starts;
        if (!h0_gone) { starts.push_back(state.hydrocarbon0_phase.composition); }
        if (!h1_gone) { starts.push_back(state.hydrocarbon1_phase.composition); }
        if (starts.empty()) {
            result.status = Sw92PhaseAssignedBoundaryStatus::source_chain_inconsistent;
            result.diagnostic = "three-phase endpoint has no surviving phase composition";
            return result;
        }
        return detail::sw92_phase_assigned_resolve_no_w_neighbor(
            c1.pressure_pa, c1.temperature_k, c1.feed, starts, model,
            c1.nacl_molality_mol_per_kg_water, options);
    }

    result.status = Sw92PhaseAssignedBoundaryStatus::no_boundary_route;
    result.diagnostic = "C2b.2 result does not request a disappearance-neighbor route";
    return result;
}

/// Boundary-aware wrapper around the existing Profile-C PT topology driver.
///
/// Interior/local-closure behavior is unchanged. If the old driver stops before
/// C2b.2 because C2b.1 reached a phase-disappearance state, this adapter
/// recovers the equation-converged disappearance point, runs C2b.2, and then
/// re-solves the requested lower topology instead of clipping the small phase.
[[nodiscard]] inline Sw92PhaseAssignedBoundaryAwareResult
solve_sw92_phase_assigned_pt_boundary_aware(
    double pressure_pa, double temperature_k, std::span<const double> feed,
    const thermodynamics::Sw92Phase<double>& model,
    double nacl_molality_mol_per_kg_water,
    Sw92PhaseAssignedPtOptions options = {}) {
    Sw92PhaseAssignedBoundaryAwareResult result;
    result.base = solve_sw92_phase_assigned_pt(
        pressure_pa, temperature_k, feed, model,
        nacl_molality_mol_per_kg_water, options);
    result.status = result.base.status;
    result.phases = result.base.phases;
    result.diagnostic = result.base.diagnostic;

    if (result.base.locally_closed_phase_candidate() ||
        result.base.status == Sw92PhaseAssignedPtStatus::higher_phase_count_or_wrong_candidate ||
        result.base.status == Sw92PhaseAssignedPtStatus::numerical_indeterminate) {
        return result;
    }
    if (!result.base.c1 || !result.base.c2a1 ||
        result.base.c2a1->status !=
            Sw92PhaseAssignedHSideWitnessStatus::additional_nonaqueous_phase_witness_found) {
        return result;
    }

    const double c1_gibbs = result.base.c1->candidate()
        ? result.base.c1->candidate()->reduced_gibbs
        : std::numeric_limits<double>::infinity();
    const double c1_guard = result.base.c1->candidate()
        ? result.base.c1->candidate()->gibbs_roundoff_guard
        : 0.0;
    if (!std::isfinite(c1_gibbs)) { return result; }

    double best_gibbs = std::numeric_limits<double>::infinity();
    for (std::size_t witness_index = 0;
         witness_index < result.base.c2a1->negative_witnesses.size(); ++witness_index) {
        const auto& witness = result.base.c2a1->negative_witnesses[witness_index];
        if (!witness.usable_h_split_seed()) { continue; }
        if (result.disappearance_attempts >= options.max_c2b1_attempts) {
            result.disappearance_attempt_limit_reached = true;
            break;
        }
        ++result.disappearance_attempts;
        auto candidate = solve_sw92_phase_assigned_c2b1_candidate(
            *result.base.c1, *result.base.c2a1, witness_index, model,
            options.c2b1);
        if (!candidate.point || !candidate.equations_converged() ||
            !detail::sw92_phase_assigned_boundary_disappearance_status(candidate.status)) {
            continue;
        }
        const auto& point = *candidate.point;
        if (point.reduced_gibbs >
            c1_gibbs + c1_guard + point.gibbs_roundoff_guard) {
            continue;
        }
        auto review = review_sw92_phase_assigned_c2b2(
            *result.base.c1, *result.base.c2a1, candidate, model,
            options.c2b2);
        if (review.status != Sw92PhaseAssignedC2b2Status::route_to_w_h &&
            review.status != Sw92PhaseAssignedC2b2Status::route_to_no_w_h0_h1 &&
            review.status !=
                Sw92PhaseAssignedC2b2Status::single_phase_endpoint_unresolved) {
            continue;
        }
        if (point.reduced_gibbs < best_gibbs) {
            best_gibbs = point.reduced_gibbs;
            result.boundary_c2b1 = std::move(candidate);
            result.boundary_c2b2 = std::move(review);
        }
    }

    if (!result.boundary_c2b1 || !result.boundary_c2b2) {
        if (result.disappearance_attempt_limit_reached) {
            result.diagnostic +=
                "; disappearance-aware C2b.1 attempt quota exhausted before a routable boundary state was retained";
        }
        return result;
    }

    result.boundary = resolve_sw92_phase_assigned_c2b2_boundary(
        *result.base.c1, *result.base.c2a1, *result.boundary_c2b1,
        *result.boundary_c2b2, model, options);
    if (!result.boundary->neighbor_locally_closed()) {
        result.status = result.boundary->status ==
                Sw92PhaseAssignedBoundaryStatus::higher_phase_count_or_wrong_candidate
            ? Sw92PhaseAssignedPtStatus::higher_phase_count_or_wrong_candidate
            : (result.boundary->status ==
                       Sw92PhaseAssignedBoundaryStatus::numerical_indeterminate
                   ? Sw92PhaseAssignedPtStatus::numerical_indeterminate
                   : Sw92PhaseAssignedPtStatus::topology_unresolved);
        result.phases.clear();
        result.diagnostic = result.boundary->diagnostic;
        return result;
    }

    result.phases = result.boundary->phases;
    switch (result.boundary->status) {
    case Sw92PhaseAssignedBoundaryStatus::resolved_to_w_h:
        result.status = Sw92PhaseAssignedPtStatus::w_h_locally_closed;
        break;
    case Sw92PhaseAssignedBoundaryStatus::resolved_to_no_w_single_h:
        result.status = Sw92PhaseAssignedPtStatus::no_w_single_h_locally_closed;
        break;
    case Sw92PhaseAssignedBoundaryStatus::resolved_to_no_w_two_h:
        result.status = Sw92PhaseAssignedPtStatus::no_w_two_h_locally_closed;
        break;
    case Sw92PhaseAssignedBoundaryStatus::no_boundary_route:
    case Sw92PhaseAssignedBoundaryStatus::neighbor_topology_not_closed:
    case Sw92PhaseAssignedBoundaryStatus::higher_phase_count_or_wrong_candidate:
    case Sw92PhaseAssignedBoundaryStatus::source_chain_inconsistent:
    case Sw92PhaseAssignedBoundaryStatus::numerical_indeterminate:
        result.status = Sw92PhaseAssignedPtStatus::topology_unresolved;
        result.phases.clear();
        break;
    }
    result.diagnostic = result.boundary->diagnostic;
    return result;
}

} // namespace mpmc::flash

#endif // MPMC_FLASH_SW92_PHASE_ASSIGNED_BOUNDARY_HPP