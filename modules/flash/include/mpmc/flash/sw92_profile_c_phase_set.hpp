#ifndef MPMC_FLASH_SW92_PROFILE_C_PHASE_SET_HPP
#define MPMC_FLASH_SW92_PROFILE_C_PHASE_SET_HPP

#include <mpmc/flash/pt_phase_set.hpp>
#include <mpmc/flash/sw92_phase_assigned_boundary.hpp>

#include <cmath>
#include <cstddef>
#include <optional>
#include <span>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

namespace mpmc::flash {

inline constexpr std::string_view sw92_profile_c_phase_set_publication_convention =
    "SW92-equilibrium/phase-assigned-aq-na-joint/authoritative-phase-set-publication/v1";

// Model-specific identity stays outside the generic PtCandidatePhase payload.
// In particular, NA family identity is not a liquid/vapor morphology label.
struct Sw92ProfileCPhaseMetadata {
    Sw92PhaseAssignedPtPhysicalRole physical_role{
        Sw92PhaseAssignedPtPhysicalRole::nonaqueous_unclassified};
    thermodynamics::SwPhaseFamily thermodynamic_family{
        thermodynamics::SwPhaseFamily::nonaqueous};
};

struct Sw92ProfileCPtPhaseSetResult {
    static constexpr bool global_stability_proven = false;
    static constexpr bool morphology_resolved = false;

    PtPhaseSetResult solution;
    double input_feed_sum{};
    double nacl_molality_mol_per_kg_water{};
    std::string dataset_id;
    std::string revision;
    std::vector<std::string> component_ids;
    std::string model_profile{thermodynamics::sw92_corrected_profile};
    std::string phase_convention{thermodynamics::sw92_pt_convention};
    std::string equilibrium_profile{sw92_phase_assigned_aq_na_joint_profile};
    std::string orchestration_convention{sw92_phase_assigned_pt_convention};
    std::string boundary_convention{sw92_phase_assigned_boundary_convention};
    std::string publication_convention{
        sw92_profile_c_phase_set_publication_convention};
    std::vector<Sw92ProfileCPhaseMetadata> phase_metadata;

    [[nodiscard]] bool accepted_phase_set_published() const noexcept {
        const auto* accepted = solution.accepted_phase_set();
        return accepted != nullptr &&
               phase_metadata.size() == accepted->phases.size();
    }
};

namespace detail {

template <class Result>
[[nodiscard]] inline bool sw92_profile_c_snapshot_matches(
    const Sw92PhaseAssignedPtResult& base, const Result& other) {
    return base.pressure_pa == other.pressure_pa &&
           base.temperature_k == other.temperature_k &&
           base.feed == other.feed &&
           base.nacl_molality_mol_per_kg_water ==
               other.nacl_molality_mol_per_kg_water &&
           base.dataset_id == other.dataset_id &&
           base.revision == other.revision &&
           base.component_ids == other.component_ids &&
           std::string_view{other.model_profile} ==
               thermodynamics::sw92_corrected_profile &&
           std::string_view{other.phase_convention} ==
               thermodynamics::sw92_pt_convention &&
           std::string_view{other.equilibrium_profile} ==
               sw92_phase_assigned_aq_na_joint_profile;
}

[[nodiscard]] inline bool sw92_profile_c_no_w_contract_matches(
    const Sw92PhaseAssignedPtResult& base,
    const Sw92PhaseAssignedNoWResult& no_w) {
    if (!sw92_profile_c_snapshot_matches(base, no_w) ||
        std::string_view{no_w.adapter_convention} !=
            sw92_phase_assigned_no_w_convention) {
        return false;
    }
    const auto& family = no_w.hydrocarbon_flash;
    const auto& solution = family.solution;
    return family.dataset_id == base.dataset_id &&
           family.revision == base.revision &&
           family.component_ids == base.component_ids &&
           family.nacl_molality_mol_per_kg_water ==
               base.nacl_molality_mol_per_kg_water &&
           family.family == thermodynamics::SwPhaseFamily::nonaqueous &&
           std::string_view{family.model_profile} ==
               thermodynamics::sw92_corrected_profile &&
           std::string_view{family.phase_convention} ==
               thermodynamics::sw92_pt_convention &&
           std::string_view{family.equilibrium_algorithm} ==
               sw92_family_vle_algorithm &&
           solution.initial_stability.pressure_pa == base.pressure_pa &&
           solution.initial_stability.temperature_k == base.temperature_k &&
           solution.initial_stability.feed == base.feed;
}

[[nodiscard]] inline bool sw92_profile_c_c1_contract_matches(
    const Sw92PhaseAssignedPtResult& base,
    const Sw92PhaseAssignedJointResult& c1) {
    return sw92_profile_c_snapshot_matches(base, c1) &&
           std::string_view{c1.primitive_convention} ==
               sw92_phase_assigned_aq_na_joint_primitive;
}

[[nodiscard]] inline bool sw92_profile_c_c2a1_contract_matches(
    const Sw92PhaseAssignedPtResult& base,
    const Sw92PhaseAssignedHSideWitnessResult& c2a1) {
    return sw92_profile_c_snapshot_matches(base, c2a1) &&
           std::string_view{c2a1.witness_convention} ==
               sw92_phase_assigned_h_side_na_witness_convention;
}

[[nodiscard]] inline bool sw92_profile_c_base_contract_valid(
    const Sw92PhaseAssignedPtResult& base) {
    if (base.feed.empty() || base.component_ids.size() != base.feed.size() ||
        std::string_view{base.model_profile} != thermodynamics::sw92_corrected_profile ||
        std::string_view{base.phase_convention} != thermodynamics::sw92_pt_convention ||
        std::string_view{base.equilibrium_profile} !=
            sw92_phase_assigned_aq_na_joint_profile ||
        std::string_view{base.orchestration_convention} !=
            sw92_phase_assigned_pt_convention ||
        !sw92_profile_c_no_w_contract_matches(base, base.no_w)) {
        return false;
    }
    try {
        (void)stability_check_composition(base.feed);
    } catch (const std::exception&) {
        return false;
    }
    return true;
}

[[nodiscard]] inline bool sw92_profile_c_activity_valid(
    const StabilityPhase& activity, std::size_t component_count) {
    if (!activity.smooth || activity.ln_phi.size() != component_count) {
        return false;
    }
    for (double value : activity.ln_phi) {
        if (!std::isfinite(value)) { return false; }
    }
    return true;
}

[[nodiscard]] inline bool sw92_profile_c_summary_matches(
    const Sw92PhaseAssignedPtPhase& summary,
    double fraction,
    const std::vector<double>& composition,
    std::optional<double> compressibility_factor,
    Sw92PhaseAssignedPtPhysicalRole role,
    thermodynamics::SwPhaseFamily family) {
    return summary.mole_phase_fraction == fraction &&
           summary.composition == composition &&
           summary.compressibility_factor == compressibility_factor &&
           summary.physical_role == role &&
           summary.thermodynamic_family == family;
}

inline void sw92_profile_c_reject_publication(
    Sw92ProfileCPtPhaseSetResult& result,
    const Sw92PhaseAssignedBoundaryAwareResult& source,
    std::string_view reason) {
    result.solution.status = PtPhaseSetStatus::indeterminate;
    result.solution.candidate_phase_set.reset();
    result.phase_metadata.clear();
    result.solution.diagnostic =
        "SW92 Profile-C phase-set publication rejected: " + std::string(reason);
    if (!source.diagnostic.empty()) {
        result.solution.diagnostic += "; source: " + source.diagnostic;
    }
}

inline void sw92_profile_c_append_phase(
    Sw92ProfileCPtPhaseSetResult& result,
    PtCandidatePhase phase,
    Sw92PhaseAssignedPtPhysicalRole role,
    thermodynamics::SwPhaseFamily family) {
    if (!result.solution.candidate_phase_set) {
        result.solution.candidate_phase_set.emplace();
    }
    result.solution.candidate_phase_set->phases.push_back(std::move(phase));
    result.phase_metadata.push_back({role, family});
}

template <class Phase>
[[nodiscard]] inline PtCandidatePhase sw92_profile_c_project_owned_phase(
    const Phase& source) {
    PtCandidatePhase phase;
    phase.mole_phase_fraction = source.mole_phase_fraction;
    phase.composition = source.composition;
    phase.activity = source.activity;
    phase.compressibility_factor = source.compressibility_factor;
    return phase;
}

[[nodiscard]] inline bool sw92_profile_c_project_no_w(
    Sw92ProfileCPtPhaseSetResult& result,
    const Sw92PhaseAssignedBoundaryAwareResult& source,
    const Sw92PhaseAssignedNoWResult& no_w) {
    const auto& base = source.base;
    if (!sw92_profile_c_no_w_contract_matches(base, no_w) ||
        !no_w.no_w_locally_closed()) {
        return false;
    }
    const auto& flash = no_w.hydrocarbon_flash.solution;
    const std::size_t count = no_w.retained_na_candidate_count();
    if (count == 0U || count > 2U || source.phases.size() != count ||
        no_w.retained_na_candidate_fractions.size() != count) {
        return false;
    }

    if (count == 1U) {
        if (source.status !=
                Sw92PhaseAssignedPtStatus::no_w_single_h_locally_closed ||
            no_w.status !=
                Sw92PhaseAssignedNoWStatus::no_w_single_h_locally_closed ||
            flash.status != PtSplitStatus::single_phase_no_instability_found ||
            !flash.initial_stability.reference ||
            no_w.retained_na_candidate_compositions.front() != base.feed ||
            no_w.retained_na_candidate_fractions.front() != 1.0 ||
            !sw92_profile_c_activity_valid(
                *flash.initial_stability.reference, base.feed.size()) ||
            !sw92_profile_c_summary_matches(
                source.phases.front(), 1.0, base.feed, std::nullopt,
                Sw92PhaseAssignedPtPhysicalRole::nonaqueous_unclassified,
                thermodynamics::SwPhaseFamily::nonaqueous)) {
            return false;
        }
        PtCandidatePhase phase;
        phase.mole_phase_fraction = 1.0;
        phase.composition = base.feed;
        phase.activity = *flash.initial_stability.reference;
        // The accepted fixed-family single phase stores activity but no selected
        // Z. A publication projection must not call the EOS merely to fill it.
        phase.compressibility_factor.reset();
        sw92_profile_c_append_phase(
            result, std::move(phase),
            Sw92PhaseAssignedPtPhysicalRole::nonaqueous_unclassified,
            thermodynamics::SwPhaseFamily::nonaqueous);
        return true;
    }

    if (source.status != Sw92PhaseAssignedPtStatus::no_w_two_h_locally_closed ||
        no_w.status != Sw92PhaseAssignedNoWStatus::no_w_two_h_locally_closed ||
        flash.status != PtSplitStatus::two_phase_no_instability_found ||
        !flash.final_stability ||
        flash.final_stability->status != StabilityStatus::no_instability_found) {
        return false;
    }
    const auto* pair = flash.candidate();
    if (pair == nullptr ||
        !sw92_profile_c_activity_valid(pair->liquid.activity, base.feed.size()) ||
        !sw92_profile_c_activity_valid(pair->vapor.activity, base.feed.size())) {
        return false;
    }
    const double fractions[2]{
        1.0 - pair->fractions.vapor_fraction,
        pair->fractions.vapor_fraction};
    const std::vector<double>* compositions[2]{
        &pair->fractions.liquid, &pair->fractions.vapor};
    const PtSplitPhase* properties[2]{&pair->liquid, &pair->vapor};
    for (std::size_t i = 0; i < 2U; ++i) {
        if (no_w.retained_na_candidate_compositions[i] != *compositions[i] ||
            no_w.retained_na_candidate_fractions[i] != fractions[i] ||
            !sw92_profile_c_summary_matches(
                source.phases[i], fractions[i], *compositions[i], properties[i]->z,
                Sw92PhaseAssignedPtPhysicalRole::nonaqueous_unclassified,
                thermodynamics::SwPhaseFamily::nonaqueous)) {
            return false;
        }
        PtCandidatePhase phase;
        phase.mole_phase_fraction = fractions[i];
        phase.composition = *compositions[i];
        phase.activity = properties[i]->activity;
        phase.compressibility_factor = properties[i]->z;
        sw92_profile_c_append_phase(
            result, std::move(phase),
            Sw92PhaseAssignedPtPhysicalRole::nonaqueous_unclassified,
            thermodynamics::SwPhaseFamily::nonaqueous);
    }
    return true;
}

[[nodiscard]] inline bool sw92_profile_c_project_w_h(
    Sw92ProfileCPtPhaseSetResult& result,
    const Sw92PhaseAssignedBoundaryAwareResult& source,
    const Sw92PhaseAssignedJointResult& c1,
    const Sw92PhaseAssignedHSideWitnessResult& c2a1) {
    const auto& base = source.base;
    const auto* point = c1.candidate();
    if (source.status != Sw92PhaseAssignedPtStatus::w_h_locally_closed ||
        source.phases.size() != 2U || point == nullptr ||
        !sw92_profile_c_c1_contract_matches(base, c1) ||
        !sw92_profile_c_c2a1_contract_matches(base, c2a1) ||
        c2a1.status != Sw92PhaseAssignedHSideWitnessStatus::
                            no_additional_nonaqueous_witness_found ||
        !sw92_profile_c_activity_valid(
            point->aqueous_phase.activity, base.feed.size()) ||
        !sw92_profile_c_activity_valid(
            point->nonaqueous_phase.activity, base.feed.size()) ||
        !sw92_profile_c_summary_matches(
            source.phases[0], point->aqueous_phase.mole_phase_fraction,
            point->aqueous_phase.composition,
            point->aqueous_phase.compressibility_factor,
            Sw92PhaseAssignedPtPhysicalRole::aqueous,
            thermodynamics::SwPhaseFamily::aqueous) ||
        !sw92_profile_c_summary_matches(
            source.phases[1], point->nonaqueous_phase.mole_phase_fraction,
            point->nonaqueous_phase.composition,
            point->nonaqueous_phase.compressibility_factor,
            Sw92PhaseAssignedPtPhysicalRole::nonaqueous_unclassified,
            thermodynamics::SwPhaseFamily::nonaqueous)) {
        return false;
    }
    sw92_profile_c_append_phase(
        result, sw92_profile_c_project_owned_phase(point->aqueous_phase),
        Sw92PhaseAssignedPtPhysicalRole::aqueous,
        thermodynamics::SwPhaseFamily::aqueous);
    sw92_profile_c_append_phase(
        result, sw92_profile_c_project_owned_phase(point->nonaqueous_phase),
        Sw92PhaseAssignedPtPhysicalRole::nonaqueous_unclassified,
        thermodynamics::SwPhaseFamily::nonaqueous);
    return true;
}

[[nodiscard]] inline bool sw92_profile_c_c2_chain_matches(
    const Sw92PhaseAssignedPtResult& base,
    const Sw92PhaseAssignedThreePhaseResult& c2b1,
    const Sw92PhaseAssignedC2b2Result& c2b2) {
    if (!base.c1 || !base.c2a1 ||
        !sw92_profile_c_c1_contract_matches(base, *base.c1) ||
        !sw92_profile_c_c2a1_contract_matches(base, *base.c2a1) ||
        !sw92_profile_c_snapshot_matches(base, c2b1) ||
        !sw92_profile_c_snapshot_matches(base, c2b2) ||
        std::string_view{c2b1.primitive_convention} !=
            sw92_phase_assigned_c2b1_three_phase_convention ||
        std::string_view{c2b2.review_convention} !=
            sw92_phase_assigned_c2b2_closure_convention ||
        !sw92_phase_assigned_c2b1_source_matches(*base.c1, *base.c2a1) ||
        !sw92_phase_assigned_c2b2_three_phase_matches(*base.c1, c2b1) ||
        !sw92_phase_assigned_c2b2_source_seed_matches(*base.c1, *base.c2a1, c2b1) ||
        base.c1->pressure_pa != c2b2.pressure_pa ||
        base.c1->temperature_k != c2b2.temperature_k ||
        base.c1->nacl_molality_mol_per_kg_water !=
            c2b2.nacl_molality_mol_per_kg_water ||
        base.c1->feed != c2b2.feed ||
        base.c1->dataset_id != c2b2.dataset_id ||
        base.c1->revision != c2b2.revision ||
        base.c1->component_ids != c2b2.component_ids ||
        c2b2.selected_c2a1_witness_index != c2b1.source_witness_index) {
        return false;
    }
    return true;
}

[[nodiscard]] inline bool sw92_profile_c_project_three_phase(
    Sw92ProfileCPtPhaseSetResult& result,
    const Sw92PhaseAssignedBoundaryAwareResult& source) {
    const auto& base = source.base;
    if (source.status != Sw92PhaseAssignedPtStatus::w_h0_h1_locally_closed ||
        source.phases.size() != 3U || !base.c2b1 || !base.c2b2 ||
        base.c2b2->status != Sw92PhaseAssignedC2b2Status::
                                w_present_h_multiplicity_locally_closed ||
        !sw92_profile_c_c2_chain_matches(base, *base.c2b1, *base.c2b2)) {
        return false;
    }
    const auto* point = base.c2b1->candidate();
    if (point == nullptr) { return false; }

    const auto append = [&](const Sw92PhaseAssignedThreePhasePhase& phase,
                            const Sw92PhaseAssignedPtPhase& summary,
                            Sw92PhaseAssignedPtPhysicalRole role,
                            thermodynamics::SwPhaseFamily family) {
        if (!sw92_profile_c_activity_valid(phase.activity, base.feed.size()) ||
            !sw92_profile_c_summary_matches(
                summary, phase.mole_phase_fraction, phase.composition,
                phase.compressibility_factor, role, family)) {
            return false;
        }
        sw92_profile_c_append_phase(
            result, sw92_profile_c_project_owned_phase(phase), role, family);
        return true;
    };
    return append(
               point->aqueous_phase, source.phases[0],
               Sw92PhaseAssignedPtPhysicalRole::aqueous,
               thermodynamics::SwPhaseFamily::aqueous) &&
           append(
               point->hydrocarbon0_phase, source.phases[1],
               Sw92PhaseAssignedPtPhysicalRole::nonaqueous_unclassified,
               thermodynamics::SwPhaseFamily::nonaqueous) &&
           append(
               point->hydrocarbon1_phase, source.phases[2],
               Sw92PhaseAssignedPtPhysicalRole::nonaqueous_unclassified,
               thermodynamics::SwPhaseFamily::nonaqueous);
}

[[nodiscard]] inline bool sw92_profile_c_boundary_chain_matches(
    const Sw92PhaseAssignedBoundaryAwareResult& source) {
    const auto& base = source.base;
    if (!source.boundary || !source.boundary_c2b1 || !source.boundary_c2b2 ||
        !base.c1 || !base.c2a1 ||
        std::string_view{source.boundary->convention} !=
            sw92_phase_assigned_boundary_convention ||
        !sw92_profile_c_c2_chain_matches(
            base, *source.boundary_c2b1, *source.boundary_c2b2)) {
        return false;
    }
    const auto review_status = source.boundary_c2b2->status;
    switch (source.boundary->status) {
    case Sw92PhaseAssignedBoundaryStatus::resolved_to_w_h:
        return review_status == Sw92PhaseAssignedC2b2Status::route_to_w_h &&
               source.boundary->neighbor_c1 && source.boundary->neighbor_c2a1 &&
               sw92_profile_c_c1_contract_matches(
                   base, *source.boundary->neighbor_c1) &&
               sw92_profile_c_c2a1_contract_matches(
                   base, *source.boundary->neighbor_c2a1);
    case Sw92PhaseAssignedBoundaryStatus::resolved_to_no_w_single_h:
    case Sw92PhaseAssignedBoundaryStatus::resolved_to_no_w_two_h:
        return (review_status ==
                    Sw92PhaseAssignedC2b2Status::route_to_no_w_h0_h1 ||
                review_status ==
                    Sw92PhaseAssignedC2b2Status::single_phase_endpoint_unresolved) &&
               source.boundary->neighbor_no_w &&
               sw92_profile_c_no_w_contract_matches(
                   base, *source.boundary->neighbor_no_w);
    case Sw92PhaseAssignedBoundaryStatus::no_boundary_route:
    case Sw92PhaseAssignedBoundaryStatus::neighbor_topology_not_closed:
    case Sw92PhaseAssignedBoundaryStatus::higher_phase_count_or_wrong_candidate:
    case Sw92PhaseAssignedBoundaryStatus::source_chain_inconsistent:
    case Sw92PhaseAssignedBoundaryStatus::numerical_indeterminate:
        return false;
    }
    return false;
}

} // namespace detail

// Pure publication projection. No EOS/provider call, root solve, stability
// search, normalization, equilibrium iteration, tolerance change, or morphology
// classification is performed here. Malformed or provenance-inconsistent source
// objects are conservatively downgraded to indeterminate.
[[nodiscard]] inline Sw92ProfileCPtPhaseSetResult
project_sw92_profile_c_pt_phase_set(
    const Sw92PhaseAssignedBoundaryAwareResult& source) {
    Sw92ProfileCPtPhaseSetResult result;
    const auto& base = source.base;
    result.solution.capability.maximum_phase_count = 3U;
    result.solution.pressure_pa = base.pressure_pa;
    result.solution.temperature_k = base.temperature_k;
    result.solution.feed = base.feed;
    result.solution.global_stability_proven = false;
    result.solution.diagnostic = source.diagnostic;
    result.input_feed_sum = base.input_feed_sum;
    result.nacl_molality_mol_per_kg_water =
        base.nacl_molality_mol_per_kg_water;
    result.dataset_id = base.dataset_id;
    result.revision = base.revision;
    result.component_ids = base.component_ids;
    result.model_profile = base.model_profile;
    result.phase_convention = base.phase_convention;
    result.equilibrium_profile = base.equilibrium_profile;
    result.orchestration_convention = base.orchestration_convention;

    if (!detail::sw92_profile_c_base_contract_valid(base)) {
        detail::sw92_profile_c_reject_publication(
            result, source, "invalid base/model/no-W provenance snapshot");
        return result;
    }

    bool projected = false;
    if (source.boundary) {
        if (!detail::sw92_profile_c_boundary_chain_matches(source)) {
            detail::sw92_profile_c_reject_publication(
                result, source, "boundary result lacks a complete consistent fresh re-solve chain");
            return result;
        }
        switch (source.boundary->status) {
        case Sw92PhaseAssignedBoundaryStatus::resolved_to_w_h:
            projected = detail::sw92_profile_c_project_w_h(
                result, source, *source.boundary->neighbor_c1,
                *source.boundary->neighbor_c2a1);
            break;
        case Sw92PhaseAssignedBoundaryStatus::resolved_to_no_w_single_h:
        case Sw92PhaseAssignedBoundaryStatus::resolved_to_no_w_two_h:
            projected = detail::sw92_profile_c_project_no_w(
                result, source, *source.boundary->neighbor_no_w);
            break;
        case Sw92PhaseAssignedBoundaryStatus::no_boundary_route:
        case Sw92PhaseAssignedBoundaryStatus::neighbor_topology_not_closed:
        case Sw92PhaseAssignedBoundaryStatus::higher_phase_count_or_wrong_candidate:
        case Sw92PhaseAssignedBoundaryStatus::source_chain_inconsistent:
        case Sw92PhaseAssignedBoundaryStatus::numerical_indeterminate:
            projected = false;
            break;
        }
    } else {
        switch (source.status) {
        case Sw92PhaseAssignedPtStatus::no_w_single_h_locally_closed:
        case Sw92PhaseAssignedPtStatus::no_w_two_h_locally_closed:
            projected = detail::sw92_profile_c_project_no_w(
                result, source, base.no_w);
            break;
        case Sw92PhaseAssignedPtStatus::w_h_locally_closed:
            projected = base.c1 && base.c2a1 &&
                detail::sw92_profile_c_project_w_h(
                    result, source, *base.c1, *base.c2a1);
            break;
        case Sw92PhaseAssignedPtStatus::w_h0_h1_locally_closed:
            projected = detail::sw92_profile_c_project_three_phase(result, source);
            break;
        case Sw92PhaseAssignedPtStatus::higher_phase_count_or_wrong_candidate:
        case Sw92PhaseAssignedPtStatus::topology_unresolved:
        case Sw92PhaseAssignedPtStatus::numerical_indeterminate:
            projected = false;
            break;
        }
    }

    if (!projected || !result.solution.candidate_phase_set ||
        result.solution.candidate_phase_set->phases.empty() ||
        result.solution.candidate_phase_set->phases.size() > 3U ||
        result.phase_metadata.size() !=
            result.solution.candidate_phase_set->phases.size()) {
        detail::sw92_profile_c_reject_publication(
            result, source,
            source.locally_closed_phase_candidate()
                ? "locally closed source failed publication integrity/provenance guards"
                : "source topology is not authoritatively closed under the declared finite review");
        return result;
    }

    result.solution.status = PtPhaseSetStatus::accepted;
    if (!result.accepted_phase_set_published()) {
        detail::sw92_profile_c_reject_publication(
            result, source, "generic accepted-phase structural guard rejected the projection");
    }
    return result;
}

// Production entry point. The boundary-aware Profile-C PT driver is the sole
// solver/orchestration source; publication is a pure owned-result projection.
[[nodiscard]] inline Sw92ProfileCPtPhaseSetResult
solve_sw92_profile_c_pt_phase_set(
    double pressure_pa, double temperature_k, std::span<const double> feed,
    const thermodynamics::Sw92Phase<double>& model,
    double nacl_molality_mol_per_kg_water,
    Sw92PhaseAssignedPtOptions options = {}) {
    const auto source = solve_sw92_phase_assigned_pt_boundary_aware(
        pressure_pa, temperature_k, feed, model,
        nacl_molality_mol_per_kg_water, std::move(options));
    return project_sw92_profile_c_pt_phase_set(source);
}

} // namespace mpmc::flash

#endif // MPMC_FLASH_SW92_PROFILE_C_PHASE_SET_HPP
