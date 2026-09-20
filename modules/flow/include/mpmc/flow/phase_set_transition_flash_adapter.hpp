#ifndef MPMC_FLOW_PHASE_SET_TRANSITION_FLASH_ADAPTER_HPP
#define MPMC_FLOW_PHASE_SET_TRANSITION_FLASH_ADAPTER_HPP

#include <mpmc/flash/pt_flash_backend.hpp>
#include <mpmc/flow/phase_set_transition.hpp>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <optional>
#include <span>
#include <stdexcept>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

namespace mpmc::flow {

inline constexpr std::string_view
    phase_set_transition_flash_adapter_convention =
        "flow/flash/phase-set-transition-adapter/v1";

namespace phase_set_transition_flash_detail {

[[nodiscard]] inline
PhaseSetTransitionTrigger
map_trigger(
    mpmc::flash::PtPhaseTransitionTrigger
        trigger) {
    using Flash =
        mpmc::flash::PtPhaseTransitionTrigger;
    switch (trigger) {
    case Flash::initial_stability_witness:
        return PhaseSetTransitionTrigger::
            stability_witness;
    case Flash::final_phase_set_instability:
        return PhaseSetTransitionTrigger::
            final_phase_set_instability;
    case Flash::phase_disappearance:
        return PhaseSetTransitionTrigger::
            phase_disappearance;
    case Flash::provider_topology_witness:
        return PhaseSetTransitionTrigger::
            provider_topology_witness;
    case Flash::provider_boundary_route:
        return PhaseSetTransitionTrigger::
            provider_boundary_route;
    }
    throw std::invalid_argument(
        "mpmc::flow: unknown flash transition trigger");
}

[[nodiscard]] inline bool
is_source_challenge(
    const mpmc::flash::
        PtPhaseTransitionEvidence&
            evidence,
    std::size_t source_phase_count) {
    using Resolution =
        mpmc::flash::
            PtPhaseTransitionResolution;
    return evidence.source_phase_count ==
               source_phase_count &&
        evidence.resolution !=
            Resolution::indeterminate &&
        evidence.resolution !=
            Resolution::target_resolve_failed;
}

[[nodiscard]] inline const
mpmc::flash::PtPhaseTransitionEvidence*
accepted_target_evidence(
    const mpmc::flash::
        PtPhaseTransitionReport& report,
    std::size_t target_phase_count) {
    using Resolution =
        mpmc::flash::
            PtPhaseTransitionResolution;
    for (const auto& evidence :
         report.evidence) {
        if (evidence.target_phase_count &&
            *evidence.target_phase_count ==
                target_phase_count &&
            evidence.resolution ==
                Resolution::accepted_target &&
            evidence.fresh_target_solve_attempted &&
            evidence.target_topology_closed) {
            return &evidence;
        }
    }
    return nullptr;
}

[[nodiscard]] inline const
mpmc::flash::PtPhaseTransitionEvidence*
direct_accepted_disappearance(
    const mpmc::flash::
        PtPhaseTransitionReport& report,
    std::size_t source_phase_count,
    std::size_t target_phase_count) {
    using Trigger =
        mpmc::flash::PtPhaseTransitionTrigger;
    using Resolution =
        mpmc::flash::
            PtPhaseTransitionResolution;
    for (const auto& evidence :
         report.evidence) {
        if (evidence.source_phase_count ==
                source_phase_count &&
            evidence.target_phase_count &&
            *evidence.target_phase_count ==
                target_phase_count &&
            evidence.resolution ==
                Resolution::accepted_target &&
            evidence.fresh_target_solve_attempted &&
            evidence.target_topology_closed &&
            (evidence.trigger ==
                 Trigger::phase_disappearance ||
             evidence.trigger ==
                 Trigger::provider_boundary_route)) {
            return &evidence;
        }
    }
    return nullptr;
}

} // namespace phase_set_transition_flash_detail

[[nodiscard]] inline std::optional<
    PhaseSetTransitionCandidate>
make_phase_set_transition_candidate_from_flash(
    std::size_t source_phase_count,
    const mpmc::flash::
        PtFlashBackendResult& result,
    std::span<const double>
        target_phase_molar_density_mol_per_m3) {
    using namespace
        phase_set_transition_flash_detail;

    if (source_phase_count == 0U ||
        source_phase_count >
            fixed_three_phase_count ||
        !result.structurally_valid()) {
        throw std::invalid_argument(
            "mpmc::flow: invalid source phase count or flash transition result");
    }

    const auto* accepted =
        result.accepted_phase_set();
    if (accepted == nullptr) {
        return std::nullopt;
    }

    const std::size_t target_phase_count =
        accepted->phases.size();
    if (target_phase_count ==
        source_phase_count) {
        return std::nullopt;
    }
    if (target_phase_count == 0U ||
        target_phase_count >
            fixed_three_phase_count ||
        target_phase_molar_density_mol_per_m3
                .size() !=
            target_phase_count) {
        throw std::invalid_argument(
            "mpmc::flow: flash transition target density cardinality mismatch");
    }

    const mpmc::flash::
        PtPhaseTransitionEvidence*
        governing_evidence = nullptr;

    if (target_phase_count <
        source_phase_count) {
        governing_evidence =
            direct_accepted_disappearance(
                result.transition_report,
                source_phase_count,
                target_phase_count);
        if (governing_evidence == nullptr) {
            throw std::invalid_argument(
                "mpmc::flow: accepted lower-phase flash target lacks fresh disappearance/boundary evidence");
        }
    } else {
        bool source_challenged = false;
        for (const auto& evidence :
             result.transition_report.evidence) {
            source_challenged =
                source_challenged ||
                is_source_challenge(
                    evidence,
                    source_phase_count);
        }
        governing_evidence =
            accepted_target_evidence(
                result.transition_report,
                target_phase_count);
        if (!source_challenged ||
            governing_evidence == nullptr) {
            throw std::invalid_argument(
                "mpmc::flow: accepted higher-phase flash target lacks source-instability and fresh-target evidence");
        }
    }

    PhaseSetTransitionCandidate candidate;
    candidate.source_phase_count =
        source_phase_count;
    candidate.target_phase_count =
        target_phase_count;
    candidate.trigger =
        map_trigger(
            governing_evidence->trigger);
    candidate.status =
        PhaseSetTransitionCandidateStatus::
            target_resolved;
    candidate.pressure_pa =
        result.solution.pressure_pa;
    candidate.temperature_k =
        result.solution.temperature_k;
    candidate.component_ids =
        result.capability.component_ids;
    candidate.evidence_profile =
        governing_evidence
            ->provider_evidence_profile;
    candidate.diagnostic =
        governing_evidence->diagnostic;

    candidate.target_phases.reserve(
        target_phase_count);
    for (std::size_t phase = 0U;
         phase < target_phase_count;
         ++phase) {
        const auto& source =
            accepted->phases[phase];
        candidate.target_phases.push_back(
            {
                source.mole_phase_fraction,
                source.composition,
                target_phase_molar_density_mol_per_m3[
                    phase]});
    }

    return candidate;
}

[[nodiscard]] inline
PhaseSetTransitionCandidate
make_unresolved_phase_set_transition_candidate(
    std::size_t source_phase_count,
    std::size_t target_phase_count,
    PhaseSetTransitionTrigger trigger,
    double pressure_pa,
    double temperature_k,
    std::vector<std::string>
        component_ids,
    std::string evidence_profile,
    std::string diagnostic) {
    PhaseSetTransitionCandidate candidate;
    candidate.source_phase_count =
        source_phase_count;
    candidate.target_phase_count =
        target_phase_count;
    candidate.trigger = trigger;
    candidate.status =
        PhaseSetTransitionCandidateStatus::
            target_resolve_required;
    candidate.pressure_pa =
        pressure_pa;
    candidate.temperature_k =
        temperature_k;
    candidate.component_ids =
        std::move(component_ids);
    candidate.evidence_profile =
        std::move(evidence_profile);
    candidate.diagnostic =
        std::move(diagnostic);

    phase_set_transition_detail::
        validate_candidate_direction(
            candidate);
    phase_set_transition_detail::
        validate_component_ids(
            candidate.component_ids);

    if (!std::isfinite(
            candidate.pressure_pa) ||
        !(candidate.pressure_pa > 0.0) ||
        !std::isfinite(
            candidate.temperature_k) ||
        !(candidate.temperature_k > 0.0) ||
        candidate.evidence_profile.empty()) {
        throw std::invalid_argument(
            "mpmc::flow: unresolved phase-transition candidate lacks valid provenance");
    }

    return candidate;
}

} // namespace mpmc::flow

#endif // MPMC_FLOW_PHASE_SET_TRANSITION_FLASH_ADAPTER_HPP
