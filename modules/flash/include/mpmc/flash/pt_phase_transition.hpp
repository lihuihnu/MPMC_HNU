#ifndef MPMC_FLASH_PT_PHASE_TRANSITION_HPP
#define MPMC_FLASH_PT_PHASE_TRANSITION_HPP

#include <algorithm>
#include <cstddef>
#include <optional>
#include <string>
#include <string_view>
#include <vector>

namespace mpmc::flash {

inline constexpr std::string_view pt_phase_transition_convention =
    "PT/phase-transition-boundary/v1";

// A backend may either only detect that a neighboring phase count must be
// reconsidered, or it may own a fresh solve/review path for that target
// topology. Detection-only support must never be interpreted as permission to
// delete/add a phase directly from the source candidate.
enum class PtPhaseTransitionSupport {
    detection_only,
    fresh_target_resolve
};

// Evidence describes why the current phase count/topology is being challenged.
// These values are algorithm-neutral; provider-specific details remain in the
// evidence profile/diagnostic and are not promoted to universal phase labels.
enum class PtPhaseTransitionTrigger {
    initial_stability_witness,
    final_phase_set_instability,
    phase_disappearance,
    provider_topology_witness,
    provider_boundary_route
};

// `accepted_target` is the only state that says the target topology itself was
// freshly solved/reviewed and may be used. All other states preserve explicit
// uncertainty or the requirement for a new target-topology solve.
enum class PtPhaseTransitionResolution {
    accepted_target,
    target_resolve_required,
    target_resolve_failed,
    candidate_not_accepted,
    broader_topology_required,
    indeterminate
};

struct PtPhaseTransitionEdgeCapability {
    std::size_t source_phase_count{};
    std::size_t target_phase_count{};
    PtPhaseTransitionSupport support{PtPhaseTransitionSupport::detection_only};
    // Phase addition/removal is never authorized by fraction clipping. When
    // true, acceptance of the target requires solving/reviewing that target
    // topology rather than projecting the source state by deleting/adding a slot.
    bool requires_fresh_target_solve{true};
};

struct PtPhaseTransitionCapability {
    static constexpr std::string_view convention = pt_phase_transition_convention;

    std::vector<PtPhaseTransitionEdgeCapability> edges;

    [[nodiscard]] const PtPhaseTransitionEdgeCapability* edge(
        std::size_t source_phase_count,
        std::size_t target_phase_count) const noexcept {
        for (const auto& candidate : edges) {
            if (candidate.source_phase_count == source_phase_count &&
                candidate.target_phase_count == target_phase_count) {
                return &candidate;
            }
        }
        return nullptr;
    }

    [[nodiscard]] bool structurally_valid(
        std::size_t maximum_phase_count) const noexcept {
        if (maximum_phase_count == 0U) { return false; }
        for (std::size_t i = 0; i < edges.size(); ++i) {
            const auto& value = edges[i];
            if (value.source_phase_count == 0U ||
                value.target_phase_count == 0U ||
                value.source_phase_count > maximum_phase_count ||
                value.target_phase_count > maximum_phase_count ||
                value.source_phase_count == value.target_phase_count) {
                return false;
            }
            for (std::size_t j = i + 1U; j < edges.size(); ++j) {
                if (edges[j].source_phase_count == value.source_phase_count &&
                    edges[j].target_phase_count == value.target_phase_count) {
                    return false;
                }
            }
            if (value.support == PtPhaseTransitionSupport::fresh_target_resolve &&
                !value.requires_fresh_target_solve) {
                return false;
            }
        }
        return true;
    }
};

// One piece of evidence generated during the current PT solve. `target_phase_count`
// is optional because a final common-tangent instability can prove that the
// current set is incomplete without proving whether the correct replacement is
// another same-count set or a higher phase count.
struct PtPhaseTransitionEvidence {
    std::size_t source_phase_count{};
    std::optional<std::size_t> target_phase_count;
    PtPhaseTransitionTrigger trigger{
        PtPhaseTransitionTrigger::provider_topology_witness};
    PtPhaseTransitionResolution resolution{
        PtPhaseTransitionResolution::indeterminate};
    bool fresh_target_solve_attempted{false};
    bool target_topology_closed{false};
    std::string provider_evidence_profile;
    std::string diagnostic;

    [[nodiscard]] bool structurally_valid(
        const PtPhaseTransitionCapability& capability,
        std::size_t maximum_phase_count) const noexcept {
        if (source_phase_count == 0U ||
            source_phase_count > maximum_phase_count ||
            provider_evidence_profile.empty()) {
            return false;
        }

        const PtPhaseTransitionEdgeCapability* declared_edge = nullptr;
        if (target_phase_count) {
            if (*target_phase_count == 0U ||
                *target_phase_count > maximum_phase_count ||
                *target_phase_count == source_phase_count) {
                return false;
            }
            declared_edge = capability.edge(source_phase_count, *target_phase_count);
            if (declared_edge == nullptr) { return false; }
        }

        switch (resolution) {
        case PtPhaseTransitionResolution::accepted_target:
            return target_phase_count.has_value() &&
                   fresh_target_solve_attempted && target_topology_closed &&
                   declared_edge != nullptr &&
                   declared_edge->support ==
                       PtPhaseTransitionSupport::fresh_target_resolve &&
                   declared_edge->requires_fresh_target_solve;
        case PtPhaseTransitionResolution::target_resolve_required:
            return target_phase_count.has_value() &&
                   !fresh_target_solve_attempted && !target_topology_closed &&
                   declared_edge != nullptr &&
                   declared_edge->requires_fresh_target_solve;
        case PtPhaseTransitionResolution::target_resolve_failed:
            return fresh_target_solve_attempted && !target_topology_closed;
        case PtPhaseTransitionResolution::candidate_not_accepted:
            return !target_topology_closed;
        case PtPhaseTransitionResolution::broader_topology_required:
            return !target_topology_closed;
        case PtPhaseTransitionResolution::indeterminate:
            return !target_topology_closed;
        }
        return false;
    }
};

struct PtPhaseTransitionReport {
    static constexpr std::string_view convention = pt_phase_transition_convention;

    std::vector<PtPhaseTransitionEvidence> evidence;

    [[nodiscard]] bool structurally_valid(
        const PtPhaseTransitionCapability& capability,
        std::size_t maximum_phase_count) const noexcept {
        if (!capability.structurally_valid(maximum_phase_count)) { return false; }
        return std::all_of(
            evidence.begin(), evidence.end(),
            [&](const PtPhaseTransitionEvidence& item) {
                return item.structurally_valid(capability, maximum_phase_count);
            });
    }
};

} // namespace mpmc::flash

#endif // MPMC_FLASH_PT_PHASE_TRANSITION_HPP
