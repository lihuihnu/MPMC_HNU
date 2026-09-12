#ifndef MPMC_FLASH_PR76_PT_FLASH_BACKEND_HPP
#define MPMC_FLASH_PR76_PT_FLASH_BACKEND_HPP

#include <mpmc/flash/pr76_phase_set.hpp>
#include <mpmc/flash/pt_flash_backend.hpp>

#include <algorithm>
#include <span>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

namespace mpmc::flash {

inline constexpr std::string_view pr76_pt_flash_backend_id =
    "PR76/PT-VLE/phase-set-backend/v1";
inline constexpr std::string_view pr76_pt_flash_backend_configuration_profile =
    "PR76/PT-VLE/backend-configuration/v1";
inline constexpr std::string_view pr76_pt_flash_backend_result_convention =
    "PR76/PT-VLE/phase-set-backend-result/v1";
inline constexpr std::string_view pr76_pt_transition_evidence_profile =
    "PR76/PT-VLE/phase-transition-evidence/v1";

struct Pr76PtFlashBackendOptions {
    PtSplitOptions split;
    std::vector<std::vector<double>> initial_starts;
    std::vector<std::vector<double>> final_starts;
};

// Pure projection of already-owned PR76/VLE solve evidence. It does not rerun
// stability, RR, EOS properties, phase split, or final review.
[[nodiscard]] inline PtPhaseTransitionReport project_pr76_pt_transition_report(
    const PtSplitResult& source) {
    PtPhaseTransitionReport report;
    const auto add = [&report](
        std::size_t from, std::optional<std::size_t> to,
        PtPhaseTransitionTrigger trigger,
        PtPhaseTransitionResolution resolution,
        bool fresh_attempted, bool target_closed,
        std::string diagnostic) {
        report.evidence.push_back({
            from, to, trigger, resolution, fresh_attempted, target_closed,
            std::string(pr76_pt_transition_evidence_profile),
            std::move(diagnostic)});
    };

    if (source.status == PtSplitStatus::two_phase_no_instability_found) {
        add(1U, 2U, PtPhaseTransitionTrigger::initial_stability_witness,
            PtPhaseTransitionResolution::accepted_target,
            true, true,
            "negative feed-stability evidence was followed by a fresh two-phase solve and final phase-set review");
    } else if (source.status == PtSplitStatus::phase_set_unstable) {
        add(1U, 2U, PtPhaseTransitionTrigger::initial_stability_witness,
            PtPhaseTransitionResolution::candidate_not_accepted,
            true, false,
            "a two-phase candidate was solved but failed the final phase-set review");
        add(2U, std::nullopt,
            PtPhaseTransitionTrigger::final_phase_set_instability,
            PtPhaseTransitionResolution::broader_topology_required,
            false, false,
            "final common-tangent review found a lower-Gibbs trial; the correct replacement phase count is not proven by the PR76 VLE backend");
    } else if (source.status == PtSplitStatus::indeterminate) {
        const bool disappearance = std::any_of(
            source.attempts.begin(), source.attempts.end(),
            [](const PtSplitAttempt& attempt) {
                return attempt.status == PtSplitAttemptStatus::phase_disappearance;
            });
        if (disappearance) {
            add(2U, 1U, PtPhaseTransitionTrigger::phase_disappearance,
                PtPhaseTransitionResolution::target_resolve_required,
                false, false,
                "a material-balanced/fugacity-converged two-phase attempt reached the phase-fraction endpoint; PR76 does not delete the phase and requires a fresh one-phase decision");
        }
        if (source.candidate() != nullptr) {
            add(1U, 2U,
                PtPhaseTransitionTrigger::initial_stability_witness,
                PtPhaseTransitionResolution::candidate_not_accepted,
                true, false,
                "a two-phase candidate exists but the current PR76 solve did not accept that topology");
            if (source.final_stability &&
                source.final_stability->status == StabilityStatus::indeterminate) {
                add(2U, std::nullopt,
                    PtPhaseTransitionTrigger::final_phase_set_instability,
                    PtPhaseTransitionResolution::indeterminate,
                    false, false,
                    "two-phase equations converged but final phase-set stability remained indeterminate");
            }
        }
    }
    return report;
}

class Pr76PtFlashBackend final : public PtFlashBackend {
public:
    explicit Pr76PtFlashBackend(
        Pr76VleEvaluator& evaluator,
        Pr76PtFlashBackendOptions options = {})
        : evaluator_(evaluator), options_(std::move(options)),
          capability_(build_capability(evaluator_)) {}

    [[nodiscard]] const PtFlashBackendCapability& capability()
        const noexcept override {
        return capability_;
    }

    [[nodiscard]] PtFlashBackendResult solve(
        const PtFlashRequest& request) override {
        // Run the established PR76 VLE path exactly once. The legacy result is
        // retained here only so the unified layer can project transition evidence;
        // phase-set publication remains the existing pure projection.
        const auto source = solve_pr76_pt_vle(
            request.pressure_pa, request.temperature_k, request.feed,
            evaluator_, options_.split,
            std::span<const std::vector<double>>{options_.initial_starts},
            std::span<const std::vector<double>>{options_.final_starts});
        const auto published = project_pr76_pt_phase_set(source);

        PtFlashBackendResult result;
        result.capability = capability_;
        result.solution = published.solution;
        result.transition_report = project_pr76_pt_transition_report(source.solution);
        result.provider_result_convention =
            std::string(pr76_pt_flash_backend_result_convention);
        result.morphology_resolved = false;

        if (source.dataset_id != capability_.dataset_id ||
            source.revision != capability_.revision ||
            source.component_ids != capability_.component_ids ||
            source.model_profile != capability_.model_profile ||
            std::string_view{source.phase_convention} !=
                thermodynamics::pr76_pt_convention) {
            reject_adapter_result(
                result, "PR76 backend: provider/model provenance mismatch");
            return result;
        }

        if (!result.structurally_valid()) {
            reject_adapter_result(
                result, "PR76 backend: generic result integrity guard failed");
        }
        return result;
    }

private:
    [[nodiscard]] static PtFlashBackendCapability build_capability(
        const Pr76VleEvaluator& evaluator) {
        PtFlashBackendCapability capability;
        capability.backend_id = std::string(pr76_pt_flash_backend_id);
        capability.model_profile = std::string(thermodynamics::pr76_profile);
        capability.algorithm_profile = PtSplitResult::convention;
        capability.publication_profile =
            std::string(PtPhaseSetResult::convention);
        capability.configuration_profile =
            std::string(pr76_pt_flash_backend_configuration_profile);
        const auto& parameters = evaluator.model().parameters();
        capability.dataset_id = parameters.dataset_id();
        capability.revision = parameters.revision();
        for (const auto& component : parameters.components().items()) {
            capability.component_ids.push_back(component.id);
        }
        capability.supported_phase_counts = {1U, 2U};
        capability.transition_capability.edges = {
            {1U, 2U, PtPhaseTransitionSupport::fresh_target_resolve, true},
            {2U, 1U, PtPhaseTransitionSupport::detection_only, true}};
        capability.performs_initial_stability_search = true;
        capability.performs_final_phase_set_review = true;
        capability.performs_boundary_neighbor_resolve = false;
        capability.global_stability_proven = false;
        return capability;
    }

    static void reject_adapter_result(
        PtFlashBackendResult& result, std::string diagnostic) {
        result.solution.status = PtPhaseSetStatus::indeterminate;
        result.solution.candidate_phase_set.reset();
        result.solution.global_stability_proven = false;
        result.phase_metadata.clear();
        result.transition_report.evidence.clear();
        result.solution.diagnostic = std::move(diagnostic);
    }

    Pr76VleEvaluator& evaluator_;
    Pr76PtFlashBackendOptions options_;
    PtFlashBackendCapability capability_;
};

} // namespace mpmc::flash

#endif // MPMC_FLASH_PR76_PT_FLASH_BACKEND_HPP
