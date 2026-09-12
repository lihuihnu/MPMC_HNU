#ifndef MPMC_FLASH_PR76_PT_FLASH_BACKEND_HPP
#define MPMC_FLASH_PR76_PT_FLASH_BACKEND_HPP

#include <mpmc/flash/pr76_max3_phase_set.hpp>
#include <mpmc/flash/pt_flash_backend.hpp>

#include <algorithm>
#include <optional>
#include <span>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

namespace mpmc::flash {

inline constexpr std::string_view pr76_pt_flash_backend_id =
    "PR76/PT/max3/phase-set-backend/v2";
inline constexpr std::string_view pr76_pt_flash_backend_configuration_profile =
    "PR76/PT/max3/backend-configuration/v2";
inline constexpr std::string_view pr76_pt_flash_backend_result_convention =
    "PR76/PT/max3/phase-set-backend-result/v2";
inline constexpr std::string_view pr76_pt_transition_evidence_profile =
    "PR76/PT/max3/phase-transition-evidence/v2";

struct Pr76PtFlashBackendOptions {
    PtSplitOptions split;
    PtThreePhaseOptions three_phase;
    StabilityOptions final_three_phase_stability;
    std::size_t max_three_phase_attempts{16};
    double new_phase_seed_fraction{0.1};
    std::vector<Pr76PtThreePhaseStart> three_phase_starts;
    std::size_t max_three_phase_starts{16};
    std::size_t max_three_phase_start_entries{12288};
    std::vector<std::vector<double>> initial_starts;
    std::vector<std::vector<double>> final_starts;
};

// Legacy VLE-only evidence projection remains source-compatible for callers that
// intentionally use solve_pr76_pt_vle directly.
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
            "final common-tangent review found lower-Gibbs evidence; VLE-only routing does not determine the replacement phase count");
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
                "a material-balanced/fugacity-converged two-phase attempt reached the phase-fraction endpoint; no phase is deleted without a fresh one-phase decision");
        }
        if (source.candidate() != nullptr) {
            add(1U, 2U,
                PtPhaseTransitionTrigger::initial_stability_witness,
                PtPhaseTransitionResolution::candidate_not_accepted,
                true, false,
                "a two-phase candidate exists but the current PR76 solve did not accept that topology");
        }
    }
    return report;
}

[[nodiscard]] inline PtPhaseTransitionReport project_pr76_pt_max3_transition_report(
    const Pr76PtMax3Result& source) {
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

    const auto& base = source.base.solution;
    if (base.status == PtSplitStatus::two_phase_no_instability_found) {
        add(1U, 2U, PtPhaseTransitionTrigger::initial_stability_witness,
            PtPhaseTransitionResolution::accepted_target,
            true, true,
            "negative feed-stability evidence was followed by a fresh two-phase solve and final review");
        return report;
    }
    if (base.status == PtSplitStatus::indeterminate) {
        return project_pr76_pt_transition_report(base);
    }
    if (base.status != PtSplitStatus::phase_set_unstable) {
        return report;
    }

    add(1U, 2U, PtPhaseTransitionTrigger::initial_stability_witness,
        PtPhaseTransitionResolution::candidate_not_accepted,
        true, false,
        "the first two-phase candidate failed its final common-tangent review");

    switch (source.status) {
    case Pr76PtMax3Status::three_phase:
        add(2U, 3U, PtPhaseTransitionTrigger::final_phase_set_instability,
            PtPhaseTransitionResolution::accepted_target,
            true, true,
            "negative final-stability evidence triggered a fresh three-phase solve and the resulting common tangent passed final stability review");
        break;
    case Pr76PtMax3Status::two_phase:
        if (source.two_phase_neighbor() != nullptr) {
            add(2U, 3U, PtPhaseTransitionTrigger::final_phase_set_instability,
                PtPhaseTransitionResolution::candidate_not_accepted,
                true, false,
                "additional-phase evidence was followed by a fresh three-phase solve that reached a phase-disappearance boundary");
            add(3U, 2U, PtPhaseTransitionTrigger::phase_disappearance,
                PtPhaseTransitionResolution::accepted_target,
                true, true,
                "the surviving pair was fresh-resolved through the full PR76 two-phase stability/split/final-review path");
        }
        break;
    case Pr76PtMax3Status::phase_boundary_unresolved:
        add(2U, 3U, PtPhaseTransitionTrigger::final_phase_set_instability,
            PtPhaseTransitionResolution::candidate_not_accepted,
            true, false,
            "additional-phase evidence reached a three-phase disappearance boundary");
        add(3U, 2U, PtPhaseTransitionTrigger::phase_disappearance,
            PtPhaseTransitionResolution::target_resolve_failed,
            true, false,
            source.diagnostic);
        break;
    case Pr76PtMax3Status::higher_phase_count_or_wrong_candidate:
        add(3U, std::nullopt,
            PtPhaseTransitionTrigger::final_phase_set_instability,
            PtPhaseTransitionResolution::broader_topology_required,
            true, false,
            source.diagnostic);
        break;
    case Pr76PtMax3Status::indeterminate:
        add(2U, 3U, PtPhaseTransitionTrigger::final_phase_set_instability,
            PtPhaseTransitionResolution::indeterminate,
            !source.attempts.empty(), false,
            source.diagnostic);
        break;
    case Pr76PtMax3Status::single_phase:
        break;
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
        Pr76PtMax3Options flash_options;
        flash_options.two_phase = options_.split;
        flash_options.three_phase = options_.three_phase;
        flash_options.final_three_phase_stability =
            options_.final_three_phase_stability;
        flash_options.max_three_phase_attempts = options_.max_three_phase_attempts;
        flash_options.new_phase_seed_fraction = options_.new_phase_seed_fraction;
        flash_options.three_phase_starts = options_.three_phase_starts;
        flash_options.max_three_phase_starts = options_.max_three_phase_starts;
        flash_options.max_three_phase_start_entries =
            options_.max_three_phase_start_entries;

        const auto source = solve_pr76_pt_max3(
            request.pressure_pa, request.temperature_k, request.feed,
            evaluator_, flash_options,
            std::span<const std::vector<double>>{options_.initial_starts},
            std::span<const std::vector<double>>{options_.final_starts});
        const auto published = project_pr76_pt_max3_phase_set(source);

        PtFlashBackendResult result;
        result.capability = capability_;
        result.solution = published.solution;
        result.transition_report =
            project_pr76_pt_max3_transition_report(source);
        result.provider_result_convention = published.publication_convention;
        result.morphology_resolved = false;

        if (published.dataset_id != capability_.dataset_id ||
            published.revision != capability_.revision ||
            published.component_ids != capability_.component_ids ||
            published.model_profile != capability_.model_profile ||
            published.orchestration_convention != capability_.algorithm_profile ||
            published.publication_convention != capability_.publication_profile ||
            std::string_view{published.phase_convention} !=
                thermodynamics::pr76_pt_convention) {
            reject_adapter_result(
                result, "PR76 max3 backend: provider/model provenance mismatch");
            return result;
        }

        if (!result.structurally_valid()) {
            reject_adapter_result(
                result, "PR76 max3 backend: generic result integrity guard failed");
        }
        return result;
    }

private:
    [[nodiscard]] static PtFlashBackendCapability build_capability(
        const Pr76VleEvaluator& evaluator) {
        PtFlashBackendCapability capability;
        capability.backend_id = std::string(pr76_pt_flash_backend_id);
        capability.model_profile = std::string(thermodynamics::pr76_profile);
        capability.algorithm_profile = pr76_pt_max3_convention;
        capability.publication_profile =
            pr76_pt_max3_phase_set_publication_convention;
        capability.configuration_profile =
            std::string(pr76_pt_flash_backend_configuration_profile);
        const auto& parameters = evaluator.model().parameters();
        capability.dataset_id = parameters.dataset_id();
        capability.revision = parameters.revision();
        for (const auto& component : parameters.components().items()) {
            capability.component_ids.push_back(component.id);
        }
        capability.supported_phase_counts = {1U, 2U, 3U};
        capability.transition_capability.edges = {
            {1U, 2U, PtPhaseTransitionSupport::fresh_target_resolve, true},
            {2U, 1U, PtPhaseTransitionSupport::detection_only, true},
            {2U, 3U, PtPhaseTransitionSupport::fresh_target_resolve, true},
            {3U, 2U, PtPhaseTransitionSupport::fresh_target_resolve, true}};
        capability.performs_initial_stability_search = true;
        capability.performs_final_phase_set_review = true;
        capability.performs_boundary_neighbor_resolve = true;
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
