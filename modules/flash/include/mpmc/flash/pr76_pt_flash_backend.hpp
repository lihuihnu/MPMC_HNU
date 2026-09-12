#ifndef MPMC_FLASH_PR76_PT_FLASH_BACKEND_HPP
#define MPMC_FLASH_PR76_PT_FLASH_BACKEND_HPP

#include <mpmc/flash/pr76_phase_set.hpp>
#include <mpmc/flash/pt_flash_backend.hpp>

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

struct Pr76PtFlashBackendOptions {
    PtSplitOptions split;
    std::vector<std::vector<double>> initial_starts;
    std::vector<std::vector<double>> final_starts;
};

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
        const auto source = solve_pr76_pt_phase_set(
            request.pressure_pa, request.temperature_k, request.feed,
            evaluator_, options_.split,
            std::span<const std::vector<double>>{options_.initial_starts},
            std::span<const std::vector<double>>{options_.final_starts});

        PtFlashBackendResult result;
        result.capability = capability_;
        result.solution = source.solution;
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
        result.solution.diagnostic = std::move(diagnostic);
    }

    Pr76VleEvaluator& evaluator_;
    Pr76PtFlashBackendOptions options_;
    PtFlashBackendCapability capability_;
};

} // namespace mpmc::flash

#endif // MPMC_FLASH_PR76_PT_FLASH_BACKEND_HPP
