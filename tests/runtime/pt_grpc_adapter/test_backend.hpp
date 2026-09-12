#ifndef MPMC_TESTS_RUNTIME_PT_GRPC_ADAPTER_TEST_BACKEND_HPP
#define MPMC_TESTS_RUNTIME_PT_GRPC_ADAPTER_TEST_BACKEND_HPP

#include <mpmc/flash/pt_flash_backend.hpp>

#include <atomic>
#include <chrono>
#include <condition_variable>
#include <cstddef>
#include <mutex>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace pt_grpc_test {

namespace fl = ::mpmc::flash;

enum class BackendOutcome {
    accepted,
    indeterminate
};

inline fl::PtFlashBackendCapability make_capability(
    std::string backend_id = "golden/backend/v1") {
    fl::PtFlashBackendCapability result;
    result.backend_id = std::move(backend_id);
    result.model_profile = "golden/model-neutral-fixture/v1";
    result.algorithm_profile = "golden/scripted-no-eos/v1";
    result.publication_profile = "golden/wire-structure-only/v1";
    result.configuration_profile = "golden/configuration/v1";
    result.dataset_id = "synthetic-wire-fixture/no-physical-data";
    result.revision = "golden-r1";
    result.component_ids = {"methane", "water"};
    result.supported_phase_counts = {1U, 2U, 3U};
    result.scalar_settings = {{"fixture-setting", 2.5, "dimensionless"}};
    result.transition_capability.edges = {
        {1U, 2U, fl::PtPhaseTransitionSupport::fresh_target_resolve, true},
        {2U, 1U, fl::PtPhaseTransitionSupport::detection_only, true},
        {2U, 3U, fl::PtPhaseTransitionSupport::fresh_target_resolve, true},
        {3U, 2U, fl::PtPhaseTransitionSupport::fresh_target_resolve, true}};
    result.performs_initial_stability_search = true;
    result.performs_final_phase_set_review = true;
    result.performs_boundary_neighbor_resolve = true;
    result.global_stability_proven = false;
    result.phase_metadata_namespace = "golden/provider-phase-metadata/v1";
    return result;
}

class Backend final : public fl::PtFlashBackend {
public:
    explicit Backend(BackendOutcome outcome = BackendOutcome::accepted)
        : outcome_(outcome), capability_(make_capability()) {}

    [[nodiscard]] const fl::PtFlashBackendCapability& capability()
        const noexcept override {
        return capability_;
    }

    [[nodiscard]] fl::PtFlashBackendResult solve(
        const fl::PtFlashRequest& request) override {
        call_count.fetch_add(1U, std::memory_order_relaxed);
        {
            std::unique_lock lock(mutex_);
            started_ = true;
            condition_.notify_all();
            condition_.wait(lock, [&] { return !blocked_ || released_; });
        }

        fl::PtFlashBackendResult result;
        result.capability = capability_;
        result.solution.capability.maximum_phase_count = 3U;
        result.solution.pressure_pa = request.pressure_pa;
        result.solution.temperature_k = request.temperature_k;
        result.solution.feed = request.feed;
        result.solution.global_stability_proven = false;
        result.provider_result_convention = "golden/provider-result/v1";
        result.morphology_resolved = false;
        if (outcome_ == BackendOutcome::accepted) {
            result.transition_report.evidence.push_back(
                {3U,
                 2U,
                 fl::PtPhaseTransitionTrigger::phase_disappearance,
                 fl::PtPhaseTransitionResolution::accepted_target,
                 true,
                 true,
                 "golden/transition-evidence/v1",
                 "synthetic accepted-target evidence for wire mapping only"});
        } else {
            result.transition_report.evidence.push_back(
                {2U,
                 3U,
                 fl::PtPhaseTransitionTrigger::final_phase_set_instability,
                 fl::PtPhaseTransitionResolution::target_resolve_failed,
                 true,
                 false,
                 "golden/transition-evidence/v1",
                 "synthetic failed-target evidence for wire mapping only"});
        }

        fl::PtCandidatePhaseSet candidates;
        for (std::size_t phase_index = 0; phase_index < 2U; ++phase_index) {
            fl::PtCandidatePhase phase;
            phase.mole_phase_fraction = phase_index == 0U ? 0.25 : 0.75;
            phase.composition = request.feed;
            phase.activity.ln_phi = {-0.125, -0.25};
            phase.activity.branch = 7U + phase_index * 2U;
            phase.activity.smooth = true;
            phase.compressibility_factor = 0.8 + 0.1 *
                static_cast<double>(phase_index);
            candidates.phases.push_back(std::move(phase));
            result.phase_metadata.push_back(
                {"opaque-role-" + std::to_string(phase_index),
                 "opaque-family"});
        }
        result.solution.candidate_phase_set = std::move(candidates);
        if (outcome_ == BackendOutcome::accepted) {
            result.solution.status = fl::PtPhaseSetStatus::accepted;
            result.solution.diagnostic =
                "synthetic accepted result for cross-language wire golden";
        } else {
            result.solution.status = fl::PtPhaseSetStatus::indeterminate;
            result.solution.diagnostic =
                "synthetic indeterminate result for cross-language wire golden";
        }
        completed_count.fetch_add(1U, std::memory_order_release);
        condition_.notify_all();
        return result;
    }

    void block() {
        std::lock_guard lock(mutex_);
        blocked_ = true;
        released_ = false;
        started_ = false;
    }

    void release() {
        std::lock_guard lock(mutex_);
        released_ = true;
        condition_.notify_all();
    }

    [[nodiscard]] bool wait_until_started(
        std::chrono::milliseconds timeout = std::chrono::seconds(2)) {
        std::unique_lock lock(mutex_);
        return condition_.wait_for(lock, timeout, [&] { return started_; });
    }

    [[nodiscard]] bool wait_until_completed(
        std::size_t expected,
        std::chrono::milliseconds timeout = std::chrono::seconds(2)) {
        std::unique_lock lock(mutex_);
        return condition_.wait_for(lock, timeout, [&] {
            return completed_count.load(std::memory_order_acquire) >= expected;
        });
    }

    std::atomic<std::size_t> call_count{0U};
    std::atomic<std::size_t> completed_count{0U};

private:
    BackendOutcome outcome_;
    fl::PtFlashBackendCapability capability_;
    std::mutex mutex_;
    std::condition_variable condition_;
    bool blocked_{false};
    bool released_{false};
    bool started_{false};
};

} // namespace pt_grpc_test

#endif // MPMC_TESTS_RUNTIME_PT_GRPC_ADAPTER_TEST_BACKEND_HPP
