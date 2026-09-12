#ifndef MPMC_FLASH_PT_FLASH_BACKEND_HPP
#define MPMC_FLASH_PT_FLASH_BACKEND_HPP

#include <mpmc/flash/pt_phase_set.hpp>

#include <algorithm>
#include <cstddef>
#include <string>
#include <string_view>
#include <vector>

namespace mpmc::flash {

inline constexpr std::string_view pt_flash_backend_convention =
    "PT/flash-backend/capability-and-dispatch/v1";

// Capability belongs to one configured backend instance. The ordered component
// snapshot and dataset/revision therefore describe the exact model instance that
// will receive PtFlashRequest::feed. This contract describes implemented search
// and publication capabilities; it is not a global thermodynamic proof.
struct PtFlashBackendCapability {
    static constexpr std::string_view convention = pt_flash_backend_convention;

    std::string backend_id;
    std::string model_profile;
    std::string algorithm_profile;
    std::string publication_profile;
    std::string dataset_id;
    std::string revision;
    std::vector<std::string> component_ids;
    std::vector<std::size_t> supported_phase_counts;

    bool performs_initial_stability_search{false};
    bool performs_final_phase_set_review{false};
    bool performs_boundary_neighbor_resolve{false};
    bool global_stability_proven{false};

    // Empty means that the generic backend result carries no provider-specific
    // phase metadata. A non-empty namespace makes role_id/family_id opaque
    // provider identifiers; generic code must not reinterpret them.
    std::string phase_metadata_namespace;

    [[nodiscard]] bool supports_phase_count(std::size_t count) const noexcept {
        return std::find(supported_phase_counts.begin(),
                         supported_phase_counts.end(), count) !=
               supported_phase_counts.end();
    }

    [[nodiscard]] std::size_t maximum_phase_count() const noexcept {
        std::size_t maximum = 0U;
        for (const std::size_t count : supported_phase_counts) {
            maximum = std::max(maximum, count);
        }
        return maximum;
    }

    [[nodiscard]] bool structurally_valid() const noexcept {
        if (backend_id.empty() || model_profile.empty() ||
            algorithm_profile.empty() || publication_profile.empty() ||
            dataset_id.empty() || revision.empty() || component_ids.empty() ||
            supported_phase_counts.empty()) {
            return false;
        }
        for (std::size_t i = 0; i < component_ids.size(); ++i) {
            if (component_ids[i].empty()) { return false; }
            for (std::size_t j = i + 1U; j < component_ids.size(); ++j) {
                if (component_ids[i] == component_ids[j]) { return false; }
            }
        }
        for (std::size_t i = 0; i < supported_phase_counts.size(); ++i) {
            if (supported_phase_counts[i] == 0U) { return false; }
            for (std::size_t j = i + 1U; j < supported_phase_counts.size(); ++j) {
                if (supported_phase_counts[i] == supported_phase_counts[j]) {
                    return false;
                }
            }
        }
        return maximum_phase_count() > 0U;
    }
};

// Model-neutral PT request. Component identity/order is owned by capability();
// the backend remains responsible for the exact validation rules of its existing
// solver. This layer never normalizes, clips or repairs the feed.
struct PtFlashRequest {
    double pressure_pa{};
    double temperature_k{};
    std::vector<double> feed;
};

// Optional provider metadata associated with the phase in the same vector slot.
// role_id/family_id have meaning only inside capability.phase_metadata_namespace.
// They are not portable liquid/vapor/aqueous classifiers.
struct PtFlashBackendPhaseMetadata {
    std::string role_id;
    std::string family_id;
};

struct PtFlashBackendResult {
    static constexpr std::string_view convention =
        "PT/flash-backend/result/v1";

    PtFlashBackendCapability capability;
    PtPhaseSetResult solution;
    std::string provider_result_convention;
    std::vector<PtFlashBackendPhaseMetadata> phase_metadata;
    bool morphology_resolved{false};

    [[nodiscard]] bool structurally_valid() const noexcept {
        if (!capability.structurally_valid() ||
            provider_result_convention.empty() ||
            solution.capability.maximum_phase_count !=
                capability.maximum_phase_count() ||
            solution.feed.size() != capability.component_ids.size() ||
            (solution.global_stability_proven &&
             !capability.global_stability_proven)) {
            return false;
        }

        if (solution.candidate_phase_set) {
            const std::size_t count = solution.candidate_phase_set->phases.size();
            if (count > capability.maximum_phase_count()) { return false; }
            if (!phase_metadata.empty() && phase_metadata.size() != count) {
                return false;
            }
        } else if (!phase_metadata.empty()) {
            return false;
        }

        if (capability.phase_metadata_namespace.empty() &&
            !phase_metadata.empty()) {
            return false;
        }

        if (const auto* accepted = solution.accepted_phase_set()) {
            if (!capability.supports_phase_count(accepted->phases.size())) {
                return false;
            }
            if (!capability.phase_metadata_namespace.empty() &&
                phase_metadata.size() != accepted->phases.size()) {
                return false;
            }
        }
        return true;
    }

    [[nodiscard]] const PtCandidatePhaseSet* accepted_phase_set() const & noexcept {
        if (!structurally_valid()) { return nullptr; }
        return solution.accepted_phase_set();
    }
    const PtCandidatePhaseSet* accepted_phase_set() const && = delete;
};

// Runtime polymorphism is intentionally confined to this coarse solve boundary.
// Concrete backends immediately enter their existing strongly typed solver;
// fugacity/property/stability inner loops remain free of virtual dispatch.
class PtFlashBackend {
public:
    virtual ~PtFlashBackend() = default;

    [[nodiscard]] virtual const PtFlashBackendCapability& capability()
        const noexcept = 0;
    [[nodiscard]] virtual PtFlashBackendResult solve(
        const PtFlashRequest& request) = 0;
};

} // namespace mpmc::flash

#endif // MPMC_FLASH_PT_FLASH_BACKEND_HPP
