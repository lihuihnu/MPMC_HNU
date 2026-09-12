#ifndef MPMC_RUNTIME_PT_SERVICE_CONTRACT_HPP
#define MPMC_RUNTIME_PT_SERVICE_CONTRACT_HPP

#include <mpmc/flash/pt_flash_backend.hpp>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <optional>
#include <string>
#include <string_view>
#include <vector>

namespace mpmc::runtime {

inline constexpr std::string_view pt_service_boundary_convention =
    "PT/service-boundary/v1";
inline constexpr std::string_view pt_service_component_inventory_convention =
    "PT/service-component-inventory/v1";
inline constexpr std::string_view pt_service_capability_convention =
    "PT/service-capability-discovery/v1";
inline constexpr std::string_view pt_service_request_convention =
    "PT/service-request/v1";
inline constexpr std::string_view pt_service_result_convention =
    "PT/service-result/v1";

namespace detail {

[[nodiscard]] inline bool pt_service_identifier_valid(
    std::string_view value) noexcept {
    return !value.empty() &&
           std::all_of(value.begin(), value.end(), [](unsigned char character) {
               return character >= 33U && character <= 126U;
           });
}

} // namespace detail

struct PtRuntimeComponent {
    std::string component_id;
    std::size_t feed_index{};
};

// This is the exact ordered component snapshot accepted by one configured
// backend. IDs are stable identities, not display labels or inferred chemistry.
struct PtRuntimeComponentInventory {
    static constexpr std::string_view convention =
        pt_service_component_inventory_convention;

    std::vector<PtRuntimeComponent> components;

    [[nodiscard]] bool structurally_valid() const noexcept {
        if (components.empty()) { return false; }
        for (std::size_t i = 0; i < components.size(); ++i) {
            if (!detail::pt_service_identifier_valid(
                    components[i].component_id) ||
                components[i].feed_index != i) {
                return false;
            }
            for (std::size_t j = i + 1U; j < components.size(); ++j) {
                if (components[i].component_id == components[j].component_id) {
                    return false;
                }
            }
        }
        return true;
    }
};

// Discovery returns an owning snapshot. The frozen model-neutral backend
// capability remains authoritative; the inventory only makes its order explicit
// for service clients and never invents names or chemical metadata.
struct PtServiceBackendDescriptor {
    static constexpr std::string_view convention =
        pt_service_capability_convention;

    // Service-local stable selector for one configured backend instance. This
    // is distinct from capability.backend_id, which identifies the adapter
    // implementation and may be shared by multiple configured instances.
    std::string configured_backend_id;
    flash::PtFlashBackendCapability capability;
    PtRuntimeComponentInventory component_inventory;

    [[nodiscard]] bool structurally_valid() const noexcept {
        if (!detail::pt_service_identifier_valid(configured_backend_id) ||
            !detail::pt_service_identifier_valid(capability.backend_id) ||
            !capability.structurally_valid() ||
            !component_inventory.structurally_valid() ||
            capability.component_ids.size() !=
                component_inventory.components.size()) {
            return false;
        }
        for (std::size_t i = 0; i < capability.component_ids.size(); ++i) {
            if (capability.component_ids[i] !=
                component_inventory.components[i].component_id) {
                return false;
            }
        }
        return true;
    }
};

struct PtServiceCompositionEntry {
    std::string component_id;
    double mole_fraction{};
};

// Request components are keyed by stable ID and may arrive in any order. The
// service maps them to the discovered backend order without normalizing,
// clipping, filling, or otherwise repairing mole fractions.
struct PtServiceRequest {
    static constexpr std::string_view convention = pt_service_request_convention;

    std::string configured_backend_id;
    double pressure_pa{};
    double temperature_k{};
    std::vector<PtServiceCompositionEntry> feed;
};

struct PtServicePhaseComponent {
    std::string component_id;
    double mole_fraction{};
    double ln_fugacity_coefficient{};
};

struct PtServicePhaseMetadata {
    std::string role_id;
    std::string family_id;
};

// Phase index and provider branch are diagnostics, not portable physical phase
// identities. Provider metadata is meaningful only in the namespace declared by
// the discovered backend capability.
struct PtServicePhase {
    std::size_t phase_index{};
    double mole_phase_fraction{};
    std::vector<PtServicePhaseComponent> components;
    std::size_t provider_branch{};
    bool provider_branch_smooth{true};
    std::optional<double> compressibility_factor;
    std::optional<PtServicePhaseMetadata> provider_metadata;
};

struct PtServiceResultProvenance {
    PtServiceBackendDescriptor backend;
    std::string backend_result_convention;
    std::string phase_set_convention;
    std::string phase_transition_convention;
    std::string provider_result_convention;

    [[nodiscard]] bool structurally_valid() const noexcept {
        return backend.structurally_valid() &&
               backend_result_convention ==
                   flash::PtFlashBackendResult::convention &&
               phase_set_convention == flash::PtPhaseSetResult::convention &&
               phase_transition_convention ==
                   flash::PtPhaseTransitionReport::convention &&
               !provider_result_convention.empty();
    }
};

namespace detail {

inline constexpr double pt_service_composition_roundoff =
    64.0 * std::numeric_limits<double>::epsilon();

inline void pt_service_add(double value, double& sum,
                           double& correction) noexcept {
    const double increment = value - correction;
    const double next = sum + increment;
    correction = (next - sum) - increment;
    sum = next;
}

template <typename Entry, typename Value>
[[nodiscard]] inline bool pt_service_normalized_entries(
    const std::vector<Entry>& entries, Value&& value) noexcept {
    if (entries.empty()) { return false; }
    double sum = 0.0;
    double correction = 0.0;
    for (const auto& entry : entries) {
        const double fraction = value(entry);
        if (!std::isfinite(fraction) || fraction < 0.0 || fraction > 1.0) {
            return false;
        }
        pt_service_add(fraction, sum, correction);
    }
    return std::abs(sum - 1.0) <= pt_service_composition_roundoff;
}

} // namespace detail

struct PtServiceComputationResult {
    PtServiceResultProvenance provenance;
    double pressure_pa{};
    double temperature_k{};
    std::vector<PtServiceCompositionEntry> feed;
    // Accepted phases only. Non-accepted backend candidates are deliberately not
    // promoted into this vector.
    std::vector<PtServicePhase> phases;
    flash::PtPhaseTransitionReport transition_report;
    bool global_stability_proven{false};
    bool morphology_resolved{false};
    std::string diagnostic;

    [[nodiscard]] bool structurally_valid() const noexcept {
        if (!provenance.structurally_valid() || !std::isfinite(pressure_pa) ||
            pressure_pa <= 0.0 || !std::isfinite(temperature_k) ||
            temperature_k <= 0.0 ||
            feed.size() !=
                provenance.backend.component_inventory.components.size() ||
            !detail::pt_service_normalized_entries(
                feed, [](const auto& entry) { return entry.mole_fraction; }) ||
            !transition_report.structurally_valid(
                provenance.backend.capability.transition_capability,
                provenance.backend.capability.maximum_phase_count()) ||
            (global_stability_proven &&
             !provenance.backend.capability.global_stability_proven) ||
            (morphology_resolved &&
             provenance.backend.capability.phase_metadata_namespace.empty())) {
            return false;
        }

        for (std::size_t i = 0; i < feed.size(); ++i) {
            if (feed[i].component_id != provenance.backend.component_inventory
                                                  .components[i]
                                                  .component_id) {
                return false;
            }
        }

        if (phases.empty()) { return true; }
        if (!provenance.backend.capability.supports_phase_count(phases.size()) ||
            !detail::pt_service_normalized_entries(
                phases,
                [](const auto& phase) { return phase.mole_phase_fraction; })) {
            return false;
        }

        const bool metadata_required =
            !provenance.backend.capability.phase_metadata_namespace.empty();
        for (std::size_t phase_index = 0; phase_index < phases.size();
             ++phase_index) {
            const auto& phase = phases[phase_index];
            if (phase.phase_index != phase_index ||
                phase.components.size() != feed.size() ||
                !detail::pt_service_normalized_entries(
                    phase.components,
                    [](const auto& component) {
                        return component.mole_fraction;
                    }) ||
                (phase.compressibility_factor.has_value() &&
                 (!std::isfinite(*phase.compressibility_factor) ||
                  *phase.compressibility_factor <= 0.0)) ||
                phase.provider_metadata.has_value() != metadata_required ||
                (phase.provider_metadata.has_value() &&
                 (phase.provider_metadata->role_id.empty() ||
                  phase.provider_metadata->family_id.empty()))) {
                return false;
            }
            for (std::size_t component_index = 0;
                 component_index < phase.components.size(); ++component_index) {
                const auto& component = phase.components[component_index];
                if (component.component_id !=
                        feed[component_index].component_id ||
                    !std::isfinite(component.ln_fugacity_coefficient)) {
                    return false;
                }
            }
        }
        return true;
    }
};

// `indeterminate` and `phase_set_unstable` are scientific computation outcomes,
// not service failures. Only `error` carries PtServiceError.
enum class PtServiceOutcome {
    accepted,
    phase_set_unstable,
    indeterminate,
    error
};

enum class PtServiceErrorCode {
    invalid_configured_backend_id,
    configured_backend_not_found,
    invalid_pressure,
    invalid_temperature,
    component_count_mismatch,
    invalid_component_id,
    duplicate_component,
    unknown_component,
    invalid_mole_fraction,
    composition_not_normalized,
    backend_rejected_request,
    backend_contract_violation,
    backend_execution_failure
};

struct PtServiceError {
    PtServiceErrorCode code{PtServiceErrorCode::backend_contract_violation};
    std::string field;
    std::string diagnostic;

    [[nodiscard]] bool structurally_valid() const noexcept {
        switch (code) {
        case PtServiceErrorCode::invalid_configured_backend_id:
        case PtServiceErrorCode::configured_backend_not_found:
        case PtServiceErrorCode::invalid_pressure:
        case PtServiceErrorCode::invalid_temperature:
        case PtServiceErrorCode::component_count_mismatch:
        case PtServiceErrorCode::invalid_component_id:
        case PtServiceErrorCode::duplicate_component:
        case PtServiceErrorCode::unknown_component:
        case PtServiceErrorCode::invalid_mole_fraction:
        case PtServiceErrorCode::composition_not_normalized:
        case PtServiceErrorCode::backend_rejected_request:
        case PtServiceErrorCode::backend_contract_violation:
        case PtServiceErrorCode::backend_execution_failure:
            return !field.empty() && !diagnostic.empty();
        }
        return false;
    }
};

struct PtServiceResponse {
    static constexpr std::string_view convention = pt_service_result_convention;

    PtServiceOutcome outcome{PtServiceOutcome::error};
    std::optional<PtServiceComputationResult> result;
    std::optional<PtServiceError> error;

    [[nodiscard]] bool structurally_valid() const noexcept {
        switch (outcome) {
        case PtServiceOutcome::error:
            return !result.has_value() && error.has_value() &&
                   error->structurally_valid();
        case PtServiceOutcome::accepted:
            return result.has_value() && !error.has_value() &&
                   result->structurally_valid() && !result->phases.empty();
        case PtServiceOutcome::phase_set_unstable:
        case PtServiceOutcome::indeterminate:
            return result.has_value() && !error.has_value() &&
                   result->structurally_valid() && result->phases.empty();
        }
        return false;
    }

    [[nodiscard]] std::size_t accepted_phase_count() const noexcept {
        return structurally_valid() && outcome == PtServiceOutcome::accepted
                   ? result->phases.size()
                   : 0U;
    }
};

} // namespace mpmc::runtime

#endif // MPMC_RUNTIME_PT_SERVICE_CONTRACT_HPP
