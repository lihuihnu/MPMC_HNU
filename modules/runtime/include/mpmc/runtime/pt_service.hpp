#ifndef MPMC_RUNTIME_PT_SERVICE_HPP
#define MPMC_RUNTIME_PT_SERVICE_HPP

#include <mpmc/runtime/pt_service_contract.hpp>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <exception>
#include <limits>
#include <new>
#include <span>
#include <stdexcept>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

namespace mpmc::runtime {

struct PtServiceLimits {
    std::size_t max_backends{64U};
    std::size_t max_components_per_backend{256U};
    std::size_t max_identifier_bytes{256U};

    [[nodiscard]] bool structurally_valid() const noexcept {
        return max_backends > 0U && max_components_per_backend > 0U &&
               max_identifier_bytes > 0U;
    }
};

struct PtServiceBackendRegistration {
    std::string configured_backend_id;
    flash::PtFlashBackend* backend{};
};

enum class PtServiceConfigurationErrorCode {
    invalid_limits,
    backend_limit,
    null_backend,
    invalid_registration_id,
    invalid_backend_capability,
    component_limit,
    identifier_limit,
    duplicate_configured_backend_id
};

class PtServiceConfigurationError final : public std::invalid_argument {
public:
    PtServiceConfigurationError(PtServiceConfigurationErrorCode code,
                                std::string message)
        : std::invalid_argument(std::move(message)), code_(code) {}

    [[nodiscard]] PtServiceConfigurationErrorCode code() const noexcept {
        return code_;
    }

private:
    PtServiceConfigurationErrorCode code_;
};

namespace detail {

[[nodiscard]] inline bool same_scalar_settings(
    const std::vector<flash::PtFlashBackendScalarSetting>& left,
    const std::vector<flash::PtFlashBackendScalarSetting>& right) noexcept {
    if (left.size() != right.size()) { return false; }
    for (std::size_t i = 0; i < left.size(); ++i) {
        if (left[i].id != right[i].id || left[i].value != right[i].value ||
            left[i].unit != right[i].unit) {
            return false;
        }
    }
    return true;
}

[[nodiscard]] inline bool same_transition_capability(
    const flash::PtPhaseTransitionCapability& left,
    const flash::PtPhaseTransitionCapability& right) noexcept {
    if (left.edges.size() != right.edges.size()) { return false; }
    for (std::size_t i = 0; i < left.edges.size(); ++i) {
        const auto& a = left.edges[i];
        const auto& b = right.edges[i];
        if (a.source_phase_count != b.source_phase_count ||
            a.target_phase_count != b.target_phase_count ||
            a.support != b.support ||
            a.requires_fresh_target_solve != b.requires_fresh_target_solve) {
            return false;
        }
    }
    return true;
}

[[nodiscard]] inline bool same_backend_capability(
    const flash::PtFlashBackendCapability& left,
    const flash::PtFlashBackendCapability& right) noexcept {
    return left.backend_id == right.backend_id &&
           left.model_profile == right.model_profile &&
           left.algorithm_profile == right.algorithm_profile &&
           left.publication_profile == right.publication_profile &&
           left.configuration_profile == right.configuration_profile &&
           left.dataset_id == right.dataset_id && left.revision == right.revision &&
           left.component_ids == right.component_ids &&
           left.supported_phase_counts == right.supported_phase_counts &&
           same_scalar_settings(left.scalar_settings, right.scalar_settings) &&
           same_transition_capability(left.transition_capability,
                                      right.transition_capability) &&
           left.performs_initial_stability_search ==
               right.performs_initial_stability_search &&
           left.performs_final_phase_set_review ==
               right.performs_final_phase_set_review &&
           left.performs_boundary_neighbor_resolve ==
               right.performs_boundary_neighbor_resolve &&
           left.global_stability_proven == right.global_stability_proven &&
           left.phase_metadata_namespace == right.phase_metadata_namespace;
}

[[nodiscard]] inline PtServiceBackendDescriptor make_descriptor(
    std::string configured_backend_id,
    const flash::PtFlashBackendCapability& capability) {
    PtServiceBackendDescriptor descriptor;
    descriptor.configured_backend_id = std::move(configured_backend_id);
    descriptor.capability = capability;
    descriptor.component_inventory.components.reserve(
        capability.component_ids.size());
    for (std::size_t i = 0; i < capability.component_ids.size(); ++i) {
        descriptor.component_inventory.components.push_back(
            {capability.component_ids[i], i});
    }
    return descriptor;
}

[[nodiscard]] inline bool same_composition_within_roundoff(
    std::span<const double> left, std::span<const double> right) noexcept {
    if (left.size() != right.size()) { return false; }
    for (std::size_t i = 0; i < left.size(); ++i) {
        const double scale = std::max({1.0, std::abs(left[i]), std::abs(right[i])});
        if (!std::isfinite(left[i]) || !std::isfinite(right[i]) ||
            std::abs(left[i] - right[i]) >
                2.0 * pt_service_composition_roundoff * scale) {
            return false;
        }
    }
    return true;
}

} // namespace detail

// The registry is immutable after construction and does not own its backends.
// Each backend and everything it references must outlive the service. Concurrent
// solve calls are permitted only when the selected configured backend is itself
// safe for that use; this boundary adds no hidden serialization or shared cache.
class PtService final {
public:
    explicit PtService(
        std::span<const PtServiceBackendRegistration> backends,
        PtServiceLimits limits = {})
        : limits_(limits) {
        if (!limits_.structurally_valid()) {
            throw PtServiceConfigurationError(
                PtServiceConfigurationErrorCode::invalid_limits,
                "PT service: invalid registry limits");
        }
        if (backends.size() > limits_.max_backends) {
            throw PtServiceConfigurationError(
                PtServiceConfigurationErrorCode::backend_limit,
                "PT service: configured backend limit exceeded");
        }

        entries_.reserve(backends.size());
        descriptors_.reserve(backends.size());
        for (const auto& registration : backends) {
            if (registration.backend == nullptr) {
                throw PtServiceConfigurationError(
                    PtServiceConfigurationErrorCode::null_backend,
                    "PT service: null configured backend");
            }
            if (registration.configured_backend_id.size() >
                limits_.max_identifier_bytes) {
                throw PtServiceConfigurationError(
                    PtServiceConfigurationErrorCode::identifier_limit,
                    "PT service: configured backend ID limit exceeded");
            }
            if (!detail::pt_service_identifier_valid(
                    registration.configured_backend_id)) {
                throw PtServiceConfigurationError(
                    PtServiceConfigurationErrorCode::invalid_registration_id,
                    "PT service: invalid configured backend ID");
            }
            const auto& capability = registration.backend->capability();
            if (!capability.structurally_valid()) {
                throw PtServiceConfigurationError(
                    PtServiceConfigurationErrorCode::invalid_backend_capability,
                    "PT service: invalid configured backend capability");
            }
            if (!detail::pt_service_identifier_valid(capability.backend_id) ||
                std::any_of(
                    capability.component_ids.begin(), capability.component_ids.end(),
                    [](const auto& id) {
                        return !detail::pt_service_identifier_valid(id);
                    })) {
                throw PtServiceConfigurationError(
                    PtServiceConfigurationErrorCode::invalid_backend_capability,
                    "PT service: backend/component IDs must be visible ASCII");
            }
            if (capability.component_ids.size() >
                limits_.max_components_per_backend) {
                throw PtServiceConfigurationError(
                    PtServiceConfigurationErrorCode::component_limit,
                    "PT service: configured component limit exceeded");
            }
            if (capability.backend_id.size() > limits_.max_identifier_bytes ||
                std::any_of(
                    capability.component_ids.begin(), capability.component_ids.end(),
                    [&](const auto& id) {
                        return id.size() > limits_.max_identifier_bytes;
                    })) {
                throw PtServiceConfigurationError(
                    PtServiceConfigurationErrorCode::identifier_limit,
                    "PT service: configured identifier limit exceeded");
            }
            if (find_descriptor(registration.configured_backend_id) != nullptr) {
                throw PtServiceConfigurationError(
                    PtServiceConfigurationErrorCode::duplicate_configured_backend_id,
                    "PT service: duplicate configured backend ID");
            }

            auto descriptor = detail::make_descriptor(
                registration.configured_backend_id, capability);
            if (!descriptor.structurally_valid()) {
                throw PtServiceConfigurationError(
                    PtServiceConfigurationErrorCode::invalid_backend_capability,
                    "PT service: component inventory projection failed");
            }
            entries_.push_back(registration.backend);
            descriptors_.push_back(std::move(descriptor));
        }
    }

    [[nodiscard]] std::span<const PtServiceBackendDescriptor>
    discover_capabilities() const noexcept {
        return descriptors_;
    }

    [[nodiscard]] const PtServiceBackendDescriptor* find_capability(
        std::string_view configured_backend_id) const noexcept {
        return find_descriptor(configured_backend_id);
    }

    [[nodiscard]] PtServiceResponse solve(const PtServiceRequest& request) {
        if (request.configured_backend_id.size() > limits_.max_identifier_bytes ||
            !detail::pt_service_identifier_valid(request.configured_backend_id)) {
            return error_response(
                PtServiceErrorCode::invalid_configured_backend_id,
                "configured_backend_id",
                "configured_backend_id must be a nonempty registered identifier "
                "within the service limit");
        }

        const std::size_t backend_index =
            find_backend_index(request.configured_backend_id);
        if (backend_index == descriptors_.size()) {
            return error_response(
                PtServiceErrorCode::configured_backend_not_found,
                "configured_backend_id",
                "requested configured PT backend is not registered");
        }
        const auto& descriptor = descriptors_[backend_index];

        if (!std::isfinite(request.pressure_pa) || request.pressure_pa <= 0.0) {
            return error_response(PtServiceErrorCode::invalid_pressure,
                                  "pressure_pa",
                                  "pressure must be finite and greater than 0 Pa");
        }
        if (!std::isfinite(request.temperature_k) ||
            request.temperature_k <= 0.0) {
            return error_response(PtServiceErrorCode::invalid_temperature,
                                  "temperature_k",
                                  "temperature must be finite and greater than 0 K");
        }

        const std::size_t component_count =
            descriptor.component_inventory.components.size();
        if (request.feed.size() != component_count) {
            return error_response(
                PtServiceErrorCode::component_count_mismatch, "feed",
                "feed must contain every discovered component exactly once");
        }

        std::vector<double> ordered_feed(component_count, 0.0);
        std::vector<bool> seen(component_count, false);
        double sum = 0.0;
        double correction = 0.0;
        for (const auto& component : request.feed) {
            if (component.component_id.size() > limits_.max_identifier_bytes ||
                !detail::pt_service_identifier_valid(component.component_id)) {
                return error_response(
                    PtServiceErrorCode::invalid_component_id,
                    "feed.component_id",
                    "component_id must be a nonempty discovered identifier "
                    "within the service limit");
            }
            const std::size_t index = component_index(
                descriptor.component_inventory, component.component_id);
            if (index == component_count) {
                return error_response(PtServiceErrorCode::unknown_component,
                                      "feed.component_id",
                                      "feed contains a component not present in "
                                      "the discovered inventory");
            }
            if (seen[index]) {
                return error_response(PtServiceErrorCode::duplicate_component,
                                      "feed.component_id",
                                      "feed contains a duplicate component ID");
            }
            if (!std::isfinite(component.mole_fraction) ||
                component.mole_fraction < 0.0 ||
                component.mole_fraction > 1.0) {
                return error_response(PtServiceErrorCode::invalid_mole_fraction,
                                      "feed.mole_fraction",
                                      "every mole fraction must be finite and in [0,1]");
            }
            seen[index] = true;
            ordered_feed[index] = component.mole_fraction;
            detail::pt_service_add(component.mole_fraction, sum, correction);
        }
        if (std::abs(sum - 1.0) > detail::pt_service_composition_roundoff) {
            return error_response(
                PtServiceErrorCode::composition_not_normalized, "feed",
                "feed mole fractions must sum to 1 within the service roundoff contract");
        }

        flash::PtFlashBackendResult backend_result;
        try {
            backend_result = entries_[backend_index]->solve(
                {request.pressure_pa, request.temperature_k, ordered_feed});
        } catch (const std::bad_alloc&) {
            throw;
        } catch (const std::invalid_argument& error) {
            return backend_error(PtServiceErrorCode::backend_rejected_request,
                                 error.what());
        } catch (const std::domain_error& error) {
            return backend_error(PtServiceErrorCode::backend_rejected_request,
                                 error.what());
        } catch (const std::length_error& error) {
            return backend_error(PtServiceErrorCode::backend_rejected_request,
                                 error.what());
        } catch (const std::logic_error& error) {
            return backend_error(PtServiceErrorCode::backend_contract_violation,
                                 error.what());
        } catch (const std::exception& error) {
            return backend_error(PtServiceErrorCode::backend_execution_failure,
                                 error.what());
        } catch (...) {
            return backend_error(PtServiceErrorCode::backend_execution_failure,
                                 "non-standard backend exception");
        }

        if (!backend_result.structurally_valid() ||
            !detail::same_backend_capability(backend_result.capability,
                                             descriptor.capability) ||
            backend_result.solution.pressure_pa != request.pressure_pa ||
            backend_result.solution.temperature_k != request.temperature_k ||
            !detail::same_composition_within_roundoff(
                backend_result.solution.feed, ordered_feed)) {
            return error_response(
                PtServiceErrorCode::backend_contract_violation, "backend_result",
                "backend result does not match the discovered capability or requested PT state");
        }

        PtServiceResponse response;
        response.result = project_result(descriptor, backend_result);
        switch (backend_result.solution.status) {
        case flash::PtPhaseSetStatus::accepted:
            response.outcome = PtServiceOutcome::accepted;
            break;
        case flash::PtPhaseSetStatus::phase_set_unstable:
            response.outcome = PtServiceOutcome::phase_set_unstable;
            break;
        case flash::PtPhaseSetStatus::indeterminate:
            response.outcome = PtServiceOutcome::indeterminate;
            break;
        }

        if (!response.structurally_valid()) {
            return error_response(
                PtServiceErrorCode::backend_contract_violation, "backend_result",
                "backend result cannot be represented by PT service boundary v1");
        }
        return response;
    }

private:
    [[nodiscard]] const PtServiceBackendDescriptor* find_descriptor(
        std::string_view configured_backend_id) const noexcept {
        for (const auto& descriptor : descriptors_) {
            if (descriptor.configured_backend_id == configured_backend_id) {
                return &descriptor;
            }
        }
        return nullptr;
    }

    [[nodiscard]] std::size_t find_backend_index(
        std::string_view configured_backend_id) const noexcept {
        for (std::size_t i = 0; i < descriptors_.size(); ++i) {
            if (descriptors_[i].configured_backend_id == configured_backend_id) {
                return i;
            }
        }
        return descriptors_.size();
    }

    [[nodiscard]] static std::size_t component_index(
        const PtRuntimeComponentInventory& inventory,
        std::string_view component_id) noexcept {
        for (std::size_t i = 0; i < inventory.components.size(); ++i) {
            if (inventory.components[i].component_id == component_id) { return i; }
        }
        return inventory.components.size();
    }

    [[nodiscard]] static PtServiceResponse error_response(
        PtServiceErrorCode code, std::string field, std::string diagnostic) {
        PtServiceResponse response;
        response.outcome = PtServiceOutcome::error;
        response.error = PtServiceError{
            code, std::move(field), std::move(diagnostic)};
        return response;
    }

    [[nodiscard]] static PtServiceResponse backend_error(
        PtServiceErrorCode code, std::string_view diagnostic) {
        std::string message = "configured PT backend call failed";
        if (!diagnostic.empty()) {
            message += ": ";
            message += diagnostic;
        }
        return error_response(code, "backend", std::move(message));
    }

    [[nodiscard]] static PtServiceComputationResult project_result(
        const PtServiceBackendDescriptor& descriptor,
        const flash::PtFlashBackendResult& source) {
        PtServiceComputationResult result;
        result.provenance.backend = descriptor;
        result.provenance.backend_result_convention =
            std::string(flash::PtFlashBackendResult::convention);
        result.provenance.phase_set_convention =
            std::string(flash::PtPhaseSetResult::convention);
        result.provenance.phase_transition_convention =
            std::string(flash::PtPhaseTransitionReport::convention);
        result.provenance.provider_result_convention =
            source.provider_result_convention;
        result.pressure_pa = source.solution.pressure_pa;
        result.temperature_k = source.solution.temperature_k;
        result.transition_report = source.transition_report;
        result.global_stability_proven =
            source.solution.global_stability_proven;
        result.morphology_resolved = source.morphology_resolved;
        result.diagnostic = source.solution.diagnostic;

        result.feed.reserve(source.solution.feed.size());
        for (std::size_t i = 0; i < source.solution.feed.size(); ++i) {
            result.feed.push_back({
                descriptor.component_inventory.components[i].component_id,
                source.solution.feed[i]});
        }

        const auto* accepted = source.accepted_phase_set();
        if (accepted == nullptr) { return result; }

        result.phases.reserve(accepted->phases.size());
        for (std::size_t phase_index = 0;
             phase_index < accepted->phases.size(); ++phase_index) {
            const auto& source_phase = accepted->phases[phase_index];
            PtServicePhase phase;
            phase.phase_index = phase_index;
            phase.mole_phase_fraction = source_phase.mole_phase_fraction;
            phase.provider_branch = source_phase.activity.branch;
            phase.provider_branch_smooth = source_phase.activity.smooth;
            phase.compressibility_factor =
                source_phase.compressibility_factor;
            if (!source.phase_metadata.empty()) {
                phase.provider_metadata = PtServicePhaseMetadata{
                    source.phase_metadata[phase_index].role_id,
                    source.phase_metadata[phase_index].family_id};
            }
            if (source_phase.composition.size() != result.feed.size() ||
                source_phase.activity.ln_phi.size() != result.feed.size()) {
                return {};
            }
            phase.components.reserve(result.feed.size());
            for (std::size_t component_index = 0;
                 component_index < result.feed.size(); ++component_index) {
                phase.components.push_back({
                    result.feed[component_index].component_id,
                    source_phase.composition[component_index],
                    source_phase.activity.ln_phi[component_index]});
            }
            result.phases.push_back(std::move(phase));
        }
        return result;
    }

    PtServiceLimits limits_;
    std::vector<flash::PtFlashBackend*> entries_;
    std::vector<PtServiceBackendDescriptor> descriptors_;
};

} // namespace mpmc::runtime

#endif // MPMC_RUNTIME_PT_SERVICE_HPP
