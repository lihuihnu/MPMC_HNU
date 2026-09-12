#include <mpmc/runtime_grpc/pt_grpc_adapter.hpp>

#include <mpmc/flash/pt_flash_backend.hpp>
#include <mpmc/flash/pt_phase_set.hpp>
#include <mpmc/flash/pt_phase_transition.hpp>
#include <mpmc/runtime/pt_service_contract.hpp>

#include <grpcpp/resource_quota.h>

#include <chrono>
#include <cstdint>
#include <exception>
#include <limits>
#include <new>
#include <stdexcept>
#include <string>
#include <string_view>
#include <utility>

namespace mpmc::runtime_grpc {
namespace {

namespace core = ::mpmc::runtime;
namespace flash = ::mpmc::flash;
namespace wire = ::mpmc::runtime::v1;

inline constexpr std::string_view pt_wire_contract =
    "mpmc.runtime.v1/PT-flash-service/v1";

[[nodiscard]] grpc::Status invalid_wire(std::string message) {
    return {grpc::StatusCode::INVALID_ARGUMENT, std::move(message)};
}

[[nodiscard]] grpc::Status resource_exhausted(std::string message) {
    return {grpc::StatusCode::RESOURCE_EXHAUSTED, std::move(message)};
}

[[nodiscard]] grpc::Status internal_error() {
    return {grpc::StatusCode::INTERNAL, "PT process adapter internal failure"};
}

[[nodiscard]] std::uint32_t checked_uint32(std::size_t value) {
    if (value > std::numeric_limits<std::uint32_t>::max()) {
        throw std::overflow_error("PT service index does not fit uint32");
    }
    return static_cast<std::uint32_t>(value);
}

[[nodiscard]] wire::PtPhaseTransitionSupport transition_support(
    flash::PtPhaseTransitionSupport value) {
    switch (value) {
    case flash::PtPhaseTransitionSupport::detection_only:
        return wire::PT_PHASE_TRANSITION_SUPPORT_DETECTION_ONLY;
    case flash::PtPhaseTransitionSupport::fresh_target_resolve:
        return wire::PT_PHASE_TRANSITION_SUPPORT_FRESH_TARGET_RESOLVE;
    }
    throw std::logic_error("unknown phase transition support");
}

[[nodiscard]] wire::PtPhaseTransitionTrigger transition_trigger(
    flash::PtPhaseTransitionTrigger value) {
    switch (value) {
    case flash::PtPhaseTransitionTrigger::initial_stability_witness:
        return wire::PT_PHASE_TRANSITION_TRIGGER_INITIAL_STABILITY_WITNESS;
    case flash::PtPhaseTransitionTrigger::final_phase_set_instability:
        return wire::PT_PHASE_TRANSITION_TRIGGER_FINAL_PHASE_SET_INSTABILITY;
    case flash::PtPhaseTransitionTrigger::phase_disappearance:
        return wire::PT_PHASE_TRANSITION_TRIGGER_PHASE_DISAPPEARANCE;
    case flash::PtPhaseTransitionTrigger::provider_topology_witness:
        return wire::PT_PHASE_TRANSITION_TRIGGER_PROVIDER_TOPOLOGY_WITNESS;
    case flash::PtPhaseTransitionTrigger::provider_boundary_route:
        return wire::PT_PHASE_TRANSITION_TRIGGER_PROVIDER_BOUNDARY_ROUTE;
    }
    throw std::logic_error("unknown phase transition trigger");
}

[[nodiscard]] wire::PtPhaseTransitionResolution transition_resolution(
    flash::PtPhaseTransitionResolution value) {
    switch (value) {
    case flash::PtPhaseTransitionResolution::accepted_target:
        return wire::PT_PHASE_TRANSITION_RESOLUTION_ACCEPTED_TARGET;
    case flash::PtPhaseTransitionResolution::target_resolve_required:
        return wire::PT_PHASE_TRANSITION_RESOLUTION_TARGET_RESOLVE_REQUIRED;
    case flash::PtPhaseTransitionResolution::target_resolve_failed:
        return wire::PT_PHASE_TRANSITION_RESOLUTION_TARGET_RESOLVE_FAILED;
    case flash::PtPhaseTransitionResolution::candidate_not_accepted:
        return wire::PT_PHASE_TRANSITION_RESOLUTION_CANDIDATE_NOT_ACCEPTED;
    case flash::PtPhaseTransitionResolution::broader_topology_required:
        return wire::PT_PHASE_TRANSITION_RESOLUTION_BROADER_TOPOLOGY_REQUIRED;
    case flash::PtPhaseTransitionResolution::indeterminate:
        return wire::PT_PHASE_TRANSITION_RESOLUTION_INDETERMINATE;
    }
    throw std::logic_error("unknown phase transition resolution");
}

[[nodiscard]] wire::PtServiceErrorCode service_error_code(
    core::PtServiceErrorCode value) {
    switch (value) {
    case core::PtServiceErrorCode::invalid_configured_backend_id:
        return wire::PT_SERVICE_ERROR_CODE_INVALID_CONFIGURED_BACKEND_ID;
    case core::PtServiceErrorCode::configured_backend_not_found:
        return wire::PT_SERVICE_ERROR_CODE_CONFIGURED_BACKEND_NOT_FOUND;
    case core::PtServiceErrorCode::invalid_pressure:
        return wire::PT_SERVICE_ERROR_CODE_INVALID_PRESSURE;
    case core::PtServiceErrorCode::invalid_temperature:
        return wire::PT_SERVICE_ERROR_CODE_INVALID_TEMPERATURE;
    case core::PtServiceErrorCode::component_count_mismatch:
        return wire::PT_SERVICE_ERROR_CODE_COMPONENT_COUNT_MISMATCH;
    case core::PtServiceErrorCode::invalid_component_id:
        return wire::PT_SERVICE_ERROR_CODE_INVALID_COMPONENT_ID;
    case core::PtServiceErrorCode::duplicate_component:
        return wire::PT_SERVICE_ERROR_CODE_DUPLICATE_COMPONENT;
    case core::PtServiceErrorCode::unknown_component:
        return wire::PT_SERVICE_ERROR_CODE_UNKNOWN_COMPONENT;
    case core::PtServiceErrorCode::invalid_mole_fraction:
        return wire::PT_SERVICE_ERROR_CODE_INVALID_MOLE_FRACTION;
    case core::PtServiceErrorCode::composition_not_normalized:
        return wire::PT_SERVICE_ERROR_CODE_COMPOSITION_NOT_NORMALIZED;
    case core::PtServiceErrorCode::backend_rejected_request:
        return wire::PT_SERVICE_ERROR_CODE_BACKEND_REJECTED_REQUEST;
    case core::PtServiceErrorCode::backend_contract_violation:
        return wire::PT_SERVICE_ERROR_CODE_BACKEND_CONTRACT_VIOLATION;
    case core::PtServiceErrorCode::backend_execution_failure:
        return wire::PT_SERVICE_ERROR_CODE_BACKEND_EXECUTION_FAILURE;
    }
    throw std::logic_error("unknown PT service error code");
}

void map_transition_capability(
    const flash::PtPhaseTransitionCapability& source,
    wire::PtPhaseTransitionCapability& target) {
    target.set_convention(flash::PtPhaseTransitionCapability::convention.data(),
                          flash::PtPhaseTransitionCapability::convention.size());
    for (const auto& edge : source.edges) {
        auto* output = target.add_edges();
        output->set_source_phase_count(checked_uint32(edge.source_phase_count));
        output->set_target_phase_count(checked_uint32(edge.target_phase_count));
        output->set_support(transition_support(edge.support));
        output->set_requires_fresh_target_solve(
            edge.requires_fresh_target_solve);
    }
}

void map_capability(const flash::PtFlashBackendCapability& source,
                    wire::PtBackendCapability& target) {
    target.set_convention(flash::PtFlashBackendCapability::convention.data(),
                          flash::PtFlashBackendCapability::convention.size());
    target.set_backend_id(source.backend_id);
    target.set_model_profile(source.model_profile);
    target.set_algorithm_profile(source.algorithm_profile);
    target.set_publication_profile(source.publication_profile);
    target.set_configuration_profile(source.configuration_profile);
    target.set_dataset_id(source.dataset_id);
    target.set_revision(source.revision);
    for (const auto& component_id : source.component_ids) {
        target.add_component_ids(component_id);
    }
    for (const std::size_t count : source.supported_phase_counts) {
        target.add_supported_phase_counts(checked_uint32(count));
    }
    for (const auto& setting : source.scalar_settings) {
        auto* output = target.add_scalar_settings();
        output->set_id(setting.id);
        output->set_value(setting.value);
        output->set_unit(setting.unit);
    }
    map_transition_capability(source.transition_capability,
                              *target.mutable_transition_capability());
    target.set_performs_initial_stability_search(
        source.performs_initial_stability_search);
    target.set_performs_final_phase_set_review(
        source.performs_final_phase_set_review);
    target.set_performs_boundary_neighbor_resolve(
        source.performs_boundary_neighbor_resolve);
    target.set_global_stability_proven(source.global_stability_proven);
    target.set_phase_metadata_namespace(source.phase_metadata_namespace);
}

void map_descriptor(const core::PtServiceBackendDescriptor& source,
                    wire::PtBackendDescriptor& target) {
    target.set_convention(core::PtServiceBackendDescriptor::convention.data(),
                          core::PtServiceBackendDescriptor::convention.size());
    target.set_configured_backend_id(source.configured_backend_id);
    map_capability(source.capability, *target.mutable_capability());
    auto* inventory = target.mutable_component_inventory();
    inventory->set_convention(core::PtRuntimeComponentInventory::convention.data(),
                              core::PtRuntimeComponentInventory::convention.size());
    for (const auto& component : source.component_inventory.components) {
        auto* output = inventory->add_components();
        output->set_component_id(component.component_id);
        output->set_feed_index(checked_uint32(component.feed_index));
    }
}

void map_transition_report(const flash::PtPhaseTransitionReport& source,
                           wire::PtPhaseTransitionReport& target) {
    target.set_convention(flash::PtPhaseTransitionReport::convention.data(),
                          flash::PtPhaseTransitionReport::convention.size());
    for (const auto& evidence : source.evidence) {
        auto* output = target.add_evidence();
        output->set_source_phase_count(
            checked_uint32(evidence.source_phase_count));
        if (evidence.target_phase_count) {
            output->set_target_phase_count(
                checked_uint32(*evidence.target_phase_count));
        }
        output->set_trigger(transition_trigger(evidence.trigger));
        output->set_resolution(transition_resolution(evidence.resolution));
        output->set_fresh_target_solve_attempted(
            evidence.fresh_target_solve_attempted);
        output->set_target_topology_closed(evidence.target_topology_closed);
        output->set_provider_evidence_profile(
            evidence.provider_evidence_profile);
        output->set_diagnostic(evidence.diagnostic);
    }
}

void map_result(const core::PtServiceComputationResult& source,
                core::PtServiceOutcome outcome,
                wire::PtComputationResult& target) {
    switch (outcome) {
    case core::PtServiceOutcome::accepted:
        target.set_outcome(wire::PT_COMPUTATION_OUTCOME_ACCEPTED);
        break;
    case core::PtServiceOutcome::phase_set_unstable:
        target.set_outcome(wire::PT_COMPUTATION_OUTCOME_PHASE_SET_UNSTABLE);
        break;
    case core::PtServiceOutcome::indeterminate:
        target.set_outcome(wire::PT_COMPUTATION_OUTCOME_INDETERMINATE);
        break;
    case core::PtServiceOutcome::error:
        throw std::logic_error("service error cannot map to result arm");
    }

    auto* provenance = target.mutable_provenance();
    map_descriptor(source.provenance.backend, *provenance->mutable_backend());
    provenance->set_backend_result_convention(
        source.provenance.backend_result_convention);
    provenance->set_phase_set_convention(source.provenance.phase_set_convention);
    provenance->set_phase_transition_convention(
        source.provenance.phase_transition_convention);
    provenance->set_provider_result_convention(
        source.provenance.provider_result_convention);

    target.set_pressure_pa(source.pressure_pa);
    target.set_temperature_k(source.temperature_k);
    for (const auto& component : source.feed) {
        auto* output = target.add_feed();
        output->set_component_id(component.component_id);
        output->set_mole_fraction(component.mole_fraction);
    }
    for (const auto& phase : source.phases) {
        auto* output = target.add_phases();
        output->set_phase_index(checked_uint32(phase.phase_index));
        output->set_mole_phase_fraction(phase.mole_phase_fraction);
        for (const auto& component : phase.components) {
            auto* wire_component = output->add_components();
            wire_component->set_component_id(component.component_id);
            wire_component->set_mole_fraction(component.mole_fraction);
            wire_component->set_ln_fugacity_coefficient(
                component.ln_fugacity_coefficient);
        }
        static_assert(sizeof(std::size_t) <= sizeof(std::uint64_t));
        output->set_provider_branch(
            static_cast<std::uint64_t>(phase.provider_branch));
        output->set_provider_branch_smooth(phase.provider_branch_smooth);
        if (phase.compressibility_factor) {
            output->set_compressibility_factor(
                *phase.compressibility_factor);
        }
        if (phase.provider_metadata) {
            auto* metadata = output->mutable_provider_metadata();
            metadata->set_role_id(phase.provider_metadata->role_id);
            metadata->set_family_id(phase.provider_metadata->family_id);
        }
    }
    map_transition_report(source.transition_report,
                          *target.mutable_transition_report());
    target.set_global_stability_proven(source.global_stability_proven);
    target.set_morphology_resolved(source.morphology_resolved);
    target.set_diagnostic(source.diagnostic);
}

void map_service_response(const core::PtServiceResponse& source,
                          wire::SolvePtFlashResponse& target) {
    if (!source.structurally_valid()) {
        throw std::logic_error("invalid PtService response");
    }
    target.set_wire_contract(pt_wire_contract.data(), pt_wire_contract.size());
    target.set_service_result_convention(core::PtServiceResponse::convention.data(),
                                         core::PtServiceResponse::convention.size());
    if (source.outcome == core::PtServiceOutcome::error) {
        const auto& source_error = *source.error;
        auto* output = target.mutable_error();
        output->set_code(service_error_code(source_error.code));
        output->set_field(source_error.field);
        output->set_diagnostic(source_error.diagnostic);
        return;
    }
    map_result(*source.result, source.outcome, *target.mutable_result());
}

[[nodiscard]] grpc::Status parse_request(
    const wire::SolvePtFlashRequest& source, std::size_t byte_limit,
    core::PtServiceRequest& target) {
    if (source.ByteSizeLong() > byte_limit) {
        return resource_exhausted(
            "serialized PT request exceeds the process adapter limit");
    }
    if (!source.has_wire_contract() || source.wire_contract() != pt_wire_contract ||
        !source.has_service_request_convention() ||
        source.service_request_convention() != core::PtServiceRequest::convention) {
        return invalid_wire("PT request wire/convention version mismatch");
    }
    if (!source.has_configured_backend_id() || !source.has_pressure_pa() ||
        !source.has_temperature_k()) {
        return invalid_wire("PT request is missing a required scalar field");
    }

    target.configured_backend_id = source.configured_backend_id();
    target.pressure_pa = source.pressure_pa();
    target.temperature_k = source.temperature_k();
    target.feed.reserve(static_cast<std::size_t>(source.feed_size()));
    for (const auto& component : source.feed()) {
        if (!component.has_component_id() || !component.has_mole_fraction()) {
            return invalid_wire(
                "PT request feed entry is missing component_id or mole_fraction");
        }
        target.feed.push_back(
            {component.component_id(), component.mole_fraction()});
    }
    return grpc::Status::OK;
}

[[nodiscard]] grpc::Status lifecycle_status(
    const grpc::ServerContext& context,
    std::chrono::milliseconds maximum_deadline) {
    const auto now = std::chrono::system_clock::now();
    const auto deadline = context.deadline();
    if (deadline == std::chrono::system_clock::time_point::max()) {
        return invalid_wire("a finite client RPC deadline is required");
    }
    if (deadline <= now) {
        return {grpc::StatusCode::DEADLINE_EXCEEDED,
                "PT RPC deadline elapsed before publication"};
    }
    if (deadline - now > maximum_deadline) {
        return invalid_wire("client RPC deadline exceeds the process policy");
    }
    if (context.IsCancelled()) {
        return {grpc::StatusCode::CANCELLED,
                "PT RPC was cancelled before publication"};
    }
    return grpc::Status::OK;
}

void clear_on_error(const grpc::Status& status,
                    wire::SolvePtFlashResponse& response) {
    if (!status.ok()) { response.Clear(); }
}

} // namespace

bool PtGrpcAdapterLimits::structurally_valid() const noexcept {
    constexpr auto max_int =
        static_cast<std::size_t>(std::numeric_limits<int>::max());
    return max_serialized_request_bytes > 0U &&
           max_serialized_request_bytes <= max_int &&
           max_serialized_response_bytes > 0U &&
           max_serialized_response_bytes <= max_int &&
           max_concurrent_solves > 0U && grpc_resource_quota_bytes > 0U &&
           grpc_resource_quota_bytes >= max_serialized_request_bytes &&
           grpc_resource_quota_bytes >= max_serialized_response_bytes &&
           max_discovery_deadline.count() > 0 &&
           max_solve_deadline.count() > 0;
}

PtGrpcServiceAdapter::PtGrpcServiceAdapter(core::PtService& service,
                                           PtGrpcAdapterLimits limits,
                                           std::shared_ptr<PtGrpcObserver> observer)
    : service_(service), limits_(limits), observer_(std::move(observer)) {
    if (!limits_.structurally_valid()) {
        throw std::invalid_argument("invalid PT gRPC process adapter limits");
    }
}

void PtGrpcServiceAdapter::observe(PtGrpcObservation observation) const noexcept {
    if (observer_) { observer_->observe(observation); }
}

grpc::Status PtGrpcServiceAdapter::DiscoverPtCapabilities(
    grpc::ServerContext* context,
    const wire::DiscoverPtCapabilitiesRequest* request,
    wire::DiscoverPtCapabilitiesResponse* response) {
    const auto started_at = std::chrono::steady_clock::now();
    const auto finish = [&](grpc::Status result, PtGrpcCompletion completion) {
        std::size_t request_bytes = 0U;
        std::size_t response_bytes = 0U;
        try {
            if (request != nullptr) { request_bytes = request->ByteSizeLong(); }
            if (response != nullptr) { response_bytes = response->ByteSizeLong(); }
        } catch (...) {
            // Observation is best-effort and must not replace the RPC result.
        }
        observe({PtGrpcRpcMethod::discover_pt_capabilities,
                 completion,
                 result.error_code(),
                 std::nullopt,
                 false,
                 request_bytes,
                 response_bytes,
                 std::chrono::duration_cast<std::chrono::microseconds>(
                     std::chrono::steady_clock::now() - started_at)});
        return result;
    };
    if (context == nullptr || request == nullptr || response == nullptr) {
        return finish(internal_error(), PtGrpcCompletion::internal_failure);
    }
    response->Clear();
    auto status = lifecycle_status(*context, limits_.max_discovery_deadline);
    if (!status.ok()) {
        return finish(status, status.error_code() == grpc::StatusCode::CANCELLED
                                  ? PtGrpcCompletion::cancelled
                                  : status.error_code() ==
                                            grpc::StatusCode::DEADLINE_EXCEEDED
                                        ? PtGrpcCompletion::deadline_exceeded
                                        : PtGrpcCompletion::invalid_request);
    }
    if (request->ByteSizeLong() > limits_.max_serialized_request_bytes) {
        return finish(
            resource_exhausted(
                "serialized discovery request exceeds the process adapter limit"),
            PtGrpcCompletion::request_limit);
    }

    PtGrpcCompletion completion = PtGrpcCompletion::completed;
    try {
        response->set_wire_contract(pt_wire_contract.data(),
                                    pt_wire_contract.size());
        response->set_service_boundary_convention(
            core::pt_service_boundary_convention.data(),
            core::pt_service_boundary_convention.size());
        response->set_capability_convention(
            core::pt_service_capability_convention.data(),
            core::pt_service_capability_convention.size());
        for (const auto& descriptor : service_.discover_capabilities()) {
            map_descriptor(descriptor, *response->add_backends());
        }
        status = lifecycle_status(*context, limits_.max_discovery_deadline);
        if (status.ok() &&
            response->ByteSizeLong() > limits_.max_serialized_response_bytes) {
            status = resource_exhausted(
                "serialized discovery response exceeds the process adapter limit");
        }
    } catch (const std::bad_alloc&) {
        status = resource_exhausted(
            "PT capability discovery exhausted process memory");
        completion = PtGrpcCompletion::memory_exhausted;
    } catch (const std::exception&) {
        status = internal_error();
        completion = PtGrpcCompletion::internal_failure;
    } catch (...) {
        status = internal_error();
        completion = PtGrpcCompletion::internal_failure;
    }
    if (!status.ok()) {
        response->Clear();
        if (completion == PtGrpcCompletion::completed &&
            status.error_code() == grpc::StatusCode::RESOURCE_EXHAUSTED) {
            completion = PtGrpcCompletion::response_limit;
        } else if (completion == PtGrpcCompletion::completed &&
                   status.error_code() == grpc::StatusCode::CANCELLED) {
            completion = PtGrpcCompletion::cancelled;
        } else if (completion == PtGrpcCompletion::completed &&
                   status.error_code() == grpc::StatusCode::DEADLINE_EXCEEDED) {
            completion = PtGrpcCompletion::deadline_exceeded;
        } else if (completion == PtGrpcCompletion::completed) {
            completion = PtGrpcCompletion::internal_failure;
        }
    }
    return finish(status, completion);
}

bool PtGrpcServiceAdapter::try_acquire_solve() noexcept {
    std::size_t current = in_flight_solves_.load(std::memory_order_relaxed);
    while (current < limits_.max_concurrent_solves) {
        if (in_flight_solves_.compare_exchange_weak(
                current, current + 1U, std::memory_order_acquire,
                std::memory_order_relaxed)) {
            return true;
        }
    }
    return false;
}

void PtGrpcServiceAdapter::release_solve() noexcept {
    in_flight_solves_.fetch_sub(1U, std::memory_order_release);
}

grpc::Status PtGrpcServiceAdapter::SolvePtFlash(
    grpc::ServerContext* context, const wire::SolvePtFlashRequest* request,
    wire::SolvePtFlashResponse* response) {
    const auto started_at = std::chrono::steady_clock::now();
    bool pt_service_called = false;
    std::optional<core::PtServiceOutcome> service_outcome;
    const auto finish = [&](grpc::Status result, PtGrpcCompletion completion) {
        std::size_t request_bytes = 0U;
        std::size_t response_bytes = 0U;
        try {
            if (request != nullptr) { request_bytes = request->ByteSizeLong(); }
            if (response != nullptr) { response_bytes = response->ByteSizeLong(); }
        } catch (...) {
            // Observation is best-effort and must not replace the RPC result.
        }
        observe({PtGrpcRpcMethod::solve_pt_flash,
                 completion,
                 result.error_code(),
                 service_outcome,
                 pt_service_called,
                 request_bytes,
                 response_bytes,
                 std::chrono::duration_cast<std::chrono::microseconds>(
                     std::chrono::steady_clock::now() - started_at)});
        return result;
    };
    if (context == nullptr || request == nullptr || response == nullptr) {
        return finish(internal_error(), PtGrpcCompletion::internal_failure);
    }
    response->Clear();
    auto status = lifecycle_status(*context, limits_.max_solve_deadline);
    if (!status.ok()) {
        return finish(status, status.error_code() == grpc::StatusCode::CANCELLED
                                  ? PtGrpcCompletion::cancelled
                                  : status.error_code() ==
                                            grpc::StatusCode::DEADLINE_EXCEEDED
                                        ? PtGrpcCompletion::deadline_exceeded
                                        : PtGrpcCompletion::invalid_request);
    }

    core::PtServiceRequest service_request;
    PtGrpcCompletion completion = PtGrpcCompletion::completed;
    try {
        status = parse_request(*request, limits_.max_serialized_request_bytes,
                               service_request);
    } catch (const std::bad_alloc&) {
        return finish(
            resource_exhausted("PT request mapping exhausted process memory"),
            PtGrpcCompletion::memory_exhausted);
    } catch (const std::exception&) {
        return finish(
            invalid_wire("PT request could not be mapped to service v1"),
            PtGrpcCompletion::invalid_request);
    }
    if (!status.ok()) {
        return finish(status,
                      status.error_code() == grpc::StatusCode::RESOURCE_EXHAUSTED
                          ? PtGrpcCompletion::request_limit
                          : PtGrpcCompletion::invalid_request);
    }
    if (!try_acquire_solve()) {
        return finish(
            resource_exhausted(
                "maximum concurrent PT solves are already in flight"),
            PtGrpcCompletion::concurrency_limit);
    }

    status = lifecycle_status(*context, limits_.max_solve_deadline);
    if (!status.ok()) {
        release_solve();
        return finish(status,
                      status.error_code() == grpc::StatusCode::CANCELLED
                          ? PtGrpcCompletion::cancelled
                          : status.error_code() ==
                                    grpc::StatusCode::DEADLINE_EXCEEDED
                                ? PtGrpcCompletion::deadline_exceeded
                                : PtGrpcCompletion::invalid_request);
    }

    core::PtServiceResponse service_response;
    try {
        // The only scientific delegation in this adapter.
        pt_service_called = true;
        service_response = service_.solve(service_request);
        service_outcome = service_response.outcome;
    } catch (const std::bad_alloc&) {
        release_solve();
        return finish(
            resource_exhausted("PT service solve exhausted process memory"),
            PtGrpcCompletion::memory_exhausted);
    } catch (const std::exception&) {
        release_solve();
        return finish(
            {grpc::StatusCode::INTERNAL,
             "PT service boundary raised an unexpected exception"},
            PtGrpcCompletion::internal_failure);
    } catch (...) {
        release_solve();
        return finish(
            {grpc::StatusCode::INTERNAL,
             "PT service boundary raised a non-standard exception"},
            PtGrpcCompletion::internal_failure);
    }
    release_solve();

    status = lifecycle_status(*context, limits_.max_solve_deadline);
    if (!status.ok()) {
        return finish(status,
                      status.error_code() == grpc::StatusCode::CANCELLED
                          ? PtGrpcCompletion::cancelled
                          : PtGrpcCompletion::deadline_exceeded);
    }
    try {
        map_service_response(service_response, *response);
        if (response->ByteSizeLong() > limits_.max_serialized_response_bytes) {
            status = resource_exhausted(
                "serialized PT response exceeds the process adapter limit");
        }
    } catch (const std::bad_alloc&) {
        status = resource_exhausted(
            "PT response mapping exhausted process memory");
        completion = PtGrpcCompletion::memory_exhausted;
    } catch (const std::exception&) {
        status = internal_error();
        completion = PtGrpcCompletion::internal_failure;
    } catch (...) {
        status = internal_error();
        completion = PtGrpcCompletion::internal_failure;
    }
    clear_on_error(status, *response);
    if (!status.ok() && completion == PtGrpcCompletion::completed) {
        completion = status.error_code() == grpc::StatusCode::RESOURCE_EXHAUSTED
                         ? PtGrpcCompletion::response_limit
                         : PtGrpcCompletion::internal_failure;
    }
    return finish(status, completion);
}

void configure_pt_grpc_server(grpc::ServerBuilder& builder,
                              PtGrpcServiceAdapter& adapter) {
    const auto& limits = adapter.limits();
    if (!limits.structurally_valid()) {
        throw std::invalid_argument("invalid PT gRPC process adapter limits");
    }
    builder.SetMaxReceiveMessageSize(
        static_cast<int>(limits.max_serialized_request_bytes));
    builder.SetMaxSendMessageSize(
        static_cast<int>(limits.max_serialized_response_bytes));
    grpc::ResourceQuota quota("mpmc-pt-grpc-process");
    quota.Resize(limits.grpc_resource_quota_bytes);
    builder.SetResourceQuota(quota);
    builder.RegisterService(&adapter);
}

} // namespace mpmc::runtime_grpc
