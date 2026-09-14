#include <mpmc/model_configuration_grpc/model_grpc_adapter.hpp>

#include <google/protobuf/descriptor.h>
#include <google/protobuf/message.h>
#include <google/protobuf/unknown_field_set.h>

#include <limits>
#include <new>
#include <utility>
#include <vector>

// Windows SDK's message-loop macro conflicts with Protobuf reflection.
// This translation unit uses the C++ method, never the Win32 message loop.
#if defined(GetMessage)
#undef GetMessage
#endif

namespace mpmc::model_configuration_grpc {
namespace {
namespace mc = ::mpmc::model_configuration;
namespace wire = ::mpmc::model_configuration::v1;
using StatusCode = grpc::StatusCode;

struct RpcFailure {
    StatusCode status;
    const char* code;
    std::string field;
};
[[noreturn]] void fail(StatusCode status, const char* code, std::string field = {}) {
    throw RpcFailure{status, code, std::move(field)};
}
grpc::Status error_status(StatusCode status, const char* code, std::string_view field = {}) {
    wire::ModelServiceError details;
    details.set_wire_contract(std::string(model_wire_contract));
    details.set_code(code);
    // Domain field paths may contain bounded component IDs; never echo what(),
    // request bodies, handles or credentials. Keep metadata small even if a host
    // permits long identifiers. Allocation failure here uses an empty fallback.
    details.set_field(field.size() <= 256 && field.find("mh1_") == std::string_view::npos
                          ? std::string(field) : "configuration");
    return model_error_status(status, "model configuration request failed", details);
}
class Admission {
public:
    Admission(std::atomic<std::size_t>& count, std::size_t limit) : count_(count) {
        auto current = count_.load(std::memory_order_relaxed);
        while (current < limit) {
            if (count_.compare_exchange_weak(current, current + 1, std::memory_order_acquire,
                                             std::memory_order_relaxed)) { return; }
        }
        fail(StatusCode::RESOURCE_EXHAUSTED, "rpc.concurrency_limit");
    }
    ~Admission() { count_.fetch_sub(1, std::memory_order_release); }
    Admission(const Admission&) = delete;
    Admission& operator=(const Admission&) = delete;
private:
    std::atomic<std::size_t>& count_;
};
void lifecycle(const grpc::ServerContext& context, std::chrono::milliseconds maximum) {
    const auto deadline = context.deadline();
    const auto now = std::chrono::system_clock::now();
    if (deadline == std::chrono::system_clock::time_point::max()) {
        fail(StatusCode::INVALID_ARGUMENT, "rpc.deadline_required");
    }
    if (deadline <= now) { fail(StatusCode::DEADLINE_EXCEEDED, "rpc.deadline_exceeded"); }
    if (deadline - now > maximum) { fail(StatusCode::INVALID_ARGUMENT, "rpc.deadline_limit"); }
    if (context.IsCancelled()) { fail(StatusCode::CANCELLED, "rpc.cancelled"); }
}
void reject_unknown(const google::protobuf::Message& message, int depth = 0) {
    if (depth > 16) { fail(StatusCode::INVALID_ARGUMENT, "wire.nesting_limit"); }
    const auto* reflection = message.GetReflection();
    if (reflection->GetUnknownFields(message).field_count() != 0) {
        fail(StatusCode::INVALID_ARGUMENT, "wire.unknown_field");
    }
    std::vector<const google::protobuf::FieldDescriptor*> fields;
    reflection->ListFields(message, &fields);
    for (const auto* field : fields) {
        if (field->cpp_type() != google::protobuf::FieldDescriptor::CPPTYPE_MESSAGE) { continue; }
        if (field->is_repeated()) {
            for (int i = 0; i < reflection->FieldSize(message, field); ++i) {
                reject_unknown(reflection->GetRepeatedMessage(message, field, i), depth + 1);
            }
        } else { reject_unknown(reflection->GetMessage(message, field), depth + 1); }
    }
}
const char* configuration_code(mc::ModelConfigurationErrorCode code) {
    switch (code) {
    case mc::ModelConfigurationErrorCode::unsupported_version: return "configuration.unsupported_version";
    case mc::ModelConfigurationErrorCode::unsupported_family: return "configuration.unsupported_family";
    case mc::ModelConfigurationErrorCode::missing_field: return "configuration.missing_field";
    case mc::ModelConfigurationErrorCode::invalid_identifier: return "configuration.invalid_identifier";
    case mc::ModelConfigurationErrorCode::duplicate_identifier: return "configuration.duplicate_identifier";
    case mc::ModelConfigurationErrorCode::unknown_component: return "configuration.unknown_component";
    case mc::ModelConfigurationErrorCode::invalid_value: return "configuration.invalid_value";
    case mc::ModelConfigurationErrorCode::invalid_unit: return "configuration.invalid_unit";
    case mc::ModelConfigurationErrorCode::invalid_source: return "configuration.invalid_source";
    case mc::ModelConfigurationErrorCode::invalid_range: return "configuration.invalid_range";
    case mc::ModelConfigurationErrorCode::duplicate_parameter: return "configuration.duplicate_parameter";
    case mc::ModelConfigurationErrorCode::missing_parameter: return "configuration.missing_parameter";
    case mc::ModelConfigurationErrorCode::invalid_pair: return "configuration.invalid_pair";
    case mc::ModelConfigurationErrorCode::resource_limit: return "configuration.resource_limit";
    case mc::ModelConfigurationErrorCode::unsupported_preset: return "configuration.unsupported_preset";
    case mc::ModelConfigurationErrorCode::invalid_settings: return "configuration.invalid_settings";
    }
    return "configuration.unknown_error";
}
grpc::Status translate_exception() {
    try { throw; }
    catch (const RpcFailure& e) { return error_status(e.status, e.code, e.field); }
    catch (const mc::ModelConfigurationError& e) {
        auto status = StatusCode::INVALID_ARGUMENT;
        switch (e.code()) {
        case mc::ModelConfigurationErrorCode::unsupported_version:
        case mc::ModelConfigurationErrorCode::unsupported_family:
        case mc::ModelConfigurationErrorCode::unsupported_preset: status = StatusCode::UNIMPLEMENTED; break;
        case mc::ModelConfigurationErrorCode::resource_limit: status = StatusCode::RESOURCE_EXHAUSTED; break;
        default: break;
        }
        return error_status(status, configuration_code(e.code()), e.field());
    }
    catch (const mc::ModelRegistryError& e) {
        switch (e.code()) {
        case mc::ModelRegistryErrorCode::capacity_exceeded:
        case mc::ModelRegistryErrorCode::handle_sequence_exhausted:
            return error_status(StatusCode::RESOURCE_EXHAUSTED, "registry.capacity_exceeded");
        case mc::ModelRegistryErrorCode::invalid_handle:
            return error_status(StatusCode::INVALID_ARGUMENT, "registry.invalid_handle", "model_handle");
        case mc::ModelRegistryErrorCode::model_not_found:
            return error_status(StatusCode::NOT_FOUND, "registry.model_not_found", "model_handle");
        case mc::ModelRegistryErrorCode::registry_closed:
            return error_status(StatusCode::FAILED_PRECONDITION, "registry.closed");
        case mc::ModelRegistryErrorCode::entropy_unavailable:
            return error_status(StatusCode::UNAVAILABLE, "registry.entropy_unavailable");
        default: return error_status(StatusCode::INTERNAL, "registry.internal_failure");
        }
    }
    catch (const mc::Pr76SolveRequestError& e) { return error_status(StatusCode::INVALID_ARGUMENT, "request.rejected", e.field()); }
    catch (const mc::Pr76ModelBusyError&) { return error_status(StatusCode::RESOURCE_EXHAUSTED, "model.busy"); }
    catch (const std::bad_alloc&) { return error_status(StatusCode::RESOURCE_EXHAUSTED, "rpc.memory_exhausted"); }
    catch (const std::invalid_argument&) { return error_status(StatusCode::INVALID_ARGUMENT, "request.rejected", "request"); }
    catch (const std::domain_error&) { return error_status(StatusCode::INVALID_ARGUMENT, "request.rejected", "request"); }
    catch (const std::length_error&) { return error_status(StatusCode::INVALID_ARGUMENT, "request.rejected", "request"); }
    catch (...) { return error_status(StatusCode::INTERNAL, "rpc.internal_failure"); }
}
template <class Request, class Response, class Work>
grpc::Status dispatch(grpc::ServerContext* context, const Request* request, Response* response,
    const ModelGrpcAuthorization& authorize, const ModelGrpcLimits& limits,
    std::atomic<std::size_t>& requests, Work&& work) {
    if (context == nullptr || request == nullptr || response == nullptr) {
        return {StatusCode::INTERNAL, "model configuration adapter null argument"};
    }
    response->Clear();
    try {
        const Admission admitted(requests, limits.max_concurrent_requests);
        if (!authorize(*context)) { fail(StatusCode::PERMISSION_DENIED, "rpc.not_authorized"); }
        lifecycle(*context, limits.max_deadline);
        if (request->ByteSizeLong() > limits.max_request_bytes) {
            fail(StatusCode::RESOURCE_EXHAUSTED, "wire.request_limit");
        }
        if (!request->has_wire_contract()) { fail(StatusCode::INVALID_ARGUMENT, "wire.missing_field", "wire_contract"); }
        if (request->wire_contract() != model_wire_contract) {
            fail(StatusCode::INVALID_ARGUMENT, "wire.unsupported_version", "wire_contract");
        }
        reject_unknown(*request);
        work();
        lifecycle(*context, limits.max_deadline);
        response->set_wire_contract(std::string(model_wire_contract));
        if (response->ByteSizeLong() > limits.max_response_bytes) {
            fail(StatusCode::RESOURCE_EXHAUSTED, "wire.response_limit");
        }
        return grpc::Status::OK;
    } catch (...) {
        response->Clear();
        try { return translate_exception(); }
        catch (...) { return {StatusCode::INTERNAL, "model configuration error mapping failed"}; }
    }
}
template <class Request>
void require_handle(const Request& request) {
    if (!request.has_model_handle()) {
        fail(StatusCode::INVALID_ARGUMENT, "wire.missing_field", "model_handle");
    }
}
} // namespace

bool ModelGrpcLimits::structurally_valid() const noexcept {
    constexpr auto max_int = static_cast<std::size_t>(std::numeric_limits<int>::max());
    return max_request_bytes > 0 && max_request_bytes <= max_int &&
           max_response_bytes > 0 && max_response_bytes <= max_int &&
           resource_quota_bytes >= max_request_bytes && resource_quota_bytes >= max_response_bytes &&
           max_concurrent_requests > 0 && max_concurrent_solves > 0 &&
           max_concurrent_solves <= max_concurrent_requests && max_deadline.count() > 0;
}
ModelGrpcServiceAdapter::ModelGrpcServiceAdapter(mc::Pr76ModelRegistry& registry,
    ModelGrpcAuthorization authorize, ModelGrpcLimits limits)
    : registry_(registry), authorize_(std::move(authorize)), limits_(limits) {
    if (!authorize_ || !limits_.structurally_valid()) {
        throw std::invalid_argument("invalid model gRPC limits or missing host authorization policy");
    }
}
grpc::Status ModelGrpcServiceAdapter::CreateModel(grpc::ServerContext* context,
    const wire::CreateModelRequest* request, wire::CreateModelResponse* response) {
    std::string created;
    auto status = dispatch(context, request, response, authorize_, limits_, requests_, [&] {
        if (!request->has_definition()) { fail(StatusCode::INVALID_ARGUMENT, "wire.missing_field", "definition"); }
        auto definition = decode_definition(request->definition(), registry_.limits().parameters);
        mc::PtSolverSettings settings;
        switch (request->solver_selection_case()) {
        case wire::CreateModelRequest::kSettings: settings = decode_settings(request->settings()); break;
        case wire::CreateModelRequest::kPresetId: settings = mc::resolve_pt_solver_preset(request->preset_id()); break;
        default: fail(StatusCode::INVALID_ARGUMENT, "wire.missing_field", "solver_selection");
        }
        lifecycle(*context, limits_.max_deadline);
        created = registry_.create(definition, settings);
        response->set_model_handle(created);
        encode_snapshot(registry_.describe(created), *response->mutable_snapshot());
    });
    // Roll back every locally detected failure after creation, including an
    // oversized response/cancel/deadline. A post-OK network loss is ambiguous;
    // the owning session must eventually close (no automatic create retry).
    if (!status.ok() && !created.empty()) {
        try { registry_.release(created); }
        catch (const mc::ModelRegistryError& e) {
            if (e.code() != mc::ModelRegistryErrorCode::registry_closed &&
                e.code() != mc::ModelRegistryErrorCode::model_not_found) {
                return {StatusCode::INTERNAL, "model creation rollback failed"};
            }
        } catch (...) { return {StatusCode::INTERNAL, "model creation rollback failed"}; }
    }
    return status;
}
grpc::Status ModelGrpcServiceAdapter::DescribeModel(grpc::ServerContext* context,
    const wire::ModelHandleRequest* request, wire::DescribeModelResponse* response) {
    return dispatch(context, request, response, authorize_, limits_, requests_, [&] {
        require_handle(*request);
        encode_snapshot(registry_.describe(request->model_handle()), *response->mutable_snapshot());
    });
}
grpc::Status ModelGrpcServiceAdapter::SolveModel(grpc::ServerContext* context,
    const wire::SolveModelRequest* request, wire::SolveModelResponse* response) {
    return dispatch(context, request, response, authorize_, limits_, requests_, [&] {
        require_handle(*request);
        if (!request->has_pressure_pa()) {
            fail(StatusCode::INVALID_ARGUMENT, "wire.missing_field", "pressure_pa");
        }
        if (!request->has_temperature_k()) {
            fail(StatusCode::INVALID_ARGUMENT, "wire.missing_field", "temperature_k");
        }
        if (static_cast<std::size_t>(request->feed_size()) > registry_.limits().solver.max_components) {
            fail(StatusCode::RESOURCE_EXHAUSTED, "wire.feed_limit", "feed");
        }
        const Admission solve_admission(solves_, limits_.max_concurrent_solves);
        const flash::PtFlashRequest input{request->pressure_pa(), request->temperature_k(),
                                         {request->feed().begin(), request->feed().end()}};
        // Exactly one registry/lease/native solve. Hints are decoded into an
        // ephemeral public DTO and are never retained by the service/registry.
        if (request->has_hints()) {
            const auto hints = decode_solve_hints(request->hints());
            encode_result(registry_.solve(request->model_handle(), input, hints),
                          *response->mutable_result());
        } else {
            encode_result(registry_.solve(request->model_handle(), input),
                          *response->mutable_result());
        }
    });
}
grpc::Status ModelGrpcServiceAdapter::ReleaseModel(grpc::ServerContext* context,
    const wire::ModelHandleRequest* request, wire::ReleaseModelResponse* response) {
    return dispatch(context, request, response, authorize_, limits_, requests_, [&] {
        require_handle(*request);
        registry_.release(request->model_handle());
        // Release is committed even if cancellation/transport failure follows it.
    });
}
void configure_model_grpc_server(grpc::ServerBuilder& builder, ModelGrpcServiceAdapter& adapter) {
    const auto& limits = adapter.limits();
    builder.SetMaxReceiveMessageSize(static_cast<int>(limits.max_request_bytes));
    builder.SetMaxSendMessageSize(static_cast<int>(limits.max_response_bytes));
    grpc::ResourceQuota quota("mpmc-model-configuration");
    quota.Resize(limits.resource_quota_bytes);
    builder.SetResourceQuota(quota);
    builder.RegisterService(&adapter);
}
} // namespace mpmc::model_configuration_grpc
