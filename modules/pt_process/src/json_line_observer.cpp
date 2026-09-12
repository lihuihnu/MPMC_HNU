#include <mpmc/pt_process/json_line_observer.hpp>

#include <mpmc/runtime/pt_service_contract.hpp>

#include <sstream>
#include <string_view>

namespace mpmc::pt_process {
namespace {

std::string_view method_name(runtime_grpc::PtGrpcRpcMethod method) noexcept {
    switch (method) {
    case runtime_grpc::PtGrpcRpcMethod::discover_pt_capabilities:
        return "discover_pt_capabilities";
    case runtime_grpc::PtGrpcRpcMethod::solve_pt_flash:
        return "solve_pt_flash";
    }
    return "unknown";
}

std::string_view completion_name(
    runtime_grpc::PtGrpcCompletion completion) noexcept {
    switch (completion) {
    case runtime_grpc::PtGrpcCompletion::completed: return "completed";
    case runtime_grpc::PtGrpcCompletion::invalid_request:
        return "invalid_request";
    case runtime_grpc::PtGrpcCompletion::request_limit: return "request_limit";
    case runtime_grpc::PtGrpcCompletion::response_limit:
        return "response_limit";
    case runtime_grpc::PtGrpcCompletion::concurrency_limit:
        return "concurrency_limit";
    case runtime_grpc::PtGrpcCompletion::deadline_exceeded:
        return "deadline_exceeded";
    case runtime_grpc::PtGrpcCompletion::cancelled: return "cancelled";
    case runtime_grpc::PtGrpcCompletion::memory_exhausted:
        return "memory_exhausted";
    case runtime_grpc::PtGrpcCompletion::internal_failure:
        return "internal_failure";
    }
    return "unknown";
}

std::string_view service_outcome_name(
    runtime::PtServiceOutcome outcome) noexcept {
    switch (outcome) {
    case runtime::PtServiceOutcome::accepted: return "accepted";
    case runtime::PtServiceOutcome::phase_set_unstable:
        return "phase_set_unstable";
    case runtime::PtServiceOutcome::indeterminate: return "indeterminate";
    case runtime::PtServiceOutcome::error: return "service_error";
    }
    return "unknown";
}

} // namespace

void PtJsonLineObserver::observe(
    const runtime_grpc::PtGrpcObservation& observation) noexcept {
    try {
        std::ostringstream line;
        line << "{\"event\":\"pt_grpc_rpc\",\"method\":\""
             << method_name(observation.method) << "\",\"completion\":\""
             << completion_name(observation.completion)
             << "\",\"grpc_status\":"
             << static_cast<int>(observation.transport_status)
             << ",\"pt_service_called\":"
             << (observation.pt_service_called ? "true" : "false")
             << ",\"service_outcome\":";
        if (observation.service_outcome) {
            line << '\"' << service_outcome_name(*observation.service_outcome)
                 << '\"';
        } else {
            line << "null";
        }
        line << ",\"request_bytes\":" << observation.request_bytes
             << ",\"response_bytes\":" << observation.response_bytes
             << ",\"elapsed_us\":" << observation.elapsed.count() << "}\n";
        const auto record = line.str();
        const std::lock_guard lock(mutex_);
        output_ << record;
        output_.flush();
    } catch (...) {
        // Observability is fail-open with respect to RPC semantics.
    }
}

} // namespace mpmc::pt_process
