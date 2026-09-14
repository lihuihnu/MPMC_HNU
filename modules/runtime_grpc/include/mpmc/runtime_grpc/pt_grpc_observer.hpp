#ifndef MPMC_RUNTIME_GRPC_PT_GRPC_OBSERVER_HPP
#define MPMC_RUNTIME_GRPC_PT_GRPC_OBSERVER_HPP

#include <mpmc/runtime/pt_service_contract.hpp>
#include <mpmc/runtime_grpc/grpc_headers.hpp>

#include <chrono>
#include <cstddef>
#include <optional>

namespace mpmc::runtime_grpc {

enum class PtGrpcRpcMethod {
    discover_pt_capabilities,
    solve_pt_flash
};

// Process-level completion only. Scientific acceptance remains exclusively in
// PtServiceOutcome and is never inferred from transport status.
enum class PtGrpcCompletion {
    completed,
    authentication_failure,
    invalid_request,
    request_limit,
    response_limit,
    concurrency_limit,
    deadline_exceeded,
    cancelled,
    memory_exhausted,
    internal_failure
};

struct PtGrpcObservation {
    PtGrpcRpcMethod method{PtGrpcRpcMethod::discover_pt_capabilities};
    PtGrpcCompletion completion{PtGrpcCompletion::internal_failure};
    grpc::StatusCode transport_status{grpc::StatusCode::UNKNOWN};
    // Present only after PtService::solve() returned. `error` is therefore a
    // service outcome, while `indeterminate` remains a scientific result.
    std::optional<::mpmc::runtime::PtServiceOutcome> service_outcome;
    bool pt_service_called{false};
    std::size_t request_bytes{};
    std::size_t response_bytes{};
    std::chrono::microseconds elapsed{};
};

// Observability must never change RPC behavior. Calls may be concurrent;
// implementations are required to synchronize their state, absorb their own
// I/O failures, and avoid logging request compositions or PEMs.
class PtGrpcObserver {
public:
    virtual ~PtGrpcObserver() = default;
    virtual void observe(const PtGrpcObservation& observation) noexcept = 0;
};

} // namespace mpmc::runtime_grpc

#endif // MPMC_RUNTIME_GRPC_PT_GRPC_OBSERVER_HPP
