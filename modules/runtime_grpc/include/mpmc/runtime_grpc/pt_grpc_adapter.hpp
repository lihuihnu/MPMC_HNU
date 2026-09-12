#ifndef MPMC_RUNTIME_GRPC_PT_GRPC_ADAPTER_HPP
#define MPMC_RUNTIME_GRPC_PT_GRPC_ADAPTER_HPP

#include <mpmc/runtime/pt_service.hpp>

#include <grpcpp/grpcpp.h>
#include <mpmc/runtime/v1/pt_service.grpc.pb.h>

#include <atomic>
#include <chrono>
#include <cstddef>

namespace mpmc::runtime_grpc {

// Process-level limits are deliberately separate from PtServiceLimits. The
// latter protects the transport-neutral registry and request semantics; these
// values bound untrusted RPC transport, scheduling, and publication resources.
struct PtGrpcAdapterLimits {
    std::size_t max_serialized_request_bytes{64U * 1024U};
    std::size_t max_serialized_response_bytes{4U * 1024U * 1024U};
    std::size_t max_concurrent_solves{1U};
    std::size_t grpc_resource_quota_bytes{64U * 1024U * 1024U};
    std::chrono::milliseconds max_discovery_deadline{10'000};
    std::chrono::milliseconds max_solve_deadline{120'000};

    [[nodiscard]] bool structurally_valid() const noexcept;
};

// Thin process adapter for the versioned wire contract. It validates only wire
// presence/version and process resource policy, converts the request without
// normalization, and calls PtService::solve exactly once. It never calls a
// model-specific backend or recomputes a scientific acceptance decision.
class PtGrpcServiceAdapter final
    : public ::mpmc::runtime::v1::PtFlashService::Service {
public:
    explicit PtGrpcServiceAdapter(
        ::mpmc::runtime::PtService& service,
        PtGrpcAdapterLimits limits = {});

    PtGrpcServiceAdapter(const PtGrpcServiceAdapter&) = delete;
    PtGrpcServiceAdapter& operator=(const PtGrpcServiceAdapter&) = delete;

    [[nodiscard]] const PtGrpcAdapterLimits& limits() const noexcept {
        return limits_;
    }

    [[nodiscard]] std::size_t in_flight_solves() const noexcept {
        return in_flight_solves_.load(std::memory_order_relaxed);
    }

    grpc::Status DiscoverPtCapabilities(
        grpc::ServerContext* context,
        const ::mpmc::runtime::v1::DiscoverPtCapabilitiesRequest* request,
        ::mpmc::runtime::v1::DiscoverPtCapabilitiesResponse* response) override;

    grpc::Status SolvePtFlash(
        grpc::ServerContext* context,
        const ::mpmc::runtime::v1::SolvePtFlashRequest* request,
        ::mpmc::runtime::v1::SolvePtFlashResponse* response) override;

private:
    [[nodiscard]] bool try_acquire_solve() noexcept;
    void release_solve() noexcept;

    ::mpmc::runtime::PtService& service_;
    PtGrpcAdapterLimits limits_;
    std::atomic<std::size_t> in_flight_solves_{0U};
};

// Applies the same message-size limits before Protobuf request deserialization,
// attaches a bounded gRPC memory quota, and registers the adapter.
// Listening address, credentials, TLS, and process lifetime remain application
// choices.
void configure_pt_grpc_server(grpc::ServerBuilder& builder,
                              PtGrpcServiceAdapter& adapter);

} // namespace mpmc::runtime_grpc

#endif // MPMC_RUNTIME_GRPC_PT_GRPC_ADAPTER_HPP
