#ifndef MPMC_MODEL_CONFIGURATION_GRPC_MODEL_GRPC_ADAPTER_HPP
#define MPMC_MODEL_CONFIGURATION_GRPC_MODEL_GRPC_ADAPTER_HPP

#include <mpmc/model_configuration_grpc/model_wire.hpp>
#include <mpmc/model_configuration/v1/model_service.grpc.pb.h>

#include <atomic>
#include <chrono>
#include <functional>

namespace mpmc::model_configuration_grpc {
struct ModelGrpcLimits {
    std::size_t max_request_bytes{1024U * 1024U};
    std::size_t max_response_bytes{4U * 1024U * 1024U};
    std::size_t resource_quota_bytes{64U * 1024U * 1024U};
    std::size_t max_concurrent_requests{4};
    std::size_t max_concurrent_solves{1};
    std::chrono::milliseconds max_deadline{120000};
    [[nodiscard]] bool structurally_valid() const noexcept;
};
// Host must authorize the caller for THIS registry/session, e.g. using existing
// authenticated transport identity. Required, thread-safe; not supplied on wire.
using ModelGrpcAuthorization = std::function<bool(const grpc::ServerContext&)>;
class ModelGrpcServiceAdapter final : public model_configuration::v1::ModelConfigurationService::Service {
public:
    ModelGrpcServiceAdapter(model_configuration::Pr76ModelRegistry& registry,
                           ModelGrpcAuthorization authorize, ModelGrpcLimits limits = {});
    ModelGrpcServiceAdapter(const ModelGrpcServiceAdapter&) = delete;
    ModelGrpcServiceAdapter& operator=(const ModelGrpcServiceAdapter&) = delete;
    [[nodiscard]] const ModelGrpcLimits& limits() const noexcept { return limits_; }
    [[nodiscard]] std::size_t in_flight_requests() const noexcept { return requests_.load(); }
    grpc::Status CreateModel(grpc::ServerContext*, const model_configuration::v1::CreateModelRequest*,
                             model_configuration::v1::CreateModelResponse*) override;
    grpc::Status DescribeModel(grpc::ServerContext*, const model_configuration::v1::ModelHandleRequest*,
                               model_configuration::v1::DescribeModelResponse*) override;
    grpc::Status SolveModel(grpc::ServerContext*, const model_configuration::v1::SolveModelRequest*,
                            model_configuration::v1::SolveModelResponse*) override;
    grpc::Status ReleaseModel(grpc::ServerContext*, const model_configuration::v1::ModelHandleRequest*,
                              model_configuration::v1::ReleaseModelResponse*) override;
private:
    model_configuration::Pr76ModelRegistry& registry_;
    const ModelGrpcAuthorization authorize_;
    const ModelGrpcLimits limits_;
    std::atomic<std::size_t> requests_{0};
    std::atomic<std::size_t> solves_{0};
};
// Configure this service's listener before BuildAndStart: receive cap applies
// before protobuf deserialization. Host owns credentials/listening address and
// must keep registry/adapter alive until server Shutdown + Wait completes.
void configure_model_grpc_server(grpc::ServerBuilder&, ModelGrpcServiceAdapter&);
} // namespace mpmc::model_configuration_grpc
#endif
