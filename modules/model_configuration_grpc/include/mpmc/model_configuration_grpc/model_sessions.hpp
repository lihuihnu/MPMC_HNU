#ifndef MPMC_MODEL_CONFIGURATION_GRPC_MODEL_SESSIONS_HPP
#define MPMC_MODEL_CONFIGURATION_GRPC_MODEL_SESSIONS_HPP
#include <mpmc/model_configuration_grpc/model_grpc_adapter.hpp>

#include <functional>
#include <memory>
#include <string>

namespace mpmc::model_configuration_grpc {
inline constexpr std::string_view model_session_contract = "mpmc.model_configuration.v1/model-session/v1";
inline constexpr char model_session_metadata[] = "x-mpmc-model-session";

struct ModelSessionLimits {
    std::size_t max_sessions{16};
    std::size_t max_models_per_session{4};
    std::chrono::milliseconds max_session_lifetime{1'800'000};
    // Also enforced globally across sessions. Open streams have their own cap.
    ModelGrpcLimits rpc{64U * 1024U};
};
struct ModelSessionStatus {
    std::size_t active_sessions{};
    std::size_t resident_sessions{}; // Includes closed sessions with retained calls.
    std::size_t in_flight_requests{};
    bool closed{};
};
// Required host policy: return a bounded stable authenticated principal, or an
// empty string to deny. Never derive identity from client-asserted metadata/peer
// socket address. The principal is kept private and never sent in errors/logs.
using ModelSessionIdentity = std::function<std::string(const grpc::ServerContext&)>;

class ModelSessionService final
    : public ::mpmc::model_configuration::v1::ModelConfigurationService::Service,
      public ::mpmc::model_configuration::v1::ModelSessionService::Service {
public:
    explicit ModelSessionService(ModelSessionIdentity identity, ModelSessionLimits limits = {},
        ::mpmc::model_configuration::ModelHandleEntropySource entropy =
            ::mpmc::model_configuration::system_model_handle_entropy);
    ~ModelSessionService();
    ModelSessionService(const ModelSessionService&) = delete;
    ModelSessionService& operator=(const ModelSessionService&) = delete;
    // Close before server Shutdown/Wait; destroy only after all handlers return.
    void close();
    [[nodiscard]] ModelSessionStatus status() const;
    [[nodiscard]] const ModelSessionLimits& limits() const noexcept;
    // Registers both services only. The host must apply compatible pre-parse
    // message/resource quotas to the listener; it owns TLS and shutdown.
    void register_services(grpc::ServerBuilder& builder);
    grpc::Status OpenModelSession(grpc::ServerContext*,
        const ::mpmc::model_configuration::v1::OpenModelSessionRequest*,
        grpc::ServerWriter<::mpmc::model_configuration::v1::ModelSessionOpened>*) override;
    grpc::Status CreateModel(grpc::ServerContext*, const ::mpmc::model_configuration::v1::CreateModelRequest*,
        ::mpmc::model_configuration::v1::CreateModelResponse*) override;
    grpc::Status DescribeModel(grpc::ServerContext*, const ::mpmc::model_configuration::v1::ModelHandleRequest*,
        ::mpmc::model_configuration::v1::DescribeModelResponse*) override;
    grpc::Status SolveModel(grpc::ServerContext*, const ::mpmc::model_configuration::v1::SolveModelRequest*,
        ::mpmc::model_configuration::v1::SolveModelResponse*) override;
    grpc::Status ReleaseModel(grpc::ServerContext*, const ::mpmc::model_configuration::v1::ModelHandleRequest*,
        ::mpmc::model_configuration::v1::ReleaseModelResponse*) override;
private:
    struct Impl;
    std::unique_ptr<Impl> impl_;
};
} // namespace mpmc::model_configuration_grpc
#endif
