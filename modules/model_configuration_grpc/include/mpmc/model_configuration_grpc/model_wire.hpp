#ifndef MPMC_MODEL_CONFIGURATION_GRPC_MODEL_WIRE_HPP
#define MPMC_MODEL_CONFIGURATION_GRPC_MODEL_WIRE_HPP

#include <mpmc/runtime_grpc/grpc_headers.hpp>
#include <mpmc/model_configuration/pr76_model_registry.hpp>
#include <mpmc/model_configuration/v1/model_service.pb.h>

namespace mpmc::model_configuration_grpc {
inline constexpr std::string_view model_wire_contract = "mpmc.model_configuration.v1/model-service/v1";
// Standard gRPC rich-error envelope containing one typed ModelServiceError.
[[nodiscard]] grpc::Status model_error_status(grpc::StatusCode code, std::string_view message,
    const model_configuration::v1::ModelServiceError& detail);

// Owning DTO codec; server calls occur only after bounded message admission.
// Preserves optional numerical presence; domain rules stay in the factory.
[[nodiscard]] model_configuration::ThermodynamicModelDefinition decode_definition(
    const model_configuration::v1::ThermodynamicModelDefinition& input,
    const model_configuration::ModelConfigurationLimits& limits);
[[nodiscard]] model_configuration::PtSolverSettings decode_settings(
    const model_configuration::v1::PtSolverSettings& input);
// Structural wire mapping only. State-dependent support/simplex/resource checks
// stay in Pr76ExecutableModel so direct C++ and RPC calls share one authority.
[[nodiscard]] model_configuration::PtSolveHints decode_solve_hints(
    const model_configuration::v1::PtSolveHints& input);
void encode_definition(const model_configuration::ThermodynamicModelDefinition& input,
                       model_configuration::v1::ThermodynamicModelDefinition& output);
void encode_settings(const model_configuration::PtSolverSettings& input,
                     model_configuration::v1::PtSolverSettings& output);
void encode_solve_hints(const model_configuration::PtSolveHints& input,
                        model_configuration::v1::PtSolveHints& output);
void encode_snapshot(const model_configuration::Pr76RegisteredModelSnapshot& input,
                     model_configuration::v1::ModelSnapshot& output);
void encode_result(const flash::PtFlashBackendResult& input,
                   model_configuration::v1::FullPtResult& output);
} // namespace mpmc::model_configuration_grpc
#endif
