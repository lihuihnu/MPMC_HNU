#include <mpmc/model_configuration_grpc/model_grpc_adapter.hpp>
#include <type_traits>
static_assert(!std::is_copy_constructible_v<mpmc::model_configuration_grpc::ModelGrpcServiceAdapter>);
static_assert(std::is_base_of_v<mpmc::model_configuration::v1::ModelConfigurationService::Service,
                               mpmc::model_configuration_grpc::ModelGrpcServiceAdapter>);
