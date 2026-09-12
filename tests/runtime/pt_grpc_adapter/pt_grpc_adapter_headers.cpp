#include <mpmc/runtime_grpc/pt_grpc_adapter.hpp>
#include <mpmc/runtime_grpc/pt_grpc_observer.hpp>

bool pt_grpc_adapter_headers() {
    using Limits = mpmc::runtime_grpc::PtGrpcAdapterLimits;
    using Observation = mpmc::runtime_grpc::PtGrpcObservation;
    return Limits{}.structurally_valid() &&
           Observation{}.completion ==
               mpmc::runtime_grpc::PtGrpcCompletion::internal_failure;
}
