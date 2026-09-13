#if !defined(pascal)
#define pascal mpmc_windows_sdk_pascal_macro_probe
#endif

#include <mpmc/runtime_grpc/pt_grpc_adapter.hpp>
#include <mpmc/runtime_grpc/pt_grpc_observer.hpp>
#include <mpmc/thermodynamics/components.hpp>

#if defined(pascal)
#error "the gRPC boundary leaked the Windows SDK pascal macro"
#endif

bool pt_grpc_adapter_headers() {
    using Limits = mpmc::runtime_grpc::PtGrpcAdapterLimits;
    using Observation = mpmc::runtime_grpc::PtGrpcObservation;
    return Limits{}.structurally_valid() &&
           Observation{}.completion ==
               mpmc::runtime_grpc::PtGrpcCompletion::internal_failure &&
           mpmc::thermodynamics::Unit::pascal !=
               mpmc::thermodynamics::Unit::unspecified;
}
