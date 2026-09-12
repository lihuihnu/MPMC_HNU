#include <mpmc/runtime_grpc/pt_grpc_adapter.hpp>

bool pt_grpc_adapter_headers() {
    using Limits = mpmc::runtime_grpc::PtGrpcAdapterLimits;
    return Limits{}.structurally_valid();
}
