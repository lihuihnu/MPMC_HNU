#ifndef MPMC_PT_PROCESS_JSON_LINE_OBSERVER_HPP
#define MPMC_PT_PROCESS_JSON_LINE_OBSERVER_HPP

#include <mpmc/runtime_grpc/pt_grpc_observer.hpp>

#include <mutex>
#include <ostream>

namespace mpmc::pt_process {

// Emits one bounded, model-neutral JSON object per completed RPC. Request state,
// compositions, parameter values, credentials, and diagnostics are deliberately
// excluded. The referenced stream must outlive this observer.
class PtJsonLineObserver final : public runtime_grpc::PtGrpcObserver {
public:
    explicit PtJsonLineObserver(std::ostream& output) noexcept
        : output_(output) {}

    void observe(
        const runtime_grpc::PtGrpcObservation& observation) noexcept override;

private:
    std::ostream& output_;
    std::mutex mutex_;
};

} // namespace mpmc::pt_process

#endif // MPMC_PT_PROCESS_JSON_LINE_OBSERVER_HPP
