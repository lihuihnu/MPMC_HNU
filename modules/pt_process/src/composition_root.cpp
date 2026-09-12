#include <mpmc/pt_process/composition_root.hpp>

#include <stdexcept>
#include <utility>

namespace mpmc::pt_process {

std::vector<runtime::PtServiceBackendRegistration>
PtCompositionRoot::registrations_for(
    const std::vector<OwnedConfiguredPtBackend>& owned_backends) {
    std::vector<runtime::PtServiceBackendRegistration> registrations;
    registrations.reserve(owned_backends.size());
    for (const auto& configured : owned_backends) {
        if (!configured.backend) {
            throw std::invalid_argument(
                "PT composition root: null configured backend owner");
        }
        registrations.push_back(
            {configured.configured_backend_id, configured.backend.get()});
    }
    return registrations;
}

PtCompositionRoot::PtCompositionRoot(
    std::vector<OwnedConfiguredPtBackend> backends,
    runtime::PtServiceLimits service_limits,
    runtime_grpc::PtGrpcAdapterLimits adapter_limits,
    std::shared_ptr<runtime_grpc::PtGrpcObserver> observer)
    : owned_backends_(std::move(backends)),
      registrations_(registrations_for(owned_backends_)),
      service_(registrations_, service_limits),
      adapter_(service_, adapter_limits, std::move(observer)) {}

} // namespace mpmc::pt_process
