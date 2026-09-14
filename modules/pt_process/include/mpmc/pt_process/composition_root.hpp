#ifndef MPMC_PT_PROCESS_COMPOSITION_ROOT_HPP
#define MPMC_PT_PROCESS_COMPOSITION_ROOT_HPP

#include <mpmc/flash/pt_flash_backend.hpp>
#include <mpmc/runtime/pt_service.hpp>
#include <mpmc/runtime_grpc/pt_grpc_adapter.hpp>

#include <cstddef>
#include <memory>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace mpmc::pt_process {

// Owns the complete lifetime behind one configured backend registration. The
// shared pointer may be an aliasing pointer into a larger model/evaluator/backend
// object graph; PtCompositionRoot keeps that graph alive for PtService.
struct OwnedConfiguredPtBackend {
    std::string configured_backend_id;
    std::shared_ptr<flash::PtFlashBackend> backend;
};

template <typename Owner>
[[nodiscard]] OwnedConfiguredPtBackend retain_configured_pt_backend(
    std::string configured_backend_id, std::shared_ptr<Owner> owner,
    flash::PtFlashBackend& backend) {
    if (!owner) {
        throw std::invalid_argument(
            "PT composition root: null backend object-graph owner");
    }
    return {std::move(configured_backend_id),
            std::shared_ptr<flash::PtFlashBackend>(std::move(owner), &backend)};
}

// Process-level composition only: immutable backend ownership -> PtService ->
// gRPC adapter. No EOS, root selection, flash retry, or result reinterpretation
// is implemented here.
class PtCompositionRoot final {
public:
    explicit PtCompositionRoot(
        std::vector<OwnedConfiguredPtBackend> backends,
        runtime::PtServiceLimits service_limits = {},
        runtime_grpc::PtGrpcAdapterLimits adapter_limits = {},
        std::shared_ptr<runtime_grpc::PtGrpcObserver> observer = {},
        runtime_grpc::PtGrpcAuthenticationOptions authentication = {});

    PtCompositionRoot(const PtCompositionRoot&) = delete;
    PtCompositionRoot& operator=(const PtCompositionRoot&) = delete;
    PtCompositionRoot(PtCompositionRoot&&) = delete;
    PtCompositionRoot& operator=(PtCompositionRoot&&) = delete;

    [[nodiscard]] runtime::PtService& service() noexcept { return service_; }
    [[nodiscard]] const runtime::PtService& service() const noexcept {
        return service_;
    }
    [[nodiscard]] runtime_grpc::PtGrpcServiceAdapter& adapter() noexcept {
        return adapter_;
    }
    [[nodiscard]] const runtime_grpc::PtGrpcServiceAdapter& adapter()
        const noexcept {
        return adapter_;
    }
    [[nodiscard]] std::size_t configured_backend_count() const noexcept {
        return owned_backends_.size();
    }

private:
    [[nodiscard]] static std::vector<runtime::PtServiceBackendRegistration>
    registrations_for(
        const std::vector<OwnedConfiguredPtBackend>& owned_backends);

    // Destruction is intentionally the reverse of this declaration order.
    std::vector<OwnedConfiguredPtBackend> owned_backends_;
    std::vector<runtime::PtServiceBackendRegistration> registrations_;
    runtime::PtService service_;
    runtime_grpc::PtGrpcServiceAdapter adapter_;
};

} // namespace mpmc::pt_process

#endif // MPMC_PT_PROCESS_COMPOSITION_ROOT_HPP
