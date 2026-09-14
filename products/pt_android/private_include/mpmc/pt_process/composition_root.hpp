#ifndef MPMC_PT_ANDROID_PRIVATE_COMPOSITION_ROOT_HPP
#define MPMC_PT_ANDROID_PRIVATE_COMPOSITION_ROOT_HPP

#include <mpmc/flash/pt_flash_backend.hpp>

#include <memory>
#include <stdexcept>
#include <string>
#include <utility>

namespace mpmc::pt_process {

// Android-private compile adapter for the two canonical repository-curated
// snapshot translation units. It intentionally mirrors only the backend
// ownership primitive from the process composition root. The desktop/process
// PtCompositionRoot and runtime_grpc adapter remain outside the Android graph.
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
            "PT Android ownership shim: null backend object-graph owner");
    }
    return {std::move(configured_backend_id),
            std::shared_ptr<flash::PtFlashBackend>(std::move(owner), &backend)};
}

} // namespace mpmc::pt_process

#endif // MPMC_PT_ANDROID_PRIVATE_COMPOSITION_ROOT_HPP
