#ifndef MPMC_PT_PROCESS_PROCESS_HOST_HPP
#define MPMC_PT_PROCESS_PROCESS_HOST_HPP

#include <mpmc/pt_process/composition_root.hpp>
#include <mpmc/runtime_grpc/grpc_headers.hpp>

#include <atomic>
#include <chrono>
#include <cstddef>
#include <filesystem>
#include <memory>
#include <string>

namespace mpmc::pt_process {

inline constexpr std::size_t pt_process_max_pem_file_bytes = 1024U * 1024U;
inline constexpr char pt_flash_grpc_health_service_name[] =
    "mpmc.runtime.v1.PtFlashService";

struct PtProcessTlsIdentity {
    std::string certificate_chain_pem;
    std::string private_key_pem;
    std::string trusted_client_ca_pem;

    [[nodiscard]] bool structurally_valid() const noexcept;
};

struct PtProcessTlsFiles {
    std::filesystem::path certificate_chain;
    std::filesystem::path private_key;
    std::filesystem::path trusted_client_ca;
    std::size_t max_file_bytes{pt_process_max_pem_file_bytes};

    [[nodiscard]] bool structurally_valid() const noexcept;
};

// Reads three bounded PEM files. File names may be reported by the caller, but
// PEM contents are never included in exceptions or observations.
[[nodiscard]] PtProcessTlsIdentity load_pt_process_tls_identity(
    const PtProcessTlsFiles& files);

struct PtProcessHostOptions {
    // Loopback is the safe default when Envoy is the public TLS/CORS edge.
    std::string listen_address{"127.0.0.1:50051"};
    PtProcessTlsIdentity tls;
    std::chrono::seconds shutdown_grace{10};
    // Explicit engineering-desktop mode. It is accepted only on an
    // OS-selected IPv4 loopback port and only when the adapter requires a
    // per-launch bearer token. Production deployment must leave this false.
    bool desktop_loopback_session{false};

    [[nodiscard]] bool structurally_valid() const noexcept;
};

// Synchronous native gRPC host. Production uses mandatory verified client
// certificates. The explicit desktop session alternative is restricted to an
// ephemeral IPv4 loopback listener plus adapter-enforced bearer authentication.
// Lifecycle methods may be called from a process-control thread; start() is a
// one-shot operation and must complete before wait()/shutdown() are invoked.
// The referenced PtCompositionRoot must outlive this host.
class PtProcessHost final {
public:
    PtProcessHost(PtCompositionRoot& root, PtProcessHostOptions options);
    ~PtProcessHost();

    PtProcessHost(const PtProcessHost&) = delete;
    PtProcessHost& operator=(const PtProcessHost&) = delete;
    PtProcessHost(PtProcessHost&&) = delete;
    PtProcessHost& operator=(PtProcessHost&&) = delete;

    void start();
    void wait();
    void shutdown() noexcept;

    [[nodiscard]] bool running() const noexcept;
    [[nodiscard]] int selected_port() const noexcept { return selected_port_; }

private:
    PtCompositionRoot& root_;
    PtProcessHostOptions options_;
    std::unique_ptr<grpc::Server> server_;
    int selected_port_{};
    std::atomic<bool> started_{false};
    std::atomic<bool> shutdown_requested_{false};
    std::atomic<bool> wait_completed_{false};
};

} // namespace mpmc::pt_process

#endif // MPMC_PT_PROCESS_PROCESS_HOST_HPP
