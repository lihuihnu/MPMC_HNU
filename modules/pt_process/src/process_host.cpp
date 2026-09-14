#include <mpmc/pt_process/process_host.hpp>

#include <mpmc/runtime_grpc/pt_grpc_adapter.hpp>

#include <algorithm>
#include <fstream>
#include <limits>
#include <mutex>
#include <stdexcept>
#include <utility>

namespace mpmc::pt_process {
namespace {

bool valid_pem_field(const std::string& value) noexcept {
    return !value.empty() && value.size() <= pt_process_max_pem_file_bytes &&
           value.find('\0') == std::string::npos;
}

std::string read_bounded_file(const std::filesystem::path& path,
                              std::size_t max_bytes) {
    std::ifstream input(path, std::ios::binary | std::ios::ate);
    if (!input) {
        throw std::runtime_error("PT process TLS: cannot open configured PEM file");
    }
    const auto end = input.tellg();
    if (end < 0 || static_cast<unsigned long long>(end) > max_bytes) {
        throw std::runtime_error(
            "PT process TLS: configured PEM file exceeds the byte limit");
    }
    std::string result(static_cast<std::size_t>(end), '\0');
    input.seekg(0, std::ios::beg);
    if (!result.empty() &&
        !input.read(result.data(), static_cast<std::streamsize>(result.size()))) {
        throw std::runtime_error("PT process TLS: cannot read configured PEM file");
    }
    if (result.empty()) {
        throw std::runtime_error("PT process TLS: configured PEM file is empty");
    }
    return result;
}

void enable_health_service_once() {
    static std::once_flag once;
    std::call_once(once, [] { grpc::EnableDefaultHealthCheckService(true); });
}

bool set_health(grpc::Server& server, bool serving) noexcept {
    try {
        if (auto* health = server.GetHealthCheckService(); health != nullptr) {
            health->SetServingStatus("", serving);
            health->SetServingStatus(pt_flash_grpc_health_service_name, serving);
            return true;
        }
    } catch (...) {
    }
    return false;
}

} // namespace

bool PtProcessTlsIdentity::structurally_valid() const noexcept {
    return valid_pem_field(certificate_chain_pem) &&
           valid_pem_field(private_key_pem) &&
           valid_pem_field(trusted_client_ca_pem);
}

bool PtProcessTlsFiles::structurally_valid() const noexcept {
    return !certificate_chain.empty() && !private_key.empty() &&
           !trusted_client_ca.empty() && max_file_bytes > 0U &&
           max_file_bytes <= pt_process_max_pem_file_bytes &&
           max_file_bytes <=
               static_cast<std::size_t>(
                   std::numeric_limits<std::streamsize>::max());
}

PtProcessTlsIdentity load_pt_process_tls_identity(
    const PtProcessTlsFiles& files) {
    if (!files.structurally_valid()) {
        throw std::invalid_argument("PT process TLS: invalid PEM file policy");
    }
    PtProcessTlsIdentity identity{
        read_bounded_file(files.certificate_chain, files.max_file_bytes),
        read_bounded_file(files.private_key, files.max_file_bytes),
        read_bounded_file(files.trusted_client_ca, files.max_file_bytes)};
    if (!identity.structurally_valid()) {
        throw std::runtime_error("PT process TLS: invalid PEM payload structure");
    }
    return identity;
}

bool PtProcessHostOptions::structurally_valid() const noexcept {
    const bool address_valid =
        !listen_address.empty() && listen_address.size() <= 1024U &&
        std::all_of(listen_address.begin(), listen_address.end(), [](char value) {
            const auto byte = static_cast<unsigned char>(value);
            return byte > 0x20U && byte <= 0x7eU;
        });
    const bool tls_empty = tls.certificate_chain_pem.empty() &&
                           tls.private_key_pem.empty() &&
                           tls.trusted_client_ca_pem.empty();
    const bool security_valid =
        desktop_loopback_session
            ? listen_address == "127.0.0.1:0" && tls_empty
            : tls.structurally_valid();
    return address_valid && security_valid && shutdown_grace.count() > 0 &&
           shutdown_grace <= std::chrono::minutes(5);
}

PtProcessHost::PtProcessHost(PtCompositionRoot& root,
                             PtProcessHostOptions options)
    : root_(root), options_(std::move(options)) {
    if (!options_.structurally_valid()) {
        throw std::invalid_argument("invalid secure PT process host options");
    }
    if (root_.configured_backend_count() == 0U) {
        throw std::invalid_argument(
            "secure PT process host requires a configured backend");
    }
    if (options_.desktop_loopback_session &&
        !root_.adapter().requires_bearer_token()) {
        throw std::invalid_argument(
            "desktop PT process host requires adapter bearer authentication");
    }
    if (!options_.desktop_loopback_session &&
        root_.adapter().requires_bearer_token()) {
        throw std::invalid_argument(
            "production PT process host cannot use desktop bearer authentication");
    }
}

PtProcessHost::~PtProcessHost() {
    shutdown();
    if (server_ && !wait_completed_.load(std::memory_order_acquire)) {
        server_->Wait();
    }
}

void PtProcessHost::start() {
    if (started_.exchange(true, std::memory_order_acq_rel)) {
        throw std::logic_error("PT process host can only be started once");
    }

    enable_health_service_once();
    grpc::ServerBuilder builder;
    if (options_.desktop_loopback_session) {
        builder.AddListeningPort(options_.listen_address,
                                 grpc::InsecureServerCredentials(),
                                 &selected_port_);
    } else {
        grpc::SslServerCredentialsOptions tls_options(
            GRPC_SSL_REQUEST_AND_REQUIRE_CLIENT_CERTIFICATE_AND_VERIFY);
        tls_options.pem_root_certs = options_.tls.trusted_client_ca_pem;
        tls_options.pem_key_cert_pairs.push_back(
            {options_.tls.private_key_pem,
             options_.tls.certificate_chain_pem});
        builder.AddListeningPort(options_.listen_address,
                                 grpc::SslServerCredentials(tls_options),
                                 &selected_port_);
    }
    runtime_grpc::configure_pt_grpc_server(builder, root_.adapter());
    server_ = builder.BuildAndStart();
    if (!server_ || selected_port_ <= 0) {
        started_.store(false, std::memory_order_release);
        throw std::runtime_error(
            "secure PT process host could not bind its configured address");
    }
    if (!set_health(*server_, true)) {
        server_->Shutdown();
        server_->Wait();
        server_.reset();
        selected_port_ = 0;
        started_.store(false, std::memory_order_release);
        throw std::runtime_error(
            "secure PT process host could not publish gRPC health status");
    }
}

void PtProcessHost::wait() {
    if (!server_) {
        throw std::logic_error("PT process host is not started");
    }
    server_->Wait();
    wait_completed_.store(true, std::memory_order_release);
}

void PtProcessHost::shutdown() noexcept {
    if (!server_ ||
        shutdown_requested_.exchange(true, std::memory_order_acq_rel)) {
        return;
    }
    (void)set_health(*server_, false);
    server_->Shutdown(std::chrono::system_clock::now() +
                      options_.shutdown_grace);
}

bool PtProcessHost::running() const noexcept {
    return server_ != nullptr &&
           !shutdown_requested_.load(std::memory_order_acquire);
}

} // namespace mpmc::pt_process
