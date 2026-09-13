#include <mpmc/pt_process/composition_root.hpp>
#include <mpmc/pt_process/json_line_observer.hpp>
#include <mpmc/pt_process/parameter_snapshot_supplier.hpp>
#include <mpmc/pt_process/process_host.hpp>

#include <atomic>
#include <charconv>
#include <chrono>
#include <csignal>
#include <filesystem>
#include <iostream>
#include <memory>
#include <stdexcept>
#include <string>
#include <string_view>
#include <system_error>
#include <thread>
#include <utility>

namespace {

namespace process = mpmc::pt_process;

constexpr std::string_view default_backend_certificate =
    "/run/secrets/mpmc-pt/backend/tls.crt";
constexpr std::string_view default_backend_private_key =
    "/run/secrets/mpmc-pt/backend/tls.key";
constexpr std::string_view default_edge_client_ca =
    "/run/secrets/mpmc-pt/backend/edge-client-ca.pem";

volatile std::sig_atomic_t stop_requested = 0;

extern "C" void request_stop(int) {
    stop_requested = 1;
}

struct CommandLine {
    std::string listen_address{"127.0.0.1:50051"};
    std::filesystem::path certificate_chain{default_backend_certificate};
    std::filesystem::path private_key{default_backend_private_key};
    std::filesystem::path trusted_edge_client_ca{default_edge_client_ca};
    std::chrono::seconds shutdown_grace{10};
    bool print_snapshot_manifest{};
    bool help{};
};

[[noreturn]] void argument_error(std::string_view message) {
    throw std::invalid_argument("PT production host arguments: " +
                                std::string(message));
}

std::string require_value(
    int argc, char** argv, int& index, std::string_view option) {
    if (index + 1 >= argc) {
        argument_error(std::string(option) + " requires a value");
    }
    ++index;
    const std::string value(argv[index]);
    if (value.empty() || value.find('\0') != std::string::npos) {
        argument_error(std::string(option) + " requires nonempty text");
    }
    return value;
}

std::chrono::seconds parse_shutdown_grace(std::string_view value) {
    unsigned int seconds{};
    const auto result = std::from_chars(
        value.data(), value.data() + value.size(), seconds);
    if (result.ec != std::errc{} ||
        result.ptr != value.data() + value.size() || seconds == 0U ||
        seconds > 300U) {
        argument_error(
            "--shutdown-grace-seconds must be an integer in [1,300]");
    }
    return std::chrono::seconds(seconds);
}

CommandLine parse_command_line(int argc, char** argv) {
    CommandLine options;
    bool listen_seen = false;
    bool certificate_seen = false;
    bool key_seen = false;
    bool client_ca_seen = false;
    bool grace_seen = false;
    for (int index = 1; index < argc; ++index) {
        const std::string_view option(argv[index]);
        if (option == "--help") {
            options.help = true;
        } else if (option == "--print-snapshot-manifest") {
            options.print_snapshot_manifest = true;
        } else if (option == "--listen-address") {
            if (listen_seen) {
                argument_error("duplicate --listen-address");
            }
            listen_seen = true;
            options.listen_address =
                require_value(argc, argv, index, option);
        } else if (option == "--tls-certificate-chain") {
            if (certificate_seen) {
                argument_error("duplicate --tls-certificate-chain");
            }
            certificate_seen = true;
            options.certificate_chain =
                require_value(argc, argv, index, option);
        } else if (option == "--tls-private-key") {
            if (key_seen) {
                argument_error("duplicate --tls-private-key");
            }
            key_seen = true;
            options.private_key = require_value(argc, argv, index, option);
        } else if (option == "--trusted-edge-client-ca") {
            if (client_ca_seen) {
                argument_error("duplicate --trusted-edge-client-ca");
            }
            client_ca_seen = true;
            options.trusted_edge_client_ca =
                require_value(argc, argv, index, option);
        } else if (option == "--shutdown-grace-seconds") {
            if (grace_seen) {
                argument_error("duplicate --shutdown-grace-seconds");
            }
            grace_seen = true;
            options.shutdown_grace = parse_shutdown_grace(
                require_value(argc, argv, index, option));
        } else {
            argument_error("unknown option: " + std::string(option));
        }
    }
    return options;
}

void validate_secret_paths(const CommandLine& options) {
    for (const auto* path : {
             &options.certificate_chain,
             &options.private_key,
             &options.trusted_edge_client_ca}) {
        if (!path->is_absolute()) {
            argument_error("TLS secret paths must be absolute");
        }
    }
}

void print_help(std::ostream& output) {
    output
        << "Usage: mpmc_pt_service_host [options]\n"
        << "  --listen-address ADDRESS\n"
        << "  --tls-certificate-chain ABSOLUTE_PATH\n"
        << "  --tls-private-key ABSOLUTE_PATH\n"
        << "  --trusted-edge-client-ca ABSOLUTE_PATH\n"
        << "  --shutdown-grace-seconds 1..300\n"
        << "  --print-snapshot-manifest\n"
        << "\nDefaults use the documented /run/secrets/mpmc-pt/backend "
           "mount contract.\n";
}

void install_signal_handlers() {
    if (std::signal(SIGINT, request_stop) == SIG_ERR ||
        std::signal(SIGTERM, request_stop) == SIG_ERR) {
        throw std::runtime_error(
            "PT production host could not install stop signal handlers");
    }
#ifdef _WIN32
    if (std::signal(SIGBREAK, request_stop) == SIG_ERR) {
        throw std::runtime_error(
            "PT production host could not install the Windows stop handler");
    }
#endif
}

void write_json_string(std::ostream& output, std::string_view value) {
    output << '"';
    for (const char raw_character : value) {
        const auto character = static_cast<unsigned char>(raw_character);
        switch (character) {
        case '"':
            output << "\\\"";
            break;
        case '\\':
            output << "\\\\";
            break;
        case '\n':
            output << "\\n";
            break;
        case '\r':
            output << "\\r";
            break;
        case '\t':
            output << "\\t";
            break;
        default:
            if (character < 0x20U) {
                output << '?';
            } else {
                output << raw_character;
            }
        }
    }
    output << '"';
}

} // namespace

int main(int argc, char** argv) {
    try {
        const auto command_line = parse_command_line(argc, argv);
        if (command_line.help) {
            print_help(std::cout);
            return 0;
        }

        auto snapshots =
            process::load_repository_curated_pt_parameter_snapshots_v1();
        if (command_line.print_snapshot_manifest) {
            process::write_pt_parameter_snapshot_manifest_json(
                std::cout, snapshots);
            return 0;
        }

        validate_secret_paths(command_line);
        auto observer =
            std::make_shared<process::PtJsonLineObserver>(std::cout);
        process::PtCompositionRoot root(
            std::move(snapshots.backends), {}, {}, std::move(observer));

        process::PtProcessHostOptions host_options;
        host_options.listen_address = command_line.listen_address;
        host_options.tls = process::load_pt_process_tls_identity({
            command_line.certificate_chain,
            command_line.private_key,
            command_line.trusted_edge_client_ca});
        host_options.shutdown_grace = command_line.shutdown_grace;

        process::PtProcessHost host(root, std::move(host_options));
        install_signal_handlers();
        host.start();
        std::cout << "{\"event\":\"pt_process_ready\","
                  << "\"snapshot_bundle\":\""
                  << process::repository_curated_pt_bundle_id
                  << "\",\"snapshot_revision\":\""
                  << process::repository_curated_pt_bundle_revision
                  << "\",\"configured_backends\":"
                  << root.configured_backend_count()
                  << ",\"selected_port\":" << host.selected_port()
                  << "}\n" << std::flush;

        std::atomic<bool> wait_completed{false};
        std::thread waiter([&host, &wait_completed] {
            host.wait();
            wait_completed.store(true, std::memory_order_release);
        });
        while (stop_requested == 0 &&
               !wait_completed.load(std::memory_order_acquire)) {
            std::this_thread::sleep_for(std::chrono::milliseconds(50));
        }
        if (stop_requested != 0) {
            host.shutdown();
        }
        waiter.join();
        std::cout << "{\"event\":\"pt_process_stopped\"}\n"
                  << std::flush;
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "{\"event\":\"pt_process_start_failed\","
                  << "\"diagnostic\":";
        write_json_string(
            std::cerr,
            "configuration or startup failure: " +
                std::string(error.what()));
        std::cerr << "}\n";
        return 1;
    }
}
