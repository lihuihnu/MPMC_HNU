#include <mpmc/pt_process/composition_root.hpp>
#include <mpmc/pt_process/configured_backends.hpp>
#include <mpmc/pt_process/json_line_observer.hpp>
#include <mpmc/pt_process/process_host.hpp>

#include "../pt_grpc_adapter/test_backend.hpp"
#include "../../flash/cpa_max3/test_support.hpp"
#include "../../flash/pr76_three_phase/synthetic_fixture.hpp"
#include "test_support.hpp"

#include <grpcpp/create_channel.h>
#include <grpcpp/security/credentials.h>
#include <grpcpp/support/channel_arguments.h>
#include <mpmc/runtime/v1/pt_service.grpc.pb.h>

#include <chrono>
#include <fstream>
#include <iostream>
#include <iterator>
#include <memory>
#include <sstream>
#include <source_location>
#include <stdexcept>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

bool pt_process_headers();

namespace {

namespace process = mpmc::pt_process;
namespace runtime = mpmc::runtime;
namespace wire = mpmc::runtime::v1;
using namespace std::chrono_literals;

void require(bool condition, std::string_view message,
             std::source_location where = std::source_location::current()) {
    if (!condition) {
        throw std::runtime_error(std::string(where.file_name()) + ":" +
                                 std::to_string(where.line()) + ": " +
                                 std::string(message));
    }
}

template <typename Error, typename Function>
void expect_error(Function&& function) {
    try {
        function();
    } catch (const Error&) {
        return;
    }
    throw std::runtime_error("expected exception missing");
}

std::string read_file(const char* path) {
    std::ifstream input(path, std::ios::binary);
    if (!input) { throw std::runtime_error("cannot read test certificate"); }
    return {std::istreambuf_iterator<char>(input),
            std::istreambuf_iterator<char>()};
}

std::shared_ptr<grpc::ChannelCredentials> client_credentials(
    const std::string& server_ca, const std::string& certificate = {},
    const std::string& private_key = {}) {
    grpc::SslCredentialsOptions options;
    options.pem_root_certs = server_ca;
    options.pem_cert_chain = certificate;
    options.pem_private_key = private_key;
    return grpc::SslCredentials(options);
}

std::unique_ptr<wire::PtFlashService::Stub> make_stub(
    int port, const std::shared_ptr<grpc::ChannelCredentials>& credentials) {
    grpc::ChannelArguments arguments;
    arguments.SetSslTargetNameOverride("pt-backend.internal");
    auto channel = grpc::CreateCustomChannel(
        "127.0.0.1:" + std::to_string(port), credentials, arguments);
    return wire::PtFlashService::NewStub(channel);
}

grpc::Status discover(wire::PtFlashService::Stub& stub,
                      wire::DiscoverPtCapabilitiesResponse& response) {
    grpc::ClientContext context;
    context.set_deadline(std::chrono::system_clock::now() + 5s);
    wire::DiscoverPtCapabilitiesRequest request;
    return stub.DiscoverPtCapabilities(&context, request, &response);
}

void composition_root() {
    std::vector<process::OwnedConfiguredPtBackend> configured;
    {
        auto pr_model = pr76_max3_test::model();
        const auto pr_parameters = pr_model.parameters();
        const auto sw_parameters =
            sw92_test::binary_parameters(sw92_test::co2);
        const auto cpa_parameters = cpa_max3_test::parameters();
        configured.push_back(process::make_pr76_configured_backend(
            "pr76.fixture", pr_parameters));
        configured.push_back(process::make_sw92_configured_backend(
            "sw92.fixture", sw_parameters));
        configured.push_back(process::make_cpa_configured_backend(
            "cpa.fixture", cpa_parameters,
            cpa_max3_test::fast_pt_options()));
    }

    process::PtCompositionRoot root(std::move(configured));
    require(root.configured_backend_count() == 3U,
            "composition root did not retain three configured backends");
    const auto capabilities = root.service().discover_capabilities();
    require(capabilities.size() == 3U,
            "PtService did not receive all configured backends");
    require(root.service().find_capability("pr76.fixture") != nullptr &&
                root.service().find_capability("sw92.fixture") != nullptr &&
                root.service().find_capability("cpa.fixture") != nullptr,
            "typed backend factory registration was not discoverable");
    for (const auto& descriptor : capabilities) {
        require(descriptor.structurally_valid() &&
                    descriptor.capability.supported_phase_counts ==
                        std::vector<std::size_t>({1U, 2U, 3U}),
                "composition root changed frozen backend capability");
    }

    std::vector<process::OwnedConfiguredPtBackend> invalid{{"null", nullptr}};
    expect_error<std::invalid_argument>([&] {
        process::PtCompositionRoot rejected(std::move(invalid));
        (void)rejected;
    });
}

void observer() {
    std::ostringstream output;
    process::PtJsonLineObserver observer(output);
    observer.observe({mpmc::runtime_grpc::PtGrpcRpcMethod::solve_pt_flash,
                      mpmc::runtime_grpc::PtGrpcCompletion::completed,
                      grpc::StatusCode::OK,
                      runtime::PtServiceOutcome::indeterminate,
                      true,
                      101U,
                      202U,
                      303us});
    observer.observe({mpmc::runtime_grpc::PtGrpcRpcMethod::solve_pt_flash,
                      mpmc::runtime_grpc::PtGrpcCompletion::completed,
                      grpc::StatusCode::OK,
                      runtime::PtServiceOutcome::error,
                      true,
                      11U,
                      22U,
                      33us});
    const auto text = output.str();
    require(text.find("\"service_outcome\":\"indeterminate\"") !=
                    std::string::npos &&
                text.find("\"service_outcome\":\"service_error\"") !=
                    std::string::npos,
            "observer conflated scientific indeterminate and service error");
    require(text.find("pressure") == std::string::npos &&
                text.find("temperature") == std::string::npos &&
                text.find("composition") == std::string::npos &&
                text.find("certificate") == std::string::npos,
            "observer exposed prohibited scientific or credential fields");
}

void mtls_host(int argc, char** argv) {
    require(argc == 10, "mTLS host test requires eight certificate paths");
    const auto edge_certificate = read_file(argv[5]);
    const auto edge_private_key = read_file(argv[6]);
    const auto backend_ca = read_file(argv[7]);
    const auto untrusted_certificate = read_file(argv[8]);
    const auto untrusted_private_key = read_file(argv[9]);

    auto backend = std::make_shared<pt_grpc_test::Backend>();
    std::vector<process::OwnedConfiguredPtBackend> configured{{
        "host.fixture", backend}};
    std::ostringstream log;
    auto observer = std::make_shared<process::PtJsonLineObserver>(log);
    process::PtCompositionRoot root(std::move(configured), {}, {}, observer);

    process::PtProcessHostOptions options;
    options.listen_address = "127.0.0.1:0";
    options.tls = process::load_pt_process_tls_identity(
        {argv[2], argv[3], argv[4]});
    options.shutdown_grace = 2s;
    process::PtProcessHost host(root, std::move(options));
    host.start();
    require(host.running() && host.selected_port() > 0,
            "secure process host did not start");

    auto unauthenticated = make_stub(
        host.selected_port(), client_credentials(backend_ca));
    wire::DiscoverPtCapabilitiesResponse rejected_response;
    const auto rejected = discover(*unauthenticated, rejected_response);
    require(!rejected.ok() && rejected_response.backends_size() == 0,
            "host accepted a client without a trusted certificate");

    auto untrusted = make_stub(
        host.selected_port(), client_credentials(
                                  backend_ca, untrusted_certificate,
                                  untrusted_private_key));
    wire::DiscoverPtCapabilitiesResponse untrusted_response;
    const auto untrusted_status = discover(*untrusted, untrusted_response);
    require(!untrusted_status.ok() && untrusted_response.backends_size() == 0,
            "host accepted a client signed by an untrusted authority");

    auto authenticated = make_stub(
        host.selected_port(),
        client_credentials(backend_ca, edge_certificate, edge_private_key));
    wire::DiscoverPtCapabilitiesResponse accepted_response;
    const auto accepted = discover(*authenticated, accepted_response);
    require(accepted.ok() && accepted_response.backends_size() == 1 &&
                accepted_response.backends(0).configured_backend_id() ==
                    "host.fixture" &&
                accepted_response.backends(0)
                        .component_inventory()
                        .components_size() == 2,
            "authenticated discovery did not traverse the composition root");

    host.shutdown();
    host.wait();
    require(!host.running(), "host remained serving after shutdown");
    require(log.str().find("\"method\":\"discover_pt_capabilities\"") !=
                std::string::npos,
            "authenticated RPC was not observed");

    process::PtProcessTlsFiles undersized{argv[2], argv[3], argv[4]};
    undersized.max_file_bytes = 1U;
    expect_error<std::runtime_error>([&] {
        (void)process::load_pt_process_tls_identity(undersized);
    });

    process::PtProcessHostOptions insecure;
    expect_error<std::invalid_argument>([&] {
        process::PtProcessHost rejected_host(root, std::move(insecure));
        (void)rejected_host;
    });

    std::vector<process::OwnedConfiguredPtBackend> none;
    process::PtCompositionRoot empty_root(std::move(none));
    process::PtProcessHostOptions empty_options;
    empty_options.tls = process::load_pt_process_tls_identity(
        {argv[2], argv[3], argv[4]});
    expect_error<std::invalid_argument>([&] {
        process::PtProcessHost empty_host(empty_root, std::move(empty_options));
        (void)empty_host;
    });
}

} // namespace

int main(int argc, char** argv) {
    try {
        if (argc < 2) {
            std::cerr << "missing PT process test case\n";
            return 2;
        }
        const std::string_view test_case = argv[1];
        if (test_case == "composition_root") {
            composition_root();
        } else if (test_case == "observer") {
            observer();
        } else if (test_case == "headers") {
            require(pt_process_headers(), "public headers failed containment test");
        } else if (test_case == "mtls_host") {
            mtls_host(argc, argv);
        } else {
            std::cerr << "unknown PT process test case\n";
            return 2;
        }
        return 0;
    } catch (const std::exception& error) {
        std::cerr << error.what() << '\n';
        return 1;
    }
}
