#include <grpcpp/create_channel.h>
#include <grpcpp/security/credentials.h>
#include <grpcpp/support/channel_arguments.h>
#include <mpmc/runtime/v1/pt_service.grpc.pb.h>

#include <array>
#include <chrono>
#include <fstream>
#include <iostream>
#include <iterator>
#include <stdexcept>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

namespace {

namespace wire = mpmc::runtime::v1;

struct ExpectedBackend {
    std::string_view configured_id;
    std::vector<std::string_view> components;
};

std::string read_file(const char* path) {
    std::ifstream input(path, std::ios::binary);
    if (!input) { throw std::runtime_error("cannot read product probe TLS file"); }
    return {std::istreambuf_iterator<char>(input),
            std::istreambuf_iterator<char>()};
}

void require(bool condition, std::string_view message) {
    if (!condition) { throw std::runtime_error(std::string(message)); }
}

void verify_discovery(const wire::DiscoverPtCapabilitiesResponse& response) {
    const std::array<ExpectedBackend, 3> expected{{
        {"pr76.methane-ethane-propane.literature-r1",
         {"methane", "ethane", "propane"}},
        {"sw92.carbon-dioxide-water.freshwater.literature-r1",
         {"carbon-dioxide", "water"}},
        {"cpa.methanol-water-333.15k.cr1.literature-r1",
         {"METHANOL", "WATER"}},
    }};
    require(response.has_wire_contract() &&
                response.wire_contract() ==
                    "mpmc.runtime.v1/PT-flash-service/v1" &&
                response.backends_size() == static_cast<int>(expected.size()),
            "staged discovery envelope changed");

    for (std::size_t index = 0; index < expected.size(); ++index) {
        const auto& actual = response.backends(static_cast<int>(index));
        const auto& wanted = expected[index];
        require(actual.configured_backend_id() == wanted.configured_id &&
                    actual.has_capability() &&
                    actual.has_component_inventory(),
                "staged backend identity or inventory is missing");
        require(actual.capability().supported_phase_counts_size() == 3 &&
                    actual.capability().supported_phase_counts(0) == 1U &&
                    actual.capability().supported_phase_counts(1) == 2U &&
                    actual.capability().supported_phase_counts(2) == 3U &&
                    actual.capability().performs_final_phase_set_review() &&
                    !actual.capability().global_stability_proven(),
                "staged backend capability contract changed");
        require(actual.component_inventory().components_size() ==
                    static_cast<int>(wanted.components.size()),
                "staged component count changed");
        for (std::size_t component = 0; component < wanted.components.size();
             ++component) {
            const auto& item = actual.component_inventory().components(
                static_cast<int>(component));
            require(item.component_id() == wanted.components[component] &&
                        item.feed_index() == component,
                    "staged ordered component inventory changed");
        }
    }
}

void verify_service_error(wire::PtFlashService::Stub& stub) {
    wire::SolvePtFlashRequest request;
    request.set_wire_contract("mpmc.runtime.v1/PT-flash-service/v1");
    request.set_service_request_convention("PT/service-request/v1");
    request.set_configured_backend_id(
        "pr76.methane-ethane-propane.literature-r1");
    request.set_pressure_pa(-1.0);
    request.set_temperature_k(300.0);
    for (const auto& [id, fraction] :
         std::array<std::pair<std::string_view, double>, 3>{{
             {"methane", 0.3}, {"ethane", 0.3}, {"propane", 0.4}}}) {
        auto* entry = request.add_feed();
        entry->set_component_id(id);
        entry->set_mole_fraction(fraction);
    }

    grpc::ClientContext context;
    context.set_deadline(std::chrono::system_clock::now() +
                         std::chrono::seconds(10));
    wire::SolvePtFlashResponse response;
    const auto status = stub.SolvePtFlash(&context, request, &response);
    require(status.ok() && response.has_error() && !response.has_result() &&
                response.error().code() ==
                    wire::PT_SERVICE_ERROR_CODE_INVALID_PRESSURE,
            "staged service error was not kept outside the result arm");
}

} // namespace

int main(int argc, char** argv) {
    if (argc != 6) {
        std::cerr << "usage: staged_probe TARGET SERVER_NAME CA CERT KEY\n";
        return 2;
    }
    try {
        grpc::SslCredentialsOptions tls;
        tls.pem_root_certs = read_file(argv[3]);
        tls.pem_cert_chain = read_file(argv[4]);
        tls.pem_private_key = read_file(argv[5]);
        grpc::ChannelArguments arguments;
        arguments.SetSslTargetNameOverride(argv[2]);
        const auto channel = grpc::CreateCustomChannel(
            argv[1], grpc::SslCredentials(tls), arguments);
        require(channel->WaitForConnected(
                    std::chrono::system_clock::now() +
                    std::chrono::seconds(10)),
                "staged product channel did not become ready");
        auto stub = wire::PtFlashService::NewStub(channel);

        grpc::ClientContext discovery_context;
        discovery_context.set_deadline(std::chrono::system_clock::now() +
                                       std::chrono::seconds(10));
        wire::DiscoverPtCapabilitiesRequest request;
        wire::DiscoverPtCapabilitiesResponse response;
        const auto status = stub->DiscoverPtCapabilities(
            &discovery_context, request, &response);
        require(status.ok(), "staged capability discovery RPC failed");
        verify_discovery(response);
        verify_service_error(*stub);
        std::cout << "STAGED_PRODUCT_OK\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << error.what() << '\n';
        return 1;
    }
}
