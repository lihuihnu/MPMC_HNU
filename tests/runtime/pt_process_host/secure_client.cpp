#include <grpcpp/create_channel.h>
#include <grpcpp/security/credentials.h>
#include <grpcpp/support/channel_arguments.h>
#include <mpmc/runtime/v1/pt_service.grpc.pb.h>

#include <chrono>
#include <fstream>
#include <iostream>
#include <iterator>
#include <memory>
#include <stdexcept>
#include <string>

namespace {

std::string read_file(const char* path) {
    std::ifstream input(path, std::ios::binary);
    if (!input) { throw std::runtime_error("cannot read client TLS file"); }
    return {std::istreambuf_iterator<char>(input),
            std::istreambuf_iterator<char>()};
}

} // namespace

int main(int argc, char** argv) {
    if (argc != 6) {
        std::cerr << "usage: secure_client TARGET SERVER_NAME CA CERT KEY\n";
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
        auto stub = mpmc::runtime::v1::PtFlashService::NewStub(channel);

        grpc::ClientContext context;
        context.set_deadline(std::chrono::system_clock::now() +
                             std::chrono::seconds(10));
        mpmc::runtime::v1::DiscoverPtCapabilitiesRequest request;
        mpmc::runtime::v1::DiscoverPtCapabilitiesResponse response;
        const auto status =
            stub->DiscoverPtCapabilities(&context, request, &response);
        if (!status.ok() || response.backends_size() != 1 ||
            response.backends(0).configured_backend_id() != "golden.accepted" ||
            response.backends(0).component_inventory().components_size() != 2) {
            std::cerr << "secure discovery failed: " << status.error_message()
                      << '\n';
            return 1;
        }
        std::cout << "DISCOVERED "
                  << response.backends(0).configured_backend_id() << '\n';
        return 0;
    } catch (const std::exception& error) {
        std::cerr << error.what() << '\n';
        return 1;
    }
}
