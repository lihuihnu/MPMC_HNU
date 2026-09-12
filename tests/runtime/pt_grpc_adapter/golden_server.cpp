#include <mpmc/runtime_grpc/pt_grpc_adapter.hpp>

#include "test_backend.hpp"

#include <grpcpp/grpcpp.h>

#include <array>
#include <iostream>
#include <memory>
#include <string>

int main(int argc, char** argv) {
    const std::string address =
        argc > 1 ? argv[1] : "127.0.0.1:50051";

    pt_grpc_test::Backend accepted;
    pt_grpc_test::Backend indeterminate(
        pt_grpc_test::BackendOutcome::indeterminate);
    const std::array<mpmc::runtime::PtServiceBackendRegistration, 2>
        registrations{{
            {"golden.accepted", &accepted},
            {"golden.indeterminate", &indeterminate}}};
    mpmc::runtime::PtService service(registrations);
    mpmc::runtime_grpc::PtGrpcServiceAdapter adapter(service);

    grpc::ServerBuilder builder;
    int selected_port = 0;
    builder.AddListeningPort(address, grpc::InsecureServerCredentials(),
                             &selected_port);
    mpmc::runtime_grpc::configure_pt_grpc_server(builder, adapter);
    std::unique_ptr<grpc::Server> server = builder.BuildAndStart();
    if (!server || selected_port <= 0) {
        std::cerr << "failed to start PT gRPC golden server\n";
        return 1;
    }

    std::cout << "READY " << selected_port << '\n' << std::flush;
    server->Wait();
    return 0;
}
