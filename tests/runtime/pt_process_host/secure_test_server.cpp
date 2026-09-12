#include <mpmc/pt_process/composition_root.hpp>
#include <mpmc/pt_process/json_line_observer.hpp>
#include <mpmc/pt_process/process_host.hpp>

#include "../pt_grpc_adapter/test_backend.hpp"

#include <iostream>
#include <memory>
#include <string>
#include <utility>
#include <vector>

int main(int argc, char** argv) {
    if (argc != 5) {
        std::cerr << "usage: secure_test_server ADDRESS CERT KEY CLIENT_CA\n";
        return 2;
    }
    try {
        auto backend = std::make_shared<pt_grpc_test::Backend>();
        std::vector<mpmc::pt_process::OwnedConfiguredPtBackend> configured{{
            "golden.accepted", backend}};
        auto observer =
            std::make_shared<mpmc::pt_process::PtJsonLineObserver>(std::clog);
        mpmc::pt_process::PtCompositionRoot root(
            std::move(configured), {}, {}, std::move(observer));

        mpmc::pt_process::PtProcessHostOptions options;
        options.listen_address = argv[1];
        options.tls = mpmc::pt_process::load_pt_process_tls_identity(
            {argv[2], argv[3], argv[4]});
        mpmc::pt_process::PtProcessHost host(root, std::move(options));
        host.start();
        std::cout << "READY " << host.selected_port() << '\n' << std::flush;
        host.wait();
        return 0;
    } catch (const std::exception& error) {
        std::cerr << error.what() << '\n';
        return 1;
    }
}
