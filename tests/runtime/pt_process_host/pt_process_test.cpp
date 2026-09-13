#include <mpmc/pt_process/composition_root.hpp>
#include <mpmc/pt_process/configured_backends.hpp>
#include <mpmc/pt_process/json_line_observer.hpp>
#include <mpmc/pt_process/parameter_snapshot_supplier.hpp>
#include <mpmc/pt_process/process_host.hpp>

#include "../pt_grpc_adapter/test_backend.hpp"
#include "../../flash/cpa_max3/test_support.hpp"
#include "../../flash/pr76_three_phase/synthetic_fixture.hpp"
#include "test_support.hpp"

#include <mpmc/runtime/v1/pt_service.grpc.pb.h>

#include <chrono>
#include <cmath>
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

void repository_snapshot_bundle() {
    const auto parameters =
        process::load_repository_curated_pt_parameters_v1();
    const auto& pr76_parameters = parameters.pr76;
    const auto pr76_tc = pr76_parameters.critical_temperatures_k();
    const auto pr76_pc = pr76_parameters.critical_pressures_pa();
    const auto pr76_omega = pr76_parameters.acentric_factors();
    require(pr76_tc.size() == 3U && pr76_tc[0] == 190.555 &&
                pr76_tc[1] == 305.4 && pr76_tc[2] == 369.825 &&
                pr76_pc.size() == 3U && pr76_pc[0] == 4.595e6 &&
                pr76_pc[1] == 4.88e6 && pr76_pc[2] == 4.248e6 &&
                pr76_omega.size() == 3U && pr76_omega[0] == 0.0 &&
                pr76_omega[1] == 0.099 &&
                pr76_omega[2] == 0.15308 &&
                pr76_parameters.binary_records().size() == 3U,
            "PR76 repository parameter payload changed");
    for (const auto& component : pr76_parameters.components().items()) {
        require(component.definition.kind ==
                    mpmc::thermodynamics::SourceKind::literature,
                "PR76 repository component lost literature provenance");
    }
    for (const auto& record : pr76_parameters.pure_records()) {
        require(record.critical_temperature.has_value() &&
                    record.critical_pressure.has_value() &&
                    record.acentric_factor.has_value() &&
                    record.critical_temperature->source.kind ==
                        mpmc::thermodynamics::SourceKind::literature &&
                    record.critical_pressure->source.kind ==
                        mpmc::thermodynamics::SourceKind::literature &&
                    record.acentric_factor->source.kind ==
                        mpmc::thermodynamics::SourceKind::literature,
                "PR76 repository pure record lost literature provenance");
    }
    for (const auto& pair : pr76_parameters.binary_records()) {
        require(pair.kij.has_value() && pair.kij->value == 0.0 &&
                    pair.kij->source.kind ==
                        mpmc::thermodynamics::SourceKind::literature,
                "PR76 repository binary record is implicit or non-literature");
    }

    const auto& sw92_parameters = parameters.sw92;
    const auto sw92_tc = sw92_parameters.critical_temperatures_k();
    const auto sw92_pc = sw92_parameters.critical_pressures_pa();
    const auto sw92_omega = sw92_parameters.acentric_factors();
    require(sw92_tc.size() == 2U && sw92_tc[0] == 304.2 &&
                sw92_tc[1] == 647.3 && sw92_pc.size() == 2U &&
                sw92_pc[0] == 7.38e6 && sw92_pc[1] == 22.12e6 &&
                sw92_omega.size() == 2U && sw92_omega[0] == 0.2273 &&
                sw92_omega[1] == 0.3434 &&
                sw92_parameters.water_nonaqueous_constant_kij(0U) == 0.1896 &&
                sw92_parameters.applicability()
                    .nacl_molality_mol_per_kg_water.has_value() &&
                sw92_parameters.applicability()
                        .nacl_molality_mol_per_kg_water->lower == 0.0 &&
                sw92_parameters.applicability()
                        .nacl_molality_mol_per_kg_water->upper == 0.0,
            "SW92 repository freshwater parameter payload changed");
    for (const auto& record : sw92_parameters.pure_records()) {
        require(record.critical_temperature.has_value() &&
                    record.critical_pressure.has_value() &&
                    record.acentric_factor.has_value() &&
                    record.critical_temperature->source.kind ==
                        mpmc::thermodynamics::SourceKind::literature &&
                    record.critical_pressure->source.kind ==
                        mpmc::thermodynamics::SourceKind::literature &&
                    record.acentric_factor->source.kind ==
                        mpmc::thermodynamics::SourceKind::literature,
                "SW92 repository pure record lost literature provenance");
    }
    require(sw92_parameters.water_binary_records().size() == 1U &&
                sw92_parameters.water_binary_records()[0]
                        .nonaqueous_kij.has_value() &&
                sw92_parameters.water_binary_records()[0]
                        .nonaqueous_kij->source.kind ==
                    mpmc::thermodynamics::SourceKind::literature,
            "SW92 repository pair record lost literature provenance");

    const auto& cpa_parameters = parameters.cpa;
    require(cpa_parameters.size() == 2U &&
                cpa_parameters.pure(0U).critical_temperature_k == 512.64 &&
                cpa_parameters.pure(0U).a0_pa_m6_per_mol2 == 0.40531 &&
                cpa_parameters.pure(0U).b_m3_per_mol == 3.0978e-5 &&
                cpa_parameters.pure(0U).c1_dimensionless == 0.43102 &&
                cpa_parameters.pure(1U).critical_temperature_k == 647.29 &&
                cpa_parameters.pure(1U).a0_pa_m6_per_mol2 == 0.12277 &&
                cpa_parameters.pure(1U).b_m3_per_mol == 1.4515e-5 &&
                cpa_parameters.pure(1U).c1_dimensionless == 0.67359 &&
                cpa_parameters.kij(0U, 1U) == -0.055 &&
                cpa_parameters.association_records().size() == 4U,
            "CPA repository parameter payload changed");
    for (const auto& record : cpa_parameters.pure_records()) {
        require(
            record.critical_temperature_k.source.kind ==
                    mpmc::thermodynamics::SourceKind::literature &&
                record.a0_pa_m6_per_mol2.source.kind ==
                    mpmc::thermodynamics::SourceKind::literature &&
                record.b_m3_per_mol.source.kind ==
                    mpmc::thermodynamics::SourceKind::literature &&
                record.c1_dimensionless.source.kind ==
                    mpmc::thermodynamics::SourceKind::literature,
            "CPA repository pure record lost literature provenance");
    }
    require(cpa_parameters.binary_records().size() == 1U &&
                cpa_parameters.binary_records()[0]
                        .kij_dimensionless.source.kind ==
                    mpmc::thermodynamics::SourceKind::literature,
            "CPA repository binary record lost literature provenance");
    for (const auto& record : cpa_parameters.association_records()) {
        require(
            record.epsilon_j_per_mol.source.kind ==
                    mpmc::thermodynamics::SourceKind::literature &&
                record.beta_dimensionless.source.kind ==
                    mpmc::thermodynamics::SourceKind::literature,
            "CPA repository association record lost literature provenance");
    }

    auto bundle =
        process::load_repository_curated_pt_parameter_snapshots_v1();
    require(bundle.structurally_valid(),
            "repository snapshot bundle is structurally invalid");
    require(bundle.convention ==
                process::pt_parameter_snapshot_bundle_convention &&
                bundle.bundle_id ==
                    process::repository_curated_pt_bundle_id &&
                bundle.revision ==
                    process::repository_curated_pt_bundle_revision &&
                bundle.snapshots.size() == 3U &&
                bundle.backends.size() == 3U,
            "repository snapshot bundle identity changed");

    const auto& pr76 = bundle.snapshots[0];
    const auto& sw92 = bundle.snapshots[1];
    const auto& cpa = bundle.snapshots[2];
    require(pr76.configured_backend_id ==
                "pr76.methane-ethane-propane.literature-r1" &&
                pr76.dataset_id ==
                    "DeitersBell-aic16730-PengRobinson1976-ternary" &&
                pr76.component_ids ==
                    std::vector<std::string>(
                        {"methane", "ethane", "propane"}) &&
                pr76.sources.size() == 2U,
            "PR76 repository snapshot identity or provenance changed");
    require(sw92.configured_backend_id ==
                "sw92.carbon-dioxide-water.freshwater.literature-r1" &&
                sw92.component_ids ==
                    std::vector<std::string>(
                        {"carbon-dioxide", "water"}) &&
                sw92.sources.size() == 1U,
            "SW92 repository snapshot identity or provenance changed");
    require(cpa.configured_backend_id ==
                "cpa.methanol-water-333.15k.cr1.literature-r1" &&
                cpa.dataset_id ==
                    "literature-cpa-water-methanol-cr1-333.15K" &&
                cpa.component_ids ==
                    std::vector<std::string>({"METHANOL", "WATER"}) &&
                cpa.sources.size() == 3U,
            "CPA repository snapshot identity or provenance changed");
    for (const auto& snapshot : bundle.snapshots) {
        require(snapshot.structurally_valid() &&
                    snapshot.snapshot_location.starts_with("builtin://"),
                "snapshot lost a bounded applicability statement");
        for (const auto& source : snapshot.sources) {
            require(source.structurally_valid(),
                    "snapshot lost a complete public source descriptor");
        }
    }

    std::ostringstream manifest;
    process::write_pt_parameter_snapshot_manifest_json(manifest, bundle);
    const auto serialized = manifest.str();
    require(serialized.find(
                "\"bundle_id\":\"MPMC/PT/repository-curated-literature-snapshots/v1\"") !=
                std::string::npos &&
                serialized.find("\"configured_backend_id\":\"pr76.") !=
                    std::string::npos &&
                serialized.find("\"configured_backend_id\":\"sw92.") !=
                    std::string::npos &&
                serialized.find("\"configured_backend_id\":\"cpa.") !=
                    std::string::npos &&
                serialized.find("\"kij\"") == std::string::npos &&
                serialized.find("\"a0\"") == std::string::npos,
            "snapshot manifest omitted identity or exposed model parameters");

    auto tampered = bundle;
    tampered.snapshots[1].revision = "tampered";
    require(!tampered.structurally_valid(),
            "snapshot/backend provenance mismatch was accepted");

    process::PtCompositionRoot root(std::move(bundle.backends));
    const auto capabilities = root.service().discover_capabilities();
    require(capabilities.size() == 3U,
            "repository snapshot supplier did not register three backends");
    for (const auto& descriptor : capabilities) {
        require(descriptor.structurally_valid() &&
                    descriptor.capability.supported_phase_counts ==
                        std::vector<std::size_t>({1U, 2U, 3U}),
                "repository snapshot changed the frozen backend contract");
    }
    const auto* sw92_capability = root.service().find_capability(
        "sw92.carbon-dioxide-water.freshwater.literature-r1");
    require(sw92_capability != nullptr &&
                sw92_capability->capability.scalar_settings.size() == 1U &&
                sw92_capability->capability.scalar_settings[0].value == 0.0,
            "SW92 repository snapshot is not explicitly fixed to fresh water");
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
        } else if (test_case == "repository_snapshot_bundle") {
            repository_snapshot_bundle();
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
