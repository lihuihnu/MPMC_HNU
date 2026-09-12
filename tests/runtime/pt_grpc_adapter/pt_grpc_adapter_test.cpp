#include <mpmc/runtime_grpc/pt_grpc_adapter.hpp>

#include "test_backend.hpp"

#include <grpcpp/grpcpp.h>
#include <mpmc/runtime/v1/pt_service.grpc.pb.h>

#include <array>
#include <chrono>
#include <cstddef>
#include <iostream>
#include <memory>
#include <source_location>
#include <stdexcept>
#include <string>
#include <string_view>
#include <thread>
#include <utility>

bool pt_grpc_adapter_headers();

namespace {

namespace rt = ::mpmc::runtime;
namespace adapter = ::mpmc::runtime_grpc;
namespace wire = ::mpmc::runtime::v1;
using namespace std::chrono_literals;

void require(bool condition, std::string_view message,
             std::source_location where = std::source_location::current()) {
    if (!condition) {
        throw std::runtime_error(std::string(where.file_name()) + ":" +
                                 std::to_string(where.line()) + ": " +
                                 std::string(message));
    }
}

class ServerHarness {
public:
    ServerHarness(rt::PtService& service,
                  adapter::PtGrpcAdapterLimits limits = {})
        : adapter_(service, limits) {
        grpc::ServerBuilder builder;
        int selected_port = 0;
        builder.AddListeningPort("127.0.0.1:0",
                                 grpc::InsecureServerCredentials(),
                                 &selected_port);
        adapter::configure_pt_grpc_server(builder, adapter_);
        server_ = builder.BuildAndStart();
        require(server_ != nullptr && selected_port > 0,
                "failed to start local PT gRPC server");
        channel_ = grpc::CreateChannel(
            "127.0.0.1:" + std::to_string(selected_port),
            grpc::InsecureChannelCredentials());
        stub_ = wire::PtFlashService::NewStub(channel_);
    }

    ~ServerHarness() {
        if (server_) {
            server_->Shutdown(std::chrono::system_clock::now() + 2s);
            server_->Wait();
        }
    }

    [[nodiscard]] wire::PtFlashService::Stub& stub() { return *stub_; }

    [[nodiscard]] adapter::PtGrpcServiceAdapter& adapter() noexcept {
        return adapter_;
    }

private:
    adapter::PtGrpcServiceAdapter adapter_;
    std::unique_ptr<grpc::Server> server_;
    std::shared_ptr<grpc::Channel> channel_;
    std::unique_ptr<wire::PtFlashService::Stub> stub_;
};

void set_deadline(grpc::ClientContext& context,
                  std::chrono::milliseconds timeout = 2s) {
    context.set_deadline(std::chrono::system_clock::now() + timeout);
}

wire::SolvePtFlashRequest request(std::string configured_backend_id) {
    wire::SolvePtFlashRequest result;
    result.set_wire_contract("mpmc.runtime.v1/PT-flash-service/v1");
    result.set_service_request_convention("PT/service-request/v1");
    result.set_configured_backend_id(std::move(configured_backend_id));
    result.set_pressure_pa(5.0e6);
    result.set_temperature_k(325.0);
    auto* water = result.add_feed();
    water->set_component_id("water");
    water->set_mole_fraction(0.6);
    auto* methane = result.add_feed();
    methane->set_component_id("methane");
    methane->set_mole_fraction(0.4);
    return result;
}

template <typename Predicate>
bool wait_until(Predicate&& predicate,
                std::chrono::milliseconds timeout = 2s) {
    const auto deadline = std::chrono::steady_clock::now() + timeout;
    while (!predicate()) {
        if (std::chrono::steady_clock::now() >= deadline) { return false; }
        std::this_thread::yield();
    }
    return true;
}

void mapping_and_outcomes() {
    pt_grpc_test::Backend accepted;
    pt_grpc_test::Backend indeterminate(
        pt_grpc_test::BackendOutcome::indeterminate);
    const std::array<rt::PtServiceBackendRegistration, 2> registrations{{
        {"golden.accepted", &accepted},
        {"golden.indeterminate", &indeterminate}}};
    rt::PtService service(registrations);
    ServerHarness server(service);

    grpc::ClientContext discovery_context;
    set_deadline(discovery_context);
    wire::DiscoverPtCapabilitiesRequest discovery_request;
    wire::DiscoverPtCapabilitiesResponse discovery;
    const auto discovery_status = server.stub().DiscoverPtCapabilities(
        &discovery_context, discovery_request, &discovery);
    require(discovery_status.ok() && discovery.has_wire_contract() &&
                discovery.wire_contract() ==
                    "mpmc.runtime.v1/PT-flash-service/v1" &&
                discovery.backends_size() == 2 &&
                discovery.backends(0).configured_backend_id() ==
                    "golden.accepted" &&
                discovery.backends(0).component_inventory().components_size() ==
                    2 &&
                discovery.backends(0)
                        .component_inventory()
                        .components(0)
                        .component_id() == "methane",
            "capability discovery did not preserve the runtime snapshot");

    grpc::ClientContext accepted_context;
    set_deadline(accepted_context);
    const auto accepted_request = request("golden.accepted");
    wire::SolvePtFlashResponse accepted_response;
    const auto accepted_status = server.stub().SolvePtFlash(
        &accepted_context, accepted_request, &accepted_response);
    require(accepted_status.ok() && accepted_response.has_result() &&
                !accepted_response.has_error() &&
                accepted_response.result().outcome() ==
                    wire::PT_COMPUTATION_OUTCOME_ACCEPTED &&
                accepted_response.result().phases_size() == 2 &&
                accepted_response.result().feed(0).component_id() == "methane" &&
                accepted_response.result().feed(0).mole_fraction() == 0.4 &&
                accepted_response.result().phases(0).provider_branch() == 7U &&
                accepted_response.result()
                        .provenance()
                        .backend()
                        .configured_backend_id() == "golden.accepted" &&
                accepted.call_count.load() == 1U,
            "accepted response/provenance mapping changed or solve was repeated");

    grpc::ClientContext indeterminate_context;
    set_deadline(indeterminate_context);
    const auto indeterminate_request = request("golden.indeterminate");
    wire::SolvePtFlashResponse indeterminate_response;
    const auto indeterminate_status = server.stub().SolvePtFlash(
        &indeterminate_context, indeterminate_request,
        &indeterminate_response);
    require(indeterminate_status.ok() &&
                indeterminate_response.has_result() &&
                !indeterminate_response.has_error() &&
                indeterminate_response.result().outcome() ==
                    wire::PT_COMPUTATION_OUTCOME_INDETERMINATE &&
                indeterminate_response.result().phases_size() == 0 &&
                indeterminate.call_count.load() == 1U,
            "indeterminate was not preserved as a phase-free result arm");

    grpc::ClientContext error_context;
    set_deadline(error_context);
    auto invalid_request = request("golden.accepted");
    invalid_request.set_pressure_pa(-1.0);
    wire::SolvePtFlashResponse error_response;
    const auto error_status = server.stub().SolvePtFlash(
        &error_context, invalid_request, &error_response);
    require(error_status.ok() && error_response.has_error() &&
                !error_response.has_result() &&
                error_response.error().code() ==
                    wire::PT_SERVICE_ERROR_CODE_INVALID_PRESSURE &&
                accepted.call_count.load() == 1U,
            "PtService error was not kept separate from scientific results");
}

void wire_and_size_limits() {
    pt_grpc_test::Backend backend;
    const std::array<rt::PtServiceBackendRegistration, 1> registrations{{
        {"golden.accepted", &backend}}};
    rt::PtService service(registrations);
    auto limits = adapter::PtGrpcAdapterLimits{};
    limits.max_serialized_request_bytes = 512U;
    ServerHarness server(service, limits);

    grpc::ClientContext missing_deadline_context;
    const auto valid_request = request("golden.accepted");
    wire::SolvePtFlashResponse missing_deadline_response;
    const auto missing_deadline = server.stub().SolvePtFlash(
        &missing_deadline_context, valid_request,
        &missing_deadline_response);
    require(missing_deadline.error_code() == grpc::StatusCode::INVALID_ARGUMENT &&
                backend.call_count.load() == 0U,
            "finite client deadline policy was not enforced before solve");

    grpc::ClientContext discovery_deadline_context;
    set_deadline(discovery_deadline_context, 11s);
    wire::DiscoverPtCapabilitiesRequest discovery_request;
    wire::DiscoverPtCapabilitiesResponse discovery_response;
    const auto discovery_deadline = server.stub().DiscoverPtCapabilities(
        &discovery_deadline_context, discovery_request, &discovery_response);
    require(discovery_deadline.error_code() ==
                grpc::StatusCode::INVALID_ARGUMENT &&
                discovery_response.backends_size() == 0,
            "over-limit discovery deadline was accepted");

    grpc::ClientContext long_deadline_context;
    set_deadline(long_deadline_context, 121s);
    wire::SolvePtFlashResponse long_deadline_response;
    const auto long_deadline = server.stub().SolvePtFlash(
        &long_deadline_context, valid_request, &long_deadline_response);
    require(long_deadline.error_code() == grpc::StatusCode::INVALID_ARGUMENT &&
                backend.call_count.load() == 0U,
            "over-limit client deadline reached PtService");

    grpc::ClientContext malformed_context;
    set_deadline(malformed_context);
    auto malformed_request = valid_request;
    malformed_request.clear_wire_contract();
    wire::SolvePtFlashResponse malformed_response;
    const auto malformed = server.stub().SolvePtFlash(
        &malformed_context, malformed_request, &malformed_response);
    require(malformed.error_code() == grpc::StatusCode::INVALID_ARGUMENT &&
                backend.call_count.load() == 0U,
            "missing wire version reached PtService");

    grpc::ClientContext oversized_context;
    set_deadline(oversized_context);
    auto oversized_request = valid_request;
    oversized_request.mutable_feed(0)->set_component_id(
        std::string(2'048U, 'x'));
    wire::SolvePtFlashResponse oversized_response;
    const auto oversized = server.stub().SolvePtFlash(
        &oversized_context, oversized_request, &oversized_response);
    require(oversized.error_code() == grpc::StatusCode::RESOURCE_EXHAUSTED &&
                backend.call_count.load() == 0U,
            "transport message-size limit did not reject request before solve");

    auto response_limits = adapter::PtGrpcAdapterLimits{};
    response_limits.max_serialized_response_bytes = 128U;
    ServerHarness response_limited_server(service, response_limits);
    grpc::ClientContext response_limited_context;
    set_deadline(response_limited_context);
    wire::SolvePtFlashResponse response_limited_response;
    const auto response_limited = response_limited_server.stub().SolvePtFlash(
        &response_limited_context, valid_request, &response_limited_response);
    require(response_limited.error_code() ==
                grpc::StatusCode::RESOURCE_EXHAUSTED &&
                response_limited_response.payload_case() ==
                    wire::SolvePtFlashResponse::PAYLOAD_NOT_SET &&
                backend.call_count.load() == 1U,
            "serialized response limit did not suppress an oversized payload");

    auto invalid_limits = adapter::PtGrpcAdapterLimits{};
    invalid_limits.max_concurrent_solves = 0U;
    bool rejected = false;
    try {
        adapter::PtGrpcServiceAdapter invalid_adapter(service, invalid_limits);
        (void)invalid_adapter;
    } catch (const std::invalid_argument&) {
        rejected = true;
    }
    require(rejected, "invalid process adapter limits were accepted");

    invalid_limits = adapter::PtGrpcAdapterLimits{};
    invalid_limits.grpc_max_threads = 0U;
    rejected = false;
    try {
        adapter::PtGrpcServiceAdapter invalid_adapter(service, invalid_limits);
        (void)invalid_adapter;
    } catch (const std::invalid_argument&) {
        rejected = true;
    }
    require(rejected, "invalid gRPC thread quota was accepted");
}

void concurrency_gate() {
    pt_grpc_test::Backend backend;
    backend.block();
    const std::array<rt::PtServiceBackendRegistration, 1> registrations{{
        {"golden.accepted", &backend}}};
    rt::PtService service(registrations);
    auto limits = adapter::PtGrpcAdapterLimits{};
    limits.max_concurrent_solves = 1U;
    ServerHarness server(service, limits);

    grpc::Status first_status;
    wire::SolvePtFlashResponse first_response;
    std::thread first([&] {
        grpc::ClientContext context;
        set_deadline(context, 5s);
        const auto value = request("golden.accepted");
        first_status = server.stub().SolvePtFlash(&context, value,
                                                  &first_response);
    });
    require(backend.wait_until_started(), "first backend solve did not start");
    require(server.adapter().in_flight_solves() == 1U,
            "in-flight solve accounting changed");

    grpc::ClientContext second_context;
    set_deadline(second_context);
    const auto second_request = request("golden.accepted");
    wire::SolvePtFlashResponse second_response;
    const auto second_status = server.stub().SolvePtFlash(
        &second_context, second_request, &second_response);
    require(second_status.error_code() ==
                grpc::StatusCode::RESOURCE_EXHAUSTED &&
                backend.call_count.load() == 1U,
            "concurrency gate started a second backend solve");

    backend.release();
    first.join();
    require(first_status.ok() && first_response.has_result() &&
                server.adapter().in_flight_solves() == 0U,
            "admitted solve did not finish or release its permit");
}

void deadline_and_cancel() {
    pt_grpc_test::Backend backend;
    const std::array<rt::PtServiceBackendRegistration, 1> registrations{{
        {"golden.accepted", &backend}}};
    rt::PtService service(registrations);
    ServerHarness server(service);

    backend.block();
    grpc::Status deadline_status;
    std::thread deadline_call([&] {
        grpc::ClientContext context;
        set_deadline(context, 75ms);
        const auto value = request("golden.accepted");
        wire::SolvePtFlashResponse response;
        deadline_status = server.stub().SolvePtFlash(&context, value, &response);
    });
    require(backend.wait_until_started(), "deadline test solve did not start");
    deadline_call.join();
    require(deadline_status.error_code() ==
                grpc::StatusCode::DEADLINE_EXCEEDED &&
                backend.call_count.load() == 1U &&
                server.adapter().in_flight_solves() == 1U,
            "elapsed deadline was published as a scientific/service result");
    backend.release();
    require(backend.wait_until_completed(1U) &&
                wait_until([&] {
                    return server.adapter().in_flight_solves() == 0U;
                }),
            "deadline-expired backend work did not drain its permit");

    backend.block();
    grpc::ClientContext cancel_context;
    set_deadline(cancel_context, 2s);
    grpc::Status cancel_status;
    std::thread cancel_call([&] {
        const auto value = request("golden.accepted");
        wire::SolvePtFlashResponse response;
        cancel_status = server.stub().SolvePtFlash(&cancel_context, value,
                                                   &response);
    });
    require(backend.wait_until_started(), "cancel test solve did not start");
    cancel_context.TryCancel();
    cancel_call.join();
    require(cancel_status.error_code() == grpc::StatusCode::CANCELLED &&
                backend.call_count.load() == 2U &&
                server.adapter().in_flight_solves() == 1U,
            "explicit cancellation was published as a scientific/service result");
    backend.release();
    require(backend.wait_until_completed(2U) &&
                wait_until([&] {
                    return server.adapter().in_flight_solves() == 0U;
                }),
            "cancelled backend work did not drain its permit");
}

void headers() {
    require(pt_grpc_adapter_headers(),
            "PT gRPC adapter public header is not self-contained");
}

} // namespace

int main(int argc, char** argv) {
    if (argc != 2) {
        std::cerr << "expected one PT gRPC adapter test case\n";
        return 2;
    }
    try {
        const std::string_view test_case(argv[1]);
        if (test_case == "mapping_and_outcomes") {
            mapping_and_outcomes();
        } else if (test_case == "wire_and_size_limits") {
            wire_and_size_limits();
        } else if (test_case == "concurrency_gate") {
            concurrency_gate();
        } else if (test_case == "deadline_and_cancel") {
            deadline_and_cancel();
        } else if (test_case == "headers") {
            headers();
        } else {
            std::cerr << "unknown PT gRPC adapter test case\n";
            return 2;
        }
        return 0;
    } catch (const std::exception& error) {
        std::cerr << error.what() << '\n';
        return 1;
    }
}
