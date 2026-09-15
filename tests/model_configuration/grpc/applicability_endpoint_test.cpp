#include "test_support.hpp"

#include <google/rpc/status.pb.h>
#include <iostream>
#include <stdexcept>
#include <string_view>

namespace {
using namespace service_test;
using SC = grpc::StatusCode;

std::string detail_field(const grpc::Status& status) {
    google::rpc::Status envelope;
    wire::ModelServiceError detail;
    require(envelope.ParseFromString(status.error_details()) &&
            envelope.details_size() == 1 && envelope.details(0).UnpackTo(&detail),
            "missing structured applicability error detail");
    return detail.field();
}

void expect_solve_rejected(Server& server, const std::string& handle,
                           const fl::PtFlashRequest& state, std::string_view field) {
    wire::SolveModelResponse out;
    const auto status = solve(server, solve_request(handle, state), out);
    error(status, SC::INVALID_ARGUMENT, "request.rejected");
    require(detail_field(status) == field, "wrong wire applicability rejection field");
    require(!out.has_result(), "rejected applicability state published a result");
}

void run() {
    mc::Pr76ModelRegistry registry({}, synthetic);
    Server server(registry);
    const auto native = pr76_max3_test::model();

    auto d = definition(native);
    d.applicability.temperature_lower_k = 250.0;
    d.applicability.pressure_upper_pa = 1.0e6;
    d.applicability.pressure_upper_exclusive = true;
    const auto request_wire = request(d);
    require(request_wire.definition().applicability().has_temperature_lower_k() &&
            !request_wire.definition().applicability().has_temperature_upper_k() &&
            request_wire.definition().applicability().pressure_upper_exclusive(),
            "public endpoint presence/open flag not encoded");

    const auto made = create(server, request_wire);
    require(made.snapshot().definition().SerializeAsString() ==
                request_wire.definition().SerializeAsString(),
            "create snapshot lost applicability endpoint semantics");
    const auto stored = registry.describe(made.model_handle());
    require(stored.definition.applicability.temperature_lower_k == 250.0 &&
            !stored.definition.applicability.temperature_upper_k &&
            stored.definition.applicability.pressure_upper_pa == 1.0e6 &&
            stored.definition.applicability.pressure_upper_exclusive,
            "registry snapshot lost applicability endpoint semantics");

    expect_solve_rejected(server, made.model_handle(),
                          {1.0e6, 250.0, single.feed}, "pressure_pa");
    expect_solve_rejected(server, made.model_handle(),
                          {0.9e6, 249.0, single.feed}, "temperature_k");
    const fl::PtFlashRequest interior_unknown{0.9e6, 250.0, single.feed};
    compare(native_result(solve(server, made.model_handle(), interior_unknown).result()),
            direct(native, interior_unknown));

    wire::DescribeModelResponse described;
    ok(describe(server, made.model_handle(), described));
    require(described.snapshot().definition().SerializeAsString() ==
                request_wire.definition().SerializeAsString(),
            "describe snapshot changed endpoint presence/open flags");
    ok(release(server, made.model_handle()));

    // Complete interval: open equality is rejected, closed equality is accepted.
    auto complete = definition(native);
    complete.applicability.temperature_lower_k = 240.0;
    complete.applicability.temperature_upper_k = 260.0;
    complete.applicability.temperature_upper_exclusive = true;
    complete.applicability.pressure_lower_pa = 0.9e6;
    complete.applicability.pressure_lower_exclusive = true;
    complete.applicability.pressure_upper_pa = 1.1e6;
    const auto complete_made = create(server, request(complete));
    expect_solve_rejected(server, complete_made.model_handle(),
                          {1.0e6, 260.0, single.feed}, "temperature_k");
    expect_solve_rejected(server, complete_made.model_handle(),
                          {0.9e6, 250.0, single.feed}, "pressure_pa");
    const fl::PtFlashRequest closed_edges{1.1e6, 240.0, single.feed};
    compare(native_result(solve(server, complete_made.model_handle(), closed_edges).result()),
            direct(native, closed_edges));
    ok(release(server, complete_made.model_handle()));

    // False/default remains the historical closed endpoint contract.
    auto closed = definition(native);
    closed.applicability.temperature_lower_k = 250.0;
    closed.applicability.temperature_upper_k = 260.0;
    const auto closed_made = create(server, request(closed));
    compare(native_result(solve(server, closed_made.model_handle(), single).result()),
            direct(native, single));
    ok(release(server, closed_made.model_handle()));

    // Exclusivity without its endpoint is a configuration error, never an
    // implicit bound. The field path remains bounded and contains no value echo.
    auto invalid = request(definition(native));
    invalid.mutable_definition()->mutable_applicability()->set_temperature_lower_exclusive(true);
    wire::CreateModelResponse create_out;
    const auto status = create(server, invalid, create_out);
    error(status, SC::INVALID_ARGUMENT, "configuration.invalid_range");
    require(detail_field(status) == "applicability.temperature_lower_exclusive",
            "exclusive-without-endpoint field location changed across wire");
    require(!create_out.has_model_handle() && registry.status().resident_models == 0U,
            "invalid applicability consumed registry capacity");
}
} // namespace

int main() {
    try {
        run();
        std::cout << "PASS grpc_applicability_endpoints\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "FAIL: " << error.what() << '\n';
        return 1;
    }
}
