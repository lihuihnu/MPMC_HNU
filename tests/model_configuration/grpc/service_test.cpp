#include "test_support.hpp"

#include <google/protobuf/unknown_field_set.h>
#include <iostream>
#include <limits>

// Match the adapter's local Win32/Protobuf reflection macro isolation.
#if defined(GetMessage)
#undef GetMessage
#endif

namespace service_test {
namespace {
using SC = grpc::StatusCode;
void roundtrip() {
    mc::Pr76ModelRegistry registry;
    Server server(registry);
    const auto native = binary_model();
    auto r = request(definition(native));
    const auto made = create(server, r);
    require(made.snapshot().definition().SerializeAsString() == r.definition().SerializeAsString(),
            "parameter/provenance roundtrip lost fields");
    require(registry.describe(made.model_handle()).settings == preset(), "preset did not resolve explicitly");
    wire::DescribeModelResponse description;
    ok(describe(server, made.model_handle(), description));
    require(description.snapshot().SerializeAsString() == made.snapshot().SerializeAsString(), "describe differs from creation");
    r.mutable_definition()->clear_components();
    compare(native_result(solve(server, made.model_handle(), binary).result()), direct(native, binary));
    ok(release(server, made.model_handle()));
    error(describe(server, made.model_handle(), description), SC::NOT_FOUND, "registry.model_not_found");
    wire::SolveModelResponse result;
    error(solve(server, solve_request(made.model_handle(), binary), result), SC::NOT_FOUND, "registry.model_not_found");
    error(release(server, made.model_handle()), SC::NOT_FOUND, "registry.model_not_found");
    require(registry.status().resident_models == 0, "release leaked model capacity");
}
void metadata_roundtrip() {
    auto d = definition(pr76_max3_test::model());
    d.components[1].kind = mc::ComponentKind::pseudo;
    d.components[1].display_name = "synthetic pseudo component";
    d.components[1].molar_mass_kg_per_mol = mc::ModelScalar{0.1, d.provenance, "kg/mol", "identity"};
    d.applicability.temperature_lower_k = 240.0; d.applicability.temperature_upper_k = 260.0;
    d.applicability.pressure_lower_pa = 0.9e6; d.applicability.pressure_upper_pa = 1.1e6;
    mc::Pr76ModelRegistry registry({}, synthetic); Server server(registry);
    const auto r = request(d); const auto made = create(server, r);
    require(made.snapshot().definition().SerializeAsString() == r.definition().SerializeAsString(),
            "pseudo/scalar/bounds metadata lost");
    require(made.snapshot().parameter_limits().max_components() == registry.limits().parameters.max_components &&
            made.snapshot().solver_limits().max_root_iterations() == registry.limits().solver.max_root_iterations,
            "host ceilings not exposed in snapshot");
    const auto stored = registry.describe(made.model_handle());
    require(stored.definition.components[1].molar_mass_kg_per_mol->value == 0.1 &&
            stored.definition.components[1].kind == mc::ComponentKind::pseudo, "native metadata mismatch");
    wire::SolveModelResponse result;
    error(solve(server, solve_request(made.model_handle(), {2e6, 250.0, single.feed}), result),
          SC::INVALID_ARGUMENT, "request.rejected");
    require(native_result(solve(server, made.model_handle()).result()).solution.accepted_phase_count() == 1,
            "valid solve after domain error failed");
}
void full_result_parity() {
    mc::Pr76ModelRegistry registry({}, synthetic); Server server(registry);
    const auto made = create(server);
    for (const auto& pt : {single, ternary}) {
        compare(native_result(solve(server, made.model_handle(), pt).result()), direct(pr76_max3_test::model(), pt));
    }
    require(solve(server, made.model_handle(), ternary).result().outcome() == old::PT_COMPUTATION_OUTCOME_INDETERMINATE,
            "cold ternary result was promoted");
}
void custom_settings() {
    mc::Pr76ModelRegistry registry({}, synthetic); Server server(registry);
    const auto a = create(server);
    auto s = custom();
    s.eos_root.max_iterations = 1;
    // Distinct final-search values ensure wire groups cannot alias each other.
    s.initial_stability.max_backtracks = 21;
    s.final_two_phase_stability.max_backtracks = 22;
    s.final_three_phase_stability.max_backtracks = 23;
    s.two_phase.max_split_attempts = 0;
    s.final_three_phase_stability.automatic_multistart = false;
    auto r = request(); api::encode_settings(s, *r.mutable_settings());
    const auto b = create(server, r);
    require(registry.describe(b.model_handle()).settings == s, "wire -> native settings mismatch");
    require(b.snapshot().settings().SerializeAsString() == r.settings().SerializeAsString(),
            "zero/false or settings presence lost on output");
    fl::Pr76PtFlashBackendOptions options;
    options.split.initial_stability.max_backtracks = 21;
    options.split.final_stability.max_backtracks = 22;
    options.final_three_phase_stability.max_backtracks = 23;
    options.split.max_split_attempts = 0;
    options.final_three_phase_stability.automatic_starts = false;
    compare(native_result(solve(server, b.model_handle()).result()), direct(pr76_max3_test::model(), single, {1}, options));
    compare(native_result(solve(server, a.model_handle()).result()), direct(pr76_max3_test::model(), single));
}
void settings_presence() {
    mc::Pr76ModelRegistry registry({}, synthetic); Server server(registry);
    auto base = request(); api::encode_settings(custom(), *base.mutable_settings());
    std::size_t checked = 0;
    for (const std::string group : {"eos_root", "initial_stability", "two_phase", "final_two_phase_stability",
                                   "three_phase", "final_three_phase_stability"}) {
        const auto* descriptor = base.settings().GetDescriptor()->FindFieldByName(group);
        const auto& message = base.settings().GetReflection()->GetMessage(base.settings(), descriptor);
        for (int i = 0; i < message.GetDescriptor()->field_count(); ++i) {
            auto r = base;
            auto* settings = r.mutable_settings();
            auto* nested = settings->GetReflection()->MutableMessage(settings, descriptor);
            nested->GetReflection()->ClearField(nested, nested->GetDescriptor()->field(i));
            wire::CreateModelResponse out;
            error(create(server, r, out), SC::INVALID_ARGUMENT, "configuration.missing_field");
            ++checked;
        }
    }
    require(checked == 57 && registry.status().resident_models == 0, "numerical wire presence coverage or rollback");
}
void wire_errors() {
    mc::Pr76ModelRegistry registry({}, synthetic); Server server(registry);
    wire::CreateModelResponse out;
    auto r = request(); r.clear_wire_contract();
    error(create(server, r, out), SC::INVALID_ARGUMENT, "wire.missing_field");
    r = request(); r.set_wire_contract("future-service/v99");
    error(create(server, r, out), SC::INVALID_ARGUMENT, "wire.unsupported_version");
    r = request(); r.clear_definition();
    error(create(server, r, out), SC::INVALID_ARGUMENT, "wire.missing_field");
    r = request(); r.clear_solver_selection();
    error(create(server, r, out), SC::INVALID_ARGUMENT, "wire.missing_field");
    r = request();
    r.mutable_definition()->GetReflection()->MutableUnknownFields(r.mutable_definition())->AddVarint(900, 1);
    error(create(server, r, out), SC::INVALID_ARGUMENT, "wire.unknown_field");
    r = request(); r.mutable_definition()->set_family(static_cast<wire::ModelFamily>(99));
    error(create(server, r, out), SC::INVALID_ARGUMENT, "configuration.invalid_value");
    r = request(); r.mutable_definition()->mutable_pr76()->mutable_pure(0)->mutable_critical_temperature_k()->clear_value();
    error(create(server, r, out), SC::INVALID_ARGUMENT, "configuration.missing_field");
    require(registry.status().resident_models == 0, "invalid wire published a model");
}
void domain_errors() {
    mc::Pr76ModelRegistry registry({}, synthetic); Server server(registry);
    wire::CreateModelResponse out;
    auto r = request(); r.mutable_definition()->set_version("future-parameters/v99");
    error(create(server, r, out), SC::UNIMPLEMENTED, "configuration.unsupported_version");
    for (const auto family : {wire::MODEL_FAMILY_SW92, wire::MODEL_FAMILY_CPA}) {
        r = request(); r.mutable_definition()->set_family(family);
        error(create(server, r, out), SC::UNIMPLEMENTED, "configuration.unsupported_family");
    }
    r = request(); r.set_preset_id("unknown-preset");
    error(create(server, r, out), SC::UNIMPLEMENTED, "configuration.unsupported_preset");
    r = request(); r.mutable_definition()->mutable_pr76()->clear_binary();
    error(create(server, r, out), SC::INVALID_ARGUMENT, "configuration.missing_parameter");
    r = request(); auto s = preset(); s.eos_root.max_iterations = 1;
    api::encode_settings(s, *r.mutable_settings());
    error(create(server, r, out), SC::INVALID_ARGUMENT, "configuration.invalid_settings");
    const auto made = create(server);
    for (const auto& pt : std::vector<fl::PtFlashRequest>{{0.0, 250.0, single.feed}, {1e6, -1.0, single.feed},
            {1e6, 250.0, {1.0}}, {1e6, 250.0, {-0.1, 0.5, 0.6}},
            {std::numeric_limits<double>::quiet_NaN(), 250.0, single.feed}}) {
        wire::SolveModelResponse response;
        error(solve(server, solve_request(made.model_handle(), pt), response), SC::INVALID_ARGUMENT, "request.rejected");
    }
    compare(native_result(solve(server, made.model_handle()).result()), direct(pr76_max3_test::model(), single));
}
void registry_errors() {
    mc::Pr76ModelRegistryLimits limits; limits.max_models = 1;
    mc::Pr76ModelRegistry registry(limits, synthetic); Server server(registry);
    const auto a = create(server); wire::CreateModelResponse out;
    error(create(server, request(), out), SC::RESOURCE_EXHAUSTED, "registry.capacity_exceeded");
    error(release(server, "bad-handle"), SC::INVALID_ARGUMENT, "registry.invalid_handle");
    error(release(server, "mh1_" + std::string(80, 'f')), SC::NOT_FOUND, "registry.model_not_found");
    ok(release(server, a.model_handle()));
    const auto b = create(server);
    require(a.model_handle() != b.model_handle(), "RPC handle reused");
    registry.close();
    error(create(server, request(), out), SC::FAILED_PRECONDITION, "registry.closed");
    error(release(server, b.model_handle()), SC::FAILED_PRECONDITION, "registry.closed");
}
void host_quota() {
    mc::Pr76ModelRegistryLimits limits; limits.parameters.max_components = 2;
    mc::Pr76ModelRegistry registry(limits, synthetic); Server server(registry);
    wire::CreateModelResponse out;
    error(create(server, request(), out), SC::RESOURCE_EXHAUSTED, "configuration.resource_limit");
    require(registry.status().resident_models == 0, "host quota allowed a model");
}
void busy_and_release() {
    mc::Pr76ModelRegistry registry({}, synthetic); Server server(registry);
    const auto made = create(server);
    auto lease = registry.acquire_solve(made.model_handle());
    wire::SolveModelResponse out;
    error(solve(server, solve_request(made.model_handle()), out), SC::RESOURCE_EXHAUSTED, "model.busy");
    ok(release(server, made.model_handle()));
    error(solve(server, solve_request(made.model_handle()), out), SC::NOT_FOUND, "registry.model_not_found");
    require(registry.status().resident_models == 1, "RPC release destroyed admitted model");
    compare(std::move(lease).solve(single), direct(pr76_max3_test::model(), single));
    require(registry.status().resident_models == 0, "lease completion did not return capacity");
}
void request_size() {
    mc::Pr76ModelRegistry registry({}, synthetic);
    auto r = request(); api::ModelGrpcLimits limits; limits.max_request_bytes = r.ByteSizeLong() + 128;
    std::atomic<int> authorized{0};
    Server server(registry, limits, [&](const grpc::ServerContext&) { authorized.fetch_add(1); return true; });
    r.mutable_definition()->set_display_name(std::string(20000, 'x'));
    wire::CreateModelResponse out;
    const auto status = create(server, r, out);
    require(status.error_code() == SC::RESOURCE_EXHAUSTED && authorized.load() == 0 &&
            registry.status().resident_models == 0, "oversized message reached typed handler/deserialization dispatch");
    (void)create(server);
}
void response_rollback() {
    mc::Pr76ModelRegistry registry({}, synthetic); api::ModelGrpcLimits limits; limits.max_response_bytes = 1;
    Server server(registry, limits); wire::CreateModelResponse out;
    error(create(server, request(), out), SC::RESOURCE_EXHAUSTED, "wire.response_limit");
    require(!out.has_model_handle() && registry.status().resident_models == 0, "undelivered handle leaked");
    const auto handle = registry.create(definition(pr76_max3_test::model()), preset());
    wire::DescribeModelResponse description;
    error(describe(server, handle, description), SC::RESOURCE_EXHAUSTED, "wire.response_limit");
    wire::SolveModelResponse solution;
    error(solve(server, solve_request(handle), solution), SC::RESOURCE_EXHAUSTED, "wire.response_limit");
    error(release(server, handle), SC::RESOURCE_EXHAUSTED, "wire.response_limit");
    require(registry.status().resident_models == 0, "release was undone by response failure");
}
void authorization_and_internal_errors() {
    mc::Pr76ModelRegistry registry({}, synthetic);
    {
        Server server(registry, {}, [](const grpc::ServerContext&) { return false; });
        wire::CreateModelResponse out;
        error(create(server, request(), out), SC::PERMISSION_DENIED, "rpc.not_authorized");
    }
    {
        Server server(registry, {}, [](const grpc::ServerContext&) -> bool { throw std::runtime_error("private-host-detail"); });
        wire::CreateModelResponse out; const auto status = create(server, request(), out);
        error(status, SC::INTERNAL, "rpc.internal_failure");
        require(status.error_message().find("private-host-detail") == std::string::npos, "internal error leaked");
    }
    mc::Pr76ModelRegistry no_entropy({}, synthetic, []() -> mc::ModelHandleEntropy {
        throw mc::ModelRegistryError(mc::ModelRegistryErrorCode::entropy_unavailable, "test entropy failure");
    });
    Server server(no_entropy); wire::CreateModelResponse out;
    error(create(server, request(), out), SC::UNAVAILABLE, "registry.entropy_unavailable");
    require(no_entropy.status().resident_models == 0 && registry.status().resident_models == 0, "failed RPC created model");
}
void deadline_policy() {
    mc::Pr76ModelRegistry registry({}, synthetic); Server server(registry);
    wire::CreateModelResponse out;
    grpc::ClientContext missing;
    error(server.stub().CreateModel(&missing, request(), &out), SC::INVALID_ARGUMENT, "rpc.deadline_required");
    grpc::ClientContext too_long; too_long.set_deadline(std::chrono::system_clock::now() + 1h);
    error(server.stub().CreateModel(&too_long, request(), &out), SC::INVALID_ARGUMENT, "rpc.deadline_limit");
    grpc::ClientContext expired; expired.set_deadline(std::chrono::system_clock::now() - 1s);
    require(server.stub().CreateModel(&expired, request(), &out).error_code() == SC::DEADLINE_EXCEEDED,
            "expired deadline accepted");
    require(registry.status().resident_models == 0, "deadline failure created a model");
}
void cancellation_and_admission() {
    Gate gate;
    mc::Pr76ModelRegistry registry({}, synthetic, [&] { gate.pause(); return mc::system_model_handle_entropy(); });
    api::ModelGrpcLimits limits; limits.max_concurrent_requests = 1;
    Server server(registry, limits);
    grpc::ClientContext context; deadline(context);
    auto task = std::async(std::launch::async, [&] {
        wire::CreateModelResponse response; return server.stub().CreateModel(&context, request(), &response);
    });
    const ResumeOnExit cleanup{gate}; gate.entered.wait();
    wire::CreateModelResponse second;
    error(create(server, request(), second), SC::RESOURCE_EXHAUSTED, "rpc.concurrency_limit");
    context.TryCancel();
    require(task.get().error_code() == SC::CANCELLED, "client cancellation was lost");
    gate.open(); server.stop(); // Wait for the server's rollback, not merely the client return.
    require(registry.status().resident_models == 0 && server.adapter().in_flight_requests() == 0,
            "cancelled creation leaked capacity/admission");
}
void legacy_coexistence() {
    const auto native = binary_model();
    auto owner = mc::make_pr76_executable_model(definition(native), preset());
    const std::vector<mpmc::runtime::PtServiceBackendRegistration> registrations{{"legacy.binary", owner.get()}};
    mpmc::runtime::PtService legacy(registrations);
    mpmc::runtime_grpc::PtGrpcServiceAdapter adapter(legacy);
    mc::Pr76ModelRegistry registry;
    api::ModelGrpcLimits limits; limits.max_request_bytes = adapter.limits().max_serialized_request_bytes;
    Server server(registry, limits, [](const grpc::ServerContext&) { return true; }, &adapter);
    auto stub = old::PtFlashService::NewStub(server.channel());
    old::SolvePtFlashRequest r;
    r.set_wire_contract("mpmc.runtime.v1/PT-flash-service/v1");
    r.set_service_request_convention(std::string(mpmc::runtime::PtServiceRequest::convention));
    r.set_configured_backend_id("legacy.binary"); r.set_pressure_pa(binary.pressure_pa); r.set_temperature_k(binary.temperature_k);
    for (std::size_t i = 0; i < binary.feed.size(); ++i) {
        auto* entry = r.add_feed(); entry->set_component_id(native.parameters().components().at(i).id);
        entry->set_mole_fraction(binary.feed[i]);
    }
    old::SolvePtFlashResponse response; grpc::ClientContext first; deadline(first);
    ok(stub->SolvePtFlash(&first, r, &response));
    require(response.has_result() && response.result().phases_size() == 2, "old v1 solve changed");
    const auto made = create(server, request(definition(native)));
    compare(native_result(solve(server, made.model_handle(), binary).result()), direct(native, binary));
    r.set_configured_backend_id(made.model_handle()); grpc::ClientContext wrong; deadline(wrong);
    ok(stub->SolvePtFlash(&wrong, r, &response));
    require(response.has_error() && response.error().code() == old::PT_SERVICE_ERROR_CODE_CONFIGURED_BACKEND_NOT_FOUND,
            "v1 backend ID was reinterpreted as a dynamic handle");
}
} // namespace
} // namespace service_test
int main(int argc, char** argv) {
    using namespace service_test;
    try {
        require(argc == 2, "expected one service test"); const std::string_view name = argv[1];
        if (name == "roundtrip") { roundtrip(); }
        else if (name == "metadata_roundtrip") { metadata_roundtrip(); }
        else if (name == "full_result_parity") { full_result_parity(); }
        else if (name == "custom_settings") { custom_settings(); }
        else if (name == "settings_presence") { settings_presence(); }
        else if (name == "wire_errors") { wire_errors(); }
        else if (name == "domain_errors") { domain_errors(); }
        else if (name == "registry_errors") { registry_errors(); }
        else if (name == "host_quota") { host_quota(); }
        else if (name == "busy_and_release") { busy_and_release(); }
        else if (name == "request_size") { request_size(); }
        else if (name == "response_rollback") { response_rollback(); }
        else if (name == "authorization_and_internal_errors") { authorization_and_internal_errors(); }
        else if (name == "deadline_policy") { deadline_policy(); }
        else if (name == "cancellation_and_admission") { cancellation_and_admission(); }
        else if (name == "legacy_coexistence") { legacy_coexistence(); }
        else { throw std::runtime_error("unknown service test"); }
        std::cout << "PASS " << name << '\n'; return 0;
    } catch (const std::exception& e) { std::cerr << "FAIL: " << e.what() << '\n'; return 1; }
}
