#include "test_support.hpp"
#include <mpmc/model_configuration_grpc/model_sessions.hpp>
#include <mpmc/pt_process/process_host.hpp>

#include <filesystem>
#include <fstream>
#include <iostream>
#include <thread>

namespace session_test {
using namespace service_test;
namespace process = mpmc::pt_process;
using SC = grpc::StatusCode;
const std::string token(43, 'a'); // Isolated test credential only.
std::string read_file(const std::filesystem::path& path) {
    std::ifstream file(path); require(file.good(), "missing generated test certificate");
    return {std::istreambuf_iterator<char>(file), std::istreambuf_iterator<char>()};
}
template <class Predicate> void eventually(Predicate predicate) {
    const auto end = std::chrono::steady_clock::now() + 8s;
    while (!predicate()) {
        require(std::chrono::steady_clock::now() < end, "session cleanup did not complete");
        std::this_thread::sleep_for(10ms);
    }
}
process::PtProcessHostOptions desktop_options() {
    process::PtProcessHostOptions options;
    options.desktop_loopback_session = true; options.listen_address = "127.0.0.1:0";
    options.enable_model_sessions = true; options.shutdown_grace = 1s;
    return options;
}
struct Host {
    std::shared_ptr<mc::Pr76ExecutableModel> model{mc::make_pr76_executable_model(definition(binary_model()), preset())};
    process::PtCompositionRoot root;
    process::PtProcessHost host;
    explicit Host(process::PtProcessHostOptions options = desktop_options())
        : root({{"legacy.binary", model}}, {}, {}, {}, {options.desktop_loopback_session ? token : ""}),
          host(root, std::move(options)) { host.start(); }
    std::string address() const { return "127.0.0.1:" + std::to_string(host.selected_port()); }
    std::shared_ptr<grpc::Channel> channel() const {
        return grpc::CreateChannel(address(), grpc::InsecureChannelCredentials());
    }
};
void metadata(grpc::ClientContext& context, const std::string& sid, const std::string& bearer = token) {
    deadline(context);
    if (!bearer.empty()) { context.AddMetadata("authorization", "Bearer " + bearer); }
    if (!sid.empty()) { context.AddMetadata(api::model_session_metadata, sid); }
}
class Session {
public:
    grpc::ClientContext context;
    std::unique_ptr<wire::ModelSessionService::Stub> stub;
    std::unique_ptr<grpc::ClientReader<wire::ModelSessionOpened>> reader;
    std::string id;
    Session(std::shared_ptr<grpc::Channel> channel, const std::string& bearer = token,
            std::chrono::milliseconds age = 10s) : stub(wire::ModelSessionService::NewStub(channel)) {
        metadata(context, "", bearer); context.set_deadline(std::chrono::system_clock::now() + age);
        wire::OpenModelSessionRequest r; r.set_wire_contract(std::string(api::model_session_contract));
        reader = stub->OpenModelSession(&context, r);
        wire::ModelSessionOpened opened;
        if (!reader->Read(&opened)) { const auto s = reader->Finish(); throw std::runtime_error("open failed: " + std::to_string(s.error_code())); }
        require(opened.wire_contract() == api::model_session_contract && opened.has_session_id(), "invalid session handshake");
        id = opened.session_id();
    }
    ~Session() { if (reader) { context.TryCancel(); (void)reader->Finish(); } }
    void cancel() { context.TryCancel(); (void)reader->Finish(); reader.reset(); }
};
wire::CreateModelResponse make(wire::ModelConfigurationService::Stub& stub, const std::string& sid,
                               const std::string& bearer = token) {
    grpc::ClientContext c; metadata(c, sid, bearer);
    wire::CreateModelResponse out; ok(stub.CreateModel(&c, request(definition(binary_model())), &out)); return out;
}
grpc::Status inspect(wire::ModelConfigurationService::Stub& stub, const std::string& sid,
                     const std::string& handle, const std::string& bearer = token) {
    grpc::ClientContext c; metadata(c, sid, bearer);
    wire::DescribeModelResponse out; return stub.DescribeModel(&c, handle_request(handle), &out);
}
void roundtrip() {
    Host h; auto channel = h.channel(); auto stub = wire::ModelConfigurationService::NewStub(channel);
    Session a(channel), b(channel);
    const auto made = make(*stub, a.id);
    ok(inspect(*stub, a.id, made.model_handle()));
    error(inspect(*stub, b.id, made.model_handle()), SC::NOT_FOUND, "registry.model_not_found");
    for (bool release_call : {false, true}) {
        grpc::ClientContext c; metadata(c, b.id);
        if (release_call) {
            wire::ReleaseModelResponse out;
            error(stub->ReleaseModel(&c, handle_request(made.model_handle()), &out), SC::NOT_FOUND, "registry.model_not_found");
        } else {
            wire::SolveModelResponse out;
            error(stub->SolveModel(&c, solve_request(made.model_handle(), binary), &out), SC::NOT_FOUND, "registry.model_not_found");
        }
    }
    grpc::ClientContext c; metadata(c, a.id); wire::SolveModelResponse out;
    ok(stub->SolveModel(&c, solve_request(made.model_handle(), binary), &out));
    compare(native_result(out.result()), direct(binary_model(), binary));
    grpc::ClientContext release_context; metadata(release_context, a.id); wire::ReleaseModelResponse released;
    ok(stub->ReleaseModel(&release_context, handle_request(made.model_handle()), &released));
    error(inspect(*stub, a.id, made.model_handle()), SC::NOT_FOUND, "registry.model_not_found");
    // Existing v1 discovery continues to use its original bearer and wire contract.
    auto old_stub = old::PtFlashService::NewStub(channel);
    grpc::ClientContext old_context; metadata(old_context, "");
    old::DiscoverPtCapabilitiesRequest discovery;
    old::DiscoverPtCapabilitiesResponse discovered;
    ok(old_stub->DiscoverPtCapabilities(&old_context, discovery, &discovered));
}
void authentication() {
    Host h; auto channel = h.channel(); Session a(channel);
    auto stub = wire::ModelConfigurationService::NewStub(channel); const auto made = make(*stub, a.id);
    for (const auto& bad : {std::string{}, std::string(43, 'b')}) {
        error(inspect(*stub, a.id, made.model_handle(), bad), SC::UNAUTHENTICATED, "session.unauthenticated");
    }
    error(inspect(*stub, "", made.model_handle()), SC::INVALID_ARGUMENT, "session.metadata_required");
    error(inspect(*stub, "bad", made.model_handle()), SC::INVALID_ARGUMENT, "session.invalid_id");
    for (bool duplicate_auth : {false, true}) {
        grpc::ClientContext c; metadata(c, a.id);
        c.AddMetadata(duplicate_auth ? "authorization" : api::model_session_metadata,
                      duplicate_auth ? "Bearer " + token : a.id);
        wire::DescribeModelResponse response;
        error(stub->DescribeModel(&c, handle_request(made.model_handle()), &response),
              duplicate_auth ? SC::UNAUTHENTICATED : SC::INVALID_ARGUMENT,
              duplicate_auth ? "session.unauthenticated" : "session.metadata_required");
    }
    ok(inspect(*stub, a.id, made.model_handle()));
    auto options = desktop_options(); options.enable_model_sessions = false; Host disabled(options);
    auto session_stub = wire::ModelSessionService::NewStub(disabled.channel()); grpc::ClientContext c; metadata(c, "");
    wire::OpenModelSessionRequest r; r.set_wire_contract(std::string(api::model_session_contract));
    auto reader = session_stub->OpenModelSession(&c, r); wire::ModelSessionOpened opened;
    require(!reader->Read(&opened) && reader->Finish().error_code() == SC::UNIMPLEMENTED, "service enabled implicitly");
}
void disconnect() {
    auto options = desktop_options(); options.model_sessions.max_sessions = 1;
    Host h(options); auto channel = h.channel(); auto stub = wire::ModelConfigurationService::NewStub(channel);
    Session a(channel); const auto old_id = a.id; const auto made = make(*stub, a.id);
    a.cancel(); eventually([&] { return h.host.model_session_status().resident_sessions == 0; });
    error(inspect(*stub, old_id, made.model_handle()), SC::NOT_FOUND, "session.not_found");
    Session b(channel); require(b.id != old_id, "session ID was reused");
    error(inspect(*stub, b.id, made.model_handle()), SC::NOT_FOUND, "registry.model_not_found");
    (void)make(*stub, b.id);
}
void open_errors() {
    auto options = desktop_options(); options.model_sessions.max_sessions = 1;
    Host h(options); auto channel = h.channel(); auto stub = wire::ModelSessionService::NewStub(channel);
    const auto reject = [&](const char* version, const std::string& bearer, bool with_deadline,
                            SC code, const char* detail) {
        grpc::ClientContext c;
        if (with_deadline) { deadline(c); }
        if (!bearer.empty()) { c.AddMetadata("authorization", "Bearer " + bearer); }
        wire::OpenModelSessionRequest r; r.set_wire_contract(version);
        auto reader = stub->OpenModelSession(&c, r); wire::ModelSessionOpened opened;
        require(!reader->Read(&opened), "invalid open published a session");
        error(reader->Finish(), code, detail);
    };
    const auto version = std::string(api::model_session_contract);
    reject(version.c_str(), "", true, SC::UNAUTHENTICATED, "session.unauthenticated");
    reject(version.c_str(), token, false, SC::INVALID_ARGUMENT, "session.deadline_required");
    reject("future", token, true, SC::INVALID_ARGUMENT, "session.invalid_request");
    require(h.host.model_session_status().resident_sessions == 0, "invalid opens consumed capacity");
    Session a(channel);
    reject(version.c_str(), token, true, SC::RESOURCE_EXHAUSTED, "session.capacity_exceeded");
    a.cancel(); eventually([&] { return h.host.model_session_status().resident_sessions == 0; });
    auto small_options = desktop_options(); small_options.model_sessions.rpc.max_response_bytes = 1;
    Host small(small_options); auto small_stub = wire::ModelSessionService::NewStub(small.channel());
    grpc::ClientContext c; metadata(c, ""); wire::OpenModelSessionRequest r; r.set_wire_contract(version);
    auto reader = small_stub->OpenModelSession(&c, r); wire::ModelSessionOpened opened;
    require(!reader->Read(&opened), "oversize handshake delivered");
    error(reader->Finish(), SC::CANCELLED, "session.publication_failed");
    eventually([&] { return small.host.model_session_status().resident_sessions == 0; });
}
void expiry_and_shutdown() {
    Host h; auto channel = h.channel(); auto stub = wire::ModelConfigurationService::NewStub(channel);
    Session a(channel, token, 2s); const auto made = make(*stub, a.id);
    eventually([&] { return h.host.model_session_status().resident_sessions == 0; });
    error(inspect(*stub, a.id, made.model_handle()), SC::NOT_FOUND, "session.not_found");
    Session b(channel); (void)make(*stub, b.id);
    h.host.shutdown(); h.host.wait();
    require(h.host.model_session_status().closed && h.host.model_session_status().resident_sessions == 0,
            "host shutdown leaked sessions");
}
std::shared_ptr<grpc::Channel> secure_channel(const Host& host, const std::filesystem::path& dir, const char* name) {
    grpc::SslCredentialsOptions tls;
    tls.pem_root_certs = read_file(dir / "ca.crt");
    if (std::string_view(name) != "none") {
        tls.pem_cert_chain = read_file(dir / (std::string(name) + ".crt"));
        tls.pem_private_key = read_file(dir / (std::string(name) + ".key"));
    }
    grpc::ChannelArguments args; args.SetSslTargetNameOverride("localhost");
    return grpc::CreateCustomChannel(host.address(), grpc::SslCredentials(tls), args);
}
void mtls(const std::filesystem::path& dir) {
    process::PtProcessHostOptions options; options.listen_address = "127.0.0.1:0"; options.enable_model_sessions = true;
    options.shutdown_grace = 1s;
    options.tls = process::load_pt_process_tls_identity({dir / "server.crt", dir / "server.key", dir / "ca.crt"});
    Host h(options); auto ca = secure_channel(h, dir, "client-a"), cb = secure_channel(h, dir, "client-b");
    Session a(ca, ""), b(cb, "");
    auto sa = wire::ModelConfigurationService::NewStub(ca), sb = wire::ModelConfigurationService::NewStub(cb);
    const auto made = make(*sa, a.id, "");
    error(inspect(*sb, a.id, made.model_handle(), ""), SC::PERMISSION_DENIED, "session.wrong_identity");
    error(inspect(*sb, b.id, made.model_handle(), ""), SC::NOT_FOUND, "registry.model_not_found");
    ok(inspect(*sa, a.id, made.model_handle(), ""));
    // Both certificates deliberately have the same CN; the leaf identity matters.
    auto no_cert = wire::ModelConfigurationService::NewStub(secure_channel(h, dir, "none"));
    require(!inspect(*no_cert, a.id, made.model_handle(), "").ok(), "missing client certificate accepted");
}
void retained_capacity(bool global_admission = false) {
    Gate gate; std::atomic<int> calls{0}; api::ModelSessionLimits limits;
    limits.max_sessions = global_admission ? 2U : 1U;
    limits.rpc.max_concurrent_requests = 1;
    api::ModelSessionService service([](const grpc::ServerContext&) { return "isolated-test-host"; }, limits, [&] {
        if (calls.fetch_add(1) == (global_admission ? 2 : 1)) { gate.pause(); }
        return mc::system_model_handle_entropy();
    });
    grpc::ServerBuilder builder; int port = 0;
    builder.AddListeningPort("127.0.0.1:0", grpc::InsecureServerCredentials(), &port);
    builder.SetMaxReceiveMessageSize(65536); service.register_services(builder);
    auto server = builder.BuildAndStart(); require(server != nullptr, "test server startup");
    const auto cleanup_server = [&] { service.close(); server->Shutdown(); server->Wait(); };
    auto channel = grpc::CreateChannel("127.0.0.1:" + std::to_string(port), grpc::InsecureChannelCredentials());
    Session a(channel); auto stub = wire::ModelConfigurationService::NewStub(channel);
    std::unique_ptr<Session> other;
    if (global_admission) { other = std::make_unique<Session>(channel); }
    auto pending = std::async(std::launch::async, [&] {
        grpc::ClientContext c; metadata(c, a.id); wire::CreateModelResponse out;
        return stub->CreateModel(&c, request(definition(binary_model())), &out);
    });
    const ResumeOnExit cleanup{gate}; gate.entered.wait();
    if (global_admission) {
        error(inspect(*stub, other->id, "unused"), SC::RESOURCE_EXHAUSTED, "session.capacity_exceeded");
        gate.open(); ok(pending.get());
        (void)make(*stub, other->id);
        cleanup_server();
        require(service.status().resident_sessions == 0, "global admission or shutdown leaked a session");
        return;
    }
    a.cancel(); eventually([&] { return service.status().active_sessions == 0; });
    require(service.status().resident_sessions == 1, "retained creation bypassed session capacity");
    grpc::ClientContext c; metadata(c, ""); wire::OpenModelSessionRequest r;
    r.set_wire_contract(std::string(api::model_session_contract));
    auto reader = a.stub->OpenModelSession(&c, r); wire::ModelSessionOpened opened;
    require(!reader->Read(&opened), "over-capacity session published");
    error(reader->Finish(), SC::RESOURCE_EXHAUSTED, "session.capacity_exceeded");
    gate.open(); error(pending.get(), SC::FAILED_PRECONDITION, "registry.closed");
    eventually([&] { return service.status().resident_sessions == 0; });
    { Session b(channel); (void)make(*stub, b.id); }
    cleanup_server(); require(service.status().resident_sessions == 0, "closed service retained capacity");
}
} // namespace session_test
int main(int argc, char** argv) {
    using namespace session_test;
    try {
        require(argc >= 2, "expected session case"); const std::string_view name(argv[1]);
        if (name == "roundtrip") { roundtrip(); }
        else if (name == "authentication") { authentication(); }
        else if (name == "disconnect") { disconnect(); }
        else if (name == "open_errors") { open_errors(); }
        else if (name == "expiry_and_shutdown") { expiry_and_shutdown(); }
        else if (name == "mtls") { require(argc == 3, "certificate directory required"); mtls(argv[2]); }
        else if (name == "retained_capacity") { retained_capacity(); }
        else if (name == "global_admission") { retained_capacity(true); }
        else { throw std::runtime_error("unknown session case"); }
        std::cout << "PASS " << name << '\n'; return 0;
    } catch (const std::exception& e) { std::cerr << "FAIL " << e.what() << '\n'; return 1; }
}
