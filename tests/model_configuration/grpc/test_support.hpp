#ifndef MPMC_MODEL_GRPC_TEST_SUPPORT_HPP
#define MPMC_MODEL_GRPC_TEST_SUPPORT_HPP
#include <mpmc/model_configuration_grpc/model_grpc_adapter.hpp>
#include <mpmc/runtime_grpc/pt_grpc_adapter.hpp>
#include <google/rpc/status.pb.h>
#include "../executable_model/test_support.hpp"

#include <chrono>
#include <future>
#include <latch>
#include <vector>

namespace service_test {
using namespace model_test;
namespace api = mpmc::model_configuration_grpc;
namespace wire = mpmc::model_configuration::v1;
namespace old = mpmc::runtime::v1;
using namespace std::chrono_literals;
constexpr auto synthetic = mc::ModelDataPolicy::allow_synthetic_tests;

class Server {
public:
    Server(mc::Pr76ModelRegistry& registry, api::ModelGrpcLimits limits = {},
           api::ModelGrpcAuthorization authorization = [](const grpc::ServerContext&) { return true; },
           mpmc::runtime_grpc::PtGrpcServiceAdapter* legacy = nullptr)
        : adapter_(registry, std::move(authorization), limits) {
        grpc::ServerBuilder builder;
        int port = 0;
        // Insecure credentials ONLY for isolated loopback test traffic.
        builder.AddListeningPort("127.0.0.1:0", grpc::InsecureServerCredentials(), &port);
        if (legacy) { mpmc::runtime_grpc::configure_pt_grpc_server(builder, *legacy); }
        api::configure_model_grpc_server(builder, adapter_);
        server_ = builder.BuildAndStart();
        require(server_ != nullptr && port > 0, "could not start model RPC test server");
        channel_ = grpc::CreateChannel("127.0.0.1:" + std::to_string(port), grpc::InsecureChannelCredentials());
        stub_ = wire::ModelConfigurationService::NewStub(channel_);
    }
    ~Server() { stop(); }
    void stop() {
        if (server_) { server_->Shutdown(std::chrono::system_clock::now() + 2s); server_->Wait(); server_.reset(); }
    }
    auto& stub() { return *stub_; }
    auto& adapter() { return adapter_; }
    auto channel() { return channel_; }
private:
    api::ModelGrpcServiceAdapter adapter_;
    std::unique_ptr<grpc::Server> server_;
    std::shared_ptr<grpc::Channel> channel_;
    std::unique_ptr<wire::ModelConfigurationService::Stub> stub_;
};
inline void deadline(grpc::ClientContext& context) {
    context.set_deadline(std::chrono::system_clock::now() + 10s);
}
template <class T> T versioned() {
    T value; value.set_wire_contract(std::string(api::model_wire_contract)); return value;
}
inline wire::CreateModelRequest request(const mc::ThermodynamicModelDefinition& definition) {
    auto r = versioned<wire::CreateModelRequest>();
    api::encode_definition(definition, *r.mutable_definition());
    r.set_preset_id(std::string(mc::mpmc_balanced_default_v1));
    return r;
}
inline auto request() { return request(definition(pr76_max3_test::model())); }
inline wire::ModelHandleRequest handle_request(const std::string& handle) {
    auto r = versioned<wire::ModelHandleRequest>(); r.set_model_handle(handle); return r;
}
inline wire::SolveModelRequest solve_request(const std::string& handle, const fl::PtFlashRequest& pt = single) {
    auto r = versioned<wire::SolveModelRequest>(); r.set_model_handle(handle);
    r.set_pressure_pa(pt.pressure_pa); r.set_temperature_k(pt.temperature_k);
    for (double z : pt.feed) { r.add_feed(z); }
    return r;
}
inline void ok(const grpc::Status& status) {
    if (!status.ok()) { throw std::runtime_error("RPC unexpectedly failed: " + status.error_message()); }
}
inline void error(const grpc::Status& status, grpc::StatusCode expected, std::string_view code) {
    require(status.error_code() == expected, "wrong gRPC status");
    wire::ModelServiceError details;
    google::rpc::Status envelope;
    require(envelope.ParseFromString(status.error_details()) && envelope.code() == static_cast<int>(expected) &&
            envelope.message() == status.error_message() && envelope.details_size() == 1 &&
            envelope.details(0).UnpackTo(&details) && details.has_code() &&
            details.wire_contract() == api::model_wire_contract && details.code() == code,
            "wrong structured service error");
    require(status.error_message().find("mh1_") == std::string::npos &&
            status.error_details().find("mh1_") == std::string::npos, "error leaks bearer handle");
}
inline grpc::Status create(Server& server, const wire::CreateModelRequest& r, wire::CreateModelResponse& out) {
    grpc::ClientContext c; deadline(c); return server.stub().CreateModel(&c, r, &out);
}
inline wire::CreateModelResponse create(Server& server, const wire::CreateModelRequest& r = request()) {
    wire::CreateModelResponse out; ok(create(server, r, out));
    require(out.has_snapshot() && out.has_model_handle() && out.wire_contract() == api::model_wire_contract,
            "missing creation result"); return out;
}
inline grpc::Status solve(Server& server, const wire::SolveModelRequest& r, wire::SolveModelResponse& out) {
    grpc::ClientContext c; deadline(c); return server.stub().SolveModel(&c, r, &out);
}
inline wire::SolveModelResponse solve(Server& server, const std::string& handle, const fl::PtFlashRequest& pt = single) {
    wire::SolveModelResponse out; ok(solve(server, solve_request(handle, pt), out));
    require(out.has_result() && out.wire_contract() == api::model_wire_contract, "missing solve envelope"); return out;
}
inline grpc::Status release(Server& server, const std::string& handle) {
    grpc::ClientContext c; deadline(c); wire::ReleaseModelResponse out;
    return server.stub().ReleaseModel(&c, handle_request(handle), &out);
}
inline grpc::Status describe(Server& server, const std::string& handle, wire::DescribeModelResponse& out) {
    grpc::ClientContext c; deadline(c); return server.stub().DescribeModel(&c, handle_request(handle), &out);
}
// Independent wire -> native reconstruction for comparison with direct PR76.
// Production encode_result is never used to manufacture the expected envelope.
inline fl::PtFlashBackendResult native_result(const wire::FullPtResult& w) {
    require(w.has_outcome() && w.has_capability() && w.has_transition_report() &&
            w.has_global_stability_proven() && w.has_morphology_resolved(), "incomplete result presence");
    require(w.backend_result_convention() == fl::PtFlashBackendResult::convention &&
            w.phase_set_convention() == fl::PtPhaseSetResult::convention, "result conventions");
    fl::PtFlashBackendResult n;
    const auto& c = w.capability(); auto& cap = n.capability;
    require(c.convention() == fl::PtFlashBackendCapability::convention, "capability convention");
#define COPY(field) cap.field = c.field()
    COPY(backend_id); COPY(model_profile); COPY(algorithm_profile); COPY(publication_profile);
    COPY(configuration_profile); COPY(dataset_id); COPY(revision); COPY(phase_metadata_namespace);
    COPY(performs_initial_stability_search); COPY(performs_final_phase_set_review);
    COPY(performs_boundary_neighbor_resolve); COPY(global_stability_proven);
#undef COPY
    cap.component_ids.assign(c.component_ids().begin(), c.component_ids().end());
    cap.supported_phase_counts.assign(c.supported_phase_counts().begin(), c.supported_phase_counts().end());
    for (const auto& item : c.scalar_settings()) { cap.scalar_settings.push_back({item.id(), item.value(), item.unit()}); }
    for (const auto& edge : c.transition_capability().edges()) {
        require(edge.support() != old::PT_PHASE_TRANSITION_SUPPORT_UNSPECIFIED, "missing transition support");
        cap.transition_capability.edges.push_back({edge.source_phase_count(), edge.target_phase_count(),
            static_cast<fl::PtPhaseTransitionSupport>(static_cast<int>(edge.support()) - 1), edge.requires_fresh_target_solve()});
    }
    n.solution.status = static_cast<fl::PtPhaseSetStatus>(static_cast<int>(w.outcome()) - 1);
    n.solution.capability.maximum_phase_count = w.maximum_phase_count();
    n.solution.pressure_pa = w.pressure_pa(); n.solution.temperature_k = w.temperature_k();
    n.solution.feed.assign(w.feed().begin(), w.feed().end());
    if (w.has_candidate_phase_set()) {
        n.solution.candidate_phase_set.emplace();
        for (const auto& p : w.candidate_phase_set().phases()) {
            require(p.has_mole_phase_fraction() && p.has_provider_branch() && p.has_provider_branch_smooth(),
                    "missing phase scalar presence");
            fl::PtCandidatePhase phase;
            phase.mole_phase_fraction = p.mole_phase_fraction();
            phase.composition.assign(p.composition().begin(), p.composition().end());
            phase.activity.ln_phi.assign(p.ln_fugacity_coefficient().begin(), p.ln_fugacity_coefficient().end());
            phase.activity.branch = static_cast<std::size_t>(p.provider_branch());
            phase.activity.smooth = p.provider_branch_smooth();
            if (p.has_compressibility_factor()) { phase.compressibility_factor = p.compressibility_factor(); }
            n.solution.candidate_phase_set->phases.push_back(std::move(phase));
        }
    }
    n.solution.global_stability_proven = w.global_stability_proven(); n.solution.diagnostic = w.diagnostic();
    require(w.transition_report().convention() == fl::PtPhaseTransitionReport::convention, "transition convention");
    for (const auto& e : w.transition_report().evidence()) {
        fl::PtPhaseTransitionEvidence evidence;
        evidence.source_phase_count = e.source_phase_count();
        if (e.has_target_phase_count()) { evidence.target_phase_count = e.target_phase_count(); }
        evidence.trigger = static_cast<fl::PtPhaseTransitionTrigger>(static_cast<int>(e.trigger()) - 1);
        evidence.resolution = static_cast<fl::PtPhaseTransitionResolution>(static_cast<int>(e.resolution()) - 1);
        evidence.fresh_target_solve_attempted = e.fresh_target_solve_attempted();
        evidence.target_topology_closed = e.target_topology_closed();
        evidence.provider_evidence_profile = e.provider_evidence_profile(); evidence.diagnostic = e.diagnostic();
        n.transition_report.evidence.push_back(std::move(evidence));
    }
    n.provider_result_convention = w.provider_result_convention();
    for (const auto& m : w.phase_metadata()) { n.phase_metadata.push_back({m.role_id(), m.family_id()}); }
    n.morphology_resolved = w.morphology_resolved();
    return n;
}
struct Gate {
    std::latch entered{1}; std::latch resume{1}; std::atomic<bool> opened{false};
    void pause() { entered.count_down(); resume.wait(); }
    void open() { if (!opened.exchange(true)) { resume.count_down(); } }
};
struct ResumeOnExit { Gate& gate; ~ResumeOnExit() { gate.open(); } };
} // namespace service_test
#endif
