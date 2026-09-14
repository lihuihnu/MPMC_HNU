#include "test_support.hpp"
#include <mpmc/model_configuration_grpc/model_sessions.hpp>
#include <google/protobuf/util/json_util.h>
#include <iostream>
#include <thread>

namespace {
using namespace service_test;

mc::PtSolveHints binary_continuation_hints() {
    auto owner = model_test::create(definition(binary_model()), preset());
    const auto cold = owner->solve(binary);
    const auto hints = mc::make_pr76_continuation_hints(cold);
    require(hints.version == mc::pt_solve_hints_v1,
            "binary continuation helper lost the public hint version");
    return hints;
}

wire::SolveModelRequest binary_hinted_request(
    const std::string& handle, const mc::PtSolveHints& hints) {
    auto request = solve_request(handle, binary);
    api::encode_solve_hints(hints, *request.mutable_hints());
    return request;
}
} // namespace

int main(int argc, char** argv) {
    using namespace service_test;
    try {
        const auto hints = binary_continuation_hints();
        if (argc == 2 && std::string_view(argv[1]) == "--fixture-json") {
            // Export the existing attributed native fixture; no duplicate TS data.
            std::string json;
            require(google::protobuf::util::MessageToJsonString(
                request(definition(binary_model())), &json).ok(), "fixture encoding failed");
            std::string solve_json;
            require(google::protobuf::util::MessageToJsonString(
                solve_request("", binary), &solve_json).ok(), "state encoding failed");
            std::string hinted_json;
            require(google::protobuf::util::MessageToJsonString(
                binary_hinted_request("", hints), &hinted_json).ok(),
                "hinted state encoding failed");
            std::cout << "{\"create\":" << json
                      << ",\"solve\":" << solve_json
                      << ",\"solveHinted\":" << hinted_json << "}\n";
            return 0;
        }
        require(argc == 3, "expected address and mode");
        std::string bearer; require(static_cast<bool>(std::getline(std::cin, bearer)), "missing test credential");
        auto channel = grpc::CreateChannel(argv[1], grpc::InsecureChannelCredentials());
        auto stub = wire::ModelConfigurationService::NewStub(channel);
        const auto configure = [&](grpc::ClientContext& c, const std::string& id) {
            deadline(c); c.AddMetadata("authorization", "Bearer " + bearer);
            if (!id.empty()) { c.AddMetadata(api::model_session_metadata, id); }
        };
        if (std::string_view(argv[2]) == "hold") {
            auto sessions = wire::ModelSessionService::NewStub(channel);
            grpc::ClientContext stream_context; configure(stream_context, "");
            wire::OpenModelSessionRequest open; open.set_wire_contract(std::string(api::model_session_contract));
            auto reader = sessions->OpenModelSession(&stream_context, open); wire::ModelSessionOpened session;
            require(reader->Read(&session), "session open failed");
            grpc::ClientContext create_context; configure(create_context, session.session_id());
            wire::CreateModelResponse made;
            ok(stub->CreateModel(&create_context, request(definition(binary_model())), &made));
            grpc::ClientContext solve_context; configure(solve_context, session.session_id());
            wire::SolveModelResponse solved;
            ok(stub->SolveModel(&solve_context, solve_request(made.model_handle(), binary), &solved));
            compare(native_result(solved.result()), direct(binary_model(), binary));

            // The process-disconnect lifecycle regression now also proves that a
            // real authenticated session can carry pt-solve-hints/v1 without
            // persisting them in the model/session after this unary call.
            grpc::ClientContext hinted_context; configure(hinted_context, session.session_id());
            wire::SolveModelResponse hinted;
            ok(stub->SolveModel(
                &hinted_context, binary_hinted_request(made.model_handle(), hints), &hinted));
            auto expected_owner = model_test::create(definition(binary_model()), preset());
            compare(native_result(hinted.result()), expected_owner->solve(binary, hints));

            // Pipe to the test parent only; these credentials are never CI logs.
            std::cout << session.session_id() << '\n' << made.model_handle() << '\n' << std::flush;
            std::string control; (void)std::getline(std::cin, control);
            stream_context.TryCancel(); (void)reader->Finish();
        } else if (std::string_view(argv[2]) == "stale") {
            std::string id, handle; require(static_cast<bool>(std::getline(std::cin, id)), "missing session");
            require(static_cast<bool>(std::getline(std::cin, handle)), "missing handle");
            const auto end = std::chrono::steady_clock::now() + 5s;
            while (true) {
                grpc::ClientContext c; configure(c, id); wire::DescribeModelResponse out;
                const auto status = stub->DescribeModel(&c, handle_request(handle), &out);
                if (status.error_code() == grpc::StatusCode::NOT_FOUND) {
                    error(status, grpc::StatusCode::NOT_FOUND, "session.not_found"); break;
                }
                ok(status); require(std::chrono::steady_clock::now() < end, "dead client session retained");
                std::this_thread::sleep_for(10ms);
            }
        } else { throw std::runtime_error("unknown process probe mode"); }
        return 0;
    } catch (const std::exception& e) { std::cerr << e.what() << '\n'; return 1; }
}
