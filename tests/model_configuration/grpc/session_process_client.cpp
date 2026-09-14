#include "test_support.hpp"
#include <mpmc/model_configuration_grpc/model_sessions.hpp>
#include <google/protobuf/util/json_util.h>
#include <iostream>
#include <thread>

int main(int argc, char** argv) {
    using namespace service_test;
    try {
        if (argc == 2 && std::string_view(argv[1]) == "--fixture-json") {
            // Export the existing attributed native fixture; no duplicate TS data.
            std::string json;
            require(google::protobuf::util::MessageToJsonString(
                request(definition(binary_model())), &json).ok(), "fixture encoding failed");
            std::string solve_json;
            require(google::protobuf::util::MessageToJsonString(
                solve_request("", binary), &solve_json).ok(), "state encoding failed");
            std::cout << "{\"create\":" << json << ",\"solve\":" << solve_json << "}\n"; return 0;
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
