#include "test_support.hpp"

#include <google/rpc/status.pb.h>
#include <iostream>

namespace {
using namespace service_test;
using SC = grpc::StatusCode;

mc::PtSolveHints structural_hints() {
    auto hints = mc::make_pt_solve_hints_v1();
    hints.initial_stability_starts = pr76_max3_test::starts();
    hints.final_two_phase_stability_starts = pr76_max3_test::starts();
    mc::PtThreePhaseContinuationHint three;
    three.compositions = pr76_max3_test::reference_phases();
    three.phase_fraction_seed = {1.0 / 3.0, 1.0 / 3.0};
    hints.three_phase_continuation_starts.push_back(std::move(three));
    return hints;
}

fl::Pr76PtFlashBackendOptions native_options(const mc::PtSolveHints& hints) {
    fl::Pr76PtFlashBackendOptions options;
    options.initial_starts = hints.initial_stability_starts;
    options.final_starts = hints.final_two_phase_stability_starts;
    for (const auto& public_hint : hints.three_phase_continuation_starts) {
        fl::Pr76PtThreePhaseStart native;
        native.compositions = public_hint.compositions;
        native.phase_fraction_seed = public_hint.phase_fraction_seed;
        options.three_phase_starts.push_back(std::move(native));
    }
    return options;
}

wire::SolveModelRequest hinted_request(
    const std::string& handle, const fl::PtFlashRequest& state,
    const mc::PtSolveHints& hints) {
    auto request = solve_request(handle, state);
    api::encode_solve_hints(hints, *request.mutable_hints());
    return request;
}

std::string detail_field(const grpc::Status& status) {
    google::rpc::Status envelope;
    wire::ModelServiceError detail;
    require(envelope.ParseFromString(status.error_details()) &&
            envelope.details_size() == 1 &&
            envelope.details(0).UnpackTo(&detail),
            "missing structured hint error detail");
    return detail.field();
}

void run() {
    mc::Pr76ModelRegistry registry({}, synthetic);
    Server server(registry);
    const auto native = pr76_max3_test::model();
    const auto made = create(server);
    const auto hints = structural_hints();
    const auto expected = direct(native, ternary, {}, native_options(hints));

    wire::SolveModelResponse response;
    ok(solve(server, hinted_request(made.model_handle(), ternary, hints), response));
    compare(native_result(response.result()), expected);
    require(response.result().outcome() == old::PT_COMPUTATION_OUTCOME_ACCEPTED &&
            response.result().candidate_phase_set().phases_size() == 3,
            "wire hint solve did not preserve established three-phase closure");

    // An omitted hints field remains exactly the old cold request.
    compare(native_result(solve(server, made.model_handle(), ternary).result()),
            direct(native, ternary));

    auto bad_fraction = hints;
    bad_fraction.three_phase_continuation_starts[0].phase_fraction_seed = {0.8, 0.8};
    response.Clear();
    auto status = solve(server,
        hinted_request(made.model_handle(), ternary, bad_fraction), response);
    error(status, SC::INVALID_ARGUMENT, "request.rejected");
    require(detail_field(status) ==
                "hints.three_phase_continuation_starts[0].phase_fraction_seed",
            "state-dependent hint field location changed across wire");

    auto bad_version = hints;
    bad_version.version = "pt-solve-hints/future";
    response.Clear();
    status = solve(server, hinted_request(made.model_handle(), ternary, bad_version), response);
    error(status, SC::UNIMPLEMENTED, "configuration.unsupported_version");
    require(detail_field(status) == "hints.version",
            "hint version field location changed across wire");

    auto malformed = hinted_request(made.model_handle(), ternary, hints);
    malformed.mutable_hints()->mutable_three_phase_continuation_starts(0)->mutable_compositions()->RemoveLast();
    response.Clear();
    status = solve(server, malformed, response);
    error(status, SC::INVALID_ARGUMENT, "configuration.invalid_value");
    require(detail_field(status) ==
                "hints.three_phase_continuation_starts[0].compositions",
            "wire-only fixed-shape rejection lost its field");

    // Rejected hinted calls release admission and do not alter the model.
    response.Clear();
    ok(solve(server, hinted_request(made.model_handle(), ternary, hints), response));
    compare(native_result(response.result()), expected);
    ok(release(server, made.model_handle()));
    require(registry.status().resident_models == 0U,
            "wire hinted solve leaked registry capacity");
}
} // namespace

int main() {
    try {
        run();
        std::cout << "PASS solve_hints_parity\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "FAIL: " << error.what() << '\n';
        return 1;
    }
}
