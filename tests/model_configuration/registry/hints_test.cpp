#include <mpmc/model_configuration/pr76_model_registry.hpp>
#include "../executable_model/test_support.hpp"

#include <iostream>

namespace {
using namespace model_test;

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

void run() {
    mc::Pr76ModelRegistryLimits limits;
    limits.max_models = 1;
    mc::Pr76ModelRegistry registry(
        limits, mc::ModelDataPolicy::allow_synthetic_tests);
    const auto native = pr76_max3_test::model();
    const auto handle = registry.create(definition(native), preset());
    const auto hints = structural_hints();
    const auto expected = direct(native, ternary, {}, native_options(hints));

    compare(registry.solve(handle, ternary, hints), expected);
    // Hints are per-call only; the next coarse solve remains the native cold path.
    compare(registry.solve(handle, ternary), direct(native, ternary));

    auto bad = hints;
    bad.version = "pt-solve-hints/future";
    bool located = false;
    try { (void)registry.solve(handle, ternary, bad); }
    catch (const mc::ModelConfigurationError& error) {
        located = error.code() == mc::ModelConfigurationErrorCode::unsupported_version &&
                  error.field() == "hints.version";
    }
    require(located, "registry changed public hint version/location semantics");
    compare(registry.solve(handle, ternary, hints), expected);

    auto lease = registry.acquire_solve(handle);
    registry.release(handle);
    const auto retained = registry.status();
    require(retained.registered_models == 0U && retained.resident_models == 1U,
            "released hinted lease did not retain exactly one resident model");
    compare(std::move(lease).solve(ternary, hints), expected);
    const auto released = registry.status();
    require(released.registered_models == 0U && released.resident_models == 0U,
            "hinted retained lease leaked registry capacity");
}
} // namespace

int main() {
    try {
        run();
        std::cout << "PASS solve_hints_lifetime\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "FAIL: " << error.what() << '\n';
        return 1;
    }
}
