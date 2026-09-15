#include "test_support.hpp"
#include "../../flash/pr76_three_phase/sour_gas_fixture.hpp"

#include <cstddef>
#include <iostream>
#include <string>
#include <utility>

namespace {
namespace mc = mpmc::model_configuration;
namespace fl = mpmc::flash;
namespace sg = pr76_sour_gas_test;
using model_test::compare;
using model_test::definition;
using model_test::direct;
using model_test::preset;
using model_test::require;

mc::PtSolveHints public_hints() {
    mc::PtSolveHints hints = mc::make_pt_solve_hints_v1();
    hints.initial_stability_starts = sg::starts();
    hints.final_two_phase_stability_starts = sg::starts();
    const auto native_start = sg::continuation_start();
    mc::PtThreePhaseContinuationHint start;
    start.compositions = native_start.compositions;
    start.phase_fraction_seed = native_start.phase_fraction_seed;
    hints.three_phase_continuation_starts.push_back(std::move(start));
    return hints;
}

fl::Pr76PtFlashBackendOptions native_options() {
    fl::Pr76PtFlashBackendOptions options;
    options.initial_starts = sg::starts();
    options.final_starts = sg::starts();
    options.three_phase_starts.push_back(sg::continuation_start());
    return options;
}

} // namespace

int main() {
    try {
        const auto native = sg::model();
        const auto draft = definition(native);
        // Ordinary product policy must accept this literature-attributed dynamic
        // definition. Synthetic-test policy is intentionally not enabled here.
        auto configured = mc::make_pr76_executable_model(draft, preset());
        const fl::PtFlashRequest request{
            pr76_sour_gas_reference::pressure_pa,
            pr76_sour_gas_reference::temperature_k,
            sg::feed()};

        const auto hints = public_hints();
        const auto actual = configured->solve(request, hints);
        const auto expected = direct(native, request, {}, native_options());
        compare(actual, expected);

        require(actual.solution.status == fl::PtPhaseSetStatus::accepted,
                "configured literature state was not accepted");
        require(actual.solution.accepted_phase_count() == 3U,
                "configured literature state did not publish three phases");
        require(!actual.solution.global_stability_proven,
                "configured path strengthened finite stability evidence");
        require(actual.transition_report.evidence.size() ==
                    expected.transition_report.evidence.size(),
                "configured path changed transition evidence");

        const auto& snapshot = configured->parameter_snapshot().definition();
        require(snapshot.dataset_id == "Li-Firoozabadi-2012-acid-gas-PR76",
                "configured model lost literature dataset identity");
        require(snapshot.components.size() == 6U,
                "configured model lost runtime component count");
        require(snapshot.components.size() == actual.capability.component_ids.size(),
                "configured capability/component snapshot diverged");
        for (std::size_t i = 0; i < snapshot.components.size(); ++i) {
            require(snapshot.components[i].component_id == actual.capability.component_ids[i],
                    "configured model changed ordered component identity");
        }

        std::cout << "PASS literature_three_phase_configured\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "FAIL: " << error.what() << '\n';
        return 1;
    }
}
