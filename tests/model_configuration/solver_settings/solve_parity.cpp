#include "test_support.hpp"
#include "../../flash/pr76_three_phase/synthetic_fixture.hpp"

namespace solver_test {
namespace {
// Reuse the existing explicitly synthetic PR76 fixture. This is configuration
// parity with the same solver, not a new independent physical validation.
fl::PtFlashBackendResult solve(th::Pr76RootOptions roots, fl::Pr76PtFlashBackendOptions options,
                              bool seeded) {
    const auto model = pr76_max3_test::model();
    fl::Pr76VleEvaluator evaluator(model, roots); // Fresh owner/workspace per call.
    auto feed = std::vector<double>{1.0, 0.0, 0.0};
    if (seeded) {
        feed = pr76_max3_test::equal_feed();
        // Hints are supplied independently and identically to both C++ paths.
        // A public hint protocol/model registry is not part of this increment.
        options.initial_starts = pr76_max3_test::starts();
        options.final_starts = pr76_max3_test::starts();
        const auto& phases = pr76_max3_test::reference_phases();
        options.three_phase_starts = {{phases, {1.0 / 3.0, 1.0 / 3.0}}};
    }
    fl::Pr76PtFlashBackend backend(evaluator, options);
    return backend.solve({1e6, 250.0, feed});
}
void compare(const fl::PtFlashBackendResult& a, const fl::PtFlashBackendResult& b) {
    require(a.structurally_valid() && b.structurally_valid(), "invalid parity result");
    const auto& x = a.solution;
    const auto& y = b.solution;
    require(x.status == y.status && x.accepted_phase_count() == y.accepted_phase_count() &&
            x.diagnostic == y.diagnostic && x.feed == y.feed &&
            x.pressure_pa == y.pressure_pa && x.temperature_k == y.temperature_k,
            "outcome, count, diagnostic or state differs");
    require(!x.global_stability_proven && !y.global_stability_proven &&
            !a.morphology_resolved && !b.morphology_resolved, "scientific claim changed");
    require(a.capability.dataset_id == b.capability.dataset_id &&
            a.capability.revision == b.capability.revision &&
            a.capability.component_ids == b.capability.component_ids &&
            a.provider_result_convention == b.provider_result_convention, "result provenance differs");
    const auto* left = x.accepted_phase_set();
    const auto* right = y.accepted_phase_set();
    require((left == nullptr) == (right == nullptr), "accepted presence differs");
    if (left) {
        for (std::size_t i = 0; i < left->phases.size(); ++i) {
            const auto& p = left->phases[i];
            const auto& q = right->phases[i];
            // Exact same-platform parity is stronger than relaxing an existing tolerance.
            require(p.mole_phase_fraction == q.mole_phase_fraction && p.composition == q.composition &&
                    p.compressibility_factor == q.compressibility_factor &&
                    p.activity.ln_phi == q.activity.ln_phi && p.activity.branch == q.activity.branch &&
                    p.activity.smooth == q.activity.smooth, "phase payload differs");
        }
    }
}
} // namespace

void solve_parity(std::string_view name) {
    auto settings = custom();
    th::Pr76RootOptions direct_roots;
    fl::Pr76PtFlashBackendOptions direct_options;
    const bool seeded = name != "single_phase_parity" && name != "root_budget_parity" &&
                        name != "initial_budget_parity";
    bool exhausted = false;
    if (name == "root_budget_parity") {
        settings.eos_root.max_iterations = 1;
        direct_roots.max_iterations = 1;
        exhausted = true;
    } else if (name == "initial_budget_parity") {
        settings.initial_stability.max_property_evaluations = 1;
        direct_options.split.initial_stability.max_evaluations = 1;
        exhausted = true;
    } else if (name == "final_two_budget_parity") {
        settings.final_two_phase_stability.max_property_evaluations = 1;
        direct_options.split.final_stability.max_evaluations = 1;
        exhausted = true;
    } else if (name == "final_three_budget_parity") {
        settings.final_three_phase_stability.max_property_evaluations = 1;
        direct_options.final_three_phase_stability.max_evaluations = 1;
        // Keep this regression focused: one attempt per initialization source.
        settings.three_phase.max_three_phase_attempts = 1;
        direct_options.max_three_phase_attempts = 1;
        exhausted = true;
    } else if (name == "single_phase_parity" || name == "three_phase_parity") {
        settings = mc::resolve_pt_solver_preset(mc::mpmc_balanced_default_v1);
    } else {
        throw std::runtime_error("unknown solve parity case");
    }
    const auto prepared = mc::Pr76SolverConfiguration::create(settings);
    const auto direct = solve(direct_roots, direct_options, seeded);
    const auto mapped = solve(prepared.root_options(), prepared.backend_options(), seeded);
    compare(mapped, direct);
    if (exhausted) {
        require(mapped.solution.status == fl::PtPhaseSetStatus::indeterminate &&
                mapped.solution.accepted_phase_count() == 0, "resource-limited result was promoted");
        require(!mapped.solution.diagnostic.empty(), "missing exhaustion diagnostic");
    } else {
        require(mapped.solution.status == fl::PtPhaseSetStatus::accepted &&
                mapped.solution.accepted_phase_count() == (seeded ? 3U : 1U), "baseline did not close");
    }
    require(prepared.settings() == settings, "solve mutated the prepared configuration");
}
} // namespace solver_test
