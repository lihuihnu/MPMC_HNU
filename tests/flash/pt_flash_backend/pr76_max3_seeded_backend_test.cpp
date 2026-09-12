#include <mpmc/flash/pr76_pt_flash_backend.hpp>

#include "../pr76_three_phase/synthetic_fixture.hpp"

#include <algorithm>
#include <stdexcept>

namespace {
namespace fl = mpmc::flash;

bool phase_equal(const fl::PtCandidatePhase& a, const fl::PtCandidatePhase& b) {
    return a.mole_phase_fraction == b.mole_phase_fraction &&
           a.composition == b.composition &&
           a.activity.branch == b.activity.branch &&
           a.activity.smooth == b.activity.smooth &&
           a.activity.ln_phi == b.activity.ln_phi &&
           a.compressibility_factor == b.compressibility_factor;
}

} // namespace

int main() {
    try {
        const auto phase_model = pr76_max3_test::model();
        fl::Pr76VleEvaluator evaluator(phase_model);
        const auto feed = pr76_max3_test::equal_feed();
        const auto starts = pr76_max3_test::starts();

        const auto direct = fl::solve_pr76_pt_max3(
            1.0e6, 250.0, feed, evaluator, {}, starts, starts);
        const auto published = fl::project_pr76_pt_max3_phase_set(direct);

        fl::Pr76PtFlashBackendOptions options;
        options.initial_starts = starts;
        options.final_starts = starts;
        fl::Pr76PtFlashBackend backend(evaluator, options);
        fl::PtFlashBackend& runtime = backend;
        const auto adapted = runtime.solve({1.0e6, 250.0, feed});

        if (!adapted.structurally_valid() ||
            !runtime.capability().supports_phase_count(3U) ||
            runtime.capability().maximum_phase_count() != 3U ||
            runtime.capability().transition_capability.edge(2U, 3U) == nullptr ||
            runtime.capability().transition_capability.edge(3U, 2U) == nullptr ||
            adapted.solution.status != published.solution.status ||
            adapted.solution.accepted_phase_count() !=
                published.solution.accepted_phase_count() ||
            adapted.solution.accepted_phase_count() != 3U ||
            adapted.solution.feed != published.solution.feed) {
            throw std::runtime_error("seeded PR76 max3 backend capability/publication mismatch");
        }
        const auto* a = adapted.solution.accepted_phase_set();
        const auto* b = published.solution.accepted_phase_set();
        if (a == nullptr || b == nullptr || a->phases.size() != b->phases.size()) {
            throw std::runtime_error("seeded PR76 max3 backend accepted set missing");
        }
        for (std::size_t i = 0; i < a->phases.size(); ++i) {
            if (!phase_equal(a->phases[i], b->phases[i])) {
                throw std::runtime_error("seeded PR76 max3 backend changed phase payload");
            }
        }
        const bool has_2_to_3 = std::any_of(
            adapted.transition_report.evidence.begin(),
            adapted.transition_report.evidence.end(),
            [](const fl::PtPhaseTransitionEvidence& evidence) {
                return evidence.source_phase_count == 2U &&
                       evidence.target_phase_count == 3U &&
                       evidence.resolution ==
                           fl::PtPhaseTransitionResolution::accepted_target &&
                       evidence.fresh_target_solve_attempted &&
                       evidence.target_topology_closed;
            });
        if (!has_2_to_3) {
            throw std::runtime_error("seeded PR76 max3 backend lost 2->3 transition evidence");
        }
        return 0;
    } catch (const std::exception&) {
        return 1;
    }
}
