#include <mpmc/flash/cpa_pt_flash_backend.hpp>

#include "../cpa_max3/test_support.hpp"

#include <algorithm>
#include <stdexcept>

namespace {
namespace fl = mpmc::flash;
namespace th = mpmc::thermodynamics;

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
        const auto parameters = cpa_max3_test::parameters();
        const auto model = th::CpaPtPhase::from_parameters(parameters);
        fl::CpaVleEvaluator evaluator(model, cpa_max3_test::fast_pt_options());
        const auto feed = cpa_max3_test::feed_from_phase_fractions(
            {1.0 / 3.0, 1.0 / 3.0, 1.0 / 3.0});
        const auto starts = cpa_max3_test::starts();
        const auto three_start = cpa_max3_test::three_phase_start();

        fl::CpaPtMax3Options direct_options;
        direct_options.two_phase.initial_stability.automatic_starts = false;
        direct_options.two_phase.final_stability.automatic_starts = false;
        direct_options.final_three_phase_stability.automatic_starts = false;
        direct_options.three_phase_starts.push_back(three_start);
        const auto direct = fl::solve_cpa_pt_max3(
            cpa_max3_test::pressure_pa, cpa_max3_test::temperature_k,
            feed, evaluator, direct_options, starts, starts);
        const auto published = fl::project_cpa_pt_max3_phase_set(direct);

        fl::CpaPtFlashBackendOptions options;
        options.split = direct_options.two_phase;
        options.three_phase = direct_options.three_phase;
        options.final_three_phase_stability =
            direct_options.final_three_phase_stability;
        options.three_phase_starts.push_back(three_start);
        options.initial_starts = starts;
        options.final_starts = starts;
        fl::CpaPtFlashBackend backend(evaluator, options);
        fl::PtFlashBackend& runtime = backend;
        const auto adapted = runtime.solve({
            cpa_max3_test::pressure_pa,
            cpa_max3_test::temperature_k,
            feed});

        const auto& capability = runtime.capability();
        if (!capability.structurally_valid() ||
            capability.backend_id != fl::cpa_pt_flash_backend_id ||
            capability.model_profile != th::cpa_profile ||
            capability.algorithm_profile != fl::cpa_pt_max3_convention ||
            capability.publication_profile !=
                fl::cpa_pt_max3_phase_set_publication_convention ||
            capability.configuration_profile !=
                fl::cpa_pt_flash_backend_configuration_profile ||
            capability.supported_phase_counts !=
                std::vector<std::size_t>({1U, 2U, 3U}) ||
            capability.maximum_phase_count() != 3U ||
            !capability.performs_initial_stability_search ||
            !capability.performs_final_phase_set_review ||
            !capability.performs_boundary_neighbor_resolve ||
            capability.global_stability_proven ||
            !capability.phase_metadata_namespace.empty() ||
            capability.transition_capability.edge(1U, 2U) == nullptr ||
            capability.transition_capability.edge(2U, 1U) == nullptr ||
            capability.transition_capability.edge(2U, 3U) == nullptr ||
            capability.transition_capability.edge(3U, 2U) == nullptr ||
            capability.transition_capability.edge(1U, 2U)->support !=
                fl::PtPhaseTransitionSupport::fresh_target_resolve ||
            capability.transition_capability.edge(2U, 1U)->support !=
                fl::PtPhaseTransitionSupport::detection_only ||
            capability.transition_capability.edge(2U, 3U)->support !=
                fl::PtPhaseTransitionSupport::fresh_target_resolve ||
            capability.transition_capability.edge(3U, 2U)->support !=
                fl::PtPhaseTransitionSupport::fresh_target_resolve) {
            throw std::runtime_error("CPA max3 backend capability freeze mismatch");
        }

        if (!adapted.structurally_valid() || adapted.morphology_resolved ||
            adapted.solution.status != published.solution.status ||
            adapted.solution.accepted_phase_count() != 3U ||
            adapted.solution.accepted_phase_count() !=
                published.solution.accepted_phase_count() ||
            adapted.solution.feed != published.solution.feed) {
            throw std::runtime_error("CPA max3 runtime publication mismatch");
        }
        const auto* actual = adapted.solution.accepted_phase_set();
        const auto* expected = published.solution.accepted_phase_set();
        if (actual == nullptr || expected == nullptr ||
            actual->phases.size() != expected->phases.size()) {
            throw std::runtime_error("CPA max3 accepted phase set missing");
        }
        for (std::size_t phase = 0; phase < actual->phases.size(); ++phase) {
            if (!phase_equal(actual->phases[phase], expected->phases[phase])) {
                throw std::runtime_error("CPA max3 backend changed phase payload");
            }
        }
        const bool has_two_to_three = std::any_of(
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
        if (!has_two_to_three) {
            throw std::runtime_error("CPA max3 backend lost accepted 2->3 evidence");
        }
        return 0;
    } catch (const std::exception&) {
        return 1;
    }
}
