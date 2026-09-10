#ifndef MPMC_FLASH_PR76_PHASE_SET_HPP
#define MPMC_FLASH_PR76_PHASE_SET_HPP

#include <mpmc/flash/pr76_split.hpp>
#include <mpmc/flash/pt_vle_phase_set.hpp>

#include <span>
#include <string>
#include <vector>

namespace mpmc::flash {

// Model metadata remains outside the model-independent phase-set payload.
struct Pr76PtPhaseSetResult {
    PtPhaseSetResult solution;
    std::string dataset_id, revision;
    std::vector<std::string> component_ids;
    std::string model_profile{thermodynamics::pr76_profile};
    std::string phase_convention{thermodynamics::pr76_pt_convention};
    thermodynamics::Pr76RootOptions root_options;
};

// Compatibility projection: preserve every PR76 metadata field while replacing
// the fixed liquid/vapor public result shape with a generic candidate phase set.
[[nodiscard]] inline Pr76PtPhaseSetResult project_pr76_pt_phase_set(
    const Pr76PtSplitResult& source) {
    Pr76PtPhaseSetResult result;
    result.solution = project_pt_vle_phase_set(source.solution);
    result.dataset_id = source.dataset_id;
    result.revision = source.revision;
    result.component_ids = source.component_ids;
    result.model_profile = source.model_profile;
    result.phase_convention = source.phase_convention;
    result.root_options = source.root_options;
    return result;
}

// New generic-result entry point. The established solve_pr76_pt_vle(...) API is
// intentionally retained unchanged for source compatibility and sensitivity use.
// This function runs exactly that validated VLE path once and only projects its
// owned result; it does not introduce a second flash algorithm.
[[nodiscard]] inline Pr76PtPhaseSetResult solve_pr76_pt_phase_set(
    double pressure_pa, double temperature_k, std::span<const double> feed,
    Pr76VleEvaluator& evaluator, PtSplitOptions options = {},
    std::span<const std::vector<double>> initial_starts = {},
    std::span<const std::vector<double>> final_starts = {}) {
    const auto legacy = solve_pr76_pt_vle(
        pressure_pa, temperature_k, feed, evaluator, options,
        initial_starts, final_starts);
    return project_pr76_pt_phase_set(legacy);
}

} // namespace mpmc::flash

#endif // MPMC_FLASH_PR76_PHASE_SET_HPP
