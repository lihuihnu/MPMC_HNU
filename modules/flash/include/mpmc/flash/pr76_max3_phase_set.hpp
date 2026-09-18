#ifndef MPMC_FLASH_PR76_MAX3_PHASE_SET_HPP
#define MPMC_FLASH_PR76_MAX3_PHASE_SET_HPP

#include <mpmc/flash/pr76_phase_set.hpp>
#include <mpmc/flash/pr76_three_phase.hpp>

#include <string>
#include <vector>

namespace mpmc::flash {

inline constexpr const char* pr76_pt_max3_phase_set_publication_convention =
    "PR76/PT/max3/authoritative-phase-set-publication/v1";

struct Pr76PtMax3PhaseSetResult {
    PtPhaseSetResult solution;
    std::string dataset_id;
    std::string revision;
    std::vector<std::string> component_ids;
    std::string model_profile{thermodynamics::pr76_profile};
    std::string phase_convention{thermodynamics::pr76_pt_convention};
    std::string orchestration_convention{pr76_pt_max3_convention};
    std::string publication_convention{
        pr76_pt_max3_phase_set_publication_convention};
    thermodynamics::Pr76RootOptions root_options;
};

[[nodiscard]] inline PtCandidatePhase project_pr76_three_phase_owned_phase(
    const PtThreePhasePhase& source) {
    return detail::project_max3_owned_phase(source);
}

[[nodiscard]] inline Pr76PtMax3PhaseSetResult project_pr76_pt_max3_phase_set(
    const Pr76PtMax3Result& source) {
    Pr76PtMax3PhaseSetResult result;
    result.dataset_id = source.base.dataset_id;
    result.revision = source.base.revision;
    result.component_ids = source.base.component_ids;
    result.model_profile = source.base.model_profile;
    result.phase_convention = source.base.phase_convention;
    result.root_options = source.base.root_options;
    result.solution = detail::project_max3_solution(source,
        [](const Pr76PtSplitResult& two_phase) {
            return project_pr76_pt_phase_set(two_phase).solution;
        }, "PR76 max3 publication: three-phase status lacks accepted candidate");
    return result;
}

} // namespace mpmc::flash

#endif // MPMC_FLASH_PR76_MAX3_PHASE_SET_HPP
