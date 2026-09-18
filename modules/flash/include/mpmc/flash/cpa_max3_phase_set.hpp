#ifndef MPMC_FLASH_CPA_MAX3_PHASE_SET_HPP
#define MPMC_FLASH_CPA_MAX3_PHASE_SET_HPP

#include <mpmc/flash/cpa_three_phase.hpp>
#include <mpmc/flash/pt_vle_phase_set.hpp>

#include <string>
#include <vector>

namespace mpmc::flash {

inline constexpr const char* cpa_pt_max3_phase_set_publication_convention =
    "CPA/PT/max3/authoritative-phase-set-publication/v1";

struct CpaPtMax3PhaseSetResult {
    PtPhaseSetResult solution;
    std::string dataset_id;
    std::string revision;
    std::vector<std::string> component_ids;
    std::string model_profile{std::string(thermodynamics::cpa_profile)};
    std::string phase_convention{std::string(thermodynamics::cpa_pt_convention)};
    std::string stability_convention{cpa_pt_stability_convention};
    std::string orchestration_convention{cpa_pt_max3_convention};
    std::string publication_convention{
        cpa_pt_max3_phase_set_publication_convention};
    thermodynamics::CpaPtOptions pt_options;
};

[[nodiscard]] inline PtCandidatePhase project_cpa_three_phase_owned_phase(
    const PtThreePhasePhase& source) {
    return detail::project_max3_owned_phase(source);
}

[[nodiscard]] inline CpaPtMax3PhaseSetResult project_cpa_pt_max3_phase_set(
    const CpaPtMax3Result& source) {
    CpaPtMax3PhaseSetResult result;
    result.dataset_id = source.base.dataset_id;
    result.revision = source.base.revision;
    result.component_ids = source.base.component_ids;
    result.model_profile = source.base.model_profile;
    result.phase_convention = source.base.phase_convention;
    result.stability_convention = source.base.stability_convention;
    result.pt_options = source.base.pt_options;
    result.solution = detail::project_max3_solution(source,
        [](const CpaPtSplitResult& two_phase) {
            return project_pt_vle_phase_set(two_phase.solution);
        }, "CPA max3 publication: three-phase status lacks accepted candidate");
    return result;
}

} // namespace mpmc::flash

#endif // MPMC_FLASH_CPA_MAX3_PHASE_SET_HPP
