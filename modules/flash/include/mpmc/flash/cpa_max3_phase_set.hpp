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
    PtCandidatePhase phase;
    phase.mole_phase_fraction = source.mole_phase_fraction;
    phase.composition = source.composition;
    phase.activity = source.activity;
    phase.compressibility_factor = source.z;
    return phase;
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

    const auto project_two_phase = [&](const CpaPtSplitResult& two_phase) {
        auto projected = project_pt_vle_phase_set(two_phase.solution);
        projected.capability.maximum_phase_count = 3U;
        return projected;
    };

    switch (source.status) {
    case CpaPtMax3Status::single_phase:
        result.solution = project_two_phase(source.base);
        break;
    case CpaPtMax3Status::two_phase:
        if (const auto* neighbor = source.two_phase_neighbor()) {
            result.solution = project_two_phase(*neighbor);
        } else {
            result.solution = project_two_phase(source.base);
        }
        break;
    case CpaPtMax3Status::three_phase: {
        result.solution.capability.maximum_phase_count = 3U;
        result.solution.pressure_pa = source.base.solution.initial_stability.pressure_pa;
        result.solution.temperature_k = source.base.solution.initial_stability.temperature_k;
        result.solution.feed = source.base.solution.initial_stability.feed;
        result.solution.global_stability_proven = false;
        const auto* candidate = source.three_phase_candidate();
        if (candidate == nullptr) {
            result.solution.status = PtPhaseSetStatus::indeterminate;
            result.solution.diagnostic =
                "CPA max3 publication: three-phase status lacks accepted candidate";
            break;
        }
        PtCandidatePhaseSet set;
        set.phases.reserve(3U);
        for (const auto& phase : candidate->phases) {
            set.phases.push_back(project_cpa_three_phase_owned_phase(phase));
        }
        result.solution.candidate_phase_set = std::move(set);
        result.solution.status = PtPhaseSetStatus::accepted;
        result.solution.diagnostic = source.diagnostic;
        break;
    }
    case CpaPtMax3Status::higher_phase_count_or_wrong_candidate:
        result.solution = project_two_phase(source.base);
        result.solution.status = PtPhaseSetStatus::phase_set_unstable;
        result.solution.diagnostic = source.diagnostic;
        break;
    case CpaPtMax3Status::phase_boundary_unresolved:
    case CpaPtMax3Status::indeterminate:
        result.solution.capability.maximum_phase_count = 3U;
        result.solution.pressure_pa = source.base.solution.initial_stability.pressure_pa;
        result.solution.temperature_k = source.base.solution.initial_stability.temperature_k;
        result.solution.feed = source.base.solution.initial_stability.feed;
        result.solution.status = PtPhaseSetStatus::indeterminate;
        result.solution.global_stability_proven = false;
        result.solution.diagnostic = source.diagnostic;
        break;
    }
    return result;
}

} // namespace mpmc::flash

#endif // MPMC_FLASH_CPA_MAX3_PHASE_SET_HPP
