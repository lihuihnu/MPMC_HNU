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
    PtCandidatePhase phase;
    phase.mole_phase_fraction = source.mole_phase_fraction;
    phase.composition = source.composition;
    phase.activity = source.activity;
    phase.compressibility_factor = source.z;
    return phase;
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

    const auto project_two_phase = [&](const Pr76PtSplitResult& two_phase) {
        auto projected = project_pr76_pt_phase_set(two_phase);
        projected.solution.capability.maximum_phase_count = 3U;
        return projected.solution;
    };

    switch (source.status) {
    case Pr76PtMax3Status::single_phase:
        result.solution = project_two_phase(source.base);
        break;
    case Pr76PtMax3Status::two_phase:
        if (const auto* neighbor = source.two_phase_neighbor()) {
            result.solution = project_two_phase(*neighbor);
        } else {
            result.solution = project_two_phase(source.base);
        }
        break;
    case Pr76PtMax3Status::three_phase: {
        result.solution.capability.maximum_phase_count = 3U;
        result.solution.pressure_pa = source.base.solution.initial_stability.pressure_pa;
        result.solution.temperature_k = source.base.solution.initial_stability.temperature_k;
        result.solution.feed = source.base.solution.initial_stability.feed;
        result.solution.global_stability_proven = false;
        const auto* candidate = source.three_phase_candidate();
        if (candidate == nullptr) {
            result.solution.status = PtPhaseSetStatus::indeterminate;
            result.solution.diagnostic =
                "PR76 max3 publication: three-phase status lacks accepted candidate";
            break;
        }
        PtCandidatePhaseSet set;
        set.phases.reserve(3U);
        for (const auto& phase : candidate->phases) {
            set.phases.push_back(project_pr76_three_phase_owned_phase(phase));
        }
        result.solution.candidate_phase_set = std::move(set);
        result.solution.status = PtPhaseSetStatus::accepted;
        result.solution.diagnostic = source.diagnostic;
        break;
    }
    case Pr76PtMax3Status::higher_phase_count_or_wrong_candidate:
        result.solution = project_two_phase(source.base);
        result.solution.status = PtPhaseSetStatus::phase_set_unstable;
        result.solution.diagnostic = source.diagnostic;
        break;
    case Pr76PtMax3Status::phase_boundary_unresolved:
    case Pr76PtMax3Status::indeterminate:
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

#endif // MPMC_FLASH_PR76_MAX3_PHASE_SET_HPP
