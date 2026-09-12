#include <mpmc/flash/pr76_max3_phase_set.hpp>
#include <mpmc/flash/pr76_pt_flash_backend.hpp>
#include <mpmc/flash/pr76_three_phase.hpp>
#include <mpmc/flash/pt_flash_backend.hpp>
#include <mpmc/flash/pt_phase_transition.hpp>
#include <mpmc/flash/pt_three_phase.hpp>
#include <mpmc/flash/sw92_profile_c_pt_flash_backend.hpp>

#include <string_view>

bool pt_flash_backend_headers() {
    mpmc::flash::PtFlashBackendCapability capability;
    mpmc::flash::PtFlashBackendResult result;
    mpmc::flash::PtPhaseTransitionCapability transitions;
    mpmc::flash::PtPhaseTransitionReport report;
    mpmc::flash::PtThreePhaseResult three_phase;
    return mpmc::flash::PtFlashBackendCapability::convention ==
               std::string_view{"PT/flash-backend/capability-and-dispatch/v1"} &&
           mpmc::flash::PtFlashBackendResult::convention ==
               std::string_view{"PT/flash-backend/result/v1"} &&
           mpmc::flash::PtPhaseTransitionCapability::convention ==
               std::string_view{"PT/phase-transition-boundary/v1"} &&
           mpmc::flash::PtPhaseTransitionReport::convention ==
               std::string_view{"PT/phase-transition-boundary/v1"} &&
           std::string_view{mpmc::flash::PtThreePhaseResult::convention} ==
               std::string_view{"PT/three-phase/logK-SSI/generalized-RR/v1"} &&
           std::string_view{mpmc::flash::Pr76PtMax3Result::convention} ==
               std::string_view{"PR76/PT/max3/TPD-VLE-additional-phase-generalized-RR/v1"} &&
           mpmc::flash::pr76_pt_flash_backend_id ==
               std::string_view{"PR76/PT/max3/phase-set-backend/v2"} &&
           mpmc::flash::sw92_profile_c_pt_flash_backend_id ==
               std::string_view{"SW92/Profile-C/PT/phase-set-backend/v1"} &&
           !capability.structurally_valid() &&
           !result.structurally_valid() &&
           transitions.structurally_valid(1U) &&
           report.structurally_valid(transitions, 1U) &&
           three_phase.candidate() == nullptr &&
           result.accepted_phase_set() == nullptr;
}
