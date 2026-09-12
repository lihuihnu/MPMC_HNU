#include <mpmc/flash/pr76_max3_phase_set.hpp>
#include <mpmc/flash/pr76_three_phase.hpp>
#include <mpmc/flash/pt_three_phase.hpp>

#include <string_view>

bool pr76_three_phase_headers() {
    mpmc::flash::PtThreePhaseResult generic;
    mpmc::flash::Pr76PtMax3Result pr;
    return std::string_view{mpmc::flash::PtThreePhaseResult::convention} ==
               "PT/three-phase/logK-SSI/generalized-RR/v1" &&
           std::string_view{mpmc::flash::Pr76PtMax3Result::convention} ==
               "PR76/PT/max3/TPD-VLE-additional-phase-generalized-RR/v1" &&
           generic.candidate() == nullptr && pr.three_phase_candidate() == nullptr;
}
