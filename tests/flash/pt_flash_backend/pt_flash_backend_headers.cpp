#include <mpmc/flash/pr76_pt_flash_backend.hpp>
#include <mpmc/flash/pt_flash_backend.hpp>
#include <mpmc/flash/sw92_profile_c_pt_flash_backend.hpp>

#include <string_view>

bool pt_flash_backend_headers() {
    mpmc::flash::PtFlashBackendCapability capability;
    mpmc::flash::PtFlashBackendResult result;
    return mpmc::flash::PtFlashBackendCapability::convention ==
               std::string_view{"PT/flash-backend/capability-and-dispatch/v1"} &&
           mpmc::flash::PtFlashBackendResult::convention ==
               std::string_view{"PT/flash-backend/result/v1"} &&
           mpmc::flash::pr76_pt_flash_backend_id ==
               std::string_view{"PR76/PT-VLE/phase-set-backend/v1"} &&
           mpmc::flash::sw92_profile_c_pt_flash_backend_id ==
               std::string_view{"SW92/Profile-C/PT/phase-set-backend/v1"} &&
           !capability.structurally_valid() &&
           !result.structurally_valid() &&
           result.accepted_phase_set() == nullptr;
}
