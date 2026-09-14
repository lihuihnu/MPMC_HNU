#include <mpmc/pt_android/portability.h>

#include <mpmc/flash/cpa_pt_flash_backend.hpp>
#include <mpmc/flash/pr76_pt_flash_backend.hpp>
#include <mpmc/flash/sw92_profile_c_pt_flash_backend.hpp>
#include <mpmc/runtime/pt_service.hpp>

#include <type_traits>

namespace {

static_assert(std::is_class_v<mpmc::runtime::PtService>);
static_assert(std::is_base_of_v<mpmc::flash::PtFlashBackend,
                                mpmc::flash::Pr76PtFlashBackend>);
static_assert(std::is_base_of_v<mpmc::flash::PtFlashBackend,
                                mpmc::flash::Sw92ProfileCPtFlashBackend>);
static_assert(std::is_base_of_v<mpmc::flash::PtFlashBackend,
                                mpmc::flash::CpaPtFlashBackend>);

inline constexpr char kConvention[] =
    "MPMC/PT/android-core-portability/v1";

} // namespace

extern "C" uint32_t mpmc_pt_android_core_abi_version(void) {
    return 1U;
}

extern "C" const char* mpmc_pt_android_core_convention(void) {
    return kConvention;
}
