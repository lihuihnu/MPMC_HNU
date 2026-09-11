#include <mpmc/flash/sw92_profile_c_sensitivity.hpp>

bool sw92_profile_c_sensitivity_header() {
    mpmc::flash::Sw92ProfileCPhaseSetSensitivityResult result;
    return result.convention ==
        mpmc::flash::sw92_profile_c_sensitivity_convention;
}
