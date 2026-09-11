#include <mpmc/flash/sw92_profile_c_sensitivity.hpp>

bool sw92_profile_c_sensitivity_header() {
    // Preserve the pre-AD public spelling while the sensitivity path may use
    // Sw92PhaseValues<Dual<...>> internally.
    mpmc::thermodynamics::Sw92PhaseValues<double> legacy_values;
    mpmc::flash::Sw92ProfileCPhaseSetSensitivityResult result;
    return legacy_values.ln_phi.empty() &&
           result.convention ==
               mpmc::flash::sw92_profile_c_sensitivity_convention;
}
