#ifndef MPMC_TEST_CPA_THERMOPACK_PARITY_THRESHOLDS_HPP
#define MPMC_TEST_CPA_THERMOPACK_PARITY_THRESHOLDS_HPP

#include <string_view>

namespace cpa_thermopack_parity_thresholds {

// Frozen numerical-parity envelope for the named ThermoPack d68c794
// methanol/water parity snapshot and the five 333.15 K TP states currently
// carried by this PR. These are cross-implementation numerical gates, not
// experimental-accuracy claims or general CPA error bounds.
inline constexpr std::string_view contract =
    "MPMC_HNU/CPA/ThermoPack-parity-thresholds/v1";

// Five-state TP-flash observables: maximum absolute error over all states.
inline constexpr double flash_max_abs_beta_vapor = 2.0e-7;
inline constexpr double flash_max_abs_liquid_methanol = 6.0e-8;
inline constexpr double flash_max_abs_vapor_methanol = 3.0e-11;

// Ten phase-state root/property observables: maximum error over all states.
inline constexpr double phase_max_relative_density = 1.0e-10;
inline constexpr double phase_max_abs_z = 1.0e-10;
inline constexpr double phase_max_abs_ln_phi = 1.0e-10;

// Common-(T,V,x) decomposition observables. Pressure thresholds are in Pa;
// chemical-potential thresholds apply to the dimensionless residual mu/(R T)
// decomposition used by the audit.
inline constexpr double phase_max_abs_pressure_physical_pa = 5.0e-7;
inline constexpr double phase_max_abs_pressure_association_pa = 5.0e-6;
inline constexpr double phase_max_abs_pressure_total_pa = 5.0e-6;
inline constexpr double phase_max_abs_mu_cubic_over_rt = 5.0e-15;
inline constexpr double phase_max_abs_mu_association_over_rt = 5.0e-11;

} // namespace cpa_thermopack_parity_thresholds

#endif // MPMC_TEST_CPA_THERMOPACK_PARITY_THRESHOLDS_HPP
