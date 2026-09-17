#ifndef MPMC_TEST_CPA_ASSOCIATION_CONTINUATION_CONTRACT_HPP
#define MPMC_TEST_CPA_ASSOCIATION_CONTINUATION_CONTRACT_HPP

#include "../cpa_physical_validation/cpa_thermopack_parity_thresholds.hpp"

#include <string_view>

namespace cpa_association_continuation_contract {

inline constexpr std::string_view contract =
    "MPMC_HNU/CPA/association-continuation-safety/v1";

// Production association convergence remains unchanged.
inline constexpr double production_site_fraction_tolerance = 1.0e-12;

// Audit-only equivalence guard for the internal association state.
inline constexpr double site_fraction_equivalence_guard = 1.0e-10;

// Observable acceptance scales are inherited from the already-frozen
// ThermoPack parity v1 phase-kernel contract. These are numerical-consistency
// scales, not new physical tolerances.
inline constexpr double pressure_error_guard_pa =
    cpa_thermopack_parity_thresholds::phase_max_abs_pressure_total_pa;
inline constexpr double ln_phi_error_guard =
    cpa_thermopack_parity_thresholds::phase_max_abs_ln_phi;

inline constexpr std::string_view cache_scope = "one_CpaPtPhase_roots_invocation";
inline constexpr std::string_view seed_policy = "nearest_accepted_density";
inline constexpr std::string_view fallback_policy = "deterministic_current_cold_solve";
inline constexpr std::string_view acceptance_policy =
    "observable_certificate_required_before_warm_accept";

} // namespace cpa_association_continuation_contract

#endif // MPMC_TEST_CPA_ASSOCIATION_CONTINUATION_CONTRACT_HPP
