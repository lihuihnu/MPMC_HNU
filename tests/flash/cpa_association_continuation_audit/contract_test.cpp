#include <mpmc/thermodynamics/cpa_association.hpp>

#include "contract.hpp"

#include <cmath>
#include <iostream>
#include <stdexcept>

namespace {
namespace th = mpmc::thermodynamics;
namespace contract = cpa_association_continuation_contract;

void require(bool condition, const char* message) {
    if (!condition) { throw std::runtime_error(message); }
}

} // namespace

int main() {
    try {
        constexpr th::CpaAssociationOptions defaults{};
        static_assert(defaults.site_fraction_tolerance ==
                      contract::production_site_fraction_tolerance);

        require(contract::contract ==
                    "MPMC_HNU/CPA/association-continuation-safety/v1",
                "continuation safety contract id drifted");
        require(contract::production_site_fraction_tolerance == 1.0e-12,
                "production association tolerance drifted from frozen v1 contract");
        require(contract::site_fraction_equivalence_guard == 1.0e-10,
                "continuation site-fraction audit guard drifted");
        require(contract::pressure_error_guard_pa == 5.0e-6,
                "continuation pressure guard drifted from ThermoPack parity v1");
        require(contract::ln_phi_error_guard == 1.0e-10,
                "continuation ln(phi) guard drifted from ThermoPack parity v1");
        require(contract::cache_scope == "one_CpaPtPhase_roots_invocation",
                "continuation cache scope drifted");
        require(contract::seed_policy == "nearest_accepted_density",
                "continuation seed policy drifted");
        require(contract::fallback_policy == "deterministic_current_cold_solve",
                "continuation fallback policy drifted");
        require(contract::acceptance_policy ==
                    "observable_certificate_required_before_warm_accept",
                "continuation acceptance policy drifted");

        std::cout
            << "CPA_ASSOC_CONTINUATION_CONTRACT_OK"
            << " contract=" << contract::contract
            << " site_tolerance=" << contract::production_site_fraction_tolerance
            << " site_guard=" << contract::site_fraction_equivalence_guard
            << " pressure_guard_pa=" << contract::pressure_error_guard_pa
            << " lnphi_guard=" << contract::ln_phi_error_guard
            << " cache_scope=" << contract::cache_scope
            << " seed_policy=" << contract::seed_policy
            << " fallback=" << contract::fallback_policy
            << " acceptance=" << contract::acceptance_policy
            << '\n';
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "[FAIL] " << error.what() << '\n';
        return 1;
    }
}
