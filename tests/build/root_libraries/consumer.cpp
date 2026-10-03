#include <mpmc/ad/dual.hpp>
#if MPMC_TEST_THERMODYNAMICS
#include <mpmc/thermodynamics/pr76_pure.hpp>
#endif
#if MPMC_TEST_FLASH
#include <mpmc/flash/rachford_rice.hpp>
#include <array>
#endif

// A build/usage-requirement probe, not new physical or numerical evidence.
int main() {
    using Number = mpmc::ad::Dual<double, 1>;
    const auto value = Number::variable(2.0, 0);
    if ((value * value).derivative(0) != 4.0) {
        return 1;
    }
#if MPMC_TEST_THERMODYNAMICS
    const mpmc::thermodynamics::Pr76PureValues<Number> coefficients{value, Number{1.0}};
    if (coefficients.a.value() != 2.0 || coefficients.b.derivative(0) != 0.0) {
        return 2;
    }
#endif
#if MPMC_TEST_FLASH
    // Equal K=1 is an algebraic degeneracy, not a phase equilibrium claim.
    const std::array feed{0.5, 0.5};
    const std::array log_k{0.0, 0.0};
    const auto split = mpmc::flash::solve_rachford_rice(feed, log_k);
    if (split.status != mpmc::flash::RachfordRiceStatus::degenerate) {
        return 3;
    }
#endif
    return 0;
}
