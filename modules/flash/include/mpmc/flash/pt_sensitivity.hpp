#ifndef MPMC_FLASH_PT_SENSITIVITY_HPP
#define MPMC_FLASH_PT_SENSITIVITY_HPP

#include <cstddef>
#include <stdexcept>
#include <string>
#include <string_view>
#include <vector>

namespace mpmc::flash {

/// Local first-order sensitivity of an already accepted interior VLE solution.
/// This is NOT a derivative of the flash iterations, TPD search, line search, or
/// phase-selection map. Columns use q=(p_Pa,T_K,z_0,...,z_{N-2}) with
/// z_{N-1}=1-sum(z_0,...,z_{N-2}); rows follow the ordered component snapshot.
enum class PtSensitivityStatus {
    success,
    solution_not_accepted,
    phase_boundary,
    unsupported_feed_support,
    property_failure,
    ill_conditioned_equilibrium,
    arithmetic_failure
};

struct PtSensitivityOptions {
    // sqrt(epsilon(double)) = 2^-26. This is a derivative-domain guard, not a
    // replacement for the flash phase-existence threshold.
    double minimum_derivative_phase_fraction{1.490116119384765625e-8};
    // Optional stricter caller requirement. Zero means use only the internal
    // floating-point invertibility guard.
    double minimum_reciprocal_condition{0.0};
    std::size_t max_components{256};
};

struct PtVleSensitivityResult {
    PtSensitivityStatus status{PtSensitivityStatus::arithmetic_failure};
    static constexpr std::string_view convention =
        "PT-VLE/implicit-logK/reduced-feed-v1";

    std::size_t component_count{};
    std::size_t input_count{}; // N+1: p, T, then N-1 independent feed fractions.
    double pressure_pa{};
    double temperature_k{};
    std::vector<double> feed; // Accepted normalized feed in snapshot order.

    // Row-major matrices: [component * input_count + column].
    std::vector<double> log_k_jacobian;
    std::vector<double> liquid_jacobian;
    std::vector<double> vapor_jacobian;
    std::vector<double> vapor_fraction_gradient;

    // Infinity-norm reciprocal condition estimate of dF/dlogK, based on the
    // computed LU inverse. It is a local numerical diagnostic, not a proof of
    // global smoothness or a forward-error bound.
    double equilibrium_jacobian_rcond{};
    // max |A*X+B| / max(1, |B|+sum|A_ij*X_j|) over all sensitivity columns.
    double linear_solve_backward_error{};
    // max |F| at the stored converged state when re-evaluated by the local
    // sensitivity equation path. This detects inconsistent/tampered result data.
    double equilibrium_residual_norm{};

    std::string diagnostic;

    [[nodiscard]] static constexpr std::size_t pressure_column() noexcept {
        return 0;
    }
    [[nodiscard]] static constexpr std::size_t temperature_column() noexcept {
        return 1;
    }
    [[nodiscard]] std::size_t dependent_feed_component() const {
        if (component_count < 2) {
            throw std::logic_error("PT sensitivity: no reduced feed chart");
        }
        return component_count - 1;
    }
    [[nodiscard]] std::size_t feed_column(std::size_t independent_component) const {
        if (component_count < 2 || independent_component + 1 >= component_count) {
            throw std::out_of_range("PT sensitivity: feed coordinate outside reduced chart");
        }
        return 2 + independent_component;
    }
    [[nodiscard]] double d_log_k(std::size_t component, std::size_t column) const {
        return log_k_jacobian.at(component * input_count + column);
    }
    [[nodiscard]] double d_liquid(std::size_t component, std::size_t column) const {
        return liquid_jacobian.at(component * input_count + column);
    }
    [[nodiscard]] double d_vapor(std::size_t component, std::size_t column) const {
        return vapor_jacobian.at(component * input_count + column);
    }
    [[nodiscard]] double d_vapor_fraction(std::size_t column) const {
        return vapor_fraction_gradient.at(column);
    }
};

} // namespace mpmc::flash

#endif // MPMC_FLASH_PT_SENSITIVITY_HPP
