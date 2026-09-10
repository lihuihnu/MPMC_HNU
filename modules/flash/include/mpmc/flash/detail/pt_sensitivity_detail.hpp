#ifndef MPMC_FLASH_DETAIL_PT_SENSITIVITY_DETAIL_HPP
#define MPMC_FLASH_DETAIL_PT_SENSITIVITY_DETAIL_HPP

#include <mpmc/flash/pt_sensitivity.hpp>
#include <mpmc/flash/pt_split.hpp>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <span>
#include <stdexcept>
#include <vector>

namespace mpmc::flash::detail {

class PtSensitivityFailure : public std::runtime_error {
public:
    PtSensitivityFailure(PtSensitivityStatus status, const char* message)
        : std::runtime_error(message), status_(status) {}

    [[nodiscard]] PtSensitivityStatus status() const noexcept { return status_; }

private:
    PtSensitivityStatus status_;
};

template <typename Number>
struct PtSensitivityRrState {
    Number vapor_fraction;
    std::vector<Number> liquid;
    std::vector<Number> vapor;
};

// Local IFT of the material-balance/K-value manifold at an accepted state.
// The production RR bisection is neither replayed nor differentiated. With
// K_i=y_i/x_i at the base point,
//
// H = sum_i (y_i-x_i)^2/z_i
// d beta = [sum_i (x_i*y_i/z_i)dlogK_i
//           +sum_i ((y_i-x_i)/z_i)dz_i] / H.
//
// The last x/y component is formed as one minus the first N-1 components so
// all phase-composition AD directions lie exactly in the simplex tangent space.
template <typename Number>
[[nodiscard]] PtSensitivityRrState<Number> pt_sensitivity_rr_linearization(
    std::span<const Number> feed, std::span<const Number> log_k,
    const PtSplitState& base) {
    const std::size_t n = feed.size();
    if (n < 2 || log_k.size() != n || base.fractions.liquid.size() != n ||
        base.fractions.vapor.size() != n) {
        throw PtSensitivityFailure(
            PtSensitivityStatus::arithmetic_failure,
            "PT sensitivity: invalid internal equilibrium dimensions");
    }

    std::vector<double> x(n);
    std::vector<double> y(n);
    double x_sum = 0.0;
    double y_sum = 0.0;
    for (std::size_t i = 0; i + 1 < n; ++i) {
        x[i] = base.fractions.liquid[i];
        y[i] = base.fractions.vapor[i];
        x_sum += x[i];
        y_sum += y[i];
    }
    x[n - 1] = 1.0 - x_sum;
    y[n - 1] = 1.0 - y_sum;

    for (std::size_t i = 0; i < n; ++i) {
        if (!(feed[i].value() > 0.0) || !(x[i] > 0.0) || !(y[i] > 0.0) ||
            !std::isfinite(x[i]) || !std::isfinite(y[i])) {
            throw PtSensitivityFailure(
                PtSensitivityStatus::phase_boundary,
                "PT sensitivity: local RR manifold requires positive compositions");
        }
    }

    double h = 0.0;
    for (std::size_t i = 0; i < n; ++i) {
        const double difference = y[i] - x[i];
        h += difference * difference / feed[i].value();
    }
    if (!std::isfinite(h) || !(h > 0.0)) {
        throw PtSensitivityFailure(
            PtSensitivityStatus::ill_conditioned_equilibrium,
            "PT sensitivity: degenerate RR implicit derivative");
    }

    typename Number::Gradient beta_gradient{};
    for (std::size_t lane = 0; lane < Number::derivative_count; ++lane) {
        double numerator = 0.0;
        for (std::size_t i = 0; i < n; ++i) {
            const double z = feed[i].value();
            numerator += (x[i] * y[i] / z) * log_k[i].derivative(lane) +
                         ((y[i] - x[i]) / z) * feed[i].derivative(lane);
        }
        beta_gradient[lane] = numerator / h;
        if (!std::isfinite(beta_gradient[lane])) {
            throw PtSensitivityFailure(
                PtSensitivityStatus::arithmetic_failure,
                "PT sensitivity: nonfinite phase-fraction derivative");
        }
    }

    PtSensitivityRrState<Number> state{
        Number{base.fractions.vapor_fraction, beta_gradient},
        std::vector<Number>(n), std::vector<Number>(n)};
    Number liquid_sum{0.0};
    Number vapor_sum{0.0};
    for (std::size_t i = 0; i + 1 < n; ++i) {
        typename Number::Gradient dx{};
        typename Number::Gradient dy{};
        const double z = feed[i].value();
        for (std::size_t lane = 0; lane < Number::derivative_count; ++lane) {
            const double d_log_x =
                feed[i].derivative(lane) / z -
                ((y[i] - x[i]) / z) * beta_gradient[lane] -
                base.fractions.vapor_fraction * (y[i] / z) *
                    log_k[i].derivative(lane);
            dx[lane] = x[i] * d_log_x;
            dy[lane] = y[i] * (d_log_x + log_k[i].derivative(lane));
            if (!std::isfinite(dx[lane]) || !std::isfinite(dy[lane])) {
                throw PtSensitivityFailure(
                    PtSensitivityStatus::arithmetic_failure,
                    "PT sensitivity: nonfinite phase-composition derivative");
            }
        }
        state.liquid[i] = Number{x[i], dx};
        state.vapor[i] = Number{y[i], dy};
        liquid_sum += state.liquid[i];
        vapor_sum += state.vapor[i];
    }
    state.liquid[n - 1] = Number{1.0} - liquid_sum;
    state.vapor[n - 1] = Number{1.0} - vapor_sum;
    if (!(state.liquid[n - 1].value() > 0.0) ||
        !(state.vapor[n - 1].value() > 0.0)) {
        throw PtSensitivityFailure(
            PtSensitivityStatus::phase_boundary,
            "PT sensitivity: dependent phase composition is not positive");
    }
    return state;
}

struct PtSensitivityLu {
    std::size_t n{};
    std::vector<double> values;
    std::vector<std::size_t> pivots;
    double norm_inf{};
};

[[nodiscard]] inline bool pt_sensitivity_factor(
    std::span<const double> matrix, std::size_t n, PtSensitivityLu& factor) {
    if (n == 0 || matrix.size() != n * n) {
        return false;
    }
    factor.n = n;
    factor.values.assign(matrix.begin(), matrix.end());
    factor.pivots.resize(n);
    factor.norm_inf = 0.0;
    for (std::size_t row = 0; row < n; ++row) {
        double row_sum = 0.0;
        for (std::size_t column = 0; column < n; ++column) {
            row_sum += std::abs(matrix[row * n + column]);
        }
        factor.norm_inf = std::max(factor.norm_inf, row_sum);
    }
    if (!(factor.norm_inf > 0.0) || !std::isfinite(factor.norm_inf)) {
        return false;
    }

    const double pivot_guard =
        64.0 * std::numeric_limits<double>::epsilon() *
        static_cast<double>(n) * factor.norm_inf;
    for (std::size_t column = 0; column < n; ++column) {
        std::size_t pivot = column;
        double magnitude = std::abs(factor.values[column * n + column]);
        for (std::size_t row = column + 1; row < n; ++row) {
            const double candidate =
                std::abs(factor.values[row * n + column]);
            if (candidate > magnitude) {
                magnitude = candidate;
                pivot = row;
            }
        }
        if (!std::isfinite(magnitude) || !(magnitude > pivot_guard)) {
            return false;
        }
        factor.pivots[column] = pivot;
        if (pivot != column) {
            for (std::size_t j = 0; j < n; ++j) {
                std::swap(factor.values[column * n + j],
                          factor.values[pivot * n + j]);
            }
        }
        const double diagonal = factor.values[column * n + column];
        for (std::size_t row = column + 1; row < n; ++row) {
            const std::size_t offset = row * n;
            const double multiplier = factor.values[offset + column] / diagonal;
            factor.values[offset + column] = multiplier;
            for (std::size_t j = column + 1; j < n; ++j) {
                factor.values[offset + j] -=
                    multiplier * factor.values[column * n + j];
            }
        }
    }
    return true;
}

[[nodiscard]] inline bool pt_sensitivity_solve(
    const PtSensitivityLu& factor, std::vector<double>& rhs) {
    const std::size_t n = factor.n;
    if (rhs.size() != n) {
        return false;
    }
    for (std::size_t column = 0; column < n; ++column) {
        if (factor.pivots[column] != column) {
            std::swap(rhs[column], rhs[factor.pivots[column]]);
        }
    }
    for (std::size_t row = 0; row < n; ++row) {
        for (std::size_t column = 0; column < row; ++column) {
            rhs[row] -= factor.values[row * n + column] * rhs[column];
        }
    }
    for (std::size_t row = n; row-- > 0;) {
        for (std::size_t column = row + 1; column < n; ++column) {
            rhs[row] -= factor.values[row * n + column] * rhs[column];
        }
        rhs[row] /= factor.values[row * n + row];
        if (!std::isfinite(rhs[row])) {
            return false;
        }
    }
    return true;
}

[[nodiscard]] inline double pt_sensitivity_rcond(
    const PtSensitivityLu& factor) {
    const std::size_t n = factor.n;
    std::vector<double> inverse_row_sums(n, 0.0);
    for (std::size_t column = 0; column < n; ++column) {
        std::vector<double> rhs(n, 0.0);
        rhs[column] = 1.0;
        if (!pt_sensitivity_solve(factor, rhs)) {
            return 0.0;
        }
        for (std::size_t row = 0; row < n; ++row) {
            inverse_row_sums[row] += std::abs(rhs[row]);
        }
    }
    const double inverse_norm =
        *std::max_element(inverse_row_sums.begin(), inverse_row_sums.end());
    const double product = factor.norm_inf * inverse_norm;
    if (!(product > 0.0) || !std::isfinite(product)) {
        return 0.0;
    }
    return 1.0 / product;
}

} // namespace mpmc::flash::detail

#endif // MPMC_FLASH_DETAIL_PT_SENSITIVITY_DETAIL_HPP
