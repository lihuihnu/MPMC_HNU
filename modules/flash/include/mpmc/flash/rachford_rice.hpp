#ifndef MPMC_FLASH_RACHFORD_RICE_HPP
#define MPMC_FLASH_RACHFORD_RICE_HPP

#include <mpmc/flash/pt_stability.hpp>
#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <numeric>
#include <span>
#include <stdexcept>
#include <vector>

namespace mpmc::flash {

// A fixed-K algebraic result, NEVER a thermodynamic phase-count decision.
enum class RachfordRiceStatus {
    interior, no_resolved_interior_root, degenerate, iteration_limit, unrepresentable
};
struct RachfordRiceOptions {
    int max_iterations{192};
    std::size_t max_components{256};
};
struct RachfordRiceResult {
    RachfordRiceStatus status{RachfordRiceStatus::unrepresentable};
    double vapor_fraction{std::numeric_limits<double>::quiet_NaN()};
    std::vector<double> liquid, vapor;
    double residual{std::numeric_limits<double>::quiet_NaN()};
    double raw_liquid_sum{}, raw_vapor_sum{};
    double mass_absolute{}, mass_relative{};
    int iterations{};
};

namespace detail {
inline void rr_check_options(const RachfordRiceOptions& o) {
    if (o.max_iterations < 0 || o.max_components == 0) {
        throw std::invalid_argument("RR: invalid resource limit");
    }
}
// log(sum z_i*exp(sign*logK_i)), using only the active feed support.
// exp(+logK) is never formed: extreme K values do not overflow the endpoints.
inline double rr_log_sum(std::span<const double> z, std::span<const double> log_k, double sign) {
    double largest = -std::numeric_limits<double>::infinity();
    for (std::size_t i = 0; i < z.size(); ++i) {
        if (z[i] > 0) { largest = std::max(largest, std::log(z[i]) + sign * log_k[i]); }
    }
    double sum = 0, correction = 0;
    for (std::size_t i = 0; i < z.size(); ++i) {
        if (z[i] > 0) { stability_add(std::exp(std::log(z[i]) + sign * log_k[i] - largest), sum, correction); }
    }
    return largest + std::log(sum);
}
inline double rr_value(double beta, std::span<const double> z, std::span<const double> log_k) {
    double value = 0, correction = 0;
    for (std::size_t i = 0; i < z.size(); ++i) {
        if (z[i] == 0) { continue; }
        const double e = std::exp(-std::abs(log_k[i]));
        const double denominator = log_k[i] >= 0 ? beta + (1-beta)*e : (1-beta) + beta*e;
        const double numerator = log_k[i] >= 0 ? -std::expm1(-log_k[i]) : std::expm1(log_k[i]);
        stability_add((z[i] / denominator) * numerator, value, correction);
    }
    return value;
}
} // namespace detail

// Gernert et al. 2014 Eqs.(11)-(14): K=y/x, F=sum z*(K-1)/(1-beta+beta*K).
// Restrict to 0<beta<1; no negative flash and no forced endpoint phase. Endpoint
// signs within 64*eps are unresolved. All logK entries must be finite, including
// unused zero-feed entries. Tiny positive compositions are never floored.
[[nodiscard]] inline RachfordRiceResult solve_rachford_rice(
    std::span<const double> feed, std::span<const double> log_k, RachfordRiceOptions options = {}) {
    detail::rr_check_options(options);
    if (feed.empty() || feed.size() > options.max_components) {
        throw std::length_error("RR: component quota exceeded or empty feed");
    }
    if (feed.size() != log_k.size()) { throw std::invalid_argument("RR: dimension mismatch"); }
    (void)detail::stability_check_composition(feed);
    double contrast = 0;
    for (std::size_t i = 0; i < feed.size(); ++i) {
        if (!std::isfinite(log_k[i])) { throw std::domain_error("RR: finite log K required"); }
        if (feed[i] > 0) { contrast = std::max(contrast, std::abs(log_k[i])); }
    }
    RachfordRiceResult result;
    if (contrast == 0) { result.status = RachfordRiceStatus::degenerate; return result; }
    const double eps = std::numeric_limits<double>::epsilon();
    const double left_sign = detail::rr_log_sum(feed, log_k, 1);
    const double right_sign = detail::rr_log_sum(feed, log_k, -1);
    if (!std::isfinite(left_sign) || !std::isfinite(right_sign)) { return result; }
    if (left_sign <= 64*eps || right_sign <= 64*eps) {
        result.status = RachfordRiceStatus::no_resolved_interior_root;
        return result;
    }
    double left = 0, right = 1, beta = 0.5, value = 0;
    bool resolved = false;
    for (int iteration = 0; iteration < options.max_iterations; ++iteration) {
        result.iterations = iteration + 1;
        beta = std::midpoint(left, right);
        if (!(beta > 0 && beta < 1)) { break; }
        value = detail::rr_value(beta, feed, log_k);
        if (!std::isfinite(value)) { return result; }
        if (value == 0 || beta == left || beta == right ||
            right-left <= 8*eps*std::min(beta, 1-beta)) {
            resolved = true;
            break;
        }
        if (value > 0) { left = beta; } else { right = beta; }
    }
    if (!resolved) { result.status = RachfordRiceStatus::iteration_limit; return result; }
    result.vapor_fraction = beta;
    result.residual = std::abs(value);
    if (result.residual > 64*eps) { return result; }
    result.liquid.resize(feed.size());
    result.vapor.resize(feed.size());
    for (std::size_t i = 0; i < feed.size(); ++i) {
        if (feed[i] == 0) { continue; }
        const double e = std::exp(-std::abs(log_k[i]));
        const double denominator = log_k[i] >= 0 ? beta + (1-beta)*e : (1-beta) + beta*e;
        const double larger = feed[i] / denominator;
        result.liquid[i] = log_k[i] >= 0 ? larger*e : larger;
        result.vapor[i] = log_k[i] >= 0 ? larger : larger*e;
        if (!std::isfinite(larger) || !(result.liquid[i] > 0 && result.vapor[i] > 0)) {
            return result;
        }
    }
    result.raw_liquid_sum = detail::stability_sum(result.liquid);
    result.raw_vapor_sum = detail::stability_sum(result.vapor);
    if (std::abs(result.raw_liquid_sum-1) > 64*eps || std::abs(result.raw_vapor_sum-1) > 64*eps) {
        return result;
    }
    // Explicit roundoff-only normalization for the EOS contract. Recheck material
    // balance on the RETURNED fractions; do not hide the correction in raw RR data.
    for (std::size_t i = 0; i < feed.size(); ++i) {
        result.liquid[i] /= result.raw_liquid_sum;
        result.vapor[i] /= result.raw_vapor_sum;
        const double recovered = std::fma(beta, result.vapor[i], (1-beta)*result.liquid[i]);
        const double error = std::abs(recovered-feed[i]);
        result.mass_absolute = std::max(result.mass_absolute, error);
        if (feed[i] > 0) {
            if (!(result.liquid[i] > 0 && result.vapor[i] > 0)) { return result; }
            result.mass_relative = std::max(result.mass_relative, error/feed[i]);
        } else if (recovered != 0) { return result; }
    }
    (void)detail::stability_check_composition(result.liquid);
    (void)detail::stability_check_composition(result.vapor);
    result.status = RachfordRiceStatus::interior;
    return result;
}

} // namespace mpmc::flash
#endif // MPMC_FLASH_RACHFORD_RICE_HPP
