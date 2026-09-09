#ifndef MPMC_THERMODYNAMICS_PR76_ROOTS_HPP
#define MPMC_THERMODYNAMICS_PR76_ROOTS_HPP

#include <algorithm>
#include <array>
#include <cmath>
#include <concepts>
#include <cstddef>
#include <limits>
#include <numeric>
#include <stdexcept>

namespace mpmc::thermodynamics {

/// This is root resolution, NOT a thermodynamic stability or phase-count result.
enum class Pr76RootStatus { success, near_multiple, iteration_limit, unrepresentable };

template <std::floating_point T>
struct Pr76Root {
    T z{};
    T free_volume_z{}; // Z-B, retained separately to avoid subtractive cancellation.
    T scaled_residual{}; // |H| / sum of absolute polynomial terms.
    T pressure_residual{}; // Relative residual of the original rational EOS.
    T slope_ratio{}; // |H'| / sum of absolute derivative terms; not an error bound.
    bool derivative_valid{false};
    int slope_sign{}; // sign(H'); positive corresponds to negative (dp/dv)_T,w.
    int iterations{};
};

template <std::floating_point T>
struct Pr76RootSet {
    std::array<Pr76Root<T>, 3> roots{}; // Increasing Z; only Z>B is admissible.
    std::size_t count{};
    T scale{1};
    Pr76RootStatus status{Pr76RootStatus::unrepresentable};
};

struct Pr76RootOptions {
    int max_iterations{2048}; // Per monotone interval, not an accuracy tolerance.
};

namespace detail {
template <typename Number>
struct Pr76CubicCoefficients { Number c2, c1, c0; };

// Substitute X=Z-B>0, then y=X/s. PR76 Eq.(5) becomes
// y^3 + (4B-1)/s*y^2 + (A+2B^2-4B)/s^2*y - 2B^2/s^3 = 0.
// The positive-domain endpoint has an exact negative constant, unlike F(B).
// Scaling s is held fixed during LOCAL implicit differentiation.
template <typename Number, std::floating_point T>
[[nodiscard]] Pr76CubicCoefficients<Number> pr76_cubic_coefficients(
    const Number& a, const Number& b, T scale) {
    const Number bs = b / scale;
    const T inverse = T{1} / scale;
    return {T{4} * bs - inverse,
            (a / scale) / scale + T{2} * bs * bs - T{4} * bs * inverse,
            ((T{-2} * bs) * bs) * inverse};
}

template <std::floating_point T>
[[nodiscard]] T pr76_poly(T y, const Pr76CubicCoefficients<T>& c) {
    return std::fma(std::fma(y + c.c2, y, c.c1), y, c.c0);
}
template <std::floating_point T>
[[nodiscard]] T pr76_poly_scale(T y, const Pr76CubicCoefficients<T>& c) {
    return ((std::abs(y) + std::abs(c.c2)) * std::abs(y) + std::abs(c.c1)) *
               std::abs(y) + std::abs(c.c0);
}
} // namespace detail

/// Resolve all positive-free-volume roots of the PR76 cubic for finite A and B>0.
/// A may be zero/negative: algebraic continuation is NOT physical validation.
/// A stationary-point sign within 64*epsilon rounding units is unresolved;
/// return near_multiple and NO guessed roots. A simple but poorly conditioned
/// root can have a usable primal with derivative_valid=false.
/// No clipping, complex-part threshold, unguarded Cardano cancellation or phase choice.
template <std::floating_point T>
[[nodiscard]] Pr76RootSet<T> pr76_roots(T a, T b, Pr76RootOptions options = {}) {
    if (!std::isfinite(a) || !std::isfinite(b) || !(b > T{0}) ||
        options.max_iterations <= 0) {
        throw std::invalid_argument("pr76_roots: finite A, B>0 and positive iteration limit required");
    }
    Pr76RootSet<T> result;
    const T eps = std::numeric_limits<T>::epsilon();
    result.scale = std::max({T{1}, b, std::sqrt(std::abs(a))});
    const auto c = detail::pr76_cubic_coefficients(a, b, result.scale);
    if (!std::isfinite(c.c2) || !std::isfinite(c.c1) || !std::isfinite(c.c0) ||
        !(c.c0 < T{0})) {
        return result; // Includes loss of the nonzero endpoint to underflow.
    }
    const T upper = T{1} + std::max({std::abs(c.c2), std::abs(c.c1), std::abs(c.c0)});
    std::array<T, 4> endpoints{};
    std::size_t endpoint_count = 1;
    const T discriminant = std::fma(c.c2, c.c2, T{-3} * c.c1);
    if (discriminant > T{0}) {
        // Stable quadratic roots for H'=0. Product of the two roots is c1/3.
        const T q = -c.c2 - std::copysign(std::sqrt(discriminant), c.c2);
        std::array<T, 2> stationary{q / T{3}, c.c1 / q};
        std::sort(stationary.begin(), stationary.end());
        for (const T point : stationary) {
            if (point > T{0} && point < upper) {
                const T value = detail::pr76_poly(point, c);
                const T magnitude = detail::pr76_poly_scale(point, c);
                if (!std::isfinite(value) || !(magnitude > T{0})) { return result; }
                if (std::abs(value) <= T{64} * eps * magnitude) {
                    result.status = Pr76RootStatus::near_multiple;
                    return result;
                }
                endpoints[endpoint_count++] = point;
            }
        }
    }
    endpoints[endpoint_count++] = upper;
    for (std::size_t interval = 1; interval < endpoint_count; ++interval) {
        T left = endpoints[interval - 1], right = endpoints[interval];
        T fl = detail::pr76_poly(left, c);
        const T fr = detail::pr76_poly(right, c);
        if (!std::isfinite(fl) || !std::isfinite(fr)) { return result; }
        if (std::signbit(fl) == std::signbit(fr)) { continue; }
        T y{};
        bool converged = false;
        int iterations = 0;
        for (; iterations < options.max_iterations; ++iterations) {
            y = std::midpoint(left, right);
            const T fy = detail::pr76_poly(y, c);
            if (fy == T{0} || y == left || y == right ||
                right - left <= T{4} * eps * std::max(y, std::numeric_limits<T>::min())) {
                converged = true;
                break;
            }
            if (std::signbit(fy) == std::signbit(fl)) {
                left = y;
                fl = fy;
            } else {
                right = y;
            }
        }
        if (!converged) {
            result.count = 0;
            result.status = Pr76RootStatus::iteration_limit;
            return result;
        }
        const T x = y * result.scale;
        const T z = x + b;
        const T magnitude = detail::pr76_poly_scale(y, c);
        const T residual = std::abs(detail::pr76_poly(y, c)) / magnitude;
        const T dh = std::fma(T{3} * y + T{2} * c.c2, y, c.c1);
        const T dh_scale = T{3} * y * y + T{2} * std::abs(c.c2) * y + std::abs(c.c1);
        const T slope_ratio = dh_scale > T{0} ? std::abs(dh) / dh_scale : T{0};
        const T bz = b / z;
        const T repulsion = T{1} / x;
        const T attraction = ((a / z) / z) / (T{1} + T{2} * bz - bz * bz);
        const T pressure_residual = std::abs(repulsion - attraction - T{1}) /
                                    (T{1} + std::abs(repulsion) + std::abs(attraction));
        if (!std::isfinite(z) || !(z > b) || !std::isfinite(residual) ||
            !std::isfinite(pressure_residual) || residual > T{128} * eps ||
            pressure_residual > T{512} * eps || result.count == result.roots.size()) {
            result.count = 0;
            result.status = Pr76RootStatus::unrepresentable;
            return result;
        }
        if (result.count != 0 && !(z > result.roots[result.count - 1].z)) {
            result.count = 0;
            result.status = Pr76RootStatus::near_multiple;
            return result;
        }
        result.roots[result.count++] = {z, x, residual, pressure_residual, slope_ratio,
            slope_ratio > T{8} * std::sqrt(eps), dh > T{0} ? 1 : (dh < T{0} ? -1 : 0),
            iterations + 1};
    }
    result.status = result.count > 0 ? Pr76RootStatus::success : Pr76RootStatus::unrepresentable;
    return result;
}

} // namespace mpmc::thermodynamics
#endif // MPMC_THERMODYNAMICS_PR76_ROOTS_HPP
