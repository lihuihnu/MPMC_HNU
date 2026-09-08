#ifndef MPMC_AD_MATH_HPP
#define MPMC_AD_MATH_HPP

#include <mpmc/ad/dual.hpp>

#include <cmath>
#include <concepts>
#include <cstddef>
#include <limits>
#include <numbers>
#include <stdexcept>
#include <type_traits>

namespace mpmc::ad {

// Checked, real-valued elementary functions. Every primal argument must be
// finite and inside the differentiable domain, even for a zero derivative seed.
// Domain violations throw std::domain_error without changing the arguments.
// Range errors follow the underlying floating-point arithmetic (Inf/NaN/zero);
// no clipping, finite-difference derivatives, or errno/fenv management is added.
// See math.md for each domain, formula, floating-point limit, and test contract.
namespace detail {

inline void require_math_domain(bool condition, const char* message) {
    if (!condition) {
        throw std::domain_error(message);
    }
}

template <std::floating_point T>
void require_finite_math_input(T value) {
    require_math_domain(std::isfinite(value), "mpmc::ad: non-finite elementary-function input");
}

// Transform each seed directly rather than always forming f'(x) first. For
// example, seed/x can remain finite when (1/x)*seed would overflow needlessly.
// The closure owns only scalar temporaries; no tape or heap allocation is used.
template <std::floating_point T, std::size_t N, typename Transform>
[[nodiscard]] Dual<T, N> map_math_derivatives(const Dual<T, N>& input, T value,
                                             Transform transform) {
    typename Dual<T, N>::Gradient gradient{};
    for (std::size_t i = 0; i < N; ++i) {
        gradient[i] = transform(input.derivatives()[i]);
    }
    return Dual<T, N>{value, gradient};
}

// Integer exponents stay integers, preserving parity even when T cannot exactly
// represent the exponent. Only unsigned shifts are used, including for INT_MIN.
template <std::floating_point T>
[[nodiscard]] T nonnegative_integer_power(T base, unsigned int exponent) {
    T result{1};
    while (exponent != 0U) {
        if ((exponent & 1U) != 0U) {
            result *= base;
        }
        exponent >>= 1U;
        if (exponent != 0U) {
            base *= base;
        }
    }
    return result;
}

} // namespace detail

/// exp(x): x in R; df = exp(x) dx.
template <std::floating_point T, std::size_t N>
[[nodiscard]] Dual<T, N> exp(const Dual<T, N>& x) {
    detail::require_finite_math_input(x.value());
    const T value = std::exp(x.value());
    return detail::map_math_derivatives(x, value, [value](T seed) { return seed * value; });
}

/// exp2(x): x in R; df = log(2) * 2^x dx.
template <std::floating_point T, std::size_t N>
[[nodiscard]] Dual<T, N> exp2(const Dual<T, N>& x) {
    detail::require_finite_math_input(x.value());
    const T value = std::exp2(x.value());
    return detail::map_math_derivatives(x, value, [value](T seed) {
        return (seed * std::numbers::ln2_v<T>) * value;
    });
}

/// expm1(x): x in R; use std::expm1 for the value and exp(x) for the slope.
template <std::floating_point T, std::size_t N>
[[nodiscard]] Dual<T, N> expm1(const Dual<T, N>& x) {
    detail::require_finite_math_input(x.value());
    const T value = std::expm1(x.value());
    // value + 1 loses a small but nonzero derivative when value rounds to -1.
    const T slope = std::exp(x.value());
    return detail::map_math_derivatives(x, value, [slope](T seed) { return seed * slope; });
}

/// log(x): x > 0; zero and negative inputs throw, including constant Duals.
template <std::floating_point T, std::size_t N>
[[nodiscard]] Dual<T, N> log(const Dual<T, N>& x) {
    const T a = x.value();
    detail::require_finite_math_input(a);
    detail::require_math_domain(a > T{0}, "mpmc::ad::log: requires x > 0");
    return detail::map_math_derivatives(x, std::log(a), [a](T seed) { return seed / a; });
}

/// log2(x): x > 0; df = dx / (x log(2)).
template <std::floating_point T, std::size_t N>
[[nodiscard]] Dual<T, N> log2(const Dual<T, N>& x) {
    const T a = x.value();
    detail::require_finite_math_input(a);
    detail::require_math_domain(a > T{0}, "mpmc::ad::log2: requires x > 0");
    return detail::map_math_derivatives(x, std::log2(a), [a](T seed) {
        return (seed / a) / std::numbers::ln2_v<T>;
    });
}

/// log10(x): x > 0; split the divisions to avoid overflowing x*log(10).
template <std::floating_point T, std::size_t N>
[[nodiscard]] Dual<T, N> log10(const Dual<T, N>& x) {
    const T a = x.value();
    detail::require_finite_math_input(a);
    detail::require_math_domain(a > T{0}, "mpmc::ad::log10: requires x > 0");
    return detail::map_math_derivatives(x, std::log10(a), [a](T seed) {
        return (seed / a) / std::numbers::ln10_v<T>;
    });
}

/// log1p(x): x > -1; df = dx/(1+x), with an accurate value near x = 0.
template <std::floating_point T, std::size_t N>
[[nodiscard]] Dual<T, N> log1p(const Dual<T, N>& x) {
    const T a = x.value();
    detail::require_finite_math_input(a);
    detail::require_math_domain(a > T{-1}, "mpmc::ad::log1p: requires x > -1");
    return detail::map_math_derivatives(x, std::log1p(a), [a](T seed) {
        return seed / (T{1} + a);
    });
}

/// sqrt(x): x > 0. The value exists at zero but has no finite derivative there.
template <std::floating_point T, std::size_t N>
[[nodiscard]] Dual<T, N> sqrt(const Dual<T, N>& x) {
    detail::require_finite_math_input(x.value());
    detail::require_math_domain(x.value() > T{0}, "mpmc::ad::sqrt: requires x > 0");
    const T value = std::sqrt(x.value());
    return detail::map_math_derivatives(x, value, [value](T seed) {
        return (seed / value) / T{2};
    });
}

/// cbrt(x): x != 0; the real cube root also accepts negative x.
template <std::floating_point T, std::size_t N>
[[nodiscard]] Dual<T, N> cbrt(const Dual<T, N>& x) {
    detail::require_finite_math_input(x.value());
    detail::require_math_domain(x.value() != T{0}, "mpmc::ad::cbrt: singular derivative at zero");
    const T value = std::cbrt(x.value());
    return detail::map_math_derivatives(x, value, [value](T seed) {
        return ((seed / value) / value) / T{3};
    });
}

/// sin(x): x in R, radians; df = cos(x) dx.
template <std::floating_point T, std::size_t N>
[[nodiscard]] Dual<T, N> sin(const Dual<T, N>& x) {
    detail::require_finite_math_input(x.value());
    const T slope = std::cos(x.value());
    return detail::map_math_derivatives(x, std::sin(x.value()), [slope](T seed) {
        return seed * slope;
    });
}

/// cos(x): x in R, radians; df = -sin(x) dx.
template <std::floating_point T, std::size_t N>
[[nodiscard]] Dual<T, N> cos(const Dual<T, N>& x) {
    detail::require_finite_math_input(x.value());
    const T slope = -std::sin(x.value());
    return detail::map_math_derivatives(x, std::cos(x.value()), [slope](T seed) {
        return seed * slope;
    });
}

/// tan(x): excludes pi/2 + k*pi. Floating approximations to pi/2 are not exact
/// poles; no epsilon exclusion band is imposed. A computed zero cosine throws.
template <std::floating_point T, std::size_t N>
[[nodiscard]] Dual<T, N> tan(const Dual<T, N>& x) {
    detail::require_finite_math_input(x.value());
    const T cosine = std::cos(x.value());
    detail::require_math_domain(cosine != T{0}, "mpmc::ad::tan: zero cosine at a pole");
    return detail::map_math_derivatives(x, std::tan(x.value()), [cosine](T seed) {
        return (seed / cosine) / cosine;
    });
}

/// asin(x): -1 < x < 1; endpoints have unbounded derivatives and throw.
template <std::floating_point T, std::size_t N>
[[nodiscard]] Dual<T, N> asin(const Dual<T, N>& x) {
    const T a = x.value();
    detail::require_finite_math_input(a);
    detail::require_math_domain(a > T{-1} && a < T{1}, "mpmc::ad::asin: requires -1 < x < 1");
    const T denominator = std::sqrt((T{1} - a) * (T{1} + a));
    return detail::map_math_derivatives(x, std::asin(a), [denominator](T seed) {
        return seed / denominator;
    });
}

/// acos(x): -1 < x < 1; df = -dx/sqrt(1-x^2).
template <std::floating_point T, std::size_t N>
[[nodiscard]] Dual<T, N> acos(const Dual<T, N>& x) {
    const T a = x.value();
    detail::require_finite_math_input(a);
    detail::require_math_domain(a > T{-1} && a < T{1}, "mpmc::ad::acos: requires -1 < x < 1");
    const T denominator = std::sqrt((T{1} - a) * (T{1} + a));
    return detail::map_math_derivatives(x, std::acos(a), [denominator](T seed) {
        return -seed / denominator;
    });
}

/// atan(x): x in R; hypot and sequential divisions avoid forming x*x.
template <std::floating_point T, std::size_t N>
[[nodiscard]] Dual<T, N> atan(const Dual<T, N>& x) {
    detail::require_finite_math_input(x.value());
    const T scale = std::hypot(T{1}, x.value());
    return detail::map_math_derivatives(x, std::atan(x.value()), [scale](T seed) {
        return (seed / scale) / scale;
    });
}

/// sinh(x): x in R; df = cosh(x) dx.
template <std::floating_point T, std::size_t N>
[[nodiscard]] Dual<T, N> sinh(const Dual<T, N>& x) {
    detail::require_finite_math_input(x.value());
    const T slope = std::cosh(x.value());
    return detail::map_math_derivatives(x, std::sinh(x.value()), [slope](T seed) {
        return seed * slope;
    });
}

/// cosh(x): x in R; df = sinh(x) dx.
template <std::floating_point T, std::size_t N>
[[nodiscard]] Dual<T, N> cosh(const Dual<T, N>& x) {
    detail::require_finite_math_input(x.value());
    const T slope = std::sinh(x.value());
    return detail::map_math_derivatives(x, std::cosh(x.value()), [slope](T seed) {
        return seed * slope;
    });
}

/// tanh(x): x in R; compute sech(x)^2 without subtracting two nearly equal
/// values after tanh(x) rounds to +/-1. Split factors preserve scaled seeds.
template <std::floating_point T, std::size_t N>
[[nodiscard]] Dual<T, N> tanh(const Dual<T, N>& x) {
    detail::require_finite_math_input(x.value());
    const T exponential = std::exp(-std::abs(x.value()));
    const T factor = exponential / (T{1} + exponential * exponential);
    return detail::map_math_derivatives(x, std::tanh(x.value()), [factor](T seed) {
        return (seed * factor) * (T{4} * factor);
    });
}

/// asinh(x): x in R; df = dx/sqrt(1+x^2).
template <std::floating_point T, std::size_t N>
[[nodiscard]] Dual<T, N> asinh(const Dual<T, N>& x) {
    detail::require_finite_math_input(x.value());
    const T denominator = std::hypot(T{1}, x.value());
    return detail::map_math_derivatives(x, std::asinh(x.value()), [denominator](T seed) {
        return seed / denominator;
    });
}

/// acosh(x): x > 1; factor x^2-1 and divide separately to avoid overflow.
template <std::floating_point T, std::size_t N>
[[nodiscard]] Dual<T, N> acosh(const Dual<T, N>& x) {
    const T a = x.value();
    detail::require_finite_math_input(a);
    detail::require_math_domain(a > T{1}, "mpmc::ad::acosh: requires x > 1");
    const T left = std::sqrt(a - T{1});
    const T right = std::sqrt(a + T{1});
    return detail::map_math_derivatives(x, std::acosh(a), [left, right](T seed) {
        return (seed / left) / right;
    });
}

/// atanh(x): -1 < x < 1; df = dx/(1-x^2).
template <std::floating_point T, std::size_t N>
[[nodiscard]] Dual<T, N> atanh(const Dual<T, N>& x) {
    const T a = x.value();
    detail::require_finite_math_input(a);
    detail::require_math_domain(a > T{-1} && a < T{1}, "mpmc::ad::atanh: requires -1 < x < 1");
    const T denominator = (T{1} - a) * (T{1} + a);
    return detail::map_math_derivatives(x, std::atanh(a), [denominator](T seed) {
        return seed / denominator;
    });
}

/// abs(x): x != 0; do not silently choose a subgradient at the kink.
template <std::floating_point T, std::size_t N>
[[nodiscard]] Dual<T, N> abs(const Dual<T, N>& x) {
    detail::require_finite_math_input(x.value());
    detail::require_math_domain(x.value() != T{0}, "mpmc::ad::abs: not differentiable at zero");
    const T sign = x.value() < T{0} ? T{-1} : T{1};
    return detail::map_math_derivatives(x, std::abs(x.value()), [sign](T seed) {
        return sign * seed;
    });
}

/// Integer power: all finite x for n >= 0, nonzero x for n < 0.
/// The polynomial convention pow(x, 0) == 1 includes x == 0, with zero gradient.
/// In contrast, floating/AD exponents below always require a positive base.
template <std::floating_point T, std::size_t N, typename Integer>
    requires std::same_as<Integer, int>
[[nodiscard]] Dual<T, N> pow(const Dual<T, N>& x, Integer exponent) {
    const T a = x.value();
    detail::require_finite_math_input(a);
    if (exponent == 0) {
        return Dual<T, N>{T{1}};
    }
    const bool negative = exponent < 0;
    detail::require_math_domain(!negative || a != T{0}, "mpmc::ad::pow: zero to negative power");
    // Unsigned subtraction handles INT_MIN without a signed negation overflow.
    const unsigned int magnitude = negative ? 0U - static_cast<unsigned int>(exponent)
                                            : static_cast<unsigned int>(exponent);
    const T base = negative ? T{1} / a : a;
    const T value = detail::nonnegative_integer_power(base, magnitude);
    const T adjacent = detail::nonnegative_integer_power(base, magnitude - 1U);
    const T n = static_cast<T>(exponent);
    return detail::map_math_derivatives(x, value, [=](T seed) {
        // A negative exponent uses two divisions, not an eagerly formed 1/x^2.
        // The positive path also keeps a derivative when x^n itself underflows.
        return negative ? (((seed / a) * adjacent) / a) * n : (seed * adjacent) * n;
    });
}

/// Constant real power: x > 0 and finite exponent; df = p*x^(p-1) dx.
/// The floating exponent is explicitly converted to T; no precision promotion.
/// Constrained overloads prevent floating exponents from converting to int.
template <std::floating_point T, std::size_t N, std::floating_point Exponent>
[[nodiscard]] Dual<T, N> pow(const Dual<T, N>& x, Exponent real_exponent) {
    detail::require_finite_math_input(real_exponent);
    if constexpr (std::numeric_limits<Exponent>::max_exponent >
                      std::numeric_limits<T>::max_exponent ||
                  (std::numeric_limits<Exponent>::max_exponent ==
                       std::numeric_limits<T>::max_exponent &&
                   std::numeric_limits<Exponent>::digits > std::numeric_limits<T>::digits)) {
        // Reject an out-of-range narrowing conversion before performing it.
        const Exponent bound = static_cast<Exponent>(std::numeric_limits<T>::max());
        detail::require_math_domain(real_exponent >= -bound && real_exponent <= bound,
                                    "mpmc::ad::pow: exponent is outside the scalar range");
    }
    const T exponent = static_cast<T>(real_exponent);
    const T a = x.value();
    detail::require_finite_math_input(a);
    detail::require_finite_math_input(exponent);
    detail::require_math_domain(a > T{0}, "mpmc::ad::pow: real exponent requires x > 0");
    if (exponent == T{0}) {
        return Dual<T, N>{T{1}};
    }
    if (exponent == T{1}) {
        return x;
    }
    const T adjacent = std::pow(a, exponent - T{1});
    return detail::map_math_derivatives(x, std::pow(a, exponent), [=](T seed) {
        return (seed * adjacent) * exponent;
    });
}

/// Variable exponent: x > 0, y in R; df = y*x^(y-1) dx + x^y*log(x) dy.
/// A zero seed does not enlarge the joint differentiable domain.
template <std::floating_point T, std::size_t N>
[[nodiscard]] Dual<T, N> pow(const Dual<T, N>& x, const Dual<T, N>& y) {
    const T a = x.value();
    const T b = y.value();
    detail::require_finite_math_input(a);
    detail::require_finite_math_input(b);
    detail::require_math_domain(a > T{0}, "mpmc::ad::pow: AD exponent requires x > 0");
    const T value = std::pow(a, b);
    // Avoid 0*pow(a,-1) when b==0, and use the exact identity when b==1.
    const T adjacent = (b == T{0} || b == T{1}) ? T{1} : std::pow(a, b - T{1});
    const T logarithm = std::log(a);
    typename Dual<T, N>::Gradient gradient{};
    for (std::size_t i = 0; i < N; ++i) {
        gradient[i] = (x.derivatives()[i] * adjacent) * b +
                      (y.derivatives()[i] * logarithm) * value;
    }
    return Dual<T, N>{value, gradient};
}

/// Constant positive base: a > 0; df = a^y*log(a) dy.
template <std::floating_point T, std::size_t N>
[[nodiscard]] Dual<T, N> pow(std::type_identity_t<T> base, const Dual<T, N>& y) {
    detail::require_finite_math_input(base);
    detail::require_finite_math_input(y.value());
    detail::require_math_domain(base > T{0}, "mpmc::ad::pow: constant base requires a > 0");
    const T value = std::pow(base, y.value());
    const T logarithm = std::log(base);
    return detail::map_math_derivatives(y, value, [=](T seed) {
        return (seed * logarithm) * value;
    });
}

} // namespace mpmc::ad

#endif // MPMC_AD_MATH_HPP
