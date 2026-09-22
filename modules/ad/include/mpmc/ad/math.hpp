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
#include <utility>

namespace mpmc::ad {

namespace detail {

inline void require_math_domain(
    bool condition,
    const char* message) {
    if (!condition) {
        throw std::domain_error(message);
    }
}

template <AdScalar Scalar>
void require_finite_math_input(
    const Scalar& value) {
    require_math_domain(
        std::isfinite(primal_value(value)),
        "mpmc::ad: non-finite elementary-function input");
}

template <
    AdScalar Scalar,
    std::size_t N,
    typename Transform>
[[nodiscard]] Dual<Scalar, N>
map_math_derivatives(
    const Dual<Scalar, N>& input,
    Scalar value,
    Transform transform) {
    typename Dual<Scalar, N>::Gradient
        gradient{};
    for (std::size_t i = 0;
         i < N;
         ++i) {
        gradient[i] =
            transform(
                input.derivatives()[i]);
    }
    return Dual<Scalar, N>{
        std::move(value),
        gradient};
}

template <AdScalar Scalar>
[[nodiscard]] Scalar
nonnegative_integer_power(
    Scalar base,
    unsigned int exponent) {
    Scalar result{
        base_scalar_t<Scalar>{1}};
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

template <AdScalar Scalar>
[[nodiscard]] Scalar hypot_one(
    const Scalar& value) {
    if constexpr (
        std::floating_point<Scalar>) {
        return std::hypot(
            Scalar{1},
            value);
    } else {
        using std::sqrt;
        return sqrt(
            Scalar{
                base_scalar_t<Scalar>{1}} +
            value * value);
    }
}

#define MPMC_AD_RECURSIVE_UNARY(name) template <AdScalar Scalar> [[nodiscard]] Scalar recursive_##name(const Scalar& x) {     using std::name;     return name(x); }

MPMC_AD_RECURSIVE_UNARY(exp)
MPMC_AD_RECURSIVE_UNARY(exp2)
MPMC_AD_RECURSIVE_UNARY(expm1)
MPMC_AD_RECURSIVE_UNARY(log)
MPMC_AD_RECURSIVE_UNARY(log2)
MPMC_AD_RECURSIVE_UNARY(log10)
MPMC_AD_RECURSIVE_UNARY(log1p)
MPMC_AD_RECURSIVE_UNARY(sqrt)
MPMC_AD_RECURSIVE_UNARY(cbrt)
MPMC_AD_RECURSIVE_UNARY(sin)
MPMC_AD_RECURSIVE_UNARY(cos)
MPMC_AD_RECURSIVE_UNARY(tan)
MPMC_AD_RECURSIVE_UNARY(asin)
MPMC_AD_RECURSIVE_UNARY(acos)
MPMC_AD_RECURSIVE_UNARY(atan)
MPMC_AD_RECURSIVE_UNARY(sinh)
MPMC_AD_RECURSIVE_UNARY(cosh)
MPMC_AD_RECURSIVE_UNARY(tanh)
MPMC_AD_RECURSIVE_UNARY(asinh)
MPMC_AD_RECURSIVE_UNARY(acosh)
MPMC_AD_RECURSIVE_UNARY(atanh)
MPMC_AD_RECURSIVE_UNARY(abs)

#undef MPMC_AD_RECURSIVE_UNARY

template <
    AdScalar Scalar,
    std::floating_point Exponent>
[[nodiscard]] Scalar
recursive_pow_real(
    const Scalar& x,
    Exponent exponent) {
    using std::pow;
    return pow(x, exponent);
}

template <AdScalar Scalar>
[[nodiscard]] Scalar
recursive_pow_ad(
    const Scalar& x,
    const Scalar& y) {
    using std::pow;
    return pow(x, y);
}

} // namespace detail

template <detail::AdScalar Scalar, std::size_t N>
[[nodiscard]] Dual<Scalar, N>
exp(const Dual<Scalar, N>& x) {
    detail::require_finite_math_input(
        x.value());
    const Scalar value =
        detail::recursive_exp(
            x.value());
    return detail::map_math_derivatives(
        x,
        value,
        [value](Scalar seed) {
            return seed * value;
        });
}

template <detail::AdScalar Scalar, std::size_t N>
[[nodiscard]] Dual<Scalar, N>
exp2(const Dual<Scalar, N>& x) {
    using Base =
        typename Dual<Scalar, N>::BaseScalar;
    detail::require_finite_math_input(
        x.value());
    const Scalar value =
        detail::recursive_exp2(
            x.value());
    const Scalar ln2{
        std::numbers::ln2_v<Base>};
    return detail::map_math_derivatives(
        x,
        value,
        [value, ln2](Scalar seed) {
            return (seed * ln2) * value;
        });
}

template <detail::AdScalar Scalar, std::size_t N>
[[nodiscard]] Dual<Scalar, N>
expm1(const Dual<Scalar, N>& x) {
    detail::require_finite_math_input(
        x.value());
    const Scalar value =
        detail::recursive_expm1(
            x.value());
    const Scalar slope =
        detail::recursive_exp(
            x.value());
    return detail::map_math_derivatives(
        x,
        value,
        [slope](Scalar seed) {
            return seed * slope;
        });
}

template <detail::AdScalar Scalar, std::size_t N>
[[nodiscard]] Dual<Scalar, N>
log(const Dual<Scalar, N>& x) {
    using Base =
        typename Dual<Scalar, N>::BaseScalar;
    const Scalar a = x.value();
    detail::require_finite_math_input(a);
    detail::require_math_domain(
        detail::primal_value(a) > Base{0},
        "mpmc::ad::log: requires x > 0");
    return detail::map_math_derivatives(
        x,
        detail::recursive_log(a),
        [a](Scalar seed) {
            return seed / a;
        });
}

template <detail::AdScalar Scalar, std::size_t N>
[[nodiscard]] Dual<Scalar, N>
log2(const Dual<Scalar, N>& x) {
    using Base =
        typename Dual<Scalar, N>::BaseScalar;
    const Scalar a = x.value();
    detail::require_finite_math_input(a);
    detail::require_math_domain(
        detail::primal_value(a) > Base{0},
        "mpmc::ad::log2: requires x > 0");
    const Scalar ln2{
        std::numbers::ln2_v<Base>};
    return detail::map_math_derivatives(
        x,
        detail::recursive_log2(a),
        [a, ln2](Scalar seed) {
            return (seed / a) / ln2;
        });
}

template <detail::AdScalar Scalar, std::size_t N>
[[nodiscard]] Dual<Scalar, N>
log10(const Dual<Scalar, N>& x) {
    using Base =
        typename Dual<Scalar, N>::BaseScalar;
    const Scalar a = x.value();
    detail::require_finite_math_input(a);
    detail::require_math_domain(
        detail::primal_value(a) > Base{0},
        "mpmc::ad::log10: requires x > 0");
    const Scalar ln10{
        std::numbers::ln10_v<Base>};
    return detail::map_math_derivatives(
        x,
        detail::recursive_log10(a),
        [a, ln10](Scalar seed) {
            return (seed / a) / ln10;
        });
}

template <detail::AdScalar Scalar, std::size_t N>
[[nodiscard]] Dual<Scalar, N>
log1p(const Dual<Scalar, N>& x) {
    using Base =
        typename Dual<Scalar, N>::BaseScalar;
    const Scalar a = x.value();
    detail::require_finite_math_input(a);
    detail::require_math_domain(
        detail::primal_value(a) > Base{-1},
        "mpmc::ad::log1p: requires x > -1");
    const Scalar one{Base{1}};
    return detail::map_math_derivatives(
        x,
        detail::recursive_log1p(a),
        [a, one](Scalar seed) {
            return seed / (one + a);
        });
}

template <detail::AdScalar Scalar, std::size_t N>
[[nodiscard]] Dual<Scalar, N>
sqrt(const Dual<Scalar, N>& x) {
    using Base =
        typename Dual<Scalar, N>::BaseScalar;
    detail::require_finite_math_input(
        x.value());
    detail::require_math_domain(
        detail::primal_value(x.value()) >
            Base{0},
        "mpmc::ad::sqrt: requires x > 0");
    const Scalar value =
        detail::recursive_sqrt(
            x.value());
    const Scalar two{Base{2}};
    return detail::map_math_derivatives(
        x,
        value,
        [value, two](Scalar seed) {
            return (seed / value) / two;
        });
}

template <detail::AdScalar Scalar, std::size_t N>
[[nodiscard]] Dual<Scalar, N>
cbrt(const Dual<Scalar, N>& x) {
    using Base =
        typename Dual<Scalar, N>::BaseScalar;
    detail::require_finite_math_input(
        x.value());
    detail::require_math_domain(
        detail::primal_value(x.value()) !=
            Base{0},
        "mpmc::ad::cbrt: singular derivative at zero");
    const Scalar value =
        detail::recursive_cbrt(
            x.value());
    const Scalar three{Base{3}};
    return detail::map_math_derivatives(
        x,
        value,
        [value, three](Scalar seed) {
            return ((seed / value) / value) /
                three;
        });
}

template <detail::AdScalar Scalar, std::size_t N>
[[nodiscard]] Dual<Scalar, N>
sin(const Dual<Scalar, N>& x) {
    detail::require_finite_math_input(
        x.value());
    const Scalar slope =
        detail::recursive_cos(
            x.value());
    return detail::map_math_derivatives(
        x,
        detail::recursive_sin(
            x.value()),
        [slope](Scalar seed) {
            return seed * slope;
        });
}

template <detail::AdScalar Scalar, std::size_t N>
[[nodiscard]] Dual<Scalar, N>
cos(const Dual<Scalar, N>& x) {
    detail::require_finite_math_input(
        x.value());
    const Scalar slope =
        -detail::recursive_sin(
            x.value());
    return detail::map_math_derivatives(
        x,
        detail::recursive_cos(
            x.value()),
        [slope](Scalar seed) {
            return seed * slope;
        });
}

template <detail::AdScalar Scalar, std::size_t N>
[[nodiscard]] Dual<Scalar, N>
tan(const Dual<Scalar, N>& x) {
    using Base =
        typename Dual<Scalar, N>::BaseScalar;
    detail::require_finite_math_input(
        x.value());
    const Scalar cosine =
        detail::recursive_cos(
            x.value());
    detail::require_math_domain(
        detail::primal_value(cosine) !=
            Base{0},
        "mpmc::ad::tan: zero cosine at a pole");
    return detail::map_math_derivatives(
        x,
        detail::recursive_tan(
            x.value()),
        [cosine](Scalar seed) {
            return (seed / cosine) /
                cosine;
        });
}

template <detail::AdScalar Scalar, std::size_t N>
[[nodiscard]] Dual<Scalar, N>
asin(const Dual<Scalar, N>& x) {
    using Base =
        typename Dual<Scalar, N>::BaseScalar;
    const Scalar a = x.value();
    detail::require_finite_math_input(a);
    const Base primal =
        detail::primal_value(a);
    detail::require_math_domain(
        primal > Base{-1} &&
            primal < Base{1},
        "mpmc::ad::asin: requires -1 < x < 1");
    const Scalar one{Base{1}};
    const Scalar denominator =
        detail::recursive_sqrt(
            (one - a) * (one + a));
    return detail::map_math_derivatives(
        x,
        detail::recursive_asin(a),
        [denominator](Scalar seed) {
            return seed / denominator;
        });
}

template <detail::AdScalar Scalar, std::size_t N>
[[nodiscard]] Dual<Scalar, N>
acos(const Dual<Scalar, N>& x) {
    using Base =
        typename Dual<Scalar, N>::BaseScalar;
    const Scalar a = x.value();
    detail::require_finite_math_input(a);
    const Base primal =
        detail::primal_value(a);
    detail::require_math_domain(
        primal > Base{-1} &&
            primal < Base{1},
        "mpmc::ad::acos: requires -1 < x < 1");
    const Scalar one{Base{1}};
    const Scalar denominator =
        detail::recursive_sqrt(
            (one - a) * (one + a));
    return detail::map_math_derivatives(
        x,
        detail::recursive_acos(a),
        [denominator](Scalar seed) {
            return -seed / denominator;
        });
}

template <detail::AdScalar Scalar, std::size_t N>
[[nodiscard]] Dual<Scalar, N>
atan(const Dual<Scalar, N>& x) {
    detail::require_finite_math_input(
        x.value());
    const Scalar scale =
        detail::hypot_one(
            x.value());
    return detail::map_math_derivatives(
        x,
        detail::recursive_atan(
            x.value()),
        [scale](Scalar seed) {
            return (seed / scale) / scale;
        });
}

template <detail::AdScalar Scalar, std::size_t N>
[[nodiscard]] Dual<Scalar, N>
sinh(const Dual<Scalar, N>& x) {
    detail::require_finite_math_input(
        x.value());
    const Scalar slope =
        detail::recursive_cosh(
            x.value());
    return detail::map_math_derivatives(
        x,
        detail::recursive_sinh(
            x.value()),
        [slope](Scalar seed) {
            return seed * slope;
        });
}

template <detail::AdScalar Scalar, std::size_t N>
[[nodiscard]] Dual<Scalar, N>
cosh(const Dual<Scalar, N>& x) {
    detail::require_finite_math_input(
        x.value());
    const Scalar slope =
        detail::recursive_sinh(
            x.value());
    return detail::map_math_derivatives(
        x,
        detail::recursive_cosh(
            x.value()),
        [slope](Scalar seed) {
            return seed * slope;
        });
}

template <detail::AdScalar Scalar, std::size_t N>
[[nodiscard]] Dual<Scalar, N>
tanh(const Dual<Scalar, N>& x) {
    using Base =
        typename Dual<Scalar, N>::BaseScalar;
    detail::require_finite_math_input(
        x.value());
    const Scalar a = x.value();
    const Scalar exponential =
        detail::primal_value(a) >= Base{0}
            ? detail::recursive_exp(-a)
            : detail::recursive_exp(a);
    const Scalar one{Base{1}};
    const Scalar factor =
        exponential /
        (one + exponential * exponential);
    const Scalar four{Base{4}};
    return detail::map_math_derivatives(
        x,
        detail::recursive_tanh(a),
        [factor, four](Scalar seed) {
            return (seed * factor) *
                (four * factor);
        });
}

template <detail::AdScalar Scalar, std::size_t N>
[[nodiscard]] Dual<Scalar, N>
asinh(const Dual<Scalar, N>& x) {
    detail::require_finite_math_input(
        x.value());
    const Scalar denominator =
        detail::hypot_one(
            x.value());
    return detail::map_math_derivatives(
        x,
        detail::recursive_asinh(
            x.value()),
        [denominator](Scalar seed) {
            return seed / denominator;
        });
}

template <detail::AdScalar Scalar, std::size_t N>
[[nodiscard]] Dual<Scalar, N>
acosh(const Dual<Scalar, N>& x) {
    using Base =
        typename Dual<Scalar, N>::BaseScalar;
    const Scalar a = x.value();
    detail::require_finite_math_input(a);
    detail::require_math_domain(
        detail::primal_value(a) > Base{1},
        "mpmc::ad::acosh: requires x > 1");
    const Scalar one{Base{1}};
    const Scalar left =
        detail::recursive_sqrt(a - one);
    const Scalar right =
        detail::recursive_sqrt(a + one);
    return detail::map_math_derivatives(
        x,
        detail::recursive_acosh(a),
        [left, right](Scalar seed) {
            return (seed / left) / right;
        });
}

template <detail::AdScalar Scalar, std::size_t N>
[[nodiscard]] Dual<Scalar, N>
atanh(const Dual<Scalar, N>& x) {
    using Base =
        typename Dual<Scalar, N>::BaseScalar;
    const Scalar a = x.value();
    detail::require_finite_math_input(a);
    const Base primal =
        detail::primal_value(a);
    detail::require_math_domain(
        primal > Base{-1} &&
            primal < Base{1},
        "mpmc::ad::atanh: requires -1 < x < 1");
    const Scalar one{Base{1}};
    const Scalar denominator =
        (one - a) * (one + a);
    return detail::map_math_derivatives(
        x,
        detail::recursive_atanh(a),
        [denominator](Scalar seed) {
            return seed / denominator;
        });
}

template <detail::AdScalar Scalar, std::size_t N>
[[nodiscard]] Dual<Scalar, N>
abs(const Dual<Scalar, N>& x) {
    using Base =
        typename Dual<Scalar, N>::BaseScalar;
    detail::require_finite_math_input(
        x.value());
    const Base primal =
        detail::primal_value(
            x.value());
    detail::require_math_domain(
        primal != Base{0},
        "mpmc::ad::abs: not differentiable at zero");
    const Scalar sign{
        primal < Base{0}
            ? Base{-1}
            : Base{1}};
    return detail::map_math_derivatives(
        x,
        detail::recursive_abs(
            x.value()),
        [sign](Scalar seed) {
            return sign * seed;
        });
}

template <
    detail::AdScalar Scalar,
    std::size_t N,
    typename Integer>
    requires std::same_as<Integer, int>
[[nodiscard]] Dual<Scalar, N>
pow(
    const Dual<Scalar, N>& x,
    Integer exponent) {
    using Base =
        typename Dual<Scalar, N>::BaseScalar;
    const Scalar a = x.value();
    detail::require_finite_math_input(a);
    if (exponent == 0) {
        return Dual<Scalar, N>{
            Scalar{Base{1}}};
    }

    const bool negative =
        exponent < 0;
    detail::require_math_domain(
        !negative ||
            detail::primal_value(a) !=
                Base{0},
        "mpmc::ad::pow: zero to negative power");

    const unsigned int magnitude =
        negative
            ? 0U -
                static_cast<unsigned int>(
                    exponent)
            : static_cast<unsigned int>(
                exponent);
    const Scalar one{Base{1}};
    const Scalar base =
        negative
            ? one / a
            : a;
    const Scalar value =
        detail::nonnegative_integer_power(
            base,
            magnitude);
    const Scalar adjacent =
        detail::nonnegative_integer_power(
            base,
            magnitude - 1U);
    const Scalar order{
        static_cast<Base>(
            exponent)};

    return detail::map_math_derivatives(
        x,
        value,
        [=](Scalar seed) {
            return negative
                ? (((seed / a) * adjacent) /
                   a) *
                      order
                : (seed * adjacent) *
                      order;
        });
}

template <
    detail::AdScalar Scalar,
    std::size_t N,
    std::floating_point Exponent>
[[nodiscard]] Dual<Scalar, N>
pow(
    const Dual<Scalar, N>& x,
    Exponent real_exponent) {
    using Base =
        typename Dual<Scalar, N>::BaseScalar;
    detail::require_finite_math_input(
        real_exponent);

    if constexpr (
        std::numeric_limits<Exponent>::
                max_exponent >
            std::numeric_limits<Base>::
                max_exponent ||
        (std::numeric_limits<Exponent>::
                 max_exponent ==
             std::numeric_limits<Base>::
                 max_exponent &&
         std::numeric_limits<Exponent>::
                 digits >
             std::numeric_limits<Base>::
                 digits)) {
        const Exponent bound =
            static_cast<Exponent>(
                std::numeric_limits<Base>::
                    max());
        detail::require_math_domain(
            real_exponent >= -bound &&
                real_exponent <= bound,
            "mpmc::ad::pow: exponent is outside the scalar range");
    }

    const Base exponent =
        static_cast<Base>(
            real_exponent);
    const Scalar a = x.value();
    detail::require_finite_math_input(a);
    detail::require_math_domain(
        detail::primal_value(a) > Base{0},
        "mpmc::ad::pow: real exponent requires x > 0");

    if (exponent == Base{0}) {
        return Dual<Scalar, N>{
            Scalar{Base{1}}};
    }
    if (exponent == Base{1}) {
        return x;
    }

    const Scalar adjacent =
        detail::recursive_pow_real(
            a,
            exponent - Base{1});
    const Scalar value =
        detail::recursive_pow_real(
            a,
            exponent);
    const Scalar exponent_scalar{
        exponent};

    return detail::map_math_derivatives(
        x,
        value,
        [=](Scalar seed) {
            return (seed * adjacent) *
                exponent_scalar;
        });
}

template <detail::AdScalar Scalar, std::size_t N>
[[nodiscard]] Dual<Scalar, N>
pow(
    const Dual<Scalar, N>& x,
    const Dual<Scalar, N>& y) {
    using Base =
        typename Dual<Scalar, N>::BaseScalar;
    const Scalar a = x.value();
    const Scalar b = y.value();
    detail::require_finite_math_input(a);
    detail::require_finite_math_input(b);
    detail::require_math_domain(
        detail::primal_value(a) > Base{0},
        "mpmc::ad::pow: AD exponent requires x > 0");

    const Scalar value =
        detail::recursive_pow_ad(
            a,
            b);

    Scalar adjacent{Base{1}};
    if constexpr (
        std::floating_point<Scalar>) {
        adjacent =
            (b == Scalar{0} ||
             b == Scalar{1})
                ? Scalar{1}
                : std::pow(
                      a,
                      b - Scalar{1});
    } else {
        adjacent =
            detail::recursive_pow_ad(
                a,
                b - Scalar{Base{1}});
    }

    const Scalar logarithm =
        detail::recursive_log(a);
    typename Dual<Scalar, N>::Gradient
        gradient{};
    for (std::size_t i = 0;
         i < N;
         ++i) {
        gradient[i] =
            (x.derivatives()[i] *
                 adjacent) *
                b +
            (y.derivatives()[i] *
                 logarithm) *
                value;
    }
    return Dual<Scalar, N>{
        value,
        gradient};
}

template <
    std::floating_point BaseInput,
    detail::AdScalar Scalar,
    std::size_t N>
[[nodiscard]] Dual<Scalar, N>
pow(
    BaseInput base_input,
    const Dual<Scalar, N>& y) {
    using Base =
        typename Dual<Scalar, N>::BaseScalar;

    detail::require_finite_math_input(
        base_input);
    const long double raw =
        static_cast<long double>(
            base_input);
    const long double limit =
        static_cast<long double>(
            std::numeric_limits<Base>::
                max());
    detail::require_math_domain(
        raw >= -limit &&
            raw <= limit,
        "mpmc::ad::pow: base is outside the scalar range");

    const Base base_value =
        static_cast<Base>(
            base_input);
    detail::require_math_domain(
        base_value > Base{0},
        "mpmc::ad::pow: constant base requires a > 0");
    detail::require_finite_math_input(
        y.value());

    const Scalar base{
        base_value};
    const Scalar value =
        detail::recursive_pow_ad(
            base,
            y.value());
    const Scalar logarithm{
        std::log(base_value)};

    return detail::map_math_derivatives(
        y,
        value,
        [=](Scalar seed) {
            return (seed * logarithm) *
                value;
        });
}

} // namespace mpmc::ad

#endif // MPMC_AD_MATH_HPP
