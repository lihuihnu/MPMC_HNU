#ifndef MPMC_AD_DUAL_HPP
#define MPMC_AD_DUAL_HPP

#include <array>
#include <concepts>
#include <cstddef>
#include <stdexcept>
#include <type_traits>

namespace mpmc::ad {

/// A first-order forward-mode scalar with N independently seeded derivative directions.
/// Owns all storage; depends only on the C++ standard library. See the module
/// README for the numerical contract and design references.
/// T is a built-in floating-point type; nested/higher-order AD is not supported.
template <std::floating_point T, std::size_t N = 1>
    requires(N > 0 && std::same_as<T, std::remove_cv_t<T>>)
class Dual {
public:
    using Scalar = T;
    using Gradient = std::array<T, N>;
    static constexpr std::size_t derivative_count = N;

    /// Zero is a constant: both its value and all derivatives are zero.
    constexpr Dual() noexcept = default;

    /// Constants are explicit to avoid accidental scalar/AD type changes.
    explicit constexpr Dual(T value) noexcept : value_(value) {}

    /// Supply arbitrary directional seeds (identity seeds produce a gradient).
    constexpr Dual(T value, const Gradient& derivatives) noexcept
        : value_(value), derivatives_(derivatives) {}

    /// Seed one independent direction. Bounds remain checked in release builds.
    [[nodiscard]] static constexpr Dual variable(T value, std::size_t index) {
        Dual result{value};
        result.derivatives_.at(index) = T{1};
        return result;
    }

    [[nodiscard]] constexpr T value() const noexcept { return value_; }

    [[nodiscard]] constexpr T derivative(std::size_t index) const {
        return derivatives_.at(index);
    }

    [[nodiscard]] constexpr const Gradient& derivatives() const & noexcept {
        return derivatives_;
    }

    /// Returning an owned array prevents a dangling view into a temporary Dual.
    [[nodiscard]] constexpr Gradient derivatives() const && noexcept {
        return derivatives_;
    }

    constexpr Dual& operator+=(const Dual& rhs) noexcept {
        value_ += rhs.value_;
        for (std::size_t i = 0; i < N; ++i) {
            derivatives_[i] += rhs.derivatives_[i];
        }
        return *this;
    }

    constexpr Dual& operator-=(const Dual& rhs) noexcept {
        value_ -= rhs.value_;
        for (std::size_t i = 0; i < N; ++i) {
            derivatives_[i] -= rhs.derivatives_[i];
        }
        return *this;
    }

    constexpr Dual& operator*=(const Dual& rhs) noexcept {
        // Preserve both primal values until all derivative lanes are updated.
        // Lanes are independent, so this is also safe for x *= x.
        const T left_value = value_;
        const T right_value = rhs.value_;
        for (std::size_t i = 0; i < N; ++i) {
            derivatives_[i] = derivatives_[i] * right_value + left_value * rhs.derivatives_[i];
        }
        value_ = left_value * right_value;
        return *this;
    }

    constexpr Dual& operator/=(const Dual& rhs) {
        const T denominator = rhs.value_;
        require_nonzero(denominator); // Check before mutation: strong exception guarantee.
        const T quotient = value_ / denominator;
        for (std::size_t i = 0; i < N; ++i) {
            // q' = (u' - q*v')/v avoids squaring v or forming 1/v.
            // Do not update value_ early: rhs may alias *this.
            derivatives_[i] = (derivatives_[i] - quotient * rhs.derivatives_[i]) / denominator;
        }
        value_ = quotient;
        return *this;
    }

    // Scalar paths avoid constructing and propagating an N-lane zero gradient.
    // Scalars follow T's ordinary conversion rules; no mixed-precision promotion.
    constexpr Dual& operator+=(T rhs) noexcept {
        value_ += rhs;
        return *this;
    }

    constexpr Dual& operator-=(T rhs) noexcept {
        value_ -= rhs;
        return *this;
    }

    constexpr Dual& operator*=(T rhs) noexcept {
        value_ *= rhs;
        for (T& derivative : derivatives_) {
            derivative *= rhs;
        }
        return *this;
    }

    constexpr Dual& operator/=(T rhs) {
        require_nonzero(rhs);
        value_ /= rhs;
        for (T& derivative : derivatives_) {
            derivative /= rhs;
        }
        return *this;
    }

    [[nodiscard]] friend constexpr Dual operator+(Dual value) noexcept { return value; }

    [[nodiscard]] friend constexpr Dual operator-(Dual value) noexcept {
        value.value_ = -value.value_;
        for (T& derivative : value.derivatives_) {
            derivative = -derivative;
        }
        return value;
    }

    [[nodiscard]] friend constexpr Dual operator+(Dual lhs, const Dual& rhs) noexcept {
        lhs += rhs;
        return lhs;
    }

    [[nodiscard]] friend constexpr Dual operator-(Dual lhs, const Dual& rhs) noexcept {
        lhs -= rhs;
        return lhs;
    }

    [[nodiscard]] friend constexpr Dual operator*(Dual lhs, const Dual& rhs) noexcept {
        lhs *= rhs;
        return lhs;
    }

    [[nodiscard]] friend constexpr Dual operator/(Dual lhs, const Dual& rhs) {
        lhs /= rhs;
        return lhs;
    }

    [[nodiscard]] friend constexpr Dual operator+(Dual lhs, T rhs) noexcept {
        lhs += rhs;
        return lhs;
    }

    [[nodiscard]] friend constexpr Dual operator+(T lhs, Dual rhs) noexcept {
        rhs += lhs;
        return rhs;
    }

    [[nodiscard]] friend constexpr Dual operator-(Dual lhs, T rhs) noexcept {
        lhs -= rhs;
        return lhs;
    }

    [[nodiscard]] friend constexpr Dual operator-(T lhs, const Dual& rhs) noexcept {
        Dual result = -rhs;
        result.value_ = lhs - rhs.value_;
        return result;
    }

    [[nodiscard]] friend constexpr Dual operator*(Dual lhs, T rhs) noexcept {
        lhs *= rhs;
        return lhs;
    }

    [[nodiscard]] friend constexpr Dual operator*(T lhs, Dual rhs) noexcept {
        rhs *= lhs;
        return rhs;
    }

    [[nodiscard]] friend constexpr Dual operator/(Dual lhs, T rhs) {
        lhs /= rhs;
        return lhs;
    }

    [[nodiscard]] friend constexpr Dual operator/(T lhs, const Dual& rhs) {
        require_nonzero(rhs.value_);
        Dual result{lhs / rhs.value_};
        for (std::size_t i = 0; i < N; ++i) {
            result.derivatives_[i] = (-result.value_ * rhs.derivatives_[i]) / rhs.value_;
        }
        return result;
    }

private:
    static constexpr void require_nonzero(T denominator) {
        if (denominator == T{0}) { // Includes both +0 and -0; never clips small nonzero values.
            throw std::domain_error("mpmc::ad::Dual: division by zero");
        }
    }

    T value_{};
    Gradient derivatives_{};
};

} // namespace mpmc::ad

#endif // MPMC_AD_DUAL_HPP
