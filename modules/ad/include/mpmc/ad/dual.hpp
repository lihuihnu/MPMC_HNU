#ifndef MPMC_AD_DUAL_HPP
#define MPMC_AD_DUAL_HPP

#include <array>
#include <concepts>
#include <cstddef>
#include <stdexcept>
#include <type_traits>
#include <utility>

namespace mpmc::ad {

namespace detail {
struct DualTag final {};

template <typename T>
concept NestedDualScalar = requires {
    typename T::MpmcDualTag;
    typename T::BaseScalar;
} && std::same_as<typename T::MpmcDualTag, DualTag> &&
     std::floating_point<typename T::BaseScalar>;

template <typename T>
concept AdScalar =
    std::same_as<T, std::remove_cv_t<T>> &&
    (std::floating_point<T> || NestedDualScalar<T>);

template <AdScalar T, bool = std::floating_point<T>>
struct BaseScalarImpl;

template <AdScalar T>
struct BaseScalarImpl<T, true> { using type = T; };

template <AdScalar T>
struct BaseScalarImpl<T, false> { using type = typename T::BaseScalar; };

template <AdScalar T>
using base_scalar_t = typename BaseScalarImpl<T>::type;

template <AdScalar T>
[[nodiscard]] constexpr base_scalar_t<T> primal_value(const T& value) noexcept {
    if constexpr (std::floating_point<T>) {
        return value;
    } else {
        return primal_value(value.value());
    }
}
} // namespace detail

/// Fixed-width forward-mode dual number.
///
/// Scalar may be a built-in floating-point type or another mpmc::ad::Dual.
/// Recursively nesting Dual types provides higher-order/mixed derivatives while
/// preserving the same value/gradient ownership model and without a tape.
template <detail::AdScalar ScalarT, std::size_t N = 1>
    requires(N > 0)
class Dual {
public:
    using MpmcDualTag = detail::DualTag;
    using Scalar = ScalarT;
    using BaseScalar = detail::base_scalar_t<Scalar>;
    using Gradient = std::array<Scalar, N>;
    static constexpr std::size_t derivative_count = N;

    constexpr Dual() noexcept(std::is_nothrow_default_constructible_v<Scalar>) = default;

    explicit constexpr Dual(const Scalar& value)
        noexcept(std::is_nothrow_copy_constructible_v<Scalar>)
        : value_(value) {}

    explicit constexpr Dual(Scalar&& value)
        noexcept(std::is_nothrow_move_constructible_v<Scalar>)
        : value_(std::move(value)) {}

    template <typename U>
        requires(!std::same_as<std::remove_cvref_t<U>, Dual> &&
                 !std::same_as<std::remove_cvref_t<U>, Scalar> &&
                 std::constructible_from<Scalar, U>)
    explicit constexpr Dual(U&& value)
        noexcept(std::is_nothrow_constructible_v<Scalar, U>)
        : value_(std::forward<U>(value)) {}

    constexpr Dual(Scalar value, const Gradient& derivatives)
        noexcept(std::is_nothrow_move_constructible_v<Scalar> &&
                 std::is_nothrow_copy_constructible_v<Gradient>)
        : value_(std::move(value)), derivatives_(derivatives) {}

    template <typename U>
        requires std::constructible_from<Scalar, U>
    [[nodiscard]] static constexpr Dual variable(U&& value, std::size_t index) {
        Dual result{Scalar{std::forward<U>(value)}};
        result.derivatives_.at(index) = Scalar{BaseScalar{1}};
        return result;
    }

    [[nodiscard]] constexpr Scalar value() const
        noexcept(std::is_nothrow_copy_constructible_v<Scalar>) {
        return value_;
    }

    [[nodiscard]] constexpr Scalar derivative(std::size_t index) const {
        return derivatives_.at(index);
    }

    [[nodiscard]] constexpr const Gradient& derivatives() const & noexcept {
        return derivatives_;
    }

    [[nodiscard]] constexpr Gradient derivatives() const &&
        noexcept(std::is_nothrow_copy_constructible_v<Gradient>) {
        return derivatives_;
    }

    constexpr Dual& operator+=(const Dual& rhs)
        noexcept(noexcept(std::declval<Scalar&>() += std::declval<const Scalar&>())) {
        value_ += rhs.value_;
        for (std::size_t i = 0; i < N; ++i) {
            derivatives_[i] += rhs.derivatives_[i];
        }
        return *this;
    }

    constexpr Dual& operator-=(const Dual& rhs)
        noexcept(noexcept(std::declval<Scalar&>() -= std::declval<const Scalar&>())) {
        value_ -= rhs.value_;
        for (std::size_t i = 0; i < N; ++i) {
            derivatives_[i] -= rhs.derivatives_[i];
        }
        return *this;
    }

    constexpr Dual& operator*=(const Dual& rhs)
        noexcept(noexcept(std::declval<Scalar>() * std::declval<Scalar>()) &&
                 noexcept(std::declval<Scalar&>() = std::declval<Scalar>())) {
        const Scalar left_value = value_;
        const Scalar right_value = rhs.value_;
        for (std::size_t i = 0; i < N; ++i) {
            derivatives_[i] =
                derivatives_[i] * right_value +
                left_value * rhs.derivatives_[i];
        }
        value_ = left_value * right_value;
        return *this;
    }

    constexpr Dual& operator/=(const Dual& rhs) {
        const Scalar denominator = rhs.value_;
        require_nonzero(denominator);
        const Scalar quotient = value_ / denominator;
        for (std::size_t i = 0; i < N; ++i) {
            derivatives_[i] =
                (derivatives_[i] -
                 quotient * rhs.derivatives_[i]) /
                denominator;
        }
        value_ = quotient;
        return *this;
    }

    template <typename U>
        requires(!std::same_as<std::remove_cvref_t<U>, Dual> &&
                 std::constructible_from<Scalar, U>)
    constexpr Dual& operator+=(U&& rhs) {
        value_ += Scalar{std::forward<U>(rhs)};
        return *this;
    }

    template <typename U>
        requires(!std::same_as<std::remove_cvref_t<U>, Dual> &&
                 std::constructible_from<Scalar, U>)
    constexpr Dual& operator-=(U&& rhs) {
        value_ -= Scalar{std::forward<U>(rhs)};
        return *this;
    }

    template <typename U>
        requires(!std::same_as<std::remove_cvref_t<U>, Dual> &&
                 std::constructible_from<Scalar, U>)
    constexpr Dual& operator*=(U&& rhs) {
        const Scalar scalar{std::forward<U>(rhs)};
        value_ *= scalar;
        for (auto& derivative : derivatives_) {
            derivative *= scalar;
        }
        return *this;
    }

    template <typename U>
        requires(!std::same_as<std::remove_cvref_t<U>, Dual> &&
                 std::constructible_from<Scalar, U>)
    constexpr Dual& operator/=(U&& rhs) {
        const Scalar scalar{std::forward<U>(rhs)};
        require_nonzero(scalar);
        value_ /= scalar;
        for (auto& derivative : derivatives_) {
            derivative /= scalar;
        }
        return *this;
    }

    [[nodiscard]] friend constexpr Dual operator+(Dual value) noexcept {
        return value;
    }

    [[nodiscard]] friend constexpr Dual operator-(Dual value)
        noexcept(noexcept(-std::declval<Scalar>())) {
        value.value_ = -value.value_;
        for (auto& derivative : value.derivatives_) {
            derivative = -derivative;
        }
        return value;
    }

    [[nodiscard]] friend constexpr Dual operator+(
        Dual lhs, const Dual& rhs)
        noexcept(noexcept(lhs += rhs)) {
        lhs += rhs;
        return lhs;
    }

    [[nodiscard]] friend constexpr Dual operator-(
        Dual lhs, const Dual& rhs)
        noexcept(noexcept(lhs -= rhs)) {
        lhs -= rhs;
        return lhs;
    }

    [[nodiscard]] friend constexpr Dual operator*(
        Dual lhs, const Dual& rhs)
        noexcept(noexcept(lhs *= rhs)) {
        lhs *= rhs;
        return lhs;
    }

    [[nodiscard]] friend constexpr Dual operator/(
        Dual lhs, const Dual& rhs) {
        lhs /= rhs;
        return lhs;
    }

    template <typename U>
        requires(!std::same_as<std::remove_cvref_t<U>, Dual> &&
                 std::constructible_from<Scalar, U>)
    [[nodiscard]] friend constexpr Dual operator+(Dual lhs, U&& rhs) {
        lhs += std::forward<U>(rhs);
        return lhs;
    }

    template <typename U>
        requires(!std::same_as<std::remove_cvref_t<U>, Dual> &&
                 std::constructible_from<Scalar, U>)
    [[nodiscard]] friend constexpr Dual operator+(U&& lhs, Dual rhs) {
        rhs += std::forward<U>(lhs);
        return rhs;
    }

    template <typename U>
        requires(!std::same_as<std::remove_cvref_t<U>, Dual> &&
                 std::constructible_from<Scalar, U>)
    [[nodiscard]] friend constexpr Dual operator-(Dual lhs, U&& rhs) {
        lhs -= std::forward<U>(rhs);
        return lhs;
    }

    template <typename U>
        requires(!std::same_as<std::remove_cvref_t<U>, Dual> &&
                 std::constructible_from<Scalar, U>)
    [[nodiscard]] friend constexpr Dual operator-(U&& lhs, const Dual& rhs) {
        return Dual{std::forward<U>(lhs)} - rhs;
    }

    template <typename U>
        requires(!std::same_as<std::remove_cvref_t<U>, Dual> &&
                 std::constructible_from<Scalar, U>)
    [[nodiscard]] friend constexpr Dual operator*(Dual lhs, U&& rhs) {
        lhs *= std::forward<U>(rhs);
        return lhs;
    }

    template <typename U>
        requires(!std::same_as<std::remove_cvref_t<U>, Dual> &&
                 std::constructible_from<Scalar, U>)
    [[nodiscard]] friend constexpr Dual operator*(U&& lhs, Dual rhs) {
        rhs *= std::forward<U>(lhs);
        return rhs;
    }

    template <typename U>
        requires(!std::same_as<std::remove_cvref_t<U>, Dual> &&
                 std::constructible_from<Scalar, U>)
    [[nodiscard]] friend constexpr Dual operator/(Dual lhs, U&& rhs) {
        lhs /= std::forward<U>(rhs);
        return lhs;
    }

    template <typename U>
        requires(!std::same_as<std::remove_cvref_t<U>, Dual> &&
                 std::constructible_from<Scalar, U>)
    [[nodiscard]] friend constexpr Dual operator/(U&& lhs, const Dual& rhs) {
        return Dual{std::forward<U>(lhs)} / rhs;
    }

private:
    static constexpr void require_nonzero(const Scalar& denominator) {
        if (detail::primal_value(denominator) == BaseScalar{0}) {
            throw std::domain_error("mpmc::ad::Dual: division by zero");
        }
    }

    Scalar value_{};
    Gradient derivatives_{};
};

} // namespace mpmc::ad

#endif // MPMC_AD_DUAL_HPP
