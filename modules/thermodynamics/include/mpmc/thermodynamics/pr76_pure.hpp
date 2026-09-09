#ifndef MPMC_THERMODYNAMICS_PR76_PURE_HPP
#define MPMC_THERMODYNAMICS_PR76_PURE_HPP

#include <mpmc/thermodynamics/pr_parameters.hpp>

#include <cmath>
#include <concepts>
#include <cstddef>
#include <limits>
#include <span>
#include <stdexcept>
#include <string>
#include <string_view>
#include <type_traits>
#include <utility>

namespace mpmc::thermodynamics {

/// Implementation convention, separate from the parameter dataset revision.
/// Printed PR76 coefficients, with modern SI R (not a value stated in the paper).
inline constexpr std::string_view pr76_pure_convention = "PR76/printed-coefficients/R-SI-2019";

/// Molar SI coefficients; derivatives have these units divided by the seed unit.
/// b is returned as a constant of the SAME number type, with zero AD derivatives.
template <typename Number>
struct Pr76PureValues {
    Number a; // Pa m^6 mol^-2, NOT the dimensionless EOS coefficient A.
    Number b; // m^3 mol^-1, NOT the dimensionless EOS coefficient B.
};

namespace detail {
// Structural scalar adapter: no dependency on AD or a third-party number library.
// Supported/tested numbers are T and mpmc::ad::Dual<T,N>. Other adapters are not
// promised; they must also supply AD-preserving arithmetic, construction and sqrt.
template <typename Number, typename T>
concept Pr76Number = std::same_as<Number, T> || requires(const Number& number) {
    typename Number::Scalar;
    requires std::same_as<typename Number::Scalar, T>;
    { number.value() } -> std::same_as<T>;
    { std::span<const T>{number.derivatives()} };
};

template <typename Number>
[[nodiscard]] auto pr76_value(const Number& number) {
    if constexpr (std::floating_point<Number>) {
        return number;
    } else {
        return number.value();
    }
}

template <typename Number>
[[nodiscard]] bool pr76_finite(const Number& number) {
    if (!std::isfinite(pr76_value(number))) {
        return false;
    }
    if constexpr (!std::floating_point<Number>) {
        for (const auto derivative : number.derivatives()) {
            if (!std::isfinite(derivative)) {
                return false;
            }
        }
    }
    return true;
}
} // namespace detail

/// Prepared coefficients for ONE entry of an existing, validated PR76 snapshot.
/// Configuration copies provenance and precomputes a(Tc), b and kappa; evaluate
/// has no allocation or mutable cache. Copy/move construction is allowed, but no
/// partial multi-field assignment. Retain a live, non-moved-from object for use.
///
/// Source: Peng & Robinson (1976), journal p.60, Eqs.(9),(10),(12),(13),(17),(18).
/// Only temperature is differentiated; Tc/Pc/omega are constant model data.
/// See pr76_pure.md for the source digest, extrapolation and numerical limits.
template <std::floating_point T = double>
    requires std::same_as<T, std::remove_cv_t<T>>
class Pr76Pure {
public:
    Pr76Pure(const Pr76Pure&) = default;
    Pr76Pure(Pr76Pure&&) noexcept = default;
    Pr76Pure& operator=(const Pr76Pure&) = delete;
    Pr76Pure& operator=(Pr76Pure&&) = delete;

    [[nodiscard]] static Pr76Pure from_parameters(const PrParameterSet& parameters,
                                                 std::size_t component_index) {
        if (component_index >= parameters.pure_records().size()) {
            throw std::out_of_range("Pr76Pure: component index outside snapshot order");
        }
        return Pr76Pure(parameters, component_index);
    }

    /// R = NA*kB in J/(mol K), exact decimal SI value rounded to the selected T.
    [[nodiscard]] static constexpr T gas_constant() noexcept {
        return static_cast<T>(8.31446261815324L);
    }
    [[nodiscard]] T critical_temperature_k() const noexcept { return tc_; }
    [[nodiscard]] T critical_attraction() const noexcept { return ac_; }
    [[nodiscard]] T covolume() const noexcept { return b_; }
    [[nodiscard]] T kappa() const noexcept { return kappa_; }
    [[nodiscard]] const PrPureRecord& source_record() const & noexcept { return record_; }
    const PrPureRecord& source_record() const && = delete;
    [[nodiscard]] const std::string& dataset_id() const & noexcept { return dataset_id_; }
    const std::string& dataset_id() const && = delete;
    [[nodiscard]] const std::string& revision() const & noexcept { return revision_; }
    const std::string& revision() const && = delete;
    [[nodiscard]] const Applicability& applicability() const & noexcept { return applicability_; }
    const Applicability& applicability() const && = delete;

    /// T_kelvin > 0, finite, with finite seeds. A declared temperature interval
    /// is enforced without clipping (inclusive endpoints); absent bounds remain
    /// UNKNOWN, not scientifically validated. Pressure is not an input here.
    ///
    /// alpha = [1+kappa*(1-sqrt(T/Tc))]^2. For T outside the paper's calibration
    /// interval this is an explicit algebraic continuation, not a validity claim;
    /// in particular a negative bracket is squared, never clipped or abs-ed.
    ///
    /// domain_error: invalid temperature/seed or known temperature-bound violation.
    /// range_error: nonfinite intermediate/result or positive attraction lost to zero.
    /// AD math exceptions propagate. No finite-difference or manual AD slope injection.
    template <typename Number>
        requires detail::Pr76Number<Number, T>
    [[nodiscard]] Pr76PureValues<Number> evaluate(const Number& temperature_k) const {
        const T temperature = detail::pr76_value(temperature_k);
        if (!detail::pr76_finite(temperature_k) || !(temperature > T{0})) {
            throw std::domain_error("Pr76Pure: temperature and seeds must be finite, T > 0 K");
        }
        const auto& bounds = applicability_.temperature_k;
        if (bounds && (static_cast<long double>(temperature) < bounds->lower ||
                       static_cast<long double>(temperature) > bounds->upper)) {
            throw std::domain_error("Pr76Pure: temperature outside declared dataset interval");
        }
        using std::sqrt; // ADL finds the caller's AD sqrt; never strip its derivatives.
        // This algebraic form avoids forming a possibly overflowing/underflowing T/Tc.
        const Number reduced_root = sqrt(temperature_k) / sqrt_tc_;
        const Number bracket = T{1} + kappa_ * (T{1} - reduced_root);
        if (!detail::pr76_finite(reduced_root) || !detail::pr76_finite(bracket)) {
            throw std::range_error("Pr76Pure: nonrepresentable temperature expression");
        }
        // Do not square the bracket first: a small ac can rescale a large bracket.
        const Number attraction = (bracket * ac_) * bracket;
        if (!detail::pr76_finite(attraction) ||
            (detail::pr76_value(attraction) == T{0} && detail::pr76_value(bracket) != T{0})) {
            throw std::range_error("Pr76Pure: nonrepresentable attraction or derivative");
        }
        return {attraction, Number{b_}}; // Eq.(13): b has no temperature dependence.
    }

private:
    [[nodiscard]] static T parameter_cast(double value) {
        const long double limit = static_cast<long double>(std::numeric_limits<T>::max());
        if (!std::isfinite(value) || static_cast<long double>(value) > limit ||
            static_cast<long double>(value) < -limit) {
            throw std::range_error("Pr76Pure: parameter outside selected scalar range");
        }
        const T converted = static_cast<T>(value);
        if (converted == T{0} && value != 0.0) {
            throw std::range_error("Pr76Pure: nonzero parameter underflows selected scalar");
        }
        return converted;
    }

    Pr76Pure(const PrParameterSet& parameters, std::size_t index)
        : record_(parameters.pure_records()[index]), dataset_id_(parameters.dataset_id()),
          revision_(parameters.revision()), applicability_(parameters.applicability()) {
        tc_ = parameter_cast(record_.critical_temperature->value);
        const T pc = parameter_cast(record_.critical_pressure->value);
        const T omega = parameter_cast(record_.acentric_factor->value);
        sqrt_tc_ = std::sqrt(tc_);
        // Eq.(18), Horner order. Printed decimal coefficients, NOT a PR78 branch.
        kappa_ = static_cast<T>(0.37464L) + omega *
                 (static_cast<T>(1.54226L) - static_cast<T>(0.26992L) * omega);
        const T rtc = gas_constant() * tc_;
        // Eq.(9),(10): split multiplication instead of eagerly forming (R*Tc)^2.
        ac_ = (static_cast<T>(0.45724L) * rtc) * (rtc / pc);
        b_ = (static_cast<T>(0.07780L) * rtc) / pc;
        if (!std::isfinite(kappa_) || !std::isfinite(ac_) || !(ac_ > T{0}) ||
            !std::isfinite(b_) || !(b_ > T{0})) {
            throw std::range_error("Pr76Pure: precomputed coefficient outside scalar range");
        }
    }

    PrPureRecord record_;
    std::string dataset_id_, revision_;
    Applicability applicability_;
    T tc_{}, sqrt_tc_{}, ac_{}, b_{}, kappa_{};
};

} // namespace mpmc::thermodynamics
#endif // MPMC_THERMODYNAMICS_PR76_PURE_HPP
