#ifndef MPMC_THERMODYNAMICS_SW92_PURE_HPP
#define MPMC_THERMODYNAMICS_SW92_PURE_HPP

#include <mpmc/thermodynamics/sw92_correlations.hpp>
#include <mpmc/thermodynamics/sw92_parameters.hpp>

#include <cmath>
#include <concepts>
#include <cstddef>
#include <limits>
#include <stdexcept>
#include <string_view>
#include <type_traits>

namespace mpmc::thermodynamics {

inline constexpr std::string_view sw92_pure_convention =
    "SW92/corrected-original/PR76-printed-coefficients/R-SI-2019/pure-v1";

template <typename Number>
struct Sw92PureValues { Number a, b, alpha; };

namespace detail {
template <typename Number, std::floating_point T>
    requires Sw92Number<Number, T>
void sw92_check_state(
    const Number& t, T molality, const Sw92Applicability& app) {
    const T temperature = sw92_value(t);
    if (!sw92_finite(t) || !(temperature > T{0}) ||
        !std::isfinite(molality) || molality < T{0})
        throw std::domain_error(
            "SW92: finite T/seeds, T>0 K and NaCl molality>=0 mol/kg H2O required");
    if (const auto& b = app.state.temperature_k; b &&
        (static_cast<long double>(temperature) < b->lower ||
         static_cast<long double>(temperature) > b->upper))
        throw std::domain_error("SW92: temperature outside declared dataset interval");
    if (const auto& b = app.nacl_molality_mol_per_kg_water; b &&
        (static_cast<long double>(molality) < b->lower ||
         static_cast<long double>(molality) > b->upper))
        throw std::domain_error("SW92: NaCl molality outside declared dataset interval");
}

template <std::floating_point T>
[[nodiscard]] T sw92_cast(double v) {
    const long double limit = static_cast<long double>(std::numeric_limits<T>::max());
    if (!std::isfinite(v) || static_cast<long double>(v) > limit ||
        static_cast<long double>(v) < -limit)
        throw std::range_error("SW92: parameter outside scalar range");
    const T out = static_cast<T>(v);
    if (out == T{0} && v != 0.0) throw std::range_error("SW92: parameter underflow");
    return out;
}
} // namespace detail

/// Prepared PR pure-component coefficients. Non-water uses PR76 alpha; water
/// uses SW92 Eq.(9). The prepared scalar type remains built-in floating point;
/// evaluate accepts the same structural AD numbers as the PR76 thermodynamic
/// kernel. Prescribed NaCl molality is intentionally not an AD coordinate.
template <std::floating_point T = double>
    requires std::same_as<T, std::remove_cv_t<T>>
class Sw92Pure {
public:
    Sw92Pure(const Sw92Pure&) = default;
    Sw92Pure(Sw92Pure&&) noexcept = default;
    Sw92Pure& operator=(const Sw92Pure&) = delete;
    Sw92Pure& operator=(Sw92Pure&&) = delete;

    [[nodiscard]] static Sw92Pure from_parameters(const Sw92ParameterSet& p,
                                                  std::size_t i) {
        if (i >= p.pure_records().size())
            throw std::out_of_range("Sw92Pure: component index outside snapshot");
        return Sw92Pure(p, i);
    }
    [[nodiscard]] static constexpr T gas_constant() noexcept {
        return static_cast<T>(8.31446261815324L);
    }
    [[nodiscard]] T critical_temperature_k() const noexcept { return tc_; }
    [[nodiscard]] T covolume() const noexcept { return b_; }
    [[nodiscard]] Sw92Species species() const noexcept { return species_; }

    template <typename Number>
        requires detail::Sw92Number<Number, T>
    [[nodiscard]] Sw92PureValues<Number> evaluate(
        const Number& temperature_k,
        T nacl_molality_mol_per_kg_water) const {
        detail::sw92_check_state(
            temperature_k, nacl_molality_mol_per_kg_water, app_);
        Number alpha{T{0}};
        if (species_ == Sw92Species::water) {
            alpha = sw92_water_alpha<Number, T>(
                temperature_k, tc_, nacl_molality_mol_per_kg_water);
        } else {
            using std::sqrt;
            const Number q = T{1} + kappa_ *
                (T{1} - sqrt(temperature_k) / sqrt_tc_);
            alpha = q * q;
            if (!detail::sw92_finite(alpha))
                throw std::range_error("Sw92Pure: non-water alpha");
        }
        const Number a = ac_ * alpha;
        if (!detail::sw92_finite(a) || detail::sw92_value(a) < T{0})
            throw std::range_error("Sw92Pure: attraction");
        return {a, Number{b_}, alpha};
    }

private:
    Sw92Pure(const Sw92ParameterSet& p, std::size_t i)
        : species_(p.species(i)), app_(p.applicability()) {
        const auto& r = p.pure_records()[i];
        tc_ = detail::sw92_cast<T>(r.critical_temperature->value);
        const T pc = detail::sw92_cast<T>(r.critical_pressure->value);
        const T omega = detail::sw92_cast<T>(r.acentric_factor->value);
        sqrt_tc_ = std::sqrt(tc_);
        const T rtc = gas_constant()*tc_;
        ac_ = (static_cast<T>(0.45724L)*rtc)*(rtc/pc);
        b_ = (static_cast<T>(0.07780L)*rtc)/pc;
        kappa_ = species_ == Sw92Species::water ? T{0} :
            static_cast<T>(0.37464L) + omega*(static_cast<T>(1.54226L) -
                                               static_cast<T>(0.26992L)*omega);
        if (!std::isfinite(sqrt_tc_) || !(sqrt_tc_ > T{0}) ||
            !std::isfinite(ac_) || !(ac_ > T{0}) || !std::isfinite(b_) ||
            !(b_ > T{0}) || !std::isfinite(kappa_))
            throw std::range_error(
                "Sw92Pure: prepared coefficient outside scalar range");
    }

    Sw92Species species_{Sw92Species::unspecified};
    Sw92Applicability app_;
    T tc_{}, sqrt_tc_{}, ac_{}, b_{}, kappa_{};
};

} // namespace mpmc::thermodynamics

#endif // MPMC_THERMODYNAMICS_SW92_PURE_HPP
