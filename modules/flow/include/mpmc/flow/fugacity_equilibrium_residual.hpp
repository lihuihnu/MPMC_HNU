#ifndef MPMC_FLOW_FUGACITY_EQUILIBRIUM_RESIDUAL_HPP
#define MPMC_FLOW_FUGACITY_EQUILIBRIUM_RESIDUAL_HPP

#include <mpmc/flow/natural_variable_cell_state.hpp>

#include <algorithm>
#include <array>
#include <cmath>
#include <concepts>
#include <cstddef>
#include <functional>
#include <limits>
#include <span>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <utility>
#include <vector>

namespace mpmc::flow {

/// Scalar-generic result returned by a configured phase fugacity-coefficient
/// evaluator.
///
/// The evaluator owns all thermodynamic-model-specific branch/family/root
/// choices. Flow only consumes ln(phi_i) for the exact supplied p, T and x.
template <typename Number>
struct PhaseLnFugacityCoefficientEvaluation {
    std::vector<Number> ln_phi;
};

/// Model-neutral callable contract for one fixed-three-phase fugacity
/// evaluation.
///
/// A concrete adapter may capture PR76 root indices, SW92 family/root metadata,
/// or another model-specific branch choice. The callable must evaluate the
/// exact phase slot at the exact pressure passed by flow and return one
/// ln(phi_i) value per component using the same scalar type.
///
/// This scalar-preserving signature is the differentiability boundary for a
/// future local/global Jacobian. It deliberately contains no finite-difference
/// fallback and performs no phase/root/family selection itself.
template <typename Evaluator, typename Number>
concept PhaseLnFugacityCoefficientEvaluator3P =
    requires(
        Evaluator& evaluator,
        PhaseSlot3 slot,
        const Number& pressure_pa,
        const Number& temperature_k,
        std::span<const Number> composition) {
        {
            std::invoke(
                evaluator,
                slot,
                pressure_pa,
                temperature_k,
                composition)
        } -> std::same_as<
            PhaseLnFugacityCoefficientEvaluation<Number>>;
    };

/// Non-owning scalar-generic local thermodynamic state used by the equilibrium
/// residual contract.
///
/// phase_pressures_pa contains the already-resolved actual pressure of every
/// phase. The residual never substitutes reference pressure for a non-reference
/// phase. The current natural-variable cell state helper produces equal
/// pressures because that earlier slice implements pc=none only; future
/// capillary-pressure code must resolve phase pressures before constructing
/// this view.
template <typename Number>
struct FugacityEquilibriumStateView3P {
    std::array<Number, 3> phase_pressures_pa;
    Number temperature_k;
    std::array<std::span<const Number>, 3>
        phase_compositions;
};

namespace fugacity_equilibrium_detail {

template <typename Number>
[[nodiscard]] inline long double
primal_value(const Number& value) {
    if constexpr (requires { value.value(); }) {
        return static_cast<long double>(
            value.value());
    } else {
        return static_cast<long double>(
            value);
    }
}

template <typename Number>
[[nodiscard]] inline long double
primal_epsilon() {
    if constexpr (requires(const Number& value) {
                      value.value();
                  }) {
        using Primal = std::remove_cvref_t<
            decltype(
                std::declval<const Number&>()
                    .value())>;
        static_assert(
            std::numeric_limits<Primal>::
                is_specialized);
        return static_cast<long double>(
            std::numeric_limits<Primal>::
                epsilon());
    } else {
        using Primal =
            std::remove_cvref_t<Number>;
        static_assert(
            std::numeric_limits<Primal>::
                is_specialized);
        return static_cast<long double>(
            std::numeric_limits<Primal>::
                epsilon());
    }
}

template <typename Number>
inline void require_finite_positive(
    const Number& value,
    const char* name) {
    const long double primal =
        primal_value(value);
    if (!std::isfinite(primal) ||
        !(primal > 0.0L)) {
        throw std::invalid_argument(
            std::string{"mpmc::flow: "} +
            name +
            " must be finite and strictly positive");
    }
}

template <typename Number>
inline void validate_state(
    const FugacityEquilibriumStateView3P<Number>& state) {
    require_finite_positive(
        state.temperature_k,
        "equilibrium temperature [K]");

    const std::size_t component_count =
        state.phase_compositions[0].size();
    if (component_count < 2U) {
        throw std::invalid_argument(
            "mpmc::flow: fugacity equilibrium requires at least two components");
    }

    const long double normalization_tolerance =
        256.0L *
        primal_epsilon<Number>() *
        static_cast<long double>(
            std::max<std::size_t>(
                component_count,
                1U));

    for (std::size_t phase = 0U;
         phase < fixed_three_phase_count;
         ++phase) {
        require_finite_positive(
            state.phase_pressures_pa[phase],
            "actual phase pressure [Pa]");
        const auto composition =
            state.phase_compositions[phase];
        if (composition.size() !=
            component_count) {
            throw std::invalid_argument(
                "mpmc::flow: phase composition sizes must match");
        }

        long double sum = 0.0L;
        for (const auto& fraction :
             composition) {
            const long double primal =
                primal_value(fraction);
            if (!std::isfinite(primal) ||
                !(primal > 0.0L)) {
                throw std::invalid_argument(
                    "mpmc::flow: fugacity equilibrium requires strictly positive finite phase compositions");
            }
            sum += primal;
        }
        if (!std::isfinite(sum) ||
            std::abs(sum - 1.0L) >
                normalization_tolerance) {
            throw std::invalid_argument(
                "mpmc::flow: fugacity equilibrium phase composition must be normalized");
        }
    }
}

template <typename Number>
inline void validate_ln_phi(
    const PhaseLnFugacityCoefficientEvaluation<Number>&
        evaluation,
    std::size_t component_count) {
    if (evaluation.ln_phi.size() !=
        component_count) {
        throw std::invalid_argument(
            "mpmc::flow: thermodynamic ln(phi) evaluator returned wrong component count");
    }
    for (const auto& value :
         evaluation.ln_phi) {
        if (!std::isfinite(
                primal_value(value))) {
            throw std::range_error(
                "mpmc::flow: thermodynamic ln(phi) evaluator returned non-finite value");
        }
    }
}

template <typename Number>
[[nodiscard]] inline auto natural_log(
    const Number& value) {
    using std::log;
    return log(value);
}

[[nodiscard]] inline std::size_t
local_row_index(
    PhaseSlot3 non_reference_phase,
    std::size_t component,
    std::size_t component_count) {
    const std::size_t phase =
        natural_variable_detail::
            checked_phase_index(
                non_reference_phase);
    if (phase == 0U) {
        throw std::invalid_argument(
            "mpmc::flow: reference phase has no fugacity-equilibrium residual block");
    }
    if (component >= component_count) {
        throw std::out_of_range(
            "mpmc::flow: fugacity-equilibrium component index out of range");
    }
    return (phase - 1U) *
               component_count +
           component;
}

} // namespace fugacity_equilibrium_detail

/// Compact 2*Nc local fugacity-equilibrium residual.
///
/// Row order is:
///   phase1-vs-phase0 components 0..Nc-1,
///   phase2-vs-phase0 components 0..Nc-1.
///
/// These local rows map directly to NaturalVariableLayout3P's final 2*Nc
/// equation rows. No conservation, mesh, flux, time-discretization or solver
/// quantity is included here.
template <typename Number>
class FugacityEquilibriumResidual3P {
public:
    FugacityEquilibriumResidual3P(
        std::size_t component_count,
        std::vector<Number> values)
        : component_count_(
              component_count),
          values_(std::move(values)) {
        if (component_count_ < 2U ||
            component_count_ >
                std::numeric_limits<std::size_t>::
                    max() /
                    2U ||
            values_.size() !=
                2U * component_count_) {
            throw std::invalid_argument(
                "mpmc::flow::FugacityEquilibriumResidual3P: invalid residual shape");
        }
    }

    [[nodiscard]] std::size_t
    component_count() const noexcept {
        return component_count_;
    }

    [[nodiscard]] std::size_t
    residual_count() const noexcept {
        return values_.size();
    }

    [[nodiscard]] std::span<const Number>
    values() const noexcept {
        return values_;
    }

    [[nodiscard]] std::size_t
    local_row_index(
        PhaseSlot3 non_reference_phase,
        std::size_t component) const {
        return fugacity_equilibrium_detail::
            local_row_index(
                non_reference_phase,
                component,
                component_count_);
    }

    [[nodiscard]] const Number&
    residual(
        PhaseSlot3 non_reference_phase,
        std::size_t component) const {
        return values_.at(
            local_row_index(
                non_reference_phase,
                component));
    }

private:
    std::size_t component_count_;
    std::vector<Number> values_;
};

/// Evaluate the fixed-three-phase positive-support equilibrium residual.
///
/// For alpha = phase1, phase2 and every component i:
///
///   R_(alpha,i) =
///       log(x_phase0,i / x_alpha,i)
///     + ln(phi_phase0,i) - ln(phi_alpha,i)
///     + log(p_phase0 / p_alpha)
///
/// which is exactly log(f_phase0,i / f_alpha,i) when
/// f_i = x_i * phi_i * p. Ratios keep logarithm arguments dimensionless.
///
/// The configured evaluator is called independently for phase0/1/2 using that
/// phase's actual pressure from state.phase_pressures_pa. Model-specific
/// root/family selection belongs inside the evaluator adapter. Exceptions from
/// the thermodynamic evaluator are intentionally propagated unchanged; flow
/// does not replace them with zero residuals, stale values or a fallback root.
///
/// Number is preserved through all arithmetic, allowing forward AD or another
/// scalar type to propagate derivatives through composition/pressure logarithms
/// and through a scalar-generic thermodynamic evaluator.
template <typename Number, typename Evaluator>
requires PhaseLnFugacityCoefficientEvaluator3P<
    Evaluator,
    Number>
[[nodiscard]] inline
FugacityEquilibriumResidual3P<Number>
evaluate_fugacity_equilibrium_residual_3p(
    const FugacityEquilibriumStateView3P<Number>&
        state,
    Evaluator&& evaluator) {
    fugacity_equilibrium_detail::
        validate_state(state);

    const std::size_t component_count =
        state.phase_compositions[0].size();
    if (component_count >
        std::numeric_limits<std::size_t>::
            max() /
            2U) {
        throw std::length_error(
            "mpmc::flow: fugacity residual size overflow");
    }

    std::array<
        PhaseLnFugacityCoefficientEvaluation<Number>,
        3>
        phase_evaluations;

    auto&& evaluator_ref = evaluator;
    for (std::size_t phase = 0U;
         phase < fixed_three_phase_count;
         ++phase) {
        const auto slot =
            static_cast<PhaseSlot3>(phase);
        phase_evaluations[phase] =
            std::invoke(
                evaluator_ref,
                slot,
                state.phase_pressures_pa[phase],
                state.temperature_k,
                state.phase_compositions[phase]);
        fugacity_equilibrium_detail::
            validate_ln_phi(
                phase_evaluations[phase],
                component_count);
    }

    std::vector<Number> values;
    values.reserve(
        2U * component_count);

    const auto reference_composition =
        state.phase_compositions[0];
    const auto& reference_ln_phi =
        phase_evaluations[0].ln_phi;
    const Number& reference_pressure =
        state.phase_pressures_pa[0];

    for (std::size_t phase = 1U;
         phase < fixed_three_phase_count;
         ++phase) {
        const auto composition =
            state.phase_compositions[phase];
        const auto& ln_phi =
            phase_evaluations[phase].ln_phi;
        const auto pressure_ratio =
            reference_pressure /
            state.phase_pressures_pa[phase];
        const auto log_pressure_ratio =
            fugacity_equilibrium_detail::
                natural_log(pressure_ratio);

        for (std::size_t component = 0U;
             component < component_count;
             ++component) {
            const auto composition_ratio =
                reference_composition[
                    component] /
                composition[component];
            values.push_back(
                fugacity_equilibrium_detail::
                    natural_log(
                        composition_ratio) +
                reference_ln_phi[component] -
                ln_phi[component] +
                log_pressure_ratio);
        }
    }

    return FugacityEquilibriumResidual3P<Number>{
        component_count,
        std::move(values)};
}

/// View the already validated natural-variable cell state as an equilibrium
/// input. Current phase pressures are equal because the current state contract
/// implements pc=none, but the residual itself does not assume equal pressure.
[[nodiscard]] inline
FugacityEquilibriumStateView3P<double>
make_fugacity_equilibrium_state_view(
    const NaturalVariableCellState3P& state) {
    return FugacityEquilibriumStateView3P<double>{
        state.phase_pressures_pa(),
        state.temperature_k(),
        {
            state.phase_composition(
                PhaseSlot3::phase0),
            state.phase_composition(
                PhaseSlot3::phase1),
            state.phase_composition(
                PhaseSlot3::phase2)}};
}

} // namespace mpmc::flow

#endif // MPMC_FLOW_FUGACITY_EQUILIBRIUM_RESIDUAL_HPP
