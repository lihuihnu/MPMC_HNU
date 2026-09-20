#ifndef MPMC_FLOW_SATURATION_CONSTITUTIVE_HPP
#define MPMC_FLOW_SATURATION_CONSTITUTIVE_HPP

#include <mpmc/flow/natural_variable_cell_state.hpp>

#include <array>
#include <cmath>
#include <concepts>
#include <cstddef>
#include <functional>
#include <stdexcept>
#include <string>
#include <string_view>
#include <type_traits>
#include <utility>

namespace mpmc::flow {

inline constexpr std::string_view
    three_phase_saturation_constitutive_convention =
        "flow/three-phase-saturation-constitutive/fixed-positive-support/v1";

/// Scalar-generic positive-support saturation chart.
///
/// Coordinates are exactly the current natural-variable saturation coordinates:
///   q_S = (S0, S1)
///   S2  = 1 - S0 - S1.
///
/// No clipping, effective-saturation transform, residual saturation or
/// hysteresis state is implied at this layer. A concrete constitutive model may
/// capture those definitions explicitly in its own configured evaluator.
template <typename Number>
struct ThreePhaseSaturationState3P {
    std::array<Number, 3> saturation;
};

template <typename Number>
struct RelativePermeabilityEvaluation3P {
    std::array<Number, 3> relative_permeability;
};

/// Capillary-pressure offsets relative to phase0.
///
/// The sign convention is fixed:
///   offset_phase1_pa = p_phase1 - p_phase0
///   offset_phase2_pa = p_phase2 - p_phase0.
///
/// phase0 is therefore the explicit pressure reference and has identically zero
/// capillary offset, including zero derivative.
template <typename Number>
struct CapillaryPressureOffsetsEvaluation3P {
    std::array<Number, 2> non_reference_offset_pa;
};

template <typename Number>
struct ThreePhaseSaturationConstitutiveEvaluation3P {
    static constexpr std::string_view convention =
        three_phase_saturation_constitutive_convention;

    ThreePhaseSaturationState3P<Number> saturation_state;
    std::array<Number, 3> relative_permeability;
    std::array<Number, 3> capillary_pressure_offset_pa;
    std::array<Number, 3> phase_pressure_pa;
};

template <typename Evaluator, typename Number>
concept RelativePermeabilityEvaluator3P =
    requires(
        Evaluator& evaluator,
        const ThreePhaseSaturationState3P<Number>& state) {
        {
            std::invoke(evaluator, state)
        } -> std::same_as<
            RelativePermeabilityEvaluation3P<Number>>;
    };

template <typename Evaluator, typename Number>
concept CapillaryPressureEvaluator3P =
    requires(
        Evaluator& evaluator,
        const ThreePhaseSaturationState3P<Number>& state) {
        {
            std::invoke(evaluator, state)
        } -> std::same_as<
            CapillaryPressureOffsetsEvaluation3P<Number>>;
    };

/// Explicit capillary-disabled model.
///
/// This is the only built-in capillary law in the first contract. It publishes
/// exactly zero pressure offsets and zero saturation derivatives. A nonzero
/// capillary model must be selected/configured explicitly.
struct NoCapillaryPressure3P {
    template <typename Number>
    [[nodiscard]]
    CapillaryPressureOffsetsEvaluation3P<Number>
    operator()(
        const ThreePhaseSaturationState3P<Number>&) const {
        return {
            std::array<Number, 2>{
                Number{},
                Number{}}};
    }
};

namespace saturation_constitutive_detail {

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
inline void require_finite_nonnegative(
    const Number& value,
    const char* name) {
    const long double primal =
        primal_value(value);
    if (!std::isfinite(primal) ||
        primal < 0.0L) {
        throw std::invalid_argument(
            std::string{"mpmc::flow: "} +
            name +
            " must be finite and nonnegative");
    }
}

template <typename Number>
inline void require_finite(
    const Number& value,
    const char* name) {
    if (!std::isfinite(
            primal_value(value))) {
        throw std::invalid_argument(
            std::string{"mpmc::flow: "} +
            name +
            " must be finite");
    }
}

} // namespace saturation_constitutive_detail

template <typename Number>
[[nodiscard]] inline
ThreePhaseSaturationState3P<Number>
make_three_phase_saturation_state(
    const Number& saturation0,
    const Number& saturation1) {
    saturation_constitutive_detail::
        require_finite_positive(
            saturation0,
            "phase0 saturation");
    saturation_constitutive_detail::
        require_finite_positive(
            saturation1,
            "phase1 saturation");

    const Number saturation2 =
        Number{1.0} -
        saturation0 -
        saturation1;
    saturation_constitutive_detail::
        require_finite_positive(
            saturation2,
            "reconstructed phase2 saturation");

    return {
        std::array<Number, 3>{
            saturation0,
            saturation1,
            saturation2}};
}

/// Evaluate model-neutral three-phase saturation constitutive data.
///
/// Concrete evaluators may capture residual saturations, endpoints, tables,
/// wettability/physical-role maps and provenance. This wrapper never infers
/// oil/gas/water identity from numerical phase slots.
///
/// The scalar type is preserved from S0/S1 through kr, capillary offsets and
/// resolved phase pressures, so an analytic scalar-generic law can propagate
/// derivatives directly with forward AD.
///
/// Generic validation requires only finite nonnegative kr. Any stricter upper
/// bound (for example kr <= 1) belongs to the configured model's declared
/// normalization/bounds contract and must not be imposed silently here.
template <
    typename Number,
    typename RelativePermeabilityEvaluator,
    typename CapillaryPressureEvaluator>
requires
    RelativePermeabilityEvaluator3P<
        RelativePermeabilityEvaluator,
        Number> &&
    CapillaryPressureEvaluator3P<
        CapillaryPressureEvaluator,
        Number>
[[nodiscard]] inline
ThreePhaseSaturationConstitutiveEvaluation3P<Number>
evaluate_three_phase_saturation_constitutive(
    const Number& reference_pressure_pa,
    const Number& saturation0,
    const Number& saturation1,
    RelativePermeabilityEvaluator&&
        relative_permeability_evaluator,
    CapillaryPressureEvaluator&&
        capillary_pressure_evaluator) {
    saturation_constitutive_detail::
        require_finite_positive(
            reference_pressure_pa,
            "reference phase pressure [Pa]");

    auto state =
        make_three_phase_saturation_state(
            saturation0,
            saturation1);

    auto&& kr_evaluator =
        relative_permeability_evaluator;
    auto&& pc_evaluator =
        capillary_pressure_evaluator;

    auto kr =
        std::invoke(
            kr_evaluator,
            state);
    auto pc =
        std::invoke(
            pc_evaluator,
            state);

    for (const auto& value :
         kr.relative_permeability) {
        saturation_constitutive_detail::
            require_finite_nonnegative(
                value,
                "relative permeability");
    }
    for (const auto& value :
         pc.non_reference_offset_pa) {
        saturation_constitutive_detail::
            require_finite(
                value,
                "capillary pressure offset [Pa]");
    }

    const std::array<Number, 3>
        capillary_offsets{
            Number{},
            pc.non_reference_offset_pa[0],
            pc.non_reference_offset_pa[1]};

    const std::array<Number, 3>
        phase_pressures{
            reference_pressure_pa,
            reference_pressure_pa +
                capillary_offsets[1],
            reference_pressure_pa +
                capillary_offsets[2]};

    for (const auto& pressure :
         phase_pressures) {
        saturation_constitutive_detail::
            require_finite_positive(
                pressure,
                "resolved phase pressure [Pa]");
    }

    return {
        std::move(state),
        std::move(kr.relative_permeability),
        capillary_offsets,
        phase_pressures};
}

/// Convenience overload for the current validated fixed-three-phase cell
/// saturation state. It does not mutate the cell state or change its existing
/// pc=none phase-pressure helper.
template <
    typename RelativePermeabilityEvaluator,
    typename CapillaryPressureEvaluator>
requires
    RelativePermeabilityEvaluator3P<
        RelativePermeabilityEvaluator,
        double> &&
    CapillaryPressureEvaluator3P<
        CapillaryPressureEvaluator,
        double>
[[nodiscard]] inline
ThreePhaseSaturationConstitutiveEvaluation3P<double>
evaluate_three_phase_saturation_constitutive(
    const NaturalVariableCellState3P& state,
    RelativePermeabilityEvaluator&&
        relative_permeability_evaluator,
    CapillaryPressureEvaluator&&
        capillary_pressure_evaluator) {
    return evaluate_three_phase_saturation_constitutive(
        state.reference_pressure_pa(),
        state.phase_saturation(
            PhaseSlot3::phase0),
        state.phase_saturation(
            PhaseSlot3::phase1),
        std::forward<
            RelativePermeabilityEvaluator>(
                relative_permeability_evaluator),
        std::forward<
            CapillaryPressureEvaluator>(
                capillary_pressure_evaluator));
}

} // namespace mpmc::flow

#endif // MPMC_FLOW_SATURATION_CONSTITUTIVE_HPP
