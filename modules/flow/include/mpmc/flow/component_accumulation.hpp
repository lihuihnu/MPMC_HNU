#ifndef MPMC_FLOW_COMPONENT_ACCUMULATION_HPP
#define MPMC_FLOW_COMPONENT_ACCUMULATION_HPP

#include <mpmc/flow/natural_variable_cell_state.hpp>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <limits>
#include <span>
#include <stdexcept>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

namespace mpmc::flow {

inline constexpr std::string_view
    pore_volume_component_accumulation_convention =
        "flow/pore-volume-component-accumulation/fixed-three-phase-positive-support/v1";

/// Phase molar-density primal and derivative payload on one frozen
/// NaturalVariableLayout3P chart.
///
/// gradient[phase][q] = dc_phase / dq, where c is [mol / phase-fluid m^3].
/// The gradient must contain every natural-variable column, including zeros.
struct PhaseMolarDensityNaturalVariableLinearization3P {
    NaturalVariableLayout3P layout;
    std::array<double, 3>
        molar_density_mol_per_m3{};
    std::array<std::vector<double>, 3>
        gradient;
};

/// Component accumulation per bulk porous-medium volume.
///
/// N_i = phi * sum_phase S_phase c_phase x_phase,i
///
/// component_accumulation_mol_per_bulk_m3 retains canonical component order.
/// This quantity is not the physics-layer total-fluid-volume inventory.
struct PoreVolumeComponentAccumulationSnapshot3P {
    static constexpr std::string_view convention =
        pore_volume_component_accumulation_convention;

    double porosity{};
    std::vector<std::string> component_ids;
    std::vector<double>
        component_accumulation_mol_per_bulk_m3;
    double total_accumulation_mol_per_bulk_m3{};

    [[nodiscard]] std::size_t
    component_count() const noexcept {
        return component_ids.size();
    }

    [[nodiscard]] double component(
        std::size_t component_index) const {
        return component_accumulation_mol_per_bulk_m3.at(
            component_index);
    }
};

/// Current-state accumulation Jacobian on the exact frozen natural-variable
/// chart used to evaluate the current accumulation.
///
/// component_jacobian[row * input_count + column]
///     = d N_row / d q_column.
///
/// Row identity is canonical component identity. Column identity is supplied by
/// layout, including the explicit per-phase composition pivot.
struct PoreVolumeComponentAccumulationLinearization3P {
    static constexpr std::string_view convention =
        pore_volume_component_accumulation_convention;

    NaturalVariableLayout3P layout;
    double porosity{};
    std::vector<std::string> component_ids;
    std::size_t input_count{};
    std::vector<double> component_jacobian;
    std::vector<double> total_accumulation_gradient;

    [[nodiscard]] double d_component(
        std::size_t component,
        std::size_t column) const {
        if (component >= component_ids.size() ||
            column >= input_count ||
            (input_count != 0U &&
             component >
                 (std::numeric_limits<std::size_t>::max() -
                  column) /
                     input_count)) {
            throw std::out_of_range(
                "mpmc::flow: component accumulation Jacobian index out of range");
        }
        return component_jacobian.at(
            component * input_count +
            column);
    }

    [[nodiscard]] double d_total(
        std::size_t column) const {
        return total_accumulation_gradient.at(
            column);
    }
};

/// Owned current/previous accumulation pair.
///
/// No time-step size, time derivative or residual is implied. The pair exists
/// only to preserve the two state snapshots needed by a later time-discrete
/// conservation law.
struct PoreVolumeComponentAccumulationPair3P {
    PoreVolumeComponentAccumulationSnapshot3P current;
    PoreVolumeComponentAccumulationSnapshot3P previous;
};

namespace component_accumulation_detail {

[[nodiscard]] inline bool near_roundoff(
    double first,
    double second) {
    if (!std::isfinite(first) ||
        !std::isfinite(second)) {
        return false;
    }
    const double scale =
        std::max(
            {1.0,
             std::abs(first),
             std::abs(second)});
    return std::abs(first - second) <=
        4096.0 *
            std::numeric_limits<double>::
                epsilon() *
            scale;
}

inline void validate_porosity(
    double porosity) {
    if (!std::isfinite(porosity) ||
        !(porosity > 0.0) ||
        porosity > 1.0) {
        throw std::invalid_argument(
            "mpmc::flow: porosity must be finite in (0, 1]");
    }
}

inline void validate_snapshot_identity(
    const PoreVolumeComponentAccumulationSnapshot3P&
        snapshot) {
    if (snapshot.component_ids.size() < 2U ||
        snapshot
                .component_accumulation_mol_per_bulk_m3
                .size() !=
            snapshot.component_ids.size() ||
        !std::isfinite(
            snapshot
                .total_accumulation_mol_per_bulk_m3) ||
        !(snapshot
              .total_accumulation_mol_per_bulk_m3 >
          0.0)) {
        throw std::invalid_argument(
            "mpmc::flow: malformed pore-volume component accumulation snapshot");
    }
    validate_porosity(snapshot.porosity);
    double sum = 0.0;
    for (std::size_t component = 0U;
         component <
         snapshot.component_ids.size();
         ++component) {
        if (snapshot.component_ids[component].empty() ||
            !std::isfinite(
                snapshot
                    .component_accumulation_mol_per_bulk_m3
                    [component]) ||
            !(snapshot
                  .component_accumulation_mol_per_bulk_m3
                  [component] >
              0.0)) {
            throw std::invalid_argument(
                "mpmc::flow: invalid component accumulation identity/value");
        }
        for (std::size_t previous = 0U;
             previous < component;
             ++previous) {
            if (snapshot.component_ids[previous] ==
                snapshot.component_ids[component]) {
                throw std::invalid_argument(
                    "mpmc::flow: component accumulation IDs must be unique and ordered");
            }
        }
        sum +=
            snapshot
                .component_accumulation_mol_per_bulk_m3
                [component];
    }
    if (!near_roundoff(
            sum,
            snapshot
                .total_accumulation_mol_per_bulk_m3)) {
        throw std::invalid_argument(
            "mpmc::flow: component accumulation does not close to total accumulation");
    }
}

[[nodiscard]] inline double
d_saturation(
    const NaturalVariableLayout3P& layout,
    std::size_t phase,
    std::size_t column) {
    const auto s0 =
        layout.independent_saturation_unknown_index(
            PhaseSlot3::phase0);
    const auto s1 =
        layout.independent_saturation_unknown_index(
            PhaseSlot3::phase1);
    if (!s0 || !s1) {
        throw std::logic_error(
            "mpmc::flow: fixed-three-phase saturation chart is malformed");
    }
    if (column == *s0) {
        if (phase == 0U) {
            return 1.0;
        }
        if (phase == 2U) {
            return -1.0;
        }
    }
    if (column == *s1) {
        if (phase == 1U) {
            return 1.0;
        }
        if (phase == 2U) {
            return -1.0;
        }
    }
    return 0.0;
}

[[nodiscard]] inline double
d_composition(
    const NaturalVariableLayout3P& layout,
    std::size_t phase,
    std::size_t component,
    std::size_t column) {
    const auto identity =
        layout.composition_unknown_identity(
            column);
    if (!identity ||
        natural_variable_detail::
                checked_phase_index(
                    identity->phase) !=
            phase) {
        return 0.0;
    }
    if (component == identity->component) {
        return 1.0;
    }
    if (component ==
        layout.dependent_composition_component(
            static_cast<PhaseSlot3>(phase))) {
        return -1.0;
    }
    return 0.0;
}

inline void validate_density_linearization(
    const NaturalVariableCellState3P& state,
    const PhaseMolarDensityNaturalVariableLinearization3P&
        density) {
    const auto& state_layout =
        state.layout();
    if (density.layout.component_count() !=
            state_layout.component_count() ||
        density.layout.composition_pivot()
                .dependent_components() !=
            state_layout.composition_pivot()
                .dependent_components() ||
        density.layout.unknown_count() !=
            state_layout.unknown_count()) {
        throw std::invalid_argument(
            "mpmc::flow: molar-density linearization chart does not match current natural-variable state");
    }

    for (std::size_t phase = 0U;
         phase < fixed_three_phase_count;
         ++phase) {
        const double expected =
            state.phase_properties(
                     static_cast<PhaseSlot3>(
                         phase))
                .molar_density_mol_per_m3;
        const double actual =
            density
                .molar_density_mol_per_m3
                [phase];
        if (!std::isfinite(actual) ||
            !(actual > 0.0) ||
            !near_roundoff(
                actual,
                expected)) {
            throw std::invalid_argument(
                "mpmc::flow: molar-density linearization primal does not match current phase property");
        }
        if (density.gradient[phase].size() !=
            state_layout.unknown_count()) {
            throw std::invalid_argument(
                "mpmc::flow: molar-density gradient shape does not match natural-variable chart");
        }
        for (double derivative :
             density.gradient[phase]) {
            if (!std::isfinite(derivative)) {
                throw std::invalid_argument(
                    "mpmc::flow: molar-density gradient contains non-finite derivative");
            }
        }
    }
}

} // namespace component_accumulation_detail

[[nodiscard]] inline
PoreVolumeComponentAccumulationSnapshot3P
build_pore_volume_component_accumulation(
    const NaturalVariableCellState3P& state,
    double porosity) {
    component_accumulation_detail::
        validate_porosity(porosity);

    const std::size_t n =
        state.layout().component_count();
    if (state.component_ids().size() != n) {
        throw std::invalid_argument(
            "mpmc::flow: component identity/state size mismatch");
    }

    PoreVolumeComponentAccumulationSnapshot3P
        result;
    result.porosity = porosity;
    result.component_ids.assign(
        state.component_ids().begin(),
        state.component_ids().end());
    result
        .component_accumulation_mol_per_bulk_m3
        .assign(n, 0.0);

    double phase_total = 0.0;
    for (std::size_t phase = 0U;
         phase < fixed_three_phase_count;
         ++phase) {
        const auto slot =
            static_cast<PhaseSlot3>(phase);
        const double saturation =
            state.phase_saturation(slot);
        const double density =
            state.phase_properties(slot)
                .molar_density_mol_per_m3;
        const auto composition =
            state.phase_composition(slot);

        if (composition.size() != n ||
            !std::isfinite(saturation) ||
            !(saturation > 0.0) ||
            !std::isfinite(density) ||
            !(density > 0.0)) {
            throw std::invalid_argument(
                "mpmc::flow: invalid phase state for pore-volume accumulation");
        }

        phase_total +=
            saturation * density;
        for (std::size_t component = 0U;
             component < n;
             ++component) {
            result
                .component_accumulation_mol_per_bulk_m3
                [component] +=
                porosity *
                saturation *
                density *
                composition[component];
        }
    }

    result.total_accumulation_mol_per_bulk_m3 =
        porosity * phase_total;
    component_accumulation_detail::
        validate_snapshot_identity(result);
    return result;
}

[[nodiscard]] inline
PoreVolumeComponentAccumulationPair3P
make_pore_volume_component_accumulation_pair(
    PoreVolumeComponentAccumulationSnapshot3P
        current,
    PoreVolumeComponentAccumulationSnapshot3P
        previous) {
    component_accumulation_detail::
        validate_snapshot_identity(current);
    component_accumulation_detail::
        validate_snapshot_identity(previous);
    if (current.component_ids !=
        previous.component_ids) {
        throw std::invalid_argument(
            "mpmc::flow: current/previous accumulation component identity/order mismatch");
    }
    if (!component_accumulation_detail::
            near_roundoff(
                current.porosity,
                previous.porosity)) {
        throw std::invalid_argument(
            "mpmc::flow: rigid-medium current/previous porosity mismatch");
    }
    return {
        std::move(current),
        std::move(previous)};
}

[[nodiscard]] inline
PoreVolumeComponentAccumulationLinearization3P
build_pore_volume_component_accumulation_linearization(
    const NaturalVariableCellState3P& current,
    double porosity,
    const PhaseMolarDensityNaturalVariableLinearization3P&
        density) {
    component_accumulation_detail::
        validate_porosity(porosity);
    component_accumulation_detail::
        validate_density_linearization(
            current,
            density);

    const auto& layout =
        current.layout();
    const std::size_t n =
        layout.component_count();
    const std::size_t q_count =
        layout.unknown_count();
    if (n >
        std::numeric_limits<std::size_t>::
                max() /
            q_count) {
        throw std::length_error(
            "mpmc::flow: component accumulation Jacobian size overflow");
    }

    PoreVolumeComponentAccumulationLinearization3P
        result{
            layout,
            porosity,
            std::vector<std::string>{
                current.component_ids().begin(),
                current.component_ids().end()},
            q_count,
            std::vector<double>(
                n * q_count,
                0.0),
            std::vector<double>(
                q_count,
                0.0)};

    for (std::size_t column = 0U;
         column < q_count;
         ++column) {
        double total_derivative = 0.0;
        for (std::size_t phase = 0U;
             phase < fixed_three_phase_count;
             ++phase) {
            const auto slot =
                static_cast<PhaseSlot3>(phase);
            const double saturation =
                current.phase_saturation(slot);
            const double molar_density =
                density
                    .molar_density_mol_per_m3
                    [phase];
            const double dc =
                density.gradient[phase][column];
            const double ds =
                component_accumulation_detail::
                    d_saturation(
                        layout,
                        phase,
                        column);
            total_derivative +=
                porosity *
                (ds * molar_density +
                 saturation * dc);
        }
        result.total_accumulation_gradient[
            column] =
            total_derivative;

        double component_derivative_sum = 0.0;
        for (std::size_t component = 0U;
             component < n;
             ++component) {
            double derivative = 0.0;
            for (std::size_t phase = 0U;
                 phase < fixed_three_phase_count;
                 ++phase) {
                const auto slot =
                    static_cast<PhaseSlot3>(
                        phase);
                const double saturation =
                    current.phase_saturation(
                        slot);
                const double molar_density =
                    density
                        .molar_density_mol_per_m3
                        [phase];
                const double dc =
                    density.gradient[phase]
                                    [column];
                const double ds =
                    component_accumulation_detail::
                        d_saturation(
                            layout,
                            phase,
                            column);
                const double dx =
                    component_accumulation_detail::
                        d_composition(
                            layout,
                            phase,
                            component,
                            column);
                const double x =
                    current
                        .phase_composition(slot)
                        [component];

                derivative +=
                    porosity *
                    (ds *
                         molar_density *
                         x +
                     saturation *
                         dc *
                         x +
                     saturation *
                         molar_density *
                         dx);
            }
            if (!std::isfinite(derivative)) {
                throw std::range_error(
                    "mpmc::flow: non-finite component accumulation derivative");
            }
            result.component_jacobian[
                component * q_count +
                column] =
                derivative;
            component_derivative_sum +=
                derivative;
        }

        if (!component_accumulation_detail::
                near_roundoff(
                    component_derivative_sum,
                    total_derivative)) {
            throw std::runtime_error(
                "mpmc::flow: differentiated component accumulation does not close to total accumulation");
        }
    }

    return result;
}

} // namespace mpmc::flow

#endif // MPMC_FLOW_COMPONENT_ACCUMULATION_HPP
