#ifndef MPMC_TEST_FLOW_FUGACITY_ADAPTER_REGRESSION_HPP
#define MPMC_TEST_FLOW_FUGACITY_ADAPTER_REGRESSION_HPP

#include <mpmc/ad/runtime_differentiate.hpp>
#include <mpmc/flow/thermodynamics_fugacity_adapters.hpp>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <limits>
#include <span>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace mpmc::test::flow_adapter {

struct JacobianTolerance {
    double pressure_absolute{5.0e-10};
    double temperature_absolute{5.0e-5};
    double composition_absolute{2.0e-4};
    double relative{5.0e-3};
};

inline void require(
    bool condition,
    const char* message) {
    if (!condition) {
        throw std::runtime_error(message);
    }
}

template <typename Adapter>
[[nodiscard]] std::vector<double>
evaluate_values(
    const flow::NaturalVariableLayout3P& layout,
    std::span<const double> inputs,
    Adapter& adapter) {
    const std::size_t n =
        layout.component_count();
    if (inputs.size() !=
        layout.unknown_count()) {
        throw std::invalid_argument(
            "flow adapter regression: natural-variable input size mismatch");
    }

    std::array<std::vector<double>, 3>
        compositions;
    for (std::size_t phase = 0U;
         phase < 3U;
         ++phase) {
        const auto slot =
            static_cast<flow::PhaseSlot3>(
                phase);
        compositions[phase].reserve(n);
        double sum = 0.0;
        for (std::size_t component = 0U;
             component + 1U < n;
             ++component) {
            const auto index =
                layout
                    .independent_composition_unknown_index(
                        slot,
                        component);
            if (!index) {
                throw std::logic_error(
                    "flow adapter regression: missing independent composition index");
            }
            const double value =
                inputs[*index];
            compositions[phase].push_back(
                value);
            sum += value;
        }
        compositions[phase].push_back(
            1.0 - sum);
    }

    const flow::
        FugacityEquilibriumStateView3P<double>
        state{
            {
                inputs[
                    layout
                        .pressure_unknown_index()],
                inputs[
                    layout
                        .pressure_unknown_index()],
                inputs[
                    layout
                        .pressure_unknown_index()]},
            inputs[
                layout
                    .temperature_unknown_index()],
            {
                std::span<const double>{
                    compositions[0]},
                std::span<const double>{
                    compositions[1]},
                std::span<const double>{
                    compositions[2]}}};

    const auto residual =
        flow::
            evaluate_fugacity_equilibrium_residual_3p(
                state,
                adapter);
    return std::vector<double>{
        residual.values().begin(),
        residual.values().end()};
}

template <typename Adapter>
void verify_accepted_state_and_jacobian(
    double pressure_pa,
    double temperature_k,
    const std::array<double, 3>&
        phase_fractions,
    const std::array<std::vector<double>, 3>&
        phase_compositions,
    Adapter& adapter,
    double residual_absolute_tolerance,
    JacobianTolerance tolerance = {}) {
    const std::size_t n =
        phase_compositions[0].size();
    require(
        n >= 2U &&
            phase_compositions[1].size() ==
                n &&
            phase_compositions[2].size() ==
                n,
        "flow adapter regression: invalid phase composition shape");

    const flow::NaturalVariableLayout3P
        layout{n};
    std::vector<double> inputs(
        layout.unknown_count(),
        0.0);
    inputs[
        layout.pressure_unknown_index()] =
        pressure_pa;
    inputs[
        layout.temperature_unknown_index()] =
        temperature_k;

    const auto s0 =
        layout
            .independent_saturation_unknown_index(
                flow::PhaseSlot3::phase0);
    const auto s1 =
        layout
            .independent_saturation_unknown_index(
                flow::PhaseSlot3::phase1);
    require(
        s0.has_value() && s1.has_value(),
        "flow adapter regression: saturation layout changed");
    inputs[*s0] = phase_fractions[0];
    inputs[*s1] = phase_fractions[1];

    for (std::size_t phase = 0U;
         phase < 3U;
         ++phase) {
        const auto slot =
            static_cast<flow::PhaseSlot3>(
                phase);
        for (std::size_t component = 0U;
             component + 1U < n;
             ++component) {
            const auto index =
                layout
                    .independent_composition_unknown_index(
                        slot,
                        component);
            require(
                index.has_value(),
                "flow adapter regression: composition layout changed");
            inputs[*index] =
                phase_compositions[phase]
                                  [component];
        }
    }

    using D =
        ad::Dual<double, 4U>;
    ad::RuntimeJacobianWorkspace<
        double,
        4U>
        workspace;

    const auto callback =
        [&](std::span<const D> active,
            std::span<D> outputs) {
            std::array<
                std::vector<D>,
                3>
                compositions;
            for (std::size_t phase = 0U;
                 phase < 3U;
                 ++phase) {
                const auto slot =
                    static_cast<
                        flow::PhaseSlot3>(
                        phase);
                D sum{};
                compositions[phase]
                    .reserve(n);
                for (std::size_t component =
                         0U;
                     component + 1U < n;
                     ++component) {
                    const auto index =
                        layout
                            .independent_composition_unknown_index(
                                slot,
                                component);
                    if (!index) {
                        throw std::logic_error(
                            "flow adapter regression: AD composition index missing");
                    }
                    compositions[phase]
                        .push_back(
                            active[*index]);
                    sum += active[*index];
                }
                compositions[phase]
                    .push_back(
                        D{1.0} - sum);
            }

            const D& pressure =
                active[
                    layout
                        .pressure_unknown_index()];
            const flow::
                FugacityEquilibriumStateView3P<
                    D>
                state{
                    {
                        pressure,
                        pressure,
                        pressure},
                    active[
                        layout
                            .temperature_unknown_index()],
                    {
                        std::span<const D>{
                            compositions[0]},
                        std::span<const D>{
                            compositions[1]},
                        std::span<const D>{
                            compositions[2]}}};

            const auto residual =
                flow::
                    evaluate_fugacity_equilibrium_residual_3p(
                        state,
                        adapter);
            if (outputs.size() !=
                residual.values().size()) {
                throw std::logic_error(
                    "flow adapter regression: AD output shape mismatch");
            }
            for (std::size_t row = 0U;
                 row < outputs.size();
                 ++row) {
                outputs[row] =
                    residual.values()[row];
            }
        };

    const auto derived =
        workspace.evaluate(
            callback,
            inputs,
            2U * n,
            {
                layout.unknown_count(),
                2U * n,
                (2U * n) *
                    layout.unknown_count()});

    require(
        derived.input_count ==
                layout.unknown_count() &&
            derived.output_count ==
                2U * n,
        "flow adapter regression: Jacobian shape changed");

    double maximum_residual = 0.0;
    for (double value :
         derived.values) {
        maximum_residual =
            std::max(
                maximum_residual,
                std::abs(value));
    }
    require(
        maximum_residual <=
            residual_absolute_tolerance,
        "flow adapter regression: accepted three-phase state no longer satisfies fugacity equilibrium");

    for (std::size_t row = 0U;
         row < derived.output_count;
         ++row) {
        require(
            derived.jacobian[
                row *
                    derived.input_count +
                *s0] == 0.0 &&
                derived.jacobian[
                    row *
                        derived.input_count +
                    *s1] == 0.0,
            "flow adapter regression: pc=none fugacity rows acquired a saturation derivative");
    }

    const auto strongest_independent =
        [&](std::size_t phase) {
            const double dependent =
                phase_compositions[phase]
                                  .back();
            std::size_t best = 0U;
            double best_margin = 0.0;
            for (std::size_t component =
                     0U;
                 component + 1U < n;
                 ++component) {
                const double margin =
                    std::min(
                        phase_compositions[
                            phase][component],
                        dependent);
                if (margin >
                    best_margin) {
                    best = component;
                    best_margin = margin;
                }
            }
            require(
                best_margin > 0.0,
                "flow adapter regression: no positive composition perturbation direction");
            return std::pair{
                best,
                best_margin};
        };

    const auto [phase1_component,
                phase1_margin] =
        strongest_independent(1U);
    const auto [phase2_component,
                phase2_margin] =
        strongest_independent(2U);

    const auto phase1_column =
        *layout
             .independent_composition_unknown_index(
                 flow::PhaseSlot3::phase1,
                 phase1_component);
    const auto phase2_column =
        *layout
             .independent_composition_unknown_index(
                 flow::PhaseSlot3::phase2,
                 phase2_component);

    const std::array<std::size_t, 4>
        checked_columns{
            layout.pressure_unknown_index(),
            layout.temperature_unknown_index(),
            phase1_column,
            phase2_column};

    const std::array<double, 4> steps{
        std::max(
            10.0,
            1.0e-6 *
                std::abs(pressure_pa)),
        std::max(
            1.0e-4,
            1.0e-6 *
                std::abs(temperature_k)),
        std::min(
            1.0e-6,
            0.05 * phase1_margin),
        std::min(
            1.0e-6,
            0.05 * phase2_margin)};

    const std::array<double, 4>
        absolute_tolerances{
            tolerance.pressure_absolute,
            tolerance.temperature_absolute,
            tolerance.composition_absolute,
            tolerance.composition_absolute};

    for (std::size_t check = 0U;
         check < checked_columns.size();
         ++check) {
        const std::size_t column =
            checked_columns[check];
        const double step =
            steps[check];
        require(
            std::isfinite(step) &&
                step > 0.0,
            "flow adapter regression: invalid finite-difference step");

        auto plus = inputs;
        auto minus = inputs;
        plus[column] += step;
        minus[column] -= step;

        const auto plus_values =
            evaluate_values(
                layout,
                plus,
                adapter);
        const auto minus_values =
            evaluate_values(
                layout,
                minus,
                adapter);

        for (std::size_t row = 0U;
             row < derived.output_count;
             ++row) {
            const double finite_difference =
                (plus_values[row] -
                 minus_values[row]) /
                (2.0 * step);
            const double automatic =
                derived.jacobian[
                    row *
                        derived.input_count +
                    column];
            const double allowed =
                absolute_tolerances[check] +
                tolerance.relative *
                    std::abs(
                        finite_difference);
            require(
                std::isfinite(automatic) &&
                    std::isfinite(
                        finite_difference) &&
                    std::abs(
                        automatic -
                        finite_difference) <=
                        allowed,
                "flow adapter regression: AD Jacobian disagrees with fixed-branch central perturbation");
        }
    }
}

} // namespace mpmc::test::flow_adapter

#endif // MPMC_TEST_FLOW_FUGACITY_ADAPTER_REGRESSION_HPP
