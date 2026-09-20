#include <mpmc/flow/component_accumulation.hpp>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <iostream>
#include <limits>
#include <optional>
#include <source_location>
#include <span>
#include <stdexcept>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

bool component_accumulation_header();

namespace {

namespace fl = mpmc::flow;

void require(
    bool condition,
    std::string_view message,
    std::source_location where =
        std::source_location::current()) {
    if (!condition) {
        throw std::runtime_error(
            std::string{where.file_name()} +
            ":" +
            std::to_string(where.line()) +
            ": " +
            std::string{message});
    }
}

void near(
    double actual,
    double expected,
    double relative = 2.0e-11,
    double absolute = 2.0e-11,
    std::source_location where =
        std::source_location::current()) {
    if (!std::isfinite(actual) ||
        !std::isfinite(expected) ||
        std::abs(actual - expected) >
            absolute +
                relative *
                    std::max(
                        std::abs(actual),
                        std::abs(expected))) {
        std::cerr
            << "actual=" << actual
            << " expected=" << expected
            << '\n';
        require(false, "numeric mismatch", where);
    }
}

template <class Function>
void expect_invalid(
    Function&& function,
    std::string_view fragment) {
    try {
        function();
    } catch (const std::invalid_argument& error) {
        require(
            std::string_view{error.what()}.find(fragment) !=
                std::string_view::npos,
            "invalid_argument diagnostic changed");
        return;
    }
    throw std::runtime_error(
        "expected std::invalid_argument");
}

fl::PhasePropertyPrerequisiteInput
properties(double molar_density) {
    return {
        molar_density,
        700.0,
        1.0e-3,
        2.5e5,
        2.0e5};
}

std::vector<double> independent_values(
    std::span<const double> full,
    std::size_t dependent) {
    std::vector<double> result;
    result.reserve(full.size() - 1U);
    for (std::size_t component = 0U;
         component < full.size();
         ++component) {
        if (component != dependent) {
            result.push_back(full[component]);
        }
    }
    return result;
}

fl::NaturalVariableCellState3P make_state(
    std::array<double, 2> independent_saturations,
    const std::array<std::vector<double>, 3>& compositions,
    const std::array<double, 3>& densities) {
    const auto pivot =
        fl::NaturalVariableCompositionPivot3P::select(
            compositions);

    fl::NaturalVariableCellStateInput3P input;
    input.component_ids = {"A", "B", "C"};
    input.reference_pressure_pa = 10.0e6;
    input.temperature_k = 350.0;
    input.independent_saturations =
        independent_saturations;
    input.composition_pivot = pivot;
    for (std::size_t phase = 0U;
         phase < 3U;
         ++phase) {
        input.independent_phase_compositions[phase] =
            independent_values(
                compositions[phase],
                pivot.dependent_component(
                    static_cast<fl::PhaseSlot3>(phase)));
        input.phase_properties[phase] =
            properties(densities[phase]);
    }
    return fl::NaturalVariableCellState3P::create(
        std::move(input));
}

const std::array<std::vector<double>, 3>&
current_compositions() {
    static const std::array<std::vector<double>, 3>
        values{
            std::vector<double>{0.10, 0.70, 0.20},
            std::vector<double>{0.60, 0.20, 0.20},
            std::vector<double>{0.20, 0.30, 0.50}};
    return values;
}

void snapshot_pair() {
    constexpr double porosity = 0.25;
    const std::array<double, 3>
        current_density{5000.0, 3000.0, 7000.0};
    const auto current = make_state(
        {0.20, 0.30},
        current_compositions(),
        current_density);

    const auto current_snapshot =
        fl::build_pore_volume_component_accumulation(
            current,
            porosity);

    const std::array<double, 3>
        saturation{0.20, 0.30, 0.50};
    std::array<double, 3> expected{};
    double expected_total = 0.0;
    for (std::size_t phase = 0U;
         phase < 3U;
         ++phase) {
        expected_total +=
            porosity *
            saturation[phase] *
            current_density[phase];
        for (std::size_t component = 0U;
             component < 3U;
             ++component) {
            expected[component] +=
                porosity *
                saturation[phase] *
                current_density[phase] *
                current_compositions()[phase][component];
        }
    }

    require(
        current_snapshot.component_ids ==
            std::vector<std::string>{"A", "B", "C"},
        "component identity/order changed");
    for (std::size_t component = 0U;
         component < 3U;
         ++component) {
        near(
            current_snapshot.component(component),
            expected[component]);
    }
    near(
        current_snapshot.total_accumulation_mol_per_bulk_m3,
        expected_total);
    near(
        current_snapshot.component(0U) +
            current_snapshot.component(1U) +
            current_snapshot.component(2U),
        current_snapshot.total_accumulation_mol_per_bulk_m3);

    const std::array<std::vector<double>, 3>
        previous_compositions{
            std::vector<double>{0.20, 0.60, 0.20},
            std::vector<double>{0.50, 0.30, 0.20},
            std::vector<double>{0.25, 0.25, 0.50}};
    const auto previous = make_state(
        {0.25, 0.25},
        previous_compositions,
        {4800.0, 3200.0, 6800.0});
    const auto previous_snapshot =
        fl::build_pore_volume_component_accumulation(
            previous,
            porosity);

    const auto pair =
        fl::make_pore_volume_component_accumulation_pair(
            current_snapshot,
            previous_snapshot);
    require(
        pair.current.component_ids ==
                pair.previous.component_ids &&
            pair.current.porosity ==
                pair.previous.porosity,
        "current/previous accumulation pair lost identity");

    auto wrong_ids = previous_snapshot;
    wrong_ids.component_ids[1] = "X";
    expect_invalid(
        [&] {
            (void)fl::
                make_pore_volume_component_accumulation_pair(
                    current_snapshot,
                    wrong_ids);
        },
        "component identity/order mismatch");

    auto wrong_porosity = previous_snapshot;
    wrong_porosity.porosity = 0.26;
    expect_invalid(
        [&] {
            (void)fl::
                make_pore_volume_component_accumulation_pair(
                    current_snapshot,
                    wrong_porosity);
        },
        "porosity mismatch");
}

struct SyntheticDensityModel {
    std::array<double, 3> base{
        1500.0, 1700.0, 1900.0};
    std::array<double, 3> pressure_coeff{
        1.0e-4, 1.2e-4, 0.8e-4};
    std::array<double, 3> temperature_coeff{
        0.7, 0.5, 0.9};
    std::array<std::array<double, 3>, 3>
        composition_coeff{{
            {120.0, 80.0, 30.0},
            {60.0, 100.0, 40.0},
            {90.0, 20.0, 110.0}}};

    double density(
        std::size_t phase,
        double pressure,
        double temperature,
        std::span<const double> composition) const {
        double result =
            base[phase] +
            pressure_coeff[phase] * pressure +
            temperature_coeff[phase] * temperature;
        for (std::size_t component = 0U;
             component < composition.size();
             ++component) {
            result +=
                composition_coeff[phase][component] *
                composition[component];
        }
        return result;
    }
};

std::vector<double> natural_inputs(
    const fl::NaturalVariableCellState3P& state) {
    const auto& layout = state.layout();
    std::vector<double> inputs(
        layout.unknown_count(),
        0.0);
    inputs[layout.pressure_unknown_index()] =
        state.reference_pressure_pa();
    inputs[layout.temperature_unknown_index()] =
        state.temperature_k();
    inputs[*layout.independent_saturation_unknown_index(
        fl::PhaseSlot3::phase0)] =
        state.phase_saturation(fl::PhaseSlot3::phase0);
    inputs[*layout.independent_saturation_unknown_index(
        fl::PhaseSlot3::phase1)] =
        state.phase_saturation(fl::PhaseSlot3::phase1);

    for (std::size_t phase = 0U;
         phase < 3U;
         ++phase) {
        const auto slot =
            static_cast<fl::PhaseSlot3>(phase);
        const auto x =
            state.phase_composition(slot);
        for (std::size_t component = 0U;
             component < x.size();
             ++component) {
            const auto column =
                layout.independent_composition_unknown_index(
                    slot,
                    component);
            if (column) {
                inputs[*column] = x[component];
            }
        }
    }
    return inputs;
}

std::array<std::vector<double>, 3>
reconstruct_compositions(
    const fl::NaturalVariableLayout3P& layout,
    std::span<const double> inputs) {
    const std::size_t n =
        layout.component_count();
    std::array<std::vector<double>, 3> result;
    for (std::size_t phase = 0U;
         phase < 3U;
         ++phase) {
        const auto slot =
            static_cast<fl::PhaseSlot3>(phase);
        result[phase].assign(n, 0.0);
        const std::size_t dependent =
            layout.dependent_composition_component(slot);
        double sum = 0.0;
        for (std::size_t component = 0U;
             component < n;
             ++component) {
            const auto column =
                layout.independent_composition_unknown_index(
                    slot,
                    component);
            if (!column) {
                continue;
            }
            result[phase][component] =
                inputs[*column];
            sum += inputs[*column];
        }
        result[phase][dependent] =
            1.0 - sum;
    }
    return result;
}

std::vector<double> evaluate_synthetic_accumulation(
    const fl::NaturalVariableLayout3P& layout,
    std::span<const double> inputs,
    double porosity,
    const SyntheticDensityModel& model) {
    const auto compositions =
        reconstruct_compositions(layout, inputs);
    const double s0 =
        inputs[*layout.independent_saturation_unknown_index(
            fl::PhaseSlot3::phase0)];
    const double s1 =
        inputs[*layout.independent_saturation_unknown_index(
            fl::PhaseSlot3::phase1)];
    const std::array<double, 3> saturation{
        s0, s1, 1.0 - s0 - s1};

    std::vector<double> result(
        layout.component_count(),
        0.0);
    for (std::size_t phase = 0U;
         phase < 3U;
         ++phase) {
        const double density =
            model.density(
                phase,
                inputs[layout.pressure_unknown_index()],
                inputs[layout.temperature_unknown_index()],
                compositions[phase]);
        for (std::size_t component = 0U;
             component < layout.component_count();
             ++component) {
            result[component] +=
                porosity *
                saturation[phase] *
                density *
                compositions[phase][component];
        }
    }
    return result;
}

void jacobian() {
    constexpr double porosity = 0.22;
    const SyntheticDensityModel model;
    const auto& compositions =
        current_compositions();

    const auto pivot =
        fl::NaturalVariableCompositionPivot3P::select(
            compositions);
    const fl::NaturalVariableLayout3P layout{
        pivot};

    std::array<double, 3> densities{};
    for (std::size_t phase = 0U;
         phase < 3U;
         ++phase) {
        densities[phase] =
            model.density(
                phase,
                10.0e6,
                350.0,
                compositions[phase]);
    }

    const auto state =
        make_state(
            {0.20, 0.30},
            compositions,
            densities);
    require(
        state.layout().composition_pivot().dependent_components() ==
            pivot.dependent_components(),
        "current state pivot changed");

    fl::PhaseMolarDensityNaturalVariableLinearization3P
        density_linearization{
            state.layout(),
            densities,
            {}};

    for (std::size_t phase = 0U;
         phase < 3U;
         ++phase) {
        auto& gradient =
            density_linearization.gradient[phase];
        gradient.assign(
            state.layout().unknown_count(),
            0.0);
        gradient[state.layout().pressure_unknown_index()] =
            model.pressure_coeff[phase];
        gradient[state.layout().temperature_unknown_index()] =
            model.temperature_coeff[phase];

        const auto slot =
            static_cast<fl::PhaseSlot3>(phase);
        const std::size_t dependent =
            state.layout().dependent_composition_component(
                slot);
        for (std::size_t component = 0U;
             component < 3U;
             ++component) {
            const auto column =
                state.layout().independent_composition_unknown_index(
                    slot,
                    component);
            if (column) {
                gradient[*column] =
                    model.composition_coeff[phase][component] -
                    model.composition_coeff[phase][dependent];
            }
        }
    }

    const auto linearization =
        fl::build_pore_volume_component_accumulation_linearization(
            state,
            porosity,
            density_linearization);

    require(
        linearization.component_ids ==
                std::vector<std::string>{"A", "B", "C"} &&
            linearization.input_count ==
                state.layout().unknown_count() &&
            linearization.layout.composition_pivot().dependent_components() ==
                std::vector<std::size_t>{
                    state.layout().composition_pivot().dependent_components().begin(),
                    state.layout().composition_pivot().dependent_components().end()},
        "accumulation Jacobian lost chart/component identity");

    const auto inputs =
        natural_inputs(state);

    for (std::size_t column = 0U;
         column < inputs.size();
         ++column) {
        double step = 1.0e-6;
        if (column ==
            state.layout().pressure_unknown_index()) {
            step = 10.0;
        } else if (
            column ==
            state.layout().temperature_unknown_index()) {
            step = 1.0e-4;
        }

        auto plus = inputs;
        auto minus = inputs;
        plus[column] += step;
        minus[column] -= step;

        const auto plus_value =
            evaluate_synthetic_accumulation(
                state.layout(),
                plus,
                porosity,
                model);
        const auto minus_value =
            evaluate_synthetic_accumulation(
                state.layout(),
                minus,
                porosity,
                model);

        double sum_jacobian = 0.0;
        for (std::size_t component = 0U;
             component < 3U;
             ++component) {
            const double finite_difference =
                (plus_value[component] -
                 minus_value[component]) /
                (2.0 * step);
            const double actual =
                linearization.d_component(
                    component,
                    column);
            const double allowed =
                2.0e-6 +
                3.0e-7 *
                    std::abs(finite_difference);
            require(
                std::isfinite(actual) &&
                    std::abs(
                        actual -
                        finite_difference) <=
                        allowed,
                "analytic accumulation Jacobian disagrees with fresh perturbation");
            sum_jacobian += actual;
        }
        near(
            sum_jacobian,
            linearization.d_total(column),
            5.0e-12,
            5.0e-9);
    }
}

void invalid_inputs() {
    const auto& compositions =
        current_compositions();
    const std::array<double, 3>
        densities{5000.0, 3000.0, 7000.0};
    const auto state =
        make_state(
            {0.20, 0.30},
            compositions,
            densities);

    expect_invalid(
        [&] {
            (void)fl::
                build_pore_volume_component_accumulation(
                    state,
                    0.0);
        },
        "porosity");

    auto wrong_layout =
        fl::NaturalVariableLayout3P{
            fl::NaturalVariableCompositionPivot3P::
                fixed_last(3U)};
    fl::PhaseMolarDensityNaturalVariableLinearization3P bad{
        wrong_layout,
        densities,
        {std::vector<double>(10U, 0.0),
         std::vector<double>(10U, 0.0),
         std::vector<double>(10U, 0.0)}};
    expect_invalid(
        [&] {
            (void)fl::
                build_pore_volume_component_accumulation_linearization(
                    state,
                    0.25,
                    bad);
        },
        "chart does not match");

    fl::PhaseMolarDensityNaturalVariableLinearization3P nonfinite{
        state.layout(),
        densities,
        {std::vector<double>(
             state.layout().unknown_count(), 0.0),
         std::vector<double>(
             state.layout().unknown_count(), 0.0),
         std::vector<double>(
             state.layout().unknown_count(), 0.0)}};
    nonfinite.gradient[1][0] =
        std::numeric_limits<double>::infinity();
    expect_invalid(
        [&] {
            (void)fl::
                build_pore_volume_component_accumulation_linearization(
                    state,
                    0.25,
                    nonfinite);
        },
        "non-finite");
}

void headers() {
    require(
        component_accumulation_header(),
        "component accumulation header probe failed");
}

using Test =
    std::pair<std::string_view, void (*)()>;

constexpr Test tests[]{
    {"snapshot_pair", snapshot_pair},
    {"jacobian", jacobian},
    {"invalid_inputs", invalid_inputs},
    {"headers", headers}};

} // namespace

int main(int argc, char** argv) {
    try {
        if (argc != 2) {
            throw std::invalid_argument(
                "one test name required");
        }
        for (const auto& [name, run] : tests) {
            if (name == argv[1]) {
                run();
                std::cout
                    << "[PASS] "
                    << name
                    << '\n';
                return 0;
            }
        }
        throw std::invalid_argument(
            "unknown test");
    } catch (const std::exception& error) {
        std::cerr
            << "[FAIL] "
            << error.what()
            << '\n';
        return 1;
    }
}
