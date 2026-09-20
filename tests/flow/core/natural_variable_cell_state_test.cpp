#include <mpmc/flow/natural_variable_cell_state.hpp>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <iostream>
#include <limits>
#include <optional>
#include <source_location>
#include <stdexcept>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

bool natural_variable_cell_state_header();

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
    double absolute = 8.0 *
        std::numeric_limits<double>::epsilon(),
    std::source_location where =
        std::source_location::current()) {
    if (!std::isfinite(actual) ||
        !std::isfinite(expected) ||
        std::abs(actual - expected) >
            absolute *
                std::max(
                    1.0,
                    std::max(
                        std::abs(actual),
                        std::abs(expected)))) {
        throw std::runtime_error(
            std::string{where.file_name()} +
            ":" +
            std::to_string(where.line()) +
            ": numeric mismatch");
    }
}

template <class Function>
void expect_invalid(
    Function&& function,
    std::string_view expected_fragment) {
    try {
        function();
    } catch (const std::invalid_argument& error) {
        require(
            std::string_view{error.what()}.find(
                expected_fragment) !=
                std::string_view::npos,
            "invalid_argument diagnostic changed");
        return;
    }
    throw std::runtime_error(
        "expected std::invalid_argument");
}

fl::PhasePropertyPrerequisiteInput
properties(double scale) {
    return fl::PhasePropertyPrerequisiteInput{
        5000.0 * scale,
        700.0 * scale,
        1.0e-3 * scale,
        2.5e5 * scale,
        2.0e5 * scale};
}

fl::NaturalVariableCellStateInput3P
valid_input() {
    fl::NaturalVariableCellStateInput3P input;
    input.component_ids = {"A", "B", "C"};
    input.reference_pressure_pa = 12.0e6;
    input.temperature_k = 360.0;
    input.independent_saturations = {0.20, 0.30};
    input.independent_phase_compositions = {
        std::vector<double>{0.20, 0.30},
        std::vector<double>{0.10, 0.60},
        std::vector<double>{0.40, 0.20}};
    input.phase_properties = {
        properties(1.0),
        properties(1.2),
        properties(1.4)};
    return input;
}

void layout() {
    const fl::NaturalVariableLayout3P layout{
        4U};

    require(
        layout.component_count() == 4U &&
            layout.phase_count() == 3U &&
            layout.unknown_count() == 13U &&
            layout.equation_count() == 13U,
        "natural-variable layout size changed");
    require(
        layout.pressure_unknown_index() == 0U &&
            layout.temperature_unknown_index() == 1U,
        "pressure/temperature unknown indices changed");
    require(
        layout
                .independent_saturation_unknown_index(
                    fl::PhaseSlot3::phase0) ==
            std::optional<std::size_t>{2U} &&
            layout
                .independent_saturation_unknown_index(
                    fl::PhaseSlot3::phase1) ==
            std::optional<std::size_t>{3U} &&
            !layout
                 .independent_saturation_unknown_index(
                     fl::PhaseSlot3::phase2)
                 .has_value(),
        "saturation independent/dependent layout changed");

    require(
        layout
                .independent_composition_unknown_index(
                    fl::PhaseSlot3::phase0,
                    0U) ==
            std::optional<std::size_t>{4U} &&
            layout
                .independent_composition_unknown_index(
                    fl::PhaseSlot3::phase0,
                    2U) ==
            std::optional<std::size_t>{6U} &&
            !layout
                 .independent_composition_unknown_index(
                     fl::PhaseSlot3::phase0,
                     3U)
                 .has_value() &&
            layout
                .independent_composition_unknown_index(
                    fl::PhaseSlot3::phase2,
                    2U) ==
            std::optional<std::size_t>{12U},
        "composition unknown block layout changed");

    require(
        layout.component_conservation_equation_index(
            3U) == 3U &&
            layout.energy_equation_index() == 4U &&
            layout.fugacity_equilibrium_equation_index(
                fl::PhaseSlot3::phase1,
                0U) == 5U &&
            layout.fugacity_equilibrium_equation_index(
                fl::PhaseSlot3::phase2,
                3U) == 12U,
        "equation row layout changed");

    expect_invalid(
        [] {
            (void)fl::NaturalVariableLayout3P{1U};
        },
        "at least two components");
    expect_invalid(
        [&] {
            (void)layout
                .fugacity_equilibrium_equation_index(
                    fl::PhaseSlot3::phase0,
                    0U);
        },
        "reference phase");
}


void composition_pivot() {
    const std::array<std::vector<double>, 3>
        compositions{
            std::vector<double>{0.10, 0.70, 0.20},
            std::vector<double>{0.60, 0.20, 0.20},
            std::vector<double>{0.20, 0.30, 0.50}};
    const auto pivot =
        fl::NaturalVariableCompositionPivot3P::
            select(compositions);

    require(
        pivot.component_count() == 3U &&
            pivot.dependent_component(
                fl::PhaseSlot3::phase0) == 1U &&
            pivot.dependent_component(
                fl::PhaseSlot3::phase1) == 0U &&
            pivot.dependent_component(
                fl::PhaseSlot3::phase2) == 2U,
        "largest-composition pivot selection changed");

    const fl::NaturalVariableLayout3P layout{
        pivot};
    require(
        layout
                .independent_composition_unknown_index(
                    fl::PhaseSlot3::phase0,
                    0U) ==
            std::optional<std::size_t>{4U} &&
            !layout
                 .independent_composition_unknown_index(
                     fl::PhaseSlot3::phase0,
                     1U)
                 .has_value() &&
            layout
                .independent_composition_unknown_index(
                    fl::PhaseSlot3::phase0,
                    2U) ==
            std::optional<std::size_t>{5U} &&
            !layout
                 .independent_composition_unknown_index(
                     fl::PhaseSlot3::phase1,
                     0U)
                 .has_value() &&
            layout
                .independent_composition_unknown_index(
                    fl::PhaseSlot3::phase1,
                    1U) ==
            std::optional<std::size_t>{6U} &&
            layout
                .independent_composition_unknown_index(
                    fl::PhaseSlot3::phase1,
                    2U) ==
            std::optional<std::size_t>{7U} &&
            layout
                .independent_composition_unknown_index(
                    fl::PhaseSlot3::phase2,
                    0U) ==
            std::optional<std::size_t>{8U} &&
            layout
                .independent_composition_unknown_index(
                    fl::PhaseSlot3::phase2,
                    1U) ==
            std::optional<std::size_t>{9U} &&
            !layout
                 .independent_composition_unknown_index(
                     fl::PhaseSlot3::phase2,
                     2U)
                 .has_value(),
        "pivoted composition-column mapping changed");

    const auto column5 =
        layout.composition_unknown_identity(5U);
    const auto column6 =
        layout.composition_unknown_identity(6U);
    const auto column9 =
        layout.composition_unknown_identity(9U);
    require(
        column5 &&
            column5->phase ==
                fl::PhaseSlot3::phase0 &&
            column5->component == 2U &&
            column6 &&
            column6->phase ==
                fl::PhaseSlot3::phase1 &&
            column6->component == 1U &&
            column9 &&
            column9->phase ==
                fl::PhaseSlot3::phase2 &&
            column9->component == 1U,
        "composition unknown identity reverse mapping changed");

    const auto phase1_row =
        layout.fugacity_equilibrium_row_identity(
            layout
                .fugacity_equilibrium_equation_index(
                    fl::PhaseSlot3::phase1,
                    2U));
    const auto phase2_row =
        layout.fugacity_equilibrium_row_identity(
            layout
                .fugacity_equilibrium_equation_index(
                    fl::PhaseSlot3::phase2,
                    0U));
    require(
        phase1_row &&
            phase1_row->non_reference_phase ==
                fl::PhaseSlot3::phase1 &&
            phase1_row->component == 2U &&
            phase2_row &&
            phase2_row->non_reference_phase ==
                fl::PhaseSlot3::phase2 &&
            phase2_row->component == 0U,
        "fugacity row identity changed under composition pivot");

    const fl::NaturalVariableLayout3P default_layout{
        3U};
    require(
        layout
                .fugacity_equilibrium_equation_index(
                    fl::PhaseSlot3::phase1,
                    2U) ==
            default_layout
                .fugacity_equilibrium_equation_index(
                    fl::PhaseSlot3::phase1,
                    2U) &&
            layout
                .fugacity_equilibrium_equation_index(
                    fl::PhaseSlot3::phase2,
                    0U) ==
            default_layout
                .fugacity_equilibrium_equation_index(
                    fl::PhaseSlot3::phase2,
                    0U),
        "residual row order depends on pivot");

    const std::array<std::vector<double>, 3>
        tied{
            std::vector<double>{0.45, 0.45, 0.10},
            std::vector<double>{0.20, 0.40, 0.40},
            std::vector<double>{0.50, 0.25, 0.25}};
    const auto tied_pivot =
        fl::NaturalVariableCompositionPivot3P::
            select(tied);
    require(
        tied_pivot.dependent_component(
            fl::PhaseSlot3::phase0) == 0U &&
            tied_pivot.dependent_component(
                fl::PhaseSlot3::phase1) == 1U,
        "pivot tie-break is not lowest component index");

    const std::array<std::vector<double>, 3>
        trace{
            std::vector<double>{
                0.9991426078018199,
                0.0008573921907259,
                7.4541155e-12},
            std::vector<double>{0.30, 0.20, 0.50},
            std::vector<double>{0.10, 0.80, 0.10}};
    const auto trace_pivot =
        fl::NaturalVariableCompositionPivot3P::
            select(trace);
    require(
        trace_pivot.dependent_component(
            fl::PhaseSlot3::phase0) == 0U,
        "trace final component was selected as dependent");

    auto input = valid_input();
    input.composition_pivot = pivot;
    input.independent_phase_compositions = {
        std::vector<double>{0.10, 0.20},
        std::vector<double>{0.20, 0.20},
        std::vector<double>{0.20, 0.30}};
    const auto state =
        fl::NaturalVariableCellState3P::create(
            std::move(input));
    require(
        state.component_ids()[0] == "A" &&
            state.component_ids()[1] == "B" &&
            state.component_ids()[2] == "C" &&
            state.layout()
                    .dependent_composition_component(
                        fl::PhaseSlot3::phase0) ==
                1U,
        "pivot changed ordered component identity");
    for (std::size_t phase = 0U;
         phase < 3U;
         ++phase) {
        const auto actual =
            state.phase_composition(
                static_cast<fl::PhaseSlot3>(
                    phase));
        for (std::size_t component = 0U;
             component < 3U;
             ++component) {
            near(
                actual[component],
                compositions[phase]
                            [component]);
        }
    }

    expect_invalid(
        [] {
            const std::array<std::vector<double>, 3>
                invalid{
                    std::vector<double>{0.5, 0.5, 0.0},
                    std::vector<double>{0.4, 0.3, 0.3},
                    std::vector<double>{0.4, 0.3, 0.3}};
            (void)fl::
                NaturalVariableCompositionPivot3P::
                    select(invalid);
        },
        "positive-support");
}

void valid_state() {
    const auto state =
        fl::NaturalVariableCellState3P::create(
            valid_input());

    require(
        state.layout().component_count() == 3U &&
            state.layout().unknown_count() == 10U &&
            state.layout().equation_count() == 10U,
        "three-component state layout changed");
    require(
        state.component_ids().size() == 3U &&
            state.component_ids()[0] == "A" &&
            state.component_ids()[1] == "B" &&
            state.component_ids()[2] == "C",
        "ordered component identity changed");
    near(
        state.reference_pressure_pa(),
        12.0e6);
    near(
        state.temperature_k(),
        360.0);

    const auto pressures =
        state.phase_pressures_pa();
    for (const double pressure : pressures) {
        near(pressure, 12.0e6);
    }
    near(
        state.phase_pressure_pa(
            fl::PhaseSlot3::phase1),
        12.0e6);

    near(
        state.phase_saturation(
            fl::PhaseSlot3::phase0),
        0.20);
    near(
        state.phase_saturation(
            fl::PhaseSlot3::phase1),
        0.30);
    near(
        state.phase_saturation(
            fl::PhaseSlot3::phase2),
        0.50);

    const auto phase0 =
        state.phase_composition(
            fl::PhaseSlot3::phase0);
    const auto phase1 =
        state.phase_composition(
            fl::PhaseSlot3::phase1);
    const auto phase2 =
        state.phase_composition(
            fl::PhaseSlot3::phase2);
    require(
        phase0.size() == 3U &&
            phase1.size() == 3U &&
            phase2.size() == 3U,
        "phase composition size changed");
    near(phase0[0], 0.20);
    near(phase0[1], 0.30);
    near(phase0[2], 0.50);
    near(phase1[2], 0.30);
    near(phase2[2], 0.40);

    const auto& p1 =
        state.phase_properties(
            fl::PhaseSlot3::phase1);
    near(
        p1.molar_density_mol_per_m3,
        6000.0);
    near(
        p1.mass_density_kg_per_m3,
        840.0);
    near(
        p1.dynamic_viscosity_pa_s,
        1.2e-3);
    near(
        p1.specific_enthalpy_j_per_kg,
        3.0e5);
    near(
        p1.specific_internal_energy_j_per_kg,
        2.4e5);
}

void invalid_state() {
    {
        auto input = valid_input();
        input.component_ids[2] = "A";
        expect_invalid(
            [&] {
                (void)fl::
                    NaturalVariableCellState3P::
                        create(std::move(input));
            },
            "component ids must be unique");
    }
    {
        auto input = valid_input();
        input.reference_pressure_pa = 0.0;
        expect_invalid(
            [&] {
                (void)fl::
                    NaturalVariableCellState3P::
                        create(std::move(input));
            },
            "reference pressure");
    }
    {
        auto input = valid_input();
        input.temperature_k =
            std::numeric_limits<double>::
                quiet_NaN();
        expect_invalid(
            [&] {
                (void)fl::
                    NaturalVariableCellState3P::
                        create(std::move(input));
            },
            "temperature");
    }
    {
        auto input = valid_input();
        input.independent_saturations = {
            0.0,
            0.30};
        expect_invalid(
            [&] {
                (void)fl::
                    NaturalVariableCellState3P::
                        create(std::move(input));
            },
            "positive-support independent saturations");
    }
    {
        auto input = valid_input();
        input.independent_saturations = {
            0.60,
            0.40};
        expect_invalid(
            [&] {
                (void)fl::
                    NaturalVariableCellState3P::
                        create(std::move(input));
            },
            "reconstructed phase2 saturation");
    }
    {
        auto input = valid_input();
        input
            .independent_phase_compositions[0]
            [0] = 0.0;
        expect_invalid(
            [&] {
                (void)fl::
                    NaturalVariableCellState3P::
                        create(std::move(input));
            },
            "positive-support composition");
    }
    {
        auto input = valid_input();
        input
            .independent_phase_compositions[1] = {
            0.50,
            0.50};
        expect_invalid(
            [&] {
                (void)fl::
                    NaturalVariableCellState3P::
                        create(std::move(input));
            },
            "reconstructed dependent composition");
    }
    {
        auto input = valid_input();
        input
            .independent_phase_compositions[2]
            .pop_back();
        expect_invalid(
            [&] {
                (void)fl::
                    NaturalVariableCellState3P::
                        create(std::move(input));
            },
            "composition size mismatch");
    }
    {
        auto input = valid_input();
        input.phase_properties[1]
            .dynamic_viscosity_pa_s.reset();
        expect_invalid(
            [&] {
                (void)fl::
                    NaturalVariableCellState3P::
                        create(std::move(input));
            },
            "unsupported required phase property 'dynamic_viscosity_pa_s'");
    }
    {
        auto input = valid_input();
        input.phase_properties[2]
            .mass_density_kg_per_m3 = 0.0;
        expect_invalid(
            [&] {
                (void)fl::
                    NaturalVariableCellState3P::
                        create(std::move(input));
            },
            "mass_density_kg_per_m3");
    }
    {
        auto input = valid_input();
        input.phase_properties[0]
            .specific_enthalpy_j_per_kg =
            std::numeric_limits<double>::
                infinity();
        expect_invalid(
            [&] {
                (void)fl::
                    NaturalVariableCellState3P::
                        create(std::move(input));
            },
            "specific_enthalpy_j_per_kg");
    }
    {
        const fl::NaturalVariableLayout3P layout{
            3U};
        expect_invalid(
            [&] {
                (void)layout
                    .independent_saturation_unknown_index(
                        static_cast<fl::PhaseSlot3>(
                            99U));
            },
            "invalid fixed-three-phase slot");
    }
}

void headers() {
    require(
        natural_variable_cell_state_header(),
        "flow public-header probe failed");
}

using Test =
    std::pair<std::string_view, void (*)()>;

constexpr Test tests[]{
    {"layout", layout},
    {"composition_pivot", composition_pivot},
    {"valid_state", valid_state},
    {"invalid_state", invalid_state},
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
