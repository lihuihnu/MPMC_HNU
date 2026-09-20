#include <mpmc/flow/energy_accumulation.hpp>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <iostream>
#include <limits>
#include <source_location>
#include <span>
#include <stdexcept>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

bool energy_accumulation_header();

namespace {

namespace flow = mpmc::flow;

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
    double relative = 2.0e-10,
    double absolute = 2.0e-8,
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

std::vector<double> independent_values(
    std::span<const double> full,
    std::size_t dependent) {
    std::vector<double> result;
    for (std::size_t component = 0U;
         component < full.size();
         ++component) {
        if (component != dependent) {
            result.push_back(
                full[component]);
        }
    }
    return result;
}

flow::NaturalVariableCellState3P
make_state() {
    const std::array<
        std::vector<double>,
        3>
        compositions{
            std::vector<double>{
                0.10, 0.70, 0.20},
            std::vector<double>{
                0.60, 0.20, 0.20},
            std::vector<double>{
                0.20, 0.30, 0.50}};
    const auto pivot =
        flow::NaturalVariableCompositionPivot3P::
            from_dependent_components(
                3U,
                {1U, 0U, 2U});

    flow::NaturalVariableCellStateInput3P input;
    input.component_ids = {"A", "B", "C"};
    input.reference_pressure_pa = 10.0e6;
    input.temperature_k = 350.0;
    input.independent_saturations =
        {0.20, 0.30};
    input.composition_pivot = pivot;

    for (std::size_t phase = 0U;
         phase < 3U;
         ++phase) {
        input.independent_phase_compositions[
            phase] =
            independent_values(
                compositions[phase],
                pivot.dependent_component(
                    static_cast<
                        flow::PhaseSlot3>(
                            phase)));
        input.phase_properties[phase] = {
            3000.0 +
                100.0 *
                    static_cast<double>(phase),
            700.0 +
                50.0 *
                    static_cast<double>(phase),
            1.0e-3 +
                1.0e-4 *
                    static_cast<double>(phase),
            2.0e5 +
                1.0e4 *
                    static_cast<double>(phase),
            1.5e5 +
                8.0e3 *
                    static_cast<double>(phase)};
    }
    return flow::NaturalVariableCellState3P::create(
        std::move(input));
}

std::vector<double> natural_inputs(
    const flow::NaturalVariableCellState3P& state) {
    const auto& layout =
        state.layout();
    std::vector<double> values(
        layout.unknown_count(),
        0.0);
    values[layout.pressure_unknown_index()] =
        state.reference_pressure_pa();
    values[layout.temperature_unknown_index()] =
        state.temperature_k();
    values[*layout.independent_saturation_unknown_index(
        flow::PhaseSlot3::phase0)] =
        state.phase_saturation(
            flow::PhaseSlot3::phase0);
    values[*layout.independent_saturation_unknown_index(
        flow::PhaseSlot3::phase1)] =
        state.phase_saturation(
            flow::PhaseSlot3::phase1);

    for (std::size_t phase = 0U;
         phase < 3U;
         ++phase) {
        const auto slot =
            static_cast<flow::PhaseSlot3>(
                phase);
        const auto composition =
            state.phase_composition(slot);
        for (std::size_t component = 0U;
             component <
             composition.size();
             ++component) {
            const auto column =
                layout.independent_composition_unknown_index(
                    slot,
                    component);
            if (column.has_value()) {
                values[*column] =
                    composition[component];
            }
        }
    }
    return values;
}

std::array<std::vector<double>, 3>
property_gradient(
    std::size_t q,
    double base_scale) {
    std::array<std::vector<double>, 3>
        result;
    for (std::size_t phase = 0U;
         phase < 3U;
         ++phase) {
        result[phase].assign(q, 0.0);
        for (std::size_t column = 0U;
             column < q;
             ++column) {
            const double column_scale =
                column == 0U
                    ? 1.0e-5
                    : (column == 1U
                           ? 0.1
                           : 0.5);
            result[phase][column] =
                base_scale *
                column_scale *
                static_cast<double>(
                    (phase + 1U) *
                    (column + 1U));
        }
    }
    return result;
}

struct Fixture {
    flow::NaturalVariableCellState3P state;
    std::vector<double> baseline_inputs;
    std::array<std::vector<double>, 3>
        rho_gradient;
    std::array<std::vector<double>, 3>
        h_gradient;
    std::array<std::vector<double>, 3>
        u_gradient;
    std::vector<double> rock_gradient;
    flow::PhaseTransportPropertyNaturalVariableLinearization3P
        transport;
    flow::PhaseCaloricPropertyNaturalVariableLinearization3P
        caloric;
    flow::StationaryRockThermalStorageLinearization3P
        rock;

    Fixture()
        : state(make_state()),
          baseline_inputs(
              natural_inputs(state)),
          rho_gradient(
              property_gradient(
                  state.layout()
                      .unknown_count(),
                  1.0)),
          h_gradient(
              property_gradient(
                  state.layout()
                      .unknown_count(),
                  600.0)),
          u_gradient(
              property_gradient(
                  state.layout()
                      .unknown_count(),
                  400.0)),
          rock_gradient(
              state.layout()
                  .unknown_count(),
              0.0),
          transport(
              flow::
                  make_phase_transport_property_linearization(
                      state,
                      rho_gradient,
                      std::array<
                          std::vector<double>,
                          3>{
                          std::vector<double>(
                              state.layout()
                                  .unknown_count(),
                              0.0),
                          std::vector<double>(
                              state.layout()
                                  .unknown_count(),
                              0.0),
                          std::vector<double>(
                              state.layout()
                                  .unknown_count(),
                              0.0)},
                      {"synthetic-rho",
                       "energy-fixture",
                       "v1"},
                      {"synthetic-mu",
                       "energy-fixture",
                       "v1"})),
          caloric(
              flow::
                  make_phase_caloric_property_linearization(
                      state,
                      h_gradient,
                      u_gradient,
                      {"synthetic-h",
                       "energy-fixture",
                       "v1"},
                      {"synthetic-u",
                       "energy-fixture",
                       "v1"})),
          rock(
              [&] {
                  rock_gradient[
                      state.layout()
                          .temperature_unknown_index()] =
                      2.0e6;
                  return flow::
                      make_stationary_rock_thermal_storage_linearization(
                          state,
                          1.0e8,
                          rock_gradient,
                          {"synthetic-rock-energy",
                           "energy-fixture",
                           "v1"});
              }()) {}
};

double fresh_total_energy(
    const Fixture& fixture,
    std::span<const double> inputs,
    double porosity) {
    const auto& layout =
        fixture.state.layout();
    const double s0 =
        inputs[
            *layout.independent_saturation_unknown_index(
                flow::PhaseSlot3::phase0)];
    const double s1 =
        inputs[
            *layout.independent_saturation_unknown_index(
                flow::PhaseSlot3::phase1)];
    const std::array<double, 3>
        saturation{
            s0,
            s1,
            1.0 - s0 - s1};

    double fluid = 0.0;
    for (std::size_t phase = 0U;
         phase < 3U;
         ++phase) {
        double rho =
            fixture.transport
                .mass_density_kg_per_m3[
                    phase];
        double u =
            fixture.caloric
                .specific_internal_energy_j_per_kg[
                    phase];
        for (std::size_t column = 0U;
             column < inputs.size();
             ++column) {
            const double delta =
                inputs[column] -
                fixture.baseline_inputs[
                    column];
            rho +=
                fixture.rho_gradient[
                    phase][column] *
                delta;
            u +=
                fixture.u_gradient[
                    phase][column] *
                delta;
        }
        fluid +=
            porosity *
            saturation[phase] *
            rho *
            u;
    }

    double rock =
        fixture.rock
            .volumetric_internal_energy_j_per_rock_m3;
    for (std::size_t column = 0U;
         column < inputs.size();
         ++column) {
        rock +=
            fixture.rock_gradient[
                column] *
            (inputs[column] -
             fixture.baseline_inputs[
                 column]);
    }

    return fluid +
        (1.0 - porosity) *
            rock;
}

void accumulation_and_backward_euler() {
    Fixture fixture;
    constexpr double porosity = 0.25;
    constexpr double dt = 20.0;

    const auto snapshot =
        flow::
            build_pore_volume_energy_accumulation_snapshot(
                fixture.state,
                porosity,
                fixture.transport,
                fixture.caloric,
                fixture.rock);
    const auto linearization =
        flow::
            build_pore_volume_energy_accumulation_linearization(
                fixture.state,
                porosity,
                fixture.transport,
                fixture.caloric,
                fixture.rock);

    double expected_fluid = 0.0;
    for (std::size_t phase = 0U;
         phase < 3U;
         ++phase) {
        expected_fluid +=
            porosity *
            fixture.state.phase_saturation(
                static_cast<
                    flow::PhaseSlot3>(
                        phase)) *
            fixture.transport
                .mass_density_kg_per_m3[
                    phase] *
            fixture.caloric
                .specific_internal_energy_j_per_kg[
                    phase];
    }
    const double expected_rock =
        (1.0 - porosity) *
        fixture.rock
            .volumetric_internal_energy_j_per_rock_m3;

    near(
        snapshot.fluid_internal_energy_j_per_bulk_m3,
        expected_fluid);
    near(
        snapshot.rock_internal_energy_j_per_bulk_m3,
        expected_rock);
    near(
        snapshot.total_internal_energy_j_per_bulk_m3,
        expected_fluid +
            expected_rock);
    near(
        linearization.total_internal_energy_j_per_bulk_m3,
        snapshot.total_internal_energy_j_per_bulk_m3);

    auto previous = snapshot;
    previous.fluid_internal_energy_j_per_bulk_m3 -=
        5000.0;
    previous.total_internal_energy_j_per_bulk_m3 -=
        5000.0;

    const auto residual =
        flow::
            build_backward_euler_energy_accumulation_residual(
                snapshot,
                linearization,
                previous,
                dt);
    near(
        residual.residual_w_per_bulk_m3,
        5000.0 / dt);
    require(
        residual.input_count ==
                fixture.state.layout()
                    .unknown_count() &&
            residual.gradient.size() ==
                residual.input_count,
        "energy backward-Euler gradient shape changed");

    for (std::size_t column = 0U;
         column < residual.input_count;
         ++column) {
        near(
            residual.d_residual(column),
            linearization
                .total_internal_energy_gradient[
                    column] /
                dt);
    }
}

void jacobian_fresh_perturbation() {
    Fixture fixture;
    constexpr double porosity = 0.25;

    const auto linearization =
        flow::
            build_pore_volume_energy_accumulation_linearization(
                fixture.state,
                porosity,
                fixture.transport,
                fixture.caloric,
                fixture.rock);

    const std::size_t q =
        fixture.state.layout()
            .unknown_count();

    for (std::size_t column = 0U;
         column < q;
         ++column) {
        const double step =
            column ==
                    fixture.state.layout()
                        .pressure_unknown_index()
                ? 10.0
                : (column ==
                           fixture.state.layout()
                               .temperature_unknown_index()
                       ? 1.0e-4
                       : 1.0e-7);

        auto plus =
            fixture.baseline_inputs;
        auto minus =
            fixture.baseline_inputs;
        plus[column] += step;
        minus[column] -= step;

        const double finite_difference =
            (fresh_total_energy(
                 fixture,
                 plus,
                 porosity) -
             fresh_total_energy(
                 fixture,
                 minus,
                 porosity)) /
            (2.0 * step);

        near(
            linearization
                .total_internal_energy_gradient[
                    column],
            finite_difference,
            2.0e-7,
            2.0e-3);
    }
}

void zero_residual_nonzero_jacobian() {
    Fixture fixture;
    constexpr double porosity = 0.25;

    const auto snapshot =
        flow::
            build_pore_volume_energy_accumulation_snapshot(
                fixture.state,
                porosity,
                fixture.transport,
                fixture.caloric,
                fixture.rock);
    const auto linearization =
        flow::
            build_pore_volume_energy_accumulation_linearization(
                fixture.state,
                porosity,
                fixture.transport,
                fixture.caloric,
                fixture.rock);
    const auto residual =
        flow::
            build_backward_euler_energy_accumulation_residual(
                snapshot,
                linearization,
                snapshot,
                10.0);

    near(
        residual.residual_w_per_bulk_m3,
        0.0,
        0.0,
        0.0);
    require(
        std::any_of(
            residual.gradient.begin(),
            residual.gradient.end(),
            [](double value) {
                return value != 0.0;
            }),
        "zero energy accumulation residual erased current-state Jacobian");
}

void invalid_inputs() {
    Fixture fixture;
    constexpr double porosity = 0.25;

    expect_invalid(
        [&] {
            (void)flow::
                build_pore_volume_energy_accumulation_snapshot(
                    fixture.state,
                    1.0,
                    fixture.transport,
                    fixture.caloric,
                    fixture.rock);
        },
        "porosity");

    expect_invalid(
        [&] {
            auto bad_gradient =
                fixture.rock_gradient;
            bad_gradient.pop_back();
            (void)flow::
                make_stationary_rock_thermal_storage_linearization(
                    fixture.state,
                    1.0e8,
                    std::move(
                        bad_gradient),
                    {"synthetic-rock-energy",
                     "energy-fixture",
                     "v1"});
        },
        "malformed");

    const auto snapshot =
        flow::
            build_pore_volume_energy_accumulation_snapshot(
                fixture.state,
                porosity,
                fixture.transport,
                fixture.caloric,
                fixture.rock);
    const auto linearization =
        flow::
            build_pore_volume_energy_accumulation_linearization(
                fixture.state,
                porosity,
                fixture.transport,
                fixture.caloric,
                fixture.rock);

    expect_invalid(
        [&] {
            auto bad_previous =
                snapshot;
            bad_previous.total_internal_energy_j_per_bulk_m3 +=
                100.0;
            (void)flow::
                build_backward_euler_energy_accumulation_residual(
                    snapshot,
                    linearization,
                    bad_previous,
                    10.0);
        },
        "does not close");

    expect_invalid(
        [&] {
            auto bad_linearization =
                linearization;
            bad_linearization
                .total_internal_energy_gradient[0] =
                std::numeric_limits<double>::infinity();
            (void)flow::
                build_backward_euler_energy_accumulation_residual(
                    snapshot,
                    bad_linearization,
                    snapshot,
                    10.0);
        },
        "non-finite");
}

void headers() {
    require(
        energy_accumulation_header(),
        "energy accumulation header probe failed");
}

using Test =
    std::pair<std::string_view, void (*)()>;

constexpr Test tests[]{
    {"accumulation_and_backward_euler",
     accumulation_and_backward_euler},
    {"jacobian_fresh_perturbation",
     jacobian_fresh_perturbation},
    {"zero_residual_nonzero_jacobian",
     zero_residual_nonzero_jacobian},
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
