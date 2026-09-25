#include <mpmc/well_discretization/component_molar_rate.hpp>

#include <algorithm>
#include <cmath>
#include <limits>
#include <source_location>
#include <stdexcept>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

namespace {

namespace flow = mpmc::flow;
namespace well = mpmc::well;
namespace wd = mpmc::well_discretization;

void require_component(
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

void near_component(
    double actual,
    double expected,
    double relative = 5.0e-12,
    double absolute = 5.0e-28,
    std::source_location where =
        std::source_location::current()) {
    const double scale =
        std::max(
            {1.0e-30,
             std::abs(actual),
             std::abs(expected)});
    require_component(
        std::isfinite(actual) &&
            std::isfinite(expected) &&
            std::abs(actual - expected) <=
                absolute + relative * scale,
        "component molar-rate numeric mismatch",
        where);
}

template <class Function>
void expect_component_invalid(
    Function&& function) {
    bool caught = false;
    try {
        std::forward<Function>(function)();
    } catch (const std::invalid_argument&) {
        caught = true;
    }
    require_component(
        caught,
        "expected invalid_argument");
}

well::PeacemanWellIndex3D
component_connection() {
    return well::make_peaceman_well_index_3d(
        {20.0, 10.0, 5.0},
        {4.0e-13, 1.0e-13, 8.0e-14},
        well::AxisAlignedWellDirection3D::z,
        0.10,
        0.0);
}

flow::LocalPhaseMobilityLinearization3P
component_mobility() {
    flow::LocalPhaseMobilityLinearization3P
        result{};
    result.state_identity.layout =
        flow::NaturalVariableLayoutDescriptor{
            2U,
            3U,
            {1U, 1U, 1U}};
    result.state_identity.component_ids =
        {"A", "B"};
    result.state_identity.reference_pressure_pa =
        1.0e6;
    result.state_identity.temperature_k =
        350.0;
    result.state_identity.saturation =
        {0.2, 0.3, 0.5};
    result.state_identity.phase_composition = {
        std::vector<double>{0.7, 0.3},
        std::vector<double>{0.5, 0.5},
        std::vector<double>{0.2, 0.8}};

    const std::size_t q =
        result.state_identity.layout
            .unknown_count();
    result.phase_pressure_pa =
        {1.20e6, 1.10e6, 0.90e6};
    result.mobility_per_pa_s =
        {1.0e-9, 2.0e-9, 3.0e-9};

    for (std::size_t phase = 0U;
         phase < flow::fixed_three_phase_count;
         ++phase) {
        result.phase_pressure_gradient[phase]
            .assign(q, 0.0);
        result.phase_pressure_gradient[phase][0] =
            1.0;
        result.phase_pressure_gradient[phase][2] =
            100.0 *
            static_cast<double>(
                phase + 1U);

        result.mobility_gradient[phase]
            .resize(q);
        for (std::size_t column = 0U;
             column < q;
             ++column) {
            result.mobility_gradient[phase][column] =
                static_cast<double>(
                    (phase + 1U) *
                    (column + 1U)) *
                1.0e-12;
        }
    }
    return result;
}

flow::PhaseMolarDensityNaturalVariableLinearization3P
component_density() {
    const flow::NaturalVariableLayout3P layout{
        flow::NaturalVariableCompositionPivot3P::
            from_dependent_components(
                2U,
                {1U, 1U, 1U})};

    flow::PhaseMolarDensityNaturalVariableLinearization3P
        result{
            layout,
            {10000.0, 12000.0, 8000.0},
            {}};

    const std::size_t q =
        layout.unknown_count();
    for (std::size_t phase = 0U;
         phase < flow::fixed_three_phase_count;
         ++phase) {
        result.gradient[phase].resize(q);
        for (std::size_t column = 0U;
             column < q;
             ++column) {
            result.gradient[phase][column] =
                5.0 *
                static_cast<double>(
                    (phase + 1U) *
                    (column + 1U));
        }
    }
    return result;
}

void mixed_phase_component_closure() {
    const auto wi =
        component_connection();
    const auto mobility =
        component_mobility();
    const auto density =
        component_density();
    const auto result =
        wd::make_connection_component_molar_rate_linearization_3p(
            wi,
            mobility,
            density,
            1.00e6);

    require_component(
        result.component_ids ==
                std::vector<std::string>{"A", "B"} &&
            result.component_count() == 2U &&
            result.input_count == 7U &&
            result.phase_volumetric_rate_m3_per_s[0] >
                0.0 &&
            result.phase_volumetric_rate_m3_per_s[1] >
                0.0 &&
            result.phase_volumetric_rate_m3_per_s[2] <
                0.0,
        "mixed production/injection phase signs were not preserved");

    double expected_total = 0.0;
    std::vector<double> expected_components(
        2U,
        0.0);
    for (std::size_t phase = 0U;
         phase < 3U;
         ++phase) {
        const auto rate =
            wd::make_pressure_drawdown_phase_rate_linearization_3p(
                wi,
                mobility,
                static_cast<flow::PhaseSlot3>(
                    phase),
                1.00e6);
        const double phase_molar =
            rate.phase_volumetric_rate_m3_per_s *
            density.molar_density_mol_per_m3[
                phase];
        expected_total +=
            phase_molar;
        for (std::size_t component = 0U;
             component < 2U;
             ++component) {
            expected_components[component] +=
                phase_molar *
                mobility
                    .state_identity
                    .phase_composition[phase]
                                      [component];
        }
    }

    near_component(
        result.total_molar_rate_mol_per_s,
        expected_total);
    near_component(
        result.component_rate(0U),
        expected_components[0]);
    near_component(
        result.component_rate(1U),
        expected_components[1]);
    near_component(
        result.component_rate(0U) +
            result.component_rate(1U),
        result.total_molar_rate_mol_per_s);
}

void reservoir_product_rule_and_composition_pivot() {
    const auto wi =
        component_connection();
    const auto mobility =
        component_mobility();
    const auto density =
        component_density();
    const auto result =
        wd::make_connection_component_molar_rate_linearization_3p(
            wi,
            mobility,
            density,
            1.00e6);

    constexpr std::size_t component = 0U;
    constexpr std::size_t pressure_column = 0U;
    double expected_pressure = 0.0;
    for (std::size_t phase = 0U;
         phase < 3U;
         ++phase) {
        const auto rate =
            wd::make_pressure_drawdown_phase_rate_linearization_3p(
                wi,
                mobility,
                static_cast<flow::PhaseSlot3>(
                    phase),
                1.00e6);
        const double c =
            density.molar_density_mol_per_m3[
                phase];
        const double x =
            mobility.state_identity
                .phase_composition[phase]
                                  [component];
        expected_pressure +=
            rate.d_phase_rate_d_reservoir_unknown(
                pressure_column) *
                c *
                x +
            rate.phase_volumetric_rate_m3_per_s *
                density.gradient[phase]
                                [pressure_column] *
                x;
    }
    near_component(
        result.d_component_rate_d_reservoir_unknown(
            component,
            pressure_column),
        expected_pressure);

    // Column 4 is phase0 component-A independent composition because B is the
    // frozen dependent component. Therefore dx_A/dq4=+1 and dx_B/dq4=-1.
    constexpr std::size_t composition_column = 4U;
    const auto phase0_rate =
        wd::make_pressure_drawdown_phase_rate_linearization_3p(
            wi,
            mobility,
            flow::PhaseSlot3::phase0,
            1.00e6);
    const double extra_pivot_term =
        phase0_rate.phase_volumetric_rate_m3_per_s *
        density.molar_density_mol_per_m3[0];

    const double component0 =
        result.d_component_rate_d_reservoir_unknown(
            0U,
            composition_column);
    const double component1 =
        result.d_component_rate_d_reservoir_unknown(
            1U,
            composition_column);
    near_component(
        component0 + component1,
        result.total_reservoir_gradient[
            composition_column]);

    // Remove the identical non-composition product-rule pieces by subtracting
    // the two-component sum closure; the pivot itself must contribute
    // +q*c to A and -q*c to B.
    double base0 = 0.0;
    double base1 = 0.0;
    for (std::size_t phase = 0U;
         phase < 3U;
         ++phase) {
        const auto rate =
            wd::make_pressure_drawdown_phase_rate_linearization_3p(
                wi,
                mobility,
                static_cast<flow::PhaseSlot3>(
                    phase),
                1.00e6);
        const double c =
            density.molar_density_mol_per_m3[
                phase];
        const double dc =
            density.gradient[phase]
                            [composition_column];
        const double q =
            rate.phase_volumetric_rate_m3_per_s;
        const double dq =
            rate.d_phase_rate_d_reservoir_unknown(
                composition_column);
        base0 +=
            dq * c *
                mobility.state_identity
                    .phase_composition[phase][0] +
            q * dc *
                mobility.state_identity
                    .phase_composition[phase][0];
        base1 +=
            dq * c *
                mobility.state_identity
                    .phase_composition[phase][1] +
            q * dc *
                mobility.state_identity
                    .phase_composition[phase][1];
    }
    near_component(
        component0 - base0,
        extra_pivot_term);
    near_component(
        component1 - base1,
        -extra_pivot_term);
}

void bhp_derivative_closure() {
    const auto wi =
        component_connection();
    const auto mobility =
        component_mobility();
    const auto density =
        component_density();
    const auto result =
        wd::make_connection_component_molar_rate_linearization_3p(
            wi,
            mobility,
            density,
            1.00e6);

    std::vector<double> expected(
        2U,
        0.0);
    double expected_total = 0.0;
    for (std::size_t phase = 0U;
         phase < 3U;
         ++phase) {
        const auto rate =
            wd::make_pressure_drawdown_phase_rate_linearization_3p(
                wi,
                mobility,
                static_cast<flow::PhaseSlot3>(
                    phase),
                1.00e6);
        const double phase_bhp =
            rate.d_phase_rate_d_bottom_hole_pressure() *
            density.molar_density_mol_per_m3[
                phase];
        expected_total +=
            phase_bhp;
        for (std::size_t component = 0U;
             component < 2U;
             ++component) {
            expected[component] +=
                phase_bhp *
                mobility.state_identity
                    .phase_composition[phase]
                                      [component];
        }
    }

    near_component(
        result.d_component_rate_d_bottom_hole_pressure(
            0U),
        expected[0]);
    near_component(
        result.d_component_rate_d_bottom_hole_pressure(
            1U),
        expected[1]);
    near_component(
        result.total_bhp_derivative_mol_per_pa_s,
        expected_total);
    near_component(
        result.d_component_rate_d_bottom_hole_pressure(
            0U) +
            result.d_component_rate_d_bottom_hole_pressure(
                1U),
        result.total_bhp_derivative_mol_per_pa_s);
}

void zero_phase_rate_preserves_component_derivative() {
    const auto wi =
        component_connection();
    auto mobility =
        component_mobility();
    const auto density =
        component_density();

    mobility.mobility_per_pa_s[0] =
        0.0;
    mobility.mobility_per_pa_s[2] =
        0.0;
    for (std::size_t column = 0U;
         column < 7U;
         ++column) {
        mobility.mobility_gradient[0][column] =
            0.0;
        mobility.mobility_gradient[2][column] =
            0.0;
    }

    const double bhp =
        mobility.phase_pressure_pa[1];
    const auto result =
        wd::make_connection_component_molar_rate_linearization_3p(
            wi,
            mobility,
            density,
            bhp);

    require_component(
        result.total_molar_rate_mol_per_s ==
            0.0 &&
        result.component_rate(0U) ==
            0.0 &&
        result.component_rate(1U) ==
            0.0,
        "zero phase drawdown did not preserve exact zero component molar rates");

    const auto phase1 =
        wd::make_pressure_drawdown_phase_rate_linearization_3p(
            wi,
            mobility,
            flow::PhaseSlot3::phase1,
            bhp);
    const double expected =
        phase1.d_phase_rate_d_reservoir_unknown(
            0U) *
        density.molar_density_mol_per_m3[1] *
        mobility.state_identity
            .phase_composition[1][0];
    near_component(
        result.d_component_rate_d_reservoir_unknown(
            0U,
            0U),
        expected);
    require_component(
        result.d_component_rate_d_reservoir_unknown(
            0U,
            0U) != 0.0,
        "zero phase rate incorrectly erased component reservoir derivative");
}

void component_invalid_inputs() {
    const auto wi =
        component_connection();
    const auto mobility =
        component_mobility();

    {
        auto density =
            component_density();
        density.layout =
            flow::NaturalVariableLayout3P{
                3U};
        expect_component_invalid(
            [&] {
                (void)wd::
                    make_connection_component_molar_rate_linearization_3p(
                        wi,
                        mobility,
                        density,
                        1.0e6);
            });
    }
    {
        auto density =
            component_density();
        density.molar_density_mol_per_m3[1] =
            0.0;
        expect_component_invalid(
            [&] {
                (void)wd::
                    make_connection_component_molar_rate_linearization_3p(
                        wi,
                        mobility,
                        density,
                        1.0e6);
            });
    }
    {
        auto density =
            component_density();
        density.gradient[1].pop_back();
        expect_component_invalid(
            [&] {
                (void)wd::
                    make_connection_component_molar_rate_linearization_3p(
                        wi,
                        mobility,
                        density,
                        1.0e6);
            });
    }
    {
        auto density =
            component_density();
        density.gradient[1][0] =
            std::numeric_limits<double>::
                quiet_NaN();
        expect_component_invalid(
            [&] {
                (void)wd::
                    make_connection_component_molar_rate_linearization_3p(
                        wi,
                        mobility,
                        density,
                        1.0e6);
            });
    }
    {
        auto bad_mobility =
            mobility;
        bad_mobility.state_identity
            .phase_composition[0] =
            {0.8, 0.3};
        expect_component_invalid(
            [&] {
                (void)wd::
                    make_connection_component_molar_rate_linearization_3p(
                        wi,
                        bad_mobility,
                        component_density(),
                        1.0e6);
            });
    }
    {
        auto bad_mobility =
            mobility;
        bad_mobility.state_identity
            .phase_composition[0]
            .pop_back();
        expect_component_invalid(
            [&] {
                (void)wd::
                    make_connection_component_molar_rate_linearization_3p(
                        wi,
                        bad_mobility,
                        component_density(),
                        1.0e6);
            });
    }
}

} // namespace

bool run_component_molar_rate_case(
    std::string_view name) {
    if (name ==
        "component_molar_rate_mixed_phase_closure") {
        mixed_phase_component_closure();
    } else if (
        name ==
        "component_molar_rate_reservoir_product_rule_and_composition_pivot") {
        reservoir_product_rule_and_composition_pivot();
    } else if (
        name ==
        "component_molar_rate_bhp_derivative_closure") {
        bhp_derivative_closure();
    } else if (
        name ==
        "component_molar_rate_zero_phase_rate_preserves_derivative") {
        zero_phase_rate_preserves_component_derivative();
    } else if (
        name ==
        "component_molar_rate_invalid_inputs") {
        component_invalid_inputs();
    } else {
        return false;
    }
    return true;
}
