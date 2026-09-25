#include <mpmc/well_discretization/cell_source_adapter.hpp>

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

namespace fd = mpmc::flow_discretization;
namespace flow = mpmc::flow;
namespace wd = mpmc::well_discretization;

void require_adapter(
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

void near_adapter(
    double actual,
    double expected,
    double relative = 5.0e-13,
    double absolute = 5.0e-25,
    std::source_location where =
        std::source_location::current()) {
    const double scale =
        std::max(
            {1.0,
             std::abs(actual),
             std::abs(expected)});
    require_adapter(
        std::isfinite(actual) &&
            std::isfinite(expected) &&
            std::abs(actual - expected) <=
                absolute + relative * scale,
        "cell-source adapter numeric mismatch",
        where);
}

template <class Function>
void expect_adapter_invalid(
    Function&& function) {
    bool caught = false;
    try {
        std::forward<Function>(function)();
    } catch (const std::invalid_argument&) {
        caught = true;
    }
    require_adapter(
        caught,
        "expected invalid_argument");
}

flow::NaturalVariableStateIdentity3P
adapter_identity() {
    flow::NaturalVariableStateIdentity3P
        result{};
    result.layout =
        flow::NaturalVariableLayoutDescriptor{
            2U,
            3U,
            {1U, 1U, 1U}};
    result.component_ids =
        {"A", "B"};
    result.reference_pressure_pa =
        1.0e6;
    result.temperature_k =
        350.0;
    result.saturation =
        {0.2, 0.3, 0.5};
    result.phase_composition = {
        std::vector<double>{0.7, 0.3},
        std::vector<double>{0.5, 0.5},
        std::vector<double>{0.2, 0.8}};
    return result;
}

wd::WellConnectionComponentMolarRateLinearization3P
adapter_component_rate() {
    const auto identity =
        adapter_identity();

    wd::WellConnectionComponentMolarRateLinearization3P
        result{};
    result.state_identity =
        identity;
    result.bottom_hole_pressure_pa =
        0.9e6;
    result.component_ids =
        identity.component_ids;
    result.input_count =
        identity.layout.unknown_count();
    result.component_molar_rate_mol_per_s =
        {10.0, -4.0};
    result.component_reservoir_jacobian.resize(
        2U * result.input_count);
    for (std::size_t index = 0U;
         index <
         result.component_reservoir_jacobian
             .size();
         ++index) {
        result.component_reservoir_jacobian[
            index] =
            static_cast<double>(
                index + 1U) *
            0.25;
    }
    result.component_bhp_derivative_mol_per_pa_s =
        {-0.5, 0.25};
    return result;
}

wd::WellConnectionAdvectiveEnergyRateLinearization3P
adapter_energy_rate() {
    const auto identity =
        adapter_identity();

    wd::WellConnectionAdvectiveEnergyRateLinearization3P
        result{};
    result.state_identity =
        identity;
    result.bottom_hole_pressure_pa =
        0.9e6;
    result.total_advective_energy_rate_w =
        50.0;
    result.reservoir_natural_variable_gradient.resize(
        identity.layout.unknown_count());
    for (std::size_t column = 0U;
         column <
         result.reservoir_natural_variable_gradient
             .size();
         ++column) {
        result.reservoir_natural_variable_gradient[
            column] =
            static_cast<double>(
                column + 1U) *
            1.5;
    }
    result.bottom_hole_pressure_derivative_w_per_pa =
        -3.0;
    return result;
}

void sign_bridge_and_reservoir_jacobian() {
    const auto component =
        adapter_component_rate();
    const auto energy =
        adapter_energy_rate();
    const auto adapted =
        wd::make_connection_cell_source_adapter_3p(
            component,
            energy,
            "well/test-connection");

    const auto& source =
        adapted.cell_source;
    require_adapter(
        source.provenance ==
                "well/test-connection" &&
            source.component_ids ==
                std::vector<std::string>{"A", "B"} &&
            source.input_count == 7U,
        "cell-source adapter metadata mismatch");

    near_adapter(
        source.component_molar_rate_mol_per_s[0],
        -10.0);
    near_adapter(
        source.component_molar_rate_mol_per_s[1],
        4.0);
    near_adapter(
        source.energy_rate_w,
        -50.0);

    for (std::size_t index = 0U;
         index <
         component.component_reservoir_jacobian
             .size();
         ++index) {
        near_adapter(
            source
                .component_molar_rate_jacobian_mol_per_s[
                    index],
            -component
                 .component_reservoir_jacobian[
                     index]);
    }
    for (std::size_t column = 0U;
         column < source.input_count;
         ++column) {
        near_adapter(
            source.energy_rate_gradient_w[
                column],
            -energy
                 .reservoir_natural_variable_gradient[
                     column]);
    }
}

void bhp_sidecar_sign_bridge() {
    const auto component =
        adapter_component_rate();
    const auto energy =
        adapter_energy_rate();
    const auto adapted =
        wd::make_connection_cell_source_adapter_3p(
            component,
            energy,
            "well/test-connection");

    near_adapter(
        adapted
            .d_component_source_d_bottom_hole_pressure(
                0U),
        0.5);
    near_adapter(
        adapted
            .d_component_source_d_bottom_hole_pressure(
                1U),
        -0.25);
    near_adapter(
        adapted
            .d_energy_source_d_bottom_hole_pressure(),
        3.0);
}

void normalized_residual_has_no_double_negation() {
    const auto component =
        adapter_component_rate();
    const auto energy =
        adapter_energy_rate();
    const auto adapted =
        wd::make_connection_cell_source_adapter_3p(
            component,
            energy,
            "well/test-connection");
    const auto normalized =
        fd::normalize_cell_source_by_bulk_volume(
            adapted.cell_source,
            2.0);

    // CellSource normalization is residual = -source/V. Since adapter source
    // is -well_rate, the final residual is +well_rate/V for production-positive
    // well convention.
    near_adapter(
        normalized
            .component_residual_mol_per_bulk_m3_s[
                0],
        5.0);
    near_adapter(
        normalized
            .component_residual_mol_per_bulk_m3_s[
                1],
        -2.0);
    near_adapter(
        normalized
            .energy_residual_w_per_bulk_m3,
        25.0);

    for (std::size_t index = 0U;
         index <
         component.component_reservoir_jacobian
             .size();
         ++index) {
        near_adapter(
            normalized
                .component_jacobian_mol_per_bulk_m3_s[
                    index],
            component
                .component_reservoir_jacobian[
                    index] /
                2.0);
    }
    for (std::size_t column = 0U;
         column <
         energy
             .reservoir_natural_variable_gradient
             .size();
         ++column) {
        near_adapter(
            normalized
                .energy_gradient_w_per_bulk_m3[
                    column],
            energy
                .reservoir_natural_variable_gradient[
                    column] /
                2.0);
    }
}

void adapter_invalid_identity_and_inputs() {
    {
        auto energy =
            adapter_energy_rate();
        energy.state_identity.temperature_k +=
            1.0;
        expect_adapter_invalid(
            [&] {
                (void)wd::
                    make_connection_cell_source_adapter_3p(
                        adapter_component_rate(),
                        energy,
                        "well/test");
            });
    }
    {
        auto energy =
            adapter_energy_rate();
        energy.bottom_hole_pressure_pa +=
            1.0;
        expect_adapter_invalid(
            [&] {
                (void)wd::
                    make_connection_cell_source_adapter_3p(
                        adapter_component_rate(),
                        energy,
                        "well/test");
            });
    }
    {
        auto component =
            adapter_component_rate();
        component.component_ids =
            {"B", "A"};
        expect_adapter_invalid(
            [&] {
                (void)wd::
                    make_connection_cell_source_adapter_3p(
                        component,
                        adapter_energy_rate(),
                        "well/test");
            });
    }
    {
        auto component =
            adapter_component_rate();
        component
            .component_reservoir_jacobian[0] =
            std::numeric_limits<double>::
                quiet_NaN();
        expect_adapter_invalid(
            [&] {
                (void)wd::
                    make_connection_cell_source_adapter_3p(
                        component,
                        adapter_energy_rate(),
                        "well/test");
            });
    }
    {
        auto energy =
            adapter_energy_rate();
        energy
            .bottom_hole_pressure_derivative_w_per_pa =
            std::numeric_limits<double>::
                quiet_NaN();
        expect_adapter_invalid(
            [&] {
                (void)wd::
                    make_connection_cell_source_adapter_3p(
                        adapter_component_rate(),
                        energy,
                        "well/test");
            });
    }
    expect_adapter_invalid(
        [&] {
            (void)wd::
                make_connection_cell_source_adapter_3p(
                    adapter_component_rate(),
                    adapter_energy_rate(),
                    "");
        });
}

} // namespace

bool run_cell_source_adapter_case(
    std::string_view name) {
    if (name ==
        "cell_source_adapter_sign_and_reservoir_jacobian") {
        sign_bridge_and_reservoir_jacobian();
    } else if (
        name ==
        "cell_source_adapter_bhp_sidecar") {
        bhp_sidecar_sign_bridge();
    } else if (
        name ==
        "cell_source_adapter_normalized_residual_sign") {
        normalized_residual_has_no_double_negation();
    } else if (
        name ==
        "cell_source_adapter_invalid_identity_and_inputs") {
        adapter_invalid_identity_and_inputs();
    } else {
        return false;
    }
    return true;
}
