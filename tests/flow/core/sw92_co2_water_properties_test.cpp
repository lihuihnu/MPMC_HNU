#include <mpmc/flow/sw92_co2_water_properties.hpp>

#include "test_support.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <span>
#include <stdexcept>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

namespace {

namespace th = mpmc::thermodynamics;
namespace flow = mpmc::flow;
namespace ad = mpmc::ad;
namespace st = sw92_test;

void require_provider(
    bool condition,
    std::string_view message) {
    if (!condition) {
        throw std::runtime_error(
            std::string{message});
    }
}

void near_provider(
    double actual,
    double expected,
    double relative = 3.0e-9,
    double absolute = 3.0e-12) {
    const double scale =
        std::max(
            {1.0,
             std::abs(actual),
             std::abs(expected)});
    require_provider(
        std::isfinite(actual) &&
            std::isfinite(expected) &&
            std::abs(actual - expected) <=
                absolute +
                    relative * scale,
        "SW92 CO2/H2O sourced property numeric mismatch");
}

th::Provenance nist_source(
    std::string locator) {
    return {
        th::SourceKind::database,
        "NIST Chemistry WebBook SRD 69",
        "accessed-2026-09-25",
        std::move(locator),
        "",
        "Public NIST Chemistry WebBook record",
        "Use as source data; no NIST code redistributed"};
}

st::Prepared sourced_input(
    bool reverse = false) {
    auto prepared =
        st::binary_input(st::co2, reverse);
    for (auto& component : prepared.catalog) {
        const double value =
            component.id == st::co2.id
                ? 0.0440095
                : 0.0180153;
        component.molar_mass =
            st::scalar(
                value,
                th::Unit::kilogram_per_mole,
                nist_source(
                    component.id +
                    std::string{
                        " molecular weight"}),
                "kg/mol",
                "identity");
    }
    return prepared;
}

th::Sw92ParameterSet sourced_parameters(
    bool reverse = false) {
    auto prepared =
        sourced_input(reverse);
    return th::Sw92ParameterSet::create(
        prepared.catalog,
        prepared.order,
        prepared.input);
}

template <class Error, class Function>
void expect_provider_error(
    Function&& function) {
    bool caught = false;
    try {
        std::forward<Function>(function)();
    } catch (const Error&) {
        caught = true;
    }
    require_provider(
        caught,
        "expected sourced provider failure missing");
}

[[nodiscard]]
std::size_t selected_na_root(
    const th::Sw92Phase<double>& model,
    double pressure_pa,
    double temperature_k,
    std::span<const double> composition) {
    th::Sw92PhaseWorkspace<double>
        workspace;
    const auto roots =
        model.roots(
            pressure_pa,
            temperature_k,
            composition,
            0.0,
            th::SwPhaseFamily::nonaqueous,
            workspace);
    require_provider(
        roots.status ==
                th::Sw92RootStatus::success &&
            roots.count > 0U,
        "sourced provider NA root fixture unavailable");
    return roots.count - 1U;
}

void sourced_formula_oracles() {
    auto parameters =
        sourced_parameters();
    auto model =
        th::Sw92Phase<double>::
            from_parameters(parameters);
    flow::Sw92Co2WaterPropertyProvider
        provider{model};

    constexpr std::array<double, 2>
        composition{0.7, 0.3};
    const double ideal =
        provider
            .ideal_gas_molar_enthalpy_j_per_mol(
                550.0,
                std::span<const double>{
                    composition});
    // Independently evaluated NIST Shomate H-H_298.15:
    // CO2=10572.267471335748 J/mol,
    // H2O=8698.84649871969 J/mol.
    near_provider(
        ideal,
        10010.24117955093,
        2.0e-12,
        2.0e-9);

    const double viscosity =
        provider.chung_dynamic_viscosity_pa_s(
            550.0,
            std::span<const double>{
                composition},
            5000.0);
    // Independent Chung/Poling equation re-evaluation using:
    // CO2 Vc=91.9 cm3/mol; H2O Vc from IAPWS rho_c=322 kg/m3;
    // H2O dipole=1.8546 D and kappa=0.076.
    near_provider(
        viscosity,
        2.6441732545837795e-05,
        3.0e-10,
        2.0e-13);

    const auto p =
        flow::Sw92Co2WaterPropertyProvider::
            provenance();
    require_provider(
        p.viscosity.model.find(
            "Chung") !=
                std::string::npos &&
            p.enthalpy.model.find(
                "NIST") !=
                std::string::npos,
        "source-backed provider provenance changed");
}

void permutation_and_source_guards() {
    auto parameters =
        sourced_parameters(false);
    auto reversed_parameters =
        sourced_parameters(true);
    auto model =
        th::Sw92Phase<double>::
            from_parameters(parameters);
    auto reversed_model =
        th::Sw92Phase<double>::
            from_parameters(
                reversed_parameters);
    flow::Sw92Co2WaterPropertyProvider
        provider{model};
    flow::Sw92Co2WaterPropertyProvider
        reversed{reversed_model};

    constexpr std::array<double, 2>
        x{0.7, 0.3};
    constexpr std::array<double, 2>
        xr{0.3, 0.7};
    near_provider(
        provider
            .ideal_gas_molar_enthalpy_j_per_mol(
                600.0,
                std::span<const double>{x}),
        reversed
            .ideal_gas_molar_enthalpy_j_per_mol(
                600.0,
                std::span<const double>{xr}),
        2.0e-12);
    near_provider(
        provider.chung_dynamic_viscosity_pa_s(
            600.0,
            std::span<const double>{x},
            4000.0),
        reversed.chung_dynamic_viscosity_pa_s(
            600.0,
            std::span<const double>{xr},
            4000.0),
        2.0e-11);

    expect_provider_error<std::domain_error>(
        [&] {
            (void)provider
                .ideal_gas_molar_enthalpy_j_per_mol(
                    499.99,
                    std::span<const double>{x});
        });
    expect_provider_error<std::domain_error>(
        [&] {
            (void)provider
                .ideal_gas_molar_enthalpy_j_per_mol(
                    1200.01,
                    std::span<const double>{x});
        });

    auto prepared =
        st::binary_input(st::co2);
    auto missing =
        th::Sw92ParameterSet::create(
            prepared.catalog,
            prepared.order,
            prepared.input);
    auto missing_model =
        th::Sw92Phase<double>::
            from_parameters(missing);
    expect_provider_error<std::invalid_argument>(
        [&] {
            (void)flow::
                Sw92Co2WaterPropertyProvider{
                    missing_model};
        });
}

void closure_and_ad_regression() {
    auto parameters =
        sourced_parameters();
    auto model =
        th::Sw92Phase<double>::
            from_parameters(parameters);

    constexpr double pressure = 5.0e6;
    constexpr double temperature = 550.0;
    constexpr std::array<double, 2>
        composition{0.7, 0.3};
    const std::size_t root =
        selected_na_root(
            model,
            pressure,
            temperature,
            std::span<const double>{
                composition});
    const th::Sw92SelectedPhase<double>
        selection{
            0.0,
            th::SwPhaseFamily::nonaqueous,
            root,
            {}};

    auto closure =
        flow::
            make_sw92_co2_water_selected_phase_property_closure(
                model,
                std::vector<
                    th::Sw92SelectedPhase<double>>{
                    selection});

    flow::NaturalVariableLayoutDescriptor
        layout{
            2U,
            1U,
            std::vector<std::size_t>{1U}};
    std::vector<double> q(
        layout.unknown_count(),
        0.0);
    q[layout.pressure_unknown_index()] =
        pressure;
    q[layout.temperature_unknown_index()] =
        temperature;
    const auto x_column =
        layout.independent_composition_unknown_index(
            flow::PhaseSlot3::phase0,
            0U);
    require_provider(
        x_column.has_value(),
        "sourced provider composition column missing");
    q[*x_column] = composition[0];

    ad::RuntimeJacobianWorkspace<double, 4U>
        workspace;
    const auto value =
        flow::
            evaluate_sw92_selected_phase_property_chart(
                layout,
                closure.component_ids(),
                q,
                closure,
                workspace);
    require_provider(
        value.phase_properties.size() == 1U &&
            value.phase_properties[0]
                    .dynamic_viscosity_pa_s >
                0.0 &&
            std::isfinite(
                value.phase_properties[0]
                    .specific_enthalpy_j_per_kg),
        "sourced provider did not publish usable flow properties");

    const std::array<std::size_t, 3>
        columns{
            layout.pressure_unknown_index(),
            layout.temperature_unknown_index(),
            *x_column};
    const std::array<double, 3>
        steps{
            10.0,
            1.0e-4,
            1.0e-6};

    for (std::size_t check = 0U;
         check < columns.size();
         ++check) {
        auto plus = q;
        auto minus = q;
        plus[columns[check]] +=
            steps[check];
        minus[columns[check]] -=
            steps[check];
        const auto pv =
            flow::
                evaluate_sw92_selected_phase_property_chart(
                    layout,
                    closure.component_ids(),
                    plus,
                    closure,
                    workspace);
        const auto mv =
            flow::
                evaluate_sw92_selected_phase_property_chart(
                    layout,
                    closure.component_ids(),
                    minus,
                    closure,
                    workspace);
        const auto central =
            [&](double p, double m) {
                return
                    (p - m) /
                    (2.0 * steps[check]);
            };
        const auto check_gradient =
            [&](double automatic,
                double finite_difference) {
                const double scale =
                    std::max(
                        {1.0,
                         std::abs(automatic),
                         std::abs(
                             finite_difference)});
                require_provider(
                    std::isfinite(automatic) &&
                        std::isfinite(
                            finite_difference) &&
                        std::abs(
                            automatic -
                            finite_difference) <=
                            5.0e-4 * scale,
                    "sourced SW92 provider AD derivative mismatch");
            };

        check_gradient(
            value.viscosity_gradient[0][
                columns[check]],
            central(
                pv.phase_properties[0]
                    .dynamic_viscosity_pa_s,
                mv.phase_properties[0]
                    .dynamic_viscosity_pa_s));
        check_gradient(
            value.enthalpy_gradient[0][
                columns[check]],
            central(
                pv.phase_properties[0]
                    .specific_enthalpy_j_per_kg,
                mv.phase_properties[0]
                    .specific_enthalpy_j_per_kg));
        check_gradient(
            value.internal_energy_gradient[0][
                columns[check]],
            central(
                pv.phase_properties[0]
                    .specific_internal_energy_j_per_kg,
                mv.phase_properties[0]
                    .specific_internal_energy_j_per_kg));
    }

    auto salted = selection;
    salted.nacl_molality_mol_per_kg_water =
        1.0;
    expect_provider_error<std::domain_error>(
        [&] {
            (void)flow::
                make_sw92_co2_water_selected_phase_property_closure(
                    model,
                    std::vector<
                        th::Sw92SelectedPhase<double>>{
                        salted});
        });
}

} // namespace

bool sw92_co2_water_properties_header();

void sw92_co2_water_provider_regression() {
    require_provider(
        sw92_co2_water_properties_header(),
        "SW92 CO2/H2O provider public header check failed");
    sourced_formula_oracles();
    permutation_and_source_guards();
    closure_and_ad_regression();
}
