#include <mpmc/flow/sw92_selected_phase_property_closure.hpp>

#include "test_support.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <iostream>
#include <limits>
#include <span>
#include <stdexcept>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

bool sw92_selected_phase_property_closure_header();
bool sw92_co2_water_properties_header();
void sw92_co2_water_provider_regression();

namespace {

namespace th = mpmc::thermodynamics;
namespace flow = mpmc::flow;
namespace ad = mpmc::ad;
namespace st = sw92_test;

void require(bool condition, std::string_view message) {
    if (!condition) {
        throw std::runtime_error(std::string{message});
    }
}

void near(
    double actual,
    double expected,
    double relative = 3.0e-10,
    double absolute = 3.0e-12) {
    const double scale =
        std::max({1.0, std::abs(actual), std::abs(expected)});
    require(
        std::isfinite(actual) &&
            std::isfinite(expected) &&
            std::abs(actual - expected) <=
                absolute + relative * scale,
        "numeric comparison failed");
}

template <class Error, class Function>
void expect_error(Function&& function) {
    bool caught = false;
    try {
        std::forward<Function>(function)();
    } catch (const Error&) {
        caught = true;
    }
    require(caught, "expected exception missing");
}

template <class Function>
void expect_pc_none_capability_error(
    Function&& function,
    std::string_view diagnostic_fragment) {
    try {
        std::forward<Function>(function)();
    } catch (const flow::
                 Sw92SelectedPhasePcNoneCapabilityError&
                     error) {
        require(
            std::string_view{error.what()}.find(
                diagnostic_fragment) !=
                std::string_view::npos,
            "SW92 pc=none capability diagnostic changed");
        return;
    }
    throw std::runtime_error(
        "expected SW92 pc=none capability rejection");
}

st::Prepared prepared_with_molar_masses(
    bool reverse = false) {
    auto prepared =
        st::binary_input(st::co2, reverse);
    const auto source =
        st::synthetic(
            "flow selected-phase mass-density molar masses");
    for (auto& component : prepared.catalog) {
        const double molar_mass =
            component.id == st::co2.id
                ? 0.0440095
                : 0.01801528;
        component.molar_mass =
            st::scalar(
                molar_mass,
                th::Unit::kilogram_per_mole,
                source,
                "kg/mol",
                "identity");
    }
    return prepared;
}

th::Sw92ParameterSet parameters(
    bool reverse = false) {
    auto prepared =
        prepared_with_molar_masses(reverse);
    return th::Sw92ParameterSet::create(
        prepared.catalog,
        prepared.order,
        prepared.input,
        th::DataPolicy::allow_synthetic_tests);
}

struct SyntheticTransportCaloricProvider {
    template <typename Number>
    [[nodiscard]]
    flow::SelectedPhaseTransportCaloricValues<Number>
    operator()(
        const th::Sw92SelectedPhase<double>& selection,
        const Number& pressure_pa,
        const Number& temperature_k,
        std::span<const Number> composition,
        const Number& mixture_molar_mass,
        const Number& molar_density,
        const Number& mass_density) const {
        Number composition_metric{0.0};
        for (const auto& value : composition) {
            composition_metric += value * value;
        }
        const double family_term =
            selection.family ==
                    th::SwPhaseFamily::aqueous
                ? 1.0
                : 2.0;
        const double selection_term =
            family_term +
            static_cast<double>(
                selection.root_index) +
            0.25 *
                selection
                    .nacl_molality_mol_per_kg_water;

        const Number viscosity =
            Number{8.0e-6 +
                   selection_term * 1.0e-7} +
            pressure_pa * 1.0e-12 +
            temperature_k * 1.0e-8 +
            composition_metric * 1.0e-6 +
            molar_density * 1.0e-12;
        const Number enthalpy =
            Number{selection_term * 100.0} +
            temperature_k * 1200.0 +
            pressure_pa * 2.0e-4 +
            composition_metric * 2000.0 +
            mixture_molar_mass * 10000.0 +
            mass_density * 0.1;
        return {viscosity, enthalpy};
    }
};

flow::SelectedPhasePropertyProvenance
provenance() {
    return {
        {
            "SW92 selected density + explicit synthetic molar masses",
            "manufactured-sw92-selected-phase-property",
            "v1"},
        {
            "synthetic-test scalar-generic viscosity",
            "manufactured-sw92-selected-phase-property",
            "v1"},
        {
            "synthetic-test scalar-generic enthalpy",
            "manufactured-sw92-selected-phase-property",
            "v1"},
        {
            "thermodynamic identity u=h-p/rho_mass",
            "manufactured-sw92-selected-phase-property",
            "v1"}};
}

[[nodiscard]]
th::Sw92SelectedPhase<double>
selected(std::size_t root) {
    return {
        0.0,
        th::SwPhaseFamily::nonaqueous,
        root,
        {}};
}

[[nodiscard]]
flow::NaturalVariableLayoutDescriptor
layout(std::size_t phase_count) {
    return {
        2U,
        phase_count,
        std::vector<std::size_t>(
            phase_count,
            1U)};
}

[[nodiscard]] std::vector<double>
natural_variables(
    const flow::NaturalVariableLayoutDescriptor& chart,
    bool reverse = false) {
    std::vector<double> q(
        chart.unknown_count(),
        0.0);
    q[chart.pressure_unknown_index()] =
        3.0e6;
    q[chart.temperature_unknown_index()] =
        340.0;

    if (chart.phase_count() >= 2U) {
        const auto s0 =
            chart.independent_saturation_unknown_index(
                flow::PhaseSlot3::phase0);
        require(s0.has_value(), "missing S0");
        q[*s0] =
            chart.phase_count() == 2U
                ? 0.4
                : 0.2;
    }
    if (chart.phase_count() == 3U) {
        const auto s1 =
            chart.independent_saturation_unknown_index(
                flow::PhaseSlot3::phase1);
        require(s1.has_value(), "missing S1");
        q[*s1] = 0.3;
    }

    const double independent =
        reverse ? 0.3 : 0.7;
    for (std::size_t phase = 0U;
         phase < chart.phase_count();
         ++phase) {
        const auto column =
            chart.independent_composition_unknown_index(
                static_cast<flow::PhaseSlot3>(phase),
                0U);
        require(
            column.has_value(),
            "missing independent composition column");
        q[*column] = independent;
    }
    return q;
}

template <typename Closure>
flow::Sw92SelectedPhasePropertyChartLinearization
evaluate(
    const flow::NaturalVariableLayoutDescriptor& chart,
    const std::vector<double>& q,
    const Closure& closure,
    ad::RuntimeJacobianWorkspace<double, 4U>&
        workspace) {
    return flow::
        evaluate_sw92_selected_phase_property_chart(
            chart,
            closure.component_ids(),
            q,
            closure,
            workspace);
}

struct LinearSaturationRelativePermeability3P {
    template <typename Number>
    [[nodiscard]]
    flow::RelativePermeabilityEvaluation3P<Number>
    operator()(
        const flow::
            ThreePhaseSaturationState3P<Number>&
                state) const {
        return {state.saturation};
    }
};

struct ConstantNonzeroCapillaryPressure3P {
    template <typename Number>
    [[nodiscard]]
    flow::CapillaryPressureOffsetsEvaluation3P<Number>
    operator()(
        const flow::
            ThreePhaseSaturationState3P<Number>&) const {
        return {
            std::array<Number, 2>{
                Number{1000.0},
                Number{-2000.0}}};
    }
};

struct ZeroValueNonzeroCapillaryJacobian3P {
    template <typename Number>
    [[nodiscard]]
    flow::CapillaryPressureOffsetsEvaluation3P<Number>
    operator()(
        const flow::
            ThreePhaseSaturationState3P<Number>&
                state) const {
        return {
            std::array<Number, 2>{
                Number{5.0e4} *
                    (state.saturation[0] -
                     Number{0.20}),
                Number{7.0e4} *
                    (state.saturation[1] -
                     Number{0.30})}};
    }
};

template <typename CapillaryPressureEvaluator>
[[nodiscard]]
flow::
    ThreePhaseSaturationConstitutiveNaturalVariableLinearization3P
make_test_saturation_linearization_from_ad(
    const flow::NaturalVariableCellState3P& state,
    CapillaryPressureEvaluator capillary_pressure) {
    using D = ad::Dual<double, 2U>;

    const double saturation0 =
        state.phase_saturation(
            flow::PhaseSlot3::phase0);
    const double saturation1 =
        state.phase_saturation(
            flow::PhaseSlot3::phase1);

    const auto primal =
        flow::
            evaluate_three_phase_saturation_constitutive(
                state,
                LinearSaturationRelativePermeability3P{},
                capillary_pressure);
    const auto differentiated =
        flow::
            evaluate_three_phase_saturation_constitutive(
                D{state.reference_pressure_pa()},
                D::variable(saturation0, 0U),
                D::variable(saturation1, 1U),
                LinearSaturationRelativePermeability3P{},
                capillary_pressure);

    flow::
        ThreePhaseSaturationCoordinateDerivatives3P
            derivatives{};
    for (std::size_t phase = 0U;
         phase < 3U;
         ++phase) {
        for (std::size_t direction = 0U;
             direction < 2U;
             ++direction) {
            derivatives
                .relative_permeability[phase][direction] =
                differentiated
                    .relative_permeability[phase]
                    .derivative(direction);
            derivatives
                .capillary_pressure_offset_pa[phase][direction] =
                differentiated
                    .capillary_pressure_offset_pa[phase]
                    .derivative(direction);
        }
    }

    return flow::
        make_saturation_constitutive_natural_variable_linearization(
            state,
            primal,
            derivatives);
}

void cardinality_and_direct_eos() {
    auto ps = parameters();
    auto model =
        th::Sw92Phase<double>::from_parameters(ps);

    auto closure1 =
        flow::make_sw92_selected_phase_property_closure(
            model,
            std::vector<th::Sw92SelectedPhase<double>>{
                selected(2U)},
            SyntheticTransportCaloricProvider{},
            provenance());
    auto closure2 =
        flow::make_sw92_selected_phase_property_closure(
            model,
            std::vector<th::Sw92SelectedPhase<double>>{
                selected(0U),
                selected(2U)},
            SyntheticTransportCaloricProvider{},
            provenance());
    auto closure3 =
        flow::make_sw92_selected_phase_property_closure(
            model,
            std::vector<th::Sw92SelectedPhase<double>>{
                selected(0U),
                selected(1U),
                selected(2U)},
            SyntheticTransportCaloricProvider{},
            provenance());

    ad::RuntimeJacobianWorkspace<double, 4U>
        workspace;
    const auto layout1 = layout(1U);
    const auto layout2 = layout(2U);
    const auto layout3 = layout(3U);
    const auto q1 = natural_variables(layout1);
    const auto q2 = natural_variables(layout2);
    const auto q3 = natural_variables(layout3);

    const auto result1 =
        evaluate(layout1, q1, closure1, workspace);
    const auto result2 =
        evaluate(layout2, q2, closure2, workspace);
    const auto result3 =
        evaluate(layout3, q3, closure3, workspace);

    require(
        result1.phase_properties.size() == 1U &&
            result1.equilibrium_residual.empty(),
        "SW92 1P property chart shape changed");
    require(
        result2.phase_properties.size() == 2U &&
            result2.equilibrium_residual.size() == 2U,
        "SW92 2P property chart shape changed");
    require(
        result3.phase_properties.size() == 3U &&
            result3.equilibrium_residual.size() == 4U,
        "SW92 3P property chart shape changed");

    constexpr std::array<double, 2> x{
        0.7, 0.3};
    const auto point =
        closure3.evaluate(
            1U,
            3.0e6,
            340.0,
            std::span<const double>{x});

    th::Sw92PhaseWorkspace<double>
        density_workspace;
    const auto direct_density =
        th::evaluate_selected_phase_molar_density(
            model,
            3.0e6,
            340.0,
            std::span<const double>{x},
            selected(1U),
            density_workspace);
    th::Sw92PhaseWorkspace<double>
        fugacity_workspace;
    const auto direct_fugacity =
        th::evaluate_selected_phase_fugacity(
            model,
            3.0e6,
            340.0,
            std::span<const double>{x},
            selected(1U),
            fugacity_workspace);

    near(
        point.molar_density_mol_per_m3,
        direct_density.molar_density_mol_per_m3);
    const double mixture_molar_mass =
        0.7 * 0.0440095 +
        0.3 * 0.01801528;
    near(
        point.mixture_molar_mass_kg_per_mol,
        mixture_molar_mass);
    near(
        point.mass_density_kg_per_m3,
        point.molar_density_mol_per_m3 *
            mixture_molar_mass);
    near(
        point.specific_internal_energy_j_per_kg,
        point.specific_enthalpy_j_per_kg -
            3.0e6 /
                point.mass_density_kg_per_m3);
    require(
        point.ln_phi.size() ==
            direct_fugacity.ln_phi.size(),
        "SW92 direct fugacity shape changed");
    for (std::size_t component = 0U;
         component < point.ln_phi.size();
         ++component) {
        near(
            point.ln_phi[component],
            direct_fugacity.ln_phi[component]);
    }

    const auto& frozen =
        closure3.selection(1U);
    require(
        frozen.family ==
                th::SwPhaseFamily::nonaqueous &&
            frozen.root_index == 1U &&
            frozen.nacl_molality_mol_per_kg_water ==
                0.0,
        "SW92 family/root/molality provenance changed");

    const auto bridge1 =
        flow::
            make_sw92_selected_phase_flow_linearization_1p(
                result1,
                q1,
                1.0,
                std::vector<double>(
                    layout1.unknown_count(),
                    0.0));
    const auto bridge2 =
        flow::
            make_sw92_selected_phase_flow_linearization_2p(
                result2,
                q2,
                std::array<double, 2>{0.4, 0.6},
                std::array<std::vector<double>, 2>{
                    std::vector<double>(
                        layout2.unknown_count(),
                        0.0),
                    std::vector<double>(
                        layout2.unknown_count(),
                        0.0)});
    const auto bridge3 =
        flow::
            make_sw92_selected_phase_flow_linearization_3p(
                result3,
                q3);

    require(
        bridge1.state.component_ids().size() == 2U &&
            bridge2.state.component_ids().size() == 2U &&
            bridge3.state.component_ids().size() == 2U,
        "SW92 natural-variable bridge component count changed");
    require(
        bridge2.fugacity.residual_count() == 2U &&
            bridge3.fugacity.residual_count() == 4U,
        "SW92 natural-variable fugacity row count changed");

    const auto s0 =
        layout3.independent_saturation_unknown_index(
            flow::PhaseSlot3::phase0);
    const auto s1 =
        layout3.independent_saturation_unknown_index(
            flow::PhaseSlot3::phase1);
    require(
        s0.has_value() && s1.has_value(),
        "SW92 3P saturation columns missing");
    for (std::size_t phase = 0U;
         phase < 3U;
         ++phase) {
        for (const auto* gradient :
             std::array<const std::vector<double>*, 5>{
                 &result3.molar_density_gradient[phase],
                 &result3.mass_density_gradient[phase],
                 &result3.viscosity_gradient[phase],
                 &result3.enthalpy_gradient[phase],
                 &result3.internal_energy_gradient[phase]}) {
            require(
                (*gradient)[*s0] == 0.0 &&
                    (*gradient)[*s1] == 0.0,
                "SW92 pc=none property acquired saturation derivative");
        }
    }
}

void derivative_oracle() {
    auto ps = parameters();
    auto model =
        th::Sw92Phase<double>::from_parameters(ps);
    auto closure =
        flow::make_sw92_selected_phase_property_closure(
            model,
            std::vector<th::Sw92SelectedPhase<double>>{
                selected(0U),
                selected(1U),
                selected(2U)},
            SyntheticTransportCaloricProvider{},
            provenance());

    const auto chart = layout(3U);
    const auto q = natural_variables(chart);
    ad::RuntimeJacobianWorkspace<double, 4U>
        workspace;
    const auto derived =
        evaluate(chart, q, closure, workspace);

    const auto composition_column =
        chart.independent_composition_unknown_index(
            flow::PhaseSlot3::phase1,
            0U);
    require(
        composition_column.has_value(),
        "SW92 composition perturbation column missing");

    const std::array<std::size_t, 3> columns{
        chart.pressure_unknown_index(),
        chart.temperature_unknown_index(),
        *composition_column};
    const std::array<double, 3> steps{
        10.0,
        1.0e-4,
        1.0e-6};

    for (std::size_t check = 0U;
         check < columns.size();
         ++check) {
        auto plus = q;
        auto minus = q;
        plus[columns[check]] += steps[check];
        minus[columns[check]] -= steps[check];

        const auto plus_values =
            evaluate(chart, plus, closure, workspace);
        const auto minus_values =
            evaluate(chart, minus, closure, workspace);

        const auto central =
            [&](double plus_value,
                double minus_value) {
                return
                    (plus_value - minus_value) /
                    (2.0 * steps[check]);
            };
        const auto compare =
            [&](double automatic,
                double finite_difference) {
                const double scale =
                    std::max(
                        {1.0,
                         std::abs(automatic),
                         std::abs(finite_difference)});
                require(
                    std::isfinite(automatic) &&
                        std::isfinite(finite_difference) &&
                        std::abs(
                            automatic -
                            finite_difference) <=
                            4.0e-4 * scale,
                    "SW92 property AD derivative disagrees with fresh central perturbation");
            };

        constexpr std::size_t phase = 1U;
        compare(
            derived.molar_density_gradient[
                phase][columns[check]],
            central(
                plus_values.phase_properties[phase]
                    .molar_density_mol_per_m3,
                minus_values.phase_properties[phase]
                    .molar_density_mol_per_m3));
        compare(
            derived.mass_density_gradient[
                phase][columns[check]],
            central(
                plus_values.phase_properties[phase]
                    .mass_density_kg_per_m3,
                minus_values.phase_properties[phase]
                    .mass_density_kg_per_m3));
        compare(
            derived.viscosity_gradient[
                phase][columns[check]],
            central(
                plus_values.phase_properties[phase]
                    .dynamic_viscosity_pa_s,
                minus_values.phase_properties[phase]
                    .dynamic_viscosity_pa_s));
        compare(
            derived.enthalpy_gradient[
                phase][columns[check]],
            central(
                plus_values.phase_properties[phase]
                    .specific_enthalpy_j_per_kg,
                minus_values.phase_properties[phase]
                    .specific_enthalpy_j_per_kg));
        compare(
            derived.internal_energy_gradient[
                phase][columns[check]],
            central(
                plus_values.phase_properties[phase]
                    .specific_internal_energy_j_per_kg,
                minus_values.phase_properties[phase]
                    .specific_internal_energy_j_per_kg));

        constexpr std::size_t component = 0U;
        compare(
            derived.ln_phi_jacobian[phase][
                component *
                    chart.unknown_count() +
                columns[check]],
            central(
                plus_values.ln_phi[phase][component],
                minus_values.ln_phi[phase][component]));

        constexpr std::size_t equilibrium_row = 0U;
        compare(
            derived.equilibrium_jacobian[
                equilibrium_row *
                    chart.unknown_count() +
                columns[check]],
            central(
                plus_values.equilibrium_residual[
                    equilibrium_row],
                minus_values.equilibrium_residual[
                    equilibrium_row]));
    }
}

void permutation_and_fixed_selection() {
    auto forward_parameters = parameters(false);
    auto reverse_parameters = parameters(true);
    auto forward_model =
        th::Sw92Phase<double>::from_parameters(
            forward_parameters);
    auto reverse_model =
        th::Sw92Phase<double>::from_parameters(
            reverse_parameters);

    auto forward =
        flow::make_sw92_selected_phase_property_closure(
            forward_model,
            std::vector<th::Sw92SelectedPhase<double>>{
                selected(2U)},
            SyntheticTransportCaloricProvider{},
            provenance());
    auto reverse =
        flow::make_sw92_selected_phase_property_closure(
            reverse_model,
            std::vector<th::Sw92SelectedPhase<double>>{
                selected(2U)},
            SyntheticTransportCaloricProvider{},
            provenance());

    constexpr std::array<double, 2> xf{
        0.7, 0.3};
    constexpr std::array<double, 2> xr{
        0.3, 0.7};
    const auto a =
        forward.evaluate(
            0U,
            3.0e6,
            340.0,
            std::span<const double>{xf});
    const auto b =
        reverse.evaluate(
            0U,
            3.0e6,
            340.0,
            std::span<const double>{xr});

    near(
        a.molar_density_mol_per_m3,
        b.molar_density_mol_per_m3,
        2.0e-11);
    near(
        a.mass_density_kg_per_m3,
        b.mass_density_kg_per_m3,
        2.0e-11);
    near(
        a.dynamic_viscosity_pa_s,
        b.dynamic_viscosity_pa_s,
        2.0e-11);
    near(
        a.specific_enthalpy_j_per_kg,
        b.specific_enthalpy_j_per_kg,
        2.0e-11);
    near(a.ln_phi[0], b.ln_phi[1], 2.0e-11);
    near(a.ln_phi[1], b.ln_phi[0], 2.0e-11);

    require(
        forward.component_ids()[0] ==
                st::co2.id &&
            reverse.component_ids()[0] ==
                st::water.id,
        "SW92 ordered component identity did not follow snapshot permutation");
    require(
        forward.selection(0U).family ==
                reverse.selection(0U).family &&
            forward.selection(0U).root_index ==
                reverse.selection(0U).root_index &&
            forward.selection(0U)
                    .nacl_molality_mol_per_kg_water ==
                reverse.selection(0U)
                    .nacl_molality_mol_per_kg_water,
        "SW92 frozen family/root/molality changed under component permutation");
}

void pc_none_capability_guard() {
    auto ps = parameters();
    auto model =
        th::Sw92Phase<double>::from_parameters(ps);
    auto closure =
        flow::make_sw92_selected_phase_property_closure(
            model,
            std::vector<th::Sw92SelectedPhase<double>>{
                selected(0U),
                selected(1U),
                selected(2U)},
            SyntheticTransportCaloricProvider{},
            provenance());

    const auto chart = layout(3U);
    const auto q = natural_variables(chart);
    ad::RuntimeJacobianWorkspace<double, 4U>
        workspace;
    const auto properties =
        evaluate(chart, q, closure, workspace);
    const auto flow_linearization =
        flow::
            make_sw92_selected_phase_flow_linearization_3p(
                properties,
                q);

    const auto valid =
        make_test_saturation_linearization_from_ad(
            flow_linearization.state,
            flow::NoCapillaryPressure3P{});
    flow::
        require_sw92_selected_phase_pc_none_capability(
            flow_linearization.state,
            valid);

    const auto nonzero_value =
        make_test_saturation_linearization_from_ad(
            flow_linearization.state,
            ConstantNonzeroCapillaryPressure3P{});
    expect_pc_none_capability_error(
        [&] {
            flow::
                require_sw92_selected_phase_pc_none_capability(
                    flow_linearization.state,
                    nonzero_value);
        },
        "nonzero capillary/phase-pressure offset");

    const auto zero_value_nonzero_jacobian =
        make_test_saturation_linearization_from_ad(
            flow_linearization.state,
            ZeroValueNonzeroCapillaryJacobian3P{});
    require(
        zero_value_nonzero_jacobian
                .capillary_pressure_offset_pa[1] ==
            0.0 &&
            zero_value_nonzero_jacobian
                .capillary_pressure_offset_pa[2] ==
            0.0,
        "SW92 zero-value/nonzero-Jacobian capillary fixture changed");
    expect_pc_none_capability_error(
        [&] {
            flow::
                require_sw92_selected_phase_pc_none_capability(
                    flow_linearization.state,
                    zero_value_nonzero_jacobian);
        },
        "nonzero capillary/phase-pressure Jacobian");
}

void invalid_inputs_and_branch_failures() {
    auto missing = st::binary_input(st::co2);
    auto missing_parameters =
        th::Sw92ParameterSet::create(
            missing.catalog,
            missing.order,
            missing.input);
    auto missing_model =
        th::Sw92Phase<double>::from_parameters(
            missing_parameters);
    expect_error<std::invalid_argument>(
        [&] {
            (void)flow::
                make_sw92_selected_phase_property_closure(
                    missing_model,
                    std::vector<
                        th::Sw92SelectedPhase<double>>{
                        selected(2U)},
                    SyntheticTransportCaloricProvider{},
                    provenance());
        });

    auto ps = parameters();
    auto model =
        th::Sw92Phase<double>::from_parameters(ps);

    auto bad_molality = selected(2U);
    bad_molality.nacl_molality_mol_per_kg_water =
        std::numeric_limits<double>::quiet_NaN();
    expect_error<std::invalid_argument>(
        [&] {
            (void)flow::
                make_sw92_selected_phase_property_closure(
                    model,
                    std::vector<
                        th::Sw92SelectedPhase<double>>{
                        bad_molality},
                    SyntheticTransportCaloricProvider{},
                    provenance());
        });

    auto bad_family = selected(2U);
    bad_family.family =
        static_cast<th::SwPhaseFamily>(99);
    expect_error<std::invalid_argument>(
        [&] {
            (void)flow::
                make_sw92_selected_phase_property_closure(
                    model,
                    std::vector<
                        th::Sw92SelectedPhase<double>>{
                        bad_family},
                    SyntheticTransportCaloricProvider{},
                    provenance());
        });

    auto invalid_root =
        flow::make_sw92_selected_phase_property_closure(
            model,
            std::vector<th::Sw92SelectedPhase<double>>{
                selected(99U)},
            SyntheticTransportCaloricProvider{},
            provenance());
    constexpr std::array<double, 2> x{
        0.7, 0.3};
    expect_error<std::out_of_range>(
        [&] {
            (void)invalid_root.evaluate(
                0U,
                3.0e6,
                340.0,
                std::span<const double>{x});
        });

    auto limited_selection = selected(2U);
    limited_selection.root_options.max_iterations = 1;
    auto limited =
        flow::make_sw92_selected_phase_property_closure(
            model,
            std::vector<th::Sw92SelectedPhase<double>>{
                limited_selection},
            SyntheticTransportCaloricProvider{},
            provenance());
    expect_error<th::Sw92PhaseError>(
        [&] {
            (void)limited.evaluate(
                0U,
                3.0e6,
                340.0,
                std::span<const double>{x});
        });
}

} // namespace

int main() {
    try {
        require(
            sw92_selected_phase_property_closure_header(),
            "SW92 property closure public header check failed");
        cardinality_and_direct_eos();
        derivative_oracle();
        permutation_and_fixed_selection();
        pc_none_capability_guard();
        invalid_inputs_and_branch_failures();
        sw92_co2_water_provider_regression();
        std::cout
            << "PASS SW92 selected-phase property -> flow natural-variable bridge\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr
            << "FAIL SW92 selected-phase property bridge: "
            << error.what() << '\n';
        return 1;
    }
}
