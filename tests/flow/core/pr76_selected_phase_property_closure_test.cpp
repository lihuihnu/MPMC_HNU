#include <mpmc/flow/pr76_selected_phase_property_closure.hpp>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <iostream>
#include <optional>
#include <span>
#include <stdexcept>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

bool pr76_selected_phase_property_closure_header();

namespace {

namespace th = mpmc::thermodynamics;
namespace flow = mpmc::flow;
namespace ad = mpmc::ad;

void require(bool condition, std::string_view message) {
    if (!condition) {
        throw std::runtime_error(std::string{message});
    }
}

void near(
    double actual,
    double expected,
    double relative = 2.0e-10,
    double absolute = 2.0e-12) {
    const double scale =
        std::max({1.0, std::abs(actual), std::abs(expected)});
    require(
        std::isfinite(actual) &&
            std::isfinite(expected) &&
            std::abs(actual - expected) <=
                absolute + relative * scale,
        "numeric comparison failed");
}

th::Provenance source(std::string locator) {
    return {
        th::SourceKind::synthetic_test,
        "test://pr76-selected-phase-property",
        "v1",
        std::move(locator),
        "Manufactured software-contract fixture; not physical property data",
        "Constructed by the flow property-closure regression",
        "No third-party data copied"};
}

th::SourcedScalar datum(
    double value,
    th::Unit unit,
    std::string locator) {
    return {
        value,
        unit,
        source(std::move(locator)),
        "specified SI",
        "identity"};
}

struct Fixture {
    std::vector<th::Component> catalog;
    th::PrParameterInput input;

    explicit Fixture(bool include_molar_mass = true) {
        constexpr std::array<double, 3> tc{
            400.0, 500.0, 600.0};
        constexpr std::array<double, 3> pc{
            4.0e6, 3.0e6, 5.0e6};
        constexpr std::array<double, 3> omega{
            -0.125, 0.25, 0.5};
        constexpr std::array<double, 3> molar_mass{
            0.016, 0.030, 0.044};
        constexpr double kij[3][3]{
            {0.0, 0.125, -0.0625},
            {0.125, 0.0, 0.0625},
            {-0.0625, 0.0625, 0.0}};

        input.model_id =
            std::string{th::pr76_profile};
        input.dataset_id =
            "manufactured-selected-phase-property";
        input.revision = "v1";
        input.applicability = {
            std::nullopt,
            std::nullopt,
            source("unknown applicability")};

        for (std::size_t i = 0U;
             i < tc.size();
             ++i) {
            const std::string id =
                "test:property:" +
                std::to_string(i);
            std::optional<th::SourcedScalar> mw;
            if (include_molar_mass) {
                mw = datum(
                    molar_mass[i],
                    th::Unit::kilogram_per_mole,
                    "molar-mass-" +
                        std::to_string(i));
            }
            catalog.push_back({
                id,
                "Artificial component",
                th::ComponentKind::pure,
                source(
                    "component-" +
                    std::to_string(i)),
                std::move(mw)});
            input.pure.push_back({
                id,
                datum(
                    tc[i],
                    th::Unit::kelvin,
                    "Tc-" +
                        std::to_string(i)),
                datum(
                    pc[i],
                    th::Unit::pascal,
                    "Pc-" +
                        std::to_string(i)),
                datum(
                    omega[i],
                    th::Unit::dimensionless,
                    "omega-" +
                        std::to_string(i))});
            for (std::size_t j = 0U;
                 j < i;
                 ++j) {
                input.binary.push_back({
                    id,
                    "test:property:" +
                        std::to_string(j),
                    datum(
                        kij[i][j],
                        th::Unit::dimensionless,
                        "kij-" +
                            std::to_string(i) +
                            "-" +
                            std::to_string(j))});
            }
        }
    }

    [[nodiscard]] th::PrParameterSet
    parameters() const {
        const std::array<std::string, 3> order{
            "test:property:0",
            "test:property:1",
            "test:property:2"};
        return th::PrParameterSet::create(
            catalog,
            order,
            input,
            th::DataPolicy::allow_synthetic_tests);
    }
};

struct SyntheticTransportCaloricProvider {
    template <typename Number>
    [[nodiscard]]
    flow::SelectedPhaseTransportCaloricValues<Number>
    operator()(
        const th::Pr76SelectedPhase& selection,
        const Number& pressure_pa,
        const Number& temperature_k,
        std::span<const Number> composition,
        const Number& mixture_molar_mass,
        const Number& molar_density,
        const Number& mass_density) const {
        Number weighted{0.0};
        for (std::size_t component = 0U;
             component < composition.size();
             ++component) {
            weighted +=
                composition[component] *
                static_cast<double>(
                    component + 1U);
        }

        const Number viscosity =
            Number{8.0e-6} +
            pressure_pa * 1.0e-12 +
            temperature_k * 1.0e-8 +
            weighted * 1.0e-6 +
            molar_density * 1.0e-12 +
            Number{
                static_cast<double>(
                    selection.root_index)} *
                1.0e-7;
        const Number enthalpy =
            temperature_k * 1200.0 +
            pressure_pa * 2.0e-4 +
            weighted * 2000.0 +
            mixture_molar_mass * 10000.0 +
            mass_density * 0.1;
        return {viscosity, enthalpy};
    }
};

flow::SelectedPhasePropertyProvenance
provenance() {
    return {
        {
            "PR76 selected density + explicit synthetic molar masses",
            "manufactured-selected-phase-property",
            "v1"},
        {
            "synthetic-test differentiable viscosity",
            "manufactured-selected-phase-property",
            "v1"},
        {
            "synthetic-test differentiable enthalpy",
            "manufactured-selected-phase-property",
            "v1"},
        {
            "thermodynamic identity u=h-p/rho_mass",
            "manufactured-selected-phase-property",
            "v1"}};
}

[[nodiscard]] flow::NaturalVariableLayoutDescriptor
layout(std::size_t phase_count) {
    return {
        3U,
        phase_count,
        std::vector<std::size_t>(
            phase_count,
            2U)};
}

[[nodiscard]] std::vector<double>
natural_variables(
    const flow::NaturalVariableLayoutDescriptor& layout) {
    std::vector<double> q(
        layout.unknown_count(),
        0.0);
    q[layout.pressure_unknown_index()] =
        1.0e6;
    q[layout.temperature_unknown_index()] =
        450.0;

    if (layout.phase_count() >= 2U) {
        const auto s0 =
            layout.independent_saturation_unknown_index(
                flow::PhaseSlot3::phase0);
        require(s0.has_value(), "missing S0");
        q[*s0] =
            layout.phase_count() == 2U
                ? 0.4
                : 0.2;
    }
    if (layout.phase_count() == 3U) {
        const auto s1 =
            layout.independent_saturation_unknown_index(
                flow::PhaseSlot3::phase1);
        require(s1.has_value(), "missing S1");
        q[*s1] = 0.3;
    }

    constexpr std::array<double, 3> composition{
        0.25, 0.50, 0.25};
    for (std::size_t phase = 0U;
         phase < layout.phase_count();
         ++phase) {
        const auto slot =
            static_cast<flow::PhaseSlot3>(phase);
        for (std::size_t component = 0U;
             component < composition.size();
             ++component) {
            const auto column =
                layout.independent_composition_unknown_index(
                    slot,
                    component);
            if (column) {
                q[*column] =
                    composition[component];
            }
        }
    }
    return q;
}

template <typename Closure>
flow::Pr76SelectedPhasePropertyChartLinearization
evaluate(
    const flow::NaturalVariableLayoutDescriptor& chart,
    const std::vector<double>& q,
    const Closure& closure,
    ad::RuntimeJacobianWorkspace<double, 4U>&
        workspace) {
    return flow::
        evaluate_pr76_selected_phase_property_chart(
            chart,
            closure.component_ids(),
            q,
            closure,
            workspace);
}

void cardinality_and_direct_eos() {
    auto parameters =
        Fixture{}.parameters();
    auto model =
        th::Pr76Phase<double>::
            from_parameters(parameters);

    auto closure1 =
        flow::make_pr76_selected_phase_property_closure(
            model,
            std::vector<th::Pr76SelectedPhase>{
                {2U, {}}},
            SyntheticTransportCaloricProvider{},
            provenance());
    auto closure2 =
        flow::make_pr76_selected_phase_property_closure(
            model,
            std::vector<th::Pr76SelectedPhase>{
                {0U, {}},
                {2U, {}}},
            SyntheticTransportCaloricProvider{},
            provenance());
    auto closure3 =
        flow::make_pr76_selected_phase_property_closure(
            model,
            std::vector<th::Pr76SelectedPhase>{
                {0U, {}},
                {1U, {}},
                {2U, {}}},
            SyntheticTransportCaloricProvider{},
            provenance());

    ad::RuntimeJacobianWorkspace<double, 4U>
        workspace;
    const auto layout1 = layout(1U);
    const auto layout2 = layout(2U);
    const auto layout3 = layout(3U);
    const auto q1 =
        natural_variables(layout1);
    const auto q2 =
        natural_variables(layout2);
    const auto q3 =
        natural_variables(layout3);

    const auto result1 =
        evaluate(
            layout1,
            q1,
            closure1,
            workspace);
    const auto result2 =
        evaluate(
            layout2,
            q2,
            closure2,
            workspace);
    const auto result3 =
        evaluate(
            layout3,
            q3,
            closure3,
            workspace);

    require(
        result1.phase_properties.size() == 1U &&
            result1.equilibrium_residual.empty(),
        "1P property chart shape changed");
    require(
        result2.phase_properties.size() == 2U &&
            result2.equilibrium_residual.size() ==
                3U,
        "2P property chart shape changed");
    require(
        result3.phase_properties.size() == 3U &&
            result3.equilibrium_residual.size() ==
                6U,
        "3P property chart shape changed");

    constexpr std::array<double, 3> x{
        0.25, 0.50, 0.25};
    const auto point =
        closure3.evaluate(
            1U,
            1.0e6,
            450.0,
            std::span<const double>{x});

    th::Pr76PhaseWorkspace<double>
        density_workspace;
    const auto direct_density =
        th::evaluate_selected_phase_molar_density(
            model,
            1.0e6,
            450.0,
            std::span<const double>{x},
            th::Pr76SelectedPhase{1U, {}},
            density_workspace);
    th::Pr76PhaseWorkspace<double>
        fugacity_workspace;
    const auto direct_fugacity =
        th::evaluate_selected_phase_fugacity(
            model,
            1.0e6,
            450.0,
            std::span<const double>{x},
            th::Pr76SelectedPhase{1U, {}},
            fugacity_workspace);

    near(
        point.molar_density_mol_per_m3,
        direct_density.molar_density_mol_per_m3);
    constexpr double mixture_molar_mass =
        0.25 * 0.016 +
        0.50 * 0.030 +
        0.25 * 0.044;
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
            1.0e6 /
                point.mass_density_kg_per_m3);
    require(
        point.ln_phi.size() ==
            direct_fugacity.ln_phi.size(),
        "direct fugacity shape changed");
    for (std::size_t component = 0U;
         component < point.ln_phi.size();
         ++component) {
        near(
            point.ln_phi[component],
            direct_fugacity.ln_phi[component]);
    }

    const auto s0 =
        layout3.independent_saturation_unknown_index(
            flow::PhaseSlot3::phase0);
    const auto s1 =
        layout3.independent_saturation_unknown_index(
            flow::PhaseSlot3::phase1);
    require(
        s0.has_value() && s1.has_value(),
        "3P saturation columns missing");
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
                "pc=none phase property acquired saturation derivative");
        }
    }
    for (std::size_t row = 0U;
         row < result3.equilibrium_residual.size();
         ++row) {
        require(
            result3.equilibrium_jacobian[
                row * layout3.unknown_count() +
                *s0] == 0.0 &&
                result3.equilibrium_jacobian[
                    row * layout3.unknown_count() +
                    *s1] == 0.0,
            "pc=none equilibrium row acquired saturation derivative");
    }

    const auto bridged =
        flow::
            make_pr76_selected_phase_flow_linearization_3p(
                result3,
                q3);
    for (std::size_t phase = 0U;
         phase < 3U;
         ++phase) {
        const auto slot =
            static_cast<flow::PhaseSlot3>(phase);
        near(
            bridged.state
                .phase_properties(slot)
                .mass_density_kg_per_m3,
            result3
                .phase_properties[phase]
                .mass_density_kg_per_m3);
        near(
            bridged.molar_density
                .molar_density_mol_per_m3[phase],
            result3
                .phase_properties[phase]
                .molar_density_mol_per_m3);
    }
    require(
        bridged.fugacity.residual_count() ==
            6U,
        "3P fugacity bridge row count changed");
}

void derivative_oracle() {
    auto parameters =
        Fixture{}.parameters();
    auto model =
        th::Pr76Phase<double>::
            from_parameters(parameters);
    auto closure =
        flow::make_pr76_selected_phase_property_closure(
            model,
            std::vector<th::Pr76SelectedPhase>{
                {0U, {}},
                {1U, {}},
                {2U, {}}},
            SyntheticTransportCaloricProvider{},
            provenance());

    const auto chart = layout(3U);
    auto q = natural_variables(chart);
    ad::RuntimeJacobianWorkspace<double, 4U>
        workspace;
    const auto derived =
        evaluate(
            chart,
            q,
            closure,
            workspace);

    const auto composition_column =
        chart.independent_composition_unknown_index(
            flow::PhaseSlot3::phase1,
            0U);
    require(
        composition_column.has_value(),
        "composition perturbation column missing");

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
        plus[columns[check]] +=
            steps[check];
        minus[columns[check]] -=
            steps[check];

        const auto plus_values =
            evaluate(
                chart,
                plus,
                closure,
                workspace);
        const auto minus_values =
            evaluate(
                chart,
                minus,
                closure,
                workspace);

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
                         std::abs(
                             finite_difference)});
                require(
                    std::isfinite(automatic) &&
                        std::isfinite(
                            finite_difference) &&
                        std::abs(
                            automatic -
                            finite_difference) <=
                            2.0e-4 * scale,
                    "PR76 property AD derivative disagrees with fresh central perturbation");
            };

        const std::size_t phase = 1U;
        compare(
            derived.mass_density_gradient[
                phase][columns[check]],
            central(
                plus_values.phase_properties[
                    phase]
                    .mass_density_kg_per_m3,
                minus_values.phase_properties[
                    phase]
                    .mass_density_kg_per_m3));
        compare(
            derived.viscosity_gradient[
                phase][columns[check]],
            central(
                plus_values.phase_properties[
                    phase]
                    .dynamic_viscosity_pa_s,
                minus_values.phase_properties[
                    phase]
                    .dynamic_viscosity_pa_s));
        compare(
            derived.enthalpy_gradient[
                phase][columns[check]],
            central(
                plus_values.phase_properties[
                    phase]
                    .specific_enthalpy_j_per_kg,
                minus_values.phase_properties[
                    phase]
                    .specific_enthalpy_j_per_kg));
        compare(
            derived.internal_energy_gradient[
                phase][columns[check]],
            central(
                plus_values.phase_properties[
                    phase]
                    .specific_internal_energy_j_per_kg,
                minus_values.phase_properties[
                    phase]
                    .specific_internal_energy_j_per_kg));

        const std::size_t component = 0U;
        compare(
            derived.ln_phi_jacobian[phase][
                component *
                    chart.unknown_count() +
                columns[check]],
            central(
                plus_values.ln_phi[phase][component],
                minus_values.ln_phi[phase][component]));

        const std::size_t equilibrium_row = 0U;
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

void invalid_inputs() {
    auto parameters =
        Fixture{}.parameters();
    auto model =
        th::Pr76Phase<double>::
            from_parameters(parameters);
    auto closure =
        flow::make_pr76_selected_phase_property_closure(
            model,
            std::vector<th::Pr76SelectedPhase>{
                {0U, {}}},
            SyntheticTransportCaloricProvider{},
            provenance());

    ad::RuntimeJacobianWorkspace<double, 4U>
        workspace;
    const auto chart2 = layout(2U);
    const auto q2 =
        natural_variables(chart2);

    bool mismatch = false;
    try {
        (void)evaluate(
            chart2,
            q2,
            closure,
            workspace);
    } catch (const std::invalid_argument&) {
        mismatch = true;
    }
    require(
        mismatch,
        "phase-cardinality mismatch was accepted");

    bool missing_molar_mass = false;
    try {
        auto missing_parameters =
            Fixture{false}.parameters();
        auto missing_model =
            th::Pr76Phase<double>::
                from_parameters(
                    missing_parameters);
        (void)flow::
            make_pr76_selected_phase_property_closure(
                missing_model,
                std::vector<
                    th::Pr76SelectedPhase>{
                    {0U, {}}},
                SyntheticTransportCaloricProvider{},
                provenance());
    } catch (const std::invalid_argument&) {
        missing_molar_mass = true;
    }
    require(
        missing_molar_mass,
        "missing component molar mass was accepted");
}

} // namespace

int main() {
    try {
        require(
            pr76_selected_phase_property_closure_header(),
            "public header probe failed");
        cardinality_and_direct_eos();
        derivative_oracle();
        invalid_inputs();
        std::cout
            << "[PASS] PR76 selected-phase property closure"
            << std::endl;
        return 0;
    } catch (const std::exception& error) {
        std::cerr
            << "[FAIL] "
            << error.what()
            << std::endl;
        return 1;
    }
}
