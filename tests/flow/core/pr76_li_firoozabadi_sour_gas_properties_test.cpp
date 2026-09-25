#include <mpmc/flow/pr76_li_firoozabadi_sour_gas_properties.hpp>

#include "../../flash/pr76_three_phase/sour_gas_fixture.hpp"
#include "../../flash/pr76_three_phase/sour_gas_references.hpp"
#include "pr76_li_firoozabadi_sour_gas_flow_reference.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <limits>
#include <span>
#include <stdexcept>
#include <string_view>
#include <vector>

namespace {

namespace ad = mpmc::ad;
namespace flow = mpmc::flow;
namespace th = mpmc::thermodynamics;
namespace fx = pr76_sour_gas_test;
namespace eq_ref = pr76_sour_gas_reference;
namespace flow_ref =
    pr76_li_firoozabadi_sour_gas_flow_reference;

void require(
    bool condition,
    std::string_view message) {
    if (!condition) {
        throw std::runtime_error(
            std::string{message});
    }
}

void near(
    double actual,
    double expected,
    double relative = 3.0e-10,
    double absolute = 3.0e-11) {
    const double scale =
        std::max(
            {1.0,
             std::abs(actual),
             std::abs(expected)});
    require(
        std::isfinite(actual) &&
            std::isfinite(expected) &&
            std::abs(actual - expected) <=
                absolute + relative * scale,
        "Li-Firoozabadi sour-gas provider numeric comparison failed");
}

std::vector<th::Pr76SelectedPhase>
reference_selections(
    const th::Pr76Phase<double>& model) {
    std::vector<th::Pr76SelectedPhase>
        selections;
    selections.reserve(3U);

    for (std::size_t phase = 0U;
         phase < 3U;
         ++phase) {
        const std::vector<double>
            composition{
                eq_ref::phases[phase].begin(),
                eq_ref::phases[phase].end()};
        th::Pr76PhaseWorkspace<double>
            workspace;
        const auto roots =
            model.roots_full(
                eq_ref::pressure_pa,
                eq_ref::temperature_k,
                composition,
                workspace);
        require(
            roots.status ==
                    th::Pr76RootStatus::success &&
                roots.count > 0U,
            "Li-Firoozabadi reference phase has no resolved PR76 root");

        std::size_t best = 0U;
        double best_error =
            std::numeric_limits<double>::
                infinity();
        for (std::size_t root = 0U;
             root < roots.count;
             ++root) {
            const double error =
                std::abs(
                    roots.roots[root].z -
                    eq_ref::
                        compressibility_factors[
                            phase]);
            if (error < best_error) {
                best_error = error;
                best = root;
            }
        }
        require(
            best_error <= 5.0e-12,
            "Li-Firoozabadi reference Z did not map to one PR76 selected root");
        selections.push_back(
            {best, {}});
    }

    return selections;
}

flow::NaturalVariableLayout3P
reference_layout() {
    const std::array<std::vector<double>, 3>
        compositions{
            std::vector<double>{
                eq_ref::phases[0].begin(),
                eq_ref::phases[0].end()},
            std::vector<double>{
                eq_ref::phases[1].begin(),
                eq_ref::phases[1].end()},
            std::vector<double>{
                eq_ref::phases[2].begin(),
                eq_ref::phases[2].end()}};

    return flow::NaturalVariableLayout3P{
        flow::
            NaturalVariableCompositionPivot3P::
                select(compositions)};
}

std::vector<double> reference_natural_variables(
    const flow::NaturalVariableLayout3P&
        layout) {
    std::vector<double>
        q(
            layout.unknown_count(),
            0.0);
    q[layout.pressure_unknown_index()] =
        eq_ref::pressure_pa;
    q[layout.temperature_unknown_index()] =
        eq_ref::temperature_k;

    const auto s0 =
        layout.independent_saturation_unknown_index(
            flow::PhaseSlot3::phase0);
    const auto s1 =
        layout.independent_saturation_unknown_index(
            flow::PhaseSlot3::phase1);
    require(
        s0.has_value() &&
            s1.has_value(),
        "Li-Firoozabadi 3P layout lost saturation columns");
    q[*s0] =
        flow_ref::phase_saturation[0];
    q[*s1] =
        flow_ref::phase_saturation[1];

    for (std::size_t phase = 0U;
         phase < 3U;
         ++phase) {
        const auto slot =
            static_cast<flow::PhaseSlot3>(
                phase);
        for (std::size_t component = 0U;
             component <
                 eq_ref::component_ids.size();
             ++component) {
            const auto column =
                layout
                    .independent_composition_unknown_index(
                        slot,
                        component);
            if (column) {
                q[*column] =
                    eq_ref::phases[
                        phase][component];
            }
        }
    }

    return q;
}

} // namespace

void pr76_li_firoozabadi_sour_gas_provider_regression() {
    const auto model =
        fx::model();
    const auto selections =
        reference_selections(model);
    auto closure =
        flow::
            make_pr76_li_firoozabadi_sour_gas_property_closure(
                model,
                selections);

    for (std::size_t phase = 0U;
         phase < 3U;
         ++phase) {
        const std::vector<double>
            composition{
                eq_ref::phases[phase].begin(),
                eq_ref::phases[phase].end()};
        const auto values =
            closure.evaluate(
                phase,
                eq_ref::pressure_pa,
                eq_ref::temperature_k,
                std::span<const double>{
                    composition});

        near(
            values.molar_density_mol_per_m3,
            flow_ref::
                molar_density_mol_per_m3[
                    phase]);
        near(
            values.mass_density_kg_per_m3,
            flow_ref::
                mass_density_kg_per_m3[
                    phase]);
        near(
            values.dynamic_viscosity_pa_s,
            flow_ref::
                dynamic_viscosity_pa_s[
                    phase],
            4.0e-10,
            2.0e-12);
        near(
            values.specific_enthalpy_j_per_kg,
            flow_ref::
                specific_enthalpy_j_per_kg[
                    phase],
            4.0e-10,
            2.0e-7);
        near(
            values.specific_internal_energy_j_per_kg,
            flow_ref::
                specific_internal_energy_j_per_kg[
                    phase],
            4.0e-10,
            2.0e-7);
        near(
            values.mixture_molar_mass_kg_per_mol,
            flow_ref::
                mixture_molar_mass_kg_per_mol[
                    phase]);
    }

    const auto layout =
        reference_layout();
    const auto q =
        reference_natural_variables(layout);
    ad::RuntimeJacobianWorkspace<double, 4U>
        workspace;
    const auto chart =
        flow::
            evaluate_pr76_selected_phase_property_chart(
                flow::NaturalVariableLayoutDescriptor{layout},
                closure.component_ids(),
                q,
                closure,
                workspace);

    require(
        chart.phase_properties.size() == 3U &&
            chart.equilibrium_residual.size() ==
                12U &&
            chart.layout.unknown_count() ==
                19U,
        "Li-Firoozabadi 3P provider chart shape changed");

    double max_equilibrium = 0.0;
    for (double residual :
         chart.equilibrium_residual) {
        max_equilibrium =
            std::max(
                max_equilibrium,
                std::abs(residual));
    }
    require(
        max_equilibrium <= 2.0e-10,
        "Li-Firoozabadi independent three-phase state does not close in production fugacity rows");

    for (std::size_t phase = 0U;
         phase < 3U;
         ++phase) {
        for (const auto* gradient :
             std::array<const std::vector<double>*, 5>{
                 &chart
                      .molar_density_gradient[
                          phase],
                 &chart
                      .mass_density_gradient[
                          phase],
                 &chart
                      .viscosity_gradient[
                          phase],
                 &chart
                      .enthalpy_gradient[
                          phase],
                 &chart
                      .internal_energy_gradient[
                          phase]}) {
            require(
                gradient->size() ==
                        layout.unknown_count() &&
                    std::all_of(
                        gradient->begin(),
                        gradient->end(),
                        [](double value) {
                            return std::isfinite(
                                value);
                        }),
                "Li-Firoozabadi provider published a malformed/non-finite AD gradient");
        }
    }

    const std::vector<double>
        probe_composition{
            eq_ref::phases[0].begin(),
            eq_ref::phases[0].end()};
    const flow::
        Pr76LiFiroozabadiSourGasPropertyProvider
            provider{model};

    bool low_rejected = false;
    try {
        (void)provider
            .ideal_gas_molar_enthalpy_j_per_mol(
                99.0,
                std::span<const double>{
                    probe_composition});
    } catch (const std::out_of_range&) {
        low_rejected = true;
    }

    bool high_rejected = false;
    try {
        (void)provider
            .ideal_gas_molar_enthalpy_j_per_mol(
                298.16,
                std::span<const double>{
                    probe_composition});
    } catch (const std::out_of_range&) {
        high_rejected = true;
    }

    require(
        low_rejected &&
            high_rejected,
        "Li-Firoozabadi low-temperature caloric provider extrapolated outside sourced [100,298.15] K data");
}
