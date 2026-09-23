#ifndef MPMC_FLOW_DISCRETIZATION_SINGLE_PHASE_TPFA_HPP
#define MPMC_FLOW_DISCRETIZATION_SINGLE_PHASE_TPFA_HPP

#include <mpmc/flow/single_phase_natural_variable.hpp>
#include <mpmc/flow_discretization/energy_face_flux.hpp>
#include <mpmc/flow_discretization/normalized_component_face_contribution.hpp>
#include <mpmc/flow_discretization/tpfa_phase_darcy_flux.hpp>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace mpmc::flow_discretization {

namespace single_phase_tpfa_detail {

[[nodiscard]] inline bool selects_owner(
    mpmc::flow::SinglePhaseUpwindCellSelection
        selection) {
    return selection !=
        mpmc::flow::
            SinglePhaseUpwindCellSelection::
                neighbour_positive_phase_potential;
}

inline void validate_side(
    const mpmc::flow::NaturalVariableCellState1P&
        state,
    const mpmc::flow::
        SinglePhaseMolarDensityNaturalVariableLinearization&
            molar_density,
    const mpmc::flow::
        SinglePhaseTransportNaturalVariableLinearization&
            transport,
    const mpmc::flow::
        SinglePhaseCaloricNaturalVariableLinearization&
            caloric) {
    const auto descriptor =
        state.layout().descriptor();
    if (!mpmc::flow::single_phase_detail::
            same_layout(
                descriptor,
                molar_density.layout) ||
        !mpmc::flow::energy_accumulation_detail::
            same_state_identity(
                mpmc::flow::single_phase_detail::
                    make_state_identity(state),
                transport.state_identity) ||
        !mpmc::flow::energy_accumulation_detail::
            same_state_identity(
                transport.state_identity,
                caloric.state_identity) ||
        molar_density.molar_density_mol_per_m3 !=
            state.phase_properties()
                .molar_density_mol_per_m3 ||
        transport.mass_density_kg_per_m3 !=
            state.phase_properties()
                .mass_density_kg_per_m3 ||
        transport.dynamic_viscosity_pa_s !=
            state.phase_properties()
                .dynamic_viscosity_pa_s ||
        caloric.specific_enthalpy_j_per_kg !=
            state.phase_properties()
                .specific_enthalpy_j_per_kg) {
        throw std::invalid_argument(
            "mpmc::flow_discretization: single-phase exact-state/property identity mismatch");
    }
}

[[nodiscard]] inline double
d_composition(
    const mpmc::flow::NaturalVariableCellState1P&
        state,
    std::size_t component,
    std::size_t column) {
    return mpmc::flow::single_phase_detail::
        d_composition(
            state.layout(),
            component,
            column);
}

} // namespace single_phase_tpfa_detail

struct SinglePhaseTpfaFaceLinearization3D {
    double volumetric_flux_m3_per_s{};
    std::vector<double> owner_flux_gradient;
    std::vector<double> neighbour_flux_gradient;
    NormalizedComponentFaceContributionLinearization3D
        component;
    NormalizedEnergyFaceContributionLinearization3D
        energy;
};

[[nodiscard]] inline
SinglePhaseTpfaFaceLinearization3D
build_single_phase_tpfa_face_linearization(
    const mpmc::discretization::
        TpfaInternalFaceTransmissibilityEntry3D&
            transmissibility,
    const mpmc::flow::
        SinglePhasePotentialUpwindLinearization3D&
            potential,
    const mpmc::flow::NaturalVariableCellState1P&
        owner_state,
    const mpmc::flow::
        SinglePhaseMolarDensityNaturalVariableLinearization&
            owner_molar_density,
    const mpmc::flow::
        SinglePhaseTransportNaturalVariableLinearization&
            owner_transport,
    const mpmc::flow::
        SinglePhaseCaloricNaturalVariableLinearization&
            owner_caloric,
    const mpmc::flow::NaturalVariableCellState1P&
        neighbour_state,
    const mpmc::flow::
        SinglePhaseMolarDensityNaturalVariableLinearization&
            neighbour_molar_density,
    const mpmc::flow::
        SinglePhaseTransportNaturalVariableLinearization&
            neighbour_transport,
    const mpmc::flow::
        SinglePhaseCaloricNaturalVariableLinearization&
            neighbour_caloric,
    TwoCellBulkVolume3D bulk_volume,
    StaticThermalFaceConductance3D
        thermal_conductance) {
    mpmc::flow_discretization::detail::
        validate_materialized_entry(
            transmissibility);
    single_phase_tpfa_detail::validate_side(
        owner_state,
        owner_molar_density,
        owner_transport,
        owner_caloric);
    single_phase_tpfa_detail::validate_side(
        neighbour_state,
        neighbour_molar_density,
        neighbour_transport,
        neighbour_caloric);

    if (owner_state.component_ids().size() !=
            neighbour_state.component_ids().size() ||
        !std::equal(
            owner_state.component_ids().begin(),
            owner_state.component_ids().end(),
            neighbour_state.component_ids().begin()) ||
        !std::isfinite(
            bulk_volume.owner_bulk_volume_m3) ||
        !(bulk_volume.owner_bulk_volume_m3 > 0.0) ||
        !std::isfinite(
            bulk_volume.neighbour_bulk_volume_m3) ||
        !(bulk_volume.neighbour_bulk_volume_m3 > 0.0) ||
        !std::isfinite(
            thermal_conductance.conductance_w_per_k) ||
        thermal_conductance.conductance_w_per_k <
            0.0) {
        throw std::invalid_argument(
            "mpmc::flow_discretization: malformed single-phase TPFA face inputs");
    }

    const auto& static_tf =
        *transmissibility.static_transmissibility;
    const double tf =
        static_tf.face_transmissibility_m3;
    const double flux =
        -tf *
        potential.upwind_mobility_per_pa_s *
        potential.phase_potential_difference_pa;
    if (!std::isfinite(flux)) {
        throw std::range_error(
            "mpmc::flow_discretization: single-phase Darcy flux is non-finite");
    }

    const std::size_t owner_q =
        owner_state.layout().unknown_count();
    const std::size_t neighbour_q =
        neighbour_state.layout().unknown_count();
    if (potential.owner_phase_potential_gradient.size() !=
            owner_q ||
        potential.neighbour_phase_potential_gradient.size() !=
            neighbour_q ||
        potential.owner_upwind_mobility_gradient.size() !=
            owner_q ||
        potential.neighbour_upwind_mobility_gradient.size() !=
            neighbour_q) {
        throw std::invalid_argument(
            "mpmc::flow_discretization: single-phase potential gradient shape mismatch");
    }

    std::vector<double> owner_flux_gradient(
        owner_q,
        0.0);
    std::vector<double> neighbour_flux_gradient(
        neighbour_q,
        0.0);
    for (std::size_t column = 0U;
         column < owner_q;
         ++column) {
        owner_flux_gradient[column] =
            -tf *
            (potential.upwind_mobility_per_pa_s *
                 potential
                     .owner_phase_potential_gradient[
                         column] +
             potential.phase_potential_difference_pa *
                 potential
                     .owner_upwind_mobility_gradient[
                         column]);
    }
    for (std::size_t column = 0U;
         column < neighbour_q;
         ++column) {
        neighbour_flux_gradient[column] =
            -tf *
            (potential.upwind_mobility_per_pa_s *
                 potential
                     .neighbour_phase_potential_gradient[
                         column] +
             potential.phase_potential_difference_pa *
                 potential
                     .neighbour_upwind_mobility_gradient[
                         column]);
    }

    const bool owner_upstream =
        single_phase_tpfa_detail::
            selects_owner(
                potential.upwind_selection);
    const auto& upstream_state =
        owner_upstream
            ? owner_state
            : neighbour_state;
    const auto& upstream_density =
        owner_upstream
            ? owner_molar_density
            : neighbour_molar_density;
    const auto upstream_composition =
        upstream_state.phase_composition();

    const std::size_t n =
        owner_state.layout().component_count();
    std::vector<double> component_rate(
        n,
        0.0);
    std::vector<double> owner_rate_jacobian(
        n * owner_q,
        0.0);
    std::vector<double> neighbour_rate_jacobian(
        n * neighbour_q,
        0.0);

    const double c =
        upstream_density
            .molar_density_mol_per_m3;
    for (std::size_t component = 0U;
         component < n;
         ++component) {
        const double content =
            c *
            upstream_composition[component];
        component_rate[component] =
            content * flux;

        for (std::size_t column = 0U;
             column < owner_q;
             ++column) {
            double dcontent = 0.0;
            if (owner_upstream) {
                dcontent =
                    owner_molar_density
                        .gradient[column] *
                        upstream_composition[component] +
                    c *
                        single_phase_tpfa_detail::
                            d_composition(
                                owner_state,
                                component,
                                column);
            }
            owner_rate_jacobian[
                component * owner_q +
                column] =
                dcontent * flux +
                content *
                    owner_flux_gradient[column];
        }

        for (std::size_t column = 0U;
             column < neighbour_q;
             ++column) {
            double dcontent = 0.0;
            if (!owner_upstream) {
                dcontent =
                    neighbour_molar_density
                        .gradient[column] *
                        upstream_composition[component] +
                    c *
                        single_phase_tpfa_detail::
                            d_composition(
                                neighbour_state,
                                component,
                                column);
            }
            neighbour_rate_jacobian[
                component * neighbour_q +
                column] =
                dcontent * flux +
                content *
                    neighbour_flux_gradient[column];
        }
    }

    double total_rate = 0.0;
    for (double value : component_rate) {
        total_rate += value;
    }
    std::vector<double> owner_total_gradient(
        owner_q,
        0.0);
    std::vector<double> neighbour_total_gradient(
        neighbour_q,
        0.0);
    for (std::size_t component = 0U;
         component < n;
         ++component) {
        for (std::size_t column = 0U;
             column < owner_q;
             ++column) {
            owner_total_gradient[column] +=
                owner_rate_jacobian[
                    component * owner_q +
                    column];
        }
        for (std::size_t column = 0U;
             column < neighbour_q;
             ++column) {
            neighbour_total_gradient[column] +=
                neighbour_rate_jacobian[
                    component * neighbour_q +
                    column];
        }
    }

    NormalizedComponentFaceContributionLinearization3D
        component{
            transmissibility.face,
            bulk_volume,
            std::vector<std::string>{
                owner_state.component_ids().begin(),
                owner_state.component_ids().end()},
            potential.owner_state_identity,
            potential.neighbour_state_identity,
            {},
            {},
            {},
            {},
            {},
            {},
            total_rate /
                bulk_volume.owner_bulk_volume_m3,
            -total_rate /
                bulk_volume.neighbour_bulk_volume_m3,
            {},
            {},
            {},
            {}};

    component.owner_component_contribution_mol_per_bulk_m3_s
        .resize(n);
    component.neighbour_component_contribution_mol_per_bulk_m3_s
        .resize(n);
    component.owner_row_owner_column_jacobian
        .resize(n * owner_q);
    component.owner_row_neighbour_column_jacobian
        .resize(n * neighbour_q);
    component.neighbour_row_owner_column_jacobian
        .resize(n * owner_q);
    component.neighbour_row_neighbour_column_jacobian
        .resize(n * neighbour_q);

    for (std::size_t i = 0U; i < n; ++i) {
        component.owner_component_contribution_mol_per_bulk_m3_s[i] =
            component_rate[i] /
            bulk_volume.owner_bulk_volume_m3;
        component.neighbour_component_contribution_mol_per_bulk_m3_s[i] =
            -component_rate[i] /
            bulk_volume.neighbour_bulk_volume_m3;
        for (std::size_t column = 0U;
             column < owner_q;
             ++column) {
            const double d =
                owner_rate_jacobian[
                    i * owner_q + column];
            component.owner_row_owner_column_jacobian[
                i * owner_q + column] =
                d /
                bulk_volume.owner_bulk_volume_m3;
            component.neighbour_row_owner_column_jacobian[
                i * owner_q + column] =
                -d /
                bulk_volume.neighbour_bulk_volume_m3;
        }
        for (std::size_t column = 0U;
             column < neighbour_q;
             ++column) {
            const double d =
                neighbour_rate_jacobian[
                    i * neighbour_q + column];
            component.owner_row_neighbour_column_jacobian[
                i * neighbour_q + column] =
                d /
                bulk_volume.owner_bulk_volume_m3;
            component.neighbour_row_neighbour_column_jacobian[
                i * neighbour_q + column] =
                -d /
                bulk_volume.neighbour_bulk_volume_m3;
        }
    }

    component.owner_total_row_owner_column_gradient
        .resize(owner_q);
    component.neighbour_total_row_owner_column_gradient
        .resize(owner_q);
    for (std::size_t column = 0U;
         column < owner_q;
         ++column) {
        component.owner_total_row_owner_column_gradient[column] =
            owner_total_gradient[column] /
            bulk_volume.owner_bulk_volume_m3;
        component.neighbour_total_row_owner_column_gradient[column] =
            -owner_total_gradient[column] /
            bulk_volume.neighbour_bulk_volume_m3;
    }
    component.owner_total_row_neighbour_column_gradient
        .resize(neighbour_q);
    component.neighbour_total_row_neighbour_column_gradient
        .resize(neighbour_q);
    for (std::size_t column = 0U;
         column < neighbour_q;
         ++column) {
        component.owner_total_row_neighbour_column_gradient[column] =
            neighbour_total_gradient[column] /
            bulk_volume.owner_bulk_volume_m3;
        component.neighbour_total_row_neighbour_column_gradient[column] =
            -neighbour_total_gradient[column] /
            bulk_volume.neighbour_bulk_volume_m3;
    }

    const auto& upstream_transport =
        owner_upstream
            ? owner_transport
            : neighbour_transport;
    const auto& upstream_caloric =
        owner_upstream
            ? owner_caloric
            : neighbour_caloric;
    const double energy_density =
        upstream_transport.mass_density_kg_per_m3 *
        upstream_caloric.specific_enthalpy_j_per_kg;
    const double advective =
        energy_density * flux;
    const double conductive =
        thermal_conductance.conductance_w_per_k *
        (owner_state.temperature_k() -
         neighbour_state.temperature_k());
    const double total_energy =
        advective + conductive;

    std::vector<double> owner_energy_gradient(
        owner_q,
        0.0);
    std::vector<double> neighbour_energy_gradient(
        neighbour_q,
        0.0);
    for (std::size_t column = 0U;
         column < owner_q;
         ++column) {
        double d_energy_density = 0.0;
        if (owner_upstream) {
            d_energy_density =
                owner_transport
                    .mass_density_gradient[column] *
                    upstream_caloric
                        .specific_enthalpy_j_per_kg +
                upstream_transport
                    .mass_density_kg_per_m3 *
                    owner_caloric
                        .specific_enthalpy_gradient[
                            column];
        }
        owner_energy_gradient[column] =
            d_energy_density * flux +
            energy_density *
                owner_flux_gradient[column];
    }
    for (std::size_t column = 0U;
         column < neighbour_q;
         ++column) {
        double d_energy_density = 0.0;
        if (!owner_upstream) {
            d_energy_density =
                neighbour_transport
                    .mass_density_gradient[column] *
                    upstream_caloric
                        .specific_enthalpy_j_per_kg +
                upstream_transport
                    .mass_density_kg_per_m3 *
                    neighbour_caloric
                        .specific_enthalpy_gradient[
                            column];
        }
        neighbour_energy_gradient[column] =
            d_energy_density * flux +
            energy_density *
                neighbour_flux_gradient[column];
    }
    owner_energy_gradient[
        owner_state.layout()
            .temperature_unknown_index()] +=
        thermal_conductance.conductance_w_per_k;
    neighbour_energy_gradient[
        neighbour_state.layout()
            .temperature_unknown_index()] -=
        thermal_conductance.conductance_w_per_k;

    NormalizedEnergyFaceContributionLinearization3D
        energy{
            transmissibility.face,
            {
                bulk_volume.owner_bulk_volume_m3,
                bulk_volume.neighbour_bulk_volume_m3},
            potential.owner_state_identity,
            potential.neighbour_state_identity,
            total_energy /
                bulk_volume.owner_bulk_volume_m3,
            -total_energy /
                bulk_volume.neighbour_bulk_volume_m3,
            {},
            {},
            {},
            {}};
    energy.owner_row_owner_column_gradient
        .resize(owner_q);
    energy.neighbour_row_owner_column_gradient
        .resize(owner_q);
    for (std::size_t column = 0U;
         column < owner_q;
         ++column) {
        energy.owner_row_owner_column_gradient[column] =
            owner_energy_gradient[column] /
            bulk_volume.owner_bulk_volume_m3;
        energy.neighbour_row_owner_column_gradient[column] =
            -owner_energy_gradient[column] /
            bulk_volume.neighbour_bulk_volume_m3;
    }
    energy.owner_row_neighbour_column_gradient
        .resize(neighbour_q);
    energy.neighbour_row_neighbour_column_gradient
        .resize(neighbour_q);
    for (std::size_t column = 0U;
         column < neighbour_q;
         ++column) {
        energy.owner_row_neighbour_column_gradient[column] =
            neighbour_energy_gradient[column] /
            bulk_volume.owner_bulk_volume_m3;
        energy.neighbour_row_neighbour_column_gradient[column] =
            -neighbour_energy_gradient[column] /
            bulk_volume.neighbour_bulk_volume_m3;
    }

    if (!std::isfinite(total_rate) ||
        !std::isfinite(total_energy)) {
        throw std::range_error(
            "mpmc::flow_discretization: single-phase face rate became non-finite");
    }

    return {
        flux,
        std::move(owner_flux_gradient),
        std::move(neighbour_flux_gradient),
        std::move(component),
        std::move(energy)};
}

} // namespace mpmc::flow_discretization

#endif // MPMC_FLOW_DISCRETIZATION_SINGLE_PHASE_TPFA_HPP
