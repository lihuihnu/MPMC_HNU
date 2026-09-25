#ifndef MPMC_FLOW_DISCRETIZATION_TWO_PHASE_TPFA_HPP
#define MPMC_FLOW_DISCRETIZATION_TWO_PHASE_TPFA_HPP

#include <mpmc/flow/two_phase_natural_variable.hpp>
#include <mpmc/flow_discretization/energy_face_flux.hpp>
#include <mpmc/flow_discretization/normalized_component_face_contribution.hpp>
#include <mpmc/flow_discretization/tpfa_phase_darcy_flux.hpp>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <limits>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace mpmc::flow_discretization {

namespace two_phase_tpfa_detail {

[[nodiscard]] inline bool selects_owner(
    mpmc::flow::TwoPhaseUpwindCellSelection
        selection) {
    return selection !=
        mpmc::flow::
            TwoPhaseUpwindCellSelection::
                neighbour_positive_phase_potential;
}

inline void validate_side(
    const mpmc::flow::NaturalVariableCellState2P&
        state,
    const mpmc::flow::
        TwoPhaseMolarDensityNaturalVariableLinearization&
            molar_density,
    const mpmc::flow::
        TwoPhaseTransportNaturalVariableLinearization&
            transport,
    const mpmc::flow::
        TwoPhaseCaloricNaturalVariableLinearization&
            caloric) {
    const auto descriptor =
        state.layout().descriptor();
    if (!mpmc::flow::two_phase_detail::
            same_layout(
                descriptor,
                molar_density.layout) ||
        !mpmc::flow::energy_accumulation_detail::
            same_state_identity(
                mpmc::flow::two_phase_detail::
                    make_state_identity(state),
                transport.state_identity) ||
        !mpmc::flow::energy_accumulation_detail::
            same_state_identity(
                transport.state_identity,
                caloric.state_identity)) {
        throw std::invalid_argument(
            "mpmc::flow_discretization: two-phase exact-state/property identity mismatch");
    }

    const std::size_t q =
        state.layout().unknown_count();
    for (std::size_t phase = 0U;
         phase < 2U;
         ++phase) {
        mpmc::flow::two_phase_detail::
            validate_gradient(
                molar_density.gradient[phase],
                q,
                "molar-density");
        mpmc::flow::two_phase_detail::
            validate_gradient(
                transport.mass_density_gradient[
                    phase],
                q,
                "mass-density");
        mpmc::flow::two_phase_detail::
            validate_gradient(
                caloric.specific_enthalpy_gradient[
                    phase],
                q,
                "enthalpy");
        if (!mpmc::flow::
                component_accumulation_detail::
                    near_roundoff(
                        molar_density
                            .molar_density_mol_per_m3[
                                phase],
                        state.phase_properties(phase)
                            .molar_density_mol_per_m3) ||
            !mpmc::flow::
                component_accumulation_detail::
                    near_roundoff(
                        transport.mass_density_kg_per_m3[
                            phase],
                        state.phase_properties(phase)
                            .mass_density_kg_per_m3) ||
            !mpmc::flow::
                component_accumulation_detail::
                    near_roundoff(
                        caloric
                            .specific_enthalpy_j_per_kg[
                                phase],
                        state.phase_properties(phase)
                            .specific_enthalpy_j_per_kg)) {
            throw std::invalid_argument(
                "mpmc::flow_discretization: two-phase property primal mismatch");
        }
    }
}

} // namespace two_phase_tpfa_detail

struct TwoPhaseTpfaFaceLinearization3D {
    std::array<double, 2>
        phase_volumetric_flux_m3_per_s{};
    std::array<std::vector<double>, 2>
        owner_phase_flux_gradient;
    std::array<std::vector<double>, 2>
        neighbour_phase_flux_gradient;
    NormalizedComponentFaceContributionLinearization3D
        component;
    NormalizedEnergyFaceContributionLinearization3D
        energy;
};

[[nodiscard]] inline
TwoPhaseTpfaFaceLinearization3D
build_two_phase_tpfa_face_linearization(
    const mpmc::discretization::
        TpfaInternalFaceTransmissibilityEntry3D&
            transmissibility,
    const mpmc::flow::
        TwoPhasePotentialUpwindLinearization3D&
            potential,
    const mpmc::flow::NaturalVariableCellState2P&
        owner_state,
    const mpmc::flow::
        TwoPhaseMolarDensityNaturalVariableLinearization&
            owner_molar_density,
    const mpmc::flow::
        TwoPhaseTransportNaturalVariableLinearization&
            owner_transport,
    const mpmc::flow::
        TwoPhaseCaloricNaturalVariableLinearization&
            owner_caloric,
    const mpmc::flow::NaturalVariableCellState2P&
        neighbour_state,
    const mpmc::flow::
        TwoPhaseMolarDensityNaturalVariableLinearization&
            neighbour_molar_density,
    const mpmc::flow::
        TwoPhaseTransportNaturalVariableLinearization&
            neighbour_transport,
    const mpmc::flow::
        TwoPhaseCaloricNaturalVariableLinearization&
            neighbour_caloric,
    TwoCellBulkVolume3D bulk_volume,
    StaticThermalFaceConductance3D
        thermal_conductance) {
    mpmc::flow_discretization::detail::
        validate_materialized_entry(
            transmissibility);
    two_phase_tpfa_detail::validate_side(
        owner_state,
        owner_molar_density,
        owner_transport,
        owner_caloric);
    two_phase_tpfa_detail::validate_side(
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
            "mpmc::flow_discretization: malformed two-phase TPFA face inputs");
    }

    const double tf =
        transmissibility
            .static_transmissibility
            ->face_transmissibility_m3;
    const std::size_t owner_q =
        owner_state.layout().unknown_count();
    const std::size_t neighbour_q =
        neighbour_state.layout().unknown_count();
    const std::size_t n =
        owner_state.layout().component_count();

    std::array<double, 2> phase_flux{};
    std::array<std::vector<double>, 2>
        owner_phase_flux_gradient;
    std::array<std::vector<double>, 2>
        neighbour_phase_flux_gradient;
    std::vector<double> component_rate(
        n,
        0.0);
    std::vector<double> owner_rate_jacobian(
        n * owner_q,
        0.0);
    std::vector<double> neighbour_rate_jacobian(
        n * neighbour_q,
        0.0);

    double advective_energy_rate = 0.0;
    std::vector<double> owner_energy_gradient(
        owner_q,
        0.0);
    std::vector<double> neighbour_energy_gradient(
        neighbour_q,
        0.0);

    for (std::size_t phase = 0U;
         phase < 2U;
         ++phase) {
        owner_phase_flux_gradient[phase]
            .assign(owner_q, 0.0);
        neighbour_phase_flux_gradient[phase]
            .assign(neighbour_q, 0.0);

        const double flux =
            -tf *
            potential.upwind_mobility_per_pa_s[
                phase] *
            potential.phase_potential_difference_pa[
                phase];
        if (!std::isfinite(flux)) {
            throw std::range_error(
                "mpmc::flow_discretization: two-phase Darcy flux is non-finite");
        }
        phase_flux[phase] = flux;

        for (std::size_t column = 0U;
             column < owner_q;
             ++column) {
            owner_phase_flux_gradient[phase][
                column] =
                -tf *
                (potential
                     .upwind_mobility_per_pa_s[
                         phase] *
                     potential
                         .owner_phase_potential_gradient[
                             phase][column] +
                 potential
                     .phase_potential_difference_pa[
                         phase] *
                     potential
                         .owner_upwind_mobility_gradient[
                             phase][column]);
        }
        for (std::size_t column = 0U;
             column < neighbour_q;
             ++column) {
            neighbour_phase_flux_gradient[phase][
                column] =
                -tf *
                (potential
                     .upwind_mobility_per_pa_s[
                         phase] *
                     potential
                         .neighbour_phase_potential_gradient[
                             phase][column] +
                 potential
                     .phase_potential_difference_pa[
                         phase] *
                     potential
                         .neighbour_upwind_mobility_gradient[
                             phase][column]);
        }

        const bool owner_upstream =
            two_phase_tpfa_detail::
                selects_owner(
                    potential.upwind_selection[
                        phase]);
        const auto& upstream_state =
            owner_upstream
                ? owner_state
                : neighbour_state;
        const auto& upstream_density =
            owner_upstream
                ? owner_molar_density
                : neighbour_molar_density;
        const auto& upstream_transport =
            owner_upstream
                ? owner_transport
                : neighbour_transport;
        const auto& upstream_caloric =
            owner_upstream
                ? owner_caloric
                : neighbour_caloric;
        const auto upstream_composition =
            upstream_state.phase_composition(
                phase);
        const double molar_density =
            upstream_density
                .molar_density_mol_per_m3[
                    phase];

        for (std::size_t component = 0U;
             component < n;
             ++component) {
            const double molar_content =
                molar_density *
                upstream_composition[component];
            component_rate[component] +=
                molar_content * flux;

            for (std::size_t column = 0U;
                 column < owner_q;
                 ++column) {
                double dcontent = 0.0;
                if (owner_upstream) {
                    dcontent =
                        owner_molar_density
                            .gradient[phase][column] *
                            upstream_composition[
                                component] +
                        molar_density *
                            mpmc::flow::
                                two_phase_detail::
                                    d_composition(
                                        owner_state
                                            .layout(),
                                        phase,
                                        component,
                                        column);
                }
                owner_rate_jacobian[
                    component * owner_q +
                    column] +=
                    dcontent * flux +
                    molar_content *
                        owner_phase_flux_gradient[
                            phase][column];
            }
            for (std::size_t column = 0U;
                 column < neighbour_q;
                 ++column) {
                double dcontent = 0.0;
                if (!owner_upstream) {
                    dcontent =
                        neighbour_molar_density
                            .gradient[phase][column] *
                            upstream_composition[
                                component] +
                        molar_density *
                            mpmc::flow::
                                two_phase_detail::
                                    d_composition(
                                        neighbour_state
                                            .layout(),
                                        phase,
                                        component,
                                        column);
                }
                neighbour_rate_jacobian[
                    component * neighbour_q +
                    column] +=
                    dcontent * flux +
                    molar_content *
                        neighbour_phase_flux_gradient[
                            phase][column];
            }
        }

        const double rho =
            upstream_transport
                .mass_density_kg_per_m3[
                    phase];
        const double h =
            upstream_caloric
                .specific_enthalpy_j_per_kg[
                    phase];
        const double energy_density =
            rho * h;
        advective_energy_rate +=
            energy_density * flux;

        for (std::size_t column = 0U;
             column < owner_q;
             ++column) {
            double d_energy_density = 0.0;
            if (owner_upstream) {
                d_energy_density =
                    owner_transport
                        .mass_density_gradient[
                            phase][column] *
                        h +
                    rho *
                        owner_caloric
                            .specific_enthalpy_gradient[
                                phase][column];
            }
            owner_energy_gradient[column] +=
                d_energy_density * flux +
                energy_density *
                    owner_phase_flux_gradient[
                        phase][column];
        }
        for (std::size_t column = 0U;
             column < neighbour_q;
             ++column) {
            double d_energy_density = 0.0;
            if (!owner_upstream) {
                d_energy_density =
                    neighbour_transport
                        .mass_density_gradient[
                            phase][column] *
                        h +
                    rho *
                        neighbour_caloric
                            .specific_enthalpy_gradient[
                                phase][column];
            }
            neighbour_energy_gradient[column] +=
                d_energy_density * flux +
                energy_density *
                    neighbour_phase_flux_gradient[
                        phase][column];
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

    for (std::size_t component_index = 0U;
         component_index < n;
         ++component_index) {
        component.owner_component_contribution_mol_per_bulk_m3_s[
            component_index] =
            component_rate[component_index] /
            bulk_volume.owner_bulk_volume_m3;
        component.neighbour_component_contribution_mol_per_bulk_m3_s[
            component_index] =
            -component_rate[component_index] /
            bulk_volume.neighbour_bulk_volume_m3;

        for (std::size_t column = 0U;
             column < owner_q;
             ++column) {
            const double derivative =
                owner_rate_jacobian[
                    component_index * owner_q +
                    column];
            component.owner_row_owner_column_jacobian[
                component_index * owner_q +
                column] =
                derivative /
                bulk_volume.owner_bulk_volume_m3;
            component.neighbour_row_owner_column_jacobian[
                component_index * owner_q +
                column] =
                -derivative /
                bulk_volume.neighbour_bulk_volume_m3;
        }
        for (std::size_t column = 0U;
             column < neighbour_q;
             ++column) {
            const double derivative =
                neighbour_rate_jacobian[
                    component_index *
                        neighbour_q +
                    column];
            component.owner_row_neighbour_column_jacobian[
                component_index *
                    neighbour_q +
                column] =
                derivative /
                bulk_volume.owner_bulk_volume_m3;
            component.neighbour_row_neighbour_column_jacobian[
                component_index *
                    neighbour_q +
                column] =
                -derivative /
                bulk_volume.neighbour_bulk_volume_m3;
        }
    }

    component.owner_total_row_owner_column_gradient =
        owner_total_gradient;
    component.neighbour_total_row_owner_column_gradient =
        owner_total_gradient;
    for (std::size_t column = 0U;
         column < owner_q;
         ++column) {
        component.owner_total_row_owner_column_gradient[
            column] /=
            bulk_volume.owner_bulk_volume_m3;
        component.neighbour_total_row_owner_column_gradient[
            column] *=
            -1.0 /
            bulk_volume.neighbour_bulk_volume_m3;
    }

    component.owner_total_row_neighbour_column_gradient =
        neighbour_total_gradient;
    component.neighbour_total_row_neighbour_column_gradient =
        neighbour_total_gradient;
    for (std::size_t column = 0U;
         column < neighbour_q;
         ++column) {
        component.owner_total_row_neighbour_column_gradient[
            column] /=
            bulk_volume.owner_bulk_volume_m3;
        component.neighbour_total_row_neighbour_column_gradient[
            column] *=
            -1.0 /
            bulk_volume.neighbour_bulk_volume_m3;
    }

    const double conductive =
        thermal_conductance.conductance_w_per_k *
        (owner_state.temperature_k() -
         neighbour_state.temperature_k());
    const double total_energy =
        advective_energy_rate +
        conductive;
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

    energy.owner_row_owner_column_gradient =
        owner_energy_gradient;
    energy.neighbour_row_owner_column_gradient =
        owner_energy_gradient;
    for (std::size_t column = 0U;
         column < owner_q;
         ++column) {
        energy.owner_row_owner_column_gradient[
            column] /=
            bulk_volume.owner_bulk_volume_m3;
        energy.neighbour_row_owner_column_gradient[
            column] *=
            -1.0 /
            bulk_volume.neighbour_bulk_volume_m3;
    }

    energy.owner_row_neighbour_column_gradient =
        neighbour_energy_gradient;
    energy.neighbour_row_neighbour_column_gradient =
        neighbour_energy_gradient;
    for (std::size_t column = 0U;
         column < neighbour_q;
         ++column) {
        energy.owner_row_neighbour_column_gradient[
            column] /=
            bulk_volume.owner_bulk_volume_m3;
        energy.neighbour_row_neighbour_column_gradient[
            column] *=
            -1.0 /
            bulk_volume.neighbour_bulk_volume_m3;
    }

    return {
        phase_flux,
        std::move(owner_phase_flux_gradient),
        std::move(neighbour_phase_flux_gradient),
        std::move(component),
        std::move(energy)};
}

} // namespace mpmc::flow_discretization

#endif // MPMC_FLOW_DISCRETIZATION_TWO_PHASE_TPFA_HPP
