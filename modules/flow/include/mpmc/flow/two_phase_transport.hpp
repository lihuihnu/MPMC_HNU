#ifndef MPMC_FLOW_TWO_PHASE_TRANSPORT_HPP
#define MPMC_FLOW_TWO_PHASE_TRANSPORT_HPP

#include <mpmc/flow/two_phase_properties.hpp>
#include <mpmc/flow/phase_potential_upwind.hpp>

#include <array>
#include <cmath>
#include <cstddef>
#include <stdexcept>
#include <utility>
#include <vector>

// Mobility and oriented face potentials on a frozen upwind branch.
namespace mpmc::flow {

struct TwoPhaseMobilityLinearization {
    NaturalVariableStateIdentity state_identity;
    std::array<double, 2>
        phase_pressure_pa{};
    std::array<double, 2>
        mass_density_kg_per_m3{};
    std::array<double, 2>
        dynamic_viscosity_pa_s{};
    std::array<double, 2>
        relative_permeability{};
    std::array<double, 2>
        mobility_per_pa_s{};
    std::array<std::vector<double>, 2>
        phase_pressure_gradient;
    std::array<std::vector<double>, 2>
        mass_density_gradient;
    std::array<std::vector<double>, 2>
        mobility_gradient;
};


[[nodiscard]] inline
TwoPhaseMobilityLinearization
build_two_phase_mobility_linearization(
    const NaturalVariableCellState2P& state,
    const TwoPhaseTransportNaturalVariableLinearization&
        transport) {
    const auto identity =
        two_phase_detail::make_state_identity(
            state);
    if (!validation_detail::same_state_identity(
            identity, transport.state_identity,
            phase_transport_detail::near_roundoff)) {
        throw std::invalid_argument(
            "mpmc::flow: two-phase mobility state identity mismatch");
    }

    TwoPhaseMobilityLinearization result;
    result.state_identity = identity;
    const std::size_t q =
        state.layout().unknown_count();
    for (std::size_t phase = 0U;
         phase < 2U;
         ++phase) {
        two_phase_detail::validate_gradient(
            transport.mass_density_gradient[
                phase],
            q,
            "mass-density");
        two_phase_detail::validate_gradient(
            transport.dynamic_viscosity_gradient[
                phase],
            q,
            "viscosity");
        two_phase_detail::validate_gradient(
            transport.relative_permeability_gradient[
                phase],
            q,
            "relative-permeability");

        const double mu =
            transport.dynamic_viscosity_pa_s[
                phase];
        const double kr =
            transport.relative_permeability[
                phase];
        if (!std::isfinite(mu) ||
            !(mu > 0.0) ||
            !std::isfinite(kr) ||
            kr < 0.0) {
            throw std::invalid_argument(
                "mpmc::flow: two-phase mobility input is invalid");
        }

        result.phase_pressure_pa[phase] =
            state.reference_pressure_pa();
        result.mass_density_kg_per_m3[phase] =
            transport.mass_density_kg_per_m3[
                phase];
        result.dynamic_viscosity_pa_s[phase] =
            mu;
        result.relative_permeability[phase] =
            kr;
        result.mobility_per_pa_s[phase] =
            kr / mu;
        result.phase_pressure_gradient[phase]
            .assign(q, 0.0);
        result.phase_pressure_gradient[phase][
            state.layout()
                .pressure_unknown_index()] =
            1.0;
        result.mass_density_gradient[phase] =
            transport.mass_density_gradient[
                phase];
        result.mobility_gradient[phase]
            .assign(q, 0.0);
        const double mu2 = mu * mu;
        for (std::size_t column = 0U;
             column < q;
             ++column) {
            result.mobility_gradient[phase][
                column] =
                (transport
                     .relative_permeability_gradient[
                         phase][column] *
                     mu -
                 kr *
                     transport
                         .dynamic_viscosity_gradient[
                             phase][column]) /
                mu2;
            if (!std::isfinite(
                    result.mobility_gradient[phase][
                        column])) {
                throw std::range_error(
                    "mpmc::flow: two-phase mobility derivative is non-finite");
            }
        }
    }
    return result;
}

enum class TwoPhaseUpwindCellSelection {
    owner_negative_phase_potential,
    neighbour_positive_phase_potential,
    owner_exact_zero_tie
};

struct TwoPhasePotentialUpwindLinearization3D {
    NaturalVariableStateIdentity
        owner_state_identity;
    NaturalVariableStateIdentity
        neighbour_state_identity;
    std::array<double, 2>
        face_mass_density_kg_per_m3{};
    double gravity_projection_m2_per_s2{};
    std::array<double, 2>
        phase_potential_difference_pa{};
    std::array<TwoPhaseUpwindCellSelection, 2>
        upwind_selection{
            TwoPhaseUpwindCellSelection::
                owner_exact_zero_tie,
            TwoPhaseUpwindCellSelection::
                owner_exact_zero_tie};
    std::array<double, 2>
        upwind_mobility_per_pa_s{};
    std::array<std::vector<double>, 2>
        owner_phase_potential_gradient;
    std::array<std::vector<double>, 2>
        neighbour_phase_potential_gradient;
    std::array<std::vector<double>, 2>
        owner_upwind_mobility_gradient;
    std::array<std::vector<double>, 2>
        neighbour_upwind_mobility_gradient;
};

[[nodiscard]] inline
TwoPhasePotentialUpwindLinearization3D
build_two_phase_potential_upwind_linearization(
    const TwoPhaseMobilityLinearization&
        owner,
    const TwoPhaseMobilityLinearization&
        neighbour,
    GravityVector3D gravity,
    OwnerToNeighbourDisplacement3D
        displacement) {
    if (owner.state_identity.component_ids !=
            neighbour.state_identity.component_ids ||
        owner.state_identity.layout.phase_count() !=
            2U ||
        neighbour.state_identity.layout.phase_count() !=
            2U ||
        owner.state_identity.layout.component_count() !=
            neighbour.state_identity.layout.component_count()) {
        throw std::invalid_argument(
            "mpmc::flow: two-phase owner/neighbour state identity mismatch");
    }
    phase_potential_upwind_detail::
        require_finite_geometry(
            gravity,
            displacement);

    TwoPhasePotentialUpwindLinearization3D result;
    result.owner_state_identity =
        owner.state_identity;
    result.neighbour_state_identity =
        neighbour.state_identity;
    result.gravity_projection_m2_per_s2 =
        phase_potential_upwind_detail::dot(
            gravity,
            displacement);

    const std::size_t owner_q =
        owner.state_identity.layout.unknown_count();
    const std::size_t neighbour_q =
        neighbour.state_identity.layout.unknown_count();

    for (std::size_t phase = 0U;
         phase < 2U;
         ++phase) {
        const double face_density =
            0.5 *
            (owner.mass_density_kg_per_m3[
                 phase] +
             neighbour.mass_density_kg_per_m3[
                 phase]);
        const double delta =
            (neighbour.phase_pressure_pa[
                 phase] -
             owner.phase_pressure_pa[phase]) -
            face_density *
                result.gravity_projection_m2_per_s2;
        if (!std::isfinite(face_density) ||
            !(face_density > 0.0) ||
            !std::isfinite(delta)) {
            throw std::range_error(
                "mpmc::flow: two-phase face density/potential is invalid");
        }

        result.face_mass_density_kg_per_m3[
            phase] =
            face_density;
        result.phase_potential_difference_pa[
            phase] =
            delta;
        result.owner_phase_potential_gradient[
            phase].assign(owner_q, 0.0);
        result.neighbour_phase_potential_gradient[
            phase].assign(neighbour_q, 0.0);

        for (std::size_t column = 0U;
             column < owner_q;
             ++column) {
            result.owner_phase_potential_gradient[
                phase][column] =
                -owner
                     .phase_pressure_gradient[phase][
                         column] -
                0.5 *
                    owner
                        .mass_density_gradient[phase][
                            column] *
                    result
                        .gravity_projection_m2_per_s2;
        }
        for (std::size_t column = 0U;
             column < neighbour_q;
             ++column) {
            result.neighbour_phase_potential_gradient[
                phase][column] =
                neighbour
                    .phase_pressure_gradient[phase][
                        column] -
                0.5 *
                    neighbour
                        .mass_density_gradient[phase][
                            column] *
                    result
                        .gravity_projection_m2_per_s2;
        }

        auto selection =
            TwoPhaseUpwindCellSelection::
                owner_exact_zero_tie;
        if (delta < 0.0) {
            selection =
                TwoPhaseUpwindCellSelection::
                    owner_negative_phase_potential;
        } else if (delta > 0.0) {
            selection =
                TwoPhaseUpwindCellSelection::
                    neighbour_positive_phase_potential;
        }
        result.upwind_selection[phase] =
            selection;

        result.owner_upwind_mobility_gradient[
            phase].assign(owner_q, 0.0);
        result.neighbour_upwind_mobility_gradient[
            phase].assign(neighbour_q, 0.0);

        if (selection ==
            TwoPhaseUpwindCellSelection::
                neighbour_positive_phase_potential) {
            result.upwind_mobility_per_pa_s[
                phase] =
                neighbour.mobility_per_pa_s[
                    phase];
            result.neighbour_upwind_mobility_gradient[
                phase] =
                neighbour.mobility_gradient[
                    phase];
        } else {
            result.upwind_mobility_per_pa_s[
                phase] =
                owner.mobility_per_pa_s[
                    phase];
            result.owner_upwind_mobility_gradient[
                phase] =
                owner.mobility_gradient[
                    phase];
        }
    }
    return result;
}

} // namespace mpmc::flow

#endif // MPMC_FLOW_TWO_PHASE_TRANSPORT_HPP
