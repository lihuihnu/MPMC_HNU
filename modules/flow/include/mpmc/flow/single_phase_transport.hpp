#ifndef MPMC_FLOW_SINGLE_PHASE_TRANSPORT_HPP
#define MPMC_FLOW_SINGLE_PHASE_TRANSPORT_HPP

#include <mpmc/flow/single_phase_properties.hpp>
#include <mpmc/flow/phase_potential_upwind.hpp>

#include <array>
#include <cmath>
#include <cstddef>
#include <stdexcept>
#include <utility>
#include <vector>

// Mobility and oriented face potentials on a frozen upwind branch.
namespace mpmc::flow {

struct SinglePhaseMobilityLinearization {
    NaturalVariableStateIdentity state_identity;
    double phase_pressure_pa{};
    double mass_density_kg_per_m3{};
    double dynamic_viscosity_pa_s{};
    double relative_permeability{};
    double mobility_per_pa_s{};
    std::vector<double> phase_pressure_gradient;
    std::vector<double> mass_density_gradient;
    std::vector<double> mobility_gradient;
};


[[nodiscard]] inline
SinglePhaseMobilityLinearization
build_single_phase_mobility_linearization(
    const NaturalVariableCellState1P& state,
    const SinglePhaseTransportNaturalVariableLinearization&
        transport) {
    const auto identity =
        single_phase_detail::
            make_state_identity(state);
    if (!validation_detail::same_state_identity(
            identity, transport.state_identity,
            phase_transport_detail::near_roundoff)) {
        throw std::invalid_argument(
            "mpmc::flow: single-phase mobility state identity mismatch");
    }
    const std::size_t q =
        state.layout().unknown_count();
    single_phase_detail::validate_gradient(
        transport.mass_density_gradient,
        q,
        "mass-density");
    single_phase_detail::validate_gradient(
        transport.dynamic_viscosity_gradient,
        q,
        "viscosity");
    single_phase_detail::validate_gradient(
        transport.relative_permeability_gradient,
        q,
        "relative-permeability");

    const double mu =
        transport.dynamic_viscosity_pa_s;
    const double kr =
        transport.relative_permeability;
    if (!std::isfinite(mu) ||
        !(mu > 0.0) ||
        !std::isfinite(kr) ||
        kr < 0.0) {
        throw std::invalid_argument(
            "mpmc::flow: single-phase mobility inputs are invalid");
    }
    const double mobility =
        kr / mu;
    std::vector<double> mobility_gradient(
        q,
        0.0);
    const double mu2 =
        mu * mu;
    for (std::size_t column = 0U;
         column < q;
         ++column) {
        mobility_gradient[column] =
            (transport
                 .relative_permeability_gradient[column] *
                 mu -
             kr *
                 transport
                     .dynamic_viscosity_gradient[column]) /
            mu2;
        if (!std::isfinite(
                mobility_gradient[column])) {
            throw std::range_error(
                "mpmc::flow: single-phase mobility derivative is non-finite");
        }
    }

    std::vector<double>
        phase_pressure_gradient(
            q,
            0.0);
    phase_pressure_gradient[
        state.layout().pressure_unknown_index()] =
        1.0;

    return {
        identity,
        state.reference_pressure_pa(),
        transport.mass_density_kg_per_m3,
        transport.dynamic_viscosity_pa_s,
        kr,
        mobility,
        std::move(
            phase_pressure_gradient),
        transport.mass_density_gradient,
        std::move(mobility_gradient)};
}

enum class SinglePhaseUpwindCellSelection {
    owner_negative_phase_potential,
    neighbour_positive_phase_potential,
    owner_exact_zero_tie
};

struct SinglePhasePotentialUpwindLinearization3D {
    NaturalVariableStateIdentity
        owner_state_identity;
    NaturalVariableStateIdentity
        neighbour_state_identity;
    double face_mass_density_kg_per_m3{};
    double gravity_projection_m2_per_s2{};
    double phase_potential_difference_pa{};
    SinglePhaseUpwindCellSelection
        upwind_selection{
            SinglePhaseUpwindCellSelection::
                owner_exact_zero_tie};
    double upwind_mobility_per_pa_s{};
    std::vector<double>
        owner_phase_potential_gradient;
    std::vector<double>
        neighbour_phase_potential_gradient;
    std::vector<double>
        owner_upwind_mobility_gradient;
    std::vector<double>
        neighbour_upwind_mobility_gradient;
};

[[nodiscard]] inline
SinglePhasePotentialUpwindLinearization3D
build_single_phase_potential_upwind_linearization(
    const SinglePhaseMobilityLinearization&
        owner,
    const SinglePhaseMobilityLinearization&
        neighbour,
    GravityVector3D gravity,
    OwnerToNeighbourDisplacement3D
        displacement) {
    if (owner.state_identity.component_ids !=
            neighbour.state_identity.component_ids ||
        owner.state_identity.layout.phase_count() !=
            1U ||
        neighbour.state_identity.layout.phase_count() !=
            1U ||
        owner.state_identity.layout.component_count() !=
            neighbour.state_identity.layout.component_count()) {
        throw std::invalid_argument(
            "mpmc::flow: single-phase owner/neighbour state identity mismatch");
    }
    phase_potential_upwind_detail::
        require_finite_geometry(
            gravity,
            displacement);
    const double gravity_projection =
        phase_potential_upwind_detail::dot(
            gravity,
            displacement);
    const double face_density =
        0.5 *
        (owner.mass_density_kg_per_m3 +
         neighbour.mass_density_kg_per_m3);
    const double delta =
        (neighbour.phase_pressure_pa -
         owner.phase_pressure_pa) -
        face_density *
            gravity_projection;
    if (!std::isfinite(face_density) ||
        !(face_density > 0.0) ||
        !std::isfinite(delta)) {
        throw std::range_error(
            "mpmc::flow: single-phase face density/potential is invalid");
    }

    const std::size_t owner_q =
        owner.state_identity.layout.unknown_count();
    const std::size_t neighbour_q =
        neighbour.state_identity.layout.unknown_count();
    std::vector<double> owner_potential(
        owner_q,
        0.0);
    std::vector<double> neighbour_potential(
        neighbour_q,
        0.0);
    for (std::size_t column = 0U;
         column < owner_q;
         ++column) {
        owner_potential[column] =
            -owner.phase_pressure_gradient[column] -
            0.5 *
                owner.mass_density_gradient[column] *
                gravity_projection;
    }
    for (std::size_t column = 0U;
         column < neighbour_q;
         ++column) {
        neighbour_potential[column] =
            neighbour.phase_pressure_gradient[column] -
            0.5 *
                neighbour.mass_density_gradient[column] *
                gravity_projection;
    }

    SinglePhaseUpwindCellSelection selection =
        SinglePhaseUpwindCellSelection::
            owner_exact_zero_tie;
    if (delta < 0.0) {
        selection =
            SinglePhaseUpwindCellSelection::
                owner_negative_phase_potential;
    } else if (delta > 0.0) {
        selection =
            SinglePhaseUpwindCellSelection::
                neighbour_positive_phase_potential;
    }

    std::vector<double> owner_upwind(
        owner_q,
        0.0);
    std::vector<double> neighbour_upwind(
        neighbour_q,
        0.0);
    double mobility = 0.0;
    if (selection ==
            SinglePhaseUpwindCellSelection::
                neighbour_positive_phase_potential) {
        mobility =
            neighbour.mobility_per_pa_s;
        neighbour_upwind =
            neighbour.mobility_gradient;
    } else {
        mobility =
            owner.mobility_per_pa_s;
        owner_upwind =
            owner.mobility_gradient;
    }
    if (!std::isfinite(mobility) ||
        mobility < 0.0) {
        throw std::range_error(
            "mpmc::flow: single-phase upwind mobility is invalid");
    }

    return {
        owner.state_identity,
        neighbour.state_identity,
        face_density,
        gravity_projection,
        delta,
        selection,
        mobility,
        std::move(owner_potential),
        std::move(neighbour_potential),
        std::move(owner_upwind),
        std::move(neighbour_upwind)};
}

} // namespace mpmc::flow

#endif // MPMC_FLOW_SINGLE_PHASE_TRANSPORT_HPP
