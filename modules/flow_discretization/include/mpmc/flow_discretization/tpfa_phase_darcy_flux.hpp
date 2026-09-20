#ifndef MPMC_FLOW_DISCRETIZATION_TPFA_PHASE_DARCY_FLUX_HPP
#define MPMC_FLOW_DISCRETIZATION_TPFA_PHASE_DARCY_FLUX_HPP

#include <mpmc/discretization/tpfa_internal_face_transmissibility_snapshot_3d.hpp>
#include <mpmc/flow/phase_potential_upwind.hpp>

#include <array>
#include <cmath>
#include <cstddef>
#include <limits>
#include <stdexcept>
#include <string_view>
#include <vector>

namespace mpmc::flow_discretization {

inline constexpr std::string_view
    materialized_tpfa_phase_flux_convention =
        "flow_discretization/materialized-tpfa-phase-volumetric-flux/v1";

struct MaterializedTpfaPhaseDarcyFluxEntry3D {
    double volumetric_flux_m3_per_s{};

    std::vector<double> owner_flux_gradient;
    std::vector<double> neighbour_flux_gradient;

    mpmc::flow::UpwindCellSelection3P upwind_selection{
        mpmc::flow::UpwindCellSelection3P::
            owner_exact_zero_tie};

    double phase_potential_difference_pa{};
    double upwind_mobility_per_pa_s{};
};

/// Phase volumetric Darcy flux on one already-materialized internal TPFA face.
///
/// Sign convention:
///   positive F_alpha  => owner -> neighbour
///   negative F_alpha  => neighbour -> owner
///
/// The static transmissibility is [m^3], mobility is [1/(Pa s)], and
/// phase-potential difference is [Pa], therefore F_alpha is [m^3/s].
///
/// T_f is treated as frozen static geometry/permeability data in this contract.
/// No derivative of T_f is included.
struct MaterializedTpfaInternalFacePhaseDarcyFluxLinearization3D {
    static constexpr std::string_view convention =
        materialized_tpfa_phase_flux_convention;

    mpmc::mesh::LocalIndex face{
        mpmc::mesh::LocalIndex::value_type{0}};
    double static_transmissibility_m3{};

    mpmc::flow::NaturalVariableStateIdentity3P
        owner_state_identity;
    mpmc::flow::NaturalVariableStateIdentity3P
        neighbour_state_identity;

    std::array<MaterializedTpfaPhaseDarcyFluxEntry3D, 3>
        phase;
};

namespace detail {

inline void validate_materialized_entry(
    const mpmc::discretization::
        TpfaInternalFaceTransmissibilityEntry3D&
            entry) {
    using mpmc::discretization::
        CombinedTransmissibilityAdmissibilityDisposition3D;
    using mpmc::discretization::
        TpfaInternalFaceTransmissibilityDisposition3D;
    using mpmc::discretization::
        TpfaStaticFaceTransmissibilityDisposition3D;

    if (entry.disposition !=
            TpfaInternalFaceTransmissibilityDisposition3D::
                materialized ||
        entry.admissibility.disposition !=
            CombinedTransmissibilityAdmissibilityDisposition3D::
                direct_normal_projection_k_orthogonal_candidate ||
        !entry.static_transmissibility.has_value()) {
        throw std::invalid_argument(
            "mpmc::flow_discretization: only admissibility-gate materialized TPFA internal faces may produce phase flux");
    }

    const auto& static_tf =
        *entry.static_transmissibility;
    if (static_tf.disposition !=
            TpfaStaticFaceTransmissibilityDisposition3D::
                positive_harmonic_combination ||
        !std::isfinite(static_tf.face_area_m2) ||
        !(static_tf.face_area_m2 > 0.0) ||
        !std::isfinite(
            static_tf.face_transmissibility_m3) ||
        !(static_tf.face_transmissibility_m3 > 0.0)) {
        throw std::invalid_argument(
            "mpmc::flow_discretization: materialized TPFA face must carry positive static transmissibility");
    }
}

inline void validate_potential_payload(
    const mpmc::flow::
        TwoCellPhasePotentialUpwindLinearization3P&
            potential) {
    const std::size_t owner_q =
        potential.owner_state_identity.layout
            .unknown_count();
    const std::size_t neighbour_q =
        potential.neighbour_state_identity.layout
            .unknown_count();

    if (owner_q == 0U ||
        neighbour_q == 0U) {
        throw std::invalid_argument(
            "mpmc::flow_discretization: phase-potential natural-variable chart is empty");
    }

    for (const auto& phase :
         potential.phase) {
        if (!std::isfinite(
                phase.phase_potential_difference_pa) ||
            !std::isfinite(
                phase.upwind_mobility_per_pa_s) ||
            phase.upwind_mobility_per_pa_s < 0.0 ||
            phase.owner_phase_potential_gradient.size() !=
                owner_q ||
            phase.neighbour_phase_potential_gradient.size() !=
                neighbour_q ||
            phase.owner_upwind_mobility_gradient.size() !=
                owner_q ||
            phase.neighbour_upwind_mobility_gradient.size() !=
                neighbour_q) {
            throw std::invalid_argument(
                "mpmc::flow_discretization: malformed phase-potential/upwind payload");
        }

        const auto require_finite =
            [](const std::vector<double>& gradient) {
                for (double value : gradient) {
                    if (!std::isfinite(value)) {
                        throw std::invalid_argument(
                            "mpmc::flow_discretization: phase-potential/upwind gradient contains non-finite derivative");
                    }
                }
            };
        require_finite(
            phase.owner_phase_potential_gradient);
        require_finite(
            phase.neighbour_phase_potential_gradient);
        require_finite(
            phase.owner_upwind_mobility_gradient);
        require_finite(
            phase.neighbour_upwind_mobility_gradient);
    }
}

} // namespace detail

[[nodiscard]] inline
MaterializedTpfaInternalFacePhaseDarcyFluxLinearization3D
build_materialized_tpfa_internal_face_phase_darcy_flux(
    const mpmc::discretization::
        TpfaInternalFaceTransmissibilityEntry3D&
            transmissibility,
    const mpmc::flow::
        TwoCellPhasePotentialUpwindLinearization3P&
            potential) {
    detail::validate_materialized_entry(
        transmissibility);
    detail::validate_potential_payload(
        potential);

    const double tf =
        transmissibility.static_transmissibility
            ->face_transmissibility_m3;
    const std::size_t owner_q =
        potential.owner_state_identity.layout
            .unknown_count();
    const std::size_t neighbour_q =
        potential.neighbour_state_identity.layout
            .unknown_count();

    MaterializedTpfaInternalFacePhaseDarcyFluxLinearization3D
        result{
            transmissibility.face,
            tf,
            potential.owner_state_identity,
            potential.neighbour_state_identity,
            {}};

    for (std::size_t phase = 0U;
         phase < result.phase.size();
         ++phase) {
        const auto& source =
            potential.phase[phase];
        auto& target =
            result.phase[phase];

        target.upwind_selection =
            source.upwind_selection;
        target.phase_potential_difference_pa =
            source.phase_potential_difference_pa;
        target.upwind_mobility_per_pa_s =
            source.upwind_mobility_per_pa_s;

        target.volumetric_flux_m3_per_s =
            -tf *
            source.upwind_mobility_per_pa_s *
            source.phase_potential_difference_pa;

        if (!std::isfinite(
                target.volumetric_flux_m3_per_s)) {
            throw std::range_error(
                "mpmc::flow_discretization: phase volumetric Darcy flux is non-finite");
        }

        target.owner_flux_gradient.assign(
            owner_q,
            0.0);
        target.neighbour_flux_gradient.assign(
            neighbour_q,
            0.0);

        for (std::size_t column = 0U;
             column < owner_q;
             ++column) {
            const double derivative =
                -tf *
                (source.upwind_mobility_per_pa_s *
                     source.owner_phase_potential_gradient[
                         column] +
                 source.phase_potential_difference_pa *
                     source.owner_upwind_mobility_gradient[
                         column]);
            if (!std::isfinite(derivative)) {
                throw std::range_error(
                    "mpmc::flow_discretization: owner phase-flux derivative is non-finite");
            }
            target.owner_flux_gradient[column] =
                derivative;
        }

        for (std::size_t column = 0U;
             column < neighbour_q;
             ++column) {
            const double derivative =
                -tf *
                (source.upwind_mobility_per_pa_s *
                     source.neighbour_phase_potential_gradient[
                         column] +
                 source.phase_potential_difference_pa *
                     source.neighbour_upwind_mobility_gradient[
                         column]);
            if (!std::isfinite(derivative)) {
                throw std::range_error(
                    "mpmc::flow_discretization: neighbour phase-flux derivative is non-finite");
            }
            target.neighbour_flux_gradient[column] =
                derivative;
        }
    }

    return result;
}

} // namespace mpmc::flow_discretization

#endif // MPMC_FLOW_DISCRETIZATION_TPFA_PHASE_DARCY_FLUX_HPP
