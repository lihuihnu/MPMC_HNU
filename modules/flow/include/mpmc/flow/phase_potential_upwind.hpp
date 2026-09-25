#ifndef MPMC_FLOW_PHASE_POTENTIAL_UPWIND_HPP
#define MPMC_FLOW_PHASE_POTENTIAL_UPWIND_HPP

#include <mpmc/flow/phase_transport.hpp>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <limits>
#include <stdexcept>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

namespace mpmc::flow {

inline constexpr std::string_view
    two_cell_phase_potential_convention =
        "flow/two-cell-phase-potential-owner-to-neighbour/v1";

enum class FacePhaseDensityPolicy3P {
    arithmetic_mean_owner_neighbour
};

/// Upstream selection is defined for the future owner->neighbour flux
///
///   F_alpha = -T_f * lambda_up * DeltaPhi_alpha,
///
/// with T_f > 0.
///
/// Therefore:
///   DeltaPhi < 0 => F > 0 => owner is upstream.
///   DeltaPhi > 0 => F < 0 => neighbour is upstream.
///
/// At DeltaPhi == 0 the upwind map is non-differentiable. The v1 contract uses
/// the owner branch deterministically and reports that tie explicitly; the
/// published Jacobian is the frozen owner-branch Jacobian, not a claim of a
/// classical derivative of the switching map.
enum class UpwindCellSelection3P {
    owner_negative_phase_potential,
    neighbour_positive_phase_potential,
    owner_exact_zero_tie
};

struct GravityVector3D {
    double x_m_per_s2{};
    double y_m_per_s2{};
    double z_m_per_s2{};
};

/// Frozen geometry from owner cell center to neighbour cell center:
///
///   d_on = x_neighbour - x_owner [m].
///
/// This model-neutral flow contract intentionally does not own mesh indices or
/// reconstruct geometry. A later mesh/flow adapter can populate this directly
/// from the existing mesh owner-to-neighbour displacement contract.
struct OwnerToNeighbourDisplacement3D {
    double x_m{};
    double y_m{};
    double z_m{};
};

struct TwoCellPhasePotentialEntry3P {
    double face_mass_density_kg_per_m3{};
    std::vector<double>
        owner_face_density_gradient;
    std::vector<double>
        neighbour_face_density_gradient;

    double gravity_projection_m2_per_s2{};
    double gravity_pressure_difference_pa{};
    double phase_potential_difference_pa{};

    std::vector<double>
        owner_phase_potential_gradient;
    std::vector<double>
        neighbour_phase_potential_gradient;

    UpwindCellSelection3P upwind_selection{
        UpwindCellSelection3P::
            owner_exact_zero_tie};

    double upwind_mobility_per_pa_s{};
    std::vector<double>
        owner_upwind_mobility_gradient;
    std::vector<double>
        neighbour_upwind_mobility_gradient;
};

/// Two-cell/one-internal-face phase-potential package.
///
/// Potential direction is strictly owner -> neighbour:
///
///   DeltaPhi_alpha =
///       (p_alpha,n - p_alpha,o)
///       - rho_alpha,* * g . (x_n - x_o).
///
/// g is the physical gravitational acceleration vector [m/s^2]. For a
/// coordinate system with z positive upward, ordinary Earth gravity therefore
/// has a negative z component.
///
/// Jacobians are split by cell so owner and neighbour may retain different
/// dependent-component pivots. No combined global column numbering is invented
/// here.
struct TwoCellPhasePotentialUpwindLinearization3P {
    static constexpr std::string_view convention =
        two_cell_phase_potential_convention;

    FacePhaseDensityPolicy3P density_policy{
        FacePhaseDensityPolicy3P::
            arithmetic_mean_owner_neighbour};

    GravityVector3D gravity;
    OwnerToNeighbourDisplacement3D
        owner_to_neighbour_displacement;

    NaturalVariableStateIdentity3P owner_state_identity;
    NaturalVariableStateIdentity3P neighbour_state_identity;

    std::array<TwoCellPhasePotentialEntry3P, 3>
        phase;
};

namespace phase_potential_upwind_detail {

inline void require_finite_geometry(
    GravityVector3D gravity,
    OwnerToNeighbourDisplacement3D displacement) {
    if (!std::isfinite(gravity.x_m_per_s2) ||
        !std::isfinite(gravity.y_m_per_s2) ||
        !std::isfinite(gravity.z_m_per_s2)) {
        throw std::invalid_argument(
            "mpmc::flow: gravity vector must be finite [m/s^2]");
    }
    if (!std::isfinite(displacement.x_m) ||
        !std::isfinite(displacement.y_m) ||
        !std::isfinite(displacement.z_m)) {
        throw std::invalid_argument(
            "mpmc::flow: owner-to-neighbour displacement must be finite [m]");
    }
    const double norm_squared =
        displacement.x_m * displacement.x_m +
        displacement.y_m * displacement.y_m +
        displacement.z_m * displacement.z_m;
    if (!std::isfinite(norm_squared) ||
        !(norm_squared > 0.0)) {
        throw std::invalid_argument(
            "mpmc::flow: owner-to-neighbour displacement must be nonzero");
    }
}

[[nodiscard]] inline double dot(
    GravityVector3D gravity,
    OwnerToNeighbourDisplacement3D displacement) {
    const double value =
        gravity.x_m_per_s2 * displacement.x_m +
        gravity.y_m_per_s2 * displacement.y_m +
        gravity.z_m_per_s2 * displacement.z_m;
    if (!std::isfinite(value)) {
        throw std::range_error(
            "mpmc::flow: gravity/displacement projection is non-finite");
    }
    return value;
}

inline void validate_local_payload(
    const LocalPhaseMobilityLinearization3P& local,
    const char* side) {
    const std::size_t q =
        local.state_identity.layout.unknown_count();
    if (q == 0U ||
        local.state_identity.component_ids.size() !=
            local.state_identity.layout.component_count()) {
        throw std::invalid_argument(
            std::string{"mpmc::flow: invalid "} +
            side +
            " local transport state identity");
    }

    const auto require_gradient =
        [q, side](
            const std::vector<double>& gradient,
            const char* quantity) {
            if (gradient.size() != q) {
                throw std::invalid_argument(
                    std::string{"mpmc::flow: "} +
                    side + " " + quantity +
                    " gradient shape does not match local natural-variable chart");
            }
            for (double value : gradient) {
                if (!std::isfinite(value)) {
                    throw std::invalid_argument(
                        std::string{"mpmc::flow: "} +
                        side + " " + quantity +
                        " gradient contains non-finite derivative");
                }
            }
        };

    for (std::size_t phase = 0U;
         phase < fixed_three_phase_count;
         ++phase) {
        if (!std::isfinite(
                local.phase_pressure_pa[phase]) ||
            !(local.phase_pressure_pa[phase] > 0.0) ||
            !std::isfinite(
                local.mass_density_kg_per_m3[phase]) ||
            !(local.mass_density_kg_per_m3[phase] > 0.0) ||
            !std::isfinite(
                local.mobility_per_pa_s[phase]) ||
            local.mobility_per_pa_s[phase] < 0.0) {
            throw std::invalid_argument(
                std::string{"mpmc::flow: invalid "} +
                side +
                " local phase pressure/density/mobility primal");
        }

        require_gradient(
            local.phase_pressure_gradient[phase],
            "phase-pressure");
        require_gradient(
            local.mass_density_gradient[phase],
            "mass-density");
        require_gradient(
            local.dynamic_viscosity_gradient[phase],
            "viscosity");
        require_gradient(
            local.relative_permeability_gradient[phase],
            "relative-permeability");
        require_gradient(
            local.mobility_gradient[phase],
            "mobility");

        if (!std::isfinite(
                local.dynamic_viscosity_pa_s[phase]) ||
            !(local.dynamic_viscosity_pa_s[phase] > 0.0) ||
            !std::isfinite(
                local.relative_permeability[phase]) ||
            local.relative_permeability[phase] < 0.0) {
            throw std::invalid_argument(
                std::string{"mpmc::flow: invalid "} +
                side +
                " local viscosity/relative-permeability primal");
        }
    }
}

inline void validate_two_cell_identity(
    const LocalPhaseMobilityLinearization3P& owner,
    const LocalPhaseMobilityLinearization3P& neighbour) {
    if (owner.state_identity.component_ids !=
        neighbour.state_identity.component_ids) {
        throw std::invalid_argument(
            "mpmc::flow: owner/neighbour component identity/order mismatch");
    }
    if (owner.state_identity.layout.component_count() !=
        neighbour.state_identity.layout.component_count()) {
        throw std::invalid_argument(
            "mpmc::flow: owner/neighbour component count mismatch");
    }
}

[[nodiscard]] inline double face_density(
    FacePhaseDensityPolicy3P policy,
    double owner_density,
    double neighbour_density) {
    switch (policy) {
    case FacePhaseDensityPolicy3P::
        arithmetic_mean_owner_neighbour:
        return 0.5 *
            (owner_density + neighbour_density);
    }
    throw std::invalid_argument(
        "mpmc::flow: unsupported face phase-density policy");
}

[[nodiscard]] inline UpwindCellSelection3P
select_upwind(double phase_potential_difference_pa) {
    if (phase_potential_difference_pa < 0.0) {
        return UpwindCellSelection3P::
            owner_negative_phase_potential;
    }
    if (phase_potential_difference_pa > 0.0) {
        return UpwindCellSelection3P::
            neighbour_positive_phase_potential;
    }
    return UpwindCellSelection3P::
        owner_exact_zero_tie;
}

[[nodiscard]] inline bool selects_owner(
    UpwindCellSelection3P selection) {
    return selection ==
               UpwindCellSelection3P::
                   owner_negative_phase_potential ||
           selection ==
               UpwindCellSelection3P::
                   owner_exact_zero_tie;
}

} // namespace phase_potential_upwind_detail

[[nodiscard]] inline
TwoCellPhasePotentialUpwindLinearization3P
build_two_cell_phase_potential_upwind_linearization(
    const LocalPhaseMobilityLinearization3P& owner,
    const LocalPhaseMobilityLinearization3P& neighbour,
    GravityVector3D gravity,
    OwnerToNeighbourDisplacement3D
        owner_to_neighbour_displacement,
    FacePhaseDensityPolicy3P density_policy =
        FacePhaseDensityPolicy3P::
            arithmetic_mean_owner_neighbour) {
    using namespace phase_potential_upwind_detail;

    validate_local_payload(owner, "owner");
    validate_local_payload(neighbour, "neighbour");
    validate_two_cell_identity(owner, neighbour);
    require_finite_geometry(
        gravity,
        owner_to_neighbour_displacement);

    const double gravity_projection =
        dot(
            gravity,
            owner_to_neighbour_displacement);

    const std::size_t owner_q =
        owner.state_identity.layout.unknown_count();
    const std::size_t neighbour_q =
        neighbour.state_identity.layout.unknown_count();

    TwoCellPhasePotentialUpwindLinearization3P
        result{
            density_policy,
            gravity,
            owner_to_neighbour_displacement,
            owner.state_identity,
            neighbour.state_identity,
            {}};

    for (std::size_t phase = 0U;
         phase < fixed_three_phase_count;
         ++phase) {
        auto& entry =
            result.phase[phase];

        entry.face_mass_density_kg_per_m3 =
            face_density(
                density_policy,
                owner.mass_density_kg_per_m3[phase],
                neighbour.mass_density_kg_per_m3[phase]);
        if (!std::isfinite(
                entry.face_mass_density_kg_per_m3) ||
            !(entry.face_mass_density_kg_per_m3 >
              0.0)) {
            throw std::range_error(
                "mpmc::flow: face phase density is non-finite or non-positive");
        }

        entry.owner_face_density_gradient.assign(
            owner_q,
            0.0);
        entry.neighbour_face_density_gradient.assign(
            neighbour_q,
            0.0);
        for (std::size_t column = 0U;
             column < owner_q;
             ++column) {
            entry.owner_face_density_gradient[column] =
                0.5 *
                owner.mass_density_gradient[phase][column];
        }
        for (std::size_t column = 0U;
             column < neighbour_q;
             ++column) {
            entry.neighbour_face_density_gradient[column] =
                0.5 *
                neighbour.mass_density_gradient[phase][column];
        }

        entry.gravity_projection_m2_per_s2 =
            gravity_projection;
        entry.gravity_pressure_difference_pa =
            entry.face_mass_density_kg_per_m3 *
            gravity_projection;

        entry.phase_potential_difference_pa =
            (neighbour.phase_pressure_pa[phase] -
             owner.phase_pressure_pa[phase]) -
            entry.gravity_pressure_difference_pa;

        if (!std::isfinite(
                entry.gravity_pressure_difference_pa) ||
            !std::isfinite(
                entry.phase_potential_difference_pa)) {
            throw std::range_error(
                "mpmc::flow: phase potential is non-finite");
        }

        entry.owner_phase_potential_gradient.assign(
            owner_q,
            0.0);
        entry.neighbour_phase_potential_gradient.assign(
            neighbour_q,
            0.0);

        for (std::size_t column = 0U;
             column < owner_q;
             ++column) {
            entry.owner_phase_potential_gradient[column] =
                -owner.phase_pressure_gradient[phase][column] -
                entry.owner_face_density_gradient[column] *
                    gravity_projection;
        }
        for (std::size_t column = 0U;
             column < neighbour_q;
             ++column) {
            entry.neighbour_phase_potential_gradient[column] =
                neighbour.phase_pressure_gradient[phase][column] -
                entry.neighbour_face_density_gradient[column] *
                    gravity_projection;
        }

        entry.upwind_selection =
            select_upwind(
                entry.phase_potential_difference_pa);

        entry.owner_upwind_mobility_gradient.assign(
            owner_q,
            0.0);
        entry.neighbour_upwind_mobility_gradient.assign(
            neighbour_q,
            0.0);

        if (selects_owner(
                entry.upwind_selection)) {
            entry.upwind_mobility_per_pa_s =
                owner.mobility_per_pa_s[phase];
            entry.owner_upwind_mobility_gradient =
                owner.mobility_gradient[phase];
        } else {
            entry.upwind_mobility_per_pa_s =
                neighbour.mobility_per_pa_s[phase];
            entry.neighbour_upwind_mobility_gradient =
                neighbour.mobility_gradient[phase];
        }

        if (!std::isfinite(
                entry.upwind_mobility_per_pa_s) ||
            entry.upwind_mobility_per_pa_s < 0.0) {
            throw std::range_error(
                "mpmc::flow: upwind phase mobility is non-finite or negative");
        }

        for (double value :
             entry.owner_phase_potential_gradient) {
            if (!std::isfinite(value)) {
                throw std::range_error(
                    "mpmc::flow: owner phase-potential gradient is non-finite");
            }
        }
        for (double value :
             entry.neighbour_phase_potential_gradient) {
            if (!std::isfinite(value)) {
                throw std::range_error(
                    "mpmc::flow: neighbour phase-potential gradient is non-finite");
            }
        }
    }

    return result;
}

} // namespace mpmc::flow

#endif // MPMC_FLOW_PHASE_POTENTIAL_UPWIND_HPP
