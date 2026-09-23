#ifndef MPMC_FLOW_DISCRETIZATION_ENERGY_FACE_FLUX_HPP
#define MPMC_FLOW_DISCRETIZATION_ENERGY_FACE_FLUX_HPP

#include <mpmc/flow/energy_accumulation.hpp>
#include <mpmc/flow_discretization/tpfa_phase_darcy_flux.hpp>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <limits>
#include <span>
#include <stdexcept>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

namespace mpmc::flow_discretization {

inline constexpr std::string_view
    internal_energy_face_rate_convention =
        "flow_discretization/internal-energy-face-rate/advection-plus-conduction/v1";
inline constexpr std::string_view
    normalized_energy_face_contribution_convention =
        "flow_discretization/rigid-cell-volume-normalized-energy-face-contribution/v1";

struct StaticThermalFaceConductance3D {
    /// G_f [W / K]. This adapter does not derive it from a conductivity tensor.
    double conductance_w_per_k{};
};

struct InternalEnergyFaceRateLinearization3D {
    static constexpr std::string_view convention =
        internal_energy_face_rate_convention;

    mpmc::mesh::LocalIndex face{
        mpmc::mesh::LocalIndex::value_type{0}};
    StaticThermalFaceConductance3D
        thermal_conductance;
    mpmc::flow::NaturalVariableStateIdentity3P
        owner_state_identity;
    mpmc::flow::NaturalVariableStateIdentity3P
        neighbour_state_identity;

    std::array<double, 3>
        phase_advective_energy_rate_w{};
    double advective_energy_rate_w{};
    double conductive_energy_rate_w{};
    double total_energy_rate_w{};

    std::vector<double> owner_gradient;
    std::vector<double> neighbour_gradient;
};

struct TwoCellEnergyBulkVolume3D {
    double owner_bulk_volume_m3{};
    double neighbour_bulk_volume_m3{};
};

struct NormalizedEnergyFaceContributionLinearization3D {
    static constexpr std::string_view convention =
        normalized_energy_face_contribution_convention;
    static constexpr bool bulk_volume_derivative_is_zero =
        true;

    mpmc::mesh::LocalIndex face{
        mpmc::mesh::LocalIndex::value_type{0}};
    TwoCellEnergyBulkVolume3D bulk_volume;
    mpmc::flow::NaturalVariableStateIdentity3P
        owner_state_identity;
    mpmc::flow::NaturalVariableStateIdentity3P
        neighbour_state_identity;

    double owner_contribution_w_per_bulk_m3{};
    double neighbour_contribution_w_per_bulk_m3{};

    std::vector<double>
        owner_row_owner_column_gradient;
    std::vector<double>
        owner_row_neighbour_column_gradient;
    std::vector<double>
        neighbour_row_owner_column_gradient;
    std::vector<double>
        neighbour_row_neighbour_column_gradient;
};

namespace energy_face_flux_detail {

[[nodiscard]] inline bool near_roundoff(
    double first,
    double second,
    double scale = 0.0) {
    if (!std::isfinite(first) ||
        !std::isfinite(second) ||
        !std::isfinite(scale) ||
        scale < 0.0) {
        return false;
    }
    const double reference =
        std::max(
            {1.0,
             std::abs(first),
             std::abs(second),
             scale});
    return std::abs(first - second) <=
        16384.0 *
            std::numeric_limits<double>::epsilon() *
            reference;
}

inline void require_gradient(
    std::span<const double> gradient,
    std::size_t expected,
    const char* message) {
    if (gradient.size() != expected) {
        throw std::invalid_argument(message);
    }
    for (double value : gradient) {
        if (!std::isfinite(value)) {
            throw std::invalid_argument(
                "mpmc::flow_discretization: energy face gradient contains non-finite derivative");
        }
    }
}

inline void validate_side_properties(
    const mpmc::flow::
        PhaseTransportPropertyNaturalVariableLinearization3P&
            transport,
    const mpmc::flow::
        PhaseCaloricPropertyNaturalVariableLinearization3P&
            caloric,
    const char* side) {
    using mpmc::flow::energy_accumulation_detail::
        same_state_identity;

    if (!same_state_identity(
            transport.state_identity,
            caloric.state_identity)) {
        throw std::invalid_argument(
            std::string{"mpmc::flow_discretization: "} +
            side +
            " transport/caloric state identity mismatch");
    }
    mpmc::flow::phase_transport_detail::
        validate_provenance(
            transport.mass_density_provenance,
            "energy-face-mass-density");
    mpmc::flow::phase_transport_detail::
        validate_provenance(
            caloric.enthalpy_provenance,
            "energy-face-specific-enthalpy");

    const std::size_t q =
        transport.state_identity.layout.unknown_count();
    if (q == 0U ||
        transport.state_identity.component_ids.size() !=
            transport.state_identity.layout.component_count()) {
        throw std::invalid_argument(
            std::string{"mpmc::flow_discretization: invalid "} +
            side +
            " energy-property state identity");
    }

    for (std::size_t phase = 0U;
         phase < mpmc::flow::fixed_three_phase_count;
         ++phase) {
        if (!std::isfinite(
                transport.mass_density_kg_per_m3[phase]) ||
            !(transport.mass_density_kg_per_m3[phase] > 0.0) ||
            !std::isfinite(
                caloric.specific_enthalpy_j_per_kg[phase])) {
            throw std::invalid_argument(
                std::string{"mpmc::flow_discretization: invalid "} +
                side +
                " density/enthalpy primal");
        }
        require_gradient(
            transport.mass_density_gradient[phase],
            q,
            "mpmc::flow_discretization: energy mass-density gradient shape mismatch");
        require_gradient(
            caloric.specific_enthalpy_gradient[phase],
            q,
            "mpmc::flow_discretization: energy enthalpy gradient shape mismatch");
    }
}

[[nodiscard]] inline bool selects_owner(
    mpmc::flow::UpwindCellSelection3P
        selection) {
    using Selection =
        mpmc::flow::UpwindCellSelection3P;
    switch (selection) {
    case Selection::owner_negative_phase_potential:
    case Selection::owner_exact_zero_tie:
        return true;
    case Selection::neighbour_positive_phase_potential:
        return false;
    }
    throw std::invalid_argument(
        "mpmc::flow_discretization: invalid phase upwind selection in energy flux");
}

inline void validate_phase_flux(
    const MaterializedTpfaInternalFacePhaseDarcyFluxLinearization3D&
        flux) {
    const std::size_t owner_q =
        flux.owner_state_identity.layout.unknown_count();
    const std::size_t neighbour_q =
        flux.neighbour_state_identity.layout.unknown_count();
    if (owner_q == 0U ||
        neighbour_q == 0U) {
        throw std::invalid_argument(
            "mpmc::flow_discretization: energy phase-flux chart is empty");
    }

    for (const auto& phase :
         flux.phase) {
        if (!std::isfinite(
                phase.volumetric_flux_m3_per_s)) {
            throw std::invalid_argument(
                "mpmc::flow_discretization: energy phase flux contains non-finite volumetric rate");
        }
        require_gradient(
            phase.owner_flux_gradient,
            owner_q,
            "mpmc::flow_discretization: owner phase-flux gradient shape mismatch");
        require_gradient(
            phase.neighbour_flux_gradient,
            neighbour_q,
            "mpmc::flow_discretization: neighbour phase-flux gradient shape mismatch");
        (void)selects_owner(
            phase.upwind_selection);
    }
}

} // namespace energy_face_flux_detail

[[nodiscard]] inline
InternalEnergyFaceRateLinearization3D
build_internal_energy_face_rate(
    const MaterializedTpfaInternalFacePhaseDarcyFluxLinearization3D&
        phase_flux,
    const mpmc::flow::
        PhaseTransportPropertyNaturalVariableLinearization3P&
            owner_transport,
    const mpmc::flow::
        PhaseCaloricPropertyNaturalVariableLinearization3P&
            owner_caloric,
    const mpmc::flow::
        PhaseTransportPropertyNaturalVariableLinearization3P&
            neighbour_transport,
    const mpmc::flow::
        PhaseCaloricPropertyNaturalVariableLinearization3P&
            neighbour_caloric,
    StaticThermalFaceConductance3D
        thermal_conductance) {
    using namespace energy_face_flux_detail;
    using mpmc::flow::energy_accumulation_detail::
        same_state_identity;

    validate_phase_flux(
        phase_flux);
    validate_side_properties(
        owner_transport,
        owner_caloric,
        "owner");
    validate_side_properties(
        neighbour_transport,
        neighbour_caloric,
        "neighbour");

    if (!same_state_identity(
            phase_flux.owner_state_identity,
            owner_transport.state_identity) ||
        !same_state_identity(
            phase_flux.neighbour_state_identity,
            neighbour_transport.state_identity) ||
        phase_flux.owner_state_identity.component_ids !=
            phase_flux.neighbour_state_identity.component_ids) {
        throw std::invalid_argument(
            "mpmc::flow_discretization: phase flux and energy-property state identities disagree");
    }

    if (!std::isfinite(
            thermal_conductance.conductance_w_per_k) ||
        thermal_conductance.conductance_w_per_k < 0.0) {
        throw std::invalid_argument(
            "mpmc::flow_discretization: thermal face conductance must be finite and non-negative [W/K]");
    }

    const std::size_t owner_q =
        phase_flux.owner_state_identity.layout.unknown_count();
    const std::size_t neighbour_q =
        phase_flux.neighbour_state_identity.layout.unknown_count();

    InternalEnergyFaceRateLinearization3D
        result{
            phase_flux.face,
            thermal_conductance,
            phase_flux.owner_state_identity,
            phase_flux.neighbour_state_identity,
            {},
            0.0,
            0.0,
            0.0,
            std::vector<double>(
                owner_q,
                0.0),
            std::vector<double>(
                neighbour_q,
                0.0)};

    for (std::size_t phase = 0U;
         phase < mpmc::flow::fixed_three_phase_count;
         ++phase) {
        const auto& flux =
            phase_flux.phase[phase];
        const bool owner_upstream =
            selects_owner(
                flux.upwind_selection);

        const auto& transport =
            owner_upstream
                ? owner_transport
                : neighbour_transport;
        const auto& caloric =
            owner_upstream
                ? owner_caloric
                : neighbour_caloric;

        const double rho =
            transport.mass_density_kg_per_m3[
                phase];
        const double h =
            caloric.specific_enthalpy_j_per_kg[
                phase];
        const double energy_density =
            rho * h;
        const double advective =
            energy_density *
            flux.volumetric_flux_m3_per_s;

        if (!std::isfinite(energy_density) ||
            !std::isfinite(advective)) {
            throw std::range_error(
                "mpmc::flow_discretization: advective phase energy rate is non-finite");
        }
        result.phase_advective_energy_rate_w[
            phase] =
            advective;
        result.advective_energy_rate_w +=
            advective;

        for (std::size_t column = 0U;
             column < owner_q;
             ++column) {
            double d_energy_density = 0.0;
            if (owner_upstream) {
                d_energy_density =
                    h *
                        owner_transport
                            .mass_density_gradient[
                                phase][column] +
                    rho *
                        owner_caloric
                            .specific_enthalpy_gradient[
                                phase][column];
            }
            result.owner_gradient[column] +=
                d_energy_density *
                    flux.volumetric_flux_m3_per_s +
                energy_density *
                    flux.owner_flux_gradient[
                        column];
        }

        for (std::size_t column = 0U;
             column < neighbour_q;
             ++column) {
            double d_energy_density = 0.0;
            if (!owner_upstream) {
                d_energy_density =
                    h *
                        neighbour_transport
                            .mass_density_gradient[
                                phase][column] +
                    rho *
                        neighbour_caloric
                            .specific_enthalpy_gradient[
                                phase][column];
            }
            result.neighbour_gradient[column] +=
                d_energy_density *
                    flux.volumetric_flux_m3_per_s +
                energy_density *
                    flux.neighbour_flux_gradient[
                        column];
        }
    }

    const double owner_temperature =
        phase_flux.owner_state_identity
            .temperature_k;
    const double neighbour_temperature =
        phase_flux.neighbour_state_identity
            .temperature_k;
    result.conductive_energy_rate_w =
        thermal_conductance.conductance_w_per_k *
        (owner_temperature -
         neighbour_temperature);
    result.total_energy_rate_w =
        result.advective_energy_rate_w +
        result.conductive_energy_rate_w;

    if (!std::isfinite(
            result.advective_energy_rate_w) ||
        !std::isfinite(
            result.conductive_energy_rate_w) ||
        !std::isfinite(
            result.total_energy_rate_w)) {
        throw std::range_error(
            "mpmc::flow_discretization: total face energy rate is non-finite");
    }

    const std::size_t owner_temperature_column =
        phase_flux.owner_state_identity.layout
            .temperature_unknown_index();
    const std::size_t neighbour_temperature_column =
        phase_flux.neighbour_state_identity.layout
            .temperature_unknown_index();
    result.owner_gradient[
        owner_temperature_column] +=
        thermal_conductance.conductance_w_per_k;
    result.neighbour_gradient[
        neighbour_temperature_column] -=
        thermal_conductance.conductance_w_per_k;

    for (double value :
         result.owner_gradient) {
        if (!std::isfinite(value)) {
            throw std::range_error(
                "mpmc::flow_discretization: owner energy-flux derivative is non-finite");
        }
    }
    for (double value :
         result.neighbour_gradient) {
        if (!std::isfinite(value)) {
            throw std::range_error(
                "mpmc::flow_discretization: neighbour energy-flux derivative is non-finite");
        }
    }

    return result;
}

[[nodiscard]] inline
NormalizedEnergyFaceContributionLinearization3D
normalize_energy_face_rate_by_bulk_volume(
    const InternalEnergyFaceRateLinearization3D&
        face_rate,
    TwoCellEnergyBulkVolume3D bulk_volume) {
    if (!std::isfinite(
            bulk_volume.owner_bulk_volume_m3) ||
        !(bulk_volume.owner_bulk_volume_m3 > 0.0) ||
        !std::isfinite(
            bulk_volume.neighbour_bulk_volume_m3) ||
        !(bulk_volume.neighbour_bulk_volume_m3 > 0.0)) {
        throw std::invalid_argument(
            "mpmc::flow_discretization: energy cell bulk volumes must be finite and strictly positive [m^3]");
    }

    const std::size_t owner_q =
        face_rate.owner_state_identity.layout
            .unknown_count();
    const std::size_t neighbour_q =
        face_rate.neighbour_state_identity.layout
            .unknown_count();
    energy_face_flux_detail::require_gradient(
        face_rate.owner_gradient,
        owner_q,
        "mpmc::flow_discretization: owner energy-rate gradient shape mismatch");
    energy_face_flux_detail::require_gradient(
        face_rate.neighbour_gradient,
        neighbour_q,
        "mpmc::flow_discretization: neighbour energy-rate gradient shape mismatch");

    NormalizedEnergyFaceContributionLinearization3D
        result{
            face_rate.face,
            bulk_volume,
            face_rate.owner_state_identity,
            face_rate.neighbour_state_identity,
            face_rate.total_energy_rate_w /
                bulk_volume.owner_bulk_volume_m3,
            -face_rate.total_energy_rate_w /
                bulk_volume.neighbour_bulk_volume_m3,
            std::vector<double>(
                owner_q,
                0.0),
            std::vector<double>(
                neighbour_q,
                0.0),
            std::vector<double>(
                owner_q,
                0.0),
            std::vector<double>(
                neighbour_q,
                0.0)};

    for (std::size_t column = 0U;
         column < owner_q;
         ++column) {
        result.owner_row_owner_column_gradient[
            column] =
            face_rate.owner_gradient[column] /
            bulk_volume.owner_bulk_volume_m3;
        result.neighbour_row_owner_column_gradient[
            column] =
            -face_rate.owner_gradient[column] /
            bulk_volume.neighbour_bulk_volume_m3;
    }

    for (std::size_t column = 0U;
         column < neighbour_q;
         ++column) {
        result.owner_row_neighbour_column_gradient[
            column] =
            face_rate.neighbour_gradient[column] /
            bulk_volume.owner_bulk_volume_m3;
        result.neighbour_row_neighbour_column_gradient[
            column] =
            -face_rate.neighbour_gradient[column] /
            bulk_volume.neighbour_bulk_volume_m3;
    }

    const double weighted_primal =
        bulk_volume.owner_bulk_volume_m3 *
            result.owner_contribution_w_per_bulk_m3 +
        bulk_volume.neighbour_bulk_volume_m3 *
            result.neighbour_contribution_w_per_bulk_m3;
    const double primal_scale =
        bulk_volume.owner_bulk_volume_m3 *
            std::abs(
                result.owner_contribution_w_per_bulk_m3) +
        bulk_volume.neighbour_bulk_volume_m3 *
            std::abs(
                result.neighbour_contribution_w_per_bulk_m3);
    if (!energy_face_flux_detail::near_roundoff(
            weighted_primal,
            0.0,
            primal_scale)) {
        throw std::runtime_error(
            "mpmc::flow_discretization: normalized energy face contribution violates volume-weighted conservation");
    }

    for (std::size_t column = 0U;
         column < owner_q;
         ++column) {
        const double weighted =
            bulk_volume.owner_bulk_volume_m3 *
                result.owner_row_owner_column_gradient[
                    column] +
            bulk_volume.neighbour_bulk_volume_m3 *
                result.neighbour_row_owner_column_gradient[
                    column];
        const double scale =
            bulk_volume.owner_bulk_volume_m3 *
                std::abs(
                    result.owner_row_owner_column_gradient[
                        column]) +
            bulk_volume.neighbour_bulk_volume_m3 *
                std::abs(
                    result.neighbour_row_owner_column_gradient[
                        column]);
        if (!energy_face_flux_detail::near_roundoff(
                weighted,
                0.0,
                scale)) {
            throw std::runtime_error(
                "mpmc::flow_discretization: normalized owner-column energy Jacobian violates volume-weighted conservation");
        }
    }
    for (std::size_t column = 0U;
         column < neighbour_q;
         ++column) {
        const double weighted =
            bulk_volume.owner_bulk_volume_m3 *
                result.owner_row_neighbour_column_gradient[
                    column] +
            bulk_volume.neighbour_bulk_volume_m3 *
                result.neighbour_row_neighbour_column_gradient[
                    column];
        const double scale =
            bulk_volume.owner_bulk_volume_m3 *
                std::abs(
                    result.owner_row_neighbour_column_gradient[
                        column]) +
            bulk_volume.neighbour_bulk_volume_m3 *
                std::abs(
                    result.neighbour_row_neighbour_column_gradient[
                        column]);
        if (!energy_face_flux_detail::near_roundoff(
                weighted,
                0.0,
                scale)) {
            throw std::runtime_error(
                "mpmc::flow_discretization: normalized neighbour-column energy Jacobian violates volume-weighted conservation");
        }
    }

    return result;
}

} // namespace mpmc::flow_discretization

#endif // MPMC_FLOW_DISCRETIZATION_ENERGY_FACE_FLUX_HPP
