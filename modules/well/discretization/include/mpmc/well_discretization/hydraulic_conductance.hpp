#ifndef MPMC_WELL_DISCRETIZATION_HYDRAULIC_CONDUCTANCE_HPP
#define MPMC_WELL_DISCRETIZATION_HYDRAULIC_CONDUCTANCE_HPP

#include <mpmc/flow/phase_transport.hpp>
#include <mpmc/well/peaceman_well_index_3d.hpp>

#include <cmath>
#include <cstddef>
#include <stdexcept>
#include <string_view>
#include <utility>
#include <vector>

namespace mpmc::well_discretization {

inline constexpr std::string_view
    hydraulic_conductance_convention =
        "well-discretization/hydraulic-conductance/frozen-wi-times-local-phase-mobility/v1";

/// One phase's cell-to-well hydraulic conductance linearized in the exact
/// natural-variable chart carried by LocalPhaseMobilityLinearization3P.
///
/// Contract:
///   C_alpha = WI * lambda_alpha                      [m3/(Pa s)]
///   dC_alpha/dq_j = WI * d(lambda_alpha)/dq_j
///
/// The Peaceman WI is geometry/property data frozen for this local nonlinear
/// linearization, so d(WI)/dq is exactly zero here. This type does not contain
/// pressure drawdown, a well rate, BHP/rate controls, well unknowns, wellbore
/// hydraulics, gravity correction, or a cell-source residual.
struct WellConnectionHydraulicConductanceLinearization3P {
    static constexpr std::string_view convention =
        hydraulic_conductance_convention;
    static constexpr bool
        well_index_derivative_is_zero = true;

    mpmc::flow::NaturalVariableStateIdentity3P
        state_identity;
    mpmc::flow::PhaseSlot3 phase{
        mpmc::flow::PhaseSlot3::phase0};

    double well_index_m3{};
    double phase_mobility_per_pa_s{};
    double hydraulic_conductance_m3_per_pa_s{};
    std::vector<double>
        hydraulic_conductance_gradient;

    [[nodiscard]] double
    d_hydraulic_conductance(
        std::size_t column) const {
        return hydraulic_conductance_gradient.at(
            column);
    }
};

namespace hydraulic_conductance_detail {

inline void validate_well_index(
    const mpmc::well::PeacemanWellIndex3D&
        connection) {
    if (!std::isfinite(connection.well_index_m3) ||
        !(connection.well_index_m3 > 0.0) ||
        !std::isfinite(
            connection.effective_radial_permeability_m2) ||
        !(connection.effective_radial_permeability_m2 > 0.0) ||
        !std::isfinite(connection.completion_length_m) ||
        !(connection.completion_length_m > 0.0) ||
        !std::isfinite(connection.equivalent_radius_m) ||
        !(connection.equivalent_radius_m > 0.0) ||
        !std::isfinite(connection.wellbore_radius_m) ||
        !(connection.wellbore_radius_m > 0.0) ||
        !std::isfinite(connection.logarithmic_denominator) ||
        !(connection.logarithmic_denominator > 0.0)) {
        throw std::invalid_argument(
            "mpmc::well_discretization::make_hydraulic_conductance_linearization_3p: invalid Peaceman connection");
    }
}

} // namespace hydraulic_conductance_detail

/// Multiply one validated Peaceman connection index by one selected phase
/// mobility and its analytic natural-variable gradient.
///
/// Zero mobility is valid. Its conductance is exactly zero while its derivative
/// may remain nonzero, which is important near a constitutive phase-mobility
/// endpoint. Negative or non-finite mobility is rejected.
[[nodiscard]] inline
WellConnectionHydraulicConductanceLinearization3P
make_hydraulic_conductance_linearization_3p(
    const mpmc::well::PeacemanWellIndex3D&
        connection,
    const mpmc::flow::LocalPhaseMobilityLinearization3P&
        mobility,
    mpmc::flow::PhaseSlot3 phase) {
    hydraulic_conductance_detail::
        validate_well_index(connection);

    const std::size_t phase_index =
        static_cast<std::size_t>(phase);
    if (phase_index >=
        mpmc::flow::fixed_three_phase_count) {
        throw std::invalid_argument(
            "mpmc::well_discretization::make_hydraulic_conductance_linearization_3p: invalid phase slot");
    }

    const auto& layout =
        mobility.state_identity.layout;
    if (layout.phase_count() !=
            mpmc::flow::fixed_three_phase_count ||
        mobility.state_identity.component_ids.size() !=
            layout.component_count()) {
        throw std::invalid_argument(
            "mpmc::well_discretization::make_hydraulic_conductance_linearization_3p: mobility state identity is not a complete three-phase natural-variable chart");
    }

    const std::size_t q =
        layout.unknown_count();
    const double phase_mobility =
        mobility.mobility_per_pa_s[
            phase_index];
    const auto& phase_gradient =
        mobility.mobility_gradient[
            phase_index];

    if (!std::isfinite(phase_mobility) ||
        phase_mobility < 0.0 ||
        phase_gradient.size() != q) {
        throw std::invalid_argument(
            "mpmc::well_discretization::make_hydraulic_conductance_linearization_3p: phase mobility value/gradient shape is invalid");
    }
    for (double derivative :
         phase_gradient) {
        if (!std::isfinite(derivative)) {
            throw std::invalid_argument(
                "mpmc::well_discretization::make_hydraulic_conductance_linearization_3p: phase mobility gradient contains non-finite derivative");
        }
    }

    const double conductance =
        connection.well_index_m3 *
        phase_mobility;
    if (!std::isfinite(conductance) ||
        conductance < 0.0) {
        throw std::range_error(
            "mpmc::well_discretization::make_hydraulic_conductance_linearization_3p: hydraulic conductance is outside representable non-negative range");
    }

    std::vector<double> gradient(
        q,
        0.0);
    for (std::size_t column = 0U;
         column < q;
         ++column) {
        gradient[column] =
            connection.well_index_m3 *
            phase_gradient[column];
        if (!std::isfinite(
                gradient[column])) {
            throw std::range_error(
                "mpmc::well_discretization::make_hydraulic_conductance_linearization_3p: hydraulic-conductance derivative is outside representable range");
        }
    }

    return {
        mobility.state_identity,
        phase,
        connection.well_index_m3,
        phase_mobility,
        conductance,
        std::move(gradient)};
}

} // namespace mpmc::well_discretization

#endif // MPMC_WELL_DISCRETIZATION_HYDRAULIC_CONDUCTANCE_HPP
