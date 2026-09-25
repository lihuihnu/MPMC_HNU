#ifndef MPMC_WELL_DISCRETIZATION_CELL_SOURCE_ADAPTER_HPP
#define MPMC_WELL_DISCRETIZATION_CELL_SOURCE_ADAPTER_HPP

#include <mpmc/flow_discretization/cell_source.hpp>
#include <mpmc/well_discretization/component_molar_rate.hpp>
#include <mpmc/well_discretization/energy_rate.hpp>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <stdexcept>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

namespace mpmc::well_discretization {

inline constexpr std::string_view
    well_connection_cell_source_adapter_convention =
        "well-discretization/connection-to-cell-source/sign-bridge/v1";

/// Adapter result between production-positive well rates and the model-neutral
/// injection-positive finite-volume cell-source contract.
///
/// The embedded CellSourceLinearization3D contains reservoir natural-variable
/// derivatives only. Explicit BHP derivatives remain a sidecar because BHP is
/// not a globally numbered reservoir/PETSc unknown in this slice.
struct WellConnectionCellSourceAdapterResult3P {
    static constexpr std::string_view convention =
        well_connection_cell_source_adapter_convention;

    mpmc::flow::NaturalVariableStateIdentity3P
        state_identity;
    double bottom_hole_pressure_pa{};

    mpmc::flow_discretization::
        CellSourceLinearization3D
            cell_source;

    std::vector<double>
        component_source_bhp_derivative_mol_per_pa_s;
    double energy_source_bhp_derivative_w_per_pa{};

    [[nodiscard]] double
    d_component_source_d_bottom_hole_pressure(
        std::size_t component) const {
        return component_source_bhp_derivative_mol_per_pa_s.at(
            component);
    }

    [[nodiscard]] double
    d_energy_source_d_bottom_hole_pressure() const noexcept {
        return energy_source_bhp_derivative_w_per_pa;
    }
};

namespace cell_source_adapter_detail {

[[nodiscard]] inline bool
same_state_identity(
    const mpmc::flow::NaturalVariableStateIdentity3P&
        first,
    const mpmc::flow::NaturalVariableStateIdentity3P&
        second) {
    return
        first.component_ids ==
            second.component_ids &&
        first.layout.component_count() ==
            second.layout.component_count() &&
        first.layout.phase_count() ==
            second.layout.phase_count() &&
        first.layout.unknown_count() ==
            second.layout.unknown_count() &&
        first.layout.composition_pivot()
                .dependent_components() ==
            second.layout.composition_pivot()
                .dependent_components() &&
        first.reference_pressure_pa ==
            second.reference_pressure_pa &&
        first.temperature_k ==
            second.temperature_k &&
        first.saturation ==
            second.saturation &&
        first.phase_composition ==
            second.phase_composition;
}

inline void validate_upstream_results(
    const WellConnectionComponentMolarRateLinearization3P&
        component_rate,
    const WellConnectionAdvectiveEnergyRateLinearization3P&
        energy_rate) {
    const std::size_t n =
        component_rate.component_ids.size();
    const std::size_t q =
        component_rate.input_count;

    if (!same_state_identity(
            component_rate.state_identity,
            energy_rate.state_identity) ||
        component_rate.bottom_hole_pressure_pa !=
            energy_rate.bottom_hole_pressure_pa ||
        !std::isfinite(
            component_rate.bottom_hole_pressure_pa) ||
        !(component_rate.bottom_hole_pressure_pa > 0.0) ||
        component_rate.component_ids !=
            component_rate.state_identity.component_ids ||
        n == 0U ||
        q == 0U ||
        q !=
            component_rate.state_identity.layout
                .unknown_count() ||
        n >
            std::numeric_limits<std::size_t>::max() /
                q ||
        component_rate
                .component_molar_rate_mol_per_s
                .size() !=
            n ||
        component_rate
                .component_reservoir_jacobian
                .size() !=
            n * q ||
        component_rate
                .component_bhp_derivative_mol_per_pa_s
                .size() !=
            n ||
        energy_rate
                .reservoir_natural_variable_gradient
                .size() !=
            q ||
        !std::isfinite(
            energy_rate.total_advective_energy_rate_w) ||
        !std::isfinite(
            energy_rate
                .bottom_hole_pressure_derivative_w_per_pa)) {
        throw std::invalid_argument(
            "mpmc::well_discretization::make_connection_cell_source_adapter_3p: upstream well-rate identity/shape mismatch");
    }

    for (double value :
         component_rate
             .component_molar_rate_mol_per_s) {
        if (!std::isfinite(value)) {
            throw std::invalid_argument(
                "mpmc::well_discretization: well component molar rate must be finite");
        }
    }
    for (double value :
         component_rate
             .component_reservoir_jacobian) {
        if (!std::isfinite(value)) {
            throw std::invalid_argument(
                "mpmc::well_discretization: well component reservoir derivative must be finite");
        }
    }
    for (double value :
         component_rate
             .component_bhp_derivative_mol_per_pa_s) {
        if (!std::isfinite(value)) {
            throw std::invalid_argument(
                "mpmc::well_discretization: well component BHP derivative must be finite");
        }
    }
    for (double value :
         energy_rate
             .reservoir_natural_variable_gradient) {
        if (!std::isfinite(value)) {
            throw std::invalid_argument(
                "mpmc::well_discretization: well energy reservoir derivative must be finite");
        }
    }
}

} // namespace cell_source_adapter_detail

/// Convert one production-positive well connection into the existing
/// injection-positive finite-volume cell-source contract.
///
/// Well convention:
///   n_dot_i^well > 0, E_dot^well > 0 : cell -> well (production)
///
/// CellSource convention:
///   q_i^source > 0, Q_E^source > 0 : injection into cell
///
/// Therefore the adapter applies exactly one sign reversal:
///   q_i^source = -n_dot_i^well
///   Q_E^source = -E_dot^well
///
/// Reservoir Jacobians and explicit BHP derivatives receive the same sign
/// reversal. Downstream CellSource normalization then contributes
///   residual = -source / V_b,
/// so production correctly appears as positive outward residual.
[[nodiscard]] inline
WellConnectionCellSourceAdapterResult3P
make_connection_cell_source_adapter_3p(
    const WellConnectionComponentMolarRateLinearization3P&
        component_rate,
    const WellConnectionAdvectiveEnergyRateLinearization3P&
        energy_rate,
    std::string provenance) {
    using namespace cell_source_adapter_detail;

    validate_upstream_results(
        component_rate,
        energy_rate);
    if (provenance.empty()) {
        throw std::invalid_argument(
            "mpmc::well_discretization::make_connection_cell_source_adapter_3p: source provenance must be nonempty");
    }

    WellConnectionCellSourceAdapterResult3P
        result;
    result.state_identity =
        component_rate.state_identity;
    result.bottom_hole_pressure_pa =
        component_rate.bottom_hole_pressure_pa;

    auto& source =
        result.cell_source;
    source.provenance =
        std::move(provenance);
    source.component_ids =
        component_rate.component_ids;
    source.input_count =
        component_rate.input_count;

    source.component_molar_rate_mol_per_s.resize(
        component_rate.component_count());
    source.component_molar_rate_jacobian_mol_per_s.resize(
        component_rate.component_reservoir_jacobian.size());
    result.component_source_bhp_derivative_mol_per_pa_s.resize(
        component_rate.component_count());

    for (std::size_t component = 0U;
         component <
         component_rate.component_count();
         ++component) {
        source.component_molar_rate_mol_per_s[
            component] =
            -component_rate
                 .component_molar_rate_mol_per_s[
                     component];
        result
            .component_source_bhp_derivative_mol_per_pa_s[
                component] =
            -component_rate
                 .component_bhp_derivative_mol_per_pa_s[
                     component];
    }

    for (std::size_t index = 0U;
         index <
         component_rate
             .component_reservoir_jacobian
             .size();
         ++index) {
        source
            .component_molar_rate_jacobian_mol_per_s[
                index] =
            -component_rate
                 .component_reservoir_jacobian[
                     index];
    }

    source.energy_rate_w =
        -energy_rate.total_advective_energy_rate_w;
    source.energy_rate_gradient_w.resize(
        component_rate.input_count);
    for (std::size_t column = 0U;
         column < component_rate.input_count;
         ++column) {
        source.energy_rate_gradient_w[
            column] =
            -energy_rate
                 .reservoir_natural_variable_gradient[
                     column];
    }
    result.energy_source_bhp_derivative_w_per_pa =
        -energy_rate
             .bottom_hole_pressure_derivative_w_per_pa;

    mpmc::flow_discretization::
        validate_cell_source_linearization(
            source,
            component_rate.component_ids,
            component_rate.input_count);

    for (double value :
         result
             .component_source_bhp_derivative_mol_per_pa_s) {
        if (!std::isfinite(value)) {
            throw std::invalid_argument(
                "mpmc::well_discretization: source BHP derivative must be finite");
        }
    }
    if (!std::isfinite(
            result
                .energy_source_bhp_derivative_w_per_pa)) {
        throw std::invalid_argument(
            "mpmc::well_discretization: energy-source BHP derivative must be finite");
    }

    return result;
}

} // namespace mpmc::well_discretization

#endif // MPMC_WELL_DISCRETIZATION_CELL_SOURCE_ADAPTER_HPP
