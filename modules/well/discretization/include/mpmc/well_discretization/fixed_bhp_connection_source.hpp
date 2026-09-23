#ifndef MPMC_WELL_DISCRETIZATION_FIXED_BHP_CONNECTION_SOURCE_HPP
#define MPMC_WELL_DISCRETIZATION_FIXED_BHP_CONNECTION_SOURCE_HPP

#include <mpmc/flow/phase_transport.hpp>
#include <mpmc/flow_discretization/cell_source.hpp>
#include <mpmc/well/peaceman_well_index_3d.hpp>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <span>
#include <stdexcept>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

namespace mpmc::well_discretization {

inline constexpr std::string_view
    fixed_bhp_variable_cardinality_connection_source_convention =
        "well-discretization/fixed-bhp-single-connection/active-phase-sum/v1";

/// Well-side injection enthalpy for the currently frozen active phase chart.
///
/// The vector length must equal the active phase count. Inactive phases are not
/// padded or assigned fictitious properties.
struct FixedBhpInjectionEnthalpy3D {
    std::string provenance;
    std::vector<double>
        specific_enthalpy_j_per_kg;
};

/// Cardinality-neutral local input required by one fixed-BHP Peaceman
/// connection. All phase arrays contain exactly layout.phase_count() entries.
struct FixedBhpActivePhaseSourceInput3D {
    mpmc::flow::NaturalVariableStateIdentity
        state_identity;

    std::vector<double>
        phase_pressure_pa;
    std::vector<std::vector<double>>
        phase_pressure_gradient;

    std::vector<double>
        mobility_per_pa_s;
    std::vector<std::vector<double>>
        mobility_gradient;

    std::vector<double>
        molar_density_mol_per_m3;
    std::vector<std::vector<double>>
        molar_density_gradient;

    std::vector<double>
        mass_density_kg_per_m3;
    std::vector<std::vector<double>>
        mass_density_gradient;

    std::vector<double>
        specific_enthalpy_j_per_kg;
    std::vector<std::vector<double>>
        specific_enthalpy_gradient;
};

/// Injection-positive cell-source result plus explicit BHP derivatives.
///
/// BHP remains external to the reservoir natural-variable numbering in this
/// contract. The CellSource payload therefore contains reservoir derivatives
/// only; BHP derivatives remain a sidecar for future well-unknown work.
struct FixedBhpConnectionCellSourceLinearization3D {
    static constexpr std::string_view convention =
        fixed_bhp_variable_cardinality_connection_source_convention;

    mpmc::flow::NaturalVariableStateIdentity
        state_identity;
    double bottom_hole_pressure_pa{};

    mpmc::flow_discretization::
        CellSourceLinearization3D
            cell_source;

    std::vector<double>
        component_source_bhp_derivative_mol_per_pa_s;
    double energy_source_bhp_derivative_w_per_pa{};

    std::vector<double>
        phase_volumetric_rate_m3_per_s;
};

namespace fixed_bhp_connection_source_detail {

[[nodiscard]] inline bool
near_roundoff(
    double first,
    double second) {
    if (!std::isfinite(first) ||
        !std::isfinite(second)) {
        return false;
    }
    const double scale =
        std::max(
            {1.0,
             std::abs(first),
             std::abs(second)});
    return std::abs(first - second) <=
        8192.0 *
            std::numeric_limits<double>::epsilon() *
            scale;
}

inline void validate_gradient(
    const std::vector<double>& gradient,
    std::size_t q,
    std::string_view name) {
    if (gradient.size() != q) {
        throw std::invalid_argument(
            "mpmc::well_discretization: fixed-BHP " +
            std::string{name} +
            " gradient shape mismatch");
    }
    for (double value : gradient) {
        if (!std::isfinite(value)) {
            throw std::invalid_argument(
                "mpmc::well_discretization: fixed-BHP " +
                std::string{name} +
                " gradient contains non-finite value");
        }
    }
}

inline void validate_input(
    const mpmc::well::PeacemanWellIndex3D&
        connection,
    const FixedBhpActivePhaseSourceInput3D&
        input,
    const FixedBhpInjectionEnthalpy3D&
        injection,
    double bottom_hole_pressure_pa) {
    const auto& layout =
        input.state_identity.layout;
    const std::size_t p =
        layout.phase_count();
    const std::size_t n =
        layout.component_count();
    const std::size_t q =
        layout.unknown_count();

    if (p == 0U || p > 3U ||
        n < 2U ||
        q == 0U ||
        input.state_identity.component_ids.size() != n ||
        input.phase_pressure_pa.size() != p ||
        input.phase_pressure_gradient.size() != p ||
        input.mobility_per_pa_s.size() != p ||
        input.mobility_gradient.size() != p ||
        input.molar_density_mol_per_m3.size() != p ||
        input.molar_density_gradient.size() != p ||
        input.mass_density_kg_per_m3.size() != p ||
        input.mass_density_gradient.size() != p ||
        input.specific_enthalpy_j_per_kg.size() != p ||
        input.specific_enthalpy_gradient.size() != p ||
        injection.provenance.empty() ||
        injection.specific_enthalpy_j_per_kg.size() != p ||
        !std::isfinite(bottom_hole_pressure_pa) ||
        !(bottom_hole_pressure_pa > 0.0) ||
        !std::isfinite(connection.well_index_m3) ||
        !(connection.well_index_m3 > 0.0)) {
        throw std::invalid_argument(
            "mpmc::well_discretization: invalid fixed-BHP variable-cardinality connection input");
    }

    for (std::size_t component = 0U;
         component < n;
         ++component) {
        if (input.state_identity
                .component_ids[component]
                .empty()) {
            throw std::invalid_argument(
                "mpmc::well_discretization: fixed-BHP component identity must be nonempty");
        }
        for (std::size_t previous = 0U;
             previous < component;
             ++previous) {
            if (input.state_identity
                    .component_ids[previous] ==
                input.state_identity
                    .component_ids[component]) {
                throw std::invalid_argument(
                    "mpmc::well_discretization: fixed-BHP component identity/order must be unique");
            }
        }
    }

    double saturation_sum = 0.0;
    for (std::size_t phase = 0U;
         phase < p;
         ++phase) {
        const double saturation =
            input.state_identity
                .saturation[phase];
        if (!std::isfinite(saturation) ||
            saturation < 0.0) {
            throw std::invalid_argument(
                "mpmc::well_discretization: fixed-BHP active saturation is invalid");
        }
        saturation_sum += saturation;

        const auto& composition =
            input.state_identity
                .phase_composition[phase];
        if (composition.size() != n) {
            throw std::invalid_argument(
                "mpmc::well_discretization: fixed-BHP phase composition shape mismatch");
        }
        double composition_sum = 0.0;
        for (double value : composition) {
            if (!std::isfinite(value) ||
                !(value > 0.0)) {
                throw std::invalid_argument(
                    "mpmc::well_discretization: fixed-BHP phase composition must retain positive support");
            }
            composition_sum += value;
        }
        if (!near_roundoff(
                composition_sum,
                1.0)) {
            throw std::invalid_argument(
                "mpmc::well_discretization: fixed-BHP phase composition does not sum to unity");
        }

        const double phase_pressure =
            input.phase_pressure_pa[phase];
        const double mobility =
            input.mobility_per_pa_s[phase];
        const double molar_density =
            input.molar_density_mol_per_m3[phase];
        const double mass_density =
            input.mass_density_kg_per_m3[phase];
        const double enthalpy =
            input.specific_enthalpy_j_per_kg[phase];
        const double injection_enthalpy =
            injection
                .specific_enthalpy_j_per_kg[
                    phase];

        if (!std::isfinite(phase_pressure) ||
            !(phase_pressure > 0.0) ||
            !std::isfinite(mobility) ||
            mobility < 0.0 ||
            !std::isfinite(molar_density) ||
            !(molar_density > 0.0) ||
            !std::isfinite(mass_density) ||
            !(mass_density > 0.0) ||
            !std::isfinite(enthalpy) ||
            !std::isfinite(injection_enthalpy)) {
            throw std::invalid_argument(
                "mpmc::well_discretization: fixed-BHP active phase property is invalid");
        }

        validate_gradient(
            input.phase_pressure_gradient[
                phase],
            q,
            "phase-pressure");
        validate_gradient(
            input.mobility_gradient[
                phase],
            q,
            "mobility");
        validate_gradient(
            input.molar_density_gradient[
                phase],
            q,
            "molar-density");
        validate_gradient(
            input.mass_density_gradient[
                phase],
            q,
            "mass-density");
        validate_gradient(
            input.specific_enthalpy_gradient[
                phase],
            q,
            "enthalpy");
    }

    if (!near_roundoff(
            saturation_sum,
            1.0)) {
        throw std::invalid_argument(
            "mpmc::well_discretization: fixed-BHP active saturations do not sum to unity");
    }
}

[[nodiscard]] inline double
phase_composition_derivative(
    const mpmc::flow::
        NaturalVariableLayoutDescriptor&
            layout,
    std::size_t phase,
    std::size_t component,
    std::size_t column) {
    const auto identity =
        layout.composition_unknown_identity(
            column);
    if (!identity ||
        static_cast<std::size_t>(
            identity->phase) != phase) {
        return 0.0;
    }

    if (component ==
        identity->component) {
        return 1.0;
    }

    if (component ==
        layout.dependent_composition_component(
            static_cast<
                mpmc::flow::PhaseSlot3>(
                    phase))) {
        return -1.0;
    }
    return 0.0;
}

} // namespace fixed_bhp_connection_source_detail

/// Build the complete fixed-BHP connection source for one frozen 1P/2P/3P
/// natural-variable chart.
///
/// Production-positive phase rate:
///   q_alpha = WI * lambda_alpha * (p_alpha - p_bhp)
///
/// Component production:
///   n_dot_i = sum_alpha q_alpha c_alpha x_alpha,i
///
/// Advective energy production:
///   E_dot = sum_alpha q_alpha rho_alpha h_selected
///
/// where production and exact-zero ties use reservoir enthalpy while injection
/// uses the explicit well-side enthalpy. The returned CellSource follows the
/// repository injection-positive convention, so all production rates and
/// reservoir derivatives receive exactly one sign reversal.
[[nodiscard]] inline
FixedBhpConnectionCellSourceLinearization3D
make_fixed_bhp_connection_cell_source_3d(
    const mpmc::well::PeacemanWellIndex3D&
        connection,
    const FixedBhpActivePhaseSourceInput3D&
        input,
    const FixedBhpInjectionEnthalpy3D&
        injection_enthalpy,
    double bottom_hole_pressure_pa,
    std::string source_provenance) {
    using namespace
        fixed_bhp_connection_source_detail;

    validate_input(
        connection,
        input,
        injection_enthalpy,
        bottom_hole_pressure_pa);
    if (source_provenance.empty()) {
        throw std::invalid_argument(
            "mpmc::well_discretization: fixed-BHP source provenance must be nonempty");
    }

    const auto& layout =
        input.state_identity.layout;
    const std::size_t p =
        layout.phase_count();
    const std::size_t n =
        layout.component_count();
    const std::size_t q =
        layout.unknown_count();

    std::vector<double>
        phase_rate(
            p,
            0.0);
    std::vector<std::vector<double>>
        phase_rate_gradient(
            p,
            std::vector<double>(
                q,
                0.0));
    std::vector<double>
        phase_bhp_derivative(
            p,
            0.0);

    for (std::size_t phase = 0U;
         phase < p;
         ++phase) {
        const double conductance =
            connection.well_index_m3 *
            input.mobility_per_pa_s[
                phase];
        const double drawdown =
            input.phase_pressure_pa[
                phase] -
            bottom_hole_pressure_pa;
        phase_rate[phase] =
            conductance * drawdown;
        phase_bhp_derivative[phase] =
            -conductance;

        if (!std::isfinite(conductance) ||
            !std::isfinite(drawdown) ||
            !std::isfinite(
                phase_rate[phase]) ||
            !std::isfinite(
                phase_bhp_derivative[phase])) {
            throw std::range_error(
                "mpmc::well_discretization: fixed-BHP phase rate is outside representable range");
        }

        for (std::size_t column = 0U;
             column < q;
             ++column) {
            phase_rate_gradient[phase][
                column] =
                connection.well_index_m3 *
                    input.mobility_gradient[
                        phase][column] *
                    drawdown +
                conductance *
                    input.phase_pressure_gradient[
                        phase][column];
            if (!std::isfinite(
                    phase_rate_gradient[
                        phase][column])) {
                throw std::range_error(
                    "mpmc::well_discretization: fixed-BHP phase-rate derivative is outside representable range");
            }
        }
    }

    std::vector<double>
        component_production_rate(
            n,
            0.0);
    std::vector<double>
        component_production_gradient(
            n * q,
            0.0);
    std::vector<double>
        component_production_bhp_derivative(
            n,
            0.0);

    for (std::size_t phase = 0U;
         phase < p;
         ++phase) {
        const double density =
            input.molar_density_mol_per_m3[
                phase];
        const auto& composition =
            input.state_identity
                .phase_composition[phase];

        for (std::size_t component = 0U;
             component < n;
             ++component) {
            const double x =
                composition[component];
            component_production_rate[
                component] +=
                phase_rate[phase] *
                density *
                x;
            component_production_bhp_derivative[
                component] +=
                phase_bhp_derivative[
                    phase] *
                density *
                x;

            for (std::size_t column = 0U;
                 column < q;
                 ++column) {
                const double dx =
                    phase_composition_derivative(
                        layout,
                        phase,
                        component,
                        column);
                component_production_gradient[
                    component * q +
                    column] +=
                    phase_rate_gradient[
                        phase][column] *
                        density *
                        x +
                    phase_rate[phase] *
                        input.molar_density_gradient[
                            phase][column] *
                        x +
                    phase_rate[phase] *
                        density *
                        dx;
            }
        }
    }

    double energy_production_rate = 0.0;
    std::vector<double>
        energy_production_gradient(
            q,
            0.0);
    double
        energy_production_bhp_derivative =
            0.0;

    for (std::size_t phase = 0U;
         phase < p;
         ++phase) {
        const bool reservoir_enthalpy =
            phase_rate[phase] >= 0.0;
        const double selected_enthalpy =
            reservoir_enthalpy
                ? input
                      .specific_enthalpy_j_per_kg[
                          phase]
                : injection_enthalpy
                      .specific_enthalpy_j_per_kg[
                          phase];
        const double mass_density =
            input.mass_density_kg_per_m3[
                phase];

        energy_production_rate +=
            phase_rate[phase] *
            mass_density *
            selected_enthalpy;
        energy_production_bhp_derivative +=
            phase_bhp_derivative[phase] *
            mass_density *
            selected_enthalpy;

        for (std::size_t column = 0U;
             column < q;
             ++column) {
            const double dh =
                reservoir_enthalpy
                    ? input
                          .specific_enthalpy_gradient[
                              phase][column]
                    : 0.0;
            energy_production_gradient[
                column] +=
                phase_rate_gradient[
                    phase][column] *
                    mass_density *
                    selected_enthalpy +
                phase_rate[phase] *
                    input.mass_density_gradient[
                        phase][column] *
                    selected_enthalpy +
                phase_rate[phase] *
                    mass_density *
                    dh;
        }
    }

    FixedBhpConnectionCellSourceLinearization3D
        result;
    result.state_identity =
        input.state_identity;
    result.bottom_hole_pressure_pa =
        bottom_hole_pressure_pa;
    result.phase_volumetric_rate_m3_per_s =
        phase_rate;

    auto& source =
        result.cell_source;
    source.provenance =
        std::move(source_provenance);
    source.component_ids =
        input.state_identity.component_ids;
    source.input_count = q;
    source.component_molar_rate_mol_per_s.resize(
        n);
    source.component_molar_rate_jacobian_mol_per_s.resize(
        n * q);
    result.component_source_bhp_derivative_mol_per_pa_s.resize(
        n);

    for (std::size_t component = 0U;
         component < n;
         ++component) {
        source.component_molar_rate_mol_per_s[
            component] =
            -component_production_rate[
                 component];
        result
            .component_source_bhp_derivative_mol_per_pa_s[
                component] =
            -component_production_bhp_derivative[
                 component];
        if (!std::isfinite(
                source
                    .component_molar_rate_mol_per_s[
                        component]) ||
            !std::isfinite(
                result
                    .component_source_bhp_derivative_mol_per_pa_s[
                        component])) {
            throw std::range_error(
                "mpmc::well_discretization: fixed-BHP component source is outside representable range");
        }
    }

    for (std::size_t index = 0U;
         index < n * q;
         ++index) {
        source.component_molar_rate_jacobian_mol_per_s[
            index] =
            -component_production_gradient[
                 index];
        if (!std::isfinite(
                source
                    .component_molar_rate_jacobian_mol_per_s[
                        index])) {
            throw std::range_error(
                "mpmc::well_discretization: fixed-BHP component source derivative is outside representable range");
        }
    }

    source.energy_rate_w =
        -energy_production_rate;
    source.energy_rate_gradient_w.resize(
        q);
    for (std::size_t column = 0U;
         column < q;
         ++column) {
        source.energy_rate_gradient_w[
            column] =
            -energy_production_gradient[
                 column];
        if (!std::isfinite(
                source
                    .energy_rate_gradient_w[
                        column])) {
            throw std::range_error(
                "mpmc::well_discretization: fixed-BHP energy source derivative is outside representable range");
        }
    }
    result.energy_source_bhp_derivative_w_per_pa =
        -energy_production_bhp_derivative;

    if (!std::isfinite(
            source.energy_rate_w) ||
        !std::isfinite(
            result
                .energy_source_bhp_derivative_w_per_pa)) {
        throw std::range_error(
            "mpmc::well_discretization: fixed-BHP energy source is outside representable range");
    }

    mpmc::flow_discretization::
        validate_cell_source_linearization(
            source,
            source.component_ids,
            q);
    return result;
}

/// Production-positive primal total for one logical fixed-BHP well.
///
/// This aggregate deliberately contains no reservoir Jacobian: connections may
/// live on different cells and therefore have different natural-variable
/// columns. Each connection source remains assembled into its own cell block.
/// The well total is only the authoritative sum of connection-level primal
/// rates at one evaluated state.
struct FixedBhpWellProductionRateAggregation3D {
    double bottom_hole_pressure_pa{};
    std::vector<std::string>
        component_ids;
    std::vector<double>
        component_molar_rate_mol_per_s;
    double energy_rate_w{};
    std::size_t authoritative_connection_count{};
};

/// Sum already-authoritative connection sources into one logical well total.
///
/// CellSourceLinearization3D is injection-positive, while this returned
/// aggregate is production-positive, so each primal source receives exactly
/// one sign reversal. Callers are responsible for supplying each physical
/// connection exactly once; distributed callers should aggregate only
/// owner-side connection evaluations before their MPI reduction.
[[nodiscard]] inline
FixedBhpWellProductionRateAggregation3D
aggregate_fixed_bhp_connection_production_rates_3d(
    std::span<
        const FixedBhpConnectionCellSourceLinearization3D>
        authoritative_connections) {
    if (authoritative_connections.empty()) {
        throw std::invalid_argument(
            "mpmc::well_discretization: fixed-BHP well aggregation requires at least one authoritative connection");
    }

    const auto& first =
        authoritative_connections.front();
    FixedBhpWellProductionRateAggregation3D
        result;
    result.bottom_hole_pressure_pa =
        first.bottom_hole_pressure_pa;
    result.component_ids =
        first.cell_source.component_ids;
    result.component_molar_rate_mol_per_s.assign(
        result.component_ids.size(),
        0.0);
    result.authoritative_connection_count =
        authoritative_connections.size();

    if (!std::isfinite(
            result.bottom_hole_pressure_pa) ||
        !(result.bottom_hole_pressure_pa > 0.0) ||
        result.component_ids.empty()) {
        throw std::invalid_argument(
            "mpmc::well_discretization: malformed fixed-BHP authoritative connection aggregate");
    }

    for (const auto& connection :
         authoritative_connections) {
        if (connection.bottom_hole_pressure_pa !=
                result.bottom_hole_pressure_pa ||
            connection.cell_source.component_ids !=
                result.component_ids ||
            connection
                    .cell_source
                    .component_molar_rate_mol_per_s
                    .size() !=
                result.component_ids.size() ||
            !std::isfinite(
                connection
                    .cell_source
                    .energy_rate_w)) {
            throw std::invalid_argument(
                "mpmc::well_discretization: fixed-BHP authoritative connections do not share BHP/component ordering");
        }

        for (std::size_t component = 0U;
             component < result.component_ids.size();
             ++component) {
            const double source_rate =
                connection
                    .cell_source
                    .component_molar_rate_mol_per_s[
                        component];
            if (!std::isfinite(source_rate)) {
                throw std::invalid_argument(
                    "mpmc::well_discretization: fixed-BHP authoritative connection contains non-finite component rate");
            }
            result.component_molar_rate_mol_per_s[
                component] -=
                source_rate;
            if (!std::isfinite(
                    result
                        .component_molar_rate_mol_per_s[
                            component])) {
                throw std::range_error(
                    "mpmc::well_discretization: fixed-BHP aggregated component production rate is outside representable range");
            }
        }

        result.energy_rate_w -=
            connection.cell_source.energy_rate_w;
        if (!std::isfinite(
                result.energy_rate_w)) {
            throw std::range_error(
                "mpmc::well_discretization: fixed-BHP aggregated energy production rate is outside representable range");
        }
    }
    return result;
}

} // namespace mpmc::well_discretization

#endif // MPMC_WELL_DISCRETIZATION_FIXED_BHP_CONNECTION_SOURCE_HPP
