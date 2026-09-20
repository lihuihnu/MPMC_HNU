#ifndef MPMC_FLOW_DISCRETIZATION_TPFA_COMPONENT_MOLAR_FLUX_HPP
#define MPMC_FLOW_DISCRETIZATION_TPFA_COMPONENT_MOLAR_FLUX_HPP

#include <mpmc/flow/component_accumulation.hpp>
#include <mpmc/flow_discretization/tpfa_phase_darcy_flux.hpp>

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

namespace mpmc::flow_discretization {

inline constexpr std::string_view
    tpfa_component_molar_flux_convention =
        "flow_discretization/upwind-phase-molar-content-component-flux/v1";

struct UpwindPhaseComponentMolarContent3P {
    mpmc::flow::UpwindCellSelection3P upwind_selection{
        mpmc::flow::UpwindCellSelection3P::
            owner_exact_zero_tie};

    double upwind_molar_density_mol_per_m3{};
    std::vector<double> upwind_mole_fraction;
    std::vector<double>
        component_molar_content_mol_per_m3;

    std::vector<double>
        owner_component_molar_content_jacobian;
    std::vector<double>
        neighbour_component_molar_content_jacobian;
};

/// Component molar flux on one already-materialized TPFA internal face.
///
/// For canonical component i:
///
///   n_dot_i^f = sum_alpha [
///       c_alpha,up * x_alpha,i,up * F_alpha
///   ].
///
/// F_alpha is the already-validated phase volumetric Darcy flux [m^3/s].
/// c*x is [mol/m^3], so the component flux is [mol/s].
///
/// Positive component flux follows the phase-flux convention:
/// owner -> neighbour.
struct MaterializedTpfaInternalFaceComponentMolarFluxLinearization3D {
    static constexpr std::string_view convention =
        tpfa_component_molar_flux_convention;

    mpmc::mesh::LocalIndex face{
        mpmc::mesh::LocalIndex::value_type{0}};

    std::vector<std::string> component_ids;

    mpmc::flow::NaturalVariableStateIdentity3P
        owner_state_identity;
    mpmc::flow::NaturalVariableStateIdentity3P
        neighbour_state_identity;

    std::array<UpwindPhaseComponentMolarContent3P, 3>
        phase_molar_content;

    std::vector<double>
        component_molar_flux_mol_per_s;
    std::vector<double>
        owner_component_flux_jacobian;
    std::vector<double>
        neighbour_component_flux_jacobian;

    double total_molar_flux_mol_per_s{};
    std::vector<double>
        owner_total_molar_flux_gradient;
    std::vector<double>
        neighbour_total_molar_flux_gradient;

    [[nodiscard]] std::size_t
    component_count() const noexcept {
        return component_ids.size();
    }

    [[nodiscard]] double component_flux(
        std::size_t component) const {
        return component_molar_flux_mol_per_s.at(
            component);
    }

    [[nodiscard]] double d_component_flux_owner(
        std::size_t component,
        std::size_t column) const {
        const std::size_t q =
            owner_state_identity.layout.unknown_count();
        if (component >= component_ids.size() ||
            column >= q ||
            (q != 0U &&
             component >
                 (std::numeric_limits<std::size_t>::max() -
                  column) /
                     q)) {
            throw std::out_of_range(
                "mpmc::flow_discretization: owner component-flux Jacobian index out of range");
        }
        return owner_component_flux_jacobian.at(
            component * q + column);
    }

    [[nodiscard]] double d_component_flux_neighbour(
        std::size_t component,
        std::size_t column) const {
        const std::size_t q =
            neighbour_state_identity.layout.unknown_count();
        if (component >= component_ids.size() ||
            column >= q ||
            (q != 0U &&
             component >
                 (std::numeric_limits<std::size_t>::max() -
                  column) /
                     q)) {
            throw std::out_of_range(
                "mpmc::flow_discretization: neighbour component-flux Jacobian index out of range");
        }
        return neighbour_component_flux_jacobian.at(
            component * q + column);
    }
};

namespace component_flux_detail {

[[nodiscard]] inline bool near_roundoff(
    double first,
    double second,
    double extra_scale = 0.0) {
    if (!std::isfinite(first) ||
        !std::isfinite(second) ||
        !std::isfinite(extra_scale) ||
        extra_scale < 0.0) {
        return false;
    }
    const double scale =
        std::max(
            {1.0,
             std::abs(first),
             std::abs(second),
             extra_scale});
    return std::abs(first - second) <=
        8192.0 *
            std::numeric_limits<double>::epsilon() *
            scale;
}

inline void validate_state_identity(
    const mpmc::flow::NaturalVariableCellState3P& state,
    const mpmc::flow::NaturalVariableStateIdentity3P&
        identity,
    const char* side) {
    const auto& layout =
        state.layout();

    if (identity.component_ids !=
            std::vector<std::string>{
                state.component_ids().begin(),
                state.component_ids().end()} ||
        identity.layout.component_count() !=
            layout.component_count() ||
        identity.layout.phase_count() !=
            mpmc::flow::fixed_three_phase_count ||
        identity.layout.composition_pivot()
                .dependent_components() !=
            std::vector<std::size_t>{
                layout.composition_pivot()
                    .dependent_components()
                    .begin(),
                layout.composition_pivot()
                    .dependent_components()
                    .end()} ||
        identity.layout.unknown_count() !=
            layout.unknown_count() ||
        identity.saturation.size() !=
            mpmc::flow::fixed_three_phase_count ||
        identity.phase_composition.size() !=
            mpmc::flow::fixed_three_phase_count ||
        !near_roundoff(
            identity.reference_pressure_pa,
            state.reference_pressure_pa()) ||
        !near_roundoff(
            identity.temperature_k,
            state.temperature_k())) {
        throw std::invalid_argument(
            std::string{"mpmc::flow_discretization: "} +
            side +
            " phase-flux state identity/chart does not match cell state");
    }

    for (std::size_t phase = 0U;
         phase < mpmc::flow::fixed_three_phase_count;
         ++phase) {
        const auto slot =
            static_cast<mpmc::flow::PhaseSlot3>(
                phase);
        if (!near_roundoff(
                identity.saturation[phase],
                state.phase_saturation(slot))) {
            throw std::invalid_argument(
                std::string{"mpmc::flow_discretization: "} +
                side +
                " phase-flux saturation identity does not match cell state");
        }

        const auto composition =
            state.phase_composition(slot);
        if (identity.phase_composition[phase].size() !=
            composition.size()) {
            throw std::invalid_argument(
                std::string{"mpmc::flow_discretization: "} +
                side +
                " phase-flux composition shape does not match cell state");
        }
        for (std::size_t component = 0U;
             component < composition.size();
             ++component) {
            if (!near_roundoff(
                    identity.phase_composition[phase][component],
                    composition[component])) {
                throw std::invalid_argument(
                    std::string{"mpmc::flow_discretization: "} +
                    side +
                    " phase-flux composition identity does not match cell state");
            }
        }
    }
}

inline void validate_density_linearization(
    const mpmc::flow::NaturalVariableCellState3P& state,
    const mpmc::flow::
        PhaseMolarDensityNaturalVariableLinearization3P&
            density,
    const char* side) {
    const auto& layout =
        state.layout();
    if (density.layout.component_count() !=
            layout.component_count() ||
        density.layout.composition_pivot()
                .dependent_components() !=
            layout.composition_pivot()
                .dependent_components() ||
        density.layout.unknown_count() !=
            layout.unknown_count()) {
        throw std::invalid_argument(
            std::string{"mpmc::flow_discretization: "} +
            side +
            " molar-density linearization chart mismatch");
    }

    for (std::size_t phase = 0U;
         phase < mpmc::flow::fixed_three_phase_count;
         ++phase) {
        const double expected =
            state.phase_properties(
                     static_cast<mpmc::flow::PhaseSlot3>(
                         phase))
                .molar_density_mol_per_m3;
        const double actual =
            density
                .molar_density_mol_per_m3[phase];

        if (!std::isfinite(actual) ||
            !(actual > 0.0) ||
            !near_roundoff(actual, expected)) {
            throw std::invalid_argument(
                std::string{"mpmc::flow_discretization: "} +
                side +
                " molar-density primal does not match exact cell state");
        }

        if (density.gradient[phase].size() !=
            layout.unknown_count()) {
            throw std::invalid_argument(
                std::string{"mpmc::flow_discretization: "} +
                side +
                " molar-density gradient shape mismatch");
        }
        for (double derivative :
             density.gradient[phase]) {
            if (!std::isfinite(derivative)) {
                throw std::invalid_argument(
                    std::string{"mpmc::flow_discretization: "} +
                    side +
                    " molar-density gradient contains non-finite derivative");
            }
        }
    }
}

inline void validate_phase_flux(
    const MaterializedTpfaInternalFacePhaseDarcyFluxLinearization3D&
        phase_flux) {
    const std::size_t owner_q =
        phase_flux.owner_state_identity.layout
            .unknown_count();
    const std::size_t neighbour_q =
        phase_flux.neighbour_state_identity.layout
            .unknown_count();

    if (!std::isfinite(
            phase_flux.static_transmissibility_m3) ||
        !(phase_flux.static_transmissibility_m3 > 0.0) ||
        owner_q == 0U ||
        neighbour_q == 0U ||
        phase_flux.owner_state_identity.component_ids !=
            phase_flux.neighbour_state_identity.component_ids) {
        throw std::invalid_argument(
            "mpmc::flow_discretization: malformed materialized phase-flux state identity");
    }

    for (const auto& phase :
         phase_flux.phase) {
        if (!std::isfinite(
                phase.volumetric_flux_m3_per_s) ||
            !std::isfinite(
                phase.phase_potential_difference_pa) ||
            !std::isfinite(
                phase.upwind_mobility_per_pa_s) ||
            phase.upwind_mobility_per_pa_s < 0.0 ||
            phase.owner_flux_gradient.size() !=
                owner_q ||
            phase.neighbour_flux_gradient.size() !=
                neighbour_q) {
            throw std::invalid_argument(
                "mpmc::flow_discretization: malformed phase Darcy-flux payload");
        }

        for (double derivative :
             phase.owner_flux_gradient) {
            if (!std::isfinite(derivative)) {
                throw std::invalid_argument(
                    "mpmc::flow_discretization: phase Darcy-flux Jacobian contains non-finite derivative");
            }
        }
        for (double derivative :
             phase.neighbour_flux_gradient) {
            if (!std::isfinite(derivative)) {
                throw std::invalid_argument(
                    "mpmc::flow_discretization: phase Darcy-flux Jacobian contains non-finite derivative");
            }
        }

        const double expected_flux =
            -phase_flux.static_transmissibility_m3 *
            phase.upwind_mobility_per_pa_s *
            phase.phase_potential_difference_pa;
        if (!near_roundoff(
                phase.volumetric_flux_m3_per_s,
                expected_flux,
                std::abs(expected_flux))) {
            throw std::invalid_argument(
                "mpmc::flow_discretization: phase Darcy-flux primal is inconsistent with T_f, upwind mobility and phase potential");
        }

        switch (phase.upwind_selection) {
        case mpmc::flow::UpwindCellSelection3P::
            owner_negative_phase_potential:
            if (!(phase.phase_potential_difference_pa < 0.0)) {
                throw std::invalid_argument(
                    "mpmc::flow_discretization: owner-upwind phase must have negative phase potential");
            }
            break;
        case mpmc::flow::UpwindCellSelection3P::
            neighbour_positive_phase_potential:
            if (!(phase.phase_potential_difference_pa > 0.0)) {
                throw std::invalid_argument(
                    "mpmc::flow_discretization: neighbour-upwind phase must have positive phase potential");
            }
            break;
        case mpmc::flow::UpwindCellSelection3P::
            owner_exact_zero_tie:
            if (phase.phase_potential_difference_pa != 0.0 ||
                phase.volumetric_flux_m3_per_s != 0.0) {
                throw std::invalid_argument(
                    "mpmc::flow_discretization: exact-zero upwind tie must have zero phase potential and zero phase flux");
            }
            break;
        default:
            throw std::invalid_argument(
                "mpmc::flow_discretization: invalid phase upwind selection");
        }
    }
}

[[nodiscard]] inline bool selects_owner(
    mpmc::flow::UpwindCellSelection3P selection) {
    return selection ==
               mpmc::flow::UpwindCellSelection3P::
                   owner_negative_phase_potential ||
           selection ==
               mpmc::flow::UpwindCellSelection3P::
                   owner_exact_zero_tie;
}

[[nodiscard]] inline double d_composition(
    const mpmc::flow::NaturalVariableLayout3P& layout,
    std::size_t phase,
    std::size_t component,
    std::size_t column) {
    const auto identity =
        layout.composition_unknown_identity(column);
    if (!identity ||
        static_cast<std::size_t>(
            identity->phase) != phase) {
        return 0.0;
    }
    if (component == identity->component) {
        return 1.0;
    }
    if (component ==
        layout.dependent_composition_component(
            static_cast<mpmc::flow::PhaseSlot3>(
                phase))) {
        return -1.0;
    }
    return 0.0;
}

} // namespace component_flux_detail

[[nodiscard]] inline
MaterializedTpfaInternalFaceComponentMolarFluxLinearization3D
build_materialized_tpfa_internal_face_component_molar_flux(
    const MaterializedTpfaInternalFacePhaseDarcyFluxLinearization3D&
        phase_flux,
    const mpmc::flow::NaturalVariableCellState3P&
        owner_state,
    const mpmc::flow::
        PhaseMolarDensityNaturalVariableLinearization3P&
            owner_molar_density,
    const mpmc::flow::NaturalVariableCellState3P&
        neighbour_state,
    const mpmc::flow::
        PhaseMolarDensityNaturalVariableLinearization3P&
            neighbour_molar_density) {
    using namespace component_flux_detail;

    validate_phase_flux(phase_flux);
    validate_state_identity(
        owner_state,
        phase_flux.owner_state_identity,
        "owner");
    validate_state_identity(
        neighbour_state,
        phase_flux.neighbour_state_identity,
        "neighbour");
    validate_density_linearization(
        owner_state,
        owner_molar_density,
        "owner");
    validate_density_linearization(
        neighbour_state,
        neighbour_molar_density,
        "neighbour");

    if (owner_state.component_ids().size() !=
            neighbour_state.component_ids().size() ||
        !std::equal(
            owner_state.component_ids().begin(),
            owner_state.component_ids().end(),
            neighbour_state.component_ids().begin())) {
        throw std::invalid_argument(
            "mpmc::flow_discretization: owner/neighbour component identity/order mismatch");
    }

    const std::size_t n =
        owner_state.component_ids().size();
    const std::size_t owner_q =
        owner_state.layout().unknown_count();
    const std::size_t neighbour_q =
        neighbour_state.layout().unknown_count();

    if (n == 0U ||
        n >
            std::numeric_limits<std::size_t>::max() /
                owner_q ||
        n >
            std::numeric_limits<std::size_t>::max() /
                neighbour_q) {
        throw std::length_error(
            "mpmc::flow_discretization: component-flux Jacobian size overflow");
    }

    MaterializedTpfaInternalFaceComponentMolarFluxLinearization3D
        result{
            phase_flux.face,
            std::vector<std::string>{
                owner_state.component_ids().begin(),
                owner_state.component_ids().end()},
            phase_flux.owner_state_identity,
            phase_flux.neighbour_state_identity,
            {},
            std::vector<double>(n, 0.0),
            std::vector<double>(
                n * owner_q,
                0.0),
            std::vector<double>(
                n * neighbour_q,
                0.0),
            0.0,
            std::vector<double>(owner_q, 0.0),
            std::vector<double>(neighbour_q, 0.0)};

    for (std::size_t phase = 0U;
         phase < mpmc::flow::fixed_three_phase_count;
         ++phase) {
        const auto& phase_flux_entry =
            phase_flux.phase[phase];
        auto& content =
            result.phase_molar_content[phase];

        content.upwind_selection =
            phase_flux_entry.upwind_selection;

        const bool owner_upwind =
            selects_owner(
                phase_flux_entry.upwind_selection);
        const auto slot =
            static_cast<mpmc::flow::PhaseSlot3>(
                phase);

        const auto& upstream_state =
            owner_upwind
                ? owner_state
                : neighbour_state;
        const auto& upstream_density =
            owner_upwind
                ? owner_molar_density
                : neighbour_molar_density;

        content.upwind_molar_density_mol_per_m3 =
            upstream_density
                .molar_density_mol_per_m3[phase];
        const auto x =
            upstream_state.phase_composition(slot);
        content.upwind_mole_fraction.assign(
            x.begin(),
            x.end());
        content
            .component_molar_content_mol_per_m3
            .assign(n, 0.0);
        content
            .owner_component_molar_content_jacobian
            .assign(
                n * owner_q,
                0.0);
        content
            .neighbour_component_molar_content_jacobian
            .assign(
                n * neighbour_q,
                0.0);

        double phase_content_sum = 0.0;
        for (std::size_t component = 0U;
             component < n;
             ++component) {
            const double molar_content =
                content
                    .upwind_molar_density_mol_per_m3 *
                x[component];
            if (!std::isfinite(molar_content) ||
                !(molar_content > 0.0)) {
                throw std::range_error(
                    "mpmc::flow_discretization: upwind component molar content is non-finite or non-positive");
            }
            content
                .component_molar_content_mol_per_m3
                [component] =
                molar_content;
            phase_content_sum +=
                molar_content;

            if (owner_upwind) {
                for (std::size_t column = 0U;
                     column < owner_q;
                     ++column) {
                    const double dx =
                        d_composition(
                            owner_state.layout(),
                            phase,
                            component,
                            column);
                    const double derivative =
                        x[component] *
                            owner_molar_density
                                .gradient[phase][column] +
                        content
                            .upwind_molar_density_mol_per_m3 *
                            dx;
                    if (!std::isfinite(derivative)) {
                        throw std::range_error(
                            "mpmc::flow_discretization: owner upwind molar-content derivative is non-finite");
                    }
                    content
                        .owner_component_molar_content_jacobian
                        [component * owner_q + column] =
                        derivative;
                }
            } else {
                for (std::size_t column = 0U;
                     column < neighbour_q;
                     ++column) {
                    const double dx =
                        d_composition(
                            neighbour_state.layout(),
                            phase,
                            component,
                            column);
                    const double derivative =
                        x[component] *
                            neighbour_molar_density
                                .gradient[phase][column] +
                        content
                            .upwind_molar_density_mol_per_m3 *
                            dx;
                    if (!std::isfinite(derivative)) {
                        throw std::range_error(
                            "mpmc::flow_discretization: neighbour upwind molar-content derivative is non-finite");
                    }
                    content
                        .neighbour_component_molar_content_jacobian
                        [component * neighbour_q + column] =
                        derivative;
                }
            }
        }

        if (!near_roundoff(
                phase_content_sum,
                content
                    .upwind_molar_density_mol_per_m3,
                content
                    .upwind_molar_density_mol_per_m3)) {
            throw std::runtime_error(
                "mpmc::flow_discretization: upwind component molar contents do not close to phase molar density");
        }

        const auto& selected_density_gradient =
            owner_upwind
                ? owner_molar_density.gradient[phase]
                : neighbour_molar_density.gradient[phase];
        const auto& selected_content_jacobian =
            owner_upwind
                ? content
                      .owner_component_molar_content_jacobian
                : content
                      .neighbour_component_molar_content_jacobian;
        const std::size_t selected_q =
            owner_upwind
                ? owner_q
                : neighbour_q;

        for (std::size_t column = 0U;
             column < selected_q;
             ++column) {
            double derivative_sum = 0.0;
            double derivative_scale =
                std::abs(
                    selected_density_gradient[column]);
            for (std::size_t component = 0U;
                 component < n;
                 ++component) {
                const double derivative =
                    selected_content_jacobian[
                        component * selected_q +
                        column];
                derivative_sum += derivative;
                derivative_scale +=
                    std::abs(derivative);
            }
            if (!near_roundoff(
                    derivative_sum,
                    selected_density_gradient[column],
                    derivative_scale)) {
                throw std::runtime_error(
                    "mpmc::flow_discretization: differentiated upwind molar contents do not close to phase molar-density derivative");
            }
        }
    }

    for (std::size_t component = 0U;
         component < n;
         ++component) {
        double flux = 0.0;

        for (std::size_t phase = 0U;
             phase < mpmc::flow::fixed_three_phase_count;
             ++phase) {
            const auto& phase_flux_entry =
                phase_flux.phase[phase];
            const auto& content =
                result.phase_molar_content[phase];
            const double molar_content =
                content
                    .component_molar_content_mol_per_m3
                    [component];

            flux +=
                molar_content *
                phase_flux_entry
                    .volumetric_flux_m3_per_s;

            for (std::size_t column = 0U;
                 column < owner_q;
                 ++column) {
                result
                    .owner_component_flux_jacobian
                    [component * owner_q +
                     column] +=
                    molar_content *
                        phase_flux_entry
                            .owner_flux_gradient[column] +
                    phase_flux_entry
                            .volumetric_flux_m3_per_s *
                        content
                            .owner_component_molar_content_jacobian
                            [component * owner_q +
                             column];
            }

            for (std::size_t column = 0U;
                 column < neighbour_q;
                 ++column) {
                result
                    .neighbour_component_flux_jacobian
                    [component * neighbour_q +
                     column] +=
                    molar_content *
                        phase_flux_entry
                            .neighbour_flux_gradient[column] +
                    phase_flux_entry
                            .volumetric_flux_m3_per_s *
                        content
                            .neighbour_component_molar_content_jacobian
                            [component * neighbour_q +
                             column];
            }
        }

        if (!std::isfinite(flux)) {
            throw std::range_error(
                "mpmc::flow_discretization: component molar face flux is non-finite");
        }
        result.component_molar_flux_mol_per_s[
            component] =
            flux;
        result.total_molar_flux_mol_per_s +=
            flux;
    }

    double expected_total = 0.0;
    for (std::size_t phase = 0U;
         phase < mpmc::flow::fixed_three_phase_count;
         ++phase) {
        expected_total +=
            result.phase_molar_content[phase]
                .upwind_molar_density_mol_per_m3 *
            phase_flux.phase[phase]
                .volumetric_flux_m3_per_s;
    }
    if (!near_roundoff(
            result.total_molar_flux_mol_per_s,
            expected_total,
            std::abs(expected_total))) {
        throw std::runtime_error(
            "mpmc::flow_discretization: component molar face fluxes do not close to total molar face flux");
    }

    for (std::size_t column = 0U;
         column < owner_q;
         ++column) {
        double component_sum = 0.0;
        double expected = 0.0;
        double scale = 0.0;
        for (std::size_t component = 0U;
             component < n;
             ++component) {
            const double derivative =
                result.owner_component_flux_jacobian[
                    component * owner_q +
                    column];
            component_sum += derivative;
            scale += std::abs(derivative);
        }

        for (std::size_t phase = 0U;
             phase < mpmc::flow::fixed_three_phase_count;
             ++phase) {
            const auto& phase_flux_entry =
                phase_flux.phase[phase];
            const auto& content =
                result.phase_molar_content[phase];
            double dc_up = 0.0;
            if (selects_owner(
                    phase_flux_entry.upwind_selection)) {
                dc_up =
                    owner_molar_density
                        .gradient[phase][column];
            }
            expected +=
                content
                    .upwind_molar_density_mol_per_m3 *
                    phase_flux_entry
                        .owner_flux_gradient[column] +
                phase_flux_entry
                        .volumetric_flux_m3_per_s *
                    dc_up;
        }
        result.owner_total_molar_flux_gradient[
            column] =
            expected;
        if (!near_roundoff(
                component_sum,
                expected,
                scale + std::abs(expected))) {
            throw std::runtime_error(
                "mpmc::flow_discretization: owner component-flux Jacobian does not close to total molar-flux derivative");
        }
    }

    for (std::size_t column = 0U;
         column < neighbour_q;
         ++column) {
        double component_sum = 0.0;
        double expected = 0.0;
        double scale = 0.0;
        for (std::size_t component = 0U;
             component < n;
             ++component) {
            const double derivative =
                result.neighbour_component_flux_jacobian[
                    component * neighbour_q +
                    column];
            component_sum += derivative;
            scale += std::abs(derivative);
        }

        for (std::size_t phase = 0U;
             phase < mpmc::flow::fixed_three_phase_count;
             ++phase) {
            const auto& phase_flux_entry =
                phase_flux.phase[phase];
            const auto& content =
                result.phase_molar_content[phase];
            double dc_up = 0.0;
            if (!selects_owner(
                    phase_flux_entry.upwind_selection)) {
                dc_up =
                    neighbour_molar_density
                        .gradient[phase][column];
            }
            expected +=
                content
                    .upwind_molar_density_mol_per_m3 *
                    phase_flux_entry
                        .neighbour_flux_gradient[column] +
                phase_flux_entry
                        .volumetric_flux_m3_per_s *
                    dc_up;
        }
        result.neighbour_total_molar_flux_gradient[
            column] =
            expected;
        if (!near_roundoff(
                component_sum,
                expected,
                scale + std::abs(expected))) {
            throw std::runtime_error(
                "mpmc::flow_discretization: neighbour component-flux Jacobian does not close to total molar-flux derivative");
        }
    }

    return result;
}

} // namespace mpmc::flow_discretization

#endif // MPMC_FLOW_DISCRETIZATION_TPFA_COMPONENT_MOLAR_FLUX_HPP
