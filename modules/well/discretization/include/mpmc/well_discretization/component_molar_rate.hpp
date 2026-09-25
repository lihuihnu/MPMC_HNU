#ifndef MPMC_WELL_DISCRETIZATION_COMPONENT_MOLAR_RATE_HPP
#define MPMC_WELL_DISCRETIZATION_COMPONENT_MOLAR_RATE_HPP

#include <mpmc/flow/component_accumulation.hpp>
#include <mpmc/well_discretization/pressure_drawdown_rate.hpp>

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

namespace mpmc::well_discretization {

inline constexpr std::string_view
    connection_component_molar_rate_convention =
        "well-discretization/connection-component-molar-rate/three-phase-sum/v1";

/// Three-phase connection component molar-rate linearization.
///
/// Production-positive phase volumetric rates from the pressure-drawdown layer
/// are converted to canonical component molar rates:
///
///   n_dot_i = sum_alpha q_alpha * c_alpha * x_alpha,i    [mol/s]
///
/// Reservoir natural-variable derivative:
///
///   d n_dot_i / d q_j = sum_alpha [
///       (d q_alpha/dq_j) c_alpha x_alpha,i
///       + q_alpha (d c_alpha/dq_j) x_alpha,i
///       + q_alpha c_alpha (d x_alpha,i/dq_j)
///   ].
///
/// Explicit BHP derivative:
///
///   d n_dot_i / d p_bhp = sum_alpha
///       (d q_alpha/dp_bhp) c_alpha x_alpha,i.
///
/// Composition derivatives are reconstructed exactly from the frozen
/// natural-variable composition pivot. No extra composition-gradient payload is
/// introduced.
struct WellConnectionComponentMolarRateLinearization3P {
    static constexpr std::string_view convention =
        connection_component_molar_rate_convention;

    mpmc::flow::NaturalVariableStateIdentity3P
        state_identity;
    double bottom_hole_pressure_pa{};

    std::vector<std::string> component_ids;
    std::size_t input_count{};

    std::array<double, 3>
        phase_volumetric_rate_m3_per_s{};
    std::array<double, 3>
        phase_molar_density_mol_per_m3{};
    std::array<double, 3>
        phase_total_molar_rate_mol_per_s{};

    std::vector<double>
        component_molar_rate_mol_per_s;
    std::vector<double>
        component_reservoir_jacobian;
    std::vector<double>
        component_bhp_derivative_mol_per_pa_s;

    double total_molar_rate_mol_per_s{};
    std::vector<double>
        total_reservoir_gradient;
    double total_bhp_derivative_mol_per_pa_s{};

    [[nodiscard]] std::size_t
    component_count() const noexcept {
        return component_ids.size();
    }

    [[nodiscard]] double
    component_rate(
        std::size_t component) const {
        return component_molar_rate_mol_per_s.at(
            component);
    }

    [[nodiscard]] double
    d_component_rate_d_reservoir_unknown(
        std::size_t component,
        std::size_t column) const {
        if (component >= component_count() ||
            column >= input_count ||
            input_count == 0U ||
            component >
                (std::numeric_limits<std::size_t>::max() -
                 column) /
                    input_count) {
            throw std::out_of_range(
                "mpmc::well_discretization: component molar-rate Jacobian index out of range");
        }
        return component_reservoir_jacobian.at(
            component * input_count +
            column);
    }

    [[nodiscard]] double
    d_component_rate_d_bottom_hole_pressure(
        std::size_t component) const {
        return component_bhp_derivative_mol_per_pa_s.at(
            component);
    }
};

namespace component_molar_rate_detail {

[[nodiscard]] inline bool near_roundoff(
    double first,
    double second) {
    if (!std::isfinite(first) ||
        !std::isfinite(second)) {
        return false;
    }
    const double scale =
        std::max(
            {1.0e-30,
             std::abs(first),
             std::abs(second)});
    return std::abs(first - second) <=
        8192.0 *
            std::numeric_limits<double>::epsilon() *
            scale;
}

inline void validate_identity_and_density(
    const mpmc::flow::NaturalVariableStateIdentity3P&
        identity,
    const mpmc::flow::
        PhaseMolarDensityNaturalVariableLinearization3P&
            density) {
    const auto& descriptor =
        identity.layout;
    const std::size_t component_count =
        descriptor.component_count();
    const std::size_t q =
        descriptor.unknown_count();

    if (descriptor.phase_count() !=
            mpmc::flow::fixed_three_phase_count ||
        component_count < 2U ||
        identity.component_ids.size() !=
            component_count ||
        density.layout.component_count() !=
            component_count ||
        density.layout.unknown_count() != q ||
        density.layout.composition_pivot()
                .dependent_components() !=
            descriptor.composition_pivot()
                .dependent_components()) {
        throw std::invalid_argument(
            "mpmc::well_discretization::make_connection_component_molar_rate_linearization_3p: state/density natural-variable chart mismatch");
    }

    for (std::size_t component = 0U;
         component < component_count;
         ++component) {
        if (identity.component_ids[component].empty()) {
            throw std::invalid_argument(
                "mpmc::well_discretization: component identity must be nonempty");
        }
        for (std::size_t previous = 0U;
             previous < component;
             ++previous) {
            if (identity.component_ids[previous] ==
                identity.component_ids[component]) {
                throw std::invalid_argument(
                    "mpmc::well_discretization: component identity/order must be unique");
            }
        }
    }

    for (std::size_t phase = 0U;
         phase < mpmc::flow::fixed_three_phase_count;
         ++phase) {
        const auto& composition =
            identity.phase_composition[phase];
        if (composition.size() !=
            component_count) {
            throw std::invalid_argument(
                "mpmc::well_discretization: phase composition shape does not match canonical component count");
        }

        double composition_sum = 0.0;
        for (double value : composition) {
            if (!std::isfinite(value) ||
                !(value > 0.0)) {
                throw std::invalid_argument(
                    "mpmc::well_discretization: phase composition must retain finite positive support");
            }
            composition_sum += value;
        }
        if (!near_roundoff(
                composition_sum,
                1.0)) {
            throw std::invalid_argument(
                "mpmc::well_discretization: phase composition does not sum to unity");
        }

        const double molar_density =
            density
                .molar_density_mol_per_m3[
                    phase];
        if (!std::isfinite(molar_density) ||
            !(molar_density > 0.0) ||
            density.gradient[phase].size() !=
                q) {
            throw std::invalid_argument(
                "mpmc::well_discretization: phase molar-density value/gradient shape is invalid");
        }
        for (double derivative :
             density.gradient[phase]) {
            if (!std::isfinite(derivative)) {
                throw std::invalid_argument(
                    "mpmc::well_discretization: phase molar-density gradient contains non-finite derivative");
            }
        }
    }
}

[[nodiscard]] inline double
d_phase_composition(
    const mpmc::flow::NaturalVariableLayoutDescriptor&
        layout,
    std::size_t phase,
    std::size_t component,
    std::size_t column) {
    const auto identity =
        layout.composition_unknown_identity(
            column);
    if (!identity ||
        static_cast<std::size_t>(
            identity->phase) !=
            phase) {
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

} // namespace component_molar_rate_detail

/// Build the canonical component molar-rate vector for one well connection by
/// summing all three phase contributions at one explicit BHP.
///
/// This remains a connection-local constitutive/discretization contract. It
/// does not convert the result into a finite-volume cell source and does not
/// create any well-control or PETSc equation.
[[nodiscard]] inline
WellConnectionComponentMolarRateLinearization3P
make_connection_component_molar_rate_linearization_3p(
    const mpmc::well::PeacemanWellIndex3D&
        connection,
    const mpmc::flow::LocalPhaseMobilityLinearization3P&
        mobility,
    const mpmc::flow::
        PhaseMolarDensityNaturalVariableLinearization3P&
            molar_density,
    double bottom_hole_pressure_pa) {
    using namespace component_molar_rate_detail;

    validate_identity_and_density(
        mobility.state_identity,
        molar_density);

    const auto& identity =
        mobility.state_identity;
    const auto& layout =
        identity.layout;
    const std::size_t component_count =
        layout.component_count();
    const std::size_t q =
        layout.unknown_count();

    if (component_count >
        std::numeric_limits<std::size_t>::max() /
            q) {
        throw std::length_error(
            "mpmc::well_discretization: component molar-rate Jacobian size overflow");
    }

    WellConnectionComponentMolarRateLinearization3P
        result;
    result.state_identity =
        identity;
    result.bottom_hole_pressure_pa =
        bottom_hole_pressure_pa;
    result.component_ids =
        identity.component_ids;
    result.input_count =
        q;
    result.phase_molar_density_mol_per_m3 =
        molar_density.molar_density_mol_per_m3;
    result.component_molar_rate_mol_per_s.assign(
        component_count,
        0.0);
    result.component_reservoir_jacobian.assign(
        component_count * q,
        0.0);
    result.component_bhp_derivative_mol_per_pa_s.assign(
        component_count,
        0.0);
    result.total_reservoir_gradient.assign(
        q,
        0.0);

    for (std::size_t phase = 0U;
         phase < mpmc::flow::fixed_three_phase_count;
         ++phase) {
        const auto slot =
            static_cast<
                mpmc::flow::PhaseSlot3>(
                    phase);
        const auto phase_rate =
            make_pressure_drawdown_phase_rate_linearization_3p(
                connection,
                mobility,
                slot,
                bottom_hole_pressure_pa);

        const double volumetric_rate =
            phase_rate
                .phase_volumetric_rate_m3_per_s;
        const double density =
            molar_density
                .molar_density_mol_per_m3[
                    phase];
        const double phase_molar_rate =
            volumetric_rate *
            density;
        if (!std::isfinite(phase_molar_rate)) {
            throw std::range_error(
                "mpmc::well_discretization: phase total molar rate is outside representable range");
        }

        result.phase_volumetric_rate_m3_per_s[
            phase] =
            volumetric_rate;
        result.phase_total_molar_rate_mol_per_s[
            phase] =
            phase_molar_rate;
        result.total_molar_rate_mol_per_s +=
            phase_molar_rate;

        const auto& composition =
            identity.phase_composition[
                phase];

        for (std::size_t component = 0U;
             component < component_count;
             ++component) {
            const double x =
                composition[component];
            const double contribution =
                phase_molar_rate *
                x;
            if (!std::isfinite(contribution)) {
                throw std::range_error(
                    "mpmc::well_discretization: component molar-rate contribution is outside representable range");
            }
            result.component_molar_rate_mol_per_s[
                component] +=
                contribution;

            const double bhp_derivative =
                phase_rate
                    .d_phase_rate_d_bottom_hole_pressure() *
                density *
                x;
            if (!std::isfinite(bhp_derivative)) {
                throw std::range_error(
                    "mpmc::well_discretization: component BHP derivative is outside representable range");
            }
            result.component_bhp_derivative_mol_per_pa_s[
                component] +=
                bhp_derivative;
        }

        const double phase_total_bhp_derivative =
            phase_rate
                .d_phase_rate_d_bottom_hole_pressure() *
            density;
        if (!std::isfinite(
                phase_total_bhp_derivative)) {
            throw std::range_error(
                "mpmc::well_discretization: total BHP derivative is outside representable range");
        }
        result.total_bhp_derivative_mol_per_pa_s +=
            phase_total_bhp_derivative;

        for (std::size_t column = 0U;
             column < q;
             ++column) {
            const double dq =
                phase_rate
                    .d_phase_rate_d_reservoir_unknown(
                        column);
            const double dc =
                molar_density.gradient[
                    phase][column];

            const double phase_total_derivative =
                dq * density +
                volumetric_rate * dc;
            if (!std::isfinite(
                    phase_total_derivative)) {
                throw std::range_error(
                    "mpmc::well_discretization: total molar-rate derivative is outside representable range");
            }
            result.total_reservoir_gradient[
                column] +=
                phase_total_derivative;

            for (std::size_t component = 0U;
                 component < component_count;
                 ++component) {
                const double x =
                    composition[component];
                const double dx =
                    d_phase_composition(
                        layout,
                        phase,
                        component,
                        column);
                const double derivative =
                    dq *
                        density *
                        x +
                    volumetric_rate *
                        dc *
                        x +
                    volumetric_rate *
                        density *
                        dx;
                if (!std::isfinite(derivative)) {
                    throw std::range_error(
                        "mpmc::well_discretization: component reservoir derivative is outside representable range");
                }
                result.component_reservoir_jacobian[
                    component * q +
                    column] +=
                    derivative;
            }
        }
    }

    double component_rate_sum = 0.0;
    double component_bhp_sum = 0.0;
    for (std::size_t component = 0U;
         component < component_count;
         ++component) {
        component_rate_sum +=
            result.component_molar_rate_mol_per_s[
                component];
        component_bhp_sum +=
            result.component_bhp_derivative_mol_per_pa_s[
                component];
    }
    if (!near_roundoff(
            component_rate_sum,
            result.total_molar_rate_mol_per_s) ||
        !near_roundoff(
            component_bhp_sum,
            result.total_bhp_derivative_mol_per_pa_s)) {
        throw std::runtime_error(
            "mpmc::well_discretization: component molar rates/BHP derivatives do not close to phase totals");
    }

    for (std::size_t column = 0U;
         column < q;
         ++column) {
        double component_gradient_sum =
            0.0;
        for (std::size_t component = 0U;
             component < component_count;
             ++component) {
            component_gradient_sum +=
                result.component_reservoir_jacobian[
                    component * q +
                    column];
        }
        if (!near_roundoff(
                component_gradient_sum,
                result.total_reservoir_gradient[
                    column])) {
            throw std::runtime_error(
                "mpmc::well_discretization: differentiated component molar rates do not close to total molar rate");
        }
    }

    return result;
}

} // namespace mpmc::well_discretization

#endif // MPMC_WELL_DISCRETIZATION_COMPONENT_MOLAR_RATE_HPP
