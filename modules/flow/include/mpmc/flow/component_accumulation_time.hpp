#ifndef MPMC_FLOW_COMPONENT_ACCUMULATION_TIME_HPP
#define MPMC_FLOW_COMPONENT_ACCUMULATION_TIME_HPP

#include <mpmc/flow/component_accumulation.hpp>

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

namespace mpmc::flow {

inline constexpr std::string_view
    backward_euler_component_accumulation_convention =
        "flow/component-accumulation/backward-euler-local/v1";

/// Local backward-Euler accumulation residual and current-state Jacobian.
///
/// Residual units:
///   mol / (bulk-m^3 s)
///
/// The Jacobian differentiates only with respect to the current
/// phase-cardinality-specific natural-variable coordinates. The previous accumulation snapshot and
/// dt are frozen history/discretization inputs, not nonlinear unknowns.
struct BackwardEulerComponentAccumulationResidual3P {
    static constexpr std::string_view convention =
        backward_euler_component_accumulation_convention;

    NaturalVariableLayoutDescriptor current_layout{
        std::size_t{2U},
        std::size_t{1U},
        std::vector<std::size_t>{1U}};
    double porosity{};
    double time_step_seconds{};
    std::vector<std::string> component_ids;
    std::vector<double>
        component_residual_mol_per_bulk_m3_s;
    double total_residual_mol_per_bulk_m3_s{};
    std::size_t input_count{};
    std::vector<double> component_jacobian;
    std::vector<double> total_residual_gradient;

    [[nodiscard]] std::size_t
    component_count() const noexcept {
        return component_ids.size();
    }

    [[nodiscard]] double residual(
        std::size_t component) const {
        return component_residual_mol_per_bulk_m3_s.at(
            component);
    }

    [[nodiscard]] std::string_view
    component_row_identity(
        std::size_t component) const {
        return component_ids.at(component);
    }

    [[nodiscard]] double d_residual(
        std::size_t component,
        std::size_t column) const {
        if (component >= component_ids.size() ||
            column >= input_count ||
            (input_count != 0U &&
             component >
                 (std::numeric_limits<std::size_t>::max() -
                  column) /
                     input_count)) {
            throw std::out_of_range(
                "mpmc::flow: backward-Euler accumulation Jacobian index out of range");
        }
        return component_jacobian.at(
            component * input_count +
            column);
    }

    [[nodiscard]] double d_total(
        std::size_t column) const {
        return total_residual_gradient.at(column);
    }
};

namespace component_accumulation_time_detail {

[[nodiscard]] inline bool near_scaled_roundoff(
    double first,
    double second,
    double scale) {
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
        4096.0 *
            std::numeric_limits<double>::
                epsilon() *
            reference;
}

inline void validate_time_step(
    double time_step_seconds) {
    if (!std::isfinite(time_step_seconds) ||
        !(time_step_seconds > 0.0)) {
        throw std::invalid_argument(
            "mpmc::flow: backward-Euler time step must be finite and strictly positive [s]");
    }
}

inline void validate_current_linearization(
    const PoreVolumeComponentAccumulationSnapshot3P&
        current,
    const PoreVolumeComponentAccumulationLinearization3P&
        linearization) {
    const std::size_t n =
        current.component_ids.size();
    const std::size_t q =
        linearization.layout.unknown_count();

    if (linearization.component_ids !=
            current.component_ids ||
        linearization.layout.component_count() != n ||
        linearization.input_count != q ||
        !component_accumulation_detail::near_roundoff(
            linearization.porosity,
            current.porosity)) {
        throw std::invalid_argument(
            "mpmc::flow: current accumulation linearization identity/porosity does not match current snapshot");
    }
    if (q == 0U ||
        n >
            std::numeric_limits<std::size_t>::max() /
                q ||
        linearization.component_jacobian.size() !=
            n * q ||
        linearization.total_accumulation_gradient.size() !=
            q) {
        throw std::invalid_argument(
            "mpmc::flow: current accumulation linearization shape is invalid");
    }

    for (std::size_t column = 0U;
         column < q;
         ++column) {
        const double total =
            linearization.total_accumulation_gradient[
                column];
        if (!std::isfinite(total)) {
            throw std::invalid_argument(
                "mpmc::flow: current accumulation linearization contains non-finite derivative");
        }

        double component_sum = 0.0;
        for (std::size_t component = 0U;
             component < n;
             ++component) {
            const double derivative =
                linearization.component_jacobian[
                    component * q + column];
            if (!std::isfinite(derivative)) {
                throw std::invalid_argument(
                    "mpmc::flow: current accumulation linearization contains non-finite derivative");
            }
            component_sum += derivative;
        }
        double derivative_scale =
            std::abs(total);
        for (std::size_t component = 0U;
             component < n;
             ++component) {
            derivative_scale +=
                std::abs(
                    linearization.component_jacobian[
                        component * q + column]);
        }
        if (!near_scaled_roundoff(
                component_sum,
                total,
                derivative_scale)) {
            throw std::invalid_argument(
                "mpmc::flow: current accumulation linearization violates differentiated component closure");
        }
    }
}

} // namespace component_accumulation_time_detail

/// Build only the local backward-Euler accumulation contribution:
///
///   R_i^acc = (N_i^{n+1} - N_i^n) / dt
///
/// and its current-state Jacobian:
///
///   d R_i^acc / d q^{n+1}
///       = (1/dt) d N_i^{n+1} / d q^{n+1}.
///
/// No spatial flux, source/well, mobility, Darcy velocity, gravity, global
/// residual assembly or nonlinear solve is implied.
[[nodiscard]] inline
BackwardEulerComponentAccumulationResidual3P
build_backward_euler_component_accumulation_residual(
    const PoreVolumeComponentAccumulationPair3P&
        accumulation,
    const PoreVolumeComponentAccumulationLinearization3P&
        current_linearization,
    double time_step_seconds) {
    component_accumulation_time_detail::
        validate_time_step(time_step_seconds);
    component_accumulation_detail::
        validate_snapshot_identity(
            accumulation.current);
    component_accumulation_detail::
        validate_snapshot_identity(
            accumulation.previous);

    if (accumulation.current.component_ids !=
        accumulation.previous.component_ids) {
        throw std::invalid_argument(
            "mpmc::flow: backward-Euler current/previous component identity/order mismatch");
    }
    if (!component_accumulation_detail::near_roundoff(
            accumulation.current.porosity,
            accumulation.previous.porosity)) {
        throw std::invalid_argument(
            "mpmc::flow: backward-Euler rigid-medium current/previous porosity mismatch");
    }

    component_accumulation_time_detail::
        validate_current_linearization(
            accumulation.current,
            current_linearization);

    const std::size_t n =
        accumulation.current.component_ids.size();
    const std::size_t q =
        current_linearization.input_count;
    const double inverse_dt =
        1.0 / time_step_seconds;
    if (!std::isfinite(inverse_dt)) {
        throw std::range_error(
            "mpmc::flow: backward-Euler inverse time step is non-finite");
    }

    BackwardEulerComponentAccumulationResidual3P
        result{
            current_linearization.layout,
            accumulation.current.porosity,
            time_step_seconds,
            accumulation.current.component_ids,
            std::vector<double>(n, 0.0),
            0.0,
            q,
            std::vector<double>(n * q, 0.0),
            std::vector<double>(q, 0.0)};

    double component_residual_sum = 0.0;
    for (std::size_t component = 0U;
         component < n;
         ++component) {
        const double value =
            (accumulation.current
                 .component_accumulation_mol_per_bulk_m3[
                     component] -
             accumulation.previous
                 .component_accumulation_mol_per_bulk_m3[
                     component]) *
            inverse_dt;
        if (!std::isfinite(value)) {
            throw std::range_error(
                "mpmc::flow: backward-Euler component accumulation residual is non-finite");
        }
        result
            .component_residual_mol_per_bulk_m3_s[
                component] =
            value;
        component_residual_sum += value;
    }

    result.total_residual_mol_per_bulk_m3_s =
        (accumulation.current
             .total_accumulation_mol_per_bulk_m3 -
         accumulation.previous
             .total_accumulation_mol_per_bulk_m3) *
        inverse_dt;
    double inventory_difference_scale =
        std::abs(
            accumulation.current
                .total_accumulation_mol_per_bulk_m3) +
        std::abs(
            accumulation.previous
                .total_accumulation_mol_per_bulk_m3);
    for (std::size_t component = 0U;
         component < n;
         ++component) {
        inventory_difference_scale +=
            std::abs(
                accumulation.current
                    .component_accumulation_mol_per_bulk_m3[
                        component]) +
            std::abs(
                accumulation.previous
                    .component_accumulation_mol_per_bulk_m3[
                        component]);
    }
    inventory_difference_scale *=
        inverse_dt;

    if (!std::isfinite(
            result.total_residual_mol_per_bulk_m3_s) ||
        !component_accumulation_time_detail::
            near_scaled_roundoff(
                component_residual_sum,
                result.total_residual_mol_per_bulk_m3_s,
                inventory_difference_scale)) {
        throw std::runtime_error(
            "mpmc::flow: backward-Euler component residuals do not close to total accumulation residual");
    }

    for (std::size_t column = 0U;
         column < q;
         ++column) {
        const double total_derivative =
            current_linearization
                .total_accumulation_gradient[column] *
            inverse_dt;
        if (!std::isfinite(total_derivative)) {
            throw std::range_error(
                "mpmc::flow: backward-Euler total accumulation derivative is non-finite");
        }
        result.total_residual_gradient[column] =
            total_derivative;

        double component_derivative_sum = 0.0;
        for (std::size_t component = 0U;
             component < n;
             ++component) {
            const double derivative =
                current_linearization
                    .component_jacobian[
                        component * q + column] *
                inverse_dt;
            if (!std::isfinite(derivative)) {
                throw std::range_error(
                    "mpmc::flow: backward-Euler component accumulation derivative is non-finite");
            }
            result.component_jacobian[
                component * q + column] =
                derivative;
            component_derivative_sum += derivative;
        }

        double derivative_scale =
            std::abs(total_derivative);
        for (std::size_t component = 0U;
             component < n;
             ++component) {
            derivative_scale +=
                std::abs(
                    result.component_jacobian[
                        component * q + column]);
        }
        if (!component_accumulation_time_detail::
                near_scaled_roundoff(
                    component_derivative_sum,
                    total_derivative,
                    derivative_scale)) {
            throw std::runtime_error(
                "mpmc::flow: backward-Euler differentiated component residuals do not close to total residual");
        }
    }

    return result;
}

} // namespace mpmc::flow

#endif // MPMC_FLOW_COMPONENT_ACCUMULATION_TIME_HPP
