#ifndef MPMC_FLOW_DISCRETIZATION_CONSERVATIVE_COMPONENT_FACE_RATE_SCATTER_HPP
#define MPMC_FLOW_DISCRETIZATION_CONSERVATIVE_COMPONENT_FACE_RATE_SCATTER_HPP

#include <mpmc/flow_discretization/tpfa_component_molar_flux.hpp>

#include <algorithm>
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
    conservative_component_face_rate_scatter_convention =
        "flow_discretization/two-cell-conservative-component-face-rate-scatter/v1";

/// Conservative two-cell scatter of one canonical component face rate.
///
/// Positive source face rate is owner -> neighbour. Under the conservation
/// convention
///
///   accumulation + net outward face rate = source,
///
/// the same internal-face rate contributes
///
///   R_owner^face     = +n_dot^f
///   R_neighbour^face = -n_dot^f.
///
/// Units remain [mol/s]. This is not yet a complete cell residual because no
/// cell bulk volume normalization or accumulation term is combined here.
struct ConservativeComponentFaceRateScatterLinearization3D {
    static constexpr std::string_view convention =
        conservative_component_face_rate_scatter_convention;

    mpmc::mesh::LocalIndex face{
        mpmc::mesh::LocalIndex::value_type{0}};

    std::vector<std::string> component_ids;

    mpmc::flow::NaturalVariableStateIdentity3P
        owner_state_identity;
    mpmc::flow::NaturalVariableStateIdentity3P
        neighbour_state_identity;

    std::vector<double>
        owner_component_face_rate_mol_per_s;
    std::vector<double>
        neighbour_component_face_rate_mol_per_s;

    /// d R_owner^face / d q_owner
    std::vector<double>
        owner_row_owner_column_jacobian;
    /// d R_owner^face / d q_neighbour
    std::vector<double>
        owner_row_neighbour_column_jacobian;
    /// d R_neighbour^face / d q_owner
    std::vector<double>
        neighbour_row_owner_column_jacobian;
    /// d R_neighbour^face / d q_neighbour
    std::vector<double>
        neighbour_row_neighbour_column_jacobian;

    double owner_total_face_rate_mol_per_s{};
    double neighbour_total_face_rate_mol_per_s{};

    std::vector<double>
        owner_total_row_owner_column_gradient;
    std::vector<double>
        owner_total_row_neighbour_column_gradient;
    std::vector<double>
        neighbour_total_row_owner_column_gradient;
    std::vector<double>
        neighbour_total_row_neighbour_column_gradient;

    [[nodiscard]] std::size_t
    component_count() const noexcept {
        return component_ids.size();
    }

    [[nodiscard]] double owner_rate(
        std::size_t component) const {
        return owner_component_face_rate_mol_per_s.at(
            component);
    }

    [[nodiscard]] double neighbour_rate(
        std::size_t component) const {
        return neighbour_component_face_rate_mol_per_s.at(
            component);
    }

    [[nodiscard]] double
    d_owner_rate_wrt_owner(
        std::size_t component,
        std::size_t column) const {
        return jacobian_entry(
            owner_row_owner_column_jacobian,
            owner_state_identity.layout.unknown_count(),
            component,
            column,
            "owner-row/owner-column");
    }

    [[nodiscard]] double
    d_owner_rate_wrt_neighbour(
        std::size_t component,
        std::size_t column) const {
        return jacobian_entry(
            owner_row_neighbour_column_jacobian,
            neighbour_state_identity.layout.unknown_count(),
            component,
            column,
            "owner-row/neighbour-column");
    }

    [[nodiscard]] double
    d_neighbour_rate_wrt_owner(
        std::size_t component,
        std::size_t column) const {
        return jacobian_entry(
            neighbour_row_owner_column_jacobian,
            owner_state_identity.layout.unknown_count(),
            component,
            column,
            "neighbour-row/owner-column");
    }

    [[nodiscard]] double
    d_neighbour_rate_wrt_neighbour(
        std::size_t component,
        std::size_t column) const {
        return jacobian_entry(
            neighbour_row_neighbour_column_jacobian,
            neighbour_state_identity.layout.unknown_count(),
            component,
            column,
            "neighbour-row/neighbour-column");
    }

private:
    [[nodiscard]] double jacobian_entry(
        const std::vector<double>& values,
        std::size_t input_count,
        std::size_t component,
        std::size_t column,
        const char* block) const {
        if (component >= component_ids.size() ||
            column >= input_count ||
            input_count == 0U ||
            component >
                (std::numeric_limits<std::size_t>::max() -
                 column) /
                    input_count) {
            throw std::out_of_range(
                std::string{
                    "mpmc::flow_discretization: "} +
                block +
                " face-rate Jacobian index out of range");
        }
        return values.at(
            component * input_count +
            column);
    }
};

namespace conservative_scatter_detail {

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

inline void validate_source(
    const MaterializedTpfaInternalFaceComponentMolarFluxLinearization3D&
        source) {
    const std::size_t n =
        source.component_ids.size();
    const std::size_t owner_q =
        source.owner_state_identity.layout
            .unknown_count();
    const std::size_t neighbour_q =
        source.neighbour_state_identity.layout
            .unknown_count();

    if (n < 2U ||
        source.owner_state_identity.component_ids !=
            source.component_ids ||
        source.neighbour_state_identity.component_ids !=
            source.component_ids ||
        source.owner_state_identity.layout
                .component_count() !=
            n ||
        source.neighbour_state_identity.layout
                .component_count() !=
            n ||
        owner_q == 0U ||
        neighbour_q == 0U) {
        throw std::invalid_argument(
            "mpmc::flow_discretization: component face-flux identity/chart is invalid for conservative scatter");
    }

    for (std::size_t component = 0U;
         component < n;
         ++component) {
        if (source.component_ids[component].empty()) {
            throw std::invalid_argument(
                "mpmc::flow_discretization: component face-flux ID is empty");
        }
        for (std::size_t previous = 0U;
             previous < component;
             ++previous) {
            if (source.component_ids[previous] ==
                source.component_ids[component]) {
                throw std::invalid_argument(
                    "mpmc::flow_discretization: component face-flux IDs must be unique and ordered");
            }
        }
    }

    if (source.component_molar_flux_mol_per_s.size() !=
            n ||
        n >
            std::numeric_limits<std::size_t>::max() /
                owner_q ||
        n >
            std::numeric_limits<std::size_t>::max() /
                neighbour_q ||
        source.owner_component_flux_jacobian.size() !=
            n * owner_q ||
        source.neighbour_component_flux_jacobian.size() !=
            n * neighbour_q ||
        source.owner_total_molar_flux_gradient.size() !=
            owner_q ||
        source.neighbour_total_molar_flux_gradient.size() !=
            neighbour_q ||
        !std::isfinite(
            source.total_molar_flux_mol_per_s)) {
        throw std::invalid_argument(
            "mpmc::flow_discretization: component face-flux value/Jacobian shape is invalid for conservative scatter");
    }

    double component_sum = 0.0;
    double component_scale =
        std::abs(
            source.total_molar_flux_mol_per_s);
    for (double value :
         source.component_molar_flux_mol_per_s) {
        if (!std::isfinite(value)) {
            throw std::invalid_argument(
                "mpmc::flow_discretization: component face flux contains non-finite value");
        }
        component_sum += value;
        component_scale +=
            std::abs(value);
    }
    if (!near_roundoff(
            component_sum,
            source.total_molar_flux_mol_per_s,
            component_scale)) {
        throw std::invalid_argument(
            "mpmc::flow_discretization: component face flux does not close to total molar face flux");
    }

    for (std::size_t column = 0U;
         column < owner_q;
         ++column) {
        double sum = 0.0;
        double scale =
            std::abs(
                source
                    .owner_total_molar_flux_gradient[
                        column]);
        for (std::size_t component = 0U;
             component < n;
             ++component) {
            const double derivative =
                source.owner_component_flux_jacobian[
                    component * owner_q +
                    column];
            if (!std::isfinite(derivative)) {
                throw std::invalid_argument(
                    "mpmc::flow_discretization: owner component face-flux Jacobian contains non-finite derivative");
            }
            sum += derivative;
            scale += std::abs(derivative);
        }
        if (!std::isfinite(
                source
                    .owner_total_molar_flux_gradient[
                        column]) ||
            !near_roundoff(
                sum,
                source
                    .owner_total_molar_flux_gradient[
                        column],
                scale)) {
            throw std::invalid_argument(
                "mpmc::flow_discretization: owner component face-flux Jacobian does not close to total molar-flux derivative");
        }
    }

    for (std::size_t column = 0U;
         column < neighbour_q;
         ++column) {
        double sum = 0.0;
        double scale =
            std::abs(
                source
                    .neighbour_total_molar_flux_gradient[
                        column]);
        for (std::size_t component = 0U;
             component < n;
             ++component) {
            const double derivative =
                source.neighbour_component_flux_jacobian[
                    component * neighbour_q +
                    column];
            if (!std::isfinite(derivative)) {
                throw std::invalid_argument(
                    "mpmc::flow_discretization: neighbour component face-flux Jacobian contains non-finite derivative");
            }
            sum += derivative;
            scale += std::abs(derivative);
        }
        if (!std::isfinite(
                source
                    .neighbour_total_molar_flux_gradient[
                        column]) ||
            !near_roundoff(
                sum,
                source
                    .neighbour_total_molar_flux_gradient[
                        column],
                scale)) {
            throw std::invalid_argument(
                "mpmc::flow_discretization: neighbour component face-flux Jacobian does not close to total molar-flux derivative");
        }
    }
}

} // namespace conservative_scatter_detail

[[nodiscard]] inline
ConservativeComponentFaceRateScatterLinearization3D
scatter_component_molar_face_flux_conservatively(
    const MaterializedTpfaInternalFaceComponentMolarFluxLinearization3D&
        source) {
    conservative_scatter_detail::
        validate_source(source);

    const std::size_t n =
        source.component_ids.size();
    const std::size_t owner_q =
        source.owner_state_identity.layout
            .unknown_count();
    const std::size_t neighbour_q =
        source.neighbour_state_identity.layout
            .unknown_count();

    ConservativeComponentFaceRateScatterLinearization3D
        result{
            source.face,
            source.component_ids,
            source.owner_state_identity,
            source.neighbour_state_identity,
            source.component_molar_flux_mol_per_s,
            {},
            source.owner_component_flux_jacobian,
            source.neighbour_component_flux_jacobian,
            {},
            {},
            source.total_molar_flux_mol_per_s,
            -source.total_molar_flux_mol_per_s,
            source.owner_total_molar_flux_gradient,
            source.neighbour_total_molar_flux_gradient,
            {},
            {}};

    result.neighbour_component_face_rate_mol_per_s
        .resize(n);
    result.neighbour_row_owner_column_jacobian
        .resize(n * owner_q);
    result.neighbour_row_neighbour_column_jacobian
        .resize(n * neighbour_q);
    result.neighbour_total_row_owner_column_gradient
        .resize(owner_q);
    result.neighbour_total_row_neighbour_column_gradient
        .resize(neighbour_q);

    for (std::size_t component = 0U;
         component < n;
         ++component) {
        result.neighbour_component_face_rate_mol_per_s[
            component] =
            -source.component_molar_flux_mol_per_s[
                component];

        if (result.owner_component_face_rate_mol_per_s[
                component] +
                result.neighbour_component_face_rate_mol_per_s[
                    component] !=
            0.0) {
            throw std::runtime_error(
                "mpmc::flow_discretization: component face-rate scatter is not exactly conservative");
        }

        for (std::size_t column = 0U;
             column < owner_q;
             ++column) {
            const std::size_t index =
                component * owner_q +
                column;
            result.neighbour_row_owner_column_jacobian[
                index] =
                -source.owner_component_flux_jacobian[
                    index];

            if (result.owner_row_owner_column_jacobian[
                    index] +
                    result.neighbour_row_owner_column_jacobian[
                        index] !=
                0.0) {
                throw std::runtime_error(
                    "mpmc::flow_discretization: owner-column component face-rate Jacobian scatter is not exactly antisymmetric");
            }
        }

        for (std::size_t column = 0U;
             column < neighbour_q;
             ++column) {
            const std::size_t index =
                component * neighbour_q +
                column;
            result.neighbour_row_neighbour_column_jacobian[
                index] =
                -source.neighbour_component_flux_jacobian[
                    index];

            if (result.owner_row_neighbour_column_jacobian[
                    index] +
                    result.neighbour_row_neighbour_column_jacobian[
                        index] !=
                0.0) {
                throw std::runtime_error(
                    "mpmc::flow_discretization: neighbour-column component face-rate Jacobian scatter is not exactly antisymmetric");
            }
        }
    }

    for (std::size_t column = 0U;
         column < owner_q;
         ++column) {
        result.neighbour_total_row_owner_column_gradient[
            column] =
            -source.owner_total_molar_flux_gradient[
                column];
        if (result
                    .owner_total_row_owner_column_gradient[
                        column] +
                result
                    .neighbour_total_row_owner_column_gradient[
                        column] !=
            0.0) {
            throw std::runtime_error(
                "mpmc::flow_discretization: owner-column total face-rate gradient scatter is not exactly antisymmetric");
        }
    }

    for (std::size_t column = 0U;
         column < neighbour_q;
         ++column) {
        result.neighbour_total_row_neighbour_column_gradient[
            column] =
            -source.neighbour_total_molar_flux_gradient[
                column];
        if (result
                    .owner_total_row_neighbour_column_gradient[
                        column] +
                result
                    .neighbour_total_row_neighbour_column_gradient[
                        column] !=
            0.0) {
            throw std::runtime_error(
                "mpmc::flow_discretization: neighbour-column total face-rate gradient scatter is not exactly antisymmetric");
        }
    }

    if (result.owner_total_face_rate_mol_per_s +
            result.neighbour_total_face_rate_mol_per_s !=
        0.0) {
        throw std::runtime_error(
            "mpmc::flow_discretization: total face-rate scatter is not exactly conservative");
    }

    return result;
}

} // namespace mpmc::flow_discretization

#endif // MPMC_FLOW_DISCRETIZATION_CONSERVATIVE_COMPONENT_FACE_RATE_SCATTER_HPP
