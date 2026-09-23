#ifndef MPMC_FLOW_DISCRETIZATION_NORMALIZED_COMPONENT_FACE_CONTRIBUTION_HPP
#define MPMC_FLOW_DISCRETIZATION_NORMALIZED_COMPONENT_FACE_CONTRIBUTION_HPP

#include <mpmc/flow_discretization/conservative_component_face_rate_scatter.hpp>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <stdexcept>
#include <string>
#include <string_view>
#include <vector>

namespace mpmc::flow_discretization {

inline constexpr std::string_view
    normalized_component_face_contribution_convention =
        "flow_discretization/rigid-cell-volume-normalized-component-face-contribution/v1";

/// Explicit geometric bulk volumes for the two cells adjacent to one internal
/// face. Units are [m^3].
///
/// v1 is a rigid-grid contract:
///   d V_b,owner = d V_b,neighbour = 0
///
/// The caller must obtain these values from authoritative mesh geometry. This
/// type does not compute or infer cell volume from porosity, face geometry or
/// transmissibility.
struct TwoCellBulkVolume3D {
    double owner_bulk_volume_m3{};
    double neighbour_bulk_volume_m3{};
};

/// Volume-normalized spatial contribution from one internal component face
/// rate.
///
/// For each canonical component i:
///
///   R_i,o^(face,V) = R_i,o^face / V_b,o
///   R_i,n^(face,V) = R_i,n^face / V_b,n
///
/// Units are [mol / (bulk-m^3 s)], matching the existing backward-Euler
/// accumulation residual units. This type still represents only one face's
/// spatial contribution; it is not a complete cell conservation residual.
///
/// Because owner and neighbour volumes may differ, the normalized values are
/// generally not direct negatives. Conservation is volume weighted:
///
///   V_b,o R_i,o^(face,V) + V_b,n R_i,n^(face,V) = 0.
struct NormalizedComponentFaceContributionLinearization3D {
    static constexpr std::string_view convention =
        normalized_component_face_contribution_convention;
    static constexpr bool bulk_volume_derivative_is_zero = true;

    mpmc::mesh::LocalIndex face{
        mpmc::mesh::LocalIndex::value_type{0}};

    TwoCellBulkVolume3D bulk_volume;

    std::vector<std::string> component_ids;

    mpmc::flow::NaturalVariableStateIdentity3P
        owner_state_identity;
    mpmc::flow::NaturalVariableStateIdentity3P
        neighbour_state_identity;

    std::vector<double>
        owner_component_contribution_mol_per_bulk_m3_s;
    std::vector<double>
        neighbour_component_contribution_mol_per_bulk_m3_s;

    /// d R_owner^(face,V) / d q_owner
    std::vector<double>
        owner_row_owner_column_jacobian;
    /// d R_owner^(face,V) / d q_neighbour
    std::vector<double>
        owner_row_neighbour_column_jacobian;
    /// d R_neighbour^(face,V) / d q_owner
    std::vector<double>
        neighbour_row_owner_column_jacobian;
    /// d R_neighbour^(face,V) / d q_neighbour
    std::vector<double>
        neighbour_row_neighbour_column_jacobian;

    double owner_total_contribution_mol_per_bulk_m3_s{};
    double neighbour_total_contribution_mol_per_bulk_m3_s{};

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

    [[nodiscard]] double owner_contribution(
        std::size_t component) const {
        return owner_component_contribution_mol_per_bulk_m3_s.at(
            component);
    }

    [[nodiscard]] double neighbour_contribution(
        std::size_t component) const {
        return neighbour_component_contribution_mol_per_bulk_m3_s.at(
            component);
    }

    [[nodiscard]] double
    d_owner_contribution_wrt_owner(
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
    d_owner_contribution_wrt_neighbour(
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
    d_neighbour_contribution_wrt_owner(
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
    d_neighbour_contribution_wrt_neighbour(
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
                " normalized face-contribution Jacobian index out of range");
        }
        return values.at(
            component * input_count +
            column);
    }
};

namespace normalized_face_contribution_detail {

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

inline void validate_volume(
    double volume_m3,
    const char* side) {
    if (!std::isfinite(volume_m3) ||
        !(volume_m3 > 0.0)) {
        throw std::invalid_argument(
            std::string{"mpmc::flow_discretization: "} +
            side +
            " cell bulk volume must be finite and strictly positive [m^3]");
    }
}

inline void validate_scatter(
    const ConservativeComponentFaceRateScatterLinearization3D&
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
            "mpmc::flow_discretization: conservative face-rate scatter identity/chart is invalid");
    }

    for (std::size_t component = 0U;
         component < n;
         ++component) {
        if (source.component_ids[component].empty()) {
            throw std::invalid_argument(
                "mpmc::flow_discretization: normalized spatial component ID must not be empty");
        }
        for (std::size_t previous = 0U;
             previous < component;
             ++previous) {
            if (source.component_ids[previous] ==
                source.component_ids[component]) {
                throw std::invalid_argument(
                    "mpmc::flow_discretization: normalized spatial component IDs must be unique and ordered");
            }
        }
    }

    if (source.owner_component_face_rate_mol_per_s.size() !=
            n ||
        source.neighbour_component_face_rate_mol_per_s.size() !=
            n ||
        n >
            std::numeric_limits<std::size_t>::max() /
                owner_q ||
        n >
            std::numeric_limits<std::size_t>::max() /
                neighbour_q ||
        source.owner_row_owner_column_jacobian.size() !=
            n * owner_q ||
        source.owner_row_neighbour_column_jacobian.size() !=
            n * neighbour_q ||
        source.neighbour_row_owner_column_jacobian.size() !=
            n * owner_q ||
        source.neighbour_row_neighbour_column_jacobian.size() !=
            n * neighbour_q ||
        source.owner_total_row_owner_column_gradient.size() !=
            owner_q ||
        source.owner_total_row_neighbour_column_gradient.size() !=
            neighbour_q ||
        source.neighbour_total_row_owner_column_gradient.size() !=
            owner_q ||
        source.neighbour_total_row_neighbour_column_gradient.size() !=
            neighbour_q ||
        !std::isfinite(
            source.owner_total_face_rate_mol_per_s) ||
        !std::isfinite(
            source.neighbour_total_face_rate_mol_per_s)) {
        throw std::invalid_argument(
            "mpmc::flow_discretization: conservative face-rate scatter value/Jacobian shape is invalid");
    }

    if (source.owner_total_face_rate_mol_per_s +
            source.neighbour_total_face_rate_mol_per_s !=
        0.0) {
        throw std::invalid_argument(
            "mpmc::flow_discretization: total face-rate scatter is not exactly conservative");
    }

    double owner_component_sum = 0.0;
    double neighbour_component_sum = 0.0;
    double primal_scale =
        std::abs(source.owner_total_face_rate_mol_per_s) +
        std::abs(source.neighbour_total_face_rate_mol_per_s);
    for (std::size_t component = 0U;
         component < n;
         ++component) {
        owner_component_sum +=
            source.owner_component_face_rate_mol_per_s[
                component];
        neighbour_component_sum +=
            source.neighbour_component_face_rate_mol_per_s[
                component];
        primal_scale +=
            std::abs(
                source.owner_component_face_rate_mol_per_s[
                    component]) +
            std::abs(
                source.neighbour_component_face_rate_mol_per_s[
                    component]);
    }
    if (!near_roundoff(
            owner_component_sum,
            source.owner_total_face_rate_mol_per_s,
            primal_scale) ||
        !near_roundoff(
            neighbour_component_sum,
            source.neighbour_total_face_rate_mol_per_s,
            primal_scale)) {
        throw std::invalid_argument(
            "mpmc::flow_discretization: component face-rate scatter does not close to total face rate");
    }

    for (std::size_t component = 0U;
         component < n;
         ++component) {
        const double owner_rate =
            source.owner_component_face_rate_mol_per_s[
                component];
        const double neighbour_rate =
            source.neighbour_component_face_rate_mol_per_s[
                component];
        if (!std::isfinite(owner_rate) ||
            !std::isfinite(neighbour_rate) ||
            owner_rate + neighbour_rate != 0.0) {
            throw std::invalid_argument(
                "mpmc::flow_discretization: component face-rate scatter is not exactly conservative");
        }

        for (std::size_t column = 0U;
             column < owner_q;
             ++column) {
            const std::size_t index =
                component * owner_q +
                column;
            const double owner_value =
                source.owner_row_owner_column_jacobian[
                    index];
            const double neighbour_value =
                source.neighbour_row_owner_column_jacobian[
                    index];
            if (!std::isfinite(owner_value) ||
                !std::isfinite(neighbour_value) ||
                owner_value + neighbour_value != 0.0) {
                throw std::invalid_argument(
                    "mpmc::flow_discretization: owner-column face-rate Jacobian scatter is not exactly antisymmetric");
            }
        }

        for (std::size_t column = 0U;
             column < neighbour_q;
             ++column) {
            const std::size_t index =
                component * neighbour_q +
                column;
            const double owner_value =
                source.owner_row_neighbour_column_jacobian[
                    index];
            const double neighbour_value =
                source.neighbour_row_neighbour_column_jacobian[
                    index];
            if (!std::isfinite(owner_value) ||
                !std::isfinite(neighbour_value) ||
                owner_value + neighbour_value != 0.0) {
                throw std::invalid_argument(
                    "mpmc::flow_discretization: neighbour-column face-rate Jacobian scatter is not exactly antisymmetric");
            }
        }
    }

    for (std::size_t column = 0U;
         column < owner_q;
         ++column) {
        const double owner_value =
            source.owner_total_row_owner_column_gradient[
                column];
        const double neighbour_value =
            source.neighbour_total_row_owner_column_gradient[
                column];
        if (!std::isfinite(owner_value) ||
            !std::isfinite(neighbour_value) ||
            owner_value + neighbour_value != 0.0) {
            throw std::invalid_argument(
                "mpmc::flow_discretization: owner-column total face-rate gradient scatter is not exactly antisymmetric");
        }

        double owner_component_derivative_sum = 0.0;
        double neighbour_component_derivative_sum = 0.0;
        double derivative_scale =
            std::abs(owner_value) +
            std::abs(neighbour_value);
        for (std::size_t component = 0U;
             component < n;
             ++component) {
            const std::size_t index =
                component * owner_q +
                column;
            const double owner_component_derivative =
                source.owner_row_owner_column_jacobian[
                    index];
            const double neighbour_component_derivative =
                source.neighbour_row_owner_column_jacobian[
                    index];
            owner_component_derivative_sum +=
                owner_component_derivative;
            neighbour_component_derivative_sum +=
                neighbour_component_derivative;
            derivative_scale +=
                std::abs(owner_component_derivative) +
                std::abs(neighbour_component_derivative);
        }
        if (!near_roundoff(
                owner_component_derivative_sum,
                owner_value,
                derivative_scale) ||
            !near_roundoff(
                neighbour_component_derivative_sum,
                neighbour_value,
                derivative_scale)) {
            throw std::invalid_argument(
                "mpmc::flow_discretization: owner-column component face-rate Jacobian does not close to total derivative");
        }
    }

    for (std::size_t column = 0U;
         column < neighbour_q;
         ++column) {
        const double owner_value =
            source.owner_total_row_neighbour_column_gradient[
                column];
        const double neighbour_value =
            source.neighbour_total_row_neighbour_column_gradient[
                column];
        if (!std::isfinite(owner_value) ||
            !std::isfinite(neighbour_value) ||
            owner_value + neighbour_value != 0.0) {
            throw std::invalid_argument(
                "mpmc::flow_discretization: neighbour-column total face-rate gradient scatter is not exactly antisymmetric");
        }

        double owner_component_derivative_sum = 0.0;
        double neighbour_component_derivative_sum = 0.0;
        double derivative_scale =
            std::abs(owner_value) +
            std::abs(neighbour_value);
        for (std::size_t component = 0U;
             component < n;
             ++component) {
            const std::size_t index =
                component * neighbour_q +
                column;
            const double owner_component_derivative =
                source.owner_row_neighbour_column_jacobian[
                    index];
            const double neighbour_component_derivative =
                source.neighbour_row_neighbour_column_jacobian[
                    index];
            owner_component_derivative_sum +=
                owner_component_derivative;
            neighbour_component_derivative_sum +=
                neighbour_component_derivative;
            derivative_scale +=
                std::abs(owner_component_derivative) +
                std::abs(neighbour_component_derivative);
        }
        if (!near_roundoff(
                owner_component_derivative_sum,
                owner_value,
                derivative_scale) ||
            !near_roundoff(
                neighbour_component_derivative_sum,
                neighbour_value,
                derivative_scale)) {
            throw std::invalid_argument(
                "mpmc::flow_discretization: neighbour-column component face-rate Jacobian does not close to total derivative");
        }
    }
}

} // namespace normalized_face_contribution_detail

[[nodiscard]] inline
NormalizedComponentFaceContributionLinearization3D
normalize_component_face_rate_by_bulk_volume(
    const ConservativeComponentFaceRateScatterLinearization3D&
        source,
    TwoCellBulkVolume3D bulk_volume) {
    using namespace normalized_face_contribution_detail;

    validate_scatter(source);
    validate_volume(
        bulk_volume.owner_bulk_volume_m3,
        "owner");
    validate_volume(
        bulk_volume.neighbour_bulk_volume_m3,
        "neighbour");

    const std::size_t n =
        source.component_ids.size();
    const std::size_t owner_q =
        source.owner_state_identity.layout
            .unknown_count();
    const std::size_t neighbour_q =
        source.neighbour_state_identity.layout
            .unknown_count();

    const double inverse_owner_volume =
        1.0 /
        bulk_volume.owner_bulk_volume_m3;
    const double inverse_neighbour_volume =
        1.0 /
        bulk_volume.neighbour_bulk_volume_m3;
    if (!std::isfinite(inverse_owner_volume) ||
        !std::isfinite(inverse_neighbour_volume)) {
        throw std::range_error(
            "mpmc::flow_discretization: inverse cell bulk volume is non-finite");
    }

    NormalizedComponentFaceContributionLinearization3D
        result{
            source.face,
            bulk_volume,
            source.component_ids,
            source.owner_state_identity,
            source.neighbour_state_identity,
            std::vector<double>(n, 0.0),
            std::vector<double>(n, 0.0),
            std::vector<double>(
                n * owner_q,
                0.0),
            std::vector<double>(
                n * neighbour_q,
                0.0),
            std::vector<double>(
                n * owner_q,
                0.0),
            std::vector<double>(
                n * neighbour_q,
                0.0),
            source.owner_total_face_rate_mol_per_s *
                inverse_owner_volume,
            source.neighbour_total_face_rate_mol_per_s *
                inverse_neighbour_volume,
            std::vector<double>(owner_q, 0.0),
            std::vector<double>(neighbour_q, 0.0),
            std::vector<double>(owner_q, 0.0),
            std::vector<double>(neighbour_q, 0.0)};

    for (std::size_t component = 0U;
         component < n;
         ++component) {
        result.owner_component_contribution_mol_per_bulk_m3_s[
            component] =
            source.owner_component_face_rate_mol_per_s[
                component] *
            inverse_owner_volume;
        result.neighbour_component_contribution_mol_per_bulk_m3_s[
            component] =
            source.neighbour_component_face_rate_mol_per_s[
                component] *
            inverse_neighbour_volume;

        const double weighted_sum =
            bulk_volume.owner_bulk_volume_m3 *
                result.owner_component_contribution_mol_per_bulk_m3_s[
                    component] +
            bulk_volume.neighbour_bulk_volume_m3 *
                result.neighbour_component_contribution_mol_per_bulk_m3_s[
                    component];
        const double weighted_scale =
            std::abs(
                source.owner_component_face_rate_mol_per_s[
                    component]) +
            std::abs(
                source.neighbour_component_face_rate_mol_per_s[
                    component]);
        if (!near_roundoff(
                weighted_sum,
                0.0,
                weighted_scale)) {
            throw std::runtime_error(
                "mpmc::flow_discretization: volume-normalized component face contributions do not preserve weighted conservation");
        }

        for (std::size_t column = 0U;
             column < owner_q;
             ++column) {
            const std::size_t index =
                component * owner_q +
                column;
            result.owner_row_owner_column_jacobian[
                index] =
                source.owner_row_owner_column_jacobian[
                    index] *
                inverse_owner_volume;
            result.neighbour_row_owner_column_jacobian[
                index] =
                source.neighbour_row_owner_column_jacobian[
                    index] *
                inverse_neighbour_volume;

            const double weighted_derivative =
                bulk_volume.owner_bulk_volume_m3 *
                    result.owner_row_owner_column_jacobian[
                        index] +
                bulk_volume.neighbour_bulk_volume_m3 *
                    result.neighbour_row_owner_column_jacobian[
                        index];
            const double derivative_scale =
                std::abs(
                    source.owner_row_owner_column_jacobian[
                        index]) +
                std::abs(
                    source.neighbour_row_owner_column_jacobian[
                        index]);
            if (!near_roundoff(
                    weighted_derivative,
                    0.0,
                    derivative_scale)) {
                throw std::runtime_error(
                    "mpmc::flow_discretization: owner-column volume-normalized component Jacobian violates weighted conservation");
            }
        }

        for (std::size_t column = 0U;
             column < neighbour_q;
             ++column) {
            const std::size_t index =
                component * neighbour_q +
                column;
            result.owner_row_neighbour_column_jacobian[
                index] =
                source.owner_row_neighbour_column_jacobian[
                    index] *
                inverse_owner_volume;
            result.neighbour_row_neighbour_column_jacobian[
                index] =
                source.neighbour_row_neighbour_column_jacobian[
                    index] *
                inverse_neighbour_volume;

            const double weighted_derivative =
                bulk_volume.owner_bulk_volume_m3 *
                    result.owner_row_neighbour_column_jacobian[
                        index] +
                bulk_volume.neighbour_bulk_volume_m3 *
                    result.neighbour_row_neighbour_column_jacobian[
                        index];
            const double derivative_scale =
                std::abs(
                    source.owner_row_neighbour_column_jacobian[
                        index]) +
                std::abs(
                    source.neighbour_row_neighbour_column_jacobian[
                        index]);
            if (!near_roundoff(
                    weighted_derivative,
                    0.0,
                    derivative_scale)) {
                throw std::runtime_error(
                    "mpmc::flow_discretization: neighbour-column volume-normalized component Jacobian violates weighted conservation");
            }
        }
    }

    for (std::size_t column = 0U;
         column < owner_q;
         ++column) {
        result.owner_total_row_owner_column_gradient[
            column] =
            source.owner_total_row_owner_column_gradient[
                column] *
            inverse_owner_volume;
        result.neighbour_total_row_owner_column_gradient[
            column] =
            source.neighbour_total_row_owner_column_gradient[
                column] *
            inverse_neighbour_volume;
    }

    for (std::size_t column = 0U;
         column < neighbour_q;
         ++column) {
        result.owner_total_row_neighbour_column_gradient[
            column] =
            source.owner_total_row_neighbour_column_gradient[
                column] *
            inverse_owner_volume;
        result.neighbour_total_row_neighbour_column_gradient[
            column] =
            source.neighbour_total_row_neighbour_column_gradient[
                column] *
            inverse_neighbour_volume;
    }

    const double total_weighted_sum =
        bulk_volume.owner_bulk_volume_m3 *
            result.owner_total_contribution_mol_per_bulk_m3_s +
        bulk_volume.neighbour_bulk_volume_m3 *
            result.neighbour_total_contribution_mol_per_bulk_m3_s;
    if (!near_roundoff(
            total_weighted_sum,
            0.0,
            std::abs(
                source.owner_total_face_rate_mol_per_s) +
                std::abs(
                    source.neighbour_total_face_rate_mol_per_s))) {
        throw std::runtime_error(
            "mpmc::flow_discretization: volume-normalized total face contribution violates weighted conservation");
    }

    for (std::size_t column = 0U;
         column < owner_q;
         ++column) {
        const double weighted_derivative =
            bulk_volume.owner_bulk_volume_m3 *
                result.owner_total_row_owner_column_gradient[
                    column] +
            bulk_volume.neighbour_bulk_volume_m3 *
                result.neighbour_total_row_owner_column_gradient[
                    column];
        if (!near_roundoff(
                weighted_derivative,
                0.0,
                std::abs(
                    source.owner_total_row_owner_column_gradient[
                        column]) +
                    std::abs(
                        source.neighbour_total_row_owner_column_gradient[
                            column]))) {
            throw std::runtime_error(
                "mpmc::flow_discretization: owner-column normalized total gradient violates weighted conservation");
        }
    }

    for (std::size_t column = 0U;
         column < neighbour_q;
         ++column) {
        const double weighted_derivative =
            bulk_volume.owner_bulk_volume_m3 *
                result.owner_total_row_neighbour_column_gradient[
                    column] +
            bulk_volume.neighbour_bulk_volume_m3 *
                result.neighbour_total_row_neighbour_column_gradient[
                    column];
        if (!near_roundoff(
                weighted_derivative,
                0.0,
                std::abs(
                    source.owner_total_row_neighbour_column_gradient[
                        column]) +
                    std::abs(
                        source.neighbour_total_row_neighbour_column_gradient[
                            column]))) {
            throw std::runtime_error(
                "mpmc::flow_discretization: neighbour-column normalized total gradient violates weighted conservation");
        }
    }

    return result;
}

} // namespace mpmc::flow_discretization

#endif // MPMC_FLOW_DISCRETIZATION_NORMALIZED_COMPONENT_FACE_CONTRIBUTION_HPP
