#ifndef MPMC_FLOW_DISCRETIZATION_LOCAL_COMPONENT_CONSERVATION_RESIDUAL_HPP
#define MPMC_FLOW_DISCRETIZATION_LOCAL_COMPONENT_CONSERVATION_RESIDUAL_HPP

#include <mpmc/flow/component_accumulation_time.hpp>
#include <mpmc/flow_discretization/normalized_component_face_contribution.hpp>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <span>
#include <stdexcept>
#include <string>
#include <string_view>
#include <vector>

namespace mpmc::flow_discretization {

inline constexpr std::string_view
    local_component_conservation_residual_convention =
        "flow_discretization/local-component-conservation/multi-internal-face/v1";

/// Which side of one already-normalized internal-face contribution belongs to
/// the local cell whose conservation row is being assembled.
enum class IncidentFaceLocalSide3D : std::uint8_t {
    owner,
    neighbour
};

/// Borrowed binding used only for the duration of
/// build_local_component_conservation_residual().
///
/// The side is explicit. The builder never infers cell identity from equal
/// pressure/temperature/composition values because two different cells may
/// legitimately have numerically identical states.
struct LocalCellNormalizedFaceContributionBinding3D {
    IncidentFaceLocalSide3D local_side{
        IncidentFaceLocalSide3D::owner};
    const NormalizedComponentFaceContributionLinearization3D*
        contribution{};
};

/// One off-diagonal natural-variable Jacobian block contributed by one
/// internal face.
///
/// The block is intentionally keyed by the stable local face identity, not by
/// a guessed neighbour-cell identity. Global topology/numbering may later map
/// this face to the authoritative neighbour cell and coalesce duplicate
/// cell-pair blocks if a topology permits more than one connecting face.
struct LocalComponentConservationNeighbourBlock3D {
    mpmc::mesh::LocalIndex face{
        mpmc::mesh::LocalIndex::value_type{0}};
    IncidentFaceLocalSide3D local_side{
        IncidentFaceLocalSide3D::owner};
    mpmc::flow::NaturalVariableStateIdentity3P
        neighbour_state_identity;
    std::vector<double> component_jacobian;
    std::vector<double> total_gradient;

    [[nodiscard]] double d_component(
        std::size_t component,
        std::size_t column) const {
        const std::size_t q =
            neighbour_state_identity.layout.unknown_count();
        if (component >=
                neighbour_state_identity.component_ids.size() ||
            column >= q ||
            q == 0U ||
            component >
                (std::numeric_limits<std::size_t>::max() -
                 column) /
                    q) {
            throw std::out_of_range(
                "mpmc::flow_discretization: local conservation neighbour-block Jacobian index out of range");
        }
        return component_jacobian.at(
            component * q + column);
    }

    [[nodiscard]] double d_total(
        std::size_t column) const {
        return total_gradient.at(column);
    }
};

/// Complete component-conservation residual for one fixed-three-phase cell,
/// restricted to backward-Euler accumulation plus already-normalized internal
/// faces:
///
///   R_i,c = R_i,c^acc + sum_{f in internal(c)} R_i,c^(face,V)
///
/// Every term is [mol / (bulk-m^3 s)].
///
/// This is still a local algebra contract. It does not add boundary/source/well
/// terms, energy/fugacity rows, global row/column numbering, PETSc values or a
/// nonlinear solve.
struct LocalComponentConservationResidualLinearization3D {
    static constexpr std::string_view convention =
        local_component_conservation_residual_convention;

    mpmc::flow::NaturalVariableStateIdentity3P
        cell_state_identity;
    double cell_bulk_volume_m3{};
    double porosity{};
    double time_step_seconds{};

    std::vector<std::string> component_ids;
    std::vector<double>
        component_residual_mol_per_bulk_m3_s;
    double total_residual_mol_per_bulk_m3_s{};

    /// d R_local / d q_local
    std::vector<double> local_component_jacobian;
    std::vector<double> local_total_gradient;

    /// One off-diagonal local-row block per incident internal face.
    std::vector<LocalComponentConservationNeighbourBlock3D>
        neighbour_blocks;

    [[nodiscard]] std::size_t component_count() const noexcept {
        return component_ids.size();
    }

    [[nodiscard]] double residual(
        std::size_t component) const {
        return component_residual_mol_per_bulk_m3_s.at(
            component);
    }

    [[nodiscard]] double d_local(
        std::size_t component,
        std::size_t column) const {
        const std::size_t q =
            cell_state_identity.layout.unknown_count();
        if (component >= component_ids.size() ||
            column >= q ||
            q == 0U ||
            component >
                (std::numeric_limits<std::size_t>::max() -
                 column) /
                    q) {
            throw std::out_of_range(
                "mpmc::flow_discretization: local conservation diagonal Jacobian index out of range");
        }
        return local_component_jacobian.at(
            component * q + column);
    }

    [[nodiscard]] double d_total_local(
        std::size_t column) const {
        return local_total_gradient.at(column);
    }
};

namespace local_component_conservation_detail {

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

[[nodiscard]] inline bool same_layout(
    const mpmc::flow::NaturalVariableLayoutDescriptor& first,
    const mpmc::flow::NaturalVariableLayoutDescriptor& second) {
    return first.component_count() ==
               second.component_count() &&
        first.phase_count() ==
            second.phase_count() &&
        first.unknown_count() ==
            second.unknown_count() &&
        first.composition_pivot()
                .dependent_components() ==
            second.composition_pivot()
                .dependent_components();
}

inline void validate_state_identity(
    const mpmc::flow::NaturalVariableStateIdentity3P&
        identity,
    const char* name) {
    const std::size_t n =
        identity.layout.component_count();
    if (n < 2U ||
        identity.component_ids.size() != n ||
        !std::isfinite(identity.reference_pressure_pa) ||
        !(identity.reference_pressure_pa > 0.0) ||
        !std::isfinite(identity.temperature_k) ||
        !(identity.temperature_k > 0.0)) {
        throw std::invalid_argument(
            std::string{"mpmc::flow_discretization: invalid "} +
            name +
            " natural-variable state identity");
    }

    for (std::size_t component = 0U;
         component < n;
         ++component) {
        if (identity.component_ids[component].empty()) {
            throw std::invalid_argument(
                std::string{"mpmc::flow_discretization: "} +
                name +
                " component identity must be nonempty");
        }
        for (std::size_t previous = 0U;
             previous < component;
             ++previous) {
            if (identity.component_ids[previous] ==
                identity.component_ids[component]) {
                throw std::invalid_argument(
                    std::string{"mpmc::flow_discretization: "} +
                    name +
                    " component identity must be unique and ordered");
            }
        }
    }

    double saturation_sum = 0.0;
    for (std::size_t phase = 0U;
         phase < identity.layout.phase_count();
         ++phase) {
        const double value =
            identity.saturation[phase];
        if (!std::isfinite(value) ||
            !(value > 0.0)) {
            throw std::invalid_argument(
                std::string{"mpmc::flow_discretization: "} +
                name +
                " active saturation identity must have positive finite support");
        }
        saturation_sum += value;
    }
    if (!near_roundoff(
            saturation_sum,
            1.0,
            saturation_sum)) {
        throw std::invalid_argument(
            std::string{"mpmc::flow_discretization: "} +
            name +
            " saturation identity must close to one");
    }

    for (std::size_t phase = 0U;
         phase < identity.layout.phase_count();
         ++phase) {
        const auto& composition =
            identity.phase_composition[phase];
        if (composition.size() != n) {
            throw std::invalid_argument(
                std::string{"mpmc::flow_discretization: "} +
                name +
                " phase-composition identity shape mismatch");
        }

        double sum = 0.0;
        for (double value : composition) {
            if (!std::isfinite(value) ||
                !(value > 0.0)) {
                throw std::invalid_argument(
                    std::string{"mpmc::flow_discretization: "} +
                    name +
                    " phase-composition identity must have positive finite support");
            }
            sum += value;
        }
        if (!near_roundoff(sum, 1.0, sum)) {
            throw std::invalid_argument(
                std::string{"mpmc::flow_discretization: "} +
                name +
                " phase-composition identity must close to one");
        }
    }
    for (std::size_t phase =
             identity.layout.phase_count();
         phase < mpmc::flow::fixed_three_phase_count;
         ++phase) {
        if (identity.saturation[phase] != 0.0 ||
            !identity.phase_composition[phase].empty()) {
            throw std::invalid_argument(
                std::string{"mpmc::flow_discretization: "} +
                name +
                " inactive phase sidecar must be zero/empty");
        }
    }
}

[[nodiscard]] inline bool same_state_identity(
    const mpmc::flow::NaturalVariableStateIdentity3P& first,
    const mpmc::flow::NaturalVariableStateIdentity3P& second) {
    if (!same_layout(first.layout, second.layout) ||
        first.component_ids != second.component_ids ||
        !near_roundoff(
            first.reference_pressure_pa,
            second.reference_pressure_pa) ||
        !near_roundoff(
            first.temperature_k,
            second.temperature_k)) {
        return false;
    }

    for (std::size_t phase = 0U;
         phase < first.layout.phase_count();
         ++phase) {
        if (!near_roundoff(
                first.saturation[phase],
                second.saturation[phase]) ||
            first.phase_composition[phase].size() !=
                second.phase_composition[phase].size()) {
            return false;
        }
        for (std::size_t component = 0U;
             component <
             first.phase_composition[phase].size();
             ++component) {
            if (!near_roundoff(
                    first.phase_composition[phase][component],
                    second.phase_composition[phase][component])) {
                return false;
            }
        }
    }
    for (std::size_t phase =
             first.layout.phase_count();
         phase < mpmc::flow::fixed_three_phase_count;
         ++phase) {
        if (first.saturation[phase] != 0.0 ||
            second.saturation[phase] != 0.0 ||
            !first.phase_composition[phase].empty() ||
            !second.phase_composition[phase].empty()) {
            return false;
        }
    }
    return true;
}

inline void validate_accumulation(
    const mpmc::flow::NaturalVariableStateIdentity3P&
        cell_identity,
    const mpmc::flow::
        BackwardEulerComponentAccumulationResidual3P&
            accumulation) {
    const std::size_t n =
        cell_identity.component_ids.size();
    const std::size_t q =
        cell_identity.layout.unknown_count();

    if (q == 0U ||
        n > std::numeric_limits<std::size_t>::max() / q ||
        accumulation.component_ids !=
            cell_identity.component_ids ||
        !same_layout(
            accumulation.current_layout,
            cell_identity.layout) ||
        accumulation.input_count != q ||
        accumulation.component_residual_mol_per_bulk_m3_s.size() !=
            n ||
        accumulation.component_jacobian.size() !=
            n * q ||
        accumulation.total_residual_gradient.size() != q ||
        !std::isfinite(accumulation.porosity) ||
        !(accumulation.porosity > 0.0) ||
        accumulation.porosity > 1.0 ||
        !std::isfinite(accumulation.time_step_seconds) ||
        !(accumulation.time_step_seconds > 0.0) ||
        !std::isfinite(
            accumulation.total_residual_mol_per_bulk_m3_s)) {
        throw std::invalid_argument(
            "mpmc::flow_discretization: backward-Euler accumulation payload does not match local conservation state/chart");
    }

    double component_sum = 0.0;
    double component_scale =
        std::abs(
            accumulation.total_residual_mol_per_bulk_m3_s);
    for (double value :
         accumulation.component_residual_mol_per_bulk_m3_s) {
        if (!std::isfinite(value)) {
            throw std::invalid_argument(
                "mpmc::flow_discretization: backward-Euler accumulation contains non-finite residual");
        }
        component_sum += value;
        component_scale += std::abs(value);
    }
    if (!near_roundoff(
            component_sum,
            accumulation.total_residual_mol_per_bulk_m3_s,
            component_scale)) {
        throw std::invalid_argument(
            "mpmc::flow_discretization: backward-Euler component residuals do not close to total");
    }

    for (std::size_t column = 0U;
         column < q;
         ++column) {
        const double total =
            accumulation.total_residual_gradient[column];
        if (!std::isfinite(total)) {
            throw std::invalid_argument(
                "mpmc::flow_discretization: backward-Euler accumulation contains non-finite total derivative");
        }

        double sum = 0.0;
        double scale = std::abs(total);
        for (std::size_t component = 0U;
             component < n;
             ++component) {
            const double derivative =
                accumulation.component_jacobian[
                    component * q + column];
            if (!std::isfinite(derivative)) {
                throw std::invalid_argument(
                    "mpmc::flow_discretization: backward-Euler accumulation contains non-finite component derivative");
            }
            sum += derivative;
            scale += std::abs(derivative);
        }
        if (!near_roundoff(sum, total, scale)) {
            throw std::invalid_argument(
                "mpmc::flow_discretization: backward-Euler accumulation Jacobian does not close to total derivative");
        }
    }
}

inline void validate_normalized_face(
    const NormalizedComponentFaceContributionLinearization3D&
        face) {
    const std::size_t n = face.component_ids.size();
    const std::size_t owner_q =
        face.owner_state_identity.layout.unknown_count();
    const std::size_t neighbour_q =
        face.neighbour_state_identity.layout.unknown_count();

    validate_state_identity(
        face.owner_state_identity,
        "owner-face");
    validate_state_identity(
        face.neighbour_state_identity,
        "neighbour-face");

    if (n < 2U ||
        owner_q == 0U ||
        neighbour_q == 0U ||
        n > std::numeric_limits<std::size_t>::max() /
                owner_q ||
        n > std::numeric_limits<std::size_t>::max() /
                neighbour_q ||
        face.owner_state_identity.component_ids !=
            face.component_ids ||
        face.neighbour_state_identity.component_ids !=
            face.component_ids ||
        face.owner_state_identity.layout.component_count() != n ||
        face.neighbour_state_identity.layout.component_count() != n ||
        !std::isfinite(
            face.bulk_volume.owner_bulk_volume_m3) ||
        !(face.bulk_volume.owner_bulk_volume_m3 > 0.0) ||
        !std::isfinite(
            face.bulk_volume.neighbour_bulk_volume_m3) ||
        !(face.bulk_volume.neighbour_bulk_volume_m3 > 0.0) ||
        face.owner_component_contribution_mol_per_bulk_m3_s.size() != n ||
        face.neighbour_component_contribution_mol_per_bulk_m3_s.size() != n ||
        face.owner_row_owner_column_jacobian.size() !=
            n * owner_q ||
        face.neighbour_row_owner_column_jacobian.size() !=
            n * owner_q ||
        face.owner_row_neighbour_column_jacobian.size() !=
            n * neighbour_q ||
        face.neighbour_row_neighbour_column_jacobian.size() !=
            n * neighbour_q ||
        face.owner_total_row_owner_column_gradient.size() !=
            owner_q ||
        face.neighbour_total_row_owner_column_gradient.size() !=
            owner_q ||
        face.owner_total_row_neighbour_column_gradient.size() !=
            neighbour_q ||
        face.neighbour_total_row_neighbour_column_gradient.size() !=
            neighbour_q ||
        !std::isfinite(
            face.owner_total_contribution_mol_per_bulk_m3_s) ||
        !std::isfinite(
            face.neighbour_total_contribution_mol_per_bulk_m3_s)) {
        throw std::invalid_argument(
            "mpmc::flow_discretization: malformed normalized internal-face contribution payload");
    }

    double owner_component_sum = 0.0;
    double neighbour_component_sum = 0.0;
    double value_scale =
        std::abs(
            face.owner_total_contribution_mol_per_bulk_m3_s) +
        std::abs(
            face.neighbour_total_contribution_mol_per_bulk_m3_s);
    for (std::size_t component = 0U;
         component < n;
         ++component) {
        const double owner =
            face.owner_component_contribution_mol_per_bulk_m3_s[
                component];
        const double neighbour =
            face.neighbour_component_contribution_mol_per_bulk_m3_s[
                component];
        if (!std::isfinite(owner) ||
            !std::isfinite(neighbour)) {
            throw std::invalid_argument(
                "mpmc::flow_discretization: normalized internal-face contribution contains non-finite component value");
        }
        owner_component_sum += owner;
        neighbour_component_sum += neighbour;
        value_scale +=
            std::abs(owner) +
            std::abs(neighbour);

        const double weighted =
            face.bulk_volume.owner_bulk_volume_m3 *
                owner +
            face.bulk_volume.neighbour_bulk_volume_m3 *
                neighbour;
        const double weighted_scale =
            face.bulk_volume.owner_bulk_volume_m3 *
                std::abs(owner) +
            face.bulk_volume.neighbour_bulk_volume_m3 *
                std::abs(neighbour);
        if (!near_roundoff(
                weighted,
                0.0,
                weighted_scale)) {
            throw std::invalid_argument(
                "mpmc::flow_discretization: normalized internal-face component value violates weighted conservation");
        }
    }

    if (!near_roundoff(
            owner_component_sum,
            face.owner_total_contribution_mol_per_bulk_m3_s,
            value_scale) ||
        !near_roundoff(
            neighbour_component_sum,
            face.neighbour_total_contribution_mol_per_bulk_m3_s,
            value_scale)) {
        throw std::invalid_argument(
            "mpmc::flow_discretization: normalized internal-face components do not close to total contribution");
    }

    const double total_weighted =
        face.bulk_volume.owner_bulk_volume_m3 *
            face.owner_total_contribution_mol_per_bulk_m3_s +
        face.bulk_volume.neighbour_bulk_volume_m3 *
            face.neighbour_total_contribution_mol_per_bulk_m3_s;
    if (!near_roundoff(
            total_weighted,
            0.0,
            value_scale *
                std::max(
                    face.bulk_volume.owner_bulk_volume_m3,
                    face.bulk_volume.neighbour_bulk_volume_m3))) {
        throw std::invalid_argument(
            "mpmc::flow_discretization: normalized internal-face total violates weighted conservation");
    }

    const auto validate_column_block =
        [&](std::size_t column_count,
            const std::vector<double>& owner_component_block,
            const std::vector<double>& neighbour_component_block,
            const std::vector<double>& owner_total_block,
            const std::vector<double>& neighbour_total_block,
            const char* column_side) {
            for (std::size_t column = 0U;
                 column < column_count;
                 ++column) {
                const double owner_total =
                    owner_total_block[column];
                const double neighbour_total =
                    neighbour_total_block[column];
                if (!std::isfinite(owner_total) ||
                    !std::isfinite(neighbour_total)) {
                    throw std::invalid_argument(
                        std::string{"mpmc::flow_discretization: normalized internal-face "} +
                        column_side +
                        " total gradient contains non-finite derivative");
                }

                double owner_sum = 0.0;
                double neighbour_sum = 0.0;
                double derivative_scale =
                    std::abs(owner_total) +
                    std::abs(neighbour_total);
                for (std::size_t component = 0U;
                     component < n;
                     ++component) {
                    const std::size_t index =
                        component * column_count +
                        column;
                    const double owner_derivative =
                        owner_component_block[index];
                    const double neighbour_derivative =
                        neighbour_component_block[index];
                    if (!std::isfinite(owner_derivative) ||
                        !std::isfinite(neighbour_derivative)) {
                        throw std::invalid_argument(
                            std::string{"mpmc::flow_discretization: normalized internal-face "} +
                            column_side +
                            " component Jacobian contains non-finite derivative");
                    }
                    owner_sum += owner_derivative;
                    neighbour_sum += neighbour_derivative;
                    derivative_scale +=
                        std::abs(owner_derivative) +
                        std::abs(neighbour_derivative);

                    const double weighted =
                        face.bulk_volume.owner_bulk_volume_m3 *
                            owner_derivative +
                        face.bulk_volume.neighbour_bulk_volume_m3 *
                            neighbour_derivative;
                    const double weighted_scale =
                        face.bulk_volume.owner_bulk_volume_m3 *
                            std::abs(owner_derivative) +
                        face.bulk_volume.neighbour_bulk_volume_m3 *
                            std::abs(neighbour_derivative);
                    if (!near_roundoff(
                            weighted,
                            0.0,
                            weighted_scale)) {
                        throw std::invalid_argument(
                            std::string{"mpmc::flow_discretization: normalized internal-face "} +
                            column_side +
                            " component Jacobian violates weighted conservation");
                    }
                }

                if (!near_roundoff(
                        owner_sum,
                        owner_total,
                        derivative_scale) ||
                    !near_roundoff(
                        neighbour_sum,
                        neighbour_total,
                        derivative_scale)) {
                    throw std::invalid_argument(
                        std::string{"mpmc::flow_discretization: normalized internal-face "} +
                        column_side +
                        " component Jacobian does not close to total gradient");
                }

                const double total_weighted_gradient =
                    face.bulk_volume.owner_bulk_volume_m3 *
                        owner_total +
                    face.bulk_volume.neighbour_bulk_volume_m3 *
                        neighbour_total;
                if (!near_roundoff(
                        total_weighted_gradient,
                        0.0,
                        face.bulk_volume.owner_bulk_volume_m3 *
                                std::abs(owner_total) +
                            face.bulk_volume.neighbour_bulk_volume_m3 *
                                std::abs(neighbour_total))) {
                    throw std::invalid_argument(
                        std::string{"mpmc::flow_discretization: normalized internal-face "} +
                        column_side +
                        " total gradient violates weighted conservation");
                }
            }
        };

    validate_column_block(
        owner_q,
        face.owner_row_owner_column_jacobian,
        face.neighbour_row_owner_column_jacobian,
        face.owner_total_row_owner_column_gradient,
        face.neighbour_total_row_owner_column_gradient,
        "owner-column");

    validate_column_block(
        neighbour_q,
        face.owner_row_neighbour_column_jacobian,
        face.neighbour_row_neighbour_column_jacobian,
        face.owner_total_row_neighbour_column_gradient,
        face.neighbour_total_row_neighbour_column_gradient,
        "neighbour-column");
}

inline void validate_final_closure(
    const LocalComponentConservationResidualLinearization3D&
        result) {
    const std::size_t n = result.component_ids.size();
    const std::size_t q =
        result.cell_state_identity.layout.unknown_count();

    double component_sum = 0.0;
    double value_scale =
        std::abs(result.total_residual_mol_per_bulk_m3_s);
    for (double value :
         result.component_residual_mol_per_bulk_m3_s) {
        if (!std::isfinite(value)) {
            throw std::range_error(
                "mpmc::flow_discretization: local component conservation residual became non-finite");
        }
        component_sum += value;
        value_scale += std::abs(value);
    }
    if (!std::isfinite(
            result.total_residual_mol_per_bulk_m3_s) ||
        !near_roundoff(
            component_sum,
            result.total_residual_mol_per_bulk_m3_s,
            value_scale)) {
        throw std::runtime_error(
            "mpmc::flow_discretization: local component conservation residuals do not close to total residual");
    }

    for (std::size_t column = 0U;
         column < q;
         ++column) {
        const double total =
            result.local_total_gradient[column];
        double sum = 0.0;
        double scale = std::abs(total);
        if (!std::isfinite(total)) {
            throw std::range_error(
                "mpmc::flow_discretization: local conservation total Jacobian became non-finite");
        }
        for (std::size_t component = 0U;
             component < n;
             ++component) {
            const double derivative =
                result.local_component_jacobian[
                    component * q + column];
            if (!std::isfinite(derivative)) {
                throw std::range_error(
                    "mpmc::flow_discretization: local conservation component Jacobian became non-finite");
            }
            sum += derivative;
            scale += std::abs(derivative);
        }
        if (!near_roundoff(sum, total, scale)) {
            throw std::runtime_error(
                "mpmc::flow_discretization: local conservation diagonal component Jacobian does not close to total gradient");
        }
    }

    for (const auto& block : result.neighbour_blocks) {
        const std::size_t neighbour_q =
            block.neighbour_state_identity.layout.unknown_count();
        if (neighbour_q == 0U ||
            n > std::numeric_limits<std::size_t>::max() /
                    neighbour_q ||
            block.component_jacobian.size() !=
                n * neighbour_q ||
            block.total_gradient.size() != neighbour_q) {
            throw std::runtime_error(
                "mpmc::flow_discretization: local conservation neighbour block shape changed during assembly");
        }
        for (std::size_t column = 0U;
             column < neighbour_q;
             ++column) {
            const double total =
                block.total_gradient[column];
            double sum = 0.0;
            double scale = std::abs(total);
            if (!std::isfinite(total)) {
                throw std::range_error(
                    "mpmc::flow_discretization: local conservation neighbour total gradient became non-finite");
            }
            for (std::size_t component = 0U;
                 component < n;
                 ++component) {
                const double derivative =
                    block.component_jacobian[
                        component * neighbour_q +
                        column];
                if (!std::isfinite(derivative)) {
                    throw std::range_error(
                        "mpmc::flow_discretization: local conservation neighbour component Jacobian became non-finite");
                }
                sum += derivative;
                scale += std::abs(derivative);
            }
            if (!near_roundoff(sum, total, scale)) {
                throw std::runtime_error(
                    "mpmc::flow_discretization: local conservation neighbour component Jacobian does not close to total gradient");
            }
        }
    }
}

} // namespace local_component_conservation_detail

/// Assemble one local component-conservation row block from the already
/// time-discretized accumulation residual and zero or more already-normalized
/// internal-face contributions.
///
/// cell_bulk_volume_m3 is repeated deliberately: it binds this local row to an
/// explicit geometric cell volume and is checked against the selected side of
/// every incident face. The current rigid-grid derivative remains dV_b = 0 in
/// the upstream normalized-face contract.
[[nodiscard]] inline
LocalComponentConservationResidualLinearization3D
build_local_component_conservation_residual(
    const mpmc::flow::NaturalVariableStateIdentity3P&
        cell_state_identity,
    double cell_bulk_volume_m3,
    const mpmc::flow::
        BackwardEulerComponentAccumulationResidual3P&
            accumulation,
    std::span<
        const LocalCellNormalizedFaceContributionBinding3D>
        incident_faces) {
    using namespace local_component_conservation_detail;

    validate_state_identity(
        cell_state_identity,
        "local-cell");
    if (!std::isfinite(cell_bulk_volume_m3) ||
        !(cell_bulk_volume_m3 > 0.0)) {
        throw std::invalid_argument(
            "mpmc::flow_discretization: local cell bulk volume must be finite and strictly positive [m^3]");
    }
    validate_accumulation(
        cell_state_identity,
        accumulation);

    const std::size_t n =
        cell_state_identity.component_ids.size();
    const std::size_t local_q =
        cell_state_identity.layout.unknown_count();

    LocalComponentConservationResidualLinearization3D
        result{
            cell_state_identity,
            cell_bulk_volume_m3,
            accumulation.porosity,
            accumulation.time_step_seconds,
            accumulation.component_ids,
            accumulation
                .component_residual_mol_per_bulk_m3_s,
            accumulation
                .total_residual_mol_per_bulk_m3_s,
            accumulation.component_jacobian,
            accumulation.total_residual_gradient,
            {}};
    result.neighbour_blocks.reserve(
        incident_faces.size());

    for (const auto& binding : incident_faces) {
        if (binding.contribution == nullptr) {
            throw std::invalid_argument(
                "mpmc::flow_discretization: incident normalized face binding must not be null");
        }
        const auto& face = *binding.contribution;
        validate_normalized_face(face);

        for (const auto& existing :
             result.neighbour_blocks) {
            if (existing.face == face.face) {
                throw std::invalid_argument(
                    "mpmc::flow_discretization: duplicate incident face in local conservation assembly");
            }
        }

        bool local_is_owner = false;
        switch (binding.local_side) {
        case IncidentFaceLocalSide3D::owner:
            local_is_owner = true;
            break;
        case IncidentFaceLocalSide3D::neighbour:
            local_is_owner = false;
            break;
        default:
            throw std::invalid_argument(
                "mpmc::flow_discretization: unsupported incident-face local side");
        }

        const auto& local_identity =
            local_is_owner
                ? face.owner_state_identity
                : face.neighbour_state_identity;
        const auto& neighbour_identity =
            local_is_owner
                ? face.neighbour_state_identity
                : face.owner_state_identity;
        const double face_local_volume =
            local_is_owner
                ? face.bulk_volume.owner_bulk_volume_m3
                : face.bulk_volume.neighbour_bulk_volume_m3;

        if (!same_state_identity(
                local_identity,
                cell_state_identity)) {
            throw std::invalid_argument(
                "mpmc::flow_discretization: incident face selected side does not match local conservation state/chart");
        }
        if (!near_roundoff(
                face_local_volume,
                cell_bulk_volume_m3,
                std::max(
                    face_local_volume,
                    cell_bulk_volume_m3))) {
            throw std::invalid_argument(
                "mpmc::flow_discretization: incident face bulk volume does not match local conservation cell bulk volume");
        }
        if (face.component_ids !=
            result.component_ids) {
            throw std::invalid_argument(
                "mpmc::flow_discretization: incident face component identity/order mismatch");
        }

        const auto& local_values =
            local_is_owner
                ? face.owner_component_contribution_mol_per_bulk_m3_s
                : face.neighbour_component_contribution_mol_per_bulk_m3_s;
        const auto& local_local_jacobian =
            local_is_owner
                ? face.owner_row_owner_column_jacobian
                : face.neighbour_row_neighbour_column_jacobian;
        const auto& local_neighbour_jacobian =
            local_is_owner
                ? face.owner_row_neighbour_column_jacobian
                : face.neighbour_row_owner_column_jacobian;
        const double local_total =
            local_is_owner
                ? face.owner_total_contribution_mol_per_bulk_m3_s
                : face.neighbour_total_contribution_mol_per_bulk_m3_s;
        const auto& local_total_local_gradient =
            local_is_owner
                ? face.owner_total_row_owner_column_gradient
                : face.neighbour_total_row_neighbour_column_gradient;
        const auto& local_total_neighbour_gradient =
            local_is_owner
                ? face.owner_total_row_neighbour_column_gradient
                : face.neighbour_total_row_owner_column_gradient;

        const std::size_t neighbour_q =
            neighbour_identity.layout.unknown_count();
        if (local_q == 0U ||
            neighbour_q == 0U ||
            n > std::numeric_limits<std::size_t>::max() /
                    local_q ||
            n > std::numeric_limits<std::size_t>::max() /
                    neighbour_q ||
            local_local_jacobian.size() !=
                n * local_q ||
            local_total_local_gradient.size() !=
                local_q ||
            local_neighbour_jacobian.size() !=
                n * neighbour_q ||
            local_total_neighbour_gradient.size() !=
                neighbour_q) {
            throw std::invalid_argument(
                "mpmc::flow_discretization: selected incident-face Jacobian shape does not match local/neighbour charts");
        }

        for (std::size_t component = 0U;
             component < n;
             ++component) {
            result
                .component_residual_mol_per_bulk_m3_s[
                    component] +=
                local_values[component];
            if (!std::isfinite(
                    result
                        .component_residual_mol_per_bulk_m3_s[
                            component])) {
                throw std::range_error(
                    "mpmc::flow_discretization: local component conservation residual overflowed/non-finite during face accumulation");
            }

            for (std::size_t column = 0U;
                 column < local_q;
                 ++column) {
                const std::size_t index =
                    component * local_q +
                    column;
                result.local_component_jacobian[index] +=
                    local_local_jacobian[index];
                if (!std::isfinite(
                        result.local_component_jacobian[
                            index])) {
                    throw std::range_error(
                        "mpmc::flow_discretization: local conservation diagonal Jacobian overflowed/non-finite during face accumulation");
                }
            }
        }

        result.total_residual_mol_per_bulk_m3_s +=
            local_total;
        if (!std::isfinite(
                result.total_residual_mol_per_bulk_m3_s)) {
            throw std::range_error(
                "mpmc::flow_discretization: local total conservation residual overflowed/non-finite during face accumulation");
        }

        for (std::size_t column = 0U;
             column < local_q;
             ++column) {
            result.local_total_gradient[column] +=
                local_total_local_gradient[column];
            if (!std::isfinite(
                    result.local_total_gradient[column])) {
                throw std::range_error(
                    "mpmc::flow_discretization: local conservation total diagonal gradient overflowed/non-finite during face accumulation");
            }
        }

        result.neighbour_blocks.push_back(
            LocalComponentConservationNeighbourBlock3D{
                face.face,
                binding.local_side,
                neighbour_identity,
                local_neighbour_jacobian,
                local_total_neighbour_gradient});
    }

    validate_final_closure(result);
    return result;
}

} // namespace mpmc::flow_discretization

#endif // MPMC_FLOW_DISCRETIZATION_LOCAL_COMPONENT_CONSERVATION_RESIDUAL_HPP
