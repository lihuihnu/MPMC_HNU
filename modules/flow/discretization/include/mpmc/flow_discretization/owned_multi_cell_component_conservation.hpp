#ifndef MPMC_FLOW_DISCRETIZATION_OWNED_MULTI_CELL_COMPONENT_CONSERVATION_HPP
#define MPMC_FLOW_DISCRETIZATION_OWNED_MULTI_CELL_COMPONENT_CONSERVATION_HPP

#include <mpmc/flow_discretization/local_component_conservation_residual.hpp>
#include <mpmc/mesh/partition_snapshot.hpp>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <optional>
#include <span>
#include <stdexcept>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

namespace mpmc::flow_discretization {

inline constexpr std::string_view
    serial_owned_multi_cell_component_conservation_convention =
        "flow_discretization/serial-closed-owned-multi-cell-component-conservation/v1";

/// PETSc-free copy of the identity/topology fields carried by one authoritative
/// assembly row. The row is produced from schedule.assembly_rows(); no ghost
/// schedule row is accepted by this adapter.
struct AuthoritativeInternalConnectionIdentity3D {
    mpmc::mesh::LocalIndex face{
        mpmc::mesh::LocalIndex::value_type{0}};
    mpmc::mesh::GlobalEntityId face_global{
        mpmc::mesh::GlobalEntityId::value_type{0}};
    mpmc::mesh::LocalIndex owner_cell{
        mpmc::mesh::LocalIndex::value_type{0}};
    mpmc::mesh::GlobalEntityId owner_cell_global{
        mpmc::mesh::GlobalEntityId::value_type{0}};
    mpmc::mesh::LocalIndex neighbour_cell{
        mpmc::mesh::LocalIndex::value_type{0}};
    mpmc::mesh::GlobalEntityId neighbour_cell_global{
        mpmc::mesh::GlobalEntityId::value_type{0}};
};

/// Copy only authoritative rows from a schedule-like object.
///
/// The production ParallelOwnedConnectionSchedule3D satisfies this interface,
/// but this header intentionally does not include discretization_petsc or PETSc.
template <class Schedule>
[[nodiscard]] std::vector<AuthoritativeInternalConnectionIdentity3D>
copy_authoritative_internal_connection_rows(
    const Schedule& schedule) {
    if (schedule.rank_count() == 0U ||
        schedule.local_rank().value() >=
            schedule.rank_count()) {
        throw std::invalid_argument(
            "mpmc::flow_discretization: authoritative connection schedule has invalid rank metadata");
    }

    std::vector<AuthoritativeInternalConnectionIdentity3D>
        rows;
    rows.reserve(schedule.assembly_rows().size());

    for (const auto& row : schedule.assembly_rows()) {
        if (row.owner_cell == row.neighbour_cell ||
            row.owner_cell_global ==
                row.neighbour_cell_global) {
            throw std::invalid_argument(
                "mpmc::flow_discretization: authoritative connection row must reference two distinct cells");
        }

        for (const auto& existing : rows) {
            if (existing.face == row.face ||
                existing.face_global ==
                    row.face_global) {
                throw std::invalid_argument(
                    "mpmc::flow_discretization: duplicate authoritative face identity");
            }
        }

        rows.push_back(
            AuthoritativeInternalConnectionIdentity3D{
                row.face,
                row.face_global,
                row.owner_cell,
                row.owner_cell_global,
                row.neighbour_cell,
                row.neighbour_cell_global});
    }
    return rows;
}

struct OwnedCellAccumulationBinding3D {
    mpmc::mesh::LocalIndex cell{
        mpmc::mesh::LocalIndex::value_type{0}};
    mpmc::mesh::GlobalEntityId cell_global{
        mpmc::mesh::GlobalEntityId::value_type{0}};
    double bulk_volume_m3{};
    mpmc::flow::NaturalVariableStateIdentity3P
        state_identity;
    const mpmc::flow::
        BackwardEulerComponentAccumulationResidual3P*
            accumulation{};
};

struct AuthoritativeNormalizedFaceContributionBinding3D {
    mpmc::mesh::LocalIndex face{
        mpmc::mesh::LocalIndex::value_type{0}};
    mpmc::mesh::GlobalEntityId face_global{
        mpmc::mesh::GlobalEntityId::value_type{0}};
    const NormalizedComponentFaceContributionLinearization3D*
        contribution{};
};

/// Coalesced off-diagonal row block keyed by the stable neighbour cell.
///
/// Multiple authoritative faces between the same pair are accumulated into one
/// cell-pair block while preserving every contributing stable face identity.
struct OwnedCellPairComponentJacobianBlock3D {
    mpmc::mesh::LocalIndex column_cell{
        mpmc::mesh::LocalIndex::value_type{0}};
    mpmc::mesh::GlobalEntityId column_cell_global{
        mpmc::mesh::GlobalEntityId::value_type{0}};
    mpmc::flow::NaturalVariableStateIdentity3P
        column_state_identity;
    std::vector<mpmc::mesh::LocalIndex>
        contributing_faces;
    std::vector<mpmc::mesh::GlobalEntityId>
        contributing_face_global_ids;
    std::vector<double> component_jacobian;
    std::vector<double> total_gradient;

    [[nodiscard]] double d_component(
        std::size_t component,
        std::size_t column) const {
        const std::size_t q =
            column_state_identity.layout.unknown_count();
        if (component >=
                column_state_identity.component_ids.size() ||
            column >= q ||
            q == 0U ||
            component >
                (std::numeric_limits<std::size_t>::max() -
                 column) /
                    q) {
            throw std::out_of_range(
                "mpmc::flow_discretization: cell-pair Jacobian index out of range");
        }
        return component_jacobian.at(
            component * q + column);
    }

    [[nodiscard]] double d_total(
        std::size_t column) const {
        return total_gradient.at(column);
    }
};

struct OwnedCellComponentConservationRow3D {
    mpmc::mesh::LocalIndex cell{
        mpmc::mesh::LocalIndex::value_type{0}};
    mpmc::mesh::GlobalEntityId cell_global{
        mpmc::mesh::GlobalEntityId::value_type{0}};
    LocalComponentConservationResidualLinearization3D
        local_residual;
    std::vector<OwnedCellPairComponentJacobianBlock3D>
        off_diagonal_cell_pair_blocks;
};

/// Complete serial closed-owned patch snapshot.
///
/// v1 is intentionally serial. In distributed execution, the rank that owns an
/// authoritative face may differ from the owner rank of either endpoint cell.
/// A complete distributed owned row therefore requires a later owner-targeted
/// residual/Jacobian exchange; this serial baseline must not be reused to claim
/// such an exchange already exists.
class SerialOwnedMultiCellComponentConservationSnapshot3D {
public:
    static constexpr std::string_view convention =
        serial_owned_multi_cell_component_conservation_convention;

    SerialOwnedMultiCellComponentConservationSnapshot3D(
        std::size_t local_cell_count,
        std::size_t authoritative_face_count,
        std::vector<std::string> component_ids,
        std::vector<OwnedCellComponentConservationRow3D>
            rows,
        std::vector<double>
            volume_weighted_spatial_component_balance_mol_per_s,
        double volume_weighted_spatial_total_balance_mol_per_s)
        : local_cell_count_(local_cell_count),
          authoritative_face_count_(
              authoritative_face_count),
          component_ids_(std::move(component_ids)),
          rows_(std::move(rows)),
          row_by_local_cell_(local_cell_count),
          volume_weighted_spatial_component_balance_mol_per_s_(
              std::move(
                  volume_weighted_spatial_component_balance_mol_per_s)),
          volume_weighted_spatial_total_balance_mol_per_s_(
              volume_weighted_spatial_total_balance_mol_per_s) {
        validate_and_index();
    }

    [[nodiscard]] std::size_t
    local_cell_count() const noexcept {
        return local_cell_count_;
    }

    [[nodiscard]] std::size_t
    authoritative_face_count() const noexcept {
        return authoritative_face_count_;
    }

    [[nodiscard]] std::span<const std::string>
    component_ids() const noexcept {
        return component_ids_;
    }

    [[nodiscard]] std::span<
        const OwnedCellComponentConservationRow3D>
    rows() const noexcept {
        return rows_;
    }

    [[nodiscard]] const OwnedCellComponentConservationRow3D&
    row(mpmc::mesh::LocalIndex cell) const {
        const std::size_t local =
            static_cast<std::size_t>(cell.value());
        if (local >= local_cell_count_) {
            throw std::out_of_range(
                "mpmc::flow_discretization: owned multi-cell row index out of range");
        }
        const auto mapped =
            row_by_local_cell_[local];
        if (!mapped.has_value()) {
            throw std::invalid_argument(
                "mpmc::flow_discretization: owned multi-cell row is absent");
        }
        return rows_.at(*mapped);
    }

    [[nodiscard]] std::span<const double>
    volume_weighted_spatial_component_balance_mol_per_s()
        const noexcept {
        return volume_weighted_spatial_component_balance_mol_per_s_;
    }

    [[nodiscard]] double
    volume_weighted_spatial_total_balance_mol_per_s()
        const noexcept {
        return volume_weighted_spatial_total_balance_mol_per_s_;
    }

private:
    void validate_and_index() {
        if (component_ids_.size() < 2U ||
            rows_.size() != local_cell_count_ ||
            volume_weighted_spatial_component_balance_mol_per_s_
                    .size() !=
                component_ids_.size() ||
            !std::isfinite(
                volume_weighted_spatial_total_balance_mol_per_s_)) {
            throw std::invalid_argument(
                "mpmc::flow_discretization: malformed serial owned multi-cell snapshot");
        }

        for (double value :
             volume_weighted_spatial_component_balance_mol_per_s_) {
            if (!std::isfinite(value)) {
                throw std::invalid_argument(
                    "mpmc::flow_discretization: non-finite spatial conservation diagnostic");
            }
        }

        for (std::size_t index = 0U;
             index < rows_.size();
             ++index) {
            const auto& row = rows_[index];
            const std::size_t local =
                static_cast<std::size_t>(
                    row.cell.value());
            if (local >= local_cell_count_ ||
                row_by_local_cell_[local].has_value() ||
                row.local_residual.component_ids !=
                    component_ids_ ||
                row.local_residual.cell_state_identity
                        .component_ids !=
                    component_ids_) {
                throw std::invalid_argument(
                    "mpmc::flow_discretization: invalid or duplicate serial owned multi-cell row");
            }
            row_by_local_cell_[local] = index;
        }
    }

    std::size_t local_cell_count_{};
    std::size_t authoritative_face_count_{};
    std::vector<std::string> component_ids_;
    std::vector<OwnedCellComponentConservationRow3D>
        rows_;
    std::vector<std::optional<std::size_t>>
        row_by_local_cell_;
    std::vector<double>
        volume_weighted_spatial_component_balance_mol_per_s_;
    double volume_weighted_spatial_total_balance_mol_per_s_{};
};

namespace owned_multi_cell_conservation_detail {

using local_component_conservation_detail::near_roundoff;
using local_component_conservation_detail::same_state_identity;

[[nodiscard]] inline std::size_t local_position(
    mpmc::mesh::LocalIndex index,
    std::size_t count,
    const char* message) {
    const std::size_t value =
        static_cast<std::size_t>(index.value());
    if (value >= count) {
        throw std::out_of_range(message);
    }
    return value;
}

inline void add_cell_pair_block(
    OwnedCellComponentConservationRow3D& row,
    const LocalComponentConservationNeighbourBlock3D&
        source,
    const AuthoritativeInternalConnectionIdentity3D&
        connection,
    mpmc::mesh::LocalIndex column_cell,
    mpmc::mesh::GlobalEntityId column_cell_global) {
    const std::size_t n =
        row.local_residual.component_count();
    const std::size_t q =
        source.neighbour_state_identity.layout
            .unknown_count();
    if (q == 0U ||
        n > std::numeric_limits<std::size_t>::max() / q ||
        source.component_jacobian.size() != n * q ||
        source.total_gradient.size() != q) {
        throw std::invalid_argument(
            "mpmc::flow_discretization: malformed per-face off-diagonal Jacobian block");
    }

    auto found =
        std::find_if(
            row.off_diagonal_cell_pair_blocks.begin(),
            row.off_diagonal_cell_pair_blocks.end(),
            [column_cell_global](
                const OwnedCellPairComponentJacobianBlock3D&
                    candidate) {
                return candidate.column_cell_global ==
                    column_cell_global;
            });

    if (found ==
        row.off_diagonal_cell_pair_blocks.end()) {
        row.off_diagonal_cell_pair_blocks.push_back(
            OwnedCellPairComponentJacobianBlock3D{
                column_cell,
                column_cell_global,
                source.neighbour_state_identity,
                {connection.face},
                {connection.face_global},
                source.component_jacobian,
                source.total_gradient});
        return;
    }

    if (found->column_cell != column_cell ||
        !same_state_identity(
            found->column_state_identity,
            source.neighbour_state_identity) ||
        found->component_jacobian.size() !=
            source.component_jacobian.size() ||
        found->total_gradient.size() !=
            source.total_gradient.size()) {
        throw std::invalid_argument(
            "mpmc::flow_discretization: repeated stable cell pair changed column-cell identity/chart");
    }

    if (std::find(
            found->contributing_face_global_ids.begin(),
            found->contributing_face_global_ids.end(),
            connection.face_global) !=
        found->contributing_face_global_ids.end()) {
        throw std::invalid_argument(
            "mpmc::flow_discretization: duplicate face in coalesced cell-pair Jacobian block");
    }

    found->contributing_faces.push_back(
        connection.face);
    found->contributing_face_global_ids.push_back(
        connection.face_global);
    for (std::size_t index = 0U;
         index < found->component_jacobian.size();
         ++index) {
        found->component_jacobian[index] +=
            source.component_jacobian[index];
        if (!std::isfinite(
                found->component_jacobian[index])) {
            throw std::range_error(
                "mpmc::flow_discretization: coalesced cell-pair component Jacobian became non-finite");
        }
    }
    for (std::size_t column = 0U;
         column < found->total_gradient.size();
         ++column) {
        found->total_gradient[column] +=
            source.total_gradient[column];
        if (!std::isfinite(
                found->total_gradient[column])) {
            throw std::range_error(
                "mpmc::flow_discretization: coalesced cell-pair total gradient became non-finite");
        }
    }
}

[[nodiscard]] inline const OwnedCellPairComponentJacobianBlock3D*
find_cell_pair_block(
    const OwnedCellComponentConservationRow3D& row,
    mpmc::mesh::GlobalEntityId column_cell_global) {
    const auto found =
        std::find_if(
            row.off_diagonal_cell_pair_blocks.begin(),
            row.off_diagonal_cell_pair_blocks.end(),
            [column_cell_global](
                const OwnedCellPairComponentJacobianBlock3D&
                    block) {
                return block.column_cell_global ==
                    column_cell_global;
            });
    return found ==
               row.off_diagonal_cell_pair_blocks.end()
        ? nullptr
        : &*found;
}

} // namespace owned_multi_cell_conservation_detail

/// Assemble a complete serial closed-owned patch.
///
/// Schedule is intentionally a template so the flow-discretization public
/// contract does not include PETSc. The actual
/// discretization_petsc::ParallelOwnedConnectionSchedule3D is compiled against
/// this interface in its integration test.
template <class Schedule>
[[nodiscard]]
SerialOwnedMultiCellComponentConservationSnapshot3D
build_serial_owned_multi_cell_component_conservation_snapshot(
    const Schedule& schedule,
    const mpmc::mesh::PartitionSnapshot& partition,
    std::span<const OwnedCellAccumulationBinding3D>
        cell_bindings,
    std::span<
        const AuthoritativeNormalizedFaceContributionBinding3D>
        face_bindings) {
    using namespace owned_multi_cell_conservation_detail;

    if (!partition.is_serial() ||
        partition.local_rank().value() != 0U ||
        partition.rank_count() != 1U ||
        schedule.local_rank() != partition.local_rank() ||
        schedule.rank_count() != partition.rank_count()) {
        throw std::invalid_argument(
            "mpmc::flow_discretization: serial owned multi-cell conservation requires matching single-rank partition/schedule metadata");
    }

    const auto connections =
        copy_authoritative_internal_connection_rows(
            schedule);

    const std::size_t cell_count =
        partition.entity_count(
            mpmc::mesh::EntityKind::cell);
    const std::size_t face_count =
        partition.entity_count(
            mpmc::mesh::EntityKind::face);

    if (cell_bindings.size() != cell_count ||
        face_bindings.size() !=
            connections.size()) {
        throw std::invalid_argument(
            "mpmc::flow_discretization: serial owned multi-cell binding counts do not match cells/authoritative faces");
    }

    std::vector<std::optional<std::size_t>>
        cell_to_binding(cell_count);
    std::vector<std::optional<std::size_t>>
        face_to_binding(face_count);
    std::vector<std::optional<std::size_t>>
        face_to_connection(face_count);

    std::vector<std::string> component_ids;

    for (std::size_t index = 0U;
         index < cell_bindings.size();
         ++index) {
        const auto& binding = cell_bindings[index];
        const std::size_t cell =
            local_position(
                binding.cell,
                cell_count,
                "mpmc::flow_discretization: cell binding local index out of range");
        if (cell_to_binding[cell].has_value() ||
            binding.accumulation == nullptr ||
            partition.global_id(
                mpmc::mesh::EntityKind::cell,
                binding.cell) !=
                binding.cell_global ||
            !partition.is_owned(
                mpmc::mesh::EntityKind::cell,
                binding.cell) ||
            !std::isfinite(binding.bulk_volume_m3) ||
            !(binding.bulk_volume_m3 > 0.0)) {
            throw std::invalid_argument(
                "mpmc::flow_discretization: invalid or duplicate owned-cell accumulation binding");
        }

        const std::vector<std::string>
            ids{
                binding.state_identity
                    .component_ids.begin(),
                binding.state_identity
                    .component_ids.end()};
        if (index == 0U) {
            component_ids = ids;
        } else if (ids != component_ids) {
            throw std::invalid_argument(
                "mpmc::flow_discretization: owned cells must use one canonical component identity/order");
        }

        const std::array<
            LocalCellNormalizedFaceContributionBinding3D,
            0>
            no_faces{};
        (void)build_local_component_conservation_residual(
            binding.state_identity,
            binding.bulk_volume_m3,
            *binding.accumulation,
            no_faces);

        cell_to_binding[cell] = index;
    }

    for (const auto mapped : cell_to_binding) {
        if (!mapped.has_value()) {
            throw std::invalid_argument(
                "mpmc::flow_discretization: every serial local cell requires one accumulation binding");
        }
    }

    for (std::size_t index = 0U;
         index < face_bindings.size();
         ++index) {
        const auto& binding = face_bindings[index];
        const std::size_t face =
            local_position(
                binding.face,
                face_count,
                "mpmc::flow_discretization: face binding local index out of range");
        if (face_to_binding[face].has_value() ||
            binding.contribution == nullptr ||
            binding.contribution->face !=
                binding.face ||
            partition.global_id(
                mpmc::mesh::EntityKind::face,
                binding.face) !=
                binding.face_global ||
            !partition.is_owned(
                mpmc::mesh::EntityKind::face,
                binding.face)) {
            throw std::invalid_argument(
                "mpmc::flow_discretization: invalid or duplicate authoritative face-contribution binding");
        }
        local_component_conservation_detail::
            validate_normalized_face(
                *binding.contribution);
        face_to_binding[face] = index;
    }

    for (std::size_t index = 0U;
         index < connections.size();
         ++index) {
        const auto& connection =
            connections[index];
        const std::size_t face =
            local_position(
                connection.face,
                face_count,
                "mpmc::flow_discretization: authoritative schedule face local index out of range");
        const std::size_t owner =
            local_position(
                connection.owner_cell,
                cell_count,
                "mpmc::flow_discretization: authoritative schedule owner cell local index out of range");
        const std::size_t neighbour =
            local_position(
                connection.neighbour_cell,
                cell_count,
                "mpmc::flow_discretization: authoritative schedule neighbour cell local index out of range");

        if (face_to_connection[face].has_value() ||
            !partition.is_owned(
                mpmc::mesh::EntityKind::face,
                connection.face) ||
            partition.global_id(
                mpmc::mesh::EntityKind::face,
                connection.face) !=
                connection.face_global ||
            partition.global_id(
                mpmc::mesh::EntityKind::cell,
                connection.owner_cell) !=
                connection.owner_cell_global ||
            partition.global_id(
                mpmc::mesh::EntityKind::cell,
                connection.neighbour_cell) !=
                connection.neighbour_cell_global ||
            owner == neighbour ||
            !face_to_binding[face].has_value()) {
            throw std::invalid_argument(
                "mpmc::flow_discretization: authoritative schedule row does not match serial partition/face binding");
        }

        const auto& face_binding =
            face_bindings[*face_to_binding[face]];
        if (face_binding.face_global !=
            connection.face_global) {
            throw std::invalid_argument(
                "mpmc::flow_discretization: authoritative stable face ID mismatch");
        }
        face_to_connection[face] = index;
    }

    for (std::size_t face = 0U;
         face < face_count;
         ++face) {
        if (face_to_binding[face].has_value() !=
            face_to_connection[face].has_value()) {
            throw std::invalid_argument(
                "mpmc::flow_discretization: face bindings must match authoritative schedule rows exactly");
        }
    }

    std::vector<
        std::vector<
            LocalCellNormalizedFaceContributionBinding3D>>
        incident_faces(cell_count);

    for (const auto& connection : connections) {
        const std::size_t face =
            static_cast<std::size_t>(
                connection.face.value());
        const auto& binding =
            face_bindings[*face_to_binding[face]];

        incident_faces[
            static_cast<std::size_t>(
                connection.owner_cell.value())]
            .push_back(
                LocalCellNormalizedFaceContributionBinding3D{
                    IncidentFaceLocalSide3D::owner,
                    binding.contribution});
        incident_faces[
            static_cast<std::size_t>(
                connection.neighbour_cell.value())]
            .push_back(
                LocalCellNormalizedFaceContributionBinding3D{
                    IncidentFaceLocalSide3D::neighbour,
                    binding.contribution});
    }

    std::vector<OwnedCellComponentConservationRow3D>
        rows;
    rows.reserve(cell_count);

    for (std::size_t cell = 0U;
         cell < cell_count;
         ++cell) {
        const auto& binding =
            cell_bindings[*cell_to_binding[cell]];
        auto local =
            build_local_component_conservation_residual(
                binding.state_identity,
                binding.bulk_volume_m3,
                *binding.accumulation,
                std::span<
                    const LocalCellNormalizedFaceContributionBinding3D>{
                        incident_faces[cell]});

        OwnedCellComponentConservationRow3D row{
            binding.cell,
            binding.cell_global,
            std::move(local),
            {}};

        for (const auto& block :
             row.local_residual.neighbour_blocks) {
            const std::size_t face =
                static_cast<std::size_t>(
                    block.face.value());
            if (face >= face_to_connection.size() ||
                !face_to_connection[face].has_value()) {
                throw std::runtime_error(
                    "mpmc::flow_discretization: local residual references a non-authoritative face");
            }
            const auto& connection =
                connections[*face_to_connection[face]];

            mpmc::mesh::LocalIndex column_cell{
                mpmc::mesh::LocalIndex::value_type{0}};
            mpmc::mesh::GlobalEntityId column_global{
                mpmc::mesh::GlobalEntityId::value_type{0}};

            if (connection.owner_cell == binding.cell) {
                if (block.local_side !=
                    IncidentFaceLocalSide3D::owner) {
                    throw std::runtime_error(
                        "mpmc::flow_discretization: owner-row face side changed during multi-cell assembly");
                }
                column_cell =
                    connection.neighbour_cell;
                column_global =
                    connection.neighbour_cell_global;
            } else if (
                connection.neighbour_cell ==
                binding.cell) {
                if (block.local_side !=
                    IncidentFaceLocalSide3D::neighbour) {
                    throw std::runtime_error(
                        "mpmc::flow_discretization: neighbour-row face side changed during multi-cell assembly");
                }
                column_cell =
                    connection.owner_cell;
                column_global =
                    connection.owner_cell_global;
            } else {
                throw std::runtime_error(
                    "mpmc::flow_discretization: local row is not incident on authoritative connection");
            }

            const auto& column_binding =
                cell_bindings[
                    *cell_to_binding[
                        static_cast<std::size_t>(
                            column_cell.value())]];
            if (column_binding.cell_global !=
                    column_global ||
                !same_state_identity(
                    column_binding.state_identity,
                    block.neighbour_state_identity)) {
                throw std::invalid_argument(
                    "mpmc::flow_discretization: authoritative neighbour cell identity/chart mismatch");
            }

            add_cell_pair_block(
                row,
                block,
                connection,
                column_cell,
                column_global);
        }

        std::sort(
            row.off_diagonal_cell_pair_blocks.begin(),
            row.off_diagonal_cell_pair_blocks.end(),
            [](const auto& left, const auto& right) {
                return left.column_cell_global <
                    right.column_cell_global;
            });

        rows.push_back(std::move(row));
    }

    std::vector<double> weighted_component_balance(
        component_ids.size(),
        0.0);
    std::vector<double> weighted_component_scale(
        component_ids.size(),
        0.0);
    double weighted_total_balance = 0.0;
    double weighted_total_scale = 0.0;

    for (const auto& row : rows) {
        const auto& binding =
            cell_bindings[
                *cell_to_binding[
                    static_cast<std::size_t>(
                        row.cell.value())]];
        for (std::size_t component = 0U;
             component < component_ids.size();
             ++component) {
            const double spatial =
                row.local_residual.residual(component) -
                binding.accumulation->residual(component);
            const double weighted =
                binding.bulk_volume_m3 *
                spatial;
            weighted_component_balance[component] +=
                weighted;
            weighted_component_scale[component] +=
                std::abs(weighted);
        }
        const double weighted_total =
            binding.bulk_volume_m3 *
            (row.local_residual
                 .total_residual_mol_per_bulk_m3_s -
             binding.accumulation
                 ->total_residual_mol_per_bulk_m3_s);
        weighted_total_balance +=
            weighted_total;
        weighted_total_scale +=
            std::abs(weighted_total);
    }

    for (std::size_t component = 0U;
         component < weighted_component_balance.size();
         ++component) {
        if (!near_roundoff(
                weighted_component_balance[component],
                0.0,
                weighted_component_scale[component])) {
            throw std::runtime_error(
                "mpmc::flow_discretization: closed serial patch violates volume-weighted component spatial conservation");
        }
    }
    if (!near_roundoff(
            weighted_total_balance,
            0.0,
            weighted_total_scale)) {
        throw std::runtime_error(
            "mpmc::flow_discretization: closed serial patch violates volume-weighted total spatial conservation");
    }

    /// Global spatial Jacobian conservation: for every source-cell natural
    /// variable column, sum all target-row spatial derivatives after restoring
    /// each target cell's bulk volume.
    for (const auto& source_binding : cell_bindings) {
        const std::size_t q =
            source_binding.state_identity.layout
                .unknown_count();
        const auto& source_row =
            rows[
                static_cast<std::size_t>(
                    source_binding.cell.value())];

        for (std::size_t column = 0U;
             column < q;
             ++column) {
            std::vector<double> component_balance(
                component_ids.size(),
                0.0);
            std::vector<double> component_scale(
                component_ids.size(),
                0.0);
            double total_balance = 0.0;
            double total_scale = 0.0;

            for (const auto& target_row : rows) {
                const auto& target_binding =
                    cell_bindings[
                        *cell_to_binding[
                            static_cast<std::size_t>(
                                target_row.cell.value())]];

                if (target_row.cell ==
                    source_binding.cell) {
                    for (std::size_t component = 0U;
                         component < component_ids.size();
                         ++component) {
                        const double weighted =
                            target_binding.bulk_volume_m3 *
                            (source_row.local_residual
                                 .d_local(
                                     component,
                                     column) -
                             source_binding.accumulation
                                 ->d_residual(
                                     component,
                                     column));
                        component_balance[component] +=
                            weighted;
                        component_scale[component] +=
                            std::abs(weighted);
                    }
                    const double weighted_total =
                        target_binding.bulk_volume_m3 *
                        (source_row.local_residual
                             .d_total_local(column) -
                         source_binding.accumulation
                             ->d_total(column));
                    total_balance +=
                        weighted_total;
                    total_scale +=
                        std::abs(weighted_total);
                    continue;
                }

                const auto* block =
                    find_cell_pair_block(
                        target_row,
                        source_binding.cell_global);
                if (block == nullptr) {
                    continue;
                }

                if (block->column_cell !=
                        source_binding.cell ||
                    !same_state_identity(
                        block->column_state_identity,
                        source_binding.state_identity)) {
                    throw std::runtime_error(
                        "mpmc::flow_discretization: cell-pair column identity changed during global conservation audit");
                }

                for (std::size_t component = 0U;
                     component < component_ids.size();
                     ++component) {
                    const double weighted =
                        target_binding.bulk_volume_m3 *
                        block->d_component(
                            component,
                            column);
                    component_balance[component] +=
                        weighted;
                    component_scale[component] +=
                        std::abs(weighted);
                }
                const double weighted_total =
                    target_binding.bulk_volume_m3 *
                    block->d_total(column);
                total_balance +=
                    weighted_total;
                total_scale +=
                    std::abs(weighted_total);
            }

            for (std::size_t component = 0U;
                 component < component_balance.size();
                 ++component) {
                if (!near_roundoff(
                        component_balance[component],
                        0.0,
                        component_scale[component])) {
                    throw std::runtime_error(
                        "mpmc::flow_discretization: closed serial patch violates volume-weighted spatial Jacobian conservation");
                }
            }
            if (!near_roundoff(
                    total_balance,
                    0.0,
                    total_scale)) {
                throw std::runtime_error(
                    "mpmc::flow_discretization: closed serial patch violates volume-weighted total spatial Jacobian conservation");
            }
        }
    }

    return SerialOwnedMultiCellComponentConservationSnapshot3D{
        cell_count,
        connections.size(),
        std::move(component_ids),
        std::move(rows),
        std::move(weighted_component_balance),
        weighted_total_balance};
}

} // namespace mpmc::flow_discretization

#endif // MPMC_FLOW_DISCRETIZATION_OWNED_MULTI_CELL_COMPONENT_CONSERVATION_HPP
