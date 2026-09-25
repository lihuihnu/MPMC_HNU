#ifndef MPMC_FLOW_DISCRETIZATION_PETSC_DISTRIBUTED_ENERGY_CONSERVATION_HPP
#define MPMC_FLOW_DISCRETIZATION_PETSC_DISTRIBUTED_ENERGY_CONSERVATION_HPP

#include <mpmc/flow_discretization/local_energy_conservation_residual.hpp>
#include <mpmc/flow_discretization_petsc/distributed_component_conservation.hpp>

#include <petscsys.h>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <numeric>
#include <optional>
#include <span>
#include <stdexcept>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

namespace mpmc::flow_discretization_petsc {

inline constexpr std::string_view
    distributed_owned_energy_conservation_convention =
        "flow_discretization_petsc/owner-targeted-energy-conservation/v1";

struct DistributedEnergyCellStateBinding3D {
    mpmc::mesh::LocalIndex cell{
        mpmc::mesh::LocalIndex::value_type{0}};
    mpmc::mesh::GlobalEntityId cell_global{
        mpmc::mesh::GlobalEntityId::value_type{0}};
    double bulk_volume_m3{};
    mpmc::flow::NaturalVariableStateIdentity3P
        state_identity;
    const mpmc::flow::
        BackwardEulerEnergyAccumulationResidual3P*
            owned_accumulation{};
};

struct AuthoritativeNormalizedEnergyFaceBinding3D {
    mpmc::mesh::LocalIndex face{
        mpmc::mesh::LocalIndex::value_type{0}};
    mpmc::mesh::GlobalEntityId face_global{
        mpmc::mesh::GlobalEntityId::value_type{0}};
    const mpmc::flow_discretization::
        NormalizedEnergyFaceContributionLinearization3D*
            contribution{};
};

struct OwnedCellPairEnergyJacobianBlock3D {
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
    std::vector<double> gradient;

    [[nodiscard]] double d_residual(
        std::size_t column) const {
        return gradient.at(column);
    }
};

struct OwnedCellEnergyConservationRow3D {
    mpmc::mesh::LocalIndex cell{
        mpmc::mesh::LocalIndex::value_type{0}};
    mpmc::mesh::GlobalEntityId cell_global{
        mpmc::mesh::GlobalEntityId::value_type{0}};
    mpmc::flow_discretization::
        LocalEnergyConservationResidualLinearization3D
            local_residual;
    std::vector<OwnedCellPairEnergyJacobianBlock3D>
        off_diagonal_cell_pair_blocks;
};

class DistributedOwnedMultiCellEnergyConservationSnapshot3D {
public:
    static constexpr std::string_view convention =
        distributed_owned_energy_conservation_convention;

    DistributedOwnedMultiCellEnergyConservationSnapshot3D(
        mpmc::mesh::PartitionRank local_rank,
        std::uint32_t rank_count,
        std::size_t local_cell_count,
        std::size_t local_authoritative_face_count,
        std::size_t received_face_side_count,
        std::vector<OwnedCellEnergyConservationRow3D>
            rows,
        double global_volume_weighted_spatial_balance_w)
        : local_rank_(local_rank),
          rank_count_(rank_count),
          local_cell_count_(local_cell_count),
          local_authoritative_face_count_(
              local_authoritative_face_count),
          received_face_side_count_(
              received_face_side_count),
          rows_(std::move(rows)),
          row_by_local_cell_(local_cell_count),
          global_volume_weighted_spatial_balance_w_(
              global_volume_weighted_spatial_balance_w) {
        validate_and_index();
    }

    [[nodiscard]] mpmc::mesh::PartitionRank
    local_rank() const noexcept {
        return local_rank_;
    }

    [[nodiscard]] std::uint32_t
    rank_count() const noexcept {
        return rank_count_;
    }

    [[nodiscard]] std::size_t
    local_cell_count() const noexcept {
        return local_cell_count_;
    }

    [[nodiscard]] std::size_t
    local_authoritative_face_count() const noexcept {
        return local_authoritative_face_count_;
    }

    [[nodiscard]] std::size_t
    received_face_side_count() const noexcept {
        return received_face_side_count_;
    }

    [[nodiscard]] std::span<
        const OwnedCellEnergyConservationRow3D>
    owned_rows() const noexcept {
        return rows_;
    }

    [[nodiscard]] const OwnedCellEnergyConservationRow3D&
    owned_row(
        mpmc::mesh::LocalIndex cell) const {
        const std::size_t local =
            static_cast<std::size_t>(
                cell.value());
        if (local >= local_cell_count_) {
            throw std::out_of_range(
                "mpmc::flow_discretization_petsc: owned energy row local cell index out of range");
        }
        const auto mapped =
            row_by_local_cell_[local];
        if (!mapped.has_value()) {
            throw std::invalid_argument(
                "mpmc::flow_discretization_petsc: requested energy cell is not locally owned");
        }
        return rows_.at(*mapped);
    }

    [[nodiscard]] double
    global_volume_weighted_spatial_balance_w()
        const noexcept {
        return global_volume_weighted_spatial_balance_w_;
    }

private:
    void validate_and_index() {
        if (rank_count_ == 0U ||
            local_rank_.value() >= rank_count_ ||
            !std::isfinite(
                global_volume_weighted_spatial_balance_w_)) {
            throw std::invalid_argument(
                "mpmc::flow_discretization_petsc: malformed distributed owned energy-conservation snapshot");
        }

        for (std::size_t index = 0U;
             index < rows_.size();
             ++index) {
            const auto& row =
                rows_[index];
            const std::size_t local =
                static_cast<std::size_t>(
                    row.cell.value());
            if (local >= local_cell_count_ ||
                row_by_local_cell_[local].has_value()) {
                throw std::invalid_argument(
                    "mpmc::flow_discretization_petsc: invalid or duplicate distributed owned energy row");
            }
            row_by_local_cell_[local] =
                index;
        }
    }

    mpmc::mesh::PartitionRank local_rank_;
    std::uint32_t rank_count_;
    std::size_t local_cell_count_{};
    std::size_t local_authoritative_face_count_{};
    std::size_t received_face_side_count_{};
    std::vector<OwnedCellEnergyConservationRow3D>
        rows_;
    std::vector<std::optional<std::size_t>>
        row_by_local_cell_;
    double global_volume_weighted_spatial_balance_w_{};
};

namespace distributed_energy_detail {

inline constexpr std::size_t
    endpoint_metadata_word_count = 6U;

struct EndpointEnergyPayload {
    std::uint64_t target_cell_global{};
    std::uint64_t column_cell_global{};
    std::uint64_t face_global{};
    std::uint64_t target_state_fingerprint{};
    std::uint64_t column_state_fingerprint{};
    std::uint64_t local_side{};

    double target_bulk_volume_m3{};
    double spatial_w_per_bulk_m3{};
    std::vector<double>
        diagonal_gradient;
    std::vector<double>
        off_diagonal_gradient;
};

[[nodiscard]] inline std::size_t payload_value_count(
    std::size_t column_count) {
    if (column_count == 0U ||
        column_count >
            (std::numeric_limits<std::size_t>::max() -
             2U) /
                2U) {
        throw std::length_error(
            "mpmc::flow_discretization_petsc: energy endpoint payload size overflow");
    }
    return 2U +
        2U * column_count;
}

[[nodiscard]] inline EndpointEnergyPayload
make_endpoint_payload(
    const mpmc::discretization_petsc::
        AssemblyReadyInternalConnectionRow3D&
            connection,
    const mpmc::flow_discretization::
        NormalizedEnergyFaceContributionLinearization3D&
            face,
    bool target_is_owner) {
    using namespace mpmc::flow_discretization;

    EndpointEnergyPayload payload;
    payload.face_global =
        connection.face_global.value();

    if (target_is_owner) {
        payload.target_cell_global =
            connection.owner_cell_global.value();
        payload.column_cell_global =
            connection.neighbour_cell_global.value();
        payload.target_state_fingerprint =
            detail::state_identity_fingerprint(
                face.owner_state_identity);
        payload.column_state_fingerprint =
            detail::state_identity_fingerprint(
                face.neighbour_state_identity);
        payload.local_side =
            static_cast<std::uint64_t>(
                EnergyIncidentFaceLocalSide3D::owner);
        payload.target_bulk_volume_m3 =
            face.bulk_volume.owner_bulk_volume_m3;
        payload.spatial_w_per_bulk_m3 =
            face.owner_contribution_w_per_bulk_m3;
        payload.diagonal_gradient =
            face.owner_row_owner_column_gradient;
        payload.off_diagonal_gradient =
            face.owner_row_neighbour_column_gradient;
    } else {
        payload.target_cell_global =
            connection.neighbour_cell_global.value();
        payload.column_cell_global =
            connection.owner_cell_global.value();
        payload.target_state_fingerprint =
            detail::state_identity_fingerprint(
                face.neighbour_state_identity);
        payload.column_state_fingerprint =
            detail::state_identity_fingerprint(
                face.owner_state_identity);
        payload.local_side =
            static_cast<std::uint64_t>(
                EnergyIncidentFaceLocalSide3D::neighbour);
        payload.target_bulk_volume_m3 =
            face.bulk_volume.neighbour_bulk_volume_m3;
        payload.spatial_w_per_bulk_m3 =
            face.neighbour_contribution_w_per_bulk_m3;
        payload.diagonal_gradient =
            face.neighbour_row_neighbour_column_gradient;
        payload.off_diagonal_gradient =
            face.neighbour_row_owner_column_gradient;
    }
    return payload;
}

inline void append_payload_values(
    const EndpointEnergyPayload& payload,
    std::vector<double>* values) {
    if (values == nullptr) {
        throw std::invalid_argument(
            "mpmc::flow_discretization_petsc: null energy payload value buffer");
    }
    values->push_back(
        payload.target_bulk_volume_m3);
    values->push_back(
        payload.spatial_w_per_bulk_m3);
    values->insert(
        values->end(),
        payload.diagonal_gradient.begin(),
        payload.diagonal_gradient.end());
    values->insert(
        values->end(),
        payload.off_diagonal_gradient.begin(),
        payload.off_diagonal_gradient.end());
}

inline void add_cell_pair_block(
    OwnedCellEnergyConservationRow3D& row,
    const mpmc::flow_discretization::
        LocalEnergyConservationNeighbourBlock3D&
            source,
    mpmc::mesh::LocalIndex column_cell,
    mpmc::mesh::GlobalEntityId column_cell_global,
    mpmc::mesh::GlobalEntityId face_global) {
    const std::size_t q =
        source.neighbour_state_identity.layout
            .unknown_count();
    if (q == 0U ||
        source.gradient.size() != q) {
        throw std::invalid_argument(
            "mpmc::flow_discretization_petsc: malformed per-face energy off-diagonal Jacobian block");
    }

    auto found =
        std::find_if(
            row.off_diagonal_cell_pair_blocks.begin(),
            row.off_diagonal_cell_pair_blocks.end(),
            [column_cell_global](
                const OwnedCellPairEnergyJacobianBlock3D&
                    candidate) {
                return candidate.column_cell_global ==
                    column_cell_global;
            });

    if (found ==
        row.off_diagonal_cell_pair_blocks.end()) {
        row.off_diagonal_cell_pair_blocks.push_back(
            OwnedCellPairEnergyJacobianBlock3D{
                column_cell,
                column_cell_global,
                source.neighbour_state_identity,
                {source.face},
                {face_global},
                source.gradient});
        return;
    }

    if (found->column_cell !=
            column_cell ||
        !mpmc::flow::
            energy_accumulation_detail::
                same_state_identity(
                    found->column_state_identity,
                    source.neighbour_state_identity) ||
        found->gradient.size() !=
            source.gradient.size() ||
        std::find(
            found->contributing_face_global_ids.begin(),
            found->contributing_face_global_ids.end(),
            face_global) !=
            found->contributing_face_global_ids.end()) {
        throw std::invalid_argument(
            "mpmc::flow_discretization_petsc: repeated energy cell pair changed identity/chart or duplicated face");
    }

    found->contributing_faces.push_back(
        source.face);
    found->contributing_face_global_ids.push_back(
        face_global);
    for (std::size_t column = 0U;
         column < q;
         ++column) {
        found->gradient[column] +=
            source.gradient[column];
        if (!std::isfinite(
                found->gradient[column])) {
            throw std::range_error(
                "mpmc::flow_discretization_petsc: coalesced energy cell-pair Jacobian became non-finite");
        }
    }
}

[[nodiscard]] inline const OwnedCellPairEnergyJacobianBlock3D*
find_cell_pair_block(
    const OwnedCellEnergyConservationRow3D& row,
    mpmc::mesh::GlobalEntityId column_cell_global) {
    const auto found =
        std::find_if(
            row.off_diagonal_cell_pair_blocks.begin(),
            row.off_diagonal_cell_pair_blocks.end(),
            [column_cell_global](
                const OwnedCellPairEnergyJacobianBlock3D&
                    block) {
                return block.column_cell_global ==
                    column_cell_global;
            });
    return found ==
               row.off_diagonal_cell_pair_blocks.end()
        ? nullptr
        : &*found;
}

} // namespace distributed_energy_detail

/// Route each authoritative internal-face energy contribution to the owner rank
/// of its endpoint cell and assemble complete owned energy rows.
///
/// This layer performs MPI exchange only. It does not assign PETSc scalar
/// equation numbers and does not create/modify Mat or Vec objects.
inline PetscErrorCode
make_distributed_owned_energy_conservation_snapshot_3d(
    MPI_Comm comm,
    const mpmc::discretization_petsc::
        ParallelOwnedConnectionSchedule3D&
            schedule,
    const mpmc::mesh::PartitionSnapshot& partition,
    std::span<
        const DistributedEnergyCellStateBinding3D>
        cell_bindings,
    std::span<
        const AuthoritativeNormalizedEnergyFaceBinding3D>
        authoritative_face_bindings,
    std::optional<
        DistributedOwnedMultiCellEnergyConservationSnapshot3D>*
            output) {
    using namespace distributed_energy_detail;
    using namespace mpmc::flow_discretization;

    int mpi_rank = -1;
    int mpi_size = -1;
    if (MPI_Comm_rank(
            comm,
            &mpi_rank) != MPI_SUCCESS ||
        MPI_Comm_size(
            comm,
            &mpi_size) != MPI_SUCCESS ||
        mpi_rank < 0 ||
        mpi_size <= 0) {
        return PETSC_ERR_MPI;
    }

    PetscErrorCode local_error =
        output == nullptr
            ? PETSC_ERR_ARG_NULL
            : PETSC_SUCCESS;
    PetscErrorCode error =
        detail::collective_error(
            comm,
            local_error);
    if (error != PETSC_SUCCESS) {
        return error;
    }
    output->reset();

    std::vector<std::optional<std::size_t>>
        cell_to_binding;
    std::vector<std::optional<std::size_t>>
        face_to_binding;
    std::vector<std::string> component_ids;
    std::size_t column_count = 0U;
    std::size_t value_width = 0U;
    std::vector<
        std::vector<EndpointEnergyPayload>>
        outbound;

    try {
        if (schedule.local_rank() !=
                partition.local_rank() ||
            schedule.rank_count() !=
                partition.rank_count() ||
            schedule.local_rank().value() !=
                static_cast<std::uint32_t>(
                    mpi_rank) ||
            schedule.rank_count() !=
                static_cast<std::uint32_t>(
                    mpi_size)) {
            throw std::invalid_argument(
                "energy exchange rank metadata mismatch");
        }

        const std::size_t local_cell_count =
            partition.entity_count(
                mpmc::mesh::EntityKind::cell);
        const std::size_t local_face_count =
            partition.entity_count(
                mpmc::mesh::EntityKind::face);
        if (local_cell_count == 0U ||
            cell_bindings.size() !=
                local_cell_count ||
            authoritative_face_bindings.size() !=
                schedule.assembly_rows().size()) {
            throw std::invalid_argument(
                "energy exchange binding cardinality mismatch");
        }

        cell_to_binding.assign(
            local_cell_count,
            std::nullopt);
        face_to_binding.assign(
            local_face_count,
            std::nullopt);

        const std::array<
            LocalCellNormalizedEnergyFaceBinding3D,
            0>
            no_faces{};

        for (std::size_t index = 0U;
             index < cell_bindings.size();
             ++index) {
            const auto& binding =
                cell_bindings[index];
            const std::size_t local =
                static_cast<std::size_t>(
                    binding.cell.value());
            if (local >= local_cell_count ||
                cell_to_binding[local].has_value() ||
                partition.global_id(
                    mpmc::mesh::EntityKind::cell,
                    binding.cell) !=
                    binding.cell_global ||
                !std::isfinite(
                    binding.bulk_volume_m3) ||
                !(binding.bulk_volume_m3 > 0.0)) {
                throw std::invalid_argument(
                    "invalid or duplicate local energy cell binding");
            }

            mpmc::flow_discretization::
                local_component_conservation_detail::
                    validate_state_identity(
                        binding.state_identity,
                        "distributed-energy-cell");

            if (index == 0U) {
                component_ids =
                    binding.state_identity.component_ids;
                column_count =
                    binding.state_identity.layout
                        .unknown_count();
                value_width =
                    payload_value_count(
                        column_count);
            } else if (
                binding.state_identity.component_ids !=
                    component_ids ||
                binding.state_identity.layout
                        .unknown_count() !=
                    column_count) {
                throw std::invalid_argument(
                    "local energy cells do not share canonical component identity/count");
            }

            const bool owned =
                partition.is_owned(
                    mpmc::mesh::EntityKind::cell,
                    binding.cell);
            if (owned) {
                if (binding.owned_accumulation ==
                        nullptr ||
                    !mpmc::flow::
                        energy_accumulation_detail::
                            same_state_identity(
                                binding.state_identity,
                                binding.owned_accumulation
                                    ->current_state_identity) ||
                    !std::isfinite(
                        binding.owned_accumulation
                            ->porosity) ||
                    !(binding.owned_accumulation
                          ->porosity > 0.0) ||
                    !(binding.owned_accumulation
                          ->porosity < 1.0) ||
                    !std::isfinite(
                        binding.owned_accumulation
                            ->time_step_seconds) ||
                    !(binding.owned_accumulation
                          ->time_step_seconds > 0.0)) {
                    throw std::invalid_argument(
                        "owned energy cell is missing/mismatching backward-Euler accumulation");
                }
                (void)build_local_energy_conservation_residual(
                    *binding.owned_accumulation,
                    binding.bulk_volume_m3,
                    no_faces);
            } else if (
                !partition.is_ghost(
                    mpmc::mesh::EntityKind::cell,
                    binding.cell) ||
                binding.owned_accumulation !=
                    nullptr) {
                throw std::invalid_argument(
                    "ghost energy cell must not carry owned accumulation");
            }

            cell_to_binding[local] =
                index;
        }

        for (const auto mapped :
             cell_to_binding) {
            if (!mapped.has_value()) {
                throw std::invalid_argument(
                    "every local energy cell requires one state binding");
            }
        }

        for (std::size_t index = 0U;
             index <
             authoritative_face_bindings.size();
             ++index) {
            const auto& binding =
                authoritative_face_bindings[index];
            const std::size_t face =
                static_cast<std::size_t>(
                    binding.face.value());
            if (face >= local_face_count ||
                face_to_binding[face].has_value() ||
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
                    "invalid or duplicate authoritative energy face binding");
            }
            local_energy_conservation_detail::
                validate_face(
                    *binding.contribution);
            face_to_binding[face] =
                index;
        }

        outbound.assign(
            static_cast<std::size_t>(
                mpi_size),
            {});

        for (const auto& connection :
             schedule.assembly_rows()) {
            const std::size_t face =
                static_cast<std::size_t>(
                    connection.face.value());
            const std::size_t owner =
                static_cast<std::size_t>(
                    connection.owner_cell.value());
            const std::size_t neighbour =
                static_cast<std::size_t>(
                    connection.neighbour_cell.value());

            if (face >= face_to_binding.size() ||
                owner >= cell_to_binding.size() ||
                neighbour >= cell_to_binding.size() ||
                !face_to_binding[face].has_value() ||
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
                !partition.is_owned(
                    mpmc::mesh::EntityKind::face,
                    connection.face)) {
                throw std::invalid_argument(
                    "authoritative energy schedule row does not match local partition/binding");
            }

            const auto& face_binding =
                authoritative_face_bindings[
                    *face_to_binding[face]];
            const auto& contribution =
                *face_binding.contribution;
            const auto& owner_binding =
                cell_bindings[
                    *cell_to_binding[owner]];
            const auto& neighbour_binding =
                cell_bindings[
                    *cell_to_binding[neighbour]];

            if (face_binding.face_global !=
                    connection.face_global ||
                owner_binding.cell_global !=
                    connection.owner_cell_global ||
                neighbour_binding.cell_global !=
                    connection.neighbour_cell_global ||
                detail::state_identity_fingerprint(
                    owner_binding.state_identity) !=
                    detail::state_identity_fingerprint(
                        contribution
                            .owner_state_identity) ||
                detail::state_identity_fingerprint(
                    neighbour_binding.state_identity) !=
                    detail::state_identity_fingerprint(
                        contribution
                            .neighbour_state_identity) ||
                !detail::near_roundoff(
                    owner_binding.bulk_volume_m3,
                    contribution.bulk_volume
                        .owner_bulk_volume_m3,
                    owner_binding.bulk_volume_m3) ||
                !detail::near_roundoff(
                    neighbour_binding.bulk_volume_m3,
                    contribution.bulk_volume
                        .neighbour_bulk_volume_m3,
                    neighbour_binding.bulk_volume_m3)) {
                throw std::invalid_argument(
                    "authoritative energy face endpoint state/chart/volume mismatch");
            }

            const auto owner_rank =
                partition.owner_rank(
                    mpmc::mesh::EntityKind::cell,
                    connection.owner_cell);
            const auto neighbour_rank =
                partition.owner_rank(
                    mpmc::mesh::EntityKind::cell,
                    connection.neighbour_cell);

            outbound[
                static_cast<std::size_t>(
                    owner_rank.value())]
                .push_back(
                    make_endpoint_payload(
                        connection,
                        contribution,
                        true));
            outbound[
                static_cast<std::size_t>(
                    neighbour_rank.value())]
                .push_back(
                    make_endpoint_payload(
                        connection,
                        contribution,
                        false));
        }
    } catch (...) {
        local_error =
            PETSC_ERR_ARG_INCOMP;
    }

    error =
        detail::collective_error(
            comm,
            local_error);
    if (error != PETSC_SUCCESS) {
        return error;
    }

    std::uint64_t local_component_hash =
        detail::component_identity_fingerprint(
            component_ids);
    std::uint64_t min_component_hash = 0U;
    std::uint64_t max_component_hash = 0U;
    std::uint64_t local_q =
        static_cast<std::uint64_t>(
            column_count);
    std::uint64_t min_q = 0U;
    std::uint64_t max_q = 0U;

    if (MPI_Allreduce(
            &local_component_hash,
            &min_component_hash,
            1,
            MPI_UINT64_T,
            MPI_MIN,
            comm) != MPI_SUCCESS ||
        MPI_Allreduce(
            &local_component_hash,
            &max_component_hash,
            1,
            MPI_UINT64_T,
            MPI_MAX,
            comm) != MPI_SUCCESS ||
        MPI_Allreduce(
            &local_q,
            &min_q,
            1,
            MPI_UINT64_T,
            MPI_MIN,
            comm) != MPI_SUCCESS ||
        MPI_Allreduce(
            &local_q,
            &max_q,
            1,
            MPI_UINT64_T,
            MPI_MAX,
            comm) != MPI_SUCCESS) {
        return PETSC_ERR_MPI;
    }
    if (min_component_hash !=
            max_component_hash ||
        min_q != max_q ||
        min_q == 0U) {
        return PETSC_ERR_ARG_INCOMP;
    }

    const std::size_t rank_count =
        static_cast<std::size_t>(
            mpi_size);
    std::vector<int> send_payload_counts(
        rank_count,
        0);
    for (std::size_t rank = 0U;
         rank < rank_count;
         ++rank) {
        if (outbound[rank].size() >
            static_cast<std::size_t>(
                std::numeric_limits<int>::max())) {
            local_error =
                PETSC_ERR_ARG_OUTOFRANGE;
            break;
        }
        send_payload_counts[rank] =
            static_cast<int>(
                outbound[rank].size());
    }

    error =
        detail::collective_error(
            comm,
            local_error);
    if (error != PETSC_SUCCESS) {
        return error;
    }

    std::vector<int> receive_payload_counts(
        rank_count,
        0);
    if (MPI_Alltoall(
            send_payload_counts.data(),
            1,
            MPI_INT,
            receive_payload_counts.data(),
            1,
            MPI_INT,
            comm) != MPI_SUCCESS) {
        return PETSC_ERR_MPI;
    }

    std::vector<int> send_metadata_counts(
        rank_count,
        0);
    std::vector<int> receive_metadata_counts(
        rank_count,
        0);
    std::vector<int> send_value_counts(
        rank_count,
        0);
    std::vector<int> receive_value_counts(
        rank_count,
        0);

    for (std::size_t rank = 0U;
         rank < rank_count;
         ++rank) {
        if (!detail::checked_int_product(
                static_cast<std::size_t>(
                    send_payload_counts[rank]),
                endpoint_metadata_word_count,
                &send_metadata_counts[rank]) ||
            !detail::checked_int_product(
                static_cast<std::size_t>(
                    receive_payload_counts[rank]),
                endpoint_metadata_word_count,
                &receive_metadata_counts[rank]) ||
            !detail::checked_int_product(
                static_cast<std::size_t>(
                    send_payload_counts[rank]),
                value_width,
                &send_value_counts[rank]) ||
            !detail::checked_int_product(
                static_cast<std::size_t>(
                    receive_payload_counts[rank]),
                value_width,
                &receive_value_counts[rank])) {
            local_error =
                PETSC_ERR_ARG_OUTOFRANGE;
            break;
        }
    }

    std::vector<int> send_metadata_displacements;
    std::vector<int> receive_metadata_displacements;
    std::vector<int> send_value_displacements;
    std::vector<int> receive_value_displacements;
    int total_send_metadata = 0;
    int total_receive_metadata = 0;
    int total_send_values = 0;
    int total_receive_values = 0;

    if (local_error == PETSC_SUCCESS &&
        (!detail::build_displacements(
             send_metadata_counts,
             &send_metadata_displacements,
             &total_send_metadata) ||
         !detail::build_displacements(
             receive_metadata_counts,
             &receive_metadata_displacements,
             &total_receive_metadata) ||
         !detail::build_displacements(
             send_value_counts,
             &send_value_displacements,
             &total_send_values) ||
         !detail::build_displacements(
             receive_value_counts,
             &receive_value_displacements,
             &total_receive_values))) {
        local_error =
            PETSC_ERR_ARG_OUTOFRANGE;
    }

    error =
        detail::collective_error(
            comm,
            local_error);
    if (error != PETSC_SUCCESS) {
        return error;
    }

    std::vector<std::uint64_t> send_metadata;
    std::vector<double> send_values;
    send_metadata.reserve(
        static_cast<std::size_t>(
            total_send_metadata));
    send_values.reserve(
        static_cast<std::size_t>(
            total_send_values));

    try {
        for (std::size_t rank = 0U;
             rank < rank_count;
             ++rank) {
            for (const auto& payload :
                 outbound[rank]) {
                if (payload.diagonal_gradient.size() !=
                        column_count ||
                    payload.off_diagonal_gradient.size() !=
                        column_count) {
                    throw std::invalid_argument(
                        "energy endpoint payload shape mismatch");
                }

                send_metadata.push_back(
                    payload.target_cell_global);
                send_metadata.push_back(
                    payload.column_cell_global);
                send_metadata.push_back(
                    payload.face_global);
                send_metadata.push_back(
                    payload.target_state_fingerprint);
                send_metadata.push_back(
                    payload.column_state_fingerprint);
                send_metadata.push_back(
                    payload.local_side);
                append_payload_values(
                    payload,
                    &send_values);
            }
        }
    } catch (...) {
        local_error =
            PETSC_ERR_ARG_INCOMP;
    }

    if (send_metadata.size() !=
            static_cast<std::size_t>(
                total_send_metadata) ||
        send_values.size() !=
            static_cast<std::size_t>(
                total_send_values)) {
        local_error =
            PETSC_ERR_PLIB;
    }

    error =
        detail::collective_error(
            comm,
            local_error);
    if (error != PETSC_SUCCESS) {
        return error;
    }

    std::vector<std::uint64_t> receive_metadata(
        static_cast<std::size_t>(
            total_receive_metadata));
    std::vector<double> receive_values(
        static_cast<std::size_t>(
            total_receive_values));

    if (MPI_Alltoallv(
            send_metadata.empty()
                ? nullptr
                : send_metadata.data(),
            send_metadata_counts.data(),
            send_metadata_displacements.data(),
            MPI_UINT64_T,
            receive_metadata.empty()
                ? nullptr
                : receive_metadata.data(),
            receive_metadata_counts.data(),
            receive_metadata_displacements.data(),
            MPI_UINT64_T,
            comm) != MPI_SUCCESS ||
        MPI_Alltoallv(
            send_values.empty()
                ? nullptr
                : send_values.data(),
            send_value_counts.data(),
            send_value_displacements.data(),
            MPI_DOUBLE,
            receive_values.empty()
                ? nullptr
                : receive_values.data(),
            receive_value_counts.data(),
            receive_value_displacements.data(),
            MPI_DOUBLE,
            comm) != MPI_SUCCESS) {
        return PETSC_ERR_MPI;
    }

    std::vector<OwnedCellEnergyConservationRow3D>
        rows;
    std::vector<std::optional<std::size_t>>
        row_by_local_cell(
            partition.entity_count(
                mpmc::mesh::EntityKind::cell));

    try {
        const std::array<
            LocalCellNormalizedEnergyFaceBinding3D,
            0>
            no_faces{};

        rows.reserve(
            partition.owned_count(
                mpmc::mesh::EntityKind::cell));

        for (const auto& binding :
             cell_bindings) {
            if (!partition.is_owned(
                    mpmc::mesh::EntityKind::cell,
                    binding.cell)) {
                continue;
            }

            auto local =
                build_local_energy_conservation_residual(
                    *binding.owned_accumulation,
                    binding.bulk_volume_m3,
                    no_faces);
            const std::size_t local_index =
                static_cast<std::size_t>(
                    binding.cell.value());
            row_by_local_cell[local_index] =
                rows.size();
            rows.push_back(
                OwnedCellEnergyConservationRow3D{
                    binding.cell,
                    binding.cell_global,
                    std::move(local),
                    {}});
        }

        std::vector<
            std::pair<std::uint64_t, std::uint64_t>>
            received_target_faces;

        std::size_t parsed_payloads = 0U;
        for (std::size_t source_rank = 0U;
             source_rank < rank_count;
             ++source_rank) {
            const int payload_count =
                receive_payload_counts[
                    source_rank];
            for (int record = 0;
                 record < payload_count;
                 ++record) {
                const std::size_t metadata_offset =
                    static_cast<std::size_t>(
                        receive_metadata_displacements[
                            source_rank]) +
                    static_cast<std::size_t>(
                        record) *
                        endpoint_metadata_word_count;
                const std::size_t value_offset =
                    static_cast<std::size_t>(
                        receive_value_displacements[
                            source_rank]) +
                    static_cast<std::size_t>(
                        record) *
                        value_width;

                const auto target_global =
                    mpmc::mesh::GlobalEntityId{
                        receive_metadata[
                            metadata_offset]};
                const auto column_global =
                    mpmc::mesh::GlobalEntityId{
                        receive_metadata[
                            metadata_offset + 1U]};
                const auto face_global =
                    mpmc::mesh::GlobalEntityId{
                        receive_metadata[
                            metadata_offset + 2U]};
                const std::uint64_t target_hash =
                    receive_metadata[
                        metadata_offset + 3U];
                const std::uint64_t column_hash =
                    receive_metadata[
                        metadata_offset + 4U];
                const std::uint64_t side_word =
                    receive_metadata[
                        metadata_offset + 5U];

                if (!partition.contains_global(
                        mpmc::mesh::EntityKind::cell,
                        target_global) ||
                    !partition.contains_global(
                        mpmc::mesh::EntityKind::cell,
                        column_global) ||
                    !partition.contains_global(
                        mpmc::mesh::EntityKind::face,
                        face_global)) {
                    throw std::invalid_argument(
                        "received energy endpoint stable identity is absent from local overlap");
                }

                const auto target_cell =
                    partition.local_index(
                        mpmc::mesh::EntityKind::cell,
                        target_global);
                const auto column_cell =
                    partition.local_index(
                        mpmc::mesh::EntityKind::cell,
                        column_global);
                const auto local_face =
                    partition.local_index(
                        mpmc::mesh::EntityKind::face,
                        face_global);

                if (!partition.is_owned(
                        mpmc::mesh::EntityKind::cell,
                        target_cell)) {
                    throw std::invalid_argument(
                        "received energy endpoint payload was routed to a non-owner rank");
                }

                const auto target_binding_index =
                    cell_to_binding[
                        static_cast<std::size_t>(
                            target_cell.value())];
                const auto column_binding_index =
                    cell_to_binding[
                        static_cast<std::size_t>(
                            column_cell.value())];
                if (!target_binding_index.has_value() ||
                    !column_binding_index.has_value()) {
                    throw std::invalid_argument(
                        "received energy endpoint cell binding is absent");
                }
                const auto& target_binding =
                    cell_bindings[
                        *target_binding_index];
                const auto& column_binding =
                    cell_bindings[
                        *column_binding_index];

                if (detail::state_identity_fingerprint(
                        target_binding.state_identity) !=
                        target_hash ||
                    detail::state_identity_fingerprint(
                        column_binding.state_identity) !=
                        column_hash) {
                    throw std::invalid_argument(
                        "received energy endpoint state/chart fingerprint mismatch");
                }

                const auto* schedule_copy =
                    detail::find_local_schedule_copy(
                        schedule,
                        face_global);
                if (schedule_copy == nullptr ||
                    schedule_copy->face !=
                        local_face) {
                    throw std::invalid_argument(
                        "received energy endpoint face is absent from local schedule overlap");
                }

                EnergyIncidentFaceLocalSide3D side{};
                if (side_word ==
                    static_cast<std::uint64_t>(
                        EnergyIncidentFaceLocalSide3D::owner)) {
                    side =
                        EnergyIncidentFaceLocalSide3D::owner;
                    if (schedule_copy
                                ->owner_cell_global !=
                            target_global ||
                        schedule_copy
                                ->neighbour_cell_global !=
                            column_global) {
                        throw std::invalid_argument(
                            "received owner-side energy endpoint orientation mismatch");
                    }
                } else if (
                    side_word ==
                    static_cast<std::uint64_t>(
                        EnergyIncidentFaceLocalSide3D::neighbour)) {
                    side =
                        EnergyIncidentFaceLocalSide3D::neighbour;
                    if (schedule_copy
                                ->neighbour_cell_global !=
                            target_global ||
                        schedule_copy
                                ->owner_cell_global !=
                            column_global) {
                        throw std::invalid_argument(
                            "received neighbour-side energy endpoint orientation mismatch");
                    }
                } else {
                    throw std::invalid_argument(
                        "received energy endpoint side is invalid");
                }

                const std::pair<
                    std::uint64_t,
                    std::uint64_t>
                    target_face_key{
                        target_global.value(),
                        face_global.value()};
                if (std::find(
                        received_target_faces.begin(),
                        received_target_faces.end(),
                        target_face_key) !=
                    received_target_faces.end()) {
                    throw std::invalid_argument(
                        "duplicate received energy endpoint face contribution");
                }
                received_target_faces.push_back(
                    target_face_key);

                std::size_t cursor =
                    value_offset;
                const double target_volume =
                    receive_values[cursor++];
                const double spatial =
                    receive_values[cursor++];

                if (!detail::near_roundoff(
                        target_volume,
                        target_binding.bulk_volume_m3,
                        std::max(
                            std::abs(target_volume),
                            std::abs(
                                target_binding
                                    .bulk_volume_m3))) ||
                    !std::isfinite(spatial)) {
                    throw std::invalid_argument(
                        "received energy endpoint bulk volume/value mismatch");
                }

                std::vector<double> diagonal(
                    column_count);
                std::vector<double> off_diagonal(
                    column_count);
                for (double& value :
                     diagonal) {
                    value =
                        receive_values[cursor++];
                }
                for (double& value :
                     off_diagonal) {
                    value =
                        receive_values[cursor++];
                }
                if (cursor !=
                    value_offset +
                        value_width) {
                    throw std::runtime_error(
                        "energy endpoint payload parser width mismatch");
                }

                for (double value : diagonal) {
                    if (!std::isfinite(value)) {
                        throw std::invalid_argument(
                            "received energy diagonal derivative is non-finite");
                    }
                }
                for (double value : off_diagonal) {
                    if (!std::isfinite(value)) {
                        throw std::invalid_argument(
                            "received energy off-diagonal derivative is non-finite");
                    }
                }

                const auto mapped_row =
                    row_by_local_cell[
                        static_cast<std::size_t>(
                            target_cell.value())];
                if (!mapped_row.has_value()) {
                    throw std::runtime_error(
                        "target owned energy row was not initialized");
                }
                auto& row =
                    rows[*mapped_row];

                row.local_residual
                    .residual_w_per_bulk_m3 +=
                    spatial;
                if (!std::isfinite(
                        row.local_residual
                            .residual_w_per_bulk_m3)) {
                    throw std::range_error(
                        "owned energy residual became non-finite during exchange");
                }

                for (std::size_t column = 0U;
                     column < column_count;
                     ++column) {
                    row.local_residual
                        .local_gradient[column] +=
                        diagonal[column];
                    if (!std::isfinite(
                            row.local_residual
                                .local_gradient[column])) {
                        throw std::range_error(
                            "owned energy diagonal Jacobian became non-finite during exchange");
                    }
                }

                LocalEnergyConservationNeighbourBlock3D
                    face_block{
                        local_face,
                        side,
                        column_binding
                            .state_identity,
                        off_diagonal};

                for (const auto& existing :
                     row.local_residual
                         .neighbour_blocks) {
                    if (existing.face ==
                        local_face) {
                        throw std::invalid_argument(
                            "duplicate per-face energy neighbour block after owner-targeted exchange");
                    }
                }
                row.local_residual
                    .neighbour_blocks
                    .push_back(
                        face_block);

                add_cell_pair_block(
                    row,
                    face_block,
                    column_cell,
                    column_global,
                    face_global);

                ++parsed_payloads;
            }
        }

        for (auto& row : rows) {
            std::size_t expected = 0U;
            const auto inspect =
                [&](auto schedule_rows) {
                    for (const auto& connection :
                         schedule_rows) {
                        if (connection
                                    .owner_cell_global ==
                                row.cell_global ||
                            connection
                                    .neighbour_cell_global ==
                                row.cell_global) {
                            ++expected;
                        }
                    }
                };
            inspect(
                schedule.assembly_rows());
            inspect(
                schedule.ghost_rows());

            if (row.local_residual
                    .neighbour_blocks.size() !=
                expected) {
                throw std::runtime_error(
                    "owned energy row did not receive every incident authoritative face contribution");
            }

            std::sort(
                row.local_residual
                    .neighbour_blocks.begin(),
                row.local_residual
                    .neighbour_blocks.end(),
                [](const auto& left,
                   const auto& right) {
                    return left.face <
                        right.face;
                });
            std::sort(
                row.off_diagonal_cell_pair_blocks.begin(),
                row.off_diagonal_cell_pair_blocks.end(),
                [](const auto& left,
                   const auto& right) {
                    return left.column_cell_global <
                        right.column_cell_global;
                });
        }

        const std::size_t expected_payloads =
            std::accumulate(
                receive_payload_counts.begin(),
                receive_payload_counts.end(),
                std::size_t{0},
                [](std::size_t sum, int count) {
                    return sum +
                        static_cast<std::size_t>(
                            count);
                });
        if (parsed_payloads !=
            expected_payloads) {
            throw std::runtime_error(
                "received energy endpoint payload count mismatch");
        }
    } catch (...) {
        local_error =
            PETSC_ERR_ARG_INCOMP;
    }

    error =
        detail::collective_error(
            comm,
            local_error);
    if (error != PETSC_SUCCESS) {
        return error;
    }

    double local_balance = 0.0;
    double local_scale = 0.0;
    for (const auto& row : rows) {
        const auto binding_index =
            cell_to_binding[
                static_cast<std::size_t>(
                    row.cell.value())];
        if (!binding_index.has_value()) {
            return PETSC_ERR_PLIB;
        }
        const auto& binding =
            cell_bindings[
                *binding_index];
        const double spatial =
            row.local_residual
                .residual_w_per_bulk_m3 -
            binding.owned_accumulation
                ->residual_w_per_bulk_m3;
        const double weighted =
            binding.bulk_volume_m3 *
            spatial;
        local_balance +=
            weighted;
        local_scale +=
            std::abs(weighted);
    }

    double global_balance = 0.0;
    double global_scale = 0.0;
    if (MPI_Allreduce(
            &local_balance,
            &global_balance,
            1,
            MPI_DOUBLE,
            MPI_SUM,
            comm) != MPI_SUCCESS ||
        MPI_Allreduce(
            &local_scale,
            &global_scale,
            1,
            MPI_DOUBLE,
            MPI_SUM,
            comm) != MPI_SUCCESS) {
        return PETSC_ERR_MPI;
    }

    if (!detail::near_roundoff(
            global_balance,
            0.0,
            global_scale)) {
        return PETSC_ERR_PLIB;
    }

    std::size_t received_face_side_count = 0U;
    for (const int count :
         receive_payload_counts) {
        received_face_side_count +=
            static_cast<std::size_t>(
                count);
    }

    try {
        output->emplace(
            partition.local_rank(),
            partition.rank_count(),
            partition.entity_count(
                mpmc::mesh::EntityKind::cell),
            schedule.authoritative_row_count(),
            received_face_side_count,
            std::move(rows),
            global_balance);
    } catch (...) {
        local_error =
            PETSC_ERR_ARG_INCOMP;
    }

    error =
        detail::collective_error(
            comm,
            local_error);
    if (error != PETSC_SUCCESS) {
        output->reset();
        return error;
    }

    return PETSC_SUCCESS;
}

} // namespace mpmc::flow_discretization_petsc

#endif // MPMC_FLOW_DISCRETIZATION_PETSC_DISTRIBUTED_ENERGY_CONSERVATION_HPP
