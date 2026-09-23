#ifndef MPMC_FLOW_DISCRETIZATION_PETSC_DISTRIBUTED_COMPONENT_CONSERVATION_HPP
#define MPMC_FLOW_DISCRETIZATION_PETSC_DISTRIBUTED_COMPONENT_CONSERVATION_HPP

#include <mpmc/discretization_petsc/adapter.hpp>
#include <mpmc/flow_discretization/owned_multi_cell_component_conservation.hpp>

#include <petscsys.h>

#include <algorithm>
#include <array>
#include <bit>
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

static_assert(
    sizeof(double) == sizeof(std::uint64_t),
    "distributed state fingerprint requires 64-bit double");

inline constexpr std::string_view
    distributed_owned_component_conservation_convention =
        "flow_discretization_petsc/owner-targeted-component-conservation/v1";

/// Local cell state required by the distributed owner-targeted exchange.
///
/// Every locally visible owned or ghost cell must have one binding. Only owned
/// cells may carry a backward-Euler accumulation residual; ghost accumulation
/// is deliberately absent because the owner rank is authoritative for that row.
struct DistributedCellStateBinding3D {
    mpmc::mesh::LocalIndex cell{
        mpmc::mesh::LocalIndex::value_type{0}};
    mpmc::mesh::GlobalEntityId cell_global{
        mpmc::mesh::GlobalEntityId::value_type{0}};
    double bulk_volume_m3{};
    mpmc::flow::NaturalVariableStateIdentity3P
        state_identity;
    const mpmc::flow::
        BackwardEulerComponentAccumulationResidual3P*
            owned_accumulation{};
};

class DistributedOwnedMultiCellComponentConservationSnapshot3D {
public:
    static constexpr std::string_view convention =
        distributed_owned_component_conservation_convention;

    DistributedOwnedMultiCellComponentConservationSnapshot3D(
        mpmc::mesh::PartitionRank local_rank,
        std::uint32_t rank_count,
        std::size_t local_cell_count,
        std::size_t local_authoritative_face_count,
        std::size_t received_face_side_count,
        std::vector<std::string> component_ids,
        std::vector<
            mpmc::flow_discretization::
                OwnedCellComponentConservationRow3D>
            rows,
        std::vector<double>
            global_volume_weighted_spatial_component_balance_mol_per_s,
        double global_volume_weighted_spatial_total_balance_mol_per_s)
        : local_rank_(local_rank),
          rank_count_(rank_count),
          local_cell_count_(local_cell_count),
          local_authoritative_face_count_(
              local_authoritative_face_count),
          received_face_side_count_(
              received_face_side_count),
          component_ids_(std::move(component_ids)),
          rows_(std::move(rows)),
          row_by_local_cell_(local_cell_count),
          global_volume_weighted_spatial_component_balance_mol_per_s_(
              std::move(
                  global_volume_weighted_spatial_component_balance_mol_per_s)),
          global_volume_weighted_spatial_total_balance_mol_per_s_(
              global_volume_weighted_spatial_total_balance_mol_per_s) {
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

    [[nodiscard]] std::span<const std::string>
    component_ids() const noexcept {
        return component_ids_;
    }

    [[nodiscard]] std::span<
        const mpmc::flow_discretization::
            OwnedCellComponentConservationRow3D>
    owned_rows() const noexcept {
        return rows_;
    }

    [[nodiscard]]
    const mpmc::flow_discretization::
        OwnedCellComponentConservationRow3D&
    owned_row(
        mpmc::mesh::LocalIndex cell) const {
        const std::size_t local =
            static_cast<std::size_t>(cell.value());
        if (local >= local_cell_count_) {
            throw std::out_of_range(
                "mpmc::flow_discretization_petsc: owned row local cell index out of range");
        }
        const auto mapped =
            row_by_local_cell_[local];
        if (!mapped.has_value()) {
            throw std::invalid_argument(
                "mpmc::flow_discretization_petsc: requested cell is not locally owned");
        }
        return rows_.at(*mapped);
    }

    [[nodiscard]] std::span<const double>
    global_volume_weighted_spatial_component_balance_mol_per_s()
        const noexcept {
        return global_volume_weighted_spatial_component_balance_mol_per_s_;
    }

    [[nodiscard]] double
    global_volume_weighted_spatial_total_balance_mol_per_s()
        const noexcept {
        return global_volume_weighted_spatial_total_balance_mol_per_s_;
    }

private:
    void validate_and_index() {
        if (rank_count_ == 0U ||
            local_rank_.value() >= rank_count_ ||
            component_ids_.size() < 2U ||
            global_volume_weighted_spatial_component_balance_mol_per_s_
                    .size() !=
                component_ids_.size() ||
            !std::isfinite(
                global_volume_weighted_spatial_total_balance_mol_per_s_)) {
            throw std::invalid_argument(
                "mpmc::flow_discretization_petsc: malformed distributed owned component-conservation snapshot");
        }

        for (double value :
             global_volume_weighted_spatial_component_balance_mol_per_s_) {
            if (!std::isfinite(value)) {
                throw std::invalid_argument(
                    "mpmc::flow_discretization_petsc: non-finite global spatial conservation diagnostic");
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
                    component_ids_) {
                throw std::invalid_argument(
                    "mpmc::flow_discretization_petsc: invalid or duplicate distributed owned row");
            }
            row_by_local_cell_[local] = index;
        }
    }

    mpmc::mesh::PartitionRank local_rank_;
    std::uint32_t rank_count_;
    std::size_t local_cell_count_{};
    std::size_t local_authoritative_face_count_{};
    std::size_t received_face_side_count_{};
    std::vector<std::string> component_ids_;
    std::vector<
        mpmc::flow_discretization::
            OwnedCellComponentConservationRow3D>
        rows_;
    std::vector<std::optional<std::size_t>>
        row_by_local_cell_;
    std::vector<double>
        global_volume_weighted_spatial_component_balance_mol_per_s_;
    double global_volume_weighted_spatial_total_balance_mol_per_s_{};
};

namespace detail {

inline constexpr std::size_t endpoint_metadata_word_count = 6U;

[[nodiscard]] inline PetscErrorCode collective_error(
    MPI_Comm comm,
    PetscErrorCode local_error) {
    int local_code =
        static_cast<int>(local_error);
    int global_code =
        static_cast<int>(PETSC_SUCCESS);
    if (MPI_Allreduce(
            &local_code,
            &global_code,
            1,
            MPI_INT,
            MPI_MAX,
            comm) != MPI_SUCCESS) {
        return PETSC_ERR_MPI;
    }
    return static_cast<PetscErrorCode>(
        global_code);
}

[[nodiscard]] inline bool near_roundoff(
    double first,
    double second,
    double scale = 0.0) {
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
        16384.0 *
            std::numeric_limits<double>::epsilon() *
            reference;
}

inline void hash_byte(
    std::uint64_t& hash,
    std::uint8_t value) noexcept {
    hash ^= static_cast<std::uint64_t>(value);
    hash *= UINT64_C(1099511628211);
}

inline void hash_u64(
    std::uint64_t& hash,
    std::uint64_t value) noexcept {
    for (unsigned int shift = 0U;
         shift < 64U;
         shift += 8U) {
        hash_byte(
            hash,
            static_cast<std::uint8_t>(
                (value >> shift) &
                UINT64_C(0xff)));
    }
}

inline void hash_double(
    std::uint64_t& hash,
    double value) noexcept {
    hash_u64(
        hash,
        std::bit_cast<std::uint64_t>(value));
}

inline void hash_string(
    std::uint64_t& hash,
    std::string_view value) noexcept {
    hash_u64(
        hash,
        static_cast<std::uint64_t>(
            value.size()));
    for (char character : value) {
        hash_byte(
            hash,
            static_cast<std::uint8_t>(
                static_cast<unsigned char>(
                    character)));
    }
}

[[nodiscard]] inline std::uint64_t
component_identity_fingerprint(
    std::span<const std::string> component_ids) noexcept {
    std::uint64_t hash =
        UINT64_C(1469598103934665603);
    hash_u64(
        hash,
        static_cast<std::uint64_t>(
            component_ids.size()));
    for (const auto& id : component_ids) {
        hash_string(hash, id);
    }
    return hash;
}

[[nodiscard]] inline std::uint64_t
state_identity_fingerprint(
    const mpmc::flow::NaturalVariableStateIdentity3P&
        identity) noexcept {
    std::uint64_t hash =
        component_identity_fingerprint(
            identity.component_ids);

    for (const auto dependent :
         identity.layout.composition_pivot()
             .dependent_components()) {
        hash_u64(
            hash,
            static_cast<std::uint64_t>(
                dependent));
    }

    hash_double(
        hash,
        identity.reference_pressure_pa);
    hash_double(
        hash,
        identity.temperature_k);
    for (double value : identity.saturation) {
        hash_double(hash, value);
    }
    for (const auto& composition :
         identity.phase_composition) {
        hash_u64(
            hash,
            static_cast<std::uint64_t>(
                composition.size()));
        for (double value : composition) {
            hash_double(hash, value);
        }
    }
    return hash;
}

[[nodiscard]] inline std::size_t
payload_value_count(
    std::size_t component_count,
    std::size_t column_count) {
    const auto maximum =
        std::numeric_limits<std::size_t>::max();
    if (column_count == 0U ||
        component_count >
            maximum / column_count) {
        throw std::length_error(
            "mpmc::flow_discretization_petsc: endpoint payload Jacobian size overflow");
    }
    const std::size_t matrix =
        component_count * column_count;
    if (column_count >
            maximum / 2U ||
        matrix >
            maximum / 2U) {
        throw std::length_error(
            "mpmc::flow_discretization_petsc: endpoint payload value size overflow");
    }
    const std::size_t twice_columns =
        2U * column_count;
    const std::size_t twice_matrix =
        2U * matrix;
    if (component_count >
            maximum - 2U ||
        twice_columns >
            maximum - component_count - 2U ||
        twice_matrix >
            maximum - component_count - 2U -
                twice_columns) {
        throw std::length_error(
            "mpmc::flow_discretization_petsc: endpoint payload value size overflow");
    }
    return component_count +
        2U +
        twice_columns +
        twice_matrix;
}

[[nodiscard]] inline bool checked_int_product(
    std::size_t left,
    std::size_t right,
    int* output) {
    if (output == nullptr ||
        (left != 0U &&
         right >
             static_cast<std::size_t>(
                 std::numeric_limits<int>::max()) /
                 left)) {
        return false;
    }
    const std::size_t value =
        left * right;
    if (value >
        static_cast<std::size_t>(
            std::numeric_limits<int>::max())) {
        return false;
    }
    *output =
        static_cast<int>(value);
    return true;
}

[[nodiscard]] inline bool build_displacements(
    std::span<const int> counts,
    std::vector<int>* displacements,
    int* total) {
    if (displacements == nullptr ||
        total == nullptr) {
        return false;
    }
    displacements->assign(
        counts.size(),
        0);
    int running = 0;
    for (std::size_t rank = 0U;
         rank < counts.size();
         ++rank) {
        const int count = counts[rank];
        if (count < 0 ||
            count >
                std::numeric_limits<int>::max() -
                    running) {
            return false;
        }
        (*displacements)[rank] =
            running;
        running += count;
    }
    *total = running;
    return true;
}

struct EndpointPayload {
    std::uint64_t target_cell_global{};
    std::uint64_t column_cell_global{};
    std::uint64_t face_global{};
    std::uint64_t target_state_fingerprint{};
    std::uint64_t column_state_fingerprint{};
    std::uint64_t local_side{};

    double target_bulk_volume_m3{};
    std::vector<double> component_spatial;
    double total_spatial{};
    std::vector<double> diagonal_component_jacobian;
    std::vector<double> diagonal_total_gradient;
    std::vector<double> off_diagonal_component_jacobian;
    std::vector<double> off_diagonal_total_gradient;
};

[[nodiscard]] inline EndpointPayload make_endpoint_payload(
    const mpmc::discretization_petsc::
        AssemblyReadyInternalConnectionRow3D&
            connection,
    const mpmc::flow_discretization::
        NormalizedComponentFaceContributionLinearization3D&
            face,
    bool target_is_owner) {
    using namespace mpmc::flow_discretization;

    EndpointPayload payload;
    if (target_is_owner) {
        payload.target_cell_global =
            connection.owner_cell_global.value();
        payload.column_cell_global =
            connection.neighbour_cell_global.value();
        payload.local_side =
            static_cast<std::uint64_t>(
                IncidentFaceLocalSide3D::owner);
        payload.target_bulk_volume_m3 =
            face.bulk_volume.owner_bulk_volume_m3;
        payload.component_spatial =
            face.owner_component_contribution_mol_per_bulk_m3_s;
        payload.total_spatial =
            face.owner_total_contribution_mol_per_bulk_m3_s;
        payload.diagonal_component_jacobian =
            face.owner_row_owner_column_jacobian;
        payload.diagonal_total_gradient =
            face.owner_total_row_owner_column_gradient;
        payload.off_diagonal_component_jacobian =
            face.owner_row_neighbour_column_jacobian;
        payload.off_diagonal_total_gradient =
            face.owner_total_row_neighbour_column_gradient;
        payload.target_state_fingerprint =
            state_identity_fingerprint(
                face.owner_state_identity);
        payload.column_state_fingerprint =
            state_identity_fingerprint(
                face.neighbour_state_identity);
    } else {
        payload.target_cell_global =
            connection.neighbour_cell_global.value();
        payload.column_cell_global =
            connection.owner_cell_global.value();
        payload.local_side =
            static_cast<std::uint64_t>(
                IncidentFaceLocalSide3D::neighbour);
        payload.target_bulk_volume_m3 =
            face.bulk_volume.neighbour_bulk_volume_m3;
        payload.component_spatial =
            face.neighbour_component_contribution_mol_per_bulk_m3_s;
        payload.total_spatial =
            face.neighbour_total_contribution_mol_per_bulk_m3_s;
        payload.diagonal_component_jacobian =
            face.neighbour_row_neighbour_column_jacobian;
        payload.diagonal_total_gradient =
            face.neighbour_total_row_neighbour_column_gradient;
        payload.off_diagonal_component_jacobian =
            face.neighbour_row_owner_column_jacobian;
        payload.off_diagonal_total_gradient =
            face.neighbour_total_row_owner_column_gradient;
        payload.target_state_fingerprint =
            state_identity_fingerprint(
                face.neighbour_state_identity);
        payload.column_state_fingerprint =
            state_identity_fingerprint(
                face.owner_state_identity);
    }
    payload.face_global =
        connection.face_global.value();
    return payload;
}

inline void append_payload_values(
    const EndpointPayload& payload,
    std::vector<double>* values) {
    if (values == nullptr) {
        throw std::invalid_argument(
            "mpmc::flow_discretization_petsc: null payload value buffer");
    }
    values->push_back(
        payload.target_bulk_volume_m3);
    values->insert(
        values->end(),
        payload.component_spatial.begin(),
        payload.component_spatial.end());
    values->push_back(
        payload.total_spatial);
    values->insert(
        values->end(),
        payload.diagonal_component_jacobian.begin(),
        payload.diagonal_component_jacobian.end());
    values->insert(
        values->end(),
        payload.diagonal_total_gradient.begin(),
        payload.diagonal_total_gradient.end());
    values->insert(
        values->end(),
        payload.off_diagonal_component_jacobian.begin(),
        payload.off_diagonal_component_jacobian.end());
    values->insert(
        values->end(),
        payload.off_diagonal_total_gradient.begin(),
        payload.off_diagonal_total_gradient.end());
}

[[nodiscard]] inline const mpmc::discretization_petsc::
    AssemblyReadyInternalConnectionRow3D*
find_local_schedule_copy(
    const mpmc::discretization_petsc::
        ParallelOwnedConnectionSchedule3D&
            schedule,
    mpmc::mesh::GlobalEntityId face_global) {
    const mpmc::discretization_petsc::
        AssemblyReadyInternalConnectionRow3D*
            found = nullptr;

    const auto inspect =
        [&](auto rows) {
            for (const auto& row : rows) {
                if (row.face_global !=
                    face_global) {
                    continue;
                }
                if (found != nullptr) {
                    throw std::invalid_argument(
                        "mpmc::flow_discretization_petsc: duplicate local schedule copy for stable face");
                }
                found = &row;
            }
        };

    inspect(schedule.assembly_rows());
    inspect(schedule.ghost_rows());
    return found;
}

inline void validate_cell_pair_block(
    const mpmc::flow_discretization::
        OwnedCellPairComponentJacobianBlock3D&
            block,
    std::size_t component_count) {
    const std::size_t q =
        block.column_state_identity.layout
            .unknown_count();
    if (q == 0U ||
        component_count >
            std::numeric_limits<std::size_t>::max() /
                q ||
        block.component_jacobian.size() !=
            component_count * q ||
        block.total_gradient.size() != q) {
        throw std::invalid_argument(
            "mpmc::flow_discretization_petsc: malformed coalesced cell-pair Jacobian block");
    }

    for (std::size_t column = 0U;
         column < q;
         ++column) {
        double sum = 0.0;
        double scale =
            std::abs(
                block.total_gradient[column]);
        if (!std::isfinite(
                block.total_gradient[column])) {
            throw std::invalid_argument(
                "mpmc::flow_discretization_petsc: non-finite cell-pair total derivative");
        }
        for (std::size_t component = 0U;
             component < component_count;
             ++component) {
            const double value =
                block.component_jacobian[
                    component * q +
                    column];
            if (!std::isfinite(value)) {
                throw std::invalid_argument(
                    "mpmc::flow_discretization_petsc: non-finite cell-pair component derivative");
            }
            sum += value;
            scale += std::abs(value);
        }
        if (!near_roundoff(
                sum,
                block.total_gradient[column],
                scale)) {
            throw std::invalid_argument(
                "mpmc::flow_discretization_petsc: cell-pair component Jacobian does not close to total gradient");
        }
    }
}

} // namespace detail

/// Route every authoritative face-side spatial contribution to the owner rank
/// of its endpoint cell and assemble complete owned component rows.
///
/// No PETSc Mat/Vec value insertion or global scalar numbering occurs here.
inline PetscErrorCode
make_distributed_owned_component_conservation_snapshot_3d(
    MPI_Comm comm,
    const mpmc::discretization_petsc::
        ParallelOwnedConnectionSchedule3D&
            schedule,
    const mpmc::mesh::PartitionSnapshot& partition,
    std::span<const DistributedCellStateBinding3D>
        cell_bindings,
    std::span<
        const mpmc::flow_discretization::
            AuthoritativeNormalizedFaceContributionBinding3D>
        authoritative_face_bindings,
    std::optional<
        DistributedOwnedMultiCellComponentConservationSnapshot3D>*
            output) {
    using namespace detail;
    using namespace mpmc::flow_discretization;

    if (output == nullptr) {
        return PETSC_ERR_ARG_NULL;
    }
    output->reset();

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
        PETSC_SUCCESS;

    std::vector<std::optional<std::size_t>>
        cell_to_binding;
    std::vector<std::optional<std::size_t>>
        face_to_binding;
    std::vector<std::string> component_ids;
    std::size_t component_count = 0U;
    std::size_t column_count = 0U;
    std::size_t value_width = 0U;
    std::vector<
        std::vector<EndpointPayload>>
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
                "rank metadata mismatch");
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
                "binding cardinality mismatch");
        }

        cell_to_binding.assign(
            local_cell_count,
            std::nullopt);
        face_to_binding.assign(
            local_face_count,
            std::nullopt);

        const std::array<
            LocalCellNormalizedFaceContributionBinding3D,
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
                    "invalid or duplicate local cell binding");
            }

            local_component_conservation_detail::
                validate_state_identity(
                    binding.state_identity,
                    "distributed-cell");

            if (index == 0U) {
                component_ids =
                    binding.state_identity
                        .component_ids;
                component_count =
                    component_ids.size();
                column_count =
                    binding.state_identity.layout
                        .unknown_count();
                value_width =
                    payload_value_count(
                        component_count,
                        column_count);
            } else if (
                binding.state_identity.component_ids !=
                    component_ids ||
                binding.state_identity.layout
                        .unknown_count() !=
                    column_count) {
                throw std::invalid_argument(
                    "local cells do not share canonical component identity/count");
            }

            const bool owned =
                partition.is_owned(
                    mpmc::mesh::EntityKind::cell,
                    binding.cell);
            if (owned) {
                if (binding.owned_accumulation ==
                    nullptr) {
                    throw std::invalid_argument(
                        "owned cell is missing backward-Euler accumulation");
                }
                (void)build_local_component_conservation_residual(
                    binding.state_identity,
                    binding.bulk_volume_m3,
                    *binding.owned_accumulation,
                    no_faces);
            } else {
                if (!partition.is_ghost(
                        mpmc::mesh::EntityKind::cell,
                        binding.cell) ||
                    binding.owned_accumulation !=
                        nullptr) {
                    throw std::invalid_argument(
                        "ghost cell must not carry owned accumulation");
                }
            }

            cell_to_binding[local] =
                index;
        }

        for (const auto mapped :
             cell_to_binding) {
            if (!mapped.has_value()) {
                throw std::invalid_argument(
                    "every local cell requires one state binding");
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
                    "invalid or duplicate authoritative face binding");
            }
            local_component_conservation_detail::
                validate_normalized_face(
                    *binding.contribution);
            if (binding.contribution
                    ->component_ids !=
                component_ids) {
                throw std::invalid_argument(
                    "authoritative face component identity mismatch");
            }
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
                    "authoritative schedule row does not match local partition/binding");
            }

            const auto& face_binding =
                authoritative_face_bindings[
                    *face_to_binding[face]];
            if (face_binding.face_global !=
                connection.face_global) {
                throw std::invalid_argument(
                    "authoritative stable face ID mismatch");
            }

            const auto& contribution =
                *face_binding.contribution;
            const auto& owner_binding =
                cell_bindings[
                    *cell_to_binding[owner]];
            const auto& neighbour_binding =
                cell_bindings[
                    *cell_to_binding[neighbour]];

            if (owner_binding.cell_global !=
                    connection.owner_cell_global ||
                neighbour_binding.cell_global !=
                    connection.neighbour_cell_global ||
                state_identity_fingerprint(
                    owner_binding.state_identity) !=
                    state_identity_fingerprint(
                        contribution
                            .owner_state_identity) ||
                state_identity_fingerprint(
                    neighbour_binding.state_identity) !=
                    state_identity_fingerprint(
                        contribution
                            .neighbour_state_identity) ||
                !near_roundoff(
                    owner_binding.bulk_volume_m3,
                    contribution.bulk_volume
                        .owner_bulk_volume_m3,
                    owner_binding.bulk_volume_m3) ||
                !near_roundoff(
                    neighbour_binding.bulk_volume_m3,
                    contribution.bulk_volume
                        .neighbour_bulk_volume_m3,
                    neighbour_binding.bulk_volume_m3)) {
                throw std::invalid_argument(
                    "authoritative face endpoint state/chart/volume mismatch");
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

    PetscErrorCode error =
        collective_error(
            comm,
            local_error);
    if (error != PETSC_SUCCESS) {
        return error;
    }

    std::uint64_t local_component_hash =
        component_identity_fingerprint(
            component_ids);
    std::uint64_t min_component_hash = 0U;
    std::uint64_t max_component_hash = 0U;
    std::uint64_t local_component_count =
        static_cast<std::uint64_t>(
            component_count);
    std::uint64_t min_component_count = 0U;
    std::uint64_t max_component_count = 0U;

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
            &local_component_count,
            &min_component_count,
            1,
            MPI_UINT64_T,
            MPI_MIN,
            comm) != MPI_SUCCESS ||
        MPI_Allreduce(
            &local_component_count,
            &max_component_count,
            1,
            MPI_UINT64_T,
            MPI_MAX,
            comm) != MPI_SUCCESS) {
        return PETSC_ERR_MPI;
    }
    if (min_component_hash !=
            max_component_hash ||
        min_component_count !=
            max_component_count ||
        min_component_count <
            2U) {
        return PETSC_ERR_ARG_INCOMP;
    }
    if (component_count >
        static_cast<std::size_t>(
            std::numeric_limits<int>::max())) {
        return PETSC_ERR_ARG_OUTOFRANGE;
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
        collective_error(
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
        if (!checked_int_product(
                static_cast<std::size_t>(
                    send_payload_counts[rank]),
                endpoint_metadata_word_count,
                &send_metadata_counts[rank]) ||
            !checked_int_product(
                static_cast<std::size_t>(
                    receive_payload_counts[rank]),
                endpoint_metadata_word_count,
                &receive_metadata_counts[rank]) ||
            !checked_int_product(
                static_cast<std::size_t>(
                    send_payload_counts[rank]),
                value_width,
                &send_value_counts[rank]) ||
            !checked_int_product(
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
        (!build_displacements(
             send_metadata_counts,
             &send_metadata_displacements,
             &total_send_metadata) ||
         !build_displacements(
             receive_metadata_counts,
             &receive_metadata_displacements,
             &total_receive_metadata) ||
         !build_displacements(
             send_value_counts,
             &send_value_displacements,
             &total_send_values) ||
         !build_displacements(
             receive_value_counts,
             &receive_value_displacements,
             &total_receive_values))) {
        local_error =
            PETSC_ERR_ARG_OUTOFRANGE;
    }

    error =
        collective_error(
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
                if (payload.component_spatial.size() !=
                        component_count ||
                    payload.diagonal_component_jacobian.size() !=
                        component_count * column_count ||
                    payload.diagonal_total_gradient.size() !=
                        column_count ||
                    payload.off_diagonal_component_jacobian.size() !=
                        component_count * column_count ||
                    payload.off_diagonal_total_gradient.size() !=
                        column_count) {
                    throw std::invalid_argument(
                        "endpoint payload shape mismatch");
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
        collective_error(
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

    std::vector<
        mpmc::flow_discretization::
            OwnedCellComponentConservationRow3D>
        rows;
    std::vector<std::optional<std::size_t>>
        row_by_local_cell(
            partition.entity_count(
                mpmc::mesh::EntityKind::cell));

    try {
        const std::array<
            LocalCellNormalizedFaceContributionBinding3D,
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
                build_local_component_conservation_residual(
                    binding.state_identity,
                    binding.bulk_volume_m3,
                    *binding.owned_accumulation,
                    no_faces);
            const std::size_t local_index =
                static_cast<std::size_t>(
                    binding.cell.value());
            row_by_local_cell[local_index] =
                rows.size();
            rows.push_back(
                OwnedCellComponentConservationRow3D{
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
                        "received endpoint stable identity is absent from local overlap");
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
                        "received endpoint payload was routed to a non-owner rank");
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
                        "received endpoint cell binding is absent");
                }
                const auto& target_binding =
                    cell_bindings[
                        *target_binding_index];
                const auto& column_binding =
                    cell_bindings[
                        *column_binding_index];

                if (state_identity_fingerprint(
                        target_binding.state_identity) !=
                        target_hash ||
                    state_identity_fingerprint(
                        column_binding.state_identity) !=
                        column_hash) {
                    throw std::invalid_argument(
                        "received endpoint state/chart fingerprint mismatch");
                }

                const auto* schedule_copy =
                    find_local_schedule_copy(
                        schedule,
                        face_global);
                if (schedule_copy == nullptr ||
                    schedule_copy->face !=
                        local_face ||
                    schedule_copy->face_global !=
                        face_global) {
                    throw std::invalid_argument(
                        "received endpoint face is absent from local schedule overlap");
                }

                IncidentFaceLocalSide3D side{};
                if (side_word ==
                    static_cast<std::uint64_t>(
                        IncidentFaceLocalSide3D::owner)) {
                    side =
                        IncidentFaceLocalSide3D::owner;
                    if (schedule_copy
                                ->owner_cell_global !=
                            target_global ||
                        schedule_copy
                                ->neighbour_cell_global !=
                            column_global) {
                        throw std::invalid_argument(
                            "received owner-side endpoint orientation mismatch");
                    }
                } else if (
                    side_word ==
                    static_cast<std::uint64_t>(
                        IncidentFaceLocalSide3D::neighbour)) {
                    side =
                        IncidentFaceLocalSide3D::neighbour;
                    if (schedule_copy
                                ->neighbour_cell_global !=
                            target_global ||
                        schedule_copy
                                ->owner_cell_global !=
                            column_global) {
                        throw std::invalid_argument(
                            "received neighbour-side endpoint orientation mismatch");
                    }
                } else {
                    throw std::invalid_argument(
                        "received endpoint side is invalid");
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
                        "duplicate received endpoint face contribution");
                }
                received_target_faces.push_back(
                    target_face_key);

                std::size_t cursor =
                    value_offset;
                const double target_volume =
                    receive_values[cursor++];
                if (!near_roundoff(
                        target_volume,
                        target_binding.bulk_volume_m3,
                        std::max(
                            std::abs(target_volume),
                            std::abs(
                                target_binding
                                    .bulk_volume_m3)))) {
                    throw std::invalid_argument(
                        "received endpoint bulk volume mismatch");
                }

                std::vector<double> component_spatial(
                    component_count);
                for (std::size_t component = 0U;
                     component < component_count;
                     ++component) {
                    component_spatial[component] =
                        receive_values[cursor++];
                }
                const double total_spatial =
                    receive_values[cursor++];

                const std::size_t matrix_count =
                    component_count *
                    column_count;
                std::vector<double> diagonal_component(
                    matrix_count);
                for (double& value :
                     diagonal_component) {
                    value =
                        receive_values[cursor++];
                }
                std::vector<double> diagonal_total(
                    column_count);
                for (double& value :
                     diagonal_total) {
                    value =
                        receive_values[cursor++];
                }
                std::vector<double> off_component(
                    matrix_count);
                for (double& value :
                     off_component) {
                    value =
                        receive_values[cursor++];
                }
                std::vector<double> off_total(
                    column_count);
                for (double& value :
                     off_total) {
                    value =
                        receive_values[cursor++];
                }

                if (cursor !=
                    value_offset +
                        value_width) {
                    throw std::runtime_error(
                        "endpoint payload parser width mismatch");
                }

                const auto mapped_row =
                    row_by_local_cell[
                        static_cast<std::size_t>(
                            target_cell.value())];
                if (!mapped_row.has_value()) {
                    throw std::runtime_error(
                        "target owned row was not initialized");
                }
                auto& row =
                    rows[*mapped_row];

                for (std::size_t component = 0U;
                     component < component_count;
                     ++component) {
                    row.local_residual
                        .component_residual_mol_per_bulk_m3_s[
                            component] +=
                        component_spatial[component];
                    for (std::size_t column = 0U;
                         column < column_count;
                         ++column) {
                        const std::size_t index =
                            component *
                                column_count +
                            column;
                        row.local_residual
                            .local_component_jacobian[
                                index] +=
                            diagonal_component[
                                index];
                    }
                }
                row.local_residual
                    .total_residual_mol_per_bulk_m3_s +=
                    total_spatial;
                for (std::size_t column = 0U;
                     column < column_count;
                     ++column) {
                    row.local_residual
                        .local_total_gradient[column] +=
                        diagonal_total[column];
                }

                LocalComponentConservationNeighbourBlock3D
                    face_block{
                        local_face,
                        side,
                        column_binding
                            .state_identity,
                        off_component,
                        off_total};

                for (const auto& existing :
                     row.local_residual
                         .neighbour_blocks) {
                    if (existing.face ==
                        local_face) {
                        throw std::invalid_argument(
                            "duplicate per-face neighbour block after owner-targeted exchange");
                    }
                }
                row.local_residual
                    .neighbour_blocks
                    .push_back(
                        face_block);

                AuthoritativeInternalConnectionIdentity3D
                    connection_identity{
                        local_face,
                        face_global,
                        side ==
                                IncidentFaceLocalSide3D::
                                    owner
                            ? target_cell
                            : column_cell,
                        side ==
                                IncidentFaceLocalSide3D::
                                    owner
                            ? target_global
                            : column_global,
                        side ==
                                IncidentFaceLocalSide3D::
                                    owner
                            ? column_cell
                            : target_cell,
                        side ==
                                IncidentFaceLocalSide3D::
                                    owner
                            ? column_global
                            : target_global};

                owned_multi_cell_conservation_detail::
                    add_cell_pair_block(
                        row,
                        face_block,
                        connection_identity,
                        column_cell,
                        column_global);

                ++parsed_payloads;
            }
        }

        for (auto& row : rows) {
            const auto expected =
                [&]() {
                    std::size_t count = 0U;
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
                                    ++count;
                                }
                            }
                        };
                    inspect(
                        schedule.assembly_rows());
                    inspect(
                        schedule.ghost_rows());
                    return count;
                }();

            if (row.local_residual
                    .neighbour_blocks.size() !=
                expected) {
                throw std::runtime_error(
                    "owned row did not receive every incident authoritative face contribution");
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

            local_component_conservation_detail::
                validate_final_closure(
                    row.local_residual);
            for (const auto& block :
                 row.off_diagonal_cell_pair_blocks) {
                validate_cell_pair_block(
                    block,
                    component_count);
            }
        }

        if (parsed_payloads !=
            std::accumulate(
                receive_payload_counts.begin(),
                receive_payload_counts.end(),
                std::size_t{0},
                [](std::size_t sum, int count) {
                    return sum +
                        static_cast<std::size_t>(
                            count);
                })) {
            throw std::runtime_error(
                "received endpoint payload count mismatch");
        }
    } catch (...) {
        local_error =
            PETSC_ERR_ARG_INCOMP;
    }

    error =
        collective_error(
            comm,
            local_error);
    if (error != PETSC_SUCCESS) {
        return error;
    }

    std::vector<double> local_balance(
        component_count,
        0.0);
    std::vector<double> local_scale(
        component_count,
        0.0);
    double local_total_balance = 0.0;
    double local_total_scale = 0.0;

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

        for (std::size_t component = 0U;
             component < component_count;
             ++component) {
            const double spatial =
                row.local_residual
                    .residual(component) -
                binding.owned_accumulation
                    ->residual(component);
            const double weighted =
                binding.bulk_volume_m3 *
                spatial;
            local_balance[component] +=
                weighted;
            local_scale[component] +=
                std::abs(weighted);
        }

        const double weighted_total =
            binding.bulk_volume_m3 *
            (row.local_residual
                 .total_residual_mol_per_bulk_m3_s -
             binding.owned_accumulation
                 ->total_residual_mol_per_bulk_m3_s);
        local_total_balance +=
            weighted_total;
        local_total_scale +=
            std::abs(weighted_total);
    }

    std::vector<double> global_balance(
        component_count,
        0.0);
    std::vector<double> global_scale(
        component_count,
        0.0);
    double global_total_balance = 0.0;
    double global_total_scale = 0.0;

    if (MPI_Allreduce(
            local_balance.data(),
            global_balance.data(),
            static_cast<int>(
                component_count),
            MPI_DOUBLE,
            MPI_SUM,
            comm) != MPI_SUCCESS ||
        MPI_Allreduce(
            local_scale.data(),
            global_scale.data(),
            static_cast<int>(
                component_count),
            MPI_DOUBLE,
            MPI_SUM,
            comm) != MPI_SUCCESS ||
        MPI_Allreduce(
            &local_total_balance,
            &global_total_balance,
            1,
            MPI_DOUBLE,
            MPI_SUM,
            comm) != MPI_SUCCESS ||
        MPI_Allreduce(
            &local_total_scale,
            &global_total_scale,
            1,
            MPI_DOUBLE,
            MPI_SUM,
            comm) != MPI_SUCCESS) {
        return PETSC_ERR_MPI;
    }

    for (std::size_t component = 0U;
         component < component_count;
         ++component) {
        if (!near_roundoff(
                global_balance[component],
                0.0,
                global_scale[component])) {
            return PETSC_ERR_PLIB;
        }
    }
    if (!near_roundoff(
            global_total_balance,
            0.0,
            global_total_scale)) {
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
            std::move(component_ids),
            std::move(rows),
            std::move(global_balance),
            global_total_balance);
    } catch (...) {
        local_error =
            PETSC_ERR_ARG_INCOMP;
    }

    error =
        collective_error(
            comm,
            local_error);
    if (error != PETSC_SUCCESS) {
        output->reset();
        return error;
    }

    return PETSC_SUCCESS;
}

} // namespace mpmc::flow_discretization_petsc

#endif // MPMC_FLOW_DISCRETIZATION_PETSC_DISTRIBUTED_COMPONENT_CONSERVATION_HPP
