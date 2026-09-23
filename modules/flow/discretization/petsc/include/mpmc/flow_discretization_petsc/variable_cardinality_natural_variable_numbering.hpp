#ifndef MPMC_FLOW_DISCRETIZATION_PETSC_VARIABLE_CARDINALITY_NATURAL_VARIABLE_NUMBERING_HPP
#define MPMC_FLOW_DISCRETIZATION_PETSC_VARIABLE_CARDINALITY_NATURAL_VARIABLE_NUMBERING_HPP

#include <mpmc/flow/natural_variable_cell_state.hpp>
#include <mpmc/mesh/partition_snapshot.hpp>

#include <petscsection.h>
#include <petscvec.h>

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <optional>
#include <span>
#include <stdexcept>
#include <string_view>
#include <utility>
#include <vector>

namespace mpmc::flow_discretization_petsc {

inline constexpr std::string_view
    variable_cardinality_natural_variable_numbering_convention =
        "flow_discretization_petsc/variable-cardinality-natural-variable-numbering/v1";

struct VariableCardinalityNaturalVariableCellDof3D {
    mpmc::mesh::LocalIndex cell{
        mpmc::mesh::LocalIndex::value_type{0}};
    mpmc::mesh::GlobalEntityId cell_global{
        mpmc::mesh::GlobalEntityId::value_type{0}};
    mpmc::mesh::PartitionRank owner_rank{
        mpmc::mesh::PartitionRank::value_type{0}};
    std::size_t phase_count{};
    std::size_t local_scalar_offset{};
    std::size_t scalar_count{};
    PetscInt petsc_global_scalar_start{};
};

class VariableCardinalityNaturalVariableNumbering3D {
public:
    VariableCardinalityNaturalVariableNumbering3D(
        mpmc::mesh::PartitionRank local_rank,
        std::uint32_t rank_count,
        std::size_t component_count,
        PetscInt petsc_owned_scalar_start,
        PetscInt petsc_owned_scalar_end,
        PetscInt petsc_global_scalar_count,
        std::size_t local_packed_scalar_count,
        std::vector<VariableCardinalityNaturalVariableCellDof3D>
            cells,
        std::vector<PetscInt>
            local_packed_scalar_global_indices)
        : local_rank_(local_rank),
          rank_count_(rank_count),
          component_count_(component_count),
          petsc_owned_scalar_start_(
              petsc_owned_scalar_start),
          petsc_owned_scalar_end_(
              petsc_owned_scalar_end),
          petsc_global_scalar_count_(
              petsc_global_scalar_count),
          local_packed_scalar_count_(
              local_packed_scalar_count),
          cells_(std::move(cells)),
          local_packed_scalar_global_indices_(
              std::move(
                  local_packed_scalar_global_indices)) {
        validate();
    }

    VariableCardinalityNaturalVariableNumbering3D(
        const VariableCardinalityNaturalVariableNumbering3D&) =
        default;
    VariableCardinalityNaturalVariableNumbering3D(
        VariableCardinalityNaturalVariableNumbering3D&&) noexcept =
        default;
    VariableCardinalityNaturalVariableNumbering3D& operator=(
        const VariableCardinalityNaturalVariableNumbering3D&) =
        delete;
    VariableCardinalityNaturalVariableNumbering3D& operator=(
        VariableCardinalityNaturalVariableNumbering3D&&) =
        delete;
    ~VariableCardinalityNaturalVariableNumbering3D() =
        default;

    [[nodiscard]] mpmc::mesh::PartitionRank
    local_rank() const noexcept {
        return local_rank_;
    }

    [[nodiscard]] std::uint32_t
    rank_count() const noexcept {
        return rank_count_;
    }

    [[nodiscard]] std::size_t
    component_count() const noexcept {
        return component_count_;
    }

    [[nodiscard]] std::size_t
    local_cell_count() const noexcept {
        return cells_.size();
    }

    [[nodiscard]] std::size_t
    local_packed_scalar_count() const noexcept {
        return local_packed_scalar_count_;
    }

    [[nodiscard]] PetscInt
    petsc_owned_scalar_start() const noexcept {
        return petsc_owned_scalar_start_;
    }

    [[nodiscard]] PetscInt
    petsc_owned_scalar_end() const noexcept {
        return petsc_owned_scalar_end_;
    }

    [[nodiscard]] PetscInt
    petsc_local_owned_scalar_count() const noexcept {
        return petsc_owned_scalar_end_ -
            petsc_owned_scalar_start_;
    }

    [[nodiscard]] PetscInt
    petsc_global_scalar_count() const noexcept {
        return petsc_global_scalar_count_;
    }

    [[nodiscard]] std::span<
        const VariableCardinalityNaturalVariableCellDof3D>
    cells() const noexcept {
        return cells_;
    }

    [[nodiscard]] const
    VariableCardinalityNaturalVariableCellDof3D&
    cell(
        mpmc::mesh::LocalIndex local_cell) const {
        const std::size_t local =
            static_cast<std::size_t>(
                local_cell.value());
        return cells_.at(local);
    }

    [[nodiscard]] bool is_owned_cell(
        mpmc::mesh::LocalIndex local_cell) const {
        return cell(local_cell).owner_rank ==
            local_rank_;
    }

    [[nodiscard]] bool is_ghost_cell(
        mpmc::mesh::LocalIndex local_cell) const {
        return !is_owned_cell(local_cell);
    }

    [[nodiscard]] PetscInt
    petsc_global_scalar(
        mpmc::mesh::LocalIndex local_cell,
        std::size_t natural_variable_slot) const {
        const auto& record =
            cell(local_cell);
        if (natural_variable_slot >=
            record.scalar_count) {
            throw std::out_of_range(
                "mpmc::flow_discretization_petsc: natural-variable slot out of range");
        }
        if (natural_variable_slot >
            static_cast<std::size_t>(
                std::numeric_limits<PetscInt>::max()) ||
            record.petsc_global_scalar_start >
                std::numeric_limits<PetscInt>::max() -
                    static_cast<PetscInt>(
                        natural_variable_slot)) {
            throw std::overflow_error(
                "mpmc::flow_discretization_petsc: PETSc scalar index overflow");
        }
        return record.petsc_global_scalar_start +
            static_cast<PetscInt>(
                natural_variable_slot);
    }

    [[nodiscard]] PetscInt
    local_packed_scalar_global_index(
        std::size_t local_scalar) const {
        return local_packed_scalar_global_indices_
            .at(local_scalar);
    }

private:
    void validate() const {
        if (rank_count_ == 0U ||
            local_rank_.value() >=
                rank_count_ ||
            component_count_ < 2U ||
            petsc_owned_scalar_start_ < 0 ||
            petsc_owned_scalar_end_ <
                petsc_owned_scalar_start_ ||
            petsc_global_scalar_count_ <
                petsc_owned_scalar_end_ ||
            local_packed_scalar_global_indices_
                    .size() !=
                local_packed_scalar_count_) {
            throw std::invalid_argument(
                "mpmc::flow_discretization_petsc: malformed variable-cardinality numbering metadata");
        }

        std::size_t expected_local_offset = 0U;
        PetscInt counted_owned = 0;
        for (std::size_t local = 0U;
             local < cells_.size();
             ++local) {
            const auto& record =
                cells_[local];
            if (static_cast<std::size_t>(
                    record.cell.value()) !=
                    local ||
                record.owner_rank.value() >=
                    rank_count_ ||
                record.phase_count == 0U ||
                record.phase_count >
                    mpmc::flow::
                        fixed_three_phase_count ||
                component_count_ >
                    (std::numeric_limits<
                         std::size_t>::max() -
                     1U) /
                        record.phase_count ||
                record.scalar_count !=
                    record.phase_count *
                        component_count_ +
                        1U ||
                record.local_scalar_offset !=
                    expected_local_offset ||
                record.petsc_global_scalar_start <
                    0 ||
                record.petsc_global_scalar_start >
                    petsc_global_scalar_count_ ||
                record.scalar_count >
                    static_cast<std::size_t>(
                        std::numeric_limits<PetscInt>::
                            max()) ||
                record.scalar_count >
                    static_cast<std::size_t>(
                        petsc_global_scalar_count_) ||
                record.petsc_global_scalar_start >
                    petsc_global_scalar_count_ -
                        static_cast<PetscInt>(
                            record.scalar_count)) {
                throw std::invalid_argument(
                    "mpmc::flow_discretization_petsc: malformed variable-cardinality cell record");
            }

            if (record.scalar_count >
                std::numeric_limits<std::size_t>::
                    max() -
                    expected_local_offset) {
                throw std::length_error(
                    "mpmc::flow_discretization_petsc: local packed scalar count overflow");
            }

            for (std::size_t slot = 0U;
                 slot < record.scalar_count;
                 ++slot) {
                const PetscInt expected =
                    record.petsc_global_scalar_start +
                    static_cast<PetscInt>(slot);
                if (local_packed_scalar_global_indices_[
                        expected_local_offset + slot] !=
                    expected) {
                    throw std::invalid_argument(
                        "mpmc::flow_discretization_petsc: local packed scalar map disagrees with cell global start");
                }
            }

            if (record.owner_rank ==
                local_rank_) {
                if (record.petsc_global_scalar_start <
                        petsc_owned_scalar_start_ ||
                    record.petsc_global_scalar_start +
                            static_cast<PetscInt>(
                                record.scalar_count) >
                        petsc_owned_scalar_end_) {
                    throw std::invalid_argument(
                        "mpmc::flow_discretization_petsc: owned variable-cardinality cell escapes PETSc ownership range");
                }
                const PetscInt record_width =
                    static_cast<PetscInt>(
                        record.scalar_count);
                if (counted_owned >
                    std::numeric_limits<PetscInt>::
                            max() -
                        record_width) {
                    throw std::length_error(
                        "mpmc::flow_discretization_petsc: owned scalar count overflow");
                }
                counted_owned +=
                    record_width;
            } else {
                const PetscInt cell_end =
                    record.petsc_global_scalar_start +
                    static_cast<PetscInt>(
                        record.scalar_count);
                const bool overlaps_local_range =
                    record.petsc_global_scalar_start <
                        petsc_owned_scalar_end_ &&
                    cell_end >
                        petsc_owned_scalar_start_;
                if (overlaps_local_range) {
                    throw std::invalid_argument(
                        "mpmc::flow_discretization_petsc: ghost scalar range overlaps local PETSc ownership");
                }
            }

            expected_local_offset +=
                record.scalar_count;
        }

        if (expected_local_offset !=
                local_packed_scalar_count_ ||
            counted_owned !=
                petsc_local_owned_scalar_count()) {
            throw std::invalid_argument(
                "mpmc::flow_discretization_petsc: variable-cardinality scalar cardinality mismatch");
        }
    }

    mpmc::mesh::PartitionRank local_rank_;
    std::uint32_t rank_count_;
    std::size_t component_count_;
    PetscInt petsc_owned_scalar_start_;
    PetscInt petsc_owned_scalar_end_;
    PetscInt petsc_global_scalar_count_;
    std::size_t local_packed_scalar_count_;
    std::vector<
        VariableCardinalityNaturalVariableCellDof3D>
        cells_;
    std::vector<PetscInt>
        local_packed_scalar_global_indices_;
};

namespace variable_cardinality_numbering_detail {

[[nodiscard]] inline PetscErrorCode
collective_error(
    MPI_Comm comm,
    PetscErrorCode local_error) {
    int local =
        static_cast<int>(local_error);
    int global = 0;
    if (MPI_Allreduce(
            &local,
            &global,
            1,
            MPI_INT,
            MPI_MAX,
            comm) != MPI_SUCCESS) {
        return PETSC_ERR_MPI;
    }
    return static_cast<PetscErrorCode>(
        global);
}

[[nodiscard]] inline bool
checked_petsc_int(
    std::size_t value,
    PetscInt* output) {
    if (output == nullptr ||
        value >
            static_cast<std::size_t>(
                std::numeric_limits<PetscInt>::
                    max())) {
        return false;
    }
    *output =
        static_cast<PetscInt>(value);
    return true;
}

struct GatheredOwnedCell {
    std::uint64_t global_id{};
    std::uint64_t global_scalar_start{};
    std::uint64_t scalar_count{};
    std::uint32_t owner_rank{};
};

} // namespace variable_cardinality_numbering_detail

inline PetscErrorCode
make_variable_cardinality_natural_variable_numbering_3d(
    MPI_Comm comm,
    const mpmc::mesh::PartitionSnapshot&
        partition,
    std::size_t component_count,
    std::span<const std::size_t>
        local_phase_counts,
    std::optional<
        VariableCardinalityNaturalVariableNumbering3D>*
            output) {
    using namespace
        variable_cardinality_numbering_detail;

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
    if (partition.local_rank().value() !=
            static_cast<std::uint32_t>(
                mpi_rank) ||
        partition.rank_count() !=
            static_cast<std::uint32_t>(
                mpi_size) ||
        component_count < 2U ||
        local_phase_counts.size() !=
            partition.entity_count(
                mpmc::mesh::EntityKind::cell)) {
        local_error =
            PETSC_ERR_ARG_INCOMP;
    }

    std::vector<std::size_t>
        local_scalar_offsets;
    std::vector<std::size_t>
        local_scalar_counts;
    std::size_t local_packed_scalar_count =
        0U;
    PetscInt local_owned_scalar_count = 0;

    if (local_error == PETSC_SUCCESS) {
        try {
            local_scalar_offsets.resize(
                local_phase_counts.size());
            local_scalar_counts.resize(
                local_phase_counts.size());

            for (std::size_t local = 0U;
                 local < local_phase_counts.size();
                 ++local) {
                const std::size_t phase_count =
                    local_phase_counts[local];
                if (phase_count == 0U ||
                    phase_count >
                        mpmc::flow::
                            fixed_three_phase_count ||
                    component_count >
                        (std::numeric_limits<
                             std::size_t>::max() -
                         1U) /
                            phase_count) {
                    throw std::invalid_argument(
                        "invalid per-cell phase count");
                }

                const std::size_t q =
                    phase_count *
                        component_count +
                    1U;
                PetscInt q_petsc = 0;
                if (!checked_petsc_int(
                        q,
                        &q_petsc) ||
                    q >
                        std::numeric_limits<
                            std::size_t>::max() -
                            local_packed_scalar_count) {
                    throw std::length_error(
                        "variable-cardinality cell width overflow");
                }

                local_scalar_offsets[local] =
                    local_packed_scalar_count;
                local_scalar_counts[local] =
                    q;
                local_packed_scalar_count +=
                    q;

                const auto cell =
                    mpmc::mesh::LocalIndex{
                        static_cast<
                            mpmc::mesh::LocalIndex::value_type>(
                                local)};
                if (partition.is_owned(
                        mpmc::mesh::EntityKind::cell,
                        cell)) {
                    if (local_owned_scalar_count >
                        std::numeric_limits<PetscInt>::
                                max() -
                            q_petsc) {
                        throw std::length_error(
                            "owned scalar count overflow");
                    }
                    local_owned_scalar_count +=
                        q_petsc;
                }
            }
        } catch (...) {
            local_error =
                PETSC_ERR_ARG_INCOMP;
        }
    }

    PetscErrorCode error =
        collective_error(
            comm,
            local_error);
    if (error != PETSC_SUCCESS) {
        return error;
    }

    const std::uint64_t local_owned_scalar_count_u64 =
        static_cast<std::uint64_t>(
            local_owned_scalar_count);
    std::uint64_t global_scalar_count_u64 = 0U;
    if (MPI_Allreduce(
            &local_owned_scalar_count_u64,
            &global_scalar_count_u64,
            1,
            MPI_UINT64_T,
            MPI_SUM,
            comm) != MPI_SUCCESS) {
        return PETSC_ERR_MPI;
    }
    if (global_scalar_count_u64 >
        static_cast<std::uint64_t>(
            std::numeric_limits<PetscInt>::max())) {
        return PETSC_ERR_ARG_OUTOFRANGE;
    }
    const PetscInt global_scalar_count =
        static_cast<PetscInt>(
            global_scalar_count_u64);

    PetscLayout layout = nullptr;
    error =
        PetscLayoutCreate(
            comm,
            &layout);
    if (error == PETSC_SUCCESS) {
        error =
            PetscLayoutSetLocalSize(
                layout,
                local_owned_scalar_count);
    }
    if (error == PETSC_SUCCESS) {
        error =
            PetscLayoutSetSize(
                layout,
                global_scalar_count);
    }
    if (error == PETSC_SUCCESS) {
        error =
            PetscLayoutSetUp(
                layout);
    }

    PetscInt owned_start = -1;
    PetscInt owned_end = -1;
    if (error == PETSC_SUCCESS) {
        error =
            PetscLayoutGetRange(
                layout,
                &owned_start,
                &owned_end);
    }
    if (error != PETSC_SUCCESS ||
        owned_start < 0 ||
        owned_end < owned_start ||
        owned_end - owned_start !=
            local_owned_scalar_count) {
        PetscLayoutDestroy(
            &layout);
        return error != PETSC_SUCCESS
            ? error
            : PETSC_ERR_PLIB;
    }

    std::vector<std::size_t>
        owned_local_cells;
    owned_local_cells.reserve(
        partition.owned_count(
            mpmc::mesh::EntityKind::cell));
    for (std::size_t local = 0U;
         local < local_phase_counts.size();
         ++local) {
        const auto cell =
            mpmc::mesh::LocalIndex{
                static_cast<
                    mpmc::mesh::LocalIndex::value_type>(
                        local)};
        if (partition.is_owned(
                mpmc::mesh::EntityKind::cell,
                cell)) {
            owned_local_cells.push_back(
                local);
        }
    }
    std::sort(
        owned_local_cells.begin(),
        owned_local_cells.end(),
        [&](std::size_t left,
            std::size_t right) {
            return partition.global_id(
                       mpmc::mesh::EntityKind::cell,
                       mpmc::mesh::LocalIndex{
                           static_cast<
                               mpmc::mesh::LocalIndex::value_type>(
                                   left)}) <
                partition.global_id(
                    mpmc::mesh::EntityKind::cell,
                    mpmc::mesh::LocalIndex{
                        static_cast<
                            mpmc::mesh::LocalIndex::value_type>(
                                right)});
        });

    std::vector<std::uint64_t>
        local_owned_global_ids;
    std::vector<std::uint64_t>
        local_owned_global_starts;
    std::vector<std::uint64_t>
        local_owned_scalar_widths;
    local_owned_global_ids.reserve(
        owned_local_cells.size());
    local_owned_global_starts.reserve(
        owned_local_cells.size());
    local_owned_scalar_widths.reserve(
        owned_local_cells.size());

    PetscInt running = owned_start;
    for (const std::size_t local :
         owned_local_cells) {
        const auto cell =
            mpmc::mesh::LocalIndex{
                static_cast<
                    mpmc::mesh::LocalIndex::value_type>(
                        local)};
        const std::size_t q =
            local_scalar_counts[local];
        local_owned_global_ids.push_back(
            partition.global_id(
                mpmc::mesh::EntityKind::cell,
                cell)
                .value());
        local_owned_global_starts.push_back(
            static_cast<std::uint64_t>(
                running));
        local_owned_scalar_widths.push_back(
            static_cast<std::uint64_t>(q));
        running +=
            static_cast<PetscInt>(q);
    }
    if (running != owned_end) {
        PetscLayoutDestroy(
            &layout);
        return PETSC_ERR_PLIB;
    }

    if (owned_local_cells.size() >
        static_cast<std::size_t>(
            std::numeric_limits<int>::max())) {
        PetscLayoutDestroy(
            &layout);
        return PETSC_ERR_ARG_OUTOFRANGE;
    }
    const int local_owned_cell_count =
        static_cast<int>(
            owned_local_cells.size());

    std::vector<int> owned_counts(
        static_cast<std::size_t>(
            mpi_size),
        0);
    if (MPI_Allgather(
            &local_owned_cell_count,
            1,
            MPI_INT,
            owned_counts.data(),
            1,
            MPI_INT,
            comm) != MPI_SUCCESS) {
        PetscLayoutDestroy(
            &layout);
        return PETSC_ERR_MPI;
    }

    std::vector<int> displacements(
        static_cast<std::size_t>(
            mpi_size),
        0);
    int total_owned_cells = 0;
    for (int rank = 0;
         rank < mpi_size;
         ++rank) {
        if (owned_counts[
                static_cast<std::size_t>(
                    rank)] < 0 ||
            owned_counts[
                static_cast<std::size_t>(
                    rank)] >
                std::numeric_limits<int>::max() -
                    total_owned_cells) {
            PetscLayoutDestroy(
                &layout);
            return PETSC_ERR_ARG_OUTOFRANGE;
        }
        displacements[
            static_cast<std::size_t>(
                rank)] =
            total_owned_cells;
        total_owned_cells +=
            owned_counts[
                static_cast<std::size_t>(
                    rank)];
    }

    std::vector<std::uint64_t>
        gathered_ids(
            static_cast<std::size_t>(
                total_owned_cells));
    std::vector<std::uint64_t>
        gathered_starts(
            static_cast<std::size_t>(
                total_owned_cells));
    std::vector<std::uint64_t>
        gathered_widths(
            static_cast<std::size_t>(
                total_owned_cells));

    const auto allgather_u64 =
        [&](const std::vector<std::uint64_t>&
                send,
            std::vector<std::uint64_t>&
                receive) {
            return MPI_Allgatherv(
                send.empty()
                    ? nullptr
                    : send.data(),
                local_owned_cell_count,
                MPI_UINT64_T,
                receive.empty()
                    ? nullptr
                    : receive.data(),
                owned_counts.data(),
                displacements.data(),
                MPI_UINT64_T,
                comm);
        };

    if (allgather_u64(
            local_owned_global_ids,
            gathered_ids) != MPI_SUCCESS ||
        allgather_u64(
            local_owned_global_starts,
            gathered_starts) != MPI_SUCCESS ||
        allgather_u64(
            local_owned_scalar_widths,
            gathered_widths) != MPI_SUCCESS) {
        PetscLayoutDestroy(
            &layout);
        return PETSC_ERR_MPI;
    }

    std::vector<GatheredOwnedCell>
        gathered;
    gathered.reserve(
        gathered_ids.size());
    for (int rank = 0;
         rank < mpi_size;
         ++rank) {
        const int begin =
            displacements[
                static_cast<std::size_t>(
                    rank)];
        const int count =
            owned_counts[
                static_cast<std::size_t>(
                    rank)];
        for (int offset = 0;
             offset < count;
             ++offset) {
            const std::size_t index =
                static_cast<std::size_t>(
                    begin + offset);
            gathered.push_back(
                {
                    gathered_ids[index],
                    gathered_starts[index],
                    gathered_widths[index],
                    static_cast<std::uint32_t>(
                        rank)});
        }
    }
    std::sort(
        gathered.begin(),
        gathered.end(),
        [](const auto& left,
           const auto& right) {
            return left.global_id <
                right.global_id;
        });

    for (std::size_t index = 0U;
         index < gathered.size();
         ++index) {
        if (gathered[index]
                .global_scalar_start >
                static_cast<std::uint64_t>(
                    global_scalar_count) ||
            gathered[index]
                .scalar_count == 0U ||
            gathered[index]
                .scalar_count >
                static_cast<std::uint64_t>(
                    global_scalar_count) ||
            gathered[index]
                .global_scalar_start >
                static_cast<std::uint64_t>(
                    global_scalar_count) -
                    gathered[index]
                        .scalar_count ||
            (index > 0U &&
             gathered[index - 1U]
                     .global_id ==
                 gathered[index]
                     .global_id)) {
            PetscLayoutDestroy(
                &layout);
            return PETSC_ERR_ARG_INCOMP;
        }
    }

    std::vector<
        VariableCardinalityNaturalVariableCellDof3D>
        cells;
    cells.reserve(
        local_phase_counts.size());
    std::vector<PetscInt>
        local_scalar_global_indices(
            local_packed_scalar_count,
            PetscInt{-1});

    local_error =
        PETSC_SUCCESS;
    try {
        for (std::size_t local = 0U;
             local < local_phase_counts.size();
             ++local) {
            const auto cell =
                mpmc::mesh::LocalIndex{
                    static_cast<
                        mpmc::mesh::LocalIndex::value_type>(
                            local)};
            const auto global =
                partition.global_id(
                    mpmc::mesh::EntityKind::cell,
                    cell);
            const auto found =
                std::lower_bound(
                    gathered.begin(),
                    gathered.end(),
                    global.value(),
                    [](const auto& record,
                       std::uint64_t target) {
                        return record.global_id <
                            target;
                    });
            if (found ==
                    gathered.end() ||
                found->global_id !=
                    global.value() ||
                found->owner_rank !=
                    partition.owner_rank(
                        mpmc::mesh::EntityKind::cell,
                        cell)
                        .value() ||
                found->scalar_count !=
                    static_cast<std::uint64_t>(
                        local_scalar_counts[local]) ||
                found->global_scalar_start >
                    static_cast<std::uint64_t>(
                        std::numeric_limits<PetscInt>::
                            max())) {
                throw std::invalid_argument(
                    "owner/ghost active-set metadata mismatch");
            }

            const PetscInt global_start =
                static_cast<PetscInt>(
                    found->global_scalar_start);
            cells.push_back(
                {
                    cell,
                    global,
                    partition.owner_rank(
                        mpmc::mesh::EntityKind::cell,
                        cell),
                    local_phase_counts[local],
                    local_scalar_offsets[local],
                    local_scalar_counts[local],
                    global_start});

            for (std::size_t slot = 0U;
                 slot <
                    local_scalar_counts[local];
                 ++slot) {
                local_scalar_global_indices[
                    local_scalar_offsets[local] +
                    slot] =
                    global_start +
                    static_cast<PetscInt>(
                        slot);
            }
        }

        output->emplace(
            partition.local_rank(),
            partition.rank_count(),
            component_count,
            owned_start,
            owned_end,
            global_scalar_count,
            local_packed_scalar_count,
            std::move(cells),
            std::move(
                local_scalar_global_indices));
    } catch (...) {
        local_error =
            PETSC_ERR_ARG_INCOMP;
    }

    PetscLayoutDestroy(
        &layout);

    error =
        collective_error(
            comm,
            local_error);
    if (error != PETSC_SUCCESS) {
        output->reset();
    }
    return error;
}

inline PetscErrorCode
create_variable_cardinality_natural_variable_section_3d(
    MPI_Comm comm,
    const VariableCardinalityNaturalVariableNumbering3D&
        numbering,
    PetscSection* output) {
    if (output == nullptr) {
        return PETSC_ERR_ARG_NULL;
    }
    if (*output != nullptr) {
        return PETSC_ERR_ARG_WRONGSTATE;
    }

    PetscInt cell_count = 0;
    if (!variable_cardinality_numbering_detail::
            checked_petsc_int(
                numbering.local_cell_count(),
                &cell_count)) {
        return PETSC_ERR_ARG_OUTOFRANGE;
    }

    PetscSection section = nullptr;
    PetscErrorCode error =
        PetscSectionCreate(
            comm,
            &section);
    if (error == PETSC_SUCCESS) {
        error =
            PetscSectionSetChart(
                section,
                0,
                cell_count);
    }
    if (error == PETSC_SUCCESS) {
        for (PetscInt point = 0;
             point < cell_count;
             ++point) {
            const auto& cell =
                numbering.cell(
                    mpmc::mesh::LocalIndex{
                        static_cast<
                            mpmc::mesh::LocalIndex::value_type>(
                                point)});
            PetscInt q = 0;
            if (!variable_cardinality_numbering_detail::
                    checked_petsc_int(
                        cell.scalar_count,
                        &q)) {
                error =
                    PETSC_ERR_ARG_OUTOFRANGE;
                break;
            }
            error =
                PetscSectionSetDof(
                    section,
                    point,
                    q);
            if (error != PETSC_SUCCESS) {
                break;
            }
        }
    }
    if (error == PETSC_SUCCESS) {
        error =
            PetscSectionSetUp(
                section);
    }
    if (error == PETSC_SUCCESS) {
        for (PetscInt point = 0;
             point < cell_count;
             ++point) {
            PetscInt offset = -1;
            error =
                PetscSectionGetOffset(
                    section,
                    point,
                    &offset);
            if (error != PETSC_SUCCESS) {
                break;
            }
            const auto& cell =
                numbering.cell(
                    mpmc::mesh::LocalIndex{
                        static_cast<
                            mpmc::mesh::LocalIndex::value_type>(
                                point)});
            if (offset < 0 ||
                static_cast<std::size_t>(
                    offset) !=
                    cell.local_scalar_offset) {
                error =
                    PETSC_ERR_ARG_INCOMP;
                break;
            }
        }
    }

    if (error != PETSC_SUCCESS) {
        PetscSectionDestroy(
            &section);
        return error;
    }

    *output = section;
    return PETSC_SUCCESS;
}

inline PetscErrorCode
create_variable_cardinality_natural_variable_vec_3d(
    MPI_Comm comm,
    const VariableCardinalityNaturalVariableNumbering3D&
        numbering,
    Vec* output) {
    if (output == nullptr) {
        return PETSC_ERR_ARG_NULL;
    }
    if (*output != nullptr) {
        return PETSC_ERR_ARG_WRONGSTATE;
    }

    Vec vector = nullptr;
    PetscErrorCode error =
        VecCreateMPI(
            comm,
            numbering
                .petsc_local_owned_scalar_count(),
            numbering
                .petsc_global_scalar_count(),
            &vector);
    if (error != PETSC_SUCCESS) {
        return error;
    }

    PetscInt local = -1;
    PetscInt global = -1;
    PetscInt start = -1;
    PetscInt end = -1;
    error =
        VecGetLocalSize(
            vector,
            &local);
    if (error == PETSC_SUCCESS) {
        error =
            VecGetSize(
                vector,
                &global);
    }
    if (error == PETSC_SUCCESS) {
        error =
            VecGetOwnershipRange(
                vector,
                &start,
                &end);
    }

    if (error != PETSC_SUCCESS ||
        local !=
            numbering
                .petsc_local_owned_scalar_count() ||
        global !=
            numbering
                .petsc_global_scalar_count() ||
        start !=
            numbering
                .petsc_owned_scalar_start() ||
        end !=
            numbering
                .petsc_owned_scalar_end()) {
        VecDestroy(
            &vector);
        return error != PETSC_SUCCESS
            ? error
            : PETSC_ERR_ARG_INCOMP;
    }

    *output = vector;
    return PETSC_SUCCESS;
}

} // namespace mpmc::flow_discretization_petsc

#endif // MPMC_FLOW_DISCRETIZATION_PETSC_VARIABLE_CARDINALITY_NATURAL_VARIABLE_NUMBERING_HPP
