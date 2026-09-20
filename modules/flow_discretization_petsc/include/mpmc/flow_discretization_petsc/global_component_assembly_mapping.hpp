#ifndef MPMC_FLOW_DISCRETIZATION_PETSC_GLOBAL_COMPONENT_ASSEMBLY_MAPPING_HPP
#define MPMC_FLOW_DISCRETIZATION_PETSC_GLOBAL_COMPONENT_ASSEMBLY_MAPPING_HPP

#include <mpmc/flow_discretization_petsc/distributed_component_conservation.hpp>
#include <mpmc/mesh/dof_numbering.hpp>
#include <mpmc/mesh_petsc/adapter.hpp>

#include <petscmat.h>

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

namespace mpmc::flow_discretization_petsc {

inline constexpr std::string_view
    component_conservation_global_assembly_mapping_convention =
        "flow_discretization_petsc/component-conservation-global-assembly-mapping/v1";

enum class ComponentJacobianCellBlockKind3D {
    diagonal_cell,
    off_diagonal_cell
};

/// One component residual scalar prepared for a future Vec insertion.
///
/// petsc_global_row follows PETSc rank-contiguous scalar ownership.
/// mesh_global_row_dof is the independent mesh-global DoF provenance from
/// DofNumberingSnapshot. The two numbering spaces are deliberately not assumed
/// to be equal.
struct AssemblyReadyComponentResidualEntry3D {
    PetscInt petsc_global_row{-1};
    mpmc::mesh::GlobalDofIndex mesh_global_row_dof{
        mpmc::mesh::GlobalDofIndex::value_type{0}};
    mpmc::mesh::GlobalEntityId row_cell_global{
        mpmc::mesh::GlobalEntityId::value_type{0}};
    std::size_t component_row{};
    double value_mol_per_bulk_m3_s{};
};

/// One dense component-row / natural-variable-column Jacobian scalar prepared
/// for a future Mat insertion.
///
/// Zero-valued entries are retained. This bridge maps numbering only; it does
/// not use numerical zero to infer matrix structure.
struct AssemblyReadyComponentJacobianEntry3D {
    PetscInt petsc_global_row{-1};
    PetscInt petsc_global_column{-1};
    mpmc::mesh::GlobalDofIndex mesh_global_row_dof{
        mpmc::mesh::GlobalDofIndex::value_type{0}};
    mpmc::mesh::GlobalDofIndex mesh_global_column_dof{
        mpmc::mesh::GlobalDofIndex::value_type{0}};
    mpmc::mesh::GlobalEntityId row_cell_global{
        mpmc::mesh::GlobalEntityId::value_type{0}};
    mpmc::mesh::GlobalEntityId column_cell_global{
        mpmc::mesh::GlobalEntityId::value_type{0}};
    std::size_t component_row{};
    std::size_t natural_variable_column{};
    ComponentJacobianCellBlockKind3D block_kind{
        ComponentJacobianCellBlockKind3D::
            diagonal_cell};
    double value{};
};

class ComponentConservationGlobalAssemblyEntries3D {
public:
    static constexpr std::string_view convention =
        component_conservation_global_assembly_mapping_convention;

    ComponentConservationGlobalAssemblyEntries3D(
        mpmc::mesh::PartitionRank local_rank,
        std::uint32_t rank_count,
        std::string natural_variable_id,
        std::size_t component_count,
        std::size_t natural_variable_count,
        PetscInt petsc_scalar_row_start,
        PetscInt petsc_scalar_row_end,
        PetscInt petsc_scalar_row_count,
        std::vector<
            AssemblyReadyComponentResidualEntry3D>
            residual_entries,
        std::vector<
            AssemblyReadyComponentJacobianEntry3D>
            jacobian_entries)
        : local_rank_(local_rank),
          rank_count_(rank_count),
          natural_variable_id_(
              std::move(natural_variable_id)),
          component_count_(component_count),
          natural_variable_count_(
              natural_variable_count),
          petsc_scalar_row_start_(
              petsc_scalar_row_start),
          petsc_scalar_row_end_(
              petsc_scalar_row_end),
          petsc_scalar_row_count_(
              petsc_scalar_row_count),
          residual_entries_(
              std::move(residual_entries)),
          jacobian_entries_(
              std::move(jacobian_entries)) {
        validate();
    }

    [[nodiscard]] mpmc::mesh::PartitionRank
    local_rank() const noexcept {
        return local_rank_;
    }

    [[nodiscard]] std::uint32_t
    rank_count() const noexcept {
        return rank_count_;
    }

    [[nodiscard]] std::string_view
    natural_variable_id() const noexcept {
        return natural_variable_id_;
    }

    [[nodiscard]] std::size_t
    component_count() const noexcept {
        return component_count_;
    }

    [[nodiscard]] std::size_t
    natural_variable_count() const noexcept {
        return natural_variable_count_;
    }

    [[nodiscard]] PetscInt
    petsc_scalar_row_start() const noexcept {
        return petsc_scalar_row_start_;
    }

    [[nodiscard]] PetscInt
    petsc_scalar_row_end() const noexcept {
        return petsc_scalar_row_end_;
    }

    [[nodiscard]] PetscInt
    petsc_scalar_row_count() const noexcept {
        return petsc_scalar_row_count_;
    }

    [[nodiscard]] std::span<
        const AssemblyReadyComponentResidualEntry3D>
    residual_entries() const noexcept {
        return residual_entries_;
    }

    [[nodiscard]] std::span<
        const AssemblyReadyComponentJacobianEntry3D>
    jacobian_entries() const noexcept {
        return jacobian_entries_;
    }

private:
    void validate() const {
        if (rank_count_ == 0U ||
            local_rank_.value() >= rank_count_ ||
            natural_variable_id_.empty() ||
            component_count_ < 2U ||
            natural_variable_count_ == 0U ||
            petsc_scalar_row_start_ < 0 ||
            petsc_scalar_row_end_ <
                petsc_scalar_row_start_ ||
            petsc_scalar_row_count_ <
                petsc_scalar_row_end_) {
            throw std::invalid_argument(
                "mpmc::flow_discretization_petsc: malformed global component assembly mapping");
        }

        PetscInt previous_residual_row = -1;
        for (const auto& entry :
             residual_entries_) {
            if (entry.petsc_global_row <
                    petsc_scalar_row_start_ ||
                entry.petsc_global_row >=
                    petsc_scalar_row_end_ ||
                entry.component_row >=
                    component_count_ ||
                !std::isfinite(
                    entry.value_mol_per_bulk_m3_s) ||
                (previous_residual_row >= 0 &&
                 entry.petsc_global_row <=
                     previous_residual_row)) {
                throw std::invalid_argument(
                    "mpmc::flow_discretization_petsc: invalid or duplicate global component residual entry");
            }
            previous_residual_row =
                entry.petsc_global_row;
        }

        PetscInt previous_row = -1;
        PetscInt previous_column = -1;
        for (const auto& entry :
             jacobian_entries_) {
            if (entry.petsc_global_row <
                    petsc_scalar_row_start_ ||
                entry.petsc_global_row >=
                    petsc_scalar_row_end_ ||
                entry.petsc_global_column < 0 ||
                entry.petsc_global_column >=
                    petsc_scalar_row_count_ ||
                entry.component_row >=
                    component_count_ ||
                entry.natural_variable_column >=
                    natural_variable_count_ ||
                !std::isfinite(entry.value) ||
                (previous_row ==
                     entry.petsc_global_row &&
                 previous_column >=
                     entry.petsc_global_column) ||
                (previous_row >
                 entry.petsc_global_row)) {
                throw std::invalid_argument(
                    "mpmc::flow_discretization_petsc: invalid, unsorted or duplicate global component Jacobian entry");
            }
            previous_row =
                entry.petsc_global_row;
            previous_column =
                entry.petsc_global_column;
        }
    }

    mpmc::mesh::PartitionRank local_rank_;
    std::uint32_t rank_count_;
    std::string natural_variable_id_;
    std::size_t component_count_{};
    std::size_t natural_variable_count_{};
    PetscInt petsc_scalar_row_start_{};
    PetscInt petsc_scalar_row_end_{};
    PetscInt petsc_scalar_row_count_{};
    std::vector<
        AssemblyReadyComponentResidualEntry3D>
        residual_entries_;
    std::vector<
        AssemblyReadyComponentJacobianEntry3D>
        jacobian_entries_;
};

namespace global_component_assembly_mapping_detail {

[[nodiscard]] inline bool checked_scalar_index(
    PetscInt cell_global_row,
    std::size_t block_width,
    std::size_t slot,
    PetscInt* output) {
    if (output == nullptr ||
        cell_global_row < 0 ||
        block_width == 0U ||
        slot >= block_width ||
        block_width >
            static_cast<std::size_t>(
                std::numeric_limits<PetscInt>::max())) {
        return false;
    }

    const PetscInt width =
        static_cast<PetscInt>(
            block_width);
    const PetscInt slot_value =
        static_cast<PetscInt>(
            slot);
    const PetscInt maximum =
        std::numeric_limits<PetscInt>::max();

    if (cell_global_row >
        (maximum - slot_value) /
            width) {
        return false;
    }

    *output =
        cell_global_row * width +
        slot_value;
    return true;
}

[[nodiscard]] inline bool checked_scalar_count(
    PetscInt cell_count,
    std::size_t block_width,
    PetscInt* output) {
    if (output == nullptr ||
        cell_count < 0 ||
        block_width == 0U ||
        block_width >
            static_cast<std::size_t>(
                std::numeric_limits<PetscInt>::max())) {
        return false;
    }
    const PetscInt width =
        static_cast<PetscInt>(
            block_width);
    if (cell_count >
        std::numeric_limits<PetscInt>::max() /
            width) {
        return false;
    }
    *output =
        cell_count * width;
    return true;
}

[[nodiscard]] inline bool contains_cell_column(
    std::span<const PetscInt> columns,
    PetscInt cell_global_row) {
    return std::binary_search(
        columns.begin(),
        columns.end(),
        cell_global_row);
}

inline void validate_cell_pair_pattern(
    const mpmc::flow_discretization::
        OwnedCellComponentConservationRow3D&
            row,
    const mpmc::mesh::PartitionSnapshot&
        partition,
    const mpmc::discretization_petsc::
        PetscMpiAijSymbolicPreallocation3D&
            cell_bridge,
    const mpmc::discretization_petsc::
        OwnedCellStructuralColumnPatternSnapshot3D&
            cell_pattern) {
    const PetscInt row_cell_global =
        cell_bridge.global_row(
            row.cell);
    const auto diagonal =
        cell_pattern.diagonal_global_columns(
            row.cell);
    const auto off_diagonal =
        cell_pattern.off_diagonal_global_columns(
            row.cell);

    if (!contains_cell_column(
            diagonal,
            row_cell_global)) {
        throw std::invalid_argument(
            "mpmc::flow_discretization_petsc: cell structural pattern is missing self column");
    }

    std::size_t owned_neighbour_count = 0U;
    std::size_t ghost_neighbour_count = 0U;

    for (const auto& block :
         row.off_diagonal_cell_pair_blocks) {
        const PetscInt column_cell_global =
            cell_bridge.global_row(
                block.column_cell);
        if (partition.global_id(
                mpmc::mesh::EntityKind::cell,
                block.column_cell) !=
                block.column_cell_global) {
            throw std::invalid_argument(
                "mpmc::flow_discretization_petsc: cell-pair block stable column cell mismatch");
        }

        if (partition.is_owned(
                mpmc::mesh::EntityKind::cell,
                block.column_cell)) {
            ++owned_neighbour_count;
            if (!contains_cell_column(
                    diagonal,
                    column_cell_global)) {
                throw std::invalid_argument(
                    "mpmc::flow_discretization_petsc: owned cell-pair block is absent from diagonal structural pattern");
            }
        } else if (
            partition.is_ghost(
                mpmc::mesh::EntityKind::cell,
                block.column_cell)) {
            ++ghost_neighbour_count;
            if (!contains_cell_column(
                    off_diagonal,
                    column_cell_global)) {
                throw std::invalid_argument(
                    "mpmc::flow_discretization_petsc: ghost cell-pair block is absent from off-diagonal structural pattern");
            }
        } else {
            throw std::invalid_argument(
                "mpmc::flow_discretization_petsc: cell-pair column cell has invalid ownership");
        }
    }

    if (diagonal.size() !=
            owned_neighbour_count + 1U ||
        off_diagonal.size() !=
            ghost_neighbour_count) {
        throw std::invalid_argument(
            "mpmc::flow_discretization_petsc: component Jacobian cell-pair blocks do not exactly match structural column pattern");
    }
}

} // namespace global_component_assembly_mapping_detail

/// Convert complete distributed owned component rows into scalar PETSc
/// row/column/value entries without creating or modifying Mat/Vec objects.
///
/// PETSc scalar numbering is derived from the existing rank-contiguous cell row
/// bridge:
///
///   scalar = cell_global_row * q + natural-variable/equation slot.
///
/// DofNumberingSnapshot remains an independent mesh-global numbering audit. Its
/// indices are stored as provenance but are never substituted for PETSc row
/// ownership.
inline PetscErrorCode
make_component_conservation_global_assembly_entries_3d(
    MPI_Comm comm,
    const DistributedOwnedMultiCellComponentConservationSnapshot3D&
        conservation,
    const mpmc::mesh::PartitionSnapshot& partition,
    const mpmc::mesh::DofLayout& dof_layout,
    const mpmc::mesh::DofNumberingSnapshot&
        dof_numbering,
    const mpmc::discretization_petsc::
        PetscMpiAijSymbolicPreallocation3D&
            cell_bridge,
    const mpmc::discretization_petsc::
        OwnedCellStructuralColumnPatternSnapshot3D&
            cell_pattern,
    std::string_view natural_variable_id,
    std::optional<
        ComponentConservationGlobalAssemblyEntries3D>*
            output) {
    using namespace global_component_assembly_mapping_detail;

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

    std::size_t component_count = 0U;
    std::size_t q = 0U;
    std::size_t variable_index = 0U;
    PetscInt scalar_row_start = -1;
    PetscInt scalar_row_end = -1;
    PetscInt scalar_row_count = -1;
    PetscInt local_scalar_rows = -1;

    try {
        if (conservation.local_rank() !=
                partition.local_rank() ||
            dof_numbering.local_rank() !=
                partition.local_rank() ||
            cell_bridge.local_rank() !=
                partition.local_rank() ||
            cell_pattern.local_rank() !=
                partition.local_rank() ||
            conservation.rank_count() !=
                partition.rank_count() ||
            dof_numbering.rank_count() !=
                partition.rank_count() ||
            cell_bridge.rank_count() !=
                partition.rank_count() ||
            cell_pattern.rank_count() !=
                partition.rank_count() ||
            partition.local_rank().value() !=
                static_cast<std::uint32_t>(
                    mpi_rank) ||
            partition.rank_count() !=
                static_cast<std::uint32_t>(
                    mpi_size) ||
            dof_numbering.local_dof_count() !=
                dof_layout.total_dof_count() ||
            dof_layout.entity_count(
                mpmc::mesh::EntityKind::cell) !=
                partition.entity_count(
                    mpmc::mesh::EntityKind::cell) ||
            conservation.local_cell_count() !=
                partition.entity_count(
                    mpmc::mesh::EntityKind::cell) ||
            cell_pattern.local_cell_count() !=
                partition.entity_count(
                    mpmc::mesh::EntityKind::cell) ||
            cell_pattern.global_row_start() !=
                cell_bridge.global_row_start() ||
            cell_pattern.global_row_end() !=
                cell_bridge.global_row_end() ||
            cell_pattern.global_row_count() !=
                cell_bridge.global_row_count() ||
            dof_numbering.global_entity_count(
                mpmc::mesh::EntityKind::cell) !=
                static_cast<std::uint64_t>(
                    cell_bridge.global_row_count()) ||
            natural_variable_id.empty()) {
            throw std::invalid_argument(
                "global mapping metadata mismatch");
        }

        component_count =
            conservation.component_ids().size();
        if (component_count < 2U) {
            throw std::invalid_argument(
                "global mapping requires at least two components");
        }

        if (!dof_layout.contains(
                natural_variable_id)) {
            throw std::invalid_argument(
                "natural-variable DoF variable is absent");
        }
        variable_index =
            dof_layout.variable_index(
                natural_variable_id);
        const auto& variable =
            dof_layout.variable(
                variable_index);
        q = variable.component_count;
        const std::size_t phase_count =
            q > 1U
                ? (q - 1U) / component_count
                : 0U;
        if (variable.location !=
                mpmc::mesh::EntityKind::cell ||
            q <= 1U ||
            (q - 1U) % component_count != 0U ||
            phase_count == 0U ||
            phase_count >
                mpmc::flow::fixed_three_phase_count ||
            dof_layout.dofs_per_entity(
                mpmc::mesh::EntityKind::cell) !=
                q) {
            throw std::invalid_argument(
                "natural-variable DoF variable must occupy the complete cell scalar block");
        }

        for (std::size_t local = 0U;
             local <
             partition.entity_count(
                 mpmc::mesh::EntityKind::cell);
             ++local) {
            const auto cell =
                mpmc::mesh::LocalIndex{
                    static_cast<
                        mpmc::mesh::LocalIndex::value_type>(
                            local)};
            const std::size_t entity_offset =
                dof_layout.entity_offset(
                    mpmc::mesh::EntityKind::cell,
                    cell);
            if (dof_layout.scalar_offset(
                    variable_index,
                    cell,
                    0U) !=
                    entity_offset ||
                dof_layout.scalar_offset(
                    variable_index,
                    cell,
                    q - 1U) !=
                    entity_offset +
                        q - 1U) {
                throw std::invalid_argument(
                    "natural-variable DoF variable does not span the cell block contiguously");
            }

            for (std::size_t column = 0U;
                 column < q;
                 ++column) {
                const std::size_t local_scalar =
                    dof_layout.scalar_offset(
                        variable_index,
                        cell,
                        column);
                if (dof_numbering.is_owned(
                        local_scalar) !=
                    partition.is_owned(
                        mpmc::mesh::EntityKind::cell,
                        cell)) {
                    throw std::invalid_argument(
                        "DofNumbering ownership disagrees with cell ownership");
                }
            }
        }

        if (!checked_scalar_count(
                cell_bridge.global_row_start(),
                q,
                &scalar_row_start) ||
            !checked_scalar_count(
                cell_bridge.global_row_end(),
                q,
                &scalar_row_end) ||
            !checked_scalar_count(
                cell_bridge.global_row_count(),
                q,
                &scalar_row_count) ||
            !checked_scalar_count(
                cell_bridge.local_owned_row_count(),
                q,
                &local_scalar_rows)) {
            throw std::length_error(
                "PETSc scalar row range overflow");
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

    PetscLayout scalar_layout = nullptr;
    error =
        PetscLayoutCreate(
            comm,
            &scalar_layout);
    if (error == PETSC_SUCCESS) {
        error =
            PetscLayoutSetLocalSize(
                scalar_layout,
                local_scalar_rows);
    }
    if (error == PETSC_SUCCESS) {
        error =
            PetscLayoutSetSize(
                scalar_layout,
                scalar_row_count);
    }
    if (error == PETSC_SUCCESS) {
        error =
            PetscLayoutSetUp(
                scalar_layout);
    }

    PetscInt actual_start = -1;
    PetscInt actual_end = -1;
    if (error == PETSC_SUCCESS) {
        error =
            PetscLayoutGetRange(
                scalar_layout,
                &actual_start,
                &actual_end);
    }
    const PetscErrorCode layout_destroy_error =
        PetscLayoutDestroy(
            &scalar_layout);
    if (error == PETSC_SUCCESS &&
        layout_destroy_error !=
            PETSC_SUCCESS) {
        error =
            layout_destroy_error;
    }
    if (error == PETSC_SUCCESS &&
        (actual_start !=
             scalar_row_start ||
         actual_end !=
             scalar_row_end)) {
        error =
            PETSC_ERR_ARG_INCOMP;
    }

    local_error = error;
    error =
        detail::collective_error(
            comm,
            local_error);
    if (error != PETSC_SUCCESS) {
        return error;
    }

    PetscSection section = nullptr;
    std::vector<PetscInt>
        local_to_mesh_global;
    local_error =
        mpmc::mesh_petsc::
            create_section_mapping(
                comm,
                dof_layout,
                dof_numbering,
                &section,
                &local_to_mesh_global);

    error =
        detail::collective_error(
            comm,
            local_error);
    if (error != PETSC_SUCCESS) {
        if (section != nullptr) {
            (void)PetscSectionDestroy(
                &section);
        }
        return error;
    }

    std::vector<
        AssemblyReadyComponentResidualEntry3D>
        residual_entries;
    std::vector<
        AssemblyReadyComponentJacobianEntry3D>
        jacobian_entries;

    try {
        const std::size_t owned_count =
            conservation.owned_rows().size();
        if (owned_count !=
                partition.owned_count(
                    mpmc::mesh::EntityKind::cell) ||
            owned_count !=
                cell_pattern.row_count()) {
            throw std::invalid_argument(
                "owned row count mismatch");
        }

        if (owned_count >
            std::numeric_limits<std::size_t>::max() /
                component_count) {
            throw std::length_error(
                "residual entry count overflow");
        }
        residual_entries.reserve(
            owned_count *
            component_count);

        const std::vector<std::string>
            canonical_component_ids{
                conservation
                    .component_ids()
                    .begin(),
                conservation
                    .component_ids()
                    .end()};

        for (const auto& row :
             conservation.owned_rows()) {
            if (!partition.is_owned(
                    mpmc::mesh::EntityKind::cell,
                    row.cell) ||
                partition.global_id(
                    mpmc::mesh::EntityKind::cell,
                    row.cell) !=
                    row.cell_global ||
                row.local_residual.component_ids !=
                    canonical_component_ids ||
                row.local_residual
                        .cell_state_identity.layout
                        .unknown_count() !=
                    q ||
                row.local_residual
                        .cell_state_identity.layout
                        .phase_count() !=
                    phase_count ||
                row.local_residual
                        .cell_state_identity.layout
                        .component_count() !=
                    component_count ||
                !cell_pattern.contains_owned_cell(
                    row.cell)) {
                throw std::invalid_argument(
                    "owned component row identity/layout mismatch");
            }

            validate_cell_pair_pattern(
                row,
                partition,
                cell_bridge,
                cell_pattern);

            const PetscInt cell_global_row =
                cell_bridge.global_row(
                    row.cell);

            for (std::size_t component = 0U;
                 component < component_count;
                 ++component) {
                const std::size_t row_slot =
                    row.local_residual
                        .cell_state_identity.layout
                        .component_conservation_equation_index(
                            component);
                PetscInt petsc_row = -1;
                if (!checked_scalar_index(
                        cell_global_row,
                        q,
                        row_slot,
                        &petsc_row) ||
                    petsc_row <
                        scalar_row_start ||
                    petsc_row >=
                        scalar_row_end) {
                    throw std::invalid_argument(
                        "component residual PETSc row is outside local scalar ownership");
                }

                const std::size_t mesh_row_local =
                    dof_layout.scalar_offset(
                        variable_index,
                        row.cell,
                        row_slot);
                if (mesh_row_local >=
                        local_to_mesh_global.size() ||
                    !dof_numbering.is_owned(
                        mesh_row_local) ||
                    local_to_mesh_global[
                        mesh_row_local] < 0 ||
                    static_cast<
                        std::uint64_t>(
                            local_to_mesh_global[
                                mesh_row_local]) !=
                        dof_numbering
                            .global_index(
                                mesh_row_local)
                            .value()) {
                    throw std::invalid_argument(
                        "component residual mesh-global row mapping mismatch");
                }

                const auto mesh_row_global =
                    dof_numbering.global_index(
                        mesh_row_local);

                residual_entries.push_back(
                    AssemblyReadyComponentResidualEntry3D{
                        petsc_row,
                        mesh_row_global,
                        row.cell_global,
                        component,
                        row.local_residual
                            .residual(
                                component)});

                for (std::size_t column = 0U;
                     column < q;
                     ++column) {
                    PetscInt petsc_column = -1;
                    if (!checked_scalar_index(
                            cell_global_row,
                            q,
                            column,
                            &petsc_column)) {
                        throw std::length_error(
                            "diagonal PETSc scalar column overflow");
                    }

                    const std::size_t
                        mesh_column_local =
                            dof_layout.scalar_offset(
                                variable_index,
                                row.cell,
                                column);
                    if (mesh_column_local >=
                            local_to_mesh_global.size() ||
                        local_to_mesh_global[
                            mesh_column_local] < 0 ||
                        static_cast<
                            std::uint64_t>(
                                local_to_mesh_global[
                                    mesh_column_local]) !=
                            dof_numbering
                                .global_index(
                                    mesh_column_local)
                                .value()) {
                        throw std::invalid_argument(
                            "diagonal mesh-global column mapping mismatch");
                    }

                    jacobian_entries.push_back(
                        AssemblyReadyComponentJacobianEntry3D{
                            petsc_row,
                            petsc_column,
                            mesh_row_global,
                            dof_numbering
                                .global_index(
                                    mesh_column_local),
                            row.cell_global,
                            row.cell_global,
                            component,
                            column,
                            ComponentJacobianCellBlockKind3D::
                                diagonal_cell,
                            row.local_residual
                                .d_local(
                                    component,
                                    column)});
                }

                for (const auto& block :
                     row.off_diagonal_cell_pair_blocks) {
                    if (block.column_state_identity
                            .layout.unknown_count() !=
                            q ||
                        block.column_state_identity
                                .component_ids !=
                            row.local_residual
                                .component_ids ||
                        partition.global_id(
                            mpmc::mesh::EntityKind::cell,
                            block.column_cell) !=
                            block.column_cell_global) {
                        throw std::invalid_argument(
                            "off-diagonal cell-pair natural-variable identity mismatch");
                    }

                    const PetscInt
                        column_cell_global_row =
                            cell_bridge.global_row(
                                block.column_cell);

                    for (std::size_t column = 0U;
                         column < q;
                         ++column) {
                        PetscInt petsc_column = -1;
                        if (!checked_scalar_index(
                                column_cell_global_row,
                                q,
                                column,
                                &petsc_column)) {
                            throw std::length_error(
                                "off-diagonal PETSc scalar column overflow");
                        }

                        const std::size_t
                            mesh_column_local =
                                dof_layout.scalar_offset(
                                    variable_index,
                                    block.column_cell,
                                    column);
                        if (mesh_column_local >=
                                local_to_mesh_global.size() ||
                            local_to_mesh_global[
                                mesh_column_local] < 0 ||
                            static_cast<
                                std::uint64_t>(
                                    local_to_mesh_global[
                                        mesh_column_local]) !=
                                dof_numbering
                                    .global_index(
                                        mesh_column_local)
                                    .value() ||
                            dof_numbering
                                    .is_owned(
                                        mesh_column_local) !=
                                partition.is_owned(
                                    mpmc::mesh::EntityKind::cell,
                                    block.column_cell)) {
                            throw std::invalid_argument(
                                "off-diagonal mesh-global column mapping/ownership mismatch");
                        }

                        jacobian_entries.push_back(
                            AssemblyReadyComponentJacobianEntry3D{
                                petsc_row,
                                petsc_column,
                                mesh_row_global,
                                dof_numbering
                                    .global_index(
                                        mesh_column_local),
                                row.cell_global,
                                block.column_cell_global,
                                component,
                                column,
                                ComponentJacobianCellBlockKind3D::
                                    off_diagonal_cell,
                                block.d_component(
                                    component,
                                    column)});
                    }
                }
            }
        }

        std::sort(
            residual_entries.begin(),
            residual_entries.end(),
            [](const auto& left,
               const auto& right) {
                return left.petsc_global_row <
                    right.petsc_global_row;
            });

        std::sort(
            jacobian_entries.begin(),
            jacobian_entries.end(),
            [](const auto& left,
               const auto& right) {
                return left.petsc_global_row <
                           right.petsc_global_row ||
                    (left.petsc_global_row ==
                         right.petsc_global_row &&
                     left.petsc_global_column <
                         right.petsc_global_column);
            });

        for (std::size_t index = 1U;
             index < jacobian_entries.size();
             ++index) {
            if (jacobian_entries[index - 1U]
                        .petsc_global_row ==
                    jacobian_entries[index]
                        .petsc_global_row &&
                jacobian_entries[index - 1U]
                        .petsc_global_column ==
                    jacobian_entries[index]
                        .petsc_global_column) {
                throw std::invalid_argument(
                    "duplicate scalar Jacobian row/column after cell-pair coalescing");
            }
        }
    } catch (...) {
        local_error =
            PETSC_ERR_ARG_INCOMP;
    }

    const PetscErrorCode
        section_destroy_error =
            PetscSectionDestroy(
                &section);
    if (local_error == PETSC_SUCCESS &&
        section_destroy_error !=
            PETSC_SUCCESS) {
        local_error =
            section_destroy_error;
    }

    error =
        detail::collective_error(
            comm,
            local_error);
    if (error != PETSC_SUCCESS) {
        return error;
    }

    try {
        output->emplace(
            partition.local_rank(),
            partition.rank_count(),
            std::string{
                natural_variable_id},
            component_count,
            q,
            scalar_row_start,
            scalar_row_end,
            scalar_row_count,
            std::move(residual_entries),
            std::move(jacobian_entries));
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

#endif // MPMC_FLOW_DISCRETIZATION_PETSC_GLOBAL_COMPONENT_ASSEMBLY_MAPPING_HPP
