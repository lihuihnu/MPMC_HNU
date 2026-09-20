#ifndef MPMC_FLOW_DISCRETIZATION_PETSC_ENERGY_GLOBAL_ASSEMBLY_MAPPING_HPP
#define MPMC_FLOW_DISCRETIZATION_PETSC_ENERGY_GLOBAL_ASSEMBLY_MAPPING_HPP

#include <mpmc/flow_discretization_petsc/distributed_energy_conservation.hpp>
#include <mpmc/flow_discretization_petsc/global_component_assembly_mapping.hpp>

#include <petscsys.h>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
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
    energy_global_assembly_mapping_convention =
        "flow_discretization_petsc/energy-global-assembly-mapping/v1";

enum class EnergyJacobianCellBlockKind3D : std::uint8_t {
    diagonal_cell,
    off_diagonal_cell
};

struct AssemblyReadyEnergyResidualEntry3D {
    PetscInt petsc_global_row{-1};
    mpmc::mesh::GlobalDofIndex mesh_global_row_dof{
        mpmc::mesh::GlobalDofIndex::value_type{0}};
    mpmc::mesh::GlobalEntityId row_cell_global{
        mpmc::mesh::GlobalEntityId::value_type{0}};
    double value_w_per_bulk_m3{};
};

struct AssemblyReadyEnergyJacobianEntry3D {
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
    std::size_t natural_variable_column{};
    EnergyJacobianCellBlockKind3D block_kind{
        EnergyJacobianCellBlockKind3D::diagonal_cell};
    double value{};
};

class EnergyConservationGlobalAssemblyEntries3D {
public:
    static constexpr std::string_view convention =
        energy_global_assembly_mapping_convention;

    EnergyConservationGlobalAssemblyEntries3D(
        mpmc::mesh::PartitionRank local_rank,
        std::uint32_t rank_count,
        std::string natural_variable_id,
        std::size_t component_count,
        std::size_t natural_variable_count,
        PetscInt petsc_scalar_row_start,
        PetscInt petsc_scalar_row_end,
        PetscInt petsc_scalar_row_count,
        std::vector<AssemblyReadyEnergyResidualEntry3D>
            residual_entries,
        std::vector<AssemblyReadyEnergyJacobianEntry3D>
            jacobian_entries)
        : local_rank_(local_rank),
          rank_count_(rank_count),
          natural_variable_id_(std::move(natural_variable_id)),
          component_count_(component_count),
          natural_variable_count_(natural_variable_count),
          petsc_scalar_row_start_(petsc_scalar_row_start),
          petsc_scalar_row_end_(petsc_scalar_row_end),
          petsc_scalar_row_count_(petsc_scalar_row_count),
          residual_entries_(std::move(residual_entries)),
          jacobian_entries_(std::move(jacobian_entries)) {
        validate();
    }

    [[nodiscard]] mpmc::mesh::PartitionRank local_rank() const noexcept {
        return local_rank_;
    }
    [[nodiscard]] std::uint32_t rank_count() const noexcept {
        return rank_count_;
    }
    [[nodiscard]] std::string_view natural_variable_id() const noexcept {
        return natural_variable_id_;
    }
    [[nodiscard]] std::size_t component_count() const noexcept {
        return component_count_;
    }
    [[nodiscard]] std::size_t natural_variable_count() const noexcept {
        return natural_variable_count_;
    }
    [[nodiscard]] PetscInt petsc_scalar_row_start() const noexcept {
        return petsc_scalar_row_start_;
    }
    [[nodiscard]] PetscInt petsc_scalar_row_end() const noexcept {
        return petsc_scalar_row_end_;
    }
    [[nodiscard]] PetscInt petsc_scalar_row_count() const noexcept {
        return petsc_scalar_row_count_;
    }
    [[nodiscard]] std::span<const AssemblyReadyEnergyResidualEntry3D>
    residual_entries() const noexcept {
        return residual_entries_;
    }
    [[nodiscard]] std::span<const AssemblyReadyEnergyJacobianEntry3D>
    jacobian_entries() const noexcept {
        return jacobian_entries_;
    }

private:
    void validate() const {
        if (rank_count_ == 0U ||
            local_rank_.value() >= rank_count_ ||
            natural_variable_id_.empty() ||
            component_count_ < 2U ||
            natural_variable_count_ !=
                mpmc::flow::fixed_three_phase_count *
                    component_count_ +
                    1U ||
            petsc_scalar_row_start_ < 0 ||
            petsc_scalar_row_end_ < petsc_scalar_row_start_ ||
            petsc_scalar_row_count_ < petsc_scalar_row_end_) {
            throw std::invalid_argument(
                "mpmc::flow_discretization_petsc: malformed energy global assembly mapping");
        }

        PetscInt previous_residual_row = -1;
        for (const auto& entry : residual_entries_) {
            if (entry.petsc_global_row < petsc_scalar_row_start_ ||
                entry.petsc_global_row >= petsc_scalar_row_end_ ||
                !std::isfinite(entry.value_w_per_bulk_m3) ||
                (previous_residual_row >= 0 &&
                 entry.petsc_global_row <= previous_residual_row)) {
                throw std::invalid_argument(
                    "mpmc::flow_discretization_petsc: invalid or duplicate global energy residual entry");
            }
            previous_residual_row = entry.petsc_global_row;
        }

        PetscInt previous_row = -1;
        PetscInt previous_column = -1;
        for (const auto& entry : jacobian_entries_) {
            if (entry.petsc_global_row < petsc_scalar_row_start_ ||
                entry.petsc_global_row >= petsc_scalar_row_end_ ||
                entry.petsc_global_column < 0 ||
                entry.petsc_global_column >= petsc_scalar_row_count_ ||
                entry.natural_variable_column >= natural_variable_count_ ||
                !std::isfinite(entry.value) ||
                previous_row > entry.petsc_global_row ||
                (previous_row == entry.petsc_global_row &&
                 previous_column >= entry.petsc_global_column)) {
                throw std::invalid_argument(
                    "mpmc::flow_discretization_petsc: invalid, unsorted or duplicate global energy Jacobian entry");
            }
            previous_row = entry.petsc_global_row;
            previous_column = entry.petsc_global_column;
        }
    }

    mpmc::mesh::PartitionRank local_rank_;
    std::uint32_t rank_count_{};
    std::string natural_variable_id_;
    std::size_t component_count_{};
    std::size_t natural_variable_count_{};
    PetscInt petsc_scalar_row_start_{};
    PetscInt petsc_scalar_row_end_{};
    PetscInt petsc_scalar_row_count_{};
    std::vector<AssemblyReadyEnergyResidualEntry3D>
        residual_entries_;
    std::vector<AssemblyReadyEnergyJacobianEntry3D>
        jacobian_entries_;
};

namespace energy_global_assembly_mapping_detail {

inline void validate_cell_pair_pattern(
    const OwnedCellEnergyConservationRow3D& row,
    const mpmc::mesh::PartitionSnapshot& partition,
    const mpmc::discretization_petsc::
        PetscMpiAijSymbolicPreallocation3D& cell_bridge,
    const mpmc::discretization_petsc::
        OwnedCellStructuralColumnPatternSnapshot3D& cell_pattern) {
    using global_component_assembly_mapping_detail::
        contains_cell_column;

    const PetscInt row_cell_global =
        cell_bridge.global_row(row.cell);
    const auto diagonal =
        cell_pattern.diagonal_global_columns(row.cell);
    const auto off_diagonal =
        cell_pattern.off_diagonal_global_columns(row.cell);

    if (!contains_cell_column(diagonal, row_cell_global)) {
        throw std::invalid_argument(
            "mpmc::flow_discretization_petsc: energy structural pattern is missing self column");
    }

    std::size_t owned_neighbour_count = 0U;
    std::size_t ghost_neighbour_count = 0U;
    for (const auto& block : row.off_diagonal_cell_pair_blocks) {
        if (partition.global_id(
                mpmc::mesh::EntityKind::cell,
                block.column_cell) !=
                block.column_cell_global ||
            block.gradient.size() !=
                block.column_state_identity.layout.unknown_count()) {
            throw std::invalid_argument(
                "mpmc::flow_discretization_petsc: malformed energy cell-pair block identity/shape");
        }

        const PetscInt column_cell_global =
            cell_bridge.global_row(block.column_cell);
        if (partition.is_owned(
                mpmc::mesh::EntityKind::cell,
                block.column_cell)) {
            ++owned_neighbour_count;
            if (!contains_cell_column(diagonal, column_cell_global)) {
                throw std::invalid_argument(
                    "mpmc::flow_discretization_petsc: owned energy neighbour is absent from diagonal structural pattern");
            }
        } else if (partition.is_ghost(
                       mpmc::mesh::EntityKind::cell,
                       block.column_cell)) {
            ++ghost_neighbour_count;
            if (!contains_cell_column(off_diagonal, column_cell_global)) {
                throw std::invalid_argument(
                    "mpmc::flow_discretization_petsc: ghost energy neighbour is absent from off-diagonal structural pattern");
            }
        } else {
            throw std::invalid_argument(
                "mpmc::flow_discretization_petsc: energy neighbour has invalid ownership");
        }
    }

    if (diagonal.size() != owned_neighbour_count + 1U ||
        off_diagonal.size() != ghost_neighbour_count) {
        throw std::invalid_argument(
            "mpmc::flow_discretization_petsc: energy cell-pair blocks do not exactly match structural column pattern");
    }
}

} // namespace energy_global_assembly_mapping_detail

/// Map complete distributed owned energy rows onto the same square
/// natural-variable scalar blocks used by component and fugacity rows.
///
/// The energy equation occupies NaturalVariableLayout3P::energy_equation_index()
/// (slot Nc). No Mat or Vec is created or modified.
inline PetscErrorCode
make_energy_conservation_global_assembly_entries_3d(
    MPI_Comm comm,
    const DistributedOwnedMultiCellEnergyConservationSnapshot3D&
        conservation,
    const mpmc::mesh::PartitionSnapshot& partition,
    const mpmc::mesh::DofLayout& dof_layout,
    const mpmc::mesh::DofNumberingSnapshot& dof_numbering,
    const mpmc::discretization_petsc::
        PetscMpiAijSymbolicPreallocation3D& cell_bridge,
    const mpmc::discretization_petsc::
        OwnedCellStructuralColumnPatternSnapshot3D& cell_pattern,
    std::string_view natural_variable_id,
    std::optional<EnergyConservationGlobalAssemblyEntries3D>*
        output) {
    using namespace global_component_assembly_mapping_detail;
    using namespace energy_global_assembly_mapping_detail;

    int mpi_rank = -1;
    int mpi_size = -1;
    if (MPI_Comm_rank(comm, &mpi_rank) != MPI_SUCCESS ||
        MPI_Comm_size(comm, &mpi_size) != MPI_SUCCESS ||
        mpi_rank < 0 ||
        mpi_size <= 0) {
        return PETSC_ERR_MPI;
    }

    PetscErrorCode local_error =
        output == nullptr ? PETSC_ERR_ARG_NULL : PETSC_SUCCESS;
    PetscErrorCode error =
        detail::collective_error(comm, local_error);
    if (error != PETSC_SUCCESS) {
        return error;
    }
    output->reset();

    std::size_t q = 0U;
    std::size_t component_count = 0U;
    std::size_t variable_index = 0U;
    PetscInt scalar_row_start = -1;
    PetscInt scalar_row_end = -1;
    PetscInt scalar_row_count = -1;
    PetscInt local_scalar_rows = -1;

    try {
        if (conservation.local_rank() != partition.local_rank() ||
            conservation.rank_count() != partition.rank_count() ||
            dof_numbering.local_rank() != partition.local_rank() ||
            dof_numbering.rank_count() != partition.rank_count() ||
            cell_bridge.local_rank() != partition.local_rank() ||
            cell_bridge.rank_count() != partition.rank_count() ||
            cell_pattern.local_rank() != partition.local_rank() ||
            cell_pattern.rank_count() != partition.rank_count() ||
            partition.local_rank().value() !=
                static_cast<std::uint32_t>(mpi_rank) ||
            partition.rank_count() !=
                static_cast<std::uint32_t>(mpi_size) ||
            conservation.local_cell_count() !=
                partition.entity_count(mpmc::mesh::EntityKind::cell) ||
            cell_pattern.local_cell_count() !=
                partition.entity_count(mpmc::mesh::EntityKind::cell) ||
            dof_numbering.local_dof_count() !=
                dof_layout.total_dof_count() ||
            dof_layout.entity_count(mpmc::mesh::EntityKind::cell) !=
                partition.entity_count(mpmc::mesh::EntityKind::cell) ||
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
                "energy global mapping metadata mismatch");
        }

        if (!dof_layout.contains(natural_variable_id)) {
            throw std::invalid_argument(
                "natural-variable DoF variable is absent");
        }
        variable_index =
            dof_layout.variable_index(natural_variable_id);
        const auto& variable =
            dof_layout.variable(variable_index);
        if (variable.location != mpmc::mesh::EntityKind::cell ||
            variable.component_count < 7U ||
            dof_layout.dofs_per_entity(
                mpmc::mesh::EntityKind::cell) !=
                variable.component_count) {
            throw std::invalid_argument(
                "energy mapping natural-variable DoF must occupy the complete cell scalar block");
        }

        q = variable.component_count;
        if ((q - 1U) %
                mpmc::flow::fixed_three_phase_count !=
            0U) {
            throw std::invalid_argument(
                "energy mapping natural-variable block width is not 3*Nc+1");
        }
        component_count =
            (q - 1U) /
            mpmc::flow::fixed_three_phase_count;
        if (component_count < 2U) {
            throw std::invalid_argument(
                "energy mapping requires at least two components");
        }

        for (std::size_t local = 0U;
             local < partition.entity_count(
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
                    variable_index, cell, 0U) !=
                    entity_offset ||
                dof_layout.scalar_offset(
                    variable_index, cell, q - 1U) !=
                    entity_offset + q - 1U) {
                throw std::invalid_argument(
                    "natural-variable DoF variable does not span the cell block contiguously");
            }

            for (std::size_t column = 0U;
                 column < q;
                 ++column) {
                const std::size_t local_scalar =
                    dof_layout.scalar_offset(
                        variable_index, cell, column);
                if (dof_numbering.is_owned(local_scalar) !=
                    partition.is_owned(
                        mpmc::mesh::EntityKind::cell,
                        cell)) {
                    throw std::invalid_argument(
                        "DofNumbering ownership disagrees with energy cell ownership");
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
                "energy PETSc scalar row range overflow");
        }
    } catch (...) {
        local_error = PETSC_ERR_ARG_INCOMP;
    }

    error = detail::collective_error(comm, local_error);
    if (error != PETSC_SUCCESS) {
        return error;
    }

    PetscLayout scalar_layout = nullptr;
    error = PetscLayoutCreate(comm, &scalar_layout);
    if (error == PETSC_SUCCESS) {
        error = PetscLayoutSetLocalSize(
            scalar_layout,
            local_scalar_rows);
    }
    if (error == PETSC_SUCCESS) {
        error = PetscLayoutSetSize(
            scalar_layout,
            scalar_row_count);
    }
    if (error == PETSC_SUCCESS) {
        error = PetscLayoutSetUp(scalar_layout);
    }

    PetscInt actual_start = -1;
    PetscInt actual_end = -1;
    if (error == PETSC_SUCCESS) {
        error = PetscLayoutGetRange(
            scalar_layout,
            &actual_start,
            &actual_end);
    }
    const PetscErrorCode layout_destroy_error =
        PetscLayoutDestroy(&scalar_layout);
    if (error == PETSC_SUCCESS &&
        layout_destroy_error != PETSC_SUCCESS) {
        error = layout_destroy_error;
    }
    if (error == PETSC_SUCCESS &&
        (actual_start != scalar_row_start ||
         actual_end != scalar_row_end)) {
        error = PETSC_ERR_ARG_INCOMP;
    }

    local_error = error;
    error = detail::collective_error(comm, local_error);
    if (error != PETSC_SUCCESS) {
        return error;
    }

    PetscSection section = nullptr;
    std::vector<PetscInt> local_to_mesh_global;
    local_error =
        mpmc::mesh_petsc::create_section_mapping(
            comm,
            dof_layout,
            dof_numbering,
            &section,
            &local_to_mesh_global);
    error = detail::collective_error(comm, local_error);
    if (error != PETSC_SUCCESS) {
        if (section != nullptr) {
            (void)PetscSectionDestroy(&section);
        }
        return error;
    }

    std::vector<AssemblyReadyEnergyResidualEntry3D>
        residual_entries;
    std::vector<AssemblyReadyEnergyJacobianEntry3D>
        jacobian_entries;

    try {
        const std::size_t owned_count =
            conservation.owned_rows().size();
        if (owned_count !=
                partition.owned_count(
                    mpmc::mesh::EntityKind::cell) ||
            owned_count != cell_pattern.row_count()) {
            throw std::invalid_argument(
                "energy owned row count mismatch");
        }

        residual_entries.reserve(owned_count);

        for (const auto& row : conservation.owned_rows()) {
            const auto& identity =
                row.local_residual.cell_state_identity;
            if (!partition.is_owned(
                    mpmc::mesh::EntityKind::cell,
                    row.cell) ||
                partition.global_id(
                    mpmc::mesh::EntityKind::cell,
                    row.cell) != row.cell_global ||
                identity.layout.component_count() !=
                    component_count ||
                identity.component_ids.size() !=
                    component_count ||
                identity.layout.unknown_count() != q ||
                row.local_residual.local_gradient.size() != q ||
                !cell_pattern.contains_owned_cell(row.cell)) {
                throw std::invalid_argument(
                    "owned energy row identity/layout mismatch");
            }

            validate_cell_pair_pattern(
                row,
                partition,
                cell_bridge,
                cell_pattern);

            const std::size_t row_slot =
                identity.layout.energy_equation_index();
            if (row_slot != component_count ||
                row_slot >= q) {
                throw std::invalid_argument(
                    "energy equation row slot mismatch");
            }

            const PetscInt cell_global_row =
                cell_bridge.global_row(row.cell);
            PetscInt petsc_row = -1;
            if (!checked_scalar_index(
                    cell_global_row,
                    q,
                    row_slot,
                    &petsc_row) ||
                petsc_row < scalar_row_start ||
                petsc_row >= scalar_row_end) {
                throw std::invalid_argument(
                    "energy PETSc row is outside local scalar ownership");
            }

            const std::size_t mesh_row_local =
                dof_layout.scalar_offset(
                    variable_index,
                    row.cell,
                    row_slot);
            if (mesh_row_local >= local_to_mesh_global.size() ||
                !dof_numbering.is_owned(mesh_row_local) ||
                local_to_mesh_global[mesh_row_local] < 0 ||
                static_cast<std::uint64_t>(
                    local_to_mesh_global[mesh_row_local]) !=
                    dof_numbering.global_index(
                        mesh_row_local).value()) {
                throw std::invalid_argument(
                    "energy mesh-global row mapping mismatch");
            }
            const auto mesh_row_global =
                dof_numbering.global_index(mesh_row_local);

            residual_entries.push_back(
                AssemblyReadyEnergyResidualEntry3D{
                    petsc_row,
                    mesh_row_global,
                    row.cell_global,
                    row.local_residual.residual_w_per_bulk_m3});

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
                        "energy diagonal PETSc scalar column overflow");
                }

                const std::size_t mesh_column_local =
                    dof_layout.scalar_offset(
                        variable_index,
                        row.cell,
                        column);
                if (mesh_column_local >= local_to_mesh_global.size() ||
                    local_to_mesh_global[mesh_column_local] < 0 ||
                    static_cast<std::uint64_t>(
                        local_to_mesh_global[mesh_column_local]) !=
                        dof_numbering.global_index(
                            mesh_column_local).value()) {
                    throw std::invalid_argument(
                        "energy diagonal mesh-global column mapping mismatch");
                }

                jacobian_entries.push_back(
                    AssemblyReadyEnergyJacobianEntry3D{
                        petsc_row,
                        petsc_column,
                        mesh_row_global,
                        dof_numbering.global_index(
                            mesh_column_local),
                        row.cell_global,
                        row.cell_global,
                        column,
                        EnergyJacobianCellBlockKind3D::diagonal_cell,
                        row.local_residual.d_local(column)});
            }

            for (const auto& block :
                 row.off_diagonal_cell_pair_blocks) {
                if (block.column_state_identity.component_ids !=
                        identity.component_ids ||
                    block.column_state_identity.layout
                            .component_count() != component_count ||
                    block.column_state_identity.layout
                            .unknown_count() != q ||
                    block.gradient.size() != q ||
                    partition.global_id(
                        mpmc::mesh::EntityKind::cell,
                        block.column_cell) !=
                        block.column_cell_global) {
                    throw std::invalid_argument(
                        "energy off-diagonal cell-pair identity/layout mismatch");
                }

                const PetscInt column_cell_global_row =
                    cell_bridge.global_row(block.column_cell);
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
                            "energy off-diagonal PETSc scalar column overflow");
                    }

                    const std::size_t mesh_column_local =
                        dof_layout.scalar_offset(
                            variable_index,
                            block.column_cell,
                            column);
                    if (mesh_column_local >= local_to_mesh_global.size() ||
                        local_to_mesh_global[mesh_column_local] < 0 ||
                        static_cast<std::uint64_t>(
                            local_to_mesh_global[mesh_column_local]) !=
                            dof_numbering.global_index(
                                mesh_column_local).value() ||
                        dof_numbering.is_owned(mesh_column_local) !=
                            partition.is_owned(
                                mpmc::mesh::EntityKind::cell,
                                block.column_cell)) {
                        throw std::invalid_argument(
                            "energy off-diagonal mesh-global column mapping/ownership mismatch");
                    }

                    jacobian_entries.push_back(
                        AssemblyReadyEnergyJacobianEntry3D{
                            petsc_row,
                            petsc_column,
                            mesh_row_global,
                            dof_numbering.global_index(
                                mesh_column_local),
                            row.cell_global,
                            block.column_cell_global,
                            column,
                            EnergyJacobianCellBlockKind3D::off_diagonal_cell,
                            block.d_residual(column)});
                }
            }
        }

        std::sort(
            residual_entries.begin(),
            residual_entries.end(),
            [](const auto& left, const auto& right) {
                return left.petsc_global_row <
                    right.petsc_global_row;
            });
        std::sort(
            jacobian_entries.begin(),
            jacobian_entries.end(),
            [](const auto& left, const auto& right) {
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
            if (jacobian_entries[index - 1U].petsc_global_row ==
                    jacobian_entries[index].petsc_global_row &&
                jacobian_entries[index - 1U].petsc_global_column ==
                    jacobian_entries[index].petsc_global_column) {
                throw std::invalid_argument(
                    "duplicate energy scalar row/column");
            }
        }
    } catch (...) {
        local_error = PETSC_ERR_ARG_INCOMP;
    }

    const PetscErrorCode section_destroy_error =
        PetscSectionDestroy(&section);
    if (local_error == PETSC_SUCCESS &&
        section_destroy_error != PETSC_SUCCESS) {
        local_error = section_destroy_error;
    }

    error = detail::collective_error(comm, local_error);
    if (error != PETSC_SUCCESS) {
        return error;
    }

    try {
        output->emplace(
            partition.local_rank(),
            partition.rank_count(),
            std::string{natural_variable_id},
            component_count,
            q,
            scalar_row_start,
            scalar_row_end,
            scalar_row_count,
            std::move(residual_entries),
            std::move(jacobian_entries));
    } catch (...) {
        local_error = PETSC_ERR_ARG_INCOMP;
    }

    error = detail::collective_error(comm, local_error);
    if (error != PETSC_SUCCESS) {
        output->reset();
        return error;
    }
    return PETSC_SUCCESS;
}

} // namespace mpmc::flow_discretization_petsc

#endif // MPMC_FLOW_DISCRETIZATION_PETSC_ENERGY_GLOBAL_ASSEMBLY_MAPPING_HPP
