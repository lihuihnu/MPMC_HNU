#ifndef MPMC_FLOW_DISCRETIZATION_PETSC_FUGACITY_EQUILIBRIUM_GLOBAL_ASSEMBLY_MAPPING_HPP
#define MPMC_FLOW_DISCRETIZATION_PETSC_FUGACITY_EQUILIBRIUM_GLOBAL_ASSEMBLY_MAPPING_HPP

#include <mpmc/flow/fugacity_equilibrium_linearization.hpp>
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
    fugacity_equilibrium_global_assembly_mapping_convention =
        "flow_discretization_petsc/fugacity-equilibrium-global-assembly-mapping/v1";

struct OwnedCellFugacityEquilibriumLinearizationBinding3D {
    mpmc::mesh::LocalIndex cell{
        mpmc::mesh::LocalIndex::value_type{0}};
    mpmc::mesh::GlobalEntityId cell_global{
        mpmc::mesh::GlobalEntityId::value_type{0}};
    const mpmc::flow::
        FugacityEquilibriumResidualLinearization3P*
            linearization{};
};

struct AssemblyReadyFugacityResidualEntry3D {
    PetscInt petsc_global_row{-1};
    mpmc::mesh::GlobalDofIndex mesh_global_row_dof{
        mpmc::mesh::GlobalDofIndex::value_type{0}};
    mpmc::mesh::GlobalEntityId row_cell_global{
        mpmc::mesh::GlobalEntityId::value_type{0}};
    mpmc::flow::PhaseSlot3 non_reference_phase{
        mpmc::flow::PhaseSlot3::phase1};
    std::size_t component{};
    double value{};
};

struct AssemblyReadyFugacityJacobianEntry3D {
    PetscInt petsc_global_row{-1};
    PetscInt petsc_global_column{-1};
    mpmc::mesh::GlobalDofIndex mesh_global_row_dof{
        mpmc::mesh::GlobalDofIndex::value_type{0}};
    mpmc::mesh::GlobalDofIndex mesh_global_column_dof{
        mpmc::mesh::GlobalDofIndex::value_type{0}};
    mpmc::mesh::GlobalEntityId row_cell_global{
        mpmc::mesh::GlobalEntityId::value_type{0}};
    mpmc::flow::PhaseSlot3 non_reference_phase{
        mpmc::flow::PhaseSlot3::phase1};
    std::size_t component{};
    std::size_t natural_variable_column{};
    double value{};
};

class FugacityEquilibriumGlobalAssemblyEntries3D {
public:
    static constexpr std::string_view convention =
        fugacity_equilibrium_global_assembly_mapping_convention;

    FugacityEquilibriumGlobalAssemblyEntries3D(
        mpmc::mesh::PartitionRank local_rank,
        std::uint32_t rank_count,
        std::string natural_variable_id,
        std::size_t component_count,
        std::size_t natural_variable_count,
        PetscInt petsc_scalar_row_start,
        PetscInt petsc_scalar_row_end,
        PetscInt petsc_scalar_row_count,
        std::vector<
            AssemblyReadyFugacityResidualEntry3D>
            residual_entries,
        std::vector<
            AssemblyReadyFugacityJacobianEntry3D>
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
        const AssemblyReadyFugacityResidualEntry3D>
    residual_entries() const noexcept {
        return residual_entries_;
    }

    [[nodiscard]] std::span<
        const AssemblyReadyFugacityJacobianEntry3D>
    jacobian_entries() const noexcept {
        return jacobian_entries_;
    }

private:
    static bool valid_phase(
        mpmc::flow::PhaseSlot3 phase) {
        return phase ==
                   mpmc::flow::PhaseSlot3::phase1 ||
            phase ==
                   mpmc::flow::PhaseSlot3::phase2;
    }

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
                "mpmc::flow_discretization_petsc: malformed fugacity global assembly mapping");
        }

        PetscInt previous_residual_row = -1;
        for (const auto& entry :
             residual_entries_) {
            if (entry.petsc_global_row <
                    petsc_scalar_row_start_ ||
                entry.petsc_global_row >=
                    petsc_scalar_row_end_ ||
                !valid_phase(
                    entry.non_reference_phase) ||
                entry.component >=
                    component_count_ ||
                !std::isfinite(
                    entry.value) ||
                (previous_residual_row >= 0 &&
                 entry.petsc_global_row <=
                     previous_residual_row)) {
                throw std::invalid_argument(
                    "mpmc::flow_discretization_petsc: invalid or duplicate fugacity residual assembly entry");
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
                !valid_phase(
                    entry.non_reference_phase) ||
                entry.component >=
                    component_count_ ||
                entry.natural_variable_column >=
                    natural_variable_count_ ||
                !std::isfinite(
                    entry.value) ||
                (previous_row ==
                     entry.petsc_global_row &&
                 previous_column >=
                     entry.petsc_global_column) ||
                previous_row >
                    entry.petsc_global_row) {
                throw std::invalid_argument(
                    "mpmc::flow_discretization_petsc: invalid, unsorted or duplicate fugacity Jacobian assembly entry");
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
        AssemblyReadyFugacityResidualEntry3D>
        residual_entries_;
    std::vector<
        AssemblyReadyFugacityJacobianEntry3D>
        jacobian_entries_;
};

namespace fugacity_global_assembly_mapping_detail {

[[nodiscard]] inline bool
same_chart(
    const mpmc::flow::NaturalVariableLayout3P&
        first,
    const mpmc::flow::NaturalVariableLayout3P&
        second) {
    return first.component_count() ==
               second.component_count() &&
        first.unknown_count() ==
            second.unknown_count() &&
        first.composition_pivot()
                .dependent_components() ==
            second.composition_pivot()
                .dependent_components();
}

} // namespace fugacity_global_assembly_mapping_detail

/// Map already-computed local fugacity-equilibrium residual/Jacobian blocks
/// onto the same square natural-variable PETSc scalar blocks used by component
/// conservation.
///
/// Only diagonal cell blocks are emitted: local thermodynamic equilibrium has
/// no spatial neighbour derivative in this contract.
inline PetscErrorCode
make_fugacity_equilibrium_global_assembly_entries_3d(
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
    std::span<
        const OwnedCellFugacityEquilibriumLinearizationBinding3D>
        fugacity_bindings,
    std::optional<
        FugacityEquilibriumGlobalAssemblyEntries3D>*
            output) {
    using namespace global_component_assembly_mapping_detail;
    using namespace fugacity_global_assembly_mapping_detail;

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
    std::vector<std::optional<std::size_t>>
        binding_by_local_cell;

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
            natural_variable_id.empty() ||
            fugacity_bindings.size() !=
                partition.owned_count(
                    mpmc::mesh::EntityKind::cell) ||
            conservation.owned_rows().size() !=
                fugacity_bindings.size()) {
            throw std::invalid_argument(
                "fugacity global mapping metadata mismatch");
        }

        component_count =
            conservation.component_ids().size();
        q =
            mpmc::flow::NaturalVariableLayout3P{
                component_count}
                .unknown_count();

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
        if (variable.location !=
                mpmc::mesh::EntityKind::cell ||
            variable.component_count != q ||
            dof_layout.dofs_per_entity(
                mpmc::mesh::EntityKind::cell) !=
                q) {
            throw std::invalid_argument(
                "natural-variable DoF variable must occupy the complete cell scalar block");
        }

        binding_by_local_cell.assign(
            partition.entity_count(
                mpmc::mesh::EntityKind::cell),
            std::nullopt);

        const std::vector<std::string>
            canonical_ids{
                conservation.component_ids().begin(),
                conservation.component_ids().end()};

        for (std::size_t index = 0U;
             index < fugacity_bindings.size();
             ++index) {
            const auto& binding =
                fugacity_bindings[index];
            const std::size_t local =
                static_cast<std::size_t>(
                    binding.cell.value());
            if (local >=
                    binding_by_local_cell.size() ||
                binding_by_local_cell[local]
                    .has_value() ||
                binding.linearization ==
                    nullptr ||
                !partition.is_owned(
                    mpmc::mesh::EntityKind::cell,
                    binding.cell) ||
                partition.global_id(
                    mpmc::mesh::EntityKind::cell,
                    binding.cell) !=
                    binding.cell_global ||
                !cell_pattern.contains_owned_cell(
                    binding.cell)) {
                throw std::invalid_argument(
                    "invalid or duplicate owned fugacity binding");
            }

            const auto& linearization =
                *binding.linearization;
            const auto& component_row =
                conservation.owned_row(
                    binding.cell);
            if (component_row.cell_global !=
                    binding.cell_global ||
                linearization.component_ids().size() !=
                    canonical_ids.size() ||
                !std::equal(
                    linearization.component_ids().begin(),
                    linearization.component_ids().end(),
                    canonical_ids.begin()) ||
                linearization.component_count() !=
                    component_count ||
                linearization.input_count() !=
                    q ||
                !same_chart(
                    linearization.layout(),
                    component_row.local_residual
                        .cell_state_identity.layout) ||
                linearization.layout()
                        .composition_pivot()
                        .dependent_components() !=
                    component_row.local_residual
                        .cell_state_identity.layout
                        .composition_pivot()
                        .dependent_components()) {
                throw std::invalid_argument(
                    "fugacity linearization does not match owned component-row chart");
            }

            const PetscInt cell_global_row =
                cell_bridge.global_row(
                    binding.cell);
            const auto diagonal_columns =
                cell_pattern.diagonal_global_columns(
                    binding.cell);
            if (!std::binary_search(
                    diagonal_columns.begin(),
                    diagonal_columns.end(),
                    cell_global_row)) {
                throw std::invalid_argument(
                    "fugacity owned cell is missing structural self column");
            }

            for (std::size_t column = 0U;
                 column < q;
                 ++column) {
                const std::size_t local_scalar =
                    dof_layout.scalar_offset(
                        variable_index,
                        binding.cell,
                        column);
                if (!dof_numbering.is_owned(
                        local_scalar)) {
                    throw std::invalid_argument(
                        "owned fugacity cell has non-owned natural-variable scalar");
                }
            }

            binding_by_local_cell[local] =
                index;
        }

        for (const auto& row :
             conservation.owned_rows()) {
            const std::size_t local =
                static_cast<std::size_t>(
                    row.cell.value());
            if (local >=
                    binding_by_local_cell.size() ||
                !binding_by_local_cell[local]
                    .has_value()) {
                throw std::invalid_argument(
                    "owned component row is missing fugacity linearization");
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
                "fugacity PETSc scalar row range overflow");
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
        AssemblyReadyFugacityResidualEntry3D>
        residual_entries;
    std::vector<
        AssemblyReadyFugacityJacobianEntry3D>
        jacobian_entries;

    try {
        const std::size_t owned_count =
            fugacity_bindings.size();
        if (owned_count >
            std::numeric_limits<std::size_t>::max() /
                (2U * component_count)) {
            throw std::length_error(
                "fugacity residual entry count overflow");
        }
        residual_entries.reserve(
            owned_count *
            2U *
            component_count);

        for (const auto& component_row :
             conservation.owned_rows()) {
            const std::size_t local =
                static_cast<std::size_t>(
                    component_row.cell.value());
            const auto& binding =
                fugacity_bindings[
                    *binding_by_local_cell[
                        local]];
            const auto& linearization =
                *binding.linearization;
            const PetscInt cell_global_row =
                cell_bridge.global_row(
                    binding.cell);

            for (std::size_t phase = 1U;
                 phase <
                 mpmc::flow::fixed_three_phase_count;
                 ++phase) {
                const auto slot =
                    static_cast<
                        mpmc::flow::PhaseSlot3>(
                            phase);
                for (std::size_t component = 0U;
                     component < component_count;
                     ++component) {
                    const std::size_t row_slot =
                        linearization.equation_index(
                            slot,
                            component);
                    const auto row_identity =
                        linearization.layout()
                            .fugacity_equilibrium_row_identity(
                                row_slot);
                    if (!row_identity.has_value() ||
                        row_identity
                                ->non_reference_phase !=
                            slot ||
                        row_identity->component !=
                            component ||
                        row_slot <= component_count ||
                        row_slot >= q) {
                        throw std::invalid_argument(
                            "fugacity equation slot identity mismatch");
                    }

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
                            "fugacity PETSc row is outside local scalar ownership");
                    }

                    const std::size_t
                        mesh_row_local =
                            dof_layout.scalar_offset(
                                variable_index,
                                binding.cell,
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
                            "fugacity mesh-global row mapping mismatch");
                    }

                    const auto mesh_row_global =
                        dof_numbering.global_index(
                            mesh_row_local);
                    residual_entries.push_back(
                        AssemblyReadyFugacityResidualEntry3D{
                            petsc_row,
                            mesh_row_global,
                            binding.cell_global,
                            slot,
                            component,
                            linearization.residual(
                                slot,
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
                                "fugacity PETSc column overflow");
                        }

                        const std::size_t
                            mesh_column_local =
                                dof_layout.scalar_offset(
                                    variable_index,
                                    binding.cell,
                                    column);
                        if (mesh_column_local >=
                                local_to_mesh_global.size() ||
                            !dof_numbering.is_owned(
                                mesh_column_local) ||
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
                                "fugacity mesh-global column mapping mismatch");
                        }

                        jacobian_entries.push_back(
                            AssemblyReadyFugacityJacobianEntry3D{
                                petsc_row,
                                petsc_column,
                                mesh_row_global,
                                dof_numbering
                                    .global_index(
                                        mesh_column_local),
                                binding.cell_global,
                                slot,
                                component,
                                column,
                                linearization.d_residual(
                                    slot,
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
                    "duplicate fugacity scalar row/column");
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

#endif // MPMC_FLOW_DISCRETIZATION_PETSC_FUGACITY_EQUILIBRIUM_GLOBAL_ASSEMBLY_MAPPING_HPP
