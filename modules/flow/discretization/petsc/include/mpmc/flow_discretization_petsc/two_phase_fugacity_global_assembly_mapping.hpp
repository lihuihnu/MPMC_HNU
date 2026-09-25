#ifndef MPMC_FLOW_DISCRETIZATION_PETSC_TWO_PHASE_FUGACITY_GLOBAL_ASSEMBLY_MAPPING_HPP
#define MPMC_FLOW_DISCRETIZATION_PETSC_TWO_PHASE_FUGACITY_GLOBAL_ASSEMBLY_MAPPING_HPP

#include <mpmc/flow/two_phase_natural_variable.hpp>
#include <mpmc/flow_discretization_petsc/fugacity_equilibrium_global_assembly_mapping.hpp>

#include <petscsys.h>

#include <algorithm>
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
    two_phase_fugacity_global_assembly_mapping_convention =
        "flow_discretization_petsc/two-phase-fugacity-global-assembly-mapping/v1";

struct OwnedCellTwoPhaseFugacityLinearizationBinding3D {
    mpmc::mesh::LocalIndex cell{
        mpmc::mesh::LocalIndex::value_type{0}};
    mpmc::mesh::GlobalEntityId cell_global{
        mpmc::mesh::GlobalEntityId::value_type{0}};
    const mpmc::flow::
        TwoPhaseFugacityEquilibriumLinearization*
            linearization{};
};

inline PetscErrorCode
make_two_phase_fugacity_global_assembly_entries_3d(
    MPI_Comm comm,
    const ComponentConservationGlobalAssemblyEntries3D&
        component_global,
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
        const OwnedCellTwoPhaseFugacityLinearizationBinding3D>
        fugacity_bindings,
    std::optional<
        FugacityEquilibriumGlobalAssemblyEntries3D>*
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

    std::size_t variable_index = 0U;
    std::size_t n = 0U;
    std::size_t q = 0U;
    std::vector<std::optional<std::size_t>>
        binding_by_local_cell;

    try {
        n = component_global.component_count();
        q = component_global.natural_variable_count();
        if (n < 2U ||
            q != 2U * n + 1U ||
            component_global.local_rank() !=
                partition.local_rank() ||
            component_global.local_rank() !=
                cell_bridge.local_rank() ||
            component_global.local_rank() !=
                cell_pattern.local_rank() ||
            component_global.rank_count() !=
                partition.rank_count() ||
            component_global.rank_count() !=
                cell_bridge.rank_count() ||
            component_global.rank_count() !=
                cell_pattern.rank_count() ||
            partition.local_rank().value() !=
                static_cast<std::uint32_t>(
                    mpi_rank) ||
            partition.rank_count() !=
                static_cast<std::uint32_t>(
                    mpi_size) ||
            dof_numbering.local_rank() !=
                partition.local_rank() ||
            dof_numbering.rank_count() !=
                partition.rank_count() ||
            dof_numbering.local_dof_count() !=
                dof_layout.total_dof_count() ||
            natural_variable_id.empty() ||
            component_global.natural_variable_id() !=
                natural_variable_id ||
            fugacity_bindings.size() !=
                partition.owned_count(
                    mpmc::mesh::EntityKind::cell) ||
            !dof_layout.contains(
                natural_variable_id)) {
            throw std::invalid_argument(
                "two-phase fugacity global mapping metadata mismatch");
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
                "two-phase fugacity DoF block must occupy the complete 2*Nc+1 cell block");
        }

        binding_by_local_cell.assign(
            partition.entity_count(
                mpmc::mesh::EntityKind::cell),
            std::nullopt);

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
                binding.linearization == nullptr ||
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
                    "invalid or duplicate two-phase fugacity binding");
            }

            const auto& linearization =
                *binding.linearization;
            if (linearization.layout.phase_count() !=
                    2U ||
                linearization.layout.component_count() !=
                    n ||
                linearization.layout.unknown_count() !=
                    q ||
                linearization.component_ids.size() !=
                    n ||
                linearization.residual_count() !=
                    n ||
                linearization.input_count !=
                    q ||
                linearization.jacobian.size() !=
                    n * q) {
                throw std::invalid_argument(
                    "two-phase fugacity linearization shape/layout mismatch");
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
                        "owned two-phase fugacity cell contains non-owned scalar");
                }
            }

            binding_by_local_cell[local] =
                index;
        }

        for (const auto& cell :
             cell_bridge
                 .owned_cells_in_petsc_row_order()) {
            const std::size_t local =
                static_cast<std::size_t>(
                    cell.value());
            if (local >=
                    binding_by_local_cell.size() ||
                !binding_by_local_cell[local]
                     .has_value()) {
                throw std::invalid_argument(
                    "owned two-phase PETSc cell lacks fugacity binding");
            }
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

    std::vector<AssemblyReadyFugacityResidualEntry3D>
        residual_entries;
    std::vector<AssemblyReadyFugacityJacobianEntry3D>
        jacobian_entries;

    try {
        const std::size_t owned_count =
            fugacity_bindings.size();
        if (owned_count >
            std::numeric_limits<std::size_t>::max() /
                n) {
            throw std::length_error(
                "two-phase fugacity residual count overflow");
        }
        residual_entries.reserve(
            owned_count * n);
        if (owned_count >
            std::numeric_limits<std::size_t>::max() /
                n ||
            owned_count * n >
            std::numeric_limits<std::size_t>::max() /
                q) {
            throw std::length_error(
                "two-phase fugacity Jacobian count overflow");
        }
        jacobian_entries.reserve(
            owned_count * n * q);

        for (const auto& cell :
             cell_bridge
                 .owned_cells_in_petsc_row_order()) {
            const std::size_t local =
                static_cast<std::size_t>(
                    cell.value());
            const auto& binding =
                fugacity_bindings[
                    *binding_by_local_cell[local]];
            const auto& linearization =
                *binding.linearization;
            const PetscInt cell_global_row =
                cell_bridge.global_row(
                    binding.cell);

            for (std::size_t component = 0U;
                 component < n;
                 ++component) {
                const std::size_t row_slot =
                    n + 1U + component;
                PetscInt petsc_row = -1;
                if (!checked_scalar_index(
                        cell_global_row,
                        q,
                        row_slot,
                        &petsc_row) ||
                    petsc_row <
                        component_global
                            .petsc_scalar_row_start() ||
                    petsc_row >=
                        component_global
                            .petsc_scalar_row_end()) {
                    throw std::invalid_argument(
                        "two-phase fugacity PETSc row outside local ownership");
                }

                const std::size_t row_local =
                    dof_layout.scalar_offset(
                        variable_index,
                        binding.cell,
                        row_slot);
                const auto mesh_row =
                    dof_numbering.global_index(
                        row_local);

                residual_entries.push_back(
                    {
                        petsc_row,
                        mesh_row,
                        binding.cell_global,
                        mpmc::flow::PhaseSlot3::phase1,
                        component,
                        linearization.residual[
                            component]});

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
                            "two-phase fugacity PETSc column overflow");
                    }
                    const std::size_t column_local =
                        dof_layout.scalar_offset(
                            variable_index,
                            binding.cell,
                            column);
                    jacobian_entries.push_back(
                        {
                            petsc_row,
                            petsc_column,
                            mesh_row,
                            dof_numbering.global_index(
                                column_local),
                            binding.cell_global,
                            mpmc::flow::PhaseSlot3::phase1,
                            component,
                            column,
                            linearization.d_residual(
                                component,
                                column)});
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

        output->emplace(
            component_global.local_rank(),
            component_global.rank_count(),
            std::string{
                component_global
                    .natural_variable_id()},
            n,
            q,
            component_global.petsc_scalar_row_start(),
            component_global.petsc_scalar_row_end(),
            component_global.petsc_scalar_row_count(),
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
    }
    return error;
}

} // namespace mpmc::flow_discretization_petsc

#endif // MPMC_FLOW_DISCRETIZATION_PETSC_TWO_PHASE_FUGACITY_GLOBAL_ASSEMBLY_MAPPING_HPP
