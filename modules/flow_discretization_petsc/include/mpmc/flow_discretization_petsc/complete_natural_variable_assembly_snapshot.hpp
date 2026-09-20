#ifndef MPMC_FLOW_DISCRETIZATION_PETSC_COMPLETE_NATURAL_VARIABLE_ASSEMBLY_SNAPSHOT_HPP
#define MPMC_FLOW_DISCRETIZATION_PETSC_COMPLETE_NATURAL_VARIABLE_ASSEMBLY_SNAPSHOT_HPP

#include <mpmc/flow_discretization_petsc/energy_global_assembly_mapping.hpp>
#include <mpmc/flow_discretization_petsc/fugacity_equilibrium_global_assembly_mapping.hpp>

#include <petscsys.h>

#include <algorithm>
#include <array>
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
    complete_natural_variable_assembly_snapshot_convention =
        "flow_discretization_petsc/complete-natural-variable-assembly-snapshot/v1";

/// Equation family determines the native residual/Jacobian row units:
/// - component_conservation: mol / (bulk-m^3 s)
/// - energy_conservation: W / bulk-m^3
/// - fugacity_equilibrium: dimensionless log-fugacity residual
enum class NaturalVariableEquationKind3D : std::uint8_t {
    component_conservation,
    energy_conservation,
    fugacity_equilibrium
};

struct CompleteNaturalVariableResidualEntry3D {
    PetscInt petsc_global_row{-1};
    mpmc::mesh::GlobalDofIndex mesh_global_row_dof{
        mpmc::mesh::GlobalDofIndex::value_type{0}};
    mpmc::mesh::GlobalEntityId row_cell_global{
        mpmc::mesh::GlobalEntityId::value_type{0}};
    std::size_t equation_slot{};
    NaturalVariableEquationKind3D equation_kind{
        NaturalVariableEquationKind3D::component_conservation};
    double native_value{};
};

struct CompleteNaturalVariableJacobianEntry3D {
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
    std::size_t equation_slot{};
    std::size_t natural_variable_column{};
    NaturalVariableEquationKind3D equation_kind{
        NaturalVariableEquationKind3D::component_conservation};
    double value{};
};

class CompleteNaturalVariableAssemblySnapshot3D {
public:
    static constexpr std::string_view convention =
        complete_natural_variable_assembly_snapshot_convention;

    CompleteNaturalVariableAssemblySnapshot3D(
        mpmc::mesh::PartitionRank local_rank,
        std::uint32_t rank_count,
        std::string natural_variable_id,
        std::size_t component_count,
        std::size_t natural_variable_count,
        PetscInt petsc_scalar_row_start,
        PetscInt petsc_scalar_row_end,
        PetscInt petsc_scalar_row_count,
        std::vector<CompleteNaturalVariableResidualEntry3D>
            residual_entries,
        std::vector<CompleteNaturalVariableJacobianEntry3D>
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
        const CompleteNaturalVariableResidualEntry3D>
    residual_entries() const noexcept {
        return residual_entries_;
    }

    [[nodiscard]] std::span<
        const CompleteNaturalVariableJacobianEntry3D>
    jacobian_entries() const noexcept {
        return jacobian_entries_;
    }

private:
    void validate() const {
        const std::size_t phase_count =
            natural_variable_count_ > 1U &&
                    component_count_ > 0U
                ? (natural_variable_count_ - 1U) /
                      component_count_
                : 0U;
        if (rank_count_ == 0U ||
            local_rank_.value() >= rank_count_ ||
            natural_variable_id_.empty() ||
            component_count_ < 2U ||
            natural_variable_count_ <= 1U ||
            (natural_variable_count_ - 1U) %
                    component_count_ !=
                0U ||
            phase_count == 0U ||
            phase_count >
                mpmc::flow::fixed_three_phase_count ||
            petsc_scalar_row_start_ < 0 ||
            petsc_scalar_row_end_ < petsc_scalar_row_start_ ||
            petsc_scalar_row_count_ < petsc_scalar_row_end_) {
            throw std::invalid_argument(
                "mpmc::flow_discretization_petsc: malformed complete natural-variable assembly snapshot");
        }

        const PetscInt local_rows =
            petsc_scalar_row_end_ -
            petsc_scalar_row_start_;
        if (local_rows < 0 ||
            static_cast<std::size_t>(local_rows) !=
                residual_entries_.size()) {
            throw std::invalid_argument(
                "mpmc::flow_discretization_petsc: complete residual set does not match local scalar ownership");
        }

        for (std::size_t index = 0U;
             index < residual_entries_.size();
             ++index) {
            const auto& entry =
                residual_entries_[index];
            const PetscInt expected_row =
                petsc_scalar_row_start_ +
                static_cast<PetscInt>(index);
            const std::size_t row_slot =
                static_cast<std::size_t>(
                    entry.petsc_global_row %
                    static_cast<PetscInt>(
                        natural_variable_count_));
            const bool kind_matches =
                (row_slot < component_count_ &&
                 entry.equation_kind ==
                     NaturalVariableEquationKind3D::
                         component_conservation) ||
                (row_slot == component_count_ &&
                 entry.equation_kind ==
                     NaturalVariableEquationKind3D::
                         energy_conservation) ||
                (row_slot > component_count_ &&
                 entry.equation_kind ==
                     NaturalVariableEquationKind3D::
                         fugacity_equilibrium);
            if (entry.petsc_global_row != expected_row ||
                entry.equation_slot != row_slot ||
                !kind_matches ||
                !std::isfinite(entry.native_value)) {
                throw std::invalid_argument(
                    "mpmc::flow_discretization_petsc: complete residual rows are not contiguous or equation semantics mismatch");
            }
        }

        PetscInt previous_row = -1;
        PetscInt previous_column = -1;
        for (const auto& entry :
             jacobian_entries_) {
            const std::size_t row_slot =
                static_cast<std::size_t>(
                    entry.petsc_global_row %
                    static_cast<PetscInt>(
                        natural_variable_count_));
            const bool kind_matches =
                (row_slot < component_count_ &&
                 entry.equation_kind ==
                     NaturalVariableEquationKind3D::
                         component_conservation) ||
                (row_slot == component_count_ &&
                 entry.equation_kind ==
                     NaturalVariableEquationKind3D::
                         energy_conservation) ||
                (row_slot > component_count_ &&
                 entry.equation_kind ==
                     NaturalVariableEquationKind3D::
                         fugacity_equilibrium);
            if (entry.petsc_global_row <
                    petsc_scalar_row_start_ ||
                entry.petsc_global_row >=
                    petsc_scalar_row_end_ ||
                entry.petsc_global_column < 0 ||
                entry.petsc_global_column >=
                    petsc_scalar_row_count_ ||
                entry.equation_slot != row_slot ||
                !kind_matches ||
                entry.natural_variable_column >=
                    natural_variable_count_ ||
                !std::isfinite(entry.value) ||
                previous_row > entry.petsc_global_row ||
                (previous_row ==
                     entry.petsc_global_row &&
                 previous_column >=
                     entry.petsc_global_column)) {
                throw std::invalid_argument(
                    "mpmc::flow_discretization_petsc: complete Jacobian triplets are invalid, unsorted or duplicated");
            }
            previous_row =
                entry.petsc_global_row;
            previous_column =
                entry.petsc_global_column;
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
    std::vector<CompleteNaturalVariableResidualEntry3D>
        residual_entries_;
    std::vector<CompleteNaturalVariableJacobianEntry3D>
        jacobian_entries_;
};

namespace complete_assembly_detail {

[[nodiscard]] inline std::size_t
fugacity_equation_slot(
    mpmc::flow::PhaseSlot3 phase,
    std::size_t component,
    std::size_t component_count) {
    const auto phase_index =
        static_cast<std::size_t>(phase);
    if ((phase_index != 1U &&
         phase_index != 2U) ||
        component >= component_count) {
        throw std::invalid_argument(
            "mpmc::flow_discretization_petsc: invalid fugacity equation identity in complete assembly");
    }
    return component_count +
        1U +
        (phase_index - 1U) *
            component_count +
        component;
}

[[nodiscard]] inline const mpmc::discretization_petsc::
    OwnedCellStructuralColumnPatternRow3D*
find_pattern_row(
    const mpmc::discretization_petsc::
        OwnedCellStructuralColumnPatternSnapshot3D&
            pattern,
    PetscInt cell_global_row) {
    for (const auto& row : pattern.rows()) {
        if (row.global_row ==
            cell_global_row) {
            return &row;
        }
    }
    return nullptr;
}

[[nodiscard]] inline std::vector<PetscInt>
expand_cell_columns(
    std::span<const PetscInt> cell_columns,
    std::size_t q) {
    std::vector<PetscInt> result;
    if (q == 0U ||
        cell_columns.size() >
            std::numeric_limits<std::size_t>::max() /
                q) {
        throw std::length_error(
            "mpmc::flow_discretization_petsc: expanded scalar column pattern size overflow");
    }
    result.reserve(
        cell_columns.size() * q);

    for (const PetscInt cell_column :
         cell_columns) {
        for (std::size_t slot = 0U;
             slot < q;
             ++slot) {
            PetscInt scalar = -1;
            if (!global_component_assembly_mapping_detail::
                    checked_scalar_index(
                        cell_column,
                        q,
                        slot,
                        &scalar)) {
                throw std::length_error(
                    "mpmc::flow_discretization_petsc: expanded scalar column index overflow");
            }
            result.push_back(scalar);
        }
    }
    return result;
}

[[nodiscard]] inline std::vector<PetscInt>
expected_columns(
    NaturalVariableEquationKind3D kind,
    const mpmc::discretization_petsc::
        OwnedCellStructuralColumnPatternSnapshot3D&
            pattern,
    mpmc::mesh::LocalIndex row_cell,
    PetscInt row_cell_global,
    std::size_t q) {
    if (kind ==
        NaturalVariableEquationKind3D::
            fugacity_equilibrium) {
        const std::array<PetscInt, 1>
            self{row_cell_global};
        return expand_cell_columns(
            self,
            q);
    }

    const auto diagonal =
        pattern.diagonal_global_columns(
            row_cell);
    const auto off_diagonal =
        pattern.off_diagonal_global_columns(
            row_cell);

    std::vector<PetscInt> cell_columns;
    cell_columns.reserve(
        diagonal.size() +
        off_diagonal.size());
    cell_columns.insert(
        cell_columns.end(),
        diagonal.begin(),
        diagonal.end());
    cell_columns.insert(
        cell_columns.end(),
        off_diagonal.begin(),
        off_diagonal.end());
    std::sort(
        cell_columns.begin(),
        cell_columns.end());
    cell_columns.erase(
        std::unique(
            cell_columns.begin(),
            cell_columns.end()),
        cell_columns.end());

    return expand_cell_columns(
        cell_columns,
        q);
}

inline void validate_mapping_metadata(
    const ComponentConservationGlobalAssemblyEntries3D&
        component,
    const EnergyConservationGlobalAssemblyEntries3D&
        energy,
    const FugacityEquilibriumGlobalAssemblyEntries3D&
        fugacity) {
    if (component.local_rank() != energy.local_rank() ||
        component.local_rank() != fugacity.local_rank() ||
        component.rank_count() != energy.rank_count() ||
        component.rank_count() != fugacity.rank_count() ||
        component.natural_variable_id() !=
            energy.natural_variable_id() ||
        component.natural_variable_id() !=
            fugacity.natural_variable_id() ||
        component.component_count() !=
            energy.component_count() ||
        component.component_count() !=
            fugacity.component_count() ||
        component.natural_variable_count() !=
            energy.natural_variable_count() ||
        component.natural_variable_count() !=
            fugacity.natural_variable_count() ||
        component.petsc_scalar_row_start() !=
            energy.petsc_scalar_row_start() ||
        component.petsc_scalar_row_start() !=
            fugacity.petsc_scalar_row_start() ||
        component.petsc_scalar_row_end() !=
            energy.petsc_scalar_row_end() ||
        component.petsc_scalar_row_end() !=
            fugacity.petsc_scalar_row_end() ||
        component.petsc_scalar_row_count() !=
            energy.petsc_scalar_row_count() ||
        component.petsc_scalar_row_count() !=
            fugacity.petsc_scalar_row_count()) {
        throw std::invalid_argument(
            "mpmc::flow_discretization_petsc: component/energy/fugacity mapping metadata disagree");
    }
}

} // namespace complete_assembly_detail

/// Combine the three already-validated equation families into one complete
/// assembly-ready Newton-system snapshot.
///
/// No PETSc Mat/Vec object is created or modified.
inline PetscErrorCode
make_complete_natural_variable_assembly_snapshot_3d(
    MPI_Comm comm,
    const ComponentConservationGlobalAssemblyEntries3D&
        component,
    const EnergyConservationGlobalAssemblyEntries3D&
        energy,
    const FugacityEquilibriumGlobalAssemblyEntries3D&
        fugacity,
    const mpmc::mesh::PartitionSnapshot& partition,
    const mpmc::discretization_petsc::
        PetscMpiAijSymbolicPreallocation3D&
            cell_bridge,
    const mpmc::discretization_petsc::
        OwnedCellStructuralColumnPatternSnapshot3D&
            cell_pattern,
    std::optional<
        CompleteNaturalVariableAssemblySnapshot3D>*
            output) {
    using namespace complete_assembly_detail;

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

    std::vector<CompleteNaturalVariableResidualEntry3D>
        residual_entries;
    std::vector<CompleteNaturalVariableJacobianEntry3D>
        jacobian_entries;

    try {
        validate_mapping_metadata(
            component,
            energy,
            fugacity);

        const std::size_t q =
            component.natural_variable_count();
        const std::size_t n =
            component.component_count();
        const std::size_t phase_count =
            q > 1U && n > 0U
                ? (q - 1U) / n
                : 0U;
        if (q <= 1U ||
            (q - 1U) % n != 0U ||
            phase_count == 0U ||
            phase_count >
                mpmc::flow::fixed_three_phase_count ||
            (phase_count == 1U &&
             (!fugacity.residual_entries().empty() ||
              !fugacity.jacobian_entries().empty())) ||
            component.local_rank() !=
                partition.local_rank() ||
            component.local_rank() !=
                cell_bridge.local_rank() ||
            component.local_rank() !=
                cell_pattern.local_rank() ||
            component.rank_count() !=
                partition.rank_count() ||
            component.rank_count() !=
                cell_bridge.rank_count() ||
            component.rank_count() !=
                cell_pattern.rank_count() ||
            partition.local_rank().value() !=
                static_cast<std::uint32_t>(
                    mpi_rank) ||
            partition.rank_count() !=
                static_cast<std::uint32_t>(
                    mpi_size) ||
            cell_pattern.local_cell_count() !=
                partition.entity_count(
                    mpmc::mesh::EntityKind::cell) ||
            cell_pattern.global_row_start() !=
                cell_bridge.global_row_start() ||
            cell_pattern.global_row_end() !=
                cell_bridge.global_row_end() ||
            cell_pattern.global_row_count() !=
                cell_bridge.global_row_count()) {
            throw std::invalid_argument(
                "complete assembly partition/symbolic metadata mismatch");
        }

        PetscInt expected_scalar_start = -1;
        PetscInt expected_scalar_end = -1;
        PetscInt expected_scalar_count = -1;
        if (!global_component_assembly_mapping_detail::
                checked_scalar_count(
                    cell_bridge.global_row_start(),
                    q,
                    &expected_scalar_start) ||
            !global_component_assembly_mapping_detail::
                checked_scalar_count(
                    cell_bridge.global_row_end(),
                    q,
                    &expected_scalar_end) ||
            !global_component_assembly_mapping_detail::
                checked_scalar_count(
                    cell_bridge.global_row_count(),
                    q,
                    &expected_scalar_count) ||
            expected_scalar_start !=
                component.petsc_scalar_row_start() ||
            expected_scalar_end !=
                component.petsc_scalar_row_end() ||
            expected_scalar_count !=
                component.petsc_scalar_row_count()) {
            throw std::invalid_argument(
                "complete assembly scalar ownership does not match cell bridge");
        }

        const PetscInt local_scalar_rows =
            expected_scalar_end -
            expected_scalar_start;
        if (local_scalar_rows < 0) {
            throw std::invalid_argument(
                "complete assembly local scalar row count is negative");
        }
        const std::size_t expected_residual_count =
            static_cast<std::size_t>(
                local_scalar_rows);

        if (component.residual_entries().size() >
                std::numeric_limits<std::size_t>::max() -
                    energy.residual_entries().size() ||
            component.residual_entries().size() +
                    energy.residual_entries().size() >
                std::numeric_limits<std::size_t>::max() -
                    fugacity.residual_entries().size()) {
            throw std::length_error(
                "complete assembly residual count overflow");
        }
        residual_entries.reserve(
            component.residual_entries().size() +
            energy.residual_entries().size() +
            fugacity.residual_entries().size());

        for (const auto& source :
             component.residual_entries()) {
            const std::size_t slot =
                source.component_row;
            residual_entries.push_back(
                CompleteNaturalVariableResidualEntry3D{
                    source.petsc_global_row,
                    source.mesh_global_row_dof,
                    source.row_cell_global,
                    slot,
                    NaturalVariableEquationKind3D::
                        component_conservation,
                    source.value_mol_per_bulk_m3_s});
        }

        for (const auto& source :
             energy.residual_entries()) {
            residual_entries.push_back(
                CompleteNaturalVariableResidualEntry3D{
                    source.petsc_global_row,
                    source.mesh_global_row_dof,
                    source.row_cell_global,
                    n,
                    NaturalVariableEquationKind3D::
                        energy_conservation,
                    source.value_w_per_bulk_m3});
        }

        for (const auto& source :
             fugacity.residual_entries()) {
            const std::size_t slot =
                fugacity_equation_slot(
                    source.non_reference_phase,
                    source.component,
                    n);
            residual_entries.push_back(
                CompleteNaturalVariableResidualEntry3D{
                    source.petsc_global_row,
                    source.mesh_global_row_dof,
                    source.row_cell_global,
                    slot,
                    NaturalVariableEquationKind3D::
                        fugacity_equilibrium,
                    source.value});
        }

        std::sort(
            residual_entries.begin(),
            residual_entries.end(),
            [](const auto& left,
               const auto& right) {
                return left.petsc_global_row <
                    right.petsc_global_row;
            });

        if (residual_entries.size() !=
                expected_residual_count) {
            throw std::invalid_argument(
                "complete assembly residual row count has a gap or collision");
        }

        for (std::size_t index = 0U;
             index < residual_entries.size();
             ++index) {
            const auto& entry =
                residual_entries[index];
            const PetscInt expected_row =
                expected_scalar_start +
                static_cast<PetscInt>(index);
            if (entry.petsc_global_row !=
                    expected_row ||
                entry.equation_slot !=
                    static_cast<std::size_t>(
                        expected_row %
                        static_cast<PetscInt>(q))) {
                throw std::invalid_argument(
                    "complete assembly residual rows do not cover the owned scalar range exactly once");
            }

            const PetscInt cell_global_row =
                expected_row /
                static_cast<PetscInt>(q);
            const auto* pattern_row =
                find_pattern_row(
                    cell_pattern,
                    cell_global_row);
            if (pattern_row == nullptr ||
                pattern_row->cell_global !=
                    entry.row_cell_global ||
                !partition.is_owned(
                    mpmc::mesh::EntityKind::cell,
                    pattern_row->cell) ||
                partition.global_id(
                    mpmc::mesh::EntityKind::cell,
                    pattern_row->cell) !=
                    pattern_row->cell_global ||
                cell_bridge.global_row(
                    pattern_row->cell) !=
                    cell_global_row) {
                throw std::invalid_argument(
                    "complete assembly residual stable cell does not match partition/PETSc cell row");
            }

            const std::size_t slot =
                entry.equation_slot;
            if ((slot < n &&
                 entry.equation_kind !=
                     NaturalVariableEquationKind3D::
                         component_conservation) ||
                (slot == n &&
                 entry.equation_kind !=
                     NaturalVariableEquationKind3D::
                         energy_conservation) ||
                (slot > n &&
                 entry.equation_kind !=
                     NaturalVariableEquationKind3D::
                         fugacity_equilibrium)) {
                throw std::invalid_argument(
                    "complete assembly equation kind does not match row slot");
            }
        }

        const std::size_t component_jacobian_count =
            component.jacobian_entries().size();
        const std::size_t energy_jacobian_count =
            energy.jacobian_entries().size();
        const std::size_t fugacity_jacobian_count =
            fugacity.jacobian_entries().size();
        if (component_jacobian_count >
                std::numeric_limits<std::size_t>::max() -
                    energy_jacobian_count) {
            throw std::length_error(
                "complete assembly Jacobian count overflow");
        }
        const std::size_t component_energy_count =
            component_jacobian_count +
            energy_jacobian_count;
        if (component_energy_count >
                std::numeric_limits<std::size_t>::max() -
                    fugacity_jacobian_count) {
            throw std::length_error(
                "complete assembly Jacobian count overflow");
        }
        const std::size_t total_jacobian_count =
            component_energy_count +
            fugacity_jacobian_count;
        jacobian_entries.reserve(
            total_jacobian_count);

        for (const auto& source :
             component.jacobian_entries()) {
            jacobian_entries.push_back(
                CompleteNaturalVariableJacobianEntry3D{
                    source.petsc_global_row,
                    source.petsc_global_column,
                    source.mesh_global_row_dof,
                    source.mesh_global_column_dof,
                    source.row_cell_global,
                    source.column_cell_global,
                    source.component_row,
                    source.natural_variable_column,
                    NaturalVariableEquationKind3D::
                        component_conservation,
                    source.value});
        }

        for (const auto& source :
             energy.jacobian_entries()) {
            jacobian_entries.push_back(
                CompleteNaturalVariableJacobianEntry3D{
                    source.petsc_global_row,
                    source.petsc_global_column,
                    source.mesh_global_row_dof,
                    source.mesh_global_column_dof,
                    source.row_cell_global,
                    source.column_cell_global,
                    n,
                    source.natural_variable_column,
                    NaturalVariableEquationKind3D::
                        energy_conservation,
                    source.value});
        }

        for (const auto& source :
             fugacity.jacobian_entries()) {
            const std::size_t slot =
                fugacity_equation_slot(
                    source.non_reference_phase,
                    source.component,
                    n);
            jacobian_entries.push_back(
                CompleteNaturalVariableJacobianEntry3D{
                    source.petsc_global_row,
                    source.petsc_global_column,
                    source.mesh_global_row_dof,
                    source.mesh_global_column_dof,
                    source.row_cell_global,
                    source.row_cell_global,
                    slot,
                    source.natural_variable_column,
                    NaturalVariableEquationKind3D::
                        fugacity_equilibrium,
                    source.value});
        }

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
                    "complete assembly contains duplicate Jacobian row/column");
            }
        }

        std::size_t jacobian_index = 0U;
        for (const auto& residual :
             residual_entries) {
            const PetscInt cell_global_row =
                residual.petsc_global_row /
                static_cast<PetscInt>(q);
            const auto* pattern_row =
                find_pattern_row(
                    cell_pattern,
                    cell_global_row);
            if (pattern_row == nullptr) {
                throw std::invalid_argument(
                    "complete assembly residual row has no symbolic cell pattern");
            }

            std::vector<PetscInt> actual_columns;
            while (jacobian_index <
                       jacobian_entries.size() &&
                   jacobian_entries[jacobian_index]
                           .petsc_global_row ==
                       residual.petsc_global_row) {
                const auto& entry =
                    jacobian_entries[
                        jacobian_index];

                if (entry.mesh_global_row_dof !=
                        residual.mesh_global_row_dof ||
                    entry.row_cell_global !=
                        residual.row_cell_global ||
                    entry.equation_slot !=
                        residual.equation_slot ||
                    entry.equation_kind !=
                        residual.equation_kind ||
                    entry.petsc_global_column < 0 ||
                    entry.petsc_global_column >=
                        expected_scalar_count ||
                    entry.natural_variable_column !=
                        static_cast<std::size_t>(
                            entry.petsc_global_column %
                            static_cast<PetscInt>(q))) {
                    throw std::invalid_argument(
                        "complete assembly Jacobian row/column provenance mismatch");
                }

                if (!partition.contains_global(
                        mpmc::mesh::EntityKind::cell,
                        entry.column_cell_global)) {
                    throw std::invalid_argument(
                        "complete assembly Jacobian column cell is absent from local overlap");
                }
                const auto local_column_cell =
                    partition.local_index(
                        mpmc::mesh::EntityKind::cell,
                        entry.column_cell_global);
                const PetscInt column_cell_global_row =
                    entry.petsc_global_column /
                    static_cast<PetscInt>(q);
                if (cell_bridge.global_row(
                        local_column_cell) !=
                    column_cell_global_row) {
                    throw std::invalid_argument(
                        "complete assembly Jacobian stable column cell does not match PETSc scalar column");
                }

                actual_columns.push_back(
                    entry.petsc_global_column);
                ++jacobian_index;
            }

            const auto expected =
                expected_columns(
                    residual.equation_kind,
                    cell_pattern,
                    pattern_row->cell,
                    pattern_row->global_row,
                    q);
            if (actual_columns != expected) {
                throw std::invalid_argument(
                    "complete assembly Jacobian scalar columns do not match equation-family structural contract");
            }
        }

        if (jacobian_index !=
            jacobian_entries.size()) {
            throw std::invalid_argument(
                "complete assembly Jacobian contains a row without residual ownership");
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

    try {
        output->emplace(
            partition.local_rank(),
            partition.rank_count(),
            std::string{
                component.natural_variable_id()},
            component.component_count(),
            component.natural_variable_count(),
            component.petsc_scalar_row_start(),
            component.petsc_scalar_row_end(),
            component.petsc_scalar_row_count(),
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

#endif // MPMC_FLOW_DISCRETIZATION_PETSC_COMPLETE_NATURAL_VARIABLE_ASSEMBLY_SNAPSHOT_HPP
