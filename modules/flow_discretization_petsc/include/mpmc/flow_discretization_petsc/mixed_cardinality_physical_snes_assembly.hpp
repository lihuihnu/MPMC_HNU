#ifndef MPMC_FLOW_DISCRETIZATION_PETSC_MIXED_CARDINALITY_PHYSICAL_SNES_ASSEMBLY_HPP
#define MPMC_FLOW_DISCRETIZATION_PETSC_MIXED_CARDINALITY_PHYSICAL_SNES_ASSEMBLY_HPP

#include <mpmc/flow/cross_cardinality_phase_identity.hpp>
#include <mpmc/flow_discretization_petsc/fixed_three_phase_snes_assembly.hpp>
#include <mpmc/flow_discretization_petsc/single_phase_snes_assembly.hpp>
#include <mpmc/flow_discretization_petsc/two_phase_snes_assembly.hpp>
#include <mpmc/flow_discretization_petsc/variable_cardinality_snes_assembly.hpp>

#include <petscmat.h>
#include <petscsnes.h>
#include <petscvec.h>

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
#include <variant>
#include <vector>

namespace mpmc::flow_discretization_petsc {

inline constexpr std::string_view
    mixed_cardinality_physical_snes_assembly_convention =
        "flow_discretization_petsc/mixed-cardinality-physical-snes-assembly/v1";

using MixedCardinalityPhysicalSnesCellInput3D =
    std::variant<
        SinglePhaseSnesCellInput3D,
        TwoPhaseSnesCellInput3D,
        FixedThreePhaseSnesCellInput3D>;

using MixedCardinalityPhysicalCurrentCellLinearization3D =
    std::variant<
        SinglePhaseCurrentCellLinearization3D,
        TwoPhaseCurrentCellLinearization3D,
        FixedThreePhaseCurrentCellLinearization3D>;

using MixedCardinalityPhysicalMobilityLinearization3D =
    std::variant<
        mpmc::flow::SinglePhaseMobilityLinearization,
        mpmc::flow::TwoPhaseMobilityLinearization,
        mpmc::flow::LocalPhaseMobilityLinearization3P>;

using MixedCardinalityPhysicalSnesAuthoritativeFaceInput3D =
    FixedThreePhaseSnesAuthoritativeFaceInput3D;

struct MixedCardinalityPhysicalCellEvaluatorBindings3D {
    SinglePhaseCurrentCellEvaluatorBinding3D single_phase;
    TwoPhaseCurrentCellEvaluatorBinding3D two_phase;
    FixedThreePhaseCurrentCellEvaluatorBinding3D three_phase;
};

struct MixedCardinalityPhysicalFaceLinearization3D {
    mpmc::flow_discretization::
        NormalizedComponentFaceContributionLinearization3D
            component;
    mpmc::flow_discretization::
        NormalizedEnergyFaceContributionLinearization3D
            energy;
};

using MixedCardinalityCrossPhaseFaceEvaluator3D =
    PetscErrorCode (*)(
        const MixedCardinalityPhysicalSnesAuthoritativeFaceInput3D&
            face_input,
        const mpmc::flow::
            CrossCardinalityFacePhaseIdentityPlan&
                phase_identity_plan,
        mpmc::flow_discretization::
            TwoCellBulkVolume3D bulk_volume,
        const MixedCardinalityPhysicalCurrentCellLinearization3D&
            owner,
        const MixedCardinalityPhysicalCurrentCellLinearization3D&
            neighbour,
        void* user_context,
        std::optional<
            MixedCardinalityPhysicalFaceLinearization3D>*
                output,
        NaturalVariableSnesEvaluationStatus3D*
            status);

struct MixedCardinalityCrossPhaseFaceEvaluatorBinding3D {
    MixedCardinalityCrossPhaseFaceEvaluator3D evaluator{};
    void* user_context{};
};

namespace mixed_cardinality_physical_detail {

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

[[nodiscard]] inline std::size_t
phase_count(
    const MixedCardinalityPhysicalSnesCellInput3D&
        input) {
    return std::visit(
        [](const auto& typed) {
            return typed.frozen_layout.phase_count();
        },
        input);
}

[[nodiscard]] inline mpmc::mesh::LocalIndex
cell(
    const MixedCardinalityPhysicalSnesCellInput3D&
        input) {
    return std::visit(
        [](const auto& typed) {
            return typed.cell;
        },
        input);
}

[[nodiscard]] inline mpmc::mesh::GlobalEntityId
cell_global(
    const MixedCardinalityPhysicalSnesCellInput3D&
        input) {
    return std::visit(
        [](const auto& typed) {
            return typed.cell_global;
        },
        input);
}

[[nodiscard]] inline double
bulk_volume(
    const MixedCardinalityPhysicalSnesCellInput3D&
        input) {
    return std::visit(
        [](const auto& typed) {
            return typed.bulk_volume_m3;
        },
        input);
}

[[nodiscard]] inline double
porosity(
    const MixedCardinalityPhysicalSnesCellInput3D&
        input) {
    return std::visit(
        [](const auto& typed) {
            return typed.porosity;
        },
        input);
}

[[nodiscard]] inline const std::vector<std::string>&
component_ids(
    const MixedCardinalityPhysicalSnesCellInput3D&
        input) {
    return std::visit(
        [](const auto& typed)
            -> const std::vector<std::string>& {
            return typed.component_ids;
        },
        input);
}

[[nodiscard]] inline const std::optional<
    mpmc::flow::PoreVolumeComponentAccumulationSnapshot3P>&
previous_component(
    const MixedCardinalityPhysicalSnesCellInput3D&
        input) {
    return std::visit(
        [](const auto& typed)
            -> const std::optional<
                mpmc::flow::
                    PoreVolumeComponentAccumulationSnapshot3P>& {
            return typed.previous_component_accumulation;
        },
        input);
}

[[nodiscard]] inline const std::optional<
    mpmc::flow::PoreVolumeEnergyAccumulationSnapshot3P>&
previous_energy(
    const MixedCardinalityPhysicalSnesCellInput3D&
        input) {
    return std::visit(
        [](const auto& typed)
            -> const std::optional<
                mpmc::flow::
                    PoreVolumeEnergyAccumulationSnapshot3P>& {
            return typed.previous_energy_accumulation;
        },
        input);
}

[[nodiscard]] inline const
mpmc::flow::NaturalVariableStateIdentity3P&
state_identity(
    const MixedCardinalityPhysicalCurrentCellLinearization3D&
        current) {
    return std::visit(
        [](const auto& typed)
            -> const mpmc::flow::
                NaturalVariableStateIdentity3P& {
            return typed.transport.state_identity;
        },
        current);
}

[[nodiscard]] inline bool
near_roundoff(
    double first,
    double second) {
    if (!std::isfinite(first) ||
        !std::isfinite(second)) {
        return false;
    }
    const double scale =
        std::max(
            {1.0,
             std::abs(first),
             std::abs(second)});
    return std::abs(first - second) <=
        8192.0 *
            std::numeric_limits<double>::
                epsilon() *
            scale;
}

[[nodiscard]] inline PetscErrorCode
unused_cell_evaluator(
    mpmc::mesh::LocalIndex,
    std::size_t,
    std::span<const double>,
    const VariableCardinalityNaturalVariableNumbering3D&,
    void*,
    VariableCardinalityNaturalVariableCellAssembly3D*,
    NaturalVariableSnesEvaluationStatus3D*) {
    return PETSC_ERR_SUP;
}

inline void add_dense_block(
    std::vector<double>* target,
    std::span<const double> source) {
    if (target == nullptr ||
        target->size() != source.size()) {
        throw std::invalid_argument(
            "mixed-cardinality physical dense block shape mismatch");
    }
    for (std::size_t index = 0U;
         index < target->size();
         ++index) {
        (*target)[index] += source[index];
        if (!std::isfinite(
                (*target)[index])) {
            throw std::range_error(
                "mixed-cardinality physical dense block became non-finite");
        }
    }
}

[[nodiscard]] inline PetscErrorCode
insert_dense_cell_assembly(
    const VariableCardinalityNaturalVariableCellDof3D&
        cell_record,
    const VariableCardinalityNaturalVariableCellAssembly3D&
        assembly,
    Vec residual,
    Mat jacobian,
    bool insert_residual,
    bool insert_jacobian) {
    const std::size_t q =
        cell_record.scalar_count;
    if (assembly.residual.size() != q ||
        assembly.jacobian_blocks.size() != 1U ||
        assembly.jacobian_blocks.front().column_cell !=
            cell_record.cell ||
        assembly.jacobian_blocks.front()
                .values_row_major.size() !=
            q * q) {
        return PETSC_ERR_ARG_SIZ;
    }

    std::vector<PetscInt> rows(q);
    for (std::size_t slot = 0U;
         slot < q;
         ++slot) {
        rows[slot] =
            cell_record.petsc_global_scalar_start +
            static_cast<PetscInt>(slot);
    }

    if (insert_residual) {
        std::vector<PetscScalar> values(q);
        for (std::size_t row = 0U;
             row < q;
             ++row) {
            values[row] =
                static_cast<PetscScalar>(
                    assembly.residual[row]);
        }
        const PetscErrorCode error =
            VecSetValues(
                residual,
                static_cast<PetscInt>(q),
                rows.data(),
                values.data(),
                ADD_VALUES);
        if (error != PETSC_SUCCESS) {
            return error;
        }
    }

    if (insert_jacobian) {
        std::vector<PetscScalar> values(
            q * q);
        for (std::size_t index = 0U;
             index < values.size();
             ++index) {
            values[index] =
                static_cast<PetscScalar>(
                    assembly.jacobian_blocks.front()
                        .values_row_major[index]);
        }
        return MatSetValues(
            jacobian,
            static_cast<PetscInt>(q),
            rows.data(),
            static_cast<PetscInt>(q),
            rows.data(),
            values.data(),
            ADD_VALUES);
    }
    return PETSC_SUCCESS;
}

[[nodiscard]] inline PetscErrorCode
insert_face_side(
    const VariableCardinalityNaturalVariableCellDof3D&
        row_cell,
    const VariableCardinalityNaturalVariableCellDof3D&
        other_cell,
    std::size_t component_count,
    std::span<const double> component_residual,
    double energy_residual,
    std::span<const double> component_diagonal,
    std::span<const double> component_off_diagonal,
    std::span<const double> energy_diagonal,
    std::span<const double> energy_off_diagonal,
    Vec residual,
    Mat jacobian,
    bool insert_residual,
    bool insert_jacobian) {
    const std::size_t conservation_rows =
        component_count + 1U;
    if (component_residual.size() !=
            component_count ||
        component_diagonal.size() !=
            component_count *
                row_cell.scalar_count ||
        component_off_diagonal.size() !=
            component_count *
                other_cell.scalar_count ||
        energy_diagonal.size() !=
            row_cell.scalar_count ||
        energy_off_diagonal.size() !=
            other_cell.scalar_count) {
        return PETSC_ERR_ARG_SIZ;
    }

    std::vector<PetscInt>
        rows(conservation_rows);
    for (std::size_t row = 0U;
         row < conservation_rows;
         ++row) {
        rows[row] =
            row_cell.petsc_global_scalar_start +
            static_cast<PetscInt>(row);
    }

    if (insert_residual) {
        std::vector<PetscScalar>
            values(conservation_rows);
        for (std::size_t component = 0U;
             component < component_count;
             ++component) {
            values[component] =
                static_cast<PetscScalar>(
                    component_residual[component]);
        }
        values[component_count] =
            static_cast<PetscScalar>(
                energy_residual);
        const PetscErrorCode error =
            VecSetValues(
                residual,
                static_cast<PetscInt>(
                    rows.size()),
                rows.data(),
                values.data(),
                ADD_VALUES);
        if (error != PETSC_SUCCESS) {
            return error;
        }
    }

    if (!insert_jacobian) {
        return PETSC_SUCCESS;
    }

    const auto make_block =
        [&](std::size_t column_count,
            std::span<const double> component,
            std::span<const double> energy) {
            std::vector<PetscScalar>
                values(
                    conservation_rows *
                    column_count,
                    PetscScalar{0.0});
            for (std::size_t row = 0U;
                 row < component_count;
                 ++row) {
                for (std::size_t column = 0U;
                     column < column_count;
                     ++column) {
                    values[
                        row * column_count +
                        column] =
                        static_cast<PetscScalar>(
                            component[
                                row * column_count +
                                column]);
                }
            }
            for (std::size_t column = 0U;
                 column < column_count;
                 ++column) {
                values[
                    component_count *
                        column_count +
                    column] =
                    static_cast<PetscScalar>(
                        energy[column]);
            }
            return values;
        };

    std::vector<PetscInt>
        diagonal_columns(
            row_cell.scalar_count);
    for (std::size_t slot = 0U;
         slot < row_cell.scalar_count;
         ++slot) {
        diagonal_columns[slot] =
            row_cell.petsc_global_scalar_start +
            static_cast<PetscInt>(slot);
    }
    auto diagonal_values =
        make_block(
            row_cell.scalar_count,
            component_diagonal,
            energy_diagonal);
    PetscErrorCode error =
        MatSetValues(
            jacobian,
            static_cast<PetscInt>(
                rows.size()),
            rows.data(),
            static_cast<PetscInt>(
                diagonal_columns.size()),
            diagonal_columns.data(),
            diagonal_values.data(),
            ADD_VALUES);
    if (error != PETSC_SUCCESS) {
        return error;
    }

    std::vector<PetscInt>
        off_columns(
            other_cell.scalar_count);
    for (std::size_t slot = 0U;
         slot < other_cell.scalar_count;
         ++slot) {
        off_columns[slot] =
            other_cell.petsc_global_scalar_start +
            static_cast<PetscInt>(slot);
    }
    auto off_values =
        make_block(
            other_cell.scalar_count,
            component_off_diagonal,
            energy_off_diagonal);
    return MatSetValues(
        jacobian,
        static_cast<PetscInt>(
            rows.size()),
        rows.data(),
        static_cast<PetscInt>(
            off_columns.size()),
        off_columns.data(),
        off_values.data(),
        ADD_VALUES);
}

} // namespace mixed_cardinality_physical_detail

class MixedCardinalityPhysicalSnesAssemblyContext3D {
public:
    static constexpr std::string_view convention =
        mixed_cardinality_physical_snes_assembly_convention;

    MixedCardinalityPhysicalSnesAssemblyContext3D(
        const MixedCardinalityPhysicalSnesAssemblyContext3D&) =
        delete;
    MixedCardinalityPhysicalSnesAssemblyContext3D& operator=(
        const MixedCardinalityPhysicalSnesAssemblyContext3D&) =
        delete;
    MixedCardinalityPhysicalSnesAssemblyContext3D& operator=(
        MixedCardinalityPhysicalSnesAssemblyContext3D&&) =
        delete;
    MixedCardinalityPhysicalSnesAssemblyContext3D(
        MixedCardinalityPhysicalSnesAssemblyContext3D&&) noexcept =
        default;
    ~MixedCardinalityPhysicalSnesAssemblyContext3D() =
        default;

    [[nodiscard]] static PetscErrorCode
    create(
        MPI_Comm comm,
        const mpmc::discretization_petsc::
            ParallelOwnedConnectionSchedule3D&
                schedule,
        const mpmc::mesh::PartitionSnapshot&
            partition,
        const VariableCardinalityNaturalVariableNumbering3D&
            numbering,
        const mpmc::discretization_petsc::
            PetscMpiAijSymbolicPreallocation3D&
                cell_bridge,
        const mpmc::discretization_petsc::
            OwnedCellStructuralColumnPatternSnapshot3D&
                cell_pattern,
        double time_step_seconds,
        std::vector<
            MixedCardinalityPhysicalSnesCellInput3D>
            cell_inputs,
        std::vector<
            mpmc::flow::FrozenActivePhaseIdentityMap>
            phase_identity_maps,
        std::vector<
            MixedCardinalityPhysicalSnesAuthoritativeFaceInput3D>
            face_inputs,
        MixedCardinalityPhysicalCellEvaluatorBindings3D
            cell_evaluators,
        MixedCardinalityCrossPhaseFaceEvaluatorBinding3D
            cross_phase_face_evaluator,
        std::optional<
            MixedCardinalityPhysicalSnesAssemblyContext3D>*
            output) {
        using namespace
            mixed_cardinality_physical_detail;

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
        std::vector<
            mpmc::flow::
                CrossCardinalityFacePhaseIdentityPlan>
            face_phase_identity_plans;
        try {
            if (!std::isfinite(
                    time_step_seconds) ||
                !(time_step_seconds > 0.0) ||
                schedule.local_rank() !=
                    partition.local_rank() ||
                schedule.rank_count() !=
                    partition.rank_count() ||
                numbering.local_rank() !=
                    partition.local_rank() ||
                numbering.rank_count() !=
                    partition.rank_count() ||
                cell_bridge.local_rank() !=
                    partition.local_rank() ||
                cell_bridge.rank_count() !=
                    partition.rank_count() ||
                cell_pattern.local_rank() !=
                    partition.local_rank() ||
                cell_pattern.rank_count() !=
                    partition.rank_count() ||
                partition.local_rank().value() !=
                    static_cast<std::uint32_t>(
                        mpi_rank) ||
                partition.rank_count() !=
                    static_cast<std::uint32_t>(
                        mpi_size) ||
                cell_inputs.size() !=
                    numbering.local_cell_count() ||
                cell_inputs.size() !=
                    partition.entity_count(
                        mpmc::mesh::EntityKind::cell) ||
                phase_identity_maps.size() !=
                    cell_inputs.size() ||
                face_inputs.size() !=
                    schedule.assembly_rows().size() ||
                cell_evaluators.single_phase.evaluator ==
                    nullptr ||
                cell_evaluators.two_phase.evaluator ==
                    nullptr ||
                cell_evaluators.three_phase.evaluator ==
                    nullptr) {
                throw std::invalid_argument(
                    "mixed-cardinality physical SNES metadata mismatch");
            }

            std::vector<std::string>
                canonical_component_ids;
            for (std::size_t local = 0U;
                 local < cell_inputs.size();
                 ++local) {
                const auto cell_index =
                    mpmc::mesh::LocalIndex{
                        static_cast<
                            mpmc::mesh::LocalIndex::value_type>(
                                local)};
                const auto& input =
                    cell_inputs[local];
                const auto& record =
                    numbering.cell(
                        cell_index);
                if (cell(input) !=
                        cell_index ||
                    cell_global(input) !=
                        record.cell_global ||
                    partition.global_id(
                        mpmc::mesh::EntityKind::cell,
                        cell_index) !=
                        record.cell_global ||
                    partition.owner_rank(
                        mpmc::mesh::EntityKind::cell,
                        cell_index) !=
                        record.owner_rank ||
                    phase_count(input) !=
                        record.phase_count ||
                    phase_identity_maps[local]
                            .phase_count() !=
                        record.phase_count ||
                    component_ids(input).size() !=
                        numbering.component_count() ||
                    record.scalar_count !=
                        record.phase_count *
                            numbering.component_count() +
                        1U ||
                    !std::isfinite(
                        bulk_volume(input)) ||
                    !(bulk_volume(input) > 0.0) ||
                    !std::isfinite(
                        porosity(input)) ||
                    !(porosity(input) > 0.0) ||
                    !(porosity(input) < 1.0)) {
                    throw std::invalid_argument(
                        "invalid mixed-cardinality physical cell binding");
                }

                const std::size_t layout_q =
                    std::visit(
                        [](const auto& typed) {
                            return typed.frozen_layout
                                .unknown_count();
                        },
                        input);
                if (layout_q !=
                    record.scalar_count) {
                    throw std::invalid_argument(
                        "mixed-cardinality frozen layout disagrees with numbering");
                }

                if (local == 0U) {
                    canonical_component_ids =
                        component_ids(input);
                } else if (
                    component_ids(input) !=
                    canonical_component_ids) {
                    throw std::invalid_argument(
                        "component identity/order is not canonical across mixed-cardinality overlap");
                }

                const bool owned =
                    record.owner_rank ==
                    numbering.local_rank();
                if (owned !=
                        previous_component(input)
                            .has_value() ||
                    owned !=
                        previous_energy(input)
                            .has_value()) {
                    throw std::invalid_argument(
                        "owned mixed-cardinality cell history is required and ghost history is forbidden");
                }
            }

            face_phase_identity_plans.reserve(
                face_inputs.size());

            for (std::size_t index = 0U;
                 index < face_inputs.size();
                 ++index) {
                const auto& row =
                    schedule.assembly_rows()[index];
                const auto& input =
                    face_inputs[index];
                if (input.face !=
                        row.face ||
                    input.face_global !=
                        row.face_global ||
                    input.transmissibility.face !=
                        row.face ||
                    !input.transmissibility
                         .static_transmissibility
                         .has_value() ||
                    !near_roundoff(
                        input.transmissibility
                            .static_transmissibility
                            ->face_transmissibility_m3,
                        row.transmissibility_m3) ||
                    !std::isfinite(
                        input.thermal_conductance
                            .conductance_w_per_k) ||
                    input.thermal_conductance
                            .conductance_w_per_k <
                        0.0) {
                    throw std::invalid_argument(
                        "invalid mixed-cardinality authoritative face binding");
                }

                const auto owner_phase_count =
                    numbering.cell(
                        row.owner_cell)
                        .phase_count;
                const auto neighbour_phase_count =
                    numbering.cell(
                        row.neighbour_cell)
                        .phase_count;
                auto phase_plan =
                    mpmc::flow::
                        make_cross_cardinality_face_phase_identity_plan(
                            phase_identity_maps.at(
                                static_cast<std::size_t>(
                                    row.owner_cell.value())),
                            phase_identity_maps.at(
                                static_cast<std::size_t>(
                                    row.neighbour_cell.value())));
                const bool direct_fixed_cardinality_tpfa =
                    owner_phase_count ==
                        neighbour_phase_count &&
                    phase_plan
                        .slot_aligned_same_active_set();
                if (!direct_fixed_cardinality_tpfa &&
                    cross_phase_face_evaluator
                            .evaluator ==
                        nullptr) {
                    local_error =
                        PETSC_ERR_SUP;
                    break;
                }
                face_phase_identity_plans.push_back(
                    std::move(phase_plan));
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

        std::optional<
            VariableCardinalityNaturalVariableSnesAssemblyContext3D>
            infrastructure;
        error =
            VariableCardinalityNaturalVariableSnesAssemblyContext3D::
                create(
                    comm,
                    numbering,
                    cell_bridge,
                    cell_pattern,
                    {
                        &unused_cell_evaluator,
                        nullptr},
                    &infrastructure);
        if (error != PETSC_SUCCESS ||
            !infrastructure.has_value()) {
            return error != PETSC_SUCCESS
                ? error
                : PETSC_ERR_PLIB;
        }

        MixedCardinalityPhysicalSnesAssemblyContext3D
            context{
                comm,
                schedule,
                partition,
                numbering,
                time_step_seconds,
                std::move(cell_inputs),
                std::move(phase_identity_maps),
                std::move(face_phase_identity_plans),
                std::move(face_inputs),
                cell_evaluators,
                cross_phase_face_evaluator,
                std::move(*infrastructure)};
        output->emplace(
            std::move(context));
        return PETSC_SUCCESS;
    }

    [[nodiscard]] NaturalVariableSnesEvaluator3D
    snes_evaluator() noexcept {
        return {
            &snes_function,
            &snes_jacobian,
            &snes_precheck,
            this};
    }

    [[nodiscard]] PetscErrorCode
    create_jacobian_structure(
        Mat* output) const {
        return infrastructure_
            .create_jacobian_structure(
                output);
    }

    [[nodiscard]] PetscErrorCode
    evaluate_local_cells_for_phase_transition(
        Vec global_state,
        std::vector<std::optional<
            MixedCardinalityPhysicalCurrentCellLinearization3D>>*
                output,
        std::vector<double>* porosities,
        NaturalVariableSnesEvaluationStatus3D*
            status) {
        using namespace
            mixed_cardinality_physical_detail;

        if (global_state == nullptr ||
            output == nullptr ||
            porosities == nullptr ||
            status == nullptr) {
            return PETSC_ERR_ARG_NULL;
        }
        output->clear();
        porosities->clear();
        *status =
            NaturalVariableSnesEvaluationStatus3D::
                success;

        std::vector<double> local_state;
        PetscErrorCode error =
            infrastructure_
                .copy_local_packed_state(
                    global_state,
                    &local_state);
        error =
            collective_error(
                comm_,
                error);
        if (error != PETSC_SUCCESS) {
            return error;
        }

        std::vector<std::optional<
            MixedCardinalityPhysicalMobilityLinearization3D>>
            unused_mobility;
        error =
            evaluate_cells(
                local_state,
                output,
                &unused_mobility,
                status);
        if (error != PETSC_SUCCESS ||
            *status !=
                NaturalVariableSnesEvaluationStatus3D::
                    success) {
            return error;
        }

        porosities->reserve(
            cell_inputs_.size());
        for (const auto& input :
             cell_inputs_) {
            porosities->push_back(
                porosity(input));
        }
        return PETSC_SUCCESS;
    }

    [[nodiscard]] PetscErrorCode
    evaluate_phase_transition_cell(
        mpmc::mesh::LocalIndex cell,
        std::span<const double> natural_variables,
        std::optional<
            MixedCardinalityPhysicalCurrentCellLinearization3D>*
                output,
        double* porosity,
        NaturalVariableSnesEvaluationStatus3D*
            status) {
        if (output == nullptr ||
            porosity == nullptr ||
            status == nullptr) {
            return PETSC_ERR_ARG_NULL;
        }
        const std::size_t local =
            static_cast<std::size_t>(
                cell.value());
        if (local >=
            numbering_->local_cell_count()) {
            return PETSC_ERR_ARG_OUTOFRANGE;
        }
        *porosity =
            mixed_cardinality_physical_detail::
                porosity(
                    cell_inputs_[local]);
        std::optional<
            MixedCardinalityPhysicalMobilityLinearization3D>
            unused_mobility;
        return evaluate_cell(
            local,
            natural_variables,
            output,
            &unused_mobility,
            status);
    }

private:
    MixedCardinalityPhysicalSnesAssemblyContext3D(
        MPI_Comm comm,
        const mpmc::discretization_petsc::
            ParallelOwnedConnectionSchedule3D&
                schedule,
        const mpmc::mesh::PartitionSnapshot&
            partition,
        const VariableCardinalityNaturalVariableNumbering3D&
            numbering,
        double time_step_seconds,
        std::vector<
            MixedCardinalityPhysicalSnesCellInput3D>
            cell_inputs,
        std::vector<
            mpmc::flow::FrozenActivePhaseIdentityMap>
            phase_identity_maps,
        std::vector<
            mpmc::flow::
                CrossCardinalityFacePhaseIdentityPlan>
            face_phase_identity_plans,
        std::vector<
            MixedCardinalityPhysicalSnesAuthoritativeFaceInput3D>
            face_inputs,
        MixedCardinalityPhysicalCellEvaluatorBindings3D
            cell_evaluators,
        MixedCardinalityCrossPhaseFaceEvaluatorBinding3D
            cross_phase_face_evaluator,
        VariableCardinalityNaturalVariableSnesAssemblyContext3D
            infrastructure)
        : comm_(comm),
          schedule_(&schedule),
          partition_(&partition),
          numbering_(&numbering),
          time_step_seconds_(
              time_step_seconds),
          cell_inputs_(
              std::move(cell_inputs)),
          phase_identity_maps_(
              std::move(phase_identity_maps)),
          face_phase_identity_plans_(
              std::move(
                  face_phase_identity_plans)),
          face_inputs_(
              std::move(face_inputs)),
          cell_evaluators_(
              cell_evaluators),
          cross_phase_face_evaluator_(
              cross_phase_face_evaluator),
          infrastructure_(
              std::move(infrastructure)) {}

    [[nodiscard]] PetscErrorCode
    evaluate_cell(
        std::size_t local,
        std::span<const double> values,
        std::optional<
            MixedCardinalityPhysicalCurrentCellLinearization3D>*
                current,
        std::optional<
            MixedCardinalityPhysicalMobilityLinearization3D>*
                mobility,
        NaturalVariableSnesEvaluationStatus3D*
            status) {
        using namespace
            mixed_cardinality_physical_detail;

        if (current == nullptr ||
            mobility == nullptr ||
            status == nullptr) {
            return PETSC_ERR_ARG_NULL;
        }
        current->reset();
        mobility->reset();
        *status =
            NaturalVariableSnesEvaluationStatus3D::
                success;

        if (local >=
            numbering_->local_cell_count()) {
            return PETSC_ERR_ARG_OUTOFRANGE;
        }
        const auto cell_index =
            mpmc::mesh::LocalIndex{
                static_cast<
                    mpmc::mesh::LocalIndex::value_type>(
                        local)};
        const auto& record =
            numbering_->cell(
                cell_index);
        const auto& input =
            cell_inputs_[local];
        if (values.size() !=
            record.scalar_count) {
            return PETSC_ERR_ARG_SIZ;
        }

        NaturalVariableSnesEvaluationStatus3D
            cell_status =
                NaturalVariableSnesEvaluationStatus3D::
                    success;
        PetscErrorCode error =
            PETSC_SUCCESS;

        if (record.phase_count == 1U) {
            const auto* typed =
                std::get_if<
                    SinglePhaseSnesCellInput3D>(
                        &input);
            if (typed == nullptr) {
                return PETSC_ERR_ARG_INCOMP;
            }
            std::optional<
                SinglePhaseCurrentCellLinearization3D>
                evaluated;
            error =
                cell_evaluators_
                    .single_phase.evaluator(
                        typed->cell,
                        typed->cell_global,
                        values,
                        typed->frozen_layout,
                        typed->component_ids,
                        cell_evaluators_
                            .single_phase
                            .user_context,
                        &evaluated,
                        &cell_status);
            if (error == PETSC_SUCCESS &&
                cell_status ==
                    NaturalVariableSnesEvaluationStatus3D::
                        success &&
                evaluated.has_value()) {
                mobility->emplace(
                    mpmc::flow::
                        build_single_phase_mobility_linearization(
                            evaluated->state,
                            evaluated->transport));
                current->emplace(
                    std::move(*evaluated));
            }
        } else if (record.phase_count == 2U) {
            const auto* typed =
                std::get_if<
                    TwoPhaseSnesCellInput3D>(
                        &input);
            if (typed == nullptr) {
                return PETSC_ERR_ARG_INCOMP;
            }
            std::optional<
                TwoPhaseCurrentCellLinearization3D>
                evaluated;
            error =
                cell_evaluators_
                    .two_phase.evaluator(
                        typed->cell,
                        typed->cell_global,
                        values,
                        typed->frozen_layout,
                        typed->component_ids,
                        cell_evaluators_
                            .two_phase
                            .user_context,
                        &evaluated,
                        &cell_status);
            if (error == PETSC_SUCCESS &&
                cell_status ==
                    NaturalVariableSnesEvaluationStatus3D::
                        success &&
                evaluated.has_value()) {
                mobility->emplace(
                    mpmc::flow::
                        build_two_phase_mobility_linearization(
                            evaluated->state,
                            evaluated->transport));
                current->emplace(
                    std::move(*evaluated));
            }
        } else if (record.phase_count == 3U) {
            const auto* typed =
                std::get_if<
                    FixedThreePhaseSnesCellInput3D>(
                        &input);
            if (typed == nullptr) {
                return PETSC_ERR_ARG_INCOMP;
            }
            std::optional<
                FixedThreePhaseCurrentCellLinearization3D>
                evaluated;
            error =
                cell_evaluators_
                    .three_phase.evaluator(
                        typed->cell,
                        typed->cell_global,
                        values,
                        typed->frozen_layout,
                        typed->component_ids,
                        cell_evaluators_
                            .three_phase
                            .user_context,
                        &evaluated,
                        &cell_status);
            if (error == PETSC_SUCCESS &&
                cell_status ==
                    NaturalVariableSnesEvaluationStatus3D::
                        success &&
                evaluated.has_value()) {
                mobility->emplace(
                    mpmc::flow::
                        build_local_phase_mobility_linearization(
                            evaluated->state,
                            evaluated->transport,
                            evaluated
                                ->saturation_constitutive));
                current->emplace(
                    std::move(*evaluated));
            }
        } else {
            return PETSC_ERR_ARG_INCOMP;
        }

        if (error != PETSC_SUCCESS) {
            return error;
        }
        if (cell_status ==
            NaturalVariableSnesEvaluationStatus3D::
                domain_error) {
            *status =
                NaturalVariableSnesEvaluationStatus3D::
                    domain_error;
            return PETSC_SUCCESS;
        }
        if (cell_status !=
                NaturalVariableSnesEvaluationStatus3D::
                    success ||
            !current->has_value() ||
            !mobility->has_value()) {
            return PETSC_ERR_ARG_INCOMP;
        }

        const auto& identity =
            state_identity(
                **current);
        if (identity.component_ids !=
                component_ids(input) ||
            identity.layout.component_count() !=
                numbering_->component_count() ||
            identity.layout.phase_count() !=
                record.phase_count ||
            identity.layout.unknown_count() !=
                record.scalar_count) {
            return PETSC_ERR_ARG_INCOMP;
        }
        return PETSC_SUCCESS;
    }

    [[nodiscard]] PetscErrorCode
    evaluate_cells(
        std::span<const double> local_state,
        std::vector<std::optional<
            MixedCardinalityPhysicalCurrentCellLinearization3D>>*
                current,
        std::vector<std::optional<
            MixedCardinalityPhysicalMobilityLinearization3D>>*
                mobility,
        NaturalVariableSnesEvaluationStatus3D*
            status) {
        if (current == nullptr ||
            mobility == nullptr ||
            status == nullptr) {
            return PETSC_ERR_ARG_NULL;
        }
        current->assign(
            numbering_->local_cell_count(),
            std::nullopt);
        mobility->assign(
            numbering_->local_cell_count(),
            std::nullopt);
        *status =
            NaturalVariableSnesEvaluationStatus3D::
                success;

        PetscErrorCode local_error =
            PETSC_SUCCESS;
        int local_domain = 0;

        for (std::size_t local = 0U;
             local < numbering_->local_cell_count();
             ++local) {
            const auto cell_index =
                mpmc::mesh::LocalIndex{
                    static_cast<
                        mpmc::mesh::LocalIndex::value_type>(
                            local)};
            const auto& record =
                numbering_->cell(
                    cell_index);
            const std::span<const double>
                values =
                    local_state.subspan(
                        record.local_scalar_offset,
                        record.scalar_count);
            NaturalVariableSnesEvaluationStatus3D
                cell_status =
                    NaturalVariableSnesEvaluationStatus3D::
                        success;
            std::optional<
                MixedCardinalityPhysicalCurrentCellLinearization3D>
                cell_current;
            std::optional<
                MixedCardinalityPhysicalMobilityLinearization3D>
                cell_mobility;
            local_error =
                evaluate_cell(
                    local,
                    values,
                    &cell_current,
                    &cell_mobility,
                    &cell_status);
            if (local_error != PETSC_SUCCESS) {
                break;
            }
            if (cell_status ==
                NaturalVariableSnesEvaluationStatus3D::
                    domain_error) {
                local_domain = 1;
                continue;
            }
            (*current)[local] =
                std::move(cell_current);
            (*mobility)[local] =
                std::move(cell_mobility);
        }

        PetscErrorCode error =
            mixed_cardinality_physical_detail::collective_error(
                comm_,
                local_error);
        if (error != PETSC_SUCCESS) {
            return error;
        }

        int global_domain = 0;
        if (MPI_Allreduce(
                &local_domain,
                &global_domain,
                1,
                MPI_INT,
                MPI_MAX,
                comm_) != MPI_SUCCESS) {
            return PETSC_ERR_MPI;
        }
        if (global_domain != 0) {
            *status =
                NaturalVariableSnesEvaluationStatus3D::
                    domain_error;
        }
        return PETSC_SUCCESS;
    }

    [[nodiscard]]
    VariableCardinalityNaturalVariableCellAssembly3D
    build_owned_cell_assembly(
        std::size_t local,
        const MixedCardinalityPhysicalCurrentCellLinearization3D&
            current) const {
        using namespace
            mixed_cardinality_physical_detail;

        const auto cell_index =
            mpmc::mesh::LocalIndex{
                static_cast<
                    mpmc::mesh::LocalIndex::value_type>(
                        local)};
        const auto& record =
            numbering_->cell(
                cell_index);
        const auto& input =
            cell_inputs_[local];
        const std::size_t q =
            record.scalar_count;
        const std::size_t nc =
            numbering_->component_count();

        VariableCardinalityNaturalVariableCellAssembly3D
            assembly;
        assembly.residual.assign(
            q,
            0.0);
        VariableCardinalityNaturalVariableJacobianBlock3D
            diagonal;
        diagonal.column_cell =
            cell_index;
        diagonal.values_row_major.assign(
            q * q,
            0.0);

        const auto add_conservation =
            [&](const mpmc::flow::
                    BackwardEulerComponentAccumulationResidual3P&
                        component,
                const mpmc::flow::
                    BackwardEulerEnergyAccumulationResidual3P&
                        energy) {
                if (component.component_count() !=
                        nc ||
                    component.input_count !=
                        q ||
                    energy.input_count !=
                        q) {
                    throw std::invalid_argument(
                        "mixed-cardinality accumulation shape mismatch");
                }
                for (std::size_t row = 0U;
                     row < nc;
                     ++row) {
                    assembly.residual[row] +=
                        component.residual(row);
                    for (std::size_t column = 0U;
                         column < q;
                         ++column) {
                        diagonal.values_row_major[
                            row * q +
                            column] +=
                            component.d_residual(
                                row,
                                column);
                    }
                }
                assembly.residual[nc] +=
                    energy.residual_w_per_bulk_m3;
                for (std::size_t column = 0U;
                     column < q;
                     ++column) {
                    diagonal.values_row_major[
                        nc * q +
                        column] +=
                        energy.d_residual(
                            column);
                }
            };

        if (record.phase_count == 1U) {
            const auto& evaluated =
                std::get<
                    SinglePhaseCurrentCellLinearization3D>(
                        current);
            const auto component_current =
                mpmc::flow::
                    build_single_phase_component_accumulation(
                        evaluated.state,
                        porosity(input));
            const auto component_linearization =
                mpmc::flow::
                    build_single_phase_component_accumulation_linearization(
                        evaluated.state,
                        porosity(input),
                        evaluated.molar_density);
            const auto pair =
                mpmc::flow::
                    make_pore_volume_component_accumulation_pair(
                        component_current,
                        *previous_component(input));
            const auto component =
                mpmc::flow::
                    build_backward_euler_component_accumulation_residual(
                        pair,
                        component_linearization,
                        time_step_seconds_);
            const auto energy_current =
                mpmc::flow::
                    build_single_phase_energy_accumulation_snapshot(
                        evaluated.state,
                        porosity(input),
                        evaluated.transport,
                        evaluated.caloric,
                        evaluated.rock);
            const auto energy_linearization =
                mpmc::flow::
                    build_single_phase_energy_accumulation_linearization(
                        evaluated.state,
                        porosity(input),
                        evaluated.transport,
                        evaluated.caloric,
                        evaluated.rock);
            const auto energy =
                mpmc::flow::
                    build_backward_euler_energy_accumulation_residual(
                        energy_current,
                        energy_linearization,
                        *previous_energy(input),
                        time_step_seconds_);
            add_conservation(
                component,
                energy);
        } else if (
            record.phase_count == 2U) {
            const auto& evaluated =
                std::get<
                    TwoPhaseCurrentCellLinearization3D>(
                        current);
            const auto component_current =
                mpmc::flow::
                    build_two_phase_component_accumulation(
                        evaluated.state,
                        porosity(input));
            const auto component_linearization =
                mpmc::flow::
                    build_two_phase_component_accumulation_linearization(
                        evaluated.state,
                        porosity(input),
                        evaluated.molar_density);
            const auto pair =
                mpmc::flow::
                    make_pore_volume_component_accumulation_pair(
                        component_current,
                        *previous_component(input));
            const auto component =
                mpmc::flow::
                    build_backward_euler_component_accumulation_residual(
                        pair,
                        component_linearization,
                        time_step_seconds_);
            const auto energy_current =
                mpmc::flow::
                    build_two_phase_energy_accumulation_snapshot(
                        evaluated.state,
                        porosity(input),
                        evaluated.transport,
                        evaluated.caloric,
                        evaluated.rock);
            const auto energy_linearization =
                mpmc::flow::
                    build_two_phase_energy_accumulation_linearization(
                        evaluated.state,
                        porosity(input),
                        evaluated.transport,
                        evaluated.caloric,
                        evaluated.rock);
            const auto energy =
                mpmc::flow::
                    build_backward_euler_energy_accumulation_residual(
                        energy_current,
                        energy_linearization,
                        *previous_energy(input),
                        time_step_seconds_);
            add_conservation(
                component,
                energy);

            if (evaluated.fugacity.residual_count() !=
                    nc ||
                evaluated.fugacity.input_count !=
                    q) {
                throw std::invalid_argument(
                    "two-phase fugacity shape mismatch in mixed dispatcher");
            }
            for (std::size_t component_index = 0U;
                 component_index < nc;
                 ++component_index) {
                const std::size_t row =
                    nc + 1U +
                    component_index;
                assembly.residual[row] =
                    evaluated.fugacity
                        .residual[
                            component_index];
                for (std::size_t column = 0U;
                     column < q;
                     ++column) {
                    diagonal.values_row_major[
                        row * q +
                        column] =
                        evaluated.fugacity
                            .d_residual(
                                component_index,
                                column);
                }
            }
        } else {
            const auto& evaluated =
                std::get<
                    FixedThreePhaseCurrentCellLinearization3D>(
                        current);
            const auto component_current =
                mpmc::flow::
                    build_pore_volume_component_accumulation(
                        evaluated.state,
                        porosity(input));
            const auto component_linearization =
                mpmc::flow::
                    build_pore_volume_component_accumulation_linearization(
                        evaluated.state,
                        porosity(input),
                        evaluated.molar_density);
            const auto pair =
                mpmc::flow::
                    make_pore_volume_component_accumulation_pair(
                        component_current,
                        *previous_component(input));
            const auto component =
                mpmc::flow::
                    build_backward_euler_component_accumulation_residual(
                        pair,
                        component_linearization,
                        time_step_seconds_);
            const auto energy_current =
                mpmc::flow::
                    build_pore_volume_energy_accumulation_snapshot(
                        evaluated.state,
                        porosity(input),
                        evaluated.transport,
                        evaluated.caloric,
                        evaluated.rock);
            const auto energy_linearization =
                mpmc::flow::
                    build_pore_volume_energy_accumulation_linearization(
                        evaluated.state,
                        porosity(input),
                        evaluated.transport,
                        evaluated.caloric,
                        evaluated.rock);
            const auto energy =
                mpmc::flow::
                    build_backward_euler_energy_accumulation_residual(
                        energy_current,
                        energy_linearization,
                        *previous_energy(input),
                        time_step_seconds_);
            add_conservation(
                component,
                energy);

            if (evaluated.fugacity.residual_count() !=
                    2U * nc ||
                evaluated.fugacity.input_count() !=
                    q) {
                throw std::invalid_argument(
                    "three-phase fugacity shape mismatch in mixed dispatcher");
            }
            const auto values =
                evaluated.fugacity.values();
            const auto jacobian =
                evaluated.fugacity.jacobian();
            for (std::size_t local_row = 0U;
                 local_row < values.size();
                 ++local_row) {
                const std::size_t row =
                    nc + 1U +
                    local_row;
                assembly.residual[row] =
                    values[local_row];
                for (std::size_t column = 0U;
                     column < q;
                     ++column) {
                    diagonal.values_row_major[
                        row * q +
                        column] =
                        jacobian[
                            local_row * q +
                            column];
                }
            }
        }

        for (double value :
             assembly.residual) {
            if (!std::isfinite(value)) {
                throw std::range_error(
                    "mixed-cardinality physical cell residual is non-finite");
            }
        }
        for (double value :
             diagonal.values_row_major) {
            if (!std::isfinite(value)) {
                throw std::range_error(
                    "mixed-cardinality physical cell Jacobian is non-finite");
            }
        }

        assembly.jacobian_blocks.push_back(
            std::move(diagonal));
        return assembly;
    }

    [[nodiscard]]
    MixedCardinalityPhysicalFaceLinearization3D
    build_same_cardinality_face(
        const mpmc::discretization_petsc::
            AssemblyReadyInternalConnectionRow3D&
                row,
        const MixedCardinalityPhysicalSnesAuthoritativeFaceInput3D&
            face_input,
        const MixedCardinalityPhysicalCurrentCellLinearization3D&
            owner,
        const MixedCardinalityPhysicalCurrentCellLinearization3D&
            neighbour,
        const MixedCardinalityPhysicalMobilityLinearization3D&
            owner_mobility,
        const MixedCardinalityPhysicalMobilityLinearization3D&
            neighbour_mobility) const {
        using namespace
            mixed_cardinality_physical_detail;

        const auto owner_phase_count =
            numbering_->cell(
                row.owner_cell)
                .phase_count;
        const auto& owner_input =
            cell_inputs_.at(
                static_cast<std::size_t>(
                    row.owner_cell.value()));
        const auto& neighbour_input =
            cell_inputs_.at(
                static_cast<std::size_t>(
                    row.neighbour_cell.value()));

        if (owner_phase_count == 1U) {
            const auto& owner_eval =
                std::get<
                    SinglePhaseCurrentCellLinearization3D>(
                        owner);
            const auto& neighbour_eval =
                std::get<
                    SinglePhaseCurrentCellLinearization3D>(
                        neighbour);
            const auto potential =
                mpmc::flow::
                    build_single_phase_potential_upwind_linearization(
                        std::get<
                            mpmc::flow::
                                SinglePhaseMobilityLinearization>(
                                    owner_mobility),
                        std::get<
                            mpmc::flow::
                                SinglePhaseMobilityLinearization>(
                                    neighbour_mobility),
                        face_input.gravity,
                        face_input
                            .owner_to_neighbour_displacement);
            auto face =
                mpmc::flow_discretization::
                    build_single_phase_tpfa_face_linearization(
                        face_input.transmissibility,
                        potential,
                        owner_eval.state,
                        owner_eval.molar_density,
                        owner_eval.transport,
                        owner_eval.caloric,
                        neighbour_eval.state,
                        neighbour_eval.molar_density,
                        neighbour_eval.transport,
                        neighbour_eval.caloric,
                        {
                            bulk_volume(owner_input),
                            bulk_volume(neighbour_input)},
                        face_input
                            .thermal_conductance);
            return {
                std::move(face.component),
                std::move(face.energy)};
        }

        if (owner_phase_count == 2U) {
            const auto& owner_eval =
                std::get<
                    TwoPhaseCurrentCellLinearization3D>(
                        owner);
            const auto& neighbour_eval =
                std::get<
                    TwoPhaseCurrentCellLinearization3D>(
                        neighbour);
            const auto potential =
                mpmc::flow::
                    build_two_phase_potential_upwind_linearization(
                        std::get<
                            mpmc::flow::
                                TwoPhaseMobilityLinearization>(
                                    owner_mobility),
                        std::get<
                            mpmc::flow::
                                TwoPhaseMobilityLinearization>(
                                    neighbour_mobility),
                        face_input.gravity,
                        face_input
                            .owner_to_neighbour_displacement);
            auto face =
                mpmc::flow_discretization::
                    build_two_phase_tpfa_face_linearization(
                        face_input.transmissibility,
                        potential,
                        owner_eval.state,
                        owner_eval.molar_density,
                        owner_eval.transport,
                        owner_eval.caloric,
                        neighbour_eval.state,
                        neighbour_eval.molar_density,
                        neighbour_eval.transport,
                        neighbour_eval.caloric,
                        {
                            bulk_volume(owner_input),
                            bulk_volume(neighbour_input)},
                        face_input
                            .thermal_conductance);
            return {
                std::move(face.component),
                std::move(face.energy)};
        }

        const auto& owner_eval =
            std::get<
                FixedThreePhaseCurrentCellLinearization3D>(
                    owner);
        const auto& neighbour_eval =
            std::get<
                FixedThreePhaseCurrentCellLinearization3D>(
                    neighbour);
        const auto potential =
            mpmc::flow::
                build_two_cell_phase_potential_upwind_linearization(
                    std::get<
                        mpmc::flow::
                            LocalPhaseMobilityLinearization3P>(
                                owner_mobility),
                    std::get<
                        mpmc::flow::
                            LocalPhaseMobilityLinearization3P>(
                                neighbour_mobility),
                    face_input.gravity,
                    face_input
                        .owner_to_neighbour_displacement);
        const auto phase_flux =
            mpmc::flow_discretization::
                build_materialized_tpfa_internal_face_phase_darcy_flux(
                    face_input.transmissibility,
                    potential);
        const auto component_flux =
            mpmc::flow_discretization::
                build_materialized_tpfa_internal_face_component_molar_flux(
                    phase_flux,
                    owner_eval.state,
                    owner_eval.molar_density,
                    neighbour_eval.state,
                    neighbour_eval.molar_density);
        const auto scatter =
            mpmc::flow_discretization::
                scatter_component_molar_face_flux_conservatively(
                    component_flux);
        auto component =
            mpmc::flow_discretization::
                normalize_component_face_rate_by_bulk_volume(
                    scatter,
                    {
                        bulk_volume(owner_input),
                        bulk_volume(neighbour_input)});
        const auto energy_rate =
            mpmc::flow_discretization::
                build_internal_energy_face_rate(
                    phase_flux,
                    owner_eval.transport,
                    owner_eval.caloric,
                    neighbour_eval.transport,
                    neighbour_eval.caloric,
                    face_input
                        .thermal_conductance);
        auto energy =
            mpmc::flow_discretization::
                normalize_energy_face_rate_by_bulk_volume(
                    energy_rate,
                    {
                        bulk_volume(owner_input),
                        bulk_volume(neighbour_input)});
        return {
            std::move(component),
            std::move(energy)};
    }

    [[nodiscard]] PetscErrorCode
    insert_face(
        const mpmc::discretization_petsc::
            AssemblyReadyInternalConnectionRow3D&
                row,
        const MixedCardinalityPhysicalFaceLinearization3D&
            face,
        Vec residual,
        Mat jacobian,
        bool insert_residual,
        bool insert_jacobian) const {
        using namespace
            mixed_cardinality_physical_detail;

        const auto& owner =
            numbering_->cell(
                row.owner_cell);
        const auto& neighbour =
            numbering_->cell(
                row.neighbour_cell);
        const std::size_t nc =
            numbering_->component_count();

        if (face.component.face !=
                row.face ||
            face.energy.face !=
                row.face ||
            face.component.component_ids.size() !=
                nc ||
            face.component.owner_state_identity
                    .layout.unknown_count() !=
                owner.scalar_count ||
            face.component.neighbour_state_identity
                    .layout.unknown_count() !=
                neighbour.scalar_count ||
            face.energy.owner_state_identity
                    .layout.unknown_count() !=
                owner.scalar_count ||
            face.energy.neighbour_state_identity
                    .layout.unknown_count() !=
                neighbour.scalar_count) {
            return PETSC_ERR_ARG_INCOMP;
        }

        PetscErrorCode error =
            insert_face_side(
                owner,
                neighbour,
                nc,
                face.component
                    .owner_component_contribution_mol_per_bulk_m3_s,
                face.energy
                    .owner_contribution_w_per_bulk_m3,
                face.component
                    .owner_row_owner_column_jacobian,
                face.component
                    .owner_row_neighbour_column_jacobian,
                face.energy
                    .owner_row_owner_column_gradient,
                face.energy
                    .owner_row_neighbour_column_gradient,
                residual,
                jacobian,
                insert_residual,
                insert_jacobian);
        if (error != PETSC_SUCCESS) {
            return error;
        }

        return insert_face_side(
            neighbour,
            owner,
            nc,
            face.component
                .neighbour_component_contribution_mol_per_bulk_m3_s,
            face.energy
                .neighbour_contribution_w_per_bulk_m3,
            face.component
                .neighbour_row_neighbour_column_jacobian,
            face.component
                .neighbour_row_owner_column_jacobian,
            face.energy
                .neighbour_row_neighbour_column_gradient,
            face.energy
                .neighbour_row_owner_column_gradient,
            residual,
            jacobian,
            insert_residual,
            insert_jacobian);
    }

    [[nodiscard]] PetscErrorCode
    evaluate_and_insert(
        Vec global_state,
        Vec residual,
        Mat jacobian,
        bool insert_residual,
        bool insert_jacobian,
        NaturalVariableSnesEvaluationStatus3D*
            status) {
        using namespace
            mixed_cardinality_physical_detail;

        if (global_state == nullptr ||
            status == nullptr ||
            (insert_residual &&
             residual == nullptr) ||
            (insert_jacobian &&
             jacobian == nullptr)) {
            return PETSC_ERR_ARG_NULL;
        }

        std::vector<double>
            local_state;
        PetscErrorCode error =
            infrastructure_
                .copy_local_packed_state(
                    global_state,
                    &local_state);
        error =
            collective_error(
                comm_,
                error);
        if (error != PETSC_SUCCESS) {
            return error;
        }

        std::vector<std::optional<
            MixedCardinalityPhysicalCurrentCellLinearization3D>>
            current;
        std::vector<std::optional<
            MixedCardinalityPhysicalMobilityLinearization3D>>
            mobility;
        error =
            evaluate_cells(
                local_state,
                &current,
                &mobility,
                status);
        if (error != PETSC_SUCCESS ||
            *status ==
                NaturalVariableSnesEvaluationStatus3D::
                    domain_error) {
            return error;
        }

        PetscErrorCode local_error =
            PETSC_SUCCESS;
        try {
            for (std::size_t local = 0U;
                 local < current.size();
                 ++local) {
                const auto cell_index =
                    mpmc::mesh::LocalIndex{
                        static_cast<
                            mpmc::mesh::LocalIndex::value_type>(
                                local)};
                const auto& record =
                    numbering_->cell(
                        cell_index);
                if (record.owner_rank !=
                    numbering_->local_rank()) {
                    continue;
                }

                auto assembly =
                    build_owned_cell_assembly(
                        local,
                        *current[local]);
                local_error =
                    insert_dense_cell_assembly(
                        record,
                        assembly,
                        residual,
                        jacobian,
                        insert_residual,
                        insert_jacobian);
                if (local_error !=
                    PETSC_SUCCESS) {
                    break;
                }
            }

            if (local_error ==
                PETSC_SUCCESS) {
                for (std::size_t index = 0U;
                     index <
                     schedule_->assembly_rows().size();
                     ++index) {
                    const auto& row =
                        schedule_
                            ->assembly_rows()[index];
                    std::optional<
                        MixedCardinalityPhysicalFaceLinearization3D>
                        face;
                    const auto& phase_plan =
                        face_phase_identity_plans_
                            .at(index);
                    if (phase_plan
                            .slot_aligned_same_active_set()) {
                        face.emplace(
                            build_same_cardinality_face(
                                row,
                                face_inputs_[index],
                                *current.at(
                                    static_cast<std::size_t>(
                                        row.owner_cell.value())),
                                *current.at(
                                    static_cast<std::size_t>(
                                        row.neighbour_cell.value())),
                                *mobility.at(
                                    static_cast<std::size_t>(
                                        row.owner_cell.value())),
                                *mobility.at(
                                    static_cast<std::size_t>(
                                        row.neighbour_cell.value()))));
                    } else {
                        if (cross_phase_face_evaluator_
                                .evaluator ==
                            nullptr) {
                            local_error =
                                PETSC_ERR_SUP;
                            break;
                        }
                        NaturalVariableSnesEvaluationStatus3D
                            face_status =
                                NaturalVariableSnesEvaluationStatus3D::
                                    success;
                        local_error =
                            cross_phase_face_evaluator_
                                .evaluator(
                                    face_inputs_[index],
                                    phase_plan,
                                    {
                                        bulk_volume(
                                            cell_inputs_.at(
                                                static_cast<std::size_t>(
                                                    row.owner_cell.value()))),
                                        bulk_volume(
                                            cell_inputs_.at(
                                                static_cast<std::size_t>(
                                                    row.neighbour_cell.value())))},
                                    *current.at(
                                        static_cast<std::size_t>(
                                            row.owner_cell.value())),
                                    *current.at(
                                        static_cast<std::size_t>(
                                            row.neighbour_cell.value())),
                                    cross_phase_face_evaluator_
                                        .user_context,
                                    &face,
                                    &face_status);
                        if (local_error !=
                            PETSC_SUCCESS) {
                            break;
                        }
                        if (face_status ==
                            NaturalVariableSnesEvaluationStatus3D::
                                domain_error) {
                            *status =
                                face_status;
                            break;
                        }
                        if (face_status !=
                                NaturalVariableSnesEvaluationStatus3D::
                                    success ||
                            !face.has_value()) {
                            local_error =
                                PETSC_ERR_ARG_INCOMP;
                            break;
                        }
                    }

                    local_error =
                        insert_face(
                            row,
                            *face,
                            residual,
                            jacobian,
                            insert_residual,
                            insert_jacobian);
                    if (local_error !=
                        PETSC_SUCCESS) {
                        break;
                    }
                }
            }
        } catch (...) {
            local_error =
                PETSC_ERR_ARG_INCOMP;
        }

        error =
            collective_error(
                comm_,
                local_error);
        if (error != PETSC_SUCCESS) {
            return error;
        }

        int local_domain =
            *status ==
                    NaturalVariableSnesEvaluationStatus3D::
                        domain_error
                ? 1
                : 0;
        int global_domain = 0;
        if (MPI_Allreduce(
                &local_domain,
                &global_domain,
                1,
                MPI_INT,
                MPI_MAX,
                comm_) != MPI_SUCCESS) {
            return PETSC_ERR_MPI;
        }
        *status =
            global_domain != 0
                ? NaturalVariableSnesEvaluationStatus3D::
                      domain_error
                : NaturalVariableSnesEvaluationStatus3D::
                      success;
        return PETSC_SUCCESS;
    }

    [[nodiscard]] PetscErrorCode
    positive_support_precheck(
        Vec state,
        Vec search_direction,
        PetscBool* changed_direction) {
        if (state == nullptr ||
            search_direction == nullptr ||
            changed_direction == nullptr) {
            return PETSC_ERR_ARG_NULL;
        }
        *changed_direction =
            PETSC_FALSE;

        PetscInt state_local = -1;
        PetscInt direction_local = -1;
        PetscErrorCode error =
            VecGetLocalSize(
                state,
                &state_local);
        if (error == PETSC_SUCCESS) {
            error =
                VecGetLocalSize(
                    search_direction,
                    &direction_local);
        }
        if (error != PETSC_SUCCESS) {
            return error;
        }
        if (state_local !=
                direction_local ||
            state_local !=
                numbering_
                    ->petsc_local_owned_scalar_count()) {
            return PETSC_ERR_ARG_SIZ;
        }

        const PetscScalar* state_values =
            nullptr;
        const PetscScalar* direction_values =
            nullptr;
        error =
            VecGetArrayRead(
                state,
                &state_values);
        if (error == PETSC_SUCCESS) {
            error =
                VecGetArrayRead(
                    search_direction,
                    &direction_values);
        }
        if (error != PETSC_SUCCESS) {
            if (state_values != nullptr) {
                (void)VecRestoreArrayRead(
                    state,
                    &state_values);
            }
            return error;
        }

        double local_scale = 1.0;
        bool invalid = false;
        for (PetscInt index = 0;
             index < state_local;
             ++index) {
            const double x =
                static_cast<double>(
                    PetscRealPart(
                        state_values[index]));
            const double y =
                static_cast<double>(
                    PetscRealPart(
                        direction_values[index]));
            if (!std::isfinite(x) ||
                !std::isfinite(y) ||
                !(x > 0.0)) {
                invalid = true;
                break;
            }
            if (y > 0.0 &&
                x - y <= 0.0) {
                local_scale =
                    std::min(
                        local_scale,
                        0.8 * x / y);
            }
        }

        const PetscErrorCode direction_restore =
            VecRestoreArrayRead(
                search_direction,
                &direction_values);
        const PetscErrorCode state_restore =
            VecRestoreArrayRead(
                state,
                &state_values);
        if (direction_restore !=
            PETSC_SUCCESS) {
            return direction_restore;
        }
        if (state_restore !=
            PETSC_SUCCESS) {
            return state_restore;
        }
        if (invalid) {
            return PETSC_ERR_ARG_OUTOFRANGE;
        }

        double global_scale = 1.0;
        if (MPI_Allreduce(
                &local_scale,
                &global_scale,
                1,
                MPI_DOUBLE,
                MPI_MIN,
                comm_) != MPI_SUCCESS) {
            return PETSC_ERR_MPI;
        }
        if (!std::isfinite(
                global_scale) ||
            !(global_scale > 0.0)) {
            return PETSC_ERR_FP;
        }
        if (global_scale < 1.0) {
            error =
                VecScale(
                    search_direction,
                    static_cast<PetscScalar>(
                        global_scale));
            if (error != PETSC_SUCCESS) {
                return error;
            }
            *changed_direction =
                PETSC_TRUE;
        }
        return PETSC_SUCCESS;
    }

    static PetscErrorCode
    snes_function(
        Vec state,
        Vec residual,
        void* raw_context,
        NaturalVariableSnesEvaluationStatus3D*
            status) {
        if (raw_context == nullptr) {
            return PETSC_ERR_ARG_NULL;
        }
        return static_cast<
            MixedCardinalityPhysicalSnesAssemblyContext3D*>(
                raw_context)
            ->evaluate_and_insert(
                state,
                residual,
                nullptr,
                true,
                false,
                status);
    }

    static PetscErrorCode
    snes_jacobian(
        Vec state,
        Mat jacobian,
        void* raw_context,
        NaturalVariableSnesEvaluationStatus3D*
            status) {
        if (raw_context == nullptr) {
            return PETSC_ERR_ARG_NULL;
        }
        return static_cast<
            MixedCardinalityPhysicalSnesAssemblyContext3D*>(
                raw_context)
            ->evaluate_and_insert(
                state,
                nullptr,
                jacobian,
                false,
                true,
                status);
    }

    static PetscErrorCode
    snes_precheck(
        Vec state,
        Vec search_direction,
        void* raw_context,
        PetscBool* changed_direction) {
        if (raw_context == nullptr) {
            return PETSC_ERR_ARG_NULL;
        }
        return static_cast<
            MixedCardinalityPhysicalSnesAssemblyContext3D*>(
                raw_context)
            ->positive_support_precheck(
                state,
                search_direction,
                changed_direction);
    }

    MPI_Comm comm_;
    const mpmc::discretization_petsc::
        ParallelOwnedConnectionSchedule3D*
            schedule_;
    const mpmc::mesh::PartitionSnapshot*
        partition_;
    const VariableCardinalityNaturalVariableNumbering3D*
        numbering_;
    double time_step_seconds_{};
    std::vector<
        MixedCardinalityPhysicalSnesCellInput3D>
        cell_inputs_;
    std::vector<
        mpmc::flow::FrozenActivePhaseIdentityMap>
        phase_identity_maps_;
    std::vector<
        mpmc::flow::
            CrossCardinalityFacePhaseIdentityPlan>
        face_phase_identity_plans_;
    std::vector<
        MixedCardinalityPhysicalSnesAuthoritativeFaceInput3D>
        face_inputs_;
    MixedCardinalityPhysicalCellEvaluatorBindings3D
        cell_evaluators_;
    MixedCardinalityCrossPhaseFaceEvaluatorBinding3D
        cross_phase_face_evaluator_;
    VariableCardinalityNaturalVariableSnesAssemblyContext3D
        infrastructure_;
};

} // namespace mpmc::flow_discretization_petsc

#endif // MPMC_FLOW_DISCRETIZATION_PETSC_MIXED_CARDINALITY_PHYSICAL_SNES_ASSEMBLY_HPP
