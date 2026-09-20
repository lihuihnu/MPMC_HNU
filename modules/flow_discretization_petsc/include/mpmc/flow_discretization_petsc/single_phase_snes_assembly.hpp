#ifndef MPMC_FLOW_DISCRETIZATION_PETSC_SINGLE_PHASE_SNES_ASSEMBLY_HPP
#define MPMC_FLOW_DISCRETIZATION_PETSC_SINGLE_PHASE_SNES_ASSEMBLY_HPP

#include <mpmc/flow/single_phase_natural_variable.hpp>
#include <mpmc/flow_discretization/single_phase_tpfa.hpp>
#include <mpmc/flow_discretization_petsc/distributed_component_conservation.hpp>
#include <mpmc/flow_discretization_petsc/distributed_energy_conservation.hpp>
#include <mpmc/flow_discretization_petsc/energy_global_assembly_mapping.hpp>
#include <mpmc/flow_discretization_petsc/fugacity_equilibrium_global_assembly_mapping.hpp>
#include <mpmc/flow_discretization_petsc/global_component_assembly_mapping.hpp>
#include <mpmc/flow_discretization_petsc/natural_variable_snes_solver.hpp>

#include <petscsf.h>

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
    single_phase_snes_assembly_convention =
        "flow_discretization_petsc/single-phase-snes-assembly/v1";

struct SinglePhaseCurrentCellLinearization3D {
    mpmc::flow::NaturalVariableCellState1P state;
    mpmc::flow::
        SinglePhaseMolarDensityNaturalVariableLinearization
            molar_density;
    mpmc::flow::
        SinglePhaseTransportNaturalVariableLinearization
            transport;
    mpmc::flow::
        SinglePhaseCaloricNaturalVariableLinearization
            caloric;
    mpmc::flow::
        SinglePhaseRockThermalStorageLinearization
            rock;
};

using SinglePhaseCurrentCellEvaluator3D =
    PetscErrorCode (*)(
        mpmc::mesh::LocalIndex cell,
        mpmc::mesh::GlobalEntityId cell_global,
        std::span<const double> natural_variables,
        const mpmc::flow::NaturalVariableLayout1P&
            frozen_layout,
        std::span<const std::string> component_ids,
        void* user_context,
        std::optional<
            SinglePhaseCurrentCellLinearization3D>* output,
        NaturalVariableSnesEvaluationStatus3D* status);

struct SinglePhaseSnesCellInput3D {
    mpmc::mesh::LocalIndex cell{
        mpmc::mesh::LocalIndex::value_type{0}};
    mpmc::mesh::GlobalEntityId cell_global{
        mpmc::mesh::GlobalEntityId::value_type{0}};
    double bulk_volume_m3{};
    double porosity{};
    mpmc::flow::NaturalVariableLayout1P
        frozen_layout{std::size_t{2U}};
    std::vector<std::string> component_ids;
    std::optional<
        mpmc::flow::
            PoreVolumeComponentAccumulationSnapshot3P>
        previous_component_accumulation;
    std::optional<
        mpmc::flow::
            PoreVolumeEnergyAccumulationSnapshot3P>
        previous_energy_accumulation;
};

struct SinglePhaseSnesAuthoritativeFaceInput3D {
    mpmc::mesh::LocalIndex face{
        mpmc::mesh::LocalIndex::value_type{0}};
    mpmc::mesh::GlobalEntityId face_global{
        mpmc::mesh::GlobalEntityId::value_type{0}};
    mpmc::discretization::
        TpfaInternalFaceTransmissibilityEntry3D
            transmissibility;
    mpmc::flow::GravityVector3D gravity;
    mpmc::flow::OwnerToNeighbourDisplacement3D
        owner_to_neighbour_displacement;
    mpmc::flow_discretization::
        StaticThermalFaceConductance3D
            thermal_conductance;
};

struct SinglePhaseCurrentCellEvaluatorBinding3D {
    SinglePhaseCurrentCellEvaluator3D evaluator{};
    void* user_context{};
};

namespace single_phase_snes_assembly_detail {

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
same_layout(
    const mpmc::flow::NaturalVariableLayout1P& first,
    const mpmc::flow::NaturalVariableLayout1P& second) {
    return first.component_count() ==
               second.component_count() &&
        first.unknown_count() ==
            second.unknown_count() &&
        first.dependent_composition_component() ==
            second.dependent_composition_component();
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
            std::numeric_limits<double>::epsilon() *
            scale;
}

} // namespace single_phase_snes_assembly_detail

class SinglePhaseSnesAssemblyContext3D {
public:
    static constexpr std::string_view convention =
        single_phase_snes_assembly_convention;

    SinglePhaseSnesAssemblyContext3D(
        const SinglePhaseSnesAssemblyContext3D&) =
        delete;
    SinglePhaseSnesAssemblyContext3D& operator=(
        const SinglePhaseSnesAssemblyContext3D&) =
        delete;
    SinglePhaseSnesAssemblyContext3D& operator=(
        SinglePhaseSnesAssemblyContext3D&&) =
        delete;

    SinglePhaseSnesAssemblyContext3D(
        SinglePhaseSnesAssemblyContext3D&& other) noexcept
        : comm_(other.comm_),
          schedule_(other.schedule_),
          partition_(other.partition_),
          dof_layout_(other.dof_layout_),
          dof_numbering_(other.dof_numbering_),
          cell_bridge_(other.cell_bridge_),
          cell_pattern_(other.cell_pattern_),
          natural_variable_id_(
              std::move(other.natural_variable_id_)),
          time_step_seconds_(
              other.time_step_seconds_),
          cell_inputs_(
              std::move(other.cell_inputs_)),
          face_inputs_(
              std::move(other.face_inputs_)),
          cell_input_by_local_(
              std::move(other.cell_input_by_local_)),
          face_input_by_local_(
              std::move(other.face_input_by_local_)),
          cell_evaluator_(other.cell_evaluator_),
          q_(other.q_),
          component_count_(
              other.component_count_),
          state_sf_(other.state_sf_),
          local_state_(other.local_state_) {
        other.state_sf_ = nullptr;
        other.local_state_ = nullptr;
    }

    ~SinglePhaseSnesAssemblyContext3D() {
        if (local_state_ != nullptr) {
            (void)VecDestroy(&local_state_);
        }
        if (state_sf_ != nullptr) {
            (void)PetscSFDestroy(&state_sf_);
        }
    }

    [[nodiscard]] static PetscErrorCode
    create(
        MPI_Comm comm,
        const mpmc::discretization_petsc::
            ParallelOwnedConnectionSchedule3D& schedule,
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
        std::string natural_variable_id,
        double time_step_seconds,
        std::vector<SinglePhaseSnesCellInput3D>
            cell_inputs,
        std::vector<
            SinglePhaseSnesAuthoritativeFaceInput3D>
            face_inputs,
        SinglePhaseCurrentCellEvaluatorBinding3D
            cell_evaluator,
        std::optional<
            SinglePhaseSnesAssemblyContext3D>*
            output) {
        using namespace
            single_phase_snes_assembly_detail;

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
        std::size_t q = 0U;
        std::size_t component_count = 0U;
        std::vector<std::optional<std::size_t>>
            cell_lookup;
        std::vector<std::optional<std::size_t>>
            face_lookup;

        try {
            if (cell_evaluator.evaluator == nullptr ||
                natural_variable_id.empty() ||
                !std::isfinite(time_step_seconds) ||
                !(time_step_seconds > 0.0) ||
                schedule.local_rank() !=
                    partition.local_rank() ||
                schedule.rank_count() !=
                    partition.rank_count() ||
                dof_numbering.local_rank() !=
                    partition.local_rank() ||
                dof_numbering.rank_count() !=
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
                    partition.entity_count(
                        mpmc::mesh::EntityKind::cell) ||
                face_inputs.size() !=
                    schedule.assembly_rows().size() ||
                cell_pattern.local_cell_count() !=
                    partition.entity_count(
                        mpmc::mesh::EntityKind::cell) ||
                dof_numbering.local_dof_count() !=
                    dof_layout.total_dof_count() ||
                !dof_layout.contains(
                    natural_variable_id)) {
                throw std::invalid_argument(
                    "single-phase SNES assembly metadata mismatch");
            }

            const std::size_t variable_index =
                dof_layout.variable_index(
                    natural_variable_id);
            const auto& variable =
                dof_layout.variable(
                    variable_index);
            q = variable.component_count;
            if (variable.location !=
                    mpmc::mesh::EntityKind::cell ||
                q < 3U ||
                dof_layout.dofs_per_entity(
                    mpmc::mesh::EntityKind::cell) !=
                    q) {
                throw std::invalid_argument(
                    "single-phase natural-variable DoF block must occupy the complete cell block");
            }
            component_count =
                q - 1U;
            if (component_count < 2U) {
                throw std::invalid_argument(
                    "single-phase SNES assembly requires at least two components");
            }

            cell_lookup.assign(
                cell_inputs.size(),
                std::nullopt);
            std::vector<std::string>
                canonical_component_ids;

            for (std::size_t index = 0U;
                 index < cell_inputs.size();
                 ++index) {
                const auto& input =
                    cell_inputs[index];
                const std::size_t local =
                    static_cast<std::size_t>(
                        input.cell.value());
                if (local >= cell_lookup.size() ||
                    cell_lookup[local].has_value() ||
                    partition.global_id(
                        mpmc::mesh::EntityKind::cell,
                        input.cell) !=
                        input.cell_global ||
                    !std::isfinite(
                        input.bulk_volume_m3) ||
                    !(input.bulk_volume_m3 > 0.0) ||
                    !std::isfinite(input.porosity) ||
                    !(input.porosity > 0.0) ||
                    !(input.porosity < 1.0) ||
                    input.component_ids.size() !=
                        component_count ||
                    input.frozen_layout
                            .component_count() !=
                        component_count ||
                    input.frozen_layout
                            .unknown_count() !=
                        q) {
                    throw std::invalid_argument(
                        "invalid single-phase cell binding");
                }

                if (index == 0U) {
                    canonical_component_ids =
                        input.component_ids;
                } else if (
                    input.component_ids !=
                    canonical_component_ids) {
                    throw std::invalid_argument(
                        "single-phase component identity/order is not canonical across local overlap");
                }

                const bool owned =
                    partition.is_owned(
                        mpmc::mesh::EntityKind::cell,
                        input.cell);
                if (owned !=
                        input.previous_component_accumulation
                            .has_value() ||
                    owned !=
                        input.previous_energy_accumulation
                            .has_value()) {
                    throw std::invalid_argument(
                        "owned single-phase cell history is required and ghost history is forbidden");
                }
                if (owned) {
                    if (input
                            .previous_component_accumulation
                            ->component_ids !=
                            input.component_ids ||
                        !near_roundoff(
                            input
                                .previous_component_accumulation
                                ->porosity,
                            input.porosity) ||
                        input
                            .previous_energy_accumulation
                            ->state_identity
                            .component_ids !=
                            input.component_ids ||
                        input
                            .previous_energy_accumulation
                            ->state_identity
                            .layout.phase_count() !=
                            1U ||
                        !near_roundoff(
                            input
                                .previous_energy_accumulation
                                ->porosity,
                            input.porosity)) {
                        throw std::invalid_argument(
                            "single-phase previous accumulation does not match frozen cell identity");
                    }
                }
                cell_lookup[local] =
                    index;
            }
            if (std::find(
                    cell_lookup.begin(),
                    cell_lookup.end(),
                    std::nullopt) !=
                cell_lookup.end()) {
                throw std::invalid_argument(
                    "single-phase local cell bindings are incomplete");
            }

            face_lookup.assign(
                partition.entity_count(
                    mpmc::mesh::EntityKind::face),
                std::nullopt);
            for (std::size_t index = 0U;
                 index < face_inputs.size();
                 ++index) {
                const auto& input =
                    face_inputs[index];
                const std::size_t local =
                    static_cast<std::size_t>(
                        input.face.value());
                if (local >= face_lookup.size() ||
                    face_lookup[local].has_value() ||
                    partition.global_id(
                        mpmc::mesh::EntityKind::face,
                        input.face) !=
                        input.face_global ||
                    input.transmissibility.face !=
                        input.face ||
                    !std::isfinite(
                        input.thermal_conductance
                            .conductance_w_per_k) ||
                    input.thermal_conductance
                            .conductance_w_per_k <
                        0.0) {
                    throw std::invalid_argument(
                        "invalid single-phase authoritative face binding");
                }
                face_lookup[local] =
                    index;
            }

            for (const auto& row :
                 schedule.assembly_rows()) {
                const std::size_t local =
                    static_cast<std::size_t>(
                        row.face.value());
                if (local >= face_lookup.size() ||
                    !face_lookup[local]
                         .has_value()) {
                    throw std::invalid_argument(
                        "single-phase authoritative schedule face lacks static binding");
                }
                const auto& input =
                    face_inputs[
                        *face_lookup[local]];
                if (input.face_global !=
                        row.face_global ||
                    !input.transmissibility
                         .static_transmissibility
                         .has_value() ||
                    !near_roundoff(
                        input.transmissibility
                            .static_transmissibility
                            ->face_transmissibility_m3,
                        row.transmissibility_m3)) {
                    throw std::invalid_argument(
                        "single-phase face transmissibility does not match schedule");
                }
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

        SinglePhaseSnesAssemblyContext3D
            context{
                comm,
                schedule,
                partition,
                dof_layout,
                dof_numbering,
                cell_bridge,
                cell_pattern,
                std::move(natural_variable_id),
                time_step_seconds,
                std::move(cell_inputs),
                std::move(face_inputs),
                std::move(cell_lookup),
                std::move(face_lookup),
                cell_evaluator,
                q,
                component_count};

        error =
            context.initialize_state_exchange();
        error =
            collective_error(
                comm,
                error);
        if (error != PETSC_SUCCESS) {
            return error;
        }

        output->emplace(
            std::move(context));
        return PETSC_SUCCESS;
    }

    [[nodiscard]]
    NaturalVariableSnesEvaluator3D
    snes_evaluator() noexcept {
        return {
            &snes_function,
            &snes_jacobian,
            &snes_precheck,
            this};
    }

    [[nodiscard]] PetscErrorCode
    evaluate_complete_assembly(
        Vec global_state,
        std::optional<
            CompleteNaturalVariableAssemblySnapshot3D>*
            output,
        NaturalVariableSnesEvaluationStatus3D*
            status) {
        using namespace
            single_phase_snes_assembly_detail;

        if (global_state == nullptr ||
            output == nullptr ||
            status == nullptr) {
            return PETSC_ERR_ARG_NULL;
        }
        output->reset();
        *status =
            NaturalVariableSnesEvaluationStatus3D::
                success;

        PetscErrorCode error =
            broadcast_state(global_state);
        error =
            collective_error(
                comm_,
                error);
        if (error != PETSC_SUCCESS) {
            return error;
        }

        const std::size_t local_cell_count =
            partition_->entity_count(
                mpmc::mesh::EntityKind::cell);
        std::vector<double> local_scalars(
            local_cell_count * q_,
            0.0);
        const PetscScalar* local_array =
            nullptr;
        error =
            VecGetArrayRead(
                local_state_,
                &local_array);
        if (error != PETSC_SUCCESS) {
            return error;
        }
        for (std::size_t index = 0U;
             index < local_scalars.size();
             ++index) {
            local_scalars[index] =
                static_cast<double>(
                    PetscRealPart(
                        local_array[index]));
        }
        const PetscErrorCode restore_error =
            VecRestoreArrayRead(
                local_state_,
                &local_array);
        if (restore_error !=
            PETSC_SUCCESS) {
            return restore_error;
        }

        std::vector<std::optional<
            SinglePhaseCurrentCellLinearization3D>>
            current(cell_inputs_.size());
        PetscErrorCode local_error =
            PETSC_SUCCESS;
        int local_domain = 0;

        for (std::size_t local = 0U;
             local < current.size();
             ++local) {
            const auto& input =
                cell_inputs_.at(
                    *cell_input_by_local_[local]);
            NaturalVariableSnesEvaluationStatus3D
                cell_status =
                    NaturalVariableSnesEvaluationStatus3D::
                        success;
            const std::span<const double> values{
                local_scalars.data() +
                    local * q_,
                q_};
            const PetscErrorCode cell_error =
                cell_evaluator_.evaluator(
                    input.cell,
                    input.cell_global,
                    values,
                    input.frozen_layout,
                    input.component_ids,
                    cell_evaluator_.user_context,
                    &current[local],
                    &cell_status);
            if (cell_error != PETSC_SUCCESS) {
                local_error =
                    cell_error;
                break;
            }
            if (cell_status ==
                NaturalVariableSnesEvaluationStatus3D::
                    domain_error) {
                local_domain = 1;
                continue;
            }
            if (cell_status !=
                    NaturalVariableSnesEvaluationStatus3D::
                        success ||
                !current[local].has_value()) {
                local_error =
                    PETSC_ERR_ARG_INCOMP;
                break;
            }
            try {
                validate_current_cell(
                    input,
                    *current[local]);
            } catch (...) {
                local_error =
                    PETSC_ERR_ARG_INCOMP;
                break;
            }
        }

        error =
            collective_error(
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
            return PETSC_SUCCESS;
        }

        std::vector<std::optional<
            mpmc::flow::
                BackwardEulerComponentAccumulationResidual3P>>
            component_accumulation(
                current.size());
        std::vector<std::optional<
            mpmc::flow::
                BackwardEulerEnergyAccumulationResidual3P>>
            energy_accumulation(
                current.size());
        std::vector<std::optional<
            mpmc::flow::
                SinglePhaseMobilityLinearization>>
            mobility(
                current.size());

        try {
            for (std::size_t local = 0U;
                 local < current.size();
                 ++local) {
                const auto& static_input =
                    cell_inputs_.at(
                        *cell_input_by_local_[local]);
                const auto& evaluation =
                    *current[local];

                mobility[local].emplace(
                    mpmc::flow::
                        build_single_phase_mobility_linearization(
                            evaluation.state,
                            evaluation.transport));

                if (!partition_->is_owned(
                        mpmc::mesh::EntityKind::cell,
                        static_input.cell)) {
                    continue;
                }

                const auto current_component =
                    mpmc::flow::
                        build_single_phase_component_accumulation(
                            evaluation.state,
                            static_input.porosity);
                const auto current_component_linearization =
                    mpmc::flow::
                        build_single_phase_component_accumulation_linearization(
                            evaluation.state,
                            static_input.porosity,
                            evaluation.molar_density);
                const auto pair =
                    mpmc::flow::
                        make_pore_volume_component_accumulation_pair(
                            current_component,
                            *static_input
                                 .previous_component_accumulation);
                component_accumulation[local].emplace(
                    mpmc::flow::
                        build_backward_euler_component_accumulation_residual(
                            pair,
                            current_component_linearization,
                            time_step_seconds_));

                const auto current_energy =
                    mpmc::flow::
                        build_single_phase_energy_accumulation_snapshot(
                            evaluation.state,
                            static_input.porosity,
                            evaluation.transport,
                            evaluation.caloric,
                            evaluation.rock);
                const auto current_energy_linearization =
                    mpmc::flow::
                        build_single_phase_energy_accumulation_linearization(
                            evaluation.state,
                            static_input.porosity,
                            evaluation.transport,
                            evaluation.caloric,
                            evaluation.rock);
                energy_accumulation[local].emplace(
                    mpmc::flow::
                        build_backward_euler_energy_accumulation_residual(
                            current_energy,
                            current_energy_linearization,
                            *static_input
                                 .previous_energy_accumulation,
                            time_step_seconds_));
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

        std::vector<
            mpmc::flow_discretization::
                NormalizedComponentFaceContributionLinearization3D>
            component_face_contributions;
        std::vector<
            mpmc::flow_discretization::
                NormalizedEnergyFaceContributionLinearization3D>
            energy_face_contributions;
        component_face_contributions.reserve(
            schedule_->assembly_rows().size());
        energy_face_contributions.reserve(
            schedule_->assembly_rows().size());

        try {
            for (const auto& row :
                 schedule_->assembly_rows()) {
                const auto& face_input =
                    face_inputs_.at(
                        *face_input_by_local_.at(
                            static_cast<std::size_t>(
                                row.face.value())));
                const std::size_t owner =
                    static_cast<std::size_t>(
                        row.owner_cell.value());
                const std::size_t neighbour =
                    static_cast<std::size_t>(
                        row.neighbour_cell.value());
                const auto& owner_eval =
                    *current.at(owner);
                const auto& neighbour_eval =
                    *current.at(neighbour);

                const auto potential =
                    mpmc::flow::
                        build_single_phase_potential_upwind_linearization(
                            *mobility.at(owner),
                            *mobility.at(neighbour),
                            face_input.gravity,
                            face_input
                                .owner_to_neighbour_displacement);

                const auto& owner_static =
                    cell_inputs_.at(
                        *cell_input_by_local_.at(
                            owner));
                const auto& neighbour_static =
                    cell_inputs_.at(
                        *cell_input_by_local_.at(
                            neighbour));

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
                                owner_static
                                    .bulk_volume_m3,
                                neighbour_static
                                    .bulk_volume_m3},
                            face_input
                                .thermal_conductance);
                component_face_contributions.push_back(
                    std::move(face.component));
                energy_face_contributions.push_back(
                    std::move(face.energy));
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

        std::vector<DistributedCellStateBinding3D>
            component_cells;
        std::vector<
            DistributedEnergyCellStateBinding3D>
            energy_cells;
        component_cells.reserve(
            current.size());
        energy_cells.reserve(
            current.size());

        for (std::size_t local = 0U;
             local < current.size();
             ++local) {
            const auto& static_input =
                cell_inputs_.at(
                    *cell_input_by_local_[local]);
            const auto& identity =
                current[local]
                    ->transport
                    .state_identity;
            component_cells.push_back(
                {
                    static_input.cell,
                    static_input.cell_global,
                    static_input.bulk_volume_m3,
                    identity,
                    component_accumulation[local]
                            .has_value()
                        ? &*component_accumulation[
                              local]
                        : nullptr});
            energy_cells.push_back(
                {
                    static_input.cell,
                    static_input.cell_global,
                    static_input.bulk_volume_m3,
                    identity,
                    energy_accumulation[local]
                            .has_value()
                        ? &*energy_accumulation[
                              local]
                        : nullptr});
        }

        std::vector<
            mpmc::flow_discretization::
                AuthoritativeNormalizedFaceContributionBinding3D>
            component_faces;
        std::vector<
            AuthoritativeNormalizedEnergyFaceBinding3D>
            energy_faces;
        component_faces.reserve(
            schedule_->assembly_rows().size());
        energy_faces.reserve(
            schedule_->assembly_rows().size());

        for (std::size_t index = 0U;
             index <
             schedule_->assembly_rows().size();
             ++index) {
            const auto& row =
                schedule_->assembly_rows()[index];
            component_faces.push_back(
                {
                    row.face,
                    row.face_global,
                    &component_face_contributions[
                        index]});
            energy_faces.push_back(
                {
                    row.face,
                    row.face_global,
                    &energy_face_contributions[
                        index]});
        }

        std::optional<
            DistributedOwnedMultiCellComponentConservationSnapshot3D>
            component_conservation;
        error =
            make_distributed_owned_component_conservation_snapshot_3d(
                comm_,
                *schedule_,
                *partition_,
                component_cells,
                component_faces,
                &component_conservation);
        if (error != PETSC_SUCCESS) {
            return error;
        }

        std::optional<
            DistributedOwnedMultiCellEnergyConservationSnapshot3D>
            energy_conservation;
        error =
            make_distributed_owned_energy_conservation_snapshot_3d(
                comm_,
                *schedule_,
                *partition_,
                energy_cells,
                energy_faces,
                &energy_conservation);
        if (error != PETSC_SUCCESS) {
            return error;
        }

        std::optional<
            ComponentConservationGlobalAssemblyEntries3D>
            component_global;
        error =
            make_component_conservation_global_assembly_entries_3d(
                comm_,
                *component_conservation,
                *partition_,
                *dof_layout_,
                *dof_numbering_,
                *cell_bridge_,
                *cell_pattern_,
                natural_variable_id_,
                &component_global);
        if (error != PETSC_SUCCESS) {
            return error;
        }

        std::optional<
            EnergyConservationGlobalAssemblyEntries3D>
            energy_global;
        error =
            make_energy_conservation_global_assembly_entries_3d(
                comm_,
                *energy_conservation,
                *partition_,
                *dof_layout_,
                *dof_numbering_,
                *cell_bridge_,
                *cell_pattern_,
                natural_variable_id_,
                &energy_global);
        if (error != PETSC_SUCCESS) {
            return error;
        }

        FugacityEquilibriumGlobalAssemblyEntries3D
            empty_fugacity{
                component_global->local_rank(),
                component_global->rank_count(),
                std::string{
                    component_global
                        ->natural_variable_id()},
                component_global->component_count(),
                component_global
                    ->natural_variable_count(),
                component_global
                    ->petsc_scalar_row_start(),
                component_global
                    ->petsc_scalar_row_end(),
                component_global
                    ->petsc_scalar_row_count(),
                {},
                {}};

        return make_complete_natural_variable_assembly_snapshot_3d(
            comm_,
            *component_global,
            *energy_global,
            empty_fugacity,
            *partition_,
            *cell_bridge_,
            *cell_pattern_,
            output);
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
        const PetscInt expected =
            cell_bridge_->local_owned_row_count() *
            static_cast<PetscInt>(q_);
        if (error != PETSC_SUCCESS) {
            return error;
        }
        if (state_local != expected ||
            direction_local != expected) {
            return PETSC_ERR_ARG_SIZ;
        }

        const PetscScalar* x = nullptr;
        const PetscScalar* y = nullptr;
        error =
            VecGetArrayRead(state, &x);
        if (error == PETSC_SUCCESS) {
            error =
                VecGetArrayRead(
                    search_direction,
                    &y);
        }
        if (error != PETSC_SUCCESS) {
            if (x != nullptr) {
                (void)VecRestoreArrayRead(
                    state,
                    &x);
            }
            return error;
        }

        double local_scale = 1.0;
        try {
            for (std::size_t block = 0U;
                 block <
                 cell_bridge_
                     ->owned_cells_in_petsc_row_order()
                     .size();
                 ++block) {
                const auto cell =
                    cell_bridge_
                        ->owned_cells_in_petsc_row_order()[
                            block];
                const auto& input =
                    cell_inputs_.at(
                        *cell_input_by_local_.at(
                            static_cast<std::size_t>(
                                cell.value())));
                const auto& layout =
                    input.frozen_layout;
                const std::size_t base =
                    block * q_;

                auto limit_positive =
                    [&](std::size_t slot) {
                        const double value =
                            static_cast<double>(
                                PetscRealPart(
                                    x[base + slot]));
                        const double direction =
                            static_cast<double>(
                                PetscRealPart(
                                    y[base + slot]));
                        if (!std::isfinite(value) ||
                            !std::isfinite(direction) ||
                            !(value > 0.0)) {
                            throw std::invalid_argument(
                                "single-phase current state is outside positive support");
                        }
                        if (direction > 0.0 &&
                            value - direction <=
                                0.0) {
                            local_scale =
                                std::min(
                                    local_scale,
                                    0.8 *
                                        value /
                                        direction);
                        }
                    };

                limit_positive(
                    layout
                        .pressure_unknown_index());
                limit_positive(
                    layout
                        .temperature_unknown_index());

                double independent_sum = 0.0;
                double direction_sum = 0.0;
                for (std::size_t component = 0U;
                     component < component_count_;
                     ++component) {
                    const auto column =
                        layout
                            .independent_composition_unknown_index(
                                component);
                    if (!column) {
                        continue;
                    }
                    limit_positive(*column);
                    independent_sum +=
                        static_cast<double>(
                            PetscRealPart(
                                x[base + *column]));
                    direction_sum +=
                        static_cast<double>(
                            PetscRealPart(
                                y[base + *column]));
                }
                const double dependent =
                    1.0 -
                    independent_sum;
                if (!std::isfinite(dependent) ||
                    !std::isfinite(direction_sum) ||
                    !(dependent > 0.0)) {
                    throw std::invalid_argument(
                        "single-phase dependent composition is outside positive support");
                }
                if (direction_sum < 0.0 &&
                    dependent + direction_sum <=
                        0.0) {
                    local_scale =
                        std::min(
                            local_scale,
                            0.8 *
                                dependent /
                                (-direction_sum));
                }
            }
        } catch (...) {
            (void)VecRestoreArrayRead(
                search_direction,
                &y);
            (void)VecRestoreArrayRead(
                state,
                &x);
            return PETSC_ERR_ARG_OUTOFRANGE;
        }

        const PetscErrorCode y_restore =
            VecRestoreArrayRead(
                search_direction,
                &y);
        const PetscErrorCode x_restore =
            VecRestoreArrayRead(
                state,
                &x);
        if (y_restore != PETSC_SUCCESS) {
            return y_restore;
        }
        if (x_restore != PETSC_SUCCESS) {
            return x_restore;
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
        if (!std::isfinite(global_scale) ||
            !(global_scale > 0.0) ||
            global_scale > 1.0) {
            return PETSC_ERR_ARG_OUTOFRANGE;
        }

        *changed_direction =
            PETSC_FALSE;
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

private:
    SinglePhaseSnesAssemblyContext3D(
        MPI_Comm comm,
        const mpmc::discretization_petsc::
            ParallelOwnedConnectionSchedule3D& schedule,
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
        std::string natural_variable_id,
        double time_step_seconds,
        std::vector<SinglePhaseSnesCellInput3D>
            cell_inputs,
        std::vector<
            SinglePhaseSnesAuthoritativeFaceInput3D>
            face_inputs,
        std::vector<std::optional<std::size_t>>
            cell_input_by_local,
        std::vector<std::optional<std::size_t>>
            face_input_by_local,
        SinglePhaseCurrentCellEvaluatorBinding3D
            cell_evaluator,
        std::size_t q,
        std::size_t component_count)
        : comm_(comm),
          schedule_(&schedule),
          partition_(&partition),
          dof_layout_(&dof_layout),
          dof_numbering_(&dof_numbering),
          cell_bridge_(&cell_bridge),
          cell_pattern_(&cell_pattern),
          natural_variable_id_(
              std::move(natural_variable_id)),
          time_step_seconds_(
              time_step_seconds),
          cell_inputs_(
              std::move(cell_inputs)),
          face_inputs_(
              std::move(face_inputs)),
          cell_input_by_local_(
              std::move(cell_input_by_local)),
          face_input_by_local_(
              std::move(face_input_by_local)),
          cell_evaluator_(cell_evaluator),
          q_(q),
          component_count_(
              component_count) {}

    [[nodiscard]] PetscErrorCode
    initialize_state_exchange() {
        const PetscInt owned_cells =
            cell_bridge_
                ->local_owned_row_count();
        if (owned_cells < 0 ||
            q_ >
                static_cast<std::size_t>(
                    std::numeric_limits<PetscInt>::max()) ||
            owned_cells >
                std::numeric_limits<PetscInt>::max() /
                    static_cast<PetscInt>(q_)) {
            return PETSC_ERR_ARG_OUTOFRANGE;
        }

        const PetscInt q_petsc =
            static_cast<PetscInt>(q_);
        const PetscInt nroots =
            owned_cells * q_petsc;
        const std::size_t local_cells =
            partition_->entity_count(
                mpmc::mesh::EntityKind::cell);
        if (local_cells >
            static_cast<std::size_t>(
                std::numeric_limits<PetscInt>::max()) /
                q_) {
            return PETSC_ERR_ARG_OUTOFRANGE;
        }
        const PetscInt nleaves =
            static_cast<PetscInt>(
                local_cells * q_);

        std::array<PetscInt, 2>
            local_range{
                cell_bridge_->global_row_start(),
                cell_bridge_->global_row_end()};
        std::vector<PetscInt> all_ranges(
            static_cast<std::size_t>(
                partition_->rank_count()) *
            2U);
        if (MPI_Allgather(
                local_range.data(),
                2,
                MPIU_INT,
                all_ranges.data(),
                2,
                MPIU_INT,
                comm_) != MPI_SUCCESS) {
            return PETSC_ERR_MPI;
        }

        std::vector<PetscSFNode> remote(
            static_cast<std::size_t>(
                nleaves));
        for (std::size_t local = 0U;
             local < local_cells;
             ++local) {
            const auto cell =
                mpmc::mesh::LocalIndex{
                    static_cast<
                        mpmc::mesh::LocalIndex::value_type>(
                            local)};
            const auto owner =
                partition_->owner_rank(
                    mpmc::mesh::EntityKind::cell,
                    cell);
            const std::size_t owner_index =
                static_cast<std::size_t>(
                    owner.value());
            if (owner_index >=
                partition_->rank_count()) {
                return PETSC_ERR_ARG_INCOMP;
            }
            const PetscInt owner_start =
                all_ranges[
                    owner_index * 2U];
            const PetscInt owner_end =
                all_ranges[
                    owner_index * 2U + 1U];
            const PetscInt cell_row =
                cell_bridge_->global_row(
                    cell);
            if (cell_row < owner_start ||
                cell_row >= owner_end) {
                return PETSC_ERR_ARG_INCOMP;
            }
            const PetscInt remote_cell =
                cell_row - owner_start;
            for (std::size_t slot = 0U;
                 slot < q_;
                 ++slot) {
                const std::size_t leaf =
                    local * q_ + slot;
                remote[leaf] =
                    PetscSFNode{
                        static_cast<PetscInt>(
                            owner.value()),
                        remote_cell * q_petsc +
                            static_cast<PetscInt>(
                                slot)};
            }
        }

        PetscErrorCode error =
            PetscSFCreate(
                comm_,
                &state_sf_);
        if (error == PETSC_SUCCESS) {
            error =
                PetscSFSetGraph(
                    state_sf_,
                    nroots,
                    nleaves,
                    nullptr,
                    PETSC_COPY_VALUES,
                    remote.data(),
                    PETSC_COPY_VALUES);
        }
        if (error == PETSC_SUCCESS) {
            error =
                PetscSFSetUp(
                    state_sf_);
        }
        if (error == PETSC_SUCCESS) {
            error =
                VecCreateSeq(
                    PETSC_COMM_SELF,
                    nleaves,
                    &local_state_);
        }
        return error;
    }

    [[nodiscard]] PetscErrorCode
    broadcast_state(
        Vec global_state) {
        if (global_state == nullptr ||
            state_sf_ == nullptr ||
            local_state_ == nullptr) {
            return PETSC_ERR_ARG_NULL;
        }

        PetscInt local_size = -1;
        PetscErrorCode error =
            VecGetLocalSize(
                global_state,
                &local_size);
        if (error != PETSC_SUCCESS) {
            return error;
        }
        PetscInt roots = -1;
        PetscInt leaves = -1;
        error =
            PetscSFGetGraph(
                state_sf_,
                &roots,
                &leaves,
                nullptr,
                nullptr);
        if (error != PETSC_SUCCESS ||
            local_size != roots) {
            return error != PETSC_SUCCESS
                ? error
                : PETSC_ERR_ARG_SIZ;
        }

        const PetscScalar* root_array =
            nullptr;
        PetscScalar* leaf_array =
            nullptr;
        error =
            VecGetArrayRead(
                global_state,
                &root_array);
        if (error == PETSC_SUCCESS) {
            error =
                VecGetArray(
                    local_state_,
                    &leaf_array);
        }
        if (error != PETSC_SUCCESS) {
            if (root_array != nullptr) {
                (void)VecRestoreArrayRead(
                    global_state,
                    &root_array);
            }
            return error;
        }

        PetscErrorCode communication =
            PetscSFBcastBegin(
                state_sf_,
                MPIU_SCALAR,
                root_array,
                leaf_array,
                MPI_REPLACE);
        if (communication ==
            PETSC_SUCCESS) {
            communication =
                PetscSFBcastEnd(
                    state_sf_,
                    MPIU_SCALAR,
                    root_array,
                    leaf_array,
                    MPI_REPLACE);
        }

        const PetscErrorCode leaf_restore =
            VecRestoreArray(
                local_state_,
                &leaf_array);
        const PetscErrorCode root_restore =
            VecRestoreArrayRead(
                global_state,
                &root_array);
        if (communication !=
            PETSC_SUCCESS) {
            return communication;
        }
        if (leaf_restore !=
            PETSC_SUCCESS) {
            return leaf_restore;
        }
        return root_restore;
    }

    static void validate_current_cell(
        const SinglePhaseSnesCellInput3D& input,
        const SinglePhaseCurrentCellLinearization3D&
            current) {
        using namespace
            single_phase_snes_assembly_detail;

        if (!same_layout(
                input.frozen_layout,
                current.state.layout()) ||
            current.state.component_ids().size() !=
                input.component_ids.size() ||
            !std::equal(
                current.state.component_ids().begin(),
                current.state.component_ids().end(),
                input.component_ids.begin()) ||
            !mpmc::flow::single_phase_detail::
                same_layout(
                    current.molar_density.layout,
                    input.frozen_layout
                        .descriptor())) {
            throw std::invalid_argument(
                "single-phase current-cell evaluator changed frozen chart/component identity");
        }

        (void)mpmc::flow::
            build_single_phase_component_accumulation_linearization(
                current.state,
                input.porosity,
                current.molar_density);
        (void)mpmc::flow::
            build_single_phase_energy_accumulation_linearization(
                current.state,
                input.porosity,
                current.transport,
                current.caloric,
                current.rock);
        (void)mpmc::flow::
            build_single_phase_mobility_linearization(
                current.state,
                current.transport);
    }

    static PetscErrorCode
    snes_function(
        Vec state,
        Vec residual,
        void* raw_context,
        NaturalVariableSnesEvaluationStatus3D*
            status) {
        if (raw_context == nullptr ||
            residual == nullptr ||
            status == nullptr) {
            return PETSC_ERR_ARG_NULL;
        }
        auto* context =
            static_cast<
                SinglePhaseSnesAssemblyContext3D*>(
                    raw_context);
        std::optional<
            CompleteNaturalVariableAssemblySnapshot3D>
            assembly;
        const PetscErrorCode error =
            context->evaluate_complete_assembly(
                state,
                &assembly,
                status);
        if (error != PETSC_SUCCESS ||
            *status ==
                NaturalVariableSnesEvaluationStatus3D::
                    domain_error) {
            return error;
        }
        if (!assembly.has_value()) {
            return PETSC_ERR_PLIB;
        }

        for (const auto& entry :
             assembly->residual_entries()) {
            const PetscInt row =
                entry.petsc_global_row;
            const PetscScalar value =
                static_cast<PetscScalar>(
                    entry.native_value);
            const PetscErrorCode set_error =
                VecSetValues(
                    residual,
                    1,
                    &row,
                    &value,
                    INSERT_VALUES);
            if (set_error != PETSC_SUCCESS) {
                return set_error;
            }
        }
        return PETSC_SUCCESS;
    }

    static PetscErrorCode
    snes_jacobian(
        Vec state,
        Mat jacobian,
        void* raw_context,
        NaturalVariableSnesEvaluationStatus3D*
            status) {
        if (raw_context == nullptr ||
            jacobian == nullptr ||
            status == nullptr) {
            return PETSC_ERR_ARG_NULL;
        }
        auto* context =
            static_cast<
                SinglePhaseSnesAssemblyContext3D*>(
                    raw_context);
        std::optional<
            CompleteNaturalVariableAssemblySnapshot3D>
            assembly;
        const PetscErrorCode error =
            context->evaluate_complete_assembly(
                state,
                &assembly,
                status);
        if (error != PETSC_SUCCESS ||
            *status ==
                NaturalVariableSnesEvaluationStatus3D::
                    domain_error) {
            return error;
        }
        if (!assembly.has_value()) {
            return PETSC_ERR_PLIB;
        }

        for (const auto& entry :
             assembly->jacobian_entries()) {
            const PetscInt row =
                entry.petsc_global_row;
            const PetscInt column =
                entry.petsc_global_column;
            const PetscScalar value =
                static_cast<PetscScalar>(
                    entry.value);
            const PetscErrorCode set_error =
                MatSetValues(
                    jacobian,
                    1,
                    &row,
                    1,
                    &column,
                    &value,
                    INSERT_VALUES);
            if (set_error != PETSC_SUCCESS) {
                return set_error;
            }
        }
        return PETSC_SUCCESS;
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
            SinglePhaseSnesAssemblyContext3D*>(
                raw_context)
            ->positive_support_precheck(
                state,
                search_direction,
                changed_direction);
    }

    MPI_Comm comm_{};
    const mpmc::discretization_petsc::
        ParallelOwnedConnectionSchedule3D*
            schedule_{};
    const mpmc::mesh::PartitionSnapshot*
        partition_{};
    const mpmc::mesh::DofLayout*
        dof_layout_{};
    const mpmc::mesh::DofNumberingSnapshot*
        dof_numbering_{};
    const mpmc::discretization_petsc::
        PetscMpiAijSymbolicPreallocation3D*
            cell_bridge_{};
    const mpmc::discretization_petsc::
        OwnedCellStructuralColumnPatternSnapshot3D*
            cell_pattern_{};

    std::string natural_variable_id_;
    double time_step_seconds_{};
    std::vector<SinglePhaseSnesCellInput3D>
        cell_inputs_;
    std::vector<
        SinglePhaseSnesAuthoritativeFaceInput3D>
        face_inputs_;
    std::vector<std::optional<std::size_t>>
        cell_input_by_local_;
    std::vector<std::optional<std::size_t>>
        face_input_by_local_;
    SinglePhaseCurrentCellEvaluatorBinding3D
        cell_evaluator_;
    std::size_t q_{};
    std::size_t component_count_{};

    PetscSF state_sf_{};
    Vec local_state_{};
};

} // namespace mpmc::flow_discretization_petsc

#endif // MPMC_FLOW_DISCRETIZATION_PETSC_SINGLE_PHASE_SNES_ASSEMBLY_HPP
