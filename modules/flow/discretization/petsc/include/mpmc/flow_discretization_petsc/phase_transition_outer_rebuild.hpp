#ifndef MPMC_FLOW_DISCRETIZATION_PETSC_PHASE_TRANSITION_OUTER_REBUILD_HPP
#define MPMC_FLOW_DISCRETIZATION_PETSC_PHASE_TRANSITION_OUTER_REBUILD_HPP

#include <mpmc/flow/single_phase_natural_variable.hpp>
#include <mpmc/flow/two_phase_natural_variable.hpp>
#include <mpmc/flow_discretization_petsc/frozen_absent_phase_coordinate_registry.hpp>

#include <petscmat.h>
#include <petscvec.h>

#include <array>
#include <cmath>
#include <cstddef>
#include <memory>
#include <optional>
#include <span>
#include <stdexcept>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

namespace mpmc::flow_discretization_petsc {

inline constexpr std::string_view
    phase_transition_outer_rebuild_convention =
        "flow_discretization_petsc/phase-transition-outer-rebuild/v1";

struct FrozenPhaseTransitionRebuildCell3D {
    mpmc::mesh::LocalIndex cell{
        mpmc::mesh::LocalIndex::value_type{0}};
    mpmc::mesh::GlobalEntityId cell_global{
        mpmc::mesh::GlobalEntityId::value_type{0}};
    double bulk_volume_m3{};
    double porosity{};
    std::vector<std::string> component_ids;

    mpmc::flow::NaturalVariableLayoutDescriptor
        target_layout{
            2U,
            1U,
            {1U}};
    std::vector<double>
        target_natural_variables;
    mpmc::flow::FrozenActivePhaseIdentityMap
        target_active_phases{
            {
                mpmc::flow::FrozenPhysicalPhaseIdentity{
                    "invalid",
                    "invalid"}}};

    std::vector<
        FrozenAbsentPhaseCoordinateEntry3D>
        frozen_absent_phases;

    std::optional<
        mpmc::flow::
            PoreVolumeComponentAccumulationSnapshot3P>
        previous_component_accumulation;
    std::optional<
        mpmc::flow::
            PoreVolumeEnergyAccumulationSnapshot3P>
        previous_energy_accumulation;

    std::string transition_evidence_profile;
};

[[nodiscard]] inline
FrozenPhaseTransitionRebuildCell3D
make_accepted_phase_transition_rebuild_cell_3d(
    mpmc::mesh::LocalIndex cell,
    mpmc::mesh::GlobalEntityId cell_global,
    double bulk_volume_m3,
    double porosity,
    const mpmc::flow::
        PhaseSetTransitionCandidate&
            candidate,
    const mpmc::flow::
        PoreVolumeComponentAccumulationSnapshot3P&
            current_component_inventory,
    const mpmc::flow::
        PoreVolumeComponentAccumulationSnapshot3P&
            previous_component_accumulation,
    const mpmc::flow::
        PoreVolumeEnergyAccumulationSnapshot3P&
            previous_energy_accumulation,
    mpmc::flow::FrozenActivePhaseIdentityMap
        source_active_phases,
    mpmc::flow::FrozenActivePhaseIdentityMap
        target_active_phases,
    std::vector<
        FrozenAbsentPhaseCoordinateEntry3D>
        frozen_absent_phases,
    mpmc::flow::
        PhaseSetTransitionProjectionOptions
            projection_options = {}) {
    const auto projection =
        mpmc::flow::
            project_phase_set_transition_candidate(
                candidate,
                current_component_inventory,
                projection_options);
    const auto continuation =
        mpmc::flow::
            make_phase_identity_continuation_snapshot(
                candidate,
                std::move(
                    source_active_phases),
                target_active_phases);
    const auto history =
        mpmc::flow::
            migrate_phase_set_transition_history(
                candidate.source_phase_count,
                candidate.target_phase_count,
                candidate.component_ids,
                porosity,
                previous_component_accumulation,
                previous_energy_accumulation);

    if (projection.layout().phase_count() !=
            target_active_phases
                .phase_count() ||
        candidate.component_ids !=
            history
                .previous_component_accumulation
                .component_ids) {
        throw std::logic_error(
            "mpmc::flow_discretization_petsc: accepted transition projection/identity/history disagree");
    }

    return {
        cell,
        cell_global,
        bulk_volume_m3,
        porosity,
        candidate.component_ids,
        projection.layout(),
        std::vector<double>{
            projection
                .natural_variables()
                .begin(),
            projection
                .natural_variables()
                .end()},
        std::move(
            target_active_phases),
        std::move(
            frozen_absent_phases),
        history
            .previous_component_accumulation,
        history
            .previous_energy_accumulation,
        std::string{
            continuation
                .evidence_profile()}};
}

namespace phase_transition_outer_rebuild_detail {

[[nodiscard]] inline PetscErrorCode
collective_error(
    MPI_Comm comm,
    PetscErrorCode local_error) {
    int local =
        static_cast<int>(
            local_error);
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

inline void validate_cell_snapshot(
    const FrozenPhaseTransitionRebuildCell3D&
        snapshot,
    std::size_t expected_local,
    const mpmc::mesh::PartitionSnapshot&
        partition,
    std::size_t component_count) {
    if (static_cast<std::size_t>(
            snapshot.cell.value()) !=
            expected_local ||
        snapshot.cell_global !=
            partition.global_id(
                mpmc::mesh::EntityKind::cell,
                snapshot.cell) ||
        !std::isfinite(
            snapshot.bulk_volume_m3) ||
        !(snapshot.bulk_volume_m3 > 0.0) ||
        !std::isfinite(
            snapshot.porosity) ||
        !(snapshot.porosity > 0.0) ||
        !(snapshot.porosity < 1.0) ||
        snapshot.component_ids.size() !=
            component_count ||
        snapshot.target_layout
                .component_count() !=
            component_count ||
        snapshot.target_layout
                .phase_count() !=
            snapshot.target_active_phases
                .phase_count() ||
        (partition.is_owned(
             mpmc::mesh::EntityKind::cell,
             snapshot.cell)
             ? snapshot.target_natural_variables
                       .size() !=
                   snapshot.target_layout
                       .unknown_count()
             : !snapshot.target_natural_variables
                       .empty() &&
                   snapshot.target_natural_variables
                           .size() !=
                       snapshot.target_layout
                           .unknown_count()) ||
        snapshot.transition_evidence_profile
            .empty()) {
        throw std::invalid_argument(
            "mpmc::flow_discretization_petsc: malformed frozen phase-transition rebuild cell");
    }

    for (std::size_t index = 0U;
         index <
             snapshot.component_ids.size();
         ++index) {
        if (snapshot.component_ids[index]
                .empty()) {
            throw std::invalid_argument(
                "mpmc::flow_discretization_petsc: rebuild component id is empty");
        }
        if (index > 0U &&
            std::find(
                snapshot.component_ids.begin(),
                snapshot.component_ids.begin() +
                    static_cast<
                        std::ptrdiff_t>(
                            index),
                snapshot.component_ids[index]) !=
                snapshot.component_ids.begin() +
                    static_cast<
                        std::ptrdiff_t>(
                            index)) {
            throw std::invalid_argument(
                "mpmc::flow_discretization_petsc: rebuild component ids are not unique");
        }
    }
    for (double value :
         snapshot.target_natural_variables) {
        if (!std::isfinite(value)) {
            throw std::invalid_argument(
                "mpmc::flow_discretization_petsc: rebuild natural-variable state contains non-finite value");
        }
    }

    for (const auto& absent :
         snapshot.frozen_absent_phases) {
        absent.reference_coordinates.validate();
        const auto& host =
            absent.reference_coordinates
                .host_state_identity;
        if (host.component_ids !=
                snapshot.component_ids ||
            host.layout.component_count() !=
                snapshot.target_layout
                    .component_count() ||
            host.layout.phase_count() !=
                snapshot.target_layout
                    .phase_count() ||
            host.layout.unknown_count() !=
                snapshot.target_layout
                    .unknown_count()) {
            throw std::invalid_argument(
                "mpmc::flow_discretization_petsc: frozen absent-phase chart does not match rebuilt target layout");
        }
    }
}

[[nodiscard]] inline
mpmc::flow::NaturalVariableLayout1P
layout_1p(
    const mpmc::flow::
        NaturalVariableLayoutDescriptor&
            descriptor) {
    if (descriptor.phase_count() != 1U) {
        throw std::invalid_argument(
            "mpmc::flow_discretization_petsc: requested 1P layout from non-1P descriptor");
    }
    return mpmc::flow::
        NaturalVariableLayout1P{
            mpmc::flow::
                NaturalVariableCompositionPivot1P::
                    from_dependent_component(
                        descriptor
                            .component_count(),
                        descriptor
                            .dependent_composition_component(
                                mpmc::flow::
                                    PhaseSlot3::
                                        phase0))};
}

[[nodiscard]] inline
mpmc::flow::NaturalVariableLayout2P
layout_2p(
    const mpmc::flow::
        NaturalVariableLayoutDescriptor&
            descriptor) {
    if (descriptor.phase_count() != 2U) {
        throw std::invalid_argument(
            "mpmc::flow_discretization_petsc: requested 2P layout from non-2P descriptor");
    }
    return mpmc::flow::
        NaturalVariableLayout2P{
            mpmc::flow::
                NaturalVariableCompositionPivot2P::
                    from_dependent_components(
                        descriptor
                            .component_count(),
                        std::array<std::size_t, 2>{
                            descriptor
                                .dependent_composition_component(
                                    mpmc::flow::
                                        PhaseSlot3::
                                            phase0),
                            descriptor
                                .dependent_composition_component(
                                    mpmc::flow::
                                        PhaseSlot3::
                                            phase1)})};
}

[[nodiscard]] inline
MixedCardinalityPhysicalSnesCellInput3D
make_cell_input(
    const FrozenPhaseTransitionRebuildCell3D&
        snapshot,
    bool owned) {
    const auto previous_component =
        owned
            ? snapshot
                  .previous_component_accumulation
            : std::nullopt;
    const auto previous_energy =
        owned
            ? snapshot
                  .previous_energy_accumulation
            : std::nullopt;
    if (owned &&
        (!previous_component.has_value() ||
         !previous_energy.has_value())) {
        throw std::invalid_argument(
            "mpmc::flow_discretization_petsc: owned rebuilt cell lacks conservation history");
    }

    switch (
        snapshot.target_layout
            .phase_count()) {
    case 1U: {
        SinglePhaseSnesCellInput3D
            input;
        input.cell = snapshot.cell;
        input.cell_global =
            snapshot.cell_global;
        input.bulk_volume_m3 =
            snapshot.bulk_volume_m3;
        input.porosity =
            snapshot.porosity;
        input.frozen_layout =
            layout_1p(
                snapshot.target_layout);
        input.component_ids =
            snapshot.component_ids;
        input.previous_component_accumulation =
            previous_component;
        input.previous_energy_accumulation =
            previous_energy;
        return input;
    }
    case 2U: {
        TwoPhaseSnesCellInput3D
            input;
        input.cell = snapshot.cell;
        input.cell_global =
            snapshot.cell_global;
        input.bulk_volume_m3 =
            snapshot.bulk_volume_m3;
        input.porosity =
            snapshot.porosity;
        input.frozen_layout =
            layout_2p(
                snapshot.target_layout);
        input.component_ids =
            snapshot.component_ids;
        input.previous_component_accumulation =
            previous_component;
        input.previous_energy_accumulation =
            previous_energy;
        return input;
    }
    case 3U: {
        FixedThreePhaseSnesCellInput3D
            input;
        input.cell = snapshot.cell;
        input.cell_global =
            snapshot.cell_global;
        input.bulk_volume_m3 =
            snapshot.bulk_volume_m3;
        input.porosity =
            snapshot.porosity;
        input.frozen_layout =
            static_cast<
                mpmc::flow::
                    NaturalVariableLayout3P>(
                        snapshot.target_layout);
        input.component_ids =
            snapshot.component_ids;
        input.previous_component_accumulation =
            previous_component;
        input.previous_energy_accumulation =
            previous_energy;
        return input;
    }
    default:
        throw std::invalid_argument(
            "mpmc::flow_discretization_petsc: rebuilt cell phase count must be 1/2/3");
    }
}

} // namespace phase_transition_outer_rebuild_detail

class PhaseTransitionRebuiltNaturalVariableSystem3D {
public:
    PhaseTransitionRebuiltNaturalVariableSystem3D(
        const PhaseTransitionRebuiltNaturalVariableSystem3D&) =
        delete;
    PhaseTransitionRebuiltNaturalVariableSystem3D& operator=(
        const PhaseTransitionRebuiltNaturalVariableSystem3D&) =
        delete;
    PhaseTransitionRebuiltNaturalVariableSystem3D(
        PhaseTransitionRebuiltNaturalVariableSystem3D&&) =
        delete;
    PhaseTransitionRebuiltNaturalVariableSystem3D& operator=(
        PhaseTransitionRebuiltNaturalVariableSystem3D&&) =
        delete;

    ~PhaseTransitionRebuiltNaturalVariableSystem3D() {
        if (jacobian_ != nullptr) {
            (void)MatDestroy(
                &jacobian_);
        }
        if (initial_state_ != nullptr) {
            (void)VecDestroy(
                &initial_state_);
        }
    }

    [[nodiscard]] const
    VariableCardinalityNaturalVariableNumbering3D&
    numbering() const noexcept {
        return *numbering_;
    }

    [[nodiscard]] const
    FrozenAbsentPhaseCoordinateRegistry3D&
    coordinate_registry() const noexcept {
        return *coordinate_registry_;
    }

    [[nodiscard]] Vec
    initial_state() const noexcept {
        return initial_state_;
    }

    [[nodiscard]] Mat
    jacobian_structure() const noexcept {
        return jacobian_;
    }

    [[nodiscard]] double
    time_step_seconds() const noexcept {
        return physical_context_ != nullptr
            ? physical_context_
                  ->time_step_seconds()
            : 0.0;
    }

    /// Update only the trial timestep used by the next nonlinear attempt.
    /// This does not advance accepted time or rewrite accepted histories/state.
    [[nodiscard]] PetscErrorCode
    set_trial_timestep_seconds(
        double timestep_seconds) noexcept {
        if (physical_context_ == nullptr) {
            return PETSC_ERR_ARG_WRONGSTATE;
        }
        return physical_context_
            ->set_trial_timestep_seconds(
                timestep_seconds);
    }

    [[nodiscard]] PetscErrorCode
    accepted_history_matches_state(
        Vec accepted_state,
        bool* matches) {
        if (physical_context_ == nullptr) {
            return PETSC_ERR_ARG_WRONGSTATE;
        }
        return physical_context_
            ->accepted_history_matches_state(
                accepted_state,
                matches);
    }

    [[nodiscard]] NaturalVariableSnesEvaluator3D
    snes_evaluator() noexcept {
        return physical_context_
            ->snes_evaluator();
    }

    [[nodiscard]] PetscErrorCode
    evaluate_local_cells_for_phase_transition(
        Vec global_state,
        std::vector<std::optional<
            MixedCardinalityPhysicalCurrentCellLinearization3D>>*
                output,
        std::vector<double>* porosities,
        NaturalVariableSnesEvaluationStatus3D*
            status) const {
        if (physical_context_ == nullptr) {
            return PETSC_ERR_ARG_WRONGSTATE;
        }
        return physical_context_
            ->evaluate_local_cells_for_phase_transition(
                global_state,
                output,
                porosities,
                status);
    }

    [[nodiscard]] PetscErrorCode
    evaluate_current_cell_for_phase_transition(
        mpmc::mesh::LocalIndex cell,
        std::span<const double> natural_variables,
        std::optional<
            MixedCardinalityPhysicalCurrentCellLinearization3D>*
                output,
        double* porosity,
        NaturalVariableSnesEvaluationStatus3D*
            status) const {
        if (physical_context_ == nullptr) {
            return PETSC_ERR_ARG_WRONGSTATE;
        }
        return physical_context_
            ->evaluate_phase_transition_cell(
                cell,
                natural_variables,
                output,
                porosity,
                status);
    }

    /// Commit one already accepted state as the exact history/state baseline
    /// for the next physical timestep. Cardinality and physical phase identity
    /// remain frozen; a topology change must still use the outer transition
    /// rebuild path. Absent-phase coordinate charts are re-anchored at the
    /// accepted host state while preserving their selected branch provenance.
    [[nodiscard]] PetscErrorCode
    rebase_accepted_timestep(
        Vec accepted_state,
        double next_time_step_seconds) {
        using namespace
            phase_transition_outer_rebuild_detail;

        if (accepted_state == nullptr ||
            physical_context_ == nullptr ||
            coordinate_registry_ == nullptr ||
            initial_state_ == nullptr ||
            !std::isfinite(
                next_time_step_seconds) ||
            !(next_time_step_seconds > 0.0)) {
            return PETSC_ERR_ARG_INCOMP;
        }

        std::vector<std::optional<
            MixedCardinalityPhysicalCurrentCellLinearization3D>>
            current;
        std::vector<double> porosities;
        NaturalVariableSnesEvaluationStatus3D
            status =
                NaturalVariableSnesEvaluationStatus3D::
                    success;
        PetscErrorCode error =
            evaluate_local_cells_for_phase_transition(
                accepted_state,
                &current,
                &porosities,
                &status);
        if (error != PETSC_SUCCESS) {
            return error;
        }
        if (status !=
                NaturalVariableSnesEvaluationStatus3D::
                    success ||
            current.size() !=
                numbering_->local_cell_count() ||
            porosities.size() !=
                current.size()) {
            return PETSC_ERR_ARG_WRONGSTATE;
        }

        PetscErrorCode local_error =
            PETSC_SUCCESS;
        std::optional<
            FrozenAbsentPhaseCoordinateRegistry3D>
            refreshed_registry;
        std::vector<
            MixedCardinalityPhysicalSnesCellInput3D>
            rebased_inputs;
        try {
            const auto old_cells =
                coordinate_registry_->cells();
            if (old_cells.size() !=
                current.size()) {
                throw std::invalid_argument(
                    "coordinate registry/cell count mismatch");
            }

            std::vector<
                FrozenAbsentPhaseCoordinateCell3D>
                refreshed_cells;
            refreshed_cells.reserve(
                old_cells.size());

            for (std::size_t local = 0U;
                 local < current.size();
                 ++local) {
                if (!current[local]
                        .has_value()) {
                    throw std::invalid_argument(
                        "accepted host cell is not evaluable");
                }
                const auto cell =
                    mpmc::mesh::LocalIndex{
                        static_cast<
                            mpmc::mesh::
                                LocalIndex::value_type>(
                                    local)};
                const auto& record =
                    numbering_->cell(
                        cell);
                const auto& old =
                    coordinate_registry_->cell(
                        record.cell_global);
                const auto host =
                    std::visit(
                        [](const auto& typed) {
                            return typed
                                .transport
                                .state_identity;
                        },
                        *current[local]);
                if (old.cell !=
                        cell ||
                    old.cell_global !=
                        record.cell_global ||
                    old.active_phases
                            .phase_count() !=
                        host.layout
                            .phase_count()) {
                    throw std::invalid_argument(
                        "accepted host topology changed without rebuild");
                }

                std::vector<
                    FrozenAbsentPhaseCoordinateEntry3D>
                    absent;
                absent.reserve(
                    old.absent_phases.size());
                for (const auto& entry :
                     old.absent_phases) {
                    auto coordinates =
                        coordinate_registry_
                            ->resolve_affine(
                                old.cell_global,
                                entry.identity,
                                host);
                    coordinates.validate();
                    absent.push_back(
                        {
                            entry.identity,
                            std::move(
                                coordinates),
                            entry
                                .selected_branch_provenance});
                }
                refreshed_cells.push_back(
                    {
                        old.cell,
                        old.cell_global,
                        old.active_phases,
                        std::move(absent)});
            }

            std::vector<
                FrozenAbsentPhaseFaceEndpoint3D>
                refreshed_faces{
                    coordinate_registry_
                        ->faces()
                        .begin(),
                    coordinate_registry_
                        ->faces()
                        .end()};
            refreshed_registry.emplace(
                std::move(
                    refreshed_cells),
                std::move(
                    refreshed_faces));

            local_error =
                physical_context_
                    ->prepare_accepted_history_rebase(
                        current,
                        next_time_step_seconds,
                        &rebased_inputs);
        } catch (const std::exception&) {
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
        if (!refreshed_registry.has_value() ||
            rebased_inputs.size() !=
                current.size()) {
            return PETSC_ERR_PLIB;
        }

        Vec next_initial = nullptr;
        error =
            VecDuplicate(
                initial_state_,
                &next_initial);
        if (error == PETSC_SUCCESS) {
            error =
                VecCopy(
                    accepted_state,
                    next_initial);
        }
        if (error != PETSC_SUCCESS) {
            if (next_initial != nullptr) {
                (void)VecDestroy(
                    &next_initial);
            }
            return error;
        }

        physical_context_
            ->commit_accepted_history_rebase(
                std::move(
                    rebased_inputs),
                next_time_step_seconds);
        *coordinate_registry_ =
            std::move(
                *refreshed_registry);

        Vec old_initial =
            initial_state_;
        initial_state_ =
            next_initial;
        next_initial = nullptr;
        return VecDestroy(
            &old_initial);
    }

    [[nodiscard]] PetscErrorCode
    solve(
        Vec* solution,
        std::optional<
            VariableCardinalityNaturalVariableSnesSolveReport3D>*
                report,
        std::optional<
            NaturalVariableSnesFailureDiagnostics3D>*
                failure_diagnostics = nullptr) {
        Vec row_scaling =
            nullptr;
        PetscErrorCode error =
            make_variable_cardinality_initial_row_equilibration_3d(
                comm_,
                *numbering_,
                initial_state_,
                jacobian_,
                physical_context_
                    ->snes_evaluator(),
                &row_scaling);
        if (error != PETSC_SUCCESS) {
            return error;
        }

        error =
            solve_variable_cardinality_natural_variable_snes_3d(
                comm_,
                *numbering_,
                initial_state_,
                jacobian_,
                physical_context_
                    ->snes_evaluator(),
                solution,
                report,
                row_scaling,
                failure_diagnostics);
        const PetscErrorCode destroy =
            VecDestroy(
                &row_scaling);
        return error != PETSC_SUCCESS
            ? error
            : destroy;
    }

private:
    PhaseTransitionRebuiltNaturalVariableSystem3D(
        MPI_Comm comm)
        : comm_(comm) {}

    friend PetscErrorCode
    rebuild_phase_transition_natural_variable_system_3d(
        MPI_Comm,
        const mpmc::discretization_petsc::
            ParallelOwnedConnectionSchedule3D&,
        const mpmc::mesh::
            PartitionSnapshot&,
        const mpmc::discretization_petsc::
            PetscMpiAijSymbolicPreallocation3D&,
        const mpmc::discretization_petsc::
            OwnedCellStructuralColumnPatternSnapshot3D&,
        double,
        std::vector<
            FrozenPhaseTransitionRebuildCell3D>,
        std::vector<
            MixedCardinalityPhysicalSnesAuthoritativeFaceInput3D>,
        MixedCardinalityPhysicalCellEvaluatorBindings3D,
        MixedCardinalityPhysicalCellSourceEvaluatorBinding3D,
        AbsentPhaseThermodynamicProviderEvaluator3D,
        void*,
        std::unique_ptr<
            PhaseTransitionRebuiltNaturalVariableSystem3D>*);

    MPI_Comm comm_;
    std::unique_ptr<
        VariableCardinalityNaturalVariableNumbering3D>
        numbering_;
    std::unique_ptr<
        FrozenAbsentPhaseCoordinateRegistry3D>
        coordinate_registry_;
    std::unique_ptr<
        ThermodynamicAbsentPhaseExtensionAdapterBinding3D>
        thermodynamic_adapter_;
    std::unique_ptr<
        CrossCardinalityTpfaBridgeBinding3D>
        bridge_binding_;
    std::unique_ptr<
        MixedCardinalityPhysicalSnesAssemblyContext3D>
        physical_context_;
    Vec initial_state_{};
    Mat jacobian_{};
};

inline PetscErrorCode
rebuild_phase_transition_natural_variable_system_3d(
    MPI_Comm comm,
    const mpmc::discretization_petsc::
        ParallelOwnedConnectionSchedule3D&
            schedule,
    const mpmc::mesh::
        PartitionSnapshot&
            partition,
    const mpmc::discretization_petsc::
        PetscMpiAijSymbolicPreallocation3D&
            cell_bridge,
    const mpmc::discretization_petsc::
        OwnedCellStructuralColumnPatternSnapshot3D&
            cell_pattern,
    double time_step_seconds,
    std::vector<
        FrozenPhaseTransitionRebuildCell3D>
        cells,
    std::vector<
        MixedCardinalityPhysicalSnesAuthoritativeFaceInput3D>
        face_inputs,
    MixedCardinalityPhysicalCellEvaluatorBindings3D
        cell_evaluators,
    MixedCardinalityPhysicalCellSourceEvaluatorBinding3D
        cell_source_evaluator,
    AbsentPhaseThermodynamicProviderEvaluator3D
        provider_evaluator,
    void* provider_context,
    std::unique_ptr<
        PhaseTransitionRebuiltNaturalVariableSystem3D>*
            output) {
    using namespace
        phase_transition_outer_rebuild_detail;

    if (output == nullptr) {
        return PETSC_ERR_ARG_NULL;
    }
    output->reset();

    const std::size_t local_cell_count =
        partition.entity_count(
            mpmc::mesh::EntityKind::cell);
    PetscErrorCode local_error =
        provider_evaluator == nullptr ||
                provider_context == nullptr ||
                cells.size() !=
                    local_cell_count ||
                face_inputs.size() !=
                    schedule
                        .assembly_rows()
                        .size() ||
                cells.empty()
            ? PETSC_ERR_ARG_INCOMP
            : PETSC_SUCCESS;

    std::size_t component_count = 0U;
    std::vector<std::size_t>
        phase_counts;
    phase_counts.reserve(
        cells.size());

    if (local_error ==
        PETSC_SUCCESS) {
        try {
            component_count =
                cells.front()
                    .component_ids
                    .size();
            if (component_count < 2U) {
                throw std::invalid_argument(
                    "component count");
            }
            const auto canonical_ids =
                cells.front()
                    .component_ids;
            for (std::size_t local = 0U;
                 local < cells.size();
                 ++local) {
                validate_cell_snapshot(
                    cells[local],
                    local,
                    partition,
                    component_count);
                if (cells[local]
                        .component_ids !=
                    canonical_ids) {
                    throw std::invalid_argument(
                        "component order");
                }
                phase_counts.push_back(
                    cells[local]
                        .target_layout
                        .phase_count());
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

    std::optional<
        VariableCardinalityNaturalVariableNumbering3D>
        numbering;
    error =
        make_variable_cardinality_natural_variable_numbering_3d(
            comm,
            partition,
            component_count,
            phase_counts,
            &numbering);
    if (error != PETSC_SUCCESS ||
        !numbering.has_value()) {
        return error != PETSC_SUCCESS
            ? error
            : PETSC_ERR_PLIB;
    }

    auto system =
        std::unique_ptr<
            PhaseTransitionRebuiltNaturalVariableSystem3D>(
                new PhaseTransitionRebuiltNaturalVariableSystem3D{
                    comm});
    system->numbering_ =
        std::make_unique<
            VariableCardinalityNaturalVariableNumbering3D>(
                std::move(*numbering));

    error =
        create_variable_cardinality_natural_variable_vec_3d(
            comm,
            *system->numbering_,
            &system->initial_state_);
    if (error != PETSC_SUCCESS) {
        return error;
    }

    local_error =
        PETSC_SUCCESS;
    try {
        for (const auto& record :
             system->numbering_->cells()) {
            if (record.owner_rank !=
                system->numbering_
                    ->local_rank()) {
                continue;
            }
            const auto& snapshot =
                cells.at(
                    static_cast<std::size_t>(
                        record.cell.value()));
            if (snapshot
                    .target_natural_variables
                    .size() !=
                record.scalar_count) {
                throw std::invalid_argument(
                    "rebuilt state width");
            }
            std::vector<PetscInt>
                indices(
                    record.scalar_count);
            std::vector<PetscScalar>
                values(
                    record.scalar_count);
            for (std::size_t slot = 0U;
                 slot <
                    record.scalar_count;
                 ++slot) {
                indices[slot] =
                    record
                        .petsc_global_scalar_start +
                    static_cast<PetscInt>(
                        slot);
                values[slot] =
                    static_cast<PetscScalar>(
                        snapshot
                            .target_natural_variables[
                                slot]);
            }
            const PetscErrorCode insert_error =
                VecSetValues(
                    system->initial_state_,
                    static_cast<PetscInt>(
                        indices.size()),
                    indices.data(),
                    values.data(),
                    INSERT_VALUES);
            if (insert_error !=
                PETSC_SUCCESS) {
                local_error =
                    insert_error;
                break;
            }
        }
    } catch (...) {
        local_error =
            PETSC_ERR_ARG_INCOMP;
    }
    error =
        collective_error(
            comm,
            local_error);
    if (error == PETSC_SUCCESS) {
        error =
            VecAssemblyBegin(
                system->initial_state_);
    }
    if (error == PETSC_SUCCESS) {
        error =
            VecAssemblyEnd(
                system->initial_state_);
    }
    if (error != PETSC_SUCCESS) {
        return error;
    }

    std::vector<
        MixedCardinalityPhysicalSnesCellInput3D>
        cell_inputs;
    std::vector<
        mpmc::flow::
            FrozenActivePhaseIdentityMap>
        phase_identity_maps;
    std::vector<
        FrozenAbsentPhaseCoordinateCell3D>
        registry_cells;
    cell_inputs.reserve(
        cells.size());
    phase_identity_maps.reserve(
        cells.size());
    registry_cells.reserve(
        cells.size());

    local_error =
        PETSC_SUCCESS;
    try {
        for (const auto& snapshot :
             cells) {
            const bool owned =
                partition.is_owned(
                    mpmc::mesh::
                        EntityKind::cell,
                    snapshot.cell);
            cell_inputs.push_back(
                make_cell_input(
                    snapshot,
                    owned));
            phase_identity_maps.push_back(
                snapshot
                    .target_active_phases);
            registry_cells.push_back(
                {
                    snapshot.cell,
                    snapshot.cell_global,
                    snapshot
                        .target_active_phases,
                    snapshot
                        .frozen_absent_phases});
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

    std::vector<
        FrozenAbsentPhaseFaceEndpoint3D>
        registry_faces;
    registry_faces.reserve(
        schedule
            .assembly_rows()
            .size());
    for (const auto& row :
         schedule.assembly_rows()) {
        registry_faces.push_back(
            {
                row.face_global,
                row.owner_cell_global,
                row.neighbour_cell_global});
    }

    local_error =
        PETSC_SUCCESS;
    try {
        system->coordinate_registry_ =
            std::make_unique<
                FrozenAbsentPhaseCoordinateRegistry3D>(
                    std::move(
                        registry_cells),
                    std::move(
                        registry_faces));
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

    system->thermodynamic_adapter_ =
        std::make_unique<
            ThermodynamicAbsentPhaseExtensionAdapterBinding3D>(
                ThermodynamicAbsentPhaseExtensionAdapterBinding3D{
                    &resolve_frozen_absent_phase_thermodynamic_coordinates_3d,
                    system
                        ->coordinate_registry_
                        .get(),
                    provider_evaluator,
                    provider_context});
    system->bridge_binding_ =
        std::make_unique<
            CrossCardinalityTpfaBridgeBinding3D>(
                CrossCardinalityTpfaBridgeBinding3D{
                    &evaluate_thermodynamic_absent_phase_extension_3d,
                    system
                        ->thermodynamic_adapter_
                        .get()});

    std::optional<
        MixedCardinalityPhysicalSnesAssemblyContext3D>
        context;
    error =
        MixedCardinalityPhysicalSnesAssemblyContext3D::
            create(
                comm,
                schedule,
                partition,
                *system->numbering_,
                cell_bridge,
                cell_pattern,
                time_step_seconds,
                std::move(
                    cell_inputs),
                std::move(
                    phase_identity_maps),
                std::move(
                    face_inputs),
                cell_evaluators,
                {
                    &evaluate_standard_cross_cardinality_tpfa_face_3d,
                    system
                        ->bridge_binding_
                        .get()},
                cell_source_evaluator,
                &context);
    if (error != PETSC_SUCCESS ||
        !context.has_value()) {
        return error != PETSC_SUCCESS
            ? error
            : PETSC_ERR_PLIB;
    }
    system->physical_context_ =
        std::make_unique<
            MixedCardinalityPhysicalSnesAssemblyContext3D>(
                std::move(*context));

    error =
        system->physical_context_
            ->create_jacobian_structure(
                &system->jacobian_);
    if (error != PETSC_SUCCESS) {
        return error;
    }

    *output =
        std::move(system);
    return PETSC_SUCCESS;
}

inline PetscErrorCode
rebuild_phase_transition_natural_variable_system_3d(
    MPI_Comm comm,
    const mpmc::discretization_petsc::
        ParallelOwnedConnectionSchedule3D& schedule,
    const mpmc::mesh::PartitionSnapshot& partition,
    const mpmc::discretization_petsc::
        PetscMpiAijSymbolicPreallocation3D& cell_bridge,
    const mpmc::discretization_petsc::
        OwnedCellStructuralColumnPatternSnapshot3D& cell_pattern,
    double time_step_seconds,
    std::vector<FrozenPhaseTransitionRebuildCell3D> cells,
    std::vector<
        MixedCardinalityPhysicalSnesAuthoritativeFaceInput3D>
        face_inputs,
    MixedCardinalityPhysicalCellEvaluatorBindings3D
        cell_evaluators,
    AbsentPhaseThermodynamicProviderEvaluator3D
        provider_evaluator,
    void* provider_context,
    std::unique_ptr<
        PhaseTransitionRebuiltNaturalVariableSystem3D>*
            output) {
    return rebuild_phase_transition_natural_variable_system_3d(
        comm,
        schedule,
        partition,
        cell_bridge,
        cell_pattern,
        time_step_seconds,
        std::move(cells),
        std::move(face_inputs),
        cell_evaluators,
        {},
        provider_evaluator,
        provider_context,
        output);
}

} // namespace mpmc::flow_discretization_petsc

#endif // MPMC_FLOW_DISCRETIZATION_PETSC_PHASE_TRANSITION_OUTER_REBUILD_HPP
