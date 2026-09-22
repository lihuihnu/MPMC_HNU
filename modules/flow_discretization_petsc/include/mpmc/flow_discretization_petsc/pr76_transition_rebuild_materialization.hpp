#ifndef MPMC_FLOW_DISCRETIZATION_PETSC_PR76_TRANSITION_REBUILD_MATERIALIZATION_HPP
#define MPMC_FLOW_DISCRETIZATION_PETSC_PR76_TRANSITION_REBUILD_MATERIALIZATION_HPP

#include <mpmc/flow_discretization_petsc/pr76_single_phase_transition_target_rebuild.hpp>

#include <petscmat.h>
#include <petscvec.h>

#include <algorithm>
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
    pr76_transition_rebuild_materialization_convention =
        "flow_discretization_petsc/pr76-transition-rebuild-materialization/v1";

template <typename Binding>
struct CellScopedCurrentCellEvaluatorEntry3D {
    mpmc::mesh::GlobalEntityId cell_global{
        mpmc::mesh::GlobalEntityId::value_type{0}};
    Binding binding;
};

class CellScopedMixedCardinalityEvaluatorDispatcher3D {
public:
    using SingleEntry =
        CellScopedCurrentCellEvaluatorEntry3D<
            SinglePhaseCurrentCellEvaluatorBinding3D>;
    using TwoEntry =
        CellScopedCurrentCellEvaluatorEntry3D<
            TwoPhaseCurrentCellEvaluatorBinding3D>;
    using ThreeEntry =
        CellScopedCurrentCellEvaluatorEntry3D<
            FixedThreePhaseCurrentCellEvaluatorBinding3D>;

    CellScopedMixedCardinalityEvaluatorDispatcher3D(
        std::vector<SingleEntry> single_phase,
        std::vector<TwoEntry> two_phase,
        std::vector<ThreeEntry> three_phase)
        : single_phase_(
              std::move(single_phase)),
          two_phase_(
              std::move(two_phase)),
          three_phase_(
              std::move(three_phase)) {
        validate_and_sort(
            &single_phase_);
        validate_and_sort(
            &two_phase_);
        validate_and_sort(
            &three_phase_);
    }

    [[nodiscard]]
    MixedCardinalityPhysicalCellEvaluatorBindings3D
    bindings() noexcept {
        return {
            {
                &evaluate_single_phase,
                this},
            {
                &evaluate_two_phase,
                this},
            {
                &evaluate_three_phase,
                this}};
    }

private:
    template <class Entry>
    static void validate_and_sort(
        std::vector<Entry>* entries) {
        std::sort(
            entries->begin(),
            entries->end(),
            [](const auto& first,
               const auto& second) {
                return first.cell_global <
                    second.cell_global;
            });
        for (std::size_t index = 0U;
             index < entries->size();
             ++index) {
            const auto& entry =
                (*entries)[index];
            if (entry.binding.evaluator ==
                    nullptr ||
                (index > 0U &&
                 (*entries)[index - 1U]
                         .cell_global ==
                     entry.cell_global)) {
                throw std::invalid_argument(
                    "cell-scoped mixed-cardinality evaluator dispatcher has invalid/duplicate binding");
            }
        }
    }

    template <class Entry>
    [[nodiscard]] static const Entry*
    find(
        const std::vector<Entry>& entries,
        mpmc::mesh::GlobalEntityId
            cell_global) noexcept {
        const auto found =
            std::lower_bound(
                entries.begin(),
                entries.end(),
                cell_global,
                [](const auto& entry,
                   const auto& id) {
                    return entry.cell_global <
                        id;
                });
        return found != entries.end() &&
                found->cell_global ==
                    cell_global
            ? &*found
            : nullptr;
    }

    static PetscErrorCode
    evaluate_single_phase(
        mpmc::mesh::LocalIndex cell,
        mpmc::mesh::GlobalEntityId cell_global,
        std::span<const double> natural_variables,
        const mpmc::flow::NaturalVariableLayout1P&
            layout,
        std::span<const std::string> component_ids,
        void* raw_context,
        std::optional<
            SinglePhaseCurrentCellLinearization3D>* output,
        NaturalVariableSnesEvaluationStatus3D*
            status) {
        if (raw_context == nullptr) {
            return PETSC_ERR_ARG_NULL;
        }
        auto* self =
            static_cast<
                CellScopedMixedCardinalityEvaluatorDispatcher3D*>(
                    raw_context);
        const auto* entry =
            find(
                self->single_phase_,
                cell_global);
        if (entry == nullptr) {
            return PETSC_ERR_ARG_INCOMP;
        }
        return entry->binding.evaluator(
            cell,
            cell_global,
            natural_variables,
            layout,
            component_ids,
            entry->binding.user_context,
            output,
            status);
    }

    static PetscErrorCode
    evaluate_two_phase(
        mpmc::mesh::LocalIndex cell,
        mpmc::mesh::GlobalEntityId cell_global,
        std::span<const double> natural_variables,
        const mpmc::flow::NaturalVariableLayout2P&
            layout,
        std::span<const std::string> component_ids,
        void* raw_context,
        std::optional<
            TwoPhaseCurrentCellLinearization3D>* output,
        NaturalVariableSnesEvaluationStatus3D*
            status) {
        if (raw_context == nullptr) {
            return PETSC_ERR_ARG_NULL;
        }
        auto* self =
            static_cast<
                CellScopedMixedCardinalityEvaluatorDispatcher3D*>(
                    raw_context);
        const auto* entry =
            find(
                self->two_phase_,
                cell_global);
        if (entry == nullptr) {
            return PETSC_ERR_ARG_INCOMP;
        }
        return entry->binding.evaluator(
            cell,
            cell_global,
            natural_variables,
            layout,
            component_ids,
            entry->binding.user_context,
            output,
            status);
    }

    static PetscErrorCode
    evaluate_three_phase(
        mpmc::mesh::LocalIndex cell,
        mpmc::mesh::GlobalEntityId cell_global,
        std::span<const double> natural_variables,
        const mpmc::flow::NaturalVariableLayout3P&
            layout,
        std::span<const std::string> component_ids,
        void* raw_context,
        std::optional<
            FixedThreePhaseCurrentCellLinearization3D>* output,
        NaturalVariableSnesEvaluationStatus3D*
            status) {
        if (raw_context == nullptr) {
            return PETSC_ERR_ARG_NULL;
        }
        auto* self =
            static_cast<
                CellScopedMixedCardinalityEvaluatorDispatcher3D*>(
                    raw_context);
        const auto* entry =
            find(
                self->three_phase_,
                cell_global);
        if (entry == nullptr) {
            return PETSC_ERR_ARG_INCOMP;
        }
        return entry->binding.evaluator(
            cell,
            cell_global,
            natural_variables,
            layout,
            component_ids,
            entry->binding.user_context,
            output,
            status);
    }

    std::vector<SingleEntry>
        single_phase_;
    std::vector<TwoEntry>
        two_phase_;
    std::vector<ThreeEntry>
        three_phase_;
};

using Pr76AbsentPhaseSelectedBranchResolver3D =
    PetscErrorCode (*)(
        mpmc::mesh::GlobalEntityId absent_cell_global,
        const mpmc::flow::FrozenPhysicalPhaseIdentity&
            identity,
        const mpmc::flow::NaturalVariableStateIdentity3P&
            absent_host_state,
        std::span<const double>
            active_reference_composition,
        const mpmc::thermodynamics::Pr76SelectedPhase&
            active_reference_selection,
        void* user_context,
        std::optional<
            mpmc::thermodynamics::Pr76SelectedPhase>*
                absent_selection);

struct Pr76AbsentPhaseSelectedBranchBinding3D {
    Pr76AbsentPhaseSelectedBranchResolver3D
        resolver{};
    void* user_context{};
};

struct Pr76TransitionRebuildMaterializedSystem3D {
    std::unique_ptr<
        CellScopedMixedCardinalityEvaluatorDispatcher3D>
        dispatcher;
    std::unique_ptr<
        mpmc::flow::
            Pr76CellScopedAbsentPhasePotentialExtensionProvider<
                double>>
        absent_provider;
    std::unique_ptr<
        PhaseTransitionRebuiltNaturalVariableSystem3D>
        system;
};

namespace pr76_transition_rebuild_materialization_detail {

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

[[nodiscard]] inline bool same_vector(
    std::span<const double> first,
    std::span<const double> second) {
    if (first.size() != second.size()) {
        return false;
    }
    for (std::size_t index = 0U;
         index < first.size();
         ++index) {
        const double scale =
            std::max(
                {1.0,
                 std::abs(first[index]),
                 std::abs(second[index])});
        if (!std::isfinite(first[index]) ||
            !std::isfinite(second[index]) ||
            std::abs(
                first[index] -
                second[index]) >
                8192.0 *
                    std::numeric_limits<double>::
                        epsilon() *
                    scale) {
            return false;
        }
    }
    return true;
}

[[nodiscard]] inline std::vector<
    MixedCardinalityPhysicalSnesAuthoritativeFaceInput3D>
mixed_faces(
    std::span<
        const SinglePhaseSnesAuthoritativeFaceInput3D>
        source) {
    std::vector<
        MixedCardinalityPhysicalSnesAuthoritativeFaceInput3D>
        result;
    result.reserve(
        source.size());
    for (const auto& face : source) {
        result.push_back(
            {
                face.face,
                face.face_global,
                face.transmissibility,
                face.gravity,
                face.owner_to_neighbour_displacement,
                face.thermal_conductance});
    }
    return result;
}

[[nodiscard]] inline const
Pr76SinglePhaseTransitionResolvedPhase3D*
resolved_phase(
    const Pr76SinglePhaseTransitionTargetRebuildPlan3D&
        plan,
    mpmc::mesh::GlobalEntityId
        cell_global,
    const mpmc::flow::
        FrozenPhysicalPhaseIdentity&
            identity) noexcept {
    const auto cell =
        std::find_if(
            plan.local_transition_cells.begin(),
            plan.local_transition_cells.end(),
            [&](const auto& entry) {
                return entry.cell_global ==
                    cell_global;
            });
    if (cell ==
        plan.local_transition_cells.end()) {
        return nullptr;
    }
    const auto phase =
        std::find_if(
            cell->target_phases.begin(),
            cell->target_phases.end(),
            [&](const auto& entry) {
                return entry.identity ==
                    identity;
            });
    return phase ==
            cell->target_phases.end()
        ? nullptr
        : &*phase;
}

[[nodiscard]] inline std::string
absent_branch_provenance(
    mpmc::mesh::GlobalEntityId
        absent_cell_global,
    const mpmc::flow::
        FrozenPhysicalPhaseIdentity&
            identity,
    std::size_t root_index) {
    return std::string{
               "PR76/PT absent-cell="} +
        std::to_string(
            absent_cell_global.value()) +
        "|phase=" +
        identity.provenance_scope +
        "/" +
        identity.opaque_phase_key +
        "|activity.branch=" +
        std::to_string(
            root_index);
}

inline PetscErrorCode
bootstrap_absent_provider(
    const mpmc::flow::
        AbsentPhaseThermodynamicCoordinateExtension&,
    void* raw_context,
    std::optional<
        mpmc::flow::
            AbsentPhasePotentialExtensionLinearization>*,
    NaturalVariableSnesEvaluationStatus3D*) {
    return raw_context == nullptr
        ? PETSC_ERR_ARG_NULL
        : PETSC_ERR_SUP;
}

} // namespace pr76_transition_rebuild_materialization_detail

[[nodiscard]] inline PetscErrorCode
materialize_pr76_single_phase_transition_target_rebuild_3d(
    MPI_Comm comm,
    const mpmc::discretization_petsc::
        ParallelOwnedConnectionSchedule3D& schedule,
    const mpmc::mesh::PartitionSnapshot& partition,
    const mpmc::discretization_petsc::
        PetscMpiAijSymbolicPreallocation3D& cell_bridge,
    const mpmc::discretization_petsc::
        OwnedCellStructuralColumnPatternSnapshot3D& cell_pattern,
    double timestep_seconds,
    Pr76SinglePhaseTransitionTargetRebuildPlan3D
        plan,
    std::span<
        const SinglePhaseSnesAuthoritativeFaceInput3D>
        source_faces,
    CellScopedMixedCardinalityEvaluatorDispatcher3D
        dispatcher,
    const mpmc::thermodynamics::Pr76Phase<double>&
        model,
    Pr76AbsentPhaseSelectedBranchBinding3D
        absent_branch,
    std::unique_ptr<
        Pr76TransitionRebuildMaterializedSystem3D>*
            output) {
    using namespace
        pr76_transition_rebuild_materialization_detail;

    if (output == nullptr) {
        return PETSC_ERR_ARG_NULL;
    }
    output->reset();
    if (!std::isfinite(timestep_seconds) ||
        !(timestep_seconds > 0.0) ||
        absent_branch.resolver == nullptr ||
        source_faces.size() !=
            schedule.assembly_rows().size()) {
        return PETSC_ERR_ARG_INCOMP;
    }

    auto runtime =
        std::make_unique<
            Pr76TransitionRebuildMaterializedSystem3D>();
    runtime->dispatcher =
        std::make_unique<
            CellScopedMixedCardinalityEvaluatorDispatcher3D>(
                std::move(dispatcher));

    auto faces =
        mixed_faces(
            source_faces);
    int bootstrap_token = 1;
    std::unique_ptr<
        PhaseTransitionRebuiltNaturalVariableSystem3D>
        bootstrap;
    PetscErrorCode error =
        rebuild_phase_transition_natural_variable_system_3d(
            comm,
            schedule,
            partition,
            cell_bridge,
            cell_pattern,
            timestep_seconds,
            plan.cells,
            faces,
            runtime->dispatcher
                ->bindings(),
            &bootstrap_absent_provider,
            &bootstrap_token,
            &bootstrap);
    if (error != PETSC_SUCCESS ||
        bootstrap == nullptr) {
        return error != PETSC_SUCCESS
            ? error
            : PETSC_ERR_PLIB;
    }

    std::vector<std::optional<
        MixedCardinalityPhysicalCurrentCellLinearization3D>>
        current;
    std::vector<double> porosities;
    NaturalVariableSnesEvaluationStatus3D
        cell_status =
            NaturalVariableSnesEvaluationStatus3D::
                success;
    error =
        bootstrap
            ->evaluate_local_cells_for_phase_transition(
                bootstrap->initial_state(),
                &current,
                &porosities,
                &cell_status);
    if (error != PETSC_SUCCESS ||
        cell_status !=
            NaturalVariableSnesEvaluationStatus3D::
                success ||
        current.size() !=
            plan.cells.size() ||
        porosities.size() !=
            plan.cells.size()) {
        return error != PETSC_SUCCESS
            ? error
            : PETSC_ERR_ARG_WRONGSTATE;
    }

    std::vector<
        mpmc::flow::
            Pr76CellScopedAbsentPhasePotentialExtensionProvider<
                double>::Binding>
        provider_bindings;
    int local_unresolved = 0;
    PetscErrorCode local_error =
        PETSC_SUCCESS;

    try {
        for (const auto& row :
             schedule.assembly_rows()) {
            const auto owner_local =
                partition.local_index(
                    mpmc::mesh::EntityKind::cell,
                    row.owner_cell_global);
            const auto neighbour_local =
                partition.local_index(
                    mpmc::mesh::EntityKind::cell,
                    row.neighbour_cell_global);
            const auto owner =
                static_cast<std::size_t>(
                    owner_local.value());
            const auto neighbour =
                static_cast<std::size_t>(
                    neighbour_local.value());

            const auto phase_plan =
                mpmc::flow::
                    make_cross_cardinality_face_phase_identity_plan(
                        plan.cells[owner]
                            .target_active_phases,
                        plan.cells[neighbour]
                            .target_active_phases);
            if (!phase_plan
                    .requires_absent_phase_extension()) {
                continue;
            }

            for (const auto& binding :
                 phase_plan.bindings()) {
                if (binding
                        .owner_active_phase_index
                        .has_value() ==
                    binding
                        .neighbour_active_phase_index
                        .has_value()) {
                    continue;
                }

                const bool owner_absent =
                    !binding
                         .owner_active_phase_index
                         .has_value();
                const std::size_t absent_local =
                    owner_absent
                    ? owner
                    : neighbour;
                const std::size_t active_local =
                    owner_absent
                    ? neighbour
                    : owner;
                const auto active_global =
                    plan.cells[active_local]
                        .cell_global;
                const auto absent_global =
                    plan.cells[absent_local]
                        .cell_global;
                const auto* seed =
                    resolved_phase(
                        plan,
                        active_global,
                        binding.identity);
                if (seed == nullptr ||
                    !current[absent_local]
                        .has_value()) {
                    local_unresolved = 1;
                    continue;
                }

                const auto& host =
                    state_identity(
                        *current[absent_local]);
                std::optional<
                    mpmc::thermodynamics::
                        Pr76SelectedPhase>
                    absent_selection;
                local_error =
                    absent_branch.resolver(
                        absent_global,
                        binding.identity,
                        host,
                        seed->composition,
                        seed->selection,
                        absent_branch.user_context,
                        &absent_selection);
                if (local_error !=
                    PETSC_SUCCESS) {
                    break;
                }
                if (!absent_selection.has_value()) {
                    local_unresolved = 1;
                    continue;
                }

                const auto branch_provenance =
                    absent_branch_provenance(
                        absent_global,
                        binding.identity,
                        absent_selection
                            ->root_index);

                mpmc::flow::
                    AbsentPhaseThermodynamicCoordinateExtension
                    coordinates;
                coordinates.identity =
                    binding.identity;
                coordinates.host_state_identity =
                    host;
                coordinates.phase_pressure_pa =
                    host.reference_pressure_pa;
                coordinates.phase_pressure_gradient
                    .assign(
                        host.layout
                            .unknown_count(),
                        0.0);
                coordinates.phase_pressure_gradient[
                    host.layout
                        .pressure_unknown_index()] =
                    1.0;
                coordinates.hypothetical_composition =
                    seed->composition;
                coordinates
                    .hypothetical_composition_jacobian
                    .assign(
                        host.layout
                                .component_count() *
                            host.layout
                                .unknown_count(),
                        0.0);
                coordinates.selected_branch_provenance =
                    branch_provenance;
                coordinates.provenance =
                    std::string{
                        "PR76/PT frozen-absent-reference active-cell="} +
                    std::to_string(
                        active_global.value()) +
                    "|evidence=" +
                    std::string{
                        plan.cells[active_local]
                            .transition_evidence_profile};
                coordinates.validate();

                auto& entries =
                    plan.cells[absent_local]
                        .frozen_absent_phases;
                const auto existing =
                    std::find_if(
                        entries.begin(),
                        entries.end(),
                        [&](const auto& entry) {
                            return entry.identity ==
                                binding.identity;
                        });
                if (existing !=
                    entries.end()) {
                    if (!same_vector(
                            existing
                                ->reference_coordinates
                                .hypothetical_composition,
                            coordinates
                                .hypothetical_composition) ||
                        existing
                            ->selected_branch_provenance !=
                            branch_provenance) {
                        local_error =
                            PETSC_ERR_ARG_INCOMP;
                        break;
                    }
                } else {
                    entries.push_back(
                        {
                            binding.identity,
                            coordinates,
                            branch_provenance});
                }

                const auto provider_existing =
                    std::find_if(
                        provider_bindings.begin(),
                        provider_bindings.end(),
                        [&](const auto& entry) {
                            return entry
                                .branch_provenance ==
                                branch_provenance;
                        });
                const std::string evidence =
                    plan.cells[active_local]
                        .transition_evidence_profile;
                if (provider_existing ==
                    provider_bindings.end()) {
                    provider_bindings.push_back(
                        {
                            binding.identity,
                            *absent_selection,
                            evidence,
                            branch_provenance});
                } else if (
                    provider_existing->identity !=
                        binding.identity ||
                    provider_existing
                            ->selection.root_index !=
                        absent_selection
                            ->root_index) {
                    local_error =
                        PETSC_ERR_ARG_INCOMP;
                    break;
                }
            }
            if (local_error !=
                PETSC_SUCCESS) {
                break;
            }
        }
    } catch (...) {
        local_error =
            PETSC_ERR_ARG_INCOMP;
    }

    int local_error_int =
        static_cast<int>(
            local_error);
    int global_error_int = 0;
    if (MPI_Allreduce(
            &local_error_int,
            &global_error_int,
            1,
            MPI_INT,
            MPI_MAX,
            comm) != MPI_SUCCESS) {
        return PETSC_ERR_MPI;
    }
    if (global_error_int != 0) {
        return static_cast<PetscErrorCode>(
            global_error_int);
    }

    int global_unresolved = 0;
    if (MPI_Allreduce(
            &local_unresolved,
            &global_unresolved,
            1,
            MPI_INT,
            MPI_MAX,
            comm) != MPI_SUCCESS) {
        return PETSC_ERR_MPI;
    }
    if (global_unresolved != 0) {
        return PETSC_SUCCESS;
    }

    runtime->absent_provider =
        std::make_unique<
            mpmc::flow::
                Pr76CellScopedAbsentPhasePotentialExtensionProvider<
                    double>>(
                        model,
                        std::move(
                            provider_bindings));

    bootstrap.reset();
    faces =
        mixed_faces(
            source_faces);
    error =
        rebuild_phase_transition_natural_variable_system_3d(
            comm,
            schedule,
            partition,
            cell_bridge,
            cell_pattern,
            timestep_seconds,
            std::move(
                plan.cells),
            std::move(faces),
            runtime->dispatcher
                ->bindings(),
            &evaluate_absent_phase_thermodynamic_provider_3d<
                mpmc::flow::
                    Pr76CellScopedAbsentPhasePotentialExtensionProvider<
                        double>>,
            runtime->absent_provider
                .get(),
            &runtime->system);
    if (error != PETSC_SUCCESS ||
        runtime->system == nullptr) {
        return error != PETSC_SUCCESS
            ? error
            : PETSC_ERR_PLIB;
    }

    output->reset(
        runtime.release());
    return PETSC_SUCCESS;
}

} // namespace mpmc::flow_discretization_petsc

#endif // MPMC_FLOW_DISCRETIZATION_PETSC_PR76_TRANSITION_REBUILD_MATERIALIZATION_HPP
