#ifndef MPMC_FLOW_DISCRETIZATION_PETSC_SW92_TRANSACTIONAL_PHASE_TRANSITION_RESTART_HPP
#define MPMC_FLOW_DISCRETIZATION_PETSC_SW92_TRANSACTIONAL_PHASE_TRANSITION_RESTART_HPP

#include <mpmc/flow_discretization_petsc/post_snes_phase_transition_controller.hpp>
#include <mpmc/flow_discretization_petsc/post_snes_sw92_profile_c_phase_transition_scanner.hpp>
#include <mpmc/flow_discretization_petsc/sw92_production_cell_evaluator.hpp>

#include <petscvec.h>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <functional>
#include <limits>
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
    sw92_transactional_phase_transition_restart_convention =
        "flow_discretization_petsc/sw92-transactional-same-dt-restart/no-cross-cardinality-face/v1";

using Sw92TransactionalTargetIdentityResolver3D =
    PetscErrorCode (*)(
        mpmc::mesh::GlobalEntityId cell_global,
        const mpmc::flow::PhaseSetTransitionCandidate& candidate,
        const Sw92AuthoritativeTargetMaterialization3D& target,
        const mpmc::flow::FrozenActivePhaseIdentityMap& source_active_phases,
        void* user_context,
        std::optional<mpmc::flow::FrozenActivePhaseIdentityMap>* target_active_phases);

struct Sw92TransactionalTargetIdentityBinding3D {
    Sw92TransactionalTargetIdentityResolver3D resolver{};
    void* user_context{};
};

struct Sw92TransactionalAcceptedCellBaseline3D {
    mpmc::mesh::LocalIndex cell{
        mpmc::mesh::LocalIndex::value_type{0}};
    mpmc::mesh::GlobalEntityId cell_global{
        mpmc::mesh::GlobalEntityId::value_type{0}};
    double bulk_volume_m3{};
    double porosity{};
    std::vector<std::string> component_ids;
    mpmc::flow::PoreVolumeComponentAccumulationSnapshot3P
        previous_component_accumulation;
    mpmc::flow::PoreVolumeEnergyAccumulationSnapshot3P
        previous_energy_accumulation;
};

struct Sw92TransactionalInitialCell3D {
    Sw92TransactionalAcceptedCellBaseline3D baseline;
    mpmc::flow::NaturalVariableLayoutDescriptor layout{
        2U, 1U, std::vector<std::size_t>{1U}};
    std::vector<double> natural_variables;
    mpmc::flow::FrozenActivePhaseIdentityMap active_phases{
        {mpmc::flow::FrozenPhysicalPhaseIdentity{
            "invalid","invalid"}}};
    std::vector<mpmc::thermodynamics::Sw92SelectedPhase<double>>
        selections;
};

template <typename Provider>
struct Sw92TransactionalGenerationRuntime3D {
    using Closure =
        mpmc::flow::Sw92SelectedPhasePropertyClosure<double, Provider>;

    std::unique_ptr<Closure> closure;
    std::unique_ptr<
        Sw92SinglePhaseProductionCellEvaluatorContext3D<Closure>>
        one_phase;
    std::unique_ptr<
        Sw92TwoPhaseProductionCellEvaluatorContext3D<Closure>>
        two_phase;
    std::unique_ptr<
        Sw92ThreePhaseProductionCellEvaluatorContext3D<Closure>>
        three_phase;

    [[nodiscard]]
    MixedCardinalityPhysicalCellEvaluatorBindings3D
    bindings() noexcept {
        return {
            {
                &evaluate_sw92_single_phase_production_cell_3d<Closure>,
                one_phase.get()},
            {
                &evaluate_sw92_two_phase_production_cell_3d<Closure>,
                two_phase.get()},
            {
                &evaluate_sw92_three_phase_production_cell_3d<Closure>,
                three_phase.get()}};
    }
};

template <typename Provider>
struct Sw92TransactionalPhaseTransitionRestartContext3D {
    const mpmc::thermodynamics::Sw92Phase<double>* model{};
    std::function<Provider()> provider_factory;
    mpmc::flow::SelectedPhasePropertyProvenance property_provenance;

    MPI_Comm comm{MPI_COMM_NULL};
    PostSnesSw92ProfileCPhaseTransitionScannerContext3D* scanner{};
    const mpmc::discretization_petsc::ParallelOwnedConnectionSchedule3D*
        schedule{};
    const mpmc::mesh::PartitionSnapshot* partition{};
    const mpmc::discretization_petsc::PetscMpiAijSymbolicPreallocation3D*
        cell_bridge{};
    const mpmc::discretization_petsc::OwnedCellStructuralColumnPatternSnapshot3D*
        cell_pattern{};

    double single_phase_relative_permeability{1.0};
    Sw92RockThermalStorageEvaluatorBinding3D rock;
    Sw92TwoPhaseRelativePermeabilityEvaluatorBinding3D two_phase_relative_permeability;
    Sw92ThreePhaseSaturationConstitutiveEvaluatorBinding3D three_phase_saturation;
    Sw92TransactionalTargetIdentityBinding3D target_identity;

    std::optional<Sw92TransactionalAcceptedCellBaseline3D>
        accepted_baseline;
    std::vector<std::unique_ptr<
        Sw92TransactionalGenerationRuntime3D<Provider>>>
        generations;
};

namespace sw92_transactional_restart_detail {

template <typename Provider>
[[nodiscard]] inline std::unique_ptr<
    Sw92TransactionalGenerationRuntime3D<Provider>>
make_generation_runtime(
    Sw92TransactionalPhaseTransitionRestartContext3D<Provider>& context,
    std::vector<mpmc::thermodynamics::Sw92SelectedPhase<double>>
        selections) {
    if (context.model == nullptr ||
        !context.provider_factory ||
        context.rock.evaluator == nullptr ||
        context.two_phase_relative_permeability.evaluator == nullptr ||
        context.three_phase_saturation.evaluator == nullptr ||
        selections.empty() ||
        selections.size() > mpmc::flow::fixed_three_phase_count) {
        throw std::invalid_argument(
            "mpmc::flow_discretization_petsc: incomplete SW92 transactional evaluator context");
    }

    auto runtime = std::make_unique<
        Sw92TransactionalGenerationRuntime3D<Provider>>();
    runtime->closure =
        std::make_unique<typename
            Sw92TransactionalGenerationRuntime3D<Provider>::Closure>(
                *context.model,
                std::move(selections),
                context.provider_factory(),
                context.property_provenance);
    runtime->one_phase =
        std::make_unique<
            Sw92SinglePhaseProductionCellEvaluatorContext3D<
                typename Sw92TransactionalGenerationRuntime3D<Provider>::Closure>>();
    runtime->one_phase->property_closure =
        runtime->closure.get();
    runtime->one_phase->relative_permeability =
        context.single_phase_relative_permeability;
    runtime->one_phase->rock = context.rock;

    runtime->two_phase =
        std::make_unique<
            Sw92TwoPhaseProductionCellEvaluatorContext3D<
                typename Sw92TransactionalGenerationRuntime3D<Provider>::Closure>>();
    runtime->two_phase->property_closure =
        runtime->closure.get();
    runtime->two_phase->relative_permeability =
        context.two_phase_relative_permeability;
    runtime->two_phase->rock = context.rock;

    runtime->three_phase =
        std::make_unique<
            Sw92ThreePhaseProductionCellEvaluatorContext3D<
                typename Sw92TransactionalGenerationRuntime3D<Provider>::Closure>>();
    runtime->three_phase->property_closure =
        runtime->closure.get();
    runtime->three_phase->saturation_constitutive =
        context.three_phase_saturation;
    runtime->three_phase->rock = context.rock;
    return runtime;
}

inline bool same_vector(
    std::span<const double> first,
    std::span<const double> second) {
    if (first.size() != second.size()) {
        return false;
    }
    for (std::size_t i = 0U; i < first.size(); ++i) {
        const double scale =
            std::max({1.0, std::abs(first[i]), std::abs(second[i])});
        if (!std::isfinite(first[i]) || !std::isfinite(second[i]) ||
            std::abs(first[i]-second[i]) >
                8192.0 * std::numeric_limits<double>::epsilon() * scale) {
            return false;
        }
    }
    return true;
}

[[nodiscard]] inline PetscErrorCode reject_absent_phase_provider(
    const mpmc::flow::AbsentPhaseThermodynamicCoordinateExtension&,
    void*,
    std::optional<mpmc::flow::AbsentPhasePotentialExtensionLinearization>*,
    NaturalVariableSnesEvaluationStatus3D*) {
    return PETSC_ERR_SUP;
}

[[nodiscard]] inline std::vector<double>
read_owned_cell_q(
    const PhaseTransitionRebuiltNaturalVariableSystem3D& system,
    Vec state,
    mpmc::mesh::LocalIndex cell) {
    const auto& record = system.numbering().cell(cell);
    if (record.owner_rank != system.numbering().local_rank()) {
        throw std::invalid_argument(
            "mpmc::flow_discretization_petsc: transactional SW92 rebuild requires an owned cell");
    }
    std::vector<PetscInt> indices(record.scalar_count);
    std::vector<PetscScalar> values(record.scalar_count);
    for (std::size_t slot=0U; slot<record.scalar_count; ++slot) {
        indices[slot] = record.petsc_global_scalar_start +
            static_cast<PetscInt>(slot);
    }
    const PetscErrorCode error =
        VecGetValues(state, static_cast<PetscInt>(indices.size()),
                     indices.data(), values.data());
    if (error != PETSC_SUCCESS) {
        throw std::runtime_error("VecGetValues failed");
    }
    std::vector<double> q(record.scalar_count);
    for (std::size_t slot=0U; slot<record.scalar_count; ++slot) {
        q[slot] = static_cast<double>(PetscRealPart(values[slot]));
    }
    return q;
}

[[nodiscard]] inline
mpmc::flow::PoreVolumeComponentAccumulationSnapshot3P
current_component_inventory(
    const MixedCardinalityPhysicalCurrentCellLinearization3D& current,
    double porosity) {
    return post_snes_pt_flash_scanner_detail::
        current_component_accumulation(current, porosity);
}

template <typename Provider>
inline void validate_context(
    const Sw92TransactionalPhaseTransitionRestartContext3D<Provider>& context) {
    if (context.model == nullptr || context.scanner == nullptr ||
        context.schedule == nullptr || context.partition == nullptr ||
        context.cell_bridge == nullptr || context.cell_pattern == nullptr ||
        context.target_identity.resolver == nullptr ||
        !context.provider_factory ||
        context.schedule->assembly_rows().size() != 0U ||
        context.partition->entity_count(mpmc::mesh::EntityKind::cell) != 1U ||
        !context.partition->is_owned(
            mpmc::mesh::EntityKind::cell,
            mpmc::mesh::LocalIndex{0U})) {
        throw std::invalid_argument(
            "mpmc::flow_discretization_petsc: SW92 transactional v1 requires one owned cell, explicit identity resolver, and zero cross-cardinality faces");
    }
}

} // namespace sw92_transactional_restart_detail

template <typename Provider>
[[nodiscard]] inline PetscErrorCode
materialize_sw92_transactional_initial_system_3d(
    MPI_Comm comm,
    double timestep_seconds,
    Sw92TransactionalPhaseTransitionRestartContext3D<Provider>* context,
    Sw92TransactionalInitialCell3D initial,
    std::unique_ptr<PhaseTransitionRebuiltNaturalVariableSystem3D>* output) {
    using namespace sw92_transactional_restart_detail;
    if (context == nullptr || output == nullptr) {
        return PETSC_ERR_ARG_NULL;
    }
    output->reset();

    try {
        validate_context(*context);
        if (comm == MPI_COMM_NULL ||
            !std::isfinite(timestep_seconds) ||
            !(timestep_seconds > 0.0) ||
            initial.baseline.cell != mpmc::mesh::LocalIndex{0U} ||
            initial.baseline.cell_global !=
                context->partition->global_id(
                    mpmc::mesh::EntityKind::cell,
                    mpmc::mesh::LocalIndex{0U}) ||
            initial.layout.phase_count() != initial.selections.size() ||
            initial.natural_variables.size() != initial.layout.unknown_count() ||
            initial.active_phases.phase_count() != initial.layout.phase_count()) {
            return PETSC_ERR_ARG_INCOMP;
        }
        if (initial.baseline.component_ids.size() !=
                initial.layout.component_count() ||
            initial.baseline.previous_component_accumulation.component_ids !=
                initial.baseline.component_ids ||
            initial.baseline.previous_energy_accumulation
                    .state_identity.component_ids !=
                initial.baseline.component_ids) {
            return PETSC_ERR_ARG_INCOMP;
        }

        auto runtime =
            make_generation_runtime(*context, initial.selections);
        FrozenPhaseTransitionRebuildCell3D cell;
        cell.cell = initial.baseline.cell;
        cell.cell_global = initial.baseline.cell_global;
        cell.bulk_volume_m3 = initial.baseline.bulk_volume_m3;
        cell.porosity = initial.baseline.porosity;
        cell.component_ids = initial.baseline.component_ids;
        cell.target_layout = initial.layout;
        cell.target_natural_variables = initial.natural_variables;
        cell.target_active_phases = initial.active_phases;
        cell.previous_component_accumulation =
            initial.baseline.previous_component_accumulation;
        cell.previous_energy_accumulation =
            initial.baseline.previous_energy_accumulation;
        cell.transition_evidence_profile =
            "SW92/transactional-initial-baseline/v1";

        std::unique_ptr<PhaseTransitionRebuiltNaturalVariableSystem3D> system;
        const PetscErrorCode error =
            rebuild_phase_transition_natural_variable_system_3d(
                comm,
                *context->schedule,
                *context->partition,
                *context->cell_bridge,
                *context->cell_pattern,
                timestep_seconds,
                std::vector<FrozenPhaseTransitionRebuildCell3D>{
                    std::move(cell)},
                {},
                runtime->bindings(),
                {},
                &reject_absent_phase_provider,
                context,
                &system);
        if (error != PETSC_SUCCESS || system == nullptr) {
            return error != PETSC_SUCCESS ? error : PETSC_ERR_PLIB;
        }
        context->comm = comm;
        context->accepted_baseline = initial.baseline;
        context->generations.push_back(std::move(runtime));
        *output = std::move(system);
        return PETSC_SUCCESS;
    } catch (const std::exception&) {
        return PETSC_ERR_ARG_INCOMP;
    }
}

template <typename Provider>
[[nodiscard]] inline PetscErrorCode
rebuild_sw92_transactional_phase_transition_3d(
    const PhaseTransitionRebuiltNaturalVariableSystem3D& current_system,
    Vec converged_state,
    const VariableCardinalityNaturalVariableSnesSolveReport3D& solve_report,
    std::span<const PostSnesPhaseTransitionProposal3D> local_owned_proposals,
    std::span<const AcceptedPhaseTransitionSummary3D> accepted_global_batch,
    void* raw_context,
    std::unique_ptr<PhaseTransitionRebuiltNaturalVariableSystem3D>* rebuilt_system) {
    using namespace sw92_transactional_restart_detail;
    if (raw_context == nullptr || rebuilt_system == nullptr ||
        converged_state == nullptr) {
        return PETSC_ERR_ARG_NULL;
    }
    rebuilt_system->reset();
    auto* context = static_cast<
        Sw92TransactionalPhaseTransitionRestartContext3D<Provider>*>(
            raw_context);

    try {
        validate_context(*context);
        if (!context->accepted_baseline.has_value() ||
            context->scanner == nullptr ||
            local_owned_proposals.size() != 1U ||
            accepted_global_batch.size() != 1U ||
            static_cast<int>(solve_report.converged_reason) <= 0) {
            return PETSC_ERR_ARG_INCOMP;
        }

        const auto& proposal = local_owned_proposals.front();
        const auto* target =
            context->scanner->find_target(proposal.cell_global);
        if (target == nullptr ||
            target->cell_global != proposal.cell_global ||
            target->component_ids != proposal.candidate.component_ids ||
            target->projection.source_phase_count() !=
                proposal.candidate.source_phase_count ||
            target->projection.target_phase_count() !=
                proposal.candidate.target_phase_count ||
            accepted_global_batch.front().cell_global != proposal.cell_global ||
            accepted_global_batch.front().source_phase_count !=
                proposal.candidate.source_phase_count ||
            accepted_global_batch.front().target_phase_count !=
                proposal.candidate.target_phase_count) {
            return PETSC_ERR_ARG_INCOMP;
        }

        const auto& record =
            current_system.numbering().cell(
                mpmc::mesh::LocalIndex{0U});
        if (record.cell_global != proposal.cell_global ||
            record.phase_count != proposal.candidate.source_phase_count) {
            return PETSC_ERR_ARG_INCOMP;
        }
        const auto q = read_owned_cell_q(
            current_system, converged_state, record.cell);

        std::optional<MixedCardinalityPhysicalCurrentCellLinearization3D>
            current;
        double porosity = 0.0;
        NaturalVariableSnesEvaluationStatus3D status =
            NaturalVariableSnesEvaluationStatus3D::success;
        PetscErrorCode error =
            current_system.evaluate_current_cell_for_phase_transition(
                record.cell, q, &current, &porosity, &status);
        if (error != PETSC_SUCCESS ||
            status != NaturalVariableSnesEvaluationStatus3D::success ||
            !current.has_value()) {
            return error != PETSC_SUCCESS ? error : PETSC_ERR_ARG_WRONGSTATE;
        }

        const auto& baseline = *context->accepted_baseline;
        if (baseline.cell != record.cell ||
            baseline.cell_global != record.cell_global ||
            baseline.component_ids != proposal.candidate.component_ids ||
            !post_snes_pt_flash_scanner_detail::near_roundoff(
                baseline.porosity, porosity) ||
            !std::isfinite(current_system.time_step_seconds()) ||
            !(current_system.time_step_seconds() > 0.0) ||
            context->comm == MPI_COMM_NULL) {
            return PETSC_ERR_ARG_INCOMP;
        }

        const auto& source_map =
            current_system.coordinate_registry()
                .cell(record.cell_global)
                .active_phases;
        std::optional<mpmc::flow::FrozenActivePhaseIdentityMap>
            target_map;
        error =
            context->target_identity.resolver(
                proposal.cell_global,
                proposal.candidate,
                *target,
                source_map,
                context->target_identity.user_context,
                &target_map);
        if (error != PETSC_SUCCESS) {
            return error;
        }
        if (!target_map.has_value()) {
            return PETSC_ERR_SUP;
        }

        const auto inventory =
            current_component_inventory(*current, porosity);
        auto rebuild_cell =
            make_accepted_phase_transition_rebuild_cell_3d(
                baseline.cell,
                baseline.cell_global,
                baseline.bulk_volume_m3,
                baseline.porosity,
                proposal.candidate,
                inventory,
                baseline.previous_component_accumulation,
                baseline.previous_energy_accumulation,
                source_map,
                *target_map,
                {});
        if (!same_vector(
                rebuild_cell.target_natural_variables,
                target->projection.natural_variables()) ||
            rebuild_cell.target_layout.phase_count() !=
                target->phases.size()) {
            return PETSC_ERR_PLIB;
        }

        std::vector<mpmc::thermodynamics::Sw92SelectedPhase<double>>
            selections;
        selections.reserve(target->phases.size());
        for (const auto& phase : target->phases) {
            selections.push_back(phase.selection);
        }
        auto runtime =
            make_generation_runtime(*context, std::move(selections));

        std::unique_ptr<PhaseTransitionRebuiltNaturalVariableSystem3D> next;
        error =
            rebuild_phase_transition_natural_variable_system_3d(
                context->comm,
                *context->schedule,
                *context->partition,
                *context->cell_bridge,
                *context->cell_pattern,
                current_system.time_step_seconds(),
                std::vector<FrozenPhaseTransitionRebuildCell3D>{
                    std::move(rebuild_cell)},
                {},
                runtime->bindings(),
                {},
                &reject_absent_phase_provider,
                context,
                &next);
        if (error != PETSC_SUCCESS || next == nullptr) {
            return error != PETSC_SUCCESS ? error : PETSC_ERR_PLIB;
        }

        context->generations.push_back(std::move(runtime));
        *rebuilt_system = std::move(next);
        return PETSC_SUCCESS;
    } catch (const std::exception&) {
        return PETSC_ERR_ARG_INCOMP;
    }
}

template <typename Provider>
[[nodiscard]] inline PetscErrorCode
solve_sw92_transactional_same_dt_3d(
    MPI_Comm comm,
    std::unique_ptr<PhaseTransitionRebuiltNaturalVariableSystem3D> initial_system,
    PostSnesSw92ProfileCPhaseTransitionScannerContext3D* scanner,
    Sw92TransactionalPhaseTransitionRestartContext3D<Provider>* rebuild,
    PostSnesPhaseTransitionControllerOptions3D options,
    std::unique_ptr<PhaseTransitionRebuiltNaturalVariableSystem3D>* final_system,
    Vec* final_state,
    std::optional<PostSnesPhaseTransitionControllerReport3D>* report) {
    if (scanner == nullptr || rebuild == nullptr ||
        scanner != rebuild->scanner) {
        return PETSC_ERR_ARG_INCOMP;
    }
    return solve_nonlinear_timestep_with_phase_transitions_3d(
        comm,
        std::move(initial_system),
        {
            &scan_post_snes_sw92_profile_c_phase_transitions_3d,
            scanner,
            &rebuild_sw92_transactional_phase_transition_3d<Provider>,
            rebuild},
        options,
        final_system,
        final_state,
        report);
}

} // namespace mpmc::flow_discretization_petsc

#endif
