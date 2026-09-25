#ifndef MPMC_FLOW_DISCRETIZATION_PETSC_SINGLE_PHASE_HANDOFF_INITIAL_SYSTEM_HPP
#define MPMC_FLOW_DISCRETIZATION_PETSC_SINGLE_PHASE_HANDOFF_INITIAL_SYSTEM_HPP

#include <mpmc/flow_discretization_petsc/phase_transition_outer_rebuild.hpp>

#include <petscsys.h>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <memory>
#include <optional>
#include <stdexcept>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

namespace mpmc::flow_discretization_petsc {

inline constexpr std::string_view
    single_phase_handoff_initial_system_convention =
        "flow_discretization_petsc/single-phase-handoff-initial-system/v1";

/// Convert one converged fixed-cardinality 1P solve into the source-topology
/// PhaseTransitionRebuiltNaturalVariableSystem3D consumed by the existing outer
/// phase-transition controller.
///
/// Only locally owned q values are read from the fixed SNES report. Ghost
/// snapshots deliberately carry no target q: the rebuilt PETSc Vec owns only
/// local rows, while ghost state is synchronized from owners by the existing
/// mixed-cardinality PetscSF path. This avoids a global solution all-gather at
/// the control-plane handoff boundary.
inline PetscErrorCode
make_single_phase_handoff_initial_system_3d(
    MPI_Comm comm,
    const mpmc::discretization_petsc::
        ParallelOwnedConnectionSchedule3D& schedule,
    const mpmc::mesh::PartitionSnapshot& partition,
    const mpmc::discretization_petsc::
        PetscMpiAijSymbolicPreallocation3D& cell_bridge,
    const mpmc::discretization_petsc::
        OwnedCellStructuralColumnPatternSnapshot3D& cell_pattern,
    double timestep_seconds,
    const std::vector<SinglePhaseSnesCellInput3D>& cells,
    const std::vector<
        SinglePhaseSnesAuthoritativeFaceInput3D>& faces,
    const std::vector<
        mpmc::flow::FrozenActivePhaseIdentityMap>&
            source_active_phases,
    const NaturalVariableSnesSolveReport3D& solve_report,
    MixedCardinalityPhysicalCellEvaluatorBindings3D
        cell_evaluators,
    AbsentPhaseThermodynamicProviderEvaluator3D
        provider_evaluator,
    void* provider_context,
    std::unique_ptr<
        PhaseTransitionRebuiltNaturalVariableSystem3D>*
            output) {
    if (output == nullptr) {
        return PETSC_ERR_ARG_NULL;
    }
    output->reset();

    if (!std::isfinite(timestep_seconds) ||
        !(timestep_seconds > 0.0) ||
        cells.empty() ||
        cells.size() !=
            partition.entity_count(
                mpmc::mesh::EntityKind::cell) ||
        source_active_phases.size() !=
            cells.size() ||
        faces.size() !=
            schedule.assembly_rows().size() ||
        solve_report.component_count() < 2U ||
        solve_report.natural_variable_count() !=
            solve_report.component_count() + 1U ||
        static_cast<int>(
            solve_report.converged_reason()) <= 0 ||
        solve_report.local_rank() !=
            partition.local_rank() ||
        solve_report.rank_count() !=
            partition.rank_count()) {
        return PETSC_ERR_ARG_INCOMP;
    }

    std::vector<FrozenPhaseTransitionRebuildCell3D>
        snapshots;
    snapshots.reserve(cells.size());

    try {
        for (std::size_t local = 0U;
             local < cells.size();
             ++local) {
            const auto cell =
                mpmc::mesh::LocalIndex{
                    static_cast<
                        mpmc::mesh::LocalIndex::value_type>(
                            local)};
            const auto& source =
                cells[local];
            if (source.cell != cell ||
                source.cell_global !=
                    partition.global_id(
                        mpmc::mesh::EntityKind::cell,
                        cell) ||
                source.component_ids.size() !=
                    solve_report.component_count() ||
                source.frozen_layout.phase_count() !=
                    1U ||
                source_active_phases[local]
                        .phase_count() !=
                    1U) {
                throw std::invalid_argument(
                    "fixed 1P handoff cell metadata");
            }

            FrozenPhaseTransitionRebuildCell3D
                snapshot;
            snapshot.cell = source.cell;
            snapshot.cell_global =
                source.cell_global;
            snapshot.bulk_volume_m3 =
                source.bulk_volume_m3;
            snapshot.porosity =
                source.porosity;
            snapshot.component_ids =
                source.component_ids;
            snapshot.target_layout =
                source.frozen_layout
                    .descriptor();
            snapshot.target_active_phases =
                source_active_phases[local];
            snapshot.previous_component_accumulation =
                source.previous_component_accumulation;
            snapshot.previous_energy_accumulation =
                source.previous_energy_accumulation;
            snapshot.transition_evidence_profile =
                "flow_discretization_petsc/fixed-1p-handoff/source-topology/v1";

            if (partition.is_owned(
                    mpmc::mesh::EntityKind::cell,
                    cell)) {
                const std::size_t q =
                    source.frozen_layout
                        .unknown_count();
                snapshot.target_natural_variables
                    .assign(q, 0.0);
                std::vector<bool>
                    seen(q, false);
                for (const auto& entry :
                     solve_report
                         .locally_owned_solution()) {
                    if (entry.cell_global !=
                        source.cell_global) {
                        continue;
                    }
                    if (entry.natural_variable_slot >=
                            q ||
                        seen[entry
                                 .natural_variable_slot] ||
                        !std::isfinite(entry.value)) {
                        throw std::invalid_argument(
                            "fixed 1P handoff solution entries");
                    }
                    snapshot
                        .target_natural_variables[
                            entry.natural_variable_slot] =
                        entry.value;
                    seen[entry
                             .natural_variable_slot] =
                        true;
                }
                if (std::find(
                        seen.begin(),
                        seen.end(),
                        false) !=
                    seen.end()) {
                    throw std::invalid_argument(
                        "fixed 1P handoff missing owned q");
                }
            }

            snapshots.push_back(
                std::move(snapshot));
        }
    } catch (...) {
        return PETSC_ERR_ARG_INCOMP;
    }

    std::vector<
        MixedCardinalityPhysicalSnesAuthoritativeFaceInput3D>
        mixed_faces;
    mixed_faces.reserve(faces.size());
    for (const auto& face : faces) {
        mixed_faces.push_back(
            {
                face.face,
                face.face_global,
                face.transmissibility,
                face.gravity,
                face.owner_to_neighbour_displacement,
                face.thermal_conductance});
    }

    return
        rebuild_phase_transition_natural_variable_system_3d(
            comm,
            schedule,
            partition,
            cell_bridge,
            cell_pattern,
            timestep_seconds,
            std::move(snapshots),
            std::move(mixed_faces),
            cell_evaluators,
            provider_evaluator,
            provider_context,
            output);
}

} // namespace mpmc::flow_discretization_petsc

#endif // MPMC_FLOW_DISCRETIZATION_PETSC_SINGLE_PHASE_HANDOFF_INITIAL_SYSTEM_HPP
