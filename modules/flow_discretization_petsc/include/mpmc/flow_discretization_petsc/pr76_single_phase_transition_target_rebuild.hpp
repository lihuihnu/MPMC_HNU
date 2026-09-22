#ifndef MPMC_FLOW_DISCRETIZATION_PETSC_PR76_SINGLE_PHASE_TRANSITION_TARGET_REBUILD_HPP
#define MPMC_FLOW_DISCRETIZATION_PETSC_PR76_SINGLE_PHASE_TRANSITION_TARGET_REBUILD_HPP

#include <mpmc/flow/thermodynamics_absent_phase_extension.hpp>
#include <mpmc/flow_discretization_petsc/post_snes_phase_transition_controller.hpp>
#include <mpmc/flow_discretization_petsc/phase_transition_outer_rebuild.hpp>
#include <mpmc/flow_discretization_petsc/single_phase_snes_assembly.hpp>
#include <mpmc/thermodynamics/selected_phase_fugacity.hpp>

#include <petscsys.h>

#include <algorithm>
#include <bit>
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
    pr76_single_phase_transition_target_rebuild_convention =
        "flow_discretization_petsc/pr76-single-phase-transition-target-rebuild/v1";

using Pr76SinglePhaseTransitionTargetIdentityResolver3D =
    PetscErrorCode (*)(
        mpmc::mesh::GlobalEntityId cell_global,
        const mpmc::flow::PhaseSetTransitionCandidate& candidate,
        const mpmc::flow::FrozenActivePhaseIdentityMap&
            source_active_phases,
        std::span<
            const mpmc::thermodynamics::Pr76SelectedPhase>
            target_selected_phases,
        void* user_context,
        std::optional<
            mpmc::flow::FrozenActivePhaseIdentityMap>*
                target_active_phases);

struct Pr76SinglePhaseTransitionTargetIdentityBinding3D {
    Pr76SinglePhaseTransitionTargetIdentityResolver3D
        resolver{};
    void* user_context{};
};

struct Pr76SinglePhaseTransitionResolvedPhase3D {
    mpmc::flow::FrozenPhysicalPhaseIdentity
        identity;
    std::vector<double> composition;
    mpmc::thermodynamics::Pr76SelectedPhase
        selection;
    std::string selected_branch_provenance;
};

struct Pr76SinglePhaseTransitionResolvedCell3D {
    mpmc::mesh::GlobalEntityId cell_global{
        mpmc::mesh::GlobalEntityId::value_type{0}};
    mpmc::flow::PhaseIdentityContinuationSnapshot
        continuation;
    mpmc::flow::FrozenSelectedPhaseBranchRegistry<
        mpmc::thermodynamics::Pr76SelectedPhase>
        branch_registry;
    std::vector<
        Pr76SinglePhaseTransitionResolvedPhase3D>
        target_phases;
};

struct Pr76SinglePhaseTransitionTargetRebuildPlan3D {
    std::vector<
        FrozenPhaseTransitionRebuildCell3D>
        cells;
    std::vector<
        Pr76SinglePhaseTransitionResolvedCell3D>
        local_transition_cells;
};

namespace pr76_single_phase_target_rebuild_detail {

inline void append_u64(
    std::vector<unsigned char>* out,
    std::uint64_t value) {
    for (int shift = 56; shift >= 0; shift -= 8) {
        out->push_back(
            static_cast<unsigned char>(
                (value >> static_cast<unsigned>(shift)) &
                UINT64_C(0xff)));
    }
}

inline void append_double(
    std::vector<unsigned char>* out,
    double value) {
    append_u64(
        out,
        std::bit_cast<std::uint64_t>(
            value));
}

inline void append_string(
    std::vector<unsigned char>* out,
    std::string_view value) {
    append_u64(
        out,
        static_cast<std::uint64_t>(
            value.size()));
    out->insert(
        out->end(),
        value.begin(),
        value.end());
}

[[nodiscard]] inline bool read_u64(
    std::span<const unsigned char> bytes,
    std::size_t* cursor,
    std::uint64_t* value) {
    if (cursor == nullptr ||
        value == nullptr ||
        *cursor > bytes.size() ||
        bytes.size() - *cursor < 8U) {
        return false;
    }
    std::uint64_t result = 0U;
    for (std::size_t i = 0U; i < 8U; ++i) {
        result =
            (result << 8U) |
            static_cast<std::uint64_t>(
                bytes[*cursor + i]);
    }
    *cursor += 8U;
    *value = result;
    return true;
}

[[nodiscard]] inline bool read_double(
    std::span<const unsigned char> bytes,
    std::size_t* cursor,
    double* value) {
    std::uint64_t bits = 0U;
    if (!read_u64(
            bytes,
            cursor,
            &bits)) {
        return false;
    }
    *value =
        std::bit_cast<double>(
            bits);
    return std::isfinite(*value);
}

[[nodiscard]] inline bool read_string(
    std::span<const unsigned char> bytes,
    std::size_t* cursor,
    std::string* value) {
    std::uint64_t size = 0U;
    if (!read_u64(bytes, cursor, &size) ||
        size >
            static_cast<std::uint64_t>(
                bytes.size() - *cursor)) {
        return false;
    }
    value->assign(
        reinterpret_cast<const char*>(
            bytes.data() + *cursor),
        static_cast<std::size_t>(size));
    *cursor +=
        static_cast<std::size_t>(size);
    return true;
}

[[nodiscard]] inline std::vector<double>
owned_q(
    const SinglePhaseSnesCellInput3D& cell,
    const NaturalVariableSnesSolveReport3D& report) {
    const std::size_t q =
        cell.frozen_layout.unknown_count();
    std::vector<double> values(
        q,
        0.0);
    std::vector<bool> seen(
        q,
        false);
    for (const auto& entry :
         report.locally_owned_solution()) {
        if (entry.cell_global !=
            cell.cell_global) {
            continue;
        }
        if (entry.natural_variable_slot >= q ||
            seen[entry.natural_variable_slot] ||
            !std::isfinite(entry.value)) {
            throw std::invalid_argument(
                "PR76 target rebuild: malformed owned fixed-1P solution entries");
        }
        values[entry.natural_variable_slot] =
            entry.value;
        seen[entry.natural_variable_slot] =
            true;
    }
    if (std::find(
            seen.begin(),
            seen.end(),
            false) !=
        seen.end()) {
        throw std::invalid_argument(
            "PR76 target rebuild: owned fixed-1P q is incomplete");
    }
    return values;
}

struct ResolvedWire {
    mpmc::mesh::GlobalEntityId cell_global{
        mpmc::mesh::GlobalEntityId::value_type{0}};
    mpmc::flow::PhaseSetTransitionTrigger
        trigger{
            mpmc::flow::PhaseSetTransitionTrigger::
                stability_witness};
    std::vector<std::size_t>
        dependent_components;
    std::vector<
        mpmc::flow::FrozenPhysicalPhaseIdentity>
        target_identities;
    std::vector<std::size_t>
        target_root_indices;
    std::vector<std::vector<double>>
        target_compositions;
    std::string evidence_profile;
    std::string diagnostic;
};

inline void serialize_wire(
    const ResolvedWire& wire,
    std::vector<unsigned char>* out) {
    append_u64(
        out,
        wire.cell_global.value());
    append_u64(
        out,
        static_cast<std::uint64_t>(
            wire.trigger));
    append_u64(
        out,
        static_cast<std::uint64_t>(
            wire.target_identities.size()));
    for (std::size_t value :
         wire.dependent_components) {
        append_u64(
            out,
            static_cast<std::uint64_t>(
                value));
    }
    for (std::size_t phase = 0U;
         phase < wire.target_identities.size();
         ++phase) {
        append_string(
            out,
            wire.target_identities[phase]
                .provenance_scope);
        append_string(
            out,
            wire.target_identities[phase]
                .opaque_phase_key);
        append_u64(
            out,
            static_cast<std::uint64_t>(
                wire.target_root_indices[phase]));
        append_u64(
            out,
            static_cast<std::uint64_t>(
                wire.target_compositions[
                    phase].size()));
        for (double value :
             wire.target_compositions[phase]) {
            append_double(
                out,
                value);
        }
    }
    append_string(
        out,
        wire.evidence_profile);
    append_string(
        out,
        wire.diagnostic);
}

[[nodiscard]] inline bool parse_wire(
    std::span<const unsigned char> bytes,
    std::size_t* cursor,
    ResolvedWire* wire) {
    std::uint64_t cell = 0U;
    std::uint64_t trigger = 0U;
    std::uint64_t count = 0U;
    if (!read_u64(bytes, cursor, &cell) ||
        !read_u64(bytes, cursor, &trigger) ||
        !read_u64(bytes, cursor, &count) ||
        count < 2U ||
        count >
            mpmc::flow::
                fixed_three_phase_count ||
        trigger >
            static_cast<std::uint64_t>(
                mpmc::flow::
                    PhaseSetTransitionTrigger::
                        provider_boundary_route)) {
        return false;
    }

    wire->cell_global =
        mpmc::mesh::GlobalEntityId{cell};
    wire->trigger =
        static_cast<
            mpmc::flow::
                PhaseSetTransitionTrigger>(
                    trigger);
    const std::size_t p =
        static_cast<std::size_t>(count);
    wire->dependent_components.resize(p);
    wire->target_identities.resize(p);
    wire->target_root_indices.resize(p);
    wire->target_compositions.resize(p);
    for (std::size_t phase = 0U;
         phase < p;
         ++phase) {
        std::uint64_t dependent = 0U;
        if (!read_u64(
                bytes,
                cursor,
                &dependent)) {
            return false;
        }
        wire->dependent_components[phase] =
            static_cast<std::size_t>(
                dependent);
    }
    for (std::size_t phase = 0U;
         phase < p;
         ++phase) {
        std::uint64_t root = 0U;
        if (!read_string(
                bytes,
                cursor,
                &wire->target_identities[phase]
                     .provenance_scope) ||
            !read_string(
                bytes,
                cursor,
                &wire->target_identities[phase]
                     .opaque_phase_key) ||
            !read_u64(
                bytes,
                cursor,
                &root)) {
            return false;
        }
        wire->target_root_indices[phase] =
            static_cast<std::size_t>(
                root);
        std::uint64_t composition_size = 0U;
        if (!read_u64(
                bytes,
                cursor,
                &composition_size) ||
            composition_size < 2U ||
            composition_size >
                UINT64_C(1024)) {
            return false;
        }
        wire->target_compositions[phase]
            .resize(
                static_cast<std::size_t>(
                    composition_size));
        for (double& value :
             wire->target_compositions[phase]) {
            if (!read_double(
                    bytes,
                    cursor,
                    &value)) {
                return false;
            }
        }
    }
    return read_string(
               bytes,
               cursor,
               &wire->evidence_profile) &&
        read_string(
               bytes,
               cursor,
               &wire->diagnostic);
}

[[nodiscard]] inline std::string
selected_branch_provenance(
    mpmc::mesh::GlobalEntityId cell_global,
    const mpmc::flow::FrozenPhysicalPhaseIdentity&
        identity,
    std::size_t root_index) {
    return std::string{
               "PR76/PT cell="} +
        std::to_string(
            cell_global.value()) +
        "|phase=" +
        identity.provenance_scope +
        "/" +
        identity.opaque_phase_key +
        "|activity.branch=" +
        std::to_string(
            root_index);
}

[[nodiscard]] inline PetscErrorCode
gather_wires(
    MPI_Comm comm,
    std::span<const ResolvedWire> local,
    std::vector<ResolvedWire>* output) {
    if (output == nullptr) {
        return PETSC_ERR_ARG_NULL;
    }
    output->clear();

    std::vector<unsigned char> payload;
    append_u64(
        &payload,
        static_cast<std::uint64_t>(
            local.size()));
    for (const auto& wire : local) {
        serialize_wire(
            wire,
            &payload);
    }
    if (payload.size() >
        static_cast<std::size_t>(
            std::numeric_limits<int>::max())) {
        return PETSC_ERR_ARG_OUTOFRANGE;
    }

    int size = 0;
    if (MPI_Comm_size(
            comm,
            &size) != MPI_SUCCESS ||
        size <= 0) {
        return PETSC_ERR_MPI;
    }
    const int local_size =
        static_cast<int>(
            payload.size());
    std::vector<int> counts(
        static_cast<std::size_t>(size),
        0);
    if (MPI_Allgather(
            &local_size,
            1,
            MPI_INT,
            counts.data(),
            1,
            MPI_INT,
            comm) != MPI_SUCCESS) {
        return PETSC_ERR_MPI;
    }

    std::vector<int> displacements(
        counts.size(),
        0);
    int total = 0;
    for (std::size_t rank = 0U;
         rank < counts.size();
         ++rank) {
        if (counts[rank] < 0 ||
            counts[rank] >
                std::numeric_limits<int>::max() -
                    total) {
            return PETSC_ERR_ARG_OUTOFRANGE;
        }
        displacements[rank] = total;
        total += counts[rank];
    }

    std::vector<unsigned char> gathered(
        static_cast<std::size_t>(total));
    if (MPI_Allgatherv(
            payload.data(),
            local_size,
            MPI_BYTE,
            gathered.data(),
            counts.data(),
            displacements.data(),
            MPI_BYTE,
            comm) != MPI_SUCCESS) {
        return PETSC_ERR_MPI;
    }

    try {
        for (std::size_t rank = 0U;
             rank < counts.size();
             ++rank) {
            const auto begin =
                static_cast<std::size_t>(
                    displacements[rank]);
            const auto length =
                static_cast<std::size_t>(
                    counts[rank]);
            std::span<const unsigned char>
                bytes{
                    gathered.data() + begin,
                    length};
            std::size_t cursor = 0U;
            std::uint64_t count = 0U;
            if (!read_u64(
                    bytes,
                    &cursor,
                    &count)) {
                return PETSC_ERR_ARG_INCOMP;
            }
            for (std::uint64_t i = 0U;
                 i < count;
                 ++i) {
                ResolvedWire wire;
                if (!parse_wire(
                        bytes,
                        &cursor,
                        &wire)) {
                    return PETSC_ERR_ARG_INCOMP;
                }
                output->push_back(
                    std::move(wire));
            }
            if (cursor != bytes.size()) {
                return PETSC_ERR_ARG_INCOMP;
            }
        }

        std::sort(
            output->begin(),
            output->end(),
            [](const auto& first,
               const auto& second) {
                return first.cell_global <
                    second.cell_global;
            });
        for (std::size_t i = 1U;
             i < output->size();
             ++i) {
            if ((*output)[i - 1U]
                    .cell_global ==
                (*output)[i]
                    .cell_global) {
                return PETSC_ERR_ARG_INCOMP;
            }
        }
    } catch (...) {
        return PETSC_ERR_ARG_INCOMP;
    }
    return PETSC_SUCCESS;
}

} // namespace pr76_single_phase_target_rebuild_detail

[[nodiscard]] inline PetscErrorCode
make_pr76_single_phase_transition_target_rebuild_plan_3d(
    MPI_Comm comm,
    const mpmc::mesh::PartitionSnapshot& partition,
    const std::vector<SinglePhaseSnesCellInput3D>& source_cells,
    const NaturalVariableSnesSolveReport3D& solve_report,
    SinglePhaseCurrentCellEvaluatorBinding3D
        source_cell_evaluator,
    std::span<
        const mpmc::flow::FrozenActivePhaseIdentityMap>
        source_active_phases,
    std::span<
        const PostSnesPhaseTransitionProposal3D>
        local_owned_proposals,
    mpmc::thermodynamics::Pr76RootOptions
        root_options,
    Pr76SinglePhaseTransitionTargetIdentityBinding3D
        identity_resolver,
    std::optional<
        Pr76SinglePhaseTransitionTargetRebuildPlan3D>*
            output) {
    using namespace
        pr76_single_phase_target_rebuild_detail;

    if (output == nullptr) {
        return PETSC_ERR_ARG_NULL;
    }
    output->reset();
    if (source_cell_evaluator.evaluator == nullptr ||
        identity_resolver.resolver == nullptr ||
        source_cells.empty() ||
        source_cells.size() !=
            partition.entity_count(
                mpmc::mesh::EntityKind::cell) ||
        source_active_phases.size() !=
            source_cells.size() ||
        solve_report.natural_variable_count() !=
            solve_report.component_count() + 1U ||
        static_cast<int>(
            solve_report.converged_reason()) <= 0) {
        return PETSC_ERR_ARG_INCOMP;
    }

    Pr76SinglePhaseTransitionTargetRebuildPlan3D
        plan;
    plan.cells.reserve(
        source_cells.size());

    try {
        for (std::size_t local = 0U;
             local < source_cells.size();
             ++local) {
            const auto& source =
                source_cells[local];
            const auto cell =
                mpmc::mesh::LocalIndex{
                    static_cast<
                        mpmc::mesh::LocalIndex::value_type>(
                            local)};
            if (source.cell != cell ||
                source.cell_global !=
                    partition.global_id(
                        mpmc::mesh::EntityKind::cell,
                        cell) ||
                source.frozen_layout.phase_count() !=
                    1U ||
                source.component_ids.size() !=
                    solve_report.component_count() ||
                source_active_phases[local]
                        .phase_count() !=
                    1U) {
                throw std::invalid_argument(
                    "PR76 target rebuild: malformed source cell");
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
                "flow_discretization_petsc/pr76-fixed-1p/no-transition/v1";
            if (partition.is_owned(
                    mpmc::mesh::EntityKind::cell,
                    cell)) {
                snapshot.target_natural_variables =
                    owned_q(
                        source,
                        solve_report);
            }
            plan.cells.push_back(
                std::move(snapshot));
        }
    } catch (...) {
        return PETSC_ERR_ARG_INCOMP;
    }

    std::vector<ResolvedWire>
        local_wires;
    local_wires.reserve(
        local_owned_proposals.size());
    int local_identity_indeterminate = 0;
    PetscErrorCode local_error =
        PETSC_SUCCESS;

    for (const auto& proposal :
         local_owned_proposals) {
        try {
            if (!partition.contains_global(
                    mpmc::mesh::EntityKind::cell,
                    proposal.cell_global)) {
                local_error =
                    PETSC_ERR_ARG_INCOMP;
                break;
            }
            const auto cell =
                partition.local_index(
                    mpmc::mesh::EntityKind::cell,
                    proposal.cell_global);
            if (!partition.is_owned(
                    mpmc::mesh::EntityKind::cell,
                    cell)) {
                local_error =
                    PETSC_ERR_ARG_INCOMP;
                break;
            }
            const auto local =
                static_cast<std::size_t>(
                    cell.value());
            const auto& source =
                source_cells.at(local);
            const auto& candidate =
                proposal.candidate;
            if (candidate.source_phase_count !=
                    1U ||
                candidate.target_phase_count <=
                    1U ||
                candidate.target_phase_count >
                    mpmc::flow::
                        fixed_three_phase_count ||
                candidate.status !=
                    mpmc::flow::
                        PhaseSetTransitionCandidateStatus::
                            target_resolved ||
                candidate.target_phases.size() !=
                    candidate.target_phase_count ||
                candidate.component_ids !=
                    source.component_ids ||
                !source.previous_component_accumulation
                    .has_value() ||
                !source.previous_energy_accumulation
                    .has_value()) {
                local_error =
                    PETSC_ERR_ARG_INCOMP;
                break;
            }

            std::vector<
                mpmc::thermodynamics::
                    Pr76SelectedPhase>
                target_selections;
            target_selections.reserve(
                candidate.target_phase_count);
            for (const auto& target :
                 candidate.target_phases) {
                if (!target
                         .provider_activity_branch
                         .has_value()) {
                    local_error =
                        PETSC_ERR_ARG_INCOMP;
                    break;
                }
                target_selections.push_back(
                    {
                        *target
                             .provider_activity_branch,
                        root_options});
            }
            if (local_error !=
                PETSC_SUCCESS) {
                break;
            }

            std::optional<
                mpmc::flow::
                    FrozenActivePhaseIdentityMap>
                target_identities;
            local_error =
                identity_resolver.resolver(
                    proposal.cell_global,
                    candidate,
                    source_active_phases[local],
                    target_selections,
                    identity_resolver
                        .user_context,
                    &target_identities);
            if (local_error !=
                PETSC_SUCCESS) {
                break;
            }
            if (!target_identities.has_value()) {
                local_identity_indeterminate = 1;
                continue;
            }

            (void)mpmc::flow::
                make_phase_identity_continuation_snapshot(
                    candidate,
                    source_active_phases[local],
                    *target_identities);

            const auto q =
                owned_q(
                    source,
                    solve_report);
            std::optional<
                SinglePhaseCurrentCellLinearization3D>
                current;
            NaturalVariableSnesEvaluationStatus3D
                status =
                    NaturalVariableSnesEvaluationStatus3D::
                        success;
            local_error =
                source_cell_evaluator.evaluator(
                    source.cell,
                    source.cell_global,
                    q,
                    source.frozen_layout,
                    source.component_ids,
                    source_cell_evaluator
                        .user_context,
                    &current,
                    &status);
            if (local_error !=
                PETSC_SUCCESS) {
                break;
            }
            if (status !=
                    NaturalVariableSnesEvaluationStatus3D::
                        success ||
                !current.has_value()) {
                local_identity_indeterminate = 1;
                continue;
            }

            const auto current_inventory =
                mpmc::flow::
                    build_single_phase_component_accumulation(
                        current->state,
                        source.porosity);
            plan.cells[local] =
                make_accepted_phase_transition_rebuild_cell_3d(
                    source.cell,
                    source.cell_global,
                    source.bulk_volume_m3,
                    source.porosity,
                    candidate,
                    current_inventory,
                    *source
                         .previous_component_accumulation,
                    *source
                         .previous_energy_accumulation,
                    source_active_phases[local],
                    *target_identities,
                    {});

            ResolvedWire wire;
            wire.cell_global =
                proposal.cell_global;
            wire.trigger =
                candidate.trigger;
            wire.evidence_profile =
                candidate.evidence_profile;
            wire.diagnostic =
                candidate.diagnostic;
            wire.target_identities.assign(
                target_identities
                    ->identities().begin(),
                target_identities
                    ->identities().end());
            wire.target_root_indices.reserve(
                target_selections.size());
            wire.target_compositions.reserve(
                target_selections.size());
            wire.dependent_components.reserve(
                candidate.target_phase_count);
            for (std::size_t phase = 0U;
                 phase <
                     candidate.target_phase_count;
                 ++phase) {
                wire.target_root_indices
                    .push_back(
                        target_selections[phase]
                            .root_index);
                wire.target_compositions
                    .push_back(
                        candidate
                            .target_phases[phase]
                            .composition);
                wire.dependent_components
                    .push_back(
                        plan.cells[local]
                            .target_layout
                            .dependent_composition_component(
                                static_cast<
                                    mpmc::flow::
                                        PhaseSlot3>(
                                            phase)));
            }
            local_wires.push_back(
                std::move(wire));
        } catch (...) {
            local_error =
                PETSC_ERR_ARG_INCOMP;
            break;
        }
    }

    int local_error_int =
        static_cast<int>(local_error);
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

    int global_identity_indeterminate = 0;
    if (MPI_Allreduce(
            &local_identity_indeterminate,
            &global_identity_indeterminate,
            1,
            MPI_INT,
            MPI_MAX,
            comm) != MPI_SUCCESS) {
        return PETSC_ERR_MPI;
    }
    if (global_identity_indeterminate != 0) {
        return PETSC_SUCCESS;
    }

    std::vector<ResolvedWire>
        global_wires;
    PetscErrorCode error =
        gather_wires(
            comm,
            local_wires,
            &global_wires);
    if (error != PETSC_SUCCESS) {
        return error;
    }

    try {
        for (const auto& wire :
             global_wires) {
            if (!partition.contains_global(
                    mpmc::mesh::EntityKind::cell,
                    wire.cell_global)) {
                continue;
            }
            const auto cell =
                partition.local_index(
                    mpmc::mesh::EntityKind::cell,
                    wire.cell_global);
            const auto local =
                static_cast<std::size_t>(
                    cell.value());
            const std::size_t target_count =
                wire.target_identities.size();

            const mpmc::flow::
                FrozenActivePhaseIdentityMap
                target_map{
                    wire.target_identities};
            mpmc::flow::
                PhaseSetTransitionCandidate
                candidate;
            candidate.source_phase_count =
                1U;
            candidate.target_phase_count =
                target_count;
            candidate.trigger =
                wire.trigger;
            candidate.status =
                mpmc::flow::
                    PhaseSetTransitionCandidateStatus::
                        target_resolved;
            candidate.component_ids =
                source_cells[local]
                    .component_ids;
            candidate.evidence_profile =
                wire.evidence_profile;
            candidate.diagnostic =
                wire.diagnostic;

            auto continuation =
                mpmc::flow::
                    make_phase_identity_continuation_snapshot(
                        candidate,
                        source_active_phases[local],
                        target_map);
            if (wire.target_compositions.size() !=
                    target_count) {
                return PETSC_ERR_ARG_INCOMP;
            }
            std::vector<
                mpmc::flow::
                    FrozenSelectedPhaseBranchBinding<
                        mpmc::thermodynamics::
                            Pr76SelectedPhase>>
                bindings;
            std::vector<
                Pr76SinglePhaseTransitionResolvedPhase3D>
                resolved_phases;
            bindings.reserve(
                target_count);
            resolved_phases.reserve(
                target_count);
            for (std::size_t phase = 0U;
                 phase < target_count;
                 ++phase) {
                if (wire.target_compositions[
                        phase].size() !=
                    source_cells[local]
                        .component_ids.size()) {
                    return PETSC_ERR_ARG_INCOMP;
                }
                const auto provenance =
                    selected_branch_provenance(
                        wire.cell_global,
                        target_map.identity(
                            phase),
                        wire
                            .target_root_indices[
                                phase]);
                const mpmc::thermodynamics::
                    Pr76SelectedPhase
                    selection{
                        wire
                            .target_root_indices[
                                phase],
                        root_options};
                bindings.push_back(
                    {
                        target_map.identity(
                            phase),
                        selection,
                        provenance});
                resolved_phases.push_back(
                    {
                        target_map.identity(
                            phase),
                        wire.target_compositions[
                            phase],
                        selection,
                        provenance});
            }
            auto registry =
                mpmc::flow::
                    make_transition_selected_phase_branch_registry(
                        continuation,
                        std::move(bindings));

            if (partition.is_ghost(
                    mpmc::mesh::EntityKind::cell,
                    cell)) {
                plan.cells[local]
                    .target_layout =
                    mpmc::flow::
                        NaturalVariableLayoutDescriptor{
                            source_cells[local]
                                .component_ids
                                .size(),
                            target_count,
                            wire
                                .dependent_components};
                plan.cells[local]
                    .target_active_phases =
                    target_map;
                plan.cells[local]
                    .target_natural_variables
                    .clear();
                plan.cells[local]
                    .transition_evidence_profile =
                    wire.evidence_profile;
            } else {
                if (plan.cells[local]
                        .target_layout
                        .phase_count() !=
                        target_count ||
                    plan.cells[local]
                        .target_active_phases
                        .identities()
                        .size() !=
                        target_count) {
                    return PETSC_ERR_PLIB;
                }
            }

            plan.local_transition_cells
                .push_back(
                    {
                        wire.cell_global,
                        std::move(
                            continuation),
                        std::move(
                            registry),
                        std::move(
                            resolved_phases)});
        }
    } catch (...) {
        return PETSC_ERR_ARG_INCOMP;
    }

    output->emplace(
        std::move(plan));
    return PETSC_SUCCESS;
}

} // namespace mpmc::flow_discretization_petsc

#endif // MPMC_FLOW_DISCRETIZATION_PETSC_PR76_SINGLE_PHASE_TRANSITION_TARGET_REBUILD_HPP
