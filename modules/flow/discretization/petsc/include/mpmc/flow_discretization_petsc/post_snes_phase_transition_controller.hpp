#ifndef MPMC_FLOW_DISCRETIZATION_PETSC_POST_SNES_PHASE_TRANSITION_CONTROLLER_HPP
#define MPMC_FLOW_DISCRETIZATION_PETSC_POST_SNES_PHASE_TRANSITION_CONTROLLER_HPP

#include <mpmc/flow_discretization_petsc/phase_transition_outer_rebuild.hpp>

#include <petscvec.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdio>
#include <cstdint>
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
    post_snes_phase_transition_controller_convention =
        "flow_discretization_petsc/post-snes-phase-transition-controller/v1";

struct PostSnesPhaseTransitionProposal3D {
    mpmc::mesh::GlobalEntityId cell_global{
        mpmc::mesh::GlobalEntityId::value_type{0}};
    mpmc::flow::PhaseSetTransitionCandidate
        candidate;
};

struct AcceptedPhaseTransitionSummary3D {
    mpmc::mesh::GlobalEntityId cell_global{
        mpmc::mesh::GlobalEntityId::value_type{0}};
    mpmc::mesh::PartitionRank owner_rank{
        mpmc::mesh::PartitionRank::value_type{0}};
    std::size_t source_phase_count{};
    std::size_t target_phase_count{};
    mpmc::flow::PhaseSetTransitionTrigger
        trigger{
            mpmc::flow::
                PhaseSetTransitionTrigger::
                    provider_topology_witness};

    friend bool operator==(
        const AcceptedPhaseTransitionSummary3D&,
        const AcceptedPhaseTransitionSummary3D&) =
        default;
};

struct GlobalPhaseSetSignatureEntry3D {
    mpmc::mesh::GlobalEntityId cell_global{
        mpmc::mesh::GlobalEntityId::value_type{0}};
    std::size_t phase_count{};
    std::vector<
        mpmc::flow::FrozenPhysicalPhaseIdentity>
        ordered_active_phases;

    friend bool operator==(
        const GlobalPhaseSetSignatureEntry3D&,
        const GlobalPhaseSetSignatureEntry3D&) =
        default;
};

enum class PostSnesPhaseTransitionOutcome3D {
    stable_phase_set,
    nonlinear_solve_diverged,
    phase_set_scan_indeterminate,
    transition_restart_budget_exhausted,
    phase_set_cycle_detected
};

enum class PostSnesPhaseTransitionScanStatus3D {
    complete,
    indeterminate
};

struct PostSnesPhaseTransitionControllerOptions3D {
    std::size_t max_transition_restarts{8U};
};

struct PostSnesPhaseTransitionGenerationReport3D {
    std::size_t generation{};
    std::vector<
        GlobalPhaseSetSignatureEntry3D>
        phase_set_before_solve;
    VariableCardinalityNaturalVariableSnesSolveReport3D
        nonlinear_solve;
    std::optional<
        PostSnesPhaseTransitionScanStatus3D>
        phase_transition_scan_status;
    std::vector<
        AcceptedPhaseTransitionSummary3D>
        accepted_transition_batch;
};

struct PostSnesPhaseTransitionControllerReport3D {
    PostSnesPhaseTransitionOutcome3D outcome{
        PostSnesPhaseTransitionOutcome3D::
            nonlinear_solve_diverged};
    std::size_t transition_restarts{};
    std::vector<
        PostSnesPhaseTransitionGenerationReport3D>
        generations;

    [[nodiscard]] bool
    timestep_accepted() const noexcept {
        return outcome ==
            PostSnesPhaseTransitionOutcome3D::
                stable_phase_set;
    }
};

using PostSnesPhaseTransitionScanner3D =
    PetscErrorCode (*)(
        const PhaseTransitionRebuiltNaturalVariableSystem3D&
            system,
        Vec converged_state,
        const VariableCardinalityNaturalVariableSnesSolveReport3D&
            solve_report,
        void* user_context,
        PostSnesPhaseTransitionScanStatus3D*
            scan_status,
        std::vector<
            PostSnesPhaseTransitionProposal3D>*
                local_owned_proposals);

using PostSnesPhaseTransitionRebuildFactory3D =
    PetscErrorCode (*)(
        const PhaseTransitionRebuiltNaturalVariableSystem3D&
            current_system,
        Vec converged_state,
        const VariableCardinalityNaturalVariableSnesSolveReport3D&
            solve_report,
        std::span<
            const PostSnesPhaseTransitionProposal3D>
            local_owned_proposals,
        std::span<
            const AcceptedPhaseTransitionSummary3D>
            accepted_global_batch,
        void* user_context,
        std::unique_ptr<
            PhaseTransitionRebuiltNaturalVariableSystem3D>*
                rebuilt_system);

struct PostSnesPhaseTransitionControllerBindings3D {
    PostSnesPhaseTransitionScanner3D scanner{};
    void* scanner_context{};
    PostSnesPhaseTransitionRebuildFactory3D
        rebuild_factory{};
    void* rebuild_context{};
};

namespace post_snes_transition_detail {

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
valid_trigger_direction(
    const mpmc::flow::
        PhaseSetTransitionCandidate&
            candidate) noexcept {
    if (candidate.source_phase_count == 0U ||
        candidate.source_phase_count >
            mpmc::flow::
                fixed_three_phase_count ||
        candidate.target_phase_count == 0U ||
        candidate.target_phase_count >
            mpmc::flow::
                fixed_three_phase_count ||
        candidate.source_phase_count ==
            candidate.target_phase_count) {
        return false;
    }
    if (candidate.target_phase_count >
            candidate.source_phase_count &&
        candidate.trigger ==
            mpmc::flow::
                PhaseSetTransitionTrigger::
                    phase_disappearance) {
        return false;
    }
    if (candidate.target_phase_count <
            candidate.source_phase_count &&
        candidate.trigger !=
            mpmc::flow::
                PhaseSetTransitionTrigger::
                    phase_disappearance &&
        candidate.trigger !=
            mpmc::flow::
                PhaseSetTransitionTrigger::
                    provider_boundary_route) {
        return false;
    }
    return true;
}

[[nodiscard]] inline const
VariableCardinalityNaturalVariableCellDof3D*
find_owned_cell(
    const VariableCardinalityNaturalVariableNumbering3D&
        numbering,
    mpmc::mesh::GlobalEntityId
        cell_global) noexcept {
    for (const auto& record :
         numbering.cells()) {
        if (record.owner_rank ==
                numbering.local_rank() &&
            record.cell_global ==
                cell_global) {
            return &record;
        }
    }
    return nullptr;
}

[[nodiscard]] inline PetscErrorCode
validate_local_proposals(
    const VariableCardinalityNaturalVariableNumbering3D&
        numbering,
    std::vector<
        PostSnesPhaseTransitionProposal3D>*
            proposals) {
    if (proposals == nullptr) {
        return PETSC_ERR_ARG_NULL;
    }
    std::sort(
        proposals->begin(),
        proposals->end(),
        [](const auto& first,
           const auto& second) {
            return first.cell_global <
                second.cell_global;
        });

    for (std::size_t index = 0U;
         index < proposals->size();
         ++index) {
        const auto& proposal =
            proposals->at(index);
        const auto* current =
            find_owned_cell(
                numbering,
                proposal.cell_global);
        if (current == nullptr ||
            proposal.candidate.status !=
                mpmc::flow::
                    PhaseSetTransitionCandidateStatus::
                        target_resolved ||
            proposal.candidate.source_phase_count !=
                current->phase_count ||
            proposal.candidate.evidence_profile.empty() ||
            !valid_trigger_direction(
                proposal.candidate) ||
            (index > 0U &&
             proposals->at(index - 1U)
                     .cell_global ==
                 proposal.cell_global)) {
            return PETSC_ERR_ARG_INCOMP;
        }
    }
    return PETSC_SUCCESS;
}

inline void append_u64_be(
    std::vector<unsigned char>* output,
    std::uint64_t value) {
    if (output == nullptr) {
        throw std::invalid_argument(
            "null phase-signature serialization output");
    }
    for (int shift = 56;
         shift >= 0;
         shift -= 8) {
        output->push_back(
            static_cast<unsigned char>(
                (value >>
                 static_cast<unsigned int>(
                     shift)) &
                UINT64_C(0xff)));
    }
}

[[nodiscard]] inline bool read_u64_be(
    std::span<const unsigned char> input,
    std::size_t* cursor,
    std::uint64_t* output) {
    if (cursor == nullptr ||
        output == nullptr ||
        *cursor > input.size() ||
        input.size() - *cursor < 8U) {
        return false;
    }
    std::uint64_t value = 0U;
    for (std::size_t byte = 0U;
         byte < 8U;
         ++byte) {
        value =
            (value << 8U) |
            static_cast<std::uint64_t>(
                input[*cursor + byte]);
    }
    *cursor += 8U;
    *output = value;
    return true;
}

inline void append_string(
    std::vector<unsigned char>* output,
    std::string_view value) {
    append_u64_be(
        output,
        static_cast<std::uint64_t>(
            value.size()));
    output->insert(
        output->end(),
        value.begin(),
        value.end());
}

[[nodiscard]] inline bool read_string(
    std::span<const unsigned char> input,
    std::size_t* cursor,
    std::string* output) {
    std::uint64_t size_u64 = 0U;
    if (output == nullptr ||
        !read_u64_be(
            input,
            cursor,
            &size_u64) ||
        size_u64 >
            static_cast<std::uint64_t>(
                input.size()) ||
        size_u64 >
            static_cast<std::uint64_t>(
                std::numeric_limits<
                    std::size_t>::max())) {
        return false;
    }
    const std::size_t size =
        static_cast<std::size_t>(
            size_u64);
    if (*cursor > input.size() ||
        size > input.size() - *cursor) {
        return false;
    }
    output->assign(
        reinterpret_cast<const char*>(
            input.data() + *cursor),
        size);
    *cursor += size;
    return true;
}

[[nodiscard]] inline PetscErrorCode
gather_phase_signature(
    MPI_Comm comm,
    const PhaseTransitionRebuiltNaturalVariableSystem3D&
        system,
    std::vector<
        GlobalPhaseSetSignatureEntry3D>*
            output) {
    if (output == nullptr) {
        return PETSC_ERR_ARG_NULL;
    }
    output->clear();

    const auto& numbering =
        system.numbering();
    std::vector<unsigned char>
        local;
    try {
        for (const auto& record :
             numbering.cells()) {
            if (record.owner_rank !=
                numbering.local_rank()) {
                continue;
            }
            const auto& registry_cell =
                system
                    .coordinate_registry()
                    .cell(
                        record.cell_global);
            if (registry_cell
                    .active_phases
                    .phase_count() !=
                record.phase_count) {
                return PETSC_ERR_ARG_INCOMP;
            }

            append_u64_be(
                &local,
                record.cell_global.value());
            append_u64_be(
                &local,
                static_cast<std::uint64_t>(
                    record.phase_count));
            append_u64_be(
                &local,
                static_cast<std::uint64_t>(
                    registry_cell
                        .active_phases
                        .phase_count()));
            for (const auto& identity :
                 registry_cell
                     .active_phases
                     .identities()) {
                append_string(
                    &local,
                    identity.provenance_scope);
                append_string(
                    &local,
                    identity.opaque_phase_key);
            }
        }
    } catch (...) {
        return PETSC_ERR_ARG_INCOMP;
    }

    if (local.size() >
        static_cast<std::size_t>(
            std::numeric_limits<int>::max())) {
        return PETSC_ERR_ARG_OUTOFRANGE;
    }
    const int local_bytes =
        static_cast<int>(
            local.size());

    int mpi_size = 0;
    if (MPI_Comm_size(
            comm,
            &mpi_size) != MPI_SUCCESS ||
        mpi_size <= 0) {
        return PETSC_ERR_MPI;
    }
    std::vector<int> counts(
        static_cast<std::size_t>(
            mpi_size),
        0);
    if (MPI_Allgather(
            &local_bytes,
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
    int total_bytes = 0;
    for (std::size_t rank = 0U;
         rank < counts.size();
         ++rank) {
        if (counts[rank] < 0 ||
            counts[rank] >
                std::numeric_limits<int>::max() -
                    total_bytes) {
            return PETSC_ERR_ARG_OUTOFRANGE;
        }
        displacements[rank] =
            total_bytes;
        total_bytes +=
            counts[rank];
    }

    std::vector<unsigned char>
        gathered(
            static_cast<std::size_t>(
                total_bytes));
    if (MPI_Allgatherv(
            local.empty()
                ? nullptr
                : local.data(),
            local_bytes,
            MPI_BYTE,
            gathered.empty()
                ? nullptr
                : gathered.data(),
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
            const std::size_t begin =
                static_cast<std::size_t>(
                    displacements[rank]);
            const std::size_t count =
                static_cast<std::size_t>(
                    counts[rank]);
            const auto payload =
                std::span<const unsigned char>{
                    gathered.data() + begin,
                    count};
            std::size_t cursor = 0U;

            while (cursor < payload.size()) {
                std::uint64_t cell_id = 0U;
                std::uint64_t phase_count_u64 = 0U;
                std::uint64_t identity_count_u64 = 0U;
                if (!read_u64_be(
                        payload,
                        &cursor,
                        &cell_id) ||
                    !read_u64_be(
                        payload,
                        &cursor,
                        &phase_count_u64) ||
                    !read_u64_be(
                        payload,
                        &cursor,
                        &identity_count_u64) ||
                    phase_count_u64 == 0U ||
                    phase_count_u64 >
                        mpmc::flow::
                            fixed_three_phase_count ||
                    identity_count_u64 !=
                        phase_count_u64) {
                    return PETSC_ERR_ARG_INCOMP;
                }

                std::vector<
                    mpmc::flow::
                        FrozenPhysicalPhaseIdentity>
                    identities;
                identities.reserve(
                    static_cast<std::size_t>(
                        identity_count_u64));
                for (std::uint64_t phase = 0U;
                     phase < identity_count_u64;
                     ++phase) {
                    std::string scope;
                    std::string key;
                    if (!read_string(
                            payload,
                            &cursor,
                            &scope) ||
                        !read_string(
                            payload,
                            &cursor,
                            &key) ||
                        scope.empty() ||
                        key.empty()) {
                        return PETSC_ERR_ARG_INCOMP;
                    }
                    identities.push_back(
                        {
                            std::move(scope),
                            std::move(key)});
                }

                mpmc::flow::
                    FrozenActivePhaseIdentityMap
                    validated{
                        identities};
                output->push_back(
                    {
                        mpmc::mesh::GlobalEntityId{
                            cell_id},
                        static_cast<std::size_t>(
                            phase_count_u64),
                        std::vector<
                            mpmc::flow::
                                FrozenPhysicalPhaseIdentity>{
                            validated
                                .identities()
                                .begin(),
                            validated
                                .identities()
                                .end()}});
            }
            if (cursor != payload.size()) {
                return PETSC_ERR_ARG_INCOMP;
            }
        }
    } catch (...) {
        return PETSC_ERR_ARG_INCOMP;
    }

    std::sort(
        output->begin(),
        output->end(),
        [](const auto& first,
           const auto& second) {
            return first.cell_global <
                second.cell_global;
        });
    for (std::size_t index = 1U;
         index < output->size();
         ++index) {
        if (output->at(index - 1U)
                .cell_global ==
            output->at(index)
                .cell_global) {
            return PETSC_ERR_ARG_INCOMP;
        }
    }
    return PETSC_SUCCESS;
}

[[nodiscard]] inline PetscErrorCode
gather_accepted_batch(
    MPI_Comm comm,
    const VariableCardinalityNaturalVariableNumbering3D&
        numbering,
    std::span<
        const PostSnesPhaseTransitionProposal3D>
        local_proposals,
    std::vector<
        AcceptedPhaseTransitionSummary3D>*
            output) {
    if (output == nullptr) {
        return PETSC_ERR_ARG_NULL;
    }
    output->clear();

    constexpr std::size_t width = 5U;
    std::vector<std::uint64_t>
        local;
    local.reserve(
        local_proposals.size() *
        width);
    for (const auto& proposal :
         local_proposals) {
        const auto* record =
            find_owned_cell(
                numbering,
                proposal.cell_global);
        if (record == nullptr) {
            return PETSC_ERR_ARG_INCOMP;
        }
        local.push_back(
            proposal.cell_global.value());
        local.push_back(
            static_cast<std::uint64_t>(
                record->owner_rank.value()));
        local.push_back(
            static_cast<std::uint64_t>(
                proposal.candidate
                    .source_phase_count));
        local.push_back(
            static_cast<std::uint64_t>(
                proposal.candidate
                    .target_phase_count));
        local.push_back(
            static_cast<std::uint64_t>(
                proposal.candidate
                    .trigger));
    }

    if (local_proposals.size() >
        static_cast<std::size_t>(
            std::numeric_limits<int>::max()) /
            width) {
        return PETSC_ERR_ARG_OUTOFRANGE;
    }
    const int local_records =
        static_cast<int>(
            local_proposals.size());

    int mpi_size = 0;
    if (MPI_Comm_size(
            comm,
            &mpi_size) != MPI_SUCCESS ||
        mpi_size <= 0) {
        return PETSC_ERR_MPI;
    }
    std::vector<int> counts(
        static_cast<std::size_t>(
            mpi_size),
        0);
    if (MPI_Allgather(
            &local_records,
            1,
            MPI_INT,
            counts.data(),
            1,
            MPI_INT,
            comm) != MPI_SUCCESS) {
        return PETSC_ERR_MPI;
    }

    std::vector<int> counts_u64(
        counts.size(),
        0);
    std::vector<int> displacements_u64(
        counts.size(),
        0);
    int total_u64 = 0;
    for (std::size_t rank = 0U;
         rank < counts.size();
         ++rank) {
        if (counts[rank] < 0 ||
            counts[rank] >
                (std::numeric_limits<int>::max() -
                 total_u64) /
                    static_cast<int>(
                        width)) {
            return PETSC_ERR_ARG_OUTOFRANGE;
        }
        counts_u64[rank] =
            counts[rank] *
            static_cast<int>(
                width);
        displacements_u64[rank] =
            total_u64;
        total_u64 +=
            counts_u64[rank];
    }

    std::vector<std::uint64_t>
        gathered(
            static_cast<std::size_t>(
                total_u64));
    if (MPI_Allgatherv(
            local.empty()
                ? nullptr
                : local.data(),
            static_cast<int>(
                local.size()),
            MPI_UINT64_T,
            gathered.empty()
                ? nullptr
                : gathered.data(),
            counts_u64.data(),
            displacements_u64.data(),
            MPI_UINT64_T,
            comm) != MPI_SUCCESS) {
        return PETSC_ERR_MPI;
    }

    if (gathered.size() % width != 0U) {
        return PETSC_ERR_PLIB;
    }
    output->reserve(
        gathered.size() /
        width);
    for (std::size_t index = 0U;
         index < gathered.size();
         index += width) {
        const auto source =
            gathered[index + 2U];
        const auto target =
            gathered[index + 3U];
        const auto trigger =
            gathered[index + 4U];
        if (source == 0U ||
            source >
                mpmc::flow::
                    fixed_three_phase_count ||
            target == 0U ||
            target >
                mpmc::flow::
                    fixed_three_phase_count ||
            source == target ||
            trigger >
                static_cast<std::uint64_t>(
                    mpmc::flow::
                        PhaseSetTransitionTrigger::
                            provider_boundary_route)) {
            return PETSC_ERR_ARG_INCOMP;
        }
        output->push_back(
            {
                mpmc::mesh::GlobalEntityId{
                    gathered[index]},
                mpmc::mesh::PartitionRank{
                    static_cast<
                        mpmc::mesh::
                            PartitionRank::
                                value_type>(
                        gathered[index + 1U])},
                static_cast<std::size_t>(
                    source),
                static_cast<std::size_t>(
                    target),
                static_cast<
                    mpmc::flow::
                        PhaseSetTransitionTrigger>(
                    trigger)});
    }
    std::sort(
        output->begin(),
        output->end(),
        [](const auto& first,
           const auto& second) {
            return first.cell_global <
                second.cell_global;
        });
    for (std::size_t index = 1U;
         index < output->size();
         ++index) {
        if (output->at(index - 1U)
                .cell_global ==
            output->at(index)
                .cell_global) {
            return PETSC_ERR_ARG_INCOMP;
        }
    }
    return PETSC_SUCCESS;
}

[[nodiscard]] inline PetscErrorCode
validate_rebuilt_signature_against_batch(
    std::span<
        const GlobalPhaseSetSignatureEntry3D>
        current,
    std::span<
        const AcceptedPhaseTransitionSummary3D>
        batch,
    std::span<
        const GlobalPhaseSetSignatureEntry3D>
        rebuilt) {
    if (current.size() !=
        rebuilt.size()) {
        return PETSC_ERR_ARG_INCOMP;
    }

    for (std::size_t index = 0U;
         index < current.size();
         ++index) {
        if (current[index].cell_global !=
            rebuilt[index].cell_global) {
            return PETSC_ERR_ARG_INCOMP;
        }
        const auto transition =
            std::lower_bound(
                batch.begin(),
                batch.end(),
                current[index].cell_global,
                [](const auto& entry,
                   mpmc::mesh::GlobalEntityId value) {
                    return entry.cell_global <
                        value;
                });
        if (transition != batch.end() &&
            transition->cell_global ==
                current[index].cell_global) {
            if (current[index].phase_count !=
                    transition
                        ->source_phase_count ||
                rebuilt[index].phase_count !=
                    transition
                        ->target_phase_count) {
                return PETSC_ERR_ARG_INCOMP;
            }
        } else if (
            current[index] !=
            rebuilt[index]) {
            return PETSC_ERR_ARG_INCOMP;
        }
    }
    return PETSC_SUCCESS;
}

[[nodiscard]] inline bool
signature_seen(
    std::span<
        const std::vector<
            GlobalPhaseSetSignatureEntry3D>>
        seen,
    std::span<
        const GlobalPhaseSetSignatureEntry3D>
        candidate) {
    return std::any_of(
        seen.begin(),
        seen.end(),
        [&](const auto& previous) {
            return previous.size() ==
                    candidate.size() &&
                std::equal(
                    previous.begin(),
                    previous.end(),
                    candidate.begin());
        });
}

} // namespace post_snes_transition_detail

inline PetscErrorCode
solve_nonlinear_timestep_with_phase_transitions_3d(
    MPI_Comm comm,
    std::unique_ptr<
        PhaseTransitionRebuiltNaturalVariableSystem3D>
        initial_system,
    PostSnesPhaseTransitionControllerBindings3D
        bindings,
    PostSnesPhaseTransitionControllerOptions3D
        options,
    std::unique_ptr<
        PhaseTransitionRebuiltNaturalVariableSystem3D>*
            final_system,
    Vec* final_state,
    std::optional<
        PostSnesPhaseTransitionControllerReport3D>*
            report) {
    using namespace
        post_snes_transition_detail;

    if (final_system == nullptr ||
        final_state == nullptr ||
        report == nullptr) {
        return PETSC_ERR_ARG_NULL;
    }
    if (*final_state != nullptr ||
        *final_system != nullptr) {
        return PETSC_ERR_ARG_WRONGSTATE;
    }
    report->reset();

    PetscErrorCode local_error =
        initial_system == nullptr ||
                bindings.scanner == nullptr ||
                bindings.rebuild_factory == nullptr
            ? PETSC_ERR_ARG_NULL
            : PETSC_SUCCESS;
    PetscErrorCode error =
        collective_error(
            comm,
            local_error);
    if (error != PETSC_SUCCESS) {
        return error;
    }

    auto current =
        std::move(initial_system);
    std::vector<
        std::vector<
            GlobalPhaseSetSignatureEntry3D>>
        seen_signatures;
    std::vector<
        GlobalPhaseSetSignatureEntry3D>
        signature;
    error =
        gather_phase_signature(
            comm,
            *current,
            &signature);
    if (error != PETSC_SUCCESS) {
        return error;
    }
    seen_signatures.push_back(
        signature);

    PostSnesPhaseTransitionControllerReport3D
        completed;
    completed.transition_restarts = 0U;

    for (std::size_t generation = 0U;;
         ++generation) {
        const char* diagnostic_stage = "solve-current-system";
        const auto diagnose = [&](PetscErrorCode code) {
            std::fprintf(stderr,
                "[phase-transition controller] generation=%zu stage=%s error=%d\n",
                generation, diagnostic_stage, static_cast<int>(code));
            std::size_t logged_cells = 0U;
            for (const auto& cell : current->numbering().cells()) {
                if (logged_cells == 8U) {
                    std::fprintf(stderr, "[phase-transition controller] remaining cells omitted\n");
                    break;
                }
                if (cell.owner_rank == current->numbering().local_rank()) {
                    std::fprintf(stderr,
                        "[phase-transition controller] owned_cell=%llu phase_count=%zu\n",
                        static_cast<unsigned long long>(cell.cell_global.value()),
                        cell.phase_count);
                    ++logged_cells;
                }
            }
        };
        Vec solved = nullptr;
        std::optional<
            VariableCardinalityNaturalVariableSnesSolveReport3D>
            solve_report;
        error =
            current->solve(
                &solved,
                &solve_report);
        if (error != PETSC_SUCCESS) {
            if (solved != nullptr) {
                (void)VecDestroy(
                    &solved);
            }
            diagnose(error);
            return error;
        }
        if (!solve_report.has_value() ||
            solved == nullptr) {
            if (solved != nullptr) {
                (void)VecDestroy(
                    &solved);
            }
            return PETSC_ERR_PLIB;
        }

        PostSnesPhaseTransitionGenerationReport3D
            generation_report;
        generation_report.generation =
            generation;
        generation_report
            .phase_set_before_solve =
            signature;
        generation_report.nonlinear_solve =
            *solve_report;

        if (static_cast<int>(
                solve_report
                    ->converged_reason) <= 0) {
            completed.outcome =
                PostSnesPhaseTransitionOutcome3D::
                    nonlinear_solve_diverged;
            completed.generations.push_back(
                std::move(
                    generation_report));
            *final_state = solved;
            *final_system =
                std::move(current);
            report->emplace(
                std::move(completed));
            return PETSC_SUCCESS;
        }

        std::vector<
            PostSnesPhaseTransitionProposal3D>
            local_proposals;
        PostSnesPhaseTransitionScanStatus3D
            local_scan_status =
                PostSnesPhaseTransitionScanStatus3D::
                    complete;
        diagnostic_stage = "scan";
        local_error =
            bindings.scanner(
                *current,
                solved,
                *solve_report,
                bindings.scanner_context,
                &local_scan_status,
                &local_proposals);
        error =
            collective_error(
                comm,
                local_error);
        if (error != PETSC_SUCCESS) {
            (void)VecDestroy(
                &solved);
            diagnose(error);
            return error;
        }

        int local_indeterminate =
            local_scan_status ==
                    PostSnesPhaseTransitionScanStatus3D::
                        indeterminate
                ? 1
                : 0;
        int global_indeterminate = 0;
        if (MPI_Allreduce(
                &local_indeterminate,
                &global_indeterminate,
                1,
                MPI_INT,
                MPI_MAX,
                comm) != MPI_SUCCESS) {
            (void)VecDestroy(
                &solved);
            return PETSC_ERR_MPI;
        }
        generation_report
            .phase_transition_scan_status =
            global_indeterminate != 0
                ? PostSnesPhaseTransitionScanStatus3D::
                      indeterminate
                : PostSnesPhaseTransitionScanStatus3D::
                      complete;
        if (global_indeterminate != 0) {
            completed.outcome =
                PostSnesPhaseTransitionOutcome3D::
                    phase_set_scan_indeterminate;
            completed.generations.push_back(
                std::move(
                    generation_report));
            *final_state = solved;
            *final_system =
                std::move(current);
            report->emplace(
                std::move(completed));
            return PETSC_SUCCESS;
        }

        diagnostic_stage = "validate-local-proposals";
        local_error =
            validate_local_proposals(
                current->numbering(),
                &local_proposals);
        error =
            collective_error(
                comm,
                local_error);
        if (error != PETSC_SUCCESS) {
            (void)VecDestroy(
                &solved);
            diagnose(error);
            return error;
        }

        std::vector<
            AcceptedPhaseTransitionSummary3D>
            accepted_batch;
        diagnostic_stage = "gather-accepted-batch";
        error =
            gather_accepted_batch(
                comm,
                current->numbering(),
                local_proposals,
                &accepted_batch);
        if (error != PETSC_SUCCESS) {
            (void)VecDestroy(
                &solved);
            diagnose(error);
            return error;
        }
        generation_report
            .accepted_transition_batch =
            accepted_batch;

        if (accepted_batch.empty()) {
            completed.outcome =
                PostSnesPhaseTransitionOutcome3D::
                    stable_phase_set;
            completed.generations.push_back(
                std::move(
                    generation_report));
            *final_state = solved;
            *final_system =
                std::move(current);
            report->emplace(
                std::move(completed));
            return PETSC_SUCCESS;
        }

        if (completed.transition_restarts >=
            options.max_transition_restarts) {
            completed.outcome =
                PostSnesPhaseTransitionOutcome3D::
                    transition_restart_budget_exhausted;
            completed.generations.push_back(
                std::move(
                    generation_report));
            *final_state = solved;
            *final_system =
                std::move(current);
            report->emplace(
                std::move(completed));
            return PETSC_SUCCESS;
        }

        std::unique_ptr<
            PhaseTransitionRebuiltNaturalVariableSystem3D>
            rebuilt;
        diagnostic_stage = "rebuild-factory";
        local_error =
            bindings.rebuild_factory(
                *current,
                solved,
                *solve_report,
                local_proposals,
                accepted_batch,
                bindings.rebuild_context,
                &rebuilt);
        error =
            collective_error(
                comm,
                local_error);
        if (error != PETSC_SUCCESS) {
            (void)VecDestroy(
                &solved);
            diagnose(error);
            return error;
        }
        if (rebuilt == nullptr) {
            (void)VecDestroy(
                &solved);
            return PETSC_ERR_PLIB;
        }

        std::vector<
            GlobalPhaseSetSignatureEntry3D>
            rebuilt_signature;
        diagnostic_stage = "gather-rebuilt-signature";
        error =
            gather_phase_signature(
                comm,
                *rebuilt,
                &rebuilt_signature);
        if (error != PETSC_SUCCESS) {
            (void)VecDestroy(
                &solved);
            diagnose(error);
            return error;
        }
        diagnostic_stage = "validate-rebuilt-signature";
        error =
            validate_rebuilt_signature_against_batch(
                signature,
                accepted_batch,
                rebuilt_signature);
        if (error != PETSC_SUCCESS) {
            (void)VecDestroy(
                &solved);
            diagnose(error);
            return error;
        }

        if (signature_seen(
                seen_signatures,
                rebuilt_signature)) {
            completed.outcome =
                PostSnesPhaseTransitionOutcome3D::
                    phase_set_cycle_detected;
            completed.generations.push_back(
                std::move(
                    generation_report));
            *final_state = solved;
            *final_system =
                std::move(current);
            report->emplace(
                std::move(completed));
            return PETSC_SUCCESS;
        }

        completed.generations.push_back(
            std::move(
                generation_report));
        ++completed.transition_restarts;
        seen_signatures.push_back(
            rebuilt_signature);
        signature =
            std::move(
                rebuilt_signature);

        diagnostic_stage = "destroy-old-state";
        error =
            VecDestroy(
                &solved);
        if (error != PETSC_SUCCESS) {
            diagnose(error);
            return error;
        }
        current =
            std::move(rebuilt);
    }
}

} // namespace mpmc::flow_discretization_petsc

#endif // MPMC_FLOW_DISCRETIZATION_PETSC_POST_SNES_PHASE_TRANSITION_CONTROLLER_HPP
