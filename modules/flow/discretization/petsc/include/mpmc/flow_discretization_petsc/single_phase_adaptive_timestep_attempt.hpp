#ifndef MPMC_FLOW_DISCRETIZATION_PETSC_SINGLE_PHASE_ADAPTIVE_TIMESTEP_ATTEMPT_HPP
#define MPMC_FLOW_DISCRETIZATION_PETSC_SINGLE_PHASE_ADAPTIVE_TIMESTEP_ATTEMPT_HPP

#include <mpmc/flow_discretization_petsc/adaptive_timestep_controller.hpp>
#include <mpmc/flow_discretization_petsc/complete_natural_variable_petsc_materialization.hpp>
#include <mpmc/flow_discretization_petsc/single_phase_snes_assembly.hpp>

#include <petscvec.h>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <optional>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

namespace mpmc::flow_discretization_petsc {

inline constexpr std::string_view
    single_phase_adaptive_timestep_attempt_convention =
        "flow_discretization_petsc/single-phase-adaptive-timestep-attempt/v1";

using SinglePhaseAdaptiveTimestepPostSnesReview3D =
    PetscErrorCode (*)(
        const AdaptiveTimestepAttemptRequest3D& request,
        Vec converged_state,
        const NaturalVariableSnesSolveReport3D& solve_report,
        void* user_context,
        AdaptiveTimestepAttemptOutcome3D* outcome,
        std::size_t* phase_transition_restarts,
        std::vector<
            PostSnesPhaseTransitionProposal3D>*
                local_owned_proposals);

struct SinglePhaseAdaptiveTimestepPostSnesReviewBinding3D {
    SinglePhaseAdaptiveTimestepPostSnesReview3D reviewer{};
    void* user_context{};
};

/// Production adapter for one fixed-cardinality single-phase timestep attempt.
///
/// The accepted state and previous accumulation snapshots remain caller-owned
/// and are read-only here. Every call rebuilds the backward-Euler assembly with
/// the request timestep, rebuilds frozen row equilibration and MPIAIJ
/// materialization, then runs the existing PETSc SNESNEWTONLS/BT + GMRES/ASM
/// solver. A rejected attempt cannot modify accepted history.
///
/// A converged SNES solve is not automatically accepted. The mandatory
/// post-SNES reviewer owns phase-set/stability policy and may publish only a
/// stable phase set or one of the recoverable phase-review rejection outcomes.
/// A stable reviewed solve is retained only as a pending solution/report. The
/// adaptive controller commit callback must explicitly take that pending result
/// and advance accepted history. Starting another attempt first destroys any
/// uncommitted pending solution from the preceding rejection.
class SinglePhaseAdaptiveTimestepAttemptContext3D {
public:
    SinglePhaseAdaptiveTimestepAttemptContext3D(
        const SinglePhaseAdaptiveTimestepAttemptContext3D&) =
        delete;
    SinglePhaseAdaptiveTimestepAttemptContext3D& operator=(
        const SinglePhaseAdaptiveTimestepAttemptContext3D&) =
        delete;
    SinglePhaseAdaptiveTimestepAttemptContext3D& operator=(
        SinglePhaseAdaptiveTimestepAttemptContext3D&&) =
        delete;

    SinglePhaseAdaptiveTimestepAttemptContext3D(
        SinglePhaseAdaptiveTimestepAttemptContext3D&& other) noexcept
        : comm_(other.comm_),
          schedule_(other.schedule_),
          partition_(other.partition_),
          dof_layout_(other.dof_layout_),
          dof_numbering_(other.dof_numbering_),
          cell_bridge_(other.cell_bridge_),
          cell_pattern_(other.cell_pattern_),
          natural_variable_id_(
              std::move(other.natural_variable_id_)),
          accepted_cell_inputs_(
              other.accepted_cell_inputs_),
          face_inputs_(other.face_inputs_),
          cell_evaluator_(other.cell_evaluator_),
          accepted_state_(other.accepted_state_),
          post_snes_review_(other.post_snes_review_),
          pending_request_(
              std::move(other.pending_request_)),
          pending_solution_(other.pending_solution_),
          pending_report_(
              std::move(other.pending_report_)),
          pending_phase_transition_(
              other.pending_phase_transition_),
          pending_local_owned_proposals_(
              std::move(
                  other.pending_local_owned_proposals_)) {
        other.pending_solution_ = nullptr;
        other.pending_request_.reset();
        other.pending_report_.reset();
        other.pending_phase_transition_ = false;
        other.pending_local_owned_proposals_.clear();
    }

    ~SinglePhaseAdaptiveTimestepAttemptContext3D() {
        (void)discard_pending();
    }

    [[nodiscard]] static PetscErrorCode
    create(
        MPI_Comm comm,
        const mpmc::discretization_petsc::
            ParallelOwnedConnectionSchedule3D* schedule,
        const mpmc::mesh::PartitionSnapshot* partition,
        const mpmc::mesh::DofLayout* dof_layout,
        const mpmc::mesh::DofNumberingSnapshot*
            dof_numbering,
        const mpmc::discretization_petsc::
            PetscMpiAijSymbolicPreallocation3D*
                cell_bridge,
        const mpmc::discretization_petsc::
            OwnedCellStructuralColumnPatternSnapshot3D*
                cell_pattern,
        std::string natural_variable_id,
        const std::vector<SinglePhaseSnesCellInput3D>*
            accepted_cell_inputs,
        const std::vector<
            SinglePhaseSnesAuthoritativeFaceInput3D>*
                face_inputs,
        SinglePhaseCurrentCellEvaluatorBinding3D
            cell_evaluator,
        Vec accepted_state,
        SinglePhaseAdaptiveTimestepPostSnesReviewBinding3D
            post_snes_review,
        std::optional<
            SinglePhaseAdaptiveTimestepAttemptContext3D>*
                output) {
        if (output == nullptr) {
            return PETSC_ERR_ARG_NULL;
        }
        output->reset();

        if (schedule == nullptr ||
            partition == nullptr ||
            dof_layout == nullptr ||
            dof_numbering == nullptr ||
            cell_bridge == nullptr ||
            cell_pattern == nullptr ||
            accepted_cell_inputs == nullptr ||
            face_inputs == nullptr ||
            accepted_state == nullptr ||
            cell_evaluator.evaluator == nullptr ||
            post_snes_review.reviewer == nullptr) {
            return PETSC_ERR_ARG_NULL;
        }
        if (natural_variable_id.empty()) {
            return PETSC_ERR_ARG_INCOMP;
        }

        output->emplace(
            SinglePhaseAdaptiveTimestepAttemptContext3D{
                comm,
                schedule,
                partition,
                dof_layout,
                dof_numbering,
                cell_bridge,
                cell_pattern,
                std::move(natural_variable_id),
                accepted_cell_inputs,
                face_inputs,
                cell_evaluator,
                accepted_state,
                post_snes_review});
        return PETSC_SUCCESS;
    }

    [[nodiscard]] bool has_pending() const noexcept {
        return pending_request_.has_value() &&
            pending_solution_ != nullptr &&
            pending_report_.has_value();
    }

    [[nodiscard]] bool
    has_pending_transition() const noexcept {
        return has_pending() &&
            pending_phase_transition_;
    }

    [[nodiscard]] PetscErrorCode
    discard_pending() noexcept {
        PetscErrorCode error = PETSC_SUCCESS;
        if (pending_solution_ != nullptr) {
            error =
                VecDestroy(
                    &pending_solution_);
        }
        pending_request_.reset();
        pending_report_.reset();
        pending_phase_transition_ = false;
        pending_local_owned_proposals_.clear();
        return error;
    }

    [[nodiscard]] PetscErrorCode
    take_pending(
        const AdaptiveTimestepAttemptRequest3D& request,
        Vec* solution,
        std::optional<NaturalVariableSnesSolveReport3D>*
            report) {
        if (solution == nullptr || report == nullptr) {
            return PETSC_ERR_ARG_NULL;
        }
        if (*solution != nullptr) {
            return PETSC_ERR_ARG_WRONGSTATE;
        }
        report->reset();
        if (!has_pending() ||
            pending_phase_transition_) {
            return PETSC_ERR_ARG_WRONGSTATE;
        }
        if (pending_request_->attempt_index !=
                request.attempt_index ||
            pending_request_->retry_index !=
                request.retry_index ||
            pending_request_->timestep_seconds !=
                request.timestep_seconds) {
            return PETSC_ERR_ARG_INCOMP;
        }

        *solution = pending_solution_;
        pending_solution_ = nullptr;
        report->emplace(
            std::move(*pending_report_));
        pending_report_.reset();
        pending_request_.reset();
        pending_phase_transition_ = false;
        pending_local_owned_proposals_.clear();
        return PETSC_SUCCESS;
    }

    [[nodiscard]] PetscErrorCode
    take_pending_transition(
        const AdaptiveTimestepAttemptRequest3D& request,
        Vec* solution,
        std::optional<NaturalVariableSnesSolveReport3D>*
            report,
        std::vector<
            PostSnesPhaseTransitionProposal3D>*
                local_owned_proposals) {
        if (solution == nullptr ||
            report == nullptr ||
            local_owned_proposals == nullptr) {
            return PETSC_ERR_ARG_NULL;
        }
        if (*solution != nullptr) {
            return PETSC_ERR_ARG_WRONGSTATE;
        }
        report->reset();
        local_owned_proposals->clear();
        if (!has_pending_transition()) {
            return PETSC_ERR_ARG_WRONGSTATE;
        }
        if (pending_request_->attempt_index !=
                request.attempt_index ||
            pending_request_->retry_index !=
                request.retry_index ||
            pending_request_->timestep_seconds !=
                request.timestep_seconds) {
            return PETSC_ERR_ARG_INCOMP;
        }

        *solution = pending_solution_;
        pending_solution_ = nullptr;
        report->emplace(
            std::move(*pending_report_));
        *local_owned_proposals =
            std::move(
                pending_local_owned_proposals_);
        pending_report_.reset();
        pending_request_.reset();
        pending_phase_transition_ = false;
        pending_local_owned_proposals_.clear();
        return PETSC_SUCCESS;
    }

    [[nodiscard]] PetscErrorCode
    evaluate_attempt(
        const AdaptiveTimestepAttemptRequest3D& request,
        AdaptiveTimestepAttemptResult3D* result) {
        if (result == nullptr) {
            return PETSC_ERR_ARG_NULL;
        }
        if (!std::isfinite(request.timestep_seconds) ||
            !(request.timestep_seconds > 0.0)) {
            return PETSC_ERR_ARG_OUTOFRANGE;
        }

        PetscErrorCode error =
            discard_pending();
        if (error != PETSC_SUCCESS) {
            return error;
        }

        std::optional<SinglePhaseSnesAssemblyContext3D>
            assembly_context;
        error =
            SinglePhaseSnesAssemblyContext3D::create(
                comm_,
                *schedule_,
                *partition_,
                *dof_layout_,
                *dof_numbering_,
                *cell_bridge_,
                *cell_pattern_,
                natural_variable_id_,
                request.timestep_seconds,
                *accepted_cell_inputs_,
                *face_inputs_,
                cell_evaluator_,
                &assembly_context);
        if (error != PETSC_SUCCESS ||
            !assembly_context.has_value()) {
            return error != PETSC_SUCCESS
                ? error
                : PETSC_ERR_PLIB;
        }

        std::optional<
            CompleteNaturalVariableAssemblySnapshot3D>
            initial_assembly;
        NaturalVariableSnesEvaluationStatus3D
            initial_status =
                NaturalVariableSnesEvaluationStatus3D::
                    success;
        error =
            assembly_context->evaluate_complete_assembly(
                accepted_state_,
                &initial_assembly,
                &initial_status);
        if (error != PETSC_SUCCESS) {
            return error;
        }
        if (initial_status !=
                NaturalVariableSnesEvaluationStatus3D::
                    success ||
            !initial_assembly.has_value()) {
            return PETSC_ERR_ARG_WRONGSTATE;
        }

        Vec row_scaling = nullptr;
        Vec structural_residual = nullptr;
        Mat jacobian_template = nullptr;
        Vec solution = nullptr;

        const auto cleanup_temporaries =
            [&]() {
                PetscErrorCode first =
                    PETSC_SUCCESS;
                if (structural_residual != nullptr) {
                    const PetscErrorCode destroy =
                        VecDestroy(
                            &structural_residual);
                    if (first == PETSC_SUCCESS &&
                        destroy != PETSC_SUCCESS) {
                        first = destroy;
                    }
                }
                if (jacobian_template != nullptr) {
                    const PetscErrorCode destroy =
                        MatDestroy(
                            &jacobian_template);
                    if (first == PETSC_SUCCESS &&
                        destroy != PETSC_SUCCESS) {
                        first = destroy;
                    }
                }
                if (row_scaling != nullptr) {
                    const PetscErrorCode destroy =
                        VecDestroy(
                            &row_scaling);
                    if (first == PETSC_SUCCESS &&
                        destroy != PETSC_SUCCESS) {
                        first = destroy;
                    }
                }
                return first;
            };

        error =
            make_natural_variable_initial_row_equilibration_3d(
                comm_,
                *initial_assembly,
                &row_scaling);
        if (error == PETSC_SUCCESS) {
            error =
                materialize_complete_natural_variable_petsc_system_3d(
                    comm_,
                    *initial_assembly,
                    *cell_bridge_,
                    &structural_residual,
                    &jacobian_template);
        }
        if (error != PETSC_SUCCESS) {
            (void)cleanup_temporaries();
            return error;
        }

        std::optional<NaturalVariableSnesSolveReport3D>
            solve_report;
        std::optional<NaturalVariableSnesFailureDiagnostics3D>
            failure_diagnostics;
        const PetscErrorCode solve_error =
            solve_natural_variable_snes_3d(
                comm_,
                *initial_assembly,
                accepted_state_,
                jacobian_template,
                assembly_context->snes_evaluator(),
                &solution,
                &solve_report,
                row_scaling,
                &failure_diagnostics);

        const PetscErrorCode cleanup_error =
            cleanup_temporaries();
        if (cleanup_error != PETSC_SUCCESS) {
            if (solution != nullptr) {
                (void)VecDestroy(&solution);
            }
            return cleanup_error;
        }

        if (solve_error == PETSC_ERR_NOT_CONVERGED) {
            if (!failure_diagnostics.has_value() ||
                solution != nullptr ||
                solve_report.has_value()) {
                if (solution != nullptr) {
                    (void)VecDestroy(&solution);
                }
                return PETSC_ERR_PLIB;
            }
            *result =
                make_adaptive_timestep_attempt_result(
                    *failure_diagnostics);
            return PETSC_SUCCESS;
        }
        if (solve_error != PETSC_SUCCESS) {
            if (solution != nullptr) {
                (void)VecDestroy(&solution);
            }
            return solve_error;
        }
        if (solution == nullptr ||
            !solve_report.has_value()) {
            if (solution != nullptr) {
                (void)VecDestroy(&solution);
            }
            return PETSC_ERR_PLIB;
        }

        AdaptiveTimestepAttemptResult3D reviewed =
            make_adaptive_timestep_attempt_result(
                *solve_report);
        AdaptiveTimestepAttemptOutcome3D review_outcome =
            AdaptiveTimestepAttemptOutcome3D::
                stable_phase_set;
        std::size_t transition_restarts = 0U;
        std::vector<
            PostSnesPhaseTransitionProposal3D>
            local_owned_proposals;
        error =
            post_snes_review_.reviewer(
                request,
                solution,
                *solve_report,
                post_snes_review_.user_context,
                &review_outcome,
                &transition_restarts,
                &local_owned_proposals);
        if (error != PETSC_SUCCESS) {
            (void)VecDestroy(&solution);
            return error;
        }
        if (!valid_post_snes_review_outcome(
                review_outcome) ||
            !valid_review_proposals_collective(
                review_outcome,
                local_owned_proposals)) {
            (void)VecDestroy(&solution);
            return PETSC_ERR_ARG_INCOMP;
        }

        reviewed.outcome = review_outcome;
        reviewed.phase_transition_restarts =
            transition_restarts;
        *result = reviewed;

        if (review_outcome !=
                AdaptiveTimestepAttemptOutcome3D::
                    stable_phase_set &&
            review_outcome !=
                AdaptiveTimestepAttemptOutcome3D::
                    phase_transition_proposed) {
            return VecDestroy(&solution);
        }

        pending_request_ = request;
        pending_solution_ = solution;
        pending_report_.emplace(
            std::move(*solve_report));
        pending_phase_transition_ =
            review_outcome ==
                AdaptiveTimestepAttemptOutcome3D::
                    phase_transition_proposed;
        pending_local_owned_proposals_ =
            std::move(
                local_owned_proposals);
        return PETSC_SUCCESS;
    }

private:
    SinglePhaseAdaptiveTimestepAttemptContext3D(
        MPI_Comm comm,
        const mpmc::discretization_petsc::
            ParallelOwnedConnectionSchedule3D* schedule,
        const mpmc::mesh::PartitionSnapshot* partition,
        const mpmc::mesh::DofLayout* dof_layout,
        const mpmc::mesh::DofNumberingSnapshot*
            dof_numbering,
        const mpmc::discretization_petsc::
            PetscMpiAijSymbolicPreallocation3D*
                cell_bridge,
        const mpmc::discretization_petsc::
            OwnedCellStructuralColumnPatternSnapshot3D*
                cell_pattern,
        std::string natural_variable_id,
        const std::vector<SinglePhaseSnesCellInput3D>*
            accepted_cell_inputs,
        const std::vector<
            SinglePhaseSnesAuthoritativeFaceInput3D>*
                face_inputs,
        SinglePhaseCurrentCellEvaluatorBinding3D
            cell_evaluator,
        Vec accepted_state,
        SinglePhaseAdaptiveTimestepPostSnesReviewBinding3D
            post_snes_review)
        : comm_(comm),
          schedule_(schedule),
          partition_(partition),
          dof_layout_(dof_layout),
          dof_numbering_(dof_numbering),
          cell_bridge_(cell_bridge),
          cell_pattern_(cell_pattern),
          natural_variable_id_(
              std::move(natural_variable_id)),
          accepted_cell_inputs_(
              accepted_cell_inputs),
          face_inputs_(face_inputs),
          cell_evaluator_(cell_evaluator),
          accepted_state_(accepted_state),
          post_snes_review_(post_snes_review) {}

    [[nodiscard]] bool
    valid_review_proposals_collective(
        AdaptiveTimestepAttemptOutcome3D outcome,
        const std::vector<
            PostSnesPhaseTransitionProposal3D>&
                proposals) const {
        int local_invalid = 0;
        std::uint64_t local_count = 0U;

        try {
            if (outcome !=
                    AdaptiveTimestepAttemptOutcome3D::
                        phase_transition_proposed &&
                !proposals.empty()) {
                local_invalid = 1;
            }

            for (std::size_t index = 0U;
                 local_invalid == 0 &&
                 index < proposals.size();
                 ++index) {
                const auto& proposal =
                    proposals[index];
                const auto cell =
                    std::find_if(
                        accepted_cell_inputs_->begin(),
                        accepted_cell_inputs_->end(),
                        [&](const auto& input) {
                            return input.cell_global ==
                                proposal.cell_global;
                        });
                if (cell ==
                        accepted_cell_inputs_->end() ||
                    !partition_->is_owned(
                        mpmc::mesh::EntityKind::cell,
                        cell->cell) ||
                    proposal.candidate.source_phase_count !=
                        1U ||
                    proposal.candidate.target_phase_count <=
                        1U ||
                    proposal.candidate.target_phase_count >
                        mpmc::flow::
                            fixed_three_phase_count ||
                    proposal.candidate.status !=
                        mpmc::flow::
                            PhaseSetTransitionCandidateStatus::
                                target_resolved ||
                    proposal.candidate.component_ids !=
                        cell->component_ids ||
                    proposal.candidate.evidence_profile.empty() ||
                    !std::isfinite(
                        proposal.candidate.pressure_pa) ||
                    !(proposal.candidate.pressure_pa > 0.0) ||
                    !std::isfinite(
                        proposal.candidate.temperature_k) ||
                    !(proposal.candidate.temperature_k > 0.0) ||
                    (index > 0U &&
                     proposals[index - 1U]
                             .cell_global ==
                         proposal.cell_global)) {
                    local_invalid = 1;
                    break;
                }
                (void)mpmc::flow::
                    phase_set_transition_overall_composition(
                        proposal.candidate);
            }
            local_count =
                static_cast<std::uint64_t>(
                    proposals.size());
        } catch (...) {
            local_invalid = 1;
        }

        int global_invalid = 0;
        std::uint64_t global_count = 0U;
        if (MPI_Allreduce(
                &local_invalid,
                &global_invalid,
                1,
                MPI_INT,
                MPI_MAX,
                comm_) != MPI_SUCCESS ||
            MPI_Allreduce(
                &local_count,
                &global_count,
                1,
                MPI_UINT64_T,
                MPI_SUM,
                comm_) != MPI_SUCCESS) {
            return false;
        }

        if (global_invalid != 0) {
            return false;
        }
        return outcome ==
                   AdaptiveTimestepAttemptOutcome3D::
                       phase_transition_proposed
            ? global_count > 0U
            : global_count == 0U;
    }

    [[nodiscard]] static bool
    valid_post_snes_review_outcome(
        AdaptiveTimestepAttemptOutcome3D
            outcome) noexcept {
        return outcome ==
                   AdaptiveTimestepAttemptOutcome3D::
                       stable_phase_set ||
            outcome ==
                AdaptiveTimestepAttemptOutcome3D::
                    phase_transition_proposed ||
            outcome ==
                AdaptiveTimestepAttemptOutcome3D::
                    phase_set_scan_indeterminate ||
            outcome ==
                AdaptiveTimestepAttemptOutcome3D::
                    transition_restart_budget_exhausted ||
            outcome ==
                AdaptiveTimestepAttemptOutcome3D::
                    phase_set_cycle_detected;
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
    const std::vector<SinglePhaseSnesCellInput3D>*
        accepted_cell_inputs_{};
    const std::vector<
        SinglePhaseSnesAuthoritativeFaceInput3D>*
            face_inputs_{};
    SinglePhaseCurrentCellEvaluatorBinding3D
        cell_evaluator_;
    Vec accepted_state_{};
    SinglePhaseAdaptiveTimestepPostSnesReviewBinding3D
        post_snes_review_;

    std::optional<AdaptiveTimestepAttemptRequest3D>
        pending_request_;
    Vec pending_solution_{};
    std::optional<NaturalVariableSnesSolveReport3D>
        pending_report_;
    bool pending_phase_transition_{};
    std::vector<
        PostSnesPhaseTransitionProposal3D>
        pending_local_owned_proposals_;
};

inline PetscErrorCode
evaluate_single_phase_adaptive_timestep_attempt_3d(
    const AdaptiveTimestepAttemptRequest3D& request,
    void* raw_context,
    AdaptiveTimestepAttemptResult3D* result) {
    if (raw_context == nullptr) {
        return PETSC_ERR_ARG_NULL;
    }
    return static_cast<
        SinglePhaseAdaptiveTimestepAttemptContext3D*>(
            raw_context)
        ->evaluate_attempt(
            request,
            result);
}

} // namespace mpmc::flow_discretization_petsc

#endif // MPMC_FLOW_DISCRETIZATION_PETSC_SINGLE_PHASE_ADAPTIVE_TIMESTEP_ATTEMPT_HPP
