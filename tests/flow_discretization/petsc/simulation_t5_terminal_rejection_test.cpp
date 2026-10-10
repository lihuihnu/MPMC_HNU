#include <mpmc/flow_discretization_petsc/physical_timestep_driver.hpp>
#include <mpmc/simulation_petsc/simulation_runner.hpp>
#include <mpmc/well_discretization_petsc/fixed_bhp_well_source_evaluator.hpp>
#include <mpmc/well_discretization_petsc/fixed_total_molar_rate_control.hpp>
#include <mpmc/well_discretization_petsc/fixed_total_molar_rate_timestep_driver.hpp>
#include <mpmc/thermodynamics/pr_parameters.hpp>

#include <petscsys.h>
#include <petscvec.h>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <iostream>
#include <limits>
#include <memory>
#include <optional>
#include <span>
#include <stdexcept>
#include <string>
#include <string_view>
#include <type_traits>
#include <utility>
#include <variant>
#include <vector>

namespace {

namespace sim = mpmc::simulation;
namespace simp = mpmc::simulation_petsc;

#include "mixed_cardinality/fixture.inc"
#include "mixed_cardinality/thermodynamic_adapter.inc"
#include "mixed_cardinality/phase_transition.inc"
#include "mixed_cardinality/well_fixture.inc"
#include "mixed_cardinality/fixed_bhp.inc"
#include "mixed_cardinality/rate_control.inc"
#include "mixed_cardinality/fixed_bhp_transition.inc"

[[nodiscard]] bool same_time(
    double first,
    double second) noexcept {
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
        64.0 *
            std::numeric_limits<double>::epsilon() *
            scale;
}

struct TerminalRejectionScanContext {
    ControllerFixture* fixture{};
    const fdp::AcceptedPhysicalTimeClock3D* clock{};
    const sim::SimulationCursor* cursor{};
    double entry_time_seconds{};
    double entry_next_timestep_seconds{};
    std::size_t entry_step_count{};
    std::size_t entry_consumed_boundaries{};
    std::size_t calls{};
    bool all_attempts_observed_uncommitted{true};
};

PetscErrorCode observe_terminal_indeterminate_scan(
    const fdp::PhaseTransitionRebuiltNaturalVariableSystem3D& system,
    Vec converged_state,
    const fdp::VariableCardinalityNaturalVariableSnesSolveReport3D&
        solve_report,
    void* raw_context,
    fdp::PostSnesPhaseTransitionScanStatus3D* scan_status,
    std::vector<fdp::PostSnesPhaseTransitionProposal3D>* proposals) {
    if (raw_context == nullptr ||
        scan_status == nullptr ||
        proposals == nullptr) {
        return PETSC_ERR_ARG_NULL;
    }
    auto* context =
        static_cast<TerminalRejectionScanContext*>(
            raw_context);
    if (context->fixture == nullptr ||
        context->clock == nullptr ||
        context->cursor == nullptr) {
        return PETSC_ERR_ARG_NULL;
    }

    const auto next_boundary =
        context->cursor->next_hard_boundary();
    const bool unchanged =
        same_time(
            context->clock->accepted_time_seconds(),
            context->entry_time_seconds) &&
        context->clock->accepted_step_count() ==
            context->entry_step_count &&
        same_time(
            context->clock->next_timestep_seconds(),
            context->entry_next_timestep_seconds) &&
        context->cursor->consumed_boundary_count() ==
            context->entry_consumed_boundaries &&
        next_boundary.has_value() &&
        next_boundary->index ==
            context->entry_consumed_boundaries &&
        same_time(
            next_boundary->time_seconds,
            1.0);
    context->all_attempts_observed_uncommitted =
        context->all_attempts_observed_uncommitted &&
        unchanged;
    ++context->calls;

    return controller_scan(
        system,
        converged_state,
        solve_report,
        context->fixture,
        scan_status,
        proposals);
}

void run_terminal_rejection_transaction() {
    int rank = -1;
    int size = -1;
    if (MPI_Comm_rank(
            PETSC_COMM_WORLD,
            &rank) != MPI_SUCCESS ||
        MPI_Comm_size(
            PETSC_COMM_WORLD,
            &size) != MPI_SUCCESS) {
        throw std::runtime_error(
            "failed to query MPI rank/size for simulation terminal rejection regression");
    }
    require_collective(
        size == 2,
        "simulation terminal rejection regression requires two MPI ranks");

    const auto partition =
        make_partition(rank);
    const auto schedule =
        make_schedule(rank, false);
    const auto bridge =
        make_cell_bridge(rank);
    const auto pattern =
        make_cell_pattern(rank);

    DispatchAudit audit;
    const auto parameters =
        thermodynamic_adapter_pr_parameters();
    const auto pr_model =
        th::Pr76Phase<double>::
            from_parameters(parameters);
    auto provider =
        make_outer_rebuild_pr_provider(
            pr_model);

    ControllerFixture fixture{
        rank,
        &schedule,
        &partition,
        &bridge,
        &pattern,
        &audit,
        &pr_model,
        &provider,
        false,
        0U,
        true,
        false};

    auto accepted_system =
        make_controller_initial_system(
            &fixture);
    auto* const entry_system_identity =
        accepted_system.get();

    bool history_matches_before = false;
    PetscErrorCode error =
        accepted_system
            ->accepted_history_matches_state(
                accepted_system->initial_state(),
                &history_matches_before);
    require_collective(
        error == PETSC_SUCCESS &&
            history_matches_before,
        "terminal rejection fixture did not begin from accepted history");

    Vec entry_state = nullptr;
    error =
        VecDuplicate(
            accepted_system->initial_state(),
            &entry_state);
    if (error == PETSC_SUCCESS) {
        error =
            VecCopy(
                accepted_system->initial_state(),
                entry_state);
    }
    require_collective(
        error == PETSC_SUCCESS &&
            entry_state != nullptr,
        "failed to snapshot terminal rejection accepted reservoir state");

    fdp::AcceptedPhysicalTimeClock3D clock{
        0.0,
        1.0};
    sim::SimulationCursor cursor{
        sim::SimulationTimeline{{1.0}},
        clock.accepted_time_seconds()};

    TerminalRejectionScanContext scan_context{
        &fixture,
        &clock,
        &cursor,
        clock.accepted_time_seconds(),
        clock.next_timestep_seconds(),
        clock.accepted_step_count(),
        cursor.consumed_boundary_count()};

    fdp::PostSnesPhaseTransitionControllerBindings3D
        transition_bindings{
            &observe_terminal_indeterminate_scan,
            &scan_context,
            &controller_rebuild,
            &fixture};

    simp::SimulationRunnerOptions3D options;
    options.physical_timestep.adaptive
        .minimum_timestep_seconds =
        0.25;
    options.physical_timestep.adaptive
        .maximum_timestep_seconds =
        2.0;
    options.physical_timestep.adaptive
        .cutback_factor =
        0.5;
    options.physical_timestep.adaptive
        .growth_factor =
        2.0;
    options.physical_timestep.adaptive
        .maximum_retries =
        1U;
    options.physical_timestep.phase_transition
        .max_transition_restarts =
        4U;

    std::optional<sim::SimulationReport> report;
    error =
        simp::advance_simulation_timeline_3d(
            PETSC_COMM_WORLD,
            &accepted_system,
            transition_bindings,
            options,
            &clock,
            &cursor,
            &report);

    PetscBool same_state = PETSC_FALSE;
    const PetscErrorCode compare_error =
        VecEqual(
            entry_state,
            accepted_system->initial_state(),
            &same_state);
    bool history_matches_after = false;
    const PetscErrorCode history_error =
        accepted_system
            ->accepted_history_matches_state(
                accepted_system->initial_state(),
                &history_matches_after);

    const auto next_boundary =
        cursor.next_hard_boundary();
    require_collective(
        error == PETSC_SUCCESS &&
            report.has_value() &&
            report->status() ==
                sim::SimulationRunStatus::
                    physical_timestep_rejected &&
            report->terminal() &&
            !report->lower_layer_error_code()
                 .has_value() &&
            report->diagnostics()
                    .physical_timestep_calls ==
                1U &&
            report->diagnostics()
                    .accepted_steps ==
                0U &&
            report->diagnostics()
                    .reached_boundaries ==
                0U &&
            report->diagnostics()
                    .accepted_retries ==
                0U &&
            report->diagnostics()
                    .accepted_phase_transition_restarts ==
                0U &&
            !report->last_accepted_step()
                 .has_value() &&
            !report->last_reached_boundary()
                 .has_value() &&
            scan_context.calls ==
                2U &&
            scan_context
                .all_attempts_observed_uncommitted &&
            accepted_system.get() ==
                entry_system_identity &&
            same_time(
                accepted_system->time_step_seconds(),
                1.0) &&
            same_time(
                clock.accepted_time_seconds(),
                0.0) &&
            clock.accepted_step_count() ==
                0U &&
            same_time(
                clock.next_timestep_seconds(),
                1.0) &&
            !cursor.complete() &&
            cursor.consumed_boundary_count() ==
                0U &&
            next_boundary.has_value() &&
            next_boundary->index ==
                0U &&
            same_time(
                next_boundary->time_seconds,
                1.0) &&
            compare_error ==
                PETSC_SUCCESS &&
            same_state == PETSC_TRUE &&
            history_error ==
                PETSC_SUCCESS &&
            history_matches_after,
        "terminal retry-budget exhaustion changed accepted reservoir/history/clock/cursor or fabricated a boundary event");

    require_collective(
        VecDestroy(&entry_state) ==
            PETSC_SUCCESS,
        "terminal rejection accepted-state snapshot cleanup failed");
}

} // namespace

int main(int argc, char** argv) {
    PetscErrorCode error =
        PetscInitialize(
            &argc,
            &argv,
            nullptr,
            nullptr);
    if (error != PETSC_SUCCESS) {
        return static_cast<int>(error);
    }

    int result = 0;
    try {
        run_terminal_rejection_transaction();
        int rank = -1;
        if (MPI_Comm_rank(
                PETSC_COMM_WORLD,
                &rank) == MPI_SUCCESS &&
            rank == 0) {
            std::cout
                << "[PASS] simulation.t5.terminal-rejection-transaction\n";
        }
    } catch (const std::exception& exception) {
        int rank = -1;
        (void)MPI_Comm_rank(
            PETSC_COMM_WORLD,
            &rank);
        std::cerr
            << "[rank "
            << rank
            << "] [FAIL] "
            << exception.what()
            << '\n';
        result = 1;
    }

    error = PetscFinalize();
    if (error != PETSC_SUCCESS) {
        return static_cast<int>(error);
    }
    return result;
}
