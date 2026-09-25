#ifndef MPMC_FLOW_DISCRETIZATION_PETSC_PHYSICAL_TIME_LOOP_HPP
#define MPMC_FLOW_DISCRETIZATION_PETSC_PHYSICAL_TIME_LOOP_HPP

#include <mpmc/flow_discretization_petsc/physical_timestep_driver.hpp>

#include <petscsys.h>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <memory>
#include <optional>
#include <string_view>
#include <utility>
#include <vector>

namespace mpmc::flow_discretization_petsc {

inline constexpr std::string_view
    physical_time_loop_convention =
        "flow_discretization_petsc/physical-time-loop/v1";

enum class PhysicalTimeLoopOutcome3D {
    target_time_reached,
    timestep_rejected
};

struct PhysicalTimeLoopReport3D {
    PhysicalTimeLoopOutcome3D outcome{
        PhysicalTimeLoopOutcome3D::
            timestep_rejected};
    double initial_time_seconds{};
    double target_time_seconds{};
    std::vector<PhysicalTimestepDriverReport3D>
        timesteps;

    [[nodiscard]] bool
    target_reached() const noexcept {
        return outcome ==
            PhysicalTimeLoopOutcome3D::
                target_time_reached;
    }
};

namespace physical_time_loop_detail {

[[nodiscard]] inline double
time_tolerance(
    double first,
    double second) noexcept {
    return 64.0 *
        std::numeric_limits<double>::epsilon() *
        std::max(
            {1.0,
             std::abs(first),
             std::abs(second)});
}

[[nodiscard]] inline bool
same_time(
    double first,
    double second) noexcept {
    return std::isfinite(first) &&
        std::isfinite(second) &&
        std::abs(first - second) <=
            time_tolerance(
                first,
                second);
}

} // namespace physical_time_loop_detail

/// Advance accepted physical time to target_time_seconds by repeatedly invoking
/// the production one-timestep driver.
///
/// The outer loop never commits state/history itself. Before each call it caps
/// only the driver's first trial dt by the remaining interval, so an accepted
/// step cannot materially overshoot the target. If the terminal remainder is
/// smaller than the configured adaptive minimum, that exact remainder is
/// permitted as the final first attempt; it becomes the retry floor for that
/// step, so a failed terminal clip is rejected rather than cut below the
/// schedule boundary.
///
/// Any terminal timestep rejection stops the interval immediately and is
/// reported without advancing that failed step. The one-step driver remains the
/// sole owner of rollback, phase-transition restart and accepted commit.
inline PetscErrorCode
advance_physical_time_to_3d(
    MPI_Comm comm,
    std::unique_ptr<
        PhaseTransitionRebuiltNaturalVariableSystem3D>*
            accepted_system,
    std::size_t
        first_step_initial_phase_transition_restarts,
    PostSnesPhaseTransitionControllerBindings3D
        transition_bindings,
    const PhysicalTimestepDriverOptions3D&
        timestep_options,
    double target_time_seconds,
    AcceptedPhysicalTimeClock3D* clock,
    std::optional<
        PhysicalTimeLoopReport3D>* report) {
    using namespace
        physical_time_loop_detail;

    if (accepted_system == nullptr ||
        clock == nullptr ||
        report == nullptr) {
        return PETSC_ERR_ARG_NULL;
    }
    report->reset();
    if (*accepted_system == nullptr ||
        !std::isfinite(target_time_seconds)) {
        return PETSC_ERR_ARG_INCOMP;
    }

    const double initial_time_seconds =
        clock->accepted_time_seconds();
    if (!std::isfinite(
            initial_time_seconds) ||
        target_time_seconds <
            initial_time_seconds -
                time_tolerance(
                    target_time_seconds,
                    initial_time_seconds)) {
        return PETSC_ERR_ARG_OUTOFRANGE;
    }

    PhysicalTimeLoopReport3D
        completed;
    completed.initial_time_seconds =
        initial_time_seconds;
    completed.target_time_seconds =
        target_time_seconds;

    if (same_time(
            initial_time_seconds,
            target_time_seconds)) {
        completed.outcome =
            PhysicalTimeLoopOutcome3D::
                target_time_reached;
        report->emplace(
            std::move(completed));
        return PETSC_SUCCESS;
    }

    bool first_step = true;
    for (;;) {
        const double time_before =
            clock->accepted_time_seconds();
        const double remaining =
            target_time_seconds -
            time_before;
        if (!std::isfinite(remaining) ||
            !(remaining > 0.0)) {
            if (same_time(
                    time_before,
                    target_time_seconds)) {
                completed.outcome =
                    PhysicalTimeLoopOutcome3D::
                        target_time_reached;
                report->emplace(
                    std::move(completed));
                return PETSC_SUCCESS;
            }
            return PETSC_ERR_FP;
        }

        auto step_options =
            timestep_options;
        step_options
            .initial_timestep_cap_seconds =
            remaining;
        if (remaining <
            step_options.adaptive
                .minimum_timestep_seconds) {
            step_options.adaptive
                .minimum_timestep_seconds =
                remaining;
        }

        std::optional<
            PhysicalTimestepDriverReport3D>
            step_report;
        const PetscErrorCode error =
            advance_one_physical_timestep_3d(
                comm,
                accepted_system,
                first_step
                    ? first_step_initial_phase_transition_restarts
                    : 0U,
                transition_bindings,
                step_options,
                clock,
                &step_report);
        first_step = false;
        if (error != PETSC_SUCCESS) {
            return error;
        }
        if (!step_report.has_value()) {
            return PETSC_ERR_PLIB;
        }

        const bool accepted =
            step_report->accepted();
        completed.timesteps.push_back(
            std::move(
                *step_report));

        if (!accepted) {
            completed.outcome =
                PhysicalTimeLoopOutcome3D::
                    timestep_rejected;
            report->emplace(
                std::move(completed));
            return PETSC_SUCCESS;
        }

        const double time_after =
            clock->accepted_time_seconds();
        const double tolerance =
            time_tolerance(
                time_after,
                target_time_seconds);
        if (!std::isfinite(time_after) ||
            !(time_after > time_before) ||
            time_after >
                target_time_seconds +
                    tolerance) {
            return PETSC_ERR_PLIB;
        }

        if (same_time(
                time_after,
                target_time_seconds)) {
            completed.outcome =
                PhysicalTimeLoopOutcome3D::
                    target_time_reached;
            report->emplace(
                std::move(completed));
            return PETSC_SUCCESS;
        }
    }
}

} // namespace mpmc::flow_discretization_petsc

#endif // MPMC_FLOW_DISCRETIZATION_PETSC_PHYSICAL_TIME_LOOP_HPP
