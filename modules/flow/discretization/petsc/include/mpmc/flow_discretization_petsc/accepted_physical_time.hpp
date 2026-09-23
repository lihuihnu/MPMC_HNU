#ifndef MPMC_FLOW_DISCRETIZATION_PETSC_ACCEPTED_PHYSICAL_TIME_HPP
#define MPMC_FLOW_DISCRETIZATION_PETSC_ACCEPTED_PHYSICAL_TIME_HPP

#include <mpmc/flow_discretization_petsc/adaptive_timestep_controller.hpp>

#include <petscvec.h>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <optional>
#include <stdexcept>
#include <string_view>

namespace mpmc::flow_discretization_petsc {

inline constexpr std::string_view
    accepted_physical_time_convention =
        "flow_discretization_petsc/accepted-physical-time/v1";

struct AcceptedPhysicalTimestepRecord3D {
    std::size_t accepted_step_index{};
    double time_n_seconds{};
    double time_np1_seconds{};
    double accepted_timestep_seconds{};
    double next_timestep_seconds{};
    std::size_t phase_transition_restarts{};
    AdaptiveTimestepDecision3D decision{
        AdaptiveTimestepDecision3D::
            accept_and_hold};
};

class AcceptedPhysicalTimeClock3D {
public:
    AcceptedPhysicalTimeClock3D(
        double initial_time_seconds,
        double initial_timestep_seconds)
        : accepted_time_seconds_(
              initial_time_seconds),
          next_timestep_seconds_(
              initial_timestep_seconds) {
        if (!std::isfinite(
                initial_time_seconds) ||
            initial_time_seconds < 0.0 ||
            !std::isfinite(
                initial_timestep_seconds) ||
            !(initial_timestep_seconds > 0.0)) {
            throw std::invalid_argument(
                "mpmc::flow_discretization_petsc: invalid accepted physical-time clock initialization");
        }
    }

    [[nodiscard]] double
    accepted_time_seconds() const noexcept {
        return accepted_time_seconds_;
    }

    [[nodiscard]] std::size_t
    accepted_step_count() const noexcept {
        return accepted_step_count_;
    }

    [[nodiscard]] double
    next_timestep_seconds() const noexcept {
        return next_timestep_seconds_;
    }

private:
    friend PetscErrorCode
    commit_accepted_physical_timestep_3d(
        PhaseTransitionRebuiltNaturalVariableSystem3D&,
        Vec,
        double,
        const AdaptiveTimestepAttemptResult3D&,
        const AdaptiveTimestepDecisionResult3D&,
        AcceptedPhysicalTimeClock3D*,
        std::optional<
            AcceptedPhysicalTimestepRecord3D>*);

    double accepted_time_seconds_{};
    std::size_t accepted_step_count_{};
    double next_timestep_seconds_{};
};

[[nodiscard]] inline PetscErrorCode
commit_accepted_physical_timestep_3d(
    PhaseTransitionRebuiltNaturalVariableSystem3D&
        accepted_system,
    Vec accepted_state,
    double accepted_timestep_seconds,
    const AdaptiveTimestepAttemptResult3D&
        attempt,
    const AdaptiveTimestepDecisionResult3D&
        decision,
    AcceptedPhysicalTimeClock3D* clock,
    std::optional<
        AcceptedPhysicalTimestepRecord3D>*
            record) {
    if (clock == nullptr ||
        record == nullptr ||
        accepted_state == nullptr) {
        return PETSC_ERR_ARG_NULL;
    }
    record->reset();

    const bool accepted_decision =
        decision.decision ==
                AdaptiveTimestepDecision3D::
                    accept_and_grow ||
        decision.decision ==
                AdaptiveTimestepDecision3D::
                    accept_and_hold;
    if (attempt.outcome !=
            AdaptiveTimestepAttemptOutcome3D::
                stable_phase_set ||
        attempt.function_domain_errors != 0 ||
        attempt.jacobian_domain_errors != 0 ||
        !accepted_decision ||
        !decision.next_timestep_seconds
             .has_value() ||
        !std::isfinite(
            accepted_timestep_seconds) ||
        !(accepted_timestep_seconds > 0.0) ||
        !std::isfinite(
            *decision
                 .next_timestep_seconds) ||
        !(*decision
               .next_timestep_seconds >
          0.0)) {
        return PETSC_ERR_ARG_INCOMP;
    }

    const double system_dt =
        accepted_system
            .time_step_seconds();
    const double dt_scale =
        std::max(
            {1.0,
             std::abs(system_dt),
             std::abs(
                 accepted_timestep_seconds)});
    if (!std::isfinite(system_dt) ||
        std::abs(
            system_dt -
            accepted_timestep_seconds) >
            64.0 *
                std::numeric_limits<double>::
                    epsilon() *
                dt_scale) {
        return PETSC_ERR_ARG_INCOMP;
    }

    const double next_time =
        clock->accepted_time_seconds_ +
        accepted_timestep_seconds;
    if (!std::isfinite(next_time) ||
        !(next_time >
          clock->accepted_time_seconds_)) {
        return PETSC_ERR_FP;
    }

    const PetscErrorCode error =
        accepted_system
            .rebase_accepted_timestep(
                accepted_state,
                *decision
                     .next_timestep_seconds);
    if (error != PETSC_SUCCESS) {
        return error;
    }

    AcceptedPhysicalTimestepRecord3D
        completed{
            clock->accepted_step_count_,
            clock->accepted_time_seconds_,
            next_time,
            accepted_timestep_seconds,
            *decision.next_timestep_seconds,
            attempt.phase_transition_restarts,
            decision.decision};

    clock->accepted_time_seconds_ =
        next_time;
    ++clock->accepted_step_count_;
    clock->next_timestep_seconds_ =
        *decision.next_timestep_seconds;
    record->emplace(
        std::move(completed));
    return PETSC_SUCCESS;
}

} // namespace mpmc::flow_discretization_petsc

#endif // MPMC_FLOW_DISCRETIZATION_PETSC_ACCEPTED_PHYSICAL_TIME_HPP
