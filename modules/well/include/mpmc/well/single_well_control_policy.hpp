#ifndef MPMC_WELL_SINGLE_WELL_CONTROL_POLICY_HPP
#define MPMC_WELL_SINGLE_WELL_CONTROL_POLICY_HPP

#include <cmath>
#include <optional>
#include <stdexcept>
#include <string_view>

namespace mpmc::well {

inline constexpr std::string_view
    single_well_control_policy_convention =
        "well/single-well-control-policy/fixed-total-molar-rate-minimum-bhp/v1";

enum class SingleWellControlMode {
    fixed_total_molar_rate,
    minimum_bottom_hole_pressure
};

enum class SingleWellControlPolicyValidationStatus {
    valid,
    inconsistent_configuration,
    value_out_of_range
};

struct SingleWellMinimumBhpConstraint {
    double minimum_bottom_hole_pressure_pa{};
    std::optional<double>
        release_rate_margin_mol_per_s;
    std::optional<double>
        release_pressure_margin_pa;

    [[nodiscard]] bool
    reactivation_enabled() const noexcept {
        return release_rate_margin_mol_per_s
                   .has_value() &&
            release_pressure_margin_pa
                .has_value();
    }
};

struct SingleWellControlPolicy {
    double target_total_molar_rate_mol_per_s{};
    std::optional<
        SingleWellMinimumBhpConstraint>
        minimum_bhp_constraint;

    [[nodiscard]] bool
    minimum_bhp_reactivation_enabled()
        const noexcept {
        return minimum_bhp_constraint
                   .has_value() &&
            minimum_bhp_constraint
                ->reactivation_enabled();
    }
};

[[nodiscard]] inline
SingleWellControlPolicyValidationStatus
validate_single_well_control_policy(
    double target_total_molar_rate_mol_per_s,
    std::optional<double>
        minimum_bottom_hole_pressure_pa,
    std::optional<double>
        minimum_bhp_release_rate_margin_mol_per_s,
    std::optional<double>
        minimum_bhp_release_pressure_margin_pa) noexcept {
    if (!std::isfinite(
            target_total_molar_rate_mol_per_s) ||
        !(target_total_molar_rate_mol_per_s >
          0.0)) {
        return
            SingleWellControlPolicyValidationStatus::
                value_out_of_range;
    }
    if (minimum_bottom_hole_pressure_pa
            .has_value() &&
        (!std::isfinite(
             *minimum_bottom_hole_pressure_pa) ||
         !(*minimum_bottom_hole_pressure_pa >
           0.0))) {
        return
            SingleWellControlPolicyValidationStatus::
                value_out_of_range;
    }

    const bool rate_margin_configured =
        minimum_bhp_release_rate_margin_mol_per_s
            .has_value();
    const bool pressure_margin_configured =
        minimum_bhp_release_pressure_margin_pa
            .has_value();
    if (rate_margin_configured !=
            pressure_margin_configured ||
        (rate_margin_configured &&
         !minimum_bottom_hole_pressure_pa
              .has_value())) {
        return
            SingleWellControlPolicyValidationStatus::
                inconsistent_configuration;
    }
    if (rate_margin_configured &&
        (!std::isfinite(
             *minimum_bhp_release_rate_margin_mol_per_s) ||
         !(*minimum_bhp_release_rate_margin_mol_per_s >
           0.0) ||
         !std::isfinite(
             *minimum_bhp_release_pressure_margin_pa) ||
         !(*minimum_bhp_release_pressure_margin_pa >
           0.0))) {
        return
            SingleWellControlPolicyValidationStatus::
                value_out_of_range;
    }
    return
        SingleWellControlPolicyValidationStatus::
            valid;
}

[[nodiscard]] inline
SingleWellControlPolicy
make_single_well_control_policy(
    double target_total_molar_rate_mol_per_s,
    std::optional<double>
        minimum_bottom_hole_pressure_pa =
            std::nullopt,
    std::optional<double>
        minimum_bhp_release_rate_margin_mol_per_s =
            std::nullopt,
    std::optional<double>
        minimum_bhp_release_pressure_margin_pa =
            std::nullopt) {
    if (validate_single_well_control_policy(
            target_total_molar_rate_mol_per_s,
            minimum_bottom_hole_pressure_pa,
            minimum_bhp_release_rate_margin_mol_per_s,
            minimum_bhp_release_pressure_margin_pa) !=
        SingleWellControlPolicyValidationStatus::
            valid) {
        throw std::invalid_argument(
            "mpmc::well: invalid single-well control policy");
    }

    SingleWellControlPolicy policy;
    policy.target_total_molar_rate_mol_per_s =
        target_total_molar_rate_mol_per_s;
    if (minimum_bottom_hole_pressure_pa
            .has_value()) {
        policy.minimum_bhp_constraint =
            SingleWellMinimumBhpConstraint{
                *minimum_bottom_hole_pressure_pa,
                minimum_bhp_release_rate_margin_mol_per_s,
                minimum_bhp_release_pressure_margin_pa};
    }
    return policy;
}

struct AcceptedSingleWellControlState {
    SingleWellControlMode control{
        SingleWellControlMode::
            fixed_total_molar_rate};
    double bottom_hole_pressure_pa{};

    [[nodiscard]] static
    AcceptedSingleWellControlState
    fixed_total_molar_rate(
        double bottom_hole_pressure_pa) {
        if (!std::isfinite(
                bottom_hole_pressure_pa) ||
            !(bottom_hole_pressure_pa > 0.0)) {
            throw std::invalid_argument(
                "mpmc::well: accepted rate-control BHP must be finite and positive");
        }
        return {
            SingleWellControlMode::
                fixed_total_molar_rate,
            bottom_hole_pressure_pa};
    }

    [[nodiscard]] static
    AcceptedSingleWellControlState
    minimum_bottom_hole_pressure(
        double bottom_hole_pressure_pa) {
        if (!std::isfinite(
                bottom_hole_pressure_pa) ||
            !(bottom_hole_pressure_pa > 0.0)) {
            throw std::invalid_argument(
                "mpmc::well: accepted minimum-BHP must be finite and positive");
        }
        return {
            SingleWellControlMode::
                minimum_bottom_hole_pressure,
            bottom_hole_pressure_pa};
    }

    [[nodiscard]] bool valid() const noexcept {
        return std::isfinite(
                   bottom_hole_pressure_pa) &&
            bottom_hole_pressure_pa > 0.0;
    }
};

enum class SingleWellRateCandidateDecision {
    accept_fixed_total_molar_rate,
    solve_minimum_bottom_hole_pressure
};

enum class SingleWellMinimumBhpProbeDecision {
    accept_minimum_bottom_hole_pressure,
    probe_fixed_total_molar_rate
};

enum class SingleWellRateReactivationDecision {
    keep_minimum_bottom_hole_pressure,
    accept_fixed_total_molar_rate
};

[[nodiscard]] inline
SingleWellRateCandidateDecision
decide_single_well_rate_candidate(
    const SingleWellControlPolicy& policy,
    double candidate_bottom_hole_pressure_pa) {
    if (!std::isfinite(
            candidate_bottom_hole_pressure_pa) ||
        !(candidate_bottom_hole_pressure_pa >
          0.0)) {
        throw std::invalid_argument(
            "mpmc::well: rate candidate BHP must be finite and positive");
    }
    if (!policy.minimum_bhp_constraint
             .has_value()) {
        return
            SingleWellRateCandidateDecision::
                accept_fixed_total_molar_rate;
    }
    return candidate_bottom_hole_pressure_pa <
            policy.minimum_bhp_constraint
                ->minimum_bottom_hole_pressure_pa
        ? SingleWellRateCandidateDecision::
              solve_minimum_bottom_hole_pressure
        : SingleWellRateCandidateDecision::
              accept_fixed_total_molar_rate;
}

[[nodiscard]] inline
SingleWellMinimumBhpProbeDecision
decide_single_well_minimum_bhp_probe(
    const SingleWellControlPolicy& policy,
    double probe_total_molar_rate_mol_per_s) {
    if (!std::isfinite(
            probe_total_molar_rate_mol_per_s)) {
        throw std::invalid_argument(
            "mpmc::well: minimum-BHP probe rate must be finite");
    }
    if (!policy.minimum_bhp_reactivation_enabled()) {
        return
            SingleWellMinimumBhpProbeDecision::
                accept_minimum_bottom_hole_pressure;
    }

    const double threshold =
        policy.target_total_molar_rate_mol_per_s +
        *policy.minimum_bhp_constraint
             ->release_rate_margin_mol_per_s;
    if (!std::isfinite(threshold)) {
        throw std::overflow_error(
            "mpmc::well: minimum-BHP release-rate threshold overflow");
    }
    return probe_total_molar_rate_mol_per_s <
            threshold
        ? SingleWellMinimumBhpProbeDecision::
              accept_minimum_bottom_hole_pressure
        : SingleWellMinimumBhpProbeDecision::
              probe_fixed_total_molar_rate;
}

[[nodiscard]] inline
SingleWellRateReactivationDecision
decide_single_well_rate_reactivation_candidate(
    const SingleWellControlPolicy& policy,
    double candidate_bottom_hole_pressure_pa) {
    if (!std::isfinite(
            candidate_bottom_hole_pressure_pa) ||
        !(candidate_bottom_hole_pressure_pa >
          0.0)) {
        throw std::invalid_argument(
            "mpmc::well: rate-reactivation BHP must be finite and positive");
    }
    if (!policy.minimum_bhp_reactivation_enabled()) {
        return
            SingleWellRateReactivationDecision::
                keep_minimum_bottom_hole_pressure;
    }

    const double threshold =
        policy.minimum_bhp_constraint
                ->minimum_bottom_hole_pressure_pa +
        *policy.minimum_bhp_constraint
             ->release_pressure_margin_pa;
    if (!std::isfinite(threshold)) {
        throw std::overflow_error(
            "mpmc::well: minimum-BHP release-pressure threshold overflow");
    }
    return candidate_bottom_hole_pressure_pa <
            threshold
        ? SingleWellRateReactivationDecision::
              keep_minimum_bottom_hole_pressure
        : SingleWellRateReactivationDecision::
              accept_fixed_total_molar_rate;
}

[[nodiscard]] inline
AcceptedSingleWellControlState
make_accepted_single_well_control_state(
    const SingleWellControlPolicy& policy,
    SingleWellControlMode accepted_control,
    double accepted_bottom_hole_pressure_pa) {
    if (accepted_control ==
        SingleWellControlMode::
            minimum_bottom_hole_pressure) {
        if (!policy.minimum_bhp_constraint
                 .has_value()) {
            throw std::invalid_argument(
                "mpmc::well: minimum-BHP accepted state requires a configured constraint");
        }
        return
            AcceptedSingleWellControlState::
                minimum_bottom_hole_pressure(
                    policy.minimum_bhp_constraint
                        ->minimum_bottom_hole_pressure_pa);
    }
    return
        AcceptedSingleWellControlState::
            fixed_total_molar_rate(
                accepted_bottom_hole_pressure_pa);
}

} // namespace mpmc::well

#endif // MPMC_WELL_SINGLE_WELL_CONTROL_POLICY_HPP
