#include <mpmc/well/single_well_control_policy.hpp>

#include <iostream>
#include <limits>
#include <source_location>
#include <stdexcept>
#include <string>
#include <string_view>
#include <utility>

namespace {

namespace well = mpmc::well;

void require(
    bool condition,
    std::string_view message,
    std::source_location where =
        std::source_location::current()) {
    if (!condition) {
        throw std::runtime_error(
            std::string{where.file_name()} +
            ":" +
            std::to_string(where.line()) +
            ": " +
            std::string{message});
    }
}

template <class Function>
void expect_invalid(Function&& function) {
    bool caught = false;
    try {
        std::forward<Function>(function)();
    } catch (const std::invalid_argument&) {
        caught = true;
    }
    require(caught, "expected invalid_argument");
}

well::SingleWellControlPolicy policy_with_hysteresis() {
    return well::make_single_well_control_policy(
        10.0,
        100.0,
        2.0,
        5.0);
}

void rate_hold_and_switch() {
    const auto policy =
        policy_with_hysteresis();
    require(
        well::decide_single_well_rate_candidate(
            policy,
            101.0) ==
            well::SingleWellRateCandidateDecision::
                accept_fixed_total_molar_rate,
        "rate candidate above p_min did not remain rate controlled");
    require(
        well::decide_single_well_rate_candidate(
            policy,
            99.0) ==
            well::SingleWellRateCandidateDecision::
                solve_minimum_bottom_hole_pressure,
        "rate candidate below p_min did not switch to minimum-BHP");
}

void rate_pressure_equality_boundary() {
    const auto policy =
        policy_with_hysteresis();
    require(
        well::decide_single_well_rate_candidate(
            policy,
            100.0) ==
            well::SingleWellRateCandidateDecision::
                accept_fixed_total_molar_rate,
        "p_rate == p_min must remain rate controlled");
}

void capacity_deadband_and_equality() {
    const auto policy =
        policy_with_hysteresis();
    require(
        well::decide_single_well_minimum_bhp_probe(
            policy,
            11.999) ==
            well::SingleWellMinimumBhpProbeDecision::
                accept_minimum_bottom_hole_pressure,
        "capacity inside release deadband reactivated rate control");
    require(
        well::decide_single_well_minimum_bhp_probe(
            policy,
            12.0) ==
            well::SingleWellMinimumBhpProbeDecision::
                probe_fixed_total_molar_rate,
        "q_bhp == target + delta_q must permit a rate probe");
}

void pressure_deadband_and_equality() {
    const auto policy =
        policy_with_hysteresis();
    require(
        well::decide_single_well_rate_reactivation_candidate(
            policy,
            104.999) ==
            well::SingleWellRateReactivationDecision::
                keep_minimum_bottom_hole_pressure,
        "candidate inside pressure deadband reactivated rate control");
    require(
        well::decide_single_well_rate_reactivation_candidate(
            policy,
            105.0) ==
            well::SingleWellRateReactivationDecision::
                accept_fixed_total_molar_rate,
        "p_rate == p_min + delta_p must reactivate rate control");
}

void disabled_reactivation() {
    const auto policy =
        well::make_single_well_control_policy(
            10.0,
            100.0);
    require(
        !policy.minimum_bhp_reactivation_enabled(),
        "missing hysteresis margins unexpectedly enabled reactivation");
    require(
        well::decide_single_well_minimum_bhp_probe(
            policy,
            1000.0) ==
            well::SingleWellMinimumBhpProbeDecision::
                accept_minimum_bottom_hole_pressure,
        "disabled reactivation did not keep minimum-BHP control");
}

void accepted_state() {
    const auto policy =
        policy_with_hysteresis();
    const auto rate =
        well::make_accepted_single_well_control_state(
            policy,
            well::SingleWellControlMode::
                fixed_total_molar_rate,
            123.0);
    require(
        rate.control ==
                well::SingleWellControlMode::
                    fixed_total_molar_rate &&
            rate.bottom_hole_pressure_pa ==
                123.0 &&
            rate.valid(),
        "accepted rate-control state mismatch");

    const auto bhp =
        well::make_accepted_single_well_control_state(
            policy,
            well::SingleWellControlMode::
                minimum_bottom_hole_pressure,
            999.0);
    require(
        bhp.control ==
                well::SingleWellControlMode::
                    minimum_bottom_hole_pressure &&
            bhp.bottom_hole_pressure_pa ==
                100.0 &&
            bhp.valid(),
        "minimum-BHP accepted state did not bind the configured limit");
}

void invalid_policy() {
    using Status =
        well::SingleWellControlPolicyValidationStatus;
    require(
        well::validate_single_well_control_policy(
            0.0,
            100.0,
            std::nullopt,
            std::nullopt) ==
            Status::value_out_of_range,
        "zero target was not rejected");
    require(
        well::validate_single_well_control_policy(
            10.0,
            100.0,
            1.0,
            std::nullopt) ==
            Status::inconsistent_configuration,
        "unpaired hysteresis margins were not rejected");
    require(
        well::validate_single_well_control_policy(
            10.0,
            std::nullopt,
            1.0,
            1.0) ==
            Status::inconsistent_configuration,
        "release margins without minimum-BHP were not rejected");
    require(
        well::validate_single_well_control_policy(
            10.0,
            100.0,
            -1.0,
            1.0) ==
            Status::value_out_of_range,
        "negative release-rate margin was not rejected");
    require(
        well::validate_single_well_control_policy(
            std::numeric_limits<double>::
                quiet_NaN(),
            100.0,
            std::nullopt,
            std::nullopt) ==
            Status::value_out_of_range,
        "NaN target was not rejected");
    expect_invalid(
        [] {
            (void)well::
                make_single_well_control_policy(
                    10.0,
                    100.0,
                    1.0,
                    std::nullopt);
        });
}

using Test =
    std::pair<std::string_view, void (*)()>;

constexpr Test tests[]{
    {"rate_hold_and_switch",
     rate_hold_and_switch},
    {"rate_pressure_equality_boundary",
     rate_pressure_equality_boundary},
    {"capacity_deadband_and_equality",
     capacity_deadband_and_equality},
    {"pressure_deadband_and_equality",
     pressure_deadband_and_equality},
    {"disabled_reactivation",
     disabled_reactivation},
    {"accepted_state",
     accepted_state},
    {"invalid_policy",
     invalid_policy}};

} // namespace

int main(int argc, char** argv) {
    if (argc != 2) {
        std::cerr
            << "expected one test-case name\n";
        return 2;
    }
    const std::string_view requested{
        argv[1]};
    for (const auto& [name, test] :
         tests) {
        if (name == requested) {
            try {
                test();
                return 0;
            } catch (const std::exception& error) {
                std::cerr
                    << "[FAIL] "
                    << name
                    << ": "
                    << error.what()
                    << '\n';
                return 1;
            }
        }
    }
    std::cerr
        << "unknown test-case: "
        << requested
        << '\n';
    return 2;
}
