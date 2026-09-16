#include <mpmc/flash/pr76_three_phase.hpp>

#include "synthetic_fixture.hpp"

#include <array>
#include <cmath>
#include <iostream>
#include <optional>
#include <stdexcept>
#include <string_view>
#include <vector>

namespace {
namespace fl = mpmc::flash;
using pr76_max3_test::Vec;

void require(bool value, const char* message) {
    if (!value) { throw std::runtime_error(message); }
}

Vec log_ratio(const Vec& numerator, const Vec& denominator) {
    require(numerator.size() == denominator.size(), "log-ratio shape mismatch");
    Vec value(numerator.size(), 0.0);
    for (std::size_t i = 0; i < value.size(); ++i) {
        value[i] = std::log(numerator[i]) - std::log(denominator[i]);
    }
    return value;
}

fl::Pr76PtThreePhaseStart reference_start(const fl::PtThreePhaseState& state) {
    fl::Pr76PtThreePhaseStart start;
    for (std::size_t phase = 0; phase < 3U; ++phase) {
        start.compositions[phase] = state.phases[phase].composition;
    }
    start.phase_fraction_seed = {
        state.phases[1].mole_phase_fraction,
        state.phases[2].mole_phase_fraction};
    return start;
}

fl::Pr76PtThreePhaseStart exact_fixture_start() {
    fl::Pr76PtThreePhaseStart start;
    const auto& phases = pr76_max3_test::reference_phases();
    start.compositions = {phases[0], phases[1], phases[2]};
    return start;
}

void automatic_cold_start_is_exercised() {
    const auto model = pr76_max3_test::model();
    fl::Pr76VleEvaluator evaluator(model);
    const auto feed = pr76_max3_test::equal_feed();
    const auto vle_starts = pr76_max3_test::starts();

    fl::Pr76PtMax3Options options;
    options.max_three_phase_attempts = 16U;
    const auto result = fl::solve_pr76_pt_max3(
        1.0e6, 250.0, feed, evaluator, options, vle_starts, vle_starts);

    require(result.base.solution.status == fl::PtSplitStatus::phase_set_unstable,
            "structural cold-start fixture no longer reaches two-phase final instability");
    bool saw_automatic = false;
    for (const auto& attempt : result.attempts) {
        if (attempt.witness_trial.has_value() && !attempt.supplied_start.has_value()) {
            saw_automatic = true;
            break;
        }
    }
    require(saw_automatic,
            "negative-TPD witness did not generate an automatic three-phase attempt");
}

void continuation_quota_does_not_starve_automatic() {
    const auto model = pr76_max3_test::model();
    fl::Pr76VleEvaluator evaluator(model);
    const auto feed = pr76_max3_test::equal_feed();
    const auto vle_starts = pr76_max3_test::starts();

    fl::Pr76PtMax3Options options;
    options.max_three_phase_attempts = 1U;
    fl::Pr76PtThreePhaseStart poor_hint;
    poor_hint.compositions = {feed, feed, feed};
    options.three_phase_starts.push_back(poor_hint);

    const auto result = fl::solve_pr76_pt_max3(
        1.0e6, 250.0, feed, evaluator, options, vle_starts, vle_starts);
    require(result.base.solution.status == fl::PtSplitStatus::phase_set_unstable,
            "budget fixture no longer has final two-phase instability evidence");
    require(result.base.solution.final_stability.has_value(),
            "budget fixture lost final two-phase TPD evidence");

    std::size_t supplied_attempts = 0U;
    std::size_t automatic_attempts = 0U;
    for (const auto& attempt : result.attempts) {
        if (attempt.supplied_start.has_value()) {
            ++supplied_attempts;
            require(!attempt.witness_trial.has_value(),
                    "supplied three-phase hint was incorrectly tagged as TPD evidence");
            continue;
        }
        if (attempt.witness_trial.has_value()) {
            ++automatic_attempts;
            const std::size_t witness = *attempt.witness_trial;
            require(witness < result.base.solution.final_stability->trials.size(),
                    "automatic attempt lost its originating final-TPD trial");
            const auto& trial = result.base.solution.final_stability->trials[witness];
            require(trial.status == fl::StabilityTrialStatus::negative_tpd &&
                        trial.point.has_value(),
                    "automatic attempt did not originate from retained negative TPD evidence");
        }
    }
    require(supplied_attempts == 1U,
            "supplied-source quota did not stop independently at one attempt");
    require(automatic_attempts == 1U,
            "supplied-source quota starved or altered the one-attempt automatic quota");
    require(result.attempt_limit_reached,
            "per-source quota exhaustion was not retained in diagnostics");
}

void accepted_three_phase_continues_three_to_three() {
    const auto model = pr76_max3_test::model();
    fl::Pr76VleEvaluator evaluator(model);
    const auto feed = pr76_max3_test::equal_feed();
    const auto vle_starts = pr76_max3_test::starts();

    fl::Pr76PtMax3Options first_options;
    first_options.three_phase_starts.push_back(exact_fixture_start());
    const auto first = fl::solve_pr76_pt_max3(
        1.0e6, 250.0, feed, evaluator,
        first_options, vle_starts, vle_starts);
    require(first.status == fl::Pr76PtMax3Status::three_phase &&
                first.three_phase_candidate() != nullptr,
            "reference seeded solve did not establish the three-phase continuation base");

    fl::Pr76PtMax3Options second_options;
    second_options.three_phase_starts.push_back(
        reference_start(*first.three_phase_candidate()));
    const auto second = fl::solve_pr76_pt_max3(
        1.0e6, 250.0, feed, evaluator,
        second_options, vle_starts, vle_starts);
    require(second.status == fl::Pr76PtMax3Status::three_phase &&
                second.three_phase_candidate() != nullptr &&
                second.selected_attempt.has_value(),
            "accepted three-phase state did not survive a fresh 3->3 continuation solve");
    require(second.attempts[*second.selected_attempt].supplied_start.has_value(),
            "3->3 continuation lost supplied-start provenance");
}

void three_to_two_boundary_requires_fresh_neighbor() {
    const auto model = pr76_max3_test::model();
    fl::Pr76VleEvaluator evaluator(model);
    const auto& phases = pr76_max3_test::reference_phases();
    constexpr double beta1 = 0.4;
    Vec feed(3U, 0.0);
    for (std::size_t i = 0; i < feed.size(); ++i) {
        feed[i] = (1.0 - beta1) * phases[0][i] + beta1 * phases[1][i];
    }

    fl::Pr76ThreePhaseEvaluator three_phase_evaluator(
        model,
        {fl::Pr76RootSide::lower_admissible,
         fl::Pr76RootSide::upper_admissible,
         fl::Pr76RootSide::lower_admissible});
    const auto boundary = fl::iterate_pt_three_phase(
        1.0e6, 250.0, feed,
        log_ratio(phases[1], phases[0]),
        log_ratio(phases[2], phases[0]),
        {beta1, 0.0}, three_phase_evaluator);
    require(boundary.status == fl::PtThreePhaseStatus::phase_disappearance &&
                boundary.disappearing_phase.has_value() &&
                *boundary.disappearing_phase == 2U && boundary.equations_converged(),
            "3->2 continuation edge did not retain explicit disappearance evidence");

    const std::vector<Vec> surviving{
        boundary.point->phases[0].composition,
        boundary.point->phases[1].composition};
    const auto neighbor = fl::solve_pr76_pt_vle(
        1.0e6, 250.0, feed, evaluator, {}, surviving, surviving);
    require(neighbor.solution.status ==
                fl::PtSplitStatus::two_phase_no_instability_found &&
                neighbor.solution.candidate() != nullptr &&
                neighbor.solution.final_stability.has_value() &&
                neighbor.solution.final_stability->status ==
                    fl::StabilityStatus::no_instability_found,
            "3->2 boundary was not independently fresh-resolved as a stable two-phase neighbor");
}

using Test = std::pair<std::string_view, void (*)()>;
constexpr Test tests[]{
    {"automatic_cold_start", automatic_cold_start_is_exercised},
    {"continuation_budget", continuation_quota_does_not_starve_automatic},
    {"continuation_3_to_3", accepted_three_phase_continues_three_to_three},
    {"continuation_3_to_2", three_to_two_boundary_requires_fresh_neighbor}};

} // namespace

int main(int argc, char** argv) {
    try {
        if (argc != 2) { throw std::invalid_argument("one test name required"); }
        for (const auto& [name, run] : tests) {
            if (name == argv[1]) {
                run();
                std::cout << "[PASS] " << name << '\n';
                return 0;
            }
        }
        throw std::invalid_argument("unknown test name");
    } catch (const std::exception& error) {
        std::cerr << "[FAIL] " << error.what() << '\n';
        return 1;
    }
}
