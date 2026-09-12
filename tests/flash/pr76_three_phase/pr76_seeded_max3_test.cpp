#include <mpmc/flash/pr76_max3_phase_set.hpp>

#include "synthetic_fixture.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <iostream>
#include <optional>
#include <stdexcept>
#include <string_view>

namespace {
namespace fl = mpmc::flash;
using pr76_max3_test::Vec;

void require(bool value, const char* message) {
    if (!value) { throw std::runtime_error(message); }
}

void require_unordered_reference_match(const fl::PtThreePhaseState& state) {
    const auto& reference = pr76_max3_test::reference_phases();
    std::array<bool, 3> used{false, false, false};
    for (const auto& phase : state.phases) {
        std::optional<std::size_t> best;
        double best_error = 1.0e100;
        for (std::size_t candidate = 0; candidate < 3U; ++candidate) {
            if (used[candidate]) { continue; }
            double error = 0.0;
            for (std::size_t i = 0; i < 3U; ++i) {
                error = std::max(
                    error,
                    std::abs(phase.composition[i] - reference[candidate][i]));
            }
            if (error < best_error) { best_error = error; best = candidate; }
        }
        require(best.has_value() && best_error < 3.0e-6,
                "seeded PR76 max3 composition left structural reference anchor");
        used[*best] = true;
    }
}

void seeded_max3() {
    const auto model = pr76_max3_test::model();
    fl::Pr76VleEvaluator evaluator(model);
    const auto feed = pr76_max3_test::equal_feed();
    const auto starts = pr76_max3_test::starts();
    const auto result = fl::solve_pr76_pt_max3(
        1.0e6, 250.0, feed, evaluator, {}, starts, starts);
    require(result.status == fl::Pr76PtMax3Status::three_phase &&
                result.three_phase_candidate() != nullptr &&
                result.selected_attempt.has_value(),
            "seeded PR76 structural fixture did not close at three phases");
    require_unordered_reference_match(*result.three_phase_candidate());
    const auto& attempt = result.attempts[*result.selected_attempt];
    require(attempt.final_stability &&
                attempt.final_stability->status ==
                    fl::StabilityStatus::no_instability_found,
            "seeded PR76 three-phase candidate lacks final stability closure");
    const auto published = fl::project_pr76_pt_max3_phase_set(result);
    require(published.solution.status == fl::PtPhaseSetStatus::accepted &&
                published.solution.accepted_phase_count() == 3U,
            "seeded PR76 max3 publication failed");
}

void seeded_component_permutation() {
    const auto first_model = pr76_max3_test::model(false);
    const auto second_model = pr76_max3_test::model(true);
    fl::Pr76VleEvaluator first_evaluator(first_model);
    fl::Pr76VleEvaluator second_evaluator(second_model);
    const auto feed = pr76_max3_test::equal_feed();
    const auto first_starts = pr76_max3_test::starts(false);
    const auto second_starts = pr76_max3_test::starts(true);
    const auto first = fl::solve_pr76_pt_max3(
        1.0e6, 250.0, feed, first_evaluator, {}, first_starts, first_starts);
    const auto second = fl::solve_pr76_pt_max3(
        1.0e6, 250.0, feed, second_evaluator, {}, second_starts, second_starts);
    require(first.three_phase_candidate() != nullptr &&
                second.three_phase_candidate() != nullptr,
            "component permutation lost seeded PR76 three-phase solution");

    std::array<Vec, 3> remapped;
    for (std::size_t p = 0; p < 3U; ++p) {
        const auto& x = second.three_phase_candidate()->phases[p].composition;
        remapped[p] = {x[0], x[2], x[1]};
    }
    std::array<bool, 3> used{false, false, false};
    for (const auto& phase : first.three_phase_candidate()->phases) {
        std::optional<std::size_t> match;
        for (std::size_t q = 0; q < 3U; ++q) {
            if (used[q]) { continue; }
            double error = 0.0;
            for (std::size_t i = 0; i < 3U; ++i) {
                error = std::max(
                    error, std::abs(phase.composition[i] - remapped[q][i]));
            }
            if (error < 3.0e-9) { match = q; break; }
        }
        require(match.has_value(),
                "seeded PR76 max3 changed under component permutation");
        used[*match] = true;
    }
}

} // namespace

int main(int argc, char** argv) {
    try {
        if (argc != 2) { throw std::invalid_argument("one test name required"); }
        const std::string_view name{argv[1]};
        if (name == "seeded_max3") seeded_max3();
        else if (name == "seeded_component_permutation") seeded_component_permutation();
        else throw std::invalid_argument("unknown test name");
        std::cout << "[PASS] " << name << '\n';
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "[FAIL] " << error.what() << '\n';
        return 1;
    }
}
