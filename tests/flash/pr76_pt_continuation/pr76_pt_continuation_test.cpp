#include <mpmc/flash/pr76_pt_continuation.hpp>

#include "../pr76_three_phase/synthetic_fixture.hpp"

#include <algorithm>
#include <array>
#include <stdexcept>
#include <string_view>

namespace {
namespace fl = mpmc::flash;

void require(bool value, const char* message) {
    if (!value) { throw std::runtime_error(message); }
}

fl::Pr76PtThreePhaseStart exact_start() {
    fl::Pr76PtThreePhaseStart start;
    const auto& phases = pr76_max3_test::reference_phases();
    start.compositions = {phases[0], phases[1], phases[2]};
    return start;
}

fl::Pr76PtContinuationOptions structural_options() {
    fl::Pr76PtContinuationOptions options;
    options.point_options.three_phase_starts.push_back(exact_start());
    options.fallback_initial_starts = pr76_max3_test::starts();
    options.fallback_final_starts = pr76_max3_test::starts();
    return options;
}

void repeated_three_phase_is_fresh_continuation() {
    const auto model = pr76_max3_test::model();
    fl::Pr76VleEvaluator evaluator(model);
    const std::array<fl::Pr76PtPathState, 3> path{{
        {1.0e6, 250.0}, {1.0e6, 250.0}, {1.0e6, 250.0}}};
    const auto result = fl::solve_pr76_pt_continuation(
        path, pr76_max3_test::equal_feed(), evaluator, structural_options());

    require(result.all_points_accepted && result.points.size() == path.size(),
            "repeated three-phase path did not accept every fresh point");
    for (const auto& point : result.points) {
        require(point.accepted_phase_count && *point.accepted_phase_count == 3U,
                "repeated structural point did not remain three-phase");
        require(point.solve.three_phase_candidate() != nullptr,
                "accepted continuation point lost owned three-phase state");
    }
    require(result.points[0].incoming_hint == fl::Pr76PtContinuationHintKind::none,
            "first point unexpectedly has continuation provenance");
    require(result.points[1].incoming_hint ==
                fl::Pr76PtContinuationHintKind::three_phase &&
                result.points[1].carried_three_phase_start &&
                result.points[1].carried_three_phase_start_consumed,
            "second point did not consume previous accepted three-phase state as hint");
    require(result.points[2].carried_three_phase_start_consumed,
            "third point did not continue from the second fresh solution");
    require(result.transition_brackets.empty(),
            "same topology created a false transition bracket");
}

void invalid_path_is_rejected_before_solve() {
    const auto model = pr76_max3_test::model();
    fl::Pr76VleEvaluator evaluator(model);
    const std::array<fl::Pr76PtPathState, 1> bad{{{-1.0, 250.0}}};
    bool threw = false;
    try {
        (void)fl::solve_pr76_pt_continuation(
            bad, pr76_max3_test::equal_feed(), evaluator);
    } catch (const std::domain_error&) {
        threw = true;
    }
    require(threw, "invalid PT path was not rejected before scientific solve");
}

void transition_bracket_contract() {
    fl::Pr76PtTransitionBracket bracket;
    bracket.left_phase_count = 2U;
    bracket.right_phase_count = 3U;
    bracket.adjacent_phase_count_step = true;
    bracket.exact_boundary_resolved = false;
    require(bracket.adjacent_phase_count_step && !bracket.exact_boundary_resolved,
            "discrete continuation bracket must not claim an exact boundary");
}

void require_fresh_state(const fl::Pr76PtContinuationPointResult& point) {
    require(point.solve.base.solution.initial_stability.pressure_pa ==
                point.state.pressure_pa &&
                point.solve.base.solution.initial_stability.temperature_k ==
                    point.state.temperature_k,
            "continuation point did not own a fresh solve at its current PT state");
}

void require_phase_sequence(
    const fl::Pr76PtContinuationResult& result,
    std::span<const std::size_t> expected) {
    require(result.all_points_accepted && result.points.size() == expected.size(),
            "ordered PT sequence contains unresolved points");
    for (std::size_t i = 0; i < expected.size(); ++i) {
        require(result.points[i].accepted_phase_count &&
                    *result.points[i].accepted_phase_count == expected[i],
                "ordered PT sequence changed accepted phase count");
        require_fresh_state(result.points[i]);
    }
    require(result.transition_brackets.size() == 2U,
            "ordered PT sequence did not expose exactly two topology brackets");
    for (const auto& bracket : result.transition_brackets) {
        require(bracket.adjacent_phase_count_step &&
                    !bracket.exact_boundary_resolved,
                "phase-count bracket invented an exact boundary or skipped a topology");
    }
}

void bidirectional_phase_sequence() {
    const auto model = pr76_max3_test::model();
    fl::Pr76VleEvaluator evaluator(model);
    const std::array<fl::Pr76PtPathState, 4> forward{{
        {1.90e6, 325.0}, // accepted 3-phase side of the unresolved 3/2 belt
        {2.80e6, 331.0}, // accepted 2-phase side of that belt
        {1.00e6, 347.0}, // accepted 2-phase side of the 2/1 temperature boundary
        {1.00e6, 348.0}  // accepted 1-phase side
    }};
    const std::array<std::size_t, 4> forward_counts{3U, 2U, 2U, 1U};
    const auto forward_result = fl::solve_pr76_pt_continuation(
        forward, pr76_max3_test::equal_feed(), evaluator, structural_options());
    require_phase_sequence(forward_result, forward_counts);
    require(forward_result.transition_brackets[0].left_phase_count == 3U &&
                forward_result.transition_brackets[0].right_phase_count == 2U &&
                forward_result.transition_brackets[1].left_phase_count == 2U &&
                forward_result.transition_brackets[1].right_phase_count == 1U,
            "forward PT scan lost 3->2->1 bracket ordering");
    require(forward_result.points[1].incoming_hint ==
                fl::Pr76PtContinuationHintKind::three_phase &&
                forward_result.points[1].carried_three_phase_start &&
                !forward_result.points[1].carried_three_phase_start_consumed,
            "3->2 point must carry, but not need to consume, the previous 3-phase hint");
    require(forward_result.points[2].incoming_hint ==
                fl::Pr76PtContinuationHintKind::two_phase &&
                !forward_result.points[2].carried_three_phase_start,
            "two-phase continuation leaked a stale three-phase hint");

    auto reverse = forward;
    std::reverse(reverse.begin(), reverse.end());
    const std::array<std::size_t, 4> reverse_counts{1U, 2U, 2U, 3U};
    const auto reverse_result = fl::solve_pr76_pt_continuation(
        reverse, pr76_max3_test::equal_feed(), evaluator, structural_options());
    require_phase_sequence(reverse_result, reverse_counts);
    require(reverse_result.transition_brackets[0].left_phase_count == 1U &&
                reverse_result.transition_brackets[0].right_phase_count == 2U &&
                reverse_result.transition_brackets[1].left_phase_count == 2U &&
                reverse_result.transition_brackets[1].right_phase_count == 3U,
            "reverse PT scan lost 1->2->3 bracket ordering");
}

using Test = std::pair<std::string_view, void (*)()>;
constexpr Test tests[]{
    {"repeated_three_phase", repeated_three_phase_is_fresh_continuation},
    {"invalid_path", invalid_path_is_rejected_before_solve},
    {"bracket_contract", transition_bracket_contract},
    {"bidirectional_phase_sequence", bidirectional_phase_sequence}};

} // namespace

int main(int argc, char** argv) {
    try {
        if (argc != 2) { throw std::invalid_argument("one test name required"); }
        for (const auto& [name, run] : tests) {
            if (name == argv[1]) {
                run();
                return 0;
            }
        }
        throw std::invalid_argument("unknown test name");
    } catch (const std::exception&) {
        return 1;
    }
}
