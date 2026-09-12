#include <mpmc/flash/pr76_pt_continuation.hpp>

#include "../pr76_three_phase/synthetic_fixture.hpp"

#include <array>
#include <iostream>
#include <stdexcept>
#include <string_view>
#include <vector>

namespace {
namespace fl = mpmc::flash;
using pr76_max3_test::Vec;

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

void diagnostic_pt_grid() {
    const auto model = pr76_max3_test::model();
    fl::Pr76VleEvaluator evaluator(model);
    constexpr std::array<double, 4> temperatures{332.0, 335.0, 338.0, 340.0};
    constexpr std::array<double, 7> pressures_mpa{0.55, 0.70, 0.85, 1.00, 1.15, 1.35, 1.60};

    std::cout << "PR76_CONTINUATION_PT_GRID";
    for (const double temperature : temperatures) {
        for (const double pressure_mpa : pressures_mpa) {
            const std::array<fl::Pr76PtPathState, 1> path{{
                {pressure_mpa * 1.0e6, temperature}}};
            const auto result = fl::solve_pr76_pt_continuation(
                path, pr76_max3_test::equal_feed(), evaluator, structural_options());
            std::cout << ' ' << temperature << 'K@' << pressure_mpa << "MPa:";
            if (result.points.front().accepted_phase_count) {
                std::cout << *result.points.front().accepted_phase_count;
            } else {
                std::cout << 'X';
            }
        }
    }
    std::cout << '\n';
}

using Test = std::pair<std::string_view, void (*)()>;
constexpr Test tests[]{
    {"repeated_three_phase", repeated_three_phase_is_fresh_continuation},
    {"invalid_path", invalid_path_is_rejected_before_solve},
    {"bracket_contract", transition_bracket_contract},
    {"diagnostic_pt_grid", diagnostic_pt_grid}};

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
