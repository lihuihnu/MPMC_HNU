#include <mpmc/flash/pr76_pt_continuation.hpp>

#include "../pr76_three_phase/sour_gas_fixture.hpp"

#include <array>
#include <cstddef>
#include <iostream>
#include <vector>

namespace fl = mpmc::flash;
namespace sg = pr76_sour_gas_test;

std::vector<double> feed_at(double z_co2) {
    auto composition = sg::feed();
    const double scale = (1.0 - z_co2) / (1.0 - composition[0]);
    composition[0] = z_co2;
    for (std::size_t i = 1U; i < composition.size(); ++i) {
        composition[i] *= scale;
    }
    return composition;
}

int main() {
    constexpr double z_co2 = 0.65;
    const auto feed = feed_at(z_co2);
    const auto model = sg::model();
    fl::Pr76VleEvaluator evaluator(model);
    fl::Pr76PtContinuationOptions options;
    options.point_options = sg::max3_options();
    options.fallback_initial_starts = sg::starts();
    options.fallback_final_starts = sg::starts();
    const std::array<fl::Pr76PtPathState, 5> path{{
        {2.00e6, 178.8}, {2.50e6, 178.8}, {3.02e6, 178.8},
        {3.50e6, 178.8}, {4.00e6, 178.8}}};

    const auto result = fl::solve_pr76_pt_continuation(path, feed, evaluator, options);
    std::cout << "all_points_accepted=" << result.all_points_accepted
              << " brackets=" << result.transition_brackets.size() << '\n';
    for (std::size_t i = 0U; i < result.points.size(); ++i) {
        const auto& point = result.points[i];
        std::cout << "index=" << i
                  << " Pbar=" << point.state.pressure_pa / 1.0e5
                  << " phase_count=";
        if (point.accepted_phase_count) {
            std::cout << *point.accepted_phase_count;
        } else {
            std::cout << "unresolved";
        }
        std::cout << " max3_status=" << static_cast<int>(point.solve.status)
                  << " base_status=" << static_cast<int>(point.solve.base.solution.status)
                  << " incoming_hint=" << static_cast<int>(point.incoming_hint)
                  << " carried_3p=" << point.carried_three_phase_start
                  << " consumed_3p=" << point.carried_three_phase_start_consumed
                  << " attempts=" << point.solve.attempts.size()
                  << " diag=" << point.solve.diagnostic << '\n';
    }
    for (const auto& bracket : result.transition_brackets) {
        std::cout << "bracket=" << bracket.left_index << "->" << bracket.right_index
                  << " counts=" << bracket.left_phase_count << "->"
                  << bracket.right_phase_count << '\n';
    }

    const auto fallback = sg::starts();
    const auto cold = fl::solve_pr76_pt_max3(
        4.00e6, 178.8, feed, evaluator, {}, fallback, fallback);
    std::cout << "cold_40bar max3_status=" << static_cast<int>(cold.status)
              << " base_status=" << static_cast<int>(cold.base.solution.status)
              << " attempts=" << cold.attempts.size()
              << " diag=" << cold.diagnostic << '\n';
    return 0;
}
