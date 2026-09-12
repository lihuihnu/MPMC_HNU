#include <mpmc/flash/cpa_split.hpp>

#include "test_support.hpp"

#include <cmath>
#include <iostream>
#include <vector>

int main() {
    namespace fl = mpmc::flash;
    namespace th = mpmc::thermodynamics;
    const auto parameters = cpa_stability_test::binary(false);
    const auto model = th::CpaPtPhase::from_parameters(parameters);
    const auto& phases = cpa_split_test::coexistence_phases();
    const auto feed = cpa_split_test::feed();
    std::vector<double> log_k(feed.size(), 0.0);
    for (std::size_t i = 0; i < feed.size(); ++i) {
        log_k[i] = std::log(phases[1][i]) - std::log(phases[0][i]);
    }

    for (int mode = 0; mode < 2; ++mode) {
        auto pt = cpa_stability_test::fast_pt_options();
        if (mode == 1) {
            pt.scan_intervals = 256U;
            pt.max_evaluations = 4096U;
            pt.pressure_absolute_tolerance_pa = 1.0e-8;
            pt.pressure_relative_tolerance = 1.0e-13;
        }
        fl::CpaVleEvaluator evaluator(model, pt);
        const auto result = fl::iterate_pt_split(
            cpa_split_test::pressure_pa,
            cpa_split_test::temperature_k,
            feed, log_k, evaluator);
        std::cout << "CPA_EXACT_SPLIT mode=" << mode
                  << " status=" << static_cast<int>(result.status)
                  << " iter=" << result.iterations
                  << " eval=" << result.evaluations;
        if (result.point) {
            std::cout << " fug=" << result.point->fugacity_norm
                      << " beta=" << result.point->fractions.vapor_fraction
                      << " xA=" << result.point->fractions.liquid[0]
                      << " yA=" << result.point->fractions.vapor[0];
        }
        std::cout << '\n';
    }
    return 0;
}
