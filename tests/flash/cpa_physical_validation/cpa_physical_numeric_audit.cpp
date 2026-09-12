#include <mpmc/flash/cpa_split.hpp>

#include "test_support.hpp"

#include <array>
#include <iostream>
#include <vector>

namespace {
namespace fl = mpmc::flash;
namespace th = mpmc::thermodynamics;

fl::PtSplitOptions split_options() {
    fl::PtSplitOptions options;
    options.initial_stability.automatic_starts = false;
    options.final_stability.automatic_starts = false;
    return options;
}

void print_result(const char* prefix, double absolute_tolerance,
                  const fl::CpaPtSplitResult& result) {
    std::cout << prefix << " abs=" << absolute_tolerance
              << " initial=" << static_cast<int>(result.solution.initial_stability.status)
              << " status=" << static_cast<int>(result.solution.status)
              << " attempts=" << result.solution.attempts.size();
    if (result.solution.initial_stability.reference_issue) {
        std::cout << " ref_issue="
                  << static_cast<int>(*result.solution.initial_stability.reference_issue);
    }
    if (const auto* point = result.solution.candidate()) {
        std::cout << " x=" << point->fractions.liquid[0]
                  << " y=" << point->fractions.vapor[0]
                  << " beta=" << point->fractions.vapor_fraction
                  << " fug=" << point->fugacity_norm;
    }
    if (result.solution.final_stability) {
        std::cout << " final="
                  << static_cast<int>(result.solution.final_stability->status);
        if (result.solution.final_stability->lowest_sampled) {
            std::cout << " tpd="
                      << result.solution.final_stability->lowest_sampled->value
                      << " stat="
                      << result.solution.final_stability->lowest_sampled->stationarity;
        }
    }
    std::cout << " diag=" << result.solution.diagnostic << '\n';
}

fl::CpaPtSplitResult solve(
    const cpa_physical_test::ExperimentalVlePoint& point,
    double absolute_tolerance, bool swapped = false) {
    const auto parameters = cpa_physical_test::parameters(swapped);
    const auto model = th::CpaPtPhase::from_parameters(parameters);
    th::CpaPtOptions pt;
    pt.scan_intervals = 256U;
    pt.max_evaluations = 4096U;
    pt.pressure_absolute_tolerance_pa = absolute_tolerance;
    pt.pressure_relative_tolerance = 1.0e-12;
    fl::CpaVleEvaluator evaluator(model, pt);
    const auto feed = cpa_physical_test::feed(point, swapped);
    const auto starts = cpa_physical_test::starts(point, swapped);
    return fl::solve_cpa_pt_vle(
        point.pressure_pa, cpa_physical_test::temperature_k,
        feed, evaluator, split_options(), starts, starts);
}

} // namespace

int main() {
    try {
        const auto& low_pressure = cpa_physical_test::points().front();
        constexpr std::array<double, 7> tolerances{
            1.0e-4, 3.0e-5, 1.0e-5, 3.0e-6,
            1.0e-6, 3.0e-7, 1.0e-7};
        for (const double tolerance : tolerances) {
            const auto result = solve(low_pressure, tolerance, false);
            print_result("CPA_ROOT_SWEEP", tolerance, result);
        }

        const auto& permutation_point = cpa_physical_test::points()[3];
        for (const double tolerance : std::array<double, 3>{1.0e-4, 1.0e-5, 1.0e-6}) {
            const auto normal = solve(permutation_point, tolerance, false);
            const auto swapped = solve(permutation_point, tolerance, true);
            print_result("CPA_PERM_NORMAL", tolerance, normal);
            print_result("CPA_PERM_SWAPPED", tolerance, swapped);
        }
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "CPA_NUMERIC_AUDIT_ERROR " << error.what() << '\n';
        return 1;
    }
}
