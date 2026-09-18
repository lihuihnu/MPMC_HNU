#include <mpmc/flash/cpa_split.hpp>

#include "cpa_thermopack_oracle_generated.hpp"
#include "cpa_thermopack_parameter_snapshot.hpp"
#include "test_support.hpp"

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <vector>

namespace {
namespace fl = mpmc::flash;
namespace th = mpmc::thermodynamics;
namespace snapshot = cpa_thermopack_snapshot;

void require(bool value, const char* message) {
    if (!value) { throw std::runtime_error(message); }
}

fl::PtSplitOptions comparison_options() {
    fl::PtSplitOptions options;
    // Keep the same deterministic, pre-existing physical-validation setup.
    // ThermoPack values are never used as MPMC_HNU iteration starts.
    options.initial_stability.automatic_starts = false;
    options.final_stability.automatic_starts = false;
    return options;
}

const fl::PtSplitState& require_closed_two_phase(
    const fl::CpaPtSplitResult& result, const char* message) {
    require(result.solution.status ==
                fl::PtSplitStatus::two_phase_no_instability_found &&
            result.solution.candidate() != nullptr &&
            result.solution.final_stability.has_value() &&
            result.solution.final_stability->status ==
                fl::StabilityStatus::no_instability_found,
            message);

    const auto& point = *result.solution.candidate();
    require(point.fugacity_norm <=
                result.solution.options.iteration.fugacity_tolerance,
            "CPA ThermoPack comparison state lost fugacity equality");
    require(point.fractions.mass_absolute <=
                result.solution.options.iteration.mass_absolute_tolerance &&
            point.fractions.mass_relative <=
                result.solution.options.iteration.mass_relative_tolerance,
            "CPA ThermoPack comparison state lost material balance");
    return point;
}

struct DeltaSummary {
    double sum_abs_beta_v{};
    double sum_abs_x{};
    double sum_abs_y{};
    double max_abs_beta_v{};
    double max_abs_x{};
    double max_abs_y{};
};

void accumulate(
    DeltaSummary& summary, double d_beta_v, double d_x, double d_y) {
    summary.sum_abs_beta_v += std::abs(d_beta_v);
    summary.sum_abs_x += std::abs(d_x);
    summary.sum_abs_y += std::abs(d_y);
    summary.max_abs_beta_v = std::max(summary.max_abs_beta_v, std::abs(d_beta_v));
    summary.max_abs_x = std::max(summary.max_abs_x, std::abs(d_x));
    summary.max_abs_y = std::max(summary.max_abs_y, std::abs(d_y));
}

} // namespace

int main() {
    try {
        const auto baseline_parameters = cpa_physical_test::parameters(false);
        const auto tc_only_parameters = snapshot::parameters(false);
        const auto baseline_model = th::CpaPtPhase::from_parameters(baseline_parameters);
        const auto tc_only_model = th::CpaPtPhase::from_parameters(tc_only_parameters);
        fl::CpaVleEvaluator baseline_evaluator(baseline_model);
        fl::CpaVleEvaluator tc_only_evaluator(tc_only_model);

        DeltaSummary baseline_summary;
        DeltaSummary tc_only_summary;
        std::size_t accepted = 0U;

        std::cout << std::setprecision(17);
        for (std::size_t index = 0U;
             index < cpa_thermopack_oracle::states.size(); ++index) {
            const auto& oracle = cpa_thermopack_oracle::states[index];
            const auto& experimental = cpa_physical_test::points()[index];

            require(oracle.phase_code == 0,
                    "frozen ThermoPack state no longer records a two-phase result");
            require(std::abs(oracle.temperature_k -
                             cpa_physical_test::temperature_k) <= 1.0e-12,
                    "frozen ThermoPack state temperature no longer matches the CPA fixture");
            require(std::abs(oracle.pressure_pa - experimental.pressure_pa) <= 1.0e-9,
                    "frozen ThermoPack pressure no longer matches the physical fixture");

            const auto feed = cpa_physical_test::composition(oracle.feed_methanol);
            const auto fixture_feed = cpa_physical_test::feed(experimental, false);
            require(std::abs(feed.at(0) - fixture_feed.at(0)) <= 1.0e-14 &&
                        std::abs(feed.at(1) - fixture_feed.at(1)) <= 1.0e-14,
                    "frozen ThermoPack feed no longer matches the pre-existing physical fixture");

            // Both A and B use only the repository's pre-existing Kurihara
            // experimental starts. ThermoPack x/y/beta remain comparison-only.
            // A is the current literature snapshot. B changes only Classic-alpha
            // critical_temperature_k to the pinned ThermoPack component values;
            // the separate parameter-identity regression proves all other CPA
            // scientific fields and their provenance remain unchanged.
            const auto starts = cpa_physical_test::starts(experimental, false);
            const std::vector<std::vector<double>> final_starts{feed};
            const auto options = comparison_options();
            const auto baseline_result = fl::solve_cpa_pt_vle(
                oracle.pressure_pa, oracle.temperature_k, feed,
                baseline_evaluator, options, starts, final_starts);
            const auto tc_only_result = fl::solve_cpa_pt_vle(
                oracle.pressure_pa, oracle.temperature_k, feed,
                tc_only_evaluator, options, starts, final_starts);

            const auto& baseline = require_closed_two_phase(
                baseline_result,
                "current literature CPA snapshot did not close a frozen ThermoPack state");
            const auto& tc_only = require_closed_two_phase(
                tc_only_result,
                "Tc-only ThermoPack-aligned snapshot did not close a frozen state");

            const double baseline_d_beta_v =
                baseline.fractions.vapor_fraction - oracle.beta_vapor;
            const double baseline_d_x =
                baseline.fractions.liquid.at(0) - oracle.liquid_methanol;
            const double baseline_d_y =
                baseline.fractions.vapor.at(0) - oracle.vapor_methanol;
            const double tc_only_d_beta_v =
                tc_only.fractions.vapor_fraction - oracle.beta_vapor;
            const double tc_only_d_x =
                tc_only.fractions.liquid.at(0) - oracle.liquid_methanol;
            const double tc_only_d_y =
                tc_only.fractions.vapor.at(0) - oracle.vapor_methanol;

            require(std::isfinite(baseline_d_beta_v) &&
                        std::isfinite(baseline_d_x) &&
                        std::isfinite(baseline_d_y) &&
                        std::isfinite(tc_only_d_beta_v) &&
                        std::isfinite(tc_only_d_x) &&
                        std::isfinite(tc_only_d_y),
                    "MPMC_HNU/ThermoPack attribution delta became nonfinite");

            accumulate(
                baseline_summary, baseline_d_beta_v, baseline_d_x, baseline_d_y);
            accumulate(
                tc_only_summary, tc_only_d_beta_v, tc_only_d_x, tc_only_d_y);
            ++accepted;

            std::cout
                << "CPA_THERMOPACK_TC_SINGLE_FACTOR"
                << " PPa=" << oracle.pressure_pa
                << " baseline_d_betaV=" << baseline_d_beta_v
                << " tc_only_d_betaV=" << tc_only_d_beta_v
                << " baseline_d_xMeOH=" << baseline_d_x
                << " tc_only_d_xMeOH=" << tc_only_d_x
                << " baseline_d_yMeOH=" << baseline_d_y
                << " tc_only_d_yMeOH=" << tc_only_d_y
                << " baseline_fug=" << baseline.fugacity_norm
                << " tc_only_fug=" << tc_only.fugacity_norm
                << '\n';
        }

        require(accepted == cpa_thermopack_oracle::states.size(),
                "not every frozen ThermoPack state completed Tc single-factor attribution");

        const double count = static_cast<double>(accepted);
        std::cout
            << "CPA_THERMOPACK_TC_SINGLE_FACTOR_SUMMARY"
            << " states=" << accepted
            << " factor=critical_temperature_k_only"
            << " Tc_MeOH_A=512.64"
            << " Tc_MeOH_B=512.6"
            << " Tc_H2O_A=647.29"
            << " Tc_H2O_B=647.3"
            << " A_mean_abs_d_betaV=" << baseline_summary.sum_abs_beta_v / count
            << " B_mean_abs_d_betaV=" << tc_only_summary.sum_abs_beta_v / count
            << " A_max_abs_d_betaV=" << baseline_summary.max_abs_beta_v
            << " B_max_abs_d_betaV=" << tc_only_summary.max_abs_beta_v
            << " A_mean_abs_d_xMeOH=" << baseline_summary.sum_abs_x / count
            << " B_mean_abs_d_xMeOH=" << tc_only_summary.sum_abs_x / count
            << " A_max_abs_d_xMeOH=" << baseline_summary.max_abs_x
            << " B_max_abs_d_xMeOH=" << tc_only_summary.max_abs_x
            << " A_mean_abs_d_yMeOH=" << baseline_summary.sum_abs_y / count
            << " B_mean_abs_d_yMeOH=" << tc_only_summary.sum_abs_y / count
            << " A_max_abs_d_yMeOH=" << baseline_summary.max_abs_y
            << " B_max_abs_d_yMeOH=" << tc_only_summary.max_abs_y
            << " magnitude_gate=none_attribution_only"
            << '\n';
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "[FAIL] " << error.what() << '\n';
        return 1;
    }
}
