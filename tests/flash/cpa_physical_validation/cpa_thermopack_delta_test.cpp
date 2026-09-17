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
namespace parity = cpa_thermopack_snapshot;

void require(bool value, const char* message) {
    if (!value) { throw std::runtime_error(message); }
}

fl::PtSplitOptions literature_initialized_options() {
    fl::PtSplitOptions options;
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

} // namespace

int main() {
    try {
        const auto literature_parameters = cpa_physical_test::parameters(false);
        const auto parity_parameters = parity::parameters(false);
        require(parity_parameters.dataset_id() == parity::parity_dataset_id &&
                    parity_parameters.revision() == parity::parity_revision,
                "ThermoPack flash parity test lost named parameter-snapshot identity");

        const auto literature_model =
            th::CpaPtPhase::from_parameters(literature_parameters);
        const auto parity_model = th::CpaPtPhase::from_parameters(parity_parameters);
        fl::CpaVleEvaluator literature_evaluator(literature_model);
        fl::CpaVleEvaluator parity_evaluator(parity_model);

        double sum_abs_dbeta_v = 0.0;
        double sum_abs_dx = 0.0;
        double sum_abs_dy = 0.0;
        double max_abs_dbeta_v = 0.0;
        double max_abs_dx = 0.0;
        double max_abs_dy = 0.0;
        double parity_sum_abs_dbeta_v = 0.0;
        double parity_sum_abs_dx = 0.0;
        double parity_sum_abs_dy = 0.0;
        double parity_max_abs_dbeta_v = 0.0;
        double parity_max_abs_dx = 0.0;
        double parity_max_abs_dy = 0.0;
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
                    "frozen ThermoPack pressure no longer matches the literature fixture");

            const auto feed = cpa_physical_test::composition(oracle.feed_methanol);
            const auto fixture_feed = cpa_physical_test::feed(experimental, false);
            require(std::abs(feed.at(0) - fixture_feed.at(0)) <= 1.0e-14 &&
                        std::abs(feed.at(1) - fixture_feed.at(1)) <= 1.0e-14,
                    "frozen ThermoPack feed no longer matches the pre-existing physical fixture");

            // ThermoPack outputs remain comparison-only. Both MPMC_HNU runs use
            // the pre-existing Kurihara experimental x/y starts; the parity run
            // differs from the literature run only by selecting the separately
            // named, provenance-locked parity parameter snapshot.
            const auto starts = cpa_physical_test::starts(experimental, false);
            const std::vector<std::vector<double>> final_starts{feed};
            const auto options = literature_initialized_options();
            const auto literature_result = fl::solve_cpa_pt_vle(
                oracle.pressure_pa, oracle.temperature_k, feed,
                literature_evaluator, options, starts, final_starts);
            const auto parity_result = fl::solve_cpa_pt_vle(
                oracle.pressure_pa, oracle.temperature_k, feed,
                parity_evaluator, options, starts, final_starts);

            const auto& point = require_closed_two_phase(
                literature_result,
                "MPMC_HNU literature snapshot did not close the frozen ThermoPack state as two phase");
            const auto& parity_point = require_closed_two_phase(
                parity_result,
                "named ThermoPack parity snapshot did not close the frozen state as two phase");

            const double beta_v_mpmc = point.fractions.vapor_fraction;
            const double beta_l_mpmc = 1.0 - beta_v_mpmc;
            const double x_mpmc = point.fractions.liquid.at(0);
            const double y_mpmc = point.fractions.vapor.at(0);
            const double dbeta_v = beta_v_mpmc - oracle.beta_vapor;
            const double dbeta_l = beta_l_mpmc - oracle.beta_liquid;
            const double dx = x_mpmc - oracle.liquid_methanol;
            const double dy = y_mpmc - oracle.vapor_methanol;

            const double parity_beta_v = parity_point.fractions.vapor_fraction;
            const double parity_beta_l = 1.0 - parity_beta_v;
            const double parity_x = parity_point.fractions.liquid.at(0);
            const double parity_y = parity_point.fractions.vapor.at(0);
            const double parity_dbeta_v = parity_beta_v - oracle.beta_vapor;
            const double parity_dbeta_l = parity_beta_l - oracle.beta_liquid;
            const double parity_dx = parity_x - oracle.liquid_methanol;
            const double parity_dy = parity_y - oracle.vapor_methanol;

            require(std::isfinite(dbeta_v) && std::isfinite(dbeta_l) &&
                        std::isfinite(dx) && std::isfinite(dy) &&
                        std::isfinite(parity_dbeta_v) &&
                        std::isfinite(parity_dbeta_l) &&
                        std::isfinite(parity_dx) && std::isfinite(parity_dy),
                    "MPMC_HNU/ThermoPack delta became nonfinite");

            sum_abs_dbeta_v += std::abs(dbeta_v);
            sum_abs_dx += std::abs(dx);
            sum_abs_dy += std::abs(dy);
            max_abs_dbeta_v = std::max(max_abs_dbeta_v, std::abs(dbeta_v));
            max_abs_dx = std::max(max_abs_dx, std::abs(dx));
            max_abs_dy = std::max(max_abs_dy, std::abs(dy));
            parity_sum_abs_dbeta_v += std::abs(parity_dbeta_v);
            parity_sum_abs_dx += std::abs(parity_dx);
            parity_sum_abs_dy += std::abs(parity_dy);
            parity_max_abs_dbeta_v = std::max(
                parity_max_abs_dbeta_v, std::abs(parity_dbeta_v));
            parity_max_abs_dx = std::max(parity_max_abs_dx, std::abs(parity_dx));
            parity_max_abs_dy = std::max(parity_max_abs_dy, std::abs(parity_dy));
            ++accepted;

            std::cout
                << "CPA_THERMOPACK_LITERATURE_SNAPSHOT_DELTA"
                << " PPa=" << oracle.pressure_pa
                << " zMeOH=" << oracle.feed_methanol
                << " betaV_mpmc=" << beta_v_mpmc
                << " betaV_thermopack=" << oracle.beta_vapor
                << " d_betaV=" << dbeta_v
                << " betaL_mpmc=" << beta_l_mpmc
                << " betaL_thermopack=" << oracle.beta_liquid
                << " d_betaL=" << dbeta_l
                << " xMeOH_mpmc=" << x_mpmc
                << " xMeOH_thermopack=" << oracle.liquid_methanol
                << " d_xMeOH=" << dx
                << " yMeOH_mpmc=" << y_mpmc
                << " yMeOH_thermopack=" << oracle.vapor_methanol
                << " d_yMeOH=" << dy
                << " fug=" << point.fugacity_norm
                << " mass_abs=" << point.fractions.mass_absolute
                << " mass_rel=" << point.fractions.mass_relative;
            if (literature_result.solution.final_stability->lowest_sampled) {
                std::cout << " tpd_min="
                          << literature_result.solution.final_stability->lowest_sampled->value;
            }
            std::cout << '\n';

            std::cout
                << "CPA_THERMOPACK_PARITY_FLASH"
                << " PPa=" << oracle.pressure_pa
                << " dataset=" << parity_parameters.dataset_id()
                << " betaV_mpmc=" << parity_beta_v
                << " betaV_thermopack=" << oracle.beta_vapor
                << " d_betaV=" << parity_dbeta_v
                << " betaL_mpmc=" << parity_beta_l
                << " betaL_thermopack=" << oracle.beta_liquid
                << " d_betaL=" << parity_dbeta_l
                << " xMeOH_mpmc=" << parity_x
                << " xMeOH_thermopack=" << oracle.liquid_methanol
                << " d_xMeOH=" << parity_dx
                << " yMeOH_mpmc=" << parity_y
                << " yMeOH_thermopack=" << oracle.vapor_methanol
                << " d_yMeOH=" << parity_dy
                << " fug=" << parity_point.fugacity_norm
                << " mass_abs=" << parity_point.fractions.mass_absolute
                << " mass_rel=" << parity_point.fractions.mass_relative;
            if (parity_result.solution.final_stability->lowest_sampled) {
                std::cout << " tpd_min="
                          << parity_result.solution.final_stability->lowest_sampled->value;
            }
            std::cout << '\n';
        }

        require(accepted == cpa_thermopack_oracle::states.size(),
                "not every frozen ThermoPack state completed the read-only delta regression");
        const double count = static_cast<double>(accepted);
        std::cout
            << "CPA_THERMOPACK_LITERATURE_SNAPSHOT_DELTA_SUMMARY"
            << " states=" << accepted
            << " initialization=preexisting_Kurihara_xy_not_ThermoPack"
            << " mean_abs_d_betaV=" << sum_abs_dbeta_v / count
            << " max_abs_d_betaV=" << max_abs_dbeta_v
            << " mean_abs_d_xMeOH=" << sum_abs_dx / count
            << " max_abs_d_xMeOH=" << max_abs_dx
            << " mean_abs_d_yMeOH=" << sum_abs_dy / count
            << " max_abs_d_yMeOH=" << max_abs_dy
            << " magnitude_gate=none_read_only"
            << '\n';
        std::cout
            << "CPA_THERMOPACK_PARITY_FLASH_SUMMARY"
            << " states=" << accepted
            << " dataset=" << parity_parameters.dataset_id()
            << " revision=" << parity_parameters.revision()
            << " initialization=preexisting_Kurihara_xy_not_ThermoPack"
            << " mean_abs_d_betaV=" << parity_sum_abs_dbeta_v / count
            << " max_abs_d_betaV=" << parity_max_abs_dbeta_v
            << " mean_abs_d_xMeOH=" << parity_sum_abs_dx / count
            << " max_abs_d_xMeOH=" << parity_max_abs_dx
            << " mean_abs_d_yMeOH=" << parity_sum_abs_dy / count
            << " max_abs_d_yMeOH=" << parity_max_abs_dy
            << " magnitude_gate=none_parity_threshold_not_frozen"
            << '\n';
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "[FAIL] " << error.what() << '\n';
        return 1;
    }
}
