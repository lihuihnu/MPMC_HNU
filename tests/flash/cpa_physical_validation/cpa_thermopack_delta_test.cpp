#include <mpmc/flash/cpa_split.hpp>

#include "cpa_thermopack_oracle_generated.hpp"
#include "cpa_thermopack_phase_kernel_generated.hpp"
#include "test_support.hpp"

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

namespace {
namespace fl = mpmc::flash;
namespace th = mpmc::thermodynamics;
namespace phase_ref = cpa_thermopack_phase_kernel;

void require(bool value, const char* message) {
    if (!value) { throw std::runtime_error(message); }
}

fl::PtSplitOptions literature_initialized_options() {
    fl::PtSplitOptions options;
    options.initial_stability.automatic_starts = false;
    options.final_stability.automatic_starts = false;
    return options;
}

th::Provenance thermopack_tc_source(std::string locator) {
    return {
        th::SourceKind::database,
        "thermotools/thermopack pinned source tree",
        "d68c794c7342bfc6938eb424a1fbb88b7780b738",
        std::move(locator),
        "Test-only critical-temperature substitution used to isolate the Classic-alpha input delta; production CPA parameters are unchanged.",
        "Read from the pinned ThermoPack component JSON used to initialize cbeos%single(i)%tc.",
        "ThermoPack source metadata used only for repository numerical audit."};
}

th::CpaParameterSet with_thermopack_alpha_tc(
    const th::CpaParameterSet& baseline) {
    th::CpaParameterInput input;
    input.dataset_id = baseline.dataset_id() + "__thermopack-alpha-tc-causal-audit";
    input.revision = baseline.revision() + "__test-only-Tc-substitution";
    input.applicability = baseline.applicability();
    input.pure.assign(baseline.pure_records().begin(), baseline.pure_records().end());
    input.binary.assign(baseline.binary_records().begin(), baseline.binary_records().end());
    input.association_pairs.assign(
        baseline.association_records().begin(), baseline.association_records().end());

    for (auto& pure : input.pure) {
        if (pure.component_id == "METHANOL") {
            pure.critical_temperature_k.value =
                phase_ref::methanol_alpha_critical_temperature_k;
            pure.critical_temperature_k.source = thermopack_tc_source(
                "fluids/Methanol.json critical.temperature = 512.6 K");
        } else if (pure.component_id == "WATER") {
            pure.critical_temperature_k.value =
                phase_ref::water_alpha_critical_temperature_k;
            pure.critical_temperature_k.source = thermopack_tc_source(
                "fluids/Water.json critical.temperature = 647.3 K");
        } else {
            throw std::runtime_error(
                "unexpected component in methanol-water ThermoPack Tc audit");
        }
    }

    std::vector<std::string> order;
    order.reserve(baseline.components().size());
    for (const auto& component : baseline.components().items()) {
        order.push_back(component.id);
    }
    return th::CpaParameterSet::create(
        baseline.components().items(), order, input);
}

const fl::PtSplitPoint& require_closed_two_phase(
    const fl::CpaPtVleResult& result, const char* message) {
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
        const auto parameters = cpa_physical_test::parameters(false);
        const auto aligned_parameters = with_thermopack_alpha_tc(parameters);
        const auto model = th::CpaPtPhase::from_parameters(parameters);
        const auto aligned_model = th::CpaPtPhase::from_parameters(aligned_parameters);
        fl::CpaVleEvaluator evaluator(model);
        fl::CpaVleEvaluator aligned_evaluator(aligned_model);

        double sum_abs_dbeta_v = 0.0;
        double sum_abs_dx = 0.0;
        double sum_abs_dy = 0.0;
        double max_abs_dbeta_v = 0.0;
        double max_abs_dx = 0.0;
        double max_abs_dy = 0.0;
        double aligned_sum_abs_dbeta_v = 0.0;
        double aligned_sum_abs_dx = 0.0;
        double aligned_sum_abs_dy = 0.0;
        double aligned_max_abs_dbeta_v = 0.0;
        double aligned_max_abs_dx = 0.0;
        double aligned_max_abs_dy = 0.0;
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

            // ThermoPack outputs remain comparison-only. MPMC_HNU is initialized
            // with the pre-existing Kurihara experimental x/y starts used by the
            // repository's physical-validation regression, never with ThermoPack
            // x/y/beta. The Tc-only counterfactual uses exactly the same starts.
            const auto starts = cpa_physical_test::starts(experimental, false);
            const std::vector<std::vector<double>> final_starts{feed};
            const auto options = literature_initialized_options();
            const auto result = fl::solve_cpa_pt_vle(
                oracle.pressure_pa, oracle.temperature_k, feed, evaluator,
                options, starts, final_starts);
            const auto aligned_result = fl::solve_cpa_pt_vle(
                oracle.pressure_pa, oracle.temperature_k, feed, aligned_evaluator,
                options, starts, final_starts);

            const auto& point = require_closed_two_phase(
                result,
                "MPMC_HNU did not close the frozen ThermoPack state as two phase");
            const auto& aligned_point = require_closed_two_phase(
                aligned_result,
                "Tc-only counterfactual did not close the frozen ThermoPack state as two phase");

            const double beta_v_mpmc = point.fractions.vapor_fraction;
            const double beta_l_mpmc = 1.0 - beta_v_mpmc;
            const double x_mpmc = point.fractions.liquid.at(0);
            const double y_mpmc = point.fractions.vapor.at(0);
            const double dbeta_v = beta_v_mpmc - oracle.beta_vapor;
            const double dbeta_l = beta_l_mpmc - oracle.beta_liquid;
            const double dx = x_mpmc - oracle.liquid_methanol;
            const double dy = y_mpmc - oracle.vapor_methanol;

            const double aligned_beta_v = aligned_point.fractions.vapor_fraction;
            const double aligned_beta_l = 1.0 - aligned_beta_v;
            const double aligned_x = aligned_point.fractions.liquid.at(0);
            const double aligned_y = aligned_point.fractions.vapor.at(0);
            const double aligned_dbeta_v = aligned_beta_v - oracle.beta_vapor;
            const double aligned_dbeta_l = aligned_beta_l - oracle.beta_liquid;
            const double aligned_dx = aligned_x - oracle.liquid_methanol;
            const double aligned_dy = aligned_y - oracle.vapor_methanol;

            require(std::isfinite(dbeta_v) && std::isfinite(dbeta_l) &&
                        std::isfinite(dx) && std::isfinite(dy) &&
                        std::isfinite(aligned_dbeta_v) &&
                        std::isfinite(aligned_dbeta_l) &&
                        std::isfinite(aligned_dx) && std::isfinite(aligned_dy),
                    "MPMC_HNU/ThermoPack delta became nonfinite");

            sum_abs_dbeta_v += std::abs(dbeta_v);
            sum_abs_dx += std::abs(dx);
            sum_abs_dy += std::abs(dy);
            max_abs_dbeta_v = std::max(max_abs_dbeta_v, std::abs(dbeta_v));
            max_abs_dx = std::max(max_abs_dx, std::abs(dx));
            max_abs_dy = std::max(max_abs_dy, std::abs(dy));
            aligned_sum_abs_dbeta_v += std::abs(aligned_dbeta_v);
            aligned_sum_abs_dx += std::abs(aligned_dx);
            aligned_sum_abs_dy += std::abs(aligned_dy);
            aligned_max_abs_dbeta_v = std::max(
                aligned_max_abs_dbeta_v, std::abs(aligned_dbeta_v));
            aligned_max_abs_dx = std::max(aligned_max_abs_dx, std::abs(aligned_dx));
            aligned_max_abs_dy = std::max(aligned_max_abs_dy, std::abs(aligned_dy));
            ++accepted;

            std::cout
                << "CPA_THERMOPACK_DELTA"
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
            if (result.solution.final_stability->lowest_sampled) {
                std::cout << " tpd_min="
                          << result.solution.final_stability->lowest_sampled->value;
            }
            std::cout << '\n';

            std::cout
                << "CPA_THERMOPACK_TC_ALIGNED_FLASH"
                << " PPa=" << oracle.pressure_pa
                << " changed_fields=critical_temperature_k_only"
                << " betaV_mpmc=" << aligned_beta_v
                << " betaV_thermopack=" << oracle.beta_vapor
                << " d_betaV=" << aligned_dbeta_v
                << " betaL_mpmc=" << aligned_beta_l
                << " betaL_thermopack=" << oracle.beta_liquid
                << " d_betaL=" << aligned_dbeta_l
                << " xMeOH_mpmc=" << aligned_x
                << " xMeOH_thermopack=" << oracle.liquid_methanol
                << " d_xMeOH=" << aligned_dx
                << " yMeOH_mpmc=" << aligned_y
                << " yMeOH_thermopack=" << oracle.vapor_methanol
                << " d_yMeOH=" << aligned_dy
                << " fug=" << aligned_point.fugacity_norm
                << " mass_abs=" << aligned_point.fractions.mass_absolute
                << " mass_rel=" << aligned_point.fractions.mass_relative;
            if (aligned_result.solution.final_stability->lowest_sampled) {
                std::cout << " tpd_min="
                          << aligned_result.solution.final_stability->lowest_sampled->value;
            }
            std::cout << '\n';
        }

        require(accepted == cpa_thermopack_oracle::states.size(),
                "not every frozen ThermoPack state completed the read-only delta regression");
        const double count = static_cast<double>(accepted);
        std::cout
            << "CPA_THERMOPACK_DELTA_SUMMARY"
            << " states=" << accepted
            << " initialization=preexisting_Kurihara_xy_not_ThermoPack"
            << " mean_abs_d_betaV=" << sum_abs_dbeta_v / count
            << " max_abs_d_betaV=" << max_abs_dbeta_v
            << " mean_abs_d_xMeOH=" << sum_abs_dx / count
            << " max_abs_d_xMeOH=" << max_abs_dx
            << " mean_abs_d_yMeOH=" << sum_abs_dy / count
            << " max_abs_d_yMeOH=" << max_abs_dy
            << " parity_magnitude_gate=none_read_only"
            << '\n';
        std::cout
            << "CPA_THERMOPACK_TC_ALIGNED_FLASH_SUMMARY"
            << " states=" << accepted
            << " changed_fields=critical_temperature_k_only"
            << " initialization=preexisting_Kurihara_xy_not_ThermoPack"
            << " mean_abs_d_betaV=" << aligned_sum_abs_dbeta_v / count
            << " max_abs_d_betaV=" << aligned_max_abs_dbeta_v
            << " mean_abs_d_xMeOH=" << aligned_sum_abs_dx / count
            << " max_abs_d_xMeOH=" << aligned_max_abs_dx
            << " mean_abs_d_yMeOH=" << aligned_sum_abs_dy / count
            << " max_abs_d_yMeOH=" << aligned_max_abs_dy
            << " magnitude_gate=none_causal_audit"
            << '\n';
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "[FAIL] " << error.what() << '\n';
        return 1;
    }
}
