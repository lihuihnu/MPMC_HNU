#include <mpmc/flash/pr76_three_phase.hpp>

#include "soria_2025_ternary_fixture.hpp"

#include <cmath>
#include <iostream>
#include <numeric>
#include <stdexcept>

namespace {
namespace fl = mpmc::flash;
namespace fx = pr76_soria_2025_ternary_test;

void require(bool condition, const char* message) {
    if (!condition) { throw std::runtime_error(message); }
}

void require_fixture_contract() {
    require(std::abs(fx::m40_feed[0] - 0.9020) <= 1.0e-15 &&
                std::abs(fx::m40_feed[1] - 0.0588) <= 1.0e-15 &&
                std::abs(fx::m40_feed[2] - 0.0392) <= 1.0e-15,
            "Soria M-40 nominal feed drifted from Tables 4/9");
    const double sum = std::accumulate(fx::m40_feed.begin(), fx::m40_feed.end(), 0.0);
    require(std::abs(sum - 1.0) <= 1.0e-15,
            "Soria M-40 nominal feed is not normalized");
    require(std::abs(fx::experimental_temperature_k - 323.15) <= 1.0e-12 &&
                std::abs(fx::experimental_transition_pressure_pa - 11.21e6) <= 1.0e-6 &&
                std::abs(fx::expanded_pressure_uncertainty_pa - 0.07e6) <= 1.0e-6,
            "Soria experimental transition datum drifted from Table 9");
    require(std::abs(fx::co2_decane_kij - 0.09670) <= 1.0e-15 &&
                std::abs(fx::co2_hexadecane_kij - 0.10127) <= 1.0e-15 &&
                fx::decane_hexadecane_kij == 0.0,
            "Soria 323.15-K kij values drifted from Table 12");
    require(std::abs(fx::source_pr_aard_fraction - 0.1539) <= 1.0e-15,
            "Soria M-40/323.15-K source AARD drifted from Table 13");

    const double lower = fx::lower_comparison_pressure_pa();
    const double upper = fx::upper_comparison_pressure_pa();
    require(std::abs(lower - 9.414781e6) <= 1.0e-6 &&
                std::abs(upper - 13.005219e6) <= 1.0e-6,
            "Soria source-derived comparison bracket changed");
    require(lower < fx::experimental_transition_pressure_pa &&
                fx::experimental_transition_pressure_pa < upper,
            "Soria comparison bracket does not contain experimental transition pressure");

    const auto model = fx::model();
    const auto& parameters = model.parameters();
    require(parameters.dataset_id() == "Soria-UFC-2025-CO2-C10-C16-PR76-transfer" &&
                parameters.revision() == "tables-6-9-12-13-m40-323K-v1",
            "Soria fixture dataset identity drifted");
    require(parameters.components().size() == 3U,
            "Soria fixture must contain exactly three components");
    require(parameters.applicability().temperature_k.has_value() &&
                parameters.applicability().temperature_k->lower == 323.15 &&
                parameters.applicability().temperature_k->upper == 323.15,
            "Soria fixture lost its 323.15-K-only applicability declaration");
    require(std::abs(parameters.kij(0U, 1U) - fx::co2_decane_kij) <= 1.0e-15 &&
                std::abs(parameters.kij(0U, 2U) - fx::co2_hexadecane_kij) <= 1.0e-15 &&
                parameters.kij(1U, 2U) == 0.0,
            "Soria fixture published the wrong kij matrix");
}

fl::Pr76PtMax3Result solve_fresh(double audit_pressure_pa) {
    const auto model = fx::model();
    fl::Pr76VleEvaluator evaluator(model);
    // No experimental phase compositions, K-values, or continuation hints are supplied.
    return fl::solve_pr76_pt_max3(
        audit_pressure_pa, fx::experimental_temperature_k,
        fx::m40_feed, evaluator);
}

void print_result(double audit_pressure_pa, const fl::Pr76PtMax3Result& result) {
    std::cerr << "P_MPa=" << audit_pressure_pa / 1.0e6
              << " max3_status=" << static_cast<int>(result.status)
              << " base_status=" << static_cast<int>(result.base.solution.status)
              << " initial_stability="
              << static_cast<int>(result.base.solution.initial_stability.status)
              << " attempts=" << result.attempts.size()
              << " diagnostic=" << result.diagnostic << '\n';
    if (const auto* pair = result.base.solution.candidate()) {
        std::cerr << "  base_beta_v=" << pair->fractions.vapor_fraction
                  << " fugacity_norm=" << pair->fugacity_norm
                  << " mass_abs=" << pair->fractions.mass_absolute
                  << " mass_rel=" << pair->fractions.mass_relative << '\n';
    }
}

const fl::Pr76PtSplitResult* stable_two_phase_result(
    const fl::Pr76PtMax3Result& result) {
    if (result.base.solution.status == fl::PtSplitStatus::two_phase_no_instability_found) {
        return &result.base;
    }
    return result.two_phase_neighbor();
}

void require_low_pressure_two_phase(const fl::Pr76PtMax3Result& result) {
    require(result.base.dataset_id == "Soria-UFC-2025-CO2-C10-C16-PR76-transfer",
            "low-pressure solve lost Soria dataset identity");
    require(result.status == fl::Pr76PtMax3Status::two_phase,
            "lower comparison pressure did not resolve as two phase");
    const auto* pair_result = stable_two_phase_result(result);
    require(pair_result != nullptr &&
                pair_result->solution.status ==
                    fl::PtSplitStatus::two_phase_no_instability_found &&
                pair_result->solution.candidate() != nullptr &&
                pair_result->solution.final_stability.has_value() &&
                pair_result->solution.final_stability->status ==
                    fl::StabilityStatus::no_instability_found,
            "lower comparison pressure lacks stable two-phase closure");
    const auto& state = *pair_result->solution.candidate();
    require(state.fugacity_norm <= 1.0e-11,
            "lower comparison pressure fugacity residual exceeds production gate");
    require(state.fractions.mass_absolute <= 1.0e-12 &&
                state.fractions.mass_relative <= 1.0e-10,
            "lower comparison pressure material balance exceeds production gate");
    require(state.fractions.vapor_fraction > 1.0e-10 &&
                state.fractions.vapor_fraction < 1.0 - 1.0e-10,
            "lower comparison pressure lies on a phase-disappearance gate");
    require(!result.global_stability_proven,
            "finite low-pressure review was promoted to global stability proof");
}

void require_high_pressure_single_phase(const fl::Pr76PtMax3Result& result) {
    require(result.base.dataset_id == "Soria-UFC-2025-CO2-C10-C16-PR76-transfer",
            "high-pressure solve lost Soria dataset identity");
    require(result.status == fl::Pr76PtMax3Status::single_phase &&
                result.base.solution.status ==
                    fl::PtSplitStatus::single_phase_no_instability_found,
            "upper comparison pressure did not resolve as single phase");
    require(result.base.solution.initial_stability.status ==
                fl::StabilityStatus::no_instability_found,
            "upper comparison pressure lacks finite single-phase stability closure");
    require(result.attempts.empty(),
            "single-phase upper comparison point unexpectedly entered three-phase attempts");
    require(!result.global_stability_proven,
            "finite high-pressure review was promoted to global stability proof");
}

void experimental_transition_pressure_gate() {
    require_fixture_contract();

    // The experiment reports L -> L+V on depressurization at 11.21 +/- 0.07 MPa.
    // The two numerical points are fixed before observing repository output. Their
    // width uses the source's 15.39% M-40/323.15-K aggregate PR AARD plus U(p).
    // AARD is only a coarse comparison scale: it is neither a pointwise bound nor
    // evidence that the source's high-omega PR alpha branch equals strict PR76.
    const double lower_pressure = fx::lower_comparison_pressure_pa();
    const double upper_pressure = fx::upper_comparison_pressure_pa();

    const auto lower = solve_fresh(lower_pressure);
    if (lower.status != fl::Pr76PtMax3Status::two_phase) {
        print_result(lower_pressure, lower);
    }
    require_low_pressure_two_phase(lower);

    const auto upper = solve_fresh(upper_pressure);
    if (upper.status != fl::Pr76PtMax3Status::single_phase) {
        print_result(upper_pressure, upper);
    }
    require_high_pressure_single_phase(upper);
}

} // namespace

int main() {
    try {
        experimental_transition_pressure_gate();
        std::cout << "[PASS] soria_2025_m40_experimental_transition_pressure\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "[FAIL] " << error.what() << '\n';
        return 1;
    }
}
