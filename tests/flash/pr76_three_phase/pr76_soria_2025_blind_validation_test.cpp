#include <mpmc/flash/pr76_three_phase.hpp>

#include "soria_2025_ternary_blind_fixture.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <iomanip>
#include <iostream>
#include <numeric>
#include <optional>
#include <stdexcept>

namespace {
namespace fl = mpmc::flash;
namespace fx = pr76_soria_2025_ternary_blind_test;

struct BubbleResidual {
    double log_sum{};
    std::array<double, 3> vapor{};
    double fugacity_norm{};
    double liquid_z{};
    double vapor_z{};
};

struct BubbleBracket {
    double left_pa{};
    double right_pa{};
    double left_residual{};
    double right_residual{};
};

void require(bool condition, const char* message) {
    if (!condition) { throw std::runtime_error(message); }
}

void require_fixture_contract() {
    require(std::abs(fx::m40_feed[0] - 0.9020) <= 1.0e-15 &&
                std::abs(fx::m40_feed[1] - 0.0588) <= 1.0e-15 &&
                std::abs(fx::m40_feed[2] - 0.0392) <= 1.0e-15,
            "M-40 nominal feed drifted from the frozen blind-validation contract");
    const double sum = std::accumulate(fx::m40_feed.begin(), fx::m40_feed.end(), 0.0);
    require(std::abs(sum - 1.0) <= 1.0e-15,
            "M-40 nominal feed is not normalized");
    require(std::abs(fx::experimental_temperature_k - 323.15) <= 1.0e-12 &&
                std::abs(fx::experimental_transition_pressure_pa - 11.21e6) <= 1.0e-6 &&
                std::abs(fx::expanded_pressure_uncertainty_pa - 0.07e6) <= 1.0e-6,
            "M-40 experimental blind target drifted from UFC Table 9");
    require(std::abs(fx::co2_decane_kij - 0.05226578047) <= 1.0e-14 &&
                std::abs(fx::co2_hexadecane_kij - 0.07884923875) <= 1.0e-14 &&
                fx::decane_hexadecane_kij == 0.0,
            "strict-PR76 binary-calibrated kij values drifted from PR #101");

    const auto model = fx::model();
    const auto& parameters = model.parameters();
    require(parameters.dataset_id() == "Soria-UFC-2025-M40-strict-PR76-blind" &&
                parameters.revision() == "table6-table9/pr101-binary-kij/frozen-v1",
            "blind-validation dataset identity drifted");
    require(parameters.components().size() == 3U,
            "blind-validation fixture must contain exactly three components");
    require(std::abs(parameters.kij(0U, 1U) - fx::co2_decane_kij) <= 1.0e-14 &&
                std::abs(parameters.kij(0U, 2U) - fx::co2_hexadecane_kij) <= 1.0e-14 &&
                parameters.kij(1U, 2U) == 0.0,
            "blind-validation fixture published the wrong kij matrix");
}

double wilson_k(double pressure_pa, double temperature_k,
                double tc, double pc, double omega) {
    return (pc / pressure_pa) *
        std::exp(5.373 * (1.0 + omega) * (1.0 - tc / temperature_k));
}

std::optional<BubbleResidual> bubble_residual(
    double pressure_pa, fl::Pr76VleEvaluator& evaluator) {
    const auto& liquid_x = fx::m40_feed;
    fl::PtSplitPhase liquid;
    try {
        liquid = evaluator(pressure_pa, fx::experimental_temperature_k,
                           liquid_x, fl::PtPhaseRole::liquid_candidate);
    } catch (const fl::StabilityPropertyError&) {
        return std::nullopt;
    }

    std::array<double, 3> vapor{};
    double raw_sum = 0.0;
    for (std::size_t i = 0U; i < vapor.size(); ++i) {
        const double k = wilson_k(
            pressure_pa, fx::experimental_temperature_k,
            fx::critical_temperatures_k[i], fx::critical_pressures_pa[i],
            fx::acentric_factors[i]);
        vapor[i] = liquid_x[i] * k;
        raw_sum += vapor[i];
    }
    if (!(raw_sum > 0.0) || !std::isfinite(raw_sum)) { return std::nullopt; }
    for (double& value : vapor) { value /= raw_sum; }

    for (int iteration = 0; iteration < 120; ++iteration) {
        fl::PtSplitPhase vapor_phase;
        try {
            vapor_phase = evaluator(pressure_pa, fx::experimental_temperature_k,
                                    vapor, fl::PtPhaseRole::vapor_candidate);
        } catch (const fl::StabilityPropertyError&) {
            return std::nullopt;
        }

        std::array<double, 3> log_raw{};
        double largest = -std::numeric_limits<double>::infinity();
        for (std::size_t i = 0U; i < vapor.size(); ++i) {
            log_raw[i] = std::log(liquid_x[i]) + liquid.activity.ln_phi[i] -
                         vapor_phase.activity.ln_phi[i];
            largest = std::max(largest, log_raw[i]);
        }
        double scaled_sum = 0.0;
        for (const double value : log_raw) {
            scaled_sum += std::exp(value - largest);
        }
        if (!(scaled_sum > 0.0) || !std::isfinite(scaled_sum)) {
            return std::nullopt;
        }
        const double log_sum = largest + std::log(scaled_sum);
        std::array<double, 3> next{};
        for (std::size_t i = 0U; i < next.size(); ++i) {
            next[i] = std::exp(log_raw[i] - log_sum);
        }

        double change = 0.0;
        for (std::size_t i = 0U; i < next.size(); ++i) {
            change = std::max(change, std::abs(next[i] - vapor[i]));
        }
        if (change <= 2.0e-12) {
            try {
                vapor_phase = evaluator(pressure_pa, fx::experimental_temperature_k,
                                        next, fl::PtPhaseRole::vapor_candidate);
            } catch (const fl::StabilityPropertyError&) {
                return std::nullopt;
            }

            std::array<double, 3> final_log_raw{};
            double final_largest = -std::numeric_limits<double>::infinity();
            for (std::size_t i = 0U; i < next.size(); ++i) {
                final_log_raw[i] = std::log(liquid_x[i]) + liquid.activity.ln_phi[i] -
                                   vapor_phase.activity.ln_phi[i];
                final_largest = std::max(final_largest, final_log_raw[i]);
            }
            double final_scaled_sum = 0.0;
            for (const double value : final_log_raw) {
                final_scaled_sum += std::exp(value - final_largest);
            }
            if (!(final_scaled_sum > 0.0) || !std::isfinite(final_scaled_sum)) {
                return std::nullopt;
            }
            const double final_log_sum = final_largest + std::log(final_scaled_sum);
            double fugacity_norm = 0.0;
            for (std::size_t i = 0U; i < next.size(); ++i) {
                const double residual =
                    std::log(liquid_x[i]) + liquid.activity.ln_phi[i] -
                    std::log(next[i]) - vapor_phase.activity.ln_phi[i];
                fugacity_norm = std::max(fugacity_norm, std::abs(residual));
            }
            const double relative_z = std::abs(vapor_phase.z - liquid.z) /
                                      std::max(vapor_phase.z, liquid.z);
            if (relative_z <= 1.0e-8 || !std::isfinite(final_log_sum) ||
                !std::isfinite(fugacity_norm)) {
                return std::nullopt;
            }
            return BubbleResidual{
                final_log_sum, next, fugacity_norm, liquid.z, vapor_phase.z};
        }

        for (std::size_t i = 0U; i < vapor.size(); ++i) {
            vapor[i] = 0.5 * vapor[i] + 0.5 * next[i];
        }
        const double sum = std::accumulate(vapor.begin(), vapor.end(), 0.0);
        for (double& value : vapor) { value /= sum; }
    }
    return std::nullopt;
}

BubbleBracket locate_unique_bubble_bracket(fl::Pr76VleEvaluator& evaluator) {
    constexpr double scan_lower_pa = 2.0e6;
    constexpr double scan_upper_pa = 20.0e6;
    constexpr int scan_intervals = 72;

    std::optional<double> previous_pressure;
    std::optional<double> previous_residual;
    std::optional<BubbleBracket> found;
    int crossings = 0;

    for (int index = 0; index <= scan_intervals; ++index) {
        const double fraction = static_cast<double>(index) /
                                static_cast<double>(scan_intervals);
        const double pressure = scan_lower_pa +
            (scan_upper_pa - scan_lower_pa) * fraction;
        const auto state = bubble_residual(pressure, evaluator);
        if (!state) {
            previous_pressure.reset();
            previous_residual.reset();
            continue;
        }
        if (previous_pressure && previous_residual &&
            ((*previous_residual <= 0.0 && state->log_sum >= 0.0) ||
             (*previous_residual >= 0.0 && state->log_sum <= 0.0))) {
            ++crossings;
            found = BubbleBracket{*previous_pressure, pressure,
                                  *previous_residual, state->log_sum};
        }
        previous_pressure = pressure;
        previous_residual = state->log_sum;
    }

    require(crossings == 1 && found.has_value(),
            "blind strict-PR76 scan did not find exactly one admissible M-40 VLE bubble root in 2-20 MPa");
    return *found;
}

std::pair<double, BubbleResidual> solve_blind_bubble_pressure(
    fl::Pr76VleEvaluator& evaluator) {
    auto bracket = locate_unique_bubble_bracket(evaluator);
    double left = bracket.left_pa;
    double right = bracket.right_pa;
    double f_left = bracket.left_residual;
    double f_right = bracket.right_residual;

    for (int iteration = 0; iteration < 16; ++iteration) {
        const double middle = 0.5 * (left + right);
        const auto state = bubble_residual(middle, evaluator);
        require(state.has_value(),
                "blind strict-PR76 bubble root lost its admissible VLE branch during bisection");
        if ((f_left <= 0.0 && state->log_sum >= 0.0) ||
            (f_left >= 0.0 && state->log_sum <= 0.0)) {
            right = middle;
            f_right = state->log_sum;
        } else {
            left = middle;
            f_left = state->log_sum;
        }
    }
    (void)f_right;
    const double pressure = 0.5 * (left + right);
    const auto final_state = bubble_residual(pressure, evaluator);
    require(final_state.has_value(),
            "blind strict-PR76 bubble root is not evaluable at its final pressure");
    return {pressure, *final_state};
}

fl::Pr76PtMax3Result solve_fresh_max3(double pressure_pa) {
    const auto model = fx::model();
    fl::Pr76VleEvaluator evaluator(model);
    return fl::solve_pr76_pt_max3(
        pressure_pa, fx::experimental_temperature_k, fx::m40_feed, evaluator);
}

const fl::Pr76PtSplitResult* stable_two_phase_result(
    const fl::Pr76PtMax3Result& result) {
    if (result.base.solution.status == fl::PtSplitStatus::two_phase_no_instability_found) {
        return &result.base;
    }
    return result.two_phase_neighbor();
}

void require_two_phase_closure(const fl::Pr76PtMax3Result& result) {
    require(result.status == fl::Pr76PtMax3Status::two_phase,
            "pressure below predicted bubble point did not resolve as two phase");
    const auto* pair_result = stable_two_phase_result(result);
    require(pair_result != nullptr && pair_result->solution.candidate() != nullptr &&
                pair_result->solution.status ==
                    fl::PtSplitStatus::two_phase_no_instability_found &&
                pair_result->solution.final_stability.has_value() &&
                pair_result->solution.final_stability->status ==
                    fl::StabilityStatus::no_instability_found,
            "pressure below predicted bubble point lacks stable two-phase closure");
    const auto& state = *pair_result->solution.candidate();
    require(state.fugacity_norm <= 1.0e-11,
            "blind-validation two-phase fugacity residual exceeds production gate");
    require(state.fractions.mass_absolute <= 1.0e-12 &&
                state.fractions.mass_relative <= 1.0e-10,
            "blind-validation two-phase material balance exceeds production gate");
    require(state.fractions.vapor_fraction > 1.0e-10 &&
                state.fractions.vapor_fraction < 1.0 - 1.0e-10,
            "blind-validation lower point lies on the phase-disappearance gate");
    require(!result.global_stability_proven,
            "finite blind-validation two-phase review was promoted to global proof");
}

void require_single_phase_closure(const fl::Pr76PtMax3Result& result) {
    require(result.status == fl::Pr76PtMax3Status::single_phase &&
                result.base.solution.status ==
                    fl::PtSplitStatus::single_phase_no_instability_found,
            "pressure above predicted bubble point did not resolve as single phase");
    require(result.base.solution.initial_stability.status ==
                fl::StabilityStatus::no_instability_found,
            "blind-validation upper point lacks finite single-phase stability closure");
    require(result.attempts.empty(),
            "blind-validation single-phase point unexpectedly entered three-phase attempts");
    require(!result.global_stability_proven,
            "finite blind-validation single-phase review was promoted to global proof");
}

void blind_m40_transition_validation() {
    require_fixture_contract();

    const auto model = fx::model();
    fl::Pr76VleEvaluator evaluator(model);
    const auto [predicted_pressure_pa, bubble] = solve_blind_bubble_pressure(evaluator);

    const double absolute_error_pa =
        std::abs(predicted_pressure_pa - fx::experimental_transition_pressure_pa);
    const double relative_error = absolute_error_pa /
                                  fx::experimental_transition_pressure_pa;

    constexpr double topology_offset_pa = 0.10e6;
    const auto lower = solve_fresh_max3(predicted_pressure_pa - topology_offset_pa);
    const auto upper = solve_fresh_max3(predicted_pressure_pa + topology_offset_pa);
    require_two_phase_closure(lower);
    require_single_phase_closure(upper);

    std::cout << std::setprecision(12)
              << "M40 strict-PR76 blind Pcalc_MPa=" << predicted_pressure_pa / 1.0e6
              << " Pexp_MPa=" << fx::experimental_transition_pressure_pa / 1.0e6
              << " Uexp_MPa=" << fx::expanded_pressure_uncertainty_pa / 1.0e6
              << " abs_error_MPa=" << absolute_error_pa / 1.0e6
              << " relative_error=" << relative_error << '\n'
              << "  bubble_log_sum=" << bubble.log_sum
              << " fugacity_norm=" << bubble.fugacity_norm
              << " liquid_Z=" << bubble.liquid_z
              << " vapor_Z=" << bubble.vapor_z
              << " yCO2=" << bubble.vapor[0]
              << " yC10=" << bubble.vapor[1]
              << " yC16=" << bubble.vapor[2] << '\n'
              << "  topology_check_MPa="
              << (predicted_pressure_pa - topology_offset_pa) / 1.0e6
              << " -> two_phase, "
              << (predicted_pressure_pa + topology_offset_pa) / 1.0e6
              << " -> single_phase\n";

    require(std::abs(bubble.log_sum) <= 2.0e-8 &&
                bubble.fugacity_norm <= 2.0e-8,
            "blind strict-PR76 ternary bubble equations did not close");
    require(absolute_error_pa <= fx::expanded_pressure_uncertainty_pa,
            "strict-PR76 independently calibrated M-40 blind prediction lies outside the experimental expanded pressure uncertainty");
}

} // namespace

int main() {
    try {
        blind_m40_transition_validation();
        std::cout << "[PASS] soria_2025_m40_strict_pr76_blind_validation\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "[FAIL] " << error.what() << '\n';
        return 1;
    }
}
