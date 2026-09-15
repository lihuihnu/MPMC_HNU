#include <mpmc/flash/pr76_three_phase.hpp>

#include "soria_2025_ternary_blind_fixture.hpp"

#include <cmath>
#include <iomanip>
#include <iostream>
#include <optional>
#include <stdexcept>
#include <utility>

namespace {
namespace fl = mpmc::flash;
namespace fx = pr76_soria_2025_ternary_blind_test;

enum class StabilityClass { unstable, stable };

struct TransitionBracket {
    double unstable_pa{};
    double stable_pa{};
};

void require(bool condition, const char* message) {
    if (!condition) { throw std::runtime_error(message); }
}

void require_fixture_contract() {
    const auto model = fx::model();
    const auto& parameters = model.parameters();
    require(parameters.dataset_id() == "Soria-UFC-2025-M40-strict-PR76-blind" &&
                parameters.revision() == "table6-table9/pr101-binary-kij/frozen-v1",
            "blind-validation dataset identity drifted");
    require(std::abs(parameters.kij(0U, 1U) - 0.05226578047) <= 1.0e-14 &&
                std::abs(parameters.kij(0U, 2U) - 0.07884923875) <= 1.0e-14 &&
                parameters.kij(1U, 2U) == 0.0,
            "blind-validation kij values drifted from PR #101");
    require(std::abs(fx::experimental_temperature_k - 323.15) <= 1.0e-12 &&
                std::abs(fx::experimental_transition_pressure_pa - 11.21e6) <= 1.0e-6 &&
                std::abs(fx::expanded_pressure_uncertainty_pa - 0.07e6) <= 1.0e-6,
            "M-40 experimental target drifted from UFC Table 9");
}

fl::Pr76StabilityResult stability_at(double pressure_pa) {
    const auto model = fx::model();
    fl::Pr76StabilityEvaluator evaluator(model);
    return fl::test_pr76_pt_stability(
        pressure_pa, fx::experimental_temperature_k, fx::m40_feed, evaluator);
}

StabilityClass classify_stability(double pressure_pa) {
    const auto result = stability_at(pressure_pa);
    if (result.search.status == fl::StabilityStatus::unstable) {
        return StabilityClass::unstable;
    }
    if (result.search.status == fl::StabilityStatus::no_instability_found) {
        return StabilityClass::stable;
    }
    throw std::runtime_error(
        "finite TPD search became indeterminate while locating the M-40 blind transition");
}

TransitionBracket locate_transition_bracket() {
    constexpr int scan_lower_mpa = 2;
    constexpr int scan_upper_mpa = 20;
    std::optional<double> previous_pressure;
    std::optional<StabilityClass> previous_class;
    std::optional<TransitionBracket> found;
    int crossings = 0;

    for (int pressure_mpa = scan_lower_mpa;
         pressure_mpa <= scan_upper_mpa; ++pressure_mpa) {
        const double pressure_pa = static_cast<double>(pressure_mpa) * 1.0e6;
        const auto current = classify_stability(pressure_pa);
        if (previous_pressure && previous_class &&
            *previous_class == StabilityClass::unstable &&
            current == StabilityClass::stable) {
            ++crossings;
            found = TransitionBracket{*previous_pressure, pressure_pa};
        }
        previous_pressure = pressure_pa;
        previous_class = current;
    }

    require(crossings == 1 && found.has_value(),
            "strict-PR76 M-40 TPD scan did not find exactly one unstable-to-stable transition in 2-20 MPa");
    return *found;
}

TransitionBracket refine_transition(TransitionBracket bracket) {
    for (int iteration = 0; iteration < 18; ++iteration) {
        const double middle = 0.5 * (bracket.unstable_pa + bracket.stable_pa);
        const auto current = classify_stability(middle);
        if (current == StabilityClass::unstable) {
            bracket.unstable_pa = middle;
        } else {
            bracket.stable_pa = middle;
        }
    }
    return bracket;
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
            "pressure below predicted TPD boundary did not resolve as two phase");
    const auto* pair_result = stable_two_phase_result(result);
    require(pair_result != nullptr && pair_result->solution.candidate() != nullptr &&
                pair_result->solution.status ==
                    fl::PtSplitStatus::two_phase_no_instability_found &&
                pair_result->solution.final_stability.has_value() &&
                pair_result->solution.final_stability->status ==
                    fl::StabilityStatus::no_instability_found,
            "pressure below predicted TPD boundary lacks stable two-phase closure");
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
            "pressure above predicted TPD boundary did not resolve as single phase");
    require(result.base.solution.initial_stability.status ==
                fl::StabilityStatus::no_instability_found,
            "blind-validation upper point lacks finite single-phase stability closure");
    require(result.attempts.empty(),
            "blind-validation upper single-phase point unexpectedly entered max3 attempts");
    require(!result.global_stability_proven,
            "finite blind-validation single-phase review was promoted to global proof");
}

double lowest_sampled_tpd(const fl::Pr76StabilityResult& result) {
    require(result.search.lowest_sampled.has_value(),
            "blind-validation stability search did not retain a sampled TPD point");
    return result.search.lowest_sampled->value;
}

void blind_m40_transition_validation() {
    require_fixture_contract();

    const auto refined = refine_transition(locate_transition_bracket());
    const double predicted_pressure_pa =
        0.5 * (refined.unstable_pa + refined.stable_pa);
    const double numerical_width_pa = refined.stable_pa - refined.unstable_pa;
    const double absolute_error_pa =
        std::abs(predicted_pressure_pa - fx::experimental_transition_pressure_pa);
    const double relative_error = absolute_error_pa /
                                  fx::experimental_transition_pressure_pa;

    constexpr double topology_offset_pa = 0.10e6;
    const double lower_pressure_pa = predicted_pressure_pa - topology_offset_pa;
    const double upper_pressure_pa = predicted_pressure_pa + topology_offset_pa;
    const auto lower_max3 = solve_fresh_max3(lower_pressure_pa);
    const auto upper_max3 = solve_fresh_max3(upper_pressure_pa);
    require_two_phase_closure(lower_max3);
    require_single_phase_closure(upper_max3);

    const auto lower_stability = stability_at(lower_pressure_pa);
    const auto upper_stability = stability_at(upper_pressure_pa);
    require(lower_stability.search.status == fl::StabilityStatus::unstable &&
                upper_stability.search.status ==
                    fl::StabilityStatus::no_instability_found,
            "TPD topology around the predicted transition is not unstable-to-stable");

    std::cout << std::setprecision(12)
              << "M40 strict-PR76 blind Pcalc_MPa=" << predicted_pressure_pa / 1.0e6
              << " Pexp_MPa=" << fx::experimental_transition_pressure_pa / 1.0e6
              << " Uexp_MPa=" << fx::expanded_pressure_uncertainty_pa / 1.0e6
              << " abs_error_MPa=" << absolute_error_pa / 1.0e6
              << " relative_error=" << relative_error
              << " numerical_bracket_MPa=[" << refined.unstable_pa / 1.0e6
              << ',' << refined.stable_pa / 1.0e6 << ']'
              << " numerical_width_MPa=" << numerical_width_pa / 1.0e6 << '\n'
              << "  topology_lower_MPa=" << lower_pressure_pa / 1.0e6
              << " max3=two tpd_min=" << lowest_sampled_tpd(lower_stability)
              << " topology_upper_MPa=" << upper_pressure_pa / 1.0e6
              << " max3=single tpd_min=" << lowest_sampled_tpd(upper_stability)
              << '\n';

    require(numerical_width_pa <= 10.0,
            "blind-validation numerical transition bracket is wider than 10 Pa");
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
