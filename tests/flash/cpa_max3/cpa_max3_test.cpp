#include <mpmc/flash/cpa_max3_phase_set.hpp>

#include "test_support.hpp"

#include <array>
#include <cmath>
#include <iostream>
#include <stdexcept>
#include <string_view>
#include <vector>

namespace {
namespace fl = mpmc::flash;
namespace th = mpmc::thermodynamics;
using Vec = std::vector<double>;

void require(bool value, const char* message) {
    if (!value) { throw std::runtime_error(message); }
}

Vec log_ratio(const Vec& numerator, const Vec& denominator) {
    require(numerator.size() == denominator.size(), "log ratio dimension mismatch");
    Vec result(numerator.size(), 0.0);
    for (std::size_t i = 0; i < result.size(); ++i) {
        result[i] = std::log(numerator[i]) - std::log(denominator[i]);
    }
    return result;
}

void fixed_three_phase_candidate() {
    const auto parameters = cpa_max3_test::parameters();
    const auto model = th::CpaPtPhase::from_parameters(parameters);
    const auto phases = cpa_max3_test::starts();
    const auto feed = cpa_max3_test::feed_from_phase_fractions(
        {1.0 / 3.0, 1.0 / 3.0, 1.0 / 3.0});
    fl::CpaThreePhaseEvaluator evaluator(
        model,
        {fl::CpaRootSide::upper_density_admissible,
         fl::CpaRootSide::lower_density_admissible,
         fl::CpaRootSide::upper_density_admissible},
        cpa_max3_test::fast_pt_options());
    const auto result = fl::iterate_pt_three_phase(
        cpa_max3_test::pressure_pa, cpa_max3_test::temperature_k,
        feed, log_ratio(phases[1], phases[0]),
        log_ratio(phases[2], phases[0]),
        {1.0 / 3.0, 1.0 / 3.0}, evaluator);
    require(result.candidate_admissible(),
            "CPA structural fixed-three-phase equations did not converge");
    require(result.point &&
                result.point->chemical_potential_norm <=
                    result.options.chemical_potential_tolerance &&
                result.point->mass_absolute <=
                    result.options.mass_absolute_tolerance,
            "CPA fixed-three-phase candidate lost equilibrium or balance closure");
}

fl::CpaPtMax3Result solve_three_phase_fraction(
    std::array<double, 3> beta) {
    const auto parameters = cpa_max3_test::parameters();
    const auto model = th::CpaPtPhase::from_parameters(parameters);
    fl::CpaVleEvaluator evaluator(model, cpa_max3_test::fast_pt_options());
    const auto starts = cpa_max3_test::starts();
    const auto feed = cpa_max3_test::feed_from_phase_fractions(beta);
    fl::CpaPtMax3Options options;
    options.two_phase.initial_stability.automatic_starts = false;
    options.two_phase.final_stability.automatic_starts = false;
    options.final_three_phase_stability.automatic_starts = false;
    options.three_phase_starts.push_back(cpa_max3_test::three_phase_start(beta));
    return fl::solve_cpa_pt_max3(
        cpa_max3_test::pressure_pa, cpa_max3_test::temperature_k,
        feed, evaluator, options, starts, starts);
}

void two_to_three_baseline() {
    constexpr std::array<double, 3> beta{
        1.0 / 3.0, 1.0 / 3.0, 1.0 / 3.0};
    const auto result = solve_three_phase_fraction(beta);
    require(result.base.solution.status == fl::PtSplitStatus::phase_set_unstable,
            "CPA max3 fixture no longer starts from rejected two-phase topology");
    require(result.status == fl::CpaPtMax3Status::three_phase &&
                result.three_phase_candidate() != nullptr &&
                result.selected_attempt.has_value(),
            "CPA max3 did not accept the structural three-phase topology");
    const auto& attempt = result.attempts[*result.selected_attempt];
    require(attempt.supplied_start && *attempt.supplied_start == 0U &&
                attempt.equilibrium.candidate_admissible() &&
                attempt.final_stability &&
                attempt.final_stability->status ==
                    fl::StabilityStatus::no_instability_found &&
                attempt.accepted_three_phase,
            "CPA 2->3 route did not close equilibrium plus final all-root stability");
    const auto published = fl::project_cpa_pt_max3_phase_set(result);
    const auto* phases = published.solution.accepted_phase_set();
    require(phases != nullptr && phases->phases.size() == 3U &&
                published.solution.capability.maximum_phase_count == 3U,
            "CPA max3 publication lost accepted three-phase state");
}

void diagnostic_disappearance_band() {
    constexpr std::array<double, 7> eps_values{
        1.0e-6, 1.0e-8, 5.0e-10, 1.0e-10,
        5.0e-11, 1.0e-11, 1.0e-12};
    std::cout << "CPA_MAX3_BOUNDARY";
    for (const double epsilon : eps_values) {
        const double remaining = 1.0 - epsilon;
        const std::array<double, 3> beta{
            0.5 * remaining, 0.5 * remaining, epsilon};
        const auto result = solve_three_phase_fraction(beta);
        std::cout << ' ' << epsilon
                  << ":base=" << static_cast<int>(result.base.solution.status)
                  << ",max3=" << static_cast<int>(result.status)
                  << ",attempts=" << result.attempts.size();
        if (result.selected_attempt) {
            const auto& attempt = result.attempts[*result.selected_attempt];
            std::cout << ",eq=" << static_cast<int>(attempt.equilibrium.status)
                      << ",neighbor="
                      << (attempt.boundary_neighbor
                          ? static_cast<int>(attempt.boundary_neighbor->solution.status)
                          : -1);
        }
    }
    std::cout << '\n';
}

using Test = std::pair<std::string_view, void (*)()>;
constexpr Test tests[]{
    {"fixed_three_phase", fixed_three_phase_candidate},
    {"two_to_three", two_to_three_baseline},
    {"diagnostic_disappearance", diagnostic_disappearance_band}};

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
