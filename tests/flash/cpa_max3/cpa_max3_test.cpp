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

const char* side_name(fl::CpaRootSide side) {
    return side == fl::CpaRootSide::lower_density_admissible ? "L" : "U";
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
    std::cout << "CPA_MAX3_FIXED status=" << static_cast<int>(result.status)
              << " iter=" << result.iterations
              << " eval=" << result.evaluations;
    if (result.point) {
        std::cout << " mu=" << result.point->chemical_potential_norm
                  << " mass=" << result.point->mass_absolute;
    }
    std::cout << '\n';
    require(result.candidate_admissible(),
            "CPA structural fixed-three-phase equations did not converge");
}

void diagnostic_equal_feed() {
    constexpr std::array<double, 3> beta{
        1.0 / 3.0, 1.0 / 3.0, 1.0 / 3.0};
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
    const auto result = fl::solve_cpa_pt_max3(
        cpa_max3_test::pressure_pa, cpa_max3_test::temperature_k,
        feed, evaluator, options, starts, starts);

    std::cout << "CPA_MAX3_EQUAL base="
              << static_cast<int>(result.base.solution.status)
              << " max3=" << static_cast<int>(result.status)
              << " attempts=" << result.attempts.size();
    if (result.base.solution.final_stability) {
        std::cout << " baseFinal="
                  << static_cast<int>(result.base.solution.final_stability->status);
    }
    for (std::size_t i = 0; i < result.attempts.size(); ++i) {
        const auto& attempt = result.attempts[i];
        std::cout << " a" << i << '='
                  << side_name(attempt.root_sides[0])
                  << side_name(attempt.root_sides[1])
                  << side_name(attempt.root_sides[2])
                  << "/eq" << static_cast<int>(attempt.equilibrium.status);
        if (attempt.equilibrium.point) {
            std::cout << "/mu" << attempt.equilibrium.point->chemical_potential_norm;
        }
        if (attempt.final_stability) {
            std::cout << "/fin" << static_cast<int>(attempt.final_stability->status);
        }
        if (attempt.supplied_start) { std::cout << "/sup" << *attempt.supplied_start; }
        if (attempt.witness_trial) { std::cout << "/wit" << *attempt.witness_trial; }
    }
    std::cout << " diag=" << result.diagnostic << '\n';

    require(result.base.solution.status == fl::PtSplitStatus::phase_set_unstable,
            "CPA max3 equal-feed fixture no longer produces additional-phase evidence");
}

using Test = std::pair<std::string_view, void (*)()>;
constexpr Test tests[]{
    {"fixed_three_phase", fixed_three_phase_candidate},
    {"diagnostic_equal_feed", diagnostic_equal_feed}};

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
