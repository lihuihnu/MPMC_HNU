#include <mpmc/flash/cpa_split.hpp>

#include "../cpa_stability/test_support.hpp"
#include "../../thermodynamics/cpa_pt_phase/test_support.hpp"

#include <algorithm>
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

void print_split_summary(const char* label, const fl::CpaPtSplitResult& result) {
    std::cout << "CPA_SPLIT_DIAGNOSTIC " << label
              << " status=" << static_cast<int>(result.solution.status)
              << " initial=" << static_cast<int>(result.solution.initial_stability.status)
              << " attempts=" << result.solution.attempts.size()
              << " selected=" << (result.solution.selected_attempt.has_value() ? 1 : 0);
    for (const auto& attempt : result.solution.attempts) {
        std::cout << " attempt=" << static_cast<int>(attempt.status)
                  << ",iter=" << attempt.iterations
                  << ",eval=" << attempt.evaluations
                  << ",reject=" << attempt.rejected_evaluations;
        if (attempt.point) {
            std::cout << ",fug=" << attempt.point->fugacity_norm;
        }
    }
    if (result.solution.final_stability) {
        std::cout << " final="
                  << static_cast<int>(result.solution.final_stability->status);
    } else {
        std::cout << " final=none";
    }
    std::cout << " diag=" << result.solution.diagnostic << '\n';
}

void explicit_density_sides() {
    const auto parameters = cpa_stability_test::pure();
    const auto model = th::CpaPtPhase::from_parameters(parameters);
    fl::CpaVleEvaluator evaluator(model);
    const Vec x{1.0};
    const auto roots = model.roots(1.0e5, 200.0, x);
    require(roots.status == th::CpaPtRootStatus::success && roots.roots.size() == 3U,
            "CPA explicit-side fixture lost its three density roots");

    const auto liquid = evaluator(
        1.0e5, 200.0, x, fl::PtPhaseRole::liquid_candidate);
    const auto vapor = evaluator(
        1.0e5, 200.0, x, fl::PtPhaseRole::vapor_candidate);
    require(liquid.activity.branch != vapor.activity.branch,
            "CPA explicit split sides collapsed onto the same root in a three-root state");
    require(vapor.z > liquid.z,
            "CPA lower-density candidate does not have larger compressibility factor");
    require(roots.roots[liquid.activity.branch].molar_density_mol_per_m3 >
                roots.roots[vapor.activity.branch].molar_density_mol_per_m3,
            "CPA split role-to-density-side mapping changed");
}

fl::CpaPtSplitResult run_binary(bool swapped) {
    const auto parameters = cpa_stability_test::binary(swapped);
    const auto model = th::CpaPtPhase::from_parameters(parameters);
    fl::CpaVleEvaluator evaluator(model, cpa_stability_test::fast_pt_options());
    const Vec feed{0.5, 0.5};
    const Vec witness = swapped ? Vec{0.22, 0.78} : Vec{0.78, 0.22};
    const std::vector<Vec> starts{witness};
    fl::PtSplitOptions options;
    options.initial_stability.automatic_starts = false;
    options.initial_stability.max_evaluations = 32U;
    options.final_stability.automatic_starts = false;
    options.final_stability.max_evaluations = 64U;
    options.iteration.max_evaluations = 100000U;
    return fl::solve_cpa_pt_vle(
        3.0e6, 180.0, feed, evaluator, options, starts, starts);
}

void binary_vle_baseline() {
    const auto result = run_binary(false);
    print_split_summary("binary", result);
    require(result.solution.initial_stability.status == fl::StabilityStatus::unstable,
            "CPA VLE fixture no longer starts from robust instability evidence");
    require(result.solution.status == fl::PtSplitStatus::two_phase_no_instability_found &&
                result.solution.candidate() != nullptr &&
                result.solution.final_stability.has_value() &&
                result.solution.final_stability->status ==
                    fl::StabilityStatus::no_instability_found,
            "CPA VLE baseline did not close through split and final common-tangent review");
}

void component_permutation() {
    const auto first = run_binary(false);
    const auto second = run_binary(true);
    print_split_summary("permutation-a", first);
    print_split_summary("permutation-b", second);
    require(first.solution.status == fl::PtSplitStatus::two_phase_no_instability_found &&
                second.solution.status == fl::PtSplitStatus::two_phase_no_instability_found &&
                first.solution.candidate() && second.solution.candidate(),
            "CPA VLE component permutation changed accepted topology");
}

void associating_phase_provider_smoke() {
    const auto parameters = cpa_pt_test::associating_binary(false);
    const auto model = th::CpaPtPhase::from_parameters(parameters);
    fl::CpaVleEvaluator evaluator(model);
    const Vec x{0.7, 0.3};
    const auto roots = model.roots(1.0e6, 330.0, x);
    require(roots.status == th::CpaPtRootStatus::success && !roots.roots.empty(),
            "associating CPA split smoke state has no resolved root");
    const auto low = evaluator(
        1.0e6, 330.0, x, fl::PtPhaseRole::vapor_candidate);
    const auto high = evaluator(
        1.0e6, 330.0, x, fl::PtPhaseRole::liquid_candidate);
    require(low.activity.ln_phi.size() == x.size() && high.activity.ln_phi.size() == x.size(),
            "associating CPA split evaluator lost fugacity dimension");
    for (double value : low.activity.ln_phi) {
        require(std::isfinite(value), "nonfinite associating vapor-side ln(phi)");
    }
    for (double value : high.activity.ln_phi) {
        require(std::isfinite(value), "nonfinite associating liquid-side ln(phi)");
    }
}

void diagnostic_vle_grid() {
    const auto parameters = cpa_stability_test::binary(false);
    const auto model = th::CpaPtPhase::from_parameters(parameters);
    fl::CpaVleEvaluator evaluator(model, cpa_stability_test::fast_pt_options());
    const Vec feed{0.5, 0.5};
    const std::vector<Vec> starts{{0.95, 0.05}, {0.78, 0.22}, {0.22, 0.78}, {0.05, 0.95}};
    constexpr std::array<double, 6> temperatures{160.0, 180.0, 200.0, 220.0, 240.0, 260.0};
    constexpr std::array<double, 6> pressures_mpa{0.10, 0.30, 0.60, 1.00, 2.00, 3.00};

    std::cout << "CPA_VLE_GRID";
    for (const double temperature : temperatures) {
        for (const double pressure_mpa : pressures_mpa) {
            fl::PtSplitOptions options;
            options.initial_stability.automatic_starts = false;
            options.initial_stability.max_evaluations = 96U;
            options.final_stability.automatic_starts = false;
            options.final_stability.max_evaluations = 128U;
            const auto result = fl::solve_cpa_pt_vle(
                pressure_mpa * 1.0e6, temperature, feed, evaluator,
                options, starts, starts);
            std::cout << ' ' << temperature << "K@" << pressure_mpa << "MPa:"
                      << static_cast<int>(result.solution.initial_stability.status)
                      << '/' << static_cast<int>(result.solution.status);
            if (result.solution.selected_attempt && result.solution.candidate()) {
                std::cout << "/C";
            }
            if (result.solution.final_stability) {
                std::cout << "/F" << static_cast<int>(result.solution.final_stability->status);
            }
        }
    }
    std::cout << '\n';
}

using Test = std::pair<std::string_view, void (*)()>;
constexpr Test tests[]{
    {"explicit_density_sides", explicit_density_sides},
    {"binary_vle", binary_vle_baseline},
    {"component_permutation", component_permutation},
    {"associating_provider", associating_phase_provider_smoke},
    {"diagnostic_vle_grid", diagnostic_vle_grid}};

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
