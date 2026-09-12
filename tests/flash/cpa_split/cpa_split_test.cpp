#include <mpmc/flash/cpa_split.hpp>

#include "test_support.hpp"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <stdexcept>
#include <string_view>
#include <vector>

namespace {
namespace fl = mpmc::flash;
namespace th = mpmc::thermodynamics;
using cpa_split_test::Vec;

void require(bool value, const char* message) {
    if (!value) { throw std::runtime_error(message); }
}

void print_attempt(const fl::PtSplitAttempt& attempt, std::string_view prefix) {
    std::cerr << prefix
              << " status=" << static_cast<int>(attempt.status)
              << " iterations=" << attempt.iterations
              << " evaluations=" << attempt.evaluations
              << " diagnostic=" << attempt.diagnostic;
    if (attempt.point) {
        std::cerr << " fugacity_norm=" << attempt.point->fugacity_norm
                  << " beta=" << attempt.point->fractions.vapor_fraction
                  << " mass_abs=" << attempt.point->fractions.mass_absolute
                  << " mass_rel=" << attempt.point->fractions.mass_relative
                  << " zL=" << attempt.point->liquid.z
                  << " zV=" << attempt.point->vapor.z;
        if (attempt.point->fractions.liquid.size() >= 2U &&
            attempt.point->fractions.vapor.size() >= 2U) {
            std::cerr << " x=(" << attempt.point->fractions.liquid[0] << ','
                      << attempt.point->fractions.liquid[1] << ')'
                      << " y=(" << attempt.point->fractions.vapor[0] << ','
                      << attempt.point->fractions.vapor[1] << ')';
        }
    }
    std::cerr << '\n';
}

void print_solution(const fl::CpaPtSplitResult& result, std::string_view prefix) {
    std::cerr << prefix
              << " solution_status=" << static_cast<int>(result.solution.status)
              << " initial_status="
              << static_cast<int>(result.solution.initial_stability.status)
              << " attempts=" << result.solution.attempts.size()
              << " selected=" << (result.solution.selected_attempt ? 1 : 0)
              << " diagnostic=" << result.solution.diagnostic;
    if (result.solution.final_stability) {
        std::cerr << " final_status="
                  << static_cast<int>(result.solution.final_stability->status);
    }
    std::cerr << '\n';
    for (std::size_t i = 0; i < result.solution.attempts.size(); ++i) {
        print_attempt(result.solution.attempts[i],
                      std::string(prefix) + ".attempt[" + std::to_string(i) + "]");
    }
}

void density_side_roles() {
    const auto parameters = cpa_split_test::pure();
    const auto model = th::CpaPtPhase::from_parameters(parameters);
    fl::CpaVleEvaluator evaluator(model);
    const Vec x{1.0};
    const auto roots = model.roots(1.0e5, 200.0, x, evaluator.pt_options());
    require(roots.status == th::CpaPtRootStatus::success &&
                roots.roots.size() == 3U,
            "CPA split role fixture no longer has three simple roots");

    const auto liquid = evaluator(
        1.0e5, 200.0, x, fl::PtPhaseRole::liquid_candidate);
    const auto vapor = evaluator(
        1.0e5, 200.0, x, fl::PtPhaseRole::vapor_candidate);
    require(liquid.activity.branch != vapor.activity.branch &&
                liquid.z < vapor.z,
            "CPA density-side candidate roles did not select distinct dense/light roots");
    require(roots.roots[liquid.activity.branch].pressure_slope_sign > 0 &&
                roots.roots[vapor.activity.branch].pressure_slope_sign > 0,
            "CPA split selected a mechanically inadmissible density root");

    double maximum_density = 0.0;
    double minimum_density = 0.0;
    bool first = true;
    for (const auto& root : roots.roots) {
        if (root.pressure_slope_sign <= 0) { continue; }
        if (first) {
            maximum_density = minimum_density = root.molar_density_mol_per_m3;
            first = false;
        } else {
            maximum_density = std::max(maximum_density, root.molar_density_mol_per_m3);
            minimum_density = std::min(minimum_density, root.molar_density_mol_per_m3);
        }
    }
    require(std::abs(
                roots.roots[liquid.activity.branch].molar_density_mol_per_m3 -
                maximum_density) < 1e-12 &&
            std::abs(
                roots.roots[vapor.activity.branch].molar_density_mol_per_m3 -
                minimum_density) < 1e-12,
            "CPA split role selection is not the requested density side");
}

void fixed_seed_reference() {
    const auto parameters = cpa_split_test::binary(false);
    const auto model = th::CpaPtPhase::from_parameters(parameters);
    fl::CpaVleEvaluator evaluator(
        model, cpa_split_test::fast_pt_options());
    const Vec feed{cpa_split_test::feed.begin(), cpa_split_test::feed.end()};
    const auto log_k = cpa_split_test::reference_log_k(false);
    const auto attempt = fl::iterate_pt_split(
        cpa_split_test::pressure_pa,
        cpa_split_test::temperature_k,
        feed, log_k, evaluator);
    if (attempt.status != fl::PtSplitAttemptStatus::converged ||
        !attempt.point.has_value()) {
        print_attempt(attempt, "fixed_seed_reference");
        throw std::runtime_error("CPA fixed-seed two-phase equations did not converge");
    }
    const auto& state = *attempt.point;
    require(state.fugacity_norm <=
                fl::PtSplitIterationOptions{}.fugacity_tolerance &&
                state.fractions.mass_absolute <= 1e-12 &&
                state.fractions.mass_relative <= 1e-10,
            "CPA fixed-seed candidate failed fugacity/material-balance invariants");
    require(std::abs(state.fractions.liquid[0] -
                     cpa_split_test::reference_liquid_a) < 3e-8 &&
                std::abs(state.fractions.vapor[0] -
                         cpa_split_test::reference_vapor_a) < 3e-8 &&
                std::abs(state.fractions.vapor_fraction -
                         cpa_split_test::reference_vapor_fraction) < 3e-8,
            "CPA fixed-seed solution left the independent SRK structural anchor");
}

fl::CpaPtSplitResult run_full(bool swapped) {
    const auto parameters = cpa_split_test::binary(swapped);
    const auto model = th::CpaPtPhase::from_parameters(parameters);
    fl::CpaVleEvaluator evaluator(
        model, cpa_split_test::fast_pt_options());
    const Vec feed{0.5, 0.5};
    const Vec witness = swapped
        ? Vec{cpa_split_test::instability_witness[1],
              cpa_split_test::instability_witness[0]}
        : Vec{cpa_split_test::instability_witness[0],
              cpa_split_test::instability_witness[1]};
    const std::vector<Vec> initial_starts{witness};

    auto options = fl::cpa_pt_vle_default_split_options();
    options.initial_stability.automatic_starts = false;
    options.initial_stability.max_evaluations = 32U;
    options.final_stability.automatic_starts = false;
    options.final_stability.max_evaluations = 64U;
    options.max_split_attempts = 4U;

    return fl::solve_cpa_pt_vle(
        cpa_split_test::pressure_pa,
        cpa_split_test::temperature_k,
        feed, evaluator, options, initial_starts, {});
}

void full_two_phase_acceptance() {
    const auto result = run_full(false);
    require(result.solution.initial_stability.status == fl::StabilityStatus::unstable,
            "CPA full split did not preserve initial negative-TPD evidence");
    if (result.solution.status !=
            fl::PtSplitStatus::two_phase_no_instability_found ||
        result.solution.candidate() == nullptr ||
        !result.solution.final_stability.has_value() ||
        result.solution.final_stability->status !=
            fl::StabilityStatus::no_instability_found) {
        print_solution(result, "full_two_phase");
        throw std::runtime_error(
            "CPA full stability->split->final-review path did not accept two phases");
    }
    require(result.solution.attempts.size() >= 2U,
            "CPA full solve did not attempt both witness role assignments");
    require(result.solution.attempts[0].witness_as_vapor,
            "CPA full solve changed deterministic first role assignment");
    require(result.solution.selected_attempt.has_value() &&
                !result.solution.attempts[*result.solution.selected_attempt].witness_as_vapor,
            "CPA full solve did not select the converged alternate role assignment");
    require(result.solution.options.max_evaluations_per_attempt == 8192U,
            "CPA split did not preserve its per-attempt fairness budget");

    const auto& state = *result.solution.candidate();
    require(state.fugacity_norm <= result.solution.options.iteration.fugacity_tolerance &&
                state.fractions.mass_absolute <=
                    result.solution.options.iteration.mass_absolute_tolerance &&
                state.fractions.mass_relative <=
                    result.solution.options.iteration.mass_relative_tolerance &&
                result.solution.gibbs_change < 0.0,
            "CPA accepted two-phase state failed equilibrium/Gibbs acceptance invariants");
    require(std::abs(state.fractions.liquid[0] -
                     cpa_split_test::reference_liquid_a) < 2e-6 &&
                std::abs(state.fractions.vapor[0] -
                         cpa_split_test::reference_vapor_a) < 2e-6 &&
                std::abs(state.fractions.vapor_fraction -
                         cpa_split_test::reference_vapor_fraction) < 2e-6,
            "CPA full two-phase solve left the independent structural anchor");
    require(result.dataset_id == "synthetic-cpa-split-binary" &&
                result.component_ids == std::vector<std::string>({"A", "B"}) &&
                result.model_profile == std::string(th::cpa_profile) &&
                result.phase_convention == std::string(th::cpa_pt_convention),
            "CPA split result lost model/provenance identity");
}

void component_permutation() {
    const auto first = run_full(false);
    const auto second = run_full(true);
    if (first.solution.status !=
            fl::PtSplitStatus::two_phase_no_instability_found ||
        second.solution.status !=
            fl::PtSplitStatus::two_phase_no_instability_found ||
        first.solution.candidate() == nullptr ||
        second.solution.candidate() == nullptr) {
        print_solution(first, "permutation.first");
        print_solution(second, "permutation.second");
        throw std::runtime_error(
            "CPA component permutation changed accepted two-phase status");
    }
    const auto& a = *first.solution.candidate();
    const auto& b = *second.solution.candidate();
    require(std::abs(a.fractions.vapor_fraction - b.fractions.vapor_fraction) < 2e-9 &&
                std::abs(a.fractions.liquid[0] - b.fractions.liquid[1]) < 2e-9 &&
                std::abs(a.fractions.liquid[1] - b.fractions.liquid[0]) < 2e-9 &&
                std::abs(a.fractions.vapor[0] - b.fractions.vapor[1]) < 2e-9 &&
                std::abs(a.fractions.vapor[1] - b.fractions.vapor[0]) < 2e-9,
            "CPA accepted two-phase state changed under runtime component permutation");
}

using Test = std::pair<std::string_view, void (*)()>;
constexpr Test tests[]{
    {"density_side_roles", density_side_roles},
    {"fixed_seed_reference", fixed_seed_reference},
    {"full_two_phase", full_two_phase_acceptance},
    {"component_permutation", component_permutation}};

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
