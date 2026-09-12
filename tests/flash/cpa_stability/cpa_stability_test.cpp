#include <mpmc/flash/cpa_stability.hpp>

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
using Vec = std::vector<double>;

void require(bool value, const char* message) {
    if (!value) { throw std::runtime_error(message); }
}

void minimum_gibbs_root_selected() {
    const auto parameters = cpa_stability_test::pure();
    const auto model = th::CpaPtPhase::from_parameters(parameters);
    const Vec x{1.0};
    const auto roots = model.roots(1.0e5, 200.0, x);
    require(roots.status == th::CpaPtRootStatus::success &&
                roots.roots.size() == 3U,
            "CPA stability root-selection fixture no longer has three simple roots");

    std::size_t expected = roots.roots.size();
    double minimum = 0.0;
    for (std::size_t k = 0; k < roots.roots.size(); ++k) {
        if (roots.roots[k].pressure_slope_sign <= 0) { continue; }
        if (expected == roots.roots.size() ||
            roots.roots[k].reduced_gibbs_offset < minimum) {
            expected = k;
            minimum = roots.roots[k].reduced_gibbs_offset;
        }
    }
    require(expected < roots.roots.size(),
            "CPA root-selection fixture has no mechanically admissible root");

    fl::CpaStabilityEvaluator evaluator(model);
    const auto selected = evaluator(1.0e5, 200.0, x);
    require(selected.branch == expected && selected.smooth,
            "CPA stability adapter did not select the minimum-Gibbs admissible root");
    require(selected.ln_phi.size() == 1U &&
                std::abs(selected.ln_phi[0] - roots.roots[expected].ln_phi[0]) < 1e-13,
            "CPA stability adapter changed the selected-root fugacity coefficient");
}

void stable_reference_path() {
    const auto parameters = cpa_stability_test::binary();
    const auto model = th::CpaPtPhase::from_parameters(parameters);
    fl::CpaStabilityEvaluator evaluator(
        model, cpa_stability_test::fast_pt_options());
    const Vec feed{0.5, 0.5};
    const std::vector<Vec> starts{feed};
    fl::StabilityOptions options;
    options.automatic_starts = false;
    options.max_evaluations = 16U;
    const auto result = fl::test_cpa_pt_stability(
        1.0e5, 250.0, feed, evaluator, options, starts);
    require(result.search.status == fl::StabilityStatus::no_instability_found &&
                result.search.reference.has_value() &&
                result.search.trials.size() == 1U &&
                result.search.trials[0].status == fl::StabilityTrialStatus::stationary,
            "CPA stable reference/feed path did not close as a stationary no-instability search");
    require(result.dataset_id == "synthetic-cpa-stability-binary" &&
                result.component_ids == std::vector<std::string>({"A", "B"}) &&
                result.model_profile == std::string(th::cpa_profile) &&
                result.phase_convention == std::string(th::cpa_pt_convention),
            "CPA stability result lost ordered model/provenance identity");
}

fl::CpaStabilityResult run_unstable(bool swapped) {
    const auto parameters = cpa_stability_test::binary(swapped);
    const auto model = th::CpaPtPhase::from_parameters(parameters);
    fl::CpaStabilityEvaluator evaluator(
        model, cpa_stability_test::fast_pt_options());
    const Vec feed{0.5, 0.5};
    const Vec witness = swapped ? Vec{0.22, 0.78} : Vec{0.78, 0.22};
    const std::vector<Vec> starts{witness};
    fl::StabilityOptions options;
    options.automatic_starts = false;
    options.max_evaluations = 16U;
    return fl::test_cpa_pt_stability(
        3.0e6, 180.0, feed, evaluator, options, starts);
}

void negative_tpd_witness() {
    const auto result = run_unstable(false);
    require(result.search.status == fl::StabilityStatus::unstable &&
                result.search.trials.size() == 1U &&
                result.search.trials[0].status == fl::StabilityTrialStatus::negative_tpd &&
                result.search.trials[0].point.has_value(),
            "CPA synthetic unstable state did not retain its negative-TPD witness");
    require(result.search.trials[0].point->value < -1.0e-3,
            "CPA synthetic negative-TPD witness lost robust negative margin");
}

void component_permutation() {
    const auto first = run_unstable(false);
    const auto second = run_unstable(true);
    require(first.search.status == fl::StabilityStatus::unstable &&
                second.search.status == fl::StabilityStatus::unstable &&
                first.search.trials[0].point && second.search.trials[0].point,
            "CPA stability component permutation changed unstable status");
    const auto& a = *first.search.trials[0].point;
    const auto& b = *second.search.trials[0].point;
    require(std::abs(a.value - b.value) < 2.0e-10 &&
                std::abs(a.composition[0] - b.composition[1]) < 2.0e-12 &&
                std::abs(a.composition[1] - b.composition[0]) < 2.0e-12,
            "CPA stability negative witness changed under runtime component permutation");
}

void root_budget_failure_is_indeterminate() {
    const auto parameters = cpa_stability_test::binary();
    const auto model = th::CpaPtPhase::from_parameters(parameters);
    th::CpaPtOptions pt_options;
    pt_options.scan_intervals = 16U;
    pt_options.max_evaluations = 16U;
    fl::CpaStabilityEvaluator evaluator(model, pt_options);
    const Vec feed{0.5, 0.5};
    const std::vector<Vec> starts{feed};
    fl::StabilityOptions options;
    options.automatic_starts = false;
    const auto result = fl::test_cpa_pt_stability(
        1.0e5, 250.0, feed, evaluator, options, starts);
    require(result.search.status == fl::StabilityStatus::indeterminate &&
                result.search.reference_issue ==
                    fl::StabilityPropertyIssue::root_iteration_limit,
            "CPA density-root budget exhaustion was not propagated as indeterminate stability");
    require(result.search.trials.size() == 1U &&
                result.search.trials[0].status ==
                    fl::StabilityTrialStatus::property_failure,
            "CPA reference root failure did not mark trial property failure");
}

using Test = std::pair<std::string_view, void (*)()>;
constexpr Test tests[]{
    {"minimum_gibbs_root", minimum_gibbs_root_selected},
    {"stable_reference", stable_reference_path},
    {"negative_tpd", negative_tpd_witness},
    {"component_permutation", component_permutation},
    {"root_budget_failure", root_budget_failure_is_indeterminate}};

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
