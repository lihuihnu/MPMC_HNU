#ifndef MPMC_TEST_CPA_SPLIT_SUPPORT_HPP
#define MPMC_TEST_CPA_SPLIT_SUPPORT_HPP

#include <mpmc/thermodynamics/cpa_pt_phase.hpp>

#include <array>
#include <cmath>
#include <optional>
#include <string>
#include <utility>
#include <vector>

namespace cpa_split_test {
namespace th = mpmc::thermodynamics;
using Vec = std::vector<double>;

inline th::Provenance source(std::string locator) {
    return {th::SourceKind::synthetic_test,
            "MPMC_HNU CPA two-phase structural regression",
            "v1", std::move(locator),
            "Synthetic values exercise CPA split equations/orchestration only; not physical validation",
            "tests/flash/cpa_split/test_support.hpp",
            "Repository structural regression"};
}

inline th::CpaSourcedValue value(double number, std::string locator) {
    return {number, source(std::move(locator))};
}

inline th::Component component(std::string id) {
    const auto provenance = source("component-" + id);
    return {id, id, th::ComponentKind::pure, provenance, std::nullopt};
}

inline th::Applicability applicability() {
    return {std::nullopt, std::nullopt, source("dataset-applicability")};
}

inline th::CpaParameterSet pure() {
    const std::vector<th::Component> catalog{component("A")};
    const std::vector<std::string> order{"A"};
    th::CpaParameterInput input;
    input.dataset_id = "synthetic-cpa-split-pure";
    input.revision = "v1";
    input.applicability = applicability();
    input.pure.push_back({
        "A", value(500.0, "Tc-A"), value(0.25, "a0-A"),
        value(4.0e-5, "b-A"), value(0.7, "c1-A"), {}});
    return th::CpaParameterSet::create(
        catalog, order, input, th::DataPolicy::allow_synthetic_tests);
}

inline th::CpaParameterSet binary(bool swapped = false) {
    const std::vector<th::Component> catalog{component("A"), component("B")};
    const std::vector<std::string> order = swapped
        ? std::vector<std::string>{"B", "A"}
        : std::vector<std::string>{"A", "B"};
    th::CpaParameterInput input;
    input.dataset_id = "synthetic-cpa-split-binary";
    input.revision = "v1";
    input.applicability = applicability();
    input.pure.push_back({
        "A", value(500.0, "Tc-A"), value(0.25, "a0-A"),
        value(4.0e-5, "b-A"), value(0.7, "c1-A"), {}});
    input.pure.push_back({
        "B", value(400.0, "Tc-B"), value(0.18, "a0-B"),
        value(5.0e-5, "b-B"), value(0.5, "c1-B"), {}});
    input.binary.push_back({"A", "B", value(0.04, "kij-A-B")});
    return th::CpaParameterSet::create(
        catalog, order, input, th::DataPolicy::allow_synthetic_tests);
}

inline th::CpaPtOptions fast_pt_options() {
    th::CpaPtOptions options;
    options.scan_intervals = 128U;
    options.max_evaluations = 2048U;
    // Match the production CPA VLE precision coordination while keeping the
    // structural scan smaller for focused test cost.
    options.pressure_absolute_tolerance_pa = 1.0e-7;
    options.pressure_relative_tolerance = 1.0e-12;
    return options;
}

inline constexpr double pressure_pa = 3.0e6;
inline constexpr double temperature_k = 180.0;
inline constexpr std::array<double, 2> feed{0.5, 0.5};
inline constexpr std::array<double, 2> instability_witness{0.78, 0.22};

// Independent stdlib/SRK structural solve of the same explicitly synthetic
// non-associating CPA limit. These are software regression anchors only.
inline constexpr double reference_liquid_a = 0.72462795705763783;
inline constexpr double reference_vapor_a = 0.46050754265946103;
inline constexpr double reference_vapor_fraction = 0.85047555892062909;

inline Vec reference_log_k(bool swapped = false) {
    const std::array<double, 2> liquid{
        reference_liquid_a, 1.0 - reference_liquid_a};
    const std::array<double, 2> vapor{
        reference_vapor_a, 1.0 - reference_vapor_a};
    Vec values(2U, 0.0);
    if (!swapped) {
        values[0] = std::log(vapor[0]) - std::log(liquid[0]);
        values[1] = std::log(vapor[1]) - std::log(liquid[1]);
    } else {
        values[0] = std::log(vapor[1]) - std::log(liquid[1]);
        values[1] = std::log(vapor[0]) - std::log(liquid[0]);
    }
    return values;
}

} // namespace cpa_split_test

#endif // MPMC_TEST_CPA_SPLIT_SUPPORT_HPP
