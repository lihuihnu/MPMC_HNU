#ifndef MPMC_TEST_CPA_STABILITY_SUPPORT_HPP
#define MPMC_TEST_CPA_STABILITY_SUPPORT_HPP

#include <mpmc/thermodynamics/cpa_pt_phase.hpp>

#include <optional>
#include <string>
#include <utility>
#include <vector>

namespace cpa_stability_test {
namespace th = mpmc::thermodynamics;

inline th::Provenance source(std::string locator) {
    return {th::SourceKind::synthetic_test,
            "MPMC_HNU CPA stability structural regression",
            "v1", std::move(locator),
            "Synthetic values exercise CPA stability/root selection only; not physical validation",
            "tests/flash/cpa_stability/test_support.hpp",
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
    input.dataset_id = "synthetic-cpa-stability-pure";
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
    input.dataset_id = "synthetic-cpa-stability-binary";
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
    return options;
}

} // namespace cpa_stability_test

#endif // MPMC_TEST_CPA_STABILITY_SUPPORT_HPP
