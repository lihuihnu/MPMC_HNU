#ifndef MPMC_TEST_CPA_PT_PHASE_SUPPORT_HPP
#define MPMC_TEST_CPA_PT_PHASE_SUPPORT_HPP

#include <mpmc/thermodynamics/cpa_parameters.hpp>

#include <optional>
#include <string>
#include <utility>
#include <vector>

namespace cpa_pt_test {
namespace th = mpmc::thermodynamics;

inline th::Provenance source(std::string locator) {
    return {th::SourceKind::synthetic_test,
            "MPMC_HNU CPA PT structural regression",
            "v1", std::move(locator),
            "Synthetic values exercise CPA PT equations/root topology only; not physical validation",
            "tests/thermodynamics/cpa_pt_phase/test_support.hpp",
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

inline th::CpaParameterSet nonassociating_pure() {
    const std::vector<th::Component> catalog{component("A")};
    const std::vector<std::string> order{"A"};
    th::CpaParameterInput input;
    input.dataset_id = "synthetic-cpa-pt-pure";
    input.revision = "v1";
    input.applicability = applicability();
    input.pure.push_back({
        "A", value(500.0, "Tc-A"), value(0.25, "a0-A"),
        value(4.0e-5, "b-A"), value(0.7, "c1-A"), {}});
    return th::CpaParameterSet::create(
        catalog, order, input, th::DataPolicy::allow_synthetic_tests);
}

inline th::CpaParameterSet associating_binary(bool swapped = false) {
    const std::vector<th::Component> catalog{component("A"), component("B")};
    const std::vector<std::string> order = swapped
        ? std::vector<std::string>{"B", "A"}
        : std::vector<std::string>{"A", "B"};
    th::CpaParameterInput input;
    input.dataset_id = "synthetic-cpa-pt-binary";
    input.revision = "v1";
    input.applicability = applicability();
    input.pure.push_back({
        "A", value(500.0, "Tc-A"), value(0.25, "a0-A"),
        value(4.0e-5, "b-A"), value(0.7, "c1-A"), {{"H", 1U}}});
    input.pure.push_back({
        "B", value(400.0, "Tc-B"), value(0.18, "a0-B"),
        value(5.0e-5, "b-B"), value(0.5, "c1-B"), {}});
    input.binary.push_back({"A", "B", value(0.04, "kij-A-B")});
    input.association_pairs.push_back({
        "A", "H", "A", "H",
        value(10000.0, "epsilon-A-H-A-H"),
        value(0.02, "beta-A-H-A-H")});
    return th::CpaParameterSet::create(
        catalog, order, input, th::DataPolicy::allow_synthetic_tests);
}

} // namespace cpa_pt_test

#endif // MPMC_TEST_CPA_PT_PHASE_SUPPORT_HPP
