#ifndef MPMC_TEST_CPA_MAX3_SUPPORT_HPP
#define MPMC_TEST_CPA_MAX3_SUPPORT_HPP

#include <mpmc/flash/cpa_three_phase.hpp>

#include <array>
#include <optional>
#include <string>
#include <utility>
#include <vector>

namespace cpa_max3_test {
namespace th = mpmc::thermodynamics;
namespace fl = mpmc::flash;
using Vec = std::vector<double>;

inline constexpr double pressure_pa = 1.0e6;
inline constexpr double temperature_k = 250.0;

inline th::Provenance source(std::string locator) {
    return {th::SourceKind::synthetic_test,
            "MPMC_HNU CPA max3 symmetric structural fixture",
            "v1", std::move(locator),
            "Artificial non-associating CPA/SRK-limit values exercise max3 topology only; not experimental validation",
            "tests/flash/cpa_max3/test_support.hpp",
            "Repository structural regression"};
}

inline th::CpaSourcedValue value(double number, std::string locator) {
    return {number, source(std::move(locator))};
}

inline th::Component component(std::string id) {
    const auto provenance = source("component-" + id);
    return {id, id, th::ComponentKind::pure, provenance, std::nullopt};
}

inline th::CpaParameterSet parameters(bool swap_heavy = false) {
    const std::vector<th::Component> catalog{
        component("A"), component("B"), component("C")};
    const std::vector<std::string> order = swap_heavy
        ? std::vector<std::string>{"A", "C", "B"}
        : std::vector<std::string>{"A", "B", "C"};

    th::CpaParameterInput input;
    input.dataset_id = "synthetic-cpa-max3-symmetric";
    input.revision = "v1";
    input.applicability = {
        std::nullopt, std::nullopt, source("dataset-applicability")};
    input.pure.push_back({
        "A", value(300.0, "Tc-A"),
        value(1.4837953514927054, "a0-A"),
        value(1.0e-5, "b-A"), value(0.5, "c1-A"), {}});
    input.pure.push_back({
        "B", value(600.0, "Tc-B"),
        value(1.0, "a0-B"),
        value(2.0e-5, "b-B"), value(0.8, "c1-B"), {}});
    input.pure.push_back({
        "C", value(600.0, "Tc-C"),
        value(1.0, "a0-C"),
        value(2.0e-5, "b-C"), value(0.8, "c1-C"), {}});
    input.binary.push_back({
        "A", "B", value(4.7164403200255292, "kij-A-B")});
    input.binary.push_back({
        "A", "C", value(4.7164403200255292, "kij-A-C")});
    input.binary.push_back({
        "B", "C", value(4.6708058982499905, "kij-B-C")});

    return th::CpaParameterSet::create(
        catalog, order, input, th::DataPolicy::allow_synthetic_tests);
}

inline th::CpaPtOptions fast_pt_options() {
    th::CpaPtOptions options;
    options.scan_intervals = 128U;
    options.max_evaluations = 2048U;
    return options;
}

// Independent offline solution of the synthetic SRK-limit common-tangent
// equations. Production CPA kernels must independently re-evaluate and accept it.
inline const std::array<Vec, 3>& reference_phases() {
    static const std::array<Vec, 3> phases{{
        {0.1, 0.8, 0.1},
        {0.8, 0.1, 0.1},
        {0.1, 0.1, 0.8}}};
    return phases;
}

inline std::vector<Vec> starts(bool swap_heavy = false) {
    std::vector<Vec> result;
    result.reserve(3U);
    for (const auto& phase : reference_phases()) {
        result.push_back(swap_heavy
            ? Vec{phase[0], phase[2], phase[1]}
            : phase);
    }
    return result;
}

inline Vec feed_from_phase_fractions(
    std::array<double, 3> beta,
    bool swap_heavy = false) {
    Vec result(3U, 0.0);
    const auto& phases = reference_phases();
    for (std::size_t phase = 0; phase < 3U; ++phase) {
        for (std::size_t i = 0; i < 3U; ++i) {
            result[i] += beta[phase] * phases[phase][i];
        }
    }
    if (swap_heavy) { std::swap(result[1], result[2]); }
    return result;
}

inline fl::CpaPtThreePhaseStart three_phase_start(
    std::array<double, 3> beta = {1.0 / 3.0, 1.0 / 3.0, 1.0 / 3.0},
    bool swap_heavy = false) {
    fl::CpaPtThreePhaseStart start;
    const auto phases = starts(swap_heavy);
    start.compositions = {phases[0], phases[1], phases[2]};
    start.phase_fraction_seed = {beta[1], beta[2]};
    return start;
}

} // namespace cpa_max3_test

#endif // MPMC_TEST_CPA_MAX3_SUPPORT_HPP
