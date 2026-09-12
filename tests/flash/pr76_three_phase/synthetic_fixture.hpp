#ifndef MPMC_TEST_PR76_THREE_PHASE_SYNTHETIC_FIXTURE_HPP
#define MPMC_TEST_PR76_THREE_PHASE_SYNTHETIC_FIXTURE_HPP

#include <mpmc/thermodynamics/pr76_phase.hpp>

#include <array>
#include <optional>
#include <string>
#include <utility>
#include <vector>

namespace pr76_max3_test {
namespace th = mpmc::thermodynamics;
using Vec = std::vector<double>;

inline th::Provenance source(std::string locator) {
    return {th::SourceKind::synthetic_test,
            "MPMC_HNU symmetric PR76 max3 structural fixture",
            "v1",
            std::move(locator),
            "Artificial PR parameters used only to exercise generic three-phase topology; not experimental validation",
            "tests/flash/pr76_three_phase/synthetic_fixture.hpp",
            "Repository structural regression"};
}

inline th::SourcedScalar scalar(double value, th::Unit unit,
                                const th::Provenance& provenance) {
    return {value, unit, provenance,
            "synthetic SI or dimensionless", "identity"};
}

inline th::Pr76Phase<double> model(bool swap_heavy = false) {
    const auto provenance = source("three-component symmetric max3 fixture");
    const std::array<std::string, 3> ids{"light", "heavy-b", "heavy-c"};
    const std::array<double, 3> tc{190.6, 500.0, 500.0};
    const std::array<double, 3> pc{4.6e6, 5.0e6, 5.0e6};
    const std::array<double, 3> omega{0.01, 0.10, 0.10};

    std::vector<th::Component> catalog;
    th::PrParameterInput input;
    input.model_id = std::string(th::pr76_profile);
    input.dataset_id = "synthetic-PR76-max3-symmetric";
    input.revision = "v1";
    input.applicability = {std::nullopt, std::nullopt, provenance};
    for (std::size_t i = 0; i < ids.size(); ++i) {
        catalog.push_back({ids[i], ids[i], th::ComponentKind::pure, provenance, {}});
        input.pure.push_back({
            ids[i], scalar(tc[i], th::Unit::kelvin, provenance),
            scalar(pc[i], th::Unit::pascal, provenance),
            scalar(omega[i], th::Unit::dimensionless, provenance)});
    }
    input.binary.push_back({
        ids[0], ids[1], scalar(0.05, th::Unit::dimensionless, provenance)});
    input.binary.push_back({
        ids[0], ids[2], scalar(0.05, th::Unit::dimensionless, provenance)});
    input.binary.push_back({
        ids[1], ids[2], scalar(0.20, th::Unit::dimensionless, provenance)});

    const std::vector<std::string> order = swap_heavy
        ? std::vector<std::string>{ids[0], ids[2], ids[1]}
        : std::vector<std::string>{ids.begin(), ids.end()};
    return th::Pr76Phase<double>::from_parameters(
        th::PrParameterSet::create(
            catalog, order, input, th::DataPolicy::allow_synthetic_tests));
}

inline const std::array<Vec, 3>& reference_phases() {
    // Independent offline root/equality solve of this explicitly synthetic PR76
    // fixture at p=1 MPa, T=250 K. Structural numerical anchors only.
    static const std::array<Vec, 3> values{{
        {0.04844395156186141, 0.9076712804842868, 0.0438847679538518},
        {0.9696550397025219, 0.015172480148739028, 0.015172480148739028},
        {0.048443951561861175, 0.04388476795385164, 0.9076712804842871}}};
    return values;
}

inline std::vector<Vec> starts(bool swap_heavy = false) {
    std::vector<Vec> values;
    values.reserve(reference_phases().size());
    for (const auto& phase : reference_phases()) {
        if (swap_heavy) {
            values.push_back({phase[0], phase[2], phase[1]});
        } else {
            values.push_back(phase);
        }
    }
    return values;
}

inline Vec equal_feed() {
    return {1.0 / 3.0, 1.0 / 3.0, 1.0 / 3.0};
}

} // namespace pr76_max3_test

#endif // MPMC_TEST_PR76_THREE_PHASE_SYNTHETIC_FIXTURE_HPP
