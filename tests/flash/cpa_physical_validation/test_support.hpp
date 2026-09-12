#ifndef MPMC_TEST_CPA_PHYSICAL_VALIDATION_SUPPORT_HPP
#define MPMC_TEST_CPA_PHYSICAL_VALIDATION_SUPPORT_HPP

#include <mpmc/thermodynamics/cpa_parameters.hpp>

#include <array>
#include <optional>
#include <string>
#include <utility>
#include <vector>

namespace cpa_physical_test {
namespace th = mpmc::thermodynamics;
using Vec = std::vector<double>;

inline constexpr double temperature_k = 333.15;
inline constexpr double kij_methanol_water = -0.055;

inline th::Provenance source(
    std::string reference, std::string revision,
    std::string locator, std::string acquisition) {
    return {th::SourceKind::literature,
            std::move(reference), std::move(revision), std::move(locator), "",
            std::move(acquisition),
            "Numerical transcription for repository validation; copyright remains with the cited source."};
}

inline th::Provenance cpa_parameter_source(std::string locator) {
    return source(
        "Kontogeorgis et al., Oil & Gas Science and Technology 63 (2008) 305-319, DOI 10.2516/ogst:2008025",
        "published 2008", std::move(locator),
        "Transcribed from Table 1a / CPA equations in the cited open-access review");
}

inline th::Provenance cr1_source(std::string locator) {
    return source(
        "Kontogeorgis et al., Oil & Gas Science and Technology 63 (2008) 305-319, DOI 10.2516/ogst:2008025",
        "published 2008", std::move(locator),
        "CR-1 cross-association rule evaluated explicitly from the cited pure parameters");
}

inline th::Provenance kij_source() {
    return source(
        "Georgios K. Folas, Modeling of Complex Mixtures Containing Hydrogen Bonding Molecules, PhD thesis, DTU",
        "Table 2.1", "water-methanol, 333.15 K, CR-1: k12=-0.055",
        "Transcribed from the thesis VLE correlation table; experimental reference is Kurihara et al. 1995");
}

inline th::Provenance experiment_source(std::string locator) {
    return source(
        "K. Kurihara, T. Minoura, K. Takeda, K. Kojima, J. Chem. Eng. Data 40 (1995) 679-684, DOI 10.1021/je00019a033",
        "published 1995", std::move(locator),
        "Experimental P-x-y value transcribed from the published 333.15 K methanol-water dataset");
}

inline th::CpaSourcedValue value(double number, th::Provenance provenance) {
    return {number, std::move(provenance)};
}

inline th::Component component(std::string id, std::string display) {
    return {std::move(id), std::move(display), th::ComponentKind::pure,
            cpa_parameter_source("Table 1a chemical identity"), std::nullopt};
}

inline th::CpaParameterSet parameters(bool swapped = false) {
    const std::vector<th::Component> catalog{
        component("METHANOL", "Methanol"), component("WATER", "Water")};
    const std::vector<std::string> order = swapped
        ? std::vector<std::string>{"WATER", "METHANOL"}
        : std::vector<std::string>{"METHANOL", "WATER"};

    th::CpaParameterInput input;
    input.dataset_id = "literature-cpa-water-methanol-cr1-333.15K";
    input.revision = "Kontogeorgis-2008-pure__Folas-CR1-kij__Kurihara-1995-VLE";
    input.applicability = {
        th::ClosedInterval{temperature_k, temperature_k},
        th::ClosedInterval{3.9223e4, 7.2832e4},
        source(
            "K. Kurihara et al., J. Chem. Eng. Data 40 (1995) 679-684, DOI 10.1021/je00019a033",
            "published 1995", "333.15 K methanol-water experimental range used by this fixture",
            "Bounds restricted to the interior experimental points selected for validation")};

    // Table 1a units: a0 [bar L^2 mol^-2] -> multiply by 0.1 for Pa m^6 mol^-2;
    // b [L mol^-1] -> multiply by 1e-3; epsilon [bar L mol^-1] -> multiply by 100 J/mol.
    input.pure.push_back({
        "METHANOL",
        value(512.64, cpa_parameter_source("Table 1a: methanol Tc=512.64 K")),
        value(0.40531, cpa_parameter_source("Table 1a: methanol a0=4.0531 bar L^2/mol^2 -> 0.40531 Pa m^6/mol^2")),
        value(3.0978e-5, cpa_parameter_source("Table 1a: methanol b=0.030978 L/mol -> 3.0978e-5 m^3/mol")),
        value(0.43102, cpa_parameter_source("Table 1a: methanol c1=0.43102")),
        {{"H", 1U}, {"e", 1U}}});
    input.pure.push_back({
        "WATER",
        value(647.29, cpa_parameter_source("Table 1a: water Tc=647.29 K")),
        value(0.12277, cpa_parameter_source("Table 1a: water a0=1.2277 bar L^2/mol^2 -> 0.12277 Pa m^6/mol^2")),
        value(1.4515e-5, cpa_parameter_source("Table 1a: water b=0.014515 L/mol -> 1.4515e-5 m^3/mol")),
        value(0.67359, cpa_parameter_source("Table 1a: water c1=0.67359")),
        {{"H", 2U}, {"e", 2U}}});

    input.binary.push_back({
        "METHANOL", "WATER", value(kij_methanol_water, kij_source())});

    // Self association, 2B methanol and 4C water.
    input.association_pairs.push_back({
        "METHANOL", "H", "METHANOL", "e",
        value(24591.0, cpa_parameter_source("Table 1a: methanol epsilon=245.91 bar L/mol -> 24591 J/mol")),
        value(0.0161, cpa_parameter_source("Table 1a: methanol beta=0.0161"))});
    input.association_pairs.push_back({
        "WATER", "H", "WATER", "e",
        value(16655.0, cpa_parameter_source("Table 1a: water epsilon=166.55 bar L/mol -> 16655 J/mol")),
        value(0.0692, cpa_parameter_source("Table 1a: water beta=0.0692"))});

    // CR-1 is represented explicitly; no hidden mixing rule is used by production code.
    // epsilon_cross=(245.91+166.55)/2=206.23 bar L/mol;
    // beta_cross=sqrt(0.0161*0.0692)=0.03337843615270194.
    for (const auto& sites : std::array<std::pair<const char*, const char*>, 2>{{
             {"H", "e"}, {"e", "H"}}}) {
        input.association_pairs.push_back({
            "METHANOL", sites.first, "WATER", sites.second,
            value(20623.0, cr1_source("CR-1 epsilon_cross arithmetic mean = 206.23 bar L/mol -> 20623 J/mol")),
            value(0.03337843615270194,
                  cr1_source("CR-1 beta_cross geometric mean sqrt(0.0161*0.0692)"))});
    }

    return th::CpaParameterSet::create(catalog, order, input);
}

struct ExperimentalVlePoint {
    double pressure_pa{};
    double liquid_methanol{};
    double vapor_methanol{};
    th::Provenance source;
};

inline const std::array<ExperimentalVlePoint, 5>& points() {
    static const std::array<ExperimentalVlePoint, 5> values{{
        {3.9223e4, 0.1686, 0.5714, experiment_source("333.15 K: P=39.223 kPa, xMeOH=0.1686, yMeOH=0.5714")},
        {4.8852e4, 0.3039, 0.6943, experiment_source("333.15 K: P=48.852 kPa, xMeOH=0.3039, yMeOH=0.6943")},
        {5.6652e4, 0.4461, 0.7742, experiment_source("333.15 K: P=56.652 kPa, xMeOH=0.4461, yMeOH=0.7742")},
        {6.3998e4, 0.6044, 0.8383, experiment_source("333.15 K: P=63.998 kPa, xMeOH=0.6044, yMeOH=0.8383")},
        {7.2832e4, 0.7776, 0.9141, experiment_source("333.15 K: P=72.832 kPa, xMeOH=0.7776, yMeOH=0.9141")}
    }};
    return values;
}

inline Vec composition(double methanol, bool swapped = false) {
    return swapped ? Vec{1.0 - methanol, methanol} : Vec{methanol, 1.0 - methanol};
}

inline Vec feed(const ExperimentalVlePoint& point, bool swapped = false) {
    constexpr double vapor_fraction = 0.1;
    const double methanol =
        (1.0 - vapor_fraction) * point.liquid_methanol +
        vapor_fraction * point.vapor_methanol;
    return composition(methanol, swapped);
}

inline std::vector<Vec> starts(const ExperimentalVlePoint& point, bool swapped = false) {
    return {composition(point.liquid_methanol, swapped),
            composition(point.vapor_methanol, swapped)};
}

} // namespace cpa_physical_test

#endif // MPMC_TEST_CPA_PHYSICAL_VALIDATION_SUPPORT_HPP
