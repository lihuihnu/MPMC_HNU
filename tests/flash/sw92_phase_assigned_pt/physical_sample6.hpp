#ifndef MPMC_TEST_SW92_PROFILE_C_PHYSICAL_SAMPLE6_HPP
#define MPMC_TEST_SW92_PROFILE_C_PHYSICAL_SAMPLE6_HPP

#include "test_support.hpp"

#include <algorithm>
#include <array>
#include <cstddef>
#include <optional>
#include <string>
#include <utility>
#include <vector>

namespace sw92_profile_c_sample6 {
namespace th = mpmc::thermodynamics;
using Vec = std::vector<double>;

inline th::Provenance source_2017(std::string locator) {
    return {th::SourceKind::literature,
            "doi:10.1016/j.fluid.2017.07.007",
            "Fluid Phase Equilibria 450 (2017) 160-174",
            std::move(locator),
            "Mortezazadeh-Rasaei synthetic gas-condensate Sample 6 physical three-phase fixture",
            "User PDF mortezazadeh2017.pdf",
            "Published properties/feed/topology only; numerical golden is independently regenerated for the exact MPMC_HNU profile"};
}

inline th::Provenance source_sw92(std::string locator) {
    return {th::SourceKind::literature,
            "doi:10.1016/0378-3812(92)85105-H",
            "Fluid Phase Equilibria 77 (1992) 217-240; supplied authors' errata",
            std::move(locator),
            "Soreide-Whitson corrected-original water/non-water interaction fixture",
            "User PDF sha256:cb5b1d5034d78d934e887449ce0c89692d43431835d1371da5606d95b2c6bf58",
            "Bibliographic/formula facts only; source PDF is not redistributed"};
}

inline th::SourcedScalar scalar(double value, th::Unit unit,
                                const th::Provenance& source,
                                std::string original = "SI or dimensionless",
                                std::string conversion = "identity") {
    return {value, unit, source, std::move(original), std::move(conversion)};
}

struct ComponentSpec {
    const char* id;
    const char* display;
    double tc_k;
    double pc_mpa;
    double omega;
    double na_water_kij;
};

// The NA constants follow Mortezazadeh-Rasaei Appendix A for the physical
// Sample-6 reproduction: C4 uses 0.5091 and normal hydrocarbons heavier than C4
// use the source's 0.5 extension. AQ values still use the corrected SW92 Eq.(12).
inline constexpr std::array<ComponentSpec, 7> hydrocarbons{{
    {"methane", "C1 methane", 190.60, 4.6042, 0.013, 0.4850},
    {"ethane", "C2 ethane", 305.43, 4.8839, 0.0986, 0.4920},
    {"propane", "C3 propane", 369.80, 4.2455, 0.1524, 0.5525},
    {"n-butane", "C4 n-butane", 419.50, 3.7470, 0.1956, 0.5091},
    {"n-pentane", "C5 n-pentane", 465.90, 3.3589, 0.2413, 0.5000},
    {"c6-pseudocomponent", "C6 pseudocomponent", 507.50, 3.0104, 0.2990, 0.5000},
    {"c7plus-pseudocomponent", "C7+ pseudocomponent", 655.04, 2.2305, 0.50879, 0.5000},
}};

inline constexpr const char* water_id = "water";

struct PhysicalGolden {
    std::array<long double, 8> feed;
    std::array<long double, 8> w;
    std::array<long double, 8> h0;
    std::array<long double, 8> h1;
    std::array<long double, 3> fractions;
    std::array<long double, 3> z;
    std::array<long double, 8> common;
};

// Independent Decimal(80) solution of the exact MPMC_HNU
// SW92/corrected-original/PR76-base profile with the physical Sample-6 input
// contract above. H0/H1 are in production canonical order (smaller Z first).
// The literature source is the topology and parameter oracle; these numerical
// values are not copied from its plots.
inline constexpr PhysicalGolden golden{
    {0.5085508550855085508550855086L,
     0.4249424942494249424942494249L,
     0.02140214021402140214021402140L,
     0.01120112011201120112011201120L,
     0.008500850085008500850085008501L,
     0.004000400040004000400040004000L,
     0.002900290029002900290029002900L,
     0.01850185018501850185018501850L},
    {0.9991426078018199301269332960849100300393L,
     0.0008014994453029587077105774192456510410457L,
     0.00004088362931127884293978478724478721262434L,
     0.00001149410256687578429662965127715657029732L,
     0.000003111202876667406662473957362369014889741L,
     0.0000003758633929285473319079278282111901738854L,
     0.00000002794727524508423864472034509314991530382L,
     0.000000000007454115499886685452786701781704380621648L},
    {0.003176089890708091250860113720892617636025L,
     0.3050473617165520640241359674743478262766L,
     0.03971246324651023390153711462312961782113L,
     0.03909922845805758424670191835065623351203L,
     0.04975088579949753885265142859770120236705L,
     0.03754497843180323156971170313929445488879L,
     0.04024213955169295852890240931246936669875L,
     0.4854268529051782976254993447825086807996L},
    {0.005529806496031613918428945505131685934088L,
     0.9042844460694589543494338630230334094172L,
     0.04359956191855830853861302314080199233566L,
     0.02136718437287454972361030761330521910263L,
     0.01461944317327349730630759077604231881256L,
     0.005752607150691583412215712685023015544716L,
     0.003131256469128431384700610867088592316622L,
     0.001715694349983061366689946390573766536550L},
    {0.5063410573883960722597623619222789213811L,
     0.03649881134708674580877254694401208320018L,
     0.4571601312645171819314650911337089954187L},
    {0.07587273597589584509672741292093732608749L,
     0.4936622569902828356674452316781931437916L,
     0.8705909142056383455851967863263326083709L},
    {-5.417285239694122061890644893950094802657L,
     -0.2041574361131346621619942807354162138240L,
     -3.541304288072553847268636839999915297573L,
     -4.497495256792170946228460352841902512891L,
     -5.103557255710785914403416316169269724613L,
     -6.275632166602889257243372031367958795463L,
     -7.138972204682732112496664293088698540380L,
     -8.870497589527808778965422264423713795614L},
};

inline std::vector<std::string> normal_order() {
    std::vector<std::string> order{water_id};
    for (const auto& component : hydrocarbons) { order.emplace_back(component.id); }
    return order;
}

inline std::vector<std::string> reversed_order() {
    auto order = normal_order();
    std::reverse(order.begin(), order.end());
    return order;
}

inline th::Sw92ParameterSet parameters(bool reverse = false) {
    const auto source = source_2017("Tables 4-5; Sample 6; non-water BIPs assumed zero");
    const auto water_source = source_2017(
        "Appendix A water-hydrocarbon NA BIPs; corrected SW92 AQ correlation is used by the model profile");

    std::vector<th::Component> catalog;
    th::Sw92ParameterInput input;
    input.model_id = std::string(th::sw92_corrected_profile);
    input.dataset_id = "MR2017-Sample6-SW92-corrected-original-PR76-base";
    input.revision = "doi-10.1016-j.fluid.2017.07.007+SW92-errata";
    input.applicability = {
        {std::nullopt, std::nullopt,
         source_2017("Sample-6 published phase envelope; selected anchor P=10 MPa, T=350 K")},
        std::nullopt};

    catalog.push_back({water_id, "Water", th::ComponentKind::pure,
                       source_2017("Table 4 H2O identity"), std::nullopt});
    input.pure.push_back({
        water_id, th::Sw92Species::water,
        scalar(647.30, th::Unit::kelvin, source, "647.30 K", "identity"),
        scalar(22.048e6, th::Unit::pascal, source, "22.048 MPa", "MPa * 1e6 -> Pa"),
        scalar(0.344, th::Unit::dimensionless, source)});

    for (const auto& component : hydrocarbons) {
        catalog.push_back({component.id, component.display, th::ComponentKind::pure,
                           source_2017("Table 4 synthetic-fluid component identity"),
                           std::nullopt});
        input.pure.push_back({
            component.id, th::Sw92Species::hydrocarbon,
            scalar(component.tc_k, th::Unit::kelvin, source, "K", "identity"),
            scalar(component.pc_mpa * 1.0e6, th::Unit::pascal, source,
                   "MPa", "MPa * 1e6 -> Pa"),
            scalar(component.omega, th::Unit::dimensionless, source)});

        th::Sw92WaterBinaryRecord water_pair;
        water_pair.component_id = component.id;
        water_pair.nonaqueous_rule = th::Sw92NonAqueousWaterRule::constant;
        water_pair.nonaqueous_kij = scalar(
            component.na_water_kij, th::Unit::dimensionless, water_source);
        input.water_binary.push_back(std::move(water_pair));
    }

    for (std::size_t i = 0; i < hydrocarbons.size(); ++i) {
        for (std::size_t j = i + 1; j < hydrocarbons.size(); ++j) {
            input.nonwater_binary.push_back({
                hydrocarbons[i].id, hydrocarbons[j].id,
                scalar(0.0, th::Unit::dimensionless, source,
                       "non-water BIPs assumed zero", "identity"),
                scalar(0.0, th::Unit::dimensionless, source,
                       "non-water BIPs assumed zero", "identity")});
        }
    }

    const auto order = reverse ? reversed_order() : normal_order();
    return th::Sw92ParameterSet::create(catalog, order, input);
}

inline th::Sw92Phase<double> model(bool reverse = false) {
    return th::Sw92Phase<double>::from_parameters(parameters(reverse));
}

inline Vec vector_from(const std::array<long double, 8>& values,
                       bool reverse = false) {
    Vec result;
    result.reserve(values.size());
    if (!reverse) {
        for (long double value : values) { result.push_back(static_cast<double>(value)); }
    } else {
        for (auto it = values.rbegin(); it != values.rend(); ++it) {
            result.push_back(static_cast<double>(*it));
        }
    }
    return result;
}

inline Vec feed(bool reverse = false) { return vector_from(golden.feed, reverse); }
inline Vec w(bool reverse = false) { return vector_from(golden.w, reverse); }
inline Vec h0(bool reverse = false) { return vector_from(golden.h0, reverse); }
inline Vec h1(bool reverse = false) { return vector_from(golden.h1, reverse); }

} // namespace sw92_profile_c_sample6

#endif // MPMC_TEST_SW92_PROFILE_C_PHYSICAL_SAMPLE6_HPP
