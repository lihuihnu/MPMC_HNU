#ifndef MPMC_TEST_SW92_PROFILE_C_PHYSICAL_SAMPLE6_HPP
#define MPMC_TEST_SW92_PROFILE_C_PHYSICAL_SAMPLE6_HPP

#include "test_support.hpp"

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

inline constexpr std::array<ComponentSpec, 7> hydrocarbons{{
    {"methane", "C1 methane", 190.60, 4.6042, 0.013, 0.4850},
    {"ethane", "C2 ethane", 305.43, 4.8839, 0.0986, 0.4920},
    {"propane", "C3 propane", 369.80, 4.2455, 0.1524, 0.5525},
    {"n-butane", "C4 n-butane", 419.50, 3.7470, 0.1956, 0.5091},
    {"n-pentane", "C5 n-pentane", 465.90, 3.3589, 0.2413, 0.5091},
    {"c6-pseudocomponent", "C6 pseudocomponent", 507.50, 3.0104, 0.2990, 0.5091},
    {"c7plus-pseudocomponent", "C7+ pseudocomponent", 655.04, 2.2305, 0.50879, 0.5091},
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
// SW92/corrected-original/PR76-base profile. H0/H1 are in production
// canonical order (smaller Z first). The literature source is the topology and
// parameter oracle; these numerical values are not copied from its plots.
inline constexpr PhysicalGolden golden{
    {0.5085508550855085508550855086L,
     0.4249424942494249424942494249L,
     0.02140214021402140214021402140L,
     0.01120112011201120112011201120L,
     0.008500850085008500850085008501L,
     0.004000400040004000400040004000L,
     0.002900290029002900290029002900L,
     0.01850185018501850185018501850L},
    {0.9989848247418118913595543401053648926216L,
     0.0009437390246448547918378251367314694930708L,
     0.00005090131750566792055103740799881763691707L,
     0.00001539112068535404654808845623374559746129L,
     0.000004502294921361228411799492988079242442326L,
     0.0000005926374988163918850182978502110622052845L,
     0.00000004884305598551627085540465878795198094257L,
     0.00000000001987606874494103569917399639430023336960L},
    {0.003431159072736137368866463407210169130387L,
     0.3047879878876917839947681909291417617457L,
     0.03967884978986410404078112552403412438011L,
     0.03907378318427884495687733496322471185180L,
     0.04972553578877213899081956960285877255227L,
     0.03753384147055697177383046571975924007300L,
     0.04023696182464144267408316726422341125872L,
     0.4855318809814585761999736825905478090080L},
    {0.006369077653896316245955656486613138123090L,
     0.9035132844657757570980508964104985421547L,
     0.04356105211026104596721953258170528172610L,
     0.02135051898098504489224847441036670205068L,
     0.01461052383168891829831565291439887212351L,
     0.005749900978503531710118972141797761945924L,
     0.003130201254861312400440558090923953166595L,
     0.001715440724028073387650256964695748709383L},
    {0.5060256031413864067026662808012660533604L,
     0.03649000924517693070739764583343665718085L,
     0.4574843876134366625899360733652972894588L},
    {0.07611469089968388060693842733302756249187L,
     0.4936978040184203266617168399363985882571L,
     0.8705435515569661648431820968661642550187L},
    {-5.272115249658086899493622706522535376994L,
     -0.2050073022742263683691499361993708681149L,
     -3.542192391926915690432347179288204344059L,
     -4.498176782098288675017853347850338437173L,
     -5.104151976819097777349954263821877945160L,
     -6.275905166135803522868553987975552454083L,
     -7.139077923509930933399810949770367699721L,
     -8.870267557929932874419276542434634442359L},
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
    const auto water_source = source_sw92("Table 5 non-aqueous water BIPs; corrected Eqs. (11)-(12) for AQ");

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
