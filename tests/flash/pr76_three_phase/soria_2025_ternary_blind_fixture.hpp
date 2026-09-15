#ifndef MPMC_TEST_PR76_SORIA_2025_TERNARY_BLIND_FIXTURE_HPP
#define MPMC_TEST_PR76_SORIA_2025_TERNARY_BLIND_FIXTURE_HPP

#include <mpmc/thermodynamics/pr76_phase.hpp>

#include <array>
#include <cstddef>
#include <optional>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

namespace pr76_soria_2025_ternary_blind_test {
namespace th = mpmc::thermodynamics;

inline constexpr std::array<std::string_view, 3> component_ids{
    "carbon-dioxide", "n-decane", "n-hexadecane"};
inline constexpr std::array<double, 3> critical_temperatures_k{
    304.21, 617.7, 717.0};
inline constexpr std::array<double, 3> critical_pressures_pa{
    7.383e6, 2.110e6, 1.489e6};
inline constexpr std::array<double, 3> acentric_factors{
    0.2236, 0.4923, 0.742};

// Frozen strict-PR76 binary calibration from MPMC_HNU PR #101. These values
// were fitted before this ternary validation and the M-40 11.21 MPa datum was
// explicitly excluded from every objective, bound, seed and assertion.
inline constexpr double co2_decane_kij = 0.05226578047;
inline constexpr double co2_hexadecane_kij = 0.07884923875;
inline constexpr double decane_hexadecane_kij = 0.0;

// Soria UFC 2025 Tables 4 and 9: M-40 = 60 mol% n-C10 + 40 mol% n-C16;
// nominal overall composition at the 323.15 K L -> L+V transition.
inline constexpr std::array<double, 3> m40_feed{
    0.9020, 0.0588, 0.0392};
inline constexpr double experimental_temperature_k = 323.15;
inline constexpr double experimental_transition_pressure_pa = 11.21e6;
inline constexpr double expanded_pressure_uncertainty_pa = 0.07e6;
inline constexpr double co2_composition_uncertainty = 0.0006;

inline th::Provenance ufc_source(std::string locator, std::string note = {}) {
    return {
        th::SourceKind::literature,
        "https://repositorio.ufc.br/handle/riufc/80725",
        "E. C. Q. Soria, UFC dissertation, 2025",
        std::move(locator),
        std::move(note),
        "Transcribed from the UFC repository dissertation",
        "Test-only numeric citation and provenance metadata"};
}

inline th::Provenance calibration_source(std::string locator) {
    return {
        th::SourceKind::assumption,
        "MPMC_HNU PR #101 / merge 57a573a676832a6aee78ef9894204702638da132",
        "strict-PR76 independent binary bubble-pressure calibration",
        std::move(locator),
        "Derived only from independent CO2+n-C10 / CO2+n-C16 binary bubble data; frozen before the M-40 ternary blind validation and never fitted to its 11.21 MPa transition",
        "Computed by the hosted-verified test merged in PR #101",
        "Repository-derived test parameter; not a literature-tabulated kij"};
}

inline th::SourcedScalar scalar(
    double value, th::Unit unit, const th::Provenance& provenance,
    std::string original_unit = "SI or dimensionless",
    std::string conversion = "identity") {
    return {value, unit, provenance, std::move(original_unit), std::move(conversion)};
}

inline std::vector<std::string> order() {
    std::vector<std::string> result;
    result.reserve(component_ids.size());
    for (const auto id : component_ids) { result.emplace_back(id); }
    return result;
}

inline th::Pr76Phase<double> model() {
    const auto pure_source = ufc_source(
        "Table 6: pure-component Tc, Pc and acentric factor",
        "Evaluated with repository strict PR76 printed-coefficient quadratic kappa for every omega");
    const auto c10_source = calibration_source(
        "PR #101 CO2+n-C10 local ~323 K fit: kij=0.05226578047");
    const auto c16_source = calibration_source(
        "PR #101 CO2+n-C16 fits at 313.2/333.2 K, linearly interpolated to 323.15 K: kij=0.07884923875");
    const auto zero_source = ufc_source(
        "Sections 4.4.2.2 / Table 12: n-C10+n-C16 kij adopted as zero",
        "The independent strict-PR76 recalibration targeted only the two CO2-heavy binary pairs required by the blind-validation plan");
    const auto applicability_source = ufc_source(
        "Table 9: M-40 transition temperature 323.15 K",
        "No pressure validity interval is inferred from one experimental transition");

    std::vector<th::Component> catalog;
    catalog.reserve(component_ids.size());
    for (const auto id_view : component_ids) {
        const std::string id{id_view};
        catalog.push_back({id, id, th::ComponentKind::pure, pure_source, {}});
    }

    th::PrParameterInput input;
    input.model_id = std::string(th::pr76_profile);
    input.dataset_id = "Soria-UFC-2025-M40-strict-PR76-blind";
    input.revision = "table6-table9/pr101-binary-kij/frozen-v1";
    input.applicability = {
        th::ClosedInterval{experimental_temperature_k, experimental_temperature_k},
        std::nullopt,
        applicability_source};

    for (std::size_t i = 0U; i < component_ids.size(); ++i) {
        const std::string id{component_ids[i]};
        input.pure.push_back({
            id,
            scalar(critical_temperatures_k[i], th::Unit::kelvin, pure_source, "K", "identity"),
            scalar(critical_pressures_pa[i], th::Unit::pascal, pure_source,
                   "MPa", "MPa * 1e6 -> Pa"),
            scalar(acentric_factors[i], th::Unit::dimensionless, pure_source,
                   "dimensionless", "identity")});
    }

    input.binary.push_back({
        std::string(component_ids[0]), std::string(component_ids[1]),
        scalar(co2_decane_kij, th::Unit::dimensionless, c10_source)});
    input.binary.push_back({
        std::string(component_ids[0]), std::string(component_ids[2]),
        scalar(co2_hexadecane_kij, th::Unit::dimensionless, c16_source)});
    input.binary.push_back({
        std::string(component_ids[1]), std::string(component_ids[2]),
        scalar(decane_hexadecane_kij, th::Unit::dimensionless, zero_source)});

    const auto component_order = order();
    return th::Pr76Phase<double>::from_parameters(
        th::PrParameterSet::create(catalog, component_order, input));
}

} // namespace pr76_soria_2025_ternary_blind_test

#endif // MPMC_TEST_PR76_SORIA_2025_TERNARY_BLIND_FIXTURE_HPP
