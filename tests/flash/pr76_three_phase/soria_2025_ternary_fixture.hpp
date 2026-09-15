#ifndef MPMC_TEST_PR76_SORIA_2025_TERNARY_FIXTURE_HPP
#define MPMC_TEST_PR76_SORIA_2025_TERNARY_FIXTURE_HPP

#include <mpmc/thermodynamics/pr76_phase.hpp>

#include <array>
#include <cstddef>
#include <optional>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

namespace pr76_soria_2025_ternary_test {
namespace th = mpmc::thermodynamics;

inline constexpr std::array<std::string_view, 3> component_ids{
    "carbon-dioxide", "n-decane", "n-hexadecane"};
inline constexpr std::array<double, 3> critical_temperatures_k{
    304.21, 617.7, 717.0};
inline constexpr std::array<double, 3> critical_pressures_pa{
    7.383e6, 2.110e6, 1.489e6};
inline constexpr std::array<double, 3> acentric_factors{
    0.2236, 0.4923, 0.742};

// Soria (UFC MSc dissertation, 2025), Table 12, T=323.15 K.
inline constexpr double co2_decane_kij = 0.09670;
inline constexpr double co2_hexadecane_kij = 0.10127;
inline constexpr double decane_hexadecane_kij = 0.0;

// Table 9, CO2 + M-40, nominal x_CO2=0.9020 at 323.15 K.
// Table 4 defines M-40 as 0.600 n-C10 + 0.400 n-C16.
inline constexpr std::array<double, 3> m40_feed{
    0.9020, 0.0588, 0.0392};
inline constexpr double experimental_temperature_k = 323.15;
inline constexpr double experimental_transition_pressure_pa = 11.21e6;
inline constexpr double expanded_pressure_uncertainty_pa = 0.07e6;
inline constexpr double co2_composition_uncertainty = 0.0006;

// Table 13 reports 15.39% AARD for the source's own PR modeling of M-40 at
// 323.15 K. This is retained only as a predeclared comparison scale. It is NOT
// an experimental uncertainty, a pointwise error bound, or proof that the
// source's PR alpha convention is identical to this repository's strict PR76.
inline constexpr double source_pr_aard_fraction = 0.1539;

[[nodiscard]] inline double lower_comparison_pressure_pa() noexcept {
    return experimental_transition_pressure_pa * (1.0 - source_pr_aard_fraction) -
           expanded_pressure_uncertainty_pa;
}

[[nodiscard]] inline double upper_comparison_pressure_pa() noexcept {
    return experimental_transition_pressure_pa * (1.0 + source_pr_aard_fraction) +
           expanded_pressure_uncertainty_pa;
}

inline th::Provenance source(std::string locator, std::string note = {}) {
    return {
        th::SourceKind::literature,
        "https://repositorio.ufc.br/handle/riufc/80725",
        "E. C. Q. Soria, UFC MSc dissertation, approved 2025-01-06",
        std::move(locator),
        std::move(note),
        "Transcribed from the uploaded/open-access UFC dissertation",
        "UFC repository record is open access; this test stores cited numeric facts and provenance only"};
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
    const auto pure_source = source(
        "Table 6: pure-component Pc, Tc and acentric factor",
        "The repository evaluates these tabulated data with strict PR76 alpha; the dissertation prints an alternate high-omega branch for omega >= 0.49, so this fixture is an experimental transfer benchmark rather than an identical-model reproduction.");
    const auto binary_source = source(
        "Table 12: binary interaction coefficients at T=323.15 K",
        "Published kij values were fitted in the dissertation's PR implementation and are transferred here without refitting to strict PR76.");
    const auto zero_source = source(
        "Sections 4.4.2.2 and Table 12: n-C10/n-C16 kij adopted as zero");
    const auto applicability_source = source(
        "Tables 6, 9 and 12: this fixture is restricted to the 323.15 K experimental gate",
        "No pressure validity interval is inferred from one transition datum.");

    std::vector<th::Component> catalog;
    catalog.reserve(component_ids.size());
    for (const auto id_view : component_ids) {
        const std::string id{id_view};
        catalog.push_back({id, id, th::ComponentKind::pure, pure_source, {}});
    }

    th::PrParameterInput input;
    input.model_id = std::string(th::pr76_profile);
    input.dataset_id = "Soria-UFC-2025-CO2-C10-C16-PR76-transfer";
    input.revision = "tables-6-9-12-13-m40-323K-v1";
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
        scalar(co2_decane_kij, th::Unit::dimensionless, binary_source)});
    input.binary.push_back({
        std::string(component_ids[0]), std::string(component_ids[2]),
        scalar(co2_hexadecane_kij, th::Unit::dimensionless, binary_source)});
    input.binary.push_back({
        std::string(component_ids[1]), std::string(component_ids[2]),
        scalar(decane_hexadecane_kij, th::Unit::dimensionless, zero_source)});

    const auto component_order = order();
    return th::Pr76Phase<double>::from_parameters(
        th::PrParameterSet::create(catalog, component_order, input));
}

} // namespace pr76_soria_2025_ternary_test

#endif // MPMC_TEST_PR76_SORIA_2025_TERNARY_FIXTURE_HPP
