#ifndef MPMC_TEST_PR76_HERINGER_2026_SOUR_GAS_FIXTURE_HPP
#define MPMC_TEST_PR76_HERINGER_2026_SOUR_GAS_FIXTURE_HPP

#include <mpmc/thermodynamics/pr76_phase.hpp>

#include <algorithm>
#include <array>
#include <cstddef>
#include <optional>
#include <stdexcept>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

namespace pr76_heringer_2026_sour_gas_test {
namespace th = mpmc::thermodynamics;

using Vec = std::vector<double>;

// Heringer et al., Fluid Phase Equilibria 604 (2026) 114653,
// Appendix B, Table B1. C1/C2/C3 are mapped to methane/ethane/propane.
inline constexpr std::array<std::string_view, 6> canonical_ids{
    "carbon-dioxide", "nitrogen", "hydrogen-sulfide",
    "methane", "ethane", "propane"};

inline constexpr std::array<double, 6> critical_temperatures_k{
    304.20, 126.20, 373.20, 190.60, 305.40, 369.80};
inline constexpr std::array<double, 6> critical_pressures_pa{
    73.90e5, 33.50e5, 89.40e5, 45.40e5, 48.20e5, 41.90e5};
inline constexpr std::array<double, 6> acentric_factors{
    0.225, 0.040, 0.081, 0.008, 0.098, 0.152};

// Preserve the published rounded Table B1 values exactly. They sum to 99.98 mol%
// because of source-table rounding; table_b1_feed() performs the explicit normalization
// required by the repository composition contract.
inline constexpr std::array<double, 6> table_b1_mole_percent{
    70.59, 7.03, 1.97, 6.86, 10.56, 2.97};

inline th::Provenance table_b1_source(std::string locator) {
    return {
        th::SourceKind::literature,
        "https://doi.org/10.1016/j.fluid.2025.114653",
        "Heringer et al., Fluid Phase Equilibria 604 (2026) 114653",
        std::move(locator),
        "PR sour-gas benchmark input from Appendix B Table B1; not experimental validation",
        "Appendix B Table B1 transcribed from the open-access article; Pc converted from bar to Pa",
        "CC BY 4.0 article; fixture stores cited tabulated numeric inputs and provenance metadata"};
}

inline th::SourcedScalar scalar(
    double value, th::Unit unit, const th::Provenance& provenance,
    std::string source_unit = "SI or dimensionless",
    std::string conversion = "identity") {
    return {value, unit, provenance, std::move(source_unit), std::move(conversion)};
}

// Table B1 is explicitly a table of non-zero BIPs. Every unlisted selected pair is
// therefore recorded as an explicit zero in the complete PR76 parameter snapshot.
inline double kij(std::size_t first, std::size_t second) {
    if (first > second) { std::swap(first, second); }
    if (first == 0U && second == 1U) return -0.020;
    if (first == 0U && second == 2U) return 0.120;
    if (first == 1U && second == 2U) return 0.200;
    if (first == 0U && second == 3U) return 0.125;
    if (first == 1U && second == 3U) return 0.031;
    if (first == 2U && second == 3U) return 0.100;
    if (first == 0U && second == 4U) return 0.135;
    if (first == 1U && second == 4U) return 0.042;
    if (first == 2U && second == 4U) return 0.080;
    if (first == 0U && second == 5U) return 0.150;
    if (first == 1U && second == 5U) return 0.091;
    if (first == 2U && second == 5U) return 0.080;
    return 0.0;
}

inline std::vector<std::string> canonical_order() {
    std::vector<std::string> order;
    order.reserve(canonical_ids.size());
    for (const auto id : canonical_ids) { order.emplace_back(id); }
    return order;
}

inline std::vector<std::string> reversed_order() {
    auto order = canonical_order();
    std::reverse(order.begin(), order.end());
    return order;
}

inline std::size_t canonical_index(std::string_view id) {
    for (std::size_t i = 0; i < canonical_ids.size(); ++i) {
        if (canonical_ids[i] == id) { return i; }
    }
    throw std::invalid_argument("unknown Heringer-2026 sour-gas component id");
}

template <typename Array>
inline Vec permute_from_canonical(
    const Array& values, const std::vector<std::string>& order) {
    Vec result;
    result.reserve(order.size());
    for (const auto& id : order) {
        result.push_back(values[canonical_index(id)]);
    }
    return result;
}

inline Vec table_b1_feed(
    const std::vector<std::string>& order = canonical_order()) {
    double total = 0.0;
    for (const double value : table_b1_mole_percent) { total += value; }
    if (!(total > 0.0)) {
        throw std::runtime_error("Heringer-2026 Table B1 composition has no positive support");
    }
    std::array<double, canonical_ids.size()> normalized{};
    for (std::size_t i = 0; i < normalized.size(); ++i) {
        normalized[i] = table_b1_mole_percent[i] / total;
    }
    return permute_from_canonical(normalized, order);
}

inline th::Pr76Phase<double> model(const std::vector<std::string>& order) {
    if (order.size() != canonical_ids.size()) {
        throw std::invalid_argument("Heringer-2026 sour-gas fixture requires six components");
    }

    const auto pure_source = table_b1_source(
        "Appendix B, Table B1: Tc, Pc, omega and component identity");
    const auto nonzero_source = table_b1_source(
        "Appendix B, Table B1: columns Ki,CO2 / Ki,N2 / Ki,H2S (non-zero BIPs)");
    const auto zero_source = table_b1_source(
        "Appendix B states non-zero BIPs; unlisted selected pairs are recorded explicitly as zero");

    std::vector<th::Component> catalog;
    catalog.reserve(canonical_ids.size());
    th::PrParameterInput input;
    input.model_id = std::string(th::pr76_profile);
    input.dataset_id = "Heringer-et-al-2026-sour-gas-PR76";
    input.revision = "appendix-b-table-b1-v1";
    input.applicability = {std::nullopt, std::nullopt, pure_source};

    for (std::size_t i = 0; i < canonical_ids.size(); ++i) {
        const std::string id{canonical_ids[i]};
        catalog.push_back({id, id, th::ComponentKind::pure, pure_source, {}});
        input.pure.push_back({
            id,
            scalar(critical_temperatures_k[i], th::Unit::kelvin, pure_source),
            scalar(critical_pressures_pa[i], th::Unit::pascal, pure_source,
                   "bar", "bar * 100000 -> Pa"),
            scalar(acentric_factors[i], th::Unit::dimensionless, pure_source)});
    }

    for (std::size_t i = 0; i < canonical_ids.size(); ++i) {
        for (std::size_t j = i + 1U; j < canonical_ids.size(); ++j) {
            const double value = kij(i, j);
            const auto& source = value == 0.0 ? zero_source : nonzero_source;
            input.binary.push_back({
                std::string(canonical_ids[i]),
                std::string(canonical_ids[j]),
                scalar(value, th::Unit::dimensionless, source)});
        }
    }

    return th::Pr76Phase<double>::from_parameters(
        th::PrParameterSet::create(catalog, order, input));
}

inline th::Pr76Phase<double> model() {
    return model(canonical_order());
}

} // namespace pr76_heringer_2026_sour_gas_test

#endif // MPMC_TEST_PR76_HERINGER_2026_SOUR_GAS_FIXTURE_HPP
