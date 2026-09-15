#ifndef MPMC_TEST_PR76_SOUR_GAS_FIXTURE_HPP
#define MPMC_TEST_PR76_SOUR_GAS_FIXTURE_HPP

#include <mpmc/flash/pr76_pt_flash_backend.hpp>
#include <mpmc/thermodynamics/pr76_phase.hpp>

#include "sour_gas_references.hpp"

#include <algorithm>
#include <array>
#include <optional>
#include <span>
#include <stdexcept>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

namespace pr76_sour_gas_test {
namespace fl = mpmc::flash;
namespace th = mpmc::thermodynamics;
namespace ref = pr76_sour_gas_reference;

using Vec = std::vector<double>;

inline constexpr std::array<std::string_view, 6> canonical_ids{
    "carbon-dioxide", "nitrogen", "hydrogen-sulfide",
    "methane", "ethane", "propane"};

inline constexpr std::array<double, 6> critical_temperatures{
    304.211, 126.2, 373.2, 190.564, 305.322, 369.825};
inline constexpr std::array<double, 6> critical_pressures_pa{
    73.819e5, 33.9e5, 89.4e5, 45.992e5, 48.718e5, 42.462e5};
inline constexpr std::array<double, 6> acentric_factors{
    0.225, 0.039, 0.081, 0.01141, 0.10574, 0.15813};

inline th::Provenance literature_source(std::string locator) {
    return {
        th::SourceKind::literature,
        "https://doi.org/10.2118/129844-PA",
        "Li and Firoozabadi, SPE Journal 17(4), 2012, 1096-1107",
        std::move(locator),
        "PR76 acid-gas three-phase model benchmark; literature-defined engineering state, not experimental validation",
        "Table values transcribed exactly; critical pressures bar * 100000 -> Pa",
        "Table 3 supplies component data/nonzero kij; Table 9 selects T=178.8 K, P=20 bar, n(CO2)=0.5"};
}

inline th::SourcedScalar scalar(
    double value, th::Unit unit, const th::Provenance& provenance,
    std::string source_unit = "SI or dimensionless",
    std::string conversion = "identity") {
    return {value, unit, provenance, std::move(source_unit), std::move(conversion)};
}

inline double kij(std::size_t first, std::size_t second) {
    if (first > second) { std::swap(first, second); }
    if (first == 0U && second == 1U) return -0.02;
    if (first == 0U && second == 2U) return 0.12;
    if (first == 1U && second == 2U) return 0.20;
    if (first == 0U && second == 3U) return 0.125;
    if (first == 1U && second == 3U) return 0.031;
    if (first == 2U && second == 3U) return 0.10;
    if (first == 0U && second == 4U) return 0.135;
    if (first == 1U && second == 4U) return 0.042;
    if (first == 2U && second == 4U) return 0.08;
    if (first == 0U && second == 5U) return 0.150;
    if (first == 1U && second == 5U) return 0.091;
    if (first == 2U && second == 5U) return 0.08;
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
    throw std::invalid_argument("unknown sour-gas component id");
}

inline th::Pr76Phase<double> model(const std::vector<std::string>& order) {
    if (order.size() != canonical_ids.size()) {
        throw std::invalid_argument("sour-gas benchmark requires six components");
    }
    const auto pure_source = literature_source(
        "Table 3: acid-gas Tc, Pc, acentric factor");
    const auto nonzero_source = literature_source(
        "Table 3: columns labelled nonzero binary interaction coefficients");
    const auto zero_source = literature_source(
        "Table 3 lists nonzero kij only; every unlisted selected pair is recorded explicitly as zero");

    std::vector<th::Component> catalog;
    catalog.reserve(canonical_ids.size());
    th::PrParameterInput input;
    input.model_id = std::string(th::pr76_profile);
    input.dataset_id = "Li-Firoozabadi-2012-acid-gas-PR76";
    input.revision = "table3-table9-reconstruction-v1";
    input.applicability = {std::nullopt, std::nullopt, pure_source};

    for (std::size_t i = 0; i < canonical_ids.size(); ++i) {
        const std::string id{canonical_ids[i]};
        catalog.push_back({
            id, id, th::ComponentKind::pure, pure_source, {}});
        input.pure.push_back({
            id,
            scalar(critical_temperatures[i], th::Unit::kelvin, pure_source),
            scalar(
                critical_pressures_pa[i], th::Unit::pascal, pure_source,
                "bar", "bar * 100000 -> Pa"),
            scalar(
                acentric_factors[i], th::Unit::dimensionless, pure_source)});
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

inline Vec feed(const std::vector<std::string>& order = canonical_order()) {
    return permute_from_canonical(ref::feed, order);
}

inline std::vector<Vec> starts(
    const std::vector<std::string>& order = canonical_order()) {
    std::vector<Vec> result;
    result.reserve(ref::phases.size());
    for (const auto& phase : ref::phases) {
        result.push_back(permute_from_canonical(phase, order));
    }
    return result;
}

inline fl::Pr76PtThreePhaseStart continuation_start(
    const std::vector<std::string>& order = canonical_order()) {
    const auto values = starts(order);
    fl::Pr76PtThreePhaseStart start;
    start.compositions = {values[0], values[1], values[2]};
    start.phase_fraction_seed = {
        ref::phase_fractions[1], ref::phase_fractions[2]};
    return start;
}

inline fl::Pr76PtMax3Options max3_options(
    const std::vector<std::string>& order = canonical_order()) {
    fl::Pr76PtMax3Options options;
    options.three_phase_starts.push_back(continuation_start(order));
    return options;
}

inline fl::Pr76PtFlashBackendOptions backend_options(
    const std::vector<std::string>& order = canonical_order()) {
    fl::Pr76PtFlashBackendOptions options;
    options.three_phase_starts.push_back(continuation_start(order));
    options.initial_starts = starts(order);
    options.final_starts = starts(order);
    return options;
}

inline Vec canonicalize(
    std::span<const double> composition,
    const std::vector<std::string>& order) {
    if (composition.size() != order.size()) {
        throw std::invalid_argument("sour-gas composition/order size mismatch");
    }
    Vec result(canonical_ids.size(), 0.0);
    for (std::size_t i = 0; i < order.size(); ++i) {
        result[canonical_index(order[i])] = composition[i];
    }
    return result;
}

} // namespace pr76_sour_gas_test

#endif // MPMC_TEST_PR76_SOUR_GAS_FIXTURE_HPP
