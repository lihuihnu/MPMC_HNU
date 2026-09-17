#ifndef MPMC_TEST_CPA_THERMOPACK_PARAMETER_SNAPSHOT_HPP
#define MPMC_TEST_CPA_THERMOPACK_PARAMETER_SNAPSHOT_HPP

#include "cpa_thermopack_phase_kernel_generated.hpp"
#include "test_support.hpp"

#include <stdexcept>
#include <string>
#include <string_view>
#include <vector>

namespace cpa_thermopack_snapshot {
namespace th = mpmc::thermodynamics;
namespace external = cpa_thermopack_phase_kernel;

inline constexpr std::string_view thermopack_commit =
    "d68c794c7342bfc6938eb424a1fbb88b7780b738";

// This is deliberately NOT named a ThermoPack-default dataset. The pinned
// external oracle explicitly injects the repository's literature a0/b/c1,
// association parameters and Folas kij=-0.055 into ThermoPack, while Classic
// alpha continues to consume the component Tc from the pinned ThermoPack
// component database. The snapshot identity records exactly that hybrid origin.
inline constexpr std::string_view parity_dataset_id =
    "thermopack-d68c794__meoh-h2o__mpmc-literature-folas-cr1-aligned";
inline constexpr std::string_view parity_revision =
    "v1__Tc-thermopack-component-db__other-cpa-fields-preserved-from-literature-fixture";

inline th::Provenance thermopack_alpha_tc_source(std::string locator) {
    return {
        th::SourceKind::database,
        "thermotools/thermopack pinned component database",
        std::string(thermopack_commit),
        std::move(locator),
        "Parity-snapshot field only: Classic-alpha Tc comes from the pinned ThermoPack component database; this does not replace the repository literature dataset.",
        "Read from the pinned ThermoPack component JSON used to initialize cbeos%single(i)%tc.",
        "ThermoPack source metadata used only for reproducible numerical parity validation."};
}

[[nodiscard]] inline th::CpaParameterSet parameters(bool swapped = false) {
    const auto baseline = cpa_physical_test::parameters(swapped);

    th::CpaParameterInput input;
    input.dataset_id = std::string(parity_dataset_id);
    input.revision = std::string(parity_revision);
    input.applicability = baseline.applicability();
    input.pure.assign(baseline.pure_records().begin(), baseline.pure_records().end());
    input.binary.assign(baseline.binary_records().begin(), baseline.binary_records().end());
    input.association_pairs.assign(
        baseline.association_records().begin(), baseline.association_records().end());

    // Only the two Classic-alpha Tc fields differ from the literature snapshot.
    // a0/b/c1, sites, epsilon/beta, cross-association records and kij retain the
    // pre-existing literature/Folas provenance verbatim.
    for (auto& pure : input.pure) {
        if (pure.component_id == "METHANOL") {
            pure.critical_temperature_k.value =
                external::methanol_alpha_critical_temperature_k;
            pure.critical_temperature_k.source = thermopack_alpha_tc_source(
                "fluids/Methanol.json: critical.temperature = 512.6 K");
        } else if (pure.component_id == "WATER") {
            pure.critical_temperature_k.value =
                external::water_alpha_critical_temperature_k;
            pure.critical_temperature_k.source = thermopack_alpha_tc_source(
                "fluids/Water.json: critical.temperature = 647.3 K");
        } else {
            throw std::runtime_error(
                "unexpected component in methanol-water ThermoPack parity snapshot");
        }
    }

    std::vector<std::string> order;
    order.reserve(baseline.components().size());
    for (const auto& component : baseline.components().items()) {
        order.push_back(component.id);
    }
    return th::CpaParameterSet::create(
        baseline.components().items(), order, input);
}

} // namespace cpa_thermopack_snapshot

#endif // MPMC_TEST_CPA_THERMOPACK_PARAMETER_SNAPSHOT_HPP
