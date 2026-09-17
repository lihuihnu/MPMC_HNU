#include <mpmc/flash/pr76_pt_flash_backend.hpp>
#include <mpmc/flash/sw92_profile_c_pt_flash_backend.hpp>
#include <mpmc/runtime/pt_service.hpp>

#include "../../thermodynamics/sw92/test_support.hpp"

#include <array>
#include <optional>
#include <stdexcept>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

namespace {
namespace fl = mpmc::flash;
namespace rt = mpmc::runtime;
namespace th = mpmc::thermodynamics;

void require(bool condition, std::string_view message) {
    if (!condition) { throw std::runtime_error(std::string(message)); }
}

th::Provenance synthetic_pr_binary_source() {
    return {th::SourceKind::synthetic_test,
            "test://runtime/pt-service/pr-sw-interchangeability",
            "v1",
            "Classic PR binary interaction fixture",
            "Zero binary interaction is used only for interface conformance; "
            "no physical-accuracy claim is made",
            "tests/runtime/pt_service/pr_sw_interchangeability_test.cpp",
            "Synthetic test datum; not experimental or fitted data"};
}

th::Pr76Phase<double> classic_pr_model_with_sw_inventory() {
    const auto pure_source = sw92_test::paper("Table 3 physical properties");
    const auto binary_source = synthetic_pr_binary_source();
    const std::array<const sw92_test::Spec*, 2> specs{
        &sw92_test::co2, &sw92_test::water};

    std::vector<th::Component> catalog;
    std::vector<std::string> order;
    th::PrParameterInput input;
    input.model_id = std::string(th::pr76_profile);
    input.dataset_id = "pt-service-pr-sw-interface-fixture";
    input.revision = "v1";
    input.applicability = {
        std::nullopt, std::nullopt,
        binary_source};

    for (const auto* spec : specs) {
        catalog.push_back({spec->id, spec->display, th::ComponentKind::pure,
                           pure_source, std::nullopt});
        order.emplace_back(spec->id);
        input.pure.push_back({
            spec->id,
            sw92_test::scalar(spec->tc, th::Unit::kelvin, pure_source,
                              "K", "identity"),
            sw92_test::scalar(spec->pc_bar * 100000.0, th::Unit::pascal,
                              pure_source, "bar", "bar * 100000 -> Pa"),
            sw92_test::scalar(spec->omega, th::Unit::dimensionless,
                              pure_source)});
    }
    input.binary.push_back({
        sw92_test::co2.id, sw92_test::water.id,
        sw92_test::scalar(0.0, th::Unit::dimensionless, binary_source,
                          "dimensionless", "identity")});

    const auto parameters = th::PrParameterSet::create(
        catalog, order, input, th::DataPolicy::allow_synthetic_tests);
    return th::Pr76Phase<double>::from_parameters(parameters);
}

void require_role_neutral_response(
    const rt::PtServiceResponse& response,
    std::string_view configured_backend_id) {
    require(response.structurally_valid(),
            "service returned a structurally invalid response");
    require(response.outcome == rt::PtServiceOutcome::accepted,
            "shared PR/SW conformance state must be accepted");
    require(response.accepted_phase_count() == 1U,
            "shared pure-CO2 conformance state must remain one phase");
    require(response.result.has_value(),
            "accepted response lost its result envelope");

    const auto& result = *response.result;
    require(result.provenance.backend.configured_backend_id ==
                configured_backend_id,
            "service did not preserve the selected configured backend ID");
    require(result.provenance.backend.capability.phase_metadata_namespace.empty(),
            "PR/SW external contract must not require provider phase identities");
    require(!result.morphology_resolved,
            "PR/SW external contract must not publish phase morphology");
    require(result.pressure_pa == 3.0e6 && result.temperature_k == 340.0,
            "service changed the common PT request");
    require(result.feed.size() == 2U &&
                result.feed[0].component_id == sw92_test::co2.id &&
                result.feed[0].mole_fraction == 1.0 &&
                result.feed[1].component_id == sw92_test::water.id &&
                result.feed[1].mole_fraction == 0.0,
            "service did not map the common component-ID feed to backend order");

    for (const auto& phase : result.phases) {
        require(!phase.provider_metadata.has_value(),
                "PR/SW public phase unexpectedly carries EOS-specific identity metadata");
        require(phase.components.size() == 2U &&
                    phase.components[0].component_id == sw92_test::co2.id &&
                    phase.components[1].component_id == sw92_test::water.id,
                "PR/SW public phase changed the shared component envelope");
    }
}

} // namespace

int main() {
    try {
        const auto pr_model = classic_pr_model_with_sw_inventory();
        fl::Pr76VleEvaluator pr_evaluator(pr_model);
        fl::Pr76PtFlashBackend pr_backend(pr_evaluator);

        const auto sw_model = th::Sw92Phase<double>::from_parameters(
            sw92_test::binary_parameters(sw92_test::co2));
        fl::Sw92ProfileCPtFlashBackend sw_backend(sw_model);

        require(pr_backend.capability().component_ids ==
                    sw_backend.capability().component_ids,
                "fixture no longer gives PR and SW the same component inventory");
        require(pr_backend.capability().phase_metadata_namespace.empty() &&
                    sw_backend.capability().phase_metadata_namespace.empty(),
                "PR and SW no longer share the role-neutral external phase contract");

        const std::array<rt::PtServiceBackendRegistration, 2> backends{{
            {"classic-pr", &pr_backend},
            {"sw92", &sw_backend}}};
        rt::PtService service(backends);

        // Deliberately submit the feed in the reverse of backend order. The exact
        // same public request payload is reused for both EOS backends; only the
        // configured_backend_id selector changes.
        rt::PtServiceRequest request;
        request.pressure_pa = 3.0e6;
        request.temperature_k = 340.0;
        request.feed = {{sw92_test::water.id, 0.0},
                        {sw92_test::co2.id, 1.0}};

        request.configured_backend_id = "classic-pr";
        const auto pr_response = service.solve(request);
        require_role_neutral_response(pr_response, "classic-pr");

        request.configured_backend_id = "sw92";
        const auto sw_response = service.solve(request);
        require_role_neutral_response(sw_response, "sw92");

        require(pr_response.result->provenance.backend.component_inventory.components ==
                    sw_response.result->provenance.backend.component_inventory.components,
                "switching EOS backend changed the public component inventory");
        return 0;
    } catch (const std::exception&) {
        return 1;
    }
}
