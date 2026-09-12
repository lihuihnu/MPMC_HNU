#include <mpmc/flash/cpa_pt_flash_backend.hpp>
#include <mpmc/flash/pr76_pt_flash_backend.hpp>
#include <mpmc/flash/sw92_profile_c_pt_flash_backend.hpp>
#include <mpmc/runtime/pt_service.hpp>

#include "../../flash/cpa_max3/test_support.hpp"
#include "../../flash/pr76_three_phase/synthetic_fixture.hpp"
#include "test_support.hpp"

#include <array>
#include <stdexcept>
#include <string_view>

namespace {
namespace fl = mpmc::flash;
namespace rt = mpmc::runtime;
namespace th = mpmc::thermodynamics;

void require(bool condition, const char* message) {
    if (!condition) { throw std::runtime_error(message); }
}

void require_descriptor(const rt::PtService& service,
                        std::string_view configured_backend_id,
                        const fl::PtFlashBackendCapability& capability) {
    const auto* descriptor = service.find_capability(configured_backend_id);
    require(descriptor != nullptr && descriptor->structurally_valid(),
            "configured backend missing from service discovery");
    require(descriptor->configured_backend_id == configured_backend_id &&
                descriptor->capability.component_ids == capability.component_ids &&
                descriptor->capability.supported_phase_counts ==
                    std::vector<std::size_t>({1U, 2U, 3U}) &&
                descriptor->component_inventory.components.size() ==
                    capability.component_ids.size(),
            "service discovery changed backend capability or inventory");
    for (std::size_t i = 0; i < capability.component_ids.size(); ++i) {
        require(descriptor->component_inventory.components[i].component_id ==
                        capability.component_ids[i] &&
                    descriptor->component_inventory.components[i].feed_index == i,
                "service component inventory lost backend order");
    }
}

} // namespace

int main() {
    try {
        const auto pr_model = pr76_max3_test::model();
        fl::Pr76VleEvaluator pr_evaluator(pr_model);
        fl::Pr76PtFlashBackend pr_backend(pr_evaluator);

        const auto cpa_parameters = cpa_max3_test::parameters();
        const auto cpa_model = th::CpaPtPhase::from_parameters(cpa_parameters);
        fl::CpaVleEvaluator cpa_evaluator(
            cpa_model, cpa_max3_test::fast_pt_options());
        fl::CpaPtFlashBackend cpa_backend(cpa_evaluator);

        const auto sw_model = th::Sw92Phase<double>::from_parameters(
            sw92_test::binary_parameters(sw92_test::co2));
        fl::Sw92ProfileCPtFlashBackend sw_backend(sw_model);

        const std::array<rt::PtServiceBackendRegistration, 3> backends{{
            {"pr76.fixture", &pr_backend},
            {"sw92.fixture", &sw_backend},
            {"cpa.fixture", &cpa_backend}}};
        const rt::PtService service(backends);
        require(service.discover_capabilities().size() == 3U,
                "three-backend service discovery count changed");
        require_descriptor(service, "pr76.fixture", pr_backend.capability());
        require_descriptor(service, "sw92.fixture", sw_backend.capability());
        require_descriptor(service, "cpa.fixture", cpa_backend.capability());

        require(service.find_capability("pr76.fixture") != nullptr &&
                    service.find_capability("sw92.fixture") != nullptr &&
                    service.find_capability("cpa.fixture") != nullptr,
                "configured backend IDs are not selectable through PT service v1");
        return 0;
    } catch (const std::exception&) {
        return 1;
    }
}
