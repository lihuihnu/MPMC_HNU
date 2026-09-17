#include <mpmc/flash/cpa_pt_flash_backend.hpp>
#include <mpmc/flash/pr76_pt_flash_backend.hpp>
#include <mpmc/flash/sw92_profile_c_pt_flash_backend.hpp>

#include "../cpa_max3/test_support.hpp"
#include "../pr76_three_phase/synthetic_fixture.hpp"
#include "test_support.hpp"

#include <stdexcept>
#include <string_view>

namespace {
namespace fl = mpmc::flash;
namespace th = mpmc::thermodynamics;

void require(bool condition, const char* message) {
    if (!condition) { throw std::runtime_error(message); }
}

void require_common_max3_contract(const fl::PtFlashBackendCapability& capability) {
    require(capability.structurally_valid(), "backend capability is structurally invalid");
    require(capability.supported_phase_counts ==
                std::vector<std::size_t>({1U, 2U, 3U}) &&
                capability.maximum_phase_count() == 3U,
            "backend phase-count capability changed");
    require(capability.performs_initial_stability_search &&
                capability.performs_final_phase_set_review &&
                capability.performs_boundary_neighbor_resolve &&
                !capability.global_stability_proven,
            "backend common search/review capability changed");

    const auto* one_to_two = capability.transition_capability.edge(1U, 2U);
    const auto* two_to_one = capability.transition_capability.edge(2U, 1U);
    const auto* two_to_three = capability.transition_capability.edge(2U, 3U);
    const auto* three_to_two = capability.transition_capability.edge(3U, 2U);
    require(one_to_two != nullptr && two_to_one != nullptr &&
                two_to_three != nullptr && three_to_two != nullptr,
            "backend lost a frozen transition edge");
    require(one_to_two->support == fl::PtPhaseTransitionSupport::fresh_target_resolve &&
                two_to_one->support == fl::PtPhaseTransitionSupport::detection_only &&
                two_to_three->support == fl::PtPhaseTransitionSupport::fresh_target_resolve &&
                three_to_two->support == fl::PtPhaseTransitionSupport::fresh_target_resolve,
            "backend transition-support semantics changed");
}

} // namespace

int main() {
    try {
        static_assert(fl::pt_flash_backend_convention ==
                      std::string_view{"PT/flash-backend/capability-and-dispatch/v1"});
        static_assert(fl::PtFlashBackendResult::convention ==
                      std::string_view{"PT/flash-backend/result/v1"});

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

        require_common_max3_contract(pr_backend.capability());
        require_common_max3_contract(sw_backend.capability());
        require_common_max3_contract(cpa_backend.capability());

        require(pr_backend.capability().phase_metadata_namespace.empty() &&
                    sw_backend.capability().phase_metadata_namespace.empty() &&
                    cpa_backend.capability().phase_metadata_namespace.empty(),
                "role-neutral backend freeze exposed provider phase identities");
        require(sw_backend.capability().publication_profile ==
                    fl::sw92_pt_flash_publication_convention,
                "SW role-neutral publication profile changed");
        require(pr_backend.capability().scalar_settings.empty() &&
                    cpa_backend.capability().scalar_settings.empty(),
                "generic freeze invented PR/CPA scalar configuration");
        require(sw_backend.capability().scalar_settings.size() == 1U &&
                    sw_backend.capability().scalar_settings.front().id ==
                        "nacl_molality_mol_per_kg_water",
                "SW configured molality provenance changed");

        require(sw_backend.capability().transition_capability.edge(3U, 1U) != nullptr &&
                    sw_backend.capability().transition_capability.edge(3U, 1U)->support ==
                        fl::PtPhaseTransitionSupport::fresh_target_resolve,
                "SW-specific 3->1 fresh-neighbor capability changed");
        require(pr_backend.capability().transition_capability.edge(3U, 1U) == nullptr &&
                    cpa_backend.capability().transition_capability.edge(3U, 1U) == nullptr,
                "generic freeze invented unsupported PR/CPA 3->1 capability");

        require(!pr_backend.capability().global_stability_proven &&
                    !sw_backend.capability().global_stability_proven &&
                    !cpa_backend.capability().global_stability_proven,
                "finite backend search was promoted to global-stability proof");
        return 0;
    } catch (const std::exception&) {
        return 1;
    }
}
