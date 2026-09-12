#include <mpmc/flash/pr76_pt_flash_backend.hpp>

#include <array>
#include <cmath>
#include <optional>
#include <stdexcept>
#include <string>
#include <vector>

namespace {
namespace fl = mpmc::flash;
namespace th = mpmc::thermodynamics;
using Vec = std::vector<double>;

th::Provenance source() {
    return {th::SourceKind::synthetic_test,
            "MPMC_HNU symmetric PR76 max3 backend fixture", "v1",
            "three-component symmetric LLV-like structural fixture",
            "Artificial values for backend conformance only; not experimental validation",
            "tests/flash/pt_flash_backend/pr76_max3_backend_test.cpp",
            "Repository structural regression"};
}

th::SourcedScalar scalar(double value, th::Unit unit,
                         const th::Provenance& provenance) {
    return {value, unit, provenance, "synthetic SI or dimensionless", "identity"};
}

th::Pr76Phase<double> model() {
    const auto provenance = source();
    const std::array<std::string, 3> ids{"light", "heavy-b", "heavy-c"};
    const std::array<double, 3> tc{190.6, 500.0, 500.0};
    const std::array<double, 3> pc{4.6e6, 5.0e6, 5.0e6};
    const std::array<double, 3> omega{0.01, 0.10, 0.10};
    std::vector<th::Component> catalog;
    th::PrParameterInput input;
    input.model_id = std::string(th::pr76_profile);
    input.dataset_id = "synthetic-PR76-max3-symmetric";
    input.revision = "v1";
    input.applicability = {std::nullopt, std::nullopt, provenance};
    for (std::size_t i = 0; i < ids.size(); ++i) {
        catalog.push_back({ids[i], ids[i], th::ComponentKind::pure, provenance, {}});
        input.pure.push_back({ids[i], scalar(tc[i], th::Unit::kelvin, provenance),
                              scalar(pc[i], th::Unit::pascal, provenance),
                              scalar(omega[i], th::Unit::dimensionless, provenance)});
    }
    input.binary.push_back({ids[0], ids[1], scalar(0.05, th::Unit::dimensionless, provenance)});
    input.binary.push_back({ids[0], ids[2], scalar(0.05, th::Unit::dimensionless, provenance)});
    input.binary.push_back({ids[1], ids[2], scalar(0.20, th::Unit::dimensionless, provenance)});
    const std::vector<std::string> order{ids.begin(), ids.end()};
    return th::Pr76Phase<double>::from_parameters(
        th::PrParameterSet::create(catalog, order, input,
                                   th::DataPolicy::allow_synthetic_tests));
}

bool phase_equal(const fl::PtCandidatePhase& a, const fl::PtCandidatePhase& b) {
    return a.mole_phase_fraction == b.mole_phase_fraction &&
           a.composition == b.composition &&
           a.activity.branch == b.activity.branch &&
           a.activity.smooth == b.activity.smooth &&
           a.activity.ln_phi == b.activity.ln_phi &&
           a.compressibility_factor == b.compressibility_factor;
}

} // namespace

int main() {
    try {
        const auto phase_model = model();
        fl::Pr76VleEvaluator evaluator(phase_model);
        const Vec feed{1.0 / 3.0, 1.0 / 3.0, 1.0 / 3.0};
        const auto direct = fl::solve_pr76_pt_max3(1.0e6, 250.0, feed, evaluator);
        const auto published = fl::project_pr76_pt_max3_phase_set(direct);
        fl::Pr76PtFlashBackend backend(evaluator);
        fl::PtFlashBackend& runtime = backend;
        const auto adapted = runtime.solve({1.0e6, 250.0, feed});

        if (!adapted.structurally_valid() ||
            !runtime.capability().supports_phase_count(3U) ||
            runtime.capability().maximum_phase_count() != 3U ||
            runtime.capability().transition_capability.edge(2U, 3U) == nullptr ||
            runtime.capability().transition_capability.edge(3U, 2U) == nullptr ||
            adapted.solution.status != published.solution.status ||
            adapted.solution.accepted_phase_count() !=
                published.solution.accepted_phase_count() ||
            adapted.solution.accepted_phase_count() != 3U ||
            adapted.solution.feed != published.solution.feed) {
            throw std::runtime_error("PR76 max3 backend capability/publication mismatch");
        }
        const auto* a = adapted.solution.accepted_phase_set();
        const auto* b = published.solution.accepted_phase_set();
        if (a == nullptr || b == nullptr || a->phases.size() != b->phases.size()) {
            throw std::runtime_error("PR76 max3 backend accepted set missing");
        }
        for (std::size_t i = 0; i < a->phases.size(); ++i) {
            if (!phase_equal(a->phases[i], b->phases[i])) {
                throw std::runtime_error("PR76 max3 backend changed phase payload");
            }
        }
        const bool has_2_to_3 = std::any_of(
            adapted.transition_report.evidence.begin(),
            adapted.transition_report.evidence.end(),
            [](const fl::PtPhaseTransitionEvidence& evidence) {
                return evidence.source_phase_count == 2U &&
                       evidence.target_phase_count == 3U &&
                       evidence.resolution ==
                           fl::PtPhaseTransitionResolution::accepted_target &&
                       evidence.fresh_target_solve_attempted &&
                       evidence.target_topology_closed;
            });
        if (!has_2_to_3) {
            throw std::runtime_error("PR76 max3 backend lost 2->3 accepted transition evidence");
        }
        return 0;
    } catch (const std::exception&) {
        return 1;
    }
}
