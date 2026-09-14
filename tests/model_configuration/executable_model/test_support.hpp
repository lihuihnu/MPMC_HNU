#ifndef MPMC_TEST_EXECUTABLE_MODEL_SUPPORT_HPP
#define MPMC_TEST_EXECUTABLE_MODEL_SUPPORT_HPP

#include <mpmc/model_configuration/pr76_executable_model.hpp>
#include "../../flash/pr76_three_phase/synthetic_fixture.hpp"

#include <algorithm>
#include <array>
#include <limits>
#include <stdexcept>
#include <string>
#include <string_view>
#include <vector>

namespace model_test {
namespace mc = mpmc::model_configuration;
namespace th = mpmc::thermodynamics;
namespace fl = mpmc::flash;

inline void require(bool value, std::string_view message) {
    if (!value) { throw std::runtime_error(std::string(message)); }
}

// Reuse the repository's attributed binary numerical fixture (not experiment):
// tests/flash/pt_flash_backend/pt_flash_backend_test.cpp::pr76_binary_model.
// Native construction is independent of the public DTO mapper under test.
inline th::Pr76Phase<double> binary_model() {
    const th::Provenance source{
        th::SourceKind::literature, "https://academicweb.nd.edu/~markst/zm97a.pdf",
        "revised October 1997", "Section 4.3 and Table 4",
        "PR numerical regression, not experimental validation",
        "Read author-hosted PDF pages",
        "Limited attributed factual parameters; no article text or code reproduced"};
    const std::array<std::string, 2> ids{"nitrogen", "ethane"};
    const std::array<std::array<double, 3>, 2> specs{{
        {126.2, 3.39e6, 0.04}, {305.4, 4.88e6, 0.098}}};
    std::vector<th::Component> catalog;
    th::PrParameterInput input;
    input.model_id = std::string(th::pr76_profile);
    input.dataset_id = "Hua1997-nitrogen-ethane";
    input.revision = "split-regression-v1";
    input.applicability = {std::nullopt, std::nullopt, source};
    const auto scalar = [&](double value, th::Unit unit, std::string original, std::string conversion) {
        return th::SourcedScalar{value, unit, source, std::move(original), std::move(conversion)};
    };
    for (std::size_t i = 0; i < ids.size(); ++i) {
        catalog.push_back({ids[i], ids[i], th::ComponentKind::pure, source, {}});
        input.pure.push_back({ids[i], scalar(specs[i][0], th::Unit::kelvin, "K", "identity"),
            scalar(specs[i][1], th::Unit::pascal, "bar", "bar * 100000 -> Pa"),
            scalar(specs[i][2], th::Unit::dimensionless, "dimensionless", "identity")});
    }
    input.binary.push_back({ids[0], ids[1],
        scalar(0.08, th::Unit::dimensionless, "dimensionless", "identity")});
    return th::Pr76Phase<double>::from_parameters(th::PrParameterSet::create(
        catalog, std::vector<std::string>{ids.begin(), ids.end()}, input));
}

// Modified synthetic fixture B; retain its artificial source classification.
inline th::Pr76Phase<double> changed_model(bool bounded = false) {
    const auto original = pr76_max3_test::model();
    const auto& p = original.parameters();
    th::PrParameterInput input;
    input.model_id = std::string(th::pr76_profile);
    input.dataset_id = p.dataset_id();
    input.revision = bounded ? "bounded-test" : "isolation-B";
    input.applicability = p.applicability();
    input.pure.assign(p.pure_records().begin(), p.pure_records().end());
    input.binary.assign(p.binary_records().begin(), p.binary_records().end());
    if (bounded) {
        input.applicability.temperature_k = th::ClosedInterval{240.0, 260.0};
        input.applicability.pressure_pa = th::ClosedInterval{0.9e6, 1.1e6};
    } else {
        input.binary[0].kij->value = 0.15;
        input.binary[0].kij->source.note += "; changed kij for model-isolation software test";
    }
    std::vector<std::string> order;
    for (const auto& c : p.components().items()) { order.push_back(c.id); }
    return th::Pr76Phase<double>::from_parameters(th::PrParameterSet::create(
        p.components().items(), order, input, th::DataPolicy::allow_synthetic_tests));
}

inline mc::ModelProvenance public_source(const th::Provenance& s) {
    require(s.kind == th::SourceKind::synthetic_test || s.kind == th::SourceKind::literature,
            "fixture source kind unsupported");
    return {s.kind == th::SourceKind::synthetic_test ? mc::SourceKind::synthetic_test
                                                  : mc::SourceKind::literature,
            s.reference, s.revision, s.locator, s.note, s.acquisition, s.usage_terms};
}
inline mc::ModelScalar public_scalar(const th::SourcedScalar& s) {
    return {s.value, public_source(s.source), s.original_unit, s.conversion};
}
inline mc::ThermodynamicModelDefinition definition(const th::Pr76Phase<double>& native) {
    const auto& p = native.parameters();
    mc::ThermodynamicModelDefinition d;
    d.version = mc::model_parameter_definition_v1;
    d.family = mc::ThermodynamicModelFamily::peng_robinson_1976;
    d.display_name = p.dataset_id();
    d.dataset_id = p.dataset_id();
    d.revision = p.revision();
    d.provenance = public_source(p.applicability().declaration);
    d.applicability.provenance = d.provenance;
    if (p.applicability().temperature_k) {
        d.applicability.temperature_lower_k = p.applicability().temperature_k->lower;
        d.applicability.temperature_upper_k = p.applicability().temperature_k->upper;
    }
    if (p.applicability().pressure_pa) {
        d.applicability.pressure_lower_pa = p.applicability().pressure_pa->lower;
        d.applicability.pressure_upper_pa = p.applicability().pressure_pa->upper;
    }
    for (const auto& component : p.components().items()) {
        d.components.push_back({component.id, component.display_name, mc::ComponentKind::pure,
                                public_source(component.definition), {}});
    }
    mc::Pr76ParameterDefinition params;
    for (const auto& r : p.pure_records()) {
        params.pure.push_back({r.component_id, public_scalar(*r.critical_temperature),
            public_scalar(*r.critical_pressure), public_scalar(*r.acentric_factor)});
    }
    for (const auto& r : p.binary_records()) {
        params.binary.push_back({r.first_id, r.second_id, public_scalar(*r.kij)});
    }
    d.parameters = std::move(params);
    return d;
}
inline auto preset() { return mc::resolve_pt_solver_preset(mc::mpmc_balanced_default_v1); }
inline auto custom() {
    auto s = preset();
    s.kind = mc::PtSolverSettingsKind::custom;
    s.preset_id.clear();
    return s;
}
inline auto create(const mc::ThermodynamicModelDefinition& d, const mc::PtSolverSettings& s) {
    return mc::make_pr76_executable_model(d, s, {}, {}, mc::ModelDataPolicy::allow_synthetic_tests);
}
inline fl::PtFlashBackendResult direct(const th::Pr76Phase<double>& model,
    const fl::PtFlashRequest& request, th::Pr76RootOptions roots = {},
    fl::Pr76PtFlashBackendOptions options = {}) {
    fl::Pr76VleEvaluator evaluator(model, roots);
    fl::Pr76PtFlashBackend backend(evaluator, std::move(options));
    return backend.solve(request);
}
inline const fl::PtFlashRequest single{1e6, 250.0, {1.0, 0.0, 0.0}};
inline const fl::PtFlashRequest ternary{1e6, 250.0, pr76_max3_test::equal_feed()};
inline const fl::PtFlashRequest binary{7.6e6, 270.0, {0.30, 0.70}};

struct RequestLocationCase { fl::PtFlashRequest request; const char* field; };
inline std::vector<RequestLocationCase> request_location_cases() {
    const double nan = std::numeric_limits<double>::quiet_NaN();
    const double inf = std::numeric_limits<double>::infinity();
    const double eps = std::numeric_limits<double>::epsilon();
    std::vector<RequestLocationCase> cases;
    for (double invalid : {0.0, -1.0, nan, inf, -inf}) {
        cases.push_back({{invalid, 250.0, single.feed}, "pressure_pa"});
        cases.push_back({{1e6, invalid, single.feed}, "temperature_k"});
    }
    for (const auto& feed : std::vector<std::vector<double>>{
             {}, {1.0}, {0.0, 0.0, 0.0}, {0.4, 0.4, 0.4}, {0.5, 0.5 + 128.0 * eps, 0.0}}) {
        cases.push_back({{1e6, 250.0, feed}, "feed"});
    }
    for (double invalid : {-0.1, 1.1, nan, inf, -inf}) {
        cases.push_back({{1e6, 250.0, {0.5, invalid, 0.5}}, "feed[1]"});
    }
    cases.push_back({{1e6, 250.0, {0.5, 0.5, nan}}, "feed[2]"});
    // Stable first-invalid-field priority, without interpolating values/IDs.
    cases.push_back({{-1.0, -1.0, single.feed}, "pressure_pa"});
    cases.push_back({{-1.0, -1.0, {1.0}}, "feed"});
    return cases;
}


// Entire existing result envelope, including NON-accepted candidates/evidence.
// Exact equality tests adapter parity, not independent thermodynamic accuracy.
inline void compare_capability(const fl::PtFlashBackendCapability& a,
                               const fl::PtFlashBackendCapability& b) {
#define SAME(field) require(a.field == b.field, "capability." #field)
    SAME(backend_id); SAME(model_profile); SAME(algorithm_profile); SAME(publication_profile);
    SAME(configuration_profile); SAME(dataset_id); SAME(revision); SAME(component_ids);
    SAME(supported_phase_counts); SAME(performs_initial_stability_search);
    SAME(performs_final_phase_set_review); SAME(performs_boundary_neighbor_resolve);
    SAME(global_stability_proven); SAME(phase_metadata_namespace);
    SAME(scalar_settings.size()); SAME(transition_capability.edges.size());
#undef SAME
    for (std::size_t i = 0; i < a.scalar_settings.size(); ++i) {
        const auto& x = a.scalar_settings[i];
        const auto& y = b.scalar_settings[i];
        require(x.id == y.id && x.value == y.value && x.unit == y.unit, "capability scalar setting");
    }
    for (std::size_t i = 0; i < a.transition_capability.edges.size(); ++i) {
        const auto& x = a.transition_capability.edges[i];
        const auto& y = b.transition_capability.edges[i];
        require(x.source_phase_count == y.source_phase_count && x.target_phase_count == y.target_phase_count &&
                x.support == y.support && x.requires_fresh_target_solve == y.requires_fresh_target_solve,
                "capability transition edge");
    }
}
inline void compare(const fl::PtFlashBackendResult& a, const fl::PtFlashBackendResult& b) {
    require(a.structurally_valid() && b.structurally_valid(), "malformed result envelope");
    compare_capability(a.capability, b.capability);
#define SAME(field) require(a.field == b.field, "result." #field)
    SAME(provider_result_convention); SAME(morphology_resolved); SAME(phase_metadata.size());
    SAME(solution.capability.maximum_phase_count); SAME(solution.status);
    SAME(solution.pressure_pa); SAME(solution.temperature_k); SAME(solution.feed);
    SAME(solution.global_stability_proven); SAME(solution.diagnostic);
    SAME(solution.candidate_phase_set.has_value()); SAME(transition_report.evidence.size());
#undef SAME
    if (a.solution.candidate_phase_set) {
        const auto& left = a.solution.candidate_phase_set->phases;
        const auto& right = b.solution.candidate_phase_set->phases;
        require(left.size() == right.size(), "candidate phase count");
        for (std::size_t i = 0; i < left.size(); ++i) {
            const auto& x = left[i];
            const auto& y = right[i];
            require(x.mole_phase_fraction == y.mole_phase_fraction && x.composition == y.composition &&
                    x.compressibility_factor == y.compressibility_factor && x.activity.ln_phi == y.activity.ln_phi &&
                    x.activity.branch == y.activity.branch && x.activity.smooth == y.activity.smooth,
                    "candidate phase payload");
        }
    }
    for (std::size_t i = 0; i < a.phase_metadata.size(); ++i) {
        require(a.phase_metadata[i].role_id == b.phase_metadata[i].role_id &&
                a.phase_metadata[i].family_id == b.phase_metadata[i].family_id, "phase metadata");
    }
    for (std::size_t i = 0; i < a.transition_report.evidence.size(); ++i) {
        const auto& x = a.transition_report.evidence[i];
        const auto& y = b.transition_report.evidence[i];
        require(x.source_phase_count == y.source_phase_count && x.target_phase_count == y.target_phase_count &&
                x.trigger == y.trigger && x.resolution == y.resolution &&
                x.fresh_target_solve_attempted == y.fresh_target_solve_attempted &&
                x.target_topology_closed == y.target_topology_closed &&
                x.provider_evidence_profile == y.provider_evidence_profile && x.diagnostic == y.diagnostic,
                "transition evidence");
    }
}
} // namespace model_test
#endif
