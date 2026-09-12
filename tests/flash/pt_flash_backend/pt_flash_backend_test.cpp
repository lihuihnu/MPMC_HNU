#include <mpmc/flash/pr76_pt_flash_backend.hpp>
#include <mpmc/flash/sw92_profile_c_pt_flash_backend.hpp>

#include "../sw92_phase_assigned_pt/physical_sample6.hpp"
#include "test_support.hpp"

#include <array>
#include <cmath>
#include <cstddef>
#include <iostream>
#include <optional>
#include <source_location>
#include <stdexcept>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

bool pt_flash_backend_headers();

namespace {
namespace fl = mpmc::flash;
namespace th = mpmc::thermodynamics;
namespace sample6 = sw92_profile_c_sample6;
using Vec = std::vector<double>;

void require(bool value, std::string_view message,
             std::source_location where = std::source_location::current()) {
    if (!value) {
        throw std::runtime_error(
            std::string(where.file_name()) + ":" + std::to_string(where.line()) +
            ": " + std::string(message));
    }
}

template <class Error, class Function>
void expect_error(Function&& function) {
    bool caught = false;
    try {
        function();
    } catch (const Error&) {
        caught = true;
    }
    require(caught, "expected exception missing");
}

th::SourcedScalar sourced_scalar(double value, th::Unit unit,
                                  const th::Provenance& source,
                                  std::string original_unit,
                                  std::string conversion) {
    return {value, unit, source, std::move(original_unit), std::move(conversion)};
}

th::Pr76Phase<double> pr76_binary_model() {
    const th::Provenance source{
        th::SourceKind::literature,
        "https://academicweb.nd.edu/~markst/zm97a.pdf",
        "revised October 1997",
        "Section 4.3 and Table 4",
        "PR numerical regression, not experimental validation",
        "Read author-hosted PDF pages",
        "Limited attributed factual parameters; no article text or code reproduced"};
    const std::array<std::string, 2> ids{"nitrogen", "ethane"};
    const std::array<std::array<double, 3>, 2> specs{{
        {126.2, 3.39e6, 0.04},
        {305.4, 4.88e6, 0.098}}};
    std::vector<th::Component> catalog;
    th::PrParameterInput input;
    input.model_id = std::string(th::pr76_profile);
    input.dataset_id = "Hua1997-nitrogen-ethane";
    input.revision = "split-regression-v1";
    input.applicability = {std::nullopt, std::nullopt, source};
    for (std::size_t i = 0; i < ids.size(); ++i) {
        catalog.push_back({ids[i], ids[i], th::ComponentKind::pure, source, {}});
        input.pure.push_back({
            ids[i],
            sourced_scalar(specs[i][0], th::Unit::kelvin, source, "K", "identity"),
            sourced_scalar(specs[i][1], th::Unit::pascal, source, "bar",
                           "bar * 100000 -> Pa"),
            sourced_scalar(specs[i][2], th::Unit::dimensionless, source,
                           "dimensionless", "identity")});
    }
    input.binary.push_back({
        ids[0], ids[1],
        sourced_scalar(0.08, th::Unit::dimensionless, source,
                       "dimensionless", "identity")});
    const std::vector<std::string> order{ids.begin(), ids.end()};
    return th::Pr76Phase<double>::from_parameters(
        th::PrParameterSet::create(catalog, order, input));
}

void require_phase_equal(const fl::PtCandidatePhase& actual,
                         const fl::PtCandidatePhase& expected) {
    require(actual.mole_phase_fraction == expected.mole_phase_fraction,
            "phase fraction changed through backend adapter");
    require(actual.composition == expected.composition,
            "phase composition changed through backend adapter");
    require(actual.compressibility_factor == expected.compressibility_factor,
            "phase Z changed through backend adapter");
    require(actual.activity.branch == expected.activity.branch &&
                actual.activity.smooth == expected.activity.smooth &&
                actual.activity.ln_phi == expected.activity.ln_phi,
            "phase activity changed through backend adapter");
}

void require_phase_set_equal(const fl::PtPhaseSetResult& actual,
                             const fl::PtPhaseSetResult& expected) {
    require(actual.status == expected.status &&
                actual.capability.maximum_phase_count ==
                    expected.capability.maximum_phase_count &&
                actual.pressure_pa == expected.pressure_pa &&
                actual.temperature_k == expected.temperature_k &&
                actual.feed == expected.feed &&
                actual.global_stability_proven ==
                    expected.global_stability_proven &&
                actual.diagnostic == expected.diagnostic,
            "generic phase-set metadata changed through backend adapter");
    require(actual.candidate_phase_set.has_value() ==
                expected.candidate_phase_set.has_value(),
            "candidate presence changed through backend adapter");
    if (!actual.candidate_phase_set) { return; }
    require(actual.candidate_phase_set->phases.size() ==
                expected.candidate_phase_set->phases.size(),
            "candidate phase count changed through backend adapter");
    for (std::size_t phase = 0;
         phase < actual.candidate_phase_set->phases.size(); ++phase) {
        require_phase_equal(actual.candidate_phase_set->phases[phase],
                            expected.candidate_phase_set->phases[phase]);
    }
}

std::string expected_role_id(fl::Sw92PhaseAssignedPtPhysicalRole role) {
    switch (role) {
    case fl::Sw92PhaseAssignedPtPhysicalRole::aqueous:
        return "aqueous";
    case fl::Sw92PhaseAssignedPtPhysicalRole::nonaqueous_unclassified:
        return "nonaqueous_unclassified";
    }
    return "unknown";
}

std::string expected_family_id(th::SwPhaseFamily family) {
    switch (family) {
    case th::SwPhaseFamily::aqueous:
        return "aqueous";
    case th::SwPhaseFamily::nonaqueous:
        return "nonaqueous";
    }
    return "unknown";
}

void pr76_direct_equivalence() {
    const auto model = pr76_binary_model();
    fl::Pr76VleEvaluator evaluator(model);
    const fl::PtFlashRequest request{7.6e6, 270.0, Vec{0.30, 0.70}};

    const auto direct = fl::solve_pr76_pt_phase_set(
        request.pressure_pa, request.temperature_k, request.feed, evaluator);
    fl::Pr76PtFlashBackend backend(evaluator);
    fl::PtFlashBackend& runtime_backend = backend;
    const auto adapted = runtime_backend.solve(request);

    const auto& capability = runtime_backend.capability();
    require(capability.structurally_valid(), "PR76 capability invalid");
    require(capability.backend_id == fl::pr76_pt_flash_backend_id &&
                capability.model_profile == th::pr76_profile &&
                capability.algorithm_profile == fl::PtSplitResult::convention &&
                capability.publication_profile ==
                    fl::PtPhaseSetResult::convention &&
                capability.supports_phase_count(1U) &&
                capability.supports_phase_count(2U) &&
                !capability.supports_phase_count(3U) &&
                capability.maximum_phase_count() == 2U &&
                capability.performs_initial_stability_search &&
                capability.performs_final_phase_set_review &&
                !capability.performs_boundary_neighbor_resolve &&
                !capability.global_stability_proven &&
                capability.phase_metadata_namespace.empty(),
            "PR76 capability semantics changed");
    require(adapted.structurally_valid(), "PR76 backend result invalid");
    require(adapted.phase_metadata.empty() && !adapted.morphology_resolved,
            "PR76 backend invented provider phase morphology");
    require_phase_set_equal(adapted.solution, direct.solution);
}

void sw92_direct_equivalence() {
    const auto model = th::Sw92Phase<double>::from_parameters(
        sw92_test::binary_parameters(sw92_test::co2));
    const fl::PtFlashRequest request{3.0e6, 340.0, Vec{0.7, 0.3}};
    constexpr double molality = 0.0;

    const auto direct = fl::solve_sw92_profile_c_pt_phase_set(
        request.pressure_pa, request.temperature_k, request.feed,
        model, molality);
    fl::Sw92ProfileCPtFlashBackendOptions options;
    options.nacl_molality_mol_per_kg_water = molality;
    fl::Sw92ProfileCPtFlashBackend backend(model, options);
    fl::PtFlashBackend& runtime_backend = backend;
    const auto adapted = runtime_backend.solve(request);

    const auto& capability = runtime_backend.capability();
    require(capability.structurally_valid(), "SW92 capability invalid");
    require(capability.backend_id == fl::sw92_profile_c_pt_flash_backend_id &&
                capability.model_profile == th::sw92_corrected_profile &&
                capability.algorithm_profile ==
                    fl::sw92_phase_assigned_pt_convention &&
                capability.publication_profile ==
                    fl::sw92_profile_c_phase_set_publication_convention &&
                capability.supports_phase_count(1U) &&
                capability.supports_phase_count(2U) &&
                capability.supports_phase_count(3U) &&
                capability.maximum_phase_count() == 3U &&
                capability.performs_initial_stability_search &&
                capability.performs_final_phase_set_review &&
                capability.performs_boundary_neighbor_resolve &&
                !capability.global_stability_proven &&
                capability.phase_metadata_namespace ==
                    fl::sw92_profile_c_phase_metadata_namespace,
            "SW92 capability semantics changed");
    require(adapted.structurally_valid(), "SW92 backend result invalid");
    require(!adapted.morphology_resolved,
            "SW92 backend invented resolved H morphology");
    require_phase_set_equal(adapted.solution, direct.solution);
    require(adapted.phase_metadata.size() == direct.phase_metadata.size(),
            "SW92 phase metadata count changed");
    for (std::size_t phase = 0; phase < direct.phase_metadata.size(); ++phase) {
        require(adapted.phase_metadata[phase].role_id ==
                    expected_role_id(direct.phase_metadata[phase].physical_role) &&
                    adapted.phase_metadata[phase].family_id ==
                    expected_family_id(
                        direct.phase_metadata[phase].thermodynamic_family),
                "SW92 opaque phase metadata changed");
    }
}

void sw92_three_phase_capability() {
    const auto model = sample6::model();
    fl::Sw92ProfileCPtFlashBackend backend(model);
    const fl::PtFlashRequest request{1.0e7, 350.0, sample6::feed()};
    const auto direct = fl::solve_sw92_profile_c_pt_phase_set(
        request.pressure_pa, request.temperature_k, request.feed,
        model, 0.0);
    const auto adapted = backend.solve(request);

    require(direct.accepted_phase_set_published() &&
                direct.solution.accepted_phase_count() == 3U,
            "Sample-6 direct fixture is no longer authoritative three phase");
    require(adapted.structurally_valid() && adapted.accepted_phase_set() != nullptr &&
                adapted.accepted_phase_set()->phases.size() == 3U,
            "unified SW backend failed to publish Sample-6 three phase");
    require_phase_set_equal(adapted.solution, direct.solution);
    require(adapted.phase_metadata.size() == 3U &&
                adapted.phase_metadata[0].role_id == "aqueous" &&
                adapted.phase_metadata[0].family_id == "aqueous" &&
                adapted.phase_metadata[1].role_id ==
                    "nonaqueous_unclassified" &&
                adapted.phase_metadata[1].family_id == "nonaqueous" &&
                adapted.phase_metadata[2].role_id ==
                    "nonaqueous_unclassified" &&
                adapted.phase_metadata[2].family_id == "nonaqueous",
            "Sample-6 backend metadata invented hydrocarbon morphology");
}

void contract_guards() {
    fl::PtFlashBackendCapability capability;
    require(!capability.structurally_valid(),
            "default capability unexpectedly valid");
    capability.backend_id = "test";
    capability.model_profile = "model";
    capability.algorithm_profile = "algorithm";
    capability.publication_profile = "publication";
    capability.dataset_id = "dataset";
    capability.revision = "revision";
    capability.component_ids = {"a", "b"};
    capability.supported_phase_counts = {1U, 2U};
    require(capability.structurally_valid(),
            "minimal capability unexpectedly invalid");
    capability.supported_phase_counts = {1U, 1U};
    require(!capability.structurally_valid(),
            "duplicate supported phase count accepted");

    const auto model = pr76_binary_model();
    fl::Pr76VleEvaluator evaluator(model);
    fl::Pr76PtFlashBackend backend(evaluator);
    expect_error<std::invalid_argument>([&] {
        (void)backend.solve(fl::PtFlashRequest{7.6e6, 270.0, Vec{1.0}});
    });
}

void headers() {
    require(pt_flash_backend_headers(), "public header probe failed");
}

using Test = std::pair<std::string_view, void (*)()>;
constexpr Test tests[]{
    {"pr76_direct_equivalence", pr76_direct_equivalence},
    {"sw92_direct_equivalence", sw92_direct_equivalence},
    {"sw92_three_phase_capability", sw92_three_phase_capability},
    {"contract_guards", contract_guards},
    {"headers", headers}};

} // namespace

int main(int argc, char** argv) {
    try {
        if (argc != 2) {
            throw std::invalid_argument("one test name required");
        }
        for (const auto& [name, run] : tests) {
            if (name == argv[1]) {
                run();
                std::cout << "[PASS] " << name << '\n';
                return 0;
            }
        }
        throw std::invalid_argument("unknown test");
    } catch (const std::exception& error) {
        std::cerr << "[FAIL] " << error.what() << '\n';
        return 1;
    }
}
