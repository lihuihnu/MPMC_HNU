#include <mpmc/flash/sw92_phase_assigned_three_phase_closure.hpp>

#include "fixture.hpp"

#include <cmath>
#include <cstddef>
#include <iostream>
#include <source_location>
#include <stdexcept>
#include <string>
#include <string_view>

bool sw92_phase_assigned_three_phase_closure_header();

namespace {
namespace fl = mpmc::flash;

void require(bool condition, std::string_view message,
             std::source_location where = std::source_location::current()) {
    if (!condition) {
        throw std::runtime_error(std::string(where.file_name()) + ":" +
                                 std::to_string(where.line()) + ": " +
                                 std::string(message));
    }
}

void near(double a, double b, double relative = 2e-8,
          double absolute = 2e-11,
          std::source_location where = std::source_location::current()) {
    if (!std::isfinite(a) || !std::isfinite(b) ||
        std::abs(a - b) > absolute + relative * std::abs(b)) {
        std::cerr << "actual=" << a << " expected=" << b << '\n';
        require(false, "numeric mismatch", where);
    }
}

void synthetic_w_present_closure() {
    const auto chain = c2b2_test::solve_chain();
    require(chain.c2b1.candidate() != nullptr,
            "C2b.1 structural source is not an interior candidate");
    const auto result = c2b2_test::review_chain(chain);
    require(result.status == fl::Sw92PhaseAssignedC2b2Status::
                                 w_present_h_multiplicity_locally_closed &&
                result.w_present_locally_closed() && result.nonaqueous_search &&
                result.usable_additional_h_witness_count() == 0U &&
                result.selected_c2a1_witness_index == chain.witness_index,
            "C2b.2 did not locally close the structural W-present H multiplicity");
    require(!result.global_stability_proven &&
                !result.accepted_phase_set_published && !result.morphology_resolved,
            "C2b.2 exceeded candidate/review publication semantics");
    require(result.revalidation_property_evaluations >= 5U,
            "C2b.2 did not fresh-evaluate retained C1/C2b.1 roots");
    require(result.candidate_c2b1_reduced_gibbs <=
                result.source_c1_reduced_gibbs + result.nested_gibbs_guard,
            "nested C2b.1 Gibbs rose above C1 beyond arithmetic guard");
    require(result.retained_h0_trivial_point && result.retained_h1_trivial_point,
            "retained H common-tangent diagnostics missing");
    require(std::abs(result.retained_h0_trivial_point->value) <=
                result.common_reference_allowance +
                    result.retained_h0_trivial_point->roundoff_guard &&
                std::abs(result.retained_h1_trivial_point->value) <=
                result.common_reference_allowance +
                    result.retained_h1_trivial_point->roundoff_guard,
            "retained H phases are not on the reconstructed final tangent");
    require(result.effective_tpd_tolerance >= result.base_tpd_tolerance,
            "C2b.2 effective TPD tolerance lost common-reference allowance");
    for (const auto& witness : result.negative_witnesses) {
        require(!witness.usable_additional_h_witness(),
                "locally closed result retained an admissible additional-H witness");
    }
}

void witness_classification_guards() {
    const auto chain = c2b2_test::solve_chain();
    const auto& state = *chain.c2b1.candidate();
    const auto& feed = chain.c1.feed;
    const std::size_t water_index = 2U;

    fl::TpdPoint w_like;
    w_like.composition = state.aqueous_phase.composition;
    w_like.value = -1.0;
    const auto rejected = fl::detail::sw92_phase_assigned_c2b2_classify_na_witness(
        0U, w_like, state.hydrocarbon0_phase.composition,
        state.hydrocarbon1_phase.composition, feed,
        state.aqueous_phase.composition[water_index], water_index, 1e-7);
    require(!rejected.water_role_admissible &&
                !rejected.usable_additional_h_witness(),
            "W-like negative NA trial passed the C2b.2 physical-role guard");

    fl::TpdPoint duplicate;
    duplicate.composition = state.hydrocarbon0_phase.composition;
    duplicate.value = -1.0;
    const auto duplicate_result =
        fl::detail::sw92_phase_assigned_c2b2_classify_na_witness(
            1U, duplicate, state.hydrocarbon0_phase.composition,
            state.hydrocarbon1_phase.composition, feed,
            state.aqueous_phase.composition[water_index], water_index, 1e-7);
    require(!duplicate_result.compositionally_distinct_from_h0 &&
                !duplicate_result.usable_additional_h_witness(),
            "retained H0 was misclassified as an additional phase");

    fl::TpdPoint distinct;
    distinct.composition = {0.30, 0.69, 0.01};
    distinct.value = -1.0;
    const auto usable = fl::detail::sw92_phase_assigned_c2b2_classify_na_witness(
        2U, distinct, state.hydrocarbon0_phase.composition,
        state.hydrocarbon1_phase.composition, feed,
        state.aqueous_phase.composition[water_index], water_index, 1e-7);
    require(usable.compositionally_distinct_from_h0 &&
                usable.compositionally_distinct_from_h1 &&
                usable.water_role_admissible &&
                usable.usable_additional_h_witness(),
            "distinct water-poor negative trial did not become additional-H evidence");
}

void disappearance_routes() {
    const auto h_chain = c2b2_test::solve_chain(0.60L, false, 0.25);
    require(h_chain.c2b1.status ==
                fl::Sw92PhaseAssignedThreePhaseStatus::hydrocarbon_phase_disappearance,
            "controlled H-disappearance source changed");
    const auto h_route = c2b2_test::review_chain(h_chain);
    require(h_route.status == fl::Sw92PhaseAssignedC2b2Status::route_to_w_h &&
                !h_route.nonaqueous_search,
            "one-H disappearance did not route back to W+H review");

    const auto w_chain = c2b2_test::solve_chain(0.80L, false, 0.25);
    require(w_chain.c2b1.status ==
                fl::Sw92PhaseAssignedThreePhaseStatus::aqueous_phase_disappearance,
            "controlled W-disappearance source changed");
    const auto w_route = c2b2_test::review_chain(w_chain);
    require(w_route.status ==
                fl::Sw92PhaseAssignedC2b2Status::route_to_no_w_h0_h1 &&
                !w_route.nonaqueous_search,
            "W disappearance did not route to the no-W H0+H1 path");

    const auto endpoint_chain = c2b2_test::solve_chain(0.73L, false, 0.28);
    require(endpoint_chain.c2b1.status ==
                fl::Sw92PhaseAssignedThreePhaseStatus::aqueous_phase_disappearance,
            "controlled one-phase endpoint source changed");
    const auto endpoint = c2b2_test::review_chain(endpoint_chain);
    require(endpoint.status ==
                fl::Sw92PhaseAssignedC2b2Status::single_phase_endpoint_unresolved &&
                !endpoint.nonaqueous_search,
            "two disappearing phases were incorrectly routed as an accepted lower topology");
}

void component_permutation() {
    const auto normal_chain = c2b2_test::solve_chain();
    const auto reverse_chain = c2b2_test::solve_chain(0.73L, true);
    const auto normal = c2b2_test::review_chain(normal_chain);
    const auto reverse = c2b2_test::review_chain(reverse_chain, true);
    require(normal.w_present_locally_closed() && reverse.w_present_locally_closed(),
            "component permutation changed C2b.2 closure status");
    require(normal.water_index == 2U && reverse.water_index == 0U,
            "water index did not follow ordered component permutation");
    near(reverse.source_c1_reduced_gibbs, normal.source_c1_reduced_gibbs);
    near(reverse.candidate_c2b1_reduced_gibbs,
         normal.candidate_c2b1_reduced_gibbs);
}

void resource_and_root_failures() {
    const auto chain = c2b2_test::solve_chain();
    fl::Sw92PhaseAssignedC2b2Options limited;
    limited.stability.automatic_starts = false;
    limited.stability.max_evaluations = 1U;
    const auto resource = c2b2_test::review_chain(chain, false, limited);
    require(resource.status == fl::Sw92PhaseAssignedC2b2Status::indeterminate &&
                resource.nonaqueous_search,
            "C2b.2 stability resource exhaustion was hidden");

    auto root_chain = c2b2_test::solve_chain();
    root_chain.c2b1.options.nonaqueous_root_options.max_iterations = 1;
    const auto root = c2b2_test::review_chain(root_chain);
    require(root.status == fl::Sw92PhaseAssignedC2b2Status::indeterminate &&
                root.property_issue.has_value(),
            "C2b.2 fresh root/property failure was not retained");
}

void source_chain_guards() {
    auto unavailable = c2b2_test::solve_chain();
    unavailable.c2b1.source_witness_index.reset();
    const auto missing = c2b2_test::review_chain(unavailable);
    require(missing.status ==
                fl::Sw92PhaseAssignedC2b2Status::source_chain_unavailable,
            "C2b.2 ran without C2b.1 source-witness provenance");

    auto inconsistent = c2b2_test::solve_chain();
    inconsistent.c2b1.dataset_id = "different-dataset";
    const auto mismatched = c2b2_test::review_chain(inconsistent);
    require(mismatched.status ==
                fl::Sw92PhaseAssignedC2b2Status::source_chain_inconsistent,
            "C2b.2 accepted mismatched model metadata");
}

void headers() {
    require(sw92_phase_assigned_three_phase_closure_header(),
            "C2b.2 public header probe failed");
}

} // namespace

int main(int argc, char** argv) {
    if (argc != 2) {
        std::cerr << "expected one test-case name\n";
        return 2;
    }
    const std::string_view name{argv[1]};
    try {
        if (name == "synthetic_locally_closed") synthetic_w_present_closure();
        else if (name == "witness_guards") witness_classification_guards();
        else if (name == "disappearance_routes") disappearance_routes();
        else if (name == "permutation") component_permutation();
        else if (name == "resource_root_failures") resource_and_root_failures();
        else if (name == "source_guards") source_chain_guards();
        else if (name == "headers") headers();
        else throw std::invalid_argument("unknown test-case name");
    } catch (const std::exception& error) {
        std::cerr << error.what() << '\n';
        return 1;
    }
    return 0;
}
