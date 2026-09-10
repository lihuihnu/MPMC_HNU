#include <mpmc/flash/sw92_asymmetric_orchestration.hpp>

#include "test_support.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <iostream>
#include <source_location>
#include <stdexcept>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

bool sw92_asymmetric_orchestration_header();

namespace {
namespace fl = mpmc::flash;
namespace th = mpmc::thermodynamics;
using Vec = std::vector<double>;

void require(bool condition, std::string_view message,
             std::source_location where = std::source_location::current()) {
    if (!condition) {
        throw std::runtime_error(std::string(where.file_name()) + ":" +
                                 std::to_string(where.line()) + ": " +
                                 std::string(message));
    }
}

void near(double actual, long double expected, long double relative = 8e-9L,
          long double absolute = 8e-12L,
          std::source_location where = std::source_location::current()) {
    if (!std::isfinite(actual) ||
        std::abs(static_cast<long double>(actual) - expected) >
            absolute + relative * std::abs(expected)) {
        std::cerr << "actual=" << actual << " expected=" << expected << '\n';
        require(false, "reference mismatch", where);
    }
}

template<class Error, class Function>
void expect_error(Function&& function) {
    bool caught = false;
    try {
        function();
    } catch (const Error&) {
        caught = true;
    }
    require(caught, "expected exception missing");
}

th::Sw92Phase<double> binary_model(bool reverse = false) {
    return th::Sw92Phase<double>::from_parameters(
        sw92_test::binary_parameters(sw92_test::co2, reverse));
}

constexpr double aq_co2_kij_340k_fresh =
    -0.058970986957265289338095088569696338832153503254491;

th::Sw92Phase<double> family_tie_model() {
    auto prepared = sw92_test::binary_input(sw92_test::co2);
    const auto source = sw92_test::synthetic(
        "NA CO2/water BIP set equal to AQ value at T=340 K,m=0 for Gate 3B.2 family-tie routing test");
    prepared.input.water_binary.front().nonaqueous_rule =
        th::Sw92NonAqueousWaterRule::constant;
    prepared.input.water_binary.front().nonaqueous_kij = sw92_test::scalar(
        aq_co2_kij_340k_fresh, th::Unit::dimensionless, source);
    return th::Sw92Phase<double>::from_parameters(
        th::Sw92ParameterSet::create(
            prepared.catalog, prepared.order, prepared.input,
            th::DataPolicy::allow_synthetic_tests));
}

constexpr long double aq_aq_phase0_gas =
    0.005936172024213170056428932801071659387598344017016019414L;
constexpr long double aq_aq_phase1_gas =
    0.9874123977699055456489128208133867393111524204260154840L;
constexpr long double aq_aq_reduced_gibbs =
    -1.496025765990406698063759588412595324322460168711074814L;

bool same_pair(const fl::Sw92AsymmetricFamilyPair& first,
               const fl::Sw92AsymmetricFamilyPair& second) {
    return first.phase0 == second.phase0 && first.phase1 == second.phase1;
}

void plan_and_same_family_alternatives() {
    const auto model = binary_model();
    auto initial = fl::test_sw92_pt_asymmetric_stability(
        3.0e6, 340.0, Vec{0.7, 0.3}, model, 0.0);
    require(initial.status == fl::StabilityStatus::unstable &&
                initial.reference_family == th::SwPhaseFamily::aqueous,
            "reference asymmetric state changed");

    bool has_aq = false;
    bool has_na = false;
    for (const auto& witness : initial.negative_witnesses) {
        has_aq = has_aq || witness.family == th::SwPhaseFamily::aqueous;
        has_na = has_na || witness.family == th::SwPhaseFamily::nonaqueous;
    }
    require(has_aq && has_na, "reference state lost one family of negative evidence");

    const auto plan = fl::build_sw92_asymmetric_pair_plan(initial);
    require(plan.size() > initial.negative_witnesses.size(),
            "opposite-family same-family alternatives were not added");
    require(plan.size() >= 2 &&
                plan[0].protected_family_representative &&
                plan[1].protected_family_representative &&
                plan[0].kind == fl::Sw92AsymmetricPairPlanKind::reference_with_witness &&
                plan[1].kind == fl::Sw92AsymmetricPairPlanKind::reference_with_witness &&
                plan[0].witness_family != plan[1].witness_family,
            "best evidence from both families is not protected before alternatives");

    for (const auto& witness : initial.negative_witnesses) {
        bool primary = false;
        bool same_family = witness.family == th::SwPhaseFamily::aqueous;
        for (const auto& entry : plan) {
            if (entry.witness_family != witness.family ||
                entry.witness_trial_index != witness.trial_index) {
                continue;
            }
            require(entry.initial_log_k.has_value(),
                    "Gate-3A robust witness failed material-balanced seed construction");
            if (entry.kind == fl::Sw92AsymmetricPairPlanKind::reference_with_witness) {
                primary = true;
                require(entry.family_pair.phase0 == th::SwPhaseFamily::aqueous &&
                            entry.family_pair.phase1 == witness.family,
                        "primary witness family pair changed");
            }
            if (entry.kind ==
                    fl::Sw92AsymmetricPairPlanKind::witness_same_family_alternative) {
                same_family = same_family ||
                    (entry.family_pair.phase0 == witness.family &&
                     entry.family_pair.phase1 == witness.family);
            }
        }
        require(primary, "negative witness lost its reference-family primary attempt");
        if (witness.family == th::SwPhaseFamily::nonaqueous) {
            require(same_family, "opposite-family witness lost NA+NA alternative");
        }
    }

    std::reverse(initial.negative_witnesses.begin(), initial.negative_witnesses.end());
    const auto reversed_plan = fl::build_sw92_asymmetric_pair_plan(initial);
    require(reversed_plan.size() == plan.size(),
            "diagnostic witness storage order changed pair-plan size");
    for (std::size_t i = 0; i < plan.size(); ++i) {
        require(plan[i].witness_family == reversed_plan[i].witness_family &&
                    plan[i].witness_trial_index == reversed_plan[i].witness_trial_index &&
                    plan[i].kind == reversed_plan[i].kind &&
                    same_pair(plan[i].family_pair, reversed_plan[i].family_pair) &&
                    plan[i].protected_family_representative ==
                        reversed_plan[i].protected_family_representative &&
                    plan[i].initial_log_k == reversed_plan[i].initial_log_k,
                "pair planning depends on Gate-3A witness collection order");
    }
}

void full_candidate_selection() {
    const auto model = binary_model();
    const auto result = fl::orchestrate_sw92_asymmetric_pair_candidates(
        3.0e6, 340.0, Vec{0.7, 0.3}, model, 0.0);
    require(result.status == fl::Sw92AsymmetricPairSelectionStatus::
                                 pair_candidate_selected_pending_final_stability &&
                !result.attempt_limit_reached && result.selected_candidate_class &&
                !result.candidate_classes.empty(),
            "Gate 3B.2 did not select a complete unique pair candidate");
    const auto* pair = result.selected_pair_candidate_pending_final_stability();
    require(pair != nullptr, "selected pending-final-stability candidate is unavailable");
    require(pair->phase0.family == th::SwPhaseFamily::aqueous &&
                pair->phase1.family == th::SwPhaseFamily::aqueous,
            "candidate selection did not retain the lower-envelope AQ+AQ assignment");
    const double low_gas = std::min(pair->phase0.composition[0], pair->phase1.composition[0]);
    const double high_gas = std::max(pair->phase0.composition[0], pair->phase1.composition[0]);
    near(low_gas, aq_aq_phase0_gas);
    near(high_gas, aq_aq_phase1_gas);
    near(pair->reduced_gibbs, aq_aq_reduced_gibbs);
    require(pair->chemical_potential_norm <= result.options.fixed_pair.chemical_potential_tolerance &&
                pair->mass_absolute <= result.options.fixed_pair.mass_absolute_tolerance &&
                pair->mass_relative <= result.options.fixed_pair.mass_relative_tolerance,
            "selected pair lost fixed-pair equation gates");
    require(result.orchestration_convention ==
                fl::sw92_xu_asymmetric_orchestration_convention &&
                result.equilibrium_profile == fl::sw92_xu_asymmetric_gibbs_profile &&
                !result.global_stability_proven,
            "Gate 3B.2 metadata/global-proof semantics changed");
}

void family_neutral_attempt_budget() {
    const auto model = binary_model();
    fl::Sw92AsymmetricPairSelectionOptions options;
    options.max_pair_attempts = 2;
    const auto result = fl::orchestrate_sw92_asymmetric_pair_candidates(
        3.0e6, 340.0, Vec{0.7, 0.3}, model, 0.0, options);
    require(result.status ==
                fl::Sw92AsymmetricPairSelectionStatus::attempt_limit_reached &&
                result.attempt_limit_reached && result.attempts.size() == 2 &&
                result.plan.size() > result.attempts.size() &&
                !result.selected_candidate_class &&
                result.selected_pair_candidate_pending_final_stability() == nullptr,
            "truncated pair plan incorrectly published a selected candidate");
    require(result.plan[0].protected_family_representative &&
                result.plan[1].protected_family_representative &&
                result.plan[0].witness_family != result.plan[1].witness_family &&
                result.attempts[0].plan_index == 0 && result.attempts[1].plan_index == 1,
            "two-slot budget did not represent the best witness from each family first");
}

void stable_single_phase_candidate() {
    const auto model = binary_model();
    const auto result = fl::orchestrate_sw92_asymmetric_pair_candidates(
        3.0e6, 340.0, Vec{0.001, 0.999}, model, 0.0);
    require(result.status == fl::Sw92AsymmetricPairSelectionStatus::
                                 single_phase_candidate_no_instability_found &&
                result.single_phase_candidate && result.plan.empty() &&
                result.attempts.empty() && !result.selected_candidate_class,
            "Gate-3A stable feed did not remain a family-aware one-phase candidate");
    require(result.initial_stability.reference_family ==
                result.single_phase_candidate->family &&
                result.single_phase_candidate->composition == result.feed &&
                !result.global_stability_proven,
            "single-phase candidate lost lower-family/global-proof metadata");
}

void initial_family_tie_is_indeterminate() {
    const auto model = family_tie_model();
    const auto result = fl::orchestrate_sw92_asymmetric_pair_candidates(
        3.0e6, 340.0, Vec{0.7, 0.3}, model, 0.0);
    require(result.status ==
                fl::Sw92AsymmetricPairSelectionStatus::initial_stability_indeterminate &&
                result.initial_stability.feed_reference_status ==
                    fl::Sw92AsymmetricFeedReferenceStatus::family_gibbs_tie &&
                result.plan.empty() && result.attempts.empty() &&
                !result.single_phase_candidate && !result.selected_candidate_class,
            "family-envelope tie was broken by Gate 3B.2 orchestration");
}

fl::Sw92AsymmetricFixedPairResult synthetic_candidate(
    th::SwPhaseFamily family0, th::SwPhaseFamily family1,
    double fraction0, Vec composition0, Vec composition1,
    double gibbs, double guard) {
    fl::Sw92AsymmetricFixedPairResult result;
    result.status = fl::Sw92AsymmetricFixedPairStatus::converged;
    result.family_pair = {family0, family1};
    fl::Sw92AsymmetricFixedPairState point;
    point.phase0.family = family0;
    point.phase1.family = family1;
    point.phase0.mole_phase_fraction = fraction0;
    point.phase1.mole_phase_fraction = 1.0 - fraction0;
    point.phase0.composition = std::move(composition0);
    point.phase1.composition = std::move(composition1);
    point.chemical_potential_norm = 0.0;
    point.mass_absolute = 0.0;
    point.mass_relative = 0.0;
    point.reduced_gibbs = gibbs;
    point.gibbs_roundoff_guard = guard;
    point.phase0_assignment.status =
        fl::Sw92AsymmetricFamilyAssignmentStatus::assigned_lower;
    point.phase1_assignment.status =
        fl::Sw92AsymmetricFamilyAssignmentStatus::assigned_lower;
    result.point = std::move(point);
    return result;
}

void slot_swap_deduplication() {
    const Vec feed{0.52, 0.48};
    std::vector<fl::Sw92AsymmetricPairAttemptRecord> attempts;
    fl::Sw92AsymmetricPairAttemptRecord first;
    first.plan_index = 0;
    first.fixed_pair = synthetic_candidate(
        th::SwPhaseFamily::aqueous, th::SwPhaseFamily::nonaqueous,
        0.4, Vec{0.1, 0.9}, Vec{0.8, 0.2}, -1.2, 1e-13);
    attempts.push_back(std::move(first));
    fl::Sw92AsymmetricPairAttemptRecord swapped;
    swapped.plan_index = 1;
    swapped.fixed_pair = synthetic_candidate(
        th::SwPhaseFamily::nonaqueous, th::SwPhaseFamily::aqueous,
        0.6, Vec{0.8, 0.2}, Vec{0.1, 0.9}, -1.2, 1e-13);
    attempts.push_back(std::move(swapped));

    fl::Sw92AsymmetricFixedPairOptions options;
    const auto classes = fl::detail::sw92_build_candidate_classes(
        attempts, feed, options);
    require(classes.size() == 1 && classes[0].equivalent_attempts.size() == 2,
            "slot-swapped representation was not deduplicated");
    const auto choice = fl::detail::sw92_choose_candidate_class(classes);
    require(choice.selected_class == 0 && !choice.distinct_gibbs_tie,
            "one deduplicated candidate class was not uniquely selectable");
}

void distinct_candidate_gibbs_tie() {
    const Vec feed{0.5, 0.5};
    std::vector<fl::Sw92AsymmetricPairAttemptRecord> attempts;
    fl::Sw92AsymmetricPairAttemptRecord first;
    first.plan_index = 0;
    first.fixed_pair = synthetic_candidate(
        th::SwPhaseFamily::aqueous, th::SwPhaseFamily::aqueous,
        0.5, Vec{0.1, 0.9}, Vec{0.9, 0.1}, -2.0, 1e-12);
    attempts.push_back(std::move(first));
    fl::Sw92AsymmetricPairAttemptRecord second;
    second.plan_index = 1;
    second.fixed_pair = synthetic_candidate(
        th::SwPhaseFamily::aqueous, th::SwPhaseFamily::nonaqueous,
        0.5, Vec{0.2, 0.8}, Vec{0.8, 0.2}, -2.0 + 5e-13, 1e-12);
    attempts.push_back(std::move(second));

    fl::Sw92AsymmetricFixedPairOptions options;
    const auto classes = fl::detail::sw92_build_candidate_classes(
        attempts, feed, options);
    require(classes.size() == 2,
            "structurally distinct synthetic candidates were incorrectly deduplicated");
    const auto choice = fl::detail::sw92_choose_candidate_class(classes);
    require(choice.distinct_gibbs_tie && !choice.selected_class,
            "distinct Gibbs near-tie was resolved by attempt order");
}

void permutation() {
    const auto model = binary_model(true);
    const auto result = fl::orchestrate_sw92_asymmetric_pair_candidates(
        3.0e6, 340.0, Vec{0.3, 0.7}, model, 0.0);
    require(result.status == fl::Sw92AsymmetricPairSelectionStatus::
                                 pair_candidate_selected_pending_final_stability,
            "component permutation changed Gate 3B.2 selection status");
    const auto* pair = result.selected_pair_candidate_pending_final_stability();
    require(pair && pair->phase0.family == th::SwPhaseFamily::aqueous &&
                pair->phase1.family == th::SwPhaseFamily::aqueous,
            "permutation changed selected family assignment");
    const double low_gas = std::min(pair->phase0.composition[1], pair->phase1.composition[1]);
    const double high_gas = std::max(pair->phase0.composition[1], pair->phase1.composition[1]);
    near(low_gas, aq_aq_phase0_gas);
    near(high_gas, aq_aq_phase1_gas);
    near(pair->reduced_gibbs, aq_aq_reduced_gibbs);
    require(result.component_ids.size() == 2 && result.component_ids[0] == "water" &&
                result.component_ids[1] == "carbon-dioxide",
            "permuted ordered-snapshot metadata changed");
}

void resource_and_contract_failures() {
    const auto model = binary_model();
    fl::Sw92AsymmetricPairSelectionOptions invalid;
    invalid.max_pair_attempts = 0;
    expect_error<std::invalid_argument>([&] {
        (void)fl::orchestrate_sw92_asymmetric_pair_candidates(
            3.0e6, 340.0, Vec{0.7, 0.3}, model, 0.0, invalid);
    });

    fl::Sw92AsymmetricPairSelectionOptions tiny;
    tiny.fixed_pair.max_evaluations = 1;
    const auto limited = fl::orchestrate_sw92_asymmetric_pair_candidates(
        3.0e6, 340.0, Vec{0.7, 0.3}, model, 0.0, tiny);
    require(limited.status ==
                fl::Sw92AsymmetricPairSelectionStatus::pair_search_indeterminate &&
                !limited.attempt_limit_reached && !limited.selected_candidate_class &&
                limited.candidate_classes.empty() && !limited.attempts.empty(),
            "per-attempt property exhaustion was incorrectly converted to a candidate");
    bool saw_limit = false;
    for (const auto& attempt : limited.attempts) {
        if (attempt.fixed_pair &&
            attempt.fixed_pair->status ==
                fl::Sw92AsymmetricFixedPairStatus::evaluation_limit) {
            saw_limit = true;
        }
    }
    require(saw_limit, "fixed-pair evaluation-limit attribution was lost");

    fl::Sw92AsymmetricPairSelectionOptions initial_limited;
    initial_limited.initial_stability.aqueous.stability.max_iterations = 0;
    const auto unresolved = fl::orchestrate_sw92_asymmetric_pair_candidates(
        3.0e6, 340.0, Vec{0.001, 0.999}, model, 0.0, initial_limited);
    require(unresolved.status ==
                fl::Sw92AsymmetricPairSelectionStatus::initial_stability_indeterminate &&
                unresolved.plan.empty() && unresolved.attempts.empty(),
            "indeterminate initial family search incorrectly entered pair planning");
}

void headers() {
    require(sw92_asymmetric_orchestration_header(),
            "public orchestration header self-containment probe failed");
}

struct Case { const char* name; void (*function)(); };
constexpr Case cases[] = {
    {"plan_and_same_family_alternatives", plan_and_same_family_alternatives},
    {"full_candidate_selection", full_candidate_selection},
    {"family_neutral_attempt_budget", family_neutral_attempt_budget},
    {"stable_single_phase_candidate", stable_single_phase_candidate},
    {"initial_family_tie_is_indeterminate", initial_family_tie_is_indeterminate},
    {"slot_swap_deduplication", slot_swap_deduplication},
    {"distinct_candidate_gibbs_tie", distinct_candidate_gibbs_tie},
    {"permutation", permutation},
    {"resource_and_contract_failures", resource_and_contract_failures},
    {"headers", headers}
};

} // namespace

int main(int argc, char** argv) {
    if (argc != 2) {
        std::cerr << "usage: sw92_asymmetric_orchestration_test <case>\n";
        return 2;
    }
    for (const auto& test : cases) {
        if (test.name == std::string_view{argv[1]}) {
            try {
                test.function();
                std::cout << "[PASS] " << test.name << '\n';
                return 0;
            } catch (const std::exception& error) {
                std::cerr << "[FAIL] " << error.what() << '\n';
                return 1;
            }
        }
    }
    std::cerr << "unknown case: " << argv[1] << '\n';
    return 2;
}
