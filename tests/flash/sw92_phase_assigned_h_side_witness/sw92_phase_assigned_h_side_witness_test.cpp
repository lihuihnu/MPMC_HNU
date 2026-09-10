#include <mpmc/flash/sw92_phase_assigned_h_side_witness.hpp>

#include "test_support.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <iostream>
#include <limits>
#include <source_location>
#include <stdexcept>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

bool sw92_phase_assigned_h_side_witness_header();

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

void near(double actual, long double expected, long double relative = 2e-7L,
          long double absolute = 2e-10L,
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

struct ReferenceCase {
    const sw92_test::Spec* gas;
    double pressure_pa;
    double temperature_k;
    double molality;
    double feed_gas;
    Vec seed;
};

ReferenceCase reference_case(std::size_t index) {
    if (index == 0U) {
        return {&sw92_test::co2, 3.0e6, 340.0, 0.0, 0.7, {5.0, -4.4}};
    }
    return {&sw92_test::methane, 1.0e7, 350.0, 1.0, 0.5, {6.5, -5.0}};
}

th::Sw92Phase<double> model_for(const ReferenceCase& reference, bool reverse = false) {
    return th::Sw92Phase<double>::from_parameters(
        sw92_test::binary_parameters(*reference.gas, reverse));
}

fl::Sw92PhaseAssignedJointResult solve_c1(
    std::size_t index, bool reverse = false,
    fl::Sw92PhaseAssignedJointOptions options = {}) {
    const auto reference = reference_case(index);
    const auto model = model_for(reference, reverse);
    Vec feed = reverse
        ? Vec{1.0 - reference.feed_gas, reference.feed_gas}
        : Vec{reference.feed_gas, 1.0 - reference.feed_gas};
    Vec seed = reference.seed;
    if (reverse) { std::swap(seed[0], seed[1]); }
    return fl::iterate_sw92_phase_assigned_aq_na_joint(
        reference.pressure_pa, reference.temperature_k, feed, seed,
        model, reference.molality, options);
}

th::Sw92Phase<double> synthetic_na_zero_model(bool reverse = false) {
    auto prepared = sw92_test::binary_input(sw92_test::co2, reverse);
    prepared.input.dataset_id = "SW92-C2a1-synthetic-NA-zero";
    prepared.input.revision = "synthetic-test-v1";
    auto source = sw92_test::synthetic(
        "C2a1 structural fixture: CO2/water nonaqueous kij=0");
    prepared.input.water_binary.at(0).nonaqueous_kij =
        sw92_test::scalar(0.0, th::Unit::dimensionless, source);
    return th::Sw92Phase<double>::from_parameters(
        th::Sw92ParameterSet::create(
            prepared.catalog, prepared.order, prepared.input,
            th::DataPolicy::allow_synthetic_tests));
}

fl::Sw92PhaseAssignedJointResult solve_synthetic_c1(bool reverse = false) {
    const auto model = synthetic_na_zero_model(reverse);
    Vec feed{0.5, 0.5};
    Vec seed{6.1739566237881105, -4.396681357954204};
    if (reverse) { std::swap(seed[0], seed[1]); }
    return fl::iterate_sw92_phase_assigned_aq_na_joint(
        3.0e6, 340.0, feed, seed, model, 0.0);
}

void check_no_witness(std::size_t index) {
    const auto reference = reference_case(index);
    const auto model = model_for(reference);
    const auto c1 = solve_c1(index);
    require(c1.candidate() != nullptr, "C1 reference candidate unavailable");
    const auto result = fl::test_sw92_phase_assigned_h_side_na_witness(c1, model);
    require(result.status == fl::Sw92PhaseAssignedHSideWitnessStatus::
                                 no_additional_nonaqueous_witness_found,
            "traceable C1 reference unexpectedly produced an H-split witness");
    require(result.nonaqueous_search.has_value(), "NA-only search result missing");
    require(result.nonaqueous_search->status == fl::StabilityStatus::no_instability_found,
            "finite NA search did not resolve to no-instability-found");
    require(result.negative_witnesses.empty() && !result.has_usable_h_split_witness(),
            "no-witness state retained a negative H-split seed");
    require(result.retained_h_trivial_point.has_value(),
            "retained H common-tangent diagnostic missing");
    require(std::abs(result.retained_h_trivial_point->value) <=
                result.common_reference_allowance +
                result.retained_h_trivial_point->roundoff_guard,
            "retained H is not on the reconstructed C1 tangent");
    require(!result.global_stability_proven && !result.accepted_phase_set_published,
            "C2a1 must not claim global stability or phase-set publication");
    require(result.witness_convention ==
                fl::sw92_phase_assigned_h_side_na_witness_convention &&
                result.equilibrium_profile == fl::sw92_phase_assigned_aq_na_joint_profile,
            "C2a1 algorithm identity changed");
}

void co2_no_additional_witness() {
    check_no_witness(0U);
}

void methane_brine_no_additional_witness() {
    check_no_witness(1U);
}

void cross_family_blocker_is_not_c2a1_rejection() {
    const auto reference = reference_case(0U);
    const auto model = model_for(reference);
    const auto c1 = solve_c1(0U);
    const auto result = fl::test_sw92_phase_assigned_h_side_na_witness(c1, model);
    require(result.status == fl::Sw92PhaseAssignedHSideWitnessStatus::
                                 no_additional_nonaqueous_witness_found,
            "C2a1 inherited cross-family rejection");

    const auto& c1_point = *c1.point;
    fl::StabilityPhase equivalent_reference;
    equivalent_reference.ln_phi.resize(c1.feed.size());
    for (std::size_t i = 0; i < c1.feed.size(); ++i) {
        equivalent_reference.ln_phi[i] = c1.feed[i] > 0.0
            ? result.common_log_activity[i] - std::log(c1.feed[i]) : 0.0;
    }
    fl::Sw92FamilyStabilityEvaluator aq_evaluator(
        model, 0.0, th::SwPhaseFamily::aqueous);
    const auto aq_at_h = aq_evaluator(
        reference.pressure_pa, reference.temperature_k,
        c1_point.nonaqueous_phase.composition);
    const auto forbidden_cross_family = fl::tangent_plane_distance(
        c1_point.nonaqueous_phase.composition, c1.feed,
        aq_at_h, equivalent_reference);
    require(forbidden_cross_family.value < -1e-4,
            "CO2 blocker reference lost negative AQ-at-H TPD");
    require(result.negative_witnesses.empty(),
            "AQ-at-H blocker leaked into NA-only witness collection");
}

void synthetic_admissible_witness() {
    const auto model = synthetic_na_zero_model();
    const auto c1 = solve_synthetic_c1();
    require(c1.candidate() != nullptr, "synthetic C1 candidate unavailable");

    fl::Sw92PhaseAssignedHSideWitnessOptions options;
    options.stability.automatic_starts = false;
    const std::vector<Vec> extra{{0.003, 0.997}};
    const auto result = fl::test_sw92_phase_assigned_h_side_na_witness(
        c1, model, options, extra);
    require(result.status == fl::Sw92PhaseAssignedHSideWitnessStatus::
                                 additional_nonaqueous_phase_witness_found,
            "synthetic distinct NA negative trial was not retained as H-split witness");
    require(result.has_usable_h_split_witness() && result.usable_witness_count() >= 1U,
            "usable synthetic H-split witness missing");

    const auto it = std::find_if(
        result.negative_witnesses.begin(), result.negative_witnesses.end(),
        [](const fl::Sw92PhaseAssignedNaNegativeWitness& witness) {
            return witness.usable_h_split_seed() &&
                   witness.point.composition.size() == 2U &&
                   std::abs(witness.point.composition[0] - 0.003) < 1e-12;
        });
    require(it != result.negative_witnesses.end(),
            "explicit 0.003 synthetic NA witness was not preserved");
    require(it->point.value < -1e-5,
            "synthetic H-split witness is not robustly negative");
    near(it->point.value, -7.643999071963782e-4L, 5e-4L, 2e-7L);
    require(it->compositionally_distinct_from_retained_h &&
                it->water_role_admissible &&
                it->retained_w_minus_trial_water_fraction >
                    it->water_role_roundoff_guard,
            "synthetic H-split witness failed required topology guards");
}

void synthetic_role_guard_rejects_w_like_negative() {
    const auto model = synthetic_na_zero_model();
    const auto c1 = solve_synthetic_c1();
    fl::Sw92PhaseAssignedHSideWitnessOptions options;
    options.stability.automatic_starts = false;
    const auto result = fl::test_sw92_phase_assigned_h_side_na_witness(
        c1, model, options);
    require(result.nonaqueous_search &&
                result.nonaqueous_search->status == fl::StabilityStatus::unstable,
            "synthetic NA-at-W negative trial was not found mathematically");
    require(result.status == fl::Sw92PhaseAssignedHSideWitnessStatus::indeterminate,
            "role-inadmissible negative trial was incorrectly promoted or called stable");
    require(!result.has_usable_h_split_witness() && !result.negative_witnesses.empty(),
            "role-inadmissible negative evidence was not retained diagnostically");
    const bool has_rejected = std::any_of(
        result.negative_witnesses.begin(), result.negative_witnesses.end(),
        [](const fl::Sw92PhaseAssignedNaNegativeWitness& witness) {
            return witness.compositionally_distinct_from_retained_h &&
                   !witness.water_role_admissible;
        });
    require(has_rejected,
            "negative W-like NA trial did not exercise the physical-role guard");
}

void resource_exhaustion_is_indeterminate() {
    const auto reference = reference_case(0U);
    const auto model = model_for(reference);
    const auto c1 = solve_c1(0U);
    fl::Sw92PhaseAssignedHSideWitnessOptions options;
    options.stability.automatic_starts = false;
    options.stability.max_evaluations = 1U;
    const auto result = fl::test_sw92_phase_assigned_h_side_na_witness(
        c1, model, options);
    require(result.status == fl::Sw92PhaseAssignedHSideWitnessStatus::indeterminate,
            "finite resource exhaustion was converted to no-witness");
    require(result.nonaqueous_search &&
                result.nonaqueous_search->status == fl::StabilityStatus::indeterminate,
            "NA resource exhaustion did not remain explicit");
    const bool limited = std::any_of(
        result.nonaqueous_search->trials.begin(), result.nonaqueous_search->trials.end(),
        [](const fl::StabilityTrial& trial) {
            return trial.status == fl::StabilityTrialStatus::evaluation_limit;
        });
    require(limited, "evaluation-limit trial diagnostic missing");
}

void root_failure_is_indeterminate() {
    const auto reference = reference_case(0U);
    const auto model = model_for(reference);
    const auto c1 = solve_c1(0U);
    fl::Sw92PhaseAssignedHSideWitnessOptions options;
    options.stability.automatic_starts = false;
    options.nonaqueous_root_options.max_iterations = 1;
    const auto result = fl::test_sw92_phase_assigned_h_side_na_witness(
        c1, model, options);
    require(result.status == fl::Sw92PhaseAssignedHSideWitnessStatus::indeterminate,
            "NA root/property failure was converted to no-witness");
    require(result.nonaqueous_search &&
                result.nonaqueous_search->status == fl::StabilityStatus::indeterminate,
            "NA root/property failure did not remain search-indeterminate");
    const bool failed = std::any_of(
        result.nonaqueous_search->trials.begin(), result.nonaqueous_search->trials.end(),
        [](const fl::StabilityTrial& trial) {
            return trial.status == fl::StabilityTrialStatus::property_failure &&
                   trial.property_issue.has_value();
        });
    require(failed, "family-specific root/property issue was not retained");
}

void source_candidate_guard() {
    const auto reference = reference_case(0U);
    const auto model = model_for(reference);
    fl::Sw92PhaseAssignedJointOptions limited;
    limited.max_evaluations = 1U;
    const auto source = solve_c1(0U, false, limited);
    require(source.candidate() == nullptr,
            "source fixture unexpectedly produced an admissible C1 candidate");
    const auto result = fl::test_sw92_phase_assigned_h_side_na_witness(source, model);
    require(result.status == fl::Sw92PhaseAssignedHSideWitnessStatus::
                                 source_candidate_unavailable &&
                !result.nonaqueous_search,
            "C2a1 searched without an admissible C1 source candidate");

    auto wrong_model_prepared = sw92_test::binary_input(sw92_test::co2);
    wrong_model_prepared.input.dataset_id = "different-dataset";
    const auto wrong_model = th::Sw92Phase<double>::from_parameters(
        th::Sw92ParameterSet::create(
            wrong_model_prepared.catalog, wrong_model_prepared.order,
            wrong_model_prepared.input));
    const auto good_source = solve_c1(0U);
    expect_error<std::invalid_argument>([&] {
        (void)fl::test_sw92_phase_assigned_h_side_na_witness(good_source, wrong_model);
    });
}

void component_permutation() {
    const auto reference = reference_case(0U);
    const auto normal_model = model_for(reference, false);
    const auto reverse_model = model_for(reference, true);
    const auto normal_source = solve_c1(0U, false);
    const auto reverse_source = solve_c1(0U, true);
    const auto normal = fl::test_sw92_phase_assigned_h_side_na_witness(
        normal_source, normal_model);
    const auto reverse = fl::test_sw92_phase_assigned_h_side_na_witness(
        reverse_source, reverse_model);
    require(normal.status == reverse.status &&
                normal.status == fl::Sw92PhaseAssignedHSideWitnessStatus::
                                     no_additional_nonaqueous_witness_found,
            "component permutation changed C2a1 topology status");
    require(normal.water_index == 1U && reverse.water_index == 0U,
            "water index did not follow ordered component permutation");
    near(normal.retained_w_water_fraction, reverse.retained_w_water_fraction);
    near(normal.retained_h_water_fraction, reverse.retained_h_water_fraction);
    near(normal.retained_h_trivial_point->value,
         reverse.retained_h_trivial_point->value, 1e-5L, 1e-12L);
}

void headers() {
    require(sw92_phase_assigned_h_side_witness_header(),
            "C2a1 public header probe failed");
}

} // namespace

int main(int argc, char** argv) {
    if (argc != 2) {
        std::cerr << "expected one test-case name\n";
        return 2;
    }
    const std::string_view name{argv[1]};
    try {
        if (name == "co2_no_witness") co2_no_additional_witness();
        else if (name == "methane_brine_no_witness") methane_brine_no_additional_witness();
        else if (name == "cross_family_separation") cross_family_blocker_is_not_c2a1_rejection();
        else if (name == "synthetic_admissible_witness") synthetic_admissible_witness();
        else if (name == "synthetic_role_guard") synthetic_role_guard_rejects_w_like_negative();
        else if (name == "resource_exhaustion") resource_exhaustion_is_indeterminate();
        else if (name == "root_failure") root_failure_is_indeterminate();
        else if (name == "source_guard") source_candidate_guard();
        else if (name == "permutation") component_permutation();
        else if (name == "headers") headers();
        else throw std::invalid_argument("unknown test-case name");
    } catch (const std::exception& error) {
        std::cerr << error.what() << '\n';
        return 1;
    }
    return 0;
}
