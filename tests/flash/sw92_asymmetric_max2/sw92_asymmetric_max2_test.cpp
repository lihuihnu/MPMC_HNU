#include <mpmc/flash/sw92_asymmetric_max2.hpp>

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

bool sw92_asymmetric_max2_header();

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

struct FinalGolden {
    long double lower_feed_gibbs;
    long double pair_gibbs;
    long double pair_minus_feed;
    long double common_d_gas;
    long double common_d_water;
    long double low_gas_na_tpd;
    long double high_gas_na_tpd;
    long double uniform_aq_tpd;
    long double uniform_na_tpd;
};

constexpr FinalGolden golden{
    -0.79037918785448424481162736483813301309149968040969608311561777113515725032080097L,
    -1.4960257659904066980637595884125953243224601687110748145260673572334134725858148L,
    -0.70564657813592245325213222357446231123096048830137873141044958609825622226501383L,
    -0.12302257740493822063922765644537909570586093753424415088069809116714384065361992L,
    -4.6996998726898331453876674296694331910945250414570130275768728644073582930088114L,
    0.017539861301589877725464986104432491292966300476862357683794437070700351502828894L,
    0.0014660970760666874933269379425177054878251130493963134442422400139919360861525528L,
    0.74779638225336205859463250543935025737418377680706966192503210804882551013601405L,
    1.2754701386896666255158362239238109660572907095872400164371934615459586170554450L
};

constexpr long double low_gas_reference =
    0.005936172024213170056428932801071659387598344017016019414L;
constexpr long double high_gas_reference =
    0.9874123977699055456489128208133867393111524204260154840L;

const fl::Sw92AsymmetricFixedPairState& selected_pair(
    const fl::Sw92AsymmetricMax2Result& result) {
    const auto* pair = result.selection.selected_pair_candidate_pending_final_stability();
    require(pair != nullptr, "selected Gate-3B.2 pair is unavailable");
    return *pair;
}

bool has_start(const fl::StabilityResult& search, const Vec& target) {
    return std::any_of(search.trials.begin(), search.trials.end(),
                       [&](const fl::StabilityTrial& trial) {
                           return trial.initial_composition == target;
                       });
}

double prescribed_tpd(double pressure_pa, double temperature_k,
                      const Vec& feed, const Vec& trial_composition,
                      const Vec& common_log_activity,
                      const th::Sw92Phase<double>& model,
                      th::SwPhaseFamily family) {
    fl::Sw92FamilyStabilityEvaluator evaluator(model, 0.0, family);
    const auto trial = evaluator(pressure_pa, temperature_k, trial_composition);
    fl::StabilityPhase equivalent_reference;
    equivalent_reference.ln_phi.resize(feed.size());
    equivalent_reference.smooth = true;
    for (std::size_t i = 0; i < feed.size(); ++i) {
        require(feed[i] > 0.0, "prescribed TPD fixture expects positive support");
        equivalent_reference.ln_phi[i] =
            common_log_activity[i] - std::log(feed[i]);
    }
    return fl::tangent_plane_distance(
        trial_composition, feed, trial, equivalent_reference).value;
}

void accepted_two_phase_reference() {
    const auto model = binary_model();
    const auto result = fl::solve_sw92_xu_asymmetric_max2(
        3.0e6, 340.0, Vec{0.7, 0.3}, model, 0.0);

    require(result.status ==
                fl::Sw92AsymmetricMax2Status::two_phase_no_instability_found &&
                result.accepted_phase_count() == 2 &&
                result.accepted_phase_set() != nullptr &&
                result.final_stability &&
                result.final_stability->status == fl::StabilityStatus::no_instability_found,
            "reference state was not accepted as a finite-search two-phase max2 state");
    require(result.selection.status ==
                fl::Sw92AsymmetricPairSelectionStatus::
                    pair_candidate_selected_pending_final_stability,
            "Gate-3B.2 selection evidence was not retained");
    require(result.max2_convention == fl::sw92_xu_asymmetric_max2_convention &&
                result.equilibrium_profile == fl::sw92_xu_asymmetric_gibbs_profile &&
                result.model_profile == th::sw92_corrected_profile &&
                result.phase_convention == th::sw92_pt_convention &&
                !result.global_stability_proven,
            "Gate-3B.3 identity/global-proof semantics changed");

    const auto& pair = selected_pair(result);
    require(pair.phase0.family == th::SwPhaseFamily::aqueous &&
                pair.phase1.family == th::SwPhaseFamily::aqueous,
            "reference selected pair lost its AQ+AQ family identity");
    const double low_gas = std::min(pair.phase0.composition[0], pair.phase1.composition[0]);
    const double high_gas = std::max(pair.phase0.composition[0], pair.phase1.composition[0]);
    near(low_gas, low_gas_reference);
    near(high_gas, high_gas_reference);
    near(pair.reduced_gibbs, golden.pair_gibbs);
    near(result.lower_feed_reduced_gibbs, golden.lower_feed_gibbs);
    near(result.pair_minus_lower_feed_reduced_gibbs, golden.pair_minus_feed);
    require(result.pair_minus_lower_feed_reduced_gibbs <
                -result.pair_feed_gibbs_combined_guard,
            "reference pair is not resolved below the lower-envelope feed Gibbs");

    const auto* accepted = result.accepted_phase_set();
    require(accepted->phases[0].compressibility_factor.has_value() &&
                accepted->phases[1].compressibility_factor.has_value() &&
                accepted->phases[0].family == th::SwPhaseFamily::aqueous &&
                accepted->phases[1].family == th::SwPhaseFamily::aqueous,
            "accepted family-aware pair lost Z/family metadata");
}

void final_common_tangent_and_required_starts() {
    const auto model = binary_model();
    const auto result = fl::solve_sw92_xu_asymmetric_max2(
        3.0e6, 340.0, Vec{0.7, 0.3}, model, 0.0);
    require(result.status ==
                fl::Sw92AsymmetricMax2Status::two_phase_no_instability_found &&
                result.final_stability,
            "reference final stability unavailable");
    const auto& pair = selected_pair(result);
    require(result.final_stability->common_log_activity == pair.common_log_activity,
            "final AQ/NA searches did not use the selected-pair midpoint tangent");
    near(pair.common_log_activity[0], golden.common_d_gas);
    near(pair.common_log_activity[1], golden.common_d_water);
    near(result.common_reference_allowance,
         0.5L * static_cast<long double>(pair.chemical_potential_norm),
         1e-12L, 1e-16L);
    near(result.final_aqueous_effective_tpd_tolerance -
             result.final_aqueous_base_tpd_tolerance,
         result.common_reference_allowance, 1e-12L, 1e-16L);
    near(result.final_nonaqueous_effective_tpd_tolerance -
             result.final_nonaqueous_base_tpd_tolerance,
         result.common_reference_allowance, 1e-12L, 1e-16L);

    require(has_start(result.final_stability->aqueous, pair.phase0.composition) &&
                has_start(result.final_stability->aqueous, pair.phase1.composition) &&
                has_start(result.final_stability->nonaqueous, pair.phase0.composition) &&
                has_start(result.final_stability->nonaqueous, pair.phase1.composition),
            "both candidate compositions were not mandatory starts in both final family searches");

    const Vec feed{0.7, 0.3};
    const Vec uniform{0.5, 0.5};
    const auto& low = pair.phase0.composition[0] < pair.phase1.composition[0]
        ? pair.phase0.composition : pair.phase1.composition;
    const auto& high = pair.phase0.composition[0] < pair.phase1.composition[0]
        ? pair.phase1.composition : pair.phase0.composition;
    require(std::abs(prescribed_tpd(3.0e6, 340.0, feed, low,
                                    pair.common_log_activity, model,
                                    th::SwPhaseFamily::aqueous)) < 2e-10 &&
                std::abs(prescribed_tpd(3.0e6, 340.0, feed, high,
                                        pair.common_log_activity, model,
                                        th::SwPhaseFamily::aqueous)) < 2e-10,
            "accepted AQ phases are not on their recorded common tangent");
    near(prescribed_tpd(3.0e6, 340.0, feed, low,
                        pair.common_log_activity, model,
                        th::SwPhaseFamily::nonaqueous),
         golden.low_gas_na_tpd);
    near(prescribed_tpd(3.0e6, 340.0, feed, high,
                        pair.common_log_activity, model,
                        th::SwPhaseFamily::nonaqueous),
         golden.high_gas_na_tpd);
    near(prescribed_tpd(3.0e6, 340.0, feed, uniform,
                        pair.common_log_activity, model,
                        th::SwPhaseFamily::aqueous),
         golden.uniform_aq_tpd);
    near(prescribed_tpd(3.0e6, 340.0, feed, uniform,
                        pair.common_log_activity, model,
                        th::SwPhaseFamily::nonaqueous),
         golden.uniform_na_tpd);
}

void accepted_single_phase() {
    const auto model = binary_model();
    const auto result = fl::solve_sw92_xu_asymmetric_max2(
        3.0e6, 340.0, Vec{0.001, 0.999}, model, 0.0);
    require(result.status ==
                fl::Sw92AsymmetricMax2Status::single_phase_no_instability_found &&
                result.selection.status ==
                    fl::Sw92AsymmetricPairSelectionStatus::
                        single_phase_candidate_no_instability_found &&
                result.accepted_phase_count() == 1 &&
                result.accepted_phase_set() != nullptr &&
                !result.final_stability && !result.global_stability_proven,
            "Gate-3A stable feed was not authoritatively published as one finite-search phase");
    const auto& phase = result.accepted_phase_set()->phases.front();
    require(result.selection.initial_stability.reference_family == phase.family &&
                phase.composition == result.feed &&
                phase.mole_phase_fraction == 1.0,
            "single-phase publication lost lower-family/feed identity");
}

void selection_indeterminate_not_published() {
    const auto model = binary_model();
    fl::Sw92AsymmetricMax2Options options;
    options.selection.max_pair_attempts = 2;
    const auto result = fl::solve_sw92_xu_asymmetric_max2(
        3.0e6, 340.0, Vec{0.7, 0.3}, model, 0.0, options);
    require(result.status == fl::Sw92AsymmetricMax2Status::indeterminate &&
                result.selection.status ==
                    fl::Sw92AsymmetricPairSelectionStatus::attempt_limit_reached &&
                result.accepted_phase_set() == nullptr &&
                !result.final_stability,
            "truncated Gate-3B.2 plan incorrectly entered final acceptance");
}

void final_indeterminate_not_published() {
    const auto model = binary_model();
    fl::Sw92AsymmetricMax2Options options;
    options.final_stability.aqueous.stability.max_iterations = 0;
    options.final_stability.nonaqueous.stability.max_iterations = 0;
    const auto result = fl::solve_sw92_xu_asymmetric_max2(
        3.0e6, 340.0, Vec{0.7, 0.3}, model, 0.0, options);
    require(result.status == fl::Sw92AsymmetricMax2Status::indeterminate &&
                result.candidate_phase_set && result.final_stability &&
                result.final_stability->status == fl::StabilityStatus::indeterminate &&
                result.accepted_phase_set() == nullptr,
            "indeterminate final family search incorrectly published the selected pair");
}

void pair_gibbs_above_feed_not_published() {
    const auto model = binary_model();
    auto selection = fl::orchestrate_sw92_asymmetric_pair_candidates(
        3.0e6, 340.0, Vec{0.7, 0.3}, model, 0.0);
    require(selection.status ==
                fl::Sw92AsymmetricPairSelectionStatus::
                    pair_candidate_selected_pending_final_stability &&
                selection.selected_candidate_class,
            "reference Gate-3B.2 selection unavailable for Gibbs-gate structural test");
    const auto class_index = *selection.selected_candidate_class;
    const auto attempt_index =
        selection.candidate_classes[class_index].representative_attempt;
    require(attempt_index < selection.attempts.size() &&
                selection.attempts[attempt_index].fixed_pair &&
                selection.attempts[attempt_index].fixed_pair->point,
            "reference pair attempt unavailable");

    // Synthetic post-selection perturbation: test only the Gate-3B.3 publication
    // guard, not a physical SW92 state or a new thermodynamic reference.
    auto& point = *selection.attempts[attempt_index].fixed_pair->point;
    const auto family = *selection.initial_stability.reference_family;
    const auto& feed_leg = family == th::SwPhaseFamily::aqueous
        ? selection.initial_stability.aqueous : selection.initial_stability.nonaqueous;
    point.reduced_gibbs = feed_leg.feed_reduced_gibbs + 0.25;

    const auto result = fl::finalize_sw92_asymmetric_max2_selection(
        std::move(selection), model);
    require(result.status == fl::Sw92AsymmetricMax2Status::pair_gibbs_above_feed &&
                result.candidate_phase_set && !result.final_stability &&
                result.accepted_phase_set() == nullptr &&
                result.pair_minus_lower_feed_reduced_gibbs >
                    result.pair_feed_gibbs_combined_guard,
            "resolved pair-above-feed Gibbs state was incorrectly accepted");
}

void publication_structural_guard() {
    const auto model = binary_model();
    auto result = fl::solve_sw92_xu_asymmetric_max2(
        3.0e6, 340.0, Vec{0.7, 0.3}, model, 0.0);
    require(result.accepted_phase_set() != nullptr && result.candidate_phase_set,
            "reference accepted phase set unavailable");
    result.candidate_phase_set->phases[0].mole_phase_fraction = 0.0;
    require(result.accepted_phase_set() == nullptr && result.accepted_phase_count() == 0,
            "accepted-phase publication helper exposed a structurally malformed phase set");
}

void component_permutation() {
    const auto model = binary_model(true);
    const auto result = fl::solve_sw92_xu_asymmetric_max2(
        3.0e6, 340.0, Vec{0.3, 0.7}, model, 0.0);
    require(result.status ==
                fl::Sw92AsymmetricMax2Status::two_phase_no_instability_found &&
                result.accepted_phase_count() == 2 && result.final_stability &&
                result.component_ids.size() == 2 &&
                result.component_ids[0] == "water" &&
                result.component_ids[1] == "carbon-dioxide",
            "component permutation changed final max2 status/metadata");
    const auto& pair = selected_pair(result);
    const double low_gas = std::min(pair.phase0.composition[1], pair.phase1.composition[1]);
    const double high_gas = std::max(pair.phase0.composition[1], pair.phase1.composition[1]);
    near(low_gas, low_gas_reference);
    near(high_gas, high_gas_reference);
    near(result.pair_minus_lower_feed_reduced_gibbs, golden.pair_minus_feed);
}

void final_start_quota_preflight() {
    const auto model = binary_model();
    fl::Sw92AsymmetricMax2Options options;
    // Binary automatic search generates four starts; Gate 3B.3 must append two
    // candidate phases to each family, so five is intentionally insufficient.
    options.final_stability.aqueous.stability.max_starts = 5;
    options.final_stability.nonaqueous.stability.max_starts = 5;
    expect_error<std::length_error>([&] {
        (void)fl::solve_sw92_xu_asymmetric_max2(
            3.0e6, 340.0, Vec{0.7, 0.3}, model, 0.0, options);
    });
}

void headers() {
    require(sw92_asymmetric_max2_header(), "public max2 header probe failed");
}

using Test = void (*)();
struct TestCase { std::string_view name; Test test; };

constexpr TestCase tests[] = {
    {"accepted_two_phase_reference", accepted_two_phase_reference},
    {"final_common_tangent_and_required_starts", final_common_tangent_and_required_starts},
    {"accepted_single_phase", accepted_single_phase},
    {"selection_indeterminate_not_published", selection_indeterminate_not_published},
    {"final_indeterminate_not_published", final_indeterminate_not_published},
    {"pair_gibbs_above_feed_not_published", pair_gibbs_above_feed_not_published},
    {"publication_structural_guard", publication_structural_guard},
    {"component_permutation", component_permutation},
    {"final_start_quota_preflight", final_start_quota_preflight},
    {"headers", headers},
};

} // namespace

int main(int argc, char** argv) {
    if (argc != 2) {
        std::cerr << "usage: sw92_asymmetric_max2_test <case>\n";
        return 2;
    }
    const std::string_view requested{argv[1]};
    for (const auto& test : tests) {
        if (test.name == requested) {
            try {
                test.test();
                return 0;
            } catch (const std::exception& error) {
                std::cerr << requested << ": " << error.what() << '\n';
                return 1;
            }
        }
    }
    std::cerr << "unknown case: " << requested << '\n';
    return 2;
}
