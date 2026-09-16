#include <mpmc/flash/pr76_max3_phase_set.hpp>
#include <mpmc/flash/pr76_pt_flash_backend.hpp>

#include "sour_gas_fixture.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <iostream>
#include <numeric>
#include <optional>
#include <span>
#include <stdexcept>
#include <string_view>
#include <vector>

namespace {
namespace fl = mpmc::flash;
namespace fx = pr76_sour_gas_test;
namespace ref = pr76_sour_gas_reference;
using Vec = std::vector<double>;

constexpr double reference_budget = 2.0e-6;
constexpr double arithmetic_budget = 5.0e-13;

void require(bool condition, const char* message) {
    if (!condition) { throw std::runtime_error(message); }
}

Vec log_ratio(std::span<const double> numerator,
              std::span<const double> denominator) {
    require(numerator.size() == denominator.size(), "log-ratio dimension mismatch");
    Vec result(numerator.size(), 0.0);
    for (std::size_t i = 0; i < result.size(); ++i) {
        require(numerator[i] > 0.0 && denominator[i] > 0.0,
                "log-ratio requires positive compositions");
        result[i] = std::log(numerator[i]) - std::log(denominator[i]);
    }
    return result;
}

void print_max3(const fl::Pr76PtMax3Result& result) {
    std::cerr << "max3_status=" << static_cast<int>(result.status)
              << " base_status=" << static_cast<int>(result.base.solution.status)
              << " base_attempts=" << result.base.solution.attempts.size()
              << " max3_attempts=" << result.attempts.size()
              << " selected=" << (result.selected_attempt ? 1 : 0)
              << " diagnostic=" << result.diagnostic << '\n';
    if (result.base.solution.final_stability) {
        std::cerr << "base_final_stability="
                  << static_cast<int>(result.base.solution.final_stability->status)
                  << " trials=" << result.base.solution.final_stability->trials.size()
                  << '\n';
    }
    for (std::size_t i = 0; i < result.attempts.size(); ++i) {
        const auto& attempt = result.attempts[i];
        std::cerr << "attempt[" << i << "]"
                  << " supplied="
                  << (attempt.supplied_start
                          ? static_cast<long long>(*attempt.supplied_start) : -1LL)
                  << " witness="
                  << (attempt.witness_trial
                          ? static_cast<long long>(*attempt.witness_trial) : -1LL)
                  << " eq_status=" << static_cast<int>(attempt.equilibrium.status)
                  << " eq_norm="
                  << (attempt.equilibrium.point
                          ? attempt.equilibrium.point->chemical_potential_norm : -1.0)
                  << " final="
                  << (attempt.final_stability
                          ? static_cast<int>(attempt.final_stability->status) : -1)
                  << " accepted3=" << attempt.accepted_three_phase
                  << " accepted2=" << attempt.accepted_two_phase_neighbor
                  << " diagnostic=" << attempt.equilibrium.diagnostic << '\n';
    }
}

struct Match {
    std::array<std::size_t, 3> reference_for_phase{};
};

Match match_reference(
    const fl::PtThreePhaseState& state,
    const std::vector<std::string>& order,
    const std::array<double, 3>& expected_fractions =
        ref::phase_fractions) {
    std::array<bool, 3> used{false, false, false};
    Match match;
    for (std::size_t phase_index = 0; phase_index < 3U; ++phase_index) {
        const auto canonical = fx::canonicalize(
            state.phases[phase_index].composition, order);
        std::optional<std::size_t> best;
        double best_error = 1.0e100;
        for (std::size_t candidate = 0; candidate < 3U; ++candidate) {
            if (used[candidate]) { continue; }
            double error = 0.0;
            for (std::size_t component = 0; component < canonical.size(); ++component) {
                error = std::max(
                    error,
                    std::abs(canonical[component] - ref::phases[candidate][component]));
            }
            if (error < best_error) {
                best_error = error;
                best = candidate;
            }
        }
        require(best.has_value() && best_error <= reference_budget,
                "PR76 sour-gas phase composition left independent reference budget");
        used[*best] = true;
        match.reference_for_phase[phase_index] = *best;
        require(
            std::abs(state.phases[phase_index].mole_phase_fraction -
                     expected_fractions[*best]) <= reference_budget,
            "PR76 sour-gas phase fraction left independent reference budget");
        require(
            std::abs(state.phases[phase_index].z -
                     ref::compressibility_factors[*best]) <= reference_budget,
            "PR76 sour-gas Z left independent reference budget");
    }
    return match;
}

void require_state_invariants(const fl::PtThreePhaseState& state) {
    require(state.chemical_potential_norm <= 1.0e-11,
            "three-phase chemical-potential residual exceeds production gate");
    require(state.generalized_rr_residual <= 2.0e-13,
            "three-phase generalized RR residual exceeds production gate");
    require(state.mass_absolute <= 1.0e-12 &&
                state.mass_relative <= 1.0e-10,
            "three-phase material balance exceeds production gate");
    double phase_sum = 0.0;
    for (const auto& phase : state.phases) {
        require(phase.mole_phase_fraction > 1.0e-10,
                "accepted benchmark phase is below production phase-fraction gate");
        phase_sum += phase.mole_phase_fraction;
        const double composition_sum =
            std::accumulate(phase.composition.begin(), phase.composition.end(), 0.0);
        require(std::abs(composition_sum - 1.0) <= arithmetic_budget,
                "three-phase composition is not normalized");
    }
    require(std::abs(phase_sum - 1.0) <= arithmetic_budget,
            "three-phase fractions do not sum to one");
}

void reference_candidate() {
    const auto order = fx::canonical_order();
    const auto model = fx::model(order);
    const auto starts = fx::starts(order);

    fl::Pr76ThreePhaseEvaluator evaluator(
        model,
        {fl::Pr76RootSide::upper_admissible,
         fl::Pr76RootSide::lower_admissible,
         fl::Pr76RootSide::lower_admissible});

    double fresh_mu_norm = 0.0;
    std::array<fl::PtThreePhaseProperty, 3> properties{
        evaluator(ref::pressure_pa, ref::temperature_k, starts[0], 0U),
        evaluator(ref::pressure_pa, ref::temperature_k, starts[1], 1U),
        evaluator(ref::pressure_pa, ref::temperature_k, starts[2], 2U)};
    for (std::size_t component = 0; component < ref::feed.size(); ++component) {
        const double mu0 =
            std::log(starts[0][component]) + properties[0].activity.ln_phi[component];
        const double mu1 =
            std::log(starts[1][component]) + properties[1].activity.ln_phi[component];
        const double mu2 =
            std::log(starts[2][component]) + properties[2].activity.ln_phi[component];
        fresh_mu_norm = std::max({
            fresh_mu_norm, std::abs(mu0 - mu1), std::abs(mu0 - mu2)});
    }
    require(fresh_mu_norm <= 5.0e-12,
            "independent sour-gas reference does not close under fresh PR76 properties");
    for (std::size_t phase = 0; phase < 3U; ++phase) {
        require(
            std::abs(properties[phase].z - ref::compressibility_factors[phase]) <=
                5.0e-12,
            "fresh PR76 property Z differs from Decimal reference");
    }

    const auto result = fl::iterate_pt_three_phase(
        ref::pressure_pa, ref::temperature_k, fx::feed(order),
        log_ratio(starts[1], starts[0]),
        log_ratio(starts[2], starts[0]),
        {ref::phase_fractions[1], ref::phase_fractions[2]},
        evaluator);
    if (result.status != fl::PtThreePhaseStatus::converged_candidate ||
        result.candidate() == nullptr) {
        std::cerr << "fixed_status=" << static_cast<int>(result.status)
                  << " iterations=" << result.iterations
                  << " evaluations=" << result.evaluations
                  << " diagnostic=" << result.diagnostic << '\n';
        throw std::runtime_error(
            "literature-defined sour-gas PR76 point did not close in fixed-three-phase primitive");
    }
    require_state_invariants(*result.candidate());
    (void)match_reference(*result.candidate(), order);

    // Audit the candidate-owned common tangent independently of max3 attempt
    // orchestration.  The three phase activities need not be bitwise equal at
    // the finite chemical-potential tolerance, so the arithmetic mean is the
    // imposed tangent and the candidate norm is a conservative componentwise
    // bound on every phase's deviation from that mean.
    const auto& candidate = *result.candidate();
    double reconstructed_norm = 0.0;
    double maximum_reference_deviation = 0.0;
    for (std::size_t component = 0; component < ref::feed.size(); ++component) {
        std::array<double, 3> mu{};
        for (std::size_t phase = 0; phase < 3U; ++phase) {
            mu[phase] = std::log(candidate.phases[phase].composition[component]) +
                        candidate.phases[phase].activity.ln_phi[component];
        }
        const double mean = (mu[0] + mu[1] + mu[2]) / 3.0;
        require(std::abs(candidate.common_log_activity[component] - mean) <=
                    arithmetic_budget,
                "three-phase common log activity is not the arithmetic mean");
        reconstructed_norm = std::max({
            reconstructed_norm, std::abs(mu[0] - mu[1]), std::abs(mu[0] - mu[2])});
        maximum_reference_deviation = std::max({
            maximum_reference_deviation,
            std::abs(mu[0] - mean), std::abs(mu[1] - mean), std::abs(mu[2] - mean)});
    }
    require(std::abs(reconstructed_norm - candidate.chemical_potential_norm) <=
                arithmetic_budget,
            "three-phase chemical-potential norm changed meaning");
    require(maximum_reference_deviation <= candidate.chemical_potential_norm +
                arithmetic_budget,
            "three-phase common-reference allowance underestimates phase mismatch");

    std::vector<Vec> phase_starts;
    phase_starts.reserve(3U);
    for (const auto& phase : candidate.phases) {
        phase_starts.push_back(phase.composition);
    }
    auto final_options = fl::StabilityOptions{};
    final_options.tpd_tolerance += candidate.chemical_potential_norm;

    fl::Pr76StabilityEvaluator final_evaluator(model);
    const auto stable_review = fl::test_pt_stability_against(
        ref::pressure_pa, ref::temperature_k, fx::feed(order),
        candidate.common_log_activity, final_evaluator,
        final_options, phase_starts);
    require(stable_review.status == fl::StabilityStatus::no_instability_found,
            "independent final common-tangent review rejected the sour-gas candidate");
    require(!stable_review.reference.has_value() &&
                stable_review.imposed_log_activity.size() == ref::feed.size(),
            "imposed three-phase tangent was replaced by a feed reference");
    require(stable_review.trials.size() == ref::feed.size() + 5U,
            "final three-phase review did not include automatic plus three phase starts");

    // Raising every imposed chemical potential by a resolvable constant shifts
    // every TPD value downward by the same constant.  The three phase starts
    // therefore provide explicit negative evidence that a wrong tangent cannot
    // be accepted merely because the fixed three-phase equations converged.
    Vec wrong_reference = candidate.common_log_activity;
    for (double& value : wrong_reference) { value += 1.0e-5; }
    fl::Pr76StabilityEvaluator wrong_evaluator(model);
    const auto wrong_review = fl::test_pt_stability_against(
        ref::pressure_pa, ref::temperature_k, fx::feed(order),
        wrong_reference, wrong_evaluator, final_options, phase_starts);
    require(wrong_review.status == fl::StabilityStatus::unstable &&
                wrong_review.lowest_sampled.has_value() &&
                wrong_review.lowest_sampled->value < -5.0e-6,
            "resolvably wrong common tangent did not produce negative TPD evidence");

    // Resource exhaustion in the final search must remain indeterminate.  It is
    // not permission to publish a three-phase state simply because the fixed
    // equilibrium equations already converged.
    auto limited_options = final_options;
    limited_options.max_evaluations = 1U;
    fl::Pr76StabilityEvaluator limited_evaluator(model);
    const auto limited_review = fl::test_pt_stability_against(
        ref::pressure_pa, ref::temperature_k, fx::feed(order),
        candidate.common_log_activity, limited_evaluator,
        limited_options, phase_starts);
    require(limited_review.status == fl::StabilityStatus::indeterminate,
            "final TPD resource exhaustion was promoted to a stability decision");
}

fl::Pr76PtMax3Result solve_max3(
    const std::vector<std::string>& order,
    const Vec& benchmark_feed,
    const std::array<double, 3>* expected_fractions = nullptr) {
    const auto model = fx::model(order);
    fl::Pr76VleEvaluator evaluator(model);
    const auto starts = fx::starts(order);
    auto options = fx::max3_options(order);
    if (expected_fractions != nullptr) {
        options.three_phase_starts.front().phase_fraction_seed = {
            (*expected_fractions)[1], (*expected_fractions)[2]};
    }
    return fl::solve_pr76_pt_max3(
        ref::pressure_pa, ref::temperature_k, benchmark_feed,
        evaluator, options, starts, starts);
}

void require_closed_max3(
    const fl::Pr76PtMax3Result& result,
    const std::vector<std::string>& order,
    const std::array<double, 3>& expected_fractions =
        ref::phase_fractions) {
    if (result.status != fl::Pr76PtMax3Status::three_phase ||
        result.three_phase_candidate() == nullptr ||
        !result.selected_attempt) {
        print_max3(result);
        throw std::runtime_error(
            "sour-gas PR76 max3 did not close the literature-defined three-phase state");
    }
    require(result.base.solution.status == fl::PtSplitStatus::phase_set_unstable,
            "sour-gas 2->3 route lacks rejected/incomplete two-phase evidence");
    require_state_invariants(*result.three_phase_candidate());
    (void)match_reference(
        *result.three_phase_candidate(), order, expected_fractions);
    const auto& accepted_attempt = result.attempts[*result.selected_attempt];
    require(accepted_attempt.final_stability &&
                accepted_attempt.final_stability->status ==
                    fl::StabilityStatus::no_instability_found,
            "sour-gas three-phase candidate lacks final common-tangent closure");
    require(accepted_attempt.accepted_three_phase,
            "sour-gas selected attempt was not marked accepted three phase");
    require(!result.global_stability_proven,
            "finite PR76 three-phase review was promoted to global proof");

    const auto published = fl::project_pr76_pt_max3_phase_set(result);
    require(published.solution.status == fl::PtPhaseSetStatus::accepted &&
                published.solution.accepted_phase_count() == 3U,
            "sour-gas authoritative phase-set projection did not publish three phases");
}

void seeded_max3() {
    const auto order = fx::canonical_order();
    const auto result = solve_max3(order, fx::feed(order));
    require_closed_max3(result, order);
}

void interior_reweighted_three_phase() {
    const auto order = fx::canonical_order();
    const auto starts = fx::starts(order);
    const std::array<double, 3> fractions{0.20, 0.35, 0.45};
    Vec reweighted(ref::feed.size(), 0.0);
    for (std::size_t component = 0; component < reweighted.size(); ++component) {
        for (std::size_t phase = 0; phase < 3U; ++phase) {
            reweighted[component] +=
                fractions[phase] * starts[phase][component];
        }
    }
    const auto result = solve_max3(order, reweighted, &fractions);
    require_closed_max3(result, order, fractions);
}

void component_permutation() {
    const auto first_order = fx::canonical_order();
    const auto second_order = fx::reversed_order();
    const auto first = solve_max3(first_order, fx::feed(first_order));
    const auto second = solve_max3(second_order, fx::feed(second_order));
    require_closed_max3(first, first_order);
    require_closed_max3(second, second_order);
}

void backend_publication() {
    const auto order = fx::canonical_order();
    const auto model = fx::model(order);
    fl::Pr76VleEvaluator evaluator(model);
    fl::Pr76PtFlashBackend backend(evaluator, fx::backend_options(order));
    fl::PtFlashBackend& runtime = backend;
    const auto result = runtime.solve({
        ref::pressure_pa, ref::temperature_k, fx::feed(order)});

    require(result.structurally_valid(),
            "sour-gas public PR76 backend result is structurally invalid");
    require(result.solution.status == fl::PtPhaseSetStatus::accepted &&
                result.solution.accepted_phase_count() == 3U,
            "sour-gas public PR76 backend did not publish three phases");
    require(!result.solution.global_stability_proven &&
                !result.capability.global_stability_proven,
            "public sour-gas backend promoted finite search to global proof");

    bool accepted_2_to_3 = false;
    for (const auto& evidence : result.transition_report.evidence) {
        if (evidence.source_phase_count == 2U &&
            evidence.target_phase_count == 3U &&
            evidence.resolution ==
                fl::PtPhaseTransitionResolution::accepted_target &&
            evidence.fresh_target_solve_attempted &&
            evidence.target_topology_closed) {
            accepted_2_to_3 = true;
        }
    }
    require(accepted_2_to_3,
            "sour-gas public backend lost fresh accepted 2->3 evidence");

    const auto* accepted = result.accepted_phase_set();
    require(accepted != nullptr && accepted->phases.size() == 3U,
            "sour-gas backend accepted helper lost phase set");

    std::array<bool, 3> used{false, false, false};
    for (const auto& phase : accepted->phases) {
        std::optional<std::size_t> best;
        double best_error = 1.0e100;
        for (std::size_t candidate = 0; candidate < 3U; ++candidate) {
            if (used[candidate]) { continue; }
            double error = 0.0;
            for (std::size_t component = 0; component < ref::feed.size(); ++component) {
                error = std::max(
                    error,
                    std::abs(
                        phase.composition[component] -
                        ref::phases[candidate][component]));
            }
            if (error < best_error) {
                best_error = error;
                best = candidate;
            }
        }
        require(best.has_value() && best_error <= reference_budget,
                "public backend phase left independent sour-gas reference");
        used[*best] = true;
        require(
            std::abs(
                phase.mole_phase_fraction -
                ref::phase_fractions[*best]) <= reference_budget,
            "public backend phase fraction left sour-gas reference");
        require(phase.compressibility_factor.has_value() &&
                    std::abs(
                        *phase.compressibility_factor -
                        ref::compressibility_factors[*best]) <= reference_budget,
                "public backend Z left sour-gas reference");
    }
}

} // namespace

int main(int argc, char** argv) {
    try {
        if (argc != 2) { throw std::invalid_argument("one test name required"); }
        const std::string_view name{argv[1]};
        if (name == "reference_candidate") reference_candidate();
        else if (name == "seeded_max3") seeded_max3();
        else if (name == "interior_reweighted") interior_reweighted_three_phase();
        else if (name == "component_permutation") component_permutation();
        else if (name == "backend_publication") backend_publication();
        else throw std::invalid_argument("unknown test name");
        std::cout << "[PASS] " << name << '\n';
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "[FAIL] " << error.what() << '\n';
        return 1;
    }
}
