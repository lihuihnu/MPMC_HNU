#include <mpmc/flash/cpa_pt_flash_backend.hpp>

#include "test_support.hpp"
#include <fugacity_adapter_regression.hpp>

#include <algorithm>
#include <array>
#include <cmath>
#include <stdexcept>
#include <string_view>
#include <vector>

namespace {
namespace fl = mpmc::flash;
namespace fo = mpmc::flow;
namespace fa = mpmc::test::flow_adapter;
namespace th = mpmc::thermodynamics;
using Vec = std::vector<double>;

void require(bool value, const char* message) {
    if (!value) { throw std::runtime_error(message); }
}

Vec log_ratio(const Vec& numerator, const Vec& denominator) {
    require(numerator.size() == denominator.size(), "log ratio dimension mismatch");
    Vec result(numerator.size(), 0.0);
    for (std::size_t i = 0; i < result.size(); ++i) {
        result[i] = std::log(numerator[i]) - std::log(denominator[i]);
    }
    return result;
}

fl::CpaPtMax3Options max3_options(
    std::array<double, 3> beta,
    bool swap_heavy = false) {
    fl::CpaPtMax3Options options;
    options.two_phase.initial_stability.automatic_starts = false;
    options.two_phase.final_stability.automatic_starts = false;
    options.final_three_phase_stability.automatic_starts = false;
    options.three_phase_starts.push_back(
        cpa_max3_test::three_phase_start(beta, swap_heavy));
    return options;
}

fl::CpaPtMax3Result solve_three_phase_fraction(
    std::array<double, 3> beta,
    bool swap_heavy = false) {
    const auto parameters = cpa_max3_test::parameters(swap_heavy);
    const auto model = th::CpaPtPhase::from_parameters(parameters);
    fl::CpaVleEvaluator evaluator(model, cpa_max3_test::fast_pt_options());
    const auto starts = cpa_max3_test::starts(swap_heavy);
    const auto feed = cpa_max3_test::feed_from_phase_fractions(beta, swap_heavy);
    auto options = max3_options(beta, swap_heavy);
    return fl::solve_cpa_pt_max3(
        cpa_max3_test::pressure_pa, cpa_max3_test::temperature_k,
        feed, evaluator, options, starts, starts);
}

bool phase_equal(const fl::PtCandidatePhase& first,
                 const fl::PtCandidatePhase& second) {
    return first.mole_phase_fraction == second.mole_phase_fraction &&
           first.composition == second.composition &&
           first.activity.branch == second.activity.branch &&
           first.activity.smooth == second.activity.smooth &&
           first.activity.ln_phi == second.activity.ln_phi &&
           first.compressibility_factor == second.compressibility_factor;
}

void fixed_three_phase_candidate() {
    const auto parameters = cpa_max3_test::parameters();
    const auto model = th::CpaPtPhase::from_parameters(parameters);
    const auto phases = cpa_max3_test::starts();
    const auto feed = cpa_max3_test::feed_from_phase_fractions(
        {1.0 / 3.0, 1.0 / 3.0, 1.0 / 3.0});
    fl::CpaThreePhaseEvaluator evaluator(
        model,
        {fl::CpaRootSide::upper_density_admissible,
         fl::CpaRootSide::lower_density_admissible,
         fl::CpaRootSide::upper_density_admissible},
        cpa_max3_test::fast_pt_options());
    const auto result = fl::iterate_pt_three_phase(
        cpa_max3_test::pressure_pa, cpa_max3_test::temperature_k,
        feed, log_ratio(phases[1], phases[0]),
        log_ratio(phases[2], phases[0]),
        {1.0 / 3.0, 1.0 / 3.0}, evaluator);
    require(result.candidate_admissible(),
            "CPA structural fixed-three-phase equations did not converge");
    require(result.point &&
                result.point->chemical_potential_norm <=
                    result.options.chemical_potential_tolerance &&
                result.point->mass_absolute <=
                    result.options.mass_absolute_tolerance,
            "CPA fixed-three-phase candidate lost equilibrium or balance closure");
}

void two_to_three_baseline() {
    constexpr std::array<double, 3> beta{
        1.0 / 3.0, 1.0 / 3.0, 1.0 / 3.0};
    const auto result = solve_three_phase_fraction(beta);
    require(result.base.solution.status == fl::PtSplitStatus::phase_set_unstable,
            "CPA max3 fixture no longer starts from rejected two-phase topology");
    require(result.status == fl::CpaPtMax3Status::three_phase &&
                result.three_phase_candidate() != nullptr &&
                result.selected_attempt.has_value(),
            "CPA max3 did not accept the structural three-phase topology");
    const auto& attempt = result.attempts[*result.selected_attempt];
    require(attempt.supplied_start && *attempt.supplied_start == 0U &&
                attempt.equilibrium.candidate_admissible() &&
                attempt.final_stability &&
                attempt.final_stability->status ==
                    fl::StabilityStatus::no_instability_found &&
                attempt.accepted_three_phase,
            "CPA 2->3 route did not close equilibrium plus final all-root stability");
    const auto published = fl::project_cpa_pt_max3_phase_set(result);
    const auto* phases = published.solution.accepted_phase_set();
    require(phases != nullptr && phases->phases.size() == 3U &&
                published.solution.capability.maximum_phase_count == 3U,
            "CPA max3 publication lost accepted three-phase state");
}


void flow_fugacity_adapter_on_accepted_three_phase() {
    constexpr std::array<double, 3> beta{
        1.0 / 3.0, 1.0 / 3.0, 1.0 / 3.0};
    const auto result = solve_three_phase_fraction(beta);
    require(result.status == fl::CpaPtMax3Status::three_phase &&
                result.three_phase_candidate() != nullptr &&
                result.selected_attempt.has_value(),
            "CPA flow adapter fixture lost accepted three-phase topology");
    const auto& attempt = result.attempts[*result.selected_attempt];
    require(attempt.accepted_three_phase &&
                attempt.final_stability &&
                attempt.final_stability->status ==
                    fl::StabilityStatus::no_instability_found,
            "CPA flow adapter fixture lost accepted/stable three-phase evidence");

    const auto& state = *result.three_phase_candidate();
    const auto parameters = cpa_max3_test::parameters();
    const auto model = th::CpaPtPhase::from_parameters(parameters);

    std::array<
        fo::CpaSelectedPhaseFugacityEvaluator3P::Selection,
        3>
        selections{};
    std::array<double, 3> fractions{};
    std::array<Vec, 3> compositions;
    for (std::size_t phase = 0U; phase < 3U; ++phase) {
        selections[phase].root_index =
            state.phases[phase].activity.branch;
        selections[phase].options =
            cpa_max3_test::fast_pt_options();
        fractions[phase] =
            state.phases[phase].mole_phase_fraction;
        compositions[phase] =
            state.phases[phase].composition;
    }

    fo::CpaSelectedPhaseFugacityEvaluator3P adapter{
        model, selections};
    fa::verify_accepted_state_and_jacobian(
        cpa_max3_test::pressure_pa,
        cpa_max3_test::temperature_k,
        fractions,
        compositions,
        adapter,
        5.0e-8,
        {2.0e-9, 5.0e-5, 2.0e-4, 5.0e-3});
}

void three_to_two_fresh_neighbor() {
    constexpr double epsilon = 1.0e-10;
    constexpr std::array<double, 3> beta{
        0.5 * (1.0 - epsilon),
        0.5 * (1.0 - epsilon),
        epsilon};
    const auto result = solve_three_phase_fraction(beta);
    require(result.base.solution.status == fl::PtSplitStatus::phase_set_unstable,
            "CPA 3->2 fixture no longer has additional-phase evidence at the base pair");
    require(result.status == fl::CpaPtMax3Status::two_phase &&
                result.selected_attempt.has_value() &&
                result.two_phase_neighbor() != nullptr,
            "CPA max3 did not fresh-resolve the disappearance boundary as two phase");
    const auto& attempt = result.attempts[*result.selected_attempt];
    require(attempt.equilibrium.status ==
                fl::PtThreePhaseStatus::phase_disappearance &&
                attempt.equilibrium.disappearing_phase.has_value() &&
                attempt.boundary_neighbor.has_value() &&
                attempt.boundary_neighbor->solution.status ==
                    fl::PtSplitStatus::two_phase_no_instability_found &&
                attempt.accepted_two_phase_neighbor,
            "CPA 3->2 route deleted a phase without a fresh stable two-phase solve");
    const auto published = fl::project_cpa_pt_max3_phase_set(result);
    const auto* phases = published.solution.accepted_phase_set();
    require(phases != nullptr && phases->phases.size() == 2U,
            "CPA 3->2 publication did not own the fresh accepted two-phase neighbor");
    const auto transition = fl::project_cpa_pt_max3_transition_report(result);
    const bool has_3_to_2 = std::any_of(
        transition.evidence.begin(), transition.evidence.end(),
        [](const fl::PtPhaseTransitionEvidence& evidence) {
            return evidence.source_phase_count == 3U &&
                   evidence.target_phase_count == 2U &&
                   evidence.resolution ==
                       fl::PtPhaseTransitionResolution::accepted_target &&
                   evidence.fresh_target_solve_attempted &&
                   evidence.target_topology_closed;
        });
    require(has_3_to_2,
            "CPA max3 transition projection lost accepted 3->2 fresh-neighbor evidence");
}

void component_permutation() {
    constexpr std::array<double, 3> beta{
        1.0 / 3.0, 1.0 / 3.0, 1.0 / 3.0};
    const auto first = solve_three_phase_fraction(beta, false);
    const auto second = solve_three_phase_fraction(beta, true);
    require(first.status == fl::CpaPtMax3Status::three_phase &&
                second.status == fl::CpaPtMax3Status::three_phase &&
                first.three_phase_candidate() && second.three_phase_candidate(),
            "CPA max3 component permutation changed accepted topology");
    const auto& a = *first.three_phase_candidate();
    const auto& b = *second.three_phase_candidate();
    for (std::size_t phase = 0; phase < 3U; ++phase) {
        require(std::abs(a.phases[phase].mole_phase_fraction -
                         b.phases[phase].mole_phase_fraction) < 2.0e-11,
                "CPA max3 component permutation changed phase fraction");
        require(std::abs(a.phases[phase].composition[0] -
                         b.phases[phase].composition[0]) < 2.0e-10 &&
                    std::abs(a.phases[phase].composition[1] -
                             b.phases[phase].composition[2]) < 2.0e-10 &&
                    std::abs(a.phases[phase].composition[2] -
                             b.phases[phase].composition[1]) < 2.0e-10,
                "CPA max3 compositions changed under B/C permutation");
    }
}

void backend_equivalence() {
    constexpr std::array<double, 3> beta{
        1.0 / 3.0, 1.0 / 3.0, 1.0 / 3.0};
    const auto parameters = cpa_max3_test::parameters();
    const auto model = th::CpaPtPhase::from_parameters(parameters);
    fl::CpaVleEvaluator evaluator(model, cpa_max3_test::fast_pt_options());
    const auto starts = cpa_max3_test::starts();
    const auto feed = cpa_max3_test::feed_from_phase_fractions(beta);

    auto direct_options = max3_options(beta);
    const auto direct = fl::solve_cpa_pt_max3(
        cpa_max3_test::pressure_pa, cpa_max3_test::temperature_k,
        feed, evaluator, direct_options, starts, starts);
    const auto published = fl::project_cpa_pt_max3_phase_set(direct);

    fl::CpaPtFlashBackendOptions options;
    options.split = direct_options.two_phase;
    options.three_phase = direct_options.three_phase;
    options.final_three_phase_stability =
        direct_options.final_three_phase_stability;
    options.three_phase_starts = direct_options.three_phase_starts;
    options.initial_starts = starts;
    options.final_starts = starts;
    fl::CpaPtFlashBackend backend(evaluator, options);
    fl::PtFlashBackend& runtime = backend;
    const auto adapted = runtime.solve({
        cpa_max3_test::pressure_pa,
        cpa_max3_test::temperature_k,
        feed});

    require(adapted.structurally_valid() &&
                runtime.capability().supports_phase_count(3U) &&
                runtime.capability().maximum_phase_count() == 3U &&
                runtime.capability().transition_capability.edge(2U, 3U) != nullptr &&
                runtime.capability().transition_capability.edge(3U, 2U) != nullptr &&
                adapted.solution.status == published.solution.status &&
                adapted.solution.accepted_phase_count() == 3U &&
                adapted.solution.feed == published.solution.feed,
            "CPA max3 runtime backend capability/publication mismatch");
    const auto* a = adapted.solution.accepted_phase_set();
    const auto* b = published.solution.accepted_phase_set();
    require(a != nullptr && b != nullptr && a->phases.size() == b->phases.size(),
            "CPA max3 runtime backend lost accepted phase set");
    for (std::size_t phase = 0; phase < a->phases.size(); ++phase) {
        require(phase_equal(a->phases[phase], b->phases[phase]),
                "CPA max3 runtime backend changed phase payload");
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
    require(has_2_to_3 && !adapted.morphology_resolved,
            "CPA max3 backend lost transition evidence or invented morphology");
}

void hints_are_not_phase_count_evidence() {
    const auto parameters = cpa_max3_test::parameters();
    const auto model = th::CpaPtPhase::from_parameters(parameters);
    fl::CpaVleEvaluator evaluator(model, cpa_max3_test::fast_pt_options());
    const auto starts = cpa_max3_test::starts();
    const auto feed = starts[0];
    fl::CpaPtMax3Options options;
    options.two_phase.initial_stability.automatic_starts = false;
    options.two_phase.final_stability.automatic_starts = false;
    options.final_three_phase_stability.automatic_starts = false;
    options.three_phase_starts.push_back(cpa_max3_test::three_phase_start());
    const auto result = fl::solve_cpa_pt_max3(
        cpa_max3_test::pressure_pa, cpa_max3_test::temperature_k,
        feed, evaluator, options, starts, starts);
    require(result.base.solution.status ==
                fl::PtSplitStatus::single_phase_no_instability_found &&
                result.status == fl::CpaPtMax3Status::single_phase &&
                result.attempts.empty(),
            "CPA max3 continuation hint created phase-count evidence at a stable feed");
}

void malformed_three_phase_start_rejected() {
    const auto parameters = cpa_max3_test::parameters();
    const auto model = th::CpaPtPhase::from_parameters(parameters);
    fl::CpaVleEvaluator evaluator(model, cpa_max3_test::fast_pt_options());
    const auto feed = cpa_max3_test::feed_from_phase_fractions(
        {1.0 / 3.0, 1.0 / 3.0, 1.0 / 3.0});
    fl::CpaPtMax3Options options;
    auto malformed = cpa_max3_test::three_phase_start();
    malformed.compositions[2].pop_back();
    options.three_phase_starts.push_back(std::move(malformed));
    bool threw = false;
    try {
        (void)fl::solve_cpa_pt_max3(
            cpa_max3_test::pressure_pa,
            cpa_max3_test::temperature_k,
            feed, evaluator, options);
    } catch (const std::invalid_argument&) {
        threw = true;
    }
    require(threw,
            "CPA max3 malformed continuation start was not rejected before scientific solve");
}

using Test = std::pair<std::string_view, void (*)()>;
constexpr Test tests[]{
    {"fixed_three_phase", fixed_three_phase_candidate},
    {"two_to_three", two_to_three_baseline},
    {"flow_fugacity_adapter", flow_fugacity_adapter_on_accepted_three_phase},
    {"three_to_two", three_to_two_fresh_neighbor},
    {"component_permutation", component_permutation},
    {"backend_equivalence", backend_equivalence},
    {"hint_not_evidence", hints_are_not_phase_count_evidence},
    {"malformed_start", malformed_three_phase_start_rejected}};

} // namespace

int main(int argc, char** argv) {
    try {
        if (argc != 2) { throw std::invalid_argument("one test name required"); }
        for (const auto& [name, run] : tests) {
            if (name == argv[1]) {
                run();
                return 0;
            }
        }
        throw std::invalid_argument("unknown test name");
    } catch (const std::exception&) {
        return 1;
    }
}
