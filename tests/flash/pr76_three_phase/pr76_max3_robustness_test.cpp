#include <mpmc/flash/pr76_max3_phase_set.hpp>

#include "synthetic_fixture.hpp"

#include <array>
#include <cmath>
#include <optional>
#include <stdexcept>
#include <string>
#include <string_view>
#include <vector>

namespace {
namespace fl = mpmc::flash;
namespace th = mpmc::thermodynamics;
using pr76_max3_test::Vec;

void require(bool value, const char* message) {
    if (!value) { throw std::runtime_error(message); }
}

Vec log_ratio(const Vec& numerator, const Vec& denominator) {
    require(numerator.size() == denominator.size(), "log-ratio shape mismatch");
    Vec value(numerator.size(), 0.0);
    for (std::size_t i = 0; i < value.size(); ++i) {
        value[i] = std::log(numerator[i]) - std::log(denominator[i]);
    }
    return value;
}

fl::Pr76PtThreePhaseStart reference_start(const fl::PtThreePhaseState& state) {
    fl::Pr76PtThreePhaseStart start;
    for (std::size_t phase = 0; phase < 3U; ++phase) {
        start.compositions[phase] = state.phases[phase].composition;
    }
    start.phase_fraction_seed = {
        state.phases[1].mole_phase_fraction,
        state.phases[2].mole_phase_fraction};
    return start;
}

fl::Pr76PtThreePhaseStart exact_fixture_start() {
    fl::Pr76PtThreePhaseStart start;
    const auto& phases = pr76_max3_test::reference_phases();
    start.compositions = {phases[0], phases[1], phases[2]};
    return start;
}

// Fixed synthetic topology fixture used ONLY to exercise the real max3 3->2
// recovery branch. These artificial kij values are structural regression data,
// not a fit, parameter recommendation, or physical validation dataset.
th::Pr76Phase<double> recovery_model() {
    const auto provenance = pr76_max3_test::source(
        "fixed three-phase disappearance recovery topology fixture");
    const std::array<std::string, 3> ids{"light", "heavy-b", "heavy-c"};
    const std::array<double, 3> tc{190.6, 500.0, 500.0};
    const std::array<double, 3> pc{4.6e6, 5.0e6, 5.0e6};
    const std::array<double, 3> omega{0.01, 0.10, 0.10};

    std::vector<th::Component> catalog;
    th::PrParameterInput input;
    input.model_id = std::string(th::pr76_profile);
    input.dataset_id = "synthetic-PR76-max3-3to2-recovery";
    input.revision = "v1";
    input.applicability = {std::nullopt, std::nullopt, provenance};
    for (std::size_t i = 0U; i < ids.size(); ++i) {
        catalog.push_back({ids[i], ids[i], th::ComponentKind::pure, provenance, {}});
        input.pure.push_back({
            ids[i], pr76_max3_test::scalar(tc[i], th::Unit::kelvin, provenance),
            pr76_max3_test::scalar(pc[i], th::Unit::pascal, provenance),
            pr76_max3_test::scalar(omega[i], th::Unit::dimensionless, provenance)});
    }
    input.binary.push_back({ids[0], ids[1],
        pr76_max3_test::scalar(0.0, th::Unit::dimensionless, provenance)});
    input.binary.push_back({ids[0], ids[2],
        pr76_max3_test::scalar(0.03, th::Unit::dimensionless, provenance)});
    input.binary.push_back({ids[1], ids[2],
        pr76_max3_test::scalar(0.36, th::Unit::dimensionless, provenance)});

    const std::vector<std::string> order{ids[0], ids[1], ids[2]};
    return th::Pr76Phase<double>::from_parameters(
        th::PrParameterSet::create(
            catalog, order, input,
            th::DataPolicy::allow_synthetic_tests));
}

const std::array<Vec, 3>& recovery_phases() {
    // Frozen numerical anchors from a real accepted strict-PR76 three-phase solve
    // of recovery_model() at p=1 MPa, T=250 K. They exist only to make this
    // orchestration regression deterministic; max3 re-solves all equations below.
    static const std::array<Vec, 3> values{{
        {0.96854260046818119, 0.015729641629537154, 0.01572775790228163},
        {0.059133645515199861, 0.93843673665865202, 0.0024296178261481098},
        {0.050218693594532689, 0.0023772256933525067, 0.94740408071211479}}};
    return values;
}

void automatic_cold_start_is_exercised() {
    const auto model = pr76_max3_test::model();
    fl::Pr76VleEvaluator evaluator(model);
    const auto feed = pr76_max3_test::equal_feed();
    const auto vle_starts = pr76_max3_test::starts();

    fl::Pr76PtMax3Options options;
    options.max_three_phase_attempts = 16U;
    const auto result = fl::solve_pr76_pt_max3(
        1.0e6, 250.0, feed, evaluator, options, vle_starts, vle_starts);

    require(result.base.solution.status == fl::PtSplitStatus::phase_set_unstable,
            "structural cold-start fixture no longer reaches two-phase final instability");
    bool saw_automatic = false;
    for (const auto& attempt : result.attempts) {
        if (attempt.witness_trial.has_value() && !attempt.supplied_start.has_value()) {
            saw_automatic = true;
            break;
        }
    }
    require(saw_automatic,
            "negative-TPD witness did not generate an automatic three-phase attempt");
}

void continuation_quota_does_not_starve_automatic() {
    const auto model = pr76_max3_test::model();
    fl::Pr76VleEvaluator evaluator(model);
    const auto feed = pr76_max3_test::equal_feed();
    const auto vle_starts = pr76_max3_test::starts();

    fl::Pr76PtMax3Options options;
    options.max_three_phase_attempts = 1U;
    fl::Pr76PtThreePhaseStart poor_hint;
    poor_hint.compositions = {feed, feed, feed};
    options.three_phase_starts.push_back(poor_hint);

    const auto result = fl::solve_pr76_pt_max3(
        1.0e6, 250.0, feed, evaluator, options, vle_starts, vle_starts);
    require(result.base.solution.status == fl::PtSplitStatus::phase_set_unstable,
            "budget fixture no longer has final two-phase instability evidence");
    require(result.base.solution.final_stability.has_value(),
            "budget fixture lost final two-phase TPD evidence");

    std::size_t supplied_attempts = 0U;
    std::size_t automatic_attempts = 0U;
    for (const auto& attempt : result.attempts) {
        if (attempt.supplied_start.has_value()) {
            ++supplied_attempts;
            require(!attempt.witness_trial.has_value(),
                    "supplied three-phase hint was incorrectly tagged as TPD evidence");
            continue;
        }
        if (attempt.witness_trial.has_value()) {
            ++automatic_attempts;
            const std::size_t witness = *attempt.witness_trial;
            require(witness < result.base.solution.final_stability->trials.size(),
                    "automatic attempt lost its originating final-TPD trial");
            const auto& trial = result.base.solution.final_stability->trials[witness];
            require(trial.status == fl::StabilityTrialStatus::negative_tpd &&
                        trial.point.has_value(),
                    "automatic attempt did not originate from retained negative TPD evidence");
        }
    }
    require(supplied_attempts == 1U,
            "supplied-source quota did not stop independently at one attempt");
    require(automatic_attempts == 1U,
            "supplied-source quota starved or altered the one-attempt automatic quota");
    require(result.attempt_limit_reached,
            "per-source quota exhaustion was not retained in diagnostics");
}

void accepted_three_phase_continues_three_to_three() {
    const auto model = pr76_max3_test::model();
    fl::Pr76VleEvaluator evaluator(model);
    const auto feed = pr76_max3_test::equal_feed();
    const auto vle_starts = pr76_max3_test::starts();

    fl::Pr76PtMax3Options first_options;
    first_options.three_phase_starts.push_back(exact_fixture_start());
    const auto first = fl::solve_pr76_pt_max3(
        1.0e6, 250.0, feed, evaluator,
        first_options, vle_starts, vle_starts);
    require(first.status == fl::Pr76PtMax3Status::three_phase &&
                first.three_phase_candidate() != nullptr,
            "reference seeded solve did not establish the three-phase continuation base");

    fl::Pr76PtMax3Options second_options;
    second_options.three_phase_starts.push_back(
        reference_start(*first.three_phase_candidate()));
    const auto second = fl::solve_pr76_pt_max3(
        1.0e6, 250.0, feed, evaluator,
        second_options, vle_starts, vle_starts);
    require(second.status == fl::Pr76PtMax3Status::three_phase &&
                second.three_phase_candidate() != nullptr &&
                second.selected_attempt.has_value(),
            "accepted three-phase state did not survive a fresh 3->3 continuation solve");
    require(second.attempts[*second.selected_attempt].supplied_start.has_value(),
            "3->3 continuation lost supplied-start provenance");
}

void three_to_two_boundary_requires_fresh_neighbor() {
    // Keep the original exact symmetric tie-edge checks: the low-level fixed
    // three-phase primitive must retain disappearance evidence, and the surviving
    // pair must independently close through a fresh full VLE solve.
    const auto model = pr76_max3_test::model();
    fl::Pr76VleEvaluator evaluator(model);
    const auto& phases = pr76_max3_test::reference_phases();
    constexpr double beta1 = 0.4;
    Vec feed(3U, 0.0);
    for (std::size_t i = 0; i < feed.size(); ++i) {
        feed[i] = (1.0 - beta1) * phases[0][i] + beta1 * phases[1][i];
    }

    fl::Pr76ThreePhaseEvaluator three_phase_evaluator(
        model,
        {fl::Pr76RootSide::lower_admissible,
         fl::Pr76RootSide::upper_admissible,
         fl::Pr76RootSide::lower_admissible});
    const auto boundary = fl::iterate_pt_three_phase(
        1.0e6, 250.0, feed,
        log_ratio(phases[1], phases[0]),
        log_ratio(phases[2], phases[0]),
        {beta1, 0.0}, three_phase_evaluator);
    require(boundary.status == fl::PtThreePhaseStatus::phase_disappearance &&
                boundary.disappearing_phase.has_value() &&
                *boundary.disappearing_phase == 2U && boundary.equations_converged(),
            "3->2 continuation edge did not retain explicit disappearance evidence");
    require(boundary.point.has_value() &&
                boundary.point->phases[2].mole_phase_fraction <=
                    boundary.options.minimum_phase_fraction,
            "disappearance evidence does not correspond to the configured phase-fraction gate");

    const std::vector<Vec> surviving{
        boundary.point->phases[0].composition,
        boundary.point->phases[1].composition};
    const auto neighbor = fl::solve_pr76_pt_vle(
        1.0e6, 250.0, feed, evaluator, {}, surviving, surviving);
    require(neighbor.solution.status ==
                fl::PtSplitStatus::two_phase_no_instability_found &&
                neighbor.solution.candidate() != nullptr &&
                neighbor.solution.final_stability.has_value() &&
                neighbor.solution.final_stability->status ==
                    fl::StabilityStatus::no_instability_found,
            "3->2 boundary was not independently fresh-resolved as a stable two-phase neighbor");
    const auto& pair = *neighbor.solution.candidate();
    require(pair.fugacity_norm <= neighbor.solution.options.iteration.fugacity_tolerance &&
                pair.fractions.mass_absolute <=
                    neighbor.solution.options.iteration.mass_absolute_tolerance &&
                pair.fractions.mass_relative <=
                    neighbor.solution.options.iteration.mass_relative_tolerance,
            "fresh two-phase neighbor did not re-close equilibrium/material-balance gates");

    // Publication must not expose a fresh neighbor until max3 explicitly accepts
    // it, and the generic projection still respects the neighbor's own status.
    fl::Pr76PtMax3Result publication_gate;
    publication_gate.status = fl::Pr76PtMax3Status::two_phase;
    publication_gate.selected_attempt = 0U;
    publication_gate.attempts.resize(1U);
    publication_gate.attempts[0].boundary_neighbor = neighbor;
    require(publication_gate.two_phase_neighbor() == nullptr,
            "unaccepted fresh-neighbor result leaked through the max3 accessor");
    require(fl::project_pr76_pt_max3_phase_set(publication_gate).solution.status !=
                fl::PtPhaseSetStatus::accepted,
            "unaccepted fresh-neighbor result was published as accepted two phase");
    publication_gate.attempts[0].accepted_two_phase_neighbor = true;
    require(publication_gate.two_phase_neighbor() != nullptr,
            "accepted fresh-neighbor result was not exposed by the max3 accessor");
    require(fl::project_pr76_pt_max3_phase_set(publication_gate).solution.status ==
                fl::PtPhaseSetStatus::accepted,
            "accepted fresh-neighbor result was not published as two phase");
    auto unresolved_neighbor = neighbor;
    unresolved_neighbor.solution.status = fl::PtSplitStatus::indeterminate;
    publication_gate.attempts[0].boundary_neighbor = std::move(unresolved_neighbor);
    require(fl::project_pr76_pt_max3_phase_set(publication_gate).solution.status !=
                fl::PtPhaseSetStatus::accepted,
            "non-accepted fresh-neighbor status was promoted by the max3 publication layer");

    // End-to-end single-call recovery regression. The fixed synthetic model has a
    // reachable metastable two-phase basin. Its final common-tangent review is
    // negative, which legitimately opens the max3 route; the supplied triple is
    // only an initialization hint. Phase 1 has exactly zero lever-rule share, so
    // the fixed-three-phase solve must report phase_disappearance and max3 must
    // fresh-resolve phases 0+2 before publishing two phase.
    auto recovery = recovery_model();
    fl::Pr76VleEvaluator recovery_evaluator(recovery);
    const auto& recovery_phase = recovery_phases();
    Vec recovery_feed(3U, 0.0);
    for (std::size_t i = 0U; i < recovery_feed.size(); ++i) {
        recovery_feed[i] = 0.25 * recovery_phase[0][i] +
                           0.75 * recovery_phase[2][i];
    }

    const Vec metastable_witness{
        0.027979967031294486,
        0.90057153296773984,
        0.071448500000965642};
    const std::vector<Vec> initial_starts{metastable_witness};
    const std::vector<Vec> final_starts{
        recovery_phase[0], recovery_phase[1], recovery_phase[2]};

    fl::Pr76PtMax3Options recovery_options;
    recovery_options.two_phase.initial_stability.automatic_starts = false;
    // The frozen metastable split is the second role assignment for this witness.
    recovery_options.two_phase.max_split_attempts = 2U;
    fl::Pr76PtThreePhaseStart recovery_start;
    recovery_start.compositions = recovery_phase;
    recovery_start.phase_fraction_seed = {0.0, 0.75};
    recovery_options.three_phase_starts.push_back(recovery_start);

    const auto routed = fl::solve_pr76_pt_max3(
        1.0e6, 250.0, recovery_feed, recovery_evaluator,
        recovery_options, initial_starts, final_starts);

    require(routed.base.solution.status == fl::PtSplitStatus::phase_set_unstable &&
                routed.base.solution.final_stability.has_value() &&
                routed.base.solution.final_stability->status ==
                    fl::StabilityStatus::unstable,
            "single-call 3->2 fixture lost the required negative two-phase final-TPD evidence");
    require(routed.status == fl::Pr76PtMax3Status::two_phase &&
                routed.selected_attempt.has_value() &&
                routed.two_phase_neighbor() != nullptr,
            "single-call max3 did not publish the fresh two-phase recovery target");

    const auto& routed_attempt = routed.attempts[*routed.selected_attempt];
    require(routed_attempt.supplied_start.has_value() &&
                *routed_attempt.supplied_start == 0U &&
                routed_attempt.equilibrium.status ==
                    fl::PtThreePhaseStatus::phase_disappearance &&
                routed_attempt.equilibrium.disappearing_phase.has_value() &&
                *routed_attempt.equilibrium.disappearing_phase == 1U &&
                routed_attempt.equilibrium.equations_converged() &&
                routed_attempt.equilibrium.point.has_value() &&
                routed_attempt.equilibrium.point->phases[1].mole_phase_fraction <=
                    routed_attempt.equilibrium.options.minimum_phase_fraction,
            "selected max3 attempt lost explicit phase-1 disappearance evidence");
    require(routed_attempt.accepted_two_phase_neighbor &&
                routed_attempt.boundary_neighbor.has_value(),
            "phase disappearance was not closed by a retained fresh two-phase neighbor");

    const auto& fresh = *routed_attempt.boundary_neighbor;
    require(fresh.solution.status == fl::PtSplitStatus::two_phase_no_instability_found &&
                fresh.solution.candidate() != nullptr &&
                fresh.solution.final_stability.has_value() &&
                fresh.solution.final_stability->status ==
                    fl::StabilityStatus::no_instability_found,
            "fresh recovery neighbor did not independently pass its own final phase-set review");
    const auto& fresh_pair = *fresh.solution.candidate();
    require(fresh_pair.fugacity_norm <= fresh.solution.options.iteration.fugacity_tolerance &&
                fresh_pair.fractions.mass_absolute <=
                    fresh.solution.options.iteration.mass_absolute_tolerance &&
                fresh_pair.fractions.mass_relative <=
                    fresh.solution.options.iteration.mass_relative_tolerance,
            "fresh recovery neighbor failed its own equilibrium/material-balance tolerances");

    const auto published = fl::project_pr76_pt_max3_phase_set(routed);
    require(published.solution.status == fl::PtPhaseSetStatus::accepted &&
                published.solution.accepted_phase_count() == 2U &&
                published.solution.accepted_phase_set() != nullptr,
            "single-call 3->2 recovery was not published as an accepted role-neutral two-phase set");
}

using Test = std::pair<std::string_view, void (*)()>;
constexpr Test tests[]{
    {"automatic_cold_start", automatic_cold_start_is_exercised},
    {"continuation_budget", continuation_quota_does_not_starve_automatic},
    {"continuation_3_to_3", accepted_three_phase_continues_three_to_three},
    {"continuation_3_to_2", three_to_two_boundary_requires_fresh_neighbor}};

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
