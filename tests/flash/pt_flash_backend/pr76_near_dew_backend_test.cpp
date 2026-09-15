#include <mpmc/flash/pr76_pt_flash_backend.hpp>

#include "../pt_split/dew_limit_references.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <exception>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <string_view>
#include <vector>

namespace {
namespace fl = mpmc::flash;
namespace th = mpmc::thermodynamics;
using Ref = dew_limit_reference::State;
using Vec = std::vector<double>;

constexpr double temperature_k = 220.0;
constexpr double feed_co2 = 0.3;
constexpr double eps = std::numeric_limits<double>::epsilon();

void require(bool condition, std::string_view message) {
    if (!condition) { throw std::runtime_error(std::string(message)); }
}

void near(double actual, long double expected, long double tolerance,
          std::string_view message) {
    if (!std::isfinite(actual) ||
        std::abs(static_cast<long double>(actual) - expected) > tolerance) {
        throw std::runtime_error(std::string(message));
    }
}

const Ref& named(std::string_view name) {
    for (const auto& state : dew_limit_reference::states) {
        if (state.name == name) { return state; }
    }
    throw std::logic_error("missing frozen near-dew reference state");
}

Vec compose(double first, bool reverse) {
    return reverse ? Vec{1.0 - first, first} : Vec{first, 1.0 - first};
}

// Exact test-side copy of the already-attributed PR11 CO2/methane parameter set.
// This regression does not fit parameters or derive expectations from production.
th::Pr76Phase<double> phase_model(bool reverse) {
    const th::Provenance source{
        th::SourceKind::literature,
        "https://academicweb.nd.edu/~markst/zm97a.pdf",
        "revised October 1997",
        "Sections 4.1/4.4 and Table 5; unchanged PR11 parameters",
        "PR76 numerical audit, not experimental validation",
        "PR11 attributed parameter contract",
        "Limited factual parameters; no third-party code reproduced"};
    const auto scalar = [&](double value, th::Unit unit) {
        return th::SourcedScalar{
            value, unit, source,
            unit == th::Unit::pascal ? "bar" : "SI or dimensionless",
            unit == th::Unit::pascal ? "bar * 100000 -> Pa" : "identity"};
    };
    const std::array<std::string, 2> ids{"carbon-dioxide", "methane"};
    const std::array<std::array<double, 3>, 2> specs{{
        {304.2, 7.38e6, 0.225},
        {190.6, 4.60e6, 0.008}}};

    th::PrParameterInput input;
    input.model_id = th::pr76_profile;
    input.dataset_id = "Hua1997-CO2-methane";
    input.revision = "dew-limit-audit-v1";
    input.applicability = {std::nullopt, std::nullopt, source};
    std::vector<th::Component> catalog;
    for (std::size_t i = 0; i < ids.size(); ++i) {
        catalog.push_back({ids[i], ids[i], th::ComponentKind::pure, source, {}});
        input.pure.push_back({
            ids[i], scalar(specs[i][0], th::Unit::kelvin),
            scalar(specs[i][1], th::Unit::pascal),
            scalar(specs[i][2], th::Unit::dimensionless)});
    }
    input.binary.push_back({
        ids[0], ids[1], scalar(0.095, th::Unit::dimensionless)});
    const std::vector<std::string> order = reverse
        ? std::vector<std::string>{ids[1], ids[0]}
        : std::vector<std::string>{ids[0], ids[1]};
    return th::Pr76Phase<double>::from_parameters(
        th::PrParameterSet::create(catalog, order, input));
}

fl::Pr76PtMax3Options max3_options(const fl::Pr76PtFlashBackendOptions& options) {
    fl::Pr76PtMax3Options result;
    result.two_phase = options.split;
    result.three_phase = options.three_phase;
    result.final_three_phase_stability = options.final_three_phase_stability;
    result.max_three_phase_attempts = options.max_three_phase_attempts;
    result.new_phase_seed_fraction = options.new_phase_seed_fraction;
    result.three_phase_starts = options.three_phase_starts;
    result.max_three_phase_starts = options.max_three_phase_starts;
    result.max_three_phase_start_entries = options.max_three_phase_start_entries;
    return result;
}

void require_transition_contract(const fl::PtFlashBackendCapability& capability) {
    const auto* edge = capability.transition_capability.edge(2U, 1U);
    require(edge != nullptr, "PR76 backend lost the declared 2->1 transition edge");
    require(edge->support == fl::PtPhaseTransitionSupport::detection_only,
            "2->1 must remain detection-only in this regression increment");
    require(edge->requires_fresh_target_solve,
            "2->1 must not authorize phase deletion without a fresh target solve");
}

bool has_transition(const fl::PtPhaseTransitionReport& report,
                    std::size_t source, std::size_t target,
                    fl::PtPhaseTransitionResolution resolution) {
    return std::any_of(report.evidence.begin(), report.evidence.end(),
        [&](const fl::PtPhaseTransitionEvidence& evidence) {
            return evidence.source_phase_count == source &&
                   evidence.target_phase_count.has_value() &&
                   *evidence.target_phase_count == target &&
                   evidence.resolution == resolution;
        });
}

void require_no_accepted_2_to_1(const fl::PtPhaseTransitionReport& report) {
    for (const auto& evidence : report.evidence) {
        if (evidence.source_phase_count == 2U &&
            evidence.target_phase_count == 1U) {
            require(evidence.resolution != fl::PtPhaseTransitionResolution::accepted_target,
                    "near-dew regression fabricated an accepted 2->1 transition");
            require(!evidence.target_topology_closed,
                    "detection-only 2->1 evidence cannot close the target topology");
        }
    }
}

void require_single_phase(const fl::PtFlashBackendResult& result,
                          std::span<const double> feed) {
    require(result.solution.status == fl::PtPhaseSetStatus::accepted,
            "single-phase near-dew state was not published as accepted");
    const auto* accepted = result.accepted_phase_set();
    require(accepted != nullptr && accepted->phases.size() == 1U,
            "single-phase publication has the wrong accepted phase count");
    const auto& phase = accepted->phases.front();
    require(phase.mole_phase_fraction == 1.0,
            "single-phase publication changed the unit phase fraction");
    require(phase.composition.size() == feed.size(),
            "single-phase composition dimension mismatch");
    for (std::size_t i = 0; i < feed.size(); ++i) {
        require(phase.composition[i] == feed[i],
                "single-phase publication changed the normalized feed");
    }
    require(result.transition_report.evidence.empty(),
            "stable single-phase state fabricated a phase-transition event");
}

void require_two_phase(const fl::PtFlashBackendResult& result,
                       const Ref& reference, std::span<const double> feed,
                       bool reverse) {
    require(result.solution.status == fl::PtPhaseSetStatus::accepted,
            "above-gate near-dew state was not published as accepted");
    const auto* accepted = result.accepted_phase_set();
    require(accepted != nullptr && accepted->phases.size() == 2U,
            "above-gate near-dew state did not publish exactly two phases");
    const auto& liquid = accepted->phases[0];
    const auto& vapor = accepted->phases[1];
    require(liquid.composition.size() == feed.size() &&
            vapor.composition.size() == feed.size(),
            "two-phase composition dimension mismatch");
    require(liquid.activity.ln_phi.size() == feed.size() &&
            vapor.activity.ln_phi.size() == feed.size(),
            "two-phase activity dimension mismatch");
    require(liquid.compressibility_factor.has_value() &&
            vapor.compressibility_factor.has_value(),
            "two-phase publication lost PR76 compressibility factors");
    require(liquid.mole_phase_fraction > 0.0 && vapor.mole_phase_fraction > 0.0,
            "accepted two-phase result contains a nonpositive phase fraction");
    require(std::abs(liquid.mole_phase_fraction + vapor.mole_phase_fraction - 1.0) <=
                64.0 * eps,
            "accepted two-phase fractions do not sum to one");

    for (std::size_t i = 0; i < feed.size(); ++i) {
        const double recovered = liquid.mole_phase_fraction * liquid.composition[i] +
                                 vapor.mole_phase_fraction * vapor.composition[i];
        require(std::abs(recovered - feed[i]) <= 1e-12,
                "public two-phase result violates absolute material balance");
        require(std::abs(recovered - feed[i]) / feed[i] <= 1e-10,
                "public two-phase result violates relative material balance");
        const double residual =
            (std::log(liquid.composition[i]) + liquid.activity.ln_phi[i]) -
            (std::log(vapor.composition[i]) + vapor.activity.ln_phi[i]);
        require(std::abs(residual) <= 1e-11 + 32.0 * eps,
                "public two-phase result violates fugacity balance");
    }

    const std::size_t co2 = reverse ? 1U : 0U;
    near(liquid.composition[co2], reference.x, 2e-9L,
         "public liquid composition moved outside the frozen independent budget");
    near(vapor.composition[co2], reference.y, 2e-9L,
         "public vapor composition moved outside the frozen independent budget");
    near(*liquid.compressibility_factor, reference.zl, 2e-9L,
         "public liquid Z moved outside the frozen independent budget");
    near(*vapor.compressibility_factor, reference.zv, 2e-9L,
         "public vapor Z moved outside the frozen independent budget");
    near(liquid.mole_phase_fraction, reference.liquid, 2e-11L,
         "public tiny liquid amount moved across the independent gate budget");
    require(has_transition(result.transition_report, 1U, 2U,
                           fl::PtPhaseTransitionResolution::accepted_target),
            "accepted near-dew two-phase result lost fresh 1->2 transition evidence");
    require_no_accepted_2_to_1(result.transition_report);
}

void require_disappearance_unresolved(const fl::Pr76PtMax3Result& direct,
                                      const fl::PtFlashBackendResult& result,
                                      const Ref& reference) {
    require(reference.liquid > 0.0L && reference.liquid < 1e-10L,
            "frozen threshold case no longer contains a positive sub-gate phase");
    require(reference.tpd < -1.1e-10L,
            "frozen threshold case no longer carries the recorded negative TPD sample");
    require(direct.status == fl::Pr76PtMax3Status::indeterminate,
            "max3 promoted the sub-gate disappearance state");
    require(direct.base.solution.status == fl::PtSplitStatus::indeterminate &&
            direct.base.solution.initial_stability.status == fl::StabilityStatus::unstable,
            "sub-gate case lost its unstable-feed/indeterminate-split classification");
    require(!direct.base.solution.final_stability.has_value(),
            "sub-gate disappearance unexpectedly ran a final two-phase review");
    require(std::any_of(direct.base.solution.attempts.begin(),
                        direct.base.solution.attempts.end(),
                        [](const fl::PtSplitAttempt& attempt) {
                            return attempt.status ==
                                   fl::PtSplitAttemptStatus::phase_disappearance;
                        }),
            "sub-gate case no longer reaches the declared phase-disappearance gate");
    require(result.solution.status == fl::PtPhaseSetStatus::indeterminate,
            "public backend promoted a sub-gate disappearance state");
    require(result.accepted_phase_set() == nullptr &&
            !result.solution.candidate_phase_set.has_value(),
            "public backend exposed a diagnostic disappearance state as a phase set");
    require(has_transition(result.transition_report, 2U, 1U,
                           fl::PtPhaseTransitionResolution::target_resolve_required),
            "public backend lost 2->1 target-resolve-required evidence");
    for (const auto& evidence : result.transition_report.evidence) {
        if (evidence.source_phase_count == 2U &&
            evidence.target_phase_count == 1U) {
            require(evidence.trigger == fl::PtPhaseTransitionTrigger::phase_disappearance,
                    "2->1 near-dew evidence changed trigger semantics");
            require(!evidence.fresh_target_solve_attempted &&
                    !evidence.target_topology_closed,
                    "detection-only 2->1 evidence falsely claims a fresh closed target");
        }
    }
}

void check_state(const Ref& reference, bool reverse) {
    const auto phase = phase_model(reverse);
    fl::Pr76VleEvaluator evaluator(phase);
    const Vec feed = compose(feed_co2, reverse);
    const fl::Pr76PtFlashBackendOptions options;
    const auto direct = fl::solve_pr76_pt_max3(
        reference.p, temperature_k, feed, evaluator, max3_options(options));
    const auto published = fl::project_pr76_pt_max3_phase_set(direct);
    fl::Pr76PtFlashBackend backend(evaluator, options);
    const auto result = backend.solve({reference.p, temperature_k, feed});

    require(result.structurally_valid(), "public PR76 backend result is structurally invalid");
    require_transition_contract(result.capability);
    require(result.solution.status == published.solution.status &&
            result.solution.accepted_phase_count() ==
                published.solution.accepted_phase_count(),
            "public backend classification differs from direct max3 publication");
    require(!result.solution.global_stability_proven &&
            !result.capability.global_stability_proven,
            "finite configured searches were promoted to global stability proof");
    const std::vector<std::string> expected_ids = reverse
        ? std::vector<std::string>{"methane", "carbon-dioxide"}
        : std::vector<std::string>{"carbon-dioxide", "methane"};
    require(result.capability.component_ids == expected_ids,
            "public backend changed runtime component order");

    if (reference.region <= 0) {
        require(direct.status == fl::Pr76PtMax3Status::single_phase,
                "outside/endpoint reference did not remain max3 single phase");
        require(direct.base.solution.status ==
                    fl::PtSplitStatus::single_phase_no_instability_found,
                "outside/endpoint reference changed initial stability classification");
        require_single_phase(result, feed);
        return;
    }
    if (reference.region == 1) {
        require_disappearance_unresolved(direct, result, reference);
        return;
    }

    require(reference.liquid > 1e-10L,
            "above-gate frozen reference no longer exceeds minimum phase fraction");
    require(direct.status == fl::Pr76PtMax3Status::two_phase &&
            direct.base.solution.status ==
                fl::PtSplitStatus::two_phase_no_instability_found,
            "above-gate reference did not remain accepted max3 two phase");
    require(direct.base.solution.final_stability.has_value() &&
            direct.base.solution.final_stability->status ==
                fl::StabilityStatus::no_instability_found,
            "above-gate reference lost final phase-set stability review");
    require_two_phase(result, reference, feed, reverse);
}

void check_reference_partition() {
    const auto& target = named("target");
    require(named("target_minus_ulp").p ==
                std::nextafter(target.p, -std::numeric_limits<double>::infinity()) &&
            named("target_plus_ulp").p ==
                std::nextafter(target.p, std::numeric_limits<double>::infinity()),
            "frozen target ULP neighborhood changed");
    for (const auto& reference : dew_limit_reference::states) {
        if (reference.region == 1) {
            require(reference.liquid > 0.0L && reference.liquid < 1e-10L,
                    "region-1 fixture no longer isolates a positive sub-gate phase");
        } else if (reference.region == 2) {
            require(reference.liquid > 1e-10L,
                    "region-2 fixture no longer lies above the phase-fraction gate");
        }
    }
}

void check_budget_failure(bool reverse) {
    const auto& reference = named("target");
    const auto phase = phase_model(reverse);
    fl::Pr76VleEvaluator evaluator(phase);
    const Vec feed = compose(feed_co2, reverse);
    fl::Pr76PtFlashBackendOptions options;
    options.split.max_split_attempts = 0U;
    const auto direct = fl::solve_pr76_pt_max3(
        reference.p, temperature_k, feed, evaluator, max3_options(options));
    fl::Pr76PtFlashBackend backend(evaluator, options);
    const auto result = backend.solve({reference.p, temperature_k, feed});

    require(direct.base.solution.initial_stability.status == fl::StabilityStatus::unstable,
            "budget regression did not reach the known unstable near-dew feed");
    require(direct.base.solution.attempt_limit_reached &&
            direct.base.solution.attempts.empty() &&
            direct.status == fl::Pr76PtMax3Status::indeterminate,
            "zero split-attempt budget did not remain explicitly indeterminate");
    require(result.structurally_valid() &&
            result.solution.status == fl::PtPhaseSetStatus::indeterminate &&
            result.accepted_phase_set() == nullptr,
            "public backend fabricated a result after split-attempt budget exhaustion");
    require(result.solution.diagnostic.find("split attempt limit reached") !=
                std::string::npos,
            "public budget failure lost the bounded-attempt diagnostic");
    for (const auto& evidence : result.transition_report.evidence) {
        require(evidence.resolution != fl::PtPhaseTransitionResolution::accepted_target &&
                !evidence.target_topology_closed,
                "budget exhaustion fabricated a closed accepted transition");
    }
    require_no_accepted_2_to_1(result.transition_report);
}

} // namespace

int main() {
    try {
        check_reference_partition();
        for (bool reverse : {false, true}) {
            for (const auto& reference : dew_limit_reference::states) {
                check_state(reference, reverse);
            }
            check_budget_failure(reverse);
        }
        std::cout << "PASS PR76 public near-dew classification: single / unresolved disappearance / two-phase / budget\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << error.what() << '\n';
        return 1;
    }
}
