#include <mpmc/flash/pr76_split.hpp>
#include "dew_limit_references.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <exception>
#include <iostream>
#include <optional>
#include <source_location>
#include <span>
#include <stdexcept>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

namespace {
namespace fl = mpmc::flash;
namespace th = mpmc::thermodynamics;
using Vec = std::vector<double>;
void require(bool condition, std::string_view message,
             std::source_location where = std::source_location::current()) {
    if (!condition) { throw std::runtime_error(std::string(where.file_name())+":"+
        std::to_string(where.line())+": "+std::string(message)); }
}
bool contains(std::string_view text, std::string_view part) {
    return text.find(part) != std::string_view::npos;
}
const dew_limit_reference::State& named(std::string_view name) {
    for (const auto& ref : dew_limit_reference::states) {
        if (ref.name == name) { return ref; }
    }
    throw std::logic_error("missing unchanged PR13 reference");
}
// Same attributed parameters as PR13, used only to test reporting on real results.
th::Pr76Phase<double> model() {
    const th::Provenance source{th::SourceKind::literature,
        "https://academicweb.nd.edu/~markst/zm97a.pdf", "revised October 1997",
        "Sections 4.1/4.4 and Table 5; unchanged PR13 parameters",
        "PR76 diagnostic regression, not experimental validation", "PR13 attributed parameter contract",
        "Limited factual parameters; no third-party code reproduced"};
    const auto scalar = [&](double value, th::Unit unit) {
        return th::SourcedScalar{value, unit, source, unit == th::Unit::pascal ? "bar" : "SI or dimensionless",
            unit == th::Unit::pascal ? "bar * 100000 -> Pa" : "identity"};
    };
    const std::array<std::string, 2> ids{"carbon-dioxide", "methane"};
    const std::array<std::array<double, 3>, 2> specs{{{304.2, 7.38e6, .225}, {190.6, 4.60e6, .008}}};
    th::PrParameterInput input;
    input.model_id = th::pr76_profile; input.dataset_id = "Hua1997-CO2-methane";
    input.revision = "diagnostic-summary-v1"; input.applicability = {std::nullopt, std::nullopt, source};
    std::vector<th::Component> catalog;
    for (std::size_t i = 0; i < 2; ++i) {
        catalog.push_back({ids[i], ids[i], th::ComponentKind::pure, source, {}});
        input.pure.push_back({ids[i], scalar(specs[i][0], th::Unit::kelvin),
            scalar(specs[i][1], th::Unit::pascal), scalar(specs[i][2], th::Unit::dimensionless)});
    }
    input.binary.push_back({ids[0], ids[1], scalar(.095, th::Unit::dimensionless)});
    return th::Pr76Phase<double>::from_parameters(th::PrParameterSet::create(catalog, ids, input));
}

void attempt_outcomes() {
    // Manufactured TERMINATION RECORDS, not manufactured fluid solutions.
    constexpr std::pair<fl::PtSplitAttemptStatus, const char*> cases[] = {
        {fl::PtSplitAttemptStatus::phase_disappearance, "phase_disappearance=1"},
        {fl::PtSplitAttemptStatus::unrepresentable_seed, "unrepresentable_seed=1"},
        {fl::PtSplitAttemptStatus::no_interior_rr_root, "no_interior_rr_root=1"},
        {fl::PtSplitAttemptStatus::rr_failure, "rr_failure=1"},
        {fl::PtSplitAttemptStatus::indistinguishable_phases, "indistinguishable_phases=1"},
        {fl::PtSplitAttemptStatus::balance_failure, "balance_failure=1"},
        {fl::PtSplitAttemptStatus::property_failure, "property_failure=1"},
        {fl::PtSplitAttemptStatus::iteration_limit, "iteration_limit=1"},
        {fl::PtSplitAttemptStatus::evaluation_limit, "evaluation_limit=1"},
        {fl::PtSplitAttemptStatus::line_search_failed, "line_search_failed=1"}
    };
    fl::PtSplitResult result;
    for (const auto& [status, label] : cases) {
        (void)label;
        fl::PtSplitAttempt attempt; attempt.status = status;
        attempt.diagnostic = std::string(4096, 'Q');
        // A rejected property call is not the TERMINAL reason for every attempt.
        attempt.property_issue = fl::StabilityPropertyIssue::root_range;
        result.attempts.push_back(std::move(attempt));
    }
    result.attempts.front().point.emplace();
    result.attempts.front().point->fractions.vapor_fraction = 5e-11;
    result.split_evaluations = 17;
    const auto text = fl::detail::split_attempt_summary(result);
    for (const auto& [status, label] : cases) { (void)status; require(contains(text, label), label); }
    require(contains(text, "balance and fugacity tolerances met") &&
        contains(text, "at or below minimum_phase_fraction") &&
        contains(text, "final phase-set stability not evaluated"), "small phase meaning missing");
    require(!contains(text, "property_failure=10") && !contains(text, "QQQQ") && text.size() < 1024,
            "history or unbounded callback text copied into terminal summary");
    require(result.status == fl::PtSplitStatus::indeterminate && !result.selected_attempt &&
        !result.candidate() && !result.equations_converged() && result.split_evaluations == 17 &&
        result.attempts.front().status == fl::PtSplitAttemptStatus::phase_disappearance &&
        result.attempts.front().point->fractions.vapor_fraction == 5e-11 &&
        result.attempts.front().diagnostic.size() == 4096, "formatter mutated a result");
    std::reverse(result.attempts.begin(), result.attempts.end());
    require(fl::detail::split_attempt_summary(result) == text, "summary depends on attempt ordering");
    result.attempts.clear(); result.split_evaluations = 0;
    result.options.max_split_attempts = 0; result.options.iteration.max_evaluations = 0;
    result.attempt_limit_reached = true;
    const auto empty = fl::detail::split_attempt_summary(result);
    require(contains(empty, "attempts=0") && contains(empty, "split attempt limit reached") &&
        contains(empty, "split property-evaluation limit reached") &&
        !contains(empty, "unrepresentable_seed=") && !contains(empty, "phase_disappearance="),
        "no work was incorrectly diagnosed as initialization/phase failure");
}
void stability_outcomes() {
    constexpr std::pair<fl::StabilityTrialStatus, const char*> cases[] = {
        {fl::StabilityTrialStatus::negative_tpd, "negative_tpd=1"},
        {fl::StabilityTrialStatus::stationary, "stationary=1"},
        {fl::StabilityTrialStatus::iteration_limit, "iteration_limit=1"},
        {fl::StabilityTrialStatus::evaluation_limit, "evaluation_limit=1"},
        {fl::StabilityTrialStatus::line_search_failed, "line_search_failed=1"},
        {fl::StabilityTrialStatus::property_failure, "property_failure=1"},
        {fl::StabilityTrialStatus::nonsmooth, "nonsmooth=1"},
        {fl::StabilityTrialStatus::unrepresentable_composition, "unrepresentable_composition=1"}
    };
    fl::StabilityResult result;
    for (const auto& [status, label] : cases) {
        (void)label;
        fl::StabilityTrial trial; trial.status = status;
        trial.property_issue = fl::StabilityPropertyIssue::root_range;
        trial.diagnostic = std::string(4096, 'Q');
        result.trials.push_back(std::move(trial));
    }
    const auto text = fl::detail::split_stability_summary(result);
    for (const auto& [status, label] : cases) { (void)status; require(contains(text, label), label); }
    require(!contains(text, "property_failure=8") && !contains(text, "QQQQ") && text.size() < 1024,
            "historical property rejections counted as terminal failures");
    std::reverse(result.trials.begin(), result.trials.end());
    require(fl::detail::split_stability_summary(result) == text, "trial ordering changed summary");
    result.reference_issue = fl::StabilityPropertyIssue::root_topology;
    const auto failed_reference = fl::detail::split_stability_summary(result);
    require(contains(failed_reference, "reference property failure") &&
        contains(failed_reference, "trial searches not started") && !contains(failed_reference, "property_failure="),
        "reference failure counted placeholders as executed searches");
    require(result.reference_issue == fl::StabilityPropertyIssue::root_topology &&
        result.trials.front().diagnostic.size() == 4096, "reference/history mutated");
    result.reference_issue.reset(); result.trials.clear(); result.imposed_log_activity = {0, 0};
    require(contains(fl::detail::split_stability_summary(result), "none started") &&
        !contains(fl::detail::split_stability_summary(result), "reference property failure"),
        "imposed reference without a physical feed phase misclassified");
}
void dew_and_neighbor() {
    const auto phase = model(); fl::Pr76VleEvaluator evaluator(phase);
    for (const char* name : {"target", "above_2e-10", "outside"}) {
        const auto& ref = named(name);
        const auto result = fl::solve_pr76_pt_vle(ref.p, 220, Vec{.3, .7}, evaluator);
        const auto& s = result.solution;
        if (ref.region == 1) {
            require(s.status == fl::PtSplitStatus::indeterminate && !s.candidate() &&
                !s.equations_converged() && !s.final_stability &&
                s.initial_stability.status == fl::StabilityStatus::unstable, "dew selection changed");
            require(contains(s.diagnostic, "phase_disappearance=1") &&
                contains(s.diagnostic, "balance and fugacity tolerances met") &&
                contains(s.diagnostic, "line_search_failed=1") &&
                contains(s.diagnostic, "final phase-set stability not evaluated"), "mixed dew summary missing");
            bool found = false;
            for (const auto& a : s.attempts) {
                if (a.status != fl::PtSplitAttemptStatus::phase_disappearance) { continue; }
                require(a.point && a.point->fugacity_norm <= 1e-11 &&
                    a.point->fractions.mass_absolute <= 1e-12 && a.point->fractions.mass_relative <= 1e-10,
                    "small-phase equations changed");
                const double liquid = 1-a.point->fractions.vapor_fraction;
                require(liquid > 0 && liquid <= 1e-10 &&
                    std::abs(static_cast<long double>(liquid)-ref.liquid) <= 2e-11L && a.diagnostic.empty(),
                    "small-phase amount or lower diagnostic changed");
                found = true;
            }
            require(found, "no actual disappearance outcome behind summary");
            std::cout << s.diagnostic << '\n';
        } else if (ref.region == 2) {
            require(s.status == fl::PtSplitStatus::two_phase_no_instability_found && s.candidate() &&
                s.final_stability && s.diagnostic.empty(), "successful neighbor reporting changed");
            require(s.candidate()->fugacity_norm <= 1e-11 &&
                s.candidate()->fractions.mass_absolute <= 1e-12 &&
                s.candidate()->fractions.mass_relative <= 1e-10 &&
                1-s.candidate()->fractions.vapor_fraction > 1e-10, "neighbor numerical contract changed");
        } else {
            require(s.status == fl::PtSplitStatus::single_phase_no_instability_found &&
                s.diagnostic.empty() && !s.candidate(), "single-phase candidate reporting changed");
        }
    }
}
void limits() {
    const auto phase = model(); fl::Pr76VleEvaluator evaluator(phase);
    for (int stage = 0; stage < 6; ++stage) {
        fl::PtSplitOptions options;
        if (stage == 0) { options.initial_stability.max_evaluations = 1; }
        if (stage == 1) { options.iteration.max_evaluations = 0; }
        if (stage == 2) { options.max_split_attempts = 0; }
        if (stage == 3) { options.iteration.max_evaluations = 1; }
        if (stage == 4) { options.iteration.max_iterations = 0; }
        if (stage == 5) { options.final_stability.max_evaluations = 1; }
        std::size_t stability_calls = 0, phase_calls = 0;
        const auto stability = [&](double p, double t, std::span<const double> w) {
            ++stability_calls; return evaluator(p, t, w);
        };
        const auto property = [&](double p, double t, std::span<const double> w, fl::PtPhaseRole role) {
            ++phase_calls; return evaluator(p, t, w, role);
        };
        const auto s = fl::solve_pt_vle(named("above_2e-10").p, 220, Vec{.3, .7}, stability, property, options);
        require(s.status == fl::PtSplitStatus::indeterminate && !s.global_stability_proven &&
            s.equations_converged() == (stage == 5), "budget changed status/candidate semantics");
        require(phase_calls == s.split_evaluations && stability_calls == s.initial_stability.evaluations+
            (s.final_stability ? s.final_stability->evaluations : 0), "summary made unaccounted property calls");
        if (stage == 0) {
            require(contains(s.diagnostic, "initial feed stability is indeterminate") &&
                contains(s.diagnostic, "evaluation_limit=") && s.attempts.empty(), "initial budget not summarized");
        } else if (stage == 1 || stage == 2) {
            const char* reason = stage == 1 ? "split property-evaluation limit reached" : "split attempt limit reached";
            require(contains(s.diagnostic, reason) && contains(s.diagnostic, "attempts=0") &&
                !contains(s.diagnostic, "unrepresentable_seed=") && s.attempts.empty(), "empty work misclassified");
        } else if (stage == 3 || stage == 4) {
            require(contains(s.diagnostic, stage == 3 ? "evaluation_limit=" : "iteration_limit=") &&
                contains(s.diagnostic, "final phase-set stability not evaluated"), "split budget not summarized");
        } else {
            require(s.candidate() && s.final_stability &&
                contains(s.diagnostic, "split equations converged; final phase-set stability is indeterminate") &&
                contains(s.diagnostic, "evaluation_limit=") &&
                !contains(s.diagnostic, "reference property failure"), "final budget/candidate misclassified");
        }
        std::cout << "stage=" << stage << " " << s.diagnostic << '\n';
    }
}
void reference_failure() {
    bool fail = true;
    std::size_t stability_calls = 0, phase_calls = 0;
    const auto stability = [&](double, double, std::span<const double> w) -> fl::StabilityPhase {
        ++stability_calls;
        if (fail) { throw fl::StabilityPropertyError(fl::StabilityPropertyIssue::root_topology, std::string(4096, 'Q')); }
        return {Vec(w.size(), 0), 0, true};
    };
    const auto property = [&](double, double, std::span<const double> w, fl::PtPhaseRole) {
        ++phase_calls; return fl::PtSplitPhase{{Vec(w.size(), 0), 0, true}, 1};
    };
    const auto failed = fl::solve_pt_vle(1e5, 300, Vec{.3, .7}, stability, property);
    require(failed.status == fl::PtSplitStatus::indeterminate && failed.attempts.empty() &&
        !failed.candidate() && stability_calls == 1 && phase_calls == 0 &&
        contains(failed.diagnostic, "initial feed stability is indeterminate; reference property failure") &&
        contains(failed.diagnostic, "trial searches not started") && !contains(failed.diagnostic, "QQQQ"),
        "reference failure summary changed execution or copied callback text");
    require(failed.initial_stability.diagnostic.size() == 4096 &&
        failed.initial_stability.reference_issue == fl::StabilityPropertyIssue::root_topology,
        "lower property diagnostic discarded");
    fail = false;
    const auto recovered = fl::solve_pt_vle(1e5, 300, Vec{.3, .7}, stability, property);
    require(recovered.status == fl::PtSplitStatus::single_phase_no_instability_found &&
        recovered.diagnostic.empty() && phase_calls == 0, "old summary leaked into successful reuse");
}
using Case = std::pair<std::string_view, void(*)()>;
constexpr Case cases[] = {{"attempt_outcomes", attempt_outcomes}, {"stability_outcomes", stability_outcomes},
    {"dew_and_neighbor", dew_and_neighbor}, {"limits", limits}, {"reference_failure", reference_failure}};
} // namespace
int main(int argc, char** argv) {
    try {
        if (argc != 2) { throw std::invalid_argument("one diagnostic-summary test name required"); }
        for (const auto& [name, run] : cases) {
            if (name == argv[1]) { run(); std::cout << "[PASS] " << name << '\n'; return 0; }
        }
        throw std::invalid_argument("unknown diagnostic-summary test");
    } catch (const std::exception& error) { std::cerr << "[FAIL] " << error.what() << '\n'; return 1; }
}
