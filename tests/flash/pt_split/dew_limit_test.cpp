#include <mpmc/flash/pr76_split.hpp>
#include "boundary_references.hpp"
#include "dew_limit_references.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <exception>
#include <iomanip>
#include <iostream>
#include <limits>
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
using Ref = dew_limit_reference::State;
constexpr double dew_eps = std::numeric_limits<double>::epsilon();
void require(bool value, std::string_view message,
             std::source_location where = std::source_location::current()) {
    if (!value) { throw std::runtime_error(std::string(where.file_name())+":"+
        std::to_string(where.line())+": "+std::string(message)); }
}
void near(double actual, long double expected, double budget, std::string_view name) {
    if (!std::isfinite(actual) || !std::isfinite(expected) ||
        std::abs(static_cast<long double>(actual)-expected) > budget) {
        std::cerr << std::setprecision(20) << name << " actual=" << actual
                  << " expected=" << expected << " budget=" << budget << '\n';
        require(false, name);
    }
}
const Ref& named(std::string_view name) {
    for (const auto& r : dew_limit_reference::states) {
        if (r.name == name) { return r; }
    }
    throw std::logic_error("missing independent near-dew reference");
}
Vec compose(double first, bool reverse = false) {
    return reverse ? Vec{1-first, first} : Vec{first, 1-first};
}
// Existing attributed CO2/methane parameters; no new fluid model or parameter fit.
th::Pr76Phase<double> model(bool reverse = false) {
    const th::Provenance source{th::SourceKind::literature,
        "https://academicweb.nd.edu/~markst/zm97a.pdf", "revised October 1997",
        "Sections 4.1/4.4 and Table 5; unchanged PR11 parameters",
        "PR76 numerical audit, not experimental validation", "PR11 attributed parameter contract",
        "Limited factual parameters; no third-party code reproduced"};
    const auto scalar = [&](double value, th::Unit unit) {
        return th::SourcedScalar{value, unit, source, unit == th::Unit::pascal ? "bar" : "SI or dimensionless",
            unit == th::Unit::pascal ? "bar * 100000 -> Pa" : "identity"};
    };
    const std::array<std::string, 2> ids{"carbon-dioxide", "methane"};
    const std::array<std::array<double, 3>, 2> specs{{{304.2, 7.38e6, .225}, {190.6, 4.60e6, .008}}};
    th::PrParameterInput input;
    input.model_id = th::pr76_profile; input.dataset_id = "Hua1997-CO2-methane";
    input.revision = "dew-limit-audit-v1"; input.applicability = {std::nullopt, std::nullopt, source};
    std::vector<th::Component> catalog;
    for (std::size_t i = 0; i < 2; ++i) {
        catalog.push_back({ids[i], ids[i], th::ComponentKind::pure, source, {}});
        input.pure.push_back({ids[i], scalar(specs[i][0], th::Unit::kelvin),
            scalar(specs[i][1], th::Unit::pascal), scalar(specs[i][2], th::Unit::dimensionless)});
    }
    input.binary.push_back({ids[0], ids[1], scalar(.095, th::Unit::dimensionless)});
    const std::vector<std::string> order = reverse ? std::vector<std::string>{ids[1], ids[0]}
                                                  : std::vector<std::string>{ids[0], ids[1]};
    return th::Pr76Phase<double>::from_parameters(th::PrParameterSet::create(catalog, order, input));
}
Vec reference_log_k(const Ref& r, bool reverse = false) {
    const Vec x = compose(static_cast<double>(r.x), reverse), y = compose(static_cast<double>(r.y), reverse);
    return {std::log(y[0])-std::log(x[0]), std::log(y[1])-std::log(x[1])};
}
void defaults(const fl::PtSplitIterationOptions& o) {
    require(o.fugacity_tolerance == 1e-11 && o.mass_absolute_tolerance == 1e-12 &&
        o.mass_relative_tolerance == 1e-10 && o.minimum_phase_fraction == 1e-10 &&
        o.max_iterations == 512 && o.max_backtracks == 24, "production acceptance or budgets changed");
}
void equations(const fl::PtSplitState& point, const Ref& ref, fl::Pr76VleEvaluator& evaluator,
               bool reverse, double liquid_budget) {
    const Vec z = compose(.3, reverse);
    const double beta = point.fractions.vapor_fraction;
    require(point.fugacity_norm <= 1e-11 && point.fractions.mass_absolute <= 1e-12 &&
        point.fractions.mass_relative <= 1e-10, "stored joint equations fail");
    const auto l = evaluator(ref.p, 220, point.fractions.liquid, fl::PtPhaseRole::liquid_candidate);
    const auto v = evaluator(ref.p, 220, point.fractions.vapor, fl::PtPhaseRole::vapor_candidate);
    for (std::size_t i = 0; i < 2; ++i) {
        const double x = point.fractions.liquid[i], y = point.fractions.vapor[i];
        const long double recovered = (1-static_cast<long double>(beta))*x+static_cast<long double>(beta)*y;
        require(std::abs(recovered-z[i]) <= 1e-12 && std::abs(recovered-z[i])/z[i] <= 1e-10,
                "fresh material balance");
        require(std::abs((std::log(x)+l.activity.ln_phi[i])-(std::log(y)+v.activity.ln_phi[i])) <=
                1e-11+32*dew_eps, "fresh fugacity balance");
    }
    const std::size_t first = reverse ? 1U : 0U;
    // Original PR11 budgets retained, with a SEPARATE gate-resolving amount check.
    near(point.fractions.liquid[first], ref.x, 2e-9, "independent x");
    near(point.fractions.vapor[first], ref.y, 2e-9, "independent y");
    near(point.liquid.z, ref.zl, 2e-9, "independent ZL");
    near(point.vapor.z, ref.zv, 2e-9, "independent ZV");
    near(beta, 1-ref.liquid, 1e-8, "original independent beta budget");
    near(1-beta, ref.liquid, liquid_budget, "gate-resolving liquid amount");
    require(point.liquid.z < point.vapor.z, "candidate density roles reversed");
}
void report(const fl::PtSplitResult& s, const Ref& ref, bool attempts) {
    std::cout << std::setprecision(17) << ref.name << " p=" << ref.p
              << " status=" << static_cast<int>(s.status)
              << " initial=" << static_cast<int>(s.initial_stability.status)
              << " attempts=" << s.attempts.size() << " diagnostic=" << s.diagnostic << '\n';
    if (!attempts) { return; }
    for (const auto& a : s.attempts) {
        std::cout << " role_as_vapor=" << a.witness_as_vapor << " attempt_status=" << static_cast<int>(a.status)
                  << " iterations=" << a.iterations;
        if (a.point) {
            std::cout << " liquid=" << 1-a.point->fractions.vapor_fraction
                      << " fugacity=" << a.point->fugacity_norm
                      << " mass=" << a.point->fractions.mass_absolute;
        }
        std::cout << " diagnostic=" << a.diagnostic << '\n';
    }
}

void reference_resolution() {
    const auto& target = named("target");
    require(named("target_minus_ulp").p == std::nextafter(target.p, -std::numeric_limits<double>::infinity()) &&
            named("target_plus_ulp").p == std::nextafter(target.p, std::numeric_limits<double>::infinity()),
            "adjacent input bit patterns changed");
    bool matched = false;
    for (const auto& old : boundary_reference::pressure_path) {
        if (old.carbon == 1 && old.coordinate == 1e-10L) {
            require(static_cast<double>(old.p) == target.p, "not the recorded PR11 input");
            require(std::abs((1-old.beta)-target.liquid) < 2e-16L, "nominal/quantized reference shift");
            matched = true;
        }
    }
    require(matched, "frozen reference missing");
    require(target.liquid > 6e-11L && target.liquid < 8e-11L && target.rr_right > 1e-11L &&
            target.sensitivity*1e-11L < (1e-10L-target.liquid)/4, "gate separation not independently resolved");
    for (const auto& r : dew_limit_reference::states) {
        const auto k = reference_log_k(r);
        const auto rr = fl::solve_rachford_rice(compose(.3), k);
        const double right = fl::detail::rr_log_sum(compose(.3), k, -1);
        near(right, r.rr_right, 1024*dew_eps, "independent RR endpoint sign");
        if (r.region <= 0) {
            require(rr.status == fl::RachfordRiceStatus::no_resolved_interior_root, "endpoint/exterior forced interior");
        } else {
            require(rr.status == fl::RachfordRiceStatus::interior && right > 4096*dew_eps,
                    "tiny phase confused with unresolved RR endpoint");
            near(1-rr.vapor_fraction, r.liquid, 1024*dew_eps, "fixed-K lever amount");
        }
    }
}
void initialization() {
    const auto phase = model(); fl::Pr76VleEvaluator evaluator(phase);
    const auto& ref = named("target");
    const auto result = fl::solve_pr76_pt_vle(ref.p, 220, compose(.3), evaluator);
    const auto& s = result.solution;
    report(s, ref, true);
    defaults(s.options.iteration);
    require(s.initial_stability.status == fl::StabilityStatus::unstable, "negative evidence lost");
    std::size_t witnessed = 0;
    for (const auto& trial : s.initial_stability.trials) {
        if (trial.status != fl::StabilityTrialStatus::negative_tpd) { continue; }
        require(trial.point && trial.point->value < -1e-10-trial.point->roundoff_guard, "invalid negative evidence");
        ++witnessed;
        for (bool as_vapor : {false, true}) {
            const auto seed = fl::detail::split_seed(s.initial_stability.feed, trial.point->composition, as_vapor);
            require(seed.has_value(), "witness initialization unrepresentable");
            const auto rr = fl::solve_rachford_rice(s.initial_stability.feed, *seed);
            require(rr.status == fl::RachfordRiceStatus::interior && rr.mass_absolute <= 1e-12 &&
                    rr.mass_relative <= 1e-10, "material-balanced witness seed failed");
            require(rr.vapor_fraction > .01 && rr.vapor_fraction < .99, "seed starts at a disappearing phase");
        }
    }
    require(witnessed > 0 && s.attempts.size() == 2*witnessed && !s.attempt_limit_reached,
            "missing witness role or exhausted initialization");
    bool disappeared = false;
    for (const auto& attempt : s.attempts) {
        if (attempt.status == fl::PtSplitAttemptStatus::phase_disappearance) {
            require(attempt.point && attempt.point->fugacity_norm <= 1e-11 &&
                    attempt.point->fractions.mass_absolute <= 1e-12, "disappearance before equations satisfied");
            const double beta = attempt.point->fractions.vapor_fraction;
            require(std::min(beta, 1-beta) <= 1e-10, "wrong disappearance gate");
            disappeared = true;
        }
    }
    require(disappeared, "recorded case was not isolated to the declared disappearance gate");
    require(s.status == fl::PtSplitStatus::indeterminate && !s.candidate() &&
            !s.equations_converged() && !s.final_stability, "disappearing phase falsely accepted");
}
void seeded_threshold() {
    for (bool reverse : {false, true}) {
        const auto phase = model(reverse); fl::Pr76VleEvaluator evaluator(phase);
        for (const auto& ref : dew_limit_reference::states) {
            if (ref.region <= 0) { continue; }
            const auto attempt = fl::iterate_pt_split(ref.p, 220, compose(.3, reverse),
                reference_log_k(ref, reverse), evaluator);
            require(attempt.point.has_value(), "reference seed was not evaluated");
            equations(*attempt.point, ref, evaluator, reverse, 1024*dew_eps);
            const auto expected = ref.region == 1 ? fl::PtSplitAttemptStatus::phase_disappearance
                                                  : fl::PtSplitAttemptStatus::converged;
            require(attempt.status == expected, "reference seed changed phase gate semantics");
            require(attempt.iterations == 0, "reference seed needs unexpected refinement");
            std::cout << ref.name << " reference-seeded status=" << static_cast<int>(attempt.status)
                      << " liquid=" << std::setprecision(17) << 1-attempt.point->fractions.vapor_fraction << '\n';
        }
    }
}
void neighbors() {
    for (bool reverse : {false, true}) {
        const auto phase = model(reverse); fl::Pr76VleEvaluator evaluator(phase);
        for (const auto& ref : dew_limit_reference::states) {
            const auto result = fl::solve_pr76_pt_vle(ref.p, 220, compose(.3, reverse), evaluator);
            const auto& s = result.solution;
            report(s, ref, false);
            defaults(s.options.iteration);
            require(!s.global_stability_proven, "finite search became a proof");
            if (ref.region <= 0) {
                require(s.status == fl::PtSplitStatus::single_phase_no_instability_found && !s.candidate(),
                        "endpoint/exterior misreported as two-phase");
            } else if (ref.region == 1) {
                require(s.initial_stability.status == fl::StabilityStatus::unstable &&
                        s.status == fl::PtSplitStatus::indeterminate && !s.candidate() && !s.final_stability,
                        "subthreshold witness falsely resolved");
            } else {
                require(s.status == fl::PtSplitStatus::two_phase_no_instability_found &&
                        s.candidate() && s.final_stability, "above-gate neighbor did not resolve");
                equations(*s.candidate(), ref, evaluator, reverse, 2e-11);
                require(1-s.candidate()->fractions.vapor_fraction > 1e-10, "accepted amount below gate");
                for (const auto& trial : s.final_stability->trials) {
                    require(trial.status == fl::StabilityTrialStatus::stationary && trial.point &&
                            trial.point->stationarity <= 1e-8, "final review silently bypassed");
                }
            }
        }
    }
}
void budgets() {
    const auto phase = model(); fl::Pr76VleEvaluator evaluator(phase);
    const auto& ref = named("target");
    fl::RachfordRiceOptions rr_options; rr_options.max_iterations = 0;
    require(fl::solve_rachford_rice(compose(.3), reference_log_k(ref), rr_options).status ==
            fl::RachfordRiceStatus::iteration_limit, "RR budget confused with physical endpoint");
    fl::PtSplitIterationOptions iteration; iteration.max_evaluations = 0;
    const auto low = fl::iterate_pt_split(ref.p, 220, compose(.3), reference_log_k(ref), evaluator, iteration);
    require(low.status == fl::PtSplitAttemptStatus::evaluation_limit && !low.point,
            "unevaluated reference seed falsely classified");
    for (int stage : {0, 1, 2}) {
        fl::PtSplitOptions options;
        if (stage == 0) { options.initial_stability.max_evaluations = 1; }
        if (stage == 1) { options.iteration.max_evaluations = 0; }
        if (stage == 2) { options.final_stability.max_evaluations = 1; }
        const double p = stage == 2 ? named("inside_1e-8").p : ref.p;
        const auto result = fl::solve_pr76_pt_vle(p, 220, compose(.3), evaluator, options);
        const auto& s = result.solution;
        require(s.status == fl::PtSplitStatus::indeterminate && !s.global_stability_proven,
                "budget exhaustion became success");
        require(s.equations_converged() == (stage == 2), "candidate invented or discarded at exhausted stage");
        if (stage == 0) { require(s.attempts.empty(), "initial failure reached split"); }
        if (stage == 1) { require(s.attempt_limit_reached, "split resource diagnostic lost"); }
        if (stage == 2) {
            require(s.candidate() && s.final_stability &&
                    s.final_stability->status == fl::StabilityStatus::indeterminate, "final candidate not retained");
        }
    }
}
using Case = std::pair<std::string_view, void(*)()>;
constexpr Case cases[] = {{"reference_resolution", reference_resolution}, {"initialization", initialization},
    {"seeded_threshold", seeded_threshold}, {"neighbors", neighbors}, {"budgets", budgets}};
} // namespace
int main(int argc, char** argv) {
    try {
        if (argc != 2) { throw std::invalid_argument("one dew-limit test name required"); }
        for (const auto& [name, run] : cases) {
            if (name == argv[1]) { run(); std::cout << "[PASS] " << name << '\n'; return 0; }
        }
        throw std::invalid_argument("unknown dew-limit test");
    } catch (const std::exception& error) { std::cerr << "[FAIL] " << error.what() << '\n'; return 1; }
}
