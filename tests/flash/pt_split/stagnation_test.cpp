#include <mpmc/flash/pr76_split.hpp>
#include "boundary_references.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <exception>
#include <iomanip>
#include <iostream>
#include <limits>
#include <numeric>
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
constexpr double stagnation_eps = std::numeric_limits<double>::epsilon();
void require(bool value, std::string_view message,
             std::source_location loc = std::source_location::current()) {
    if (!value) { throw std::runtime_error(std::string(loc.file_name())+":"+
        std::to_string(loc.line())+": "+std::string(message)); }
}
void near(double actual, long double expected, double budget, std::string_view name) {
    if (!std::isfinite(actual) || !std::isfinite(expected) ||
        std::abs(static_cast<long double>(actual)-expected) > budget) {
        std::cerr << name << " actual=" << std::setprecision(20) << actual
                  << " expected=" << expected << " budget=" << budget << '\n';
        require(false, name);
    }
}
struct RegularSolution {
    double chi;
    fl::StabilityPhase operator()(double, double, std::span<const double> w) const {
        require(w.size() == 2, "binary analytic fixture");
        return {{chi*w[1]*w[1], chi*w[0]*w[0]}, 0, true};
    }
};

// Same attributed PR76 parameter contract as PR11; no new fluid data or oracle.
th::Pr76Phase<double> co2_model(bool reverse = false) {
    const th::Provenance source{th::SourceKind::literature,
        "https://academicweb.nd.edu/~markst/zm97a.pdf", "revised October 1997",
        "Sections 4.1/4.4 and Table 5; unchanged PR11 fixture",
        "PR76 numerical regression, not experimental validation", "PR11 attributed parameter contract",
        "Limited factual parameters; no third-party code reproduced"};
    const auto scalar = [&](double v, th::Unit u) {
        return th::SourcedScalar{v, u, source, u == th::Unit::pascal ? "bar" : "SI or dimensionless",
                                u == th::Unit::pascal ? "bar * 100000 -> Pa" : "identity"};
    };
    const std::array<std::string, 2> ids{"carbon-dioxide", "methane"};
    const std::array<std::array<double, 3>, 2> specs{{{304.2, 7.38e6, .225}, {190.6, 4.60e6, .008}}};
    th::PrParameterInput input;
    input.model_id = th::pr76_profile;
    input.dataset_id = "Hua1997-CO2-methane";
    input.revision = "tpd-stagnation-regression-v1";
    input.applicability = {std::nullopt, std::nullopt, source};
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
const boundary_reference::State& state_at(long double coordinate) {
    for (const auto& state : boundary_reference::pressure_path) {
        if (state.carbon == 1 && state.coordinate == coordinate) { return state; }
    }
    throw std::logic_error("missing frozen PR11 reference state");
}
Vec composition(double first, bool reverse = false) {
    return reverse ? Vec{1-first, first} : Vec{first, 1-first};
}

// Bounded diagnostic only: one full log-descent step from each incomplete trial.
// Extra provider calls here are NOT solver calls and are not added to its budget.
// In the red run this exposes the actual old C++ iterate, rather than pretending
// an independent Python calculation is the production iteration trajectory.
template <typename Provider>
void report_search(const fl::StabilityResult& search, Provider& provider) {
    std::cout << std::setprecision(17) << "TPD p=" << search.pressure_pa
              << " status=" << static_cast<int>(search.status)
              << " evaluations=" << search.evaluations << '\n';
    fl::StabilityPhase reference;
    if (search.reference) { reference = *search.reference; }
    else if (!search.imposed_log_activity.empty()) {
        reference.ln_phi.resize(search.feed.size());
        for (std::size_t i = 0; i < search.feed.size(); ++i) {
            reference.ln_phi[i] = search.feed[i] > 0 ?
                search.imposed_log_activity[i]-std::log(search.feed[i]) : 0;
        }
    } else { return; }
    for (const auto& trial : search.trials) {
        if (trial.status == fl::StabilityTrialStatus::stationary ||
            trial.status == fl::StabilityTrialStatus::negative_tpd || !trial.point) { continue; }
        const auto& point = *trial.point;
        const double scale = std::max(1.0, point.stationarity/search.options.log_step_limit);
        Vec w(point.composition.size());
        double largest = -std::numeric_limits<double>::infinity();
        for (std::size_t i = 0; i < w.size(); ++i) {
            if (search.feed[i] > 0) {
                w[i] = std::log(point.composition[i])-point.residual[i]/scale;
                largest = std::max(largest, w[i]);
            }
        }
        double sum = 0, predicted = 0;
        for (std::size_t i = 0; i < w.size(); ++i) {
            w[i] = search.feed[i] > 0 ? std::exp(w[i]-largest) : 0;
            sum += w[i];
            predicted += point.composition[i]*point.residual[i]*point.residual[i]/scale;
        }
        for (double& value : w) { value /= sum; }
        std::cout << " trial=" << static_cast<int>(trial.status)
                  << " iterations=" << trial.iterations << " w0=" << point.composition[0]
                  << " D=" << point.value << " R=" << point.stationarity;
        try {
            const auto next = fl::tangent_plane_distance(w, search.feed,
                provider(search.pressure_pa, search.temperature_k, w), reference);
            std::cout << " next_R=" << next.stationarity << " predicted=" << predicted
                      << " guard_sum=" << point.roundoff_guard+next.roundoff_guard
                      << " scaled_contraction=" << scale*(point.stationarity-next.stationarity)/point.stationarity
                      << " step_accepted=" << fl::detail::stability_accept_step(
                            point, next, predicted, search.options.armijo, 1/scale);
        } catch (const fl::StabilityPropertyError& error) { std::cout << " diagnostic=" << error.what(); }
        std::cout << '\n';
    }
}
void stationary(const fl::StabilityResult& s) {
    require(s.status == fl::StabilityStatus::no_instability_found, "TPD search still unresolved");
    require(s.options.stationarity_tolerance == 1e-8 && !s.global_stability_proven, "stationarity contract changed");
    for (const auto& trial : s.trials) {
        require(trial.status == fl::StabilityTrialStatus::stationary && trial.point &&
                trial.point->smooth && trial.point->stationarity <= 1e-8, "unconverged trial accepted");
    }
}

void progress_guard() {
    fl::TpdPoint point, next;
    point.value = 1; point.roundoff_guard = 1e-12; point.stationarity = 1e-7;
    next = point;
    require(!fl::detail::stability_accept_step(point, next, 1e-16, 1e-4, 1), "no-op accepted");
    next.stationarity = 1.1e-7;
    require(!fl::detail::stability_accept_step(point, next, 1e-16, 1e-4, 1), "growing residual accepted");
    next.stationarity = .975e-7;
    require(fl::detail::stability_accept_step(point, next, 1e-16, 1e-4, 1), "flat but resolved residual progress rejected");
    next.value = point.value+4e-12;
    require(!fl::detail::stability_accept_step(point, next, 1e-16, 1e-4, 1), "objective guard bypassed");
    next.value = .99; next.stationarity = 2e-7;
    require(fl::detail::stability_accept_step(point, next, .01, 1e-4, 1), "resolved Armijo changed");
}
void flat_convex() {
    // g=x ln x+(1-x)ln(1-x)+chi*x*(1-x), g'' >= 4-2chi > 0.
    // In logit coordinates the local contraction at x=.5 is 1-chi/2.
    // A 0.1 contraction lower bound is false even for these strictly convex models.
    fl::StabilityOptions options; options.automatic_starts = false;
    const std::vector<Vec> starts{{.500001, .499999}};
    bool complete = true;
    for (double chi : {1.8, 1.9, 1.95}) {
        RegularSolution provider{chi};
        const auto result = fl::test_pt_stability(1e5, 300, Vec{.5, .5}, provider, options, starts);
        report_search(result, provider);
        complete = complete && result.status == fl::StabilityStatus::no_instability_found;
        if (result.status == fl::StabilityStatus::no_instability_found) {
            stationary(result);
            near(result.trials[0].point->composition[0], .5L, 3e-7, "convex minimizer");
        }
    }
    require(complete, "strictly convex low-curvature search stalled");
}
void initial_co2() {
    const auto model = co2_model();
    fl::Pr76VleEvaluator evaluator(model);
    bool complete = true;
    for (long double coordinate : {1.0L, 1.000001L, 1.001L}) {
        const auto& ref = state_at(coordinate);
        const auto search = fl::test_pt_stability(static_cast<double>(ref.p), 220, Vec{.3, .7}, evaluator);
        report_search(search, evaluator);
        complete = complete && search.status == fl::StabilityStatus::no_instability_found;
        if (search.status == fl::StabilityStatus::no_instability_found) { stationary(search); }
    }
    require(complete, "recorded initial CO2 TPD stagnation remains");
}
void common_co2() {
    const auto model = co2_model();
    fl::Pr76VleEvaluator evaluator(model);
    bool complete = true;
    for (long double coordinate : {.99L, .9999L, .999999L}) {
        const auto& ref = state_at(coordinate);
        const double p = static_cast<double>(ref.p);
        const Vec x = composition(static_cast<double>(ref.x)), y = composition(static_cast<double>(ref.y));
        const auto l = evaluator(p, 220, x, fl::PtPhaseRole::liquid_candidate);
        const auto v = evaluator(p, 220, y, fl::PtPhaseRole::vapor_candidate);
        Vec common(2);
        for (std::size_t i = 0; i < 2; ++i) {
            common[i] = std::midpoint(std::log(x[i])+l.activity.ln_phi[i], std::log(y[i])+v.activity.ln_phi[i]);
        }
        const std::vector<Vec> starts{x, y};
        const auto search = fl::test_pt_stability_against(p, 220, Vec{.3, .7}, common, evaluator, {}, starts);
        report_search(search, evaluator);
        complete = complete && search.status == fl::StabilityStatus::no_instability_found;
        if (search.status == fl::StabilityStatus::no_instability_found) {
            stationary(search);
            require(!search.reference && search.imposed_log_activity == common, "common reference changed");
        }
    }
    require(complete, "recorded common-tangent CO2 TPD stagnation remains");
}
void flash_co2() {
    bool complete = true;
    for (bool reverse : {false, true}) {
        const auto model = co2_model(reverse);
        fl::Pr76VleEvaluator evaluator(model);
        for (long double coordinate : {.99L, .9999L, .999999L}) {
            const auto& ref = state_at(coordinate);
            const double p = static_cast<double>(ref.p);
            const Vec z = composition(.3, reverse);
            // Default starts only: no independent equilibrium seed is supplied.
            const auto result = fl::solve_pr76_pt_vle(p, 220, z, evaluator);
            const auto& s = result.solution;
            require(s.equations_converged() && s.candidate() && s.final_stability, "lost converged pair");
            const auto& c = *s.candidate();
            const auto& o = s.options.iteration;
            require(o.fugacity_tolerance == 1e-11 && o.mass_absolute_tolerance == 1e-12 &&
                    o.mass_relative_tolerance == 1e-10 && o.minimum_phase_fraction == 1e-10, "defaults changed");
            require(c.fugacity_norm <= 1e-11 && c.fractions.mass_absolute <= 1e-12 &&
                    c.fractions.mass_relative <= 1e-10, "stored equations fail");
            const double beta = c.fractions.vapor_fraction;
            require(beta > 1e-10 && 1-beta > 1e-10 && c.liquid.z < c.vapor.z, "invalid two-phase pair");
            const auto l = evaluator(p, 220, c.fractions.liquid, fl::PtPhaseRole::liquid_candidate);
            const auto v = evaluator(p, 220, c.fractions.vapor, fl::PtPhaseRole::vapor_candidate);
            for (std::size_t i = 0; i < 2; ++i) {
                const double x = c.fractions.liquid[i], y = c.fractions.vapor[i];
                const long double recovered = (1-static_cast<long double>(beta))*x+static_cast<long double>(beta)*y;
                require(std::abs(recovered-z[i]) <= 1e-12 && std::abs(recovered-z[i])/z[i] <= 1e-10,
                        "fresh material balance");
                require(std::abs((std::log(x)+l.activity.ln_phi[i])-(std::log(y)+v.activity.ln_phi[i])) <=
                        1e-11+32*stagnation_eps, "fresh fugacity balance");
            }
            const std::size_t first = reverse ? 1U : 0U;
            // Exactly PR11 comparison budgets; frozen header is never rewritten.
            near(c.fractions.liquid[first], ref.x, 2e-9, "reference x");
            near(c.fractions.vapor[first], ref.y, 2e-9, "reference y");
            near(beta, ref.beta, 1e-8, "reference beta");
            near(c.liquid.z, ref.zl, 2e-9, "reference ZL");
            near(c.vapor.z, ref.zv, 2e-9, "reference ZV");
            require(result.component_ids[first] == "carbon-dioxide", "component order");
            require(s.initial_stability.status == fl::StabilityStatus::unstable &&
                    !s.global_stability_proven && !s.final_stability->reference, "stability contract");
            report_search(*s.final_stability, evaluator);
            std::cout << "full flash coordinate=" << coordinate << " beta=" << beta
                      << " fugacity=" << c.fugacity_norm << " status=" << static_cast<int>(s.status) << '\n';
            complete = complete && s.status == fl::PtSplitStatus::two_phase_no_instability_found;
            if (s.status == fl::PtSplitStatus::two_phase_no_instability_found) { stationary(*s.final_stability); }
        }
    }
    require(complete, "recorded converged pair still fails final TPD review");
}
void budgets() {
    RegularSolution provider{1.95};
    fl::StabilityOptions options; options.automatic_starts = false; options.max_iterations = 0;
    const std::vector<Vec> starts{{.500001, .499999}};
    const auto exhausted = fl::test_pt_stability(1e5, 300, Vec{.5, .5}, provider, options, starts);
    require(exhausted.status == fl::StabilityStatus::indeterminate &&
            exhausted.trials[0].status == fl::StabilityTrialStatus::iteration_limit, "iteration exhaustion accepted");
    options.max_iterations = 512; options.max_evaluations = 1;
    const auto limited = fl::test_pt_stability(1e5, 300, Vec{.5, .5}, provider, options, starts);
    require(limited.status == fl::StabilityStatus::indeterminate && limited.evaluations == 1, "evaluation budget changed");
    const auto model = co2_model();
    fl::Pr76VleEvaluator evaluator(model);
    const double p = static_cast<double>(state_at(.99L).p);
    for (int stage : {0, 1, 2}) {
        fl::PtSplitOptions o;
        if (stage == 0) { o.initial_stability.max_evaluations = 1; }
        if (stage == 1) { o.iteration.max_evaluations = 0; }
        if (stage == 2) { o.final_stability.max_evaluations = 1; }
        const auto result = fl::solve_pr76_pt_vle(p, 220, Vec{.3, .7}, evaluator, o);
        const auto& s = result.solution;
        require(s.status == fl::PtSplitStatus::indeterminate && !s.global_stability_proven, "budget exhaustion false success");
        require(s.equations_converged() == (stage == 2), "candidate lost or invented at exhausted stage");
        if (stage == 2) {
            require(s.candidate() && s.candidate()->fugacity_norm <= 1e-11 && s.final_stability &&
                    s.final_stability->status == fl::StabilityStatus::indeterminate, "final budget did not preserve candidate");
        }
    }
}
using Case = std::pair<std::string_view, void(*)()>;
constexpr Case cases[] = {{"progress_guard", progress_guard}, {"flat_convex", flat_convex},
    {"initial_co2", initial_co2}, {"common_co2", common_co2}, {"flash_co2", flash_co2}, {"budgets", budgets}};
} // namespace
int main(int argc, char** argv) {
    try {
        if (argc != 2) { throw std::invalid_argument("one stagnation test name required"); }
        for (const auto& [name, run] : cases) {
            if (name == argv[1]) { run(); std::cout << "[PASS] " << name << '\n'; return 0; }
        }
        throw std::invalid_argument("unknown stagnation test");
    } catch (const std::exception& error) { std::cerr << "[FAIL] " << error.what() << '\n'; return 1; }
}
