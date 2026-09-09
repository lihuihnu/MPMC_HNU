#include <mpmc/flash/pr76_split.hpp>
#include "boundary_references.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <exception>
#include <iomanip>
#include <iostream>
#include <iterator>
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
namespace ref = boundary_reference;
using Vec = std::vector<double>;
constexpr double boundary_eps = std::numeric_limits<double>::epsilon();
// Error budgets declared before production runs: the outer 1e-11 chemical
// residual is not an x/y/beta forward-error bound near saturation. References
// retain 40 decimal digits; these comparison budgets do not replace acceptance.
constexpr double reference_budget = 2e-9;
constexpr double beta_budget = 1e-8; // Independent sampled |d(beta)/d(residual)| < 500.
constexpr double bulk_budget = 2e-8;

void require(bool value, std::string_view message,
             std::source_location loc = std::source_location::current()) {
    if (!value) {
        throw std::runtime_error(std::string(loc.file_name())+":"+
            std::to_string(loc.line())+": "+std::string(message));
    }
}
void near(double actual, long double expected, double tolerance, std::string_view name) {
    if (!std::isfinite(actual) || !std::isfinite(expected) ||
        std::abs(static_cast<long double>(actual)-expected) > tolerance) {
        std::cerr << std::setprecision(20) << name << " actual=" << actual
                  << " expected=" << expected << " budget=" << tolerance << '\n';
        require(false, name);
    }
}

th::Pr76Phase<double> model(bool carbon, bool reverse = false) {
    const th::Provenance source{th::SourceKind::literature,
        "https://academicweb.nd.edu/~markst/zm97a.pdf", "revised October 1997",
        carbon ? "Sections 4.1/4.4 and Table 5" : "Section 4.3 and Table 4",
        "PR76 numerical boundary regression, not experimental validation", "Read author-hosted PDF pages",
        "Limited attributed factual parameters; no third-party code reproduced"};
    const auto scalar = [&](double v, th::Unit u) {
        return th::SourcedScalar{v, u, source, u == th::Unit::pascal ? "bar" : "SI or dimensionless",
                                u == th::Unit::pascal ? "bar * 100000 -> Pa" : "identity"};
    };
    const std::array<std::string, 2> ids = carbon ? std::array<std::string, 2>{"carbon-dioxide", "methane"}
                                                 : std::array<std::string, 2>{"nitrogen", "ethane"};
    const std::array<std::array<double, 3>, 2> specs = carbon
        ? std::array<std::array<double, 3>, 2>{{{304.2, 7.38e6, .225}, {190.6, 4.60e6, .008}}}
        : std::array<std::array<double, 3>, 2>{{{126.2, 3.39e6, .04}, {305.4, 4.88e6, .098}}};
    th::PrParameterInput input;
    input.model_id = th::pr76_profile;
    input.dataset_id = carbon ? "Hua1997-CO2-methane" : "Hua1997-nitrogen-ethane";
    input.revision = "boundary-regression-v1";
    input.applicability = {std::nullopt, std::nullopt, source};
    std::vector<th::Component> catalog;
    for (std::size_t i = 0; i < 2; ++i) {
        catalog.push_back({ids[i], ids[i], th::ComponentKind::pure, source, {}});
        input.pure.push_back({ids[i], scalar(specs[i][0], th::Unit::kelvin),
            scalar(specs[i][1], th::Unit::pascal), scalar(specs[i][2], th::Unit::dimensionless)});
    }
    input.binary.push_back({ids[0], ids[1], scalar(carbon ? .095 : .08, th::Unit::dimensionless)});
    const std::vector<std::string> order = reverse ? std::vector<std::string>{ids[1], ids[0]}
                                                  : std::vector<std::string>{ids[0], ids[1]};
    return th::Pr76Phase<double>::from_parameters(th::PrParameterSet::create(catalog, order, input));
}

const ref::Tie& tie(bool carbon) { return ref::ties[carbon ? 1 : 0]; }
Vec feed_at(const ref::Tie& a, long double beta) {
    const double z = static_cast<double>((1-beta)*a.x+beta*a.y);
    return {z, 1-z};
}
Vec reference_log_k(const ref::Tie& a) {
    return {static_cast<double>(std::log(a.y/a.x)),
            static_cast<double>(std::log((1-a.y)/(1-a.x)))};
}
std::vector<Vec> reference_starts(const ref::Tie& a) {
    const double x = static_cast<double>(a.x), y = static_cast<double>(a.y);
    return {{x, 1-x}, {y, 1-y}};
}

// Re-evaluate actual returned states through fresh PR calls, independently of
// the stored norms/flags. Both the production criterion and fresh residual are
// checked; 32eps only accounts for the different grouping of the fresh arithmetic.
void inspect_pair(const fl::PtSplitState& c, std::span<const double> z,
                  double p, double t, fl::Pr76VleEvaluator& eval,
                  const fl::PtSplitIterationOptions& o, bool accepted) {
    const double beta = c.fractions.vapor_fraction;
    require(c.fractions.status == fl::RachfordRiceStatus::interior, "candidate RR not interior");
    require(std::isfinite(beta) && beta > 0 && beta < 1, "invalid nonzero fraction");
    require(c.fractions.liquid.size() == 2 && c.fractions.vapor.size() == 2, "dimension");
    near(std::accumulate(c.fractions.liquid.begin(), c.fractions.liquid.end(), 0.0), 1, 64*boundary_eps, "sum x");
    near(std::accumulate(c.fractions.vapor.begin(), c.fractions.vapor.end(), 0.0), 1, 64*boundary_eps, "sum y");
    require(c.fugacity_norm <= o.fugacity_tolerance, "stored fugacity acceptance");
    require(c.fractions.mass_absolute <= o.mass_absolute_tolerance &&
            c.fractions.mass_relative <= o.mass_relative_tolerance, "stored mass acceptance");
    const auto liquid = eval(p, t, c.fractions.liquid, fl::PtPhaseRole::liquid_candidate);
    const auto vapor = eval(p, t, c.fractions.vapor, fl::PtPhaseRole::vapor_candidate);
    near(liquid.z, c.liquid.z, 32*boundary_eps, "fresh liquid Z");
    near(vapor.z, c.vapor.z, 32*boundary_eps, "fresh vapor Z");
    double contrast = 0;
    for (std::size_t i = 0; i < 2; ++i) {
        const double x = c.fractions.liquid[i], y = c.fractions.vapor[i];
        require(std::isfinite(x) && std::isfinite(y) && x > 0 && y > 0 && x <= 1 && y <= 1, "composition range");
        const long double recovered = (1-static_cast<long double>(beta))*x+static_cast<long double>(beta)*y;
        const long double error = std::abs(recovered-z[i]);
        require(error <= o.mass_absolute_tolerance && error/z[i] <= o.mass_relative_tolerance,
                "fresh independent material balance");
        const double mu_l = std::log(x)+liquid.activity.ln_phi[i];
        const double mu_v = std::log(y)+vapor.activity.ln_phi[i];
        require(std::isfinite(mu_l-mu_v) && std::abs(mu_l-mu_v) <= o.fugacity_tolerance+32*boundary_eps,
                "fresh unweighted chemical residual");
        contrast = std::max(contrast, std::abs(std::log(y)-std::log(x)));
    }
    if (accepted) {
        require(beta > o.minimum_phase_fraction && 1-beta > o.minimum_phase_fraction, "disappearing phase accepted");
        require(contrast > o.log_k_separation && vapor.z-liquid.z >
                o.relative_z_separation*std::max(vapor.z, liquid.z), "trivial pair accepted");
    }
}

// One failed trial per stage is enough to locate an unresolved path without
// logging every inner iteration. These diagnostics do not alter acceptance.
void report_unresolved_stage(std::string_view name, const fl::StabilityResult& result) {
    if (result.status != fl::StabilityStatus::indeterminate) { return; }
    std::cout << "unresolved " << name << " evaluations=" << result.evaluations
              << " reason=" << result.diagnostic << '\n';
    for (const auto& trial : result.trials) {
        if (trial.status == fl::StabilityTrialStatus::stationary ||
            trial.status == fl::StabilityTrialStatus::negative_tpd) { continue; }
        std::cout << "  trial_status=" << static_cast<int>(trial.status)
                  << " iterations=" << trial.iterations << " residual="
                  << (trial.point ? trial.point->stationarity : -1.0)
                  << " reason=" << trial.diagnostic << '\n';
        break;
    }
}

struct Observation {
    fl::PtSplitStatus status;
    bool equations;
    double beta{std::numeric_limits<double>::quiet_NaN()};
    double bulk_z{std::numeric_limits<double>::quiet_NaN()};
};
Observation inspect(const fl::Pr76PtSplitResult& r, double p, double t, const Vec& z,
                    fl::Pr76VleEvaluator& eval, long double expected_beta,
                    long double x, long double y, long double zl, long double zv,
                    bool must_two, bool must_not_two) {
    const auto& s = r.solution;
    const auto& o = s.options.iteration;
    require(!s.global_stability_proven && !s.initial_stability.global_stability_proven, "global proof invented");
    require(o.fugacity_tolerance == 1e-11 && o.mass_absolute_tolerance == 1e-12 &&
            o.mass_relative_tolerance == 1e-10 && o.minimum_phase_fraction == 1e-10,
            "baseline accuracy or disappearance defaults changed");
    require(s.initial_stability.feed.size() == z.size(), "missing original feed");
    for (std::size_t i = 0; i < z.size(); ++i) {
        near(s.initial_stability.feed[i], z[i], 64*boundary_eps, "feed snapshot");
    }
    require(s.equations_converged() == (s.candidate() != nullptr), "candidate/convergence disagreement");
    if (must_two) { require(s.status == fl::PtSplitStatus::two_phase_no_instability_found, "resolved interior not accepted"); }
    if (must_not_two) { require(s.status != fl::PtSplitStatus::two_phase_no_instability_found, "exterior/endpoint false two-phase success"); }
    require(s.status != fl::PtSplitStatus::phase_set_unstable, "unexpected lower branch on checked binary reference path");
    Observation observation{s.status, s.equations_converged()};
    if (const auto* c = s.candidate()) {
        inspect_pair(*c, s.initial_stability.feed, p, t, eval, o, true);
        near(c->fractions.liquid[0], x, reference_budget, "oracle x");
        near(c->fractions.vapor[0], y, reference_budget, "oracle y");
        near(c->fractions.vapor_fraction, expected_beta, beta_budget, "oracle beta");
        near(c->liquid.z, zl, reference_budget, "oracle ZL");
        near(c->vapor.z, zv, reference_budget, "oracle ZV");
        require(s.initial_stability.status == fl::StabilityStatus::unstable, "pair without initial negative evidence");
        require(s.final_stability.has_value(), "converged pair without final review");
        require(!s.final_stability->reference && s.final_stability->imposed_log_activity.size() == 2 &&
                !s.final_stability->global_stability_proven, "not a common-tangent review");
        near(s.common_reference_allowance, .5*c->fugacity_norm, 4*boundary_eps, "reference allowance");
        near(s.final_stability->options.tpd_tolerance,
             s.options.final_stability.tpd_tolerance+s.common_reference_allowance, 4*boundary_eps, "effective TPD tolerance");
        observation.beta = c->fractions.vapor_fraction;
        observation.bulk_z = (1-observation.beta)*c->liquid.z+observation.beta*c->vapor.z;
    }
    switch (s.status) {
    case fl::PtSplitStatus::two_phase_no_instability_found:
        require(s.equations_converged() && s.final_stability &&
                s.final_stability->status == fl::StabilityStatus::no_instability_found, "false final success");
        break;
    case fl::PtSplitStatus::single_phase_no_instability_found: {
        require(s.initial_stability.status == fl::StabilityStatus::no_instability_found &&
                !s.equations_converged() && s.attempts.empty() && !s.final_stability,
                "single candidate contains stale split state");
        require(s.initial_stability.reference.has_value(), "single candidate lacks reference");
        th::Pr76PhaseWorkspace<double> workspace;
        const auto values = eval.model().evaluate_full(p, t, s.initial_stability.feed,
            s.initial_stability.reference->branch, workspace);
        observation.bulk_z = values.z;
        break;
    }
    case fl::PtSplitStatus::indeterminate:
        require(!s.diagnostic.empty(), "unresolved result lacks diagnostic");
        std::cout << "unresolved split p=" << std::setprecision(17) << p
                  << " T=" << t << " reason=" << s.diagnostic << '\n';
        report_unresolved_stage("initial TPD", s.initial_stability);
        if (s.final_stability) { report_unresolved_stage("final TPD", *s.final_stability); }
        if (s.candidate()) {
            require(s.final_stability->status == fl::StabilityStatus::indeterminate,
                    "unexpected unresolved reason on reference path");
        }
        break;
    case fl::PtSplitStatus::phase_set_unstable: break;
    }
    return observation;
}

void same_observation(const Observation& a, const Observation& b) {
    require(a.status == b.status && a.equations == b.equations, "sequential reuse changed status");
    if (a.equations) { near(a.beta, b.beta, 1e-12, "reverse scan beta"); }
    if (std::isfinite(a.bulk_z)) { near(a.bulk_z, b.bulk_z, 1e-12, "reverse scan bulk Z"); }
}
void print_observation(std::string_view kind, bool carbon, long double coordinate, const Observation& o) {
    std::cout << std::setprecision(12) << kind << " carbon=" << carbon << " coordinate=" << coordinate
              << " status=" << static_cast<int>(o.status) << " equations=" << o.equations
              << " beta=" << o.beta << '\n';
}

void composition_path(bool carbon) {
    const auto& a = tie(carbon);
    const auto pr = model(carbon);
    fl::Pr76VleEvaluator eval(pr);
    // Fixed p/T; coordinate is the exact lever-rule extension, not a forced phase.
    constexpr long double coordinates[] = {-1e-3L, -1e-6L, -1e-10L, 0.0L, 1e-10L,
        1e-8L, 1e-6L, 1e-4L, .01L, .1L, .25L, .5L, .75L, .9L, .99L,
        .9999L, .999999L, .99999999L, .9999999999L, 1.0L, 1.0000000001L, 1.000001L, 1.001L};
    std::vector<Observation> saved;
    const auto run = [&](long double q) {
        const Vec z = feed_at(a, q);
        const auto r = fl::solve_pr76_pt_vle(static_cast<double>(a.p), static_cast<double>(a.t), z, eval);
        const long double beta = (static_cast<long double>(z[0])-a.x)/(a.y-a.x);
        const auto o = inspect(r, static_cast<double>(a.p), static_cast<double>(a.t), z, eval,
            beta, a.x, a.y, a.zl, a.zv, q >= .01L && q <= .99L, q <= 0 || q >= 1);
        // For unresolved near-boundary single candidates this only checks the
        // continuity budget, never upgrades their finite stability conclusion.
        if (q >= 0 && q <= 1 && std::isfinite(o.bulk_z) && (o.equations || q <= 1e-8L || q >= .99999999L)) {
            near(o.bulk_z, (1-beta)*a.zl+beta*a.zv, bulk_budget, "composition bulk-volume limit");
        }
        return o;
    };
    for (long double q : coordinates) {
        saved.push_back(run(q));
        print_observation("composition", carbon, q, saved.back());
    }
    for (std::size_t j = std::size(coordinates); j-- > 0;) { same_observation(run(coordinates[j]), saved[j]); }
    double previous = -1;
    for (const auto& observation : saved) {
        if (observation.equations) {
            require(observation.beta+2*beta_budget >= previous, "nonmonotone composition lever path");
            previous = observation.beta;
        }
    }
}

void pressure_path(bool carbon) {
    const auto pr = model(carbon);
    fl::Pr76VleEvaluator eval(pr);
    std::vector<std::pair<const ref::State*, Observation>> saved;
    const Vec z{.3, .7};
    const auto run = [&](const ref::State& a) {
        const double p = static_cast<double>(a.p), t = static_cast<double>(a.t);
        const auto r = fl::solve_pr76_pt_vle(p, t, z, eval);
        const auto o = inspect(r, p, t, z, eval, a.beta, a.x, a.y, a.zl, a.zv,
            a.coordinate >= .1L && a.coordinate <= .9L, a.coordinate <= 0 || a.coordinate >= 1);
        const long double bulk = a.beta > 0 && a.beta < 1 ? (1-a.beta)*a.zl+a.beta*a.zv : a.feed_z;
        if (std::isfinite(o.bulk_z) && (o.equations || a.coordinate <= 1e-10L || a.coordinate >= .9999999999L)) {
            near(o.bulk_z, bulk, bulk_budget, "pressure bulk-volume limit");
        }
        return o;
    };
    for (const auto& a : ref::pressure_path) {
        if (a.carbon != static_cast<int>(carbon)) { continue; }
        saved.emplace_back(&a, run(a));
        print_observation("pressure", carbon, a.coordinate, saved.back().second);
    }
    for (std::size_t j = saved.size(); j-- > 0;) { same_observation(run(*saved[j].first), saved[j].second); }
    double previous = 2;
    for (const auto& [unused, o] : saved) {
        (void)unused;
        if (o.equations) { require(o.beta <= previous+2*beta_budget, "nonmonotone pressure quality"); previous = o.beta; }
    }
}

void disappearance(bool carbon) {
    const auto& a = tie(carbon);
    const auto pr = model(carbon);
    fl::Pr76VleEvaluator eval(pr);
    const Vec lk = reference_log_k(a);
    constexpr long double small[] = {2.5e-11L, 5e-11L, 2e-10L, 4e-10L, 1e-8L, 1e-6L, 1e-3L};
    for (bool dew : {false, true}) {
        for (long double amount : small) {
            const long double beta = dew ? 1-amount : amount;
            const Vec z = feed_at(a, beta);
            const auto r = fl::iterate_pt_split(static_cast<double>(a.p), static_cast<double>(a.t), z, lk, eval);
            require(r.point.has_value(), "oracle-seeded small phase lacks iterate");
            inspect_pair(*r.point, z, static_cast<double>(a.p), static_cast<double>(a.t), eval, {}, amount > 1e-10L);
            near(r.point->fractions.vapor_fraction, beta, 256*boundary_eps, "small-phase RR beta");
            require(r.status == (amount < 1e-10L ? fl::PtSplitAttemptStatus::phase_disappearance
                                               : fl::PtSplitAttemptStatus::converged), "disappearance gate misreported");
            std::cout << "disappearance carbon=" << carbon << " dew=" << dew << " amount=" << amount
                      << " status=" << static_cast<int>(r.status) << '\n';
        }
    }
    for (long double beta : {-1e-6L, 0.0L, 1.0L, 1.000001L}) {
        const auto r = fl::iterate_pt_split(static_cast<double>(a.p), static_cast<double>(a.t), feed_at(a, beta), lk, eval);
        require(r.status == fl::PtSplitAttemptStatus::no_interior_rr_root && !r.point, "endpoint forced into two phases");
    }
}

void assisted_neighborhoods() {
    for (bool carbon : {false, true}) {
        const auto pr = model(carbon);
        fl::Pr76VleEvaluator eval(pr);
        for (const auto& a : ref::pressure_path) {
            if (a.carbon != static_cast<int>(carbon) || !(a.coordinate > 0 && a.coordinate < 1)) { continue; }
            const double x = static_cast<double>(a.x), y = static_cast<double>(a.y);
            const std::vector<Vec> starts{{x, 1-x}, {y, 1-y}};
            const Vec z{.3, .7};
            const auto r = fl::solve_pr76_pt_vle(static_cast<double>(a.p), static_cast<double>(a.t), z, eval, {}, starts);
            const auto o = inspect(r, static_cast<double>(a.p), static_cast<double>(a.t), z, eval,
                a.beta, a.x, a.y, a.zl, a.zv, a.coordinate >= .1L && a.coordinate <= .9L, false);
            // A supplied independently negative witness well beyond arithmetic
            // guards must not be discarded as an allegedly single-phase result.
            if (a.witness_tpd < -1e-8L) {
                require(r.solution.initial_stability.status == fl::StabilityStatus::unstable &&
                        o.status != fl::PtSplitStatus::single_phase_no_instability_found, "supplied resolved witness lost");
            }
            print_observation("reference-assisted", carbon, a.coordinate, o);
        }
    }
}

void ulp_endpoints() {
    for (bool carbon : {false, true}) {
        const auto pr = model(carbon);
        fl::Pr76VleEvaluator eval(pr);
        const auto& a = tie(carbon);
        for (long double endpoint : {a.x, a.y}) {
            const double center = static_cast<double>(endpoint);
            for (double z0 : {std::nextafter(center, 0.0), center, std::nextafter(center, 1.0)}) {
                const auto r = fl::solve_pr76_pt_vle(static_cast<double>(a.p), static_cast<double>(a.t), Vec{z0, 1-z0}, eval);
                require(r.solution.status != fl::PtSplitStatus::two_phase_no_instability_found &&
                        !r.solution.equations_converged() && !r.solution.final_stability,
                        "binary64 composition endpoint claims resolved pair");
                require(!r.solution.global_stability_proven, "endpoint global proof");
            }
        }
        for (const auto& endpoint : ref::saturation) {
            if (endpoint.carbon != static_cast<int>(carbon)) { continue; }
            const double center = static_cast<double>(endpoint.p);
            for (double p : {std::nextafter(center, 0.0), center, std::nextafter(center, std::numeric_limits<double>::infinity())}) {
                const auto r = fl::solve_pr76_pt_vle(p, static_cast<double>(endpoint.t), Vec{.3, .7}, eval);
                require(r.solution.status != fl::PtSplitStatus::two_phase_no_instability_found &&
                        !r.solution.equations_converged(), "binary64 pressure endpoint claims resolved pair");
                require(!r.solution.global_stability_proven, "endpoint global proof");
            }
        }
    }
}

void budget_paths() {
    for (bool carbon : {false, true}) {
        const auto pr = model(carbon);
        fl::Pr76VleEvaluator eval(pr);
        const auto& a = tie(carbon);
        const auto starts = reference_starts(a);
        for (long double beta : {.001L, .01L, .99L, .999L}) {
            const Vec z = feed_at(a, beta);
            for (int stage : {0, 1, 2}) {
                fl::PtSplitOptions options;
                if (stage == 0) { options.initial_stability.max_evaluations = 1; }
                if (stage == 1) { options.iteration.max_evaluations = 0; }
                if (stage == 2) { options.final_stability.max_evaluations = 1; }
                const auto r = fl::solve_pr76_pt_vle(static_cast<double>(a.p), static_cast<double>(a.t), z, eval, options, starts);
                const auto& s = r.solution;
                require(s.status == fl::PtSplitStatus::indeterminate && !s.global_stability_proven && !s.diagnostic.empty(),
                        "budget failure became success");
                if (stage == 0) {
                    require(s.initial_stability.status == fl::StabilityStatus::indeterminate &&
                            s.attempts.empty() && !s.candidate() && !s.final_stability, "initial budget semantics");
                } else if (stage == 1) {
                    require(s.initial_stability.status == fl::StabilityStatus::unstable && s.attempt_limit_reached &&
                            !s.candidate() && !s.final_stability, "split budget semantics");
                } else {
                    require(s.equations_converged() && s.candidate() && s.final_stability &&
                            s.final_stability->status == fl::StabilityStatus::indeterminate,
                            "final budget lost candidate or faked stability");
                    inspect_pair(*s.candidate(), s.initial_stability.feed, static_cast<double>(a.p), static_cast<double>(a.t), eval, options.iteration, true);
                    near(s.candidate()->fractions.vapor_fraction, beta, reference_budget, "budget candidate beta");
                }
            }
        }
    }
}

void component_order() {
    for (bool carbon : {false, true}) {
        const auto pr = model(carbon), reversed = model(carbon, true);
        fl::Pr76VleEvaluator first(pr), second(reversed);
        const auto& a = tie(carbon);
        for (long double beta : {.01L, .1L, .5L, .9L, .99L}) {
            const Vec z = feed_at(a, beta), zr{z[1], z[0]};
            const auto r = fl::solve_pr76_pt_vle(static_cast<double>(a.p), static_cast<double>(a.t), z, first);
            const auto q = fl::solve_pr76_pt_vle(static_cast<double>(a.p), static_cast<double>(a.t), zr, second);
            require(r.solution.status == fl::PtSplitStatus::two_phase_no_instability_found &&
                    q.solution.status == r.solution.status, "order changed resolved state");
            require(r.component_ids.size() == 2 && q.component_ids[0] == r.component_ids[1], "ordered metadata");
            const auto& c = *q.solution.candidate();
            inspect_pair(c, zr, static_cast<double>(a.p), static_cast<double>(a.t), second, {}, true);
            near(c.fractions.vapor_fraction, beta, reference_budget, "permuted beta");
            near(c.fractions.liquid[1], a.x, reference_budget, "permuted liquid");
            near(c.fractions.vapor[1], a.y, reference_budget, "permuted vapor");
        }
    }
}
} // namespace

int main(int argc, char** argv) {
    try {
        if (argc != 2) { throw std::invalid_argument("one boundary test name required"); }
        const std::string_view name = argv[1];
        if (name == "composition_nitrogen") { composition_path(false); }
        else if (name == "composition_carbon") { composition_path(true); }
        else if (name == "pressure_nitrogen") { pressure_path(false); }
        else if (name == "pressure_carbon") { pressure_path(true); }
        else if (name == "disappearance_nitrogen") { disappearance(false); }
        else if (name == "disappearance_carbon") { disappearance(true); }
        else if (name == "assisted_neighborhoods") { assisted_neighborhoods(); }
        else if (name == "ulp_endpoints") { ulp_endpoints(); }
        else if (name == "budget_paths") { budget_paths(); }
        else if (name == "component_order") { component_order(); }
        else { throw std::invalid_argument("unknown boundary test"); }
        std::cout << "[PASS] " << name << '\n';
        return 0;
    } catch (const std::exception& e) {
        std::cerr << "[FAIL] " << e.what() << '\n';
        return 1;
    }
}
