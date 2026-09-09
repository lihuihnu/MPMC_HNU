#include <mpmc/flash/pr76_stability.hpp>
#include <mpmc/ad/math.hpp>

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
#include <sstream>
#include <stdexcept>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

bool stability_generic_header();
bool stability_pr_header(const mpmc::thermodynamics::Pr76Phase<double>&);

namespace {
namespace fl = mpmc::flash;
namespace th = mpmc::thermodynamics;
namespace ad = mpmc::ad;
using Vec = std::vector<double>;

void require(bool condition, std::string_view message,
             std::source_location where = std::source_location::current()) {
    if (!condition) {
        throw std::runtime_error(std::string(where.file_name()) + ":" +
                                 std::to_string(where.line()) + ": " + std::string(message));
    }
}
void near(double actual, long double expected, double relative = 2e-11, double absolute = 2e-13,
          std::source_location where = std::source_location::current()) {
    if (!std::isfinite(actual) || !std::isfinite(expected) ||
        std::abs(static_cast<long double>(actual)-expected) >
        absolute+relative*std::abs(expected)) {
        std::ostringstream out;
        out << std::setprecision(24) << "actual=" << actual << " expected=" << expected;
        require(false, out.str(), where);
    }
}
template <typename Error, typename Function>
void expect_error(Function&& fn) {
    bool caught = false;
    try { fn(); } catch (const Error&) { caught = true; }
    require(caught, "expected exception missing");
}
struct Ideal {
    fl::StabilityPhase operator()(double, double, std::span<const double> w) const {
        return {Vec(w.size(), 0.0), 0, true};
    }
};
struct RegularSolution {
    double interaction;
    fl::StabilityPhase operator()(double, double, std::span<const double> w) const {
        require(w.size() == 2, "binary regular-solution fixture");
        return {{interaction*w[1]*w[1], interaction*w[0]*w[0]}, 0, true};
    }
};
long double binary_g(long double x, long double interaction) {
    const auto xlogx = [](long double a) { return a == 0 ? 0 : a*std::log(a); };
    return xlogx(x)+xlogx(1-x)+interaction*x*(1-x);
}
long double regular_tpd(long double w, long double z, long double interaction) {
    return binary_g(w,interaction)-binary_g(z,interaction)-
        (std::log(z/(1-z))+interaction*(1-2*z))*(w-z);
}
void check_witness(const fl::StabilityResult& result) {
    require(result.status == fl::StabilityStatus::unstable, "missing instability");
    require(result.lowest_sampled.has_value(), "missing witness");
    require(result.lowest_sampled->value < -result.options.tpd_tolerance-
            result.lowest_sampled->roundoff_guard, "witness does not clear error guard");
    near(std::accumulate(result.lowest_sampled->composition.begin(),
                         result.lowest_sampled->composition.end(),0.0),1.0);
    require(!result.global_stability_proven, "finite search cannot prove global stability");
}
struct Spec { double tc, pc, omega; };
// Artificial data only, shared in meaning (not implementation) with the PT fixtures.
constexpr std::array<Spec,4> artificial_specs{{{400,4e6,-.125},{500,3e6,.25},
                                              {600,5e6,.5},{350,2.5e6,.125}}};
constexpr double artificial_pairs[4][4]={{0,.125,-.0625,.03125},{.125,0,.0625,-.125},
                                        {-.0625,.0625,0,.25},{.03125,-.125,.25,0}};
th::PrParameterSet fixture(const std::vector<Spec>& specs, const Vec& interactions,
                           const std::vector<std::size_t>& order, bool literature = false,
                           std::string dataset = "manufactured-stability",
                           const std::vector<std::string>& names = {}) {
    const th::Provenance source{
        literature ? th::SourceKind::literature : th::SourceKind::synthetic_test,
        literature ? "https://academicweb.nd.edu/~markst/zm97a.pdf" : "test://pt-stability",
        literature ? "author manuscript revised October 1997" : "v1",
        literature ? "Sections 4.3/4.4, Tables 4/5 (see dataset id)" : "artificial coefficients",
        literature ? "Numerical PR regression, not experimental validation" : "Not real fluid data",
        literature ? "Read author-hosted PDF and rendered pages" : "Constructed in tests",
        literature ? "Attributed limited factual parameters; no article text reproduced" : "Original fixture"};
    const auto datum = [&](double value, th::Unit unit) {
        return th::SourcedScalar{value,unit,source,
            literature && unit == th::Unit::pascal ? "bar" : "SI or dimensionless",
            literature && unit == th::Unit::pascal ? "bar * 100000 -> Pa" : "identity"};
    };
    std::vector<th::Component> catalog;
    std::vector<std::string> ids, catalog_ids;
    if (literature) { require(names.size()==specs.size(),"literature identity metadata"); }
    th::PrParameterInput input;
    input.model_id=th::pr76_profile; input.dataset_id=std::move(dataset); input.revision="v1";
    input.applicability={std::nullopt,std::nullopt,source};
    for (std::size_t i=0;i<specs.size();++i) {
        const std::string id=literature ? "literature:"+names[i] : "fixture:"+std::to_string(i);
        catalog_ids.push_back(id);
        catalog.push_back({id,literature ? names[i] : "Artificial",
                           th::ComponentKind::pure,source,{}});
        input.pure.push_back({id,datum(specs[i].tc,th::Unit::kelvin),
            datum(specs[i].pc,th::Unit::pascal),datum(specs[i].omega,th::Unit::dimensionless)});
        for (std::size_t j=0;j<i;++j) {
            input.binary.push_back({id,catalog_ids[j],
                                    datum(interactions[i*specs.size()+j],th::Unit::dimensionless)});
        }
    }
    for (std::size_t i : order) { ids.push_back(catalog_ids.at(i)); }
    return th::PrParameterSet::create(catalog,ids,input,
        literature ? th::DataPolicy::ordinary : th::DataPolicy::allow_synthetic_tests);
}
th::Pr76Phase<double> artificial(const std::vector<std::size_t>& order={0,1,2}) {
    Vec interactions;
    for (const auto& row : artificial_pairs) { interactions.insert(interactions.end(),std::begin(row),std::end(row)); }
    return th::Pr76Phase<double>::from_parameters(fixture(
        {artificial_specs.begin(),artificial_specs.end()},interactions,order));
}
th::Pr76Phase<double> literature_model(bool carbon) {
    const std::vector<Spec> specs=carbon ? std::vector<Spec>{{304.2,7.38e6,.225},{190.6,4.60e6,.008}}
                                        : std::vector<Spec>{{126.2,3.39e6,.04},{305.4,4.88e6,.098}};
    const double interaction=carbon ? .095 : .08;
    return th::Pr76Phase<double>::from_parameters(fixture(specs,{0,interaction,interaction,0},{0,1},true,
        carbon ? "Hua1997-Table5-CO2-methane" : "Hua1997-Table4-nitrogen-ethane",
        carbon ? std::vector<std::string>{"carbon-dioxide","methane"} :
                 std::vector<std::string>{"nitrogen","ethane"}));
}

void ideal_distance() {
    const Vec z{.25,.5,.25}, w{.5,.25,.25};
    const auto point=fl::tangent_plane_distance(w,z,Ideal{}(0,0,w),Ideal{}(0,0,z));
    near(point.value,.25L*std::log(2.0L));
    const auto zero=fl::tangent_plane_distance(z,z,Ideal{}(0,0,z),Ideal{}(0,0,z));
    require(zero.value==0 && zero.stationarity==0,"trivial TPD must be exactly zero");
    const Vec boundary{1,0,0};
    const auto endpoint=fl::tangent_plane_distance(boundary,z,Ideal{}(0,0,z),Ideal{}(0,0,z));
    near(endpoint.value,std::log(4.0L));
    require(std::isinf(endpoint.stationarity),"boundary is not smooth stationary point");
}
void ideal_search() {
    for (const Vec& z : {Vec{.2,.3,.5},Vec{1.0-1e-12,1e-12},Vec{1.0,1e-250}}) {
        const auto result=fl::test_pt_stability(1e5,300,z,Ideal{});
        require(result.status==fl::StabilityStatus::no_instability_found,"ideal search incomplete");
        require(!result.global_stability_proven,"global proof flag");
        for (const auto& trial : result.trials) {
            require(trial.point && trial.point->stationarity<=result.options.stationarity_tolerance,
                    "unweighted stationarity not achieved");
        }
    }
}
void regular_distance() {
    const Vec z{.4,.6},w{.1,.9};
    const RegularSolution model{3};
    const auto point=fl::tangent_plane_distance(w,z,model(0,0,w),model(0,0,z));
    near(point.value,regular_tpd(.1L,.4L,3));
    near(point.residual[0]-point.residual[1],
         std::log(.1L/.9L)+3*(1-.2L)-(std::log(.4L/.6L)+3*(1-.8L)));
}
void regular_search() {
    const Vec z{.5,.5};
    const auto unstable=fl::test_pt_stability(1e5,300,z,RegularSolution{3});
    check_witness(unstable);
    near(unstable.lowest_sampled->value,regular_tpd(unstable.lowest_sampled->composition[0],.5L,3));
    const auto stable=fl::test_pt_stability(1e5,300,z,RegularSolution{1});
    require(stable.status==fl::StabilityStatus::no_instability_found,"convex regular solution");
}
void trivial_is_not_proof() {
    fl::StabilityOptions options; options.automatic_starts=false;
    const std::vector<Vec> starts{{.5,.5}};
    const auto result=fl::test_pt_stability(1e5,300,starts[0],RegularSolution{3},options,starts);
    require(result.status==fl::StabilityStatus::no_instability_found,"trivial stationary point");
    require(!result.global_stability_proven,"trivial saddle must not be a proof");
    // Independent explicit negative witness demonstrates the limited-search counterexample.
    require(regular_tpd(.05L,.5L,3)<-.1L,"counterexample lost");
}
void zero_trace() {
    const Vec z{.7,0,.3};
    const auto result=fl::test_pt_stability(1e5,300,z,Ideal{});
    require(result.status==fl::StabilityStatus::no_instability_found,"zero-support search");
    for (const auto& trial : result.trials) {
        require(trial.point && trial.point->composition[1]==0,"absent component was introduced");
    }
    const Vec tiny{1,1e-250}, wrong{1,1e-200};
    const auto point=fl::tangent_plane_distance(wrong,tiny,Ideal{}(0,0,tiny),Ideal{}(0,0,tiny));
    require(point.stationarity>100,"trace residual was hidden by mole-fraction weighting");
    expect_error<std::domain_error>([&]{
        (void)fl::tangent_plane_distance(Vec{.7,.01,.29},z,Ideal{}(0,0,z),Ideal{}(0,0,z));
    });
}
void input_domains() {
    std::size_t calls=0;
    const auto provider=[&](double p,double t,std::span<const double> w) { ++calls; return Ideal{}(p,t,w); };
    for (const Vec& invalid : {Vec{-.1,1.1},Vec{.4,.4},Vec{std::numeric_limits<double>::quiet_NaN(),1}}) {
        expect_error<std::domain_error>([&]{ (void)fl::test_pt_stability(1e5,300,invalid,provider); });
    }
    expect_error<std::domain_error>([&]{ (void)fl::test_pt_stability(0,300,Vec{.5,.5},provider); });
    expect_error<std::domain_error>([&]{ (void)fl::test_pt_stability(1e5,-1,Vec{.5,.5},provider); });
    const std::vector<Vec> bad{{.5,.5},{1,0}};
    expect_error<std::domain_error>([&]{ (void)fl::test_pt_stability(1e5,300,Vec{.5,.5},provider,{},bad); });
    fl::StabilityOptions options; options.armijo=1;
    expect_error<std::invalid_argument>([&]{ (void)fl::test_pt_stability(1e5,300,Vec{.5,.5},provider,options); });
    require(calls==0,"invalid input reached provider");
    const Vec rounded{.5,.5+4*std::numeric_limits<double>::epsilon()};
    const auto normalized=fl::test_pt_stability(1e5,300,rounded,provider);
    require(normalized.input_feed_sum!=1,"input correction not recorded");
    near(std::accumulate(normalized.feed.begin(),normalized.feed.end(),0.0),1);
}
void budgets() {
    fl::StabilityOptions o; o.max_iterations=0;
    auto result=fl::test_pt_stability(1e5,300,Vec{.2,.8},Ideal{},o);
    require(result.status==fl::StabilityStatus::indeterminate,"iteration exhaustion claimed stability");
    o.max_evaluations=1;
    result=fl::test_pt_stability(1e5,300,Vec{.2,.8},Ideal{},o);
    require(result.evaluations==1 && result.status==fl::StabilityStatus::indeterminate,"global budget");
    o.max_starts=1;
    expect_error<std::length_error>([&]{ (void)fl::test_pt_stability(1e5,300,Vec{.2,.8},Ideal{},o); });
    o={}; o.max_start_entries=2;
    expect_error<std::length_error>([&]{ (void)fl::test_pt_stability(1e5,300,Vec{.2,.8},Ideal{},o); });
    o={}; o.automatic_starts=false;
    expect_error<std::length_error>([&]{ (void)fl::test_pt_stability(1e5,300,Vec{.2,.8},Ideal{},o); });
}
void failure_semantics() {
    const auto fails=[](double,double,std::span<const double>)->fl::StabilityPhase {
        throw fl::StabilityPropertyError(fl::StabilityPropertyIssue::root_topology,"injected topology failure");
    };
    const auto failed=fl::test_pt_stability(1e5,300,Vec{.5,.5},fails);
    require(!failed.reference && failed.reference_issue==fl::StabilityPropertyIssue::root_topology &&
            failed.status==fl::StabilityStatus::indeterminate,"reference failure lost");
    const auto mixed=[](double p,double t,std::span<const double> w) {
        if (w[0]>.9) { throw fl::StabilityPropertyError(fl::StabilityPropertyIssue::root_range,"injected"); }
        return RegularSolution{3}(p,t,w);
    };
    const auto negative=fl::test_pt_stability(1e5,300,Vec{.5,.5},mixed);
    check_witness(negative);
    require(std::any_of(negative.trials.begin(),negative.trials.end(),[](const auto& trial){
        return trial.status==fl::StabilityTrialStatus::property_failure;
    }),"negative evidence must coexist with other failed starts");
    const auto bad_shape=[](double,double,std::span<const double>){ return fl::StabilityPhase{{0},0,true}; };
    expect_error<std::logic_error>([&]{ (void)fl::test_pt_stability(1e5,300,Vec{.5,.5},bad_shape); });
    const auto programmer=[](double,double,std::span<const double>)->fl::StabilityPhase {
        throw std::logic_error("must propagate");
    };
    expect_error<std::logic_error>([&]{ (void)fl::test_pt_stability(1e5,300,Vec{.5,.5},programmer); });
}
void backtracking_recovery() {
    fl::StabilityOptions options; options.automatic_starts=false;
    const std::vector<Vec> starts{{.6,.4}};
    const auto provider=[](double p,double t,std::span<const double> w) {
        if (w[0]<.2) {
            throw fl::StabilityPropertyError(fl::StabilityPropertyIssue::root_range,"injected rejected step");
        }
        return RegularSolution{-10}(p,t,w);
    };
    const auto result=fl::test_pt_stability(1e5,300,Vec{.5,.5},provider,options,starts);
    require(result.status==fl::StabilityStatus::no_instability_found,"failed to recover by backtracking");
    require(result.trials[0].backtracks>=2 && result.trials[0].rejected_property_evaluations>=1,
            "backtracking path not exercised");
    require(result.trials[0].point->stationarity<=options.stationarity_tolerance,"premature stopping");
}
void underflow_guard() {
    const Vec z{1,std::numeric_limits<double>::denorm_min()};
    const auto result=fl::test_pt_stability(1e5,300,z,Ideal{});
    require(result.status==fl::StabilityStatus::indeterminate,"lost trace component claimed stability");
    require(std::any_of(result.trials.begin(),result.trials.end(),[](const auto& trial){
        return trial.status==fl::StabilityTrialStatus::unrepresentable_composition;
    }),"generated-start underflow diagnostic missing");
}
void nonsmooth() {
    const auto provider=[](double p,double t,std::span<const double> w) {
        auto phase=Ideal{}(p,t,w); phase.smooth=false; return phase;
    };
    const auto result=fl::test_pt_stability(1e5,300,Vec{.5,.5},provider);
    require(result.status==fl::StabilityStatus::indeterminate,"nonsmooth stationarity misreported");
    for (const auto& trial : result.trials) { require(trial.status==fl::StabilityTrialStatus::nonsmooth,"tie diagnostic"); }
}
void ownership_recovery() {
    const auto model=artificial();
    fl::Pr76StabilityEvaluator evaluator(model);
    const auto first=fl::test_pr76_pt_stability(1e5,700,Vec{.25,.5,.25},evaluator);
    const Vec saved=first.search.reference->ln_phi;
    expect_error<std::domain_error>([&]{ (void)fl::test_pr76_pt_stability(-1,700,Vec{.25,.5,.25},evaluator); });
    const auto second=fl::test_pr76_pt_stability(1e5,700,Vec{.5,.25,.25},evaluator);
    require(first.search.reference->ln_phi==saved,"results alias workspace");
    require(second.search.reference && first.dataset_id=="manufactured-stability" &&
            first.component_ids.size()==3,"metadata or recovery");
}
void minimum_gibbs_root() {
    const auto model=artificial({0});
    fl::Pr76StabilityEvaluator evaluator(model);
    th::Pr76PhaseWorkspace<double> work;
    for (double p : {1e5,1e6}) {
        const auto roots=model.roots_full(p,300,Vec{1},work);
        require(roots.status==th::Pr76RootStatus::success && roots.count==3,"fixture needs three roots");
        const auto phase=evaluator(p,300,Vec{1});
        require(phase.branch==(p==1e5 ? 2U : 0U),"always choosing largest root is incorrect");
        for (std::size_t k : {std::size_t{0},std::size_t{2}}) {
            const auto candidate=model.evaluate_full(p,300.0,Vec{1},k,work);
            require(phase.ln_phi[0]<=candidate.ln_phi[0]+2e-13,"wrong Gibbs branch");
        }
    }
}
void pr_gradient() {
    const auto model=artificial();
    fl::Pr76StabilityEvaluator evaluator(model);
    const Vec z{.25,.5,.25},w{.4,.2,.4};
    const auto reference=evaluator(1e6,450,z),phase=evaluator(1e6,450,w);
    require(phase.smooth,"gradient test away from root switch");
    const auto point=fl::tangent_plane_distance(w,z,phase,reference);
    using D=ad::Dual<double,2>;
    const D x=D::variable(w[0],0), y=D::variable(w[1],1);
    const std::vector<D> composition{x,y,1.0-x-y};
    th::Pr76PhaseWorkspace<D> work;
    const auto values=model.evaluate_full(D{1e6},D{450},composition,phase.branch,work);
    D distance{0};
    for (std::size_t i=0;i<3;++i) {
        distance+=composition[i]*(ad::log(composition[i])+values.ln_phi[i]-
                                   std::log(z[i])-reference.ln_phi[i]);
    }
    near(distance.value(),point.value);
    for (std::size_t i=0;i<2;++i) { near(distance.derivative(i),point.residual[i]-point.residual[2],2e-10,5e-12); }
}
void pr_failure() {
    const auto model=artificial();
    fl::Pr76StabilityEvaluator evaluator(model,th::Pr76RootOptions{1});
    const auto result=fl::test_pr76_pt_stability(1e6,450,Vec{.25,.5,.25},evaluator);
    require(result.search.status==fl::StabilityStatus::indeterminate &&
            result.search.reference_issue==fl::StabilityPropertyIssue::root_iteration_limit,"root failure lost");
    expect_error<std::invalid_argument>([&]{ (void)fl::test_pr76_pt_stability(1e6,450,Vec{.5,.5},evaluator); });
}
void pr_permutations_shapes() {
    std::vector<std::size_t> order{0,1,2,3};
    const Vec w{.125,.25,.375,.25};
    const auto base=artificial(order);
    fl::Pr76StabilityEvaluator initial(base);
    const auto reference=initial(1e5,700,w);
    do {
        const auto model=artificial(order);
        fl::Pr76StabilityEvaluator evaluator(model);
        Vec permuted;
        for (std::size_t i : order) { permuted.push_back(w[i]); }
        const auto result=fl::test_pr76_pt_stability(1e5,700,permuted,evaluator);
        require(result.search.status==fl::StabilityStatus::no_instability_found,"permuted stable-state search");
        for (std::size_t i=0;i<4;++i) { near(result.search.reference->ln_phi[i],reference.ln_phi[order[i]]); }
    } while (std::next_permutation(order.begin(),order.end()));
    for (std::size_t n : {1U,4U,2U,3U,1U}) {
        order.resize(n); std::iota(order.begin(),order.end(),std::size_t{0});
        const auto model=artificial(order);
        fl::Pr76StabilityEvaluator evaluator(model);
        const auto result=fl::test_pr76_pt_stability(1e5,700,Vec(n,1.0/static_cast<double>(n)),evaluator);
        require(result.search.status==fl::StabilityStatus::no_instability_found &&
                result.component_ids.size()==n,"runtime component count");
    }
}
void literature_cases(bool carbon) {
    const auto model=literature_model(carbon);
    fl::Pr76StabilityEvaluator evaluator(model);
    const Vec feeds=carbon ? Vec{.1,.2,.3,.43,.6} : Vec{.1,.18,.3,.44,.6};
    // Table 4's last feed is .60; prose says .65. Deliberately follow TABLE 4.
    for (std::size_t i=0;i<feeds.size();++i) {
        const auto result=fl::test_pr76_pt_stability(carbon ? 6.08e6 : 7.6e6,
            carbon ? 220 : 270,Vec{feeds[i],1-feeds[i]},evaluator);
        if (i==0 || i+1==feeds.size()) {
            require(result.search.status==fl::StabilityStatus::no_instability_found,"literature stable-state regression");
        } else { check_witness(result.search); }
        std::cout << "feed=" << feeds[i] << " status=" << static_cast<int>(result.search.status)
                  << " evaluations=" << result.search.evaluations << '\n';
    }
}
void literature_nitrogen() { literature_cases(false); }
void literature_carbon() { literature_cases(true); }
void decimal_anchors() {
    // Independently generated from raw data + original Z equation in Decimal(80).
    constexpr long double tpd_golden[] = {
        -0.009721854065246847410657016954526935061095L,
        -0.007389434256492901608922287314822617578491L,
         0.1886752804453729352603140947670558866353L
    };
    for (std::size_t k=0;k<3;++k) {
        const auto model=k<2 ? literature_model(k==1) : artificial();
        fl::Pr76StabilityEvaluator evaluator(model);
        const double p=k==0 ? 7.6e6 : (k==1 ? 6.08e6 : 1e6);
        const double t=k==0 ? 270 : (k==1 ? 220 : 450);
        const Vec z=k==0 ? Vec{.18,.82} : (k==1 ? Vec{.2,.8} : Vec{.25,.5,.25});
        const Vec w=k==0 ? Vec{.4943,.5057} : (k==1 ? Vec{.4972,.5028} : Vec{.4,.2,.4});
        const auto point=fl::tangent_plane_distance(w,z,evaluator(p,t,w),evaluator(p,t,z));
        near(point.value,tpd_golden[k]);
    }
}
void headers() {
    const auto model=artificial({0});
    require(stability_generic_header() && stability_pr_header(model),"standalone public headers");
}

using Test=std::pair<std::string_view,void(*)()>;
constexpr Test cases[]={{"ideal_distance",ideal_distance},{"ideal_search",ideal_search},
    {"regular_distance",regular_distance},{"regular_search",regular_search},
    {"trivial_is_not_proof",trivial_is_not_proof},{"zero_trace",zero_trace},
    {"input_domains",input_domains},{"budgets",budgets},{"failure_semantics",failure_semantics},
    {"backtracking_recovery",backtracking_recovery},{"underflow_guard",underflow_guard},
    {"nonsmooth",nonsmooth},{"ownership_recovery",ownership_recovery},
    {"minimum_gibbs_root",minimum_gibbs_root},{"pr_gradient",pr_gradient},{"pr_failure",pr_failure},
    {"pr_permutations_shapes",pr_permutations_shapes},{"literature_nitrogen",literature_nitrogen},
    {"literature_carbon",literature_carbon},{"decimal_anchors",decimal_anchors},{"headers",headers}};
} // namespace
int main(int argc,char** argv) {
    try {
        if (argc!=2) { throw std::invalid_argument("exactly one test name required"); }
        for (const auto& [name,run] : cases) {
            if (name==argv[1]) { run(); std::cout << "[PASS] " << name << '\n'; return 0; }
        }
        throw std::invalid_argument("unknown test");
    } catch (const std::exception& error) { std::cerr << "[FAIL] " << error.what() << '\n'; return 1; }
}
