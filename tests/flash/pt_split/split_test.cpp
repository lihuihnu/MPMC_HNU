#include <mpmc/flash/pr76_split.hpp>
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

bool pt_split_headers();
namespace {
namespace fl = mpmc::flash;
namespace th = mpmc::thermodynamics;
using Vec = std::vector<double>;
void require(bool value, std::string_view message,
             std::source_location where = std::source_location::current()) {
    if (!value) { throw std::runtime_error(std::string(where.file_name())+":"+
        std::to_string(where.line())+": "+std::string(message)); }
}
void near(double a, long double b, double tolerance=2e-12) {
    require(std::isfinite(a) && std::abs(static_cast<long double>(a)-b)<=tolerance*(1+std::abs(b)), "reference mismatch");
}
template <typename Error, typename F> void expect_error(F&& f) {
    bool caught=false; try { f(); } catch (const Error&) { caught=true; }
    require(caught,"expected exception missing");
}
Vec log_k() { return {std::log(2.0),std::log(.25)}; }

// Manufactured, Gibbs-consistent ideal branches. Constant ln(phi_L)=(ln2,ln1/4)
// and ln(phi_V)=0 give exact x_A=3/7, y_A=6/7, beta=.4 for z_A=.6.
// These are mathematical branches, NOT experimental fluids.
struct IdealBranches {
    bool third{false};
    fl::PtSplitPhase operator()(double, double, std::span<const double> w, fl::PtPhaseRole role) const {
        Vec phi(w.size());
        if (role==fl::PtPhaseRole::liquid_candidate) { phi[0]=std::log(2.0); phi[1]=std::log(.25); }
        return {{std::move(phi),role==fl::PtPhaseRole::liquid_candidate ? 0U : 1U,true},
                role==fl::PtPhaseRole::liquid_candidate ? .1 : 1.0};
    }
    fl::StabilityPhase operator()(double p, double t, std::span<const double> w) const {
        auto liquid=(*this)(p,t,w,fl::PtPhaseRole::liquid_candidate).activity;
        auto vapor=(*this)(p,t,w,fl::PtPhaseRole::vapor_candidate).activity;
        const double gl=w[0]*liquid.ln_phi[0]+w[1]*liquid.ln_phi[1];
        auto best=gl<0 ? liquid : vapor;
        best.smooth=std::abs(gl)>1e-13;
        if (third) {
            const double gt=-w[0]+2*w[1];
            if (gt<std::min(gl,0.0)) { best={Vec(w.size()),2,true}; best.ln_phi[0]=-1; best.ln_phi[1]=2; }
        }
        return best;
    }
};
struct Ideal {
    fl::StabilityPhase operator()(double,double,std::span<const double> w) const { return {Vec(w.size()),0,true}; }
};
void inspect(const fl::PtSplitResult& result, std::span<const double> z) {
    if (!result.candidate()) {
        std::cerr<<"status="<<static_cast<int>(result.status)<<" "<<result.diagnostic<<'\n';
        for (const auto& a:result.attempts) {
            std::cerr<<"attempt="<<static_cast<int>(a.status)<<" n="<<a.iterations
                     <<" r="<<(a.point ? a.point->fugacity_norm : -1)<<" "<<a.diagnostic<<'\n';
        }
    }
    require(result.equations_converged() && result.candidate(),"missing converged candidate");
    const auto& c=*result.candidate();
    require(c.fugacity_norm<=result.options.iteration.fugacity_tolerance,"fugacity norm");
    require(c.fractions.mass_absolute<=result.options.iteration.mass_absolute_tolerance,"absolute mass");
    require(c.fractions.mass_relative<=result.options.iteration.mass_relative_tolerance,"relative mass");
    require(c.fractions.vapor_fraction>0 && c.fractions.vapor_fraction<1 && c.liquid.z<c.vapor.z,"distinct positive phases");
    near(fl::detail::stability_sum(c.fractions.liquid),1);
    near(fl::detail::stability_sum(c.fractions.vapor),1);
    for (std::size_t i=0;i<z.size();++i) {
        const double recovered=(1-c.fractions.vapor_fraction)*c.fractions.liquid[i]+c.fractions.vapor_fraction*c.fractions.vapor[i];
        if (z[i]==0) { require(c.fractions.liquid[i]==0 && c.fractions.vapor[i]==0,"zero support"); }
        else { require(std::abs(recovered-z[i])/z[i]<1e-10,"independent relative balance"); }
    }
    require(!result.global_stability_proven,"no global proof");
    require(result.final_stability && !result.final_stability->reference &&
        result.final_stability->imposed_log_activity.size()==z.size(),"common tangent not recorded");
    near(result.final_stability->options.tpd_tolerance,
         result.options.final_stability.tpd_tolerance+result.common_reference_allowance,1e-14);
}

th::Pr76Phase<double> literature_model(bool carbon, bool reverse=false) {
    const th::Provenance source{th::SourceKind::literature,
        "https://academicweb.nd.edu/~markst/zm97a.pdf", "revised October 1997",
        carbon ? "Sections 4.1/4.4 and Table 5" : "Section 4.3 and Table 4",
        "PR numerical regression, not experimental validation", "Read author-hosted PDF pages",
        "Limited attributed factual parameters; no article text or code reproduced"};
    const auto scalar=[&](double v,th::Unit u) {
        return th::SourcedScalar{v,u,source,u==th::Unit::pascal ? "bar" : "SI or dimensionless",
                                u==th::Unit::pascal ? "bar * 100000 -> Pa" : "identity"};
    };
    const std::array<std::string,2> ids=carbon ? std::array<std::string,2>{"carbon-dioxide","methane"}
                                             : std::array<std::string,2>{"nitrogen","ethane"};
    const std::array<std::array<double,3>,2> specs=carbon
        ? std::array<std::array<double,3>,2>{{{304.2,7.38e6,.225},{190.6,4.60e6,.008}}}
        : std::array<std::array<double,3>,2>{{{126.2,3.39e6,.04},{305.4,4.88e6,.098}}};
    std::vector<th::Component> catalog;
    th::PrParameterInput input;
    input.model_id=th::pr76_profile;
    input.dataset_id=carbon ? "Hua1997-CO2-methane" : "Hua1997-nitrogen-ethane";
    input.revision="split-regression-v1";
    input.applicability={std::nullopt,std::nullopt,source};
    for (std::size_t i=0;i<2;++i) {
        catalog.push_back({ids[i],ids[i],th::ComponentKind::pure,source,{}});
        input.pure.push_back({ids[i],scalar(specs[i][0],th::Unit::kelvin),
                             scalar(specs[i][1],th::Unit::pascal),scalar(specs[i][2],th::Unit::dimensionless)});
    }
    input.binary.push_back({ids[0],ids[1],scalar(carbon ? .095 : .08,th::Unit::dimensionless)});
    const std::vector<std::string> order=reverse ? std::vector<std::string>{ids[1],ids[0]} :
                                                   std::vector<std::string>{ids[0],ids[1]};
    return th::Pr76Phase<double>::from_parameters(th::PrParameterSet::create(catalog,order,input));
}

void rr_binary() {
    const auto r=fl::solve_rachford_rice(Vec{.6,.4},log_k());
    require(r.status==fl::RachfordRiceStatus::interior,"RR interior");
    near(r.vapor_fraction,.4L); near(r.liquid[0],3.0L/7); near(r.vapor[0],6.0L/7);
}
void rr_zero_trace() {
    for (double trace : {0.0,1e-250}) {
        const auto r=fl::solve_rachford_rice(Vec{.6,.4,trace},Vec{std::log(2.0),std::log(.25),0});
        require(r.status==fl::RachfordRiceStatus::interior,"trace RR");
        near(r.vapor_fraction,.4L);
        require(r.liquid[2]==trace && r.vapor[2]==trace,"trace floored or lost");
        require(r.mass_relative<1e-12,"trace relative balance");
    }
}
void rr_extreme() {
    const auto good=fl::solve_rachford_rice(Vec{.5,.5},Vec{700,-700});
    require(good.status==fl::RachfordRiceStatus::interior,"extreme finite K overflowed");
    near(good.vapor_fraction,.5);
    require(good.liquid[0]>0 && good.vapor[1]>0,"small phase fraction lost");
    const auto lost=fl::solve_rachford_rice(Vec{.5,.5},Vec{1000,-1000});
    require(lost.status==fl::RachfordRiceStatus::unrepresentable,"underflow hidden");
}
void rr_boundary() {
    require(fl::solve_rachford_rice(Vec{.5,.5},Vec{1,1}).status==
        fl::RachfordRiceStatus::no_resolved_interior_root,"one-sided K");
    require(fl::solve_rachford_rice(Vec{.5,.5},Vec{0,0}).status==fl::RachfordRiceStatus::degenerate,"K=1 ambiguity");
    require(fl::solve_rachford_rice(Vec{.5,.5},Vec{1e-12,-1e-12}).status==
        fl::RachfordRiceStatus::no_resolved_interior_root,"unresolved endpoint signs");
}
void rr_domains_budget() {
    expect_error<std::domain_error>([]{(void)fl::solve_rachford_rice(Vec{.4,.4},log_k());});
    expect_error<std::domain_error>([]{(void)fl::solve_rachford_rice(Vec{.6,.4},Vec{0,std::numeric_limits<double>::infinity()});});
    expect_error<std::invalid_argument>([]{(void)fl::solve_rachford_rice(Vec{1},log_k());});
    fl::RachfordRiceOptions o; o.max_iterations=0;
    require(fl::solve_rachford_rice(Vec{.6,.4},log_k(),o).status==fl::RachfordRiceStatus::iteration_limit,"RR budget");
    o.max_components=1;
    expect_error<std::length_error>([&]{(void)fl::solve_rachford_rice(Vec{.6,.4},log_k(),o);});
}
void imposed_equivalent() {
    const Vec z{.2,.8},mu{std::log(.2),std::log(.8)};
    const auto a=fl::test_pt_stability(1e5,300,z,Ideal{});
    const auto b=fl::test_pt_stability_against(1e5,300,z,mu,Ideal{});
    require(a.status==b.status && a.status==fl::StabilityStatus::no_instability_found,"reference equivalence");
    require(a.evaluations==b.evaluations+1 && !b.reference && b.imposed_log_activity==mu,"imposed accounting");
}
void imposed_common_tangent() {
    const Vec z{.6,.4},mu{std::log(6.0/7),std::log(1.0/7)};
    const auto feed=fl::test_pt_stability(1e5,300,z,IdealBranches{});
    const auto final=fl::test_pt_stability_against(1e5,300,z,mu,IdealBranches{});
    require(feed.status==fl::StabilityStatus::unstable,"original feed needs split");
    require(final.status==fl::StabilityStatus::no_instability_found,"used feed tangent instead of common tangent");
}
void imposed_domains_budget() {
    std::size_t calls=0;
    const auto provider=[&](double p,double t,std::span<const double> w){++calls; return Ideal{}(p,t,w);};
    expect_error<std::invalid_argument>([&]{(void)fl::test_pt_stability_against(1e5,300,Vec{.5,.5},Vec{},provider);});
    expect_error<std::invalid_argument>([&]{(void)fl::test_pt_stability_against(1e5,300,Vec{.5,.5},Vec{0},provider);});
    expect_error<std::domain_error>([&]{(void)fl::test_pt_stability_against(1e5,300,Vec{.5,.5},Vec{0,std::numeric_limits<double>::quiet_NaN()},provider);});
    require(calls==0,"invalid imposed reference reached provider");
    fl::StabilityOptions o; o.max_evaluations=1;
    const auto r=fl::test_pt_stability_against(1e5,300,Vec{.5,.5},Vec{std::log(.5),std::log(.5)},provider,o);
    require(r.status==fl::StabilityStatus::indeterminate && r.evaluations==1,"imposed budget");
}
void balanced_seed() {
    const Vec z{.6,.4},w{.9,.1};
    for (bool as_vapor : {true,false}) {
        const auto k=fl::detail::split_seed(z,w,as_vapor);
        require(k.has_value(),"balanced seed unavailable");
        const auto r=fl::solve_rachford_rice(z,*k);
        require(r.status==fl::RachfordRiceStatus::interior,"balanced seed RR");
        near(r.vapor_fraction,as_vapor ? .1L : .9L);
        near((as_vapor ? r.vapor : r.liquid)[0],.9L);
    }
}
void iterate_exact() {
    const auto r=fl::iterate_pt_split(1e5,300,Vec{.6,.4},log_k(),IdealBranches{});
    require(r.status==fl::PtSplitAttemptStatus::converged && r.iterations==0,"exact equilibrium");
    near(r.point->fractions.vapor_fraction,.4L);
}
void iterate_refine() {
    Vec k=log_k(); for (auto& value:k) { value*=.5; }
    const auto r=fl::iterate_pt_split(1e5,300,Vec{.6,.4},k,IdealBranches{});
    require(r.status==fl::PtSplitAttemptStatus::converged && r.iterations>0,"SSI failed to refine");
    near(r.point->fractions.vapor_fraction,.4L);
    const auto outside=fl::iterate_pt_split(1e5,300,Vec{.6,.4},Vec{2,1},IdealBranches{});
    require(outside.status==fl::PtSplitAttemptStatus::no_interior_rr_root && outside.evaluations==0,"RR failure called single phase");
}
void iterate_limits() {
    Vec k=log_k(); for (auto& value:k) { value*=.5; }
    fl::PtSplitIterationOptions o; o.max_iterations=0;
    const auto r=fl::iterate_pt_split(1e5,300,Vec{.6,.4},k,IdealBranches{},o);
    require(r.status==fl::PtSplitAttemptStatus::iteration_limit && r.point,"iteration limit lost iterate");
    o.max_evaluations=1;
    const auto e=fl::iterate_pt_split(1e5,300,Vec{.6,.4},k,IdealBranches{},o);
    require(e.status==fl::PtSplitAttemptStatus::evaluation_limit && e.evaluations==1,"property budget");
}
void iterate_callbacks() {
    const auto wrong=[](double,double,std::span<const double>,fl::PtPhaseRole){return fl::PtSplitPhase{{{0},0,true},1};};
    expect_error<std::logic_error>([&]{(void)fl::iterate_pt_split(1e5,300,Vec{.6,.4},log_k(),wrong);});
    const auto fails=[](double,double,std::span<const double>,fl::PtPhaseRole)->fl::PtSplitPhase {
        throw fl::StabilityPropertyError(fl::StabilityPropertyIssue::root_topology,"injected");
    };
    const auto r=fl::iterate_pt_split(1e5,300,Vec{.6,.4},log_k(),fails);
    require(r.status==fl::PtSplitAttemptStatus::property_failure &&
        r.property_issue==fl::StabilityPropertyIssue::root_topology,"numeric error lost");
    const auto bug=[](double,double,std::span<const double>,fl::PtPhaseRole)->fl::PtSplitPhase {throw std::domain_error("input error");};
    expect_error<std::domain_error>([&]{(void)fl::iterate_pt_split(1e5,300,Vec{.6,.4},log_k(),bug);});
}
void split_ideal() {
    const Vec z{.6,.4}; IdealBranches model;
    const auto r=fl::solve_pt_vle(1e5,300,z,model,model);
    inspect(r,z);
    require(r.status==fl::PtSplitStatus::two_phase_no_instability_found,"ideal two-phase classification");
    near(r.candidate()->fractions.vapor_fraction,.4L);
    require(r.gibbs_change<0,"split does not lower Gibbs");
    for (const auto& a:r.attempts) { require(a.witness_trial.has_value(),"seed provenance missing"); }
}
void split_zero_trace() {
    IdealBranches model;
    for (double trace : {0.0,1e-250}) {
        const Vec z{.6,.4,trace};
        const auto r=fl::solve_pt_vle(1e5,300,z,model,model);
        inspect(r,z);
        require(r.status==fl::PtSplitStatus::two_phase_no_instability_found,"trace split classification");
        if (trace>0) { require(r.candidate()->fractions.liquid[2]>0 && r.candidate()->fractions.vapor[2]>0,"trace lost"); }
    }
}
void split_single() {
    IdealBranches model;
    const auto r=fl::solve_pt_vle(1e5,300,Vec{.05,.95},model,model);
    require(r.status==fl::PtSplitStatus::single_phase_no_instability_found && r.attempts.empty() &&
        !r.equations_converged() && !r.final_stability,"forced nonzero second phase");
}
void split_initial_unknown() {
    IdealBranches model; fl::PtSplitOptions o; o.initial_stability.max_evaluations=1;
    const auto r=fl::solve_pt_vle(1e5,300,Vec{.6,.4},model,model,o);
    require(r.status==fl::PtSplitStatus::indeterminate && r.attempts.empty(),"initial unknown became single phase");
    o={}; o.max_split_attempts=0;
    const auto a=fl::solve_pt_vle(1e5,300,Vec{.6,.4},model,model,o);
    require(a.status==fl::PtSplitStatus::indeterminate && a.attempt_limit_reached,"attempt budget");
}
void split_final_unknown() {
    IdealBranches model; fl::PtSplitOptions o; o.final_stability.max_evaluations=1;
    const Vec z{.6,.4}; const auto r=fl::solve_pt_vle(1e5,300,z,model,model,o);
    inspect(r,z);
    require(r.status==fl::PtSplitStatus::indeterminate &&
        r.final_stability->status==fl::StabilityStatus::indeterminate,"final unknown lost converged candidate");
}
void split_third_phase() {
    IdealBranches model{true}; const Vec z{.6,.4};
    const auto r=fl::solve_pt_vle(1e5,300,z,model,model);
    inspect(r,z);
    require(r.status==fl::PtSplitStatus::phase_set_unstable &&
        r.final_stability->status==fl::StabilityStatus::unstable,"missed lower third branch");
}
void split_disappearance() {
    const double z=3.0/7+1e-8*(3.0/7);
    fl::PtSplitIterationOptions o; o.minimum_phase_fraction=1e-6;
    const auto r=fl::iterate_pt_split(1e5,300,Vec{z,1-z},log_k(),IdealBranches{},o);
    require(r.status==fl::PtSplitAttemptStatus::phase_disappearance && r.point,"vanishing phase clipped");
    require(r.point->fractions.vapor_fraction>0 && r.point->fractions.vapor_fraction<1e-6,"lost vanishing fraction");
}
void split_domains() {
    IdealBranches model;
    expect_error<std::domain_error>([&]{(void)fl::solve_pt_vle(0,300,Vec{.6,.4},model,model);});
    fl::PtSplitOptions o; o.iteration.fugacity_tolerance=-1;
    expect_error<std::invalid_argument>([&]{(void)fl::solve_pt_vle(1e5,300,Vec{.6,.4},model,model,o);});
    o={}; o.final_stability.max_starts=1;
    expect_error<std::length_error>([&]{(void)fl::solve_pt_vle(1e5,300,Vec{.6,.4},model,model,o);});
    const std::vector<Vec> wrong{{1,0}};
    expect_error<std::domain_error>([&]{(void)fl::solve_pt_vle(1e5,300,Vec{.6,.4},model,model,{}, {}, wrong);});
}

// Independently generated by reference_split_decimal.py, not by production SSI.
constexpr long double golden[6][5] = {
    {0.1693257626108686289478276421222509076982L,0.4763119129749667607176185293170945103570L,0.03477107151730229156295424196056042350219L,0.2616014849057086564083344263812769752280L,0.6570958110320588707234614569898785442645L},
    {0.1693257626108686289478276421222509076982L,0.4763119129749667607176185293170945103570L,0.4256681848159807216281115892714050176037L,0.2616014849057086564083344263812769752280L,0.6570958110320588707234614569898785442645L},
    {0.1693257626108686289478276421222509076982L,0.4763119129749667607176185293170945103570L,0.8817148169977722233707951611340570440555L,0.2616014849057086564083344263812769752280L,0.6570958110320588707234614569898785442645L},
    {0.4421758532513619927110949047694758656995L,0.1925394685819822170424616895370635970132L,0.9701144068886490066089960173847518001253L,0.1704076964930220960977537969820239082426L,0.4588312140428790114288087939755051237850L},
    {0.4421758532513619927110949047694758656995L,0.1925394685819822170424616895370635970132L,0.5695317749440279545264642717985890698883L,0.1704076964930220960977537969820239082426L,0.4588312140428790114288087939755051237850L},
    {0.4421758532513619927110949047694758656995L,0.1925394685819822170424616895370635970132L,0.04877435341602058681917300253657752058011L,0.1704076964930220960977537969820239082426L,0.4588312140428790114288087939755051237850L}
};
void pr_cases(bool carbon) {
    const auto model=literature_model(carbon); fl::Pr76VleEvaluator evaluator(model);
    const Vec feeds=carbon ? Vec{.2,.3,.43} : Vec{.18,.3,.44};
    for (std::size_t k=0;k<feeds.size();++k) {
        const Vec z{feeds[k],1-feeds[k]};
        const auto result=fl::solve_pr76_pt_vle(carbon ? 6.08e6 : 7.6e6,carbon ? 220 : 270,z,evaluator);
        const auto& r=result.solution; inspect(r,z);
        std::cout<<std::setprecision(17)<<"feed="<<feeds[k]<<" status="<<static_cast<int>(r.status)
                 <<" beta="<<r.candidate()->fractions.vapor_fraction<<" residual="<<r.candidate()->fugacity_norm
                 <<" evaluations="<<r.split_evaluations<<'\n';
        require(r.status==fl::PtSplitStatus::two_phase_no_instability_found,"PR final stability");
        const auto& c=*r.candidate(); const auto& g=golden[(carbon ? 3U : 0U)+k];
        near(c.fractions.liquid[0],g[0],2e-9); near(c.fractions.vapor[0],g[1],2e-9);
        near(c.fractions.vapor_fraction,g[2],2e-9); near(c.liquid.z,g[3],2e-9); near(c.vapor.z,g[4],2e-9);
        require(r.gibbs_change<0,"PR split Gibbs decrease");
    }
}
void pr_nitrogen() {pr_cases(false);}
void pr_carbon() {pr_cases(true);}
void pr_order() {
    const auto a=literature_model(false), b=literature_model(false,true);
    fl::Pr76VleEvaluator ea(a), eb(b);
    const auto r=fl::solve_pr76_pt_vle(7.6e6,270,Vec{.3,.7},ea);
    const auto s=fl::solve_pr76_pt_vle(7.6e6,270,Vec{.7,.3},eb);
    inspect(r.solution,Vec{.3,.7}); inspect(s.solution,Vec{.7,.3});
    near(r.solution.candidate()->fractions.liquid[0],s.solution.candidate()->fractions.liquid[1],2e-9);
    require(r.component_ids[0]==s.component_ids[1],"ordered snapshot identity");
}
void pr_ownership() {
    const auto model=literature_model(false); fl::Pr76VleEvaluator evaluator(model);
    const auto first=fl::solve_pr76_pt_vle(7.6e6,270,Vec{.3,.7},evaluator);
    inspect(first.solution,Vec{.3,.7}); const Vec saved=first.solution.candidate()->fractions.liquid;
    expect_error<std::domain_error>([&]{(void)fl::solve_pr76_pt_vle(-1,270,Vec{.3,.7},evaluator);});
    const auto other=fl::solve_pr76_pt_vle(7.6e6,270,Vec{.44,.56},evaluator);
    inspect(other.solution,Vec{.44,.56});
    require(first.solution.candidate()->fractions.liquid==saved && first.dataset_id=="Hua1997-nitrogen-ethane", "result aliases scratch");
}
void pr_failure() {
    const auto model=literature_model(false); fl::Pr76VleEvaluator evaluator(model,th::Pr76RootOptions{1});
    const auto r=fl::solve_pr76_pt_vle(7.6e6,270,Vec{.3,.7},evaluator);
    require(r.solution.status==fl::PtSplitStatus::indeterminate &&
        r.solution.initial_stability.reference_issue==fl::StabilityPropertyIssue::root_iteration_limit,"root failure lost");
    expect_error<std::invalid_argument>([&]{(void)fl::solve_pr76_pt_vle(7.6e6,270,Vec{1},evaluator);});
}
void headers() {require(pt_split_headers(),"self-contained headers");}
using Test=std::pair<std::string_view,void(*)()>;
constexpr Test tests[]={
    {"rr_binary",rr_binary},{"rr_zero_trace",rr_zero_trace},{"rr_extreme",rr_extreme},
    {"rr_boundary",rr_boundary},{"rr_domains_budget",rr_domains_budget},
    {"imposed_equivalent",imposed_equivalent},{"imposed_common_tangent",imposed_common_tangent},
    {"imposed_domains_budget",imposed_domains_budget},{"balanced_seed",balanced_seed},
    {"iterate_exact",iterate_exact},{"iterate_refine",iterate_refine},{"iterate_limits",iterate_limits},
    {"iterate_callbacks",iterate_callbacks},{"split_ideal",split_ideal},{"split_zero_trace",split_zero_trace},
    {"split_single",split_single},{"split_initial_unknown",split_initial_unknown},
    {"split_final_unknown",split_final_unknown},{"split_third_phase",split_third_phase},
    {"split_disappearance",split_disappearance},{"split_domains",split_domains},
    {"pr_nitrogen",pr_nitrogen},{"pr_carbon",pr_carbon},{"pr_order",pr_order},
    {"pr_ownership",pr_ownership},{"pr_failure",pr_failure},{"headers",headers}
};
} // namespace
int main(int argc,char** argv) {
    try {
        if (argc!=2) {throw std::invalid_argument("one test name required");}
        for (const auto& [name,run]:tests) {
            if (name==argv[1]) {run(); std::cout<<"[PASS] "<<name<<'\n'; return 0;}
        }
        throw std::invalid_argument("unknown test");
    } catch(const std::exception& e) {std::cerr<<"[FAIL] "<<e.what()<<'\n'; return 1;}
}
