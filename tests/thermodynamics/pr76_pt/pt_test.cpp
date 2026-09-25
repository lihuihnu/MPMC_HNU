#include "test_support.hpp"
#include <mpmc/thermodynamics/selected_phase_fugacity.hpp>
#include <mpmc/thermodynamics/selected_phase_density.hpp>
#include <mpmc/ad/runtime_differentiate.hpp>
#include <iostream>
#include <numeric>
#include <span>
#include <type_traits>

double pt_plain_header(const mpmc::thermodynamics::Pr76Phase<double>& model);
bool roots_plain_header();

namespace {
using namespace test;

template <typename T>
void roots_factored() {
    // X=Z-B: (X-1/4)*[(X-1/4)^2-1/32], with A=51/128, B=1/16.
    const auto roots = th::pr76_roots(T{51}/T{128}, T{1}/T{16});
    require(roots.status == th::Pr76RootStatus::success && roots.count == 3, "three exact roots");
    const long double d = std::sqrt(2.0L)/8;
    const std::array<long double,3> expected{5.0L/16-d, 5.0L/16, 5.0L/16+d};
    for (std::size_t i = 0; i < 3; ++i) {
        near(roots.roots[i].z, expected[i]);
        require(roots.roots[i].derivative_valid, "simple manufactured root");
        require(roots.roots[i].slope_sign == (i == 1 ? -1 : 1), "mechanical slope sign");
    }
    const auto zero = th::pr76_roots(T{0}, T{1}/T{16});
    require(zero.status == th::Pr76RootStatus::success && zero.count == 1, "A=0 physical root");
    near(zero.roots[0].z, 17.0L/16);
}
template <typename T>
void root_topology() {
    // EXACT binary double root: H=(X-1/8)^2*(X-1/2), not a fitted fixture.
    const T a = T{49}/T{128}, b = T{1}/T{16};
    const auto multiple = th::pr76_roots(a,b);
    require(multiple.status == th::Pr76RootStatus::near_multiple && multiple.count == 0,
            "unresolved multiplicity must not invent roots");
    const T delta = T{1}/T{10000};
    require(th::pr76_roots(a+delta,b).count == 3, "split double root");
    require(th::pr76_roots(a-delta,b).count == 1, "one-root side of spinodal");
    const auto stopped = th::pr76_roots(T{51}/T{128}, b, {1});
    require(stopped.status == th::Pr76RootStatus::iteration_limit && stopped.count == 0,
            "no result after iteration limit");
    const T ac = static_cast<T>(.457235528921382189383460196225183789L);
    const T bc = static_cast<T>(.077796073903888455971844710037333184L);
    // These constants solve (1-4B)^3=54B^2 and A=3[(1-4B)/3]^2-2B^2+4B.
    const auto close = th::pr76_roots(ac+T{16}*std::numeric_limits<T>::epsilon(),bc);
    require(close.status == th::Pr76RootStatus::success && close.count == 1,
            "simple but near-critical root still has a primal");
    require(!close.roots[0].derivative_valid, "near-critical derivative must be flagged");
    expect_error<std::invalid_argument>([&] { (void)th::pr76_roots(a,T{0}); });
    expect_error<std::invalid_argument>([&] { (void)th::pr76_roots(a,b,{0}); });
    expect_error<std::invalid_argument>([&] {
        (void)th::pr76_roots(std::numeric_limits<T>::quiet_NaN(),b);
    });
}
template <typename T>
void root_scale() {
    for (const auto& [a,b] : std::array<std::array<T,2>,4>{{
             {T{0},T{1}/T{1000000000}}, {T{-1},T{1}/T{10}},
             {T{1000000},T{100}}, {T{1}/T{10000000},T{1}/T{100000000}}}}) {
        const auto roots = th::pr76_roots(a,b);
        require(roots.status == th::Pr76RootStatus::success && roots.count > 0, "scaled roots");
        for (std::size_t i=0; i<roots.count; ++i) {
            require(roots.roots[i].z>b && roots.roots[i].scaled_residual<=T{128}*std::numeric_limits<T>::epsilon(),
                    "root domain/backward residual");
        }
    }
    const auto tiny = th::pr76_roots(T{1},std::numeric_limits<T>::denorm_min());
    require(tiny.status == th::Pr76RootStatus::unrepresentable, "endpoint underflow must be diagnosed");
}
template <typename T>
void phase_values() {
    auto model=kernel<T>();
    th::Pr76PhaseWorkspace<T> work;
    const std::vector<T> w{T{1}/T{4},T{1}/T{2},T{1}/T{4}};
    for (const auto& [p,t] : std::array<std::array<T,2>,4>{{
             {T{1000000},T{450}}, {T{200000},T{350}}, {T{5000000},T{650}}, {T{10000000},T{300}}}}) {
        const auto ref=reference(p,t,w,{0,1,2});
        const auto roots=model.roots_full(p,t,w,work);
        require(roots.status==th::Pr76RootStatus::success && roots.count==ref.size(), "phase root count");
        for(std::size_t k=0;k<roots.count;++k) {
            const auto result=model.evaluate_full(p,t,w,k,work);
            near(result.z,ref[k].z);
            for(std::size_t i=0;i<w.size();++i) { near(result.ln_phi[i],ref[k].ln_phi[i]); }
        }
    }
    // Decimal(80) direct-Z Newton/log references; regeneration script is independent.
    constexpr long double golden[3][4] = {
        {.046044017300579607027738136516236690L,1.337140171088663018018446468398940819L,
         .180548816123331146991749684999174121L,-.609684880855792126056105516526145202L},
        {.087033790295512388734717465744650001L,1.037812206671159998731746601060865614L,
         .181878902503568685239764483271857615L,-.173911599152081002983672080004092468L},
        {.843005896107611707941248101442817013L,-.045959304113661024760838233711854714L,
         -.173187163838049554415087111179631632L,-.202518429356737116471648864036772263L}};
    for(std::size_t k=0;k<3;++k) {
        const auto result=model.evaluate_full(T{1000000},T{450},w,k,work);
        near(result.z,golden[k][0]);
        for(std::size_t i=0;i<3;++i) { near(result.ln_phi[i],golden[k][i+1]); }
    }
}
template <typename T>
void full_derivatives() {
    using D=ad::Dual<T,5>;
    auto model=kernel<T>();
    th::Pr76PhaseWorkspace<D> work;
    const std::vector<T> w{T{1}/T{4},T{1}/T{2},T{1}/T{4}};
    const auto ref=reference(T{1000000},T{450},w,{0,1,2});
    const std::vector<D> x{D::variable(w[0],2),D::variable(w[1],3),D::variable(w[2],4)};
    for(std::size_t k=0;k<ref.size();++k) {
        const auto result=model.evaluate_full(D::variable(T{1000000},0),D::variable(T{450},1),x,k,work);
        near(result.z.value(),ref[k].z);
        for(std::size_t c=0;c<5;++c) {
            near(result.z.derivative(c),ref[k].dz[c]);
            for(std::size_t i=0;i<3;++i) { near(result.ln_phi[i].derivative(c),ref[k].dln[i][c]); }
        }
    }
}
template <typename T>
void reduced_derivatives() {
    using D=ad::Dual<T,4>;
    auto model=kernel<T>();
    th::Pr76PhaseWorkspace<D> work;
    const std::vector<T> w{T{1}/T{4},T{1}/T{2},T{1}/T{4}};
    // Decimal(80), original Z IFT and directly differentiated Eq.(19).
    // Keep the original 4096*epsilon relative checks; do not assume long double
    // provides extra reference precision on MSVC. Inputs and all 48 checks are unchanged.
    constexpr long double reduced_golden[3][4][4] = {
        {
            {4.350357709924831983135725769138258721e-8L, 3.103581657790207862742935261560588777e-4L, 7.227691622190335238245502720444835137e-2L, 6.297587295236221273588651911130909930e-2L},
            {-9.312362320091739850533638526360815966e-7L, 1.055792512849407185698216151598953972e-3L, -8.382126288182792092465454009256223835e-1L, -6.545765224006856386566210420246526337e-2L},
            {-9.405372752787151246999323607292208487e-7L, 1.077486805379385879478566965525377425e-2L, -2.689247244605726914533949052978115007e-1L, -7.570357488599300779377029542070631794e-1L},
            {-1.003513148231077337435818879840529948e-6L, 2.335989216019409097835001331925226238e-2L, 1.376062077739424592153335211521245385e+0L, 1.579529149959928719741068012616591622e+0L},
        },
        {
            {1.097776349208289083936756798015315009e-7L, -1.467930029154822116274765512897460838e-3L, -2.521274081623674534483079021422591531e-1L, -1.082083986877770711719859514718324500e-1L},
            {-1.047957566482374665765520485126128139e-6L, 1.108357361857154893289010360824635155e-2L, 1.341570812652056546075118977164829842e+0L, 9.893663603374980817092439774319572441e-1L},
            {-9.040385570077842834891985344557014355e-7L, 8.053172606651301768450904443722985635e-3L, -4.773184757622566848374166371093645958e-3L, -3.257407311479036589442457384061827920e-1L},
            {-7.958301583200072123172125829838689855e-7L, 3.300494680615739680406485610129189961e-3L, -1.332024443136811412378370644422642550e+0L, -3.378848980416907638207525006195916602e-1L},
        },
        {
            {-1.771975083163735245213292337892103844e-7L, 1.210719188478681988436685978510949285e-3L, 1.833082697182418788436306527155885795e-1L, 3.716437758726671028795128421237520260e-2L},
            {-3.809509039734023807000455112667915368e-8L, 1.480320952572682344977012944923857755e-4L, -4.880284840622855100235414484856261370e-2L, 4.304382097686431202560219862171895979e-2L},
            {-1.842389825283154066256839196298925306e-7L, 1.259469929560761669036912893901961725e-3L, 2.164941607082419829152386921706868244e-2L, -3.497112107374229008949346722039724299e-2L},
            {-2.214033601155821169136352038422677332e-7L, 1.502529118263404061133616685546260636e-3L, 5.504016264580154419306406414425248832e-3L, 2.689842117062026815338473581907552618e-2L},
        },
    };
    const std::vector<D> x{D::variable(w[0],2),D::variable(w[1],3)};
    for(std::size_t k=0;k<3;++k) {
        const auto result=model.evaluate_reduced(D::variable(T{1000000},0),D::variable(T{450},1),x,k,work);
        for(std::size_t c=0;c<4;++c) {
            near(result.z.derivative(c),reduced_golden[k][0][c]);
            for(std::size_t i=0;i<3;++i) {
                near(result.ln_phi[i].derivative(c),reduced_golden[k][i+1][c]);
            }
        }
    }
}
template <typename T>
void gibbs_duhem() {
    using D=ad::Dual<T,2>;
    auto model=kernel<T>();
    th::Pr76PhaseWorkspace<D> work;
    const std::array<T,3> w{T{1}/T{4},T{1}/T{2},T{1}/T{4}};
    const std::vector<D> x{D::variable(w[0],0),D::variable(w[1],1)};
    for(std::size_t k=0;k<3;++k) {
        const auto result=model.evaluate_reduced(D{T{1000000}},D{T{450}},x,k,work);
        for(std::size_t c=0;c<2;++c) {
            long double sum=0,scale=0;
            for(std::size_t i=0;i<3;++i) {
                const long double term=static_cast<long double>(w[i])*result.ln_phi[i].derivative(c);
                sum+=term; scale+=std::abs(term);
            }
            require(std::abs(sum)<=8192.0L*std::numeric_limits<T>::epsilon()*scale,
                    "constrained Gibbs-Duhem identity");
        }
    }
}
template <typename T>
void pure_limit() {
    using D=ad::Dual<T,2>;
    auto model=kernel<T>({0});
    th::Pr76PhaseWorkspace<T> plain;
    th::Pr76PhaseWorkspace<D> work;
    const std::vector<T> w{T{1}};
    const auto ref=reference(T{1000000},T{350},w,{0});
    const auto roots=model.roots_reduced(T{1000000},T{350},{},plain);
    require(roots.count==ref.size(), "pure empty reduced composition");
    for(std::size_t k=0;k<roots.count;++k) {
        const auto result=model.evaluate_reduced(D::variable(T{1000000},0),D::variable(T{350},1),{},k,work);
        near(result.z.value(),ref[k].z); near(result.ln_phi[0].value(),ref[k].ln_phi[0]);
        for(std::size_t c=0;c<2;++c) {
            near(result.z.derivative(c),ref[k].dz[c]); near(result.ln_phi[0].derivative(c),ref[k].dln[0][c]);
        }
    }
}
template <typename T>
void low_pressure() {
    using D=ad::Dual<T,1>;
    auto model=kernel<T>();
    th::Pr76PhaseWorkspace<T> plain;
    th::Pr76PhaseWorkspace<D> work;
    const T p=T{1000000}*std::numeric_limits<T>::epsilon(), t=T{450};
    const std::vector<T> w{T{1}/T{4},T{1}/T{2},T{1}/T{4}};
    const auto roots=model.roots_full(p,t,w,plain);
    require(roots.status==th::Pr76RootStatus::success, "sub-epsilon gas departure");
    std::vector<D> x; for(const T v:w) { x.emplace_back(v); }
    const auto result=model.evaluate_full(D::variable(p,0),D{t},x,roots.count-1,work);
    std::array<PureReference,3> pure{};
    for(std::size_t i=0;i<3;++i) { pure[i]=pure_reference(t,specs[i]); }
    std::array<long double,3> rows{};
    long double a=0,b=0;
    for(std::size_t i=0;i<3;++i) {
        b+=static_cast<long double>(w[i])*pure[i].b;
        for(std::size_t j=0;j<3;++j) {
            rows[i]+=static_cast<long double>(w[j])*(1-static_cast<long double>(fixture_kij[i][j]))*
                     std::sqrt(pure[i].a)*std::sqrt(pure[j].a);
        }
        a+=static_cast<long double>(w[i])*rows[i];
    }
    const long double rt=8.31446261815324L*t;
    near(result.z.derivative(0),b/rt-a/(rt*rt));
    for(std::size_t i=0;i<3;++i) {
        // First virial coefficient of ln(phi); omitted O(p^2) is below this
        // test's relative rounding tolerance at p=1e6*epsilon Pa.
        const long double slope=pure[i].b/rt-(2*rows[i]-a)/(rt*rt);
        near(result.ln_phi[i].value()/p,slope);
        near(result.ln_phi[i].derivative(0),slope);
        require(result.ln_phi[i].value()!=T{0}, "small nonzero ln(phi) must survive Z rounding to 1");
    }
}
template <typename T>
void log_ratio() {
    using D=ad::Dual<T,1>;
    for(const T u:std::array<T,6>{T{0},std::numeric_limits<T>::epsilon(),
          T{1}/T{8}-T{1}/T{100000},T{1}/T{8}+T{1}/T{100000},T{1}/T{2},T{4}}) {
        const auto actual=th::detail::pr76_log1p_over_x<T>(D::variable(u,0));
        long double value=0,derivative=0,power=1;
        const long double v=u;
        if(v<T{1}/T{4}) {
            for(int k=0;k<80;++k) {
                value+=(k%2==0?1.0L:-1.0L)*power/(k+1);
                if(k>0) { derivative+=(k%2==0?1.0L:-1.0L)*k*std::pow(v,k-1)/(k+1); }
                power*=v;
            }
        } else {
            value=std::log1p(v)/v;
            derivative=(1/(1+v)-value)/v;
        }
        near(actual.value(),value); near(actual.derivative(0),derivative);
    }
}
template <typename T>
void zero_trace() {
    using D=ad::Dual<T,5>;
    auto model=kernel<T>();
    th::Pr76PhaseWorkspace<D> work;
    for(const auto& w:std::array<std::vector<T>,2>{{{T{1},T{0},T{0}},
           {T{1}-std::numeric_limits<T>::epsilon(),std::numeric_limits<T>::epsilon(),T{0}}}}) {
        const auto ref=reference(T{1000000},T{450},w,{0,1,2});
        std::vector<D> x; for(std::size_t i=0;i<3;++i) { x.push_back(D::variable(w[i],i+2)); }
        const auto k=ref.size()-1;
        const auto result=model.evaluate_full(D::variable(T{1000000},0),D::variable(T{450},1),x,k,work);
        for(std::size_t i=0;i<3;++i) {
            near(result.ln_phi[i].value(),ref[k].ln_phi[i]);
            for(std::size_t c=2;c<5;++c) {
                const T actual = result.ln_phi[i].derivative(c);
                const long double expected = ref[k].dln[i][c];
                if(i == 0) {
                    // At w=(1,0,0), differentiating Eq.(19) gives exactly zero
                    // for the present component in EVERY full-composition column:
                    // dg*e+g*de=0, and the remaining Z/B terms cancel by the EOS.
                    // With an epsilon trace the derivative is O(epsilon), while
                    // those individual terms remain O(1). A relative error divided
                    // by the cancelled sum (or reference roundoff) is undefined.
                    // Use the SAME 4096*epsilon budget on the independently
                    // derived sum of absolute analytic terms; no arbitrary floor.
                    const long double scale = ref[k].dln_absolute_terms[i][c];
                    require(std::isfinite(actual) && scale > 0 && std::isfinite(scale),
                            "finite boundary derivative and analytic scale");
                    require(std::abs(static_cast<long double>(actual)-expected) <=
                            4096.0L*std::numeric_limits<T>::epsilon()*scale,
                            "boundary derivative cancellation residual");
                    if(w[1] == T{0}) {
                        require(std::abs(expected) <=
                                4096.0L*std::numeric_limits<long double>::epsilon()*scale,
                                "independent reference must satisfy the exact pure-vertex identity");
                    }
                } else {
                    near(actual, expected); // Nonzero solute insertion derivatives.
                }
            }
        }
    }
}
template <typename T>
void permutations() {
    using D=ad::Dual<T,6>;
    const std::vector<T> original{T{1}/T{16},T{3}/T{16},T{1}/T{4},T{1}/T{2}};
    const auto ref=reference(T{1000000},T{450},original,{0,1,2,3});
    std::vector<std::size_t> order{0,1,2,3};
    th::Pr76PhaseWorkspace<D> work;
    int count=0;
    do {
        auto model=kernel<T>(order);
        std::vector<D> x;
        for(std::size_t i=0;i<4;++i) { x.push_back(D::variable(original[order[i]],i+2)); }
        for(std::size_t k=0;k<ref.size();++k) {
            const auto value=model.evaluate_full(D::variable(T{1000000},0),D::variable(T{450},1),x,k,work);
            near(value.z.value(),ref[k].z);
            for(std::size_t i=0;i<4;++i) {
                near(value.ln_phi[i].value(),ref[k].ln_phi[order[i]]);
                for(std::size_t j=0;j<4;++j) {
                    near(value.ln_phi[i].derivative(j+2),ref[k].dln[order[i]][order[j]+2]);
                }
            }
        }
        ++count;
    } while(std::next_permutation(order.begin(),order.end()));
    require(count==24,"all 24 snapshot permutations");
}
template <typename T>
void runtime_shapes() {
    using D=ad::Dual<T,2>;
    ad::RuntimeJacobianWorkspace<T,2> jacobian;
    th::Pr76PhaseWorkspace<D> work;
    th::Pr76PhaseWorkspace<T> plain;
    for(const std::size_t n:std::array<std::size_t,5>{1,4,2,3,1}) {
        std::vector<std::size_t> order(n); std::iota(order.begin(),order.end(),std::size_t{0});
        auto model=kernel<T>(order);
        std::vector<T> w(n,T{1}/static_cast<T>(n));
        T sum=0; for(std::size_t i=0;i+1<n;++i) { sum+=w[i]; } w.back()=T{1}-sum;
        std::vector<T> inputs{T{1000000},T{450}};
        inputs.insert(inputs.end(),w.begin(),w.end()-1);
        const auto roots=model.roots_full(inputs[0],inputs[1],w,plain);
        require(roots.status==th::Pr76RootStatus::success,"runtime roots");
        const auto k=roots.count-1;
        const auto callback=[&](std::span<const D> in,std::span<D> out) {
            const auto result=model.evaluate_reduced(in[0],in[1],in.subspan(2),k,work);
            out[0]=result.z;
            for(std::size_t i=0;i<n;++i) { out[i+1]=result.ln_phi[i]; }
        };
        const auto result=jacobian.evaluate(callback,inputs,n+1);
        const auto ref=reference(inputs[0],inputs[1],w,order);
        require(result.input_count==n+1 && result.output_count==n+1,"runtime rectangular contract");
        near(result.values[0],ref[k].z);
        for(std::size_t c=0;c<n+1;++c) {
            near(result.jacobian[c],ref[k].dz[c]-(c>=2?ref[k].dz[n+1]:0));
            for(std::size_t i=0;i<n;++i) {
                near(result.jacobian[(i+1)*(n+1)+c],ref[k].dln[i][c]-(c>=2?ref[k].dln[i][n+1]:0));
            }
        }
    }
}
template <typename T>
void input_domains() {
    auto model=kernel<T>(); th::Pr76PhaseWorkspace<T> work;
    const std::vector<T> w{T{1}/T{4},T{1}/T{2},T{1}/T{4}};
    for(const T bad:std::array<T,4>{T{0},T{-1},std::numeric_limits<T>::infinity(),std::numeric_limits<T>::quiet_NaN()}) {
        expect_error<std::domain_error>([&]{(void)model.evaluate_full(bad,T{450},w,0,work);});
        expect_error<std::domain_error>([&]{(void)model.evaluate_full(T{1000000},bad,w,0,work);});
    }
    expect_error<std::out_of_range>([&]{(void)model.evaluate_full(T{1000000},T{450},w,3,work);});
    expect_error<std::invalid_argument>([&]{(void)model.evaluate_full(T{1000000},T{450},{},0,work);});
    for(const std::vector<T>& bad:std::array<std::vector<T>,3>{{{T{1},T{1},T{0}},
            {T{1},-std::numeric_limits<T>::epsilon(),T{0}}, {T{1},std::numeric_limits<T>::quiet_NaN(),T{0}}}}) {
        expect_error<std::domain_error>([&]{(void)model.evaluate_full(T{1000000},T{450},bad,0,work);});
    }
    Fixture f; f.input.applicability.pressure_pa=th::ClosedInterval{2e5,2e6};
    f.input.applicability.temperature_k=th::ClosedInterval{300,600};
    auto bounded=th::Pr76Phase<T>::from_parameters(f.select({0,1,2}));
    const auto low_pressure=bounded.roots_full(T{100000},T{450},w,work);
    const auto high_temperature=bounded.roots_full(T{1000000},T{650},w,work);
    require(low_pressure.status==th::Pr76RootStatus::success && low_pressure.count>0,
            "declared pressure extrapolation must remain computationally available");
    require(high_temperature.status==th::Pr76RootStatus::success && high_temperature.count>0,
            "declared temperature extrapolation must remain computationally available");
    const auto& applicability=bounded.parameters().applicability();
    require(applicability.assess(450,100000)==th::RangeAssessment::outside_declared_bounds,
            "pressure extrapolation lost advisory outside status");
    require(applicability.assess(650,1000000)==th::RangeAssessment::outside_declared_bounds,
            "temperature extrapolation lost advisory outside status");
    require(applicability.assess(300,200000)==th::RangeAssessment::inside_declared_bounds &&
                applicability.assess(600,2000000)==th::RangeAssessment::inside_declared_bounds,
            "inclusive declared T/p bounds changed");
    require(bounded.roots_full(T{200000},T{300},w,work).status==th::Pr76RootStatus::success,
            "inclusive bounds remain computationally available");
    using D=ad::Dual<T,1>; th::Pr76PhaseWorkspace<D> dw;
    std::vector<D> dx; for(const T v:w) {dx.emplace_back(v);}
    const D bad_seed{T{1000000},{std::numeric_limits<T>::infinity()}};
    expect_error<std::domain_error>([&]{(void)model.evaluate_full(bad_seed,D{T{450}},dx,0,dw);});
    dx[0]=D{w[0],{std::numeric_limits<T>::infinity()}};
    expect_error<std::domain_error>([&]{(void)model.evaluate_full(D{T{1000000}},D{T{450}},dx,0,dw);});
}
template <typename T>
void ownership_recovery() {
    auto model=kernel<T>(); auto copy=model;
    th::Pr76PhaseWorkspace<T> work,other;
    static_assert(!std::is_copy_constructible_v<th::Pr76PhaseWorkspace<T>>);
    static_assert(!std::is_move_constructible_v<th::Pr76PhaseWorkspace<T>>);
    const std::vector<T> w{T{1}/T{4},T{1}/T{2},T{1}/T{4}};
    const auto first=model.evaluate_full(T{1000000},T{450},w,2,work);
    const auto saved=first.ln_phi; const auto capacity=work.capacity();
    expect_error<th::Pr76PhaseError>([&]{(void)model.evaluate_full(T{1000000},T{450},w,2,work,{1});});
    const auto again=model.evaluate_full(T{1000000},T{450},w,2,work);
    require(saved==first.ln_phi && saved==again.ln_phi && work.capacity()==capacity,"result ownership and recovery");
    const auto distinct=copy.evaluate_full(T{200000},T{350},w,2,other);
    require(distinct.ln_phi!=first.ln_phi && first.ln_phi==saved,"independent workspaces");
    auto moved=std::move(copy);
    expect_error<std::invalid_argument>([&]{(void)copy.roots_full(T{1000000},T{450},w,other);});
    near(moved.evaluate_full(T{1000000},T{450},w,2,other).z,first.z);
}
template <typename T>
void attraction_cancellation() {
    Fixture f;
    f.input.pure[1].critical_temperature->value=specs[0].tc;
    f.input.pure[1].critical_pressure->value=specs[0].pc;
    f.input.pure[1].acentric_factor->value=specs[0].omega;
    th::Pr76PhaseWorkspace<T> work;
    const std::vector<T> w{T{1}/T{2},T{1}/T{2}};
    for(const double interaction:std::array<double,2>{2,4}) {
        f.input.binary[0].kij->value=interaction;
        auto model=th::Pr76Phase<T>::from_parameters(f.select({0,1}));
        const auto roots=model.roots_full(T{1000000},T{400},w,work);
        require(roots.status==th::Pr76RootStatus::success && roots.count==1,"zero/negative attraction continuation");
        const auto result=model.evaluate_full(T{1000000},T{400},w,0,work);
        const auto pure=pure_reference(400,specs[0]);
        const long double rt=8.31446261815324L*400;
        const long double B=pure.b*1e6L/rt, A=(1-interaction/2)*pure.a*1e6L/(rt*rt);
        const auto zr=reference_roots(A,B).back();
        const long double expected=(zr-1)-std::log(zr-B)-A/(2*std::sqrt(2.0L)*B)*
            std::log((zr+(1+std::sqrt(2.0L))*B)/(zr+(1-std::sqrt(2.0L))*B));
        near(result.ln_phi[0],expected); near(result.ln_phi[1],expected);
    }
}
void derivative_guard() {
    Fixture f;
    f.input.pure[1].critical_temperature->value=specs[0].tc;
    f.input.pure[1].critical_pressure->value=specs[0].pc;
    f.input.pure[1].acentric_factor->value=specs[0].omega;
    const double A=.457235528921382189383460196225183789+16*std::numeric_limits<double>::epsilon();
    const double B=.077796073903888455971844710037333184;
    f.input.binary[0].kij->value=2*(1-(A/B)*.07780/.45724);
    const double p=B*specs[0].pc/.07780;
    auto model=th::Pr76Phase<double>::from_parameters(f.select({0,1}));
    th::Pr76PhaseWorkspace<double> plain;
    const std::vector<double> w{.5,.5};
    const auto roots=model.roots_full(p,400,w,plain);
    require(roots.status==th::Pr76RootStatus::success && roots.count==1 && !roots.roots[0].derivative_valid,
            "ill-conditioned simple primal distinguished from derivative");
    require(std::isfinite(model.evaluate_full(p,400.0,w,0,plain).z),"primal available");
    using D=ad::Dual<double,1>; th::Pr76PhaseWorkspace<D> work;
    const std::vector<D> x{D{.5},D{.5}};
    bool caught=false;
    try { (void)model.evaluate_full(D::variable(p,0),D{400},x,0,work); }
    catch(const th::Pr76PhaseError& error) { caught=error.code()==th::Pr76PhaseErrorCode::ill_conditioned_derivative; }
    require(caught,"AD must reject unreliable root slopes");
}

template <typename T>
void nested_mixture_temperature_contract() {
    using Outer = ad::Dual<T, 2>;
    using Nested = ad::Dual<Outer, 1>;

    const auto parameters =
        Fixture{}.select({0, 1, 2});
    auto mixture =
        th::Pr76Mixture<T>::from_parameters(
            parameters);

    const Outer temperature =
        Outer::variable(T{450}, 0);
    const std::vector<Outer> outer_composition{
        Outer{T{0.25}, {T{0}, T{1}}},
        Outer{T{0.50}},
        Outer{T{0.25}, {T{0}, T{-1}}}};

    std::vector<Nested> nested_composition;
    nested_composition.reserve(
        outer_composition.size());
    for (const auto& fraction :
         outer_composition) {
        nested_composition.emplace_back(
            fraction);
    }

    th::Pr76MixtureWorkspace<Nested>
        nested_workspace;
    const auto nested =
        mixture.evaluate_full(
            Nested::variable(
                temperature,
                0U),
            std::span<const Nested>{
                nested_composition},
            nested_workspace);

    const Outer da_dtemperature =
        nested.a.derivative(0U);

    using First = ad::Dual<T, 1>;
    const auto first_temperature_derivative =
        [&](T temperature_value,
            T x0,
            T x2) {
            const std::vector<First> composition{
                First{x0},
                First{T{0.50}},
                First{x2}};
            th::Pr76MixtureWorkspace<First>
                workspace;
            const auto value =
                mixture.evaluate_full(
                    First::variable(
                        temperature_value,
                        0U),
                    std::span<const First>{
                        composition},
                    workspace);
            return value.a.derivative(0U);
        };

    const T direct =
        first_temperature_derivative(
            T{450},
            T{0.25},
            T{0.25});
    near(
        da_dtemperature.value(),
        static_cast<long double>(direct));

    const T temperature_step =
        T{1} / T{100};
    const T d2a_dtemperature2 =
        (first_temperature_derivative(
             T{450} + temperature_step,
             T{0.25},
             T{0.25}) -
         first_temperature_derivative(
             T{450} - temperature_step,
             T{0.25},
             T{0.25})) /
        (T{2} * temperature_step);

    const T composition_step =
        T{1} / T{100000};
    const T mixed_tangent =
        (first_temperature_derivative(
             T{450},
             T{0.25} + composition_step,
             T{0.25} - composition_step) -
         first_temperature_derivative(
             T{450},
             T{0.25} - composition_step,
             T{0.25} + composition_step)) /
        (T{2} * composition_step);

    const auto finite_difference_close =
        [](T actual, T expected) {
            const T scale =
                std::max(
                    {T{1},
                     std::abs(actual),
                     std::abs(expected)});
            return std::isfinite(actual) &&
                std::isfinite(expected) &&
                std::abs(actual - expected) <=
                    (T{2} / T{100000}) * scale;
        };

    require(
        finite_difference_close(
            da_dtemperature.derivative(0U),
            d2a_dtemperature2),
        "nested PR76 mixture d2a/dT2 disagrees with fresh central perturbation");
    require(
        finite_difference_close(
            da_dtemperature.derivative(1U),
            mixed_tangent),
        "nested PR76 mixture d(da/dT)/dx tangent disagrees with fresh central perturbation");

    using InvalidNested =
        ad::Dual<ad::Dual<T, 1>, 1>;
    static_assert(
        !th::detail::Pr76PhaseNumber<
            InvalidNested,
            T>);
    static_assert(
        th::detail::Pr76Number<
            InvalidNested,
            T>);
}

template <typename T>
void selected_phase_fugacity_contract() {
    auto model = kernel<T>();
    th::Pr76PhaseWorkspace<T> direct_workspace;
    th::Pr76PhaseWorkspace<T> facade_workspace;
    const std::vector<T> composition{T{0.25}, T{0.50}, T{0.25}};
    const T pressure = T{1.0e6};
    const T temperature = T{450.0};
    const auto roots = model.roots_full(
        pressure, temperature, composition, direct_workspace);
    require(
        roots.status == th::Pr76RootStatus::success && roots.count > 0U,
        "PR76 selected-phase fugacity fixture has no resolved root");
    const std::size_t root_index = roots.count - 1U;
    const auto direct = model.evaluate_full(
        pressure, temperature, composition, root_index, direct_workspace);
    const auto wrapped = th::evaluate_selected_phase_fugacity(
        model, pressure, temperature,
        std::span<const T>{composition},
        th::Pr76SelectedPhase{root_index, {}},
        facade_workspace);
    require(
        wrapped.ln_phi.size() == direct.ln_phi.size(),
        "PR76 selected-phase fugacity size changed");
    for (std::size_t i = 0U; i < direct.ln_phi.size(); ++i) {
        near(wrapped.ln_phi[i], static_cast<long double>(direct.ln_phi[i]));
    }
    static_assert(
        th::SelectedPhaseFugacityCapabilities<th::Pr76Phase<T>>::derivative_support ==
        th::SelectedPhaseFugacityDerivativeSupport::scalar_generic_first_order);
}

template <typename T>
void selected_phase_density_contract() {
    auto model = kernel<T>();
    th::Pr76PhaseWorkspace<T> plain_workspace;
    const std::vector<T> composition{
        T{0.25}, T{0.50}, T{0.25}};
    const T pressure = T{1.0e6};
    const T temperature = T{450.0};
    const auto roots = model.roots_full(
        pressure, temperature, composition,
        plain_workspace);
    require(
        roots.status == th::Pr76RootStatus::success &&
            roots.count > 0U,
        "PR76 selected-phase density fixture has no root");
    const std::size_t root_index =
        roots.count - 1U;
    const auto direct = model.evaluate_full(
        pressure,
        temperature,
        composition,
        root_index,
        plain_workspace);
    th::Pr76PhaseWorkspace<T> density_workspace;
    const auto wrapped =
        th::evaluate_selected_phase_molar_density(
            model,
            pressure,
            temperature,
            std::span<const T>{composition},
            th::Pr76SelectedPhase{
                root_index, {}},
            density_workspace);
    const long double expected =
        static_cast<long double>(pressure) /
        (static_cast<long double>(direct.z) *
         8.31446261815324L *
         static_cast<long double>(temperature));
    near(
        wrapped.molar_density_mol_per_m3,
        expected);

    using D = ad::Dual<T, 2>;
    std::vector<D> x;
    x.reserve(composition.size());
    for (const T value : composition) {
        x.emplace_back(value);
    }
    th::Pr76PhaseWorkspace<D> derivative_workspace;
    const auto derivative =
        th::evaluate_selected_phase_molar_density(
            model,
            D::variable(pressure, 0U),
            D::variable(temperature, 1U),
            std::span<const D>{x},
            th::Pr76SelectedPhase{
                root_index, {}},
            derivative_workspace);
    th::Pr76PhaseWorkspace<D> phase_workspace;
    const auto phase =
        model.evaluate_full(
            D::variable(pressure, 0U),
            D::variable(temperature, 1U),
            std::span<const D>{x},
            root_index,
            phase_workspace);
    const T rho =
        derivative
            .molar_density_mol_per_m3
            .value();
    for (std::size_t lane = 0U;
         lane < 2U;
         ++lane) {
        const T dp =
            lane == 0U ? T{1} : T{0};
        const T dt =
            lane == 1U ? T{1} : T{0};
        const T expected_derivative =
            rho *
            (dp / pressure -
             phase.z.derivative(lane) /
                 phase.z.value() -
             dt / temperature);
        near(
            derivative
                .molar_density_mol_per_m3
                .derivative(lane),
            static_cast<long double>(
                expected_derivative));
    }
}

void headers() {
    auto model=kernel<double>();
    require(std::isfinite(pt_plain_header(model)),"plain PT header translation unit");
    require(roots_plain_header(),"standalone root header translation unit");
}
template <typename T>
void run_typed(std::string_view name) {
    if(name=="roots_factored") {roots_factored<T>();}
    else if(name=="root_topology") {root_topology<T>();}
    else if(name=="root_scale") {root_scale<T>();}
    else if(name=="phase_values") {phase_values<T>();}
    else if(name=="full_derivatives") {full_derivatives<T>();}
    else if(name=="reduced_derivatives") {reduced_derivatives<T>();}
    else if(name=="gibbs_duhem") {gibbs_duhem<T>();}
    else if(name=="pure_limit") {pure_limit<T>();}
    else if(name=="low_pressure") {low_pressure<T>();}
    else if(name=="log_ratio") {log_ratio<T>();}
    else if(name=="zero_trace") {zero_trace<T>();}
    else if(name=="permutations") {permutations<T>();}
    else if(name=="runtime_shapes") {runtime_shapes<T>();}
    else if(name=="input_domains") {input_domains<T>();}
    else if(name=="ownership_recovery") {ownership_recovery<T>();}
    else if(name=="attraction_cancellation") {attraction_cancellation<T>();}
    else if(name=="nested_mixture_temperature") {nested_mixture_temperature_contract<T>();}
    else if(name=="selected_phase_fugacity") {selected_phase_fugacity_contract<T>();}
    else if(name=="selected_phase_density") {selected_phase_density_contract<T>();}
    else {throw std::invalid_argument("unknown case");}
}
} // namespace
int main(int argc,char** argv) {
    try {
        if(argc!=2) {throw std::invalid_argument("one test case required");}
        const std::string_view name=argv[1];
        if(name=="headers") {headers();}
        else if(name=="derivative_guard") {derivative_guard();}
        else {run_typed<float>(name);run_typed<double>(name);run_typed<long double>(name);}
        std::cout<<"[PASS] "<<name<<'\n';
        return 0;
    } catch(const std::exception& error) {
        std::cerr<<"[FAIL] "<<error.what()<<'\n'; return 1;
    }
}
