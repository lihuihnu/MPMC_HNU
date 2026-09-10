#ifndef MPMC_TEST_SW92_SUPPORT_HPP
#define MPMC_TEST_SW92_SUPPORT_HPP

#include <mpmc/thermodynamics/sw92_phase.hpp>

#include <cmath>
#include <cstddef>
#include <optional>
#include <source_location>
#include <stdexcept>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

namespace sw92_test {
namespace th = mpmc::thermodynamics;
using Vec = std::vector<double>;

inline void require(bool value, std::string_view message,
                    std::source_location where = std::source_location::current()) {
    if (!value) throw std::runtime_error(std::string(where.file_name()) + ":" +
        std::to_string(where.line()) + ": " + std::string(message));
}
inline void near(double actual, long double expected, long double tol=3e-12L) {
    require(std::isfinite(actual), "nonfinite result");
    require(std::abs(static_cast<long double>(actual)-expected) <=
            tol*(1.0L+std::abs(expected)), "reference mismatch");
}
template<class Error,class F> inline void expect_error(F&& f) {
    bool caught=false; try { f(); } catch(const Error&) { caught=true; }
    require(caught,"expected exception missing");
}

inline th::Provenance paper(std::string locator) {
    return {th::SourceKind::literature,
        "doi:10.1016/0378-3812(92)85105-H",
        "Fluid Phase Equilibria 77 (1992) 217-240; authors' errata in supplied PDF",
        std::move(locator),
        "SW92 corrected-original regression fixture; not experimental validation",
        "User PDF sha256:cb5b1d5034d78d934e887449ce0c89692d43431835d1371da5606d95b2c6bf58",
        "Bibliographic/formula facts only; source PDF is not redistributed"};
}
inline th::Provenance synthetic(std::string locator) {
    return {th::SourceKind::synthetic_test,"MPMC_HNU SW92 structural fixture","v1",
        std::move(locator),"Artificial gas/gas BIP only for routing tests",
        "tests/thermodynamics/sw92/test_support.hpp","Repository test data"};
}
inline th::SourcedScalar scalar(double value, th::Unit unit, const th::Provenance& source,
                                std::string original="SI or dimensionless",
                                std::string conversion="identity") {
    return {value,unit,source,std::move(original),std::move(conversion)};
}

struct Spec { const char* id; const char* display; th::Sw92Species species;
              double tc,pc_bar,omega,na; bool eq17; };
inline constexpr Spec methane{"methane","Methane",th::Sw92Species::hydrocarbon,190.6,46.0,.0108,.4850,false};
inline constexpr Spec nitrogen{"nitrogen","Nitrogen",th::Sw92Species::nitrogen,126.1,34.0,.0403,.4778,false};
inline constexpr Spec co2{"carbon-dioxide","Carbon dioxide",th::Sw92Species::carbon_dioxide,304.2,73.8,.2273,.1896,false};
inline constexpr Spec h2s{"hydrogen-sulfide","Hydrogen sulfide",th::Sw92Species::hydrogen_sulfide,373.2,89.4,.1081,0,true};
inline constexpr Spec water{"water","Water",th::Sw92Species::water,647.3,221.2,.3434,0,false};

struct Prepared { std::vector<th::Component> catalog; th::Sw92ParameterInput input;
                  std::vector<std::string> order; };
inline void add_pure(Prepared& p,const Spec& s) {
    p.catalog.push_back({s.id,s.display,th::ComponentKind::pure,paper("Table 3 component identity"),std::nullopt});
    auto src=paper("Table 3 physical properties");
    p.input.pure.push_back({s.id,s.species,
        scalar(s.tc,th::Unit::kelvin,src,"K","identity"),
        scalar(s.pc_bar*100000.0,th::Unit::pascal,src,"bar","bar * 100000 -> Pa"),
        scalar(s.omega,th::Unit::dimensionless,src)});
}
inline Prepared binary_input(const Spec& gas,bool reverse=false) {
    Prepared p; p.input.model_id=std::string(th::sw92_corrected_profile);
    p.input.dataset_id="SW92-Table3-Table5-corrected";
    p.input.revision="user-PDF-sha256-cb5b1d50";
    p.input.applicability={{std::nullopt,std::nullopt,paper("Eqs. (9)-(17), Tables 2-5; no global bounds inferred")},std::nullopt};
    add_pure(p,gas); add_pure(p,water);
    th::Sw92WaterBinaryRecord wb; wb.component_id=gas.id;
    if(gas.eq17) wb.nonaqueous_rule=th::Sw92NonAqueousWaterRule::hydrogen_sulfide_eq17;
    else { wb.nonaqueous_rule=th::Sw92NonAqueousWaterRule::constant;
           wb.nonaqueous_kij=scalar(gas.na,th::Unit::dimensionless,paper("Table 5 non-aqueous water BIP")); }
    p.input.water_binary.push_back(std::move(wb));
    p.order=reverse?std::vector<std::string>{water.id,gas.id}:std::vector<std::string>{gas.id,water.id};
    return p;
}
inline th::Sw92ParameterSet binary_parameters(const Spec& s,bool reverse=false) {
    auto p=binary_input(s,reverse); return th::Sw92ParameterSet::create(p.catalog,p.order,p.input);
}

inline constexpr long double alpha_golden[]={
1.4468437971786006357582097503026120589044361383121L,
1.4533912933678893669426945675202156346403030019293L,
1.4855093864815170593420136191463497302792515365434L};
inline constexpr long double kij_golden[4][2]={
{-0.050087331231713200135077384503025312225024009212560L,.4850L},
{-0.31585475849345757335448057097541633624107850911975L,.4778L},
{0.026447245451780916487237455263752592949016222261592L,.1896L},
{0.032329439978563772775991425509110396570203644158628L,0.13002865621650589496248660235798499464094319399786L}};
inline constexpr long double methane_mix_golden[2][2]={
{0.86745438746364740129193172097955612789815088749160L,0.000018952907751928144677141936001257960531488324553818L},
{0.86619990096009627570063138316944955797585931781058L,0.000018952907751928144677141936001257960531488324553818L}};
inline constexpr long double co2_na_phase_golden[3][3]={
{0.057180176767981325431846968271992646995581273293750L,1.0208909593399308307491967248213269787130187979653L,-1.3825890900379509448216125538859566769979889994945L},
{0.078212323419495588475436214780373737724354591089947L,0.84558993181719901783378266277409243108980981854591L,-0.96019712374354804789208960161513467098981239359815L},
{0.83877383250027555663944121834758847218423364534672L,-0.095622586793950066191409126737705562686980512832323L,-0.28634773253331448691253791693978729483642277716321L}};
inline constexpr long double co2_aq_phase_golden[3]={
0.13538851508960205544139583262231128087805764877944L,
3.6869408752143027087094016483296241185811218859165L,
-3.6649702467280267477565463531227039969075629099210L};

} // namespace sw92_test
#endif
