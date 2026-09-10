#include "test_support.hpp"
#include <iostream>
#include <limits>
#include <string_view>
#include <utility>

bool sw92_headers();
namespace st=sw92_test; namespace th=mpmc::thermodynamics;

namespace {
void contract(){
 auto p=st::binary_input(st::methane); auto ps=th::Sw92ParameterSet::create(p.catalog,p.order,p.input);
 st::require(ps.components().size()==2&&ps.water_index()==1,"parameter order/water");
 st::require(ps.species(0)==th::Sw92Species::hydrocarbon,"species routing");
 auto bad=p; bad.input.model_id="SW"; st::expect_error<th::ContractError>([&]{(void)th::Sw92ParameterSet::create(bad.catalog,bad.order,bad.input);});
 bad=p; bad.input.water_binary.clear(); st::expect_error<th::ContractError>([&]{(void)th::Sw92ParameterSet::create(bad.catalog,bad.order,bad.input);});
 bad=st::binary_input(st::h2s); bad.input.water_binary[0].nonaqueous_rule=th::Sw92NonAqueousWaterRule::constant;
 bad.input.water_binary[0].nonaqueous_kij=st::scalar(.1,th::Unit::dimensionless,st::paper("invalid H2S override"));
 st::expect_error<th::ContractError>([&]{(void)th::Sw92ParameterSet::create(bad.catalog,bad.order,bad.input);});
}
void water_alpha(){
 st::near(th::sw92_water_alpha(377.15,647.3,0.0),st::alpha_golden[0]);
 st::near(th::sw92_water_alpha(377.15,647.3,1.0),st::alpha_golden[1]);
 st::near(th::sw92_water_alpha(377.15,647.3,5.0),st::alpha_golden[2]);
 st::expect_error<std::domain_error>([]{(void)th::sw92_water_alpha(377.15,647.3,-1.0);});
}
void bip_correlations(){
 const st::Spec* gases[]={&st::methane,&st::nitrogen,&st::co2,&st::h2s}; const double t[]={377.15,376.15,423.15,377.15};
 for(std::size_t i=0;i<4;++i){auto ps=st::binary_parameters(*gases[i]);auto m=th::Sw92Mixture<double>::from_parameters(ps);
  st::near(m.kij(t[i],1.0,th::SwPhaseFamily::aqueous,0,1),st::kij_golden[i][0]);
  st::near(m.kij(t[i],1.0,th::SwPhaseFamily::nonaqueous,0,1),st::kij_golden[i][1]);
  st::near(m.kij(t[i],1.0,th::SwPhaseFamily::aqueous,1,0),st::kij_golden[i][0]);}
}
void mixing_reference(){auto ps=st::binary_parameters(st::methane);auto m=th::Sw92Mixture<double>::from_parameters(ps);th::Sw92MixtureWorkspace<double>w;st::Vec x{.003,.997};
 auto aq=m.evaluate(377.15,x,1.0,th::SwPhaseFamily::aqueous,w);auto na=m.evaluate(377.15,x,1.0,th::SwPhaseFamily::nonaqueous,w);
 st::near(aq.a,st::methane_mix_golden[0][0]);st::near(aq.b,st::methane_mix_golden[0][1]);st::near(na.a,st::methane_mix_golden[1][0]);st::near(na.b,st::methane_mix_golden[1][1]);}
void phase_three_roots(){auto ps=st::binary_parameters(st::co2);auto p=th::Sw92Phase<double>::from_parameters(ps);th::Sw92PhaseWorkspace<double>w;st::Vec x{.7,.3};auto roots=p.roots(3e6,340.0,x,0.0,th::SwPhaseFamily::nonaqueous,w);st::require(roots.status==th::Sw92RootStatus::success&&roots.count==3,"three-root topology");
 for(std::size_t r=0;r<3;++r){auto v=p.evaluate(3e6,340.0,x,0.0,th::SwPhaseFamily::nonaqueous,r,w);st::near(v.z,st::co2_na_phase_golden[r][0],5e-12L);st::near(v.ln_phi.at(0),st::co2_na_phase_golden[r][1],7e-12L);st::near(v.ln_phi.at(1),st::co2_na_phase_golden[r][2],7e-12L);st::require(v.family==th::SwPhaseFamily::nonaqueous&&v.root_index==r,"phase/root diagnostic");}}
void phase_aqueous_reference(){auto ps=st::binary_parameters(st::co2);auto p=th::Sw92Phase<double>::from_parameters(ps);th::Sw92PhaseWorkspace<double>w;st::Vec x{.02,.98};auto roots=p.roots(20e6,423.15,x,1.0,th::SwPhaseFamily::aqueous,w);st::require(roots.status==th::Sw92RootStatus::success&&roots.count==1,"AQ root topology");auto v=p.evaluate(20e6,423.15,x,1.0,th::SwPhaseFamily::aqueous,0,w);st::near(v.z,st::co2_aq_phase_golden[0],7e-12L);st::near(v.ln_phi.at(0),st::co2_aq_phase_golden[1],1e-11L);st::near(v.ln_phi.at(1),st::co2_aq_phase_golden[2],1e-11L);}
void permutations(){auto a=th::Sw92Phase<double>::from_parameters(st::binary_parameters(st::co2,false));auto b=th::Sw92Phase<double>::from_parameters(st::binary_parameters(st::co2,true));th::Sw92PhaseWorkspace<double>wa,wb;auto va=a.evaluate(3e6,340.0,st::Vec{.7,.3},0.0,th::SwPhaseFamily::nonaqueous,2,wa);auto vb=b.evaluate(3e6,340.0,st::Vec{.3,.7},0.0,th::SwPhaseFamily::nonaqueous,2,wb);st::near(va.z,vb.z,1e-13L);st::near(va.ln_phi[0],vb.ln_phi[1],1e-13L);st::near(va.ln_phi[1],vb.ln_phi[0],1e-13L);}
void background_pairs(){st::Prepared p;p.input.model_id=std::string(th::sw92_corrected_profile);p.input.dataset_id="SW92-background-test";p.input.revision="v1";p.input.applicability={{std::nullopt,std::nullopt,st::paper("model basis; synthetic gas/gas BIP")},std::nullopt};st::add_pure(p,st::methane);st::add_pure(p,st::nitrogen);st::add_pure(p,st::water);p.input.water_binary=st::binary_input(st::methane).input.water_binary;p.input.water_binary.push_back(st::binary_input(st::nitrogen).input.water_binary.front());auto src=st::synthetic("methane/nitrogen family BIPs");p.input.nonwater_binary.push_back({st::methane.id,st::nitrogen.id,st::scalar(.125,th::Unit::dimensionless,src),st::scalar(-.0625,th::Unit::dimensionless,src)});p.order={st::methane.id,st::nitrogen.id,st::water.id};auto ps=th::Sw92ParameterSet::create(p.catalog,p.order,p.input,th::DataPolicy::allow_synthetic_tests);auto m=th::Sw92Mixture<double>::from_parameters(ps);st::near(m.kij(350,1,th::SwPhaseFamily::aqueous,0,1),.125L);st::near(m.kij(350,1,th::SwPhaseFamily::nonaqueous,0,1),-.0625L);auto missing=p;missing.input.nonwater_binary.clear();st::expect_error<th::ContractError>([&]{(void)th::Sw92ParameterSet::create(missing.catalog,missing.order,missing.input,th::DataPolicy::allow_synthetic_tests);});}
void input_domains(){auto ps=st::binary_parameters(st::methane);auto m=th::Sw92Mixture<double>::from_parameters(ps);th::Sw92MixtureWorkspace<double>w;st::expect_error<std::domain_error>([&]{(void)m.evaluate(377.15,st::Vec{.4,.5},1,th::SwPhaseFamily::aqueous,w);});st::expect_error<std::domain_error>([&]{(void)m.evaluate(377.15,st::Vec{.5,.5},-1,th::SwPhaseFamily::aqueous,w);});st::expect_error<std::invalid_argument>([&]{(void)m.kij(377.15,1,static_cast<th::SwPhaseFamily>(99),0,1);});st::expect_error<std::out_of_range>([&]{(void)m.kij(377.15,1,th::SwPhaseFamily::aqueous,0,2);});auto p=th::Sw92Phase<double>::from_parameters(ps);th::Sw92PhaseWorkspace<double>pw;st::expect_error<std::domain_error>([&]{(void)p.roots(0,377.15,st::Vec{.5,.5},1,th::SwPhaseFamily::aqueous,pw);});st::expect_error<std::out_of_range>([&]{(void)p.evaluate(1e6,377.15,st::Vec{.5,.5},1,th::SwPhaseFamily::aqueous,3,pw);});}
void headers(){st::require(sw92_headers(),"self-contained headers");}
using Test=std::pair<std::string_view,void(*)()>;constexpr Test tests[]={{"contract",contract},{"water_alpha",water_alpha},{"bip_correlations",bip_correlations},{"mixing_reference",mixing_reference},{"phase_three_roots",phase_three_roots},{"phase_aqueous_reference",phase_aqueous_reference},{"permutations",permutations},{"background_pairs",background_pairs},{"input_domains",input_domains},{"headers",headers}};
}
int main(int argc,char**argv){try{if(argc!=2)throw std::invalid_argument("one test name required");for(auto [n,f]:tests)if(n==argv[1]){f();std::cout<<"[PASS] "<<n<<'\n';return 0;}throw std::invalid_argument("unknown test");}catch(const std::exception&e){std::cerr<<"[FAIL] "<<e.what()<<'\n';return 1;}}
