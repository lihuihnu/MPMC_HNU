#include <mpmc/flow_discretization_petsc/sw92_transactional_phase_transition_restart.hpp>

#include "test_support.hpp"

#include <petscsys.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <optional>
#include <span>
#include <stdexcept>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

bool sw92_transactional_phase_transition_restart_header();

namespace {
namespace th = mpmc::thermodynamics;
namespace flow = mpmc::flow;
namespace mesh = mpmc::mesh;
namespace dp = mpmc::discretization_petsc;
namespace fdp = mpmc::flow_discretization_petsc;
namespace st = sw92_test;

void require_tx(bool condition, std::string_view message) {
    if (!condition) throw std::runtime_error(std::string{message});
}

th::Provenance nist_source(std::string locator) {
    return {th::SourceKind::database,
            "NIST Chemistry WebBook SRD 69",
            "accessed-2026-09-25",
            std::move(locator),"",
            "Public NIST Chemistry WebBook record",
            "Use as source data; no NIST code redistributed"};
}

th::Sw92Phase<double> binary_model() {
    auto prepared=st::binary_input(st::co2);
    for(auto& component:prepared.catalog) {
        const double value=component.id==st::co2.id?0.0440095:0.0180153;
        component.molar_mass=st::scalar(
            value,th::Unit::kilogram_per_mole,
            nist_source(component.id+" molecular weight"),
            "kg/mol","identity");
    }
    auto parameters=th::Sw92ParameterSet::create(
        prepared.catalog,prepared.order,prepared.input);
    return th::Sw92Phase<double>::from_parameters(parameters);
}

mesh::PartitionSnapshot one_cell_partition() {
    mesh::Topology topology{{{},{},{},{mesh::GlobalEntityId{UINT64_C(10)}}},{}};
    mesh::EntityOwnerRanks owners;
    owners.cells={mesh::PartitionRank{0U}};
    return mesh::PartitionSnapshot::create(
        topology,mesh::PartitionRank{0U},1U,std::move(owners));
}
dp::ParallelOwnedConnectionSchedule3D empty_schedule() {
    return {mesh::PartitionRank{0U},1U,{},{}};
}
dp::PetscMpiAijSymbolicPreallocation3D one_cell_bridge() {
    return {mesh::PartitionRank{0U},1U,0,1,1,
            {mesh::LocalIndex{0U}},
            {mesh::GlobalEntityId{UINT64_C(10)}},
            {0},{1},{0},{0}};
}
dp::OwnedCellStructuralColumnPatternSnapshot3D one_cell_pattern() {
    return {mesh::PartitionRank{0U},1U,1U,0,1,1,
            {{mesh::LocalIndex{0U},mesh::GlobalEntityId{UINT64_C(10)},
              0,0U,1U,0U,0U}},
            {0},{}};
}

PetscErrorCode zero_rock(
    const flow::NaturalVariableLayoutDescriptor& layout,
    std::span<const double>,void*,
    std::optional<fdp::Sw92RockThermalStorageLinearization3D>* output,
    fdp::NaturalVariableSnesEvaluationStatus3D* status) {
    if(!output||!status) return PETSC_ERR_ARG_NULL;
    output->emplace(fdp::Sw92RockThermalStorageLinearization3D{
        0.0,std::vector<double>(layout.unknown_count(),0.0),
        {"zero rock","sw92-transactional-test","v1"}});
    *status=fdp::NaturalVariableSnesEvaluationStatus3D::success;
    return PETSC_SUCCESS;
}
PetscErrorCode linear_kr(
    const flow::NaturalVariableCellState2P& state,void*,
    std::optional<fdp::Sw92TwoPhaseRelativePermeabilityLinearization3D>* output,
    fdp::NaturalVariableSnesEvaluationStatus3D* status) {
    if(!output||!status) return PETSC_ERR_ARG_NULL;
    const auto q=state.layout().unknown_count();
    std::array<std::vector<double>,2> g{
        std::vector<double>(q,0.0),std::vector<double>(q,0.0)};
    const auto s=state.layout().independent_saturation_unknown_index();
    g[0][s]=1.0; g[1][s]=-1.0;
    output->emplace(fdp::Sw92TwoPhaseRelativePermeabilityLinearization3D{
        {state.phase_saturation(0U),state.phase_saturation(1U)},std::move(g)});
    *status=fdp::NaturalVariableSnesEvaluationStatus3D::success;
    return PETSC_SUCCESS;
}
PetscErrorCode no_three_phase(
    const flow::NaturalVariableCellState3P&,void*,
    std::optional<flow::ThreePhaseSaturationConstitutiveNaturalVariableLinearization3P>*,
    fdp::NaturalVariableSnesEvaluationStatus3D*) {
    return PETSC_ERR_SUP;
}

th::Sw92SelectedPhase<double> na_selection(
    const th::Sw92Phase<double>& model,double p,double t,
    std::span<const double> x) {
    th::Sw92PhaseWorkspace<double> workspace;
    const auto roots=model.roots(
        p,t,x,0.0,th::SwPhaseFamily::nonaqueous,workspace);
    require_tx(roots.status==th::Sw92RootStatus::success && roots.count>0U,
               "SW92 NA source root unavailable");
    return {0.0,th::SwPhaseFamily::nonaqueous,roots.count-1U,{}};
}

struct IdentityResolverContext {
    flow::FrozenPhysicalPhaseIdentity aqueous{
        "test/sw92-transaction","aqueous"};
    flow::FrozenPhysicalPhaseIdentity nonaqueous{
        "test/sw92-transaction","nonaqueous"};
    flow::FrozenPhysicalPhaseIdentity extra{
        "test/sw92-transaction","extra"};
};

PetscErrorCode resolve_identity(
    mesh::GlobalEntityId,
    const flow::PhaseSetTransitionCandidate& candidate,
    const fdp::Sw92AuthoritativeTargetMaterialization3D& target,
    const flow::FrozenActivePhaseIdentityMap& source,
    void* raw,
    std::optional<flow::FrozenActivePhaseIdentityMap>* output) {
    if(!raw||!output) return PETSC_ERR_ARG_NULL;
    output->reset();
    auto* ctx=static_cast<IdentityResolverContext*>(raw);
    std::vector<flow::FrozenPhysicalPhaseIdentity> ids;
    ids.reserve(target.phases.size());
    if(candidate.source_phase_count==1U && candidate.target_phase_count==2U) {
        for(const auto& phase:target.phases) {
            ids.push_back(
                phase.metadata.thermodynamic_family==th::SwPhaseFamily::aqueous
                ? ctx->aqueous : source.identity(0U));
        }
    } else if(candidate.source_phase_count==2U &&
              candidate.target_phase_count==1U) {
        ids.push_back(ctx->nonaqueous);
        if(!source.find(ctx->nonaqueous).has_value()) return PETSC_ERR_ARG_INCOMP;
    } else {
        return PETSC_ERR_SUP;
    }
    try {
        output->emplace(std::move(ids));
        (void)flow::make_phase_identity_continuation_snapshot(
            candidate,source,**output);
    } catch(...) {
        output->reset();
        return PETSC_ERR_ARG_INCOMP;
    }
    return PETSC_SUCCESS;
}

template<class Closure>
std::pair<
    flow::PoreVolumeComponentAccumulationSnapshot3P,
    flow::PoreVolumeEnergyAccumulationSnapshot3P>
history_1p(
    const Closure& closure,
    const flow::NaturalVariableLayout1P& layout,
    std::span<const std::string> ids,
    std::span<const double> q,
    double porosity) {
    fdp::Sw92SinglePhaseProductionCellEvaluatorContext3D<Closure> ctx{
        &closure,1.0,{&zero_rock,nullptr},{}};
    std::optional<fdp::SinglePhaseCurrentCellLinearization3D> current;
    fdp::NaturalVariableSnesEvaluationStatus3D status=
        fdp::NaturalVariableSnesEvaluationStatus3D::success;
    require_tx(
        fdp::evaluate_sw92_single_phase_production_cell_3d<Closure>(
            mesh::LocalIndex{0U},mesh::GlobalEntityId{UINT64_C(10)},q,layout,ids,
            &ctx,&current,&status)==PETSC_SUCCESS &&
        status==fdp::NaturalVariableSnesEvaluationStatus3D::success &&
        current.has_value(),"failed to build SW92 1P baseline history");
    return {
        flow::build_single_phase_component_accumulation(current->state,porosity),
        flow::build_single_phase_energy_accumulation_snapshot(
            current->state,porosity,current->transport,current->caloric,current->rock)};
}

template<class Closure>
std::pair<
    flow::PoreVolumeComponentAccumulationSnapshot3P,
    flow::PoreVolumeEnergyAccumulationSnapshot3P>
history_2p(
    const Closure& closure,
    const flow::NaturalVariableLayout2P& layout,
    std::span<const std::string> ids,
    std::span<const double> q,
    double porosity) {
    fdp::Sw92TwoPhaseProductionCellEvaluatorContext3D<Closure> ctx{
        &closure,{&linear_kr,nullptr},{&zero_rock,nullptr},{}};
    std::optional<fdp::TwoPhaseCurrentCellLinearization3D> current;
    fdp::NaturalVariableSnesEvaluationStatus3D status=
        fdp::NaturalVariableSnesEvaluationStatus3D::success;
    require_tx(
        fdp::evaluate_sw92_two_phase_production_cell_3d<Closure>(
            mesh::LocalIndex{0U},mesh::GlobalEntityId{UINT64_C(10)},q,layout,ids,
            &ctx,&current,&status)==PETSC_SUCCESS &&
        status==fdp::NaturalVariableSnesEvaluationStatus3D::success &&
        current.has_value(),"failed to build SW92 2P baseline history");
    return {
        flow::build_two_phase_component_accumulation(current->state,porosity),
        flow::build_two_phase_energy_accumulation_snapshot(
            current->state,porosity,current->transport,current->caloric,current->rock)};
}

template<typename Provider>
void verify_transaction(
    fdp::Sw92TransactionalPhaseTransitionRestartContext3D<Provider>& rebuild,
    fdp::PostSnesSw92ProfileCPhaseTransitionScannerContext3D& scanner,
    fdp::Sw92TransactionalInitialCell3D initial,
    std::size_t expected_final_count) {
    std::unique_ptr<fdp::PhaseTransitionRebuiltNaturalVariableSystem3D> system;
    require_tx(
        fdp::materialize_sw92_transactional_initial_system_3d(
            PETSC_COMM_SELF,&rebuild,std::move(initial),&system)==PETSC_SUCCESS &&
        system!=nullptr,"failed to materialize SW92 transactional initial system");

    std::unique_ptr<fdp::PhaseTransitionRebuiltNaturalVariableSystem3D> final_system;
    Vec final_state=nullptr;
    std::optional<fdp::PostSnesPhaseTransitionControllerReport3D> report;
    fdp::PostSnesPhaseTransitionControllerOptions3D options;
    options.max_transition_restarts=3U;
    const PetscErrorCode error=
        fdp::solve_sw92_transactional_same_dt_3d(
            PETSC_COMM_SELF,std::move(system),&scanner,&rebuild,options,
            &final_system,&final_state,&report);
    require_tx(
        error==PETSC_SUCCESS && final_system && final_state &&
        report.has_value() &&
        report->outcome==fdp::PostSnesPhaseTransitionOutcome3D::stable_phase_set &&
        report->transition_restarts==1U &&
        final_system->numbering().cell(mesh::LocalIndex{0U}).phase_count==
            expected_final_count &&
        final_system->time_step_seconds()==1.0,
        "SW92 same-dt transactional controller did not reach the expected stable cardinality");

    bool history_matches=false;
    require_tx(
        final_system->accepted_history_matches_state(final_state,&history_matches)==
            PETSC_SUCCESS && history_matches,
        "SW92 same-dt restart changed the accepted BE baseline");
    (void)VecDestroy(&final_state);
}

void sourced_one_to_two() {
    auto model=binary_model();
    using Provider=flow::Sw92Co2WaterPropertyProvider;
    const double p=3.0e6,t=340.0;
    const std::array<double,2> z{0.7,0.3};
    const auto selection=na_selection(model,p,t,z);
    auto closure=flow::make_sw92_co2_water_selected_phase_property_closure(
        model,std::vector<th::Sw92SelectedPhase<double>>{selection});

    flow::NaturalVariableLayout1P layout{
        flow::NaturalVariableCompositionPivot1P::from_dependent_component(2U,0U)};
    std::vector<double> q(layout.unknown_count(),0.0);
    q[layout.pressure_unknown_index()]=p;
    q[layout.temperature_unknown_index()]=t;
    q[*layout.independent_composition_unknown_index(1U)]=z[1];
    const std::vector<std::string> ids{st::co2.id,st::water.id};
    const double phi=0.2;
    auto history=history_1p(closure,layout,ids,q,phi);

    auto partition=one_cell_partition();
    auto schedule=empty_schedule();
    auto bridge=one_cell_bridge();
    auto pattern=one_cell_pattern();
    fdp::PostSnesSw92ProfileCPhaseTransitionScannerContext3D scanner{&model,{}};
    IdentityResolverContext identity;
    fdp::Sw92TransactionalPhaseTransitionRestartContext3D<Provider> rebuild;
    rebuild.model=&model;
    rebuild.provider_factory=[&model](){return Provider{model};};
    rebuild.property_provenance=Provider::provenance();
    rebuild.scanner=&scanner;
    rebuild.schedule=&schedule; rebuild.partition=&partition;
    rebuild.cell_bridge=&bridge; rebuild.cell_pattern=&pattern;
    rebuild.rock={&zero_rock,nullptr};
    rebuild.two_phase_relative_permeability={&linear_kr,nullptr};
    rebuild.three_phase_saturation={&no_three_phase,nullptr};
    rebuild.target_identity={&resolve_identity,&identity};

    fdp::Sw92TransactionalInitialCell3D initial{
        {mesh::LocalIndex{0U},mesh::GlobalEntityId{UINT64_C(10)},1.0,phi,ids,
         history.first,history.second},
        layout.descriptor(),q,
        flow::FrozenActivePhaseIdentityMap{{identity.nonaqueous}},
        {selection}};
    verify_transaction(rebuild,scanner,std::move(initial),2U);
}

void sourced_two_to_one() {
    auto model=binary_model();
    using Provider=flow::Sw92Co2WaterPropertyProvider;
    const double p=3.0e6,t=340.0;
    const std::array<double,2> z{0.995,0.005};
    const auto selection=na_selection(model,p,t,z);
    auto closure=flow::make_sw92_co2_water_selected_phase_property_closure(
        model,std::vector<th::Sw92SelectedPhase<double>>{selection,selection});

    flow::NaturalVariableLayout2P layout{
        flow::NaturalVariableCompositionPivot2P::from_dependent_components(
            2U,{0U,0U})};
    std::vector<double> q(layout.unknown_count(),0.0);
    q[layout.pressure_unknown_index()]=p;
    q[layout.temperature_unknown_index()]=t;
    q[layout.independent_saturation_unknown_index()]=0.5;
    q[*layout.independent_composition_unknown_index(0U,1U)]=z[1];
    q[*layout.independent_composition_unknown_index(1U,1U)]=z[1];
    const std::vector<std::string> ids{st::co2.id,st::water.id};
    const double phi=0.2;
    auto history=history_2p(closure,layout,ids,q,phi);

    auto partition=one_cell_partition();
    auto schedule=empty_schedule();
    auto bridge=one_cell_bridge();
    auto pattern=one_cell_pattern();
    fdp::PostSnesSw92ProfileCPhaseTransitionScannerContext3D scanner{&model,{}};
    IdentityResolverContext identity;
    fdp::Sw92TransactionalPhaseTransitionRestartContext3D<Provider> rebuild;
    rebuild.model=&model;
    rebuild.provider_factory=[&model](){return Provider{model};};
    rebuild.property_provenance=Provider::provenance();
    rebuild.scanner=&scanner;
    rebuild.schedule=&schedule; rebuild.partition=&partition;
    rebuild.cell_bridge=&bridge; rebuild.cell_pattern=&pattern;
    rebuild.rock={&zero_rock,nullptr};
    rebuild.two_phase_relative_permeability={&linear_kr,nullptr};
    rebuild.three_phase_saturation={&no_three_phase,nullptr};
    rebuild.target_identity={&resolve_identity,&identity};

    fdp::Sw92TransactionalInitialCell3D initial{
        {mesh::LocalIndex{0U},mesh::GlobalEntityId{UINT64_C(10)},1.0,phi,ids,
         history.first,history.second},
        layout.descriptor(),q,
        flow::FrozenActivePhaseIdentityMap{{identity.extra,identity.nonaqueous}},
        {selection,selection}};
    verify_transaction(rebuild,scanner,std::move(initial),1U);
}

} // namespace

void sw92_transactional_phase_transition_restart_test() {
    require_tx(sw92_transactional_phase_transition_restart_header(),
               "SW92 transactional restart public header probe failed");
    sourced_one_to_two();
    sourced_two_to_one();
}
