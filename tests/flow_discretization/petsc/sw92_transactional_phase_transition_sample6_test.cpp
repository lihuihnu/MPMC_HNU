#include <mpmc/flow_discretization_petsc/sw92_transactional_phase_transition_restart.hpp>

#include "../../flash/sw92_phase_assigned_pt/physical_sample6.hpp"

#include <petscsys.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <memory>
#include <optional>
#include <span>
#include <stdexcept>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

namespace {

namespace dp = mpmc::discretization_petsc;
namespace fdp = mpmc::flow_discretization_petsc;
namespace flow = mpmc::flow;
namespace fl = mpmc::flash;
namespace mesh = mpmc::mesh;
namespace sample6 = sw92_profile_c_sample6;
namespace th = mpmc::thermodynamics;

void require_sample6_transaction(
    bool condition,
    std::string_view message) {
    if (!condition) {
        throw std::runtime_error(
            std::string{message});
    }
}

th::Provenance synthetic_mass_source(
    std::string locator) {
    return {
        th::SourceKind::synthetic_test,
        "MPMC_HNU manufactured SW92 Sample-6 transactional fixture",
        "v1",
        std::move(locator),
        "Nonphysical positive molar mass used only to exercise selected-phase mass density and PETSc transaction mechanics.",
        "Constructed inside the regression from the published Sample-6 thermodynamic snapshot.",
        "Test-only generated fixture"};
}

th::Sw92ParameterSet
sample6_transaction_parameters() {
    const auto source =
        sample6::parameters();

    std::vector<th::Component> catalog{
        source.components().items().begin(),
        source.components().items().end()};
    std::vector<std::string> order;
    order.reserve(catalog.size());
    for (std::size_t i = 0U;
         i < catalog.size();
         ++i) {
        order.push_back(catalog[i].id);
        catalog[i].molar_mass =
            th::SourcedScalar{
                0.020 +
                    0.005 *
                        static_cast<double>(i),
                th::Unit::kilogram_per_mole,
                synthetic_mass_source(
                    catalog[i].id +
                    std::string{
                        " manufactured molar mass"}),
                "manufactured kg/mol",
                "identity"};
    }

    th::Sw92ParameterInput input;
    input.model_id =
        std::string{
            th::sw92_corrected_profile};
    input.dataset_id =
        source.dataset_id();
    input.revision =
        source.revision();
    input.applicability =
        source.applicability();
    input.pure.assign(
        source.pure_records().begin(),
        source.pure_records().end());
    input.water_binary.assign(
        source.water_binary_records().begin(),
        source.water_binary_records().end());
    input.nonwater_binary.assign(
        source.nonwater_binary_records().begin(),
        source.nonwater_binary_records().end());

    return th::Sw92ParameterSet::create(
        catalog,
        order,
        input,
        th::DataPolicy::allow_synthetic_tests);
}

struct ManufacturedSample6TransportCaloricProvider {
    template <typename Number>
    [[nodiscard]]
    flow::SelectedPhaseTransportCaloricValues<Number>
    operator()(
        const th::Sw92SelectedPhase<double>&,
        const Number& pressure_pa,
        const Number&,
        std::span<const Number>,
        const Number&,
        const Number&,
        const Number& mass_density_kg_per_m3) const {
        // Test-only transport/caloric closure.  The production selected-phase
        // closure computes u=h-p/rho, so h=u0+p/rho makes u exactly topology-
        // independent.  This isolates the transactional rebuild from missing
        // eight-component transport/caloric source data.
        return {
            Number{1.0e-5},
            Number{2.5e5} +
                pressure_pa /
                    mass_density_kg_per_m3};
    }
};

flow::SelectedPhasePropertyProvenance
manufactured_sample6_provenance() {
    return {
        {
            "SW92 selected molar density x manufactured test molar mass",
            "sample6-transaction-manufactured",
            "v1"},
        {
            "manufactured constant viscosity",
            "sample6-transaction-manufactured",
            "v1"},
        {
            "manufactured h=u0+p/rho",
            "sample6-transaction-manufactured",
            "v1"},
        {
            "production identity u=h-p/rho",
            "sample6-transaction-manufactured",
            "v1"}};
}

mesh::PartitionSnapshot
sample6_partition() {
    mesh::Topology topology{
        {
            {},
            {},
            {},
            {mesh::GlobalEntityId{
                UINT64_C(9201)}}},
        {}};
    mesh::EntityOwnerRanks owners;
    owners.cells = {
        mesh::PartitionRank{0U}};
    return mesh::PartitionSnapshot::create(
        topology,
        mesh::PartitionRank{0U},
        1U,
        std::move(owners));
}

dp::ParallelOwnedConnectionSchedule3D
sample6_schedule() {
    return {
        mesh::PartitionRank{0U},
        1U,
        {},
        {}};
}

dp::PetscMpiAijSymbolicPreallocation3D
sample6_bridge() {
    return {
        mesh::PartitionRank{0U},
        1U,
        0,
        1,
        1,
        {mesh::LocalIndex{0U}},
        {mesh::GlobalEntityId{
            UINT64_C(9201)}},
        {0},
        {1},
        {0},
        {0}};
}

dp::OwnedCellStructuralColumnPatternSnapshot3D
sample6_pattern() {
    return {
        mesh::PartitionRank{0U},
        1U,
        1U,
        0,
        1,
        1,
        {{
            mesh::LocalIndex{0U},
            mesh::GlobalEntityId{
                UINT64_C(9201)},
            0,
            0U,
            1U,
            0U,
            0U}},
        {0},
        {}};
}

PetscErrorCode sample6_zero_rock(
    const flow::NaturalVariableLayoutDescriptor& layout,
    std::span<const double>,
    void*,
    std::optional<
        fdp::Sw92RockThermalStorageLinearization3D>* output,
    fdp::NaturalVariableSnesEvaluationStatus3D* status) {
    if (output == nullptr ||
        status == nullptr) {
        return PETSC_ERR_ARG_NULL;
    }
    output->emplace(
        fdp::Sw92RockThermalStorageLinearization3D{
            0.0,
            std::vector<double>(
                layout.unknown_count(),
                0.0),
            {
                "zero rock storage",
                "sample6-transaction-manufactured",
                "v1"}});
    *status =
        fdp::NaturalVariableSnesEvaluationStatus3D::
            success;
    return PETSC_SUCCESS;
}

PetscErrorCode sample6_two_phase_kr(
    const flow::NaturalVariableCellState2P& state,
    void*,
    std::optional<
        fdp::Sw92TwoPhaseRelativePermeabilityLinearization3D>* output,
    fdp::NaturalVariableSnesEvaluationStatus3D* status) {
    if (output == nullptr ||
        status == nullptr) {
        return PETSC_ERR_ARG_NULL;
    }
    const std::size_t q =
        state.layout().unknown_count();
    output->emplace(
        fdp::Sw92TwoPhaseRelativePermeabilityLinearization3D{
            {1.0, 1.0},
            {
                std::vector<double>(q, 0.0),
                std::vector<double>(q, 0.0)}});
    *status =
        fdp::NaturalVariableSnesEvaluationStatus3D::
            success;
    return PETSC_SUCCESS;
}

struct Sample6LinearThreePhaseKr {
    template <typename Number>
    [[nodiscard]]
    flow::RelativePermeabilityEvaluation3P<Number>
    operator()(
        const flow::ThreePhaseSaturationState3P<Number>& state)
        const {
        return {state.saturation};
    }
};

PetscErrorCode sample6_three_phase_saturation(
    const flow::NaturalVariableCellState3P& state,
    void*,
    std::optional<
        flow::
            ThreePhaseSaturationConstitutiveNaturalVariableLinearization3P>*
                output,
    fdp::NaturalVariableSnesEvaluationStatus3D* status) {
    if (output == nullptr ||
        status == nullptr) {
        return PETSC_ERR_ARG_NULL;
    }
    try {
        const auto primal =
            flow::
                evaluate_three_phase_saturation_constitutive(
                    state,
                    Sample6LinearThreePhaseKr{},
                    flow::NoCapillaryPressure3P{});
        flow::
            ThreePhaseSaturationCoordinateDerivatives3P
                derivatives{};
        derivatives.relative_permeability[0] =
            {1.0, 0.0};
        derivatives.relative_permeability[1] =
            {0.0, 1.0};
        derivatives.relative_permeability[2] =
            {-1.0, -1.0};
        output->emplace(
            flow::
                make_saturation_constitutive_natural_variable_linearization(
                    state,
                    primal,
                    derivatives));
        *status =
            fdp::NaturalVariableSnesEvaluationStatus3D::
                success;
    } catch (const std::exception&) {
        *status =
            fdp::NaturalVariableSnesEvaluationStatus3D::
                domain_error;
    }
    return PETSC_SUCCESS;
}

PetscErrorCode sample6_no_absent_provider(
    const flow::
        AbsentPhaseThermodynamicCoordinateExtension&,
    void* user_context,
    std::optional<
        flow::
            AbsentPhasePotentialExtensionLinearization>*,
    fdp::NaturalVariableSnesEvaluationStatus3D*) {
    return user_context == nullptr
        ? PETSC_ERR_ARG_NULL
        : PETSC_ERR_SUP;
}

flow::FrozenPhysicalPhaseIdentity
sample6_identity(
    std::string key) {
    return {
        "SW92/Profile-C/sample6-transaction-test",
        std::move(key)};
}

PetscErrorCode sample6_target_identity(
    mesh::GlobalEntityId,
    const fdp::Sw92AuthoritativeTargetMaterialization3D& target,
    const flow::FrozenActivePhaseIdentityMap& source,
    void*,
    std::optional<
        flow::FrozenActivePhaseIdentityMap>* output) {
    if (output == nullptr) {
        return PETSC_ERR_ARG_NULL;
    }
    output->reset();

    const auto aq =
        sample6_identity("aqueous");
    const auto na =
        sample6_identity("nonaqueous");
    const auto extra =
        sample6_identity("extra");

    if (source.phase_count() == 2U &&
        target.phases.size() == 3U &&
        target.phases[0].metadata.thermodynamic_family ==
            th::SwPhaseFamily::aqueous &&
        target.phases[1].metadata.thermodynamic_family ==
            th::SwPhaseFamily::nonaqueous &&
        target.phases[2].metadata.thermodynamic_family ==
            th::SwPhaseFamily::nonaqueous &&
        source.find(aq).has_value() &&
        source.find(na).has_value()) {
        output->emplace(
            std::vector<
                flow::FrozenPhysicalPhaseIdentity>{
                aq,
                na,
                extra});
        return PETSC_SUCCESS;
    }

    if (source.phase_count() == 3U &&
        target.phases.size() == 2U &&
        target.phases[0].metadata.thermodynamic_family ==
            th::SwPhaseFamily::aqueous &&
        target.phases[1].metadata.thermodynamic_family ==
            th::SwPhaseFamily::nonaqueous &&
        source.find(aq).has_value() &&
        source.find(na).has_value()) {
        output->emplace(
            std::vector<
                flow::FrozenPhysicalPhaseIdentity>{
                aq,
                na});
        return PETSC_SUCCESS;
    }

    return PETSC_ERR_ARG_INCOMP;
}

std::size_t dependent_largest(
    std::span<const double> composition) {
    const auto found =
        std::max_element(
            composition.begin(),
            composition.end());
    require_sample6_transaction(
        found != composition.end() &&
            std::isfinite(*found) &&
            *found > 0.0,
        "invalid Sample-6 composition");
    return static_cast<std::size_t>(
        std::distance(
            composition.begin(),
            found));
}

std::vector<double> sample6_edge_feed(
    std::span<const double> first,
    std::span<const double> second,
    double second_fraction) {
    require_sample6_transaction(
        first.size() == second.size() &&
            second_fraction > 0.0 &&
            second_fraction < 1.0,
        "invalid Sample-6 edge feed");
    std::vector<double> result(
        first.size(),
        0.0);
    for (std::size_t i = 0U;
         i < result.size();
         ++i) {
        result[i] =
            (1.0 - second_fraction) *
                first[i] +
            second_fraction *
                second[i];
    }
    return result;
}

template <typename Closure>
std::vector<double> two_phase_q(
    const flow::NaturalVariableLayout2P& layout,
    double pressure_pa,
    double temperature_k,
    const std::array<std::vector<double>, 2>& compositions,
    const std::array<double, 2>& mole_fractions,
    const Closure& closure) {
    std::vector<double> q(
        layout.unknown_count(),
        0.0);
    q[layout.pressure_unknown_index()] =
        pressure_pa;
    q[layout.temperature_unknown_index()] =
        temperature_k;

    std::array<double, 2> weight{};
    for (std::size_t phase = 0U;
         phase < 2U;
         ++phase) {
        const auto properties =
            closure.evaluate(
                phase,
                pressure_pa,
                temperature_k,
                std::span<const double>{
                    compositions[phase]});
        weight[phase] =
            mole_fractions[phase] /
            properties
                .molar_density_mol_per_m3;
    }
    q[layout.independent_saturation_unknown_index()] =
        weight[0] /
        (weight[0] + weight[1]);

    for (std::size_t phase = 0U;
         phase < 2U;
         ++phase) {
        for (std::size_t component = 0U;
             component <
                 compositions[phase].size();
             ++component) {
            const auto column =
                layout
                    .independent_composition_unknown_index(
                        phase,
                        component);
            if (column.has_value()) {
                q[*column] =
                    compositions[phase][component];
            }
        }
    }
    return q;
}

template <typename Closure>
std::vector<double> three_phase_q(
    const flow::NaturalVariableLayout3P& layout,
    double pressure_pa,
    double temperature_k,
    const std::array<std::vector<double>, 3>& compositions,
    const std::array<double, 3>& mole_fractions,
    const Closure& closure) {
    std::vector<double> q(
        layout.unknown_count(),
        0.0);
    q[layout.pressure_unknown_index()] =
        pressure_pa;
    q[layout.temperature_unknown_index()] =
        temperature_k;

    std::array<double, 3> weight{};
    double sum = 0.0;
    for (std::size_t phase = 0U;
         phase < 3U;
         ++phase) {
        const auto properties =
            closure.evaluate(
                phase,
                pressure_pa,
                temperature_k,
                std::span<const double>{
                    compositions[phase]});
        weight[phase] =
            mole_fractions[phase] /
            properties
                .molar_density_mol_per_m3;
        sum += weight[phase];
    }
    for (std::size_t phase = 0U;
         phase < 2U;
         ++phase) {
        const auto column =
            layout
                .independent_saturation_unknown_index(
                    static_cast<
                        flow::PhaseSlot3>(
                            phase));
        require_sample6_transaction(
            column.has_value(),
            "missing Sample-6 3P saturation column");
        q[*column] =
            weight[phase] / sum;
    }

    for (std::size_t phase = 0U;
         phase < 3U;
         ++phase) {
        const auto slot =
            static_cast<
                flow::PhaseSlot3>(
                    phase);
        for (std::size_t component = 0U;
             component <
                 compositions[phase].size();
             ++component) {
            const auto column =
                layout
                    .independent_composition_unknown_index(
                        slot,
                        component);
            if (column.has_value()) {
                q[*column] =
                    compositions[phase][component];
            }
        }
    }
    return q;
}

template <typename Provider>
void run_controller(
    std::unique_ptr<
        fdp::PhaseTransitionRebuiltNaturalVariableSystem3D>
        initial_system,
    fdp::PostSnesSw92ProfileCPhaseTransitionScannerContext3D*
        scanner,
    fdp::Sw92TransactionalPhaseTransitionRebuildContext3D<
        Provider>* rebuild,
    std::size_t source_count,
    std::size_t target_count) {
    std::unique_ptr<
        fdp::PhaseTransitionRebuiltNaturalVariableSystem3D>
        final_system;
    Vec final_state = nullptr;
    std::optional<
        fdp::PostSnesPhaseTransitionControllerReport3D>
        report;

    const PetscErrorCode error =
        fdp::
            solve_nonlinear_timestep_with_phase_transitions_3d(
                PETSC_COMM_SELF,
                std::move(initial_system),
                {
                    &fdp::
                        scan_post_snes_sw92_profile_c_phase_transitions_3d,
                    scanner,
                    &fdp::
                        rebuild_sw92_transactional_phase_transition_system_3d<
                            Provider>,
                    rebuild},
                {.max_transition_restarts = 4U},
                &final_system,
                &final_state,
                &report);

    require_sample6_transaction(
        error == PETSC_SUCCESS &&
            final_system != nullptr &&
            final_state != nullptr &&
            report.has_value() &&
            report->outcome ==
                fdp::
                    PostSnesPhaseTransitionOutcome3D::
                        stable_phase_set &&
            report->transition_restarts == 1U &&
            report->generations.size() == 2U &&
            report->generations[0]
                    .accepted_transition_batch
                    .size() == 1U &&
            report->generations[0]
                    .accepted_transition_batch[0]
                    .source_phase_count ==
                source_count &&
            report->generations[0]
                    .accepted_transition_batch[0]
                    .target_phase_count ==
                target_count &&
            report->generations[1]
                    .accepted_transition_batch
                    .empty() &&
            final_system->numbering()
                    .cell(mesh::LocalIndex{0U})
                    .phase_count ==
                target_count &&
            final_system->time_step_seconds() == 1.0,
        "Sample-6 SW92 same-dt transactional lifecycle failed");

    bool history_matches = false;
    require_sample6_transaction(
        final_system
                ->accepted_history_matches_state(
                    final_state,
                    &history_matches) ==
                PETSC_SUCCESS &&
            history_matches,
        "Sample-6 SW92 restart changed frozen component/energy history");

    (void)VecDestroy(&final_state);
}

void sample6_two_to_three() {
    auto parameters =
        sample6_transaction_parameters();
    auto model =
        th::Sw92Phase<double>::
            from_parameters(parameters);

    constexpr double pressure_pa = 1.0e7;
    constexpr double temperature_k = 350.0;
    const auto feed =
        sample6::feed();
    const auto owned =
        fl::solve_sw92_phase_assigned_pt_boundary_aware(
            pressure_pa,
            temperature_k,
            feed,
            model,
            0.0);
    require_sample6_transaction(
        owned.base.c1.has_value() &&
            owned.base.c1->candidate() != nullptr,
        "Sample-6 C1 2P source is unavailable");
    const auto& c1 =
        *owned.base.c1->candidate();

    std::array<std::vector<double>, 2>
        compositions{
            c1.aqueous_phase.composition,
            c1.nonaqueous_phase.composition};
    std::array<double, 2>
        fractions{
            c1.aqueous_phase.mole_phase_fraction,
            c1.nonaqueous_phase.mole_phase_fraction};
    std::vector<
        th::Sw92SelectedPhase<double>>
        selections{
            {
                0.0,
                th::SwPhaseFamily::aqueous,
                c1.aqueous_phase.activity.branch,
                {}},
            {
                0.0,
                th::SwPhaseFamily::nonaqueous,
                c1.nonaqueous_phase.activity.branch,
                {}}};

    ManufacturedSample6TransportCaloricProvider
        provider;
    auto closure =
        flow::make_sw92_selected_phase_property_closure(
            model,
            selections,
            provider,
            manufactured_sample6_provenance());
    using Closure = decltype(closure);

    const flow::NaturalVariableLayout2P layout{
        flow::NaturalVariableCompositionPivot2P::
            from_dependent_components(
                feed.size(),
                {
                    dependent_largest(
                        compositions[0]),
                    dependent_largest(
                        compositions[1])})};
    const auto q =
        two_phase_q(
            layout,
            pressure_pa,
            temperature_k,
            compositions,
            fractions,
            closure);

    fdp::Sw92TwoPhaseProductionCellEvaluatorContext3D<
        Closure>
        source_context{
            &closure,
            {&sample6_two_phase_kr, nullptr},
            {&sample6_zero_rock, nullptr},
            {}};
    std::optional<
        fdp::TwoPhaseCurrentCellLinearization3D>
        source_current;
    fdp::NaturalVariableSnesEvaluationStatus3D
        source_status =
            fdp::NaturalVariableSnesEvaluationStatus3D::
                success;
    require_sample6_transaction(
        fdp::
            evaluate_sw92_two_phase_production_cell_3d<
                Closure>(
                mesh::LocalIndex{0U},
                mesh::GlobalEntityId{UINT64_C(9201)},
                q,
                layout,
                closure.component_ids(),
                &source_context,
                &source_current,
                &source_status) ==
                PETSC_SUCCESS &&
            source_status ==
                fdp::
                    NaturalVariableSnesEvaluationStatus3D::
                        success &&
            source_current.has_value(),
        "Sample-6 2P source evaluator failed");

    constexpr double porosity = 0.20;
    const auto previous_component =
        flow::build_two_phase_component_accumulation(
            source_current->state,
            porosity);
    const auto previous_energy =
        flow::build_two_phase_energy_accumulation_snapshot(
            source_current->state,
            porosity,
            source_current->transport,
            source_current->caloric,
            source_current->rock);

    auto partition =
        sample6_partition();
    auto schedule =
        sample6_schedule();
    auto bridge =
        sample6_bridge();
    auto pattern =
        sample6_pattern();

    fdp::CellScopedMixedCardinalityEvaluatorDispatcher3D
        initial_dispatcher{
            {},
            {{
                mesh::GlobalEntityId{UINT64_C(9201)},
                {
                    &fdp::
                        evaluate_sw92_two_phase_production_cell_3d<
                            Closure>,
                    &source_context}}},
            {}};

    fdp::FrozenPhaseTransitionRebuildCell3D
        cell;
    cell.cell = mesh::LocalIndex{0U};
    cell.cell_global =
        mesh::GlobalEntityId{UINT64_C(9201)};
    cell.bulk_volume_m3 = 1.0;
    cell.porosity = porosity;
    cell.component_ids.assign(
        closure.component_ids().begin(),
        closure.component_ids().end());
    cell.target_layout = layout.descriptor();
    cell.target_natural_variables = q;
    cell.target_active_phases =
        flow::FrozenActivePhaseIdentityMap{
            {
                sample6_identity("aqueous"),
                sample6_identity("nonaqueous")}};
    cell.previous_component_accumulation =
        previous_component;
    cell.previous_energy_accumulation =
        previous_energy;
    cell.transition_evidence_profile =
        "SW92/sample6/transaction-source-2P/v1";

    int absent_token = 1;
    std::unique_ptr<
        fdp::PhaseTransitionRebuiltNaturalVariableSystem3D>
        initial_system;
    require_sample6_transaction(
        fdp::
            rebuild_phase_transition_natural_variable_system_3d(
                PETSC_COMM_SELF,
                schedule,
                partition,
                bridge,
                pattern,
                1.0,
                {cell},
                {},
                initial_dispatcher.bindings(),
                {},
                &sample6_no_absent_provider,
                &absent_token,
                &initial_system) ==
                PETSC_SUCCESS &&
            initial_system != nullptr,
        "failed to build Sample-6 2P source system");

    fdp::PostSnesSw92ProfileCPhaseTransitionScannerContext3D
        scanner{&model, {}, {}};

    using Provider =
        ManufacturedSample6TransportCaloricProvider;
    fdp::Sw92TransactionalPhaseTransitionRebuildContext3D<
        Provider>
        rebuild;
    rebuild.comm = PETSC_COMM_SELF;
    rebuild.schedule = &schedule;
    rebuild.partition = &partition;
    rebuild.cell_bridge = &bridge;
    rebuild.cell_pattern = &pattern;
    rebuild.model = &model;
    rebuild.provider = provider;
    rebuild.provenance =
        manufactured_sample6_provenance();
    rebuild.scanner_context = &scanner;
    rebuild.target_identity = {
        &sample6_target_identity,
        nullptr};
    rebuild.initial_evaluators =
        initial_dispatcher.bindings();
    rebuild.two_phase_relative_permeability = {
        &sample6_two_phase_kr,
        nullptr};
    rebuild.three_phase_saturation = {
        &sample6_three_phase_saturation,
        nullptr};
    rebuild.rock_storage = {
        &sample6_zero_rock,
        nullptr};
    rebuild.baseline_cells.push_back(
        {
            mesh::LocalIndex{0U},
            mesh::GlobalEntityId{UINT64_C(9201)},
            1.0,
            porosity,
            cell.component_ids,
            previous_component,
            previous_energy});

    run_controller(
        std::move(initial_system),
        &scanner,
        &rebuild,
        2U,
        3U);
}

void sample6_three_to_two() {
    auto parameters =
        sample6_transaction_parameters();
    auto model =
        th::Sw92Phase<double>::
            from_parameters(parameters);

    constexpr double pressure_pa = 1.0e7;
    constexpr double temperature_k = 350.0;
    const auto feed =
        sample6_edge_feed(
            sample6::w(),
            sample6::h0(),
            0.45);
    const auto authoritative =
        fl::solve_sw92_profile_c_pt_phase_set(
            pressure_pa,
            temperature_k,
            feed,
            model,
            0.0);
    require_sample6_transaction(
        authoritative.solution.status ==
                fl::PtPhaseSetStatus::accepted &&
            authoritative
                .accepted_phase_set_published() &&
            authoritative.solution
                    .accepted_phase_count() == 2U &&
            authoritative.phase_metadata.size() == 2U,
        "Sample-6 W+H target is not authoritative 2P");
    const auto* accepted =
        authoritative.solution
            .accepted_phase_set();
    require_sample6_transaction(
        accepted != nullptr &&
            accepted->phases.size() == 2U,
        "Sample-6 W+H accepted phases missing");

    const auto& aq =
        accepted->phases[0];
    const auto& na =
        accepted->phases[1];
    std::array<std::vector<double>, 3>
        compositions{
            aq.composition,
            na.composition,
            na.composition};
    std::array<double, 3>
        fractions{
            aq.mole_phase_fraction,
            0.5 * na.mole_phase_fraction,
            0.5 * na.mole_phase_fraction};
    std::vector<
        th::Sw92SelectedPhase<double>>
        selections{
            {
                0.0,
                th::SwPhaseFamily::aqueous,
                aq.activity.branch,
                {}},
            {
                0.0,
                th::SwPhaseFamily::nonaqueous,
                na.activity.branch,
                {}},
            {
                0.0,
                th::SwPhaseFamily::nonaqueous,
                na.activity.branch,
                {}}};

    ManufacturedSample6TransportCaloricProvider
        provider;
    auto closure =
        flow::make_sw92_selected_phase_property_closure(
            model,
            selections,
            provider,
            manufactured_sample6_provenance());
    using Closure = decltype(closure);

    const flow::NaturalVariableLayout3P layout{
        flow::NaturalVariableCompositionPivot3P::
            from_dependent_components(
                feed.size(),
                {
                    dependent_largest(
                        compositions[0]),
                    dependent_largest(
                        compositions[1]),
                    dependent_largest(
                        compositions[2])})};
    const auto q =
        three_phase_q(
            layout,
            pressure_pa,
            temperature_k,
            compositions,
            fractions,
            closure);

    fdp::Sw92ThreePhaseProductionCellEvaluatorContext3D<
        Closure>
        source_context{
            &closure,
            {&sample6_three_phase_saturation, nullptr},
            {&sample6_zero_rock, nullptr},
            {}};
    std::optional<
        fdp::FixedThreePhaseCurrentCellLinearization3D>
        source_current;
    fdp::NaturalVariableSnesEvaluationStatus3D
        source_status =
            fdp::NaturalVariableSnesEvaluationStatus3D::
                success;
    require_sample6_transaction(
        fdp::
            evaluate_sw92_three_phase_production_cell_3d<
                Closure>(
                mesh::LocalIndex{0U},
                mesh::GlobalEntityId{UINT64_C(9201)},
                q,
                layout,
                closure.component_ids(),
                &source_context,
                &source_current,
                &source_status) ==
                PETSC_SUCCESS &&
            source_status ==
                fdp::
                    NaturalVariableSnesEvaluationStatus3D::
                        success &&
            source_current.has_value(),
        "Sample-6 stale 3P source evaluator failed");

    constexpr double porosity = 0.20;
    const auto previous_component =
        flow::build_pore_volume_component_accumulation(
            source_current->state,
            porosity);
    const auto previous_energy =
        flow::build_pore_volume_energy_accumulation_snapshot(
            source_current->state,
            porosity,
            source_current->transport,
            source_current->caloric,
            source_current->rock);

    auto partition =
        sample6_partition();
    auto schedule =
        sample6_schedule();
    auto bridge =
        sample6_bridge();
    auto pattern =
        sample6_pattern();

    fdp::CellScopedMixedCardinalityEvaluatorDispatcher3D
        initial_dispatcher{
            {},
            {},
            {{
                mesh::GlobalEntityId{UINT64_C(9201)},
                {
                    &fdp::
                        evaluate_sw92_three_phase_production_cell_3d<
                            Closure>,
                    &source_context}}}};

    fdp::FrozenPhaseTransitionRebuildCell3D
        cell;
    cell.cell = mesh::LocalIndex{0U};
    cell.cell_global =
        mesh::GlobalEntityId{UINT64_C(9201)};
    cell.bulk_volume_m3 = 1.0;
    cell.porosity = porosity;
    cell.component_ids.assign(
        closure.component_ids().begin(),
        closure.component_ids().end());
    cell.target_layout =
        flow::NaturalVariableLayoutDescriptor{
            layout};
    cell.target_natural_variables = q;
    cell.target_active_phases =
        flow::FrozenActivePhaseIdentityMap{
            {
                sample6_identity("aqueous"),
                sample6_identity("nonaqueous"),
                sample6_identity("extra")}};
    cell.previous_component_accumulation =
        previous_component;
    cell.previous_energy_accumulation =
        previous_energy;
    cell.transition_evidence_profile =
        "SW92/sample6/transaction-source-3P/v1";

    int absent_token = 1;
    std::unique_ptr<
        fdp::PhaseTransitionRebuiltNaturalVariableSystem3D>
        initial_system;
    require_sample6_transaction(
        fdp::
            rebuild_phase_transition_natural_variable_system_3d(
                PETSC_COMM_SELF,
                schedule,
                partition,
                bridge,
                pattern,
                1.0,
                {cell},
                {},
                initial_dispatcher.bindings(),
                {},
                &sample6_no_absent_provider,
                &absent_token,
                &initial_system) ==
                PETSC_SUCCESS &&
            initial_system != nullptr,
        "failed to build Sample-6 stale 3P source system");

    fdp::PostSnesSw92ProfileCPhaseTransitionScannerContext3D
        scanner{&model, {}, {}};

    using Provider =
        ManufacturedSample6TransportCaloricProvider;
    fdp::Sw92TransactionalPhaseTransitionRebuildContext3D<
        Provider>
        rebuild;
    rebuild.comm = PETSC_COMM_SELF;
    rebuild.schedule = &schedule;
    rebuild.partition = &partition;
    rebuild.cell_bridge = &bridge;
    rebuild.cell_pattern = &pattern;
    rebuild.model = &model;
    rebuild.provider = provider;
    rebuild.provenance =
        manufactured_sample6_provenance();
    rebuild.scanner_context = &scanner;
    rebuild.target_identity = {
        &sample6_target_identity,
        nullptr};
    rebuild.initial_evaluators =
        initial_dispatcher.bindings();
    rebuild.two_phase_relative_permeability = {
        &sample6_two_phase_kr,
        nullptr};
    rebuild.three_phase_saturation = {
        &sample6_three_phase_saturation,
        nullptr};
    rebuild.rock_storage = {
        &sample6_zero_rock,
        nullptr};
    rebuild.baseline_cells.push_back(
        {
            mesh::LocalIndex{0U},
            mesh::GlobalEntityId{UINT64_C(9201)},
            1.0,
            porosity,
            cell.component_ids,
            previous_component,
            previous_energy});

    run_controller(
        std::move(initial_system),
        &scanner,
        &rebuild,
        3U,
        2U);
}

} // namespace

void sw92_transactional_phase_transition_sample6_test() {
    sample6_two_to_three();
    sample6_three_to_two();
}
