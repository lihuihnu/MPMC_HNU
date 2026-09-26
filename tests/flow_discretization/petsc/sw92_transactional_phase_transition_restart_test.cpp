#include <mpmc/flow_discretization_petsc/sw92_transactional_phase_transition_restart.hpp>
#include <mpmc/flow/sw92_co2_water_properties.hpp>

#include "test_support.hpp"

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
#include <type_traits>
#include <utility>
#include <vector>

bool sw92_transactional_phase_transition_restart_header();

void sw92_transactional_phase_transition_sample6_test();

namespace {

namespace dp = mpmc::discretization_petsc;
namespace fdp = mpmc::flow_discretization_petsc;
namespace flow = mpmc::flow;
namespace mesh = mpmc::mesh;
namespace st = sw92_test;
namespace th = mpmc::thermodynamics;

void require_transaction(
    bool condition,
    std::string_view message) {
    if (!condition) {
        throw std::runtime_error(
            std::string{message});
    }
}

// Analytic affine-simplex bounds, including a non-last composition pivot.
// Owned by the existing PETSc executable/Gate; these are software invariants.
void conservation_anchor_feasible_step_test() {
    using fdp::sw92_transactional_restart_detail::conservation_anchor_feasible_damping;
    for (std::size_t phases = 1U; phases <= 3U; ++phases) {
        const flow::NaturalVariableLayoutDescriptor layout{
            3U, phases, std::vector<std::size_t>(phases, 0U)};
        std::vector<double> q(layout.unknown_count(), 0.2);
        q[0] = 1.0e7;
        q[1] = 350.0;
        std::vector<double> step(q.size(), 0.0);
        require_transaction(conservation_anchor_feasible_damping(layout, q, step) == 1.0,
            "zero direction must retain unit damping");
        const auto check = [&](double expected) {
            const double alpha = conservation_anchor_feasible_damping(layout, q, step);
            require_transaction(std::abs(alpha / expected - 1.0) < 1.0e-12,
                "incorrect analytic conservation-anchor feasibility bound");
            auto trial = q;
            for (std::size_t i = 0U; i < q.size(); ++i) {
                trial[i] += alpha * step[i];
            }
            const std::vector<double> zero(q.size(), 0.0);
            require_transaction(conservation_anchor_feasible_damping(layout, trial, zero) == 1.0,
                "bounded trial must retain strict-positive independent and dependent coordinates");
        };
        const auto composition = *layout.independent_composition_unknown_index(
            flow::PhaseSlot3::phase0, 1U);
        q[composition] = 7.454115499883185e-12;
        step[composition] = -0.013885768807024454;
        check(0.99 * q[composition] / -step[composition]);
        q[composition] = 0.2;
        step[composition] = 2.0;
        check(0.99 * 0.6 / 2.0); // Reconstructed pivot reaches zero first.
        step[composition] = 0.0;
        if (phases > 1U) {
            step[2] = -2.0;
            check(0.99 * 0.2 / 2.0);
            step[2] = 2.0;
            check(0.99 * (1.0 - 0.2 * static_cast<double>(phases - 1U)) / 2.0);
            step[2] = 0.0;
        }
        step[0] = -2.0e7;
        check(0.495);
        step[0] = 0.0;
        step[1] = -700.0;
        check(0.495);
        step[1] = 0.0;
        q[composition] = 0.0;
        bool rejected = false;
        try {
            (void)conservation_anchor_feasible_damping(layout, q, step);
        } catch (const std::range_error&) {
            rejected = true;
        }
        require_transaction(rejected, "zero composition must not be clipped into the interior");
    }
}

th::Provenance nist_transaction_source(
    std::string locator) {
    return {
        th::SourceKind::database,
        "NIST Chemistry WebBook SRD 69",
        "accessed-2026-09-25",
        std::move(locator),
        "",
        "Public NIST Chemistry WebBook record",
        "Use as source data; no NIST code redistributed"};
}

st::Prepared sourced_transaction_input() {
    auto prepared =
        st::binary_input(st::co2);
    for (auto& component : prepared.catalog) {
        const double value =
            component.id == st::co2.id
                ? 0.0440095
                : 0.0180153;
        component.molar_mass =
            st::scalar(
                value,
                th::Unit::kilogram_per_mole,
                nist_transaction_source(
                    component.id +
                    std::string{
                        " molecular weight"}),
                "kg/mol",
                "identity");
    }
    return prepared;
}

th::Sw92ParameterSet
sourced_transaction_parameters() {
    auto prepared =
        sourced_transaction_input();
    return th::Sw92ParameterSet::create(
        prepared.catalog,
        prepared.order,
        prepared.input);
}

mesh::PartitionSnapshot
transaction_partition() {
    mesh::Topology topology{
        {
            {},
            {},
            {},
            {mesh::GlobalEntityId{UINT64_C(9101)}}},
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
transaction_schedule() {
    return {
        mesh::PartitionRank{0U},
        1U,
        {},
        {}};
}

dp::PetscMpiAijSymbolicPreallocation3D
transaction_bridge() {
    return {
        mesh::PartitionRank{0U},
        1U,
        0,
        1,
        1,
        {mesh::LocalIndex{0U}},
        {mesh::GlobalEntityId{UINT64_C(9101)}},
        {0},
        {1},
        {0},
        {0}};
}

dp::OwnedCellStructuralColumnPatternSnapshot3D
transaction_pattern() {
    return {
        mesh::PartitionRank{0U},
        1U,
        1U,
        0,
        1,
        1,
        {{
            mesh::LocalIndex{0U},
            mesh::GlobalEntityId{UINT64_C(9101)},
            0,
            0U,
            1U,
            0U,
            0U}},
        {0},
        {}};
}

PetscErrorCode transaction_zero_rock(
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
                "sw92-transaction-restart",
                "v1"}});
    *status =
        fdp::NaturalVariableSnesEvaluationStatus3D::
            success;
    return PETSC_SUCCESS;
}

PetscErrorCode transaction_two_phase_kr(
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

PetscErrorCode transaction_three_phase_saturation(
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
            flow::evaluate_three_phase_saturation_constitutive(
                state,
                [](const auto& saturation) {
                    return flow::
                        RelativePermeabilityEvaluation3P<
                            std::decay_t<
                                decltype(
                                    saturation.saturation[0])>>{
                            saturation.saturation};
                },
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

PetscErrorCode transaction_target_identity(
    mesh::GlobalEntityId,
    const fdp::Sw92AuthoritativeTargetMaterialization3D&
        target,
    const flow::FrozenActivePhaseIdentityMap&
        source_active,
    void*,
    std::optional<
        flow::FrozenActivePhaseIdentityMap>*
            output) {
    if (output == nullptr) {
        return PETSC_ERR_ARG_NULL;
    }
    output->reset();
    if (source_active.phase_count() == 1U &&
        target.phases.size() == 2U &&
        target.phases[0]
                .metadata
                .thermodynamic_family ==
            th::SwPhaseFamily::aqueous &&
        target.phases[1]
                .metadata
                .thermodynamic_family ==
            th::SwPhaseFamily::nonaqueous) {
        output->emplace(
            std::vector<
                flow::FrozenPhysicalPhaseIdentity>{
                {
                    "SW92/Profile-C/transaction-test",
                    "aqueous"},
                source_active.identity(0U)});
        return PETSC_SUCCESS;
    }
    if (source_active.phase_count() == 2U &&
        target.phases.size() == 1U &&
        target.phases[0]
                .metadata
                .thermodynamic_family ==
            th::SwPhaseFamily::nonaqueous) {
        const flow::FrozenPhysicalPhaseIdentity
            survivor{
                "SW92/Profile-C/transaction-test",
                "nonaqueous"};
        if (!source_active.find(survivor).has_value()) {
            return PETSC_ERR_ARG_INCOMP;
        }
        output->emplace(
            std::vector<
                flow::FrozenPhysicalPhaseIdentity>{
                survivor});
        return PETSC_SUCCESS;
    }
    return PETSC_ERR_ARG_INCOMP;
}

PetscErrorCode transaction_no_absent_provider(
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

void sw92_same_dt_one_to_two_transaction() {
    auto parameters =
        sourced_transaction_parameters();
    auto model =
        th::Sw92Phase<double>::
            from_parameters(parameters);

    constexpr double pressure_pa = 3.0e6;
    constexpr double temperature_k = 340.0;
    constexpr std::array<double, 2>
        composition{0.70, 0.30};

    th::Sw92PhaseWorkspace<double>
        root_workspace;
    const auto roots =
        model.roots(
            pressure_pa,
            temperature_k,
            std::span<const double>{
                composition},
            0.0,
            th::SwPhaseFamily::nonaqueous,
            root_workspace);
    require_transaction(
        roots.status ==
                th::Sw92RootStatus::success &&
            roots.count > 0U,
        "SW92 transactional 1P source root is unavailable");

    const th::Sw92SelectedPhase<double>
        source_selection{
            0.0,
            th::SwPhaseFamily::nonaqueous,
            roots.count - 1U,
            {}};
    auto source_closure =
        flow::
            make_sw92_co2_water_selected_phase_property_closure(
                model,
                std::vector<
                    th::Sw92SelectedPhase<double>>{
                    source_selection});
    using Closure =
        decltype(source_closure);

    const flow::NaturalVariableLayout1P
        source_layout{
            flow::
                NaturalVariableCompositionPivot1P::
                    from_dependent_component(
                        2U,
                        0U)};
    std::vector<double> source_q{
        pressure_pa,
        temperature_k,
        composition[1]};

    fdp::Sw92SinglePhaseProductionCellEvaluatorContext3D<
        Closure>
        source_context{
            &source_closure,
            1.0,
            {&transaction_zero_rock, nullptr},
            {}};

    std::optional<
        fdp::SinglePhaseCurrentCellLinearization3D>
        source_current;
    fdp::NaturalVariableSnesEvaluationStatus3D
        source_status =
            fdp::NaturalVariableSnesEvaluationStatus3D::
                success;
    require_transaction(
        fdp::
            evaluate_sw92_single_phase_production_cell_3d<
                Closure>(
                mesh::LocalIndex{0U},
                mesh::GlobalEntityId{UINT64_C(9101)},
                source_q,
                source_layout,
                source_closure.component_ids(),
                &source_context,
                &source_current,
                &source_status) ==
                PETSC_SUCCESS &&
            source_status ==
                fdp::
                    NaturalVariableSnesEvaluationStatus3D::
                        success &&
            source_current.has_value(),
        "SW92 transactional source evaluator failed");

    constexpr double porosity = 0.20;
    const auto previous_component =
        flow::
            build_single_phase_component_accumulation(
                source_current->state,
                porosity);
    const auto previous_energy =
        flow::
            build_single_phase_energy_accumulation_snapshot(
                source_current->state,
                porosity,
                source_current->transport,
                source_current->caloric,
                source_current->rock);

    const flow::FrozenActivePhaseIdentityMap
        source_identity{
            {
                {
                    "SW92/Profile-C/transaction-test",
                    "nonaqueous"}}};

    auto partition =
        transaction_partition();
    auto schedule =
        transaction_schedule();
    auto bridge =
        transaction_bridge();
    auto pattern =
        transaction_pattern();

    fdp::CellScopedMixedCardinalityEvaluatorDispatcher3D
        initial_dispatcher{
            {
                {
                    mesh::GlobalEntityId{UINT64_C(9101)},
                    {
                        &fdp::
                            evaluate_sw92_single_phase_production_cell_3d<
                                Closure>,
                        &source_context}}
            },
            {},
            {}};

    fdp::FrozenPhaseTransitionRebuildCell3D
        initial_cell;
    initial_cell.cell =
        mesh::LocalIndex{0U};
    initial_cell.cell_global =
        mesh::GlobalEntityId{UINT64_C(9101)};
    initial_cell.bulk_volume_m3 = 1.0;
    initial_cell.porosity = porosity;
    initial_cell.component_ids.assign(
        source_closure
            .component_ids()
            .begin(),
        source_closure
            .component_ids()
            .end());
    initial_cell.target_layout =
        source_layout.descriptor();
    initial_cell.target_natural_variables =
        source_q;
    initial_cell.target_active_phases =
        source_identity;
    initial_cell.previous_component_accumulation =
        previous_component;
    initial_cell.previous_energy_accumulation =
        previous_energy;
    initial_cell.transition_evidence_profile =
        "SW92/transaction/source-1P/v1";

    int absent_token = 1;
    std::unique_ptr<
        fdp::PhaseTransitionRebuiltNaturalVariableSystem3D>
        initial_system;
    require_transaction(
        fdp::
            rebuild_phase_transition_natural_variable_system_3d(
                PETSC_COMM_SELF,
                schedule,
                partition,
                bridge,
                pattern,
                1.0,
                {initial_cell},
                {},
                initial_dispatcher.bindings(),
                {},
                &transaction_no_absent_provider,
                &absent_token,
                &initial_system) ==
                PETSC_SUCCESS &&
            initial_system != nullptr,
        "failed to build SW92 transactional 1P source system");

    fdp::PostSnesSw92ProfileCPhaseTransitionScannerContext3D
        scanner{
            &model,
            {},
            {}};

    using Provider =
        flow::Sw92Co2WaterPropertyProvider;
    fdp::Sw92TransactionalPhaseTransitionRebuildContext3D<
        Provider>
        rebuild_context{
            PETSC_COMM_SELF,
            &schedule,
            &partition,
            &bridge,
            &pattern,
            &model,
            Provider{model},
            Provider::provenance(),
            &scanner,
            {&transaction_target_identity, nullptr},
            initial_dispatcher.bindings(),
            1.0,
            {&transaction_two_phase_kr, nullptr},
            {&transaction_three_phase_saturation, nullptr},
            {&transaction_zero_rock, nullptr},
            {},
            {{
                mesh::LocalIndex{0U},
                mesh::GlobalEntityId{UINT64_C(9101)},
                1.0,
                porosity,
                initial_cell.component_ids,
                previous_component,
                previous_energy}},
            {},
            1};

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
                    &scanner,
                    &fdp::
                        rebuild_sw92_transactional_phase_transition_system_3d<
                            Provider>,
                    &rebuild_context},
                {.max_transition_restarts = 4U},
                &final_system,
                &final_state,
                &report);
    require_transaction(
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
                    .source_phase_count == 1U &&
            report->generations[0]
                    .accepted_transition_batch[0]
                    .target_phase_count == 2U &&
            report->generations[1]
                    .accepted_transition_batch
                    .empty(),
        "SW92 transactional 1P->2P same-dt controller lifecycle failed");

    require_transaction(
        final_system
                ->numbering()
                .cell(
                    mesh::LocalIndex{0U})
                .phase_count == 2U &&
            final_system
                    ->time_step_seconds() ==
                1.0 &&
            rebuild_context.runtimes.size() ==
                1U &&
            scanner.local_owned_targets.empty(),
        "SW92 transactional restart did not finish on the rebuilt 2P same-dt system");

    bool history_matches = false;
    require_transaction(
        final_system
                ->accepted_history_matches_state(
                    final_state,
                    &history_matches) ==
                PETSC_SUCCESS &&
            history_matches,
        "SW92 same-dt restarted solution no longer matches the frozen source component/energy history");

    std::vector<std::optional<
        fdp::MixedCardinalityPhysicalCurrentCellLinearization3D>>
        current;
    std::vector<double> porosities;
    fdp::NaturalVariableSnesEvaluationStatus3D
        status =
            fdp::NaturalVariableSnesEvaluationStatus3D::
                success;
    require_transaction(
        final_system
                ->evaluate_local_cells_for_phase_transition(
                    final_state,
                    &current,
                    &porosities,
                    &status) ==
                PETSC_SUCCESS &&
            status ==
                fdp::
                    NaturalVariableSnesEvaluationStatus3D::
                        success &&
            current.size() == 1U &&
            current[0].has_value(),
        "failed to evaluate final SW92 transactional 2P state");
    const auto final_inventory =
        fdp::
            post_snes_pt_flash_scanner_detail::
                current_component_accumulation(
                    *current[0],
                    porosities[0]);
    const double total =
        final_inventory
            .total_accumulation_mol_per_bulk_m3;
    require_transaction(
        total > 0.0,
        "final SW92 transactional total component inventory is invalid");
    for (std::size_t component = 0U;
         component < 2U;
         ++component) {
        const double fraction =
            final_inventory
                .component_accumulation_mol_per_bulk_m3[
                    component] /
            total;
        require_transaction(
            std::abs(
                fraction -
                composition[component]) <
                2.0e-8,
            "SW92 transactional same-dt restart changed overall composition");
    }

    (void)VecDestroy(&final_state);
}

void sw92_same_dt_two_to_one_transaction() {
    auto parameters =
        sourced_transaction_parameters();
    auto model =
        th::Sw92Phase<double>::
            from_parameters(parameters);

    constexpr double pressure_pa = 3.0e6;
    constexpr double temperature_k = 340.0;
    constexpr std::array<double, 2>
        composition{0.995, 0.005};

    th::Sw92PhaseWorkspace<double>
        root_workspace;
    const auto roots =
        model.roots(
            pressure_pa,
            temperature_k,
            std::span<const double>{
                composition},
            0.0,
            th::SwPhaseFamily::nonaqueous,
            root_workspace);
    require_transaction(
        roots.status ==
                th::Sw92RootStatus::success &&
            roots.count > 0U,
        "SW92 transactional 2P source root is unavailable");

    const th::Sw92SelectedPhase<double>
        source_selection{
            0.0,
            th::SwPhaseFamily::nonaqueous,
            roots.count - 1U,
            {}};
    auto source_closure =
        flow::
            make_sw92_co2_water_selected_phase_property_closure(
                model,
                std::vector<
                    th::Sw92SelectedPhase<double>>{
                    source_selection,
                    source_selection});
    using Closure =
        decltype(source_closure);

    const flow::NaturalVariableLayout2P
        source_layout{
            flow::
                NaturalVariableCompositionPivot2P::
                    from_dependent_components(
                        2U,
                        {0U, 0U})};
    std::vector<double> source_q(
        source_layout.unknown_count(),
        0.0);
    source_q[
        source_layout.pressure_unknown_index()] =
        pressure_pa;
    source_q[
        source_layout.temperature_unknown_index()] =
        temperature_k;
    source_q[
        source_layout.independent_saturation_unknown_index()] =
        0.5;
    source_q[
        *source_layout
             .independent_composition_unknown_index(
                 0U,
                 1U)] =
        composition[1];
    source_q[
        *source_layout
             .independent_composition_unknown_index(
                 1U,
                 1U)] =
        composition[1];

    fdp::Sw92TwoPhaseProductionCellEvaluatorContext3D<
        Closure>
        source_context{
            &source_closure,
            {&transaction_two_phase_kr, nullptr},
            {&transaction_zero_rock, nullptr},
            {}};

    std::optional<
        fdp::TwoPhaseCurrentCellLinearization3D>
        source_current;
    fdp::NaturalVariableSnesEvaluationStatus3D
        source_status =
            fdp::NaturalVariableSnesEvaluationStatus3D::
                success;
    require_transaction(
        fdp::
            evaluate_sw92_two_phase_production_cell_3d<
                Closure>(
                mesh::LocalIndex{0U},
                mesh::GlobalEntityId{UINT64_C(9101)},
                source_q,
                source_layout,
                source_closure.component_ids(),
                &source_context,
                &source_current,
                &source_status) ==
                PETSC_SUCCESS &&
            source_status ==
                fdp::
                    NaturalVariableSnesEvaluationStatus3D::
                        success &&
            source_current.has_value(),
        "SW92 transactional 2P source evaluator failed");

    constexpr double porosity = 0.20;
    const auto previous_component =
        flow::
            build_two_phase_component_accumulation(
                source_current->state,
                porosity);
    const auto previous_energy =
        flow::
            build_two_phase_energy_accumulation_snapshot(
                source_current->state,
                porosity,
                source_current->transport,
                source_current->caloric,
                source_current->rock);

    const flow::FrozenActivePhaseIdentityMap
        source_identity{
            {
                {
                    "SW92/Profile-C/transaction-test",
                    "extra"},
                {
                    "SW92/Profile-C/transaction-test",
                    "nonaqueous"}}};

    auto partition =
        transaction_partition();
    auto schedule =
        transaction_schedule();
    auto bridge =
        transaction_bridge();
    auto pattern =
        transaction_pattern();

    fdp::CellScopedMixedCardinalityEvaluatorDispatcher3D
        initial_dispatcher{
            {},
            {
                {
                    mesh::GlobalEntityId{UINT64_C(9101)},
                    {
                        &fdp::
                            evaluate_sw92_two_phase_production_cell_3d<
                                Closure>,
                        &source_context}}
            },
            {}};

    fdp::FrozenPhaseTransitionRebuildCell3D
        initial_cell;
    initial_cell.cell =
        mesh::LocalIndex{0U};
    initial_cell.cell_global =
        mesh::GlobalEntityId{UINT64_C(9101)};
    initial_cell.bulk_volume_m3 = 1.0;
    initial_cell.porosity = porosity;
    initial_cell.component_ids.assign(
        source_closure
            .component_ids()
            .begin(),
        source_closure
            .component_ids()
            .end());
    initial_cell.target_layout =
        source_layout.descriptor();
    initial_cell.target_natural_variables =
        source_q;
    initial_cell.target_active_phases =
        source_identity;
    initial_cell.previous_component_accumulation =
        previous_component;
    initial_cell.previous_energy_accumulation =
        previous_energy;
    initial_cell.transition_evidence_profile =
        "SW92/transaction/source-2P/v1";

    int absent_token = 1;
    std::unique_ptr<
        fdp::PhaseTransitionRebuiltNaturalVariableSystem3D>
        initial_system;
    require_transaction(
        fdp::
            rebuild_phase_transition_natural_variable_system_3d(
                PETSC_COMM_SELF,
                schedule,
                partition,
                bridge,
                pattern,
                1.0,
                {initial_cell},
                {},
                initial_dispatcher.bindings(),
                {},
                &transaction_no_absent_provider,
                &absent_token,
                &initial_system) ==
                PETSC_SUCCESS &&
            initial_system != nullptr,
        "failed to build SW92 transactional 2P source system");

    fdp::PostSnesSw92ProfileCPhaseTransitionScannerContext3D
        scanner{
            &model,
            {},
            {}};

    using Provider =
        flow::Sw92Co2WaterPropertyProvider;
    fdp::Sw92TransactionalPhaseTransitionRebuildContext3D<
        Provider>
        rebuild_context{
            PETSC_COMM_SELF,
            &schedule,
            &partition,
            &bridge,
            &pattern,
            &model,
            Provider{model},
            Provider::provenance(),
            &scanner,
            {&transaction_target_identity, nullptr},
            initial_dispatcher.bindings(),
            1.0,
            {&transaction_two_phase_kr, nullptr},
            {&transaction_three_phase_saturation, nullptr},
            {&transaction_zero_rock, nullptr},
            {},
            {{
                mesh::LocalIndex{0U},
                mesh::GlobalEntityId{UINT64_C(9101)},
                1.0,
                porosity,
                initial_cell.component_ids,
                previous_component,
                previous_energy}},
            {},
            1};

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
                    &scanner,
                    &fdp::
                        rebuild_sw92_transactional_phase_transition_system_3d<
                            Provider>,
                    &rebuild_context},
                {.max_transition_restarts = 4U},
                &final_system,
                &final_state,
                &report);
    require_transaction(
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
                    .source_phase_count == 2U &&
            report->generations[0]
                    .accepted_transition_batch[0]
                    .target_phase_count == 1U &&
            report->generations[1]
                    .accepted_transition_batch
                    .empty(),
        "SW92 transactional 2P->1P same-dt controller lifecycle failed");

    require_transaction(
        final_system
                ->numbering()
                .cell(
                    mesh::LocalIndex{0U})
                .phase_count == 1U &&
            final_system
                    ->time_step_seconds() ==
                1.0 &&
            rebuild_context.runtimes.size() ==
                1U &&
            scanner.local_owned_targets.empty(),
        "SW92 transactional restart did not finish on the rebuilt 1P same-dt system");

    bool history_matches = false;
    require_transaction(
        final_system
                ->accepted_history_matches_state(
                    final_state,
                    &history_matches) ==
                PETSC_SUCCESS &&
            history_matches,
        "SW92 2P->1P same-dt restarted solution no longer matches the frozen source component/energy history");

    (void)VecDestroy(&final_state);
}

} // namespace

void sw92_transactional_phase_transition_restart_test() {
    require_transaction(
        sw92_transactional_phase_transition_restart_header(),
        "SW92 transactional restart public header probe failed");
    conservation_anchor_feasible_step_test();
    sw92_same_dt_one_to_two_transaction();
    sw92_same_dt_two_to_one_transaction();
    sw92_transactional_phase_transition_sample6_test();
}
