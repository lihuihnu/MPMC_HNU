#include <mpmc/discretization_petsc/adapter.hpp>
#include <mpmc/flash/sw92_profile_c_phase_set.hpp>
#include <mpmc/flow/sw92_co2_water_properties.hpp>
#include <mpmc/flow_discretization_petsc/mixed_cardinality_physical_snes_assembly.hpp>
#include <mpmc/flow_discretization_petsc/sw92_production_cell_evaluator.hpp>
#include <mpmc/flow_discretization_petsc/variable_cardinality_natural_variable_numbering.hpp>
#include <mpmc/flow_discretization_petsc/variable_cardinality_snes_assembly.hpp>
#include <mpmc/mesh/partition_snapshot.hpp>
#include <mpmc/mesh/topology.hpp>

#include "test_support.hpp"

#include <petscmat.h>
#include <petscvec.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <iostream>
#include <optional>
#include <span>
#include <stdexcept>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

bool selected_phase_production_cell_evaluator_header();
bool sw92_production_cell_evaluator_header();

namespace {

namespace th = mpmc::thermodynamics;
namespace fl = mpmc::flash;
namespace flow = mpmc::flow;
namespace mesh = mpmc::mesh;
namespace dp = mpmc::discretization_petsc;
namespace fdp = mpmc::flow_discretization_petsc;
namespace st = sw92_test;

void require_sw92_petsc(
    bool condition,
    std::string_view message) {
    if (!condition) {
        throw std::runtime_error(
            std::string{message});
    }
}

void near_sw92_petsc(
    double actual,
    double expected,
    double relative = 5.0e-9,
    double absolute = 5.0e-11) {
    const double scale =
        std::max(
            {1.0,
             std::abs(actual),
             std::abs(expected)});
    require_sw92_petsc(
        std::isfinite(actual) &&
            std::isfinite(expected) &&
            std::abs(actual - expected) <=
                absolute +
                    relative * scale,
        "SW92 production PETSc numeric mismatch");
}

th::Provenance nist_source(
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

st::Prepared sourced_input() {
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
                nist_source(
                    component.id +
                    std::string{
                        " molecular weight"}),
                "kg/mol",
                "identity");
    }
    return prepared;
}

th::Sw92ParameterSet sourced_parameters() {
    auto prepared = sourced_input();
    return th::Sw92ParameterSet::create(
        prepared.catalog,
        prepared.order,
        prepared.input);
}

mesh::PartitionSnapshot one_cell_partition() {
    mesh::Topology topology{
        {
            {},
            {},
            {},
            {mesh::GlobalEntityId{UINT64_C(10)}}},
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
empty_schedule() {
    return {
        mesh::PartitionRank{0U},
        1U,
        {},
        {}};
}

dp::PetscMpiAijSymbolicPreallocation3D
one_cell_bridge() {
    return {
        mesh::PartitionRank{0U},
        1U,
        0,
        1,
        1,
        {mesh::LocalIndex{0U}},
        {mesh::GlobalEntityId{UINT64_C(10)}},
        {0},
        {1},
        {0},
        {0}};
}

dp::OwnedCellStructuralColumnPatternSnapshot3D
one_cell_pattern() {
    return {
        mesh::PartitionRank{0U},
        1U,
        1U,
        0,
        1,
        1,
        {{
            mesh::LocalIndex{0U},
            mesh::GlobalEntityId{UINT64_C(10)},
            0,
            0U,
            1U,
            0U,
            0U}},
        {0},
        {}};
}

PetscErrorCode zero_rock_storage(
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
                "stationary zero-rock thermal storage",
                "sw92-production-short-step",
                "v1"}});
    *status =
        fdp::NaturalVariableSnesEvaluationStatus3D::
            success;
    return PETSC_SUCCESS;
}

void set_state(
    Vec state,
    const fdp::VariableCardinalityNaturalVariableNumbering3D&
        numbering,
    std::span<const double> q) {
    for (std::size_t slot = 0U;
         slot < q.size();
         ++slot) {
        const PetscInt row =
            numbering.petsc_global_scalar(
                mesh::LocalIndex{0U},
                slot);
        const PetscScalar value =
            static_cast<PetscScalar>(
                q[slot]);
        require_sw92_petsc(
            VecSetValue(
                state,
                row,
                value,
                INSERT_VALUES) ==
                PETSC_SUCCESS,
            "failed to set SW92 initial state");
    }
    require_sw92_petsc(
        VecAssemblyBegin(state) ==
                PETSC_SUCCESS &&
            VecAssemblyEnd(state) ==
                PETSC_SUCCESS,
        "failed to assemble SW92 initial state");
}

std::vector<double> residual_values(
    fdp::NaturalVariableSnesEvaluator3D evaluator,
    Vec state,
    std::size_t q) {
    Vec residual = nullptr;
    require_sw92_petsc(
        VecDuplicate(state, &residual) ==
            PETSC_SUCCESS,
        "failed to allocate SW92 residual");
    const auto cleanup =
        [&]() {
            if (residual != nullptr) {
                (void)VecDestroy(&residual);
            }
        };
    require_sw92_petsc(
        VecSet(residual, 0.0) ==
            PETSC_SUCCESS,
        "failed to zero SW92 residual");
    fdp::NaturalVariableSnesEvaluationStatus3D
        status =
            fdp::NaturalVariableSnesEvaluationStatus3D::
                success;
    const PetscErrorCode error =
        evaluator.function(
            state,
            residual,
            evaluator.user_context,
            &status);
    if (error != PETSC_SUCCESS ||
        status !=
            fdp::NaturalVariableSnesEvaluationStatus3D::
                success ||
        VecAssemblyBegin(residual) !=
            PETSC_SUCCESS ||
        VecAssemblyEnd(residual) !=
            PETSC_SUCCESS) {
        cleanup();
        throw std::runtime_error(
            "failed to evaluate SW92 residual");
    }
    const PetscScalar* data = nullptr;
    require_sw92_petsc(
        VecGetArrayRead(
            residual,
            &data) ==
            PETSC_SUCCESS,
        "failed to read SW92 residual");
    std::vector<double> values(q);
    for (std::size_t row = 0U;
         row < q;
         ++row) {
        values[row] =
            static_cast<double>(
                PetscRealPart(data[row]));
    }
    require_sw92_petsc(
        VecRestoreArrayRead(
            residual,
            &data) ==
            PETSC_SUCCESS,
        "failed to restore SW92 residual");
    cleanup();
    return values;
}

void sw92_stationary_short_step() {
    auto parameters =
        sourced_parameters();
    auto model =
        th::Sw92Phase<double>::
            from_parameters(parameters);

    const auto source =
        fl::solve_sw92_profile_c_pt_phase_set(
            3.0e6,
            550.0,
            std::vector<double>{0.995, 0.005},
            model,
            0.0);
    require_sw92_petsc(
        source.solution.status ==
                fl::PtPhaseSetStatus::accepted &&
            source.accepted_phase_set_published() &&
            source.solution.accepted_phase_count() ==
                1U,
        "SW92 sourced 550 K fixture is not authoritative single phase");

    auto materialized =
        fdp::
            materialize_sw92_co2_water_profile_c_frozen_cell_3d(
                source,
                model);
    const auto* layout =
        std::get_if<
            flow::NaturalVariableLayout1P>(
                &materialized.frozen_layout);
    require_sw92_petsc(
        layout != nullptr &&
            materialized.natural_variables.size() ==
                layout->unknown_count() &&
            materialized.component_ids.size() ==
                2U,
        "SW92 Profile-C materializer did not produce a 1P natural-variable cell");

    using Closure =
        fdp::Sw92Co2WaterSelectedPhasePropertyClosure3D;
    fdp::Sw92SinglePhaseProductionCellEvaluatorContext3D<
        Closure>
        evaluator_context{
            &materialized.property_closure,
            1.0,
            {&zero_rock_storage, nullptr},
            {}};

    std::optional<
        fdp::SinglePhaseCurrentCellLinearization3D>
        initial_current;
    fdp::NaturalVariableSnesEvaluationStatus3D
        evaluation_status =
            fdp::NaturalVariableSnesEvaluationStatus3D::
                success;
    require_sw92_petsc(
        fdp::
            evaluate_sw92_single_phase_production_cell_3d<
                Closure>(
                mesh::LocalIndex{0U},
                mesh::GlobalEntityId{
                    UINT64_C(10)},
                materialized.natural_variables,
                *layout,
                materialized.component_ids,
                &evaluator_context,
                &initial_current,
                &evaluation_status) ==
                PETSC_SUCCESS &&
            evaluation_status ==
                fdp::NaturalVariableSnesEvaluationStatus3D::
                    success &&
            initial_current.has_value(),
        "SW92 production cell evaluator did not materialize the initial state");

    constexpr double porosity = 0.20;
    fdp::SinglePhaseSnesCellInput3D
        cell;
    cell.cell = mesh::LocalIndex{0U};
    cell.cell_global =
        mesh::GlobalEntityId{UINT64_C(10)};
    cell.bulk_volume_m3 = 1.0;
    cell.porosity = porosity;
    cell.frozen_layout = *layout;
    cell.component_ids =
        materialized.component_ids;
    cell.previous_component_accumulation =
        flow::
            build_single_phase_component_accumulation(
                initial_current->state,
                porosity);
    cell.previous_energy_accumulation =
        flow::
            build_single_phase_energy_accumulation_snapshot(
                initial_current->state,
                porosity,
                initial_current->transport,
                initial_current->caloric,
                initial_current->rock);

    auto partition =
        one_cell_partition();
    auto schedule =
        empty_schedule();
    auto cell_bridge =
        one_cell_bridge();
    auto cell_pattern =
        one_cell_pattern();

    std::optional<
        fdp::VariableCardinalityNaturalVariableNumbering3D>
        numbering;
    const std::array<std::size_t, 1>
        phase_counts{1U};
    require_sw92_petsc(
        fdp::
            make_variable_cardinality_natural_variable_numbering_3d(
                PETSC_COMM_SELF,
                partition,
                2U,
                phase_counts,
                &numbering) ==
                PETSC_SUCCESS &&
            numbering.has_value(),
        "failed SW92 variable-cardinality numbering");

    fdp::MixedCardinalityPhysicalCellEvaluatorBindings3D
        bindings;
    bindings.single_phase = {
        &fdp::
            evaluate_sw92_single_phase_production_cell_3d<
                Closure>,
        &evaluator_context};

    std::vector<
        fdp::MixedCardinalityPhysicalSnesCellInput3D>
        cells;
    cells.emplace_back(
        std::move(cell));
    std::vector<
        flow::FrozenActivePhaseIdentityMap>
        identities;
    identities.emplace_back(
        std::vector<
            flow::FrozenPhysicalPhaseIdentity>{
            {
                "SW92/Profile-C/zero-salinity-CO2-H2O",
                "authoritative-single-phase"}});

    std::optional<
        fdp::MixedCardinalityPhysicalSnesAssemblyContext3D>
        assembly;
    require_sw92_petsc(
        fdp::
            MixedCardinalityPhysicalSnesAssemblyContext3D::
                create(
                    PETSC_COMM_SELF,
                    schedule,
                    partition,
                    *numbering,
                    cell_bridge,
                    cell_pattern,
                    1.0,
                    std::move(cells),
                    std::move(identities),
                    {},
                    bindings,
                    {},
                    &assembly) ==
                PETSC_SUCCESS &&
            assembly.has_value(),
        "failed to create SW92 mixed-cardinality physical assembly");

    Vec initial_state = nullptr;
    Mat jacobian = nullptr;
    require_sw92_petsc(
        fdp::
            create_variable_cardinality_natural_variable_vec_3d(
                PETSC_COMM_SELF,
                *numbering,
                &initial_state) ==
                PETSC_SUCCESS,
        "failed to allocate SW92 natural-variable Vec");
    set_state(
        initial_state,
        *numbering,
        materialized.natural_variables);
    require_sw92_petsc(
        assembly->create_jacobian_structure(
            &jacobian) ==
                PETSC_SUCCESS,
        "failed to create SW92 Jacobian structure");

    const auto evaluator =
        assembly->snes_evaluator();
    const std::size_t q =
        layout->unknown_count();
    const auto residual =
        residual_values(
            evaluator,
            initial_state,
            q);
    for (double value : residual) {
        near_sw92_petsc(
            value,
            0.0,
            0.0,
            2.0e-8);
    }

    require_sw92_petsc(
        MatZeroEntries(jacobian) ==
            PETSC_SUCCESS,
        "failed to zero SW92 Jacobian");
    evaluation_status =
        fdp::NaturalVariableSnesEvaluationStatus3D::
            success;
    require_sw92_petsc(
        evaluator.jacobian(
            initial_state,
            jacobian,
            evaluator.user_context,
            &evaluation_status) ==
                PETSC_SUCCESS &&
            evaluation_status ==
                fdp::NaturalVariableSnesEvaluationStatus3D::
                    success &&
            MatAssemblyBegin(
                jacobian,
                MAT_FINAL_ASSEMBLY) ==
                PETSC_SUCCESS &&
            MatAssemblyEnd(
                jacobian,
                MAT_FINAL_ASSEMBLY) ==
                PETSC_SUCCESS,
        "failed SW92 analytic Jacobian assembly");

    const std::array<double, 3>
        steps{10.0, 1.0e-4, 1.0e-6};
    for (std::size_t column = 0U;
         column < q;
         ++column) {
        Vec plus = nullptr;
        Vec minus = nullptr;
        require_sw92_petsc(
            VecDuplicate(initial_state, &plus) ==
                    PETSC_SUCCESS &&
                VecDuplicate(initial_state, &minus) ==
                    PETSC_SUCCESS &&
                VecCopy(initial_state, plus) ==
                    PETSC_SUCCESS &&
                VecCopy(initial_state, minus) ==
                    PETSC_SUCCESS,
            "failed to allocate SW92 finite-perturbation states");
        const PetscInt global_column =
            numbering->petsc_global_scalar(
                mesh::LocalIndex{0U},
                column);
        require_sw92_petsc(
            VecSetValue(
                plus,
                global_column,
                materialized.natural_variables[column] +
                    steps[column],
                INSERT_VALUES) ==
                    PETSC_SUCCESS &&
                VecSetValue(
                    minus,
                    global_column,
                    materialized.natural_variables[column] -
                        steps[column],
                    INSERT_VALUES) ==
                    PETSC_SUCCESS &&
                VecAssemblyBegin(plus) ==
                    PETSC_SUCCESS &&
                VecAssemblyEnd(plus) ==
                    PETSC_SUCCESS &&
                VecAssemblyBegin(minus) ==
                    PETSC_SUCCESS &&
                VecAssemblyEnd(minus) ==
                    PETSC_SUCCESS,
            "failed to assemble SW92 finite-perturbation states");

        const auto rp =
            residual_values(
                evaluator,
                plus,
                q);
        const auto rm =
            residual_values(
                evaluator,
                minus,
                q);
        for (std::size_t row = 0U;
             row < q;
             ++row) {
            const PetscInt global_row =
                numbering->petsc_global_scalar(
                    mesh::LocalIndex{0U},
                    row);
            PetscScalar analytic = 0.0;
            require_sw92_petsc(
                MatGetValues(
                    jacobian,
                    1,
                    &global_row,
                    1,
                    &global_column,
                    &analytic) ==
                    PETSC_SUCCESS,
                "failed to read SW92 analytic Jacobian");
            const double finite =
                (rp[row] - rm[row]) /
                (2.0 * steps[column]);
            near_sw92_petsc(
                static_cast<double>(
                    PetscRealPart(analytic)),
                finite,
                8.0e-4,
                2.0e-7);
        }
        (void)VecDestroy(&plus);
        (void)VecDestroy(&minus);
    }

    Vec solution = nullptr;
    std::optional<
        fdp::
            VariableCardinalityNaturalVariableSnesSolveReport3D>
        report;
    require_sw92_petsc(
        fdp::
            solve_variable_cardinality_natural_variable_snes_3d(
                PETSC_COMM_SELF,
                *numbering,
                initial_state,
                jacobian,
                evaluator,
                &solution,
                &report) ==
                PETSC_SUCCESS &&
            solution != nullptr &&
            report.has_value() &&
            static_cast<int>(
                report->converged_reason) >
                0,
        "SW92 stationary physical short-step did not converge");

    const PetscScalar* solved = nullptr;
    require_sw92_petsc(
        VecGetArrayRead(
            solution,
            &solved) ==
            PETSC_SUCCESS,
        "failed to read SW92 solved state");
    for (std::size_t slot = 0U;
         slot < q;
         ++slot) {
        near_sw92_petsc(
            static_cast<double>(
                PetscRealPart(
                    solved[slot])),
            materialized.natural_variables[slot],
            2.0e-10,
            2.0e-9);
    }
    require_sw92_petsc(
        VecRestoreArrayRead(
            solution,
            &solved) ==
            PETSC_SUCCESS,
        "failed to restore SW92 solved state");

    bool history_matches = false;
    require_sw92_petsc(
        assembly->accepted_history_matches_state(
            solution,
            &history_matches) ==
                PETSC_SUCCESS &&
            history_matches,
        "SW92 stationary short-step changed component or energy inventory");

    (void)VecDestroy(&solution);
    (void)MatDestroy(&jacobian);
    (void)VecDestroy(&initial_state);
}

} // namespace

void sw92_production_fully_implicit_short_step_test() {
    require_sw92_petsc(
        selected_phase_production_cell_evaluator_header(),
        "generic selected-phase production evaluator header probe failed");
    require_sw92_petsc(
        sw92_production_cell_evaluator_header(),
        "SW92 production evaluator header probe failed");
    sw92_stationary_short_step();
}
