#include <mpmc/flow_discretization_petsc/physical_timestep_driver.hpp>
#include <mpmc/flow_discretization_petsc/post_snes_phase_transition_controller.hpp>
#include <mpmc/flow_discretization_petsc/post_snes_phase_transition_handoff_scanner.hpp>
#include <mpmc/well_discretization_petsc/fixed_bhp_well_source_evaluator.hpp>
#include <mpmc/well_discretization_petsc/fixed_total_molar_rate_control.hpp>
#include <mpmc/well_discretization_petsc/fixed_total_molar_rate_timestep_driver.hpp>
#include <mpmc/thermodynamics/pr_parameters.hpp>

#include <petscmat.h>
#include <petscvec.h>

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
#include <type_traits>
#include <utility>
#include <variant>
#include <vector>

namespace {

#include "mixed_cardinality/fixture.inc"
#include "mixed_cardinality/thermodynamic_adapter.inc"
#include "mixed_cardinality/phase_transition.inc"
#include "mixed_cardinality/well_fixture.inc"
#include "mixed_cardinality/fixed_bhp.inc"
#include "mixed_cardinality/rate_control.inc"
#include "mixed_cardinality/fixed_bhp_transition.inc"
} // namespace

void mixed_cardinality_physical_snes_assembly_test() {
    int rank = -1;
    int size = -1;
    if (MPI_Comm_rank(
            PETSC_COMM_WORLD,
            &rank) != MPI_SUCCESS ||
        MPI_Comm_size(
            PETSC_COMM_WORLD,
            &size) != MPI_SUCCESS) {
        throw std::runtime_error(
            "failed to query MPI rank/size for mixed physical dispatcher test");
    }
    require_collective(
        size == 2,
        "mixed physical dispatcher regression requires two ranks");

    const auto partition =
        make_partition(rank);
    const std::vector<std::size_t>
        phase_counts{
            1U, 1U, 2U, 2U, 3U, 3U};
    std::optional<
        fdp::
            VariableCardinalityNaturalVariableNumbering3D>
        numbering;
    PetscErrorCode error =
        fdp::
            make_variable_cardinality_natural_variable_numbering_3d(
                PETSC_COMM_WORLD,
                partition,
                3U,
                phase_counts,
                &numbering);
    require_collective(
        error == PETSC_SUCCESS &&
            numbering.has_value() &&
            numbering->petsc_global_scalar_count() ==
                42,
        "failed to build 1P/2P/3P mixed physical numbering");

    const auto bridge =
        make_cell_bridge(rank);
    const auto pattern =
        make_cell_pattern(rank);
    const auto schedule =
        make_schedule(
            rank,
            false);

    DispatchAudit audit;

    auto fixed_bhp_well_context =
        wdp::
            FixedBhpThreePhasePeacemanWellSourceEvaluatorContext3D::
                create(
                    mesh::GlobalEntityId{
                        UINT64_C(60)},
                    well::make_peaceman_well_index_3d(
                        {10.0, 10.0, 5.0},
                        {
                            1.0e-12,
                            1.0e-12,
                            1.0e-12},
                        well::
                            AxisAlignedWellDirection3D::z,
                        0.10,
                        0.0),
                    25.0,
                    wd::InjectionPhaseSpecificEnthalpy3P{
                        "fixture/fixed-bhp-injection-enthalpy/v1",
                        {1000.0, 2000.0, 3000.0}},
                    "fixture/fixed-bhp-cell60-source/v1");

    // The first bridge intentionally supports only a frozen 3P target. A 1P
    // target is rejected explicitly rather than padded with fictitious phases.
    {
        auto unsupported_context =
            wdp::
                FixedBhpThreePhasePeacemanWellSourceEvaluatorContext3D::
                    create(
                        mesh::GlobalEntityId{
                            UINT64_C(20)},
                        fixed_bhp_well_context
                            .connection(),
                        5.0,
                        wd::InjectionPhaseSpecificEnthalpy3P{
                            "fixture/fixed-bhp-unsupported/v1",
                            {1000.0, 2000.0, 3000.0}},
                        "fixture/fixed-bhp-unsupported-source/v1");
        auto current_1p =
            evaluate_target(
                UINT64_C(20),
                &audit);
        std::optional<
            fd::CellSourceLinearization3D>
            unsupported_source;
        fdp::NaturalVariableSnesEvaluationStatus3D
            unsupported_status =
                fdp::
                    NaturalVariableSnesEvaluationStatus3D::
                        success;
        const auto q1 =
            target_1p();
        const PetscErrorCode unsupported_error =
            wdp::
                evaluate_fixed_bhp_three_phase_peaceman_well_source_3d(
                    mesh::LocalIndex{1U},
                    mesh::GlobalEntityId{
                        UINT64_C(20)},
                    q1,
                    current_1p,
                    &unsupported_context,
                    &unsupported_source,
                    &unsupported_status);
        require_collective(
            unsupported_error ==
                    PETSC_ERR_SUP &&
                !unsupported_source.has_value(),
            "fixed-BHP well bridge fabricated inactive phases for a non-3P target");
    }

    auto cell_inputs =
        make_cell_inputs(
            rank,
            &audit);
    auto face_inputs =
        make_face_inputs(
            rank,
            false);

    std::optional<
        fdp::
            MixedCardinalityPhysicalSnesAssemblyContext3D>
        context;
    error =
        fdp::
            MixedCardinalityPhysicalSnesAssemblyContext3D::
                create(
                    PETSC_COMM_WORLD,
                    schedule,
                    partition,
                    *numbering,
                    bridge,
                    pattern,
                    1.0,
                    std::move(cell_inputs),
                    make_phase_identity_maps(),
                    std::move(face_inputs),
                    {
                        {&evaluate_1p, &audit},
                        {&evaluate_2p, &audit},
                        {&evaluate_3p, &audit}},
                    {},
                    {
                        &evaluate_explicit_cell_source,
                        &audit},
                    &context);
    require_collective(
        error == PETSC_SUCCESS &&
            context.has_value(),
        "failed to create mixed physical production dispatcher");

    Mat jacobian = nullptr;
    error =
        context->create_jacobian_structure(
            &jacobian);
    require_collective(
        error == PETSC_SUCCESS &&
            jacobian != nullptr,
        "failed to create mixed physical ragged MPIAIJ structure");

    Vec state = nullptr;
    error =
        fdp::
            create_variable_cardinality_natural_variable_vec_3d(
                PETSC_COMM_WORLD,
                *numbering,
                &state);
    require_collective(
        error == PETSC_SUCCESS &&
            state != nullptr,
        "failed to create mixed physical state vector");
    insert_target(
        state,
        *numbering);

    Vec residual = nullptr;
    error =
        VecDuplicate(
            state,
            &residual);
    if (error == PETSC_SUCCESS) {
        error =
            VecSet(
                residual,
                PetscScalar{0.0});
    }
    auto evaluator =
        context->snes_evaluator();
    fdp::NaturalVariableSnesEvaluationStatus3D
        status =
            fdp::
                NaturalVariableSnesEvaluationStatus3D::
                    success;
    if (error == PETSC_SUCCESS) {
        error =
            evaluator.function(
                state,
                residual,
                evaluator.user_context,
                &status);
    }
    if (error == PETSC_SUCCESS) {
        error =
            VecAssemblyBegin(
                residual);
    }
    if (error == PETSC_SUCCESS) {
        error =
            VecAssemblyEnd(
                residual);
    }
    PetscReal residual_norm = -1.0;
    if (error == PETSC_SUCCESS) {
        error =
            VecNorm(
                residual,
                NORM_2,
                &residual_norm);
    }
    require_collective(
        error == PETSC_SUCCESS &&
            status ==
                fdp::
                    NaturalVariableSnesEvaluationStatus3D::
                        success &&
            residual_norm < 1.0e-12,
        "target mixed physical residual is not zero");

    error =
        MatZeroEntries(
            jacobian);
    status =
        fdp::
            NaturalVariableSnesEvaluationStatus3D::
                success;
    if (error == PETSC_SUCCESS) {
        error =
            evaluator.jacobian(
                state,
                jacobian,
                evaluator.user_context,
                &status);
    }
    if (error == PETSC_SUCCESS) {
        error =
            MatAssemblyBegin(
                jacobian,
                MAT_FINAL_ASSEMBLY);
    }
    if (error == PETSC_SUCCESS) {
        error =
            MatAssemblyEnd(
                jacobian,
                MAT_FINAL_ASSEMBLY);
    }
    require_collective(
        error == PETSC_SUCCESS &&
            status ==
                fdp::
                    NaturalVariableSnesEvaluationStatus3D::
                        success,
        "failed to assemble mixed physical Jacobian");

    bool local_spatial_blocks_ok = true;
    if (rank == 0) {
        const std::array<
            std::pair<PetscInt, PetscInt>,
            3>
            probes{{
                {0, 21},
                {4, 25},
                {11, 32}}};
        for (const auto& [row, column] :
             probes) {
            PetscScalar value = 0.0;
            const auto get_error =
                MatGetValues(
                    jacobian,
                    1,
                    &row,
                    1,
                    &column,
                    &value);
            local_spatial_blocks_ok =
                local_spatial_blocks_ok &&
                get_error ==
                    PETSC_SUCCESS &&
                std::isfinite(
                    static_cast<double>(
                        PetscRealPart(
                            value))) &&
                std::abs(
                    static_cast<double>(
                        PetscRealPart(
                            value))) >
                    0.0;
        }
    }
    require_collective(
        local_spatial_blocks_ok,
        "1P/2P/3P real TPFA off-diagonal blocks were not assembled");

    require_collective(
        audit.single_calls > 0U &&
            audit.two_calls > 0U &&
            audit.three_calls > 0U,
        "mixed physical dispatcher did not exercise all three cell closures");

    PetscScalar baseline_component_d = 0.0;
    PetscScalar baseline_energy_d = 0.0;
    if (rank == 1) {
        const auto& source_record =
            numbering->cell(
                mesh::LocalIndex{1U});
        const PetscInt row0 =
            source_record.petsc_global_scalar_start;
        const PetscInt energy_row =
            row0 + 3;
        const PetscInt column = row0;
        require_collective(
            MatGetValues(
                jacobian,
                1,
                &row0,
                1,
                &column,
                &baseline_component_d) ==
                PETSC_SUCCESS &&
            MatGetValues(
                jacobian,
                1,
                &energy_row,
                1,
                &column,
                &baseline_energy_d) ==
                PETSC_SUCCESS,
            "failed to capture source-disabled diagonal Jacobian");
    } else {
        require_collective(
            true,
            "source-disabled diagonal Jacobian");
    }

    PetscScalar baseline_well_component_d = 0.0;
    PetscScalar baseline_well_energy_d = 0.0;
    if (rank == 1) {
        const auto& well_record =
            numbering->cell(
                mesh::LocalIndex{5U});
        const PetscInt row0 =
            well_record.petsc_global_scalar_start;
        const PetscInt energy_row =
            row0 + 3;
        const PetscInt column =
            row0;
        require_collective(
            MatGetValues(
                jacobian,
                1,
                &row0,
                1,
                &column,
                &baseline_well_component_d) ==
                PETSC_SUCCESS &&
            MatGetValues(
                jacobian,
                1,
                &energy_row,
                1,
                &column,
                &baseline_well_energy_d) ==
                PETSC_SUCCESS,
            "failed to capture fixed-BHP well target baseline Jacobian");
    } else {
        require_collective(
            true,
            "fixed-BHP well target baseline Jacobian");
    }

    audit.source_enabled = true;
    error = VecSet(
        residual,
        PetscScalar{0.0});
    status =
        fdp::NaturalVariableSnesEvaluationStatus3D::
            success;
    if (error == PETSC_SUCCESS) {
        error =
            evaluator.function(
                state,
                residual,
                evaluator.user_context,
                &status);
    }
    if (error == PETSC_SUCCESS) {
        error = VecAssemblyBegin(residual);
    }
    if (error == PETSC_SUCCESS) {
        error = VecAssemblyEnd(residual);
    }
    if (error == PETSC_SUCCESS) {
        error = MatZeroEntries(jacobian);
    }
    if (error == PETSC_SUCCESS) {
        error =
            evaluator.jacobian(
                state,
                jacobian,
                evaluator.user_context,
                &status);
    }
    if (error == PETSC_SUCCESS) {
        error =
            MatAssemblyBegin(
                jacobian,
                MAT_FINAL_ASSEMBLY);
    }
    if (error == PETSC_SUCCESS) {
        error =
            MatAssemblyEnd(
                jacobian,
                MAT_FINAL_ASSEMBLY);
    }

    bool local_source_ok = true;
    if (rank == 1) {
        const auto& source_record =
            numbering->cell(
                mesh::LocalIndex{1U});
        const PetscInt row0 =
            source_record.petsc_global_scalar_start;
        const PetscInt energy_row =
            row0 + 3;
        const PetscInt column = row0;
        std::array<PetscInt, 4> rows{
            row0,
            row0 + 1,
            row0 + 2,
            energy_row};
        std::array<PetscScalar, 4> values{};
        local_source_ok =
            VecGetValues(
                residual,
                static_cast<PetscInt>(rows.size()),
                rows.data(),
                values.data()) ==
            PETSC_SUCCESS;
        const std::array<double, 4> expected{
            -2.0 / 3.0,
            1.0 / 3.0,
            -0.5 / 3.0,
            -20.0 / 3.0};
        for (std::size_t i = 0U;
             i < expected.size();
             ++i) {
            local_source_ok =
                local_source_ok &&
                std::abs(
                    static_cast<double>(
                        PetscRealPart(values[i])) -
                    expected[i]) <
                    1.0e-11;
        }

        PetscScalar sourced_component_d = 0.0;
        PetscScalar sourced_energy_d = 0.0;
        local_source_ok =
            local_source_ok &&
            MatGetValues(
                jacobian,
                1,
                &row0,
                1,
                &column,
                &sourced_component_d) ==
                PETSC_SUCCESS &&
            MatGetValues(
                jacobian,
                1,
                &energy_row,
                1,
                &column,
                &sourced_energy_d) ==
                PETSC_SUCCESS &&
            std::abs(
                static_cast<double>(
                    PetscRealPart(
                        sourced_component_d -
                        baseline_component_d)) +
                0.1 / 3.0) <
                1.0e-11 &&
            std::abs(
                static_cast<double>(
                    PetscRealPart(
                        sourced_energy_d -
                        baseline_energy_d)) +
                0.2 / 3.0) <
                1.0e-11;
    }
    require_collective(
        error == PETSC_SUCCESS &&
            status ==
                fdp::NaturalVariableSnesEvaluationStatus3D::
                    success &&
            local_source_ok &&
            (rank == 0
                 ? audit.source_cell20_calls == 0U
                 : audit.source_cell20_calls > 0U),
        "explicit cell source did not map to owner component/energy residual and diagonal Jacobian");

    // The actual fixed-BHP well bridge now enters the same owner-only source
    // callback and therefore the fully implicit mixed-cardinality residual and
    // diagonal Jacobian. The connection is attached to stable cell60, owned by
    // rank1 and ghosted on rank0.
    FixedBhpWellSourceAudit
        fixed_bhp_source_audit{
            &fixed_bhp_well_context};
    auto well_cells =
        make_cell_inputs(
            rank,
            &audit);
    auto well_faces =
        make_face_inputs(
            rank,
            false);
    std::optional<
        fdp::
            MixedCardinalityPhysicalSnesAssemblyContext3D>
        well_assembly_context;
    error =
        fdp::
            MixedCardinalityPhysicalSnesAssemblyContext3D::
                create(
                    PETSC_COMM_WORLD,
                    schedule,
                    partition,
                    *numbering,
                    bridge,
                    pattern,
                    1.0,
                    std::move(well_cells),
                    make_phase_identity_maps(),
                    std::move(well_faces),
                    {
                        {&evaluate_1p, &audit},
                        {&evaluate_2p, &audit},
                        {&evaluate_3p, &audit}},
                    {},
                    {
                        &evaluate_audited_fixed_bhp_well_source,
                        &fixed_bhp_source_audit},
                    &well_assembly_context);
    require_collective(
        error == PETSC_SUCCESS &&
            well_assembly_context.has_value(),
        "failed to create fixed-BHP well mixed physical assembly context");

    auto direct_current =
        evaluate_target(
            UINT64_C(60),
            &audit);
    const auto& direct_three_phase =
        std::get<
            fdp::
                FixedThreePhaseCurrentCellLinearization3D>(
                    direct_current);
    const auto expected_well_source =
        wdp::
            build_fixed_bhp_three_phase_peaceman_well_source_3d(
                fixed_bhp_well_context,
                direct_three_phase);
    const auto expected_well_normalized =
        fd::normalize_cell_source_by_bulk_volume(
            expected_well_source.cell_source,
            7.0);

    Mat well_jacobian = nullptr;
    Vec well_residual = nullptr;
    error =
        well_assembly_context
            ->create_jacobian_structure(
                &well_jacobian);
    if (error == PETSC_SUCCESS) {
        error =
            VecDuplicate(
                state,
                &well_residual);
    }
    if (error == PETSC_SUCCESS) {
        error =
            VecSet(
                well_residual,
                PetscScalar{0.0});
    }
    auto well_evaluator =
        well_assembly_context
            ->snes_evaluator();
    status =
        fdp::
            NaturalVariableSnesEvaluationStatus3D::
                success;
    if (error == PETSC_SUCCESS) {
        error =
            well_evaluator.function(
                state,
                well_residual,
                well_evaluator.user_context,
                &status);
    }
    if (error == PETSC_SUCCESS) {
        error =
            VecAssemblyBegin(
                well_residual);
    }
    if (error == PETSC_SUCCESS) {
        error =
            VecAssemblyEnd(
                well_residual);
    }
    if (error == PETSC_SUCCESS) {
        error =
            MatZeroEntries(
                well_jacobian);
    }
    if (error == PETSC_SUCCESS) {
        error =
            well_evaluator.jacobian(
                state,
                well_jacobian,
                well_evaluator.user_context,
                &status);
    }
    if (error == PETSC_SUCCESS) {
        error =
            MatAssemblyBegin(
                well_jacobian,
                MAT_FINAL_ASSEMBLY);
    }
    if (error == PETSC_SUCCESS) {
        error =
            MatAssemblyEnd(
                well_jacobian,
                MAT_FINAL_ASSEMBLY);
    }

    bool local_well_source_ok =
        error == PETSC_SUCCESS &&
        status ==
            fdp::
                NaturalVariableSnesEvaluationStatus3D::
                    success;
    if (rank == 1 &&
        local_well_source_ok) {
        const auto& well_record =
            numbering->cell(
                mesh::LocalIndex{5U});
        const PetscInt row0 =
            well_record.petsc_global_scalar_start;
        const PetscInt energy_row =
            row0 + 3;
        const PetscInt column =
            row0;
        std::array<PetscInt, 4>
            rows{
                row0,
                row0 + 1,
                row0 + 2,
                energy_row};
        std::array<PetscScalar, 4>
            values{};
        local_well_source_ok =
            VecGetValues(
                well_residual,
                static_cast<PetscInt>(
                    rows.size()),
                rows.data(),
                values.data()) ==
            PETSC_SUCCESS;

        for (std::size_t component = 0U;
             component < 3U;
             ++component) {
            const double actual =
                static_cast<double>(
                    PetscRealPart(
                        values[component]));
            const double expected =
                expected_well_normalized
                    .component_residual_mol_per_bulk_m3_s[
                        component];
            local_well_source_ok =
                local_well_source_ok &&
                std::abs(
                    actual -
                    expected) <=
                    1.0e-10 *
                        std::max(
                            1.0e-30,
                            std::abs(expected)) +
                    1.0e-20;
        }
        {
            const double actual =
                static_cast<double>(
                    PetscRealPart(
                        values[3]));
            const double expected =
                expected_well_normalized
                    .energy_residual_w_per_bulk_m3;
            local_well_source_ok =
                local_well_source_ok &&
                std::abs(
                    actual -
                    expected) <=
                    1.0e-10 *
                        std::max(
                            1.0e-30,
                            std::abs(expected)) +
                    1.0e-20;
        }

        PetscScalar well_component_d = 0.0;
        PetscScalar well_energy_d = 0.0;
        local_well_source_ok =
            local_well_source_ok &&
            MatGetValues(
                well_jacobian,
                1,
                &row0,
                1,
                &column,
                &well_component_d) ==
                PETSC_SUCCESS &&
            MatGetValues(
                well_jacobian,
                1,
                &energy_row,
                1,
                &column,
                &well_energy_d) ==
                PETSC_SUCCESS;

        const double component_delta =
            static_cast<double>(
                PetscRealPart(
                    well_component_d -
                    baseline_well_component_d));
        const double energy_delta =
            static_cast<double>(
                PetscRealPart(
                    well_energy_d -
                    baseline_well_energy_d));
        const double expected_component_delta =
            expected_well_normalized
                .d_component_residual(
                    0U,
                    0U);
        const double expected_energy_delta =
            expected_well_normalized
                .d_energy_residual(
                    0U);
        local_well_source_ok =
            local_well_source_ok &&
            std::abs(
                component_delta -
                expected_component_delta) <=
                1.0e-10 *
                    std::max(
                        1.0e-30,
                        std::abs(
                            expected_component_delta)) +
                1.0e-20 &&
            std::abs(
                energy_delta -
                expected_energy_delta) <=
                1.0e-10 *
                    std::max(
                        1.0e-30,
                        std::abs(
                            expected_energy_delta)) +
                1.0e-20;
    }

    require_collective(
        local_well_source_ok &&
            (rank == 0
                 ? fixed_bhp_source_audit
                           .target_calls ==
                       0U
                 : fixed_bhp_source_audit
                           .target_calls ==
                       2U),
        "fixed-BHP production well source did not enter owner-only residual/Jacobian correctly");

    require_collective(
        (well_residual == nullptr ||
         VecDestroy(
             &well_residual) ==
             PETSC_SUCCESS) &&
            (well_jacobian == nullptr ||
             MatDestroy(
                 &well_jacobian) ==
                 PETSC_SUCCESS),
        "fixed-BHP well source assembly cleanup failed");

    audit.source_enabled = false;


    // Exercise the same production fixed-BHP source/driver/conservation path
    // on one frozen target cell of each supported phase cardinality. The mixed
    // reservoir around the completion is unchanged; only the well target chart
    // and its exact active-phase injection-enthalpy payload differ.
    run_frozen_fixed_bhp_timestep_case(
        rank,
        UINT64_C(20),
        1U,
        1U,
        5.0,
        schedule,
        partition,
        bridge,
        pattern,
        &audit);
    run_frozen_fixed_bhp_timestep_case(
        rank,
        UINT64_C(40),
        3U,
        2U,
        15.0,
        schedule,
        partition,
        bridge,
        pattern,
        &audit);
    run_frozen_fixed_bhp_timestep_case(
        rank,
        UINT64_C(60),
        5U,
        3U,
        25.0,
        schedule,
        partition,
        bridge,
        pattern,
        &audit);

    run_multi_connection_fixed_bhp_timestep_case(
        rank,
        schedule,
        partition,
        bridge,
        pattern,
        &audit);

    run_fixed_total_molar_rate_control_case(
        rank,
        schedule,
        partition,
        bridge,
        pattern,
        &audit);

    Vec solution = nullptr;
    std::optional<
        fdp::
            VariableCardinalityNaturalVariableSnesSolveReport3D>
        report;
    error =
        fdp::
            solve_variable_cardinality_natural_variable_snes_3d(
                PETSC_COMM_WORLD,
                *numbering,
                state,
                jacobian,
                context->snes_evaluator(),
                &solution,
                &report);
    require_collective(
        error == PETSC_SUCCESS &&
            solution != nullptr &&
            report.has_value() &&
            static_cast<int>(
                report->converged_reason) >
                0 &&
            report->snes_type ==
                SNESNEWTONLS &&
            report->ksp_type ==
                KSPGMRES &&
            report->pc_type ==
                PCASM &&
            local_solution_matches_target(
                *report),
        "mixed physical system did not pass PETSc SNES production solve");

    // A cross-cardinality face must never be guessed from the fixed-cardinality
    // TPFA kernels. Without an explicit bridge evaluator, construction is
    // collectively rejected before nonlinear callbacks are installed.
    const auto cross_schedule =
        make_schedule(
            rank,
            true);
    auto cross_cells =
        make_cell_inputs(
            rank,
            &audit);
    auto cross_faces =
        make_face_inputs(
            rank,
            true);
    std::optional<
        fdp::
            MixedCardinalityPhysicalSnesAssemblyContext3D>
        rejected;
    error =
        fdp::
            MixedCardinalityPhysicalSnesAssemblyContext3D::
                create(
                    PETSC_COMM_WORLD,
                    cross_schedule,
                    partition,
                    *numbering,
                    bridge,
                    pattern,
                    1.0,
                    std::move(cross_cells),
                    make_phase_identity_maps(),
                    std::move(cross_faces),
                    {
                        {&evaluate_1p, &audit},
                        {&evaluate_2p, &audit},
                        {&evaluate_3p, &audit}},
                    {},
                    &rejected);
    require_collective(
        error == PETSC_ERR_SUP &&
            !rejected.has_value(),
        "cross-cardinality TPFA face was not rejected without an explicit bridge");

    // Explicit potential extensions enable the standard physical bridge for
    // all three mixed-cardinality pair classes. The bridge computes actual
    // TPFA component/energy flux from active upwind payloads and never requests
    // absent-side composition or enthalpy.
    check_direct_cross_cardinality_bridge(
        UINT64_C(20),
        UINT64_C(30),
        3U,
        &audit);
    check_direct_cross_cardinality_bridge(
        UINT64_C(20),
        UINT64_C(50),
        4U,
        &audit);
    check_direct_cross_cardinality_bridge(
        UINT64_C(40),
        UINT64_C(50),
        5U,
        &audit);

    // The 1P<->2P face is now admitted through the standard bridge into the
    // real distributed dispatcher. Rank 0 owns the authoritative face while
    // cell20 is owned by rank 1, so this exercises off-process residual and
    // rectangular Jacobian insertion, not just a direct local bridge call.
    auto bridged_cells =
        make_cell_inputs(
            rank,
            &audit);
    auto bridged_faces =
        make_face_inputs(
            rank,
            true);
    const auto pr_parameters =
        thermodynamic_adapter_pr_parameters();
    const auto pr_model =
        th::Pr76Phase<double>::
            from_parameters(
                pr_parameters);
    auto pr_provider =
        make_thermodynamic_adapter_pr_provider(
            pr_model);
    ThermodynamicCoordinateResolverAudit
        coordinate_audit;
    fdp::
        ThermodynamicAbsentPhaseExtensionAdapterBinding3D
        thermodynamic_adapter{
            &resolve_absent_phase_thermodynamic_coordinates,
            &coordinate_audit,
            &fdp::
                evaluate_absent_phase_thermodynamic_provider_3d<
                    flow::
                        Pr76AbsentPhasePotentialExtensionProvider<
                            double>>,
            &pr_provider};
    fdp::CrossCardinalityTpfaBridgeBinding3D
        standard_bridge{
            &fdp::
                evaluate_thermodynamic_absent_phase_extension_3d,
            &thermodynamic_adapter};
    std::optional<
        fdp::
            MixedCardinalityPhysicalSnesAssemblyContext3D>
        bridged_context;
    error =
        fdp::
            MixedCardinalityPhysicalSnesAssemblyContext3D::
                create(
                    PETSC_COMM_WORLD,
                    cross_schedule,
                    partition,
                    *numbering,
                    bridge,
                    pattern,
                    1.0,
                    std::move(bridged_cells),
                    make_phase_identity_maps(),
                    std::move(bridged_faces),
                    {
                        {&evaluate_1p, &audit},
                        {&evaluate_2p, &audit},
                        {&evaluate_3p, &audit}},
                    {
                        &fdp::
                            evaluate_standard_cross_cardinality_tpfa_face_3d,
                        &standard_bridge},
                    &bridged_context);
    require_collective(
        error == PETSC_SUCCESS &&
            bridged_context.has_value(),
        "standard cross-cardinality bridge was not admitted by distributed dispatcher");

    Mat bridged_jacobian = nullptr;
    error =
        bridged_context
            ->create_jacobian_structure(
                &bridged_jacobian);
    require_collective(
        error == PETSC_SUCCESS &&
            bridged_jacobian != nullptr,
        "failed to create bridged cross-cardinality MPIAIJ structure");

    Vec bridged_state = nullptr;
    error =
        fdp::
            create_variable_cardinality_natural_variable_vec_3d(
                PETSC_COMM_WORLD,
                *numbering,
                &bridged_state);
    require_collective(
        error == PETSC_SUCCESS &&
            bridged_state != nullptr,
        "failed to create bridged cross-cardinality state");
    insert_target(
        bridged_state,
        *numbering);

    Vec bridged_residual = nullptr;
    error =
        VecDuplicate(
            bridged_state,
            &bridged_residual);
    if (error == PETSC_SUCCESS) {
        error =
            VecSet(
                bridged_residual,
                PetscScalar{0.0});
    }
    auto bridged_evaluator =
        bridged_context->snes_evaluator();
    status =
        fdp::
            NaturalVariableSnesEvaluationStatus3D::
                success;
    if (error == PETSC_SUCCESS) {
        error =
            bridged_evaluator.function(
                bridged_state,
                bridged_residual,
                bridged_evaluator.user_context,
                &status);
    }
    if (error == PETSC_SUCCESS) {
        error =
            VecAssemblyBegin(
                bridged_residual);
    }
    if (error == PETSC_SUCCESS) {
        error =
            VecAssemblyEnd(
                bridged_residual);
    }
    PetscReal bridged_norm = 0.0;
    if (error == PETSC_SUCCESS) {
        error =
            VecNorm(
                bridged_residual,
                NORM_2,
                &bridged_norm);
    }
    require_collective(
        error == PETSC_SUCCESS &&
            status ==
                fdp::
                    NaturalVariableSnesEvaluationStatus3D::
                        success &&
            bridged_norm > 0.0,
        "distributed cross-cardinality face did not contribute a physical residual");

    error =
        MatZeroEntries(
            bridged_jacobian);
    status =
        fdp::
            NaturalVariableSnesEvaluationStatus3D::
                success;
    if (error == PETSC_SUCCESS) {
        error =
            bridged_evaluator.jacobian(
                bridged_state,
                bridged_jacobian,
                bridged_evaluator.user_context,
                &status);
    }
    if (error == PETSC_SUCCESS) {
        error =
            MatAssemblyBegin(
                bridged_jacobian,
                MAT_FINAL_ASSEMBLY);
    }
    if (error == PETSC_SUCCESS) {
        error =
            MatAssemblyEnd(
                bridged_jacobian,
                MAT_FINAL_ASSEMBLY);
    }
    require_collective(
        error == PETSC_SUCCESS &&
            status ==
                fdp::
                    NaturalVariableSnesEvaluationStatus3D::
                        success,
        "distributed cross-cardinality Jacobian assembly failed");

    std::array<double, 4>
        local_cross_rate{};
    if (rank == 0) {
        const auto& cell30 =
            numbering->cell(
                mesh::LocalIndex{2U});
        for (std::size_t row = 0U;
             row < local_cross_rate.size();
             ++row) {
            const PetscInt index =
                cell30.petsc_global_scalar_start +
                static_cast<PetscInt>(row);
            PetscScalar value = 0.0;
            if (VecGetValues(
                    bridged_residual,
                    1,
                    &index,
                    &value) !=
                PETSC_SUCCESS) {
                throw std::runtime_error(
                    "failed to read cell30 bridged residual");
            }
            local_cross_rate[row] =
                static_cast<double>(
                    PetscRealPart(value)) *
                4.0;
        }
    } else {
        const auto& cell20 =
            numbering->cell(
                mesh::LocalIndex{1U});
        for (std::size_t row = 0U;
             row < local_cross_rate.size();
             ++row) {
            const PetscInt index =
                cell20.petsc_global_scalar_start +
                static_cast<PetscInt>(row);
            PetscScalar value = 0.0;
            if (VecGetValues(
                    bridged_residual,
                    1,
                    &index,
                    &value) !=
                PETSC_SUCCESS) {
                throw std::runtime_error(
                    "failed to read cell20 bridged residual");
            }
            local_cross_rate[row] =
                static_cast<double>(
                    PetscRealPart(value)) *
                3.0;
        }
    }

    std::array<double, 4>
        global_cross_rate{};
    require_collective(
        MPI_Allreduce(
            local_cross_rate.data(),
            global_cross_rate.data(),
            static_cast<int>(
                global_cross_rate.size()),
            MPI_DOUBLE,
            MPI_SUM,
            PETSC_COMM_WORLD) ==
            MPI_SUCCESS,
        "failed to reduce bridged cross-face conservation rows");
    bool conserved = true;
    for (const double value :
         global_cross_rate) {
        conserved =
            conserved &&
            std::abs(value) <
                1.0e-10;
    }
    require_collective(
        conserved,
        "distributed cross-cardinality component/energy face rows are not conservative");

    bool cross_block_nonzero = true;
    {
        const auto& row_cell =
            numbering->cell(
                rank == 0
                    ? mesh::LocalIndex{2U}
                    : mesh::LocalIndex{1U});
        const auto& column_cell =
            numbering->cell(
                rank == 0
                    ? mesh::LocalIndex{1U}
                    : mesh::LocalIndex{2U});
        const PetscInt row =
            row_cell
                .petsc_global_scalar_start;
        const PetscInt column =
            column_cell
                .petsc_global_scalar_start;
        PetscScalar value = 0.0;
        cross_block_nonzero =
            MatGetValues(
                bridged_jacobian,
                1,
                &row,
                1,
                &column,
                &value) ==
                PETSC_SUCCESS &&
            std::isfinite(
                static_cast<double>(
                    PetscRealPart(value))) &&
            std::abs(
                static_cast<double>(
                    PetscRealPart(value))) >
                0.0;
    }
    require_collective(
        cross_block_nonzero,
        "distributed 4x7/7x4 cross-cardinality Jacobian block is missing");
    std::uint64_t global_coordinate_calls = 0U;
    require_collective(
        MPI_Allreduce(
            &coordinate_audit.calls,
            &global_coordinate_calls,
            1,
            MPI_UINT64_T,
            MPI_SUM,
            PETSC_COMM_WORLD) ==
            MPI_SUCCESS,
        "failed to reduce thermodynamic coordinate resolver call count");
    require_collective(
        global_coordinate_calls > 0U,
        "authoritative distributed bridge did not invoke the thermodynamic coordinate resolver");

    require_collective(
        VecDestroy(&bridged_residual) ==
                PETSC_SUCCESS &&
            VecDestroy(&bridged_state) ==
                PETSC_SUCCESS &&
            MatDestroy(&bridged_jacobian) ==
                PETSC_SUCCESS,
        "bridged cross-cardinality fixture cleanup failed");

    // Accepted phase transitions are applied only outside SNES.  Rebuild the
    // complete q-ragged production system after cell30 changes 2P -> 3P.
    // This must recreate global numbering/state, the PetscSF-backed physical
    // context and MPIAIJ structure while preserving BE histories.
    auto rebuild_cells =
        make_outer_rebuild_cells(
            rank,
            &audit);
    auto rebuild_faces =
        make_face_inputs(
            rank,
            true);
    auto outer_provider =
        make_outer_rebuild_pr_provider(
            pr_model);

    // A real fixed-BHP completion now follows stable cell30 through the
    // existing 2P->3P outer restart, resolves its well-side enthalpy by stable
    // phase identity, and then advances one post-rebuild physical timestep.
    run_phase_transition_rebound_fixed_bhp_case(
        rank,
        schedule,
        partition,
        bridge,
        pattern,
        &audit,
        &outer_provider);

    run_multi_connection_local_transition_rebind_case(
        rank,
        schedule,
        partition,
        bridge,
        pattern,
        &audit,
        &outer_provider);

    std::unique_ptr<
        fdp::
            PhaseTransitionRebuiltNaturalVariableSystem3D>
        rebuilt_system;
    error =
        fdp::
            rebuild_phase_transition_natural_variable_system_3d(
                PETSC_COMM_WORLD,
                cross_schedule,
                partition,
                bridge,
                pattern,
                1.0,
                std::move(rebuild_cells),
                std::move(rebuild_faces),
                {
                    {&evaluate_1p, &audit},
                    {&evaluate_2p, &audit},
                    {&evaluate_3p, &audit}},
                &fdp::
                    evaluate_absent_phase_thermodynamic_provider_3d<
                        flow::
                            Pr76AbsentPhasePotentialExtensionProvider<
                                double>>,
                &outer_provider,
                &rebuilt_system);
    require_collective(
        error == PETSC_SUCCESS &&
            rebuilt_system != nullptr &&
            rebuilt_system
                    ->numbering()
                    .petsc_global_scalar_count() ==
                45 &&
            rebuilt_system
                    ->numbering()
                    .cell(
                        mesh::LocalIndex{2U})
                    .phase_count ==
                3U &&
            rebuilt_system
                    ->numbering()
                    .cell(
                        mesh::LocalIndex{2U})
                    .scalar_count ==
                10U,
        "accepted 2P->3P transition did not rebuild q-ragged numbering");

    Vec rebuild_residual = nullptr;
    error =
        VecDuplicate(
            rebuilt_system
                ->initial_state(),
            &rebuild_residual);
    if (error == PETSC_SUCCESS) {
        error =
            VecSet(
                rebuild_residual,
                PetscScalar{0.0});
    }
    auto rebuild_evaluator =
        rebuilt_system
            ->snes_evaluator();
    status =
        fdp::
            NaturalVariableSnesEvaluationStatus3D::
                success;
    if (error == PETSC_SUCCESS) {
        error =
            rebuild_evaluator.function(
                rebuilt_system
                    ->initial_state(),
                rebuild_residual,
                rebuild_evaluator.user_context,
                &status);
    }
    if (error == PETSC_SUCCESS) {
        error =
            VecAssemblyBegin(
                rebuild_residual);
    }
    if (error == PETSC_SUCCESS) {
        error =
            VecAssemblyEnd(
                rebuild_residual);
    }
    require_collective(
        error == PETSC_SUCCESS &&
            status ==
                fdp::
                    NaturalVariableSnesEvaluationStatus3D::
                        success,
        "rebuilt phase-transition residual callback failed");

    error =
        MatZeroEntries(
            rebuilt_system
                ->jacobian_structure());
    status =
        fdp::
            NaturalVariableSnesEvaluationStatus3D::
                success;
    if (error == PETSC_SUCCESS) {
        error =
            rebuild_evaluator.jacobian(
                rebuilt_system
                    ->initial_state(),
                rebuilt_system
                    ->jacobian_structure(),
                rebuild_evaluator.user_context,
                &status);
    }
    if (error == PETSC_SUCCESS) {
        error =
            MatAssemblyBegin(
                rebuilt_system
                    ->jacobian_structure(),
                MAT_FINAL_ASSEMBLY);
    }
    if (error == PETSC_SUCCESS) {
        error =
            MatAssemblyEnd(
                rebuilt_system
                    ->jacobian_structure(),
                MAT_FINAL_ASSEMBLY);
    }
    require_collective(
        error == PETSC_SUCCESS &&
            status ==
                fdp::
                    NaturalVariableSnesEvaluationStatus3D::
                        success,
        "rebuilt phase-transition Jacobian callback failed");

    require_collective(
        rebuilt_system
                ->coordinate_registry()
                .reference_entry(
                    mesh::GlobalEntityId{
                        UINT64_C(20)},
                    mixed_physical_phase_identity(
                        "hydrocarbon-1"))
                .selected_branch_provenance ==
            "fixture/hydrocarbon1-root0",
        "outer rebuild did not retain frozen absent-phase branch provenance");

    require_collective(
        VecDestroy(&rebuild_residual) ==
            PETSC_SUCCESS,
        "phase-transition rebuild residual cleanup failed");
    rebuilt_system.reset();

    // Full post-SNES controller: stable after one accepted restart.
    ControllerFixture stable_controller{
        rank,
        &schedule,
        &partition,
        &bridge,
        &pattern,
        &audit,
        &pr_model,
        &outer_provider,
        false,
        0U};
    run_controller_case(
        &stable_controller,
        4U,
        fdp::
            PostSnesPhaseTransitionOutcome3D::
                stable_phase_set,
        1U);

    // Preserved proposal handoff: the first outer-controller scan is not
    // recomputed. It replays the exact local-owned batch captured by the
    // preceding fixed-cardinality solve, then delegates subsequent generations
    // to the normal controller scanner.
    stable_controller.rebuild_calls = 0U;
    auto handoff_initial =
        make_controller_initial_system(
            &stable_controller);
    fdp::PostSnesPhaseTransitionHandoffScannerContext3D
        handoff_scanner;
    if (rank == 0) {
        handoff_scanner
            .initial_local_owned_proposals
            .push_back(
                controller_two_to_three_proposal());
    }
    handoff_scanner.subsequent_scanner =
        &controller_scan;
    handoff_scanner.subsequent_scanner_context =
        &stable_controller;

    std::unique_ptr<
        fdp::
            PhaseTransitionRebuiltNaturalVariableSystem3D>
        handoff_final_system;
    Vec handoff_final_state = nullptr;
    std::optional<
        fdp::
            PostSnesPhaseTransitionControllerReport3D>
        handoff_controller_report;
    error =
        fdp::
            solve_nonlinear_timestep_with_phase_transitions_3d(
                PETSC_COMM_WORLD,
                std::move(
                    handoff_initial),
                {
                    &fdp::
                        scan_post_snes_phase_transition_handoff_3d,
                    &handoff_scanner,
                    &controller_rebuild,
                    &stable_controller},
                {4U},
                &handoff_final_system,
                &handoff_final_state,
                &handoff_controller_report);
    require_collective(
        error == PETSC_SUCCESS &&
            handoff_scanner.initial_batch_consumed &&
            handoff_controller_report.has_value() &&
            handoff_controller_report->outcome ==
                fdp::
                    PostSnesPhaseTransitionOutcome3D::
                        stable_phase_set &&
            handoff_controller_report
                    ->transition_restarts ==
                1U &&
            handoff_controller_report
                    ->generations.size() ==
                2U &&
            handoff_controller_report
                    ->generations.front()
                    .accepted_transition_batch
                    .size() ==
                1U &&
            stable_controller.rebuild_calls ==
                1U &&
            handoff_final_system != nullptr &&
            handoff_final_system
                    ->numbering()
                    .cell(
                        mesh::LocalIndex{2U})
                    .phase_count ==
                3U,
        "preserved proposal did not drive one outer-controller topology restart");
    require_collective(
        VecDestroy(
            &handoff_final_state) ==
            PETSC_SUCCESS,
        "handoff outer-controller final state cleanup failed");

    // Reverse 3P->2P proposal returns to the previously visited global
    // signature and must be rejected as a phase-set cycle before rebuilding.
    ControllerFixture cycle_controller{
        rank,
        &schedule,
        &partition,
        &bridge,
        &pattern,
        &audit,
        &pr_model,
        &outer_provider,
        true,
        0U};
    run_controller_case(
        &cycle_controller,
        4U,
        fdp::
            PostSnesPhaseTransitionOutcome3D::
                phase_set_cycle_detected,
        1U);

    // A zero restart budget permits the converged scan but never mutates the
    // system; the timestep therefore remains unaccepted.
    ControllerFixture budget_controller{
        rank,
        &schedule,
        &partition,
        &bridge,
        &pattern,
        &audit,
        &pr_model,
        &outer_provider,
        false,
        0U};
    run_controller_case(
        &budget_controller,
        0U,
        fdp::
            PostSnesPhaseTransitionOutcome3D::
                transition_restart_budget_exhausted,
        0U);

    ControllerFixture indeterminate_controller{
        rank,
        &schedule,
        &partition,
        &bridge,
        &pattern,
        &audit,
        &pr_model,
        &outer_provider,
        false,
        0U};
    indeterminate_controller.scan_indeterminate =
        true;
    run_controller_case(
        &indeterminate_controller,
        4U,
        fdp::
            PostSnesPhaseTransitionOutcome3D::
                phase_set_scan_indeterminate,
        0U);

    require_collective(
        VecDestroy(&solution) ==
                PETSC_SUCCESS &&
            VecDestroy(&residual) ==
                PETSC_SUCCESS &&
            VecDestroy(&state) ==
                PETSC_SUCCESS &&
            MatDestroy(&jacobian) ==
                PETSC_SUCCESS,
        "mixed physical dispatcher fixture cleanup failed");
}
