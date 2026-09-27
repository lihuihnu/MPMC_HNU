struct ControlledWellPhaseTransitionFixture {
    int rank{};
    const dp::
        ParallelOwnedConnectionSchedule3D*
            schedule{};
    const mesh::PartitionSnapshot*
        partition{};
    const dp::
        PetscMpiAijSymbolicPreallocation3D*
            bridge{};
    const dp::
        OwnedCellStructuralColumnPatternSnapshot3D*
            pattern{};
    DispatchAudit* audit{};
    flow::
        Pr76AbsentPhasePotentialExtensionProvider<
            double>*
        provider{};
    wdp::
        FixedTotalMolarRateWellSourceEvaluatorContext3D*
        source_context{};
    wdp::
        AcceptedFixedTotalMolarRateWellControlState3D*
        accepted_control{};
    fdp::AcceptedPhysicalTimeClock3D*
        clock{};

    wdp::
        FixedTotalMolarRatePhysicalTimestepControlMode3D
        expected_entry_control{
            wdp::
                FixedTotalMolarRatePhysicalTimestepControlMode3D::
                    fixed_total_molar_rate};
    double expected_entry_bhp_pa{};

    std::size_t scans{};
    std::size_t rebuild_calls{};
    bool accepted_anchor_unchanged_at_rebuild{
        true};
    bool rebuilt_previous_totals_captured{};
    std::array<double, 4>
        rebuilt_previous_totals{};
    bool accepted_anchor_unchanged_at_scan{
        true};
    std::vector<double>
        scanned_source_bottom_hole_pressure_pa;
};

PetscErrorCode
controlled_well_transition_scan(
    const fdp::
        PhaseTransitionRebuiltNaturalVariableSystem3D&
            system,
    Vec,
    const fdp::
        VariableCardinalityNaturalVariableSnesSolveReport3D&
            solve_report,
    void* raw_context,
    fdp::PostSnesPhaseTransitionScanStatus3D*
        scan_status,
    std::vector<
        fdp::PostSnesPhaseTransitionProposal3D>*
            output) {
    if (raw_context == nullptr ||
        scan_status == nullptr ||
        output == nullptr) {
        return PETSC_ERR_ARG_NULL;
    }
    auto* fixture =
        static_cast<
            ControlledWellPhaseTransitionFixture*>(
                raw_context);
    *scan_status =
        fdp::
            PostSnesPhaseTransitionScanStatus3D::
                complete;
    output->clear();

    if (static_cast<int>(
            solve_report.converged_reason) <=
            0 ||
        solve_report
                .locally_owned_solution
                .size() !=
            static_cast<std::size_t>(
                system
                    .numbering()
                    .petsc_local_owned_scalar_count())) {
        return PETSC_ERR_ARG_INCOMP;
    }

    if (fixture->source_context == nullptr ||
        fixture->accepted_control == nullptr ||
        fixture->clock == nullptr) {
        return PETSC_ERR_ARG_INCOMP;
    }
    fixture
        ->accepted_anchor_unchanged_at_scan =
        fixture
            ->accepted_anchor_unchanged_at_scan &&
        fixture
                ->accepted_control
                ->control ==
            fixture
                ->expected_entry_control &&
        fixture
                ->accepted_control
                ->bottom_hole_pressure_pa ==
            fixture
                ->expected_entry_bhp_pa &&
        fixture
                ->clock
                ->accepted_time_seconds() ==
            0.0 &&
        fixture
                ->clock
                ->accepted_step_count() ==
            0U;
    fixture
        ->scanned_source_bottom_hole_pressure_pa
        .push_back(
            fixture
                ->source_context
                ->current_bottom_hole_pressure_pa());

    ++fixture->scans;
    const auto& cell30 =
        system.numbering().cell(
            mesh::LocalIndex{2U});
    if (cell30.owner_rank !=
        system.numbering().local_rank()) {
        return PETSC_SUCCESS;
    }
    if (cell30.cell_global !=
            mesh::GlobalEntityId{
                UINT64_C(30)}) {
        return PETSC_ERR_ARG_INCOMP;
    }

    if (cell30.phase_count == 2U) {
        output->push_back(
            controller_two_to_three_proposal());
    } else if (cell30.phase_count != 3U) {
        return PETSC_ERR_ARG_INCOMP;
    }
    return PETSC_SUCCESS;
}

PetscErrorCode
controlled_well_transition_rebuild(
    const fdp::
        PhaseTransitionRebuiltNaturalVariableSystem3D&
            current_system,
    Vec,
    const fdp::
        VariableCardinalityNaturalVariableSnesSolveReport3D&,
    std::span<
        const fdp::
            PostSnesPhaseTransitionProposal3D>
        local_owned_proposals,
    std::span<
        const fdp::
            AcceptedPhaseTransitionSummary3D>
        accepted_global_batch,
    void* raw_context,
    std::unique_ptr<
        fdp::
            PhaseTransitionRebuiltNaturalVariableSystem3D>*
        output) {
    if (raw_context == nullptr ||
        output == nullptr) {
        return PETSC_ERR_ARG_NULL;
    }
    output->reset();
    auto* fixture =
        static_cast<
            ControlledWellPhaseTransitionFixture*>(
                raw_context);
    if (fixture->schedule == nullptr ||
        fixture->partition == nullptr ||
        fixture->bridge == nullptr ||
        fixture->pattern == nullptr ||
        fixture->audit == nullptr ||
        fixture->provider == nullptr ||
        fixture->source_context == nullptr ||
        fixture->accepted_control == nullptr ||
        fixture->clock == nullptr ||
        accepted_global_batch.size() !=
            1U ||
        accepted_global_batch.front()
                .cell_global !=
            mesh::GlobalEntityId{
                UINT64_C(30)} ||
        accepted_global_batch.front()
                .source_phase_count !=
            2U ||
        accepted_global_batch.front()
                .target_phase_count !=
            3U) {
        return PETSC_ERR_ARG_INCOMP;
    }

    if (fixture->rank == 0) {
        if (local_owned_proposals.size() !=
                1U ||
            local_owned_proposals.front()
                    .cell_global !=
                mesh::GlobalEntityId{
                    UINT64_C(30)}) {
            return PETSC_ERR_ARG_INCOMP;
        }
    } else if (
        !local_owned_proposals.empty()) {
        return PETSC_ERR_ARG_INCOMP;
    }

    fixture
        ->accepted_anchor_unchanged_at_rebuild =
        fixture
            ->accepted_anchor_unchanged_at_rebuild &&
        fixture
                ->accepted_control
                ->control ==
            fixture
                ->expected_entry_control &&
        fixture
                ->accepted_control
                ->bottom_hole_pressure_pa ==
            fixture
                ->expected_entry_bhp_pa &&
        fixture
                ->clock
                ->accepted_time_seconds() ==
            0.0 &&
        fixture
                ->clock
                ->accepted_step_count() ==
            0U;

    ++fixture->rebuild_calls;
    auto cells =
        make_controller_cells(
            fixture->rank,
            true,
            fixture->audit);
    auto faces =
        make_face_inputs(
            fixture->rank,
            false);

    PetscErrorCode error =
        fdp::
            rebuild_phase_transition_natural_variable_system_3d(
                PETSC_COMM_WORLD,
                *fixture->schedule,
                *fixture->partition,
                *fixture->bridge,
                *fixture->pattern,
                current_system
                    .time_step_seconds(),
                std::move(cells),
                std::move(faces),
                {
                    {&evaluate_1p,
                     fixture->audit},
                    {&evaluate_2p,
                     fixture->audit},
                    {&evaluate_3p,
                     fixture->audit}},
                wdp::
                    fixed_total_molar_rate_well_source_binding_3d(
                        fixture
                            ->source_context),
                &fdp::
                    evaluate_absent_phase_thermodynamic_provider_3d<
                        flow::
                            Pr76AbsentPhasePotentialExtensionProvider<
                                double>>,
                fixture->provider,
                output);
    if (error != PETSC_SUCCESS ||
        *output == nullptr) {
        return error != PETSC_SUCCESS
            ? error
            : PETSC_ERR_PLIB;
    }

    const auto local_previous =
        owned_conserved_totals(
            **output,
            (*output)->initial_state());
    if (MPI_Allreduce(
            local_previous.data(),
            fixture
                ->rebuilt_previous_totals
                .data(),
            4,
            MPI_DOUBLE,
            MPI_SUM,
            PETSC_COMM_WORLD) !=
        MPI_SUCCESS) {
        output->reset();
        return PETSC_ERR_MPI;
    }
    fixture
        ->rebuilt_previous_totals_captured =
        true;
    return PETSC_SUCCESS;
}

void run_fixed_total_molar_rate_control_case(
    int rank,
    const dp::
        ParallelOwnedConnectionSchedule3D&
            schedule,
    const mesh::PartitionSnapshot&
        partition,
    const dp::
        PetscMpiAijSymbolicPreallocation3D&
            bridge,
    const dp::
        OwnedCellStructuralColumnPatternSnapshot3D&
            pattern,
    DispatchAudit* audit) {
    if (audit == nullptr) {
        throw std::invalid_argument(
            "invalid fixed-total-molar-rate control fixture");
    }

    audit->well_timestep_compressibility =
        true;

    constexpr double initial_bhp_pa = 5.0;
    constexpr double target_reference_bhp_pa = 4.0;

    auto connection30 =
        wdp::
            FixedBhpPeacemanWellSourceEvaluatorContext3D::
                create(
                    mesh::GlobalEntityId{
                        UINT64_C(30)},
                    well::make_peaceman_well_index_3d(
                        {10.0, 10.0, 5.0},
                        {
                            1.0e-8,
                            1.0e-8,
                            1.0e-8},
                        well::
                            AxisAlignedWellDirection3D::z,
                        0.10,
                        0.0),
                    initial_bhp_pa,
                    wd::FixedBhpInjectionEnthalpy3D{
                        "fixture/rate-control-cell30-enthalpy/v1",
                        {1000.0, 2000.0}},
                    "fixture/rate-control-cell30-source/v1");
    auto connection60 =
        wdp::
            FixedBhpPeacemanWellSourceEvaluatorContext3D::
                create(
                    mesh::GlobalEntityId{
                        UINT64_C(60)},
                    well::make_peaceman_well_index_3d(
                        {12.0, 9.0, 4.0},
                        {
                            2.0e-8,
                            1.5e-8,
                            2.5e-8},
                        well::
                            AxisAlignedWellDirection3D::z,
                        0.12,
                        0.1),
                    initial_bhp_pa,
                    wd::FixedBhpInjectionEnthalpy3D{
                        "fixture/rate-control-cell60-enthalpy/v1",
                        {1100.0, 2200.0, 3300.0}},
                    "fixture/rate-control-cell60-source/v1");

    auto multi_context =
        wdp::
            FixedBhpMultiConnectionWellSourceEvaluatorContext3D::
                create(
                    "fixture/rate-controlled-two-connection-well/v1",
                    {
                        connection60,
                        connection30});

    auto initial30 =
        evaluate_target(
            UINT64_C(30),
            audit);
    auto initial60 =
        evaluate_target(
            UINT64_C(60),
            audit);

    const auto target_source30 =
        wdp::
            build_fixed_bhp_peaceman_well_source_3d(
                multi_context
                    .find_connection(
                        mesh::GlobalEntityId{
                            UINT64_C(30)})
                    ->with_bottom_hole_pressure(
                        target_reference_bhp_pa),
                initial30);
    const auto target_source60 =
        wdp::
            build_fixed_bhp_peaceman_well_source_3d(
                multi_context
                    .find_connection(
                        mesh::GlobalEntityId{
                            UINT64_C(60)})
                    ->with_bottom_hole_pressure(
                        target_reference_bhp_pa),
                initial60);

    const auto initial_source30 =
        wdp::
            build_fixed_bhp_peaceman_well_source_3d(
                *multi_context.find_connection(
                    mesh::GlobalEntityId{
                        UINT64_C(30)}),
                initial30);
    const auto initial_source60 =
        wdp::
            build_fixed_bhp_peaceman_well_source_3d(
                *multi_context.find_connection(
                    mesh::GlobalEntityId{
                        UINT64_C(60)}),
                initial60);

    auto total_molar_production =
        [](const auto& first,
           const auto& second) {
            double total = 0.0;
            for (double value :
                 first
                     .cell_source
                     .component_molar_rate_mol_per_s) {
                total -= value;
            }
            for (double value :
                 second
                     .cell_source
                     .component_molar_rate_mol_per_s) {
                total -= value;
            }
            return total;
        };

    const double target_rate =
        total_molar_production(
            target_source30,
            target_source60);
    const double initial_rate =
        total_molar_production(
            initial_source30,
            initial_source60);
    require_collective(
        std::isfinite(target_rate) &&
            target_rate > 0.0 &&
            std::isfinite(initial_rate) &&
            initial_rate > 0.0 &&
            std::abs(
                target_rate -
                initial_rate) >
                1.0e-8 *
                    std::max(
                        target_rate,
                        initial_rate),
        "fixed-total-molar-rate target does not require a nontrivial BHP correction");

    auto source_context =
        wdp::
            FixedTotalMolarRateWellSourceEvaluatorContext3D::
                create(
                    multi_context,
                    initial_bhp_pa);

    const auto pr_parameters =
        thermodynamic_adapter_pr_parameters();
    const auto pr_model =
        th::Pr76Phase<double>::
            from_parameters(
                pr_parameters);
    auto provider =
        make_outer_rebuild_pr_provider(
            pr_model);

    auto reservoir_system =
        make_fixed_total_molar_rate_reservoir_system(
            rank,
            schedule,
            partition,
            bridge,
            pattern,
            audit,
            &source_context,
            &provider);

    require_collective(
        reservoir_system != nullptr &&
            reservoir_system
                    ->numbering()
                    .petsc_global_scalar_count() ==
                42 &&
            reservoir_system
                    ->numbering()
                    .cell(
                        mesh::LocalIndex{2U})
                    .phase_count ==
                2U &&
            reservoir_system
                    ->numbering()
                    .cell(
                        mesh::LocalIndex{5U})
                    .phase_count ==
                3U,
        "rate-control reservoir system changed frozen phase cardinality");

    const auto local_previous_total =
        owned_conserved_totals(
            *reservoir_system,
            reservoir_system
                ->initial_state());
    std::array<double, 4>
        global_previous_total{};
    require_collective(
        MPI_Allreduce(
            local_previous_total.data(),
            global_previous_total.data(),
            4,
            MPI_DOUBLE,
            MPI_SUM,
            PETSC_COMM_WORLD) ==
            MPI_SUCCESS,
        "failed to reduce rate-control previous conserved totals");

    std::unique_ptr<
        wdp::
            FixedTotalMolarRateWellControlSystem3D>
        control_system;
    PetscErrorCode error =
        wdp::
            FixedTotalMolarRateWellControlSystem3D::
                create(
                    PETSC_COMM_WORLD,
                    reservoir_system.get(),
                    &source_context,
                    initial_bhp_pa,
                    target_rate,
                    &control_system);
    require_collective(
        error == PETSC_SUCCESS &&
            control_system != nullptr &&
            control_system
                    ->reservoir_global_scalar_count() ==
                42 &&
            control_system
                    ->global_scalar_count() ==
                43 &&
            control_system
                    ->well_global_scalar() ==
                42 &&
            control_system
                    ->well_owner_rank() ==
                1,
        "rate-control augmented numbering did not add exactly one authoritative BHP scalar");

    Vec solution = nullptr;
    std::optional<
        wdp::
            FixedTotalMolarRateWellControlSolveReport3D>
        report;
    error =
        control_system->solve(
            &solution,
            &report);

    require_collective(
        error == PETSC_SUCCESS &&
            solution != nullptr &&
            report.has_value() &&
            report->converged() &&
            report->global_scalar_count ==
                43 &&
            report->well_global_scalar ==
                42 &&
            report->well_owner_rank ==
                1 &&
            report->snes_type ==
                SNESNEWTONLS &&
            report->ksp_type ==
                KSPGMRES &&
            report->pc_type ==
                PCASM &&
            std::isfinite(
                report
                    ->bottom_hole_pressure_pa) &&
            report->bottom_hole_pressure_pa >
                0.0 &&
            std::abs(
                report
                    ->bottom_hole_pressure_pa -
                initial_bhp_pa) >
                1.0e-7 &&
            std::isfinite(
                report
                    ->achieved_total_molar_rate_mol_per_s) &&
            std::abs(
                report
                    ->total_molar_rate_residual_mol_per_s()) <=
                1.0e-8 *
                    std::max(
                        1.0,
                        std::abs(
                            target_rate)),
        "fixed-total-molar-rate SNES did not converge to a nontrivial BHP/control root");

    Mat controlled_jacobian = nullptr;
    error =
        MatDuplicate(
            control_system
                ->jacobian_structure(),
            MAT_DO_NOT_COPY_VALUES,
            &controlled_jacobian);
    fdp::NaturalVariableSnesEvaluationStatus3D
        controlled_status =
            fdp::
                NaturalVariableSnesEvaluationStatus3D::
                    success;
    auto controlled_evaluator =
        control_system
            ->snes_evaluator();
    if (error == PETSC_SUCCESS) {
        error =
            controlled_evaluator.jacobian(
                solution,
                controlled_jacobian,
                controlled_evaluator
                    .user_context,
                &controlled_status);
    }
    if (error == PETSC_SUCCESS) {
        error =
            MatAssemblyBegin(
                controlled_jacobian,
                MAT_FINAL_ASSEMBLY);
    }
    if (error == PETSC_SUCCESS) {
        error =
            MatAssemblyEnd(
                controlled_jacobian,
                MAT_FINAL_ASSEMBLY);
    }

    bool local_four_blocks_ok =
        error == PETSC_SUCCESS &&
        controlled_status ==
            fdp::
                NaturalVariableSnesEvaluationStatus3D::
                    success;
    const PetscInt well_row =
        control_system
            ->well_global_scalar();

    if (local_four_blocks_ok &&
        rank == 0) {
        const auto& cell30 =
            reservoir_system
                ->numbering()
                .cell(
                    mesh::LocalIndex{2U});
        const PetscInt row =
            cell30
                .petsc_global_scalar_start;
        const PetscInt column =
            row;
        PetscScalar jrr = 0.0;
        PetscScalar jrw = 0.0;
        if (MatGetValues(
                controlled_jacobian,
                1,
                &row,
                1,
                &column,
                &jrr) !=
                PETSC_SUCCESS ||
            MatGetValues(
                controlled_jacobian,
                1,
                &row,
                1,
                &well_row,
                &jrw) !=
                PETSC_SUCCESS ||
            std::abs(
                static_cast<double>(
                    PetscRealPart(
                        jrr))) <=
                0.0 ||
            std::abs(
                static_cast<double>(
                    PetscRealPart(
                        jrw))) <=
                0.0) {
            local_four_blocks_ok =
                false;
        }
    }

    if (local_four_blocks_ok &&
        rank == 1) {
        const auto& cell30 =
            reservoir_system
                ->numbering()
                .cell(
                    mesh::LocalIndex{2U});
        const auto& cell60 =
            reservoir_system
                ->numbering()
                .cell(
                    mesh::LocalIndex{5U});
        const PetscInt cell30_column =
            cell30
                .petsc_global_scalar_start;
        const PetscInt cell60_column =
            cell60
                .petsc_global_scalar_start;
        const PetscInt cell60_row =
            cell60
                .petsc_global_scalar_start;
        PetscScalar jrr = 0.0;
        PetscScalar jrw = 0.0;
        PetscScalar jwr30 = 0.0;
        PetscScalar jwr60 = 0.0;
        PetscScalar jww = 0.0;
        if (MatGetValues(
                controlled_jacobian,
                1,
                &cell60_row,
                1,
                &cell60_column,
                &jrr) !=
                PETSC_SUCCESS ||
            MatGetValues(
                controlled_jacobian,
                1,
                &cell60_row,
                1,
                &well_row,
                &jrw) !=
                PETSC_SUCCESS ||
            MatGetValues(
                controlled_jacobian,
                1,
                &well_row,
                1,
                &cell30_column,
                &jwr30) !=
                PETSC_SUCCESS ||
            MatGetValues(
                controlled_jacobian,
                1,
                &well_row,
                1,
                &cell60_column,
                &jwr60) !=
                PETSC_SUCCESS ||
            MatGetValues(
                controlled_jacobian,
                1,
                &well_row,
                1,
                &well_row,
                &jww) !=
                PETSC_SUCCESS ||
            std::abs(
                static_cast<double>(
                    PetscRealPart(
                        jrr))) <=
                0.0 ||
            std::abs(
                static_cast<double>(
                    PetscRealPart(
                        jrw))) <=
                0.0 ||
            std::abs(
                static_cast<double>(
                    PetscRealPart(
                        jwr30))) <=
                0.0 ||
            std::abs(
                static_cast<double>(
                    PetscRealPart(
                        jwr60))) <=
                0.0 ||
            std::abs(
                static_cast<double>(
                    PetscRealPart(
                        jww))) <=
                0.0) {
            local_four_blocks_ok =
                false;
        }
    }
    require_collective(
        local_four_blocks_ok,
        "rate-control Jacobian did not materialize analytic Jrr/Jrw/Jwr/Jww blocks");

    Vec reservoir_solution = nullptr;
    error =
        control_system
            ->copy_reservoir_state(
                solution,
                &reservoir_solution);
    require_collective(
        error == PETSC_SUCCESS &&
            reservoir_solution != nullptr,
        "failed to extract reservoir state from rate-controlled solution");

    const auto local_final_total =
        owned_conserved_totals(
            *reservoir_system,
            reservoir_solution);
    std::array<double, 4>
        global_final_total{};
    require_collective(
        MPI_Allreduce(
            local_final_total.data(),
            global_final_total.data(),
            4,
            MPI_DOUBLE,
            MPI_SUM,
            PETSC_COMM_WORLD) ==
            MPI_SUCCESS,
        "failed to reduce rate-control final conserved totals");

    std::vector<std::optional<
        fdp::
            MixedCardinalityPhysicalCurrentCellLinearization3D>>
        final_current;
    std::vector<double>
        final_porosity;
    fdp::NaturalVariableSnesEvaluationStatus3D
        final_status =
            fdp::
                NaturalVariableSnesEvaluationStatus3D::
                    success;
    require_collective(
        reservoir_system
                ->evaluate_local_cells_for_phase_transition(
                    reservoir_solution,
                    &final_current,
                    &final_porosity,
                    &final_status) ==
            PETSC_SUCCESS &&
            final_status ==
                fdp::
                    NaturalVariableSnesEvaluationStatus3D::
                        success &&
            final_current.size() > 5U &&
            final_current[2U]
                .has_value() &&
            final_current[5U]
                .has_value(),
        "failed to evaluate rate-controlled reservoir solution");

    std::array<double, 4>
        local_well_rate{};
    std::uint64_t
        local_authoritative_count = 0U;

    const auto& record30 =
        reservoir_system
            ->numbering()
            .cell(
                mesh::LocalIndex{2U});
    if (record30.owner_rank ==
        reservoir_system
            ->numbering()
            .local_rank()) {
        const auto source =
            wdp::
                build_fixed_bhp_peaceman_well_source_3d(
                    multi_context
                        .find_connection(
                            record30
                                .cell_global)
                        ->with_bottom_hole_pressure(
                            report
                                ->bottom_hole_pressure_pa),
                    *final_current[2U]);
        for (std::size_t component = 0U;
             component < 3U;
             ++component) {
            local_well_rate[
                component] -=
                source
                    .cell_source
                    .component_molar_rate_mol_per_s[
                        component];
        }
        local_well_rate[3] -=
            source
                .cell_source
                .energy_rate_w;
        ++local_authoritative_count;
    }

    const auto& record60 =
        reservoir_system
            ->numbering()
            .cell(
                mesh::LocalIndex{5U});
    if (record60.owner_rank ==
        reservoir_system
            ->numbering()
            .local_rank()) {
        const auto source =
            wdp::
                build_fixed_bhp_peaceman_well_source_3d(
                    multi_context
                        .find_connection(
                            record60
                                .cell_global)
                        ->with_bottom_hole_pressure(
                            report
                                ->bottom_hole_pressure_pa),
                    *final_current[5U]);
        for (std::size_t component = 0U;
             component < 3U;
             ++component) {
            local_well_rate[
                component] -=
                source
                    .cell_source
                    .component_molar_rate_mol_per_s[
                        component];
        }
        local_well_rate[3] -=
            source
                .cell_source
                .energy_rate_w;
        ++local_authoritative_count;
    }

    std::array<double, 4>
        global_well_rate{};
    std::uint64_t
        global_authoritative_count = 0U;
    require_collective(
        MPI_Allreduce(
            local_well_rate.data(),
            global_well_rate.data(),
            4,
            MPI_DOUBLE,
            MPI_SUM,
            PETSC_COMM_WORLD) ==
                MPI_SUCCESS &&
            MPI_Allreduce(
                &local_authoritative_count,
                &global_authoritative_count,
                1,
                MPI_UINT64_T,
                MPI_SUM,
                PETSC_COMM_WORLD) ==
                MPI_SUCCESS &&
            global_authoritative_count ==
                2U,
        "rate-controlled well did not retain two unique authoritative connections");

    const double summed_molar_rate =
        global_well_rate[0] +
        global_well_rate[1] +
        global_well_rate[2];
    near_collective(
        summed_molar_rate,
        target_rate,
        1.0e-8,
        1.0e-10);

    const double dt =
        reservoir_system
            ->time_step_seconds();
    for (std::size_t quantity = 0U;
         quantity < 4U;
         ++quantity) {
        near_collective(
            global_final_total[
                quantity],
            global_previous_total[
                quantity] -
                dt *
                    global_well_rate[
                        quantity],
            3.0e-7,
            quantity < 3U
                ? 3.0e-8
                : 3.0e-7);
    }

    // Re-run the same physical step through the production rate-control
    // timestep bridge. The earlier direct solve remains an independent
    // algebra/Jacobian/conservation oracle; only this bridge is allowed to
    // advance accepted reservoir history/time and the accepted BHP initial
    // guess.
    source_context.begin_evaluation(
        initial_bhp_pa);
    double accepted_bhp_pa =
        initial_bhp_pa;
    fdp::AcceptedPhysicalTimeClock3D
        rate_control_clock{
            0.0,
            reservoir_system
                ->time_step_seconds()};
    wdp::
        FixedTotalMolarRatePhysicalTimestepDriverOptions3D
        timestep_options;
    timestep_options.adaptive
        .minimum_timestep_seconds =
        0.25;
    timestep_options.adaptive
        .maximum_timestep_seconds =
        4.0;
    timestep_options.adaptive
        .maximum_retries =
        2U;

    std::optional<
        wdp::
            FixedTotalMolarRatePhysicalTimestepDriverReport3D>
        timestep_report;
    error =
        wdp::
            advance_fixed_total_molar_rate_controlled_physical_timestep_3d(
                PETSC_COMM_WORLD,
                reservoir_system.get(),
                &source_context,
                target_rate,
                timestep_options,
                &rate_control_clock,
                &accepted_bhp_pa,
                &timestep_report);

    require_collective(
        error == PETSC_SUCCESS &&
            timestep_report.has_value() &&
            timestep_report->accepted() &&
            timestep_report
                    ->accepted_record
                    ->accepted_step_index ==
                0U &&
            timestep_report
                    ->accepted_record
                    ->time_n_seconds ==
                0.0 &&
            timestep_report
                    ->accepted_record
                    ->time_np1_seconds ==
                1.0 &&
            timestep_report
                    ->accepted_record
                    ->accepted_timestep_seconds ==
                1.0 &&
            timestep_report
                    ->accepted_solve
                    ->global_scalar_count ==
                43 &&
            timestep_report
                    ->accepted_solve
                    ->well_global_scalar ==
                42 &&
            timestep_report
                    ->accepted_solve
                    ->well_owner_rank ==
                1 &&
            std::isfinite(
                accepted_bhp_pa) &&
            accepted_bhp_pa > 0.0 &&
            std::abs(
                accepted_bhp_pa -
                initial_bhp_pa) >
                1.0e-7 &&
            std::abs(
                timestep_report
                    ->accepted_solve
                    ->total_molar_rate_residual_mol_per_s()) <=
                1.0e-8 *
                    std::max(
                        1.0,
                        std::abs(
                            target_rate)) &&
            rate_control_clock
                    .accepted_time_seconds() ==
                1.0 &&
            rate_control_clock
                    .accepted_step_count() ==
                1U &&
            source_context
                    .current_bottom_hole_pressure_pa() ==
                accepted_bhp_pa,
        "fixed-total-molar-rate physical-timestep driver did not commit reservoir time/history and BHP exactly once");

    near_collective(
        accepted_bhp_pa,
        report->bottom_hole_pressure_pa,
        2.0e-8,
        2.0e-10);

    Vec committed_difference = nullptr;
    error =
        VecDuplicate(
            reservoir_system
                ->initial_state(),
            &committed_difference);
    if (error == PETSC_SUCCESS) {
        error =
            VecCopy(
                reservoir_system
                    ->initial_state(),
                committed_difference);
    }
    if (error == PETSC_SUCCESS) {
        error =
            VecAXPY(
                committed_difference,
                PetscScalar{-1.0},
                reservoir_solution);
    }
    PetscReal committed_difference_norm =
        -1.0;
    if (error == PETSC_SUCCESS) {
        error =
            VecNorm(
                committed_difference,
                NORM_2,
                &committed_difference_norm);
    }
    require_collective(
        error == PETSC_SUCCESS &&
            committed_difference_norm <=
                2.0e-8,
        "rate-control timestep driver committed a different reservoir root from the independently verified augmented solve");

    bool history_matches = false;
    require_collective(
        reservoir_system
                ->accepted_history_matches_state(
                    reservoir_system
                        ->initial_state(),
                    &history_matches) ==
            PETSC_SUCCESS &&
            history_matches,
        "rate-controlled accepted reservoir history did not rebase through the physical-timestep driver");

    Vec accepted_snapshot = nullptr;
    error =
        VecDuplicate(
            reservoir_system
                ->initial_state(),
            &accepted_snapshot);
    if (error == PETSC_SUCCESS) {
        error =
            VecCopy(
                reservoir_system
                    ->initial_state(),
                accepted_snapshot);
    }
    require_collective(
        error == PETSC_SUCCESS &&
            accepted_snapshot != nullptr,
        "failed to snapshot accepted rate-control state for rollback regression");

    const double accepted_time_before_rejection =
        rate_control_clock
            .accepted_time_seconds();
    const std::size_t
        accepted_steps_before_rejection =
            rate_control_clock
                .accepted_step_count();
    const double
        next_dt_before_rejection =
            rate_control_clock
                .next_timestep_seconds();
    const double
        bhp_before_rejection =
            accepted_bhp_pa;

    auto invalid_timestep_options =
        timestep_options;
    invalid_timestep_options.adaptive
        .minimum_timestep_seconds =
        next_dt_before_rejection *
        2.0;
    invalid_timestep_options.adaptive
        .maximum_timestep_seconds =
        next_dt_before_rejection *
        4.0;

    std::optional<
        wdp::
            FixedTotalMolarRatePhysicalTimestepDriverReport3D>
        rejected_timestep_report;
    const PetscErrorCode rejection_error =
        wdp::
            advance_fixed_total_molar_rate_controlled_physical_timestep_3d(
                PETSC_COMM_WORLD,
                reservoir_system.get(),
                &source_context,
                target_rate,
                invalid_timestep_options,
                &rate_control_clock,
                &accepted_bhp_pa,
                &rejected_timestep_report);

    Vec rollback_difference = nullptr;
    error =
        VecDuplicate(
            reservoir_system
                ->initial_state(),
            &rollback_difference);
    if (error == PETSC_SUCCESS) {
        error =
            VecCopy(
                reservoir_system
                    ->initial_state(),
                rollback_difference);
    }
    if (error == PETSC_SUCCESS) {
        error =
            VecAXPY(
                rollback_difference,
                PetscScalar{-1.0},
                accepted_snapshot);
    }
    PetscReal rollback_difference_norm =
        -1.0;
    if (error == PETSC_SUCCESS) {
        error =
            VecNorm(
                rollback_difference,
                NORM_2,
                &rollback_difference_norm);
    }
    require_collective(
        rejection_error ==
                PETSC_ERR_ARG_OUTOFRANGE &&
            !rejected_timestep_report
                 .has_value() &&
            error == PETSC_SUCCESS &&
            rollback_difference_norm <=
                1.0e-14 &&
            rate_control_clock
                    .accepted_time_seconds() ==
                accepted_time_before_rejection &&
            rate_control_clock
                    .accepted_step_count() ==
                accepted_steps_before_rejection &&
            rate_control_clock
                    .next_timestep_seconds() ==
                next_dt_before_rejection &&
            reservoir_system
                    ->time_step_seconds() ==
                next_dt_before_rejection &&
            accepted_bhp_pa ==
                bhp_before_rejection &&
            source_context
                    .current_bottom_hole_pressure_pa() ==
                bhp_before_rejection,
        "rejected rate-control timestep entry mutated accepted state/time/dt/BHP ownership");

    require_collective(
        VecDestroy(
            &rollback_difference) ==
                PETSC_SUCCESS &&
            VecDestroy(
                &accepted_snapshot) ==
                PETSC_SUCCESS &&
            VecDestroy(
                &committed_difference) ==
                PETSC_SUCCESS,
        "rate-control timestep bridge regression cleanup failed");

    // Exercise a real nonlinear retry rather than a test hook.  With a fixed
    // positive whole-well production target, an enormous physical timestep
    // would require negative global molar inventory.  The frozen natural-
    // variable domains cannot represent that state, so the augmented SNES
    // must reject the first trial.  A single cutback returns exactly to the
    // independently validated 1 s problem above.
    auto cutback_source_context =
        wdp::
            FixedTotalMolarRateWellSourceEvaluatorContext3D::
                create(
                    multi_context,
                    initial_bhp_pa);
    auto cutback_reservoir_system =
        make_fixed_total_molar_rate_reservoir_system(
            rank,
            schedule,
            partition,
            bridge,
            pattern,
            audit,
            &cutback_source_context,
            &provider);

    constexpr double impossible_first_dt_seconds =
        1.0e9;
    require_collective(
        cutback_reservoir_system != nullptr &&
            cutback_reservoir_system
                    ->set_trial_timestep_seconds(
                        impossible_first_dt_seconds) ==
                PETSC_SUCCESS,
        "failed to prepare rate-control cutback baseline");

    fdp::AcceptedPhysicalTimeClock3D
        cutback_clock{
            0.0,
            impossible_first_dt_seconds};
    wdp::
        FixedTotalMolarRatePhysicalTimestepDriverOptions3D
        cutback_options;
    cutback_options.adaptive
        .minimum_timestep_seconds =
        1.0;
    cutback_options.adaptive
        .maximum_timestep_seconds =
        impossible_first_dt_seconds;
    cutback_options.adaptive
        .cutback_factor =
        1.0e-9;
    cutback_options.adaptive
        .growth_factor =
        2.0;
    cutback_options.adaptive
        .maximum_retries =
        1U;
    cutback_options.adaptive
        .growth_nonlinear_iteration_limit =
        30;
    cutback_options.adaptive
        .growth_line_search_direction_change_limit =
        30;

    double cutback_accepted_bhp_pa =
        initial_bhp_pa;
    std::optional<
        wdp::
            FixedTotalMolarRatePhysicalTimestepDriverReport3D>
        cutback_report;
    error =
        wdp::
            advance_fixed_total_molar_rate_controlled_physical_timestep_3d(
                PETSC_COMM_WORLD,
                cutback_reservoir_system.get(),
                &cutback_source_context,
                target_rate,
                cutback_options,
                &cutback_clock,
                &cutback_accepted_bhp_pa,
                &cutback_report);

    require_collective(
        error == PETSC_SUCCESS &&
            cutback_report.has_value() &&
            cutback_report->accepted() &&
            cutback_report
                    ->adaptive
                    .attempts
                    .size() ==
                2U &&
            cutback_report
                    ->adaptive
                    .retries ==
                1U &&
            cutback_report
                    ->adaptive
                    .accepted_timestep_seconds
                    .has_value() &&
            *cutback_report
                  ->adaptive
                  .accepted_timestep_seconds ==
                1.0 &&
            cutback_report
                    ->accepted_record
                    ->accepted_timestep_seconds ==
                1.0 &&
            cutback_clock
                    .accepted_time_seconds() ==
                1.0 &&
            cutback_clock
                    .accepted_step_count() ==
                1U,
        "rate-control physical-timestep driver did not cut back one real failed augmented solve to the accepted 1 s retry");

    const auto& failed_attempt =
        cutback_report
            ->adaptive
            .attempts[0U];
    const auto& accepted_retry =
        cutback_report
            ->adaptive
            .attempts[1U];
    require_collective(
        failed_attempt
                .request
                .attempt_index ==
            0U &&
            failed_attempt
                .request
                .retry_index ==
            0U &&
            failed_attempt
                .request
                .timestep_seconds ==
            impossible_first_dt_seconds &&
            failed_attempt
                .result
                .outcome ==
            fdp::
                AdaptiveTimestepAttemptOutcome3D::
                    nonlinear_solve_diverged &&
            failed_attempt
                .decision ==
            fdp::
                AdaptiveTimestepDecision3D::
                    reject_and_cutback &&
            failed_attempt
                .next_attempt_timestep_seconds
                .has_value() &&
            *failed_attempt
                  .next_attempt_timestep_seconds ==
                1.0 &&
            accepted_retry
                .request
                .attempt_index ==
            1U &&
            accepted_retry
                .request
                .retry_index ==
            1U &&
            accepted_retry
                .request
                .timestep_seconds ==
            1.0 &&
            accepted_retry
                .result
                .outcome ==
            fdp::
                AdaptiveTimestepAttemptOutcome3D::
                    stable_phase_set,
        "rate-control adaptive report did not preserve failed-first-attempt/cutback/retry provenance");

    near_collective(
        cutback_accepted_bhp_pa,
        report->bottom_hole_pressure_pa,
        2.0e-8,
        2.0e-10);

    Vec cutback_root_difference = nullptr;
    error =
        VecDuplicate(
            cutback_reservoir_system
                ->initial_state(),
            &cutback_root_difference);
    if (error == PETSC_SUCCESS) {
        error =
            VecCopy(
                cutback_reservoir_system
                    ->initial_state(),
                cutback_root_difference);
    }
    if (error == PETSC_SUCCESS) {
        error =
            VecAXPY(
                cutback_root_difference,
                PetscScalar{-1.0},
                reservoir_solution);
    }
    PetscReal cutback_root_difference_norm =
        -1.0;
    if (error == PETSC_SUCCESS) {
        error =
            VecNorm(
                cutback_root_difference,
                NORM_2,
                &cutback_root_difference_norm);
    }
    bool cutback_history_matches = false;
    if (error == PETSC_SUCCESS) {
        error =
            cutback_reservoir_system
                ->accepted_history_matches_state(
                    cutback_reservoir_system
                        ->initial_state(),
                    &cutback_history_matches);
    }
    require_collective(
        error == PETSC_SUCCESS &&
            cutback_root_difference_norm <=
                2.0e-8 &&
            cutback_history_matches &&
            cutback_source_context
                    .current_bottom_hole_pressure_pa() ==
                cutback_accepted_bhp_pa &&
            std::abs(
                cutback_report
                    ->accepted_solve
                    ->total_molar_rate_residual_mol_per_s()) <=
                1.0e-8 *
                    std::max(
                        1.0,
                        std::abs(
                            target_rate)),
        "rate-control retry reused failed-trial state/BHP instead of rebuilding the validated 1 s accepted root");

    require_collective(
        VecDestroy(
            &cutback_root_difference) ==
            PETSC_SUCCESS,
        "rate-control nonlinear cutback regression cleanup failed");

    // Bind a true producer minimum-BHP constraint between the unconstrained
    // rate-control root and the entry BHP. The fixed-rate candidate must be
    // discarded, then the same 1 s physical timestep must be re-solved from
    // the original accepted reservoir history at exactly p_min.
    require_collective(
        report->bottom_hole_pressure_pa <
            initial_bhp_pa,
        "rate-control switching fixture did not produce a lower-BHP unconstrained rate root");
    const double minimum_bhp_pa =
        0.5 *
        (report->bottom_hole_pressure_pa +
         initial_bhp_pa);
    require_collective(
        minimum_bhp_pa >
                report
                    ->bottom_hole_pressure_pa &&
            minimum_bhp_pa <
                initial_bhp_pa,
        "minimum-BHP switching fixture did not place p_min strictly above the rate-control candidate");

    auto switch_source_context =
        wdp::
            FixedTotalMolarRateWellSourceEvaluatorContext3D::
                create(
                    multi_context,
                    initial_bhp_pa);
    auto switch_reservoir_system =
        make_fixed_total_molar_rate_reservoir_system(
            rank,
            schedule,
            partition,
            bridge,
            pattern,
            audit,
            &switch_source_context,
            &provider);

    const auto switch_local_previous_total =
        owned_conserved_totals(
            *switch_reservoir_system,
            switch_reservoir_system
                ->initial_state());
    std::array<double, 4>
        switch_global_previous_total{};
    require_collective(
        MPI_Allreduce(
            switch_local_previous_total.data(),
            switch_global_previous_total.data(),
            4,
            MPI_DOUBLE,
            MPI_SUM,
            PETSC_COMM_WORLD) ==
            MPI_SUCCESS,
        "failed to reduce minimum-BHP switching previous conserved totals");

    fdp::AcceptedPhysicalTimeClock3D
        switch_clock{
            0.0,
            1.0};
    wdp::
        FixedTotalMolarRatePhysicalTimestepDriverOptions3D
        switch_options;
    switch_options.adaptive
        .minimum_timestep_seconds =
        0.25;
    switch_options.adaptive
        .maximum_timestep_seconds =
        4.0;
    switch_options.adaptive
        .maximum_retries =
        2U;
    switch_options.minimum_bottom_hole_pressure_pa =
        minimum_bhp_pa;

    auto switch_control_state =
        wdp::
            AcceptedFixedTotalMolarRateWellControlState3D::
                fixed_total_molar_rate(
                    initial_bhp_pa);
    std::optional<
        wdp::
            FixedTotalMolarRatePhysicalTimestepDriverReport3D>
        switch_report;
    error =
        wdp::
            advance_fixed_total_molar_rate_controlled_physical_timestep_3d(
                PETSC_COMM_WORLD,
                switch_reservoir_system.get(),
                &switch_source_context,
                target_rate,
                switch_options,
                &switch_clock,
                &switch_control_state,
                &switch_report);

    require_collective(
        error == PETSC_SUCCESS &&
            switch_report.has_value() &&
            switch_report->accepted() &&
            switch_report
                    ->rate_to_bhp_switch_triggered &&
            switch_report
                    ->accepted_control ==
                wdp::
                    FixedTotalMolarRatePhysicalTimestepControlMode3D::
                        minimum_bottom_hole_pressure &&
            !switch_report
                 ->accepted_solve
                 .has_value() &&
            switch_report
                    ->accepted_fixed_bhp_solve
                    .has_value() &&
            switch_report
                    ->discarded_rate_control_candidate
                    .has_value() &&
            switch_report
                    ->discarded_rate_control_candidate
                    ->bottom_hole_pressure_pa <
                minimum_bhp_pa &&
            switch_report
                    ->accepted_record
                    ->accepted_timestep_seconds ==
                1.0 &&
            switch_clock
                    .accepted_time_seconds() ==
                1.0 &&
            switch_clock
                    .accepted_step_count() ==
                1U &&
            switch_control_state
                    .control ==
                wdp::
                    FixedTotalMolarRatePhysicalTimestepControlMode3D::
                        minimum_bottom_hole_pressure &&
            switch_control_state
                    .bottom_hole_pressure_pa ==
                minimum_bhp_pa &&
            switch_source_context
                    .current_bottom_hole_pressure_pa() ==
                minimum_bhp_pa,
        "minimum-BHP constraint did not discard the violating rate candidate and accept a same-timestep fixed-BHP re-solve");

    near_collective(
        switch_report
            ->discarded_rate_control_candidate
            ->bottom_hole_pressure_pa,
        report->bottom_hole_pressure_pa,
        2.0e-8,
        2.0e-10);

    std::vector<std::optional<
        fdp::
            MixedCardinalityPhysicalCurrentCellLinearization3D>>
        switch_final_current;
    std::vector<double>
        switch_final_porosity;
    fdp::NaturalVariableSnesEvaluationStatus3D
        switch_final_status =
            fdp::
                NaturalVariableSnesEvaluationStatus3D::
                    success;
    require_collective(
        switch_reservoir_system
                ->evaluate_local_cells_for_phase_transition(
                    switch_reservoir_system
                        ->initial_state(),
                    &switch_final_current,
                    &switch_final_porosity,
                    &switch_final_status) ==
            PETSC_SUCCESS &&
            switch_final_status ==
                fdp::
                    NaturalVariableSnesEvaluationStatus3D::
                        success &&
            switch_final_current.size() > 5U &&
            switch_final_current[2U]
                .has_value() &&
            switch_final_current[5U]
                .has_value(),
        "failed to evaluate accepted minimum-BHP reservoir state");

    std::array<double, 4>
        switch_local_well_rate{};
    std::uint64_t
        switch_local_authoritative_count = 0U;
    const auto& switch_record30 =
        switch_reservoir_system
            ->numbering()
            .cell(
                mesh::LocalIndex{2U});
    if (switch_record30.owner_rank ==
        switch_reservoir_system
            ->numbering()
            .local_rank()) {
        const auto source =
            wdp::
                build_fixed_bhp_peaceman_well_source_3d(
                    multi_context
                        .find_connection(
                            switch_record30
                                .cell_global)
                        ->with_bottom_hole_pressure(
                            minimum_bhp_pa),
                    *switch_final_current[2U]);
        for (std::size_t component = 0U;
             component < 3U;
             ++component) {
            switch_local_well_rate[
                component] -=
                source
                    .cell_source
                    .component_molar_rate_mol_per_s[
                        component];
        }
        switch_local_well_rate[3] -=
            source
                .cell_source
                .energy_rate_w;
        ++switch_local_authoritative_count;
    }

    const auto& switch_record60 =
        switch_reservoir_system
            ->numbering()
            .cell(
                mesh::LocalIndex{5U});
    if (switch_record60.owner_rank ==
        switch_reservoir_system
            ->numbering()
            .local_rank()) {
        const auto source =
            wdp::
                build_fixed_bhp_peaceman_well_source_3d(
                    multi_context
                        .find_connection(
                            switch_record60
                                .cell_global)
                        ->with_bottom_hole_pressure(
                            minimum_bhp_pa),
                    *switch_final_current[5U]);
        for (std::size_t component = 0U;
             component < 3U;
             ++component) {
            switch_local_well_rate[
                component] -=
                source
                    .cell_source
                    .component_molar_rate_mol_per_s[
                        component];
        }
        switch_local_well_rate[3] -=
            source
                .cell_source
                .energy_rate_w;
        ++switch_local_authoritative_count;
    }

    std::array<double, 4>
        switch_global_well_rate{};
    std::uint64_t
        switch_global_authoritative_count = 0U;
    require_collective(
        MPI_Allreduce(
            switch_local_well_rate.data(),
            switch_global_well_rate.data(),
            4,
            MPI_DOUBLE,
            MPI_SUM,
            PETSC_COMM_WORLD) ==
                MPI_SUCCESS &&
            MPI_Allreduce(
                &switch_local_authoritative_count,
                &switch_global_authoritative_count,
                1,
                MPI_UINT64_T,
                MPI_SUM,
                PETSC_COMM_WORLD) ==
                MPI_SUCCESS &&
            switch_global_authoritative_count ==
                2U,
        "minimum-BHP switched well lost owner-only multi-connection aggregation");

    const double switch_total_molar_rate =
        switch_global_well_rate[0] +
        switch_global_well_rate[1] +
        switch_global_well_rate[2];
    require_collective(
        std::isfinite(
            switch_total_molar_rate) &&
            switch_total_molar_rate <
                target_rate -
                    1.0e-10 *
                        std::max(
                            1.0,
                            std::abs(
                                target_rate)),
        "minimum-BHP constrained solution still forced the infeasible fixed-rate target");

    const auto switch_local_final_total =
        owned_conserved_totals(
            *switch_reservoir_system,
            switch_reservoir_system
                ->initial_state());
    std::array<double, 4>
        switch_global_final_total{};
    require_collective(
        MPI_Allreduce(
            switch_local_final_total.data(),
            switch_global_final_total.data(),
            4,
            MPI_DOUBLE,
            MPI_SUM,
            PETSC_COMM_WORLD) ==
            MPI_SUCCESS,
        "failed to reduce minimum-BHP switching final conserved totals");

    for (std::size_t quantity = 0U;
         quantity < 4U;
         ++quantity) {
        near_collective(
            switch_global_final_total[
                quantity],
            switch_global_previous_total[
                quantity] -
                switch_global_well_rate[
                    quantity],
            3.0e-7,
            quantity < 3U
                ? 3.0e-8
                : 3.0e-7);
    }

    bool switch_history_matches = false;
    require_collective(
        switch_reservoir_system
                ->accepted_history_matches_state(
                    switch_reservoir_system
                        ->initial_state(),
                    &switch_history_matches) ==
                PETSC_SUCCESS &&
            switch_history_matches,
        "minimum-BHP switch did not commit exactly the final fixed-BHP reservoir state/history");

    // Advance a second accepted physical timestep with the same accepted
    // well-control state. It must enter directly in minimum-BHP mode: no
    // augmented rate solve, no discarded rate candidate and no new rate->BHP
    // transition are allowed.
    const double second_step_dt =
        switch_clock
            .next_timestep_seconds();
    const double second_step_time_n =
        switch_clock
            .accepted_time_seconds();
    const auto second_step_local_previous_total =
        owned_conserved_totals(
            *switch_reservoir_system,
            switch_reservoir_system
                ->initial_state());
    std::array<double, 4>
        second_step_global_previous_total{};
    require_collective(
        MPI_Allreduce(
            second_step_local_previous_total.data(),
            second_step_global_previous_total.data(),
            4,
            MPI_DOUBLE,
            MPI_SUM,
            PETSC_COMM_WORLD) ==
            MPI_SUCCESS,
        "failed to reduce persisted minimum-BHP second-step previous totals");

    std::optional<
        wdp::
            FixedTotalMolarRatePhysicalTimestepDriverReport3D>
        persisted_bhp_report;
    error =
        wdp::
            advance_fixed_total_molar_rate_controlled_physical_timestep_3d(
                PETSC_COMM_WORLD,
                switch_reservoir_system.get(),
                &switch_source_context,
                target_rate,
                switch_options,
                &switch_clock,
                &switch_control_state,
                &persisted_bhp_report);

    require_collective(
        error == PETSC_SUCCESS &&
            persisted_bhp_report.has_value() &&
            persisted_bhp_report->accepted() &&
            persisted_bhp_report
                    ->entry_control ==
                wdp::
                    FixedTotalMolarRatePhysicalTimestepControlMode3D::
                        minimum_bottom_hole_pressure &&
            persisted_bhp_report
                    ->accepted_control ==
                wdp::
                    FixedTotalMolarRatePhysicalTimestepControlMode3D::
                        minimum_bottom_hole_pressure &&
            !persisted_bhp_report
                 ->rate_to_bhp_switch_triggered &&
            !persisted_bhp_report
                 ->discarded_rate_control_candidate
                 .has_value() &&
            !persisted_bhp_report
                 ->accepted_solve
                 .has_value() &&
            persisted_bhp_report
                    ->accepted_fixed_bhp_solve
                    .has_value() &&
            persisted_bhp_report
                    ->accepted_record
                    ->accepted_timestep_seconds ==
                second_step_dt &&
            switch_clock
                    .accepted_time_seconds() ==
                second_step_time_n +
                    second_step_dt &&
            switch_clock
                    .accepted_step_count() ==
                2U &&
            switch_control_state
                    .control ==
                wdp::
                    FixedTotalMolarRatePhysicalTimestepControlMode3D::
                        minimum_bottom_hole_pressure &&
            switch_control_state
                    .bottom_hole_pressure_pa ==
                minimum_bhp_pa &&
            switch_source_context
                    .current_bottom_hole_pressure_pa() ==
                minimum_bhp_pa,
        "accepted minimum-BHP control state did not persist directly into the next physical timestep");

    std::vector<std::optional<
        fdp::
            MixedCardinalityPhysicalCurrentCellLinearization3D>>
        second_step_final_current;
    std::vector<double>
        second_step_final_porosity;
    fdp::NaturalVariableSnesEvaluationStatus3D
        second_step_final_status =
            fdp::
                NaturalVariableSnesEvaluationStatus3D::
                    success;
    require_collective(
        switch_reservoir_system
                ->evaluate_local_cells_for_phase_transition(
                    switch_reservoir_system
                        ->initial_state(),
                    &second_step_final_current,
                    &second_step_final_porosity,
                    &second_step_final_status) ==
            PETSC_SUCCESS &&
            second_step_final_status ==
                fdp::
                    NaturalVariableSnesEvaluationStatus3D::
                        success &&
            second_step_final_current.size() > 5U &&
            second_step_final_current[2U]
                .has_value() &&
            second_step_final_current[5U]
                .has_value(),
        "failed to evaluate persisted minimum-BHP second-step state");

    std::array<double, 4>
        second_step_local_well_rate{};
    std::uint64_t
        second_step_local_authoritative_count = 0U;
    if (switch_record30.owner_rank ==
        switch_reservoir_system
            ->numbering()
            .local_rank()) {
        const auto source =
            wdp::
                build_fixed_bhp_peaceman_well_source_3d(
                    multi_context
                        .find_connection(
                            switch_record30
                                .cell_global)
                        ->with_bottom_hole_pressure(
                            minimum_bhp_pa),
                    *second_step_final_current[2U]);
        for (std::size_t component = 0U;
             component < 3U;
             ++component) {
            second_step_local_well_rate[
                component] -=
                source
                    .cell_source
                    .component_molar_rate_mol_per_s[
                        component];
        }
        second_step_local_well_rate[3] -=
            source
                .cell_source
                .energy_rate_w;
        ++second_step_local_authoritative_count;
    }
    if (switch_record60.owner_rank ==
        switch_reservoir_system
            ->numbering()
            .local_rank()) {
        const auto source =
            wdp::
                build_fixed_bhp_peaceman_well_source_3d(
                    multi_context
                        .find_connection(
                            switch_record60
                                .cell_global)
                        ->with_bottom_hole_pressure(
                            minimum_bhp_pa),
                    *second_step_final_current[5U]);
        for (std::size_t component = 0U;
             component < 3U;
             ++component) {
            second_step_local_well_rate[
                component] -=
                source
                    .cell_source
                    .component_molar_rate_mol_per_s[
                        component];
        }
        second_step_local_well_rate[3] -=
            source
                .cell_source
                .energy_rate_w;
        ++second_step_local_authoritative_count;
    }

    std::array<double, 4>
        second_step_global_well_rate{};
    std::uint64_t
        second_step_global_authoritative_count = 0U;
    require_collective(
        MPI_Allreduce(
            second_step_local_well_rate.data(),
            second_step_global_well_rate.data(),
            4,
            MPI_DOUBLE,
            MPI_SUM,
            PETSC_COMM_WORLD) ==
                MPI_SUCCESS &&
            MPI_Allreduce(
                &second_step_local_authoritative_count,
                &second_step_global_authoritative_count,
                1,
                MPI_UINT64_T,
                MPI_SUM,
                PETSC_COMM_WORLD) ==
                MPI_SUCCESS &&
            second_step_global_authoritative_count ==
                2U,
        "persisted minimum-BHP second step lost owner-only well aggregation");

    const auto second_step_local_final_total =
        owned_conserved_totals(
            *switch_reservoir_system,
            switch_reservoir_system
                ->initial_state());
    std::array<double, 4>
        second_step_global_final_total{};
    require_collective(
        MPI_Allreduce(
            second_step_local_final_total.data(),
            second_step_global_final_total.data(),
            4,
            MPI_DOUBLE,
            MPI_SUM,
            PETSC_COMM_WORLD) ==
            MPI_SUCCESS,
        "failed to reduce persisted minimum-BHP second-step final totals");

    for (std::size_t quantity = 0U;
         quantity < 4U;
         ++quantity) {
        near_collective(
            second_step_global_final_total[
                quantity],
            second_step_global_previous_total[
                quantity] -
                second_step_dt *
                    second_step_global_well_rate[
                        quantity],
            3.0e-7,
            quantity < 3U
                ? 3.0e-8
                : 3.0e-7);
    }

    // A rejected third entry must not mutate the already accepted BHP control
    // state. This checks control-mode rollback together with the existing
    // reservoir history/time ownership.
    const auto control_before_rejection =
        switch_control_state;
    const double time_before_control_rejection =
        switch_clock
            .accepted_time_seconds();
    const std::size_t steps_before_control_rejection =
        switch_clock
            .accepted_step_count();
    const double dt_before_control_rejection =
        switch_clock
            .next_timestep_seconds();

    auto invalid_persisted_options =
        switch_options;
    invalid_persisted_options.adaptive
        .minimum_timestep_seconds =
        dt_before_control_rejection *
        2.0;
    invalid_persisted_options.adaptive
        .maximum_timestep_seconds =
        dt_before_control_rejection *
        4.0;

    std::optional<
        wdp::
            FixedTotalMolarRatePhysicalTimestepDriverReport3D>
        rejected_persisted_report;
    const PetscErrorCode persisted_rejection_error =
        wdp::
            advance_fixed_total_molar_rate_controlled_physical_timestep_3d(
                PETSC_COMM_WORLD,
                switch_reservoir_system.get(),
                &switch_source_context,
                target_rate,
                invalid_persisted_options,
                &switch_clock,
                &switch_control_state,
                &rejected_persisted_report);

    require_collective(
        persisted_rejection_error ==
                PETSC_ERR_ARG_OUTOFRANGE &&
            !rejected_persisted_report
                 .has_value() &&
            switch_control_state.control ==
                control_before_rejection
                    .control &&
            switch_control_state
                    .bottom_hole_pressure_pa ==
                control_before_rejection
                    .bottom_hole_pressure_pa &&
            switch_clock
                    .accepted_time_seconds() ==
                time_before_control_rejection &&
            switch_clock
                    .accepted_step_count() ==
                steps_before_control_rejection &&
            switch_clock
                    .next_timestep_seconds() ==
                dt_before_control_rejection &&
            switch_reservoir_system
                    ->time_step_seconds() ==
                dt_before_control_rejection &&
            switch_source_context
                    .current_bottom_hole_pressure_pa() ==
                minimum_bhp_pa,
        "rejected persisted minimum-BHP timestep mutated accepted control state/time/history");

    // Build an independent accepted minimum-BHP baseline and measure its
    // p_min capacity without committing it. The guarded reactivation target is
    // then chosen strictly below that capacity, so both release gates can be
    // exercised without introducing schedule logic into the production API.
    auto reactivation_source_context =
        wdp::
            FixedTotalMolarRateWellSourceEvaluatorContext3D::
                create(
                    multi_context,
                    minimum_bhp_pa);
    auto reactivation_reservoir_system =
        make_fixed_total_molar_rate_reservoir_system(
            rank,
            schedule,
            partition,
            bridge,
            pattern,
            audit,
            &reactivation_source_context,
            &provider);

    constexpr double reactivation_dt_seconds =
        0.25;
    require_collective(
        reactivation_reservoir_system !=
                nullptr &&
            reactivation_reservoir_system
                    ->set_trial_timestep_seconds(
                        reactivation_dt_seconds) ==
                PETSC_SUCCESS,
        "failed to prepare minimum-BHP reactivation baseline");

    reactivation_source_context
        .begin_fixed_bhp_reservoir_evaluation(
            minimum_bhp_pa);
    Vec reactivation_probe_state = nullptr;
    std::optional<
        fdp::
            VariableCardinalityNaturalVariableSnesSolveReport3D>
        reactivation_probe_solve;
    error =
        reactivation_reservoir_system
            ->solve(
                &reactivation_probe_state,
                &reactivation_probe_solve);
    require_collective(
        error == PETSC_SUCCESS &&
            reactivation_probe_state !=
                nullptr &&
            reactivation_probe_solve
                .has_value(),
        "failed to solve independent minimum-BHP capacity oracle");

    std::uint64_t
        reactivation_probe_authoritative_count = 0U;
    const auto reactivation_probe_rate =
        global_fixed_bhp_well_production_rate(
            *reactivation_reservoir_system,
            reactivation_probe_state,
            multi_context,
            minimum_bhp_pa,
            &reactivation_probe_authoritative_count);
    const double
        reactivation_probe_total_molar_rate =
            reactivation_probe_rate[0] +
            reactivation_probe_rate[1] +
            reactivation_probe_rate[2];
    require_collective(
        reactivation_probe_authoritative_count ==
                2U &&
            std::isfinite(
                reactivation_probe_total_molar_rate) &&
            reactivation_probe_total_molar_rate >
                0.0,
        "minimum-BHP capacity oracle is not a positive two-connection production rate");

    const double reactivation_target_rate =
        0.50 *
        reactivation_probe_total_molar_rate;
    const double reactivation_rate_margin =
        0.10 *
        reactivation_probe_total_molar_rate;
    constexpr double
        reactivation_pressure_margin_pa =
            1.0e-3;

    const auto reactivation_local_previous_total =
        owned_conserved_totals(
            *reactivation_reservoir_system,
            reactivation_reservoir_system
                ->initial_state());
    std::array<double, 4>
        reactivation_global_previous_total{};
    require_collective(
        MPI_Allreduce(
            reactivation_local_previous_total.data(),
            reactivation_global_previous_total.data(),
            4,
            MPI_DOUBLE,
            MPI_SUM,
            PETSC_COMM_WORLD) ==
            MPI_SUCCESS,
        "failed to reduce reactivation previous conserved totals");

    auto reactivation_control_state =
        wdp::
            AcceptedFixedTotalMolarRateWellControlState3D::
                minimum_bottom_hole_pressure(
                    minimum_bhp_pa);
    fdp::AcceptedPhysicalTimeClock3D
        reactivation_clock{
            0.0,
            reactivation_dt_seconds};
    wdp::
        FixedTotalMolarRatePhysicalTimestepDriverOptions3D
        reactivation_options;
    reactivation_options.adaptive
        .minimum_timestep_seconds =
        reactivation_dt_seconds;
    reactivation_options.adaptive
        .maximum_timestep_seconds =
        reactivation_dt_seconds;
    reactivation_options.adaptive
        .maximum_retries =
        1U;
    reactivation_options
        .minimum_bottom_hole_pressure_pa =
        minimum_bhp_pa;
    reactivation_options
        .minimum_bhp_release_rate_margin_mol_per_s =
        reactivation_rate_margin;
    reactivation_options
        .minimum_bhp_release_pressure_margin_pa =
        reactivation_pressure_margin_pa;

    std::optional<
        wdp::
            FixedTotalMolarRatePhysicalTimestepDriverReport3D>
        reactivation_report;
    error =
        wdp::
            advance_fixed_total_molar_rate_controlled_physical_timestep_3d(
                PETSC_COMM_WORLD,
                reactivation_reservoir_system.get(),
                &reactivation_source_context,
                reactivation_target_rate,
                reactivation_options,
                &reactivation_clock,
                &reactivation_control_state,
                &reactivation_report);

    require_collective(
        error == PETSC_SUCCESS &&
            reactivation_report.has_value() &&
            reactivation_report->accepted() &&
            reactivation_report
                    ->entry_control ==
                wdp::
                    FixedTotalMolarRatePhysicalTimestepControlMode3D::
                        minimum_bottom_hole_pressure &&
            reactivation_report
                    ->accepted_control ==
                wdp::
                    FixedTotalMolarRatePhysicalTimestepControlMode3D::
                        fixed_total_molar_rate &&
            reactivation_report
                    ->minimum_bhp_probe_total_molar_rate_mol_per_s
                    .has_value() &&
            reactivation_report
                    ->minimum_bhp_to_rate_reactivation_attempted &&
            reactivation_report
                    ->minimum_bhp_to_rate_reactivation_accepted &&
            reactivation_report
                    ->discarded_minimum_bhp_candidate
                    .has_value() &&
            !reactivation_report
                 ->discarded_rate_reactivation_candidate
                 .has_value() &&
            reactivation_report
                    ->accepted_solve
                    .has_value() &&
            !reactivation_report
                 ->accepted_fixed_bhp_solve
                 .has_value() &&
            reactivation_report
                    ->accepted_solve
                    ->bottom_hole_pressure_pa >=
                minimum_bhp_pa +
                    reactivation_pressure_margin_pa &&
            std::abs(
                reactivation_report
                    ->accepted_solve
                    ->total_molar_rate_residual_mol_per_s()) <=
                1.0e-8 *
                    std::max(
                        1.0,
                        std::abs(
                            reactivation_target_rate)) &&
            reactivation_control_state
                    .control ==
                wdp::
                    FixedTotalMolarRatePhysicalTimestepControlMode3D::
                        fixed_total_molar_rate &&
            reactivation_control_state
                    .bottom_hole_pressure_pa ==
                reactivation_report
                    ->accepted_solve
                    ->bottom_hole_pressure_pa &&
            reactivation_clock
                    .accepted_time_seconds() ==
                reactivation_dt_seconds &&
            reactivation_clock
                    .accepted_step_count() ==
                1U,
        "minimum-BHP feasibility probe did not reactivate fixed-total-molar-rate control through both hysteresis guards");

    near_collective(
        *reactivation_report
             ->minimum_bhp_probe_total_molar_rate_mol_per_s,
        reactivation_probe_total_molar_rate,
        2.0e-8,
        2.0e-10);

    std::uint64_t
        reactivation_final_authoritative_count = 0U;
    const auto reactivation_final_rate =
        global_fixed_bhp_well_production_rate(
            *reactivation_reservoir_system,
            reactivation_reservoir_system
                ->initial_state(),
            multi_context,
            reactivation_control_state
                .bottom_hole_pressure_pa,
            &reactivation_final_authoritative_count);
    const double
        reactivation_final_total_molar_rate =
            reactivation_final_rate[0] +
            reactivation_final_rate[1] +
            reactivation_final_rate[2];
    require_collective(
        reactivation_final_authoritative_count ==
                2U,
        "reactivated rate-control state lost owner-only two-connection aggregation");
    near_collective(
        reactivation_final_total_molar_rate,
        reactivation_target_rate,
        1.0e-8,
        1.0e-10);

    const auto reactivation_local_final_total =
        owned_conserved_totals(
            *reactivation_reservoir_system,
            reactivation_reservoir_system
                ->initial_state());
    std::array<double, 4>
        reactivation_global_final_total{};
    require_collective(
        MPI_Allreduce(
            reactivation_local_final_total.data(),
            reactivation_global_final_total.data(),
            4,
            MPI_DOUBLE,
            MPI_SUM,
            PETSC_COMM_WORLD) ==
            MPI_SUCCESS,
        "failed to reduce reactivated rate-control final conserved totals");
    for (std::size_t quantity = 0U;
         quantity < 4U;
         ++quantity) {
        near_collective(
            reactivation_global_final_total[
                quantity],
            reactivation_global_previous_total[
                quantity] -
                reactivation_dt_seconds *
                    reactivation_final_rate[
                        quantity],
            3.0e-7,
            quantity < 3U
                ? 3.0e-8
                : 3.0e-7);
    }

    Vec reactivation_probe_difference = nullptr;
    error =
        VecDuplicate(
            reactivation_reservoir_system
                ->initial_state(),
            &reactivation_probe_difference);
    if (error == PETSC_SUCCESS) {
        error =
            VecCopy(
                reactivation_reservoir_system
                    ->initial_state(),
                reactivation_probe_difference);
    }
    if (error == PETSC_SUCCESS) {
        error =
            VecAXPY(
                reactivation_probe_difference,
                PetscScalar{-1.0},
                reactivation_probe_state);
    }
    PetscReal reactivation_probe_difference_norm =
        -1.0;
    if (error == PETSC_SUCCESS) {
        error =
            VecNorm(
                reactivation_probe_difference,
                NORM_2,
                &reactivation_probe_difference_norm);
    }
    require_collective(
        error == PETSC_SUCCESS &&
            reactivation_probe_difference_norm >
                1.0e-9,
        "reactivation committed the discarded minimum-BHP feasibility probe");

    // The accepted reactivated state must persist: the next step enters rate
    // control directly and does not run another fixed-BHP feasibility probe.
    const double reactivated_second_time_n =
        reactivation_clock
            .accepted_time_seconds();
    std::optional<
        wdp::
            FixedTotalMolarRatePhysicalTimestepDriverReport3D>
        reactivated_second_report;
    error =
        wdp::
            advance_fixed_total_molar_rate_controlled_physical_timestep_3d(
                PETSC_COMM_WORLD,
                reactivation_reservoir_system.get(),
                &reactivation_source_context,
                reactivation_target_rate,
                reactivation_options,
                &reactivation_clock,
                &reactivation_control_state,
                &reactivated_second_report);
    require_collective(
        error == PETSC_SUCCESS &&
            reactivated_second_report
                .has_value() &&
            reactivated_second_report
                    ->accepted() &&
            reactivated_second_report
                    ->entry_control ==
                wdp::
                    FixedTotalMolarRatePhysicalTimestepControlMode3D::
                        fixed_total_molar_rate &&
            reactivated_second_report
                    ->accepted_control ==
                wdp::
                    FixedTotalMolarRatePhysicalTimestepControlMode3D::
                        fixed_total_molar_rate &&
            !reactivated_second_report
                 ->minimum_bhp_probe_total_molar_rate_mol_per_s
                 .has_value() &&
            !reactivated_second_report
                 ->minimum_bhp_to_rate_reactivation_attempted &&
            !reactivated_second_report
                 ->minimum_bhp_to_rate_reactivation_accepted &&
            reactivation_clock
                    .accepted_time_seconds() ==
                reactivated_second_time_n +
                    reactivation_dt_seconds &&
            reactivation_clock
                    .accepted_step_count() ==
                2U,
        "reactivated accepted rate control did not persist directly into the next timestep");

    // Rate-capacity deadband: p_min capacity does not clear target + Δq, so
    // no augmented rate candidate may be created.
    auto capacity_deadband_source_context =
        wdp::
            FixedTotalMolarRateWellSourceEvaluatorContext3D::
                create(
                    multi_context,
                    minimum_bhp_pa);
    auto capacity_deadband_system =
        make_fixed_total_molar_rate_reservoir_system(
            rank,
            schedule,
            partition,
            bridge,
            pattern,
            audit,
            &capacity_deadband_source_context,
            &provider);
    require_collective(
        capacity_deadband_system
                    ->set_trial_timestep_seconds(
                        reactivation_dt_seconds) ==
                PETSC_SUCCESS,
        "failed to prepare rate-capacity deadband system");

    auto capacity_deadband_control =
        wdp::
            AcceptedFixedTotalMolarRateWellControlState3D::
                minimum_bottom_hole_pressure(
                    minimum_bhp_pa);
    fdp::AcceptedPhysicalTimeClock3D
        capacity_deadband_clock{
            0.0,
            reactivation_dt_seconds};
    auto capacity_deadband_options =
        reactivation_options;
    const double capacity_deadband_target =
        0.95 *
        reactivation_probe_total_molar_rate;
    capacity_deadband_options
        .minimum_bhp_release_rate_margin_mol_per_s =
        0.10 *
        reactivation_probe_total_molar_rate;

    std::optional<
        wdp::
            FixedTotalMolarRatePhysicalTimestepDriverReport3D>
        capacity_deadband_report;
    error =
        wdp::
            advance_fixed_total_molar_rate_controlled_physical_timestep_3d(
                PETSC_COMM_WORLD,
                capacity_deadband_system.get(),
                &capacity_deadband_source_context,
                capacity_deadband_target,
                capacity_deadband_options,
                &capacity_deadband_clock,
                &capacity_deadband_control,
                &capacity_deadband_report);
    require_collective(
        error == PETSC_SUCCESS &&
            capacity_deadband_report
                .has_value() &&
            capacity_deadband_report
                    ->accepted() &&
            capacity_deadband_report
                    ->accepted_control ==
                wdp::
                    FixedTotalMolarRatePhysicalTimestepControlMode3D::
                        minimum_bottom_hole_pressure &&
            capacity_deadband_report
                    ->minimum_bhp_probe_total_molar_rate_mol_per_s
                    .has_value() &&
            !capacity_deadband_report
                 ->minimum_bhp_to_rate_reactivation_attempted &&
            !capacity_deadband_report
                 ->minimum_bhp_to_rate_reactivation_accepted &&
            !capacity_deadband_report
                 ->discarded_minimum_bhp_candidate
                 .has_value() &&
            !capacity_deadband_report
                 ->discarded_rate_reactivation_candidate
                 .has_value() &&
            !capacity_deadband_report
                 ->accepted_solve
                 .has_value() &&
            capacity_deadband_report
                    ->accepted_fixed_bhp_solve
                    .has_value() &&
            capacity_deadband_control
                    .control ==
                wdp::
                    FixedTotalMolarRatePhysicalTimestepControlMode3D::
                        minimum_bottom_hole_pressure,
        "minimum-BHP rate-capacity hysteresis gate reactivated rate control inside the deadband");

    // Pressure deadband: capacity clears Δq and the augmented rate root
    // converges, but an intentionally wide Δp guard rejects that candidate and
    // keeps the already converged fixed-BHP probe as the accepted state.
    auto pressure_deadband_source_context =
        wdp::
            FixedTotalMolarRateWellSourceEvaluatorContext3D::
                create(
                    multi_context,
                    minimum_bhp_pa);
    auto pressure_deadband_system =
        make_fixed_total_molar_rate_reservoir_system(
            rank,
            schedule,
            partition,
            bridge,
            pattern,
            audit,
            &pressure_deadband_source_context,
            &provider);
    require_collective(
        pressure_deadband_system
                    ->set_trial_timestep_seconds(
                        reactivation_dt_seconds) ==
                PETSC_SUCCESS,
        "failed to prepare pressure-deadband reactivation system");

    auto pressure_deadband_control =
        wdp::
            AcceptedFixedTotalMolarRateWellControlState3D::
                minimum_bottom_hole_pressure(
                    minimum_bhp_pa);
    fdp::AcceptedPhysicalTimeClock3D
        pressure_deadband_clock{
            0.0,
            reactivation_dt_seconds};
    auto pressure_deadband_options =
        reactivation_options;
    pressure_deadband_options
        .minimum_bhp_release_pressure_margin_pa =
        100.0;

    std::optional<
        wdp::
            FixedTotalMolarRatePhysicalTimestepDriverReport3D>
        pressure_deadband_report;
    error =
        wdp::
            advance_fixed_total_molar_rate_controlled_physical_timestep_3d(
                PETSC_COMM_WORLD,
                pressure_deadband_system.get(),
                &pressure_deadband_source_context,
                reactivation_target_rate,
                pressure_deadband_options,
                &pressure_deadband_clock,
                &pressure_deadband_control,
                &pressure_deadband_report);
    require_collective(
        error == PETSC_SUCCESS &&
            pressure_deadband_report
                .has_value() &&
            pressure_deadband_report
                    ->accepted() &&
            pressure_deadband_report
                    ->accepted_control ==
                wdp::
                    FixedTotalMolarRatePhysicalTimestepControlMode3D::
                        minimum_bottom_hole_pressure &&
            pressure_deadband_report
                    ->minimum_bhp_to_rate_reactivation_attempted &&
            !pressure_deadband_report
                 ->minimum_bhp_to_rate_reactivation_accepted &&
            pressure_deadband_report
                    ->discarded_rate_reactivation_candidate
                    .has_value() &&
            !pressure_deadband_report
                 ->discarded_minimum_bhp_candidate
                 .has_value() &&
            !pressure_deadband_report
                 ->accepted_solve
                 .has_value() &&
            pressure_deadband_report
                    ->accepted_fixed_bhp_solve
                    .has_value() &&
            pressure_deadband_report
                    ->discarded_rate_reactivation_candidate
                    ->bottom_hole_pressure_pa <
                minimum_bhp_pa +
                    *pressure_deadband_options
                         .minimum_bhp_release_pressure_margin_pa &&
            pressure_deadband_control
                    .control ==
                wdp::
                    FixedTotalMolarRatePhysicalTimestepControlMode3D::
                        minimum_bottom_hole_pressure &&
            pressure_deadband_control
                    .bottom_hole_pressure_pa ==
                minimum_bhp_pa,
        "minimum-BHP pressure hysteresis gate accepted a rate candidate inside the release deadband");

    require_collective(
        VecDestroy(
            &reactivation_probe_difference) ==
                PETSC_SUCCESS &&
            VecDestroy(
                &reactivation_probe_state) ==
                PETSC_SUCCESS,
        "minimum-BHP reactivation regression cleanup failed");

    Vec discarded_rate_difference = nullptr;
    error =
        VecDuplicate(
            switch_reservoir_system
                ->initial_state(),
            &discarded_rate_difference);
    if (error == PETSC_SUCCESS) {
        error =
            VecCopy(
                switch_reservoir_system
                    ->initial_state(),
                discarded_rate_difference);
    }
    if (error == PETSC_SUCCESS) {
        error =
            VecAXPY(
                discarded_rate_difference,
                PetscScalar{-1.0},
                reservoir_solution);
    }
    PetscReal discarded_rate_difference_norm =
        -1.0;
    if (error == PETSC_SUCCESS) {
        error =
            VecNorm(
                discarded_rate_difference,
                NORM_2,
                &discarded_rate_difference_norm);
    }
    require_collective(
        error == PETSC_SUCCESS &&
            discarded_rate_difference_norm >
                1.0e-9,
        "minimum-BHP switch appears to have committed the discarded rate-control reservoir candidate");

    require_collective(
        VecDestroy(
            &discarded_rate_difference) ==
            PETSC_SUCCESS,
        "minimum-BHP switching regression cleanup failed");

    // Controlled post-SNES phase-transition restart: cell30 starts 2P,
    // cell60 stays frozen 3P. The first selected rate candidate is scanned,
    // discarded, rebuilt to 3P by stable identity, and the same 1 s physical
    // timestep is solved again from the entry accepted rate-control state.
    const auto aqueous_identity =
        mixed_physical_phase_identity(
            "aqueous");
    const auto hc0_identity =
        mixed_physical_phase_identity(
            "hydrocarbon-0");
    const auto hc1_identity =
        mixed_physical_phase_identity(
            "hydrocarbon-1");
    const wdp::
        FixedBhpPhaseIdentityInjectionEnthalpy3D
        controlled_phase_enthalpies{
            "fixture/controlled-transition-enthalpy/v1",
            {
                {aqueous_identity, 1000.0},
                {hc0_identity, 2000.0},
                {hc1_identity, 3000.0}}};

    const auto make_transition_well =
        [&](double bhp_pa) {
            auto reboundable30 =
                wdp::
                    FixedBhpPeacemanWellSourceEvaluatorContext3D::
                        create_phase_identity_bound(
                            mesh::GlobalEntityId{
                                UINT64_C(30)},
                            well::
                                make_peaceman_well_index_3d(
                                    {10.0, 10.0, 5.0},
                                    {
                                        1.0e-8,
                                        1.0e-8,
                                        1.0e-8},
                                    well::
                                        AxisAlignedWellDirection3D::z,
                                    0.10,
                                    0.0),
                            bhp_pa,
                            controlled_phase_enthalpies,
                            phase_identity_map_for_stable(
                                UINT64_C(30)),
                            "fixture/controlled-transition-cell30/v1");
            auto frozen60 =
                wdp::
                    FixedBhpPeacemanWellSourceEvaluatorContext3D::
                        create(
                            mesh::GlobalEntityId{
                                UINT64_C(60)},
                            well::
                                make_peaceman_well_index_3d(
                                    {12.0, 9.0, 4.0},
                                    {
                                        2.0e-8,
                                        1.5e-8,
                                        2.5e-8},
                                    well::
                                        AxisAlignedWellDirection3D::z,
                                    0.12,
                                    0.1),
                            bhp_pa,
                            wd::
                                FixedBhpInjectionEnthalpy3D{
                                    "fixture/controlled-transition-cell60-enthalpy/v1",
                                    {
                                        1100.0,
                                        2200.0,
                                        3300.0}},
                            "fixture/controlled-transition-cell60/v1");
            return wdp::
                FixedBhpMultiConnectionWellSourceEvaluatorContext3D::
                    create(
                        "fixture/controlled-transition-well/v1",
                        {
                            frozen60,
                            reboundable30});
        };

    const flow::FrozenActivePhaseIdentityMap
        controlled_transition_three_phase_map{
            {
                aqueous_identity,
                hc0_identity,
                hc1_identity}};

    // Construct a transition-specific control target from the 2P chart, then
    // prove directly with independent augmented solves that the same target
    // has a valid root on both frozen phase charts. Transactional restart does
    // not require those roots to occur in the same BHP interval.
    auto transition_fixed_bhp_rate =
        [&](double bhp_pa,
            bool cell30_is_three_phase)
            -> std::optional<double> {
            auto reference_source =
                wdp::
                    FixedTotalMolarRateWellSourceEvaluatorContext3D::
                        create(
                            make_transition_well(
                                bhp_pa),
                            bhp_pa);
            if (cell30_is_three_phase) {
                reference_source
                    .rebind_connection(
                        mesh::GlobalEntityId{
                            UINT64_C(30)},
                        controlled_transition_three_phase_map);
            }
            reference_source
                .begin_fixed_bhp_reservoir_evaluation(
                    bhp_pa);
            auto reference_system =
                make_fixed_total_molar_rate_reservoir_system(
                    rank,
                    schedule,
                    partition,
                    bridge,
                    pattern,
                    audit,
                    &reference_source,
                    &provider,
                    cell30_is_three_phase);

            Vec reference_state = nullptr;
            std::optional<
                fdp::
                    VariableCardinalityNaturalVariableSnesSolveReport3D>
                reference_solve;
            const PetscErrorCode reference_error =
                reference_system
                    ->solve(
                        &reference_state,
                        &reference_solve);
            const bool solved =
                reference_error ==
                    PETSC_SUCCESS &&
                reference_state !=
                    nullptr &&
                reference_solve
                    .has_value() &&
                static_cast<int>(
                    reference_solve
                        ->converged_reason) >
                    0;
            int local_solved =
                solved ? 1 : 0;
            int global_solved = 0;
            if (MPI_Allreduce(
                    &local_solved,
                    &global_solved,
                    1,
                    MPI_INT,
                    MPI_MIN,
                    PETSC_COMM_WORLD) !=
                MPI_SUCCESS) {
                if (reference_state !=
                    nullptr) {
                    (void)VecDestroy(
                        &reference_state);
                }
                throw std::runtime_error(
                    "failed to reduce transition fixed-BHP solve status");
            }
            if (global_solved == 0) {
                if (reference_state !=
                    nullptr) {
                    (void)VecDestroy(
                        &reference_state);
                }
                return std::nullopt;
            }

            std::uint64_t
                authoritative_count = 0U;
            const auto rate =
                global_fixed_bhp_well_production_rate(
                    *reference_system,
                    reference_state,
                    reference_source.well(),
                    bhp_pa,
                    &authoritative_count);
            const double total_molar_rate =
                rate[0] +
                rate[1] +
                rate[2];
            const PetscErrorCode destroy_error =
                VecDestroy(
                    &reference_state);
            if (destroy_error !=
                PETSC_SUCCESS ||
                authoritative_count !=
                    2U ||
                !std::isfinite(
                    total_molar_rate) ||
                !(total_molar_rate > 0.0)) {
                return std::nullopt;
            }
            return total_molar_rate;
        };

    auto transition_rate_root =
        [&](double target,
            bool cell30_is_three_phase)
            -> std::optional<double> {
            auto oracle_source =
                wdp::
                    FixedTotalMolarRateWellSourceEvaluatorContext3D::
                        create(
                            make_transition_well(
                                initial_bhp_pa),
                            initial_bhp_pa);
            if (cell30_is_three_phase) {
                oracle_source
                    .rebind_connection(
                        mesh::GlobalEntityId{
                            UINT64_C(30)},
                        controlled_transition_three_phase_map);
            }
            auto oracle_system =
                make_fixed_total_molar_rate_reservoir_system(
                    rank,
                    schedule,
                    partition,
                    bridge,
                    pattern,
                    audit,
                    &oracle_source,
                    &provider,
                    cell30_is_three_phase);
            std::unique_ptr<
                wdp::
                    FixedTotalMolarRateWellControlSystem3D>
                oracle_control;
            PetscErrorCode oracle_error =
                wdp::
                    FixedTotalMolarRateWellControlSystem3D::
                        create(
                            PETSC_COMM_WORLD,
                            oracle_system.get(),
                            &oracle_source,
                            initial_bhp_pa,
                            target,
                            &oracle_control);

            Vec oracle_state = nullptr;
            std::optional<
                wdp::
                    FixedTotalMolarRateWellControlSolveReport3D>
                oracle_report;
            if (oracle_error ==
                PETSC_SUCCESS) {
                oracle_error =
                    oracle_control
                        ->solve(
                            &oracle_state,
                            &oracle_report);
            }

            bool root_valid =
                oracle_error ==
                    PETSC_SUCCESS &&
                oracle_state !=
                    nullptr &&
                oracle_report
                    .has_value() &&
                oracle_report
                    ->converged() &&
                std::isfinite(
                    oracle_report
                        ->bottom_hole_pressure_pa) &&
                oracle_report
                        ->bottom_hole_pressure_pa >
                    0.0 &&
                std::abs(
                    oracle_report
                        ->total_molar_rate_residual_mol_per_s()) <=
                    1.0e-8 *
                        std::max(
                            1.0,
                            std::abs(
                                target));
            int local_root_valid =
                root_valid ? 1 : 0;
            int global_root_valid = 0;
            if (MPI_Allreduce(
                    &local_root_valid,
                    &global_root_valid,
                    1,
                    MPI_INT,
                    MPI_MIN,
                    PETSC_COMM_WORLD) !=
                MPI_SUCCESS) {
                if (oracle_state !=
                    nullptr) {
                    (void)VecDestroy(
                        &oracle_state);
                }
                throw std::runtime_error(
                    "failed to reduce transition rate-root status");
            }

            double root_bhp = 0.0;
            if (global_root_valid != 0) {
                root_bhp =
                    oracle_report
                        ->bottom_hole_pressure_pa;
            }
            if (oracle_state !=
                nullptr) {
                if (VecDestroy(
                        &oracle_state) !=
                    PETSC_SUCCESS) {
                    throw std::runtime_error(
                        "transition rate-root cleanup failed");
                }
            }
            if (global_root_valid == 0) {
                return std::nullopt;
            }
            return root_bhp;
        };

    const auto
        transition_two_phase_rate_at_entry =
            transition_fixed_bhp_rate(
                initial_bhp_pa,
                false);
    const auto
        transition_two_phase_rate_at_reference =
            transition_fixed_bhp_rate(
                target_reference_bhp_pa,
                false);
    require_collective(
        transition_two_phase_rate_at_entry
                .has_value() &&
            transition_two_phase_rate_at_reference
                .has_value() &&
            target_reference_bhp_pa <
                initial_bhp_pa &&
            *transition_two_phase_rate_at_reference >
                *transition_two_phase_rate_at_entry,
        "2P transition chart did not provide a lower-BHP production interval");

    constexpr std::array<double, 6>
        transition_target_fractions{
            0.01,
            0.025,
            0.05,
            0.10,
            0.25,
            0.50};

    std::optional<double>
        transition_target_rate_candidate;
    std::optional<double>
        transition_two_phase_rate_root_bhp;
    std::optional<double>
        transition_three_phase_rate_root_bhp;

    const double transition_rate_span =
        *transition_two_phase_rate_at_reference -
        *transition_two_phase_rate_at_entry;
    for (double fraction :
         transition_target_fractions) {
        const double candidate =
            *transition_two_phase_rate_at_entry +
            fraction *
                transition_rate_span;
        const auto two_phase_root =
            transition_rate_root(
                candidate,
                false);
        if (!two_phase_root.has_value() ||
            !(*two_phase_root <
              initial_bhp_pa)) {
            continue;
        }
        const auto three_phase_root =
            transition_rate_root(
                candidate,
                true);
        if (!three_phase_root.has_value()) {
            continue;
        }

        transition_target_rate_candidate =
            candidate;
        transition_two_phase_rate_root_bhp =
            two_phase_root;
        transition_three_phase_rate_root_bhp =
            three_phase_root;
        break;
    }

    require_collective(
        transition_target_rate_candidate
                .has_value() &&
            transition_two_phase_rate_root_bhp
                .has_value() &&
            transition_three_phase_rate_root_bhp
                .has_value(),
        "no transition-specific target has valid frozen 2P and 3P augmented rate roots");

    const double transition_target_rate =
        *transition_target_rate_candidate;
    const double
        transition_two_phase_rate_root_bhp_pa =
            *transition_two_phase_rate_root_bhp;

    auto transition_rate_source =
        wdp::
            FixedTotalMolarRateWellSourceEvaluatorContext3D::
                create(
                    make_transition_well(
                        initial_bhp_pa),
                    initial_bhp_pa);
    auto transition_rate_system =
        make_fixed_total_molar_rate_reservoir_system(
            rank,
            schedule,
            partition,
            bridge,
            pattern,
            audit,
            &transition_rate_source,
            &provider);

    auto transition_rate_control =
        wdp::
            AcceptedFixedTotalMolarRateWellControlState3D::
                fixed_total_molar_rate(
                    initial_bhp_pa);
    fdp::AcceptedPhysicalTimeClock3D
        transition_rate_clock{
            0.0,
            1.0};
    ControlledWellPhaseTransitionFixture
        transition_rate_fixture{
            rank,
            &schedule,
            &partition,
            &bridge,
            &pattern,
            audit,
            &provider,
            &transition_rate_source,
            &transition_rate_control,
            &transition_rate_clock,
            wdp::
                FixedTotalMolarRatePhysicalTimestepControlMode3D::
                    fixed_total_molar_rate,
            initial_bhp_pa,
            0U,
            0U,
            true,
            false,
            {},
            true,
            {}};

    wdp::
        FixedTotalMolarRatePhysicalTimestepDriverOptions3D
        transition_rate_options;
    transition_rate_options.adaptive
        .minimum_timestep_seconds =
        1.0;
    transition_rate_options.adaptive
        .maximum_timestep_seconds =
        1.0;
    transition_rate_options.adaptive
        .maximum_retries =
        1U;
    transition_rate_options.phase_transition
        .max_transition_restarts =
        2U;

    std::optional<
        wdp::
            FixedTotalMolarRatePhysicalTimestepDriverReport3D>
        transition_rate_report;
    error =
        wdp::
            advance_fixed_total_molar_rate_controlled_physical_timestep_with_phase_transitions_3d(
                PETSC_COMM_WORLD,
                &transition_rate_system,
                &transition_rate_source,
                transition_target_rate,
                {
                    &controlled_well_transition_scan,
                    &transition_rate_fixture,
                    &controlled_well_transition_rebuild,
                    &transition_rate_fixture},
                transition_rate_options,
                &transition_rate_clock,
                &transition_rate_control,
                &transition_rate_report);

    const auto* rebound_rate_cell30 =
        transition_rate_source
            .well()
            .find_connection(
                mesh::GlobalEntityId{
                    UINT64_C(30)});
    const bool transition_rate_report_present =
        transition_rate_report.has_value();
    const bool transition_rate_accepted =
        transition_rate_report_present &&
        transition_rate_report->accepted();
    const bool transition_rate_record_present =
        transition_rate_report_present &&
        transition_rate_report
            ->accepted_record
            .has_value();
    const bool transition_rate_solve_present =
        transition_rate_report_present &&
        transition_rate_report
            ->accepted_solve
            .has_value();
    const bool transition_rate_reservoir_report_present =
        transition_rate_solve_present &&
        transition_rate_report
            ->accepted_solve
            ->reservoir_nonlinear_solve
            .has_value();
    const std::string transition_rate_restart_diagnostic =
        std::string{"error="} +
        std::to_string(
            static_cast<int>(error)) +
        ", report=" +
        std::to_string(
            transition_rate_report_present ? 1 : 0) +
        ", accepted=" +
        std::to_string(
            transition_rate_accepted ? 1 : 0) +
        ", adaptive_outcome=" +
        std::to_string(
            transition_rate_report_present
                ? static_cast<int>(
                      transition_rate_report
                          ->adaptive
                          .outcome)
                : -999) +
        ", adaptive_retries=" +
        std::to_string(
            transition_rate_report_present
                ? transition_rate_report
                      ->adaptive
                      .retries
                : 999U) +
        ", adaptive_attempts=" +
        std::to_string(
            transition_rate_report_present
                ? transition_rate_report
                      ->adaptive
                      .attempts
                      .size()
                : 999U) +
        ", attempt_outcome=" +
        std::to_string(
            transition_rate_report_present &&
                    !transition_rate_report
                         ->adaptive
                         .attempts
                         .empty()
                ? static_cast<int>(
                      transition_rate_report
                          ->adaptive
                          .attempts
                          .back()
                          .result
                          .outcome)
                : -999) +
        ", attempt_decision=" +
        std::to_string(
            transition_rate_report_present &&
                    !transition_rate_report
                         ->adaptive
                         .attempts
                         .empty()
                ? static_cast<int>(
                      transition_rate_report
                          ->adaptive
                          .attempts
                          .back()
                          .decision)
                : -999) +
        ", attempt_nonlinear_iterations=" +
        std::to_string(
            transition_rate_report_present &&
                    !transition_rate_report
                         ->adaptive
                         .attempts
                         .empty()
                ? transition_rate_report
                      ->adaptive
                      .attempts
                      .back()
                      .result
                      .nonlinear_iterations
                : -999) +
        ", attempt_function_domain_errors=" +
        std::to_string(
            transition_rate_report_present &&
                    !transition_rate_report
                         ->adaptive
                         .attempts
                         .empty()
                ? transition_rate_report
                      ->adaptive
                      .attempts
                      .back()
                      .result
                      .function_domain_errors
                : -999) +
        ", attempt_jacobian_domain_errors=" +
        std::to_string(
            transition_rate_report_present &&
                    !transition_rate_report
                         ->adaptive
                         .attempts
                         .empty()
                ? transition_rate_report
                      ->adaptive
                      .attempts
                      .back()
                      .result
                      .jacobian_domain_errors
                : -999) +
        ", attempt_line_search_changes=" +
        std::to_string(
            transition_rate_report_present &&
                    !transition_rate_report
                         ->adaptive
                         .attempts
                         .empty()
                ? transition_rate_report
                      ->adaptive
                      .attempts
                      .back()
                      .result
                      .line_search_direction_changes
                : -999) +
        ", attempt_transition_restarts=" +
        std::to_string(
            transition_rate_report_present &&
                    !transition_rate_report
                         ->adaptive
                         .attempts
                         .empty()
                ? transition_rate_report
                      ->adaptive
                      .attempts
                      .back()
                      .result
                      .phase_transition_restarts
                : 999U) +
        ", report_restarts=" +
        std::to_string(
            transition_rate_report_present
                ? transition_rate_report
                      ->phase_transition_restarts
                : 999U) +
        ", record=" +
        std::to_string(
            transition_rate_record_present ? 1 : 0) +
        ", record_restarts=" +
        std::to_string(
            transition_rate_record_present
                ? transition_rate_report
                      ->accepted_record
                      ->phase_transition_restarts
                : 999U) +
        ", scans=" +
        std::to_string(
            transition_rate_fixture.scans) +
        ", rebuilds=" +
        std::to_string(
            transition_rate_fixture.rebuild_calls) +
        ", anchor_rebuild=" +
        std::to_string(
            transition_rate_fixture
                    .accepted_anchor_unchanged_at_rebuild
                ? 1
                : 0) +
        ", anchor_scan=" +
        std::to_string(
            transition_rate_fixture
                    .accepted_anchor_unchanged_at_scan
                ? 1
                : 0) +
        ", previous_totals=" +
        std::to_string(
            transition_rate_fixture
                    .rebuilt_previous_totals_captured
                ? 1
                : 0) +
        ", cell30_phases=" +
        std::to_string(
            transition_rate_system
                ->numbering()
                .cell(
                    mesh::LocalIndex{2U})
                .phase_count) +
        ", cell60_phases=" +
        std::to_string(
            transition_rate_system
                ->numbering()
                .cell(
                    mesh::LocalIndex{5U})
                .phase_count) +
        ", control=" +
        std::to_string(
            static_cast<int>(
                transition_rate_control
                    .control)) +
        ", solve=" +
        std::to_string(
            transition_rate_solve_present ? 1 : 0) +
        ", reservoir_report=" +
        std::to_string(
            transition_rate_reservoir_report_present
                ? 1
                : 0) +
        ", clock=" +
        std::to_string(
            transition_rate_clock
                .accepted_time_seconds()) +
        ", steps=" +
        std::to_string(
            transition_rate_clock
                .accepted_step_count());

    require_collective(
        error == PETSC_SUCCESS &&
            transition_rate_report_present &&
            transition_rate_accepted &&
            transition_rate_report
                    ->phase_transition_restarts ==
                1U &&
            transition_rate_record_present &&
            transition_rate_report
                    ->accepted_record
                    ->phase_transition_restarts ==
                1U &&
            transition_rate_fixture.scans ==
                2U &&
            transition_rate_fixture.rebuild_calls ==
                1U &&
            transition_rate_fixture
                .accepted_anchor_unchanged_at_rebuild &&
            transition_rate_fixture
                .accepted_anchor_unchanged_at_scan &&
            transition_rate_fixture
                .rebuilt_previous_totals_captured &&
            transition_rate_system
                    ->numbering()
                    .cell(
                        mesh::LocalIndex{2U})
                    .phase_count ==
                3U &&
            transition_rate_system
                    ->numbering()
                    .cell(
                        mesh::LocalIndex{5U})
                    .phase_count ==
                3U &&
            transition_rate_control
                    .control ==
                wdp::
                    FixedTotalMolarRatePhysicalTimestepControlMode3D::
                        fixed_total_molar_rate &&
            transition_rate_solve_present &&
            transition_rate_reservoir_report_present &&
            transition_rate_report
                    ->accepted_solve
                    ->global_scalar_count ==
                transition_rate_system
                    ->numbering()
                    .petsc_global_scalar_count() +
                    1 &&
            transition_rate_clock
                    .accepted_time_seconds() ==
                1.0 &&
            transition_rate_clock
                    .accepted_step_count() ==
                1U &&
            rebound_rate_cell30 !=
                nullptr &&
            rebound_rate_cell30
                    ->active_phase_identities()
                    .has_value() &&
            rebound_rate_cell30
                    ->active_phase_identities()
                    ->phase_count() ==
                3U,
        std::string{
            "rate-controlled post-SNES transition invariant failed: "} +
            transition_rate_restart_diagnostic);

    std::uint64_t
        transition_rate_authoritative_count = 0U;
    const auto transition_rate_well_rate =
        global_fixed_bhp_well_production_rate(
            *transition_rate_system,
            transition_rate_system
                ->initial_state(),
            transition_rate_source.well(),
            transition_rate_control
                .bottom_hole_pressure_pa,
            &transition_rate_authoritative_count);
    const double
        transition_rate_total_molar =
            transition_rate_well_rate[0] +
            transition_rate_well_rate[1] +
            transition_rate_well_rate[2];
    require_collective(
        transition_rate_authoritative_count ==
            2U,
        "rate-controlled transition lost owner-only two-connection aggregation");
    const double
        transition_rate_residual_tolerance =
            1.0e-8 *
            std::max(
                1.0,
                std::abs(
                    transition_target_rate));
    const double
        transition_report_rate_residual =
            transition_rate_report
                ->accepted_solve
                ->total_molar_rate_residual_mol_per_s();
    require_collective(
        std::abs(
            transition_report_rate_residual) <=
            transition_rate_residual_tolerance,
        std::string{
            "transitioned augmented rate solve violated control residual: residual="} +
            std::to_string(
                transition_report_rate_residual) +
            ", target=" +
            std::to_string(
                transition_target_rate));
    require_collective(
        std::abs(
            transition_rate_total_molar -
            transition_target_rate) <=
            transition_rate_residual_tolerance,
        std::string{
            "transitioned committed whole-well rate drifted from control root: actual="} +
            std::to_string(
                transition_rate_total_molar) +
            ", target=" +
            std::to_string(
                transition_target_rate) +
            ", accepted_bhp=" +
            std::to_string(
                transition_rate_control
                    .bottom_hole_pressure_pa));

    const auto transition_rate_local_final =
        owned_conserved_totals(
            *transition_rate_system,
            transition_rate_system
                ->initial_state());
    std::array<double, 4>
        transition_rate_global_final{};
    require_collective(
        MPI_Allreduce(
            transition_rate_local_final.data(),
            transition_rate_global_final.data(),
            4,
            MPI_DOUBLE,
            MPI_SUM,
            PETSC_COMM_WORLD) ==
            MPI_SUCCESS,
        "failed to reduce transitioned rate-control final totals");
    for (std::size_t quantity = 0U;
         quantity < 4U;
         ++quantity) {
        near_collective(
            transition_rate_global_final[
                quantity],
            transition_rate_fixture
                    .rebuilt_previous_totals[
                        quantity] -
                transition_rate_well_rate[
                    quantity],
            4.0e-7,
            quantity < 3U
                ? 4.0e-8
                : 4.0e-7);
    }

    bool transition_rate_history_matches =
        false;
    require_collective(
        transition_rate_system
                ->accepted_history_matches_state(
                    transition_rate_system
                        ->initial_state(),
                    &transition_rate_history_matches) ==
                PETSC_SUCCESS &&
            transition_rate_history_matches,
        "transitioned rate-control accepted history did not commit exactly once");

    // Cross the two state machines in one physical timestep. Entry control is
    // accepted RATE, but the first 2P rate candidate is known to violate the
    // already validated p_min and therefore selects a tentative fixed-BHP
    // candidate. That selected BHP candidate triggers 2P -> 3P. The phase
    // restart must discard the complete trial control generation and restart
    // the rebuilt system from the accepted RATE/BHP anchor, not from p_min.
    auto transition_switch_source =
        wdp::
            FixedTotalMolarRateWellSourceEvaluatorContext3D::
                create(
                    make_transition_well(
                        initial_bhp_pa),
                    initial_bhp_pa);
    auto transition_switch_system =
        make_fixed_total_molar_rate_reservoir_system(
            rank,
            schedule,
            partition,
            bridge,
            pattern,
            audit,
            &transition_switch_source,
            &provider);
    auto transition_switch_control =
        wdp::
            AcceptedFixedTotalMolarRateWellControlState3D::
                fixed_total_molar_rate(
                    initial_bhp_pa);
    fdp::AcceptedPhysicalTimeClock3D
        transition_switch_clock{
            0.0,
            1.0};
    ControlledWellPhaseTransitionFixture
        transition_switch_fixture{
            rank,
            &schedule,
            &partition,
            &bridge,
            &pattern,
            audit,
            &provider,
            &transition_switch_source,
            &transition_switch_control,
            &transition_switch_clock,
            wdp::
                FixedTotalMolarRatePhysicalTimestepControlMode3D::
                    fixed_total_molar_rate,
            initial_bhp_pa,
            0U,
            0U,
            true,
            false,
            {},
            true,
            {}};

    const double
        transition_switch_minimum_bhp_pa =
            0.5 *
            (transition_two_phase_rate_root_bhp_pa +
             initial_bhp_pa);
    require_collective(
        transition_switch_minimum_bhp_pa >
                transition_two_phase_rate_root_bhp_pa &&
            transition_switch_minimum_bhp_pa <
                initial_bhp_pa,
        "transactional switch p_min is not strictly above the feasible 2P rate root");

    auto transition_switch_options =
        transition_rate_options;
    transition_switch_options
        .minimum_bottom_hole_pressure_pa =
        transition_switch_minimum_bhp_pa;

    std::optional<
        wdp::
            FixedTotalMolarRatePhysicalTimestepDriverReport3D>
        transition_switch_report;
    error =
        wdp::
            advance_fixed_total_molar_rate_controlled_physical_timestep_with_phase_transitions_3d(
                PETSC_COMM_WORLD,
                &transition_switch_system,
                &transition_switch_source,
                transition_target_rate,
                {
                    &controlled_well_transition_scan,
                    &transition_switch_fixture,
                    &controlled_well_transition_rebuild,
                    &transition_switch_fixture},
                transition_switch_options,
                &transition_switch_clock,
                &transition_switch_control,
                &transition_switch_report);

    const wdp::
        FixedTotalMolarRateWellControlSolveReport3D*
        restarted_rate_candidate =
            nullptr;
    if (transition_switch_report
            .has_value()) {
        if (transition_switch_report
                ->accepted_solve
                .has_value()) {
            restarted_rate_candidate =
                &*transition_switch_report
                      ->accepted_solve;
        } else if (
            transition_switch_report
                ->discarded_rate_control_candidate
                .has_value()) {
            restarted_rate_candidate =
                &*transition_switch_report
                      ->discarded_rate_control_candidate;
        }
    }

    const auto* rebound_switch_cell30 =
        transition_switch_source
            .well()
            .find_connection(
                mesh::GlobalEntityId{
                    UINT64_C(30)});
    require_collective(
        error == PETSC_SUCCESS &&
            transition_switch_report
                .has_value() &&
            transition_switch_report
                    ->accepted() &&
            transition_switch_report
                    ->phase_transition_restarts ==
                1U &&
            transition_switch_report
                    ->accepted_record
                    ->phase_transition_restarts ==
                1U &&
            transition_switch_fixture.scans ==
                2U &&
            transition_switch_fixture
                    .scanned_source_bottom_hole_pressure_pa
                    .size() ==
                2U &&
            transition_switch_fixture
                    .scanned_source_bottom_hole_pressure_pa[
                        0U] ==
                transition_switch_minimum_bhp_pa &&
            transition_switch_fixture
                .accepted_anchor_unchanged_at_rebuild &&
            transition_switch_fixture
                .accepted_anchor_unchanged_at_scan &&
            restarted_rate_candidate !=
                nullptr &&
            restarted_rate_candidate
                    ->initial_bottom_hole_pressure_pa ==
                initial_bhp_pa &&
            transition_switch_clock
                    .accepted_time_seconds() ==
                1.0 &&
            transition_switch_clock
                    .accepted_step_count() ==
                1U &&
            rebound_switch_cell30 !=
                nullptr &&
            rebound_switch_cell30
                    ->active_phase_identities()
                    .has_value() &&
            rebound_switch_cell30
                    ->active_phase_identities()
                    ->phase_count() ==
                3U,
        "rate-to-BHP trial control generation leaked across the phase-transition restart");

    std::uint64_t
        transition_switch_authoritative_count =
            0U;
    const auto transition_switch_well_rate =
        global_fixed_bhp_well_production_rate(
            *transition_switch_system,
            transition_switch_system
                ->initial_state(),
            transition_switch_source.well(),
            transition_switch_control
                .bottom_hole_pressure_pa,
            &transition_switch_authoritative_count);
    require_collective(
        transition_switch_authoritative_count ==
            2U,
        "control/phase transactional restart lost owner-only two-connection aggregation");

    const double
        transition_switch_total_molar_rate =
            transition_switch_well_rate[0] +
            transition_switch_well_rate[1] +
            transition_switch_well_rate[2];
    if (transition_switch_control.control ==
        wdp::
            FixedTotalMolarRatePhysicalTimestepControlMode3D::
                fixed_total_molar_rate) {
        require_collective(
            transition_switch_report
                    ->accepted_solve
                    .has_value() &&
                std::abs(
                    transition_switch_total_molar_rate -
                    transition_target_rate) <=
                    transition_rate_residual_tolerance &&
                std::abs(
                    transition_switch_report
                        ->accepted_solve
                        ->total_molar_rate_residual_mol_per_s()) <=
                    transition_rate_residual_tolerance,
            "restarted RATE generation did not satisfy the rate residual contract");
    } else {
        require_collective(
            transition_switch_control.control ==
                    wdp::
                        FixedTotalMolarRatePhysicalTimestepControlMode3D::
                            minimum_bottom_hole_pressure &&
                transition_switch_control
                        .bottom_hole_pressure_pa ==
                    transition_switch_minimum_bhp_pa &&
                transition_switch_report
                        ->accepted_fixed_bhp_solve
                        .has_value(),
            "restarted arbitration produced an invalid final minimum-BHP control state");
    }

    const auto transition_switch_local_final =
        owned_conserved_totals(
            *transition_switch_system,
            transition_switch_system
                ->initial_state());
    std::array<double, 4>
        transition_switch_global_final{};
    require_collective(
        MPI_Allreduce(
            transition_switch_local_final.data(),
            transition_switch_global_final.data(),
            4,
            MPI_DOUBLE,
            MPI_SUM,
            PETSC_COMM_WORLD) ==
            MPI_SUCCESS,
        "failed to reduce control/phase transactional final totals");
    for (std::size_t quantity = 0U;
         quantity < 4U;
         ++quantity) {
        near_collective(
            transition_switch_global_final[
                quantity],
            transition_switch_fixture
                    .rebuilt_previous_totals[
                        quantity] -
                transition_switch_well_rate[
                    quantity],
            4.0e-7,
            quantity < 3U
                ? 4.0e-8
                : 4.0e-7);
    }

    bool transition_switch_history_matches =
        false;
    require_collective(
        transition_switch_system
                ->accepted_history_matches_state(
                    transition_switch_system
                        ->initial_state(),
                    &transition_switch_history_matches) ==
                PETSC_SUCCESS &&
            transition_switch_history_matches,
        "control/phase transactional restart did not commit exactly one final history");

    // Repeat from an already accepted minimum-BHP mode. The 2P fixed-BHP
    // candidate is discarded by the phase transition, the completion is
    // rebound to 3P, and the same timestep restarts from the accepted
    // minimum-BHP mode rather than carrying any trial control decision.
    constexpr double
        transition_minimum_bhp_pa =
            5.0;
    auto transition_bhp_source =
        wdp::
            FixedTotalMolarRateWellSourceEvaluatorContext3D::
                create(
                    make_transition_well(
                        transition_minimum_bhp_pa),
                    transition_minimum_bhp_pa);
    auto transition_bhp_system =
        make_fixed_total_molar_rate_reservoir_system(
            rank,
            schedule,
            partition,
            bridge,
            pattern,
            audit,
            &transition_bhp_source,
            &provider);
    auto transition_bhp_control =
        wdp::
            AcceptedFixedTotalMolarRateWellControlState3D::
                minimum_bottom_hole_pressure(
                    transition_minimum_bhp_pa);
    fdp::AcceptedPhysicalTimeClock3D
        transition_bhp_clock{
            0.0,
            1.0};
    ControlledWellPhaseTransitionFixture
        transition_bhp_fixture{
            rank,
            &schedule,
            &partition,
            &bridge,
            &pattern,
            audit,
            &provider,
            &transition_bhp_source,
            &transition_bhp_control,
            &transition_bhp_clock,
            wdp::
                FixedTotalMolarRatePhysicalTimestepControlMode3D::
                    minimum_bottom_hole_pressure,
            transition_minimum_bhp_pa,
            0U,
            0U,
            true,
            false,
            {},
            true,
            {}};

    auto transition_bhp_options =
        transition_rate_options;
    transition_bhp_options
        .minimum_bottom_hole_pressure_pa =
        transition_minimum_bhp_pa;

    std::optional<
        wdp::
            FixedTotalMolarRatePhysicalTimestepDriverReport3D>
        transition_bhp_report;
    error =
        wdp::
            advance_fixed_total_molar_rate_controlled_physical_timestep_with_phase_transitions_3d(
                PETSC_COMM_WORLD,
                &transition_bhp_system,
                &transition_bhp_source,
                transition_target_rate,
                {
                    &controlled_well_transition_scan,
                    &transition_bhp_fixture,
                    &controlled_well_transition_rebuild,
                    &transition_bhp_fixture},
                transition_bhp_options,
                &transition_bhp_clock,
                &transition_bhp_control,
                &transition_bhp_report);

    const auto* rebound_bhp_cell30 =
        transition_bhp_source
            .well()
            .find_connection(
                mesh::GlobalEntityId{
                    UINT64_C(30)});
    require_collective(
        error == PETSC_SUCCESS &&
            transition_bhp_report
                .has_value() &&
            transition_bhp_report
                    ->accepted() &&
            transition_bhp_report
                    ->phase_transition_restarts ==
                1U &&
            transition_bhp_report
                    ->accepted_record
                    ->phase_transition_restarts ==
                1U &&
            transition_bhp_fixture.scans ==
                2U &&
            transition_bhp_fixture.rebuild_calls ==
                1U &&
            transition_bhp_fixture
                .accepted_anchor_unchanged_at_rebuild &&
            transition_bhp_fixture
                .accepted_anchor_unchanged_at_scan &&
            transition_bhp_control
                    .control ==
                wdp::
                    FixedTotalMolarRatePhysicalTimestepControlMode3D::
                        minimum_bottom_hole_pressure &&
            transition_bhp_control
                    .bottom_hole_pressure_pa ==
                transition_minimum_bhp_pa &&
            !transition_bhp_report
                 ->accepted_solve
                 .has_value() &&
            transition_bhp_report
                    ->accepted_fixed_bhp_solve
                    .has_value() &&
            transition_bhp_system
                    ->numbering()
                    .cell(
                        mesh::LocalIndex{2U})
                    .phase_count ==
                3U &&
            transition_bhp_clock
                    .accepted_time_seconds() ==
                1.0 &&
            transition_bhp_clock
                    .accepted_step_count() ==
                1U &&
            rebound_bhp_cell30 !=
                nullptr &&
            rebound_bhp_cell30
                    ->active_phase_identities()
                    .has_value() &&
            rebound_bhp_cell30
                    ->active_phase_identities()
                    ->phase_count() ==
                3U,
        "minimum-BHP post-SNES transition did not restart from accepted BHP control after rebind");

    std::uint64_t
        transition_bhp_authoritative_count = 0U;
    const auto transition_bhp_well_rate =
        global_fixed_bhp_well_production_rate(
            *transition_bhp_system,
            transition_bhp_system
                ->initial_state(),
            transition_bhp_source.well(),
            transition_minimum_bhp_pa,
            &transition_bhp_authoritative_count);
    require_collective(
        transition_bhp_authoritative_count ==
            2U,
        "minimum-BHP transition lost owner-only two-connection aggregation");

    const auto transition_bhp_local_final =
        owned_conserved_totals(
            *transition_bhp_system,
            transition_bhp_system
                ->initial_state());
    std::array<double, 4>
        transition_bhp_global_final{};
    require_collective(
        MPI_Allreduce(
            transition_bhp_local_final.data(),
            transition_bhp_global_final.data(),
            4,
            MPI_DOUBLE,
            MPI_SUM,
            PETSC_COMM_WORLD) ==
            MPI_SUCCESS,
        "failed to reduce transitioned minimum-BHP final totals");
    for (std::size_t quantity = 0U;
         quantity < 4U;
         ++quantity) {
        near_collective(
            transition_bhp_global_final[
                quantity],
            transition_bhp_fixture
                    .rebuilt_previous_totals[
                        quantity] -
                transition_bhp_well_rate[
                    quantity],
            4.0e-7,
            quantity < 3U
                ? 4.0e-8
                : 4.0e-7);
    }

    bool transition_bhp_history_matches =
        false;
    require_collective(
        transition_bhp_system
                ->accepted_history_matches_state(
                    transition_bhp_system
                        ->initial_state(),
                    &transition_bhp_history_matches) ==
                PETSC_SUCCESS &&
            transition_bhp_history_matches,
        "transitioned minimum-BHP accepted history did not commit exactly once");

    require_collective(
        MatDestroy(
            &controlled_jacobian) ==
                PETSC_SUCCESS &&
            VecDestroy(
                &reservoir_solution) ==
                PETSC_SUCCESS &&
            VecDestroy(
                &solution) ==
                PETSC_SUCCESS,
        "rate-control regression cleanup failed");

    audit->well_timestep_compressibility =
        false;
}

