void run_rebound_fixed_bhp_physical_timestep(
    int rank,
    std::unique_ptr<
        fdp::
            PhaseTransitionRebuiltNaturalVariableSystem3D>*
                system,
    Vec* stable_state,
    ControllerFixture* fixture,
    std::size_t expected_phase_count,
    std::size_t expected_scalar_count) {
    if (system == nullptr ||
        *system == nullptr ||
        stable_state == nullptr ||
        *stable_state == nullptr ||
        fixture == nullptr ||
        fixture->well_context == nullptr ||
        expected_phase_count == 0U ||
        expected_phase_count > 3U) {
        throw std::invalid_argument(
            "invalid rebound fixed-BHP physical timestep fixture");
    }

    const auto& target_record =
        (*system)
            ->numbering()
            .cell(
                mesh::LocalIndex{2U});
    require_collective(
        target_record.cell_global ==
                mesh::GlobalEntityId{
                    UINT64_C(30)} &&
            target_record.phase_count ==
                expected_phase_count &&
            target_record.scalar_count ==
                expected_scalar_count,
        "rebound fixed-BHP target chart does not match expected cardinality");

    const bool local_target_owner =
        target_record.owner_rank ==
        (*system)
            ->numbering()
            .local_rank();

    require_collective(
        (*system)
                ->rebase_accepted_timestep(
                    *stable_state,
                    1.0) ==
            PETSC_SUCCESS,
        "failed to rebase stable rebound fixed-BHP state");
    require_collective(
        VecDestroy(
            stable_state) ==
            PETSC_SUCCESS,
        "failed to destroy stable rebound fixed-BHP controller state");

    bool rebased_history_matches = false;
    require_collective(
        (*system)
                ->accepted_history_matches_state(
                    (*system)
                        ->initial_state(),
                    &rebased_history_matches) ==
            PETSC_SUCCESS &&
            rebased_history_matches,
        "rebound fixed-BHP accepted history did not match stable state");

    const auto local_previous_total =
        owned_conserved_totals(
            **system,
            (*system)
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
        "failed to reduce rebound fixed-BHP previous conserved totals");

    Vec previous_state = nullptr;
    require_collective(
        VecDuplicate(
            (*system)
                ->initial_state(),
            &previous_state) ==
                PETSC_SUCCESS &&
            VecCopy(
                (*system)
                    ->initial_state(),
                previous_state) ==
                PETSC_SUCCESS,
        "failed to preserve rebound fixed-BHP previous state");

    fixture->well_source_audit.evaluator_calls =
        0U;
    fixture->well_source_audit.target_calls =
        0U;

    fdp::PhysicalTimestepDriverOptions3D
        timestep_options;
    timestep_options.adaptive
        .minimum_timestep_seconds =
        0.125;
    timestep_options.adaptive
        .maximum_timestep_seconds =
        1.0;
    timestep_options.adaptive
        .cutback_factor =
        0.5;
    timestep_options.adaptive
        .growth_factor =
        2.0;
    timestep_options.adaptive
        .maximum_retries =
        4U;
    timestep_options.adaptive
        .growth_nonlinear_iteration_limit =
        20;
    timestep_options.adaptive
        .growth_line_search_direction_change_limit =
        4;
    timestep_options.adaptive
        .growth_transition_restart_limit =
        0U;
    timestep_options.phase_transition
        .max_transition_restarts =
        0U;

    fdp::AcceptedPhysicalTimeClock3D
        timestep_clock{
            0.0,
            1.0};
    FrozenWellTimestepControlAudit
        frozen_scan{
            mesh::LocalIndex{2U},
            mesh::GlobalEntityId{
                UINT64_C(30)},
            expected_phase_count,
            expected_scalar_count,
            0U,
            0U};
    std::optional<
        fdp::PhysicalTimestepDriverReport3D>
        timestep_report;
    const PetscErrorCode error =
        fdp::advance_one_physical_timestep_3d(
            PETSC_COMM_WORLD,
            system,
            0U,
            {
                &frozen_well_timestep_scan,
                &frozen_scan,
                &unexpected_frozen_well_timestep_rebuild,
                &frozen_scan},
            timestep_options,
            &timestep_clock,
            &timestep_report);

    require_collective(
        error == PETSC_SUCCESS &&
            timestep_report.has_value() &&
            timestep_report->accepted() &&
            timestep_report
                    ->adaptive
                    .accepted_timestep_seconds
                    .has_value() &&
            timestep_report
                    ->accepted_record
                    .has_value() &&
            timestep_report
                    ->accepted_record
                    ->phase_transition_restarts ==
                0U &&
            timestep_clock
                    .accepted_step_count() ==
                1U &&
            frozen_scan.scans > 0U &&
            frozen_scan.rebuild_calls == 0U &&
            (*system)
                    ->numbering()
                    .cell(
                        mesh::LocalIndex{2U})
                    .phase_count ==
                expected_phase_count &&
            (local_target_owner
                 ? fixture
                           ->well_source_audit
                           .target_calls >
                       0U
                 : fixture
                           ->well_source_audit
                           .target_calls ==
                       0U),
        "rebound fixed-BHP completion did not complete an owner-only accepted physical timestep");

    const double accepted_dt =
        *timestep_report
             ->adaptive
             .accepted_timestep_seconds;

    Vec state_delta = nullptr;
    PetscReal state_delta_norm = 0.0;
    require_collective(
        VecDuplicate(
            (*system)
                ->initial_state(),
            &state_delta) ==
                PETSC_SUCCESS &&
            VecCopy(
                (*system)
                    ->initial_state(),
                state_delta) ==
                PETSC_SUCCESS &&
            VecAXPY(
                state_delta,
                PetscScalar{-1.0},
                previous_state) ==
                PETSC_SUCCESS &&
            VecNorm(
                state_delta,
                NORM_2,
                &state_delta_norm) ==
                PETSC_SUCCESS &&
            static_cast<double>(
                state_delta_norm) >
                1.0e-10,
        "rebound fixed-BHP physical timestep did not change reservoir state");

    const auto local_final_total =
        owned_conserved_totals(
            **system,
            (*system)
                ->initial_state());
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
        "failed to reduce rebound fixed-BHP final conserved totals");

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
        (*system)
                ->evaluate_local_cells_for_phase_transition(
                    (*system)
                        ->initial_state(),
                    &final_current,
                    &final_porosity,
                    &final_status) ==
            PETSC_SUCCESS &&
            final_status ==
                fdp::
                    NaturalVariableSnesEvaluationStatus3D::
                        success &&
            final_current.size() >
                2U &&
            final_current[2U]
                .has_value(),
        "failed to evaluate accepted rebound fixed-BHP well state");

    std::array<double, 4>
        local_well_production_rate{};
    bool local_final_source_shape_ok =
        true;
    if (local_target_owner) {
        const auto final_well =
            wdp::
                build_fixed_bhp_peaceman_well_source_3d(
                    *fixture->well_context,
                    *final_current[2U]);
        local_final_source_shape_ok =
            final_well
                    .cell_source
                    .input_count ==
                expected_scalar_count;
        for (std::size_t component = 0U;
             component < 3U;
             ++component) {
            local_well_production_rate[
                component] =
                -final_well
                     .cell_source
                     .component_molar_rate_mol_per_s[
                         component];
        }
        local_well_production_rate[3] =
            -final_well
                 .cell_source
                 .energy_rate_w;
    }
    require_collective(
        local_final_source_shape_ok,
        "accepted rebound fixed-BHP source lost expected Jacobian cardinality");

    std::array<double, 4>
        global_well_production_rate{};
    require_collective(
        MPI_Allreduce(
            local_well_production_rate.data(),
            global_well_production_rate.data(),
            4,
            MPI_DOUBLE,
            MPI_SUM,
            PETSC_COMM_WORLD) ==
            MPI_SUCCESS,
        "failed to reduce accepted rebound fixed-BHP well rates");

    for (std::size_t component = 0U;
         component < 3U;
         ++component) {
        require_collective(
            std::isfinite(
                global_well_production_rate[
                    component]) &&
                std::abs(
                    global_well_production_rate[
                        component]) >
                    0.0,
            "accepted rebound fixed-BHP component rate vanished or became non-finite");
        near_collective(
            global_final_total[
                component],
            global_previous_total[
                component] -
                accepted_dt *
                    global_well_production_rate[
                        component],
            2.0e-7,
            2.0e-8);
    }
    require_collective(
        std::isfinite(
            global_well_production_rate[3]) &&
            std::abs(
                global_well_production_rate[3]) >
                0.0,
        "accepted rebound fixed-BHP energy rate vanished or became non-finite");
    near_collective(
        global_final_total[3],
        global_previous_total[3] -
            accepted_dt *
                global_well_production_rate[3],
        2.0e-7,
        2.0e-7);

    bool accepted_history_matches = false;
    require_collective(
        (*system)
                ->accepted_history_matches_state(
                    (*system)
                        ->initial_state(),
                    &accepted_history_matches) ==
            PETSC_SUCCESS &&
            accepted_history_matches,
        "rebound fixed-BHP accepted timestep did not rebase history");

    require_collective(
        VecDestroy(
            &state_delta) ==
                PETSC_SUCCESS &&
            VecDestroy(
                &previous_state) ==
                PETSC_SUCCESS,
        "rebound fixed-BHP physical timestep cleanup failed");

    (void)rank;
}

void run_rebound_multi_connection_fixed_bhp_physical_timestep(
    int rank,
    std::unique_ptr<
        fdp::
            PhaseTransitionRebuiltNaturalVariableSystem3D>*
                system,
    Vec* stable_state,
    ControllerFixture* fixture,
    std::size_t expected_cell30_phase_count,
    std::size_t expected_cell30_scalar_count,
    PetscInt expected_global_scalar_count) {
    if (system == nullptr ||
        *system == nullptr ||
        stable_state == nullptr ||
        *stable_state == nullptr ||
        fixture == nullptr ||
        fixture->multi_well_context == nullptr) {
        throw std::invalid_argument(
            "invalid rebound multi-connection fixed-BHP physical timestep fixture");
    }

    const auto& record30 =
        (*system)->numbering().cell(
            mesh::LocalIndex{2U});
    const auto& record60 =
        (*system)->numbering().cell(
            mesh::LocalIndex{5U});
    require_collective(
        (*system)->numbering()
                .petsc_global_scalar_count() ==
            expected_global_scalar_count &&
            record30.cell_global ==
                mesh::GlobalEntityId{
                    UINT64_C(30)} &&
            record30.phase_count ==
                expected_cell30_phase_count &&
            record30.scalar_count ==
                expected_cell30_scalar_count &&
            record30.owner_rank.value() == 0U &&
            record60.cell_global ==
                mesh::GlobalEntityId{
                    UINT64_C(60)} &&
            record60.phase_count == 3U &&
            record60.scalar_count == 10U &&
            record60.owner_rank.value() == 1U,
        "rebound multi-connection fixed-BHP target charts or ownership changed unexpectedly");

    require_collective(
        (*system)
                ->rebase_accepted_timestep(
                    *stable_state,
                    1.0) ==
            PETSC_SUCCESS &&
            VecDestroy(
                stable_state) ==
            PETSC_SUCCESS,
        "failed to rebase stable multi-connection fixed-BHP state");

    bool rebased_history_matches = false;
    require_collective(
        (*system)
                ->accepted_history_matches_state(
                    (*system)
                        ->initial_state(),
                    &rebased_history_matches) ==
            PETSC_SUCCESS &&
            rebased_history_matches,
        "rebound multi-connection fixed-BHP history did not match stable state");

    const auto local_previous_total =
        owned_conserved_totals(
            **system,
            (*system)
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
        "failed to reduce rebound multi-connection previous conserved totals");

    Vec previous_state = nullptr;
    require_collective(
        VecDuplicate(
            (*system)
                ->initial_state(),
            &previous_state) ==
                PETSC_SUCCESS &&
            VecCopy(
                (*system)
                    ->initial_state(),
                previous_state) ==
                PETSC_SUCCESS,
        "failed to preserve rebound multi-connection previous state");

    fixture
        ->multi_well_source_audit
        .evaluator_calls = 0U;
    fixture
        ->multi_well_source_audit
        .cell30_calls = 0U;
    fixture
        ->multi_well_source_audit
        .cell60_calls = 0U;

    fdp::PhysicalTimestepDriverOptions3D
        timestep_options;
    timestep_options.adaptive
        .minimum_timestep_seconds =
        0.125;
    timestep_options.adaptive
        .maximum_timestep_seconds =
        1.0;
    timestep_options.adaptive
        .cutback_factor =
        0.5;
    timestep_options.adaptive
        .growth_factor =
        2.0;
    timestep_options.adaptive
        .maximum_retries =
        4U;
    timestep_options.adaptive
        .growth_nonlinear_iteration_limit =
        20;
    timestep_options.adaptive
        .growth_line_search_direction_change_limit =
        4;
    timestep_options.adaptive
        .growth_transition_restart_limit =
        0U;
    timestep_options.phase_transition
        .max_transition_restarts =
        0U;

    FrozenMultiWellTimestepControlAudit
        control_audit;
    control_audit
        .expected_cell30_phase_count =
        expected_cell30_phase_count;
    control_audit
        .expected_cell30_scalar_count =
        expected_cell30_scalar_count;
    control_audit
        .expected_global_scalar_count =
        expected_global_scalar_count;

    fdp::AcceptedPhysicalTimeClock3D
        timestep_clock{
            0.0,
            1.0};
    std::optional<
        fdp::PhysicalTimestepDriverReport3D>
        timestep_report;
    const PetscErrorCode error =
        fdp::advance_one_physical_timestep_3d(
            PETSC_COMM_WORLD,
            system,
            0U,
            {
                &frozen_multi_well_timestep_scan,
                &control_audit,
                &unexpected_frozen_multi_well_timestep_rebuild,
                &control_audit},
            timestep_options,
            &timestep_clock,
            &timestep_report);

    require_collective(
        error == PETSC_SUCCESS &&
            timestep_report.has_value() &&
            timestep_report->accepted() &&
            timestep_report
                    ->adaptive
                    .accepted_timestep_seconds
                    .has_value() &&
            timestep_clock
                    .accepted_step_count() ==
                1U &&
            control_audit.scans > 0U &&
            control_audit.rebuild_calls == 0U &&
            fixture
                    ->multi_well_source_audit
                    .evaluator_calls >
                0U &&
            (rank == 0
                 ? fixture
                           ->multi_well_source_audit
                           .cell30_calls >
                       0U &&
                       fixture
                               ->multi_well_source_audit
                               .cell60_calls ==
                           0U
                 : fixture
                           ->multi_well_source_audit
                           .cell30_calls ==
                       0U &&
                       fixture
                               ->multi_well_source_audit
                               .cell60_calls >
                           0U),
        "rebound multi-connection sources were not assembled exactly on authoritative owners");

    const double accepted_dt =
        *timestep_report
             ->adaptive
             .accepted_timestep_seconds;

    Vec state_delta = nullptr;
    PetscReal state_delta_norm = 0.0;
    require_collective(
        VecDuplicate(
            (*system)
                ->initial_state(),
            &state_delta) ==
                PETSC_SUCCESS &&
            VecCopy(
                (*system)
                    ->initial_state(),
                state_delta) ==
                PETSC_SUCCESS &&
            VecAXPY(
                state_delta,
                PetscScalar{-1.0},
                previous_state) ==
                PETSC_SUCCESS &&
            VecNorm(
                state_delta,
                NORM_2,
                &state_delta_norm) ==
                PETSC_SUCCESS &&
            static_cast<double>(
                state_delta_norm) >
                1.0e-10,
        "rebound multi-connection accepted timestep did not change reservoir state");

    const auto local_final_total =
        owned_conserved_totals(
            **system,
            (*system)
                ->initial_state());
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
        "failed to reduce rebound multi-connection final conserved totals");

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
        (*system)
                ->evaluate_local_cells_for_phase_transition(
                    (*system)
                        ->initial_state(),
                    &final_current,
                    &final_porosity,
                    &final_status) ==
            PETSC_SUCCESS &&
            final_status ==
                fdp::
                    NaturalVariableSnesEvaluationStatus3D::
                        success &&
            final_current.size() > 5U &&
            final_current[2U].has_value() &&
            final_current[5U].has_value(),
        "failed to evaluate accepted rebound multi-connection states");

    std::vector<
        wd::
            FixedBhpConnectionCellSourceLinearization3D>
        local_authoritative_sources;
    std::array<double, 4>
        local_connection30_rate{};
    std::array<double, 4>
        local_connection60_rate{};
    bool local_connection_shapes_ok =
        true;

    if (record30.owner_rank ==
        (*system)->numbering().local_rank()) {
        const auto* connection =
            fixture->multi_well_context
                ->find_connection(
                    record30.cell_global);
        if (connection == nullptr) {
            local_connection_shapes_ok =
                false;
        } else {
            const auto source30 =
                wdp::
                    build_fixed_bhp_peaceman_well_source_3d(
                        *connection,
                        *final_current[2U]);
            local_connection_shapes_ok =
                local_connection_shapes_ok &&
                source30.cell_source.input_count ==
                    expected_cell30_scalar_count;
            for (std::size_t component = 0U;
                 component < 3U;
                 ++component) {
                local_connection30_rate[
                    component] =
                    -source30
                         .cell_source
                         .component_molar_rate_mol_per_s[
                             component];
            }
            local_connection30_rate[3] =
                -source30
                     .cell_source
                     .energy_rate_w;
            local_authoritative_sources
                .push_back(source30);
        }
    }

    if (record60.owner_rank ==
        (*system)->numbering().local_rank()) {
        const auto* connection =
            fixture->multi_well_context
                ->find_connection(
                    record60.cell_global);
        if (connection == nullptr) {
            local_connection_shapes_ok =
                false;
        } else {
            const auto source60 =
                wdp::
                    build_fixed_bhp_peaceman_well_source_3d(
                        *connection,
                        *final_current[5U]);
            local_connection_shapes_ok =
                local_connection_shapes_ok &&
                source60.cell_source.input_count ==
                    10U;
            for (std::size_t component = 0U;
                 component < 3U;
                 ++component) {
                local_connection60_rate[
                    component] =
                    -source60
                         .cell_source
                         .component_molar_rate_mol_per_s[
                             component];
            }
            local_connection60_rate[3] =
                -source60
                     .cell_source
                     .energy_rate_w;
            local_authoritative_sources
                .push_back(source60);
        }
    }

    require_collective(
        local_connection_shapes_ok &&
            local_authoritative_sources.size() ==
                1U,
        "rebound multi-connection owner did not expose exactly one authoritative connection");

    const auto local_well_total =
        wd::
            aggregate_fixed_bhp_connection_production_rates_3d(
                local_authoritative_sources);
    std::array<double, 4>
        local_well_rate{};
    for (std::size_t component = 0U;
         component < 3U;
         ++component) {
        local_well_rate[component] =
            local_well_total
                .component_molar_rate_mol_per_s[
                    component];
    }
    local_well_rate[3] =
        local_well_total.energy_rate_w;

    std::array<double, 4>
        global_well_rate{};
    std::array<double, 4>
        global_connection30_rate{};
    std::array<double, 4>
        global_connection60_rate{};
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
                local_connection30_rate.data(),
                global_connection30_rate.data(),
                4,
                MPI_DOUBLE,
                MPI_SUM,
                PETSC_COMM_WORLD) ==
                MPI_SUCCESS &&
            MPI_Allreduce(
                local_connection60_rate.data(),
                global_connection60_rate.data(),
                4,
                MPI_DOUBLE,
                MPI_SUM,
                PETSC_COMM_WORLD) ==
                MPI_SUCCESS,
        "failed to reduce rebound multi-connection authoritative rates");

    std::uint64_t local_authoritative_count =
        static_cast<std::uint64_t>(
            local_well_total
                .authoritative_connection_count);
    std::uint64_t global_authoritative_count =
        0U;
    require_collective(
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
        "rebound multi-connection well did not retain exactly two authoritative connections");

    for (std::size_t quantity = 0U;
         quantity < 4U;
         ++quantity) {
        require_collective(
            std::isfinite(
                global_connection30_rate[
                    quantity]) &&
                std::isfinite(
                    global_connection60_rate[
                        quantity]) &&
                std::abs(
                    global_connection30_rate[
                        quantity]) >
                    0.0 &&
                std::abs(
                    global_connection60_rate[
                        quantity]) >
                    0.0,
            "rebound multi-connection accepted rate vanished or became non-finite");
        near_collective(
            global_well_rate[quantity],
            global_connection30_rate[
                quantity] +
                global_connection60_rate[
                    quantity],
            2.0e-12,
            2.0e-12);
        near_collective(
            global_final_total[quantity],
            global_previous_total[quantity] -
                accepted_dt *
                    global_well_rate[quantity],
            2.0e-7,
            quantity < 3U
                ? 2.0e-8
                : 2.0e-7);
    }

    bool history_matches = false;
    require_collective(
        (*system)
                ->accepted_history_matches_state(
                    (*system)
                        ->initial_state(),
                    &history_matches) ==
            PETSC_SUCCESS &&
            history_matches,
        "rebound multi-connection accepted timestep did not rebase history");

    require_collective(
        VecDestroy(
            &state_delta) ==
                PETSC_SUCCESS &&
            VecDestroy(
                &previous_state) ==
                PETSC_SUCCESS,
        "rebound multi-connection physical timestep cleanup failed");
}

void run_frozen_fixed_bhp_timestep_case(
    int rank,
    std::uint64_t target_stable,
    std::size_t target_local,
    std::size_t expected_phase_count,
    double bottom_hole_pressure_pa,
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
    if (audit == nullptr ||
        expected_phase_count == 0U ||
        expected_phase_count > 3U ||
        target_local >= 6U) {
        throw std::invalid_argument(
            "invalid fixed-BHP cardinality timestep case");
    }

    audit->well_timestep_compressibility =
        true;

    std::vector<double>
        injection_enthalpy;
    injection_enthalpy.reserve(
        expected_phase_count);
    for (std::size_t phase = 0U;
         phase < expected_phase_count;
         ++phase) {
        injection_enthalpy.push_back(
            1000.0 *
            static_cast<double>(
                phase + 1U));
    }

    auto timestep_well_context =
        wdp::
            FixedBhpPeacemanWellSourceEvaluatorContext3D::
                create(
                    mesh::GlobalEntityId{
                        target_stable},
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
                    bottom_hole_pressure_pa,
                    wd::FixedBhpInjectionEnthalpy3D{
                        "fixture/fixed-bhp-variable-cardinality-injection-enthalpy/v1",
                        std::move(
                            injection_enthalpy)},
                    "fixture/fixed-bhp-variable-cardinality-source/v1");

    const auto initial_target_current =
        evaluate_target(
            target_stable,
            audit);
    const auto initial_well =
        wdp::
            build_fixed_bhp_peaceman_well_source_3d(
                timestep_well_context,
                initial_target_current);
    require_collective(
        initial_well
                .phase_volumetric_rate_m3_per_s
                .size() ==
            expected_phase_count,
        "fixed-BHP initial production source changed active phase cardinality");
    for (double phase_rate :
         initial_well
             .phase_volumetric_rate_m3_per_s) {
        require_collective(
            phase_rate > 0.0,
            "fixed-BHP regression was not configured production-positive at the frozen initial state");
    }

    const auto timestep_pr_parameters =
        thermodynamic_adapter_pr_parameters();
    const auto timestep_pr_model =
        th::Pr76Phase<double>::
            from_parameters(
                timestep_pr_parameters);
    auto timestep_provider =
        make_outer_rebuild_pr_provider(
            timestep_pr_model);

    FixedBhpWellSourceAudit
        timestep_source_audit{
            &timestep_well_context};
    FixedBhpWellSourceAudit
        verification_source_audit{
            &timestep_well_context};

    auto timestep_system =
        make_frozen_fixed_bhp_timestep_system(
            rank,
            schedule,
            partition,
            bridge,
            pattern,
            audit,
            &timestep_source_audit,
            &timestep_provider);
    auto verification_system =
        make_frozen_fixed_bhp_timestep_system(
            rank,
            schedule,
            partition,
            bridge,
            pattern,
            audit,
            &verification_source_audit,
            &timestep_provider);

    const auto target_cell =
        mesh::LocalIndex{
            static_cast<
                mesh::LocalIndex::value_type>(
                    target_local)};
    const auto& target_record =
        timestep_system
            ->numbering()
            .cell(target_cell);
    const std::size_t expected_scalar_count =
        expected_phase_count * 3U +
        1U;
    require_collective(
        timestep_system
                ->numbering()
                .petsc_global_scalar_count() ==
            42 &&
            target_record.cell_global ==
                mesh::GlobalEntityId{
                    target_stable} &&
            target_record.phase_count ==
                expected_phase_count &&
            target_record.scalar_count ==
                expected_scalar_count,
        "fixed-BHP timestep fixture changed variable-cardinality reservoir numbering");

    const bool local_target_owner =
        target_record.owner_rank ==
        timestep_system
            ->numbering()
            .local_rank();

    const auto local_previous_total =
        owned_conserved_totals(
            *timestep_system,
            timestep_system
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
        "failed to reduce variable-cardinality fixed-BHP previous conserved totals");

    Vec previous_state = nullptr;
    require_collective(
        VecDuplicate(
            timestep_system
                ->initial_state(),
            &previous_state) ==
                PETSC_SUCCESS &&
            VecCopy(
                timestep_system
                    ->initial_state(),
                previous_state) ==
                PETSC_SUCCESS,
        "failed to preserve variable-cardinality fixed-BHP initial state");

    fdp::PhysicalTimestepDriverOptions3D
        timestep_options;
    timestep_options.adaptive
        .minimum_timestep_seconds =
        0.125;
    timestep_options.adaptive
        .maximum_timestep_seconds =
        1.0;
    timestep_options.adaptive
        .cutback_factor =
        0.5;
    timestep_options.adaptive
        .growth_factor =
        2.0;
    timestep_options.adaptive
        .maximum_retries =
        4U;
    timestep_options.adaptive
        .growth_nonlinear_iteration_limit =
        20;
    timestep_options.adaptive
        .growth_line_search_direction_change_limit =
        4;
    timestep_options.adaptive
        .growth_transition_restart_limit =
        0U;
    timestep_options.phase_transition
        .max_transition_restarts =
        0U;

    fdp::AcceptedPhysicalTimeClock3D
        timestep_clock{
            0.0,
            1.0};
    FrozenWellTimestepControlAudit
        timestep_control_audit{
            target_cell,
            mesh::GlobalEntityId{
                target_stable},
            expected_phase_count,
            expected_scalar_count,
            0U,
            0U};
    std::optional<
        fdp::PhysicalTimestepDriverReport3D>
        timestep_report;
    PetscErrorCode error =
        fdp::advance_one_physical_timestep_3d(
            PETSC_COMM_WORLD,
            &timestep_system,
            0U,
            {
                &frozen_well_timestep_scan,
                &timestep_control_audit,
                &unexpected_frozen_well_timestep_rebuild,
                &timestep_control_audit},
            timestep_options,
            &timestep_clock,
            &timestep_report);

    require_collective(
        error == PETSC_SUCCESS &&
            timestep_report.has_value() &&
            timestep_report->accepted() &&
            timestep_report
                ->adaptive
                .accepted_timestep_seconds
                .has_value() &&
            timestep_report
                ->accepted_record
                .has_value() &&
            timestep_report
                ->accepted_record
                ->phase_transition_restarts ==
            0U &&
            timestep_clock
                .accepted_step_count() ==
            1U &&
            timestep_control_audit.scans >
            0U &&
            timestep_control_audit
                .rebuild_calls ==
            0U &&
            timestep_system
                ->numbering()
                .cell(target_cell)
                .phase_count ==
            expected_phase_count &&
            (local_target_owner
                 ? timestep_source_audit
                       .target_calls >
                   0U
                 : timestep_source_audit
                       .target_calls ==
                   0U),
        "PhysicalTimestepDriver did not accept the frozen variable-cardinality fixed-BHP well step owner-only");

    const double accepted_dt =
        *timestep_report
             ->adaptive
             .accepted_timestep_seconds;

    bool history_matches = false;
    require_collective(
        timestep_system
                ->accepted_history_matches_state(
                    timestep_system
                        ->initial_state(),
                    &history_matches) ==
            PETSC_SUCCESS &&
            history_matches,
        "variable-cardinality fixed-BHP accepted timestep did not rebase component/energy history");

    Vec state_delta = nullptr;
    PetscReal state_delta_norm = 0.0;
    require_collective(
        VecDuplicate(
            timestep_system
                ->initial_state(),
            &state_delta) ==
                PETSC_SUCCESS &&
            VecCopy(
                timestep_system
                    ->initial_state(),
                state_delta) ==
                PETSC_SUCCESS &&
            VecAXPY(
                state_delta,
                PetscScalar{-1.0},
                previous_state) ==
                PETSC_SUCCESS &&
            VecNorm(
                state_delta,
                NORM_2,
                &state_delta_norm) ==
                PETSC_SUCCESS &&
            static_cast<double>(
                state_delta_norm) >
                1.0e-8,
        "variable-cardinality fixed-BHP well did not change the accepted reservoir state");

    require_collective(
        verification_system
                ->set_trial_timestep_seconds(
                    accepted_dt) ==
            PETSC_SUCCESS,
        "failed to align independent variable-cardinality fixed-BHP verification timestep");

    Vec verification_residual = nullptr;
    require_collective(
        VecDuplicate(
            timestep_system
                ->initial_state(),
            &verification_residual) ==
                PETSC_SUCCESS &&
            VecSet(
                verification_residual,
                PetscScalar{0.0}) ==
                PETSC_SUCCESS,
        "failed to allocate independent variable-cardinality fixed-BHP residual");

    auto verification_evaluator =
        verification_system
            ->snes_evaluator();
    fdp::NaturalVariableSnesEvaluationStatus3D
        verification_status =
            fdp::
                NaturalVariableSnesEvaluationStatus3D::
                    success;
    error =
        verification_evaluator.function(
            timestep_system
                ->initial_state(),
            verification_residual,
            verification_evaluator
                .user_context,
            &verification_status);
    if (error == PETSC_SUCCESS) {
        error =
            VecAssemblyBegin(
                verification_residual);
    }
    if (error == PETSC_SUCCESS) {
        error =
            VecAssemblyEnd(
                verification_residual);
    }
    PetscReal verification_residual_norm =
        0.0;
    if (error == PETSC_SUCCESS) {
        error =
            VecNorm(
                verification_residual,
                NORM_2,
                &verification_residual_norm);
    }
    require_collective(
        error == PETSC_SUCCESS &&
            verification_status ==
                fdp::
                    NaturalVariableSnesEvaluationStatus3D::
                        success &&
            static_cast<double>(
                verification_residual_norm) <=
                1.0e-6,
        "accepted variable-cardinality fixed-BHP timestep failed independent residual reassembly");

    const auto local_final_total =
        owned_conserved_totals(
            *verification_system,
            timestep_system
                ->initial_state());
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
        "failed to reduce variable-cardinality fixed-BHP final conserved totals");

    std::vector<std::optional<
        fdp::
            MixedCardinalityPhysicalCurrentCellLinearization3D>>
        final_current;
    std::vector<double>
        final_porosities;
    fdp::NaturalVariableSnesEvaluationStatus3D
        final_status =
            fdp::
                NaturalVariableSnesEvaluationStatus3D::
                    success;
    const PetscErrorCode final_error =
        verification_system
            ->evaluate_local_cells_for_phase_transition(
                timestep_system
                    ->initial_state(),
                &final_current,
                &final_porosities,
                &final_status);
    require_collective(
        final_error == PETSC_SUCCESS &&
            final_status ==
                fdp::
                    NaturalVariableSnesEvaluationStatus3D::
                        success &&
            final_current.size() >
                target_local &&
            final_porosities.size() ==
                final_current.size() &&
            final_current[target_local]
                .has_value() &&
            std::visit(
                [](const auto& typed) {
                    return typed
                        .transport
                        .state_identity
                        .layout
                        .phase_count();
                },
                *final_current[
                    target_local]) ==
                expected_phase_count,
        "failed to collectively recover accepted variable-cardinality well-cell state");

    std::array<double, 4>
        local_well_production_rate{};
    if (local_target_owner) {
        const auto final_well =
            wdp::
                build_fixed_bhp_peaceman_well_source_3d(
                    timestep_well_context,
                    *final_current[
                        target_local]);
        for (std::size_t component = 0U;
             component < 3U;
             ++component) {
            local_well_production_rate[
                component] =
                -final_well
                     .cell_source
                     .component_molar_rate_mol_per_s[
                         component];
        }
        local_well_production_rate[3] =
            -final_well
                 .cell_source
                 .energy_rate_w;
    }

    std::array<double, 4>
        global_well_production_rate{};
    require_collective(
        MPI_Allreduce(
            local_well_production_rate.data(),
            global_well_production_rate.data(),
            4,
            MPI_DOUBLE,
            MPI_SUM,
            PETSC_COMM_WORLD) ==
            MPI_SUCCESS,
        "failed to reduce variable-cardinality fixed-BHP production rates");

    for (std::size_t component = 0U;
         component < 3U;
         ++component) {
        require_collective(
            std::isfinite(
                global_well_production_rate[
                    component]) &&
                std::abs(
                    global_well_production_rate[
                        component]) >
                    0.0,
            "accepted variable-cardinality fixed-BHP component rate vanished or became non-finite");
        near_collective(
            global_final_total[
                component],
            global_previous_total[
                component] -
                accepted_dt *
                    global_well_production_rate[
                        component],
            2.0e-7,
            2.0e-8);
    }
    require_collective(
        std::isfinite(
            global_well_production_rate[3]) &&
            std::abs(
                global_well_production_rate[3]) >
                0.0,
        "accepted variable-cardinality fixed-BHP energy rate vanished or became non-finite");
    near_collective(
        global_final_total[3],
        global_previous_total[3] -
            accepted_dt *
                global_well_production_rate[
                    3],
        2.0e-7,
        2.0e-7);

    require_collective(
        VecDestroy(
            &verification_residual) ==
                PETSC_SUCCESS &&
            VecDestroy(
                &state_delta) ==
                PETSC_SUCCESS &&
            VecDestroy(
                &previous_state) ==
                PETSC_SUCCESS,
        "variable-cardinality fixed-BHP timestep regression cleanup failed");

    audit->well_timestep_compressibility =
        false;
}

void run_multi_connection_fixed_bhp_timestep_case(
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
            "invalid multi-connection fixed-BHP regression context");
    }

    audit->well_timestep_compressibility =
        true;

    constexpr double shared_bhp_pa = 5.0;
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
                    shared_bhp_pa,
                    wd::FixedBhpInjectionEnthalpy3D{
                        "fixture/multi-well-cell30-injection-enthalpy/v1",
                        {1000.0, 2000.0}},
                    "fixture/multi-well-cell30-source/v1");
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
                    shared_bhp_pa,
                    wd::FixedBhpInjectionEnthalpy3D{
                        "fixture/multi-well-cell60-injection-enthalpy/v1",
                        {1100.0, 2200.0, 3300.0}},
                    "fixture/multi-well-cell60-source/v1");

    bool duplicate_rejected = false;
    try {
        (void)wdp::
            FixedBhpMultiConnectionWellSourceEvaluatorContext3D::
                create(
                    "fixture/duplicate-completion",
                    {connection30, connection30});
    } catch (const std::invalid_argument&) {
        duplicate_rejected = true;
    }

    bool bhp_mismatch_rejected = false;
    try {
        auto mismatched =
            wdp::
                FixedBhpPeacemanWellSourceEvaluatorContext3D::
                    create(
                        mesh::GlobalEntityId{
                            UINT64_C(60)},
                        connection60.connection(),
                        shared_bhp_pa + 1.0,
                        connection60.injection_enthalpy(),
                        "fixture/multi-well-bhp-mismatch/v1");
        (void)wdp::
            FixedBhpMultiConnectionWellSourceEvaluatorContext3D::
                create(
                    "fixture/mismatched-bhp",
                    {connection30, mismatched});
    } catch (const std::invalid_argument&) {
        bhp_mismatch_rejected = true;
    }

    require_collective(
        duplicate_rejected &&
            bhp_mismatch_rejected,
        "multi-connection fixed-BHP context accepted duplicate stable cell or inconsistent BHP");

    auto multi_context =
        wdp::
            FixedBhpMultiConnectionWellSourceEvaluatorContext3D::
                create(
                    "fixture/single-well-two-connections/v1",
                    {connection60, connection30});

    require_collective(
        multi_context.connection_count() ==
                2U &&
            multi_context
                    .bottom_hole_pressure_pa() ==
                shared_bhp_pa &&
            multi_context
                    .connections()[0U]
                    .target_cell_global() ==
                mesh::GlobalEntityId{
                    UINT64_C(30)} &&
            multi_context
                    .connections()[1U]
                    .target_cell_global() ==
                mesh::GlobalEntityId{
                    UINT64_C(60)} &&
            multi_context.find_connection(
                mesh::GlobalEntityId{
                    UINT64_C(40)}) ==
                nullptr,
        "multi-connection fixed-BHP context did not retain deterministic unique stable-cell lookup");

    auto initial30 =
        evaluate_target(
            UINT64_C(30),
            audit);
    auto initial60 =
        evaluate_target(
            UINT64_C(60),
            audit);
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
    const std::vector<
        wd::
            FixedBhpConnectionCellSourceLinearization3D>
        direct_connection_sources{
            initial_source30,
            initial_source60};
    const auto direct_total =
        wd::
            aggregate_fixed_bhp_connection_production_rates_3d(
                direct_connection_sources);

    bool initial_production_positive =
        initial_source30
                .cell_source
                .input_count ==
            7U &&
        initial_source60
                .cell_source
                .input_count ==
            10U &&
        direct_total
                .authoritative_connection_count ==
            2U &&
        direct_total.bottom_hole_pressure_pa ==
            shared_bhp_pa;
    for (double rate :
         initial_source30
             .phase_volumetric_rate_m3_per_s) {
        initial_production_positive =
            initial_production_positive &&
            rate > 0.0;
    }
    for (double rate :
         initial_source60
             .phase_volumetric_rate_m3_per_s) {
        initial_production_positive =
            initial_production_positive &&
            rate > 0.0;
    }
    require_collective(
        initial_production_positive,
        "multi-connection fixed-BHP fixture was not production-positive with q=7/q=10 connection Jacobians");

    for (std::size_t component = 0U;
         component < 3U;
         ++component) {
        near_collective(
            direct_total
                .component_molar_rate_mol_per_s[
                    component],
            -initial_source30
                 .cell_source
                 .component_molar_rate_mol_per_s[
                     component] -
                initial_source60
                    .cell_source
                    .component_molar_rate_mol_per_s[
                        component],
            1.0e-12,
            1.0e-12);
    }
    near_collective(
        direct_total.energy_rate_w,
        -initial_source30
             .cell_source
             .energy_rate_w -
            initial_source60
                .cell_source
                .energy_rate_w,
        1.0e-12,
        1.0e-12);

    const auto pr_parameters =
        thermodynamic_adapter_pr_parameters();
    const auto pr_model =
        th::Pr76Phase<double>::
            from_parameters(
                pr_parameters);
    auto provider =
        make_outer_rebuild_pr_provider(
            pr_model);

    FixedBhpMultiConnectionWellSourceAudit
        source_audit{
            &multi_context};
    auto system =
        make_frozen_multi_connection_fixed_bhp_timestep_system(
            rank,
            schedule,
            partition,
            bridge,
            pattern,
            audit,
            &source_audit,
            &provider);

    const auto& record30 =
        system->numbering().cell(
            mesh::LocalIndex{2U});
    const auto& record60 =
        system->numbering().cell(
            mesh::LocalIndex{5U});
    require_collective(
        system->numbering()
                .petsc_global_scalar_count() ==
            42 &&
            record30.phase_count == 2U &&
            record30.scalar_count == 7U &&
            record30.owner_rank.value() == 0U &&
            record60.phase_count == 3U &&
            record60.scalar_count == 10U &&
            record60.owner_rank.value() == 1U,
        "multi-connection fixed-BHP well changed reservoir global numbering or target ownership");

    const auto local_previous_total =
        owned_conserved_totals(
            *system,
            system->initial_state());
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
        "failed to reduce multi-connection fixed-BHP previous conserved totals");

    Vec previous_state = nullptr;
    require_collective(
        VecDuplicate(
            system->initial_state(),
            &previous_state) ==
                PETSC_SUCCESS &&
            VecCopy(
                system->initial_state(),
                previous_state) ==
                PETSC_SUCCESS,
        "failed to preserve multi-connection fixed-BHP previous state");

    fdp::PhysicalTimestepDriverOptions3D
        timestep_options;
    timestep_options.adaptive
        .minimum_timestep_seconds =
        0.125;
    timestep_options.adaptive
        .maximum_timestep_seconds =
        1.0;
    timestep_options.adaptive
        .cutback_factor =
        0.5;
    timestep_options.adaptive
        .growth_factor =
        2.0;
    timestep_options.adaptive
        .maximum_retries =
        4U;
    timestep_options.adaptive
        .growth_nonlinear_iteration_limit =
        20;
    timestep_options.adaptive
        .growth_line_search_direction_change_limit =
        4;
    timestep_options.adaptive
        .growth_transition_restart_limit =
        0U;
    timestep_options.phase_transition
        .max_transition_restarts =
        0U;

    FrozenMultiWellTimestepControlAudit
        control_audit;
    fdp::AcceptedPhysicalTimeClock3D
        timestep_clock{
            0.0,
            1.0};
    std::optional<
        fdp::PhysicalTimestepDriverReport3D>
        timestep_report;
    PetscErrorCode error =
        fdp::advance_one_physical_timestep_3d(
            PETSC_COMM_WORLD,
            &system,
            0U,
            {
                &frozen_multi_well_timestep_scan,
                &control_audit,
                &unexpected_frozen_multi_well_timestep_rebuild,
                &control_audit},
            timestep_options,
            &timestep_clock,
            &timestep_report);

    require_collective(
        error == PETSC_SUCCESS &&
            timestep_report.has_value() &&
            timestep_report->accepted() &&
            timestep_report
                    ->adaptive
                    .accepted_timestep_seconds
                    .has_value() &&
            timestep_clock
                    .accepted_step_count() ==
                1U &&
            control_audit.scans > 0U &&
            control_audit.rebuild_calls == 0U &&
            source_audit.evaluator_calls > 0U &&
            (rank == 0
                 ? source_audit.cell30_calls > 0U &&
                       source_audit.cell60_calls == 0U
                 : source_audit.cell30_calls == 0U &&
                       source_audit.cell60_calls > 0U),
        "multi-connection fixed-BHP sources were not assembled exactly on their authoritative owners");

    const double accepted_dt =
        *timestep_report
             ->adaptive
             .accepted_timestep_seconds;

    Vec state_delta = nullptr;
    PetscReal state_delta_norm = 0.0;
    require_collective(
        VecDuplicate(
            system->initial_state(),
            &state_delta) ==
                PETSC_SUCCESS &&
            VecCopy(
                system->initial_state(),
                state_delta) ==
                PETSC_SUCCESS &&
            VecAXPY(
                state_delta,
                PetscScalar{-1.0},
                previous_state) ==
                PETSC_SUCCESS &&
            VecNorm(
                state_delta,
                NORM_2,
                &state_delta_norm) ==
                PETSC_SUCCESS &&
            static_cast<double>(
                state_delta_norm) >
                1.0e-10,
        "multi-connection fixed-BHP accepted timestep did not change reservoir state");

    const auto local_final_total =
        owned_conserved_totals(
            *system,
            system->initial_state());
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
        "failed to reduce multi-connection fixed-BHP final conserved totals");

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
        system->evaluate_local_cells_for_phase_transition(
                system->initial_state(),
                &final_current,
                &final_porosity,
                &final_status) ==
            PETSC_SUCCESS &&
            final_status ==
                fdp::
                    NaturalVariableSnesEvaluationStatus3D::
                        success &&
            final_current.size() > 5U &&
            final_current[2U].has_value() &&
            final_current[5U].has_value(),
        "failed to evaluate accepted multi-connection fixed-BHP states");

    std::vector<
        wd::
            FixedBhpConnectionCellSourceLinearization3D>
        local_authoritative_sources;
    std::array<double, 4>
        local_connection30_rate{};
    std::array<double, 4>
        local_connection60_rate{};
    bool local_connection_shapes_ok =
        true;

    if (record30.owner_rank ==
        system->numbering().local_rank()) {
        const auto source30 =
            wdp::
                build_fixed_bhp_peaceman_well_source_3d(
                    *multi_context.find_connection(
                        record30.cell_global),
                    *final_current[2U]);
        local_connection_shapes_ok =
            local_connection_shapes_ok &&
            source30.cell_source.input_count ==
                7U;
        for (std::size_t component = 0U;
             component < 3U;
             ++component) {
            local_connection30_rate[
                component] =
                -source30
                     .cell_source
                     .component_molar_rate_mol_per_s[
                         component];
        }
        local_connection30_rate[3] =
            -source30.cell_source.energy_rate_w;
        local_authoritative_sources.push_back(
            source30);
    }

    if (record60.owner_rank ==
        system->numbering().local_rank()) {
        const auto source60 =
            wdp::
                build_fixed_bhp_peaceman_well_source_3d(
                    *multi_context.find_connection(
                        record60.cell_global),
                    *final_current[5U]);
        local_connection_shapes_ok =
            local_connection_shapes_ok &&
            source60.cell_source.input_count ==
                10U;
        for (std::size_t component = 0U;
             component < 3U;
             ++component) {
            local_connection60_rate[
                component] =
                -source60
                     .cell_source
                     .component_molar_rate_mol_per_s[
                         component];
        }
        local_connection60_rate[3] =
            -source60.cell_source.energy_rate_w;
        local_authoritative_sources.push_back(
            source60);
    }

    require_collective(
        local_connection_shapes_ok &&
            local_authoritative_sources.size() ==
                1U,
        "multi-connection fixed-BHP authoritative owner did not expose exactly one local connection");

    const auto local_well_total =
        wd::
            aggregate_fixed_bhp_connection_production_rates_3d(
                local_authoritative_sources);
    std::array<double, 4>
        local_well_rate{};
    for (std::size_t component = 0U;
         component < 3U;
         ++component) {
        local_well_rate[component] =
            local_well_total
                .component_molar_rate_mol_per_s[
                    component];
    }
    local_well_rate[3] =
        local_well_total.energy_rate_w;

    std::array<double, 4>
        global_well_rate{};
    std::array<double, 4>
        global_connection30_rate{};
    std::array<double, 4>
        global_connection60_rate{};
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
                local_connection30_rate.data(),
                global_connection30_rate.data(),
                4,
                MPI_DOUBLE,
                MPI_SUM,
                PETSC_COMM_WORLD) ==
                MPI_SUCCESS &&
            MPI_Allreduce(
                local_connection60_rate.data(),
                global_connection60_rate.data(),
                4,
                MPI_DOUBLE,
                MPI_SUM,
                PETSC_COMM_WORLD) ==
                MPI_SUCCESS,
        "failed to reduce multi-connection fixed-BHP authoritative rates");

    std::uint64_t local_authoritative_count =
        static_cast<std::uint64_t>(
            local_well_total
                .authoritative_connection_count);
    std::uint64_t global_authoritative_count =
        0U;
    require_collective(
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
        "multi-connection fixed-BHP well did not aggregate exactly two authoritative connections");

    for (std::size_t quantity = 0U;
         quantity < 4U;
         ++quantity) {
        require_collective(
            std::isfinite(
                global_connection30_rate[
                    quantity]) &&
                std::isfinite(
                    global_connection60_rate[
                        quantity]) &&
                std::abs(
                    global_connection30_rate[
                        quantity]) >
                    0.0 &&
                std::abs(
                    global_connection60_rate[
                        quantity]) >
                    0.0,
            "multi-connection fixed-BHP accepted connection rate vanished or became non-finite");
        near_collective(
            global_well_rate[quantity],
            global_connection30_rate[
                quantity] +
                global_connection60_rate[
                    quantity],
            2.0e-12,
            2.0e-12);
        near_collective(
            global_final_total[quantity],
            global_previous_total[quantity] -
                accepted_dt *
                    global_well_rate[quantity],
            2.0e-7,
            quantity < 3U
                ? 2.0e-8
                : 2.0e-7);
    }

    bool history_matches = false;
    require_collective(
        system->accepted_history_matches_state(
                system->initial_state(),
                &history_matches) ==
            PETSC_SUCCESS &&
            history_matches,
        "multi-connection fixed-BHP accepted timestep did not rebase history");

    require_collective(
        VecDestroy(
            &state_delta) ==
                PETSC_SUCCESS &&
            VecDestroy(
                &previous_state) ==
                PETSC_SUCCESS,
        "multi-connection fixed-BHP regression cleanup failed");

    audit->well_timestep_compressibility =
        false;
}

void run_multi_connection_local_transition_rebind_case(
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
    DispatchAudit* audit,
    flow::
        Pr76AbsentPhasePotentialExtensionProvider<
            double>* provider) {
    if (audit == nullptr ||
        provider == nullptr) {
        throw std::invalid_argument(
            "invalid local-transition multi-connection fixed-BHP fixture");
    }

    audit->well_timestep_compressibility =
        true;

    const auto aqueous =
        mixed_physical_phase_identity(
            "aqueous");
    const auto hydrocarbon0 =
        mixed_physical_phase_identity(
            "hydrocarbon-0");
    const auto hydrocarbon1 =
        mixed_physical_phase_identity(
            "hydrocarbon-1");

    wdp::FixedBhpPhaseIdentityInjectionEnthalpy3D
        cell30_phase_registry{
            "fixture/multi-transition-cell30-phase-enthalpy/v1",
            {
                {aqueous, 1000.0},
                {hydrocarbon0, 2000.0},
                {hydrocarbon1, 3000.0}}};

    constexpr double shared_bhp_pa =
        5.0;
    auto connection30 =
        wdp::
            FixedBhpPeacemanWellSourceEvaluatorContext3D::
                create_phase_identity_bound(
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
                    shared_bhp_pa,
                    cell30_phase_registry,
                    phase_identity_map_for_stable(
                        UINT64_C(30)),
                    "fixture/multi-transition-cell30-source/v1");

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
                    shared_bhp_pa,
                    wd::FixedBhpInjectionEnthalpy3D{
                        "fixture/multi-transition-cell60-enthalpy/v1",
                        {1100.0, 2200.0, 3300.0}},
                    "fixture/multi-transition-cell60-source/v1");

    const double frozen60_wi =
        connection60.connection()
            .well_index_m3;
    const auto frozen60_enthalpy =
        connection60
            .injection_enthalpy()
            .specific_enthalpy_j_per_kg;
    const std::string frozen60_source =
        connection60.source_provenance();

    auto multi_context =
        wdp::
            FixedBhpMultiConnectionWellSourceEvaluatorContext3D::
                create(
                    "fixture/multi-transition-single-well/v1",
                    {
                        connection60,
                        connection30});

    ControllerFixture fixture{
        rank,
        &schedule,
        &partition,
        &bridge,
        &pattern,
        audit,
        nullptr,
        provider,
        false,
        0U,
        false,
        false};
    fixture.multi_well_context =
        &multi_context;

    auto initial_system =
        make_controller_initial_system(
            &fixture);
    require_collective(
        initial_system != nullptr &&
            initial_system
                    ->numbering()
                    .petsc_global_scalar_count() ==
                42 &&
            initial_system
                    ->numbering()
                    .cell(
                        mesh::LocalIndex{2U})
                    .phase_count ==
                2U &&
            initial_system
                    ->numbering()
                    .cell(
                        mesh::LocalIndex{2U})
                    .scalar_count ==
                7U &&
            initial_system
                    ->numbering()
                    .cell(
                        mesh::LocalIndex{5U})
                    .phase_count ==
                3U &&
            initial_system
                    ->numbering()
                    .cell(
                        mesh::LocalIndex{5U})
                    .scalar_count ==
                10U,
        "local-transition multi-well fixture did not start on 2P(cell30)+3P(cell60)");

    std::unique_ptr<
        fdp::
            PhaseTransitionRebuiltNaturalVariableSystem3D>
        appeared_system;
    Vec appeared_state = nullptr;
    std::optional<
        fdp::
            PostSnesPhaseTransitionControllerReport3D>
        appearance_report;
    PetscErrorCode error =
        fdp::
            solve_nonlinear_timestep_with_phase_transitions_3d(
                PETSC_COMM_WORLD,
                std::move(initial_system),
                {
                    &controller_scan,
                    &fixture,
                    &controller_rebuild,
                    &fixture},
                {4U},
                &appeared_system,
                &appeared_state,
                &appearance_report);

    require_collective(
        error == PETSC_SUCCESS &&
            appeared_system != nullptr &&
            appeared_state != nullptr &&
            appearance_report.has_value() &&
            appearance_report->outcome ==
                fdp::
                    PostSnesPhaseTransitionOutcome3D::
                        stable_phase_set &&
            appearance_report
                    ->transition_restarts ==
                1U &&
            fixture.rebuild_calls ==
                1U &&
            fixture.multi_well_rebound &&
            fixture.multi_well_context !=
                nullptr &&
            appeared_system
                    ->numbering()
                    .petsc_global_scalar_count() ==
                45 &&
            appeared_system
                    ->numbering()
                    .cell(
                        mesh::LocalIndex{2U})
                    .phase_count ==
                3U &&
            appeared_system
                    ->numbering()
                    .cell(
                        mesh::LocalIndex{2U})
                    .scalar_count ==
                10U &&
            appeared_system
                    ->numbering()
                    .cell(
                        mesh::LocalIndex{5U})
                    .phase_count ==
                3U &&
            appeared_system
                    ->numbering()
                    .cell(
                        mesh::LocalIndex{5U})
                    .scalar_count ==
                10U,
        "multi-connection well did not survive local cell30 2P-to-3P restart");

    const auto* appeared30 =
        fixture.multi_well_context
            ->find_connection(
                mesh::GlobalEntityId{
                    UINT64_C(30)});
    const auto* frozen60 =
        fixture.multi_well_context
            ->find_connection(
                mesh::GlobalEntityId{
                    UINT64_C(60)});
    require_collective(
        appeared30 != nullptr &&
            frozen60 != nullptr &&
            appeared30
                    ->active_phase_identities()
                    .has_value() &&
            appeared30
                    ->active_phase_identities()
                    ->phase_count() ==
                3U &&
            appeared30
                    ->injection_enthalpy()
                    .specific_enthalpy_j_per_kg ==
                std::vector<double>{
                    1000.0,
                    2000.0,
                    3000.0} &&
            frozen60
                    ->target_cell_global() ==
                mesh::GlobalEntityId{
                    UINT64_C(60)} &&
            std::abs(
                frozen60
                        ->connection()
                        .well_index_m3 -
                    frozen60_wi) <=
                1.0e-14 *
                    std::max(
                        1.0,
                        std::abs(
                            frozen60_wi)) &&
            frozen60
                    ->bottom_hole_pressure_pa() ==
                shared_bhp_pa &&
            frozen60
                    ->injection_enthalpy()
                    .specific_enthalpy_j_per_kg ==
                frozen60_enthalpy &&
            frozen60
                    ->source_provenance() ==
                frozen60_source &&
            !frozen60
                 ->phase_identity_rebindable(),
        "cell30 local appearance incorrectly mutated frozen cell60 connection");

    run_rebound_multi_connection_fixed_bhp_physical_timestep(
        rank,
        &appeared_system,
        &appeared_state,
        &fixture,
        3U,
        10U,
        45);

    // A second independent controller session removes only cell30's third
    // phase. Cell60 remains a frozen 3P connection throughout the round trip.
    fixture.phase_disappearance_only =
        true;
    fixture.oscillate_after_restart =
        false;
    fixture.rebuild_calls =
        0U;
    fixture.multi_well_rebound =
        false;
    fixture
        .multi_well_source_audit
        .evaluator_calls = 0U;
    fixture
        .multi_well_source_audit
        .cell30_calls = 0U;
    fixture
        .multi_well_source_audit
        .cell60_calls = 0U;

    std::unique_ptr<
        fdp::
            PhaseTransitionRebuiltNaturalVariableSystem3D>
        disappeared_system;
    Vec disappeared_state = nullptr;
    std::optional<
        fdp::
            PostSnesPhaseTransitionControllerReport3D>
        disappearance_report;
    error =
        fdp::
            solve_nonlinear_timestep_with_phase_transitions_3d(
                PETSC_COMM_WORLD,
                std::move(appeared_system),
                {
                    &controller_scan,
                    &fixture,
                    &controller_rebuild,
                    &fixture},
                {4U},
                &disappeared_system,
                &disappeared_state,
                &disappearance_report);

    require_collective(
        error == PETSC_SUCCESS &&
            disappeared_system != nullptr &&
            disappeared_state != nullptr &&
            disappearance_report.has_value() &&
            disappearance_report->outcome ==
                fdp::
                    PostSnesPhaseTransitionOutcome3D::
                        stable_phase_set &&
            disappearance_report
                    ->transition_restarts ==
                1U &&
            fixture.rebuild_calls ==
                1U &&
            fixture.multi_well_rebound &&
            fixture.multi_well_context !=
                nullptr &&
            disappeared_system
                    ->numbering()
                    .petsc_global_scalar_count() ==
                42 &&
            disappeared_system
                    ->numbering()
                    .cell(
                        mesh::LocalIndex{2U})
                    .phase_count ==
                2U &&
            disappeared_system
                    ->numbering()
                    .cell(
                        mesh::LocalIndex{2U})
                    .scalar_count ==
                7U &&
            disappeared_system
                    ->numbering()
                    .cell(
                        mesh::LocalIndex{5U})
                    .phase_count ==
                3U &&
            disappeared_system
                    ->numbering()
                    .cell(
                        mesh::LocalIndex{5U})
                    .scalar_count ==
                10U,
        "multi-connection well did not survive local cell30 3P-to-2P restart");

    const auto* disappeared30 =
        fixture.multi_well_context
            ->find_connection(
                mesh::GlobalEntityId{
                    UINT64_C(30)});
    frozen60 =
        fixture.multi_well_context
            ->find_connection(
                mesh::GlobalEntityId{
                    UINT64_C(60)});
    const auto initial30_map =
        phase_identity_map_for_stable(
            UINT64_C(30));
    bool cell30_roundtrip_identity_ok =
        disappeared30 != nullptr &&
        disappeared30
                ->active_phase_identities()
                .has_value() &&
        disappeared30
                ->active_phase_identities()
                ->phase_count() ==
            initial30_map.phase_count();
    if (cell30_roundtrip_identity_ok) {
        for (std::size_t phase = 0U;
             phase < initial30_map.phase_count();
             ++phase) {
            cell30_roundtrip_identity_ok =
                cell30_roundtrip_identity_ok &&
                disappeared30
                        ->active_phase_identities()
                        ->identity(phase) ==
                    initial30_map.identity(
                        phase);
        }
    }

    const auto& retained_registry =
        disappeared30 != nullptr
            ? disappeared30
                  ->phase_identity_injection_enthalpy()
            : std::optional<
                  wdp::
                      FixedBhpPhaseIdentityInjectionEnthalpy3D>{};
    const bool disappeared_phase_retained =
        retained_registry.has_value() &&
        retained_registry->phases.size() ==
            3U &&
        std::any_of(
            retained_registry
                ->phases.begin(),
            retained_registry
                ->phases.end(),
            [&](const auto& entry) {
                return entry.identity ==
                           hydrocarbon1 &&
                    entry
                            .specific_enthalpy_j_per_kg ==
                        3000.0;
            });

    require_collective(
        cell30_roundtrip_identity_ok &&
            disappeared_phase_retained &&
            disappeared30
                    ->injection_enthalpy()
                    .specific_enthalpy_j_per_kg ==
                std::vector<double>{
                    1000.0,
                    2000.0} &&
            frozen60 != nullptr &&
            frozen60
                    ->target_cell_global() ==
                mesh::GlobalEntityId{
                    UINT64_C(60)} &&
            std::abs(
                frozen60
                        ->connection()
                        .well_index_m3 -
                    frozen60_wi) <=
                1.0e-14 *
                    std::max(
                        1.0,
                        std::abs(
                            frozen60_wi)) &&
            frozen60
                    ->bottom_hole_pressure_pa() ==
                shared_bhp_pa &&
            frozen60
                    ->injection_enthalpy()
                    .specific_enthalpy_j_per_kg ==
                frozen60_enthalpy &&
            frozen60
                    ->source_provenance() ==
                frozen60_source &&
            !frozen60
                 ->phase_identity_rebindable(),
        "cell30 local disappearance drifted identity or mutated frozen cell60 connection");

    run_rebound_multi_connection_fixed_bhp_physical_timestep(
        rank,
        &disappeared_system,
        &disappeared_state,
        &fixture,
        2U,
        7U,
        42);

    audit->well_timestep_compressibility =
        false;
}

