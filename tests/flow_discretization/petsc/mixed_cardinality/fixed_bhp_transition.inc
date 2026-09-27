void run_phase_transition_rebound_fixed_bhp_case(
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
            "invalid phase-transition fixed-BHP well fixture");
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
        phase_enthalpy_registry{
            "fixture/fixed-bhp-phase-identity-enthalpy/v1",
            {
                {aqueous, 1100.0},
                {hydrocarbon0, 2200.0},
                {hydrocarbon1, 3300.0}}};

    // Prove the registry is identity-driven rather than slot-driven before it
    // enters the nonlinear system: a deliberately reordered active map must
    // reorder the resolved enthalpy payload accordingly.
    const auto reordered_enthalpy =
        wdp::
            resolve_fixed_bhp_injection_enthalpy_by_phase_identity_3d(
                phase_enthalpy_registry,
                flow::FrozenActivePhaseIdentityMap{
                    {hydrocarbon0, aqueous}});
    require_collective(
        reordered_enthalpy
                .specific_enthalpy_j_per_kg
                .size() ==
            2U &&
            reordered_enthalpy
                    .specific_enthalpy_j_per_kg[
                        0U] ==
                2200.0 &&
            reordered_enthalpy
                    .specific_enthalpy_j_per_kg[
                        1U] ==
                1100.0,
        "fixed-BHP injection enthalpy followed slot order instead of stable phase identity");

    auto well_context =
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
                    15.0,
                    phase_enthalpy_registry,
                    phase_identity_map_for_stable(
                        UINT64_C(30)),
                    "fixture/fixed-bhp-transition-cell30/v1");

    require_collective(
        well_context.phase_identity_rebindable() &&
            well_context
                    .active_phase_identities()
                    .has_value() &&
            well_context
                    .active_phase_identities()
                    ->phase_count() ==
                2U &&
            well_context
                    .injection_enthalpy()
                    .specific_enthalpy_j_per_kg ==
                std::vector<double>{
                    1100.0,
                    2200.0},
        "fixed-BHP initial 2P phase-identity binding is malformed");

    const double initial_well_index =
        well_context.connection()
            .well_index_m3;
    const double initial_bhp =
        well_context
            .bottom_hole_pressure_pa();

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
        false};
    fixture.well_context =
        &well_context;

    auto initial_system =
        make_controller_initial_system(
            &fixture);
    require_collective(
        initial_system != nullptr &&
            initial_system
                    ->numbering()
                    .cell(
                        mesh::LocalIndex{2U})
                    .cell_global ==
                mesh::GlobalEntityId{
                    UINT64_C(30)} &&
            initial_system
                    ->numbering()
                    .cell(
                        mesh::LocalIndex{2U})
                    .phase_count ==
                2U &&
            initial_system
                    ->coordinate_registry()
                    .cell(
                        mesh::GlobalEntityId{
                            UINT64_C(30)})
                    .active_phases
                    .phase_count() ==
                2U,
        "fixed-BHP transition fixture did not start on the stable 2P completion");

    std::unique_ptr<
        fdp::
            PhaseTransitionRebuiltNaturalVariableSystem3D>
        stable_system;
    Vec stable_state = nullptr;
    std::optional<
        fdp::
            PostSnesPhaseTransitionControllerReport3D>
        transition_report;
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
                &stable_system,
                &stable_state,
                &transition_report);

    require_collective(
        error == PETSC_SUCCESS &&
            stable_system != nullptr &&
            stable_state != nullptr &&
            transition_report.has_value() &&
            transition_report->outcome ==
                fdp::
                    PostSnesPhaseTransitionOutcome3D::
                        stable_phase_set &&
            transition_report
                    ->transition_restarts ==
                1U &&
            transition_report
                    ->generations.size() ==
                2U &&
            fixture.rebuild_calls ==
                1U &&
            fixture.well_rebound &&
            fixture.well_context !=
                nullptr &&
            stable_system
                    ->numbering()
                    .cell(
                        mesh::LocalIndex{2U})
                    .phase_count ==
                3U &&
            stable_system
                    ->numbering()
                    .cell(
                        mesh::LocalIndex{2U})
                    .scalar_count ==
                10U,
        "fixed-BHP completion did not survive the accepted 2P-to-3P rebuild and SNES restart");

    const auto& rebound_context =
        *fixture.well_context;
    const auto& rebuilt_phase_map =
        stable_system
            ->coordinate_registry()
            .cell(
                mesh::GlobalEntityId{
                    UINT64_C(30)})
            .active_phases;
    bool rebound_identity_matches =
        rebound_context
                .active_phase_identities()
                .has_value() &&
        rebound_context
                .active_phase_identities()
                ->phase_count() ==
            rebuilt_phase_map.phase_count();
    if (rebound_identity_matches) {
        for (std::size_t phase = 0U;
             phase < rebuilt_phase_map.phase_count();
             ++phase) {
            rebound_identity_matches =
                rebound_identity_matches &&
                rebound_context
                        .active_phase_identities()
                        ->identity(phase) ==
                    rebuilt_phase_map.identity(
                        phase);
        }
    }

    require_collective(
        rebound_identity_matches &&
            rebound_context
                    .target_cell_global() ==
                mesh::GlobalEntityId{
                    UINT64_C(30)} &&
            std::abs(
                rebound_context
                        .connection()
                        .well_index_m3 -
                    initial_well_index) <=
                1.0e-14 *
                    std::max(
                        1.0,
                        std::abs(
                            initial_well_index)) &&
            rebound_context
                    .bottom_hole_pressure_pa() ==
                initial_bhp &&
            rebound_context
                    .injection_enthalpy()
                    .specific_enthalpy_j_per_kg ==
                std::vector<double>{
                    1100.0,
                    2200.0,
                    3300.0},
        "fixed-BHP completion rebinding changed stable cell/WI/BHP or failed stable phase-identity enthalpy mapping");

    std::vector<std::optional<
        fdp::
            MixedCardinalityPhysicalCurrentCellLinearization3D>>
        rebound_current;
    std::vector<double>
        rebound_porosity;
    fdp::NaturalVariableSnesEvaluationStatus3D
        rebound_status =
            fdp::
                NaturalVariableSnesEvaluationStatus3D::
                    success;
    require_collective(
        stable_system
                ->evaluate_local_cells_for_phase_transition(
                    stable_state,
                    &rebound_current,
                    &rebound_porosity,
                    &rebound_status) ==
            PETSC_SUCCESS &&
            rebound_status ==
                fdp::
                    NaturalVariableSnesEvaluationStatus3D::
                        success &&
            rebound_current.size() >
                2U &&
            rebound_current[2U]
                .has_value(),
        "failed to evaluate rebuilt fixed-BHP completion state");

    const auto rebound_well =
        wdp::
            build_fixed_bhp_peaceman_well_source_3d(
                rebound_context,
                *rebound_current[2U]);
    require_collective(
        rebound_well
                .cell_source
                .input_count ==
            10U,
        "fixed-BHP rebound source did not rebuild its Jacobian on q=3Nc+1");

    // First prove the post-appearance 3P completion can run a real accepted
    // physical timestep with owner-only well source and conservation.
    run_rebound_fixed_bhp_physical_timestep(
        rank,
        &stable_system,
        &stable_state,
        &fixture,
        3U,
        10U);

    // Start a second, independent controller session from the already accepted
    // 3P system. This preserves the controller's cycle-detection contract: the
    // reverse 3P->2P transition is not a rollback inside the first session.
    fixture.phase_disappearance_only =
        true;
    fixture.oscillate_after_restart =
        false;
    fixture.rebuild_calls =
        0U;
    fixture.well_rebound =
        false;
    fixture.well_source_audit.evaluator_calls =
        0U;
    fixture.well_source_audit.target_calls =
        0U;

    std::unique_ptr<
        fdp::
            PhaseTransitionRebuiltNaturalVariableSystem3D>
        roundtrip_system;
    Vec roundtrip_state = nullptr;
    std::optional<
        fdp::
            PostSnesPhaseTransitionControllerReport3D>
        disappearance_report;
    error =
        fdp::
            solve_nonlinear_timestep_with_phase_transitions_3d(
                PETSC_COMM_WORLD,
                std::move(stable_system),
                {
                    &controller_scan,
                    &fixture,
                    &controller_rebuild,
                    &fixture},
                {4U},
                &roundtrip_system,
                &roundtrip_state,
                &disappearance_report);

    require_collective(
        error == PETSC_SUCCESS &&
            roundtrip_system != nullptr &&
            roundtrip_state != nullptr &&
            disappearance_report.has_value() &&
            disappearance_report->outcome ==
                fdp::
                    PostSnesPhaseTransitionOutcome3D::
                        stable_phase_set &&
            disappearance_report
                    ->transition_restarts ==
                1U &&
            disappearance_report
                    ->generations.size() ==
                2U &&
            fixture.rebuild_calls ==
                1U &&
            fixture.well_rebound &&
            fixture.well_context !=
                nullptr &&
            roundtrip_system
                    ->numbering()
                    .cell(
                        mesh::LocalIndex{2U})
                    .cell_global ==
                mesh::GlobalEntityId{
                    UINT64_C(30)} &&
            roundtrip_system
                    ->numbering()
                    .cell(
                        mesh::LocalIndex{2U})
                    .phase_count ==
                2U &&
            roundtrip_system
                    ->numbering()
                    .cell(
                        mesh::LocalIndex{2U})
                    .scalar_count ==
                7U &&
            roundtrip_system
                    ->numbering()
                    .petsc_global_scalar_count() ==
                42,
        "fixed-BHP completion did not survive the independent 3P-to-2P phase-disappearance restart");

    const auto initial_two_phase_map =
        phase_identity_map_for_stable(
            UINT64_C(30));
    const auto& roundtrip_phase_map =
        roundtrip_system
            ->coordinate_registry()
            .cell(
                mesh::GlobalEntityId{
                    UINT64_C(30)})
            .active_phases;
    bool roundtrip_identity_matches =
        roundtrip_phase_map.phase_count() ==
            initial_two_phase_map.phase_count() &&
        fixture.well_context
                ->active_phase_identities()
                .has_value() &&
        fixture.well_context
                ->active_phase_identities()
                ->phase_count() ==
            initial_two_phase_map.phase_count();
    if (roundtrip_identity_matches) {
        for (std::size_t phase = 0U;
             phase <
             initial_two_phase_map.phase_count();
             ++phase) {
            roundtrip_identity_matches =
                roundtrip_identity_matches &&
                roundtrip_phase_map.identity(
                    phase) ==
                    initial_two_phase_map.identity(
                        phase) &&
                fixture.well_context
                        ->active_phase_identities()
                        ->identity(phase) ==
                    initial_two_phase_map.identity(
                        phase);
        }
    }

    const auto& retained_registry =
        fixture.well_context
            ->phase_identity_injection_enthalpy();
    const bool retained_disappeared_phase =
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
                        3300.0;
            });

    require_collective(
        roundtrip_identity_matches &&
            retained_disappeared_phase &&
            fixture.well_context
                    ->target_cell_global() ==
                mesh::GlobalEntityId{
                    UINT64_C(30)} &&
            std::abs(
                fixture.well_context
                        ->connection()
                        .well_index_m3 -
                    initial_well_index) <=
                1.0e-14 *
                    std::max(
                        1.0,
                        std::abs(
                            initial_well_index)) &&
            fixture.well_context
                    ->bottom_hole_pressure_pa() ==
                initial_bhp &&
            fixture.well_context
                    ->injection_enthalpy()
                    .specific_enthalpy_j_per_kg ==
                std::vector<double>{
                    1100.0,
                    2200.0},
        "fixed-BHP 2P->3P->2P round-trip drifted stable phase identity or discarded inactive enthalpy registry data");

    std::vector<std::optional<
        fdp::
            MixedCardinalityPhysicalCurrentCellLinearization3D>>
        roundtrip_current;
    std::vector<double>
        roundtrip_porosity;
    fdp::NaturalVariableSnesEvaluationStatus3D
        roundtrip_status =
            fdp::
                NaturalVariableSnesEvaluationStatus3D::
                    success;
    require_collective(
        roundtrip_system
                ->evaluate_local_cells_for_phase_transition(
                    roundtrip_state,
                    &roundtrip_current,
                    &roundtrip_porosity,
                    &roundtrip_status) ==
            PETSC_SUCCESS &&
            roundtrip_status ==
                fdp::
                    NaturalVariableSnesEvaluationStatus3D::
                        success &&
            roundtrip_current.size() >
                2U &&
            roundtrip_current[2U]
                .has_value(),
        "failed to evaluate round-trip 2P fixed-BHP completion state");

    const auto roundtrip_well =
        wdp::
            build_fixed_bhp_peaceman_well_source_3d(
                *fixture.well_context,
                *roundtrip_current[2U]);
    require_collective(
        roundtrip_well
                .cell_source
                .input_count ==
            7U,
        "phase-disappearance rebinding did not shrink fixed-BHP Jacobian back to q=2Nc+1");

    run_rebound_fixed_bhp_physical_timestep(
        rank,
        &roundtrip_system,
        &roundtrip_state,
        &fixture,
        2U,
        7U);

    audit->well_timestep_compressibility =
        false;
}

