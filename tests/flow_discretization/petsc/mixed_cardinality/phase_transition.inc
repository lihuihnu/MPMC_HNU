std::vector<double>
controller_transitioned_3p_state() {
    auto q =
        target_3p();
    const auto layout =
        layout_3p();
    q[layout.pressure_unknown_index()] =
        20.0;
    q[layout.temperature_unknown_index()] =
        9.0;
    return q;
}

std::vector<
    fdp::FrozenPhaseTransitionRebuildCell3D>
make_controller_cells(
    int rank,
    bool cell30_is_three_phase,
    DispatchAudit* audit) {
    auto base_inputs =
        make_cell_inputs(
            rank,
            audit);
    auto phase_maps =
        make_phase_identity_maps();

    std::vector<
        fdp::
            FrozenPhaseTransitionRebuildCell3D>
        result;
    result.reserve(
        base_inputs.size());

    for (std::size_t local = 0U;
         local < base_inputs.size();
         ++local) {
        const std::uint64_t stable =
            static_cast<std::uint64_t>(
                (local + 1U) * 10U);
        fdp::FrozenPhaseTransitionRebuildCell3D
            snapshot;
        std::visit(
            [&](const auto& typed) {
                using Typed =
                    std::remove_cvref_t<
                        decltype(typed)>;
                snapshot.cell =
                    typed.cell;
                snapshot.cell_global =
                    typed.cell_global;
                snapshot.bulk_volume_m3 =
                    typed.bulk_volume_m3;
                snapshot.porosity =
                    typed.porosity;
                snapshot.component_ids =
                    typed.component_ids;
                if constexpr (
                    std::is_same_v<
                        Typed,
                        fdp::
                            FixedThreePhaseSnesCellInput3D>) {
                    snapshot.target_layout =
                        flow::
                            NaturalVariableLayoutDescriptor{
                                typed.frozen_layout};
                } else {
                    snapshot.target_layout =
                        typed.frozen_layout
                            .descriptor();
                }
                snapshot.previous_component_accumulation =
                    typed.previous_component_accumulation;
                snapshot.previous_energy_accumulation =
                    typed.previous_energy_accumulation;
            },
            base_inputs[local]);

        snapshot.target_natural_variables =
            target_state(stable);
        snapshot.target_active_phases =
            phase_maps[local];
        snapshot.transition_evidence_profile =
            "fixture/controller-frozen-active-set/v1";

        if (stable == UINT64_C(40)) {
            snapshot.frozen_absent_phases.push_back(
                make_frozen_absent_coordinate_entry(
                    stable,
                    mixed_physical_phase_identity(
                        "hydrocarbon-1"),
                    {0.10, 0.50, 0.40},
                    "fixture/hydrocarbon1-root0",
                    audit));
        }

        result.push_back(
            std::move(snapshot));
    }

    if (!cell30_is_three_phase) {
        return result;
    }

    constexpr std::size_t local = 2U;
    const auto q =
        controller_transitioned_3p_state();
    const auto layout =
        layout_3p();

    // This controller-policy fixture isolates outer orchestration from the
    // already-covered BE history-migration contract.  Seed the rebuilt 3P
    // previous state from the same controlled 3P closure so the restarted
    // SNES is an exact zero-residual solve; make_accepted_phase_transition_
    // rebuild_cell_3d above separately verifies preservation of real 2P
    // histories through an accepted transition.
    std::optional<
        fdp::
            FixedThreePhaseCurrentCellLinearization3D>
        controlled_three_phase;
    fdp::NaturalVariableSnesEvaluationStatus3D
        controlled_status =
            fdp::
                NaturalVariableSnesEvaluationStatus3D::
                    success;
    const std::vector<std::string>
        controlled_ids{
            "A", "B", "C"};
    const PetscErrorCode controlled_error =
        evaluate_3p(
            mesh::LocalIndex{
                static_cast<
                    mesh::LocalIndex::value_type>(
                        local)},
            mesh::GlobalEntityId{
                UINT64_C(30)},
            q,
            layout,
            controlled_ids,
            audit,
            &controlled_three_phase,
            &controlled_status);
    if (controlled_error != PETSC_SUCCESS ||
        controlled_status !=
            fdp::
                NaturalVariableSnesEvaluationStatus3D::
                    success ||
        !controlled_three_phase.has_value()) {
        throw std::runtime_error(
            "failed to build controlled 3P controller state");
    }

    fdp::FrozenPhaseTransitionRebuildCell3D
        transitioned;
    transitioned.cell =
        mesh::LocalIndex{
            static_cast<
                mesh::LocalIndex::value_type>(
                    local)};
    transitioned.cell_global =
        mesh::GlobalEntityId{
            UINT64_C(30)};
    transitioned.bulk_volume_m3 =
        4.0;
    transitioned.porosity =
        0.26;
    transitioned.component_ids =
        {"A", "B", "C"};
    transitioned.target_layout =
        flow::NaturalVariableLayoutDescriptor{
            layout};
    transitioned.target_natural_variables =
        q;
    transitioned.target_active_phases =
        flow::FrozenActivePhaseIdentityMap{
            {
                mixed_physical_phase_identity(
                    "aqueous"),
                mixed_physical_phase_identity(
                    "hydrocarbon-0"),
                mixed_physical_phase_identity(
                    "hydrocarbon-1")}};
    transitioned.transition_evidence_profile =
        "fixture/controller-accepted-2p-to-3p/v1";
    if (rank == 0) {
        transitioned
            .previous_component_accumulation =
            flow::
                build_pore_volume_component_accumulation(
                    controlled_three_phase->state,
                    transitioned.porosity);
        transitioned
            .previous_energy_accumulation =
            flow::
                build_pore_volume_energy_accumulation_snapshot(
                    controlled_three_phase->state,
                    transitioned.porosity,
                    controlled_three_phase->transport,
                    controlled_three_phase->caloric,
                    controlled_three_phase->rock);
    }

    result[local] =
        std::move(
            transitioned);
    return result;
}

struct ControllerFixture {
    ControllerFixture(
        int input_rank,
        const dp::
            ParallelOwnedConnectionSchedule3D*
                input_schedule,
        const mesh::PartitionSnapshot*
            input_partition,
        const dp::
            PetscMpiAijSymbolicPreallocation3D*
                input_bridge,
        const dp::
            OwnedCellStructuralColumnPatternSnapshot3D*
                input_pattern,
        DispatchAudit* input_audit,
        const th::Pr76Phase<double>*
            input_pr_model,
        flow::
            Pr76AbsentPhasePotentialExtensionProvider<
                double>*
                input_provider,
        bool input_oscillate_after_restart,
        std::size_t input_rebuild_calls,
        bool input_scan_indeterminate =
            false,
        bool input_phase_disappearance_only =
            false)
        : rank(input_rank),
          schedule(input_schedule),
          partition(input_partition),
          bridge(input_bridge),
          pattern(input_pattern),
          audit(input_audit),
          pr_model(input_pr_model),
          provider(input_provider),
          oscillate_after_restart(
              input_oscillate_after_restart),
          rebuild_calls(
              input_rebuild_calls),
          scan_indeterminate(
              input_scan_indeterminate),
          phase_disappearance_only(
              input_phase_disappearance_only) {}

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
    const th::Pr76Phase<double>* pr_model{};
    flow::
        Pr76AbsentPhasePotentialExtensionProvider<
            double>*
            provider{};
    bool oscillate_after_restart{};
    std::size_t rebuild_calls{};
    bool scan_indeterminate{};
    bool phase_disappearance_only{};

    wdp::
        FixedBhpPeacemanWellSourceEvaluatorContext3D*
            well_context{};
    std::optional<
        wdp::
            FixedBhpPeacemanWellSourceEvaluatorContext3D>
        rebound_well_context;
    FixedBhpWellSourceAudit
        well_source_audit{};
    bool well_rebound{};

    wdp::
        FixedBhpMultiConnectionWellSourceEvaluatorContext3D*
            multi_well_context{};
    std::optional<
        wdp::
            FixedBhpMultiConnectionWellSourceEvaluatorContext3D>
        rebound_multi_well_context;
    FixedBhpMultiConnectionWellSourceAudit
        multi_well_source_audit{};
    bool multi_well_rebound{};
};

fdp::PostSnesPhaseTransitionProposal3D
controller_two_to_three_proposal() {
    flow::PhaseSetTransitionCandidate
        candidate;
    candidate.status =
        flow::
            PhaseSetTransitionCandidateStatus::
                target_resolved;
    candidate.component_ids =
        {"A", "B", "C"};
    candidate.pressure_pa = 20.0;
    candidate.temperature_k = 9.0;
    candidate.evidence_profile =
        "fixture/post-snes-scan/v1";
    candidate.source_phase_count = 2U;
    candidate.target_phase_count = 3U;
    candidate.trigger =
        flow::
            PhaseSetTransitionTrigger::
                stability_witness;
    candidate.target_phases = {
        {0.3, {0.3, 0.3, 0.4}, 2.0},
        {0.3, {0.3, 0.3, 0.4}, 4.0},
        {0.4, {0.3, 0.3, 0.4}, 7.0}};
    return {
        mesh::GlobalEntityId{
            UINT64_C(30)},
        std::move(candidate)};
}

PetscErrorCode
controller_scan(
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
    *scan_status =
        fdp::PostSnesPhaseTransitionScanStatus3D::
            complete;
    output->clear();
    auto* fixture =
        static_cast<
            ControllerFixture*>(
                raw_context);

    if (static_cast<int>(
            solve_report.converged_reason) <=
        0) {
        return PETSC_ERR_ARG_WRONGSTATE;
    }
    if (fixture->scan_indeterminate) {
        *scan_status =
            fdp::PostSnesPhaseTransitionScanStatus3D::
                indeterminate;
        return PETSC_SUCCESS;
    }

    const auto& record =
        system.numbering().cell(
            mesh::LocalIndex{2U});
    if (record.owner_rank !=
        system.numbering()
            .local_rank()) {
        return PETSC_SUCCESS;
    }

    if (record.phase_count == 2U) {
        if (!fixture
                 ->phase_disappearance_only) {
            output->push_back(
                controller_two_to_three_proposal());
        }
        return PETSC_SUCCESS;
    }

    if (record.phase_count == 3U &&
        (fixture
             ->oscillate_after_restart ||
         fixture
             ->phase_disappearance_only)) {
        flow::PhaseSetTransitionCandidate
            candidate;
        candidate.status =
            flow::
                PhaseSetTransitionCandidateStatus::
                    target_resolved;
        candidate.component_ids =
            {"A", "B", "C"};
        candidate.pressure_pa = 20.0;
        candidate.temperature_k = 9.0;
        candidate.evidence_profile =
            "fixture/post-snes-scan/v1";
        candidate.source_phase_count = 3U;
        candidate.target_phase_count = 2U;
        candidate.trigger =
            flow::
                PhaseSetTransitionTrigger::
                    phase_disappearance;
        candidate.target_phases = {
            {0.5, {0.3, 0.3, 0.4}, 4.0},
            {0.5, {0.3, 0.3, 0.4}, 4.0}};
        output->push_back(
            {
                record.cell_global,
                std::move(candidate)});
    }
    return PETSC_SUCCESS;
}

PetscErrorCode
controller_rebuild(
    const fdp::
        PhaseTransitionRebuiltNaturalVariableSystem3D&,
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
            ControllerFixture*>(
                raw_context);

    if (accepted_global_batch.size() !=
            1U ||
        accepted_global_batch[0]
                .cell_global !=
            mesh::GlobalEntityId{
                UINT64_C(30)}) {
        return PETSC_ERR_ARG_INCOMP;
    }
    const bool target_three_phase =
        accepted_global_batch[0]
                .source_phase_count ==
            2U &&
        accepted_global_batch[0]
                .target_phase_count ==
            3U;
    const bool target_two_phase =
        accepted_global_batch[0]
                .source_phase_count ==
            3U &&
        accepted_global_batch[0]
                .target_phase_count ==
            2U;
    if (!target_three_phase &&
        !target_two_phase) {
        return PETSC_ERR_ARG_INCOMP;
    }

    if (fixture->rank == 0) {
        if (local_owned_proposals.size() !=
                1U ||
            local_owned_proposals[0]
                    .cell_global !=
                mesh::GlobalEntityId{
                    UINT64_C(30)} ||
            local_owned_proposals[0]
                    .candidate
                    .source_phase_count !=
                accepted_global_batch[0]
                    .source_phase_count ||
            local_owned_proposals[0]
                    .candidate
                    .target_phase_count !=
                accepted_global_batch[0]
                    .target_phase_count) {
            return PETSC_ERR_ARG_INCOMP;
        }
    } else if (
        !local_owned_proposals.empty()) {
        return PETSC_ERR_ARG_INCOMP;
    }

    ++fixture->rebuild_calls;
    auto cells =
        make_controller_cells(
            fixture->rank,
            target_three_phase,
            fixture->audit);
    auto faces =
        make_face_inputs(
            fixture->rank,
            false);

    fdp::
        MixedCardinalityPhysicalCellSourceEvaluatorBinding3D
        source_binding{};
    if (fixture->well_context != nullptr) {
        const auto target =
            std::find_if(
                cells.begin(),
                cells.end(),
                [&](const auto& cell) {
                    return cell.cell_global ==
                        fixture
                            ->well_context
                            ->target_cell_global();
                });
        if (target == cells.end()) {
            return PETSC_ERR_ARG_INCOMP;
        }
        try {
            fixture->rebound_well_context.emplace(
                wdp::
                    rebind_fixed_bhp_peaceman_well_source_context_3d(
                        *fixture->well_context,
                        target->cell_global,
                        target->target_active_phases));
        } catch (...) {
            return PETSC_ERR_ARG_INCOMP;
        }
        fixture->well_context =
            &*fixture->rebound_well_context;
        fixture->well_source_audit.context =
            fixture->well_context;
        fixture->well_rebound = true;
        source_binding = {
            &evaluate_audited_fixed_bhp_well_source,
            &fixture->well_source_audit};
    } else if (
        fixture->multi_well_context != nullptr) {
        const auto transitioned =
            accepted_global_batch.front()
                .cell_global;
        const auto target =
            std::find_if(
                cells.begin(),
                cells.end(),
                [&](const auto& cell) {
                    return cell.cell_global ==
                        transitioned;
                });
        if (target == cells.end()) {
            return PETSC_ERR_ARG_INCOMP;
        }
        try {
            fixture
                ->rebound_multi_well_context
                .emplace(
                    wdp::
                        rebind_fixed_bhp_multi_connection_well_source_context_3d(
                            *fixture
                                 ->multi_well_context,
                            transitioned,
                            target
                                ->target_active_phases));
        } catch (...) {
            return PETSC_ERR_ARG_INCOMP;
        }
        fixture->multi_well_context =
            &*fixture
                 ->rebound_multi_well_context;
        fixture
            ->multi_well_source_audit
            .context =
            fixture->multi_well_context;
        fixture->multi_well_rebound =
            true;
        source_binding = {
            &evaluate_audited_fixed_bhp_multi_connection_well_source,
            &fixture
                 ->multi_well_source_audit};
    }

    return fdp::
        rebuild_phase_transition_natural_variable_system_3d(
            PETSC_COMM_WORLD,
            *fixture->schedule,
            *fixture->partition,
            *fixture->bridge,
            *fixture->pattern,
            1.0,
            std::move(cells),
            std::move(faces),
            {
                {&evaluate_1p, fixture->audit},
                {&evaluate_2p, fixture->audit},
                {&evaluate_3p, fixture->audit}},
            source_binding,
            &fdp::
                evaluate_absent_phase_thermodynamic_provider_3d<
                    flow::
                        Pr76AbsentPhasePotentialExtensionProvider<
                            double>>,
            fixture->provider,
            output);
}

std::unique_ptr<
    fdp::PhaseTransitionRebuiltNaturalVariableSystem3D>
make_controller_initial_system(
    ControllerFixture* fixture) {
    auto cells =
        make_controller_cells(
            fixture->rank,
            false,
            fixture->audit);
    auto faces =
        make_face_inputs(
            fixture->rank,
            false);
    fdp::
        MixedCardinalityPhysicalCellSourceEvaluatorBinding3D
        source_binding{};
    if (fixture->well_context != nullptr) {
        fixture->well_source_audit.context =
            fixture->well_context;
        source_binding = {
            &evaluate_audited_fixed_bhp_well_source,
            &fixture->well_source_audit};
    } else if (
        fixture->multi_well_context != nullptr) {
        fixture
            ->multi_well_source_audit
            .context =
            fixture->multi_well_context;
        source_binding = {
            &evaluate_audited_fixed_bhp_multi_connection_well_source,
            &fixture
                 ->multi_well_source_audit};
    }

    std::unique_ptr<
        fdp::
            PhaseTransitionRebuiltNaturalVariableSystem3D>
        system;
    const PetscErrorCode error =
        fdp::
            rebuild_phase_transition_natural_variable_system_3d(
                PETSC_COMM_WORLD,
                *fixture->schedule,
                *fixture->partition,
                *fixture->bridge,
                *fixture->pattern,
                1.0,
                std::move(cells),
                std::move(faces),
                {
                    {&evaluate_1p, fixture->audit},
                    {&evaluate_2p, fixture->audit},
                    {&evaluate_3p, fixture->audit}},
                source_binding,
                &fdp::
                    evaluate_absent_phase_thermodynamic_provider_3d<
                        flow::
                            Pr76AbsentPhasePotentialExtensionProvider<
                                double>>,
                fixture->provider,
                &system);
    if (error != PETSC_SUCCESS ||
        system == nullptr) {
        throw std::runtime_error(
            "failed to build controller initial system");
    }
    return system;
}

void run_controller_case(
    ControllerFixture* fixture,
    std::size_t max_restarts,
    fdp::PostSnesPhaseTransitionOutcome3D
        expected_outcome,
    std::size_t expected_restarts) {
    fixture->rebuild_calls = 0U;
    auto initial =
        make_controller_initial_system(
            fixture);

    std::unique_ptr<
        fdp::
            PhaseTransitionRebuiltNaturalVariableSystem3D>
        final_system;
    Vec final_state = nullptr;
    std::optional<
        fdp::
            PostSnesPhaseTransitionControllerReport3D>
        controller_report;
    const PetscErrorCode error =
        fdp::
            solve_nonlinear_timestep_with_phase_transitions_3d(
                PETSC_COMM_WORLD,
                std::move(initial),
                {
                    &controller_scan,
                    fixture,
                    &controller_rebuild,
                    fixture},
                {max_restarts},
                &final_system,
                &final_state,
                &controller_report);

    require_collective(
        error == PETSC_SUCCESS &&
            final_system != nullptr &&
            final_state != nullptr &&
            controller_report.has_value() &&
            controller_report->outcome ==
                expected_outcome &&
            controller_report
                    ->transition_restarts ==
                expected_restarts,
        "post-SNES phase-transition controller outcome mismatch");

    const bool expected_accepted =
        expected_outcome ==
        fdp::
            PostSnesPhaseTransitionOutcome3D::
                stable_phase_set;
    require_collective(
        controller_report
                ->timestep_accepted() ==
            expected_accepted,
        "post-SNES controller timestep acceptance semantics mismatch");

    if (expected_outcome ==
        fdp::
            PostSnesPhaseTransitionOutcome3D::
                stable_phase_set) {
        require_collective(
            controller_report
                    ->generations.size() ==
                2U &&
                final_system
                        ->numbering()
                        .cell(
                            mesh::LocalIndex{2U})
                        .phase_count ==
                    3U &&
                final_system
                        ->numbering()
                        .petsc_global_scalar_count() ==
                    45,
            "stable controller path did not finish on rebuilt 3P phase set");
    } else if (
        expected_outcome ==
        fdp::
            PostSnesPhaseTransitionOutcome3D::
                phase_set_cycle_detected) {
        require_collective(
            controller_report
                    ->generations.size() ==
                2U &&
                controller_report
                    ->generations[1]
                    .accepted_transition_batch
                    .size() ==
                1U &&
                fixture->rebuild_calls ==
                    2U,
            "physical phase-set cycle was not detected after the reverse rebuild");
    } else if (
        expected_outcome ==
        fdp::
            PostSnesPhaseTransitionOutcome3D::
                phase_set_scan_indeterminate) {
        require_collective(
            controller_report
                    ->generations.size() ==
                1U &&
                controller_report
                    ->generations.front()
                    .phase_transition_scan_status ==
                    fdp::
                        PostSnesPhaseTransitionScanStatus3D::
                            indeterminate &&
                fixture->rebuild_calls ==
                    0U,
            "indeterminate phase-set scan was incorrectly accepted or rebuilt");
    } else if (
        expected_outcome ==
        fdp::
            PostSnesPhaseTransitionOutcome3D::
                transition_restart_budget_exhausted) {
        require_collective(
            controller_report
                    ->generations.size() ==
                1U &&
                fixture->rebuild_calls ==
                    0U,
            "transition budget exhausted path rebuilt the system");
    }

    require_collective(
        VecDestroy(&final_state) ==
            PETSC_SUCCESS,
        "post-SNES controller final state cleanup failed");
}

[[nodiscard]] bool
local_solution_matches_target(
    const fdp::
        VariableCardinalityNaturalVariableSnesSolveReport3D&
            report) {
    for (const auto& entry :
         report.locally_owned_solution) {
        const auto target =
            target_state(
                entry.cell_global.value());
        if (entry.natural_variable_slot >=
                target.size() ||
            std::abs(
                entry.value -
                target[
                    entry.natural_variable_slot]) >
                1.0e-10 *
                    std::max(
                        1.0,
                        std::abs(
                            target[
                                entry.natural_variable_slot]))) {
            return false;
        }
    }
    return true;
}


