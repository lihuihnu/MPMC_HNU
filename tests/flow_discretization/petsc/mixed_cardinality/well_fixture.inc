struct FrozenWellTimestepControlAudit {
    mesh::LocalIndex target_cell{
        mesh::LocalIndex::value_type{0}};
    mesh::GlobalEntityId target_cell_global{
        mesh::GlobalEntityId::value_type{0}};
    std::size_t expected_phase_count{};
    std::size_t expected_scalar_count{};
    std::size_t scans{};
    std::size_t rebuild_calls{};
};

PetscErrorCode
frozen_well_timestep_scan(
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
    auto* audit =
        static_cast<
            FrozenWellTimestepControlAudit*>(
                raw_context);
    *scan_status =
        fdp::PostSnesPhaseTransitionScanStatus3D::
            complete;
    output->clear();

    if (static_cast<int>(
            solve_report.converged_reason) <=
        0) {
        return PETSC_ERR_ARG_WRONGSTATE;
    }

    const auto& well_cell =
        system.numbering().cell(
            audit->target_cell);
    if (well_cell.cell_global !=
            audit->target_cell_global ||
        well_cell.phase_count !=
            audit->expected_phase_count ||
        well_cell.scalar_count !=
            audit->expected_scalar_count) {
        return PETSC_ERR_ARG_INCOMP;
    }

    ++audit->scans;
    return PETSC_SUCCESS;
}

struct FrozenMultiWellTimestepControlAudit {
    std::size_t expected_cell30_phase_count{2U};
    std::size_t expected_cell30_scalar_count{7U};
    PetscInt expected_global_scalar_count{42};
    std::size_t scans{};
    std::size_t rebuild_calls{};
};

PetscErrorCode
frozen_multi_well_timestep_scan(
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
    auto* audit =
        static_cast<
            FrozenMultiWellTimestepControlAudit*>(
                raw_context);
    *scan_status =
        fdp::PostSnesPhaseTransitionScanStatus3D::
            complete;
    output->clear();

    if (static_cast<int>(
            solve_report.converged_reason) <=
        0) {
        return PETSC_ERR_ARG_WRONGSTATE;
    }

    const auto& cell30 =
        system.numbering().cell(
            mesh::LocalIndex{2U});
    const auto& cell60 =
        system.numbering().cell(
            mesh::LocalIndex{5U});
    if (system.numbering()
                .petsc_global_scalar_count() !=
            audit->expected_global_scalar_count ||
        cell30.cell_global !=
            mesh::GlobalEntityId{
                UINT64_C(30)} ||
        cell30.phase_count !=
            audit->expected_cell30_phase_count ||
        cell30.scalar_count !=
            audit->expected_cell30_scalar_count ||
        cell60.cell_global !=
            mesh::GlobalEntityId{
                UINT64_C(60)} ||
        cell60.phase_count != 3U ||
        cell60.scalar_count != 10U) {
        return PETSC_ERR_ARG_INCOMP;
    }

    ++audit->scans;
    return PETSC_SUCCESS;
}

PetscErrorCode
unexpected_frozen_multi_well_timestep_rebuild(
    const fdp::
        PhaseTransitionRebuiltNaturalVariableSystem3D&,
    Vec,
    const fdp::
        VariableCardinalityNaturalVariableSnesSolveReport3D&,
    std::span<
        const fdp::
            PostSnesPhaseTransitionProposal3D>,
    std::span<
        const fdp::
            AcceptedPhaseTransitionSummary3D>,
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
    auto* audit =
        static_cast<
            FrozenMultiWellTimestepControlAudit*>(
                raw_context);
    ++audit->rebuild_calls;
    return PETSC_ERR_PLIB;
}

PetscErrorCode
unexpected_frozen_well_timestep_rebuild(
    const fdp::
        PhaseTransitionRebuiltNaturalVariableSystem3D&,
    Vec,
    const fdp::
        VariableCardinalityNaturalVariableSnesSolveReport3D&,
    std::span<
        const fdp::
            PostSnesPhaseTransitionProposal3D>,
    std::span<
        const fdp::
            AcceptedPhaseTransitionSummary3D>,
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
    auto* audit =
        static_cast<
            FrozenWellTimestepControlAudit*>(
                raw_context);
    ++audit->rebuild_calls;
    return PETSC_ERR_PLIB;
}

std::unique_ptr<
    fdp::PhaseTransitionRebuiltNaturalVariableSystem3D>
make_frozen_fixed_bhp_timestep_system(
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
    DispatchAudit* dispatch_audit,
    FixedBhpWellSourceAudit* source_audit,
    flow::
        Pr76AbsentPhasePotentialExtensionProvider<
            double>* provider) {
    if (dispatch_audit == nullptr ||
        source_audit == nullptr ||
        source_audit->context == nullptr ||
        provider == nullptr) {
        throw std::invalid_argument(
            "invalid frozen fixed-BHP timestep fixture context");
    }

    auto cells =
        make_controller_cells(
            rank,
            false,
            dispatch_audit);
    auto faces =
        make_face_inputs(
            rank,
            false);

    std::unique_ptr<
        fdp::
            PhaseTransitionRebuiltNaturalVariableSystem3D>
        system;
    const PetscErrorCode error =
        fdp::
            rebuild_phase_transition_natural_variable_system_3d(
                PETSC_COMM_WORLD,
                schedule,
                partition,
                bridge,
                pattern,
                1.0,
                std::move(cells),
                std::move(faces),
                {
                    {&evaluate_1p, dispatch_audit},
                    {&evaluate_2p, dispatch_audit},
                    {&evaluate_3p, dispatch_audit}},
                {
                    &evaluate_audited_fixed_bhp_well_source,
                    source_audit},
                &fdp::
                    evaluate_absent_phase_thermodynamic_provider_3d<
                        flow::
                            Pr76AbsentPhasePotentialExtensionProvider<
                                double>>,
                provider,
                &system);
    if (error != PETSC_SUCCESS ||
        system == nullptr) {
        throw std::runtime_error(
            "failed to build frozen fixed-BHP physical timestep system");
    }
    return system;
}

std::unique_ptr<
    fdp::PhaseTransitionRebuiltNaturalVariableSystem3D>
make_frozen_multi_connection_fixed_bhp_timestep_system(
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
    DispatchAudit* dispatch_audit,
    FixedBhpMultiConnectionWellSourceAudit*
        source_audit,
    flow::
        Pr76AbsentPhasePotentialExtensionProvider<
            double>* provider) {
    if (dispatch_audit == nullptr ||
        source_audit == nullptr ||
        source_audit->context == nullptr ||
        provider == nullptr) {
        throw std::invalid_argument(
            "invalid multi-connection fixed-BHP timestep fixture context");
    }

    auto cells =
        make_controller_cells(
            rank,
            false,
            dispatch_audit);
    auto faces =
        make_face_inputs(
            rank,
            false);

    std::unique_ptr<
        fdp::
            PhaseTransitionRebuiltNaturalVariableSystem3D>
        system;
    const PetscErrorCode error =
        fdp::
            rebuild_phase_transition_natural_variable_system_3d(
                PETSC_COMM_WORLD,
                schedule,
                partition,
                bridge,
                pattern,
                1.0,
                std::move(cells),
                std::move(faces),
                {
                    {&evaluate_1p, dispatch_audit},
                    {&evaluate_2p, dispatch_audit},
                    {&evaluate_3p, dispatch_audit}},
                {
                    &evaluate_audited_fixed_bhp_multi_connection_well_source,
                    source_audit},
                &fdp::
                    evaluate_absent_phase_thermodynamic_provider_3d<
                        flow::
                            Pr76AbsentPhasePotentialExtensionProvider<
                                double>>,
                provider,
                &system);
    if (error != PETSC_SUCCESS ||
        system == nullptr) {
        throw std::runtime_error(
            "failed to build multi-connection fixed-BHP physical timestep system");
    }
    return system;
}

std::unique_ptr<
    fdp::PhaseTransitionRebuiltNaturalVariableSystem3D>
make_fixed_total_molar_rate_reservoir_system(
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
    DispatchAudit* dispatch_audit,
    wdp::
        FixedTotalMolarRateWellSourceEvaluatorContext3D*
            source_context,
    flow::
        Pr76AbsentPhasePotentialExtensionProvider<
            double>* provider,
    bool cell30_is_three_phase = false) {
    if (dispatch_audit == nullptr ||
        source_context == nullptr ||
        provider == nullptr) {
        throw std::invalid_argument(
            "invalid fixed-total-molar-rate reservoir fixture context");
    }

    auto cells =
        make_controller_cells(
            rank,
            cell30_is_three_phase,
            dispatch_audit);
    auto faces =
        make_face_inputs(
            rank,
            false);

    std::unique_ptr<
        fdp::
            PhaseTransitionRebuiltNaturalVariableSystem3D>
        system;
    const PetscErrorCode error =
        fdp::
            rebuild_phase_transition_natural_variable_system_3d(
                PETSC_COMM_WORLD,
                schedule,
                partition,
                bridge,
                pattern,
                1.0,
                std::move(cells),
                std::move(faces),
                {
                    {&evaluate_1p, dispatch_audit},
                    {&evaluate_2p, dispatch_audit},
                    {&evaluate_3p, dispatch_audit}},
                wdp::
                    fixed_total_molar_rate_well_source_binding_3d(
                        source_context),
                &fdp::
                    evaluate_absent_phase_thermodynamic_provider_3d<
                        flow::
                            Pr76AbsentPhasePotentialExtensionProvider<
                                double>>,
                provider,
                &system);
    if (error != PETSC_SUCCESS ||
        system == nullptr) {
        throw std::runtime_error(
            "failed to build fixed-total-molar-rate reservoir system");
    }
    return system;
}

std::array<double, 4>
owned_conserved_totals(
    const fdp::
        PhaseTransitionRebuiltNaturalVariableSystem3D&
            system,
    Vec state) {
    std::vector<std::optional<
        fdp::
            MixedCardinalityPhysicalCurrentCellLinearization3D>>
        current;
    std::vector<double>
        porosities;
    fdp::NaturalVariableSnesEvaluationStatus3D
        status =
            fdp::
                NaturalVariableSnesEvaluationStatus3D::
                    success;
    const PetscErrorCode error =
        system.evaluate_local_cells_for_phase_transition(
            state,
            &current,
            &porosities,
            &status);
    if (error != PETSC_SUCCESS ||
        status !=
            fdp::
                NaturalVariableSnesEvaluationStatus3D::
                    success ||
        current.size() != porosities.size()) {
        throw std::runtime_error(
            "failed to evaluate conserved totals for frozen well timestep");
    }

    std::array<double, 4>
        totals{};
    for (std::size_t local = 0U;
         local < current.size();
         ++local) {
        const auto& record =
            system.numbering().cell(
                mesh::LocalIndex{
                    static_cast<
                        mesh::LocalIndex::value_type>(
                            local)});
        if (record.owner_rank !=
            system.numbering().local_rank()) {
            continue;
        }
        if (!current[local].has_value()) {
            throw std::runtime_error(
                "owned frozen well timestep cell has no current state");
        }

        const double bulk_volume_m3 =
            2.0 +
            static_cast<double>(local);
        std::visit(
            [&](const auto& typed) {
                using Typed =
                    std::remove_cvref_t<
                        decltype(typed)>;

                if constexpr (
                    std::is_same_v<
                        Typed,
                        fdp::
                            SinglePhaseCurrentCellLinearization3D>) {
                    const auto component =
                        flow::
                            build_single_phase_component_accumulation(
                                typed.state,
                                porosities[local]);
                    const auto energy =
                        flow::
                            build_single_phase_energy_accumulation_snapshot(
                                typed.state,
                                porosities[local],
                                typed.transport,
                                typed.caloric,
                                typed.rock);
                    for (std::size_t i = 0U;
                         i < 3U;
                         ++i) {
                        totals[i] +=
                            component
                                .component_accumulation_mol_per_bulk_m3[i] *
                            bulk_volume_m3;
                    }
                    totals[3] +=
                        energy
                            .total_internal_energy_j_per_bulk_m3 *
                        bulk_volume_m3;
                } else if constexpr (
                    std::is_same_v<
                        Typed,
                        fdp::
                            TwoPhaseCurrentCellLinearization3D>) {
                    const auto component =
                        flow::
                            build_two_phase_component_accumulation(
                                typed.state,
                                porosities[local]);
                    const auto energy =
                        flow::
                            build_two_phase_energy_accumulation_snapshot(
                                typed.state,
                                porosities[local],
                                typed.transport,
                                typed.caloric,
                                typed.rock);
                    for (std::size_t i = 0U;
                         i < 3U;
                         ++i) {
                        totals[i] +=
                            component
                                .component_accumulation_mol_per_bulk_m3[i] *
                            bulk_volume_m3;
                    }
                    totals[3] +=
                        energy
                            .total_internal_energy_j_per_bulk_m3 *
                        bulk_volume_m3;
                } else {
                    const auto component =
                        flow::
                            build_pore_volume_component_accumulation(
                                typed.state,
                                porosities[local]);
                    const auto energy =
                        flow::
                            build_pore_volume_energy_accumulation_snapshot(
                                typed.state,
                                porosities[local],
                                typed.transport,
                                typed.caloric,
                                typed.rock);
                    for (std::size_t i = 0U;
                         i < 3U;
                         ++i) {
                        totals[i] +=
                            component
                                .component_accumulation_mol_per_bulk_m3[i] *
                            bulk_volume_m3;
                    }
                    totals[3] +=
                        energy
                            .total_internal_energy_j_per_bulk_m3 *
                        bulk_volume_m3;
                }
            },
            *current[local]);
    }
    return totals;
}

std::array<double, 4>
global_fixed_bhp_well_production_rate(
    const fdp::
        PhaseTransitionRebuiltNaturalVariableSystem3D&
            system,
    Vec state,
    const wdp::
        FixedBhpMultiConnectionWellSourceEvaluatorContext3D&
            well_context,
    double bottom_hole_pressure_pa,
    std::uint64_t*
        authoritative_connection_count) {
    if (state == nullptr ||
        authoritative_connection_count ==
            nullptr ||
        !std::isfinite(
            bottom_hole_pressure_pa) ||
        !(bottom_hole_pressure_pa > 0.0)) {
        throw std::invalid_argument(
            "invalid fixed-BHP well-rate probe input");
    }

    std::vector<std::optional<
        fdp::
            MixedCardinalityPhysicalCurrentCellLinearization3D>>
        current;
    std::vector<double>
        porosities;
    fdp::NaturalVariableSnesEvaluationStatus3D
        status =
            fdp::
                NaturalVariableSnesEvaluationStatus3D::
                    success;
    if (system
            .evaluate_local_cells_for_phase_transition(
                state,
                &current,
                &porosities,
                &status) !=
            PETSC_SUCCESS ||
        status !=
            fdp::
                NaturalVariableSnesEvaluationStatus3D::
                    success ||
        current.size() !=
            porosities.size()) {
        throw std::runtime_error(
            "failed to evaluate fixed-BHP well-rate probe state");
    }

    std::array<double, 4>
        local_rate{};
    std::uint64_t
        local_count = 0U;
    for (std::size_t local = 0U;
         local < current.size();
         ++local) {
        const auto& record =
            system.numbering().cell(
                mesh::LocalIndex{
                    static_cast<
                        mesh::LocalIndex::value_type>(
                            local)});
        if (record.owner_rank !=
            system.numbering().local_rank()) {
            continue;
        }
        const auto* connection =
            well_context.find_connection(
                record.cell_global);
        if (connection == nullptr) {
            continue;
        }
        if (!current[local].has_value()) {
            throw std::runtime_error(
                "authoritative fixed-BHP connection has no current cell state");
        }

        const auto source =
            wdp::
                build_fixed_bhp_peaceman_well_source_3d(
                    connection
                        ->with_bottom_hole_pressure(
                            bottom_hole_pressure_pa),
                    *current[local]);
        for (std::size_t component = 0U;
             component < 3U;
             ++component) {
            local_rate[component] -=
                source
                    .cell_source
                    .component_molar_rate_mol_per_s[
                        component];
        }
        local_rate[3] -=
            source
                .cell_source
                .energy_rate_w;
        ++local_count;
    }

    std::array<double, 4>
        global_rate{};
    std::uint64_t
        global_count = 0U;
    if (MPI_Allreduce(
            local_rate.data(),
            global_rate.data(),
            4,
            MPI_DOUBLE,
            MPI_SUM,
            PETSC_COMM_WORLD) !=
            MPI_SUCCESS ||
        MPI_Allreduce(
            &local_count,
            &global_count,
            1,
            MPI_UINT64_T,
            MPI_SUM,
            PETSC_COMM_WORLD) !=
            MPI_SUCCESS) {
        throw std::runtime_error(
            "failed to reduce fixed-BHP well-rate probe");
    }

    *authoritative_connection_count =
        global_count;
    return global_rate;
}


