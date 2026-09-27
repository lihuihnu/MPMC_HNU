th::Provenance thermodynamic_adapter_source(
    std::string locator) {
    return {
        th::SourceKind::synthetic_test,
        "MPMC_HNU thermodynamic cross-cardinality adapter regression",
        "v1",
        std::move(locator),
        "Manufactured values verify identity-to-selected-branch TPFA integration only",
        "tests/flow_discretization/petsc/mixed_cardinality_physical_snes_assembly_test.cpp",
        "Repository structural regression"};
}

th::SourcedScalar thermodynamic_adapter_scalar(
    double value,
    th::Unit unit,
    std::string locator) {
    return {
        value,
        unit,
        thermodynamic_adapter_source(
            std::move(locator)),
        "SI",
        "identity"};
}

th::Component thermodynamic_adapter_component(
    std::string id,
    double molar_mass_kg_per_mol) {
    const auto definition =
        thermodynamic_adapter_source(
            "component-" + id);
    return {
        id,
        id,
        th::ComponentKind::pure,
        definition,
        thermodynamic_adapter_scalar(
            molar_mass_kg_per_mol,
            th::Unit::kilogram_per_mole,
            "molar-mass-" + id)};
}

th::PrParameterSet
thermodynamic_adapter_pr_parameters() {
    std::vector<th::Component> catalog{
        thermodynamic_adapter_component(
            "A", 0.020),
        thermodynamic_adapter_component(
            "B", 0.030),
        thermodynamic_adapter_component(
            "C", 0.040)};
    const std::vector<std::string>
        order{"A", "B", "C"};

    th::PrParameterInput input;
    input.model_id =
        std::string{th::pr76_profile};
    input.dataset_id =
        "synthetic-cross-cardinality-pr76";
    input.revision = "v1";
    input.applicability = {
        std::nullopt,
        std::nullopt,
        thermodynamic_adapter_source(
            "applicability")};

    const std::array<double, 3>
        critical_temperature{
            20.0, 25.0, 30.0};
    const std::array<double, 3>
        critical_pressure{
            1.0e6, 1.2e6, 1.5e6};
    const std::array<double, 3>
        acentric_factor{
            0.0, 0.1, 0.2};
    for (std::size_t index = 0U;
         index < order.size();
         ++index) {
        input.pure.push_back({
            order[index],
            thermodynamic_adapter_scalar(
                critical_temperature[index],
                th::Unit::kelvin,
                order[index] + "-Tc"),
            thermodynamic_adapter_scalar(
                critical_pressure[index],
                th::Unit::pascal,
                order[index] + "-Pc"),
            thermodynamic_adapter_scalar(
                acentric_factor[index],
                th::Unit::dimensionless,
                order[index] + "-omega")});
    }
    for (std::size_t first = 0U;
         first < order.size();
         ++first) {
        for (std::size_t second =
                 first + 1U;
             second < order.size();
             ++second) {
            input.binary.push_back({
                order[first],
                order[second],
                thermodynamic_adapter_scalar(
                    0.0,
                    th::Unit::dimensionless,
                    order[first] + "-" +
                        order[second] +
                        "-kij")});
        }
    }
    return th::PrParameterSet::create(
        catalog,
        order,
        input,
        th::DataPolicy::
            allow_synthetic_tests);
}

flow::FrozenPhysicalPhaseIdentity
mixed_physical_phase_identity(
    std::string key) {
    return {
        "fixture/mixed-physical-dispatch",
        std::move(key)};
}

flow::PhaseIdentityContinuationSnapshot
thermodynamic_adapter_continuation() {
    flow::PhaseSetTransitionCandidate
        candidate;
    candidate.source_phase_count = 1U;
    candidate.target_phase_count = 2U;
    candidate.trigger =
        flow::
            PhaseSetTransitionTrigger::
                stability_witness;
    candidate.status =
        flow::
            PhaseSetTransitionCandidateStatus::
                target_resolved;
    candidate.evidence_profile =
        "fixture/cross-cardinality-fresh-target/v1";
    candidate.diagnostic =
        "explicit 1P->2P phase identity continuation";

    return flow::
        make_phase_identity_continuation_snapshot(
            candidate,
            flow::FrozenActivePhaseIdentityMap{
                {
                    mixed_physical_phase_identity(
                        "aqueous")}},
            flow::FrozenActivePhaseIdentityMap{
                {
                    mixed_physical_phase_identity(
                        "aqueous"),
                    mixed_physical_phase_identity(
                        "hydrocarbon-0")}});
}

flow::Pr76AbsentPhasePotentialExtensionProvider<double>
make_thermodynamic_adapter_pr_provider(
    const th::Pr76Phase<double>& model) {
    const auto continuation =
        thermodynamic_adapter_continuation();
    auto registry =
        flow::
            make_transition_selected_phase_branch_registry(
                continuation,
                std::vector<
                    flow::
                        FrozenSelectedPhaseBranchBinding<
                            th::Pr76SelectedPhase>>{
                    {
                        mixed_physical_phase_identity(
                            "aqueous"),
                        {0U, {}},
                        "fixture/aqueous-root0"},
                    {
                        mixed_physical_phase_identity(
                            "hydrocarbon-0"),
                        {0U, {}},
                        "fixture/hydrocarbon0-root0"}});
    return {
        model,
        std::move(registry)};
}

struct ThermodynamicCoordinateResolverAudit {
    std::uint64_t calls{};
};

PetscErrorCode
resolve_absent_phase_thermodynamic_coordinates(
    const fdp::
        MixedCardinalityPhysicalSnesAuthoritativeFaceInput3D&,
    const flow::
        CrossCardinalityFacePhaseIdentityBinding&
            phase_binding,
    fdp::CrossCardinalityAbsentPhaseSide3D
        absent_side,
    const fdp::
        MixedCardinalityPhysicalCurrentCellLinearization3D&
            owner,
    const fdp::
        MixedCardinalityPhysicalCurrentCellLinearization3D&
            neighbour,
    void* raw_context,
    std::optional<
        flow::
            AbsentPhaseThermodynamicCoordinateExtension>*
                output,
    fdp::NaturalVariableSnesEvaluationStatus3D*
        status) {
    if (raw_context == nullptr ||
        output == nullptr ||
        status == nullptr) {
        return PETSC_ERR_ARG_NULL;
    }
    output->reset();
    auto* audit =
        static_cast<
            ThermodynamicCoordinateResolverAudit*>(
                raw_context);
    ++audit->calls;

    try {
        if (phase_binding
                .identity.opaque_phase_key !=
                "hydrocarbon-0" ||
            absent_side !=
                fdp::
                    CrossCardinalityAbsentPhaseSide3D::
                        owner) {
            return PETSC_ERR_ARG_INCOMP;
        }

        const auto host_identity =
            std::visit(
                [](const auto& typed) {
                    return typed.transport
                        .state_identity;
                },
                owner);
        (void)neighbour;

        const std::size_t q =
            host_identity.layout
                .unknown_count();
        const std::size_t n =
            host_identity.layout
                .component_count();
        if (n != 3U) {
            return PETSC_ERR_ARG_INCOMP;
        }

        flow::
            AbsentPhaseThermodynamicCoordinateExtension
            coordinates;
        coordinates.identity =
            phase_binding.identity;
        coordinates.host_state_identity =
            host_identity;
        coordinates.phase_pressure_pa =
            host_identity
                .reference_pressure_pa;
        coordinates.phase_pressure_gradient
            .assign(q, 0.0);
        coordinates.phase_pressure_gradient[
            host_identity.layout
                .pressure_unknown_index()] =
            1.0;
        coordinates.hypothetical_composition =
            {0.40, 0.20, 0.40};
        coordinates
            .hypothetical_composition_jacobian
            .assign(n * q, 0.0);
        coordinates.provenance =
            "fixture/stability-hypothetical-composition/v1";
        coordinates.validate();

        output->emplace(
            std::move(coordinates));
        *status =
            fdp::
                NaturalVariableSnesEvaluationStatus3D::
                    success;
        return PETSC_SUCCESS;
    } catch (const std::invalid_argument&) {
        return PETSC_ERR_ARG_INCOMP;
    } catch (...) {
        return PETSC_ERR_LIB;
    }
}

PetscErrorCode
evaluate_absent_phase_extension(
    const fdp::
        MixedCardinalityPhysicalSnesAuthoritativeFaceInput3D&,
    const flow::
        CrossCardinalityFacePhaseIdentityBinding&
            phase_binding,
    fdp::CrossCardinalityAbsentPhaseSide3D
        absent_side,
    const fdp::
        MixedCardinalityPhysicalCurrentCellLinearization3D&
            owner,
    const fdp::
        MixedCardinalityPhysicalCurrentCellLinearization3D&
            neighbour,
    void*,
    std::optional<
        flow::
            AbsentPhasePotentialExtensionLinearization>*
                output,
    fdp::NaturalVariableSnesEvaluationStatus3D*
        status) {
    if (output == nullptr ||
        status == nullptr) {
        return PETSC_ERR_ARG_NULL;
    }
    output->reset();

    try {
        const auto& absent_cell =
            absent_side ==
                    fdp::
                        CrossCardinalityAbsentPhaseSide3D::
                            owner
                ? owner
                : neighbour;
        const auto identity =
            std::visit(
                [](const auto& typed) {
                    return typed.transport.state_identity;
                },
                absent_cell);
        const std::size_t q =
            identity.layout.unknown_count();

        std::vector<double>
            pressure_gradient(
                q,
                0.0);
        pressure_gradient[
            identity.layout
                .pressure_unknown_index()] =
            1.0;
        std::vector<double>
            density_gradient(
                q,
                0.0);

        double density =
            1.0;
        if (phase_binding
                .identity.opaque_phase_key ==
            "hydrocarbon-0") {
            density = 1.5;
        } else if (
            phase_binding
                .identity.opaque_phase_key ==
            "hydrocarbon-1") {
            density = 1.8;
        }

        output->emplace(
            phase_binding.identity,
            identity.reference_pressure_pa,
            density,
            q,
            std::move(
                pressure_gradient),
            std::move(
                density_gradient),
            "fixture/cross-cardinality-potential-extension/v1");
        *status =
            fdp::
                NaturalVariableSnesEvaluationStatus3D::
                    success;
        return PETSC_SUCCESS;
    } catch (const std::invalid_argument&) {
        *status =
            fdp::
                NaturalVariableSnesEvaluationStatus3D::
                    domain_error;
        return PETSC_SUCCESS;
    } catch (...) {
        return PETSC_ERR_ARG_INCOMP;
    }
}

void check_direct_cross_cardinality_bridge(
    std::uint64_t owner_stable,
    std::uint64_t neighbour_stable,
    std::size_t face_index,
    DispatchAudit* audit) {
    const auto owner =
        evaluate_target(
            owner_stable,
            audit);
    const auto neighbour =
        evaluate_target(
            neighbour_stable,
            audit);
    const auto plan =
        flow::
            make_cross_cardinality_face_phase_identity_plan(
                phase_identity_map_for_stable(
                    owner_stable),
                phase_identity_map_for_stable(
                    neighbour_stable));

    fdp::CrossCardinalityTpfaBridgeBinding3D
        bridge_binding{
            &evaluate_absent_phase_extension,
            nullptr};
    std::optional<
        fdp::
            MixedCardinalityPhysicalFaceLinearization3D>
        face;
    fdp::NaturalVariableSnesEvaluationStatus3D
        status =
            fdp::
                NaturalVariableSnesEvaluationStatus3D::
                    success;

    const auto face_input =
        make_face_input(
            face_index,
            UINT64_C(1000) +
                static_cast<std::uint64_t>(
                    face_index));
    const PetscErrorCode error =
        fdp::
            evaluate_standard_cross_cardinality_tpfa_face_3d(
                face_input,
                plan,
                {2.0, 3.0},
                owner,
                neighbour,
                &bridge_binding,
                &face,
                &status);

    require_collective(
        error == PETSC_SUCCESS &&
            status ==
                fdp::
                    NaturalVariableSnesEvaluationStatus3D::
                        success &&
            face.has_value(),
        "standard cross-cardinality TPFA bridge evaluation failed");

    const auto& component =
        face->component;
    const auto& energy =
        face->energy;
    require_collective(
        component.owner_state_identity
                    .layout.phase_count() ==
                phase_count_for(
                    owner_stable) &&
            component.neighbour_state_identity
                    .layout.phase_count() ==
                phase_count_for(
                    neighbour_stable) &&
            component
                    .owner_row_neighbour_column_jacobian
                    .size() ==
                component.component_count() *
                    component
                        .neighbour_state_identity
                        .layout.unknown_count() &&
            component
                    .neighbour_row_owner_column_jacobian
                    .size() ==
                component.component_count() *
                    component
                        .owner_state_identity
                        .layout.unknown_count(),
        "cross-cardinality TPFA rectangular Jacobian shape mismatch");

    double component_norm = 0.0;
    for (const double value :
         component
             .owner_component_contribution_mol_per_bulk_m3_s) {
        component_norm +=
            std::abs(value);
    }
    require_collective(
        component_norm > 0.0 &&
            std::abs(
                energy
                    .owner_contribution_w_per_bulk_m3) >
                0.0,
        "cross-cardinality bridge did not produce physical component/energy flux");

    for (std::size_t component_index = 0U;
         component_index <
            component.component_count();
         ++component_index) {
        const double owner_rate =
            component
                .owner_component_contribution_mol_per_bulk_m3_s[
                    component_index] *
            component.bulk_volume
                .owner_bulk_volume_m3;
        const double neighbour_rate =
            component
                .neighbour_component_contribution_mol_per_bulk_m3_s[
                    component_index] *
            component.bulk_volume
                .neighbour_bulk_volume_m3;
        require_collective(
            std::abs(
                owner_rate +
                neighbour_rate) <
                1.0e-12 *
                    std::max(
                        1.0,
                        std::abs(
                            owner_rate)),
            "cross-cardinality component face rate is not conservative");
    }

    const double owner_energy =
        energy
            .owner_contribution_w_per_bulk_m3 *
        energy.bulk_volume
            .owner_bulk_volume_m3;
    const double neighbour_energy =
        energy
            .neighbour_contribution_w_per_bulk_m3 *
        energy.bulk_volume
            .neighbour_bulk_volume_m3;
    require_collective(
        std::abs(
            owner_energy +
            neighbour_energy) <
            1.0e-12 *
                std::max(
                    1.0,
                    std::abs(
                        owner_energy)),
        "cross-cardinality energy face rate is not conservative");
}

void insert_target(
    Vec state,
    const fdp::
        VariableCardinalityNaturalVariableNumbering3D&
            numbering) {
    bool local_ok = true;
    for (const auto& cell :
         numbering.cells()) {
        if (cell.owner_rank !=
            numbering.local_rank()) {
            continue;
        }
        const auto target =
            target_state(
                cell.cell_global.value());
        local_ok =
            local_ok &&
            target.size() ==
                cell.scalar_count;
        for (std::size_t slot = 0U;
             slot < cell.scalar_count;
             ++slot) {
            const PetscInt index =
                cell.petsc_global_scalar_start +
                static_cast<PetscInt>(
                    slot);
            const PetscScalar value =
                static_cast<PetscScalar>(
                    target[slot]);
            local_ok =
                local_ok &&
                VecSetValues(
                    state,
                    1,
                    &index,
                    &value,
                    INSERT_VALUES) ==
                    PETSC_SUCCESS;
        }
    }
    require_collective(
        local_ok,
        "failed to insert mixed physical target state");
    require_collective(
        VecAssemblyBegin(state) ==
                PETSC_SUCCESS &&
            VecAssemblyEnd(state) ==
                PETSC_SUCCESS,
        "failed to assemble mixed physical target state");
}


fdp::FrozenAbsentPhaseCoordinateEntry3D
make_frozen_absent_coordinate_entry(
    std::uint64_t stable,
    flow::FrozenPhysicalPhaseIdentity identity,
    std::vector<double> composition,
    std::string branch_provenance,
    DispatchAudit* audit) {
    const auto current =
        evaluate_target(
            stable,
            audit);
    const auto host_identity =
        std::visit(
            [](const auto& typed) {
                return typed.transport
                    .state_identity;
            },
            current);
    const std::size_t q =
        host_identity.layout
            .unknown_count();
    const std::size_t n =
        host_identity.layout
            .component_count();
    if (composition.size() != n) {
        throw std::invalid_argument(
            "frozen absent coordinate composition size mismatch");
    }

    flow::
        AbsentPhaseThermodynamicCoordinateExtension
        coordinates;
    coordinates.identity =
        identity;
    coordinates.host_state_identity =
        host_identity;
    coordinates.phase_pressure_pa =
        host_identity
            .reference_pressure_pa;
    coordinates.phase_pressure_gradient
        .assign(q, 0.0);
    coordinates.phase_pressure_gradient[
        host_identity.layout
            .pressure_unknown_index()] =
        1.0;
    coordinates.hypothetical_composition =
        std::move(composition);
    coordinates
        .hypothetical_composition_jacobian
        .assign(n * q, 0.0);
    coordinates.provenance =
        "fixture/frozen-absent-coordinate-reference/v1";
    coordinates.validate();

    return {
        std::move(identity),
        std::move(coordinates),
        std::move(branch_provenance)};
}

flow::Pr76AbsentPhasePotentialExtensionProvider<double>
make_outer_rebuild_pr_provider(
    const th::Pr76Phase<double>& model) {
    flow::PhaseSetTransitionCandidate
        candidate;
    candidate.source_phase_count = 1U;
    candidate.target_phase_count = 3U;
    candidate.status =
        flow::
            PhaseSetTransitionCandidateStatus::
                target_resolved;
    candidate.trigger =
        flow::
            PhaseSetTransitionTrigger::
                stability_witness;
    candidate.evidence_profile =
        "fixture/outer-rebuild-phase-universe/v1";

    const auto continuation =
        flow::
            make_phase_identity_continuation_snapshot(
                candidate,
                flow::FrozenActivePhaseIdentityMap{
                    {
                        mixed_physical_phase_identity(
                            "aqueous")}},
                flow::FrozenActivePhaseIdentityMap{
                    {
                        mixed_physical_phase_identity(
                            "aqueous"),
                        mixed_physical_phase_identity(
                            "hydrocarbon-0"),
                        mixed_physical_phase_identity(
                            "hydrocarbon-1")}});

    auto registry =
        flow::
            make_transition_selected_phase_branch_registry(
                continuation,
                std::vector<
                    flow::
                        FrozenSelectedPhaseBranchBinding<
                            th::Pr76SelectedPhase>>{
                    {
                        mixed_physical_phase_identity(
                            "aqueous"),
                        {0U, {}},
                        "fixture/aqueous-root0"},
                    {
                        mixed_physical_phase_identity(
                            "hydrocarbon-0"),
                        {0U, {}},
                        "fixture/hydrocarbon0-root0"},
                    {
                        mixed_physical_phase_identity(
                            "hydrocarbon-1"),
                        {0U, {}},
                        "fixture/hydrocarbon1-root0"}});

    return {
        model,
        std::move(registry)};
}

std::vector<
    fdp::FrozenPhaseTransitionRebuildCell3D>
make_outer_rebuild_cells(
    int rank,
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
                using Typed =
                    std::remove_cvref_t<
                        decltype(typed)>;
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
            "fixture/frozen-active-set/no-transition/v1";

        if (stable == UINT64_C(20)) {
            snapshot.frozen_absent_phases.push_back(
                make_frozen_absent_coordinate_entry(
                    stable,
                    mixed_physical_phase_identity(
                        "hydrocarbon-0"),
                    {0.40, 0.20, 0.40},
                    "fixture/hydrocarbon0-root0",
                    audit));
            snapshot.frozen_absent_phases.push_back(
                make_frozen_absent_coordinate_entry(
                    stable,
                    mixed_physical_phase_identity(
                        "hydrocarbon-1"),
                    {0.10, 0.50, 0.40},
                    "fixture/hydrocarbon1-root0",
                    audit));
        } else if (
            stable == UINT64_C(40)) {
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

    // Accept a real 2P -> 3P transition on stable cell30.  The target phase
    // compositions are identical to the source overall inventory so the
    // transition projection passes the existing material-balance gate.
    constexpr std::size_t transitioned_local = 2U;
    constexpr double porosity = 0.26;
    const auto source_current =
        evaluate_target(
            UINT64_C(30),
            audit);
    const auto& source_two_phase =
        std::get<
            fdp::
                TwoPhaseCurrentCellLinearization3D>(
                    source_current);
    const auto source_component =
        flow::
            build_two_phase_component_accumulation(
                source_two_phase.state,
                porosity);
    const auto source_energy =
        flow::
            build_two_phase_energy_accumulation_snapshot(
                source_two_phase.state,
                porosity,
                source_two_phase.transport,
                source_two_phase.caloric,
                source_two_phase.rock);

    double inventory_total = 0.0;
    for (const double value :
         source_component
             .component_accumulation_mol_per_bulk_m3) {
        inventory_total += value;
    }
    std::vector<double>
        overall;
    overall.reserve(3U);
    for (const double value :
         source_component
             .component_accumulation_mol_per_bulk_m3) {
        overall.push_back(
            value / inventory_total);
    }

    flow::PhaseSetTransitionCandidate
        candidate;
    candidate.source_phase_count = 2U;
    candidate.target_phase_count = 3U;
    candidate.trigger =
        flow::
            PhaseSetTransitionTrigger::
                stability_witness;
    candidate.status =
        flow::
            PhaseSetTransitionCandidateStatus::
                target_resolved;
    candidate.pressure_pa = 20.0;
    candidate.temperature_k = 9.0;
    candidate.component_ids =
        {"A", "B", "C"};
    candidate.target_phases = {
        {0.20, overall, 2.0},
        {0.30, overall, 4.0},
        {0.50, overall, 7.0}};
    candidate.evidence_profile =
        "fixture/accepted-2p-to-3p/v1";
    candidate.diagnostic =
        "controlled outer rebuild transition";

    result[transitioned_local] =
        fdp::
            make_accepted_phase_transition_rebuild_cell_3d(
                mesh::LocalIndex{
                    static_cast<
                        mesh::LocalIndex::value_type>(
                            transitioned_local)},
                mesh::GlobalEntityId{
                    UINT64_C(30)},
                4.0,
                porosity,
                candidate,
                source_component,
                source_component,
                source_energy,
                phase_maps[transitioned_local],
                flow::FrozenActivePhaseIdentityMap{
                    {
                        mixed_physical_phase_identity(
                            "aqueous"),
                        mixed_physical_phase_identity(
                            "hydrocarbon-0"),
                        mixed_physical_phase_identity(
                            "hydrocarbon-1")}},
                {});

    return result;
}


