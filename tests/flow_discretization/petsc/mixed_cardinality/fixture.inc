namespace mesh = mpmc::mesh;
namespace flow = mpmc::flow;
namespace disc = mpmc::discretization;
namespace dp = mpmc::discretization_petsc;
namespace fd = mpmc::flow_discretization;
namespace fdp = mpmc::flow_discretization_petsc;
namespace wd = mpmc::well_discretization;
namespace wdp = mpmc::well_discretization_petsc;
namespace well = mpmc::well;
namespace th = mpmc::thermodynamics;

void require_collective(
    bool condition,
    std::string_view message) {
    int local = condition ? 1 : 0;
    int global = 0;
    if (MPI_Allreduce(
            &local,
            &global,
            1,
            MPI_INT,
            MPI_MIN,
            PETSC_COMM_WORLD) !=
        MPI_SUCCESS) {
        throw std::runtime_error(
            "MPI_Allreduce failed in mixed physical dispatcher regression");
    }
    if (global == 0) {
        throw std::runtime_error(
            std::string{message});
    }
}

void near_collective(
    double actual,
    double expected,
    double relative,
    double absolute) {
    require_collective(
        std::isfinite(actual) &&
            std::isfinite(expected) &&
            std::abs(actual - expected) <=
                absolute +
                    relative *
                        std::max(
                            std::abs(actual),
                            std::abs(expected)),
        "numeric mismatch in mixed physical dispatcher regression");
}

[[nodiscard]] std::size_t
phase_count_for(
    std::uint64_t stable) {
    if (stable == UINT64_C(10) ||
        stable == UINT64_C(20)) {
        return 1U;
    }
    if (stable == UINT64_C(30) ||
        stable == UINT64_C(40)) {
        return 2U;
    }
    if (stable == UINT64_C(50) ||
        stable == UINT64_C(60)) {
        return 3U;
    }
    throw std::invalid_argument(
        "unknown mixed physical stable cell");
}

[[nodiscard]] flow::NaturalVariableLayout1P
layout_1p() {
    return flow::NaturalVariableLayout1P{
        flow::NaturalVariableCompositionPivot1P::
            fixed_last(3U)};
}

[[nodiscard]] flow::NaturalVariableLayout2P
layout_2p() {
    return flow::NaturalVariableLayout2P{
        flow::NaturalVariableCompositionPivot2P::
            fixed_last(3U)};
}

[[nodiscard]] flow::NaturalVariableLayout3P
layout_3p() {
    return flow::NaturalVariableLayout3P{
        flow::NaturalVariableCompositionPivot3P::
            fixed_last(3U)};
}

[[nodiscard]] std::vector<double>
target_1p() {
    auto layout = layout_1p();
    std::vector<double> q(
        layout.unknown_count(),
        0.0);
    q[layout.pressure_unknown_index()] =
        10.0;
    q[layout.temperature_unknown_index()] =
        8.0;
    q[*layout.independent_composition_unknown_index(
        0U)] = 0.20;
    q[*layout.independent_composition_unknown_index(
        1U)] = 0.30;
    return q;
}

[[nodiscard]] std::vector<double>
target_2p() {
    auto layout = layout_2p();
    std::vector<double> q(
        layout.unknown_count(),
        0.0);
    q[layout.pressure_unknown_index()] =
        20.0;
    q[layout.temperature_unknown_index()] =
        9.0;
    q[layout.independent_saturation_unknown_index()] =
        0.55;
    const std::array<
        std::array<double, 3>,
        2>
        x{{
            {{0.20, 0.30, 0.50}},
            {{0.40, 0.20, 0.40}}
        }};
    for (std::size_t phase = 0U;
         phase < 2U;
         ++phase) {
        for (std::size_t component = 0U;
             component < 3U;
             ++component) {
            const auto column =
                layout
                    .independent_composition_unknown_index(
                        phase,
                        component);
            if (column) {
                q[*column] =
                    x[phase][component];
            }
        }
    }
    return q;
}

[[nodiscard]] std::vector<double>
target_3p() {
    auto layout = layout_3p();
    std::vector<double> q(
        layout.unknown_count(),
        0.0);
    q[layout.pressure_unknown_index()] =
        30.0;
    q[layout.temperature_unknown_index()] =
        10.0;
    const auto s0 =
        layout.independent_saturation_unknown_index(
            flow::PhaseSlot3::phase0);
    const auto s1 =
        layout.independent_saturation_unknown_index(
            flow::PhaseSlot3::phase1);
    if (!s0 || !s1) {
        throw std::logic_error(
            "three-phase target layout lacks saturation coordinates");
    }
    q[*s0] = 0.20;
    q[*s1] = 0.30;

    const std::array<
        std::array<double, 3>,
        3>
        x{{
            {{0.20, 0.30, 0.50}},
            {{0.40, 0.20, 0.40}},
            {{0.10, 0.50, 0.40}}
        }};
    for (std::size_t phase = 0U;
         phase < 3U;
         ++phase) {
        const auto slot =
            static_cast<flow::PhaseSlot3>(
                phase);
        for (std::size_t component = 0U;
             component < 3U;
             ++component) {
            const auto column =
                layout
                    .independent_composition_unknown_index(
                        slot,
                        component);
            if (column) {
                q[*column] =
                    x[phase][component];
            }
        }
    }
    return q;
}

[[nodiscard]] std::vector<double>
target_state(
    std::uint64_t stable) {
    switch (phase_count_for(stable)) {
    case 1U:
        return target_1p();
    case 2U:
        return target_2p();
    case 3U:
        return target_3p();
    default:
        throw std::logic_error(
            "invalid target phase count");
    }
}

struct DispatchAudit {
    std::uint64_t single_calls{};
    std::uint64_t two_calls{};
    std::uint64_t three_calls{};
    bool source_enabled{};
    bool well_timestep_compressibility{};
    std::uint64_t source_evaluator_calls{};
    std::uint64_t source_cell20_calls{};
};

PetscErrorCode evaluate_explicit_cell_source(
    mesh::LocalIndex,
    mesh::GlobalEntityId cell_global,
    std::span<const double> natural_variables,
    const fdp::MixedCardinalityPhysicalCurrentCellLinearization3D&
        current,
    void* raw_context,
    std::optional<fd::CellSourceLinearization3D>* output,
    fdp::NaturalVariableSnesEvaluationStatus3D* status) {
    if (raw_context == nullptr ||
        output == nullptr ||
        status == nullptr) {
        return PETSC_ERR_ARG_NULL;
    }
    output->reset();
    auto* audit =
        static_cast<DispatchAudit*>(raw_context);
    ++audit->source_evaluator_calls;
    *status =
        fdp::NaturalVariableSnesEvaluationStatus3D::
            success;
    if (!audit->source_enabled ||
        cell_global !=
            mesh::GlobalEntityId{UINT64_C(20)}) {
        return PETSC_SUCCESS;
    }
    ++audit->source_cell20_calls;

    const auto identity =
        std::visit(
            [](const auto& typed) {
                return typed.transport.state_identity;
            },
            current);
    const std::size_t q =
        natural_variables.size();
    if (identity.component_ids.size() != 3U ||
        q != identity.layout.unknown_count()) {
        return PETSC_ERR_ARG_INCOMP;
    }

    fd::CellSourceLinearization3D source;
    source.provenance =
        "test/mixed-cell-source/v1";
    source.component_ids =
        identity.component_ids;
    source.input_count = q;
    source.component_molar_rate_mol_per_s =
        {2.0, -1.0, 0.5};
    source.component_molar_rate_jacobian_mol_per_s.assign(
        3U * q,
        0.0);
    for (std::size_t component = 0U;
         component < 3U;
         ++component) {
        for (std::size_t column = 0U;
             column < q;
             ++column) {
            source.component_molar_rate_jacobian_mol_per_s[
                component * q + column] =
                0.1 *
                static_cast<double>(
                    (component + 1U) *
                    (column + 1U));
        }
    }
    source.energy_rate_w = 20.0;
    source.energy_rate_gradient_w.resize(q);
    for (std::size_t column = 0U;
         column < q;
         ++column) {
        source.energy_rate_gradient_w[column] =
            0.2 *
            static_cast<double>(
                column + 1U);
    }
    output->emplace(std::move(source));
    return PETSC_SUCCESS;
}

struct FixedBhpWellSourceAudit {
    wdp::
        FixedBhpPeacemanWellSourceEvaluatorContext3D*
            context{};
    std::uint64_t evaluator_calls{};
    std::uint64_t target_calls{};
};

PetscErrorCode
evaluate_audited_fixed_bhp_well_source(
    mesh::LocalIndex cell,
    mesh::GlobalEntityId cell_global,
    std::span<const double> natural_variables,
    const fdp::
        MixedCardinalityPhysicalCurrentCellLinearization3D&
            current,
    void* raw_context,
    std::optional<fd::CellSourceLinearization3D>* output,
    fdp::NaturalVariableSnesEvaluationStatus3D* status) {
    if (raw_context == nullptr) {
        return PETSC_ERR_ARG_NULL;
    }
    auto* audit =
        static_cast<FixedBhpWellSourceAudit*>(
            raw_context);
    if (audit->context == nullptr) {
        return PETSC_ERR_ARG_NULL;
    }
    ++audit->evaluator_calls;
    if (cell_global ==
        audit->context->target_cell_global()) {
        ++audit->target_calls;
    }
    return wdp::
        evaluate_fixed_bhp_peaceman_well_source_3d(
            cell,
            cell_global,
            natural_variables,
            current,
            audit->context,
            output,
            status);
}

struct FixedBhpMultiConnectionWellSourceAudit {
    wdp::
        FixedBhpMultiConnectionWellSourceEvaluatorContext3D*
            context{};
    std::uint64_t evaluator_calls{};
    std::uint64_t cell30_calls{};
    std::uint64_t cell60_calls{};
};

PetscErrorCode
evaluate_audited_fixed_bhp_multi_connection_well_source(
    mesh::LocalIndex cell,
    mesh::GlobalEntityId cell_global,
    std::span<const double> natural_variables,
    const fdp::
        MixedCardinalityPhysicalCurrentCellLinearization3D&
            current,
    void* raw_context,
    std::optional<fd::CellSourceLinearization3D>* output,
    fdp::NaturalVariableSnesEvaluationStatus3D* status) {
    if (raw_context == nullptr) {
        return PETSC_ERR_ARG_NULL;
    }
    auto* audit =
        static_cast<
            FixedBhpMultiConnectionWellSourceAudit*>(
                raw_context);
    if (audit->context == nullptr) {
        return PETSC_ERR_ARG_NULL;
    }

    ++audit->evaluator_calls;
    if (audit->context
            ->find_connection(cell_global) !=
        nullptr) {
        if (cell_global ==
            mesh::GlobalEntityId{
                UINT64_C(30)}) {
            ++audit->cell30_calls;
        } else if (
            cell_global ==
            mesh::GlobalEntityId{
                UINT64_C(60)}) {
            ++audit->cell60_calls;
        }
    }

    return wdp::
        evaluate_fixed_bhp_multi_connection_well_source_3d(
            cell,
            cell_global,
            natural_variables,
            current,
            audit->context,
            output,
            status);
}

PetscErrorCode
evaluate_1p(
    mesh::LocalIndex,
    mesh::GlobalEntityId,
    std::span<const double> natural_variables,
    const flow::NaturalVariableLayout1P&
        frozen_layout,
    std::span<const std::string> component_ids,
    void* raw_context,
    std::optional<
        fdp::SinglePhaseCurrentCellLinearization3D>*
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
        static_cast<DispatchAudit*>(
            raw_context);
    ++audit->single_calls;

    try {
        if (natural_variables.size() !=
                frozen_layout.unknown_count() ||
            component_ids.size() !=
                frozen_layout.component_count()) {
            return PETSC_ERR_ARG_INCOMP;
        }
        const double p =
            natural_variables[
                frozen_layout
                    .pressure_unknown_index()];
        const double temperature =
            natural_variables[
                frozen_layout
                    .temperature_unknown_index()];
        if (!std::isfinite(p) ||
            !(p > 0.0) ||
            !std::isfinite(temperature) ||
            !(temperature > 0.0)) {
            *status =
                fdp::NaturalVariableSnesEvaluationStatus3D::
                    domain_error;
            return PETSC_SUCCESS;
        }

        flow::NaturalVariableCellStateInput1P
            input;
        input.component_ids.assign(
            component_ids.begin(),
            component_ids.end());
        input.reference_pressure_pa = p;
        input.temperature_k = temperature;
        input.composition_pivot =
            frozen_layout.composition_pivot();
        for (std::size_t rank = 0U;
             rank < 2U;
             ++rank) {
            const auto component =
                frozen_layout
                    .independent_composition_component(
                        rank);
            const auto column =
                frozen_layout
                    .independent_composition_unknown_index(
                        component);
            if (!column) {
                return PETSC_ERR_PLIB;
            }
            input.independent_composition
                .push_back(
                    natural_variables[*column]);
        }
        const bool synthetic_compressibility =
            audit
                ->well_timestep_compressibility;
        const double pressure_delta =
            p - 10.0;
        const double molar_density =
            5.0 +
            (synthetic_compressibility
                 ? 0.05 * pressure_delta
                 : 0.0);
        const double mass_density =
            2.0 +
            (synthetic_compressibility
                 ? 0.02 * pressure_delta
                 : 0.0);
        if (!(molar_density > 0.0) ||
            !(mass_density > 0.0)) {
            *status =
                fdp::NaturalVariableSnesEvaluationStatus3D::
                    domain_error;
            return PETSC_SUCCESS;
        }
        input.phase_properties =
            flow::PhasePropertyPrerequisiteInput{
                molar_density,
                mass_density,
                1.0,
                temperature,
                0.5 * temperature};

        auto state =
            flow::NaturalVariableCellState1P::
                create(std::move(input));
        const std::size_t q =
            frozen_layout.unknown_count();
        std::vector<double> zero(
            q,
            0.0);
        std::vector<double> molar_density_gradient(
            q,
            0.0);
        std::vector<double> mass_density_gradient(
            q,
            0.0);
        if (synthetic_compressibility) {
            molar_density_gradient[
                frozen_layout
                    .pressure_unknown_index()] =
                0.05;
            mass_density_gradient[
                frozen_layout
                    .pressure_unknown_index()] =
                0.02;
        }
        auto molar =
            flow::
                make_single_phase_molar_density_linearization(
                    state,
                    molar_density_gradient);
        const flow::TransportPropertyProvenance
            provenance{
                "mixed-physical-dispatch",
                "controlled-regression",
                "v1"};
        auto transport =
            flow::
                make_single_phase_transport_linearization(
                    state,
                    mass_density_gradient,
                    zero,
                    1.0,
                    zero,
                    provenance,
                    provenance);
        std::vector<double> dh(
            q,
            0.0);
        std::vector<double> du(
            q,
            0.0);
        dh[frozen_layout
               .temperature_unknown_index()] =
            1.0;
        du[frozen_layout
               .temperature_unknown_index()] =
            0.5;
        auto caloric =
            flow::
                make_single_phase_caloric_linearization(
                    state,
                    dh,
                    du,
                    provenance,
                    provenance);
        std::vector<double> rock_gradient(
            q,
            0.0);
        rock_gradient[
            frozen_layout
                .temperature_unknown_index()] =
            2.0;
        auto rock =
            flow::
                make_single_phase_rock_thermal_storage_linearization(
                    state,
                    2.0 * temperature,
                    rock_gradient,
                    provenance);

        output->emplace(
            fdp::
                SinglePhaseCurrentCellLinearization3D{
                    std::move(state),
                    std::move(molar),
                    std::move(transport),
                    std::move(caloric),
                    std::move(rock)});
        *status =
            fdp::NaturalVariableSnesEvaluationStatus3D::
                success;
        return PETSC_SUCCESS;
    } catch (const std::invalid_argument&) {
        *status =
            fdp::NaturalVariableSnesEvaluationStatus3D::
                domain_error;
        return PETSC_SUCCESS;
    } catch (...) {
        return PETSC_ERR_ARG_INCOMP;
    }
}

PetscErrorCode
evaluate_2p(
    mesh::LocalIndex,
    mesh::GlobalEntityId,
    std::span<const double> natural_variables,
    const flow::NaturalVariableLayout2P&
        frozen_layout,
    std::span<const std::string> component_ids,
    void* raw_context,
    std::optional<
        fdp::TwoPhaseCurrentCellLinearization3D>*
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
        static_cast<DispatchAudit*>(
            raw_context);
    ++audit->two_calls;

    try {
        if (natural_variables.size() !=
                frozen_layout.unknown_count() ||
            component_ids.size() !=
                frozen_layout.component_count()) {
            return PETSC_ERR_ARG_INCOMP;
        }

        const double p =
            natural_variables[
                frozen_layout
                    .pressure_unknown_index()];
        const double temperature =
            natural_variables[
                frozen_layout
                    .temperature_unknown_index()];
        const double s0 =
            natural_variables[
                frozen_layout
                    .independent_saturation_unknown_index()];
        if (!std::isfinite(p) ||
            !(p > 0.0) ||
            !std::isfinite(temperature) ||
            !(temperature > 0.0) ||
            !std::isfinite(s0) ||
            !(s0 > 0.0) ||
            !(s0 < 1.0)) {
            *status =
                fdp::NaturalVariableSnesEvaluationStatus3D::
                    domain_error;
            return PETSC_SUCCESS;
        }

        flow::NaturalVariableCellStateInput2P
            input;
        input.component_ids.assign(
            component_ids.begin(),
            component_ids.end());
        input.reference_pressure_pa = p;
        input.temperature_k = temperature;
        input.independent_saturation = s0;
        input.composition_pivot =
            frozen_layout.composition_pivot();

        for (std::size_t phase = 0U;
             phase < 2U;
             ++phase) {
            for (std::size_t rank = 0U;
                 rank < 2U;
                 ++rank) {
                const auto component =
                    frozen_layout
                        .independent_composition_component(
                            phase,
                            rank);
                const auto column =
                    frozen_layout
                        .independent_composition_unknown_index(
                            phase,
                            component);
                if (!column) {
                    return PETSC_ERR_PLIB;
                }
                input.independent_phase_compositions[
                    phase]
                    .push_back(
                        natural_variables[
                            *column]);
            }
        }

        const bool synthetic_compressibility =
            audit->well_timestep_compressibility;
        const double pressure_delta =
            p - 20.0;
        const std::array<double, 2>
            molar_density{
                6.0 +
                    (synthetic_compressibility
                         ? 0.06 * pressure_delta
                         : 0.0),
                3.0 +
                    (synthetic_compressibility
                         ? 0.03 * pressure_delta
                         : 0.0)};
        const std::array<double, 2>
            mass_density{
                2.0 +
                    (synthetic_compressibility
                         ? 0.02 * pressure_delta
                         : 0.0),
                1.5 +
                    (synthetic_compressibility
                         ? 0.015 * pressure_delta
                         : 0.0)};
        for (std::size_t phase = 0U;
             phase < 2U;
             ++phase) {
            if (!(molar_density[phase] > 0.0) ||
                !(mass_density[phase] > 0.0)) {
                *status =
                    fdp::NaturalVariableSnesEvaluationStatus3D::
                        domain_error;
                return PETSC_SUCCESS;
            }
        }

        input.phase_properties[0] =
            flow::PhasePropertyPrerequisiteInput{
                molar_density[0],
                mass_density[0],
                1.0,
                temperature,
                0.5 * temperature};
        input.phase_properties[1] =
            flow::PhasePropertyPrerequisiteInput{
                molar_density[1],
                mass_density[1],
                1.0,
                1.2 * temperature,
                0.8 * temperature};

        auto state =
            flow::NaturalVariableCellState2P::
                create(std::move(input));
        const std::size_t q =
            frozen_layout.unknown_count();
        std::array<std::vector<double>, 2>
            zero{
                std::vector<double>(q, 0.0),
                std::vector<double>(q, 0.0)};
        std::array<std::vector<double>, 2>
            molar_density_gradient = zero;
        std::array<std::vector<double>, 2>
            mass_density_gradient = zero;
        if (synthetic_compressibility) {
            const auto pressure_column =
                frozen_layout
                    .pressure_unknown_index();
            molar_density_gradient[0][
                pressure_column] =
                0.06;
            molar_density_gradient[1][
                pressure_column] =
                0.03;
            mass_density_gradient[0][
                pressure_column] =
                0.02;
            mass_density_gradient[1][
                pressure_column] =
                0.015;
        }
        auto molar =
            flow::
                make_two_phase_molar_density_linearization(
                    state,
                    molar_density_gradient);
        const flow::TransportPropertyProvenance
            provenance{
                "mixed-physical-dispatch",
                "controlled-regression",
                "v1"};
        auto transport =
            flow::
                make_two_phase_transport_linearization(
                    state,
                    mass_density_gradient,
                    zero,
                    {0.60, 0.40},
                    zero,
                    provenance,
                    provenance);
        std::array<std::vector<double>, 2>
            dh{
                std::vector<double>(q, 0.0),
                std::vector<double>(q, 0.0)};
        std::array<std::vector<double>, 2>
            du{
                std::vector<double>(q, 0.0),
                std::vector<double>(q, 0.0)};
        dh[0][frozen_layout
                  .temperature_unknown_index()] =
            1.0;
        dh[1][frozen_layout
                  .temperature_unknown_index()] =
            1.2;
        du[0][frozen_layout
                  .temperature_unknown_index()] =
            0.5;
        du[1][frozen_layout
                  .temperature_unknown_index()] =
            0.8;
        auto caloric =
            flow::
                make_two_phase_caloric_linearization(
                    state,
                    dh,
                    du,
                    provenance,
                    provenance);
        std::vector<double> rock_gradient(
            q,
            0.0);
        rock_gradient[
            frozen_layout
                .temperature_unknown_index()] =
            2.0;
        auto rock =
            flow::
                make_two_phase_rock_thermal_storage_linearization(
                    state,
                    2.0 * temperature,
                    rock_gradient,
                    provenance);

        flow::TwoPhaseFugacityEquilibriumLinearization
            fugacity;
        fugacity.layout =
            frozen_layout.descriptor();
        fugacity.component_ids.assign(
            component_ids.begin(),
            component_ids.end());
        fugacity.residual.assign(
            3U,
            0.0);
        fugacity.input_count = q;
        fugacity.jacobian.assign(
            3U * q,
            0.0);
        for (std::size_t row = 0U;
             row < 3U;
             ++row) {
            fugacity.jacobian[
                row * q +
                (3U + row)] =
                1.0;
        }

        output->emplace(
            fdp::
                TwoPhaseCurrentCellLinearization3D{
                    std::move(state),
                    std::move(molar),
                    std::move(transport),
                    std::move(caloric),
                    std::move(rock),
                    std::move(fugacity)});
        *status =
            fdp::NaturalVariableSnesEvaluationStatus3D::
                success;
        return PETSC_SUCCESS;
    } catch (const std::invalid_argument&) {
        *status =
            fdp::NaturalVariableSnesEvaluationStatus3D::
                domain_error;
        return PETSC_SUCCESS;
    } catch (...) {
        return PETSC_ERR_ARG_INCOMP;
    }
}

struct ConstantRelativePermeability3P {
    template <typename Number>
    flow::RelativePermeabilityEvaluation3P<Number>
    operator()(
        const flow::ThreePhaseSaturationState3P<Number>&)
        const {
        return {
            std::array<Number, 3>{
                Number{0.50},
                Number{0.30},
                Number{0.20}}};
    }
};

PetscErrorCode
evaluate_3p(
    mesh::LocalIndex,
    mesh::GlobalEntityId,
    std::span<const double> natural_variables,
    const flow::NaturalVariableLayout3P&
        frozen_layout,
    std::span<const std::string> component_ids,
    void* raw_context,
    std::optional<
        fdp::FixedThreePhaseCurrentCellLinearization3D>*
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
        static_cast<DispatchAudit*>(
            raw_context);
    ++audit->three_calls;

    try {
        if (natural_variables.size() !=
                frozen_layout.unknown_count() ||
            component_ids.size() !=
                frozen_layout.component_count()) {
            return PETSC_ERR_ARG_INCOMP;
        }

        const std::size_t q =
            frozen_layout.unknown_count();
        const double p =
            natural_variables[
                frozen_layout
                    .pressure_unknown_index()];
        const double temperature =
            natural_variables[
                frozen_layout
                    .temperature_unknown_index()];
        const auto s0 =
            frozen_layout
                .independent_saturation_unknown_index(
                    flow::PhaseSlot3::phase0);
        const auto s1 =
            frozen_layout
                .independent_saturation_unknown_index(
                    flow::PhaseSlot3::phase1);
        if (!s0 || !s1 ||
            !std::isfinite(p) ||
            !(p > 0.0) ||
            !std::isfinite(temperature) ||
            !(temperature > 0.0)) {
            *status =
                fdp::NaturalVariableSnesEvaluationStatus3D::
                    domain_error;
            return PETSC_SUCCESS;
        }

        flow::NaturalVariableCellStateInput3P
            input;
        input.component_ids.assign(
            component_ids.begin(),
            component_ids.end());
        input.reference_pressure_pa = p;
        input.temperature_k = temperature;
        input.independent_saturations = {
            natural_variables[*s0],
            natural_variables[*s1]};
        input.composition_pivot =
            frozen_layout.composition_pivot();

        for (std::size_t phase = 0U;
             phase < 3U;
             ++phase) {
            const auto slot =
                static_cast<flow::PhaseSlot3>(
                    phase);
            for (std::size_t rank = 0U;
                 rank < 2U;
                 ++rank) {
                const auto component =
                    frozen_layout
                        .independent_composition_component(
                            slot,
                            rank);
                const auto column =
                    frozen_layout
                        .independent_composition_unknown_index(
                            slot,
                            component);
                if (!column) {
                    return PETSC_ERR_PLIB;
                }
                input.independent_phase_compositions[
                    phase]
                    .push_back(
                        natural_variables[
                            *column]);
            }
        }

        const bool synthetic_compressibility =
            audit->well_timestep_compressibility;
        const double pressure_delta =
            p - 30.0;
        const std::array<double, 3>
            molar_density{
                2.0 +
                    (synthetic_compressibility
                         ? 0.02 * pressure_delta
                         : 0.0),
                4.0 +
                    (synthetic_compressibility
                         ? 0.04 * pressure_delta
                         : 0.0),
                7.0 +
                    (synthetic_compressibility
                         ? 0.07 * pressure_delta
                         : 0.0)};
        const std::array<double, 3>
            mass_density{
                1.0 +
                    (synthetic_compressibility
                         ? 0.01 * pressure_delta
                         : 0.0),
                2.0 +
                    (synthetic_compressibility
                         ? 0.02 * pressure_delta
                         : 0.0),
                3.0 +
                    (synthetic_compressibility
                         ? 0.03 * pressure_delta
                         : 0.0)};
        const std::array<double, 3>
            viscosity{
                1.0, 1.2, 1.4};
        for (std::size_t phase = 0U;
             phase < 3U;
             ++phase) {
            if (!(molar_density[phase] > 0.0) ||
                !(mass_density[phase] > 0.0)) {
                *status =
                    fdp::NaturalVariableSnesEvaluationStatus3D::
                        domain_error;
                return PETSC_SUCCESS;
            }
        }
        for (std::size_t phase = 0U;
             phase < 3U;
             ++phase) {
            const double phase_scale =
                static_cast<double>(
                    phase + 1U);
            input.phase_properties[phase] =
                flow::PhasePropertyPrerequisiteInput{
                    molar_density[phase],
                    mass_density[phase],
                    viscosity[phase],
                    phase_scale *
                        temperature,
                    0.5 *
                        phase_scale *
                        temperature};
        }

        auto state =
            flow::NaturalVariableCellState3P::
                create(std::move(input));

        std::array<std::vector<double>, 3>
            zero;
        std::array<std::vector<double>, 3>
            molar_density_gradient;
        std::array<std::vector<double>, 3>
            mass_density_gradient;
        std::array<std::vector<double>, 3>
            dh;
        std::array<std::vector<double>, 3>
            du;
        for (std::size_t phase = 0U;
             phase < 3U;
             ++phase) {
            zero[phase].assign(
                q,
                0.0);
            molar_density_gradient[phase].assign(
                q,
                0.0);
            mass_density_gradient[phase].assign(
                q,
                0.0);
            if (synthetic_compressibility) {
                const auto pressure_column =
                    frozen_layout
                        .pressure_unknown_index();
                const double phase_scale =
                    static_cast<double>(
                        phase + 1U);
                molar_density_gradient[phase][
                    pressure_column] =
                    0.02 * phase_scale;
                mass_density_gradient[phase][
                    pressure_column] =
                    0.01 * phase_scale;
            }
            dh[phase].assign(
                q,
                0.0);
            du[phase].assign(
                q,
                0.0);
            const double phase_scale =
                static_cast<double>(
                    phase + 1U);
            dh[phase][
                frozen_layout
                    .temperature_unknown_index()] =
                phase_scale;
            du[phase][
                frozen_layout
                    .temperature_unknown_index()] =
                0.5 * phase_scale;
        }

        flow::
            PhaseMolarDensityNaturalVariableLinearization3P
            molar{
                frozen_layout,
                molar_density,
                molar_density_gradient};
        const flow::TransportPropertyProvenance
            provenance{
                "mixed-physical-dispatch",
                "controlled-regression",
                "v1"};
        auto transport =
            flow::
                make_phase_transport_property_linearization(
                    state,
                    mass_density_gradient,
                    zero,
                    provenance,
                    provenance);
        auto caloric =
            flow::
                make_phase_caloric_property_linearization(
                    state,
                    dh,
                    du,
                    provenance,
                    provenance);
        std::vector<double> rock_gradient(
            q,
            0.0);
        rock_gradient[
            frozen_layout
                .temperature_unknown_index()] =
            2.0;
        auto rock =
            flow::
                make_stationary_rock_thermal_storage_linearization(
                    state,
                    2.0 * temperature,
                    rock_gradient,
                    provenance);

        const auto saturation =
            flow::
                evaluate_three_phase_saturation_constitutive(
                    p,
                    natural_variables[*s0],
                    natural_variables[*s1],
                    ConstantRelativePermeability3P{},
                    flow::NoCapillaryPressure3P{});
        flow::
            ThreePhaseSaturationCoordinateDerivatives3P
            saturation_derivatives{};
        auto saturation_linearization =
            flow::
                make_saturation_constitutive_natural_variable_linearization(
                    state,
                    saturation,
                    saturation_derivatives);

        std::vector<double> fugacity_values(
            6U,
            0.0);
        std::vector<double> fugacity_jacobian(
            6U * q,
            0.0);
        for (std::size_t row = 0U;
             row < 6U;
             ++row) {
            fugacity_jacobian[
                row * q +
                (4U + row)] =
                1.0;
        }
        flow::
            FugacityEquilibriumResidualLinearization3P
            fugacity{
                frozen_layout,
                std::vector<std::string>{
                    component_ids.begin(),
                    component_ids.end()},
                flow::
                    FugacityEquilibriumResidual3P<double>{
                        3U,
                        std::move(
                            fugacity_values)},
                q,
                std::move(
                    fugacity_jacobian)};

        output->emplace(
            fdp::
                FixedThreePhaseCurrentCellLinearization3D{
                    std::move(state),
                    std::move(molar),
                    std::move(transport),
                    std::move(caloric),
                    std::move(rock),
                    std::move(
                        saturation_linearization),
                    std::move(fugacity)});
        *status =
            fdp::NaturalVariableSnesEvaluationStatus3D::
                success;
        return PETSC_SUCCESS;
    } catch (const std::invalid_argument&) {
        *status =
            fdp::NaturalVariableSnesEvaluationStatus3D::
                domain_error;
        return PETSC_SUCCESS;
    } catch (...) {
        return PETSC_ERR_ARG_INCOMP;
    }
}

[[nodiscard]] mesh::Topology
make_topology() {
    mesh::Topology::EntityIds ids;
    ids.faces = {
        mesh::GlobalEntityId{UINT64_C(100)},
        mesh::GlobalEntityId{UINT64_C(200)},
        mesh::GlobalEntityId{UINT64_C(300)},
        mesh::GlobalEntityId{UINT64_C(400)}};
    ids.cells = {
        mesh::GlobalEntityId{UINT64_C(10)},
        mesh::GlobalEntityId{UINT64_C(20)},
        mesh::GlobalEntityId{UINT64_C(30)},
        mesh::GlobalEntityId{UINT64_C(40)},
        mesh::GlobalEntityId{UINT64_C(50)},
        mesh::GlobalEntityId{UINT64_C(60)}};
    return {
        std::move(ids),
        {}};
}

[[nodiscard]] mesh::PartitionSnapshot
make_partition(
    int rank) {
    auto topology =
        make_topology();
    mesh::EntityOwnerRanks owners;
    owners.faces = {
        mesh::PartitionRank{0U},
        mesh::PartitionRank{0U},
        mesh::PartitionRank{0U},
        mesh::PartitionRank{0U}};
    owners.cells = {
        mesh::PartitionRank{0U},
        mesh::PartitionRank{1U},
        mesh::PartitionRank{0U},
        mesh::PartitionRank{1U},
        mesh::PartitionRank{0U},
        mesh::PartitionRank{1U}};
    return mesh::PartitionSnapshot::create(
        topology,
        mesh::PartitionRank{
            static_cast<
                mesh::PartitionRank::value_type>(
                    rank)},
        2U,
        std::move(owners));
}

[[nodiscard]]
dp::AssemblyReadyInternalConnectionRow3D
connection(
    std::size_t face,
    std::size_t owner,
    std::size_t neighbour,
    std::uint64_t face_global,
    std::uint64_t owner_global,
    std::uint64_t neighbour_global) {
    return {
        mesh::LocalIndex{
            static_cast<
                mesh::LocalIndex::value_type>(
                    face)},
        mesh::GlobalEntityId{
            face_global},
        mesh::LocalIndex{
            static_cast<
                mesh::LocalIndex::value_type>(
                    owner)},
        mesh::GlobalEntityId{
            owner_global},
        mesh::LocalIndex{
            static_cast<
                mesh::LocalIndex::value_type>(
                    neighbour)},
        mesh::GlobalEntityId{
            neighbour_global},
        2.0e-12};
}

[[nodiscard]]
dp::ParallelOwnedConnectionSchedule3D
make_schedule(
    int rank,
    bool include_cross_cardinality) {
    std::vector<
        dp::AssemblyReadyInternalConnectionRow3D>
        rows{
            connection(
                0U, 0U, 1U,
                UINT64_C(100),
                UINT64_C(10),
                UINT64_C(20)),
            connection(
                1U, 2U, 3U,
                UINT64_C(200),
                UINT64_C(30),
                UINT64_C(40)),
            connection(
                2U, 4U, 5U,
                UINT64_C(300),
                UINT64_C(50),
                UINT64_C(60))};
    if (include_cross_cardinality) {
        rows.push_back(
            connection(
                3U, 1U, 2U,
                UINT64_C(400),
                UINT64_C(20),
                UINT64_C(30)));
    }

    if (rank == 0) {
        return {
            mesh::PartitionRank{0U},
            2U,
            std::move(rows),
            {}};
    }
    return {
        mesh::PartitionRank{1U},
        2U,
        {},
        std::move(rows)};
}

[[nodiscard]]
dp::PetscMpiAijSymbolicPreallocation3D
make_cell_bridge(
    int rank) {
    const std::vector<PetscInt>
        local_cell_rows{
            0, 3, 1, 4, 2, 5};
    if (rank == 0) {
        return {
            mesh::PartitionRank{0U},
            2U,
            0,
            3,
            6,
            {
                mesh::LocalIndex{0U},
                mesh::LocalIndex{2U},
                mesh::LocalIndex{4U}},
            {
                mesh::GlobalEntityId{UINT64_C(10)},
                mesh::GlobalEntityId{UINT64_C(30)},
                mesh::GlobalEntityId{UINT64_C(50)}},
            {0, 1, 2},
            {1, 1, 1},
            {1, 2, 1},
            local_cell_rows};
    }
    return {
        mesh::PartitionRank{1U},
        2U,
        3,
        6,
        6,
        {
            mesh::LocalIndex{1U},
            mesh::LocalIndex{3U},
            mesh::LocalIndex{5U}},
        {
            mesh::GlobalEntityId{UINT64_C(20)},
            mesh::GlobalEntityId{UINT64_C(40)},
            mesh::GlobalEntityId{UINT64_C(60)}},
        {3, 4, 5},
        {1, 1, 1},
        {2, 1, 1},
        local_cell_rows};
}

[[nodiscard]]
dp::OwnedCellStructuralColumnPatternSnapshot3D
make_cell_pattern(
    int rank) {
    if (rank == 0) {
        return {
            mesh::PartitionRank{0U},
            2U,
            6U,
            0,
            3,
            6,
            {
                {
                    mesh::LocalIndex{0U},
                    mesh::GlobalEntityId{UINT64_C(10)},
                    0, 0U, 1U, 0U, 1U},
                {
                    mesh::LocalIndex{2U},
                    mesh::GlobalEntityId{UINT64_C(30)},
                    1, 1U, 1U, 1U, 2U},
                {
                    mesh::LocalIndex{4U},
                    mesh::GlobalEntityId{UINT64_C(50)},
                    2, 2U, 1U, 3U, 1U}},
            {0, 1, 2},
            {3, 3, 4, 5}};
    }
    return {
        mesh::PartitionRank{1U},
        2U,
        6U,
        3,
        6,
        6,
        {
            {
                mesh::LocalIndex{1U},
                mesh::GlobalEntityId{UINT64_C(20)},
                3, 0U, 1U, 0U, 2U},
            {
                mesh::LocalIndex{3U},
                mesh::GlobalEntityId{UINT64_C(40)},
                4, 1U, 1U, 2U, 1U},
            {
                mesh::LocalIndex{5U},
                mesh::GlobalEntityId{UINT64_C(60)},
                5, 2U, 1U, 3U, 1U}},
        {3, 4, 5},
        {0, 1, 1, 2}};
}

[[nodiscard]]
disc::CombinedTransmissibilityAdmissibility3D
direct_admissibility() {
    return {
        disc::
            CombinedTransmissibilityAdmissibilityDisposition3D::
                direct_normal_projection_k_orthogonal_candidate,
        {
            disc::
                TransmissibilityGeometryDisposition3D::
                    direct_normal_projection_allowed,
            0.0,
            0.05},
        {
            disc::
                KOrthogonalityDisposition3D::
                    k_orthogonal_within_policy,
            std::nullopt,
            std::nullopt,
            {
                std::nullopt,
                std::nullopt},
            0.05}};
}

[[nodiscard]]
disc::TpfaInternalFaceTransmissibilityEntry3D
materialized_entry(
    std::size_t face) {
    return {
        mesh::LocalIndex{
            static_cast<
                mesh::LocalIndex::value_type>(
                    face)},
        disc::
            TpfaInternalFaceTransmissibilityDisposition3D::
                materialized,
        direct_admissibility(),
        disc::TpfaStaticFaceTransmissibility3D{
            disc::
                TpfaStaticFaceTransmissibilityDisposition3D::
                    positive_harmonic_combination,
            1.0,
            2.0e-12}};
}

[[nodiscard]]
fdp::MixedCardinalityPhysicalSnesAuthoritativeFaceInput3D
make_face_input(
    std::size_t face,
    std::uint64_t face_global) {
    return {
        mesh::LocalIndex{
            static_cast<
                mesh::LocalIndex::value_type>(
                    face)},
        mesh::GlobalEntityId{
            face_global},
        materialized_entry(face),
        flow::GravityVector3D{
            0.0, 0.0, 0.0},
        flow::OwnerToNeighbourDisplacement3D{
            1.0, 0.0, 0.0},
        mpmc::flow_discretization::
            StaticThermalFaceConductance3D{
                0.0}};
}

[[nodiscard]]
std::vector<
    fdp::
        MixedCardinalityPhysicalSnesAuthoritativeFaceInput3D>
make_face_inputs(
    int rank,
    bool include_cross_cardinality) {
    if (rank != 0) {
        return {};
    }
    std::vector<
        fdp::
            MixedCardinalityPhysicalSnesAuthoritativeFaceInput3D>
        result{
            make_face_input(
                0U,
                UINT64_C(100)),
            make_face_input(
                1U,
                UINT64_C(200)),
            make_face_input(
                2U,
                UINT64_C(300))};
    if (include_cross_cardinality) {
        result.push_back(
            make_face_input(
                3U,
                UINT64_C(400)));
    }
    return result;
}

[[nodiscard]]
fdp::MixedCardinalityPhysicalCurrentCellLinearization3D
evaluate_target(
    std::uint64_t stable,
    DispatchAudit* audit) {
    const std::vector<std::string>
        ids{"A", "B", "C"};
    const auto q =
        target_state(stable);
    fdp::NaturalVariableSnesEvaluationStatus3D
        status =
            fdp::
                NaturalVariableSnesEvaluationStatus3D::
                    success;

    if (phase_count_for(stable) == 1U) {
        std::optional<
            fdp::SinglePhaseCurrentCellLinearization3D>
            output;
        auto layout = layout_1p();
        const auto error =
            evaluate_1p(
                mesh::LocalIndex{0U},
                mesh::GlobalEntityId{stable},
                q,
                layout,
                ids,
                audit,
                &output,
                &status);
        if (error != PETSC_SUCCESS ||
            status !=
                fdp::
                    NaturalVariableSnesEvaluationStatus3D::
                        success ||
            !output.has_value()) {
            throw std::runtime_error(
                "failed direct 1P target evaluation");
        }
        return std::move(*output);
    }

    if (phase_count_for(stable) == 2U) {
        std::optional<
            fdp::TwoPhaseCurrentCellLinearization3D>
            output;
        auto layout = layout_2p();
        const auto error =
            evaluate_2p(
                mesh::LocalIndex{0U},
                mesh::GlobalEntityId{stable},
                q,
                layout,
                ids,
                audit,
                &output,
                &status);
        if (error != PETSC_SUCCESS ||
            status !=
                fdp::
                    NaturalVariableSnesEvaluationStatus3D::
                        success ||
            !output.has_value()) {
            throw std::runtime_error(
                "failed direct 2P target evaluation");
        }
        return std::move(*output);
    }

    std::optional<
        fdp::FixedThreePhaseCurrentCellLinearization3D>
        output;
    auto layout = layout_3p();
    const auto error =
        evaluate_3p(
            mesh::LocalIndex{0U},
            mesh::GlobalEntityId{stable},
            q,
            layout,
            ids,
            audit,
            &output,
            &status);
    if (error != PETSC_SUCCESS ||
        status !=
            fdp::
                NaturalVariableSnesEvaluationStatus3D::
                    success ||
        !output.has_value()) {
        throw std::runtime_error(
            "failed direct 3P target evaluation");
    }
    return std::move(*output);
}

[[nodiscard]]
std::vector<
    fdp::MixedCardinalityPhysicalSnesCellInput3D>
make_cell_inputs(
    int rank,
    DispatchAudit* audit) {
    std::vector<
        fdp::MixedCardinalityPhysicalSnesCellInput3D>
        result;
    result.reserve(6U);

    for (std::size_t local = 0U;
         local < 6U;
         ++local) {
        const std::uint64_t stable =
            static_cast<std::uint64_t>(
                (local + 1U) * 10U);
        const bool owned =
            (rank == 0 &&
             local % 2U == 0U) ||
            (rank == 1 &&
             local % 2U == 1U);
        const double bulk_volume =
            2.0 +
            static_cast<double>(local);
        const double porosity =
            0.25 +
            0.005 *
                static_cast<double>(
                    local);
        auto current =
            evaluate_target(
                stable,
                audit);

        if (phase_count_for(stable) == 1U) {
            fdp::SinglePhaseSnesCellInput3D
                input;
            input.cell =
                mesh::LocalIndex{
                    static_cast<
                        mesh::LocalIndex::value_type>(
                            local)};
            input.cell_global =
                mesh::GlobalEntityId{stable};
            input.bulk_volume_m3 =
                bulk_volume;
            input.porosity =
                porosity;
            input.frozen_layout =
                layout_1p();
            input.component_ids =
                {"A", "B", "C"};
            if (owned) {
                const auto& typed =
                    std::get<
                        fdp::
                            SinglePhaseCurrentCellLinearization3D>(
                                current);
                input.previous_component_accumulation =
                    flow::
                        build_single_phase_component_accumulation(
                            typed.state,
                            porosity);
                input.previous_energy_accumulation =
                    flow::
                        build_single_phase_energy_accumulation_snapshot(
                            typed.state,
                            porosity,
                            typed.transport,
                            typed.caloric,
                            typed.rock);
            }
            result.emplace_back(
                std::move(input));
            continue;
        }

        if (phase_count_for(stable) == 2U) {
            fdp::TwoPhaseSnesCellInput3D
                input;
            input.cell =
                mesh::LocalIndex{
                    static_cast<
                        mesh::LocalIndex::value_type>(
                            local)};
            input.cell_global =
                mesh::GlobalEntityId{stable};
            input.bulk_volume_m3 =
                bulk_volume;
            input.porosity =
                porosity;
            input.frozen_layout =
                layout_2p();
            input.component_ids =
                {"A", "B", "C"};
            if (owned) {
                const auto& typed =
                    std::get<
                        fdp::
                            TwoPhaseCurrentCellLinearization3D>(
                                current);
                input.previous_component_accumulation =
                    flow::
                        build_two_phase_component_accumulation(
                            typed.state,
                            porosity);
                input.previous_energy_accumulation =
                    flow::
                        build_two_phase_energy_accumulation_snapshot(
                            typed.state,
                            porosity,
                            typed.transport,
                            typed.caloric,
                            typed.rock);
            }
            result.emplace_back(
                std::move(input));
            continue;
        }

        fdp::FixedThreePhaseSnesCellInput3D
            input;
        input.cell =
            mesh::LocalIndex{
                static_cast<
                    mesh::LocalIndex::value_type>(
                        local)};
        input.cell_global =
            mesh::GlobalEntityId{stable};
        input.bulk_volume_m3 =
            bulk_volume;
        input.porosity =
            porosity;
        input.frozen_layout =
            layout_3p();
        input.component_ids =
            {"A", "B", "C"};
        if (owned) {
            const auto& typed =
                std::get<
                    fdp::
                        FixedThreePhaseCurrentCellLinearization3D>(
                            current);
            input.previous_component_accumulation =
                flow::
                    build_pore_volume_component_accumulation(
                        typed.state,
                        porosity);
            input.previous_energy_accumulation =
                flow::
                    build_pore_volume_energy_accumulation_snapshot(
                        typed.state,
                        porosity,
                        typed.transport,
                        typed.caloric,
                        typed.rock);
        }
        result.emplace_back(
            std::move(input));
    }
    return result;
}

[[nodiscard]]
std::vector<
    flow::FrozenActivePhaseIdentityMap>
make_phase_identity_maps() {
    const flow::FrozenPhysicalPhaseIdentity
        aqueous{
            "fixture/mixed-physical-dispatch",
            "aqueous"};
    const flow::FrozenPhysicalPhaseIdentity
        hydrocarbon0{
            "fixture/mixed-physical-dispatch",
            "hydrocarbon-0"};
    const flow::FrozenPhysicalPhaseIdentity
        hydrocarbon1{
            "fixture/mixed-physical-dispatch",
            "hydrocarbon-1"};

    return {
        flow::FrozenActivePhaseIdentityMap{
            {aqueous}},
        flow::FrozenActivePhaseIdentityMap{
            {aqueous}},
        flow::FrozenActivePhaseIdentityMap{
            {
                aqueous,
                hydrocarbon0}},
        flow::FrozenActivePhaseIdentityMap{
            {
                aqueous,
                hydrocarbon0}},
        flow::FrozenActivePhaseIdentityMap{
            {
                aqueous,
                hydrocarbon0,
                hydrocarbon1}},
        flow::FrozenActivePhaseIdentityMap{
            {
                aqueous,
                hydrocarbon0,
                hydrocarbon1}}};
}


[[nodiscard]]
flow::FrozenActivePhaseIdentityMap
phase_identity_map_for_stable(
    std::uint64_t stable) {
    const auto maps =
        make_phase_identity_maps();
    const std::size_t index =
        static_cast<std::size_t>(
            stable / UINT64_C(10) -
            UINT64_C(1));
    return maps.at(index);
}


