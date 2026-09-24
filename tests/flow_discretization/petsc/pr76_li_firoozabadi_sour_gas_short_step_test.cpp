#include <mpmc/flow/pr76_li_firoozabadi_sour_gas_properties.hpp>
#include <mpmc/flow_discretization_petsc/mixed_cardinality_physical_snes_assembly.hpp>
#include <mpmc/flow_discretization_petsc/pr76_production_cell_evaluator.hpp>
#include <mpmc/well/peaceman_well_index_3d.hpp>
#include <mpmc/well_discretization_petsc/fixed_bhp_well_source_evaluator.hpp>

#include "../../flash/pr76_three_phase/sour_gas_fixture.hpp"
#include "../../flash/pr76_three_phase/sour_gas_references.hpp"
#include "../../flow/core/pr76_li_firoozabadi_sour_gas_flow_reference.hpp"

#include <petscmat.h>
#include <petscsnes.h>
#include <petscvec.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <optional>
#include <span>
#include <stdexcept>
#include <string>
#include <string_view>
#include <utility>
#include <variant>
#include <vector>

namespace {

namespace dp = mpmc::discretization_petsc;
namespace fdp = mpmc::flow_discretization_petsc;
namespace flow = mpmc::flow;
namespace mesh = mpmc::mesh;
namespace th = mpmc::thermodynamics;
namespace wd = mpmc::well_discretization;
namespace wdp = mpmc::well_discretization_petsc;
namespace well = mpmc::well;
namespace fx = pr76_sour_gas_test;
namespace eq_ref = pr76_sour_gas_reference;
namespace ref =
    pr76_li_firoozabadi_sour_gas_flow_reference;

void require(
    bool condition,
    std::string_view message) {
    if (!condition) {
        throw std::runtime_error(
            std::string{message});
    }
}

void near(
    double actual,
    double expected,
    double relative,
    double absolute,
    std::string_view message) {
    const double scale =
        std::max(
            {1.0,
             std::abs(actual),
             std::abs(expected)});
    require(
        std::isfinite(actual) &&
            std::isfinite(expected) &&
            std::abs(actual - expected) <=
                absolute + relative * scale,
        message);
}

std::vector<th::Pr76SelectedPhase>
reference_selections(
    const th::Pr76Phase<double>& model) {
    std::vector<th::Pr76SelectedPhase>
        result;
    result.reserve(3U);

    for (std::size_t phase = 0U;
         phase < 3U;
         ++phase) {
        const std::vector<double>
            composition{
                eq_ref::phases[phase].begin(),
                eq_ref::phases[phase].end()};
        th::Pr76PhaseWorkspace<double>
            workspace;
        const auto roots =
            model.roots_full(
                ref::pressure_pa,
                ref::temperature_k,
                composition,
                workspace);
        require(
            roots.status ==
                    th::Pr76RootStatus::success &&
                roots.count > 0U,
            "Li-Firoozabadi reference phase has no resolved PR76 roots");

        std::size_t best = 0U;
        double best_error =
            std::numeric_limits<double>::
                infinity();
        for (std::size_t root = 0U;
             root < roots.count;
             ++root) {
            const double error =
                std::abs(
                    roots.roots[root].z -
                    eq_ref::
                        compressibility_factors[
                            phase]);
            if (error < best_error) {
                best_error = error;
                best = root;
            }
        }
        require(
            best_error <= 5.0e-12,
            "Li-Firoozabadi reference Z does not identify a production PR76 root");
        result.push_back(
            {best, {}});
    }
    return result;
}

flow::NaturalVariableLayout3P
reference_layout() {
    const std::array<std::vector<double>, 3>
        compositions{
            std::vector<double>{
                eq_ref::phases[0].begin(),
                eq_ref::phases[0].end()},
            std::vector<double>{
                eq_ref::phases[1].begin(),
                eq_ref::phases[1].end()},
            std::vector<double>{
                eq_ref::phases[2].begin(),
                eq_ref::phases[2].end()}};
    return flow::NaturalVariableLayout3P{
        flow::
            NaturalVariableCompositionPivot3P::
                select(compositions)};
}

std::vector<double> reference_q(
    const flow::NaturalVariableLayout3P&
        layout) {
    std::vector<double>
        q(
            layout.unknown_count(),
            0.0);
    q[layout.pressure_unknown_index()] =
        ref::pressure_pa;
    q[layout.temperature_unknown_index()] =
        ref::temperature_k;

    const auto s0 =
        layout.independent_saturation_unknown_index(
            flow::PhaseSlot3::phase0);
    const auto s1 =
        layout.independent_saturation_unknown_index(
            flow::PhaseSlot3::phase1);
    require(
        s0.has_value() &&
            s1.has_value(),
        "Li-Firoozabadi 3P layout lost saturation columns");
    q[*s0] =
        ref::phase_saturation[0];
    q[*s1] =
        ref::phase_saturation[1];

    for (std::size_t phase = 0U;
         phase < 3U;
         ++phase) {
        const auto slot =
            static_cast<flow::PhaseSlot3>(
                phase);
        for (std::size_t component = 0U;
             component < 6U;
             ++component) {
            const auto column =
                layout
                    .independent_composition_unknown_index(
                        slot,
                        component);
            if (column.has_value()) {
                q[*column] =
                    eq_ref::phases[
                        phase][component];
            }
        }
    }
    return q;
}

/// Software-structure-only carrier.  This is deliberately NOT part of the
/// external V-L1-L2 physical validation.  With no internal faces and
/// p_bhp == p_phase, every validated rate is independent of kr.
struct StructuralOnlyLinearRelativePermeability3P {
    template <typename Number>
    [[nodiscard]]
    flow::RelativePermeabilityEvaluation3P<Number>
    operator()(
        const flow::
            ThreePhaseSaturationState3P<Number>&
                state) const {
        return {
            state.saturation};
    }
};

PetscErrorCode evaluate_structural_saturation(
    const flow::NaturalVariableCellState3P&
        state,
    void*,
    std::optional<
        flow::
            ThreePhaseSaturationConstitutiveNaturalVariableLinearization3P>*
                output,
    fdp::NaturalVariableSnesEvaluationStatus3D*
        status) {
    if (output == nullptr ||
        status == nullptr) {
        return PETSC_ERR_ARG_NULL;
    }

    try {
        const auto primal =
            flow::
                evaluate_three_phase_saturation_constitutive(
                    state,
                    StructuralOnlyLinearRelativePermeability3P{},
                    flow::NoCapillaryPressure3P{});

        flow::
            ThreePhaseSaturationCoordinateDerivatives3P
                derivatives;
        derivatives.relative_permeability[0] =
            {1.0, 0.0};
        derivatives.relative_permeability[1] =
            {0.0, 1.0};
        derivatives.relative_permeability[2] =
            {-1.0, -1.0};
        derivatives.capillary_pressure_offset_pa[0] =
            {0.0, 0.0};
        derivatives.capillary_pressure_offset_pa[1] =
            {0.0, 0.0};
        derivatives.capillary_pressure_offset_pa[2] =
            {0.0, 0.0};

        output->emplace(
            flow::
                make_saturation_constitutive_natural_variable_linearization(
                    state,
                    primal,
                    derivatives));
        *status =
            fdp::
                NaturalVariableSnesEvaluationStatus3D::
                    success;
        return PETSC_SUCCESS;
    } catch (const std::invalid_argument&) {
        output->reset();
        *status =
            fdp::
                NaturalVariableSnesEvaluationStatus3D::
                    domain_error;
        return PETSC_SUCCESS;
    }
}

PetscErrorCode evaluate_quartz_rock_storage(
    const flow::NaturalVariableLayoutDescriptor&
        layout,
    std::span<const double>
        natural_variables,
    void*,
    std::optional<
        fdp::
            Pr76RockThermalStorageLinearization3D>*
                output,
    fdp::NaturalVariableSnesEvaluationStatus3D*
        status) {
    if (output == nullptr ||
        status == nullptr) {
        return PETSC_ERR_ARG_NULL;
    }
    if (natural_variables.size() !=
        layout.unknown_count()) {
        return PETSC_ERR_ARG_SIZ;
    }

    constexpr double t0 = 169.1;
    constexpr double t1 = 184.8;
    constexpr double cp0 = 469.3;
    constexpr double cp1 = 510.6;
    constexpr double rho_rock = 2637.8;
    constexpr double t_reference =
        ref::temperature_k;

    const double temperature =
        natural_variables[
            layout.temperature_unknown_index()];
    if (!std::isfinite(temperature) ||
        temperature < t0 ||
        temperature > t1) {
        output->reset();
        *status =
            fdp::
                NaturalVariableSnesEvaluationStatus3D::
                    domain_error;
        return PETSC_SUCCESS;
    }

    const double cp_slope =
        (cp1 - cp0) /
        (t1 - t0);
    const double cp_reference =
        cp0 +
        cp_slope *
            (t_reference - t0);
    const double delta =
        temperature -
        t_reference;
    const double specific_energy =
        cp_reference * delta +
        0.5 * cp_slope *
            delta * delta;

    std::vector<double>
        gradient(
            layout.unknown_count(),
            0.0);
    gradient[
        layout.temperature_unknown_index()] =
        rho_rock *
        (cp_reference +
         cp_slope * delta);

    output->emplace(
        fdp::
            Pr76RockThermalStorageLinearization3D{
                rho_rock *
                    specific_energy,
                std::move(gradient),
                {
                    "NIST SRD 30 quartz low-temperature heat capacity",
                    "NIST-SRD30-Z00788-Anderson-1936",
                    "rho=2.6378g/cm3__Cp-linear-169.1K-184.8K__u(178.8K)=0-v1"}});
    *status =
        fdp::
            NaturalVariableSnesEvaluationStatus3D::
                success;
    return PETSC_SUCCESS;
}

PetscErrorCode unsupported_single_phase(
    mesh::LocalIndex,
    mesh::GlobalEntityId,
    std::span<const double>,
    const flow::NaturalVariableLayout1P&,
    std::span<const std::string>,
    void*,
    std::optional<
        fdp::SinglePhaseCurrentCellLinearization3D>*
            output,
    fdp::NaturalVariableSnesEvaluationStatus3D*
        status) {
    if (output == nullptr ||
        status == nullptr) {
        return PETSC_ERR_ARG_NULL;
    }
    output->reset();
    *status =
        fdp::
            NaturalVariableSnesEvaluationStatus3D::
                success;
    return PETSC_ERR_SUP;
}

PetscErrorCode unsupported_two_phase(
    mesh::LocalIndex,
    mesh::GlobalEntityId,
    std::span<const double>,
    const flow::NaturalVariableLayout2P&,
    std::span<const std::string>,
    void*,
    std::optional<
        fdp::TwoPhaseCurrentCellLinearization3D>*
            output,
    fdp::NaturalVariableSnesEvaluationStatus3D*
        status) {
    if (output == nullptr ||
        status == nullptr) {
        return PETSC_ERR_ARG_NULL;
    }
    output->reset();
    *status =
        fdp::
            NaturalVariableSnesEvaluationStatus3D::
                success;
    return PETSC_ERR_SUP;
}

mesh::Topology make_topology() {
    mesh::Topology::EntityIds ids;
    ids.cells = {
        mesh::GlobalEntityId{
            UINT64_C(10)}};
    return {
        std::move(ids),
        {}};
}

mesh::PartitionSnapshot make_partition() {
    const auto topology =
        make_topology();
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
make_schedule() {
    return {
        mesh::PartitionRank{0U},
        1U,
        {},
        {}};
}

dp::PetscMpiAijSymbolicPreallocation3D
make_cell_bridge() {
    return {
        mesh::PartitionRank{0U},
        1U,
        0,
        1,
        1,
        {mesh::LocalIndex{0U}},
        {
            mesh::GlobalEntityId{
                UINT64_C(10)}},
        {0},
        {1},
        {0},
        {0}};
}

dp::OwnedCellStructuralColumnPatternSnapshot3D
make_cell_pattern() {
    return {
        mesh::PartitionRank{0U},
        1U,
        1U,
        0,
        1,
        1,
        {
            {
                mesh::LocalIndex{0U},
                mesh::GlobalEntityId{
                    UINT64_C(10)},
                0,
                0U,
                1U,
                0U,
                0U}
        },
        {0},
        {}};
}

flow::FrozenActivePhaseIdentityMap
make_phase_identities() {
    return flow::FrozenActivePhaseIdentityMap{
        std::vector<flow::FrozenPhysicalPhaseIdentity>{
            {
                "Li-Firoozabadi-2012-Table9",
                "equilibrium-phase-0"},
            {
                "Li-Firoozabadi-2012-Table9",
                "equilibrium-phase-1"},
            {
                "Li-Firoozabadi-2012-Table9",
                "equilibrium-phase-2"}}};
}

template <typename Closure>
fdp::FixedThreePhaseCurrentCellLinearization3D
evaluate_current(
    std::span<const double> q,
    const flow::NaturalVariableLayout3P&
        layout,
    std::span<const std::string>
        component_ids,
    fdp::
        Pr76ThreePhaseProductionCellEvaluatorContext3D<
            Closure>* context) {
    std::optional<
        fdp::
            FixedThreePhaseCurrentCellLinearization3D>
        output;
    fdp::NaturalVariableSnesEvaluationStatus3D
        status =
            fdp::
                NaturalVariableSnesEvaluationStatus3D::
                    success;
    const PetscErrorCode error =
        fdp::
            evaluate_pr76_three_phase_production_cell_3d<
                Closure>(
                mesh::LocalIndex{0U},
                mesh::GlobalEntityId{
                    UINT64_C(10)},
                q,
                layout,
                component_ids,
                context,
                &output,
                &status);
    require(
        error == PETSC_SUCCESS &&
            status ==
                fdp::
                    NaturalVariableSnesEvaluationStatus3D::
                        success &&
            output.has_value(),
        "Li-Firoozabadi production cell evaluation failed");
    return std::move(*output);
}

void insert_state(
    Vec state,
    const fdp::
        VariableCardinalityNaturalVariableNumbering3D&
            numbering,
    std::span<const double> q) {
    const auto& cell =
        numbering.cell(
            mesh::LocalIndex{0U});
    require(
        cell.scalar_count == q.size(),
        "Li-Firoozabadi numbering width mismatch");

    for (std::size_t slot = 0U;
         slot < q.size();
         ++slot) {
        const PetscInt index =
            cell.petsc_global_scalar_start +
            static_cast<PetscInt>(slot);
        const PetscScalar value =
            static_cast<PetscScalar>(
                q[slot]);
        require(
            VecSetValues(
                state,
                1,
                &index,
                &value,
                INSERT_VALUES) ==
                PETSC_SUCCESS,
            "failed to insert Li-Firoozabadi initial state");
    }
    require(
        VecAssemblyBegin(state) ==
                PETSC_SUCCESS &&
            VecAssemblyEnd(state) ==
                PETSC_SUCCESS,
        "failed to assemble Li-Firoozabadi initial state");
}

std::vector<double> read_state(
    Vec state,
    const fdp::
        VariableCardinalityNaturalVariableNumbering3D&
            numbering) {
    const auto& cell =
        numbering.cell(
            mesh::LocalIndex{0U});
    std::vector<double>
        q(
            cell.scalar_count,
            0.0);

    for (std::size_t slot = 0U;
         slot < q.size();
         ++slot) {
        const PetscInt index =
            cell.petsc_global_scalar_start +
            static_cast<PetscInt>(slot);
        PetscScalar value{};
        require(
            VecGetValues(
                state,
                1,
                &index,
                &value) ==
                PETSC_SUCCESS,
            "failed to read Li-Firoozabadi solved state");
        q[slot] =
            static_cast<double>(
                PetscRealPart(value));
    }
    return q;
}

} // namespace

void pr76_li_firoozabadi_sour_gas_short_step_test() {
    int rank = -1;
    int size = -1;
    require(
        MPI_Comm_rank(
            PETSC_COMM_SELF,
            &rank) ==
                MPI_SUCCESS &&
            MPI_Comm_size(
                PETSC_COMM_SELF,
                &size) ==
                MPI_SUCCESS &&
            rank == 0 &&
            size == 1,
        "Li-Firoozabadi short-step benchmark requires one MPI rank");

    const auto model =
        fx::model();
    const auto selections =
        reference_selections(model);
    auto closure =
        flow::
            make_pr76_li_firoozabadi_sour_gas_property_closure(
                model,
                selections);
    using Closure =
        decltype(closure);

    fdp::
        Pr76ThreePhaseProductionCellEvaluatorContext3D<
            Closure>
        evaluator_context{
            &closure,
            {
                &evaluate_structural_saturation,
                nullptr},
            {
                &evaluate_quartz_rock_storage,
                nullptr},
            {}};

    const auto layout =
        reference_layout();
    const auto q =
        reference_q(layout);
    const auto component_ids =
        fx::canonical_order();
    auto initial_current =
        evaluate_current(
            q,
            layout,
            component_ids,
            &evaluator_context);

    for (std::size_t phase = 0U;
         phase < 3U;
         ++phase) {
        const auto slot =
            static_cast<flow::PhaseSlot3>(
                phase);
        near(
            initial_current
                .state
                .phase_saturation(slot),
            ref::phase_saturation[phase],
            2.0e-12,
            2.0e-13,
            "independent phase-volume saturation mismatch");
        near(
            initial_current
                .molar_density
                .molar_density_mol_per_m3[
                    phase],
            ref::
                molar_density_mol_per_m3[
                    phase],
            4.0e-10,
            2.0e-8,
            "independent molar-density mismatch");
        near(
            initial_current
                .transport
                .mass_density_kg_per_m3[
                    phase],
            ref::
                mass_density_kg_per_m3[
                    phase],
            4.0e-10,
            2.0e-8,
            "independent mass-density mismatch");
        near(
            initial_current
                .transport
                .dynamic_viscosity_pa_s[
                    phase],
            ref::
                dynamic_viscosity_pa_s[
                    phase],
            4.0e-10,
            2.0e-12,
            "independent LBC viscosity mismatch");
        near(
            initial_current
                .caloric
                .specific_enthalpy_j_per_kg[
                    phase],
            ref::
                specific_enthalpy_j_per_kg[
                    phase],
            4.0e-10,
            2.0e-7,
            "independent enthalpy mismatch");
        near(
            initial_current
                .caloric
                .specific_internal_energy_j_per_kg[
                    phase],
            ref::
                specific_internal_energy_j_per_kg[
                    phase],
            4.0e-10,
            2.0e-7,
            "independent internal-energy mismatch");
    }

    near(
        initial_current
            .rock
            .volumetric_internal_energy_j_per_rock_m3,
        0.0,
        0.0,
        1.0e-12,
        "quartz rock-energy reference changed");
    near(
        initial_current
            .rock
            .volumetric_internal_energy_gradient[
                layout
                    .temperature_unknown_index()],
        ref::
            quartz_volumetric_heat_capacity_j_per_m3_k,
        3.0e-12,
        3.0e-6,
        "NIST quartz heat-capacity interpolation mismatch");

    const auto previous_component =
        flow::
            build_pore_volume_component_accumulation(
                initial_current.state,
                ref::porosity);
    const auto previous_energy =
        flow::
            build_pore_volume_energy_accumulation_snapshot(
                initial_current.state,
                ref::porosity,
                initial_current.transport,
                initial_current.caloric,
                initial_current.rock);

    for (std::size_t component = 0U;
         component < 6U;
         ++component) {
        near(
            previous_component
                .component_accumulation_mol_per_bulk_m3[
                    component],
            ref::
                component_accumulation_mol_per_bulk_m3[
                    component],
            5.0e-10,
            2.0e-7,
            "independent component-accumulation mismatch");
    }
    near(
        previous_energy
            .total_internal_energy_j_per_bulk_m3,
        ref::
            total_internal_energy_j_per_bulk_m3,
        5.0e-10,
        5.0e-3,
        "independent energy-accumulation mismatch");

    const auto partition =
        make_partition();
    std::optional<
        fdp::
            VariableCardinalityNaturalVariableNumbering3D>
        numbering;
    require(
        fdp::
            make_variable_cardinality_natural_variable_numbering_3d(
                PETSC_COMM_SELF,
                partition,
                6U,
                std::vector<std::size_t>{3U},
                &numbering) ==
                PETSC_SUCCESS &&
            numbering.has_value() &&
            numbering
                    ->petsc_global_scalar_count() ==
                19,
        "failed to build serial six-component 3P numbering");

    const auto connection =
        well::make_peaceman_well_index_3d(
            {
                ref::dx_m,
                ref::dy_m,
                ref::dz_m},
            {
                ref::permeability_m2,
                ref::permeability_m2,
                ref::permeability_m2},
            well::
                AxisAlignedWellDirection3D::z,
            ref::wellbore_radius_m,
            0.0);

    auto well_context =
        wdp::
            FixedBhpPeacemanWellSourceEvaluatorContext3D::
                create(
                    mesh::GlobalEntityId{
                        UINT64_C(10)},
                    connection,
                    ref::
                        bottom_hole_pressure_pa,
                    wd::
                        FixedBhpInjectionEnthalpy3D{
                            "Li-Firoozabadi stationary benchmark injection fallback; unused because q_phase=0",
                            {
                                ref::
                                    specific_enthalpy_j_per_kg
                                        .begin(),
                                ref::
                                    specific_enthalpy_j_per_kg
                                        .end()}},
                    "Li-Firoozabadi-2012 stationary fixed-BHP Peaceman benchmark");

    const auto initial_well =
        wdp::
            build_fixed_bhp_peaceman_well_source_3d(
                well_context,
                fdp::
                    MixedCardinalityPhysicalCurrentCellLinearization3D{
                        initial_current});
    require(
        initial_well
                .phase_volumetric_rate_m3_per_s
                .size() ==
            3U &&
        initial_well
                .cell_source
                .component_molar_rate_mol_per_s
                .size() ==
            6U,
        "stationary fixed-BHP well source shape mismatch");
    for (double rate :
         initial_well
             .phase_volumetric_rate_m3_per_s) {
        near(
            rate,
            0.0,
            0.0,
            2.0e-18,
            "independent stationary Peaceman phase rate is not zero");
    }
    for (double source :
         initial_well
             .cell_source
             .component_molar_rate_mol_per_s) {
        near(
            source,
            0.0,
            0.0,
            2.0e-12,
            "independent stationary component source is not zero");
    }
    near(
        initial_well
            .cell_source
            .energy_rate_w,
        0.0,
        0.0,
        2.0e-6,
        "independent stationary energy source is not zero");

    fdp::FixedThreePhaseSnesCellInput3D
        cell_input;
    cell_input.cell =
        mesh::LocalIndex{0U};
    cell_input.cell_global =
        mesh::GlobalEntityId{
            UINT64_C(10)};
    cell_input.bulk_volume_m3 =
        ref::bulk_volume_m3;
    cell_input.porosity =
        ref::porosity;
    cell_input.frozen_layout =
        layout;
    cell_input.component_ids =
        component_ids;
    cell_input.previous_component_accumulation =
        previous_component;
    cell_input.previous_energy_accumulation =
        previous_energy;

    auto schedule =
        make_schedule();
    auto bridge =
        make_cell_bridge();
    auto pattern =
        make_cell_pattern();
    std::optional<
        fdp::
            MixedCardinalityPhysicalSnesAssemblyContext3D>
        assembly;
    require(
        fdp::
            MixedCardinalityPhysicalSnesAssemblyContext3D::
                create(
                    PETSC_COMM_SELF,
                    schedule,
                    partition,
                    *numbering,
                    bridge,
                    pattern,
                    ref::timestep_seconds,
                    {
                        fdp::
                            MixedCardinalityPhysicalSnesCellInput3D{
                                std::move(
                                    cell_input)}
                    },
                    {
                        make_phase_identities()
                    },
                    {},
                    {
                        {
                            &unsupported_single_phase,
                            nullptr},
                        {
                            &unsupported_two_phase,
                            nullptr},
                        {
                            &fdp::
                                evaluate_pr76_three_phase_production_cell_3d<
                                    Closure>,
                            &evaluator_context}},
                    {},
                    wdp::
                        fixed_bhp_peaceman_well_source_binding_3d(
                            &well_context),
                    &assembly) ==
                PETSC_SUCCESS &&
            assembly.has_value(),
        "failed to create source-complete serial sour-gas assembly");

    Mat jacobian = nullptr;
    require(
        assembly
                ->create_jacobian_structure(
                    &jacobian) ==
                PETSC_SUCCESS &&
            jacobian != nullptr,
        "failed to create serial sour-gas Jacobian");

    Vec initial = nullptr;
    require(
        fdp::
            create_variable_cardinality_natural_variable_vec_3d(
                PETSC_COMM_SELF,
                *numbering,
                &initial) ==
                PETSC_SUCCESS &&
            initial != nullptr,
        "failed to create serial sour-gas state Vec");
    insert_state(
        initial,
        *numbering,
        q);

    Vec solution = nullptr;
    std::optional<
        fdp::
            VariableCardinalityNaturalVariableSnesSolveReport3D>
        report;
    require(
        fdp::
            solve_variable_cardinality_natural_variable_snes_3d(
                PETSC_COMM_SELF,
                *numbering,
                initial,
                jacobian,
                assembly->snes_evaluator(),
                &solution,
                &report) ==
                PETSC_SUCCESS &&
            solution != nullptr &&
            report.has_value() &&
            static_cast<int>(
                report
                    ->converged_reason) >
                0 &&
            report->snes_type ==
                SNESNEWTONLS &&
            report->ksp_type ==
                KSPGMRES &&
            report->pc_type ==
                PCASM &&
            report
                    ->final_function_l2_norm <=
                2.0e-8,
        "source-complete serial sour-gas production SNES did not converge");

    const auto solved =
        read_state(
            solution,
            *numbering);
    require(
        solved.size() ==
            q.size(),
        "serial sour-gas solved-state width changed");
    for (std::size_t slot = 0U;
         slot < q.size();
         ++slot) {
        near(
            solved[slot],
            q[slot],
            3.0e-7,
            slot < 2U
                ? 2.0e-5
                : 2.0e-7,
            "serial sour-gas state left independent stationary reference");
    }

    auto final_current =
        evaluate_current(
            solved,
            layout,
            component_ids,
            &evaluator_context);
    const auto final_component =
        flow::
            build_pore_volume_component_accumulation(
                final_current.state,
                ref::porosity);
    const auto final_energy =
        flow::
            build_pore_volume_energy_accumulation_snapshot(
                final_current.state,
                ref::porosity,
                final_current.transport,
                final_current.caloric,
                final_current.rock);

    for (std::size_t component = 0U;
         component < 6U;
         ++component) {
        near(
            final_component
                .component_accumulation_mol_per_bulk_m3[
                    component] *
                ref::bulk_volume_m3,
            ref::
                cell_component_inventory_mol[
                    component],
            8.0e-7,
            5.0e-2,
            "serial sour-gas component inventory is not conserved");
    }
    near(
        final_energy
                .total_internal_energy_j_per_bulk_m3 *
            ref::bulk_volume_m3,
        ref::
            cell_total_internal_energy_j,
        8.0e-7,
        2.0e7,
        "serial sour-gas energy inventory is not conserved");

    const auto final_well =
        wdp::
            build_fixed_bhp_peaceman_well_source_3d(
                well_context,
                fdp::
                    MixedCardinalityPhysicalCurrentCellLinearization3D{
                        final_current});
    for (double rate :
         final_well
             .phase_volumetric_rate_m3_per_s) {
        near(
            rate,
            0.0,
            0.0,
            2.0e-10,
            "accepted stationary Peaceman phase rate drifted from zero");
    }
    for (double source :
         final_well
             .cell_source
             .component_molar_rate_mol_per_s) {
        near(
            source,
            0.0,
            0.0,
            2.0e-5,
            "accepted stationary Peaceman component rate drifted from zero");
    }
    near(
        final_well
            .cell_source
            .energy_rate_w,
        0.0,
        0.0,
        20.0,
        "accepted stationary Peaceman energy rate drifted from zero");

    Vec residual = nullptr;
    require(
        VecDuplicate(
            solution,
            &residual) ==
                PETSC_SUCCESS &&
            VecSet(
                residual,
                PetscScalar{0.0}) ==
                PETSC_SUCCESS,
        "failed to allocate final sour-gas residual");

    auto evaluator =
        assembly->snes_evaluator();
    fdp::NaturalVariableSnesEvaluationStatus3D
        final_status =
            fdp::
                NaturalVariableSnesEvaluationStatus3D::
                    success;
    require(
        evaluator.function(
            solution,
            residual,
            evaluator.user_context,
            &final_status) ==
                PETSC_SUCCESS &&
            final_status ==
                fdp::
                    NaturalVariableSnesEvaluationStatus3D::
                        success &&
            VecAssemblyBegin(
                residual) ==
                PETSC_SUCCESS &&
            VecAssemblyEnd(
                residual) ==
                PETSC_SUCCESS,
        "final serial sour-gas residual evaluation failed");

    PetscReal final_norm = 0.0;
    require(
        VecNorm(
            residual,
            NORM_2,
            &final_norm) ==
                PETSC_SUCCESS &&
            static_cast<double>(
                final_norm) <=
                2.0e-8,
        "final serial sour-gas residual is not closed");

    require(
        VecDestroy(
            &residual) ==
                PETSC_SUCCESS &&
            VecDestroy(
                &solution) ==
                PETSC_SUCCESS &&
            VecDestroy(
                &initial) ==
                PETSC_SUCCESS &&
            MatDestroy(
                &jacobian) ==
                PETSC_SUCCESS,
        "serial sour-gas PETSc cleanup failed");
}
