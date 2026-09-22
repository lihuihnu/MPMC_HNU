#include <mpmc/flow/pr76_methane_ethane_propane_properties.hpp>
#include <mpmc/flash/pr76_pt_flash_backend.hpp>
#include <mpmc/flow_discretization/single_phase_tpfa.hpp>
#include <mpmc/flow_discretization_petsc/complete_natural_variable_petsc_materialization.hpp>
#include <mpmc/flow_discretization_petsc/adaptive_timestep_controller.hpp>
#include <mpmc/flow_discretization_petsc/pr76_production_cell_evaluator.hpp>
#include <mpmc/flow_discretization_petsc/post_snes_pt_flash_phase_transition_scanner.hpp>
#include <mpmc/flow_discretization_petsc/single_phase_adaptive_timestep_attempt.hpp>

#include <petscsys.h>

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
#include <utility>
#include <vector>

namespace {

namespace disc = mpmc::discretization;
namespace dp = mpmc::discretization_petsc;
namespace fd = mpmc::flow_discretization;
namespace fdp = mpmc::flow_discretization_petsc;
namespace fl = mpmc::flash;
namespace flow = mpmc::flow;
namespace mesh = mpmc::mesh;
namespace th = mpmc::thermodynamics;

void require_real_collective(
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
            PETSC_COMM_WORLD) != MPI_SUCCESS) {
        throw std::runtime_error(
            "MPI_Allreduce failed in real PR76 transient regression");
    }
    if (global == 0) {
        throw std::runtime_error(
            std::string{message});
    }
}

void near_real_collective(
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
    require_real_collective(
        std::isfinite(actual) &&
            std::isfinite(expected) &&
            std::abs(actual - expected) <=
                absolute + relative * scale,
        message);
}

th::Provenance real_fixture_source(
    std::string locator) {
    return {
        th::SourceKind::synthetic_test,
        "test://pr76-production-transient",
        "v1",
        std::move(locator),
        "Software fixture reproducing the repository-curated CH4/C2H6/C3H8 PR76 numeric identity",
        "Constructed by the PETSc production-bridge regression",
        "No third-party source text redistributed"};
}

th::SourcedScalar real_fixture_datum(
    double value,
    th::Unit unit,
    std::string locator) {
    return {
        value,
        unit,
        real_fixture_source(
            std::move(locator)),
        "specified SI",
        "identity"};
}

th::PrParameterSet real_pr76_parameters() {
    constexpr std::array<std::string_view, 3>
        ids{"methane", "ethane", "propane"};
    constexpr std::array<double, 3>
        tc{190.555, 305.4, 369.825};
    constexpr std::array<double, 3>
        pc{4.595e6, 4.88e6, 4.248e6};
    constexpr std::array<double, 3>
        omega{0.0, 0.099, 0.15308};
    constexpr std::array<double, 3>
        mw{0.0160425, 0.0300690, 0.0440956};

    std::vector<th::Component> catalog;
    th::PrParameterInput input;
    input.model_id =
        std::string{th::pr76_profile};
    input.dataset_id =
        "DeitersBell-aic16730-PengRobinson1976-ternary";
    input.revision =
        "production-transient-regression-v1";
    input.applicability = {
        std::nullopt,
        std::nullopt,
        real_fixture_source(
            "unknown applicability")};

    for (std::size_t i = 0U;
         i < ids.size();
         ++i) {
        const std::string id{ids[i]};
        catalog.push_back({
            id,
            id,
            th::ComponentKind::pure,
            real_fixture_source(
                "component-" +
                std::to_string(i)),
            real_fixture_datum(
                mw[i],
                th::Unit::kilogram_per_mole,
                "molar-mass-" +
                    std::to_string(i))});
        input.pure.push_back({
            id,
            real_fixture_datum(
                tc[i],
                th::Unit::kelvin,
                "Tc-" +
                    std::to_string(i)),
            real_fixture_datum(
                pc[i],
                th::Unit::pascal,
                "Pc-" +
                    std::to_string(i)),
            real_fixture_datum(
                omega[i],
                th::Unit::dimensionless,
                "omega-" +
                    std::to_string(i))});
    }

    const auto zero =
        real_fixture_datum(
            0.0,
            th::Unit::dimensionless,
            "zero-kij");
    input.binary.push_back(
        {"methane", "ethane", zero});
    input.binary.push_back(
        {"methane", "propane", zero});
    input.binary.push_back(
        {"ethane", "propane", zero});

    const std::array<std::string, 3> order{
        "methane",
        "ethane",
        "propane"};
    return th::PrParameterSet::create(
        catalog,
        order,
        input,
        th::DataPolicy::allow_synthetic_tests);
}

const std::vector<std::string>& real_component_ids() {
    static const std::vector<std::string> ids{
        "methane",
        "ethane",
        "propane"};
    return ids;
}

flow::NaturalVariableLayout1P real_layout_1p() {
    return flow::NaturalVariableLayout1P{
        flow::NaturalVariableCompositionPivot1P::
            from_dependent_component(
                3U,
                1U)};
}

std::vector<double>
real_previous_state(std::uint64_t stable) {
    const auto layout =
        real_layout_1p();
    std::vector<double> q(
        layout.unknown_count(),
        0.0);
    q[layout.pressure_unknown_index()] =
        stable == UINT64_C(10)
            ? 8.0e6
            : 7.99e6;
    q[layout.temperature_unknown_index()] =
        stable == UINT64_C(10)
            ? 450.0
            : 449.0;
    q[*layout.independent_composition_unknown_index(
        0U)] =
        stable == UINT64_C(10)
            ? 0.50
            : 0.49;
    q[*layout.independent_composition_unknown_index(
        2U)] =
        0.20;
    return q;
}

PetscErrorCode disabled_rock_storage(
    const flow::NaturalVariableLayoutDescriptor& layout,
    std::span<const double>,
    void*,
    std::optional<
        fdp::Pr76RockThermalStorageLinearization3D>* output,
    fdp::NaturalVariableSnesEvaluationStatus3D*
        status) {
    if (output == nullptr || status == nullptr) {
        return PETSC_ERR_ARG_NULL;
    }
    output->emplace(
        fdp::Pr76RockThermalStorageLinearization3D{
            0.0,
            std::vector<double>(
                layout.unknown_count(),
                0.0),
            {
                "disabled-stationary-rock-storage",
                "pr76-production-transient",
                "v1"}});
    *status =
        fdp::NaturalVariableSnesEvaluationStatus3D::
            success;
    return PETSC_SUCCESS;
}

PetscErrorCode linear_two_phase_kr(
    const flow::NaturalVariableCellState2P& state,
    void*,
    std::optional<
        fdp::Pr76TwoPhaseRelativePermeabilityLinearization3D>* output,
    fdp::NaturalVariableSnesEvaluationStatus3D*
        status) {
    if (output == nullptr || status == nullptr) {
        return PETSC_ERR_ARG_NULL;
    }
    const std::size_t q =
        state.layout().unknown_count();
    std::array<std::vector<double>, 2>
        gradient{
            std::vector<double>(q, 0.0),
            std::vector<double>(q, 0.0)};
    const std::size_t s =
        state.layout()
            .independent_saturation_unknown_index();
    gradient[0][s] = 1.0;
    gradient[1][s] = -1.0;
    output->emplace(
        fdp::Pr76TwoPhaseRelativePermeabilityLinearization3D{
            {
                state.phase_saturation(0U),
                state.phase_saturation(1U)},
            std::move(gradient)});
    *status =
        fdp::NaturalVariableSnesEvaluationStatus3D::
            success;
    return PETSC_SUCCESS;
}

struct SaturationKr3P {
    template <typename Number>
    flow::RelativePermeabilityEvaluation3P<Number>
    operator()(
        const flow::ThreePhaseSaturationState3P<Number>&
            state) const {
        return {state.saturation};
    }
};

PetscErrorCode linear_three_phase_constitutive(
    const flow::NaturalVariableCellState3P& state,
    void*,
    std::optional<
        flow::
            ThreePhaseSaturationConstitutiveNaturalVariableLinearization3P>*
                output,
    fdp::NaturalVariableSnesEvaluationStatus3D*
        status) {
    if (output == nullptr || status == nullptr) {
        return PETSC_ERR_ARG_NULL;
    }
    try {
        const auto primal =
            flow::evaluate_three_phase_saturation_constitutive(
                state,
                SaturationKr3P{},
                flow::NoCapillaryPressure3P{});
        flow::
            ThreePhaseSaturationCoordinateDerivatives3P
            derivatives{};
        derivatives.relative_permeability[0] =
            {1.0, 0.0};
        derivatives.relative_permeability[1] =
            {0.0, 1.0};
        derivatives.relative_permeability[2] =
            {-1.0, -1.0};
        output->emplace(
            flow::
                make_saturation_constitutive_natural_variable_linearization(
                    state,
                    primal,
                    derivatives));
        *status =
            fdp::NaturalVariableSnesEvaluationStatus3D::
                success;
        return PETSC_SUCCESS;
    } catch (const std::exception&) {
        *status =
            fdp::NaturalVariableSnesEvaluationStatus3D::
                domain_error;
        return PETSC_SUCCESS;
    }
}

mesh::Topology real_topology(int rank) {
    mesh::Topology::EntityIds ids;
    ids.faces = {
        mesh::GlobalEntityId{UINT64_C(100)}};
    ids.cells =
        rank == 0
            ? std::vector<mesh::GlobalEntityId>{
                  mesh::GlobalEntityId{UINT64_C(10)},
                  mesh::GlobalEntityId{UINT64_C(20)}}
            : std::vector<mesh::GlobalEntityId>{
                  mesh::GlobalEntityId{UINT64_C(20)},
                  mesh::GlobalEntityId{UINT64_C(10)}};
    return {std::move(ids), {}};
}

mesh::PartitionSnapshot
real_partition(int rank) {
    auto topology =
        real_topology(rank);
    mesh::EntityOwnerRanks owners;
    owners.faces = {
        mesh::PartitionRank{0U}};
    owners.cells =
        rank == 0
            ? std::vector<mesh::PartitionRank>{
                  mesh::PartitionRank{0U},
                  mesh::PartitionRank{1U}}
            : std::vector<mesh::PartitionRank>{
                  mesh::PartitionRank{1U},
                  mesh::PartitionRank{0U}};
    return mesh::PartitionSnapshot::create(
        topology,
        mesh::PartitionRank{
            static_cast<
                mesh::PartitionRank::value_type>(
                    rank)},
        2U,
        std::move(owners));
}

dp::ParallelOwnedConnectionSchedule3D
real_schedule(int rank) {
    constexpr double transmissibility =
        1.0e-16;
    const dp::AssemblyReadyInternalConnectionRow3D
        row{
            mesh::LocalIndex{0U},
            mesh::GlobalEntityId{UINT64_C(100)},
            rank == 0
                ? mesh::LocalIndex{0U}
                : mesh::LocalIndex{1U},
            mesh::GlobalEntityId{UINT64_C(10)},
            rank == 0
                ? mesh::LocalIndex{1U}
                : mesh::LocalIndex{0U},
            mesh::GlobalEntityId{UINT64_C(20)},
            transmissibility};
    if (rank == 0) {
        return {
            mesh::PartitionRank{0U},
            2U,
            {row},
            {}};
    }
    return {
        mesh::PartitionRank{1U},
        2U,
        {},
        {row}};
}

mesh::DofLayout real_dof_layout(int rank) {
    return mesh::DofLayout::create(
        real_topology(rank),
        {
            mesh::DofVariable{
                "real_pr76_natural_state_1p",
                mesh::EntityKind::cell,
                4U}
        });
}

mesh::DofNumberingSnapshot
real_numbering(
    int rank,
    const mesh::DofLayout& layout,
    const mesh::PartitionSnapshot& partition) {
    mesh::GlobalEntityNumberingInput input;
    input.global_cell_count = 2U;
    input.global_face_count = 1U;
    input.faces = {{
        mesh::GlobalEntityId{UINT64_C(100)},
        mesh::GlobalEntityOrdinal{UINT64_C(0)}}};
    input.cells =
        rank == 0
            ? std::vector<
                  mesh::GlobalEntityOrdinalRecord>{
                  {
                      mesh::GlobalEntityId{UINT64_C(10)},
                      mesh::GlobalEntityOrdinal{UINT64_C(1)}},
                  {
                      mesh::GlobalEntityId{UINT64_C(20)},
                      mesh::GlobalEntityOrdinal{UINT64_C(0)}}}
            : std::vector<
                  mesh::GlobalEntityOrdinalRecord>{
                  {
                      mesh::GlobalEntityId{UINT64_C(20)},
                      mesh::GlobalEntityOrdinal{UINT64_C(0)}},
                  {
                      mesh::GlobalEntityId{UINT64_C(10)},
                      mesh::GlobalEntityOrdinal{UINT64_C(1)}}};
    return mesh::DofNumberingSnapshot::create_local(
        layout,
        partition,
        std::move(input));
}

dp::PetscMpiAijSymbolicPreallocation3D
real_cell_bridge(int rank) {
    const PetscInt row =
        static_cast<PetscInt>(rank);
    return {
        mesh::PartitionRank{
            static_cast<
                mesh::PartitionRank::value_type>(
                    rank)},
        2U,
        row,
        row + 1,
        2,
        {mesh::LocalIndex{0U}},
        {
            mesh::GlobalEntityId{
                rank == 0
                    ? UINT64_C(10)
                    : UINT64_C(20)}
        },
        {row},
        {1},
        {1},
        rank == 0
            ? std::vector<PetscInt>{0, 1}
            : std::vector<PetscInt>{1, 0}};
}

dp::OwnedCellStructuralColumnPatternSnapshot3D
real_cell_pattern(int rank) {
    const PetscInt row =
        static_cast<PetscInt>(rank);
    return {
        mesh::PartitionRank{
            static_cast<
                mesh::PartitionRank::value_type>(
                    rank)},
        2U,
        2U,
        row,
        row + 1,
        2,
        {
            dp::OwnedCellStructuralColumnPatternRow3D{
                mesh::LocalIndex{0U},
                mesh::GlobalEntityId{
                    rank == 0
                        ? UINT64_C(10)
                        : UINT64_C(20)},
                row,
                0U,
                1U,
                0U,
                1U}
        },
        {row},
        {
            static_cast<PetscInt>(
                rank == 0 ? 1 : 0)}
    };
}

disc::CombinedTransmissibilityAdmissibility3D
real_admissibility() {
    return {
        disc::
            CombinedTransmissibilityAdmissibilityDisposition3D::
                direct_normal_projection_k_orthogonal_candidate,
        {
            disc::TransmissibilityGeometryDisposition3D::
                direct_normal_projection_allowed,
            0.0,
            0.05},
        {
            disc::KOrthogonalityDisposition3D::
                k_orthogonal_within_policy,
            std::nullopt,
            std::nullopt,
            {
                std::nullopt,
                std::nullopt},
            0.05}};
}

disc::TpfaInternalFaceTransmissibilityEntry3D
real_transmissibility() {
    return {
        mesh::LocalIndex{0U},
        disc::TpfaInternalFaceTransmissibilityDisposition3D::
            materialized,
        real_admissibility(),
        disc::TpfaStaticFaceTransmissibility3D{
            disc::TpfaStaticFaceTransmissibilityDisposition3D::
                positive_harmonic_combination,
            1.5,
            1.0e-16}};
}

template <typename Closure>
std::optional<
    fdp::SinglePhaseCurrentCellLinearization3D>
real_direct_cell_typed(
    std::uint64_t stable,
    std::span<const double> q,
    fdp::Pr76SinglePhaseProductionCellEvaluatorContext3D<
        Closure>* context) {
    std::optional<
        fdp::SinglePhaseCurrentCellLinearization3D>
        output;
    auto layout =
        real_layout_1p();
    fdp::NaturalVariableSnesEvaluationStatus3D
        status =
            fdp::NaturalVariableSnesEvaluationStatus3D::
                success;
    const PetscErrorCode error =
        fdp::
            evaluate_pr76_single_phase_production_cell_3d<
                Closure>(
                mesh::LocalIndex{0U},
                mesh::GlobalEntityId{stable},
                q,
                layout,
                real_component_ids(),
                context,
                &output,
                &status);
    if (error != PETSC_SUCCESS ||
        status !=
            fdp::NaturalVariableSnesEvaluationStatus3D::
                success ||
        !output) {
        throw std::runtime_error(
            "real PR76 direct 1P production evaluation failed");
    }
    return output;
}

template <typename Closure>
std::vector<fdp::SinglePhaseSnesCellInput3D>
real_cells(
    int rank,
    fdp::Pr76SinglePhaseProductionCellEvaluatorContext3D<
        Closure>* context) {
    std::vector<fdp::SinglePhaseSnesCellInput3D>
        cells;
    cells.reserve(2U);
    for (std::size_t local = 0U;
         local < 2U;
         ++local) {
        const std::uint64_t stable =
            rank == 0
                ? (local == 0U
                       ? UINT64_C(10)
                       : UINT64_C(20))
                : (local == 0U
                       ? UINT64_C(20)
                       : UINT64_C(10));
        const bool owned =
            (rank == 0 &&
             stable == UINT64_C(10)) ||
            (rank == 1 &&
             stable == UINT64_C(20));

        const auto previous =
            real_previous_state(stable);
        auto evaluated =
            real_direct_cell_typed(
                stable,
                previous,
                context);

        fdp::SinglePhaseSnesCellInput3D cell;
        cell.cell =
            mesh::LocalIndex{
                static_cast<
                    mesh::LocalIndex::value_type>(
                        local)};
        cell.cell_global =
            mesh::GlobalEntityId{stable};
        cell.bulk_volume_m3 =
            stable == UINT64_C(10)
                ? 2.0
                : 5.0;
        cell.porosity =
            stable == UINT64_C(10)
                ? 0.25
                : 0.30;
        cell.frozen_layout =
            real_layout_1p();
        cell.component_ids =
            real_component_ids();

        if (owned) {
            cell.previous_component_accumulation =
                flow::
                    build_single_phase_component_accumulation(
                        evaluated->state,
                        cell.porosity);
            cell.previous_energy_accumulation =
                flow::
                    build_single_phase_energy_accumulation_snapshot(
                        evaluated->state,
                        cell.porosity,
                        evaluated->transport,
                        evaluated->caloric,
                        evaluated->rock);
        }
        cells.push_back(
            std::move(cell));
    }
    return cells;
}

std::vector<
    fdp::SinglePhaseSnesAuthoritativeFaceInput3D>
real_faces(int rank) {
    if (rank != 0) {
        return {};
    }
    return {{
        mesh::LocalIndex{0U},
        mesh::GlobalEntityId{UINT64_C(100)},
        real_transmissibility(),
        flow::GravityVector3D{
            0.0, 0.0, -9.81},
        flow::OwnerToNeighbourDisplacement3D{
            0.0, 0.0, 1.0},
        fd::StaticThermalFaceConductance3D{
            0.05}}};
}

const fdp::CompleteNaturalVariableResidualEntry3D&
find_real_residual_entry(
    const fdp::CompleteNaturalVariableAssemblySnapshot3D& snapshot,
    mesh::GlobalEntityId row_cell,
    std::size_t equation_slot,
    fdp::NaturalVariableEquationKind3D equation_kind) {
    const auto it =
        std::find_if(
            snapshot.residual_entries().begin(),
            snapshot.residual_entries().end(),
            [&](const auto& entry) {
                return entry.row_cell_global == row_cell &&
                    entry.equation_slot == equation_slot &&
                    entry.equation_kind == equation_kind;
            });
    if (it == snapshot.residual_entries().end()) {
        throw std::runtime_error(
            "real PR76 residual entry not found");
    }
    return *it;
}

const fdp::CompleteNaturalVariableJacobianEntry3D&
find_real_jacobian_entry(
    const fdp::CompleteNaturalVariableAssemblySnapshot3D& snapshot,
    mesh::GlobalEntityId row_cell,
    std::size_t equation_slot,
    fdp::NaturalVariableEquationKind3D equation_kind,
    mesh::GlobalEntityId column_cell,
    std::size_t natural_variable_column) {
    const auto it =
        std::find_if(
            snapshot.jacobian_entries().begin(),
            snapshot.jacobian_entries().end(),
            [&](const auto& entry) {
                return entry.row_cell_global == row_cell &&
                    entry.equation_slot == equation_slot &&
                    entry.equation_kind == equation_kind &&
                    entry.column_cell_global == column_cell &&
                    entry.natural_variable_column ==
                        natural_variable_column;
            });
    if (it == snapshot.jacobian_entries().end()) {
        throw std::runtime_error(
            "real PR76 Jacobian entry not found");
    }
    return *it;
}

void set_real_previous_state(
    Vec state,
    int rank) {
    const auto values =
        real_previous_state(
            rank == 0
                ? UINT64_C(10)
                : UINT64_C(20));
    const PetscInt start =
        static_cast<PetscInt>(
            rank * 4);
    for (std::size_t slot = 0U;
         slot < values.size();
         ++slot) {
        const PetscInt index =
            start +
            static_cast<PetscInt>(
                slot);
        const PetscScalar value =
            static_cast<PetscScalar>(
                values[slot]);
        require_real_collective(
            VecSetValues(
                state,
                1,
                &index,
                &value,
                INSERT_VALUES) ==
                PETSC_SUCCESS,
            "failed to seed real PR76 previous state");
    }
    require_real_collective(
        VecAssemblyBegin(state) ==
                PETSC_SUCCESS &&
            VecAssemblyEnd(state) ==
                PETSC_SUCCESS,
        "failed to assemble real PR76 previous state");
}

template <typename Closure>
void check_real_face_flux(
    fdp::Pr76SinglePhaseProductionCellEvaluatorContext3D<
        Closure>* context) {
    const auto q10 =
        real_previous_state(UINT64_C(10));
    const auto q20 =
        real_previous_state(UINT64_C(20));
    auto owner =
        real_direct_cell_typed(
            UINT64_C(10),
            q10,
            context);
    auto neighbour =
        real_direct_cell_typed(
            UINT64_C(20),
            q20,
            context);

    const auto owner_mobility =
        flow::build_single_phase_mobility_linearization(
            owner->state,
            owner->transport);
    const auto neighbour_mobility =
        flow::build_single_phase_mobility_linearization(
            neighbour->state,
            neighbour->transport);
    const auto potential =
        flow::build_single_phase_potential_upwind_linearization(
            owner_mobility,
            neighbour_mobility,
            {0.0, 0.0, -9.81},
            {0.0, 0.0, 1.0});

    const auto advective_only =
        fd::build_single_phase_tpfa_face_linearization(
            real_transmissibility(),
            potential,
            owner->state,
            owner->molar_density,
            owner->transport,
            owner->caloric,
            neighbour->state,
            neighbour->molar_density,
            neighbour->transport,
            neighbour->caloric,
            {2.0, 5.0},
            {0.0});
    const auto full =
        fd::build_single_phase_tpfa_face_linearization(
            real_transmissibility(),
            potential,
            owner->state,
            owner->molar_density,
            owner->transport,
            owner->caloric,
            neighbour->state,
            neighbour->molar_density,
            neighbour->transport,
            neighbour->caloric,
            {2.0, 5.0},
            {0.05});

    require_real_collective(
        std::abs(full.volumetric_flux_m3_per_s) >
            0.0,
        "real PR76 TPFA Darcy flux is zero");

    double component_norm = 0.0;
    for (const double value :
         full.component
             .owner_component_contribution_mol_per_bulk_m3_s) {
        component_norm +=
            std::abs(value);
    }
    require_real_collective(
        component_norm > 0.0,
        "real PR76 component molar face flux is zero");

    const double advective_rate =
        advective_only.energy
            .owner_contribution_w_per_bulk_m3 *
        2.0;
    const double full_rate =
        full.energy
            .owner_contribution_w_per_bulk_m3 *
        2.0;
    require_real_collective(
        std::abs(advective_rate) > 0.0,
        "real PR76 advective enthalpy flux is zero");
    require_real_collective(
        std::abs(full_rate - advective_rate) >
            0.0,
        "real PR76 conductive heat flux is zero");

    for (std::size_t component = 0U;
         component < 3U;
         ++component) {
        const double balance =
            full.component
                .owner_component_contribution_mol_per_bulk_m3_s[
                    component] *
                2.0 +
            full.component
                .neighbour_component_contribution_mol_per_bulk_m3_s[
                    component] *
                5.0;
        near_real_collective(
            balance,
            0.0,
            0.0,
            1.0e-12,
            "real PR76 face component rate is not conservative");
    }

    const double energy_balance =
        full.energy
            .owner_contribution_w_per_bulk_m3 *
            2.0 +
        full.energy
            .neighbour_contribution_w_per_bulk_m3 *
            5.0;
    near_real_collective(
        energy_balance,
        0.0,
        0.0,
        1.0e-10,
        "real PR76 face energy rate is not conservative");
}

template <typename Closure1>
void check_real_cardinality_bridges(
    const th::Pr76Phase<double>& model,
    Closure1* one_phase_closure) {
    fdp::Pr76SinglePhaseProductionCellEvaluatorContext3D<
        Closure1>
        one_context{
            one_phase_closure,
            1.0,
            {&disabled_rock_storage, nullptr},
            {}};
    auto one =
        real_direct_cell_typed(
            UINT64_C(10),
            real_previous_state(UINT64_C(10)),
            &one_context);
    require_real_collective(
        one->transport.dynamic_viscosity_pa_s > 0.0 &&
            std::isfinite(
                one->caloric
                    .specific_enthalpy_j_per_kg),
        "real PR76 1P production bridge did not publish real properties");

    auto two_closure =
        flow::
            make_pr76_methane_ethane_propane_property_closure(
                model,
                std::vector<
                    th::Pr76SelectedPhase>{
                    {0U, {}},
                    {0U, {}}});
    using Closure2 =
        decltype(two_closure);
    fdp::Pr76TwoPhaseProductionCellEvaluatorContext3D<
        Closure2>
        two_context{
            &two_closure,
            {&linear_two_phase_kr, nullptr},
            {&disabled_rock_storage, nullptr},
            {}};
    const auto two_layout =
        flow::NaturalVariableLayout2P{
            flow::NaturalVariableCompositionPivot2P::
                from_dependent_components(
                    3U,
                    {1U, 1U})};
    std::vector<double> q2(
        two_layout.unknown_count(),
        0.0);
    q2[0] = 7.5e6;
    q2[1] = 450.0;
    q2[2] = 0.55;
    q2[3] = 0.50;
    q2[4] = 0.20;
    q2[5] = 0.40;
    q2[6] = 0.20;
    std::optional<
        fdp::TwoPhaseCurrentCellLinearization3D>
        two;
    fdp::NaturalVariableSnesEvaluationStatus3D
        two_status =
            fdp::NaturalVariableSnesEvaluationStatus3D::
                success;
    require_real_collective(
        fdp::
            evaluate_pr76_two_phase_production_cell_3d<
                Closure2>(
                mesh::LocalIndex{0U},
                mesh::GlobalEntityId{UINT64_C(10)},
                q2,
                two_layout,
                real_component_ids(),
                &two_context,
                &two,
                &two_status) ==
                PETSC_SUCCESS &&
            two_status ==
                fdp::NaturalVariableSnesEvaluationStatus3D::
                    success &&
            two.has_value() &&
            two->fugacity.residual_count() == 3U &&
            two->transport.dynamic_viscosity_pa_s[0] > 0.0 &&
            two->transport.dynamic_viscosity_pa_s[1] > 0.0,
        "real PR76 2P production bridge failed");

    auto three_closure =
        flow::
            make_pr76_methane_ethane_propane_property_closure(
                model,
                std::vector<
                    th::Pr76SelectedPhase>{
                    {0U, {}},
                    {0U, {}},
                    {0U, {}}});
    using Closure3 =
        decltype(three_closure);
    fdp::Pr76ThreePhaseProductionCellEvaluatorContext3D<
        Closure3>
        three_context{
            &three_closure,
            {&linear_three_phase_constitutive, nullptr},
            {&disabled_rock_storage, nullptr},
            {}};
    const auto three_layout =
        flow::NaturalVariableLayout3P{
            flow::NaturalVariableCompositionPivot3P::
                from_dependent_components(
                    3U,
                    {1U, 1U, 1U})};
    std::vector<double> q3(
        three_layout.unknown_count(),
        0.0);
    q3[0] = 7.5e6;
    q3[1] = 450.0;
    q3[2] = 0.30;
    q3[3] = 0.35;
    const std::array<
        std::array<double, 2>,
        3>
        independent{{
            {{0.50, 0.20}},
            {{0.40, 0.20}},
            {{0.35, 0.25}}
        }};
    for (std::size_t phase = 0U;
         phase < 3U;
         ++phase) {
        q3[4U + phase * 2U] =
            independent[phase][0];
        q3[5U + phase * 2U] =
            independent[phase][1];
    }

    std::optional<
        fdp::FixedThreePhaseCurrentCellLinearization3D>
        three;
    fdp::NaturalVariableSnesEvaluationStatus3D
        three_status =
            fdp::NaturalVariableSnesEvaluationStatus3D::
                success;
    require_real_collective(
        fdp::
            evaluate_pr76_three_phase_production_cell_3d<
                Closure3>(
                mesh::LocalIndex{0U},
                mesh::GlobalEntityId{UINT64_C(10)},
                q3,
                three_layout,
                real_component_ids(),
                &three_context,
                &three,
                &three_status) ==
                PETSC_SUCCESS &&
            three_status ==
                fdp::NaturalVariableSnesEvaluationStatus3D::
                    success &&
            three.has_value() &&
            three->fugacity.residual_count() == 6U &&
            three->saturation_constitutive
                    .relative_permeability[0] > 0.0,
        "real PR76 3P production bridge failed");
}

std::array<double, 4>
read_owned_real_state(
    Vec state,
    int rank) {
    std::array<double, 4> result{};
    const PetscInt start =
        static_cast<PetscInt>(
            rank * 4);
    for (std::size_t slot = 0U;
         slot < result.size();
         ++slot) {
        const PetscInt index =
            start +
            static_cast<PetscInt>(slot);
        PetscScalar value{};
        require_real_collective(
            VecGetValues(
                state,
                1,
                &index,
                &value) ==
                PETSC_SUCCESS,
            "failed to read real PR76 converged state");
        result[slot] =
            static_cast<double>(
                PetscRealPart(value));
    }
    return result;
}

std::array<double, 4>
real_owned_history_signature(
    const std::vector<
        fdp::SinglePhaseSnesCellInput3D>& cells,
    int rank) {
    const mesh::GlobalEntityId stable{
        rank == 0
            ? UINT64_C(10)
            : UINT64_C(20)};
    const auto it =
        std::find_if(
            cells.begin(),
            cells.end(),
            [&](const auto& cell) {
                return cell.cell_global == stable;
            });
    require_real_collective(
        it != cells.end() &&
            it->previous_component_accumulation
                .has_value() &&
            it->previous_energy_accumulation
                .has_value() &&
            it->previous_component_accumulation
                    ->component_accumulation_mol_per_bulk_m3
                    .size() == 3U,
        "real PR76 accepted history signature is incomplete");

    std::array<double, 4> result{};
    for (std::size_t component = 0U;
         component < 3U;
         ++component) {
        result[component] =
            it->previous_component_accumulation
                ->component_accumulation_mol_per_bulk_m3[
                    component];
    }
    result[3] =
        it->previous_energy_accumulation
            ->total_internal_energy_j_per_bulk_m3;
    return result;
}

template <typename Closure>
struct RealPr76PostSnesPtReviewContext {
    fl::Pr76PtFlashBackend* backend{};
    fdp::Pr76SinglePhaseProductionCellEvaluatorContext3D<
        Closure>* evaluator_context{};
    int rank{-1};
    std::size_t equal_cardinality_stable_scans{};
    std::vector<
        fdp::PostSnesPtFlashSourceCellSnapshot3D>
        scanned_sources;
};

template <typename Closure>
PetscErrorCode real_pr76_post_snes_pt_review(
    const fdp::AdaptiveTimestepAttemptRequest3D&,
    Vec converged_state,
    const fdp::NaturalVariableSnesSolveReport3D&
        solve_report,
    void* raw_context,
    fdp::AdaptiveTimestepAttemptOutcome3D* outcome,
    std::size_t* phase_transition_restarts,
    std::vector<
        fdp::PostSnesPhaseTransitionProposal3D>*
            local_owned_proposals) {
    if (converged_state == nullptr ||
        raw_context == nullptr ||
        outcome == nullptr ||
        phase_transition_restarts == nullptr ||
        local_owned_proposals == nullptr) {
        return PETSC_ERR_ARG_NULL;
    }
    if (static_cast<int>(
            solve_report.converged_reason()) <= 0 ||
        solve_report.function_domain_errors() != 0 ||
        solve_report.jacobian_domain_errors() != 0) {
        return PETSC_ERR_ARG_INCOMP;
    }

    local_owned_proposals->clear();
    auto* context =
        static_cast<
            RealPr76PostSnesPtReviewContext<Closure>*>(
                raw_context);
    if (context->backend == nullptr ||
        context->evaluator_context == nullptr ||
        context->rank < 0) {
        return PETSC_ERR_ARG_INCOMP;
    }

    fdp::PostSnesPtFlashSourceCellSnapshot3D
        source;
    try {
        const std::uint64_t stable =
            context->rank == 0
                ? UINT64_C(10)
                : UINT64_C(20);
        const auto converged =
            read_owned_real_state(
                converged_state,
                context->rank);
        auto evaluated =
            real_direct_cell_typed(
                stable,
                std::span<const double>{
                    converged.data(),
                    converged.size()},
                context->evaluator_context);

        source.cell_global =
            mesh::GlobalEntityId{stable};
        source.source_phase_count = 1U;
        source.pressure_pa =
            evaluated->state.reference_pressure_pa();
        source.temperature_k =
            evaluated->state.temperature_k();
        source.component_ids.assign(
            evaluated->state.component_ids().begin(),
            evaluated->state.component_ids().end());
        source.overall_composition.assign(
            evaluated->state.phase_composition().begin(),
            evaluated->state.phase_composition().end());
    } catch (const std::exception&) {
        return PETSC_ERR_LIB;
    }

    std::optional<
        fdp::PostSnesPhaseTransitionProposal3D>
        proposal;
    fdp::PostSnesPhaseTransitionScanStatus3D
        scan_status =
            fdp::PostSnesPhaseTransitionScanStatus3D::
                indeterminate;
    const PetscErrorCode scan_error =
        fdp::scan_post_snes_pt_flash_source_cell_3d(
            source,
            *context->backend,
            {},
            &proposal,
            &scan_status);
    if (scan_error != PETSC_SUCCESS) {
        return scan_error;
    }

    const int local_indeterminate =
        scan_status ==
                fdp::
                    PostSnesPhaseTransitionScanStatus3D::
                        indeterminate
        ? 1
        : 0;
    int global_indeterminate = 0;
    if (MPI_Allreduce(
            &local_indeterminate,
            &global_indeterminate,
            1,
            MPI_INT,
            MPI_MAX,
            PETSC_COMM_WORLD) != MPI_SUCCESS) {
        return PETSC_ERR_MPI;
    }

    if (proposal.has_value()) {
        local_owned_proposals->push_back(
            std::move(*proposal));
    }
    const int local_transition =
        local_owned_proposals->empty()
        ? 0
        : 1;
    int global_transition = 0;
    if (MPI_Allreduce(
            &local_transition,
            &global_transition,
            1,
            MPI_INT,
            MPI_MAX,
            PETSC_COMM_WORLD) != MPI_SUCCESS) {
        return PETSC_ERR_MPI;
    }

    *phase_transition_restarts = 0U;
    if (global_indeterminate != 0) {
        local_owned_proposals->clear();
        *outcome =
            fdp::AdaptiveTimestepAttemptOutcome3D::
                phase_set_scan_indeterminate;
        return PETSC_SUCCESS;
    }
    if (global_transition != 0) {
        *outcome =
            fdp::AdaptiveTimestepAttemptOutcome3D::
                phase_transition_proposed;
        return PETSC_SUCCESS;
    }

    try {
        context->scanned_sources.push_back(
            std::move(source));
    } catch (const std::exception&) {
        return PETSC_ERR_MEM;
    }
    ++context->equal_cardinality_stable_scans;
    *outcome =
        fdp::AdaptiveTimestepAttemptOutcome3D::
            stable_phase_set;
    return PETSC_SUCCESS;
}

template <typename Closure>
struct RealControlledTransitionReviewContext {
    fdp::Pr76SinglePhaseProductionCellEvaluatorContext3D<
        Closure>* evaluator_context{};
    int rank{-1};
};

template <typename Closure>
PetscErrorCode real_controlled_transition_review(
    const fdp::AdaptiveTimestepAttemptRequest3D&,
    Vec converged_state,
    const fdp::NaturalVariableSnesSolveReport3D&
        solve_report,
    void* raw_context,
    fdp::AdaptiveTimestepAttemptOutcome3D* outcome,
    std::size_t* phase_transition_restarts,
    std::vector<
        fdp::PostSnesPhaseTransitionProposal3D>*
            local_owned_proposals) {
    if (converged_state == nullptr ||
        raw_context == nullptr ||
        outcome == nullptr ||
        phase_transition_restarts == nullptr ||
        local_owned_proposals == nullptr ||
        static_cast<int>(
            solve_report.converged_reason()) <= 0) {
        return PETSC_ERR_ARG_INCOMP;
    }
    local_owned_proposals->clear();

    auto* context =
        static_cast<
            RealControlledTransitionReviewContext<
                Closure>*>(raw_context);
    if (context->evaluator_context == nullptr ||
        context->rank < 0) {
        return PETSC_ERR_ARG_INCOMP;
    }

    try {
        const std::uint64_t stable =
            context->rank == 0
                ? UINT64_C(10)
                : UINT64_C(20);
        const auto converged =
            read_owned_real_state(
                converged_state,
                context->rank);
        auto evaluated =
            real_direct_cell_typed(
                stable,
                std::span<const double>{
                    converged.data(),
                    converged.size()},
                context->evaluator_context);
        const auto z =
            evaluated->state.phase_composition();
        if (z.size() != 3U) {
            return PETSC_ERR_ARG_INCOMP;
        }

        const double delta =
            0.1 *
            std::min(z[0], z[1]);
        std::vector<double> first{
            z.begin(), z.end()};
        std::vector<double> second{
            z.begin(), z.end()};
        first[0] += delta;
        first[1] -= delta;
        second[0] -=
            (2.0 / 3.0) * delta;
        second[1] +=
            (2.0 / 3.0) * delta;

        mpmc::flow::
            PhaseSetTransitionCandidate
            candidate;
        candidate.source_phase_count = 1U;
        candidate.target_phase_count = 2U;
        candidate.trigger =
            mpmc::flow::
                PhaseSetTransitionTrigger::
                    stability_witness;
        candidate.status =
            mpmc::flow::
                PhaseSetTransitionCandidateStatus::
                    target_resolved;
        candidate.pressure_pa =
            evaluated->state
                .reference_pressure_pa();
        candidate.temperature_k =
            evaluated->state.temperature_k();
        candidate.component_ids.assign(
            evaluated->state.component_ids().begin(),
            evaluated->state.component_ids().end());
        candidate.evidence_profile =
            "test/adaptive-transition-handoff/v1";
        candidate.diagnostic =
            "controlled resolved two-phase handoff";
        const double density =
            evaluated->molar_density
                .molar_density_mol_per_m3;
        candidate.target_phases = {
            {
                0.4,
                std::move(first),
                0.9 * density},
            {
                0.6,
                std::move(second),
                1.1 * density}
        };

        (void)mpmc::flow::
            phase_set_transition_overall_composition(
                candidate);
        local_owned_proposals->push_back(
            {
                mesh::GlobalEntityId{stable},
                std::move(candidate)});
    } catch (const std::exception&) {
        return PETSC_ERR_LIB;
    }

    *phase_transition_restarts = 0U;
    *outcome =
        fdp::AdaptiveTimestepAttemptOutcome3D::
            phase_transition_proposed;
    return PETSC_SUCCESS;
}

PetscErrorCode real_unexpected_transition_commit(
    const fdp::AdaptiveTimestepAttemptRequest3D&,
    const fdp::AdaptiveTimestepAttemptResult3D&,
    void* raw_context) {
    if (raw_context == nullptr) {
        return PETSC_ERR_ARG_NULL;
    }
    auto* count =
        static_cast<std::size_t*>(
            raw_context);
    ++(*count);
    return PETSC_ERR_PLIB;
}

template <typename Closure>
struct RealAdaptiveTimestepHarness {
    fdp::SinglePhaseAdaptiveTimestepAttemptContext3D*
        production{};
    Vec accepted_state{};
    std::vector<
        fdp::SinglePhaseSnesCellInput3D>*
            accepted_cells{};
    fdp::Pr76SinglePhaseProductionCellEvaluatorContext3D<
        Closure>* evaluator_context{};
    int rank{-1};

    std::array<double, 4>
        initial_history{};
    std::array<double, 4>
        initial_state{};
    std::vector<double>
        attempted_dt;
    std::vector<std::array<double, 4>>
        history_before_attempt;
    std::vector<std::array<double, 4>>
        state_before_attempt;

    bool forced_recoverable_rejection{};
    std::size_t commits{};
    std::optional<
        fdp::NaturalVariableSnesSolveReport3D>
        committed_report;
};

template <typename Closure>
PetscErrorCode real_adaptive_attempt(
    const fdp::AdaptiveTimestepAttemptRequest3D&
        request,
    void* raw_context,
    fdp::AdaptiveTimestepAttemptResult3D* result) {
    if (raw_context == nullptr || result == nullptr) {
        return PETSC_ERR_ARG_NULL;
    }
    auto* context =
        static_cast<
            RealAdaptiveTimestepHarness<Closure>*>(
                raw_context);
    if (context->production == nullptr ||
        context->accepted_state == nullptr ||
        context->accepted_cells == nullptr) {
        return PETSC_ERR_ARG_NULL;
    }

    try {
        context->attempted_dt.push_back(
            request.timestep_seconds);
        context->history_before_attempt.push_back(
            real_owned_history_signature(
                *context->accepted_cells,
                context->rank));
        context->state_before_attempt.push_back(
            read_owned_real_state(
                context->accepted_state,
                context->rank));
    } catch (const std::exception&) {
        return PETSC_ERR_LIB;
    }

    PetscErrorCode error =
        fdp::
            evaluate_single_phase_adaptive_timestep_attempt_3d(
                request,
                context->production,
                result);
    if (error != PETSC_SUCCESS) {
        return error;
    }

    if (request.attempt_index == 0U) {
        if (result->outcome !=
                fdp::AdaptiveTimestepAttemptOutcome3D::
                    stable_phase_set ||
            !context->production->has_pending()) {
            return PETSC_ERR_PLIB;
        }
        error =
            context->production->discard_pending();
        if (error != PETSC_SUCCESS) {
            return error;
        }
        result->outcome =
            fdp::AdaptiveTimestepAttemptOutcome3D::
                phase_set_scan_indeterminate;
        result->phase_transition_restarts = 0U;
        context->forced_recoverable_rejection =
            true;
    }
    return PETSC_SUCCESS;
}

template <typename Closure>
PetscErrorCode real_adaptive_commit(
    const fdp::AdaptiveTimestepAttemptRequest3D&
        request,
    const fdp::AdaptiveTimestepAttemptResult3D&
        result,
    void* raw_context) {
    if (raw_context == nullptr) {
        return PETSC_ERR_ARG_NULL;
    }
    auto* context =
        static_cast<
            RealAdaptiveTimestepHarness<Closure>*>(
                raw_context);
    if (context->production == nullptr ||
        context->accepted_state == nullptr ||
        context->accepted_cells == nullptr ||
        context->evaluator_context == nullptr ||
        result.outcome !=
            fdp::AdaptiveTimestepAttemptOutcome3D::
                stable_phase_set) {
        return PETSC_ERR_ARG_INCOMP;
    }

    Vec pending_solution = nullptr;
    std::optional<
        fdp::NaturalVariableSnesSolveReport3D>
        pending_report;
    PetscErrorCode error =
        context->production->take_pending(
            request,
            &pending_solution,
            &pending_report);
    if (error != PETSC_SUCCESS ||
        pending_solution == nullptr ||
        !pending_report.has_value()) {
        if (pending_solution != nullptr) {
            (void)VecDestroy(
                &pending_solution);
        }
        return error != PETSC_SUCCESS
            ? error
            : PETSC_ERR_PLIB;
    }

    try {
        const std::uint64_t stable =
            context->rank == 0
                ? UINT64_C(10)
                : UINT64_C(20);
        const auto converged =
            read_owned_real_state(
                pending_solution,
                context->rank);
        auto evaluated =
            real_direct_cell_typed(
                stable,
                std::span<const double>{
                    converged.data(),
                    converged.size()},
                context->evaluator_context);

        auto it =
            std::find_if(
                context->accepted_cells->begin(),
                context->accepted_cells->end(),
                [&](const auto& cell) {
                    return cell.cell_global ==
                        mesh::GlobalEntityId{stable};
                });
        if (it ==
                context->accepted_cells->end() ||
            !it->previous_component_accumulation
                .has_value() ||
            !it->previous_energy_accumulation
                .has_value()) {
            (void)VecDestroy(
                &pending_solution);
            return PETSC_ERR_ARG_INCOMP;
        }

        const auto next_component =
            flow::
                build_single_phase_component_accumulation(
                    evaluated->state,
                    it->porosity);
        const auto next_energy =
            flow::
                build_single_phase_energy_accumulation_snapshot(
                    evaluated->state,
                    it->porosity,
                    evaluated->transport,
                    evaluated->caloric,
                    evaluated->rock);

        error =
            VecCopy(
                pending_solution,
                context->accepted_state);
        if (error != PETSC_SUCCESS) {
            (void)VecDestroy(
                &pending_solution);
            return error;
        }

        it->previous_component_accumulation =
            next_component;
        it->previous_energy_accumulation =
            next_energy;
        context->committed_report =
            std::move(pending_report);
        ++context->commits;
    } catch (const std::exception&) {
        (void)VecDestroy(
            &pending_solution);
        return PETSC_ERR_LIB;
    }

    return VecDestroy(
        &pending_solution);
}

} // namespace

void pr76_production_fully_implicit_transient_test() {
    int rank = -1;
    int size = -1;
    if (MPI_Comm_rank(
            PETSC_COMM_WORLD,
            &rank) != MPI_SUCCESS ||
        MPI_Comm_size(
            PETSC_COMM_WORLD,
            &size) != MPI_SUCCESS) {
        throw std::runtime_error(
            "MPI rank/size query failed for real PR76 transient regression");
    }
    require_real_collective(
        size == 2,
        "real PR76 transient regression requires two MPI ranks");

    const auto parameters =
        real_pr76_parameters();
    const auto model =
        th::Pr76Phase<double>::
            from_parameters(parameters);
    auto closure =
        flow::
            make_pr76_methane_ethane_propane_property_closure(
                model,
                std::vector<
                    th::Pr76SelectedPhase>{
                    {0U, {}}});
    using Closure =
        decltype(closure);

    check_real_cardinality_bridges(
        model,
        &closure);

    fdp::Pr76SinglePhaseProductionCellEvaluatorContext3D<
        Closure>
        evaluator_context{
            &closure,
            1.0,
            {&disabled_rock_storage, nullptr},
            {}};

    fl::Pr76VleEvaluator pt_evaluator(
        model);
    fl::Pr76PtFlashBackend pt_backend(
        pt_evaluator);
    RealPr76PostSnesPtReviewContext<Closure>
        pt_review_context{
            &pt_backend,
            &evaluator_context,
            rank,
            0U,
            {}};

    check_real_face_flux(
        &evaluator_context);

    auto partition =
        real_partition(rank);
    auto schedule =
        real_schedule(rank);
    auto dof_layout =
        real_dof_layout(rank);
    auto dof_numbering =
        real_numbering(
            rank,
            dof_layout,
            partition);
    auto cell_bridge =
        real_cell_bridge(rank);
    auto cell_pattern =
        real_cell_pattern(rank);

    auto previous_cells =
        real_cells(
            rank,
            &evaluator_context);
    auto accepted_cells =
        previous_cells;
    auto faces =
        real_faces(rank);

    std::array<double, 4>
        local_previous_total{};
    {
        const auto stable =
            rank == 0
                ? UINT64_C(10)
                : UINT64_C(20);
        const auto previous =
            real_previous_state(stable);
        auto evaluated =
            real_direct_cell_typed(
                stable,
                previous,
                &evaluator_context);
        const double volume =
            stable == UINT64_C(10)
                ? 2.0
                : 5.0;
        const double porosity =
            stable == UINT64_C(10)
                ? 0.25
                : 0.30;
        const auto component =
            flow::
                build_single_phase_component_accumulation(
                    evaluated->state,
                    porosity);
        const auto energy =
            flow::
                build_single_phase_energy_accumulation_snapshot(
                    evaluated->state,
                    porosity,
                    evaluated->transport,
                    evaluated->caloric,
                    evaluated->rock);
        for (std::size_t i = 0U;
             i < 3U;
             ++i) {
            local_previous_total[i] =
                component
                    .component_accumulation_mol_per_bulk_m3[i] *
                volume;
        }
        local_previous_total[3] =
            energy.total_internal_energy_j_per_bulk_m3 *
            volume;
    }

    std::optional<
        fdp::SinglePhaseSnesAssemblyContext3D>
        context;
    constexpr double dt_seconds = 0.1;
    PetscErrorCode error =
        fdp::SinglePhaseSnesAssemblyContext3D::
            create(
                PETSC_COMM_WORLD,
                schedule,
                partition,
                dof_layout,
                dof_numbering,
                cell_bridge,
                cell_pattern,
                "real_pr76_natural_state_1p",
                dt_seconds,
                previous_cells,
                faces,
                {
                    &fdp::
                        evaluate_pr76_single_phase_production_cell_3d<
                            Closure>,
                    &evaluator_context},
                &context);
    require_real_collective(
        error == PETSC_SUCCESS &&
            context.has_value(),
        "failed to create real PR76 single-phase production context");

    Vec initial = nullptr;
    require_real_collective(
        VecCreateMPI(
            PETSC_COMM_WORLD,
            4,
            8,
            &initial) ==
            PETSC_SUCCESS,
        "failed to create real PR76 initial vector");
    set_real_previous_state(
        initial,
        rank);

    std::optional<
        fdp::CompleteNaturalVariableAssemblySnapshot3D>
        initial_assembly;
    fdp::NaturalVariableSnesEvaluationStatus3D
        initial_status =
            fdp::NaturalVariableSnesEvaluationStatus3D::
                success;
    error =
        context->evaluate_complete_assembly(
            initial,
            &initial_assembly,
            &initial_status);
    require_real_collective(
        error == PETSC_SUCCESS &&
            initial_status ==
                fdp::NaturalVariableSnesEvaluationStatus3D::
                    success &&
            initial_assembly.has_value(),
        "real PR76 initial fully implicit assembly failed");

    double local_initial_residual2 = 0.0;
    for (const auto& entry :
         initial_assembly->residual_entries()) {
        local_initial_residual2 +=
            entry.native_value *
            entry.native_value;
    }
    double global_initial_residual2 = 0.0;
    require_real_collective(
        MPI_Allreduce(
            &local_initial_residual2,
            &global_initial_residual2,
            1,
            MPI_DOUBLE,
            MPI_SUM,
            PETSC_COMM_WORLD) ==
            MPI_SUCCESS &&
            global_initial_residual2 > 0.0,
        "real PR76 transient initial state did not contain nonzero internal-flow residual");

    // Fresh central perturbation of the complete distributed residual.
    // Only cell10 (owned by rank0) is perturbed; rank1 therefore checks the
    // off-diagonal response of cell20 through the same authoritative face.
    const std::array<std::size_t, 3>
        perturb_columns{0U, 1U, 2U};
    const std::array<double, 3>
        perturb_steps{100.0, 1.0e-3, 1.0e-6};

    for (std::size_t probe = 0U;
         probe < perturb_columns.size();
         ++probe) {
        Vec plus = nullptr;
        Vec minus = nullptr;
        require_real_collective(
            VecDuplicate(initial, &plus) ==
                    PETSC_SUCCESS &&
                VecDuplicate(initial, &minus) ==
                    PETSC_SUCCESS &&
                VecCopy(initial, plus) ==
                    PETSC_SUCCESS &&
                VecCopy(initial, minus) ==
                    PETSC_SUCCESS,
            "failed to create real PR76 finite-difference states");

        bool local_perturb_ok = true;
        if (rank == 0) {
            const PetscInt index =
                static_cast<PetscInt>(
                    perturb_columns[probe]);
            PetscScalar base{};
            local_perturb_ok =
                VecGetValues(
                    initial,
                    1,
                    &index,
                    &base) ==
                PETSC_SUCCESS;
            if (local_perturb_ok) {
                const PetscScalar positive =
                    base + perturb_steps[probe];
                const PetscScalar negative =
                    base - perturb_steps[probe];
                local_perturb_ok =
                    VecSetValues(
                        plus,
                        1,
                        &index,
                        &positive,
                        INSERT_VALUES) ==
                            PETSC_SUCCESS &&
                    VecSetValues(
                        minus,
                        1,
                        &index,
                        &negative,
                        INSERT_VALUES) ==
                            PETSC_SUCCESS;
            }
        }
        require_real_collective(
            local_perturb_ok,
            "failed to insert real PR76 finite-difference perturbation");
        require_real_collective(
            VecAssemblyBegin(plus) ==
                    PETSC_SUCCESS &&
                VecAssemblyEnd(plus) ==
                    PETSC_SUCCESS &&
                VecAssemblyBegin(minus) ==
                    PETSC_SUCCESS &&
                VecAssemblyEnd(minus) ==
                    PETSC_SUCCESS,
            "failed to assemble real PR76 finite-difference states");

        std::optional<
            fdp::CompleteNaturalVariableAssemblySnapshot3D>
            plus_assembly;
        std::optional<
            fdp::CompleteNaturalVariableAssemblySnapshot3D>
            minus_assembly;
        auto plus_status =
            fdp::NaturalVariableSnesEvaluationStatus3D::
                success;
        auto minus_status =
            fdp::NaturalVariableSnesEvaluationStatus3D::
                success;
        require_real_collective(
            context->evaluate_complete_assembly(
                plus,
                &plus_assembly,
                &plus_status) ==
                    PETSC_SUCCESS &&
                context->evaluate_complete_assembly(
                    minus,
                    &minus_assembly,
                    &minus_status) ==
                    PETSC_SUCCESS &&
                plus_status ==
                    fdp::NaturalVariableSnesEvaluationStatus3D::
                        success &&
                minus_status ==
                    fdp::NaturalVariableSnesEvaluationStatus3D::
                        success &&
                plus_assembly.has_value() &&
                minus_assembly.has_value(),
            "real PR76 finite-difference residual evaluation failed");

        const mesh::GlobalEntityId local_cell{
            rank == 0
                ? UINT64_C(10)
                : UINT64_C(20)};
        for (std::size_t equation_slot = 0U;
             equation_slot < 4U;
             ++equation_slot) {
            const auto kind =
                equation_slot < 3U
                    ? fdp::NaturalVariableEquationKind3D::
                          component_conservation
                    : fdp::NaturalVariableEquationKind3D::
                          energy_conservation;
            const auto& positive =
                find_real_residual_entry(
                    *plus_assembly,
                    local_cell,
                    equation_slot,
                    kind);
            const auto& negative =
                find_real_residual_entry(
                    *minus_assembly,
                    local_cell,
                    equation_slot,
                    kind);
            const double finite_difference =
                (positive.native_value -
                 negative.native_value) /
                (2.0 * perturb_steps[probe]);
            const auto& analytic =
                find_real_jacobian_entry(
                    *initial_assembly,
                    local_cell,
                    equation_slot,
                    kind,
                    mesh::GlobalEntityId{
                        UINT64_C(10)},
                    perturb_columns[probe]);
            near_real_collective(
                analytic.value,
                finite_difference,
                8.0e-4,
                2.0e-8,
                "real PR76 complete residual Jacobian disagrees with central perturbation");
        }

        require_real_collective(
            VecDestroy(&plus) ==
                    PETSC_SUCCESS &&
                VecDestroy(&minus) ==
                    PETSC_SUCCESS,
            "real PR76 finite-difference state cleanup failed");
    }

    {
        RealControlledTransitionReviewContext<Closure>
            controlled_review{
                &evaluator_context,
                rank};
        std::optional<
            fdp::SinglePhaseAdaptiveTimestepAttemptContext3D>
            handoff_attempt;
        error =
            fdp::SinglePhaseAdaptiveTimestepAttemptContext3D::
                create(
                    PETSC_COMM_WORLD,
                    &schedule,
                    &partition,
                    &dof_layout,
                    &dof_numbering,
                    &cell_bridge,
                    &cell_pattern,
                    "real_pr76_natural_state_1p",
                    &accepted_cells,
                    &faces,
                    {
                        &fdp::
                            evaluate_pr76_single_phase_production_cell_3d<
                                Closure>,
                        &evaluator_context},
                    initial,
                    {
                        &real_controlled_transition_review<
                            Closure>,
                        &controlled_review},
                    &handoff_attempt);
        require_real_collective(
            error == PETSC_SUCCESS &&
                handoff_attempt.has_value(),
            "failed to create controlled transition handoff attempt");

        fdp::AdaptiveTimestepControllerOptions3D
            handoff_options;
        handoff_options.minimum_timestep_seconds =
            1.0e-4;
        handoff_options.maximum_timestep_seconds =
            1.0;
        handoff_options.maximum_retries = 4U;

        std::size_t unexpected_commits = 0U;
        std::optional<
            fdp::AdaptiveTimestepControllerReport3D>
            handoff_report;
        error =
            fdp::solve_adaptive_timestep_3d(
                0.1,
                handoff_options,
                {
                    &fdp::
                        evaluate_single_phase_adaptive_timestep_attempt_3d,
                    &*handoff_attempt,
                    &real_unexpected_transition_commit,
                    &unexpected_commits},
                &handoff_report);
        require_real_collective(
            error == PETSC_SUCCESS &&
                handoff_report.has_value() &&
                handoff_report->outcome ==
                    fdp::
                        AdaptiveTimestepControllerOutcome3D::
                            phase_transition_handoff_required &&
                handoff_report->attempts.size() == 1U &&
                handoff_report->attempts.front()
                        .decision ==
                    fdp::AdaptiveTimestepDecision3D::
                        handoff_phase_transition &&
                unexpected_commits == 0U &&
                handoff_attempt
                    ->has_pending_transition(),
            "controlled phase proposal was not preserved for outer handoff");

        Vec handoff_state = nullptr;
        std::optional<
            fdp::NaturalVariableSnesSolveReport3D>
            handoff_solve_report;
        std::vector<
            fdp::PostSnesPhaseTransitionProposal3D>
            handoff_proposals;
        error =
            handoff_attempt
                ->take_pending_transition(
                    handoff_report->attempts.front()
                        .request,
                    &handoff_state,
                    &handoff_solve_report,
                    &handoff_proposals);
        require_real_collective(
            error == PETSC_SUCCESS &&
                handoff_state != nullptr &&
                handoff_solve_report.has_value() &&
                static_cast<int>(
                    handoff_solve_report
                        ->converged_reason()) > 0 &&
                handoff_proposals.size() == 1U &&
                handoff_proposals.front()
                        .candidate
                        .source_phase_count == 1U &&
                handoff_proposals.front()
                        .candidate
                        .target_phase_count == 2U &&
                handoff_proposals.front()
                        .candidate
                        .evidence_profile ==
                    "test/adaptive-transition-handoff/v1" &&
                !handoff_attempt
                    ->has_pending_transition(),
            "controlled transition proposal/state/report handoff changed");

        const auto handoff_state_values =
            read_owned_real_state(
                handoff_state,
                rank);
        const auto projected_feed =
            mpmc::flow::
                phase_set_transition_overall_composition(
                    handoff_proposals.front()
                        .candidate);
        auto handoff_eval =
            real_direct_cell_typed(
                rank == 0
                    ? UINT64_C(10)
                    : UINT64_C(20),
                std::span<const double>{
                    handoff_state_values.data(),
                    handoff_state_values.size()},
                &evaluator_context);
        for (std::size_t component = 0U;
             component < projected_feed.size();
             ++component) {
            near_real_collective(
                projected_feed[component],
                handoff_eval->state
                    .phase_composition()[component],
                1.0e-14,
                1.0e-14,
                "transition handoff changed candidate material balance");
        }
        require_real_collective(
            VecDestroy(
                &handoff_state) ==
                    PETSC_SUCCESS &&
                real_owned_history_signature(
                    accepted_cells,
                    rank) ==
                    real_owned_history_signature(
                        previous_cells,
                        rank) &&
                read_owned_real_state(
                    initial,
                    rank) ==
                    real_previous_state(
                        rank == 0
                            ? UINT64_C(10)
                            : UINT64_C(20)),
            "transition handoff advanced accepted state/history");
    }

    std::optional<
        fdp::SinglePhaseAdaptiveTimestepAttemptContext3D>
        production_attempt;
    error =
        fdp::SinglePhaseAdaptiveTimestepAttemptContext3D::
            create(
                PETSC_COMM_WORLD,
                &schedule,
                &partition,
                &dof_layout,
                &dof_numbering,
                &cell_bridge,
                &cell_pattern,
                "real_pr76_natural_state_1p",
                &accepted_cells,
                &faces,
                {
                    &fdp::
                        evaluate_pr76_single_phase_production_cell_3d<
                            Closure>,
                    &evaluator_context},
                initial,
                {
                    &real_pr76_post_snes_pt_review<
                        Closure>,
                    &pt_review_context},
                &production_attempt);
    require_real_collective(
        error == PETSC_SUCCESS &&
            production_attempt.has_value(),
        "failed to create real PR76 adaptive production attempt bridge");

    RealAdaptiveTimestepHarness<Closure>
        harness;
    harness.production =
        &*production_attempt;
    harness.accepted_state =
        initial;
    harness.accepted_cells =
        &accepted_cells;
    harness.evaluator_context =
        &evaluator_context;
    harness.rank =
        rank;
    harness.initial_history =
        real_owned_history_signature(
            accepted_cells,
            rank);
    harness.initial_state =
        read_owned_real_state(
            initial,
            rank);

    fdp::AdaptiveTimestepControllerOptions3D
        adaptive_options;
    adaptive_options.minimum_timestep_seconds =
        1.0e-4;
    adaptive_options.maximum_timestep_seconds =
        1.0;
    adaptive_options.cutback_factor = 0.5;
    adaptive_options.growth_factor = 2.0;
    adaptive_options.maximum_retries = 4U;
    adaptive_options.growth_nonlinear_iteration_limit =
        8;
    adaptive_options
        .growth_line_search_direction_change_limit =
        1;
    adaptive_options.growth_transition_restart_limit =
        0U;

    std::optional<
        fdp::AdaptiveTimestepControllerReport3D>
        adaptive_report;
    error =
        fdp::solve_adaptive_timestep_3d(
            dt_seconds,
            adaptive_options,
            {
                &real_adaptive_attempt<Closure>,
                &harness,
                &real_adaptive_commit<Closure>,
                &harness},
            &adaptive_report);
    require_real_collective(
        error == PETSC_SUCCESS &&
            adaptive_report.has_value() &&
            adaptive_report->accepted() &&
            adaptive_report->retries == 1U &&
            adaptive_report->attempts.size() == 2U &&
            harness.forced_recoverable_rejection &&
            harness.commits == 1U &&
            harness.committed_report.has_value() &&
            pt_review_context
                    .equal_cardinality_stable_scans ==
                2U &&
            pt_review_context.scanned_sources.size() ==
                2U &&
            harness.attempted_dt.size() == 2U &&
            harness.history_before_attempt.size() ==
                2U &&
            harness.state_before_attempt.size() ==
                2U &&
            harness.history_before_attempt[0] ==
                harness.initial_history &&
            harness.history_before_attempt[1] ==
                harness.initial_history &&
            harness.state_before_attempt[0] ==
                harness.initial_state &&
            harness.state_before_attempt[1] ==
                harness.initial_state &&
            adaptive_report->attempts[0].decision ==
                fdp::AdaptiveTimestepDecision3D::
                    reject_and_cutback &&
            (adaptive_report->attempts[1].decision ==
                 fdp::AdaptiveTimestepDecision3D::
                     accept_and_grow ||
             adaptive_report->attempts[1].decision ==
                 fdp::AdaptiveTimestepDecision3D::
                     accept_and_hold),
        "real PR76 adaptive retry/commit lifecycle failed");

    near_real_collective(
        harness.attempted_dt[0],
        dt_seconds,
        0.0,
        0.0,
        "real PR76 first adaptive attempt used wrong dt");
    near_real_collective(
        harness.attempted_dt[1],
        0.5 * dt_seconds,
        0.0,
        0.0,
        "real PR76 retry did not use cutback dt");
    near_real_collective(
        *adaptive_report->accepted_timestep_seconds,
        0.5 * dt_seconds,
        0.0,
        0.0,
        "real PR76 adaptive controller accepted wrong dt");

    const auto& committed_report =
        *harness.committed_report;
    require_real_collective(
        static_cast<int>(
            committed_report.converged_reason()) > 0 &&
            committed_report.function_domain_errors() ==
                0 &&
            committed_report.jacobian_domain_errors() ==
                0 &&
            committed_report.snes_type() ==
                std::string_view{SNESNEWTONLS} &&
            committed_report.line_search_type() ==
                std::string_view{SNESLINESEARCHBT} &&
            committed_report.ksp_type() ==
                std::string_view{KSPGMRES} &&
            committed_report.pc_type() ==
                std::string_view{PCASM},
        "real PR76 adaptive retry did not retain the production PETSc solver contract");

    const double accepted_dt_seconds =
        *adaptive_report->accepted_timestep_seconds;

    for (const auto& source :
         pt_review_context.scanned_sources) {
        require_real_collective(
            source.source_phase_count == 1U &&
                source.component_ids ==
                    real_component_ids() &&
                source.overall_composition.size() ==
                    real_component_ids().size() &&
                source.pressure_pa > 0.0 &&
                source.temperature_k > 0.0,
            "real PR76 post-SNES PT scanner source snapshot changed");
    }

    const auto converged =
        read_owned_real_state(
            initial,
            rank);
    const auto previous =
        real_previous_state(
            rank == 0
                ? UINT64_C(10)
                : UINT64_C(20));
    double local_change2 = 0.0;
    for (std::size_t i = 0U;
         i < converged.size();
         ++i) {
        const double delta =
            converged[i] -
            previous[i];
        local_change2 +=
            delta * delta;
    }
    double global_change2 = 0.0;
    require_real_collective(
        MPI_Allreduce(
            &local_change2,
            &global_change2,
            1,
            MPI_DOUBLE,
            MPI_SUM,
            PETSC_COMM_WORLD) ==
            MPI_SUCCESS &&
            global_change2 > 0.0,
        "real PR76 fully implicit timestep did not change state");

    std::optional<
        fdp::SinglePhaseSnesAssemblyContext3D>
        accepted_step_context;
    error =
        fdp::SinglePhaseSnesAssemblyContext3D::
            create(
                PETSC_COMM_WORLD,
                schedule,
                partition,
                dof_layout,
                dof_numbering,
                cell_bridge,
                cell_pattern,
                "real_pr76_natural_state_1p",
                accepted_dt_seconds,
                previous_cells,
                faces,
                {
                    &fdp::
                        evaluate_pr76_single_phase_production_cell_3d<
                            Closure>,
                    &evaluator_context},
                &accepted_step_context);
    require_real_collective(
        error == PETSC_SUCCESS &&
            accepted_step_context.has_value(),
        "failed to rebuild accepted-dt real PR76 context");

    std::optional<
        fdp::CompleteNaturalVariableAssemblySnapshot3D>
        final_assembly;
    fdp::NaturalVariableSnesEvaluationStatus3D
        final_status =
            fdp::NaturalVariableSnesEvaluationStatus3D::
                success;
    error =
        accepted_step_context->evaluate_complete_assembly(
            initial,
            &final_assembly,
            &final_status);
    require_real_collective(
        error == PETSC_SUCCESS &&
            final_status ==
                fdp::NaturalVariableSnesEvaluationStatus3D::
                    success &&
            final_assembly.has_value(),
        "independent real PR76 final assembly failed");

    double local_final_residual2 = 0.0;
    for (const auto& entry :
         final_assembly->residual_entries()) {
        local_final_residual2 +=
            entry.native_value *
            entry.native_value;
    }
    double global_final_residual2 = 0.0;
    require_real_collective(
        MPI_Allreduce(
            &local_final_residual2,
            &global_final_residual2,
            1,
            MPI_DOUBLE,
            MPI_SUM,
            PETSC_COMM_WORLD) ==
            MPI_SUCCESS,
        "failed to reduce real PR76 final residual");
    require_real_collective(
        std::sqrt(global_final_residual2) <=
            1.0e-6 &&
            committed_report.final_function_l2_norm() <=
                1.0e-6,
        "real PR76 fully implicit timestep did not converge on independent reassembly");

    const std::uint64_t stable =
        rank == 0
            ? UINT64_C(10)
            : UINT64_C(20);
    auto final_eval =
        real_direct_cell_typed(
            stable,
            std::span<const double>{
                converged.data(),
                converged.size()},
            &evaluator_context);

    const auto& accepted_scan_source =
        pt_review_context.scanned_sources.back();
    near_real_collective(
        accepted_scan_source.pressure_pa,
        final_eval->state.reference_pressure_pa(),
        0.0,
        0.0,
        "accepted PR76 PT scan did not use the converged pressure");
    near_real_collective(
        accepted_scan_source.temperature_k,
        final_eval->state.temperature_k(),
        0.0,
        0.0,
        "accepted PR76 PT scan did not use the converged temperature");
    for (std::size_t component = 0U;
         component <
             accepted_scan_source
                 .overall_composition.size();
         ++component) {
        near_real_collective(
            accepted_scan_source
                .overall_composition[component],
            final_eval->state
                .phase_composition()[component],
            0.0,
            0.0,
            "accepted PR76 PT scan did not use the converged overall composition");
    }
    const double volume =
        stable == UINT64_C(10)
            ? 2.0
            : 5.0;
    const double porosity =
        stable == UINT64_C(10)
            ? 0.25
            : 0.30;
    const auto final_component =
        flow::
            build_single_phase_component_accumulation(
                final_eval->state,
                porosity);
    const auto final_energy =
        flow::
            build_single_phase_energy_accumulation_snapshot(
                final_eval->state,
                porosity,
                final_eval->transport,
                final_eval->caloric,
                final_eval->rock);

    std::array<double, 4>
        local_final_total{};
    for (std::size_t i = 0U;
         i < 3U;
         ++i) {
        local_final_total[i] =
            final_component
                .component_accumulation_mol_per_bulk_m3[i] *
            volume;
    }
    local_final_total[3] =
        final_energy
            .total_internal_energy_j_per_bulk_m3 *
        volume;

    const auto committed_history =
        real_owned_history_signature(
            accepted_cells,
            rank);
    for (std::size_t component = 0U;
         component < 3U;
         ++component) {
        near_real_collective(
            committed_history[component],
            final_component
                .component_accumulation_mol_per_bulk_m3[
                    component],
            0.0,
            0.0,
            "real PR76 commit did not advance component history to the accepted state");
    }
    near_real_collective(
        committed_history[3],
        final_energy
            .total_internal_energy_j_per_bulk_m3,
        0.0,
        0.0,
        "real PR76 commit did not advance energy history to the accepted state");

    std::array<double, 4>
        global_previous_total{};
    std::array<double, 4>
        global_final_total{};
    require_real_collective(
        MPI_Allreduce(
            local_previous_total.data(),
            global_previous_total.data(),
            4,
            MPI_DOUBLE,
            MPI_SUM,
            PETSC_COMM_WORLD) ==
                MPI_SUCCESS &&
            MPI_Allreduce(
                local_final_total.data(),
                global_final_total.data(),
                4,
                MPI_DOUBLE,
                MPI_SUM,
                PETSC_COMM_WORLD) ==
                MPI_SUCCESS,
        "failed to reduce real PR76 conserved inventories");

    for (std::size_t i = 0U;
         i < 3U;
         ++i) {
        near_real_collective(
            global_final_total[i],
            global_previous_total[i],
            3.0e-8,
            1.0e-8,
            "real PR76 fully implicit component inventory is not conserved");
    }
    near_real_collective(
        global_final_total[3],
        global_previous_total[3],
        3.0e-8,
        1.0e-3,
        "real PR76 fully implicit energy inventory is not conserved");

    require_real_collective(
        VecDestroy(&initial) ==
                PETSC_SUCCESS,
        "real PR76 transient fixture cleanup failed");
}
