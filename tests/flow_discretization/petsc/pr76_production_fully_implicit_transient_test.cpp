#include <mpmc/flow/pr76_methane_ethane_propane_properties.hpp>
#include <mpmc/flow_discretization/single_phase_tpfa.hpp>
#include <mpmc/flow_discretization_petsc/complete_natural_variable_petsc_materialization.hpp>
#include <mpmc/flow_discretization_petsc/adaptive_timestep_controller.hpp>
#include <mpmc/flow_discretization_petsc/pr76_production_cell_evaluator.hpp>

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

    auto cells =
        real_cells(
            rank,
            &evaluator_context);
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
                std::move(cells),
                std::move(faces),
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

    Vec row_scaling = nullptr;
    require_real_collective(
        fdp::
            make_natural_variable_initial_row_equilibration_3d(
                PETSC_COMM_WORLD,
                *initial_assembly,
                &row_scaling) ==
                PETSC_SUCCESS &&
            row_scaling != nullptr,
        "failed to build frozen analytic row equilibration for real PR76 solve");

    Vec unused_residual = nullptr;
    Mat jacobian_template = nullptr;
    error =
        fdp::
            materialize_complete_natural_variable_petsc_system_3d(
                PETSC_COMM_WORLD,
                *initial_assembly,
                cell_bridge,
                &unused_residual,
                &jacobian_template);
    require_real_collective(
        error == PETSC_SUCCESS &&
            unused_residual != nullptr &&
            jacobian_template != nullptr,
        "failed to materialize real PR76 PETSc system");

    Vec solution = nullptr;
    std::optional<
        fdp::NaturalVariableSnesSolveReport3D>
        report;
    std::optional<
        fdp::NaturalVariableSnesFailureDiagnostics3D>
        failure_diagnostics;
    error =
        fdp::solve_natural_variable_snes_3d(
            PETSC_COMM_WORLD,
            *initial_assembly,
            initial,
            jacobian_template,
            context->snes_evaluator(),
            &solution,
            &report,
            row_scaling,
            &failure_diagnostics);
    const bool solve_ok =
        error == PETSC_SUCCESS &&
        solution != nullptr &&
        report.has_value() &&
        static_cast<int>(
            report->converged_reason()) >
            0 &&
        report->function_domain_errors() ==
            0 &&
        report->jacobian_domain_errors() ==
            0 &&
        report->snes_type() ==
            std::string_view{SNESNEWTONLS} &&
        report->line_search_type() ==
            std::string_view{SNESLINESEARCHBT} &&
        report->ksp_type() ==
            std::string_view{KSPGMRES} &&
        report->pc_type() ==
            std::string_view{PCASM};
    std::string solve_message =
        "real PR76 fully implicit PETSc solve failed";
    if (!solve_ok &&
        failure_diagnostics.has_value()) {
        solve_message +=
            " snes_reason=" +
            std::to_string(
                static_cast<int>(
                    failure_diagnostics
                        ->snes_reason)) +
            " ksp_reason=" +
            std::to_string(
                static_cast<int>(
                    failure_diagnostics
                        ->ksp_reason)) +
            " pc_failed_reason=" +
            std::to_string(
                failure_diagnostics
                    ->pc_failed_reason) +
            " asm_sub_ksp_reason=" +
            std::to_string(
                static_cast<int>(
                    failure_diagnostics
                        ->asm_sub_ksp_reason)) +
            " asm_sub_pc_failed_reason=" +
            std::to_string(
                failure_diagnostics
                    ->asm_sub_pc_failed_reason) +
            " nonlinear_iterations=" +
            std::to_string(
                failure_diagnostics
                    ->nonlinear_iterations) +
            " function_evaluations=" +
            std::to_string(
                failure_diagnostics
                    ->function_evaluations) +
            " jacobian_evaluations=" +
            std::to_string(
                failure_diagnostics
                    ->jacobian_evaluations) +
            " function_domain_errors=" +
            std::to_string(
                failure_diagnostics
                    ->function_domain_errors) +
            " jacobian_domain_errors=" +
            std::to_string(
                failure_diagnostics
                    ->jacobian_domain_errors) +
            " line_search_prechecks=" +
            std::to_string(
                failure_diagnostics
                    ->line_search_prechecks) +
            " line_search_direction_changes=" +
            std::to_string(
                failure_diagnostics
                    ->line_search_direction_changes) +
            " function_l2_norm=" +
            std::to_string(
                failure_diagnostics
                    ->function_l2_norm);
    } else if (!solve_ok) {
        solve_message +=
            " petsc_error=" +
            std::to_string(
                static_cast<int>(error));
    }
    require_real_collective(
        solve_ok,
        solve_message);

    const auto adaptive_result =
        fdp::make_adaptive_timestep_attempt_result(
            *report);
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
    const auto adaptive_decision =
        fdp::decide_adaptive_timestep_3d(
            dt_seconds,
            0U,
            adaptive_result,
            adaptive_options);
    require_real_collective(
        adaptive_result.outcome ==
                fdp::AdaptiveTimestepAttemptOutcome3D::
                    stable_phase_set &&
            (adaptive_decision.decision ==
                 fdp::AdaptiveTimestepDecision3D::
                     accept_and_grow ||
             adaptive_decision.decision ==
                 fdp::AdaptiveTimestepDecision3D::
                     accept_and_hold) &&
            adaptive_decision
                .next_timestep_seconds
                .has_value() &&
            *adaptive_decision
                 .next_timestep_seconds >=
                dt_seconds,
        "real PR76 converged SNES report was not accepted by adaptive timestep policy");

    const auto converged =
        read_owned_real_state(
            solution,
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
        fdp::CompleteNaturalVariableAssemblySnapshot3D>
        final_assembly;
    fdp::NaturalVariableSnesEvaluationStatus3D
        final_status =
            fdp::NaturalVariableSnesEvaluationStatus3D::
                success;
    error =
        context->evaluate_complete_assembly(
            solution,
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
            report->final_function_l2_norm() <=
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
        VecDestroy(&solution) ==
                PETSC_SUCCESS &&
            VecDestroy(&unused_residual) ==
                PETSC_SUCCESS &&
            MatDestroy(&jacobian_template) ==
                PETSC_SUCCESS &&
            VecDestroy(&row_scaling) ==
                PETSC_SUCCESS &&
            VecDestroy(&initial) ==
                PETSC_SUCCESS,
        "real PR76 transient fixture cleanup failed");
}
