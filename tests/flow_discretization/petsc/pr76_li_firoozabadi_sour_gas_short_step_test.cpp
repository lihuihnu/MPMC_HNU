#include <mpmc/flow/pr76_li_firoozabadi_sour_gas_properties.hpp>
#include <mpmc/flow/phase_potential_upwind.hpp>
#include <mpmc/flow_discretization/energy_face_flux.hpp>
#include <mpmc/flow_discretization/tpfa_component_molar_flux.hpp>
#include <mpmc/flow_discretization/tpfa_phase_darcy_flux.hpp>
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

namespace disc = mpmc::discretization;
namespace dp = mpmc::discretization_petsc;
namespace fd = mpmc::flow_discretization;
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



void require_comm(
    MPI_Comm comm,
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
            comm) != MPI_SUCCESS) {
        throw std::runtime_error(
            "Li-Firoozabadi decomposition-invariance MPI reduction failed");
    }
    if (global == 0) {
        throw std::runtime_error(
            std::string{message});
    }
}

void near_comm(
    MPI_Comm comm,
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
    require_comm(
        comm,
        std::isfinite(actual) &&
            std::isfinite(expected) &&
            std::abs(actual - expected) <=
                absolute + relative * scale,
        message);
}

[[nodiscard]] double
two_cell_face_area_m2() {
    return ref::dy_m * ref::dz_m;
}

[[nodiscard]] double
two_cell_transmissibility_m3() {
    return ref::permeability_m2 *
        two_cell_face_area_m2() /
        ref::dx_m;
}

[[nodiscard]]
disc::CombinedTransmissibilityAdmissibility3D
two_cell_direct_admissibility() {
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
two_cell_materialized_transmissibility() {
    return {
        mesh::LocalIndex{0U},
        disc::
            TpfaInternalFaceTransmissibilityDisposition3D::
                materialized,
        two_cell_direct_admissibility(),
        disc::TpfaStaticFaceTransmissibility3D{
            disc::
                TpfaStaticFaceTransmissibilityDisposition3D::
                    positive_harmonic_combination,
            two_cell_face_area_m2(),
            two_cell_transmissibility_m3()}};
}

[[nodiscard]]
dp::AssemblyReadyInternalConnectionRow3D
two_cell_connection_row() {
    return {
        mesh::LocalIndex{0U},
        mesh::GlobalEntityId{
            UINT64_C(100)},
        mesh::LocalIndex{0U},
        mesh::GlobalEntityId{
            UINT64_C(10)},
        mesh::LocalIndex{1U},
        mesh::GlobalEntityId{
            UINT64_C(20)},
        two_cell_transmissibility_m3()};
}

[[nodiscard]]
fdp::MixedCardinalityPhysicalSnesAuthoritativeFaceInput3D
two_cell_face_input() {
    return {
        mesh::LocalIndex{0U},
        mesh::GlobalEntityId{
            UINT64_C(100)},
        two_cell_materialized_transmissibility(),
        flow::GravityVector3D{
            0.0, 0.0, 0.0},
        flow::OwnerToNeighbourDisplacement3D{
            ref::dx_m, 0.0, 0.0},
        fd::StaticThermalFaceConductance3D{
            0.0}};
}

[[nodiscard]] mesh::Topology
make_two_cell_topology() {
    mesh::Topology::EntityIds ids;
    ids.faces = {
        mesh::GlobalEntityId{
            UINT64_C(100)}};
    ids.cells = {
        mesh::GlobalEntityId{
            UINT64_C(10)},
        mesh::GlobalEntityId{
            UINT64_C(20)}};
    return {
        std::move(ids),
        {}};
}

[[nodiscard]] mesh::PartitionSnapshot
make_two_cell_partition(
    MPI_Comm comm,
    bool distributed) {
    int rank = -1;
    int size = -1;
    if (MPI_Comm_rank(
            comm,
            &rank) != MPI_SUCCESS ||
        MPI_Comm_size(
            comm,
            &size) != MPI_SUCCESS) {
        throw std::runtime_error(
            "failed to query two-cell partition communicator");
    }
    if ((!distributed &&
         (rank != 0 || size != 1)) ||
        (distributed &&
         (rank < 0 || rank > 1 ||
          size != 2))) {
        throw std::invalid_argument(
            "invalid communicator for Li-Firoozabadi two-cell decomposition fixture");
    }

    const auto topology =
        make_two_cell_topology();
    mesh::EntityOwnerRanks owners;
    owners.faces = {
        mesh::PartitionRank{0U}};
    owners.cells = distributed
        ? std::vector<mesh::PartitionRank>{
              mesh::PartitionRank{0U},
              mesh::PartitionRank{1U}}
        : std::vector<mesh::PartitionRank>{
              mesh::PartitionRank{0U},
              mesh::PartitionRank{0U}};

    return mesh::PartitionSnapshot::create(
        topology,
        mesh::PartitionRank{
            static_cast<
                mesh::PartitionRank::value_type>(
                    rank)},
        static_cast<std::uint32_t>(
            size),
        std::move(owners));
}

[[nodiscard]]
dp::ParallelOwnedConnectionSchedule3D
make_two_cell_schedule(
    int rank,
    bool distributed) {
    const auto row =
        two_cell_connection_row();
    if (!distributed ||
        rank == 0) {
        return {
            mesh::PartitionRank{
                static_cast<
                    mesh::PartitionRank::value_type>(
                        rank)},
            distributed ? 2U : 1U,
            {row},
            {}};
    }
    return {
        mesh::PartitionRank{1U},
        2U,
        {},
        {row}};
}

[[nodiscard]]
dp::PetscMpiAijSymbolicPreallocation3D
make_two_cell_bridge(
    int rank,
    bool distributed) {
    if (!distributed) {
        return {
            mesh::PartitionRank{0U},
            1U,
            0,
            2,
            2,
            {
                mesh::LocalIndex{0U},
                mesh::LocalIndex{1U}},
            {
                mesh::GlobalEntityId{
                    UINT64_C(10)},
                mesh::GlobalEntityId{
                    UINT64_C(20)}},
            {0, 1},
            {2, 2},
            {0, 0},
            {0, 1}};
    }

    if (rank == 0) {
        return {
            mesh::PartitionRank{0U},
            2U,
            0,
            1,
            2,
            {mesh::LocalIndex{0U}},
            {
                mesh::GlobalEntityId{
                    UINT64_C(10)}},
            {0},
            {1},
            {1},
            {0, 1}};
    }

    return {
        mesh::PartitionRank{1U},
        2U,
        1,
        2,
        2,
        {mesh::LocalIndex{1U}},
        {
            mesh::GlobalEntityId{
                UINT64_C(20)}},
        {1},
        {1},
        {1},
        {0, 1}};
}

[[nodiscard]]
dp::OwnedCellStructuralColumnPatternSnapshot3D
make_two_cell_pattern(
    int rank,
    bool distributed) {
    if (!distributed) {
        return {
            mesh::PartitionRank{0U},
            1U,
            2U,
            0,
            2,
            2,
            {
                {
                    mesh::LocalIndex{0U},
                    mesh::GlobalEntityId{
                        UINT64_C(10)},
                    0,
                    0U,
                    2U,
                    0U,
                    0U},
                {
                    mesh::LocalIndex{1U},
                    mesh::GlobalEntityId{
                        UINT64_C(20)},
                    1,
                    2U,
                    2U,
                    0U,
                    0U}},
            {0, 1, 0, 1},
            {}};
    }

    if (rank == 0) {
        return {
            mesh::PartitionRank{0U},
            2U,
            2U,
            0,
            1,
            2,
            {
                {
                    mesh::LocalIndex{0U},
                    mesh::GlobalEntityId{
                        UINT64_C(10)},
                    0,
                    0U,
                    1U,
                    0U,
                    1U}},
            {0},
            {1}};
    }

    return {
        mesh::PartitionRank{1U},
        2U,
        2U,
        1,
        2,
        2,
        {
            {
                mesh::LocalIndex{1U},
                mesh::GlobalEntityId{
                    UINT64_C(20)},
                1,
                0U,
                1U,
                0U,
                1U}},
        {1},
        {0}};
}

[[nodiscard]]
std::array<std::string, 3>
two_cell_phase_identity_keys() {
    const auto map =
        make_phase_identities();
    std::array<std::string, 3>
        result;
    for (std::size_t phase = 0U;
         phase < result.size();
         ++phase) {
        result[phase] =
            map.identity(phase)
                .opaque_phase_key;
    }
    return result;
}

[[nodiscard]]
std::vector<double>
natural_variables_from_state(
    const flow::NaturalVariableCellState3P&
        state) {
    const auto& layout =
        state.layout();
    std::vector<double>
        q(
            layout.unknown_count(),
            0.0);
    q[layout.pressure_unknown_index()] =
        state.reference_pressure_pa();
    q[layout.temperature_unknown_index()] =
        state.temperature_k();

    for (std::size_t phase = 0U;
         phase < 3U;
         ++phase) {
        const auto slot =
            static_cast<
                flow::PhaseSlot3>(
                    phase);
        if (const auto saturation =
                layout
                    .independent_saturation_unknown_index(
                        slot);
            saturation.has_value()) {
            q[*saturation] =
                state.phase_saturation(
                    slot);
        }

        const auto composition =
            state.phase_composition(
                slot);
        for (std::size_t component = 0U;
             component <
                 composition.size();
             ++component) {
            if (const auto column =
                    layout
                        .independent_composition_unknown_index(
                            slot,
                            component);
                column.has_value()) {
                q[*column] =
                    composition[component];
            }
        }
    }
    return q;
}

void insert_two_cell_owned_state(
    MPI_Comm comm,
    Vec state,
    const fdp::
        VariableCardinalityNaturalVariableNumbering3D&
            numbering,
    std::span<const double> q) {
    bool local_ok = true;
    for (const auto& cell :
         numbering.cells()) {
        if (cell.owner_rank !=
            numbering.local_rank()) {
            continue;
        }
        local_ok =
            local_ok &&
            cell.scalar_count ==
                q.size();
        for (std::size_t slot = 0U;
             slot < q.size();
             ++slot) {
            const PetscInt index =
                cell.petsc_global_scalar_start +
                static_cast<PetscInt>(
                    slot);
            const PetscScalar value =
                static_cast<PetscScalar>(
                    q[slot]);
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

    require_comm(
        comm,
        local_ok,
        "failed to insert two-cell owned natural-variable state");
    require_comm(
        comm,
        VecAssemblyBegin(
            state) ==
                PETSC_SUCCESS &&
            VecAssemblyEnd(
                state) ==
                PETSC_SUCCESS,
        "failed to assemble two-cell natural-variable state");
}

struct TwoCellDecompositionSnapshot {
    std::array<
        std::vector<double>,
        2>
        natural_variables;
    std::array<std::string, 3>
        phase_identity_keys;
    std::array<double, 3>
        face_phase_volumetric_rate_m3_per_s{};
    std::array<double, 6>
        face_component_molar_rate_mol_per_s{};
    double face_energy_rate_w{};
    std::array<double, 3>
        well_phase_volumetric_rate_m3_per_s{};
    std::array<double, 6>
        well_component_source_mol_per_s{};
    double well_energy_source_w{};
    std::array<double, 6>
        global_component_inventory_mol{};
    double global_total_internal_energy_j{};
    double final_residual_l2{};
    double spatial_cross_block_abs_sum{};
};

template <typename Closure>
[[nodiscard]]
TwoCellDecompositionSnapshot
solve_two_cell_decomposition_case(
    MPI_Comm comm,
    bool distributed,
    fdp::
        Pr76ThreePhaseProductionCellEvaluatorContext3D<
            Closure>* evaluator_context,
    const flow::NaturalVariableLayout3P&
        layout,
    std::span<const double> q,
    const std::vector<std::string>&
        component_ids) {
    if (evaluator_context == nullptr) {
        throw std::invalid_argument(
            "null Li-Firoozabadi two-cell evaluator context");
    }

    int rank = -1;
    int size = -1;
    if (MPI_Comm_rank(
            comm,
            &rank) != MPI_SUCCESS ||
        MPI_Comm_size(
            comm,
            &size) != MPI_SUCCESS) {
        throw std::runtime_error(
            "failed to query two-cell solve communicator");
    }
    require_comm(
        comm,
        (!distributed &&
         rank == 0 &&
         size == 1) ||
            (distributed &&
             size == 2 &&
             rank >= 0 &&
             rank < 2),
        "Li-Firoozabadi two-cell decomposition communicator mismatch");

    const auto reference_current =
        evaluate_current(
            q,
            layout,
            component_ids,
            evaluator_context);
    const auto previous_component =
        flow::
            build_pore_volume_component_accumulation(
                reference_current.state,
                ref::porosity);
    const auto previous_energy =
        flow::
            build_pore_volume_energy_accumulation_snapshot(
                reference_current.state,
                ref::porosity,
                reference_current.transport,
                reference_current.caloric,
                reference_current.rock);

    auto partition =
        make_two_cell_partition(
            comm,
            distributed);
    auto schedule =
        make_two_cell_schedule(
            rank,
            distributed);
    auto bridge =
        make_two_cell_bridge(
            rank,
            distributed);
    auto pattern =
        make_two_cell_pattern(
            rank,
            distributed);

    std::optional<
        fdp::
            VariableCardinalityNaturalVariableNumbering3D>
        numbering;
    require_comm(
        comm,
        fdp::
            make_variable_cardinality_natural_variable_numbering_3d(
                comm,
                partition,
                6U,
                std::vector<std::size_t>{
                    3U, 3U},
                &numbering) ==
                PETSC_SUCCESS &&
            numbering.has_value() &&
            numbering
                    ->petsc_global_scalar_count() ==
                38,
        "failed to build Li-Firoozabadi two-cell natural-variable numbering");

    std::vector<
        fdp::
            MixedCardinalityPhysicalSnesCellInput3D>
        cell_inputs;
    cell_inputs.reserve(2U);
    for (std::size_t local = 0U;
         local < 2U;
         ++local) {
        fdp::FixedThreePhaseSnesCellInput3D
            input;
        input.cell =
            mesh::LocalIndex{
                static_cast<
                    mesh::LocalIndex::value_type>(
                        local)};
        input.cell_global =
            mesh::GlobalEntityId{
                local == 0U
                    ? UINT64_C(10)
                    : UINT64_C(20)};
        input.bulk_volume_m3 =
            ref::bulk_volume_m3;
        input.porosity =
            ref::porosity;
        input.frozen_layout =
            layout;
        input.component_ids =
            component_ids;

        const bool owned =
            numbering
                ->is_owned_cell(
                    input.cell);
        if (owned) {
            input.previous_component_accumulation =
                previous_component;
            input.previous_energy_accumulation =
                previous_energy;
        }
        cell_inputs.emplace_back(
            std::move(input));
    }

    std::vector<
        flow::
            FrozenActivePhaseIdentityMap>
        phase_maps;
    phase_maps.reserve(2U);
    phase_maps.push_back(
        make_phase_identities());
    phase_maps.push_back(
        make_phase_identities());

    std::vector<
        fdp::
            MixedCardinalityPhysicalSnesAuthoritativeFaceInput3D>
        face_inputs;
    if (!distributed ||
        rank == 0) {
        face_inputs.push_back(
            two_cell_face_input());
    }

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
                            "Li-Firoozabadi two-cell stationary decomposition benchmark injection fallback; unused because q_phase=0",
                            {
                                ref::
                                    specific_enthalpy_j_per_kg
                                        .begin(),
                                ref::
                                    specific_enthalpy_j_per_kg
                                        .end()}},
                    "Li-Firoozabadi-2012 two-cell stationary fixed-BHP decomposition benchmark");

    std::optional<
        fdp::
            MixedCardinalityPhysicalSnesAssemblyContext3D>
        assembly;
    require_comm(
        comm,
        fdp::
            MixedCardinalityPhysicalSnesAssemblyContext3D::
                create(
                    comm,
                    schedule,
                    partition,
                    *numbering,
                    bridge,
                    pattern,
                    ref::timestep_seconds,
                    std::move(cell_inputs),
                    std::move(phase_maps),
                    std::move(face_inputs),
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
                            evaluator_context}},
                    {},
                    wdp::
                        fixed_bhp_peaceman_well_source_binding_3d(
                            &well_context),
                    &assembly) ==
                PETSC_SUCCESS &&
            assembly.has_value(),
        "failed to create Li-Firoozabadi two-cell production assembly");

    Mat jacobian = nullptr;
    require_comm(
        comm,
        assembly
                ->create_jacobian_structure(
                    &jacobian) ==
                PETSC_SUCCESS &&
            jacobian != nullptr,
        "failed to create Li-Firoozabadi two-cell Jacobian");

    Vec initial = nullptr;
    require_comm(
        comm,
        fdp::
            create_variable_cardinality_natural_variable_vec_3d(
                comm,
                *numbering,
                &initial) ==
                PETSC_SUCCESS &&
            initial != nullptr,
        "failed to create Li-Firoozabadi two-cell state Vec");
    insert_two_cell_owned_state(
        comm,
        initial,
        *numbering,
        q);

    auto evaluator =
        assembly->snes_evaluator();
    fdp::NaturalVariableSnesEvaluationStatus3D
        jacobian_status =
            fdp::
                NaturalVariableSnesEvaluationStatus3D::
                    success;
    PetscErrorCode error =
        MatZeroEntries(
            jacobian);
    if (error == PETSC_SUCCESS) {
        error =
            evaluator.jacobian(
                initial,
                jacobian,
                evaluator.user_context,
                &jacobian_status);
    }
    if (error == PETSC_SUCCESS) {
        error =
            MatAssemblyBegin(
                jacobian,
                MAT_FINAL_ASSEMBLY);
    }
    if (error == PETSC_SUCCESS) {
        error =
            MatAssemblyEnd(
                jacobian,
                MAT_FINAL_ASSEMBLY);
    }
    require_comm(
        comm,
        error == PETSC_SUCCESS &&
            jacobian_status ==
                fdp::
                    NaturalVariableSnesEvaluationStatus3D::
                        success,
        "Li-Firoozabadi two-cell spatial Jacobian assembly failed");

    double local_cross_abs_sum = 0.0;
    for (std::size_t local = 0U;
         local < 2U;
         ++local) {
        const auto cell =
            mesh::LocalIndex{
                static_cast<
                    mesh::LocalIndex::value_type>(
                        local)};
        if (!numbering
                 ->is_owned_cell(
                     cell)) {
            continue;
        }
        const auto other =
            mesh::LocalIndex{
                static_cast<
                    mesh::LocalIndex::value_type>(
                        1U - local)};
        const auto& row_record =
            numbering->cell(
                cell);
        const auto& column_record =
            numbering->cell(
                other);

        for (std::size_t row_slot = 0U;
             row_slot <
                 row_record.scalar_count;
             ++row_slot) {
            const PetscInt row =
                row_record
                    .petsc_global_scalar_start +
                static_cast<PetscInt>(
                    row_slot);
            for (std::size_t column_slot = 0U;
                 column_slot <
                     column_record.scalar_count;
                 ++column_slot) {
                const PetscInt column =
                    column_record
                        .petsc_global_scalar_start +
                    static_cast<PetscInt>(
                        column_slot);
                PetscScalar value = 0.0;
                require_comm(
                    comm,
                    MatGetValues(
                        jacobian,
                        1,
                        &row,
                        1,
                        &column,
                        &value) ==
                        PETSC_SUCCESS,
                    "failed to inspect Li-Firoozabadi cross-cell Jacobian block");
                local_cross_abs_sum +=
                    std::abs(
                        static_cast<double>(
                            PetscRealPart(
                                value)));
            }
        }
    }

    double global_cross_abs_sum = 0.0;
    require_comm(
        comm,
        MPI_Allreduce(
            &local_cross_abs_sum,
            &global_cross_abs_sum,
            1,
            MPI_DOUBLE,
            MPI_SUM,
            comm) ==
            MPI_SUCCESS &&
            std::isfinite(
                global_cross_abs_sum) &&
            global_cross_abs_sum >
                0.0,
        "Li-Firoozabadi internal face did not create a cross-cell Jacobian coupling");

    Vec solution = nullptr;
    std::optional<
        fdp::
            VariableCardinalityNaturalVariableSnesSolveReport3D>
        report;
    require_comm(
        comm,
        fdp::
            solve_variable_cardinality_natural_variable_snes_3d(
                comm,
                *numbering,
                initial,
                jacobian,
                evaluator,
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
                PCASM,
        "Li-Firoozabadi two-cell production SNES did not converge with NewtonLS/GMRES/ASM");

    Vec residual = nullptr;
    require_comm(
        comm,
        VecDuplicate(
            solution,
            &residual) ==
                PETSC_SUCCESS &&
            VecSet(
                residual,
                PetscScalar{0.0}) ==
                PETSC_SUCCESS,
        "failed to allocate Li-Firoozabadi two-cell final residual");

    fdp::NaturalVariableSnesEvaluationStatus3D
        residual_status =
            fdp::
                NaturalVariableSnesEvaluationStatus3D::
                    success;
    error =
        evaluator.function(
            solution,
            residual,
            evaluator.user_context,
            &residual_status);
    if (error == PETSC_SUCCESS) {
        error =
            VecAssemblyBegin(
                residual);
    }
    if (error == PETSC_SUCCESS) {
        error =
            VecAssemblyEnd(
                residual);
    }
    PetscReal residual_norm = 0.0;
    if (error == PETSC_SUCCESS) {
        error =
            VecNorm(
                residual,
                NORM_2,
                &residual_norm);
    }
    require_comm(
        comm,
        error == PETSC_SUCCESS &&
            residual_status ==
                fdp::
                    NaturalVariableSnesEvaluationStatus3D::
                        success &&
            static_cast<double>(
                residual_norm) <=
                3.0e-8,
        "Li-Firoozabadi two-cell final residual is not closed");

    std::vector<std::optional<
        fdp::
            MixedCardinalityPhysicalCurrentCellLinearization3D>>
        current;
    std::vector<double>
        porosities;
    fdp::NaturalVariableSnesEvaluationStatus3D
        current_status =
            fdp::
                NaturalVariableSnesEvaluationStatus3D::
                    success;
    require_comm(
        comm,
        assembly
                ->evaluate_local_cells_for_phase_transition(
                    solution,
                    &current,
                    &porosities,
                    &current_status) ==
                PETSC_SUCCESS &&
            current_status ==
                fdp::
                    NaturalVariableSnesEvaluationStatus3D::
                        success &&
            current.size() ==
                2U &&
            porosities.size() ==
                2U &&
            current[0].has_value() &&
            current[1].has_value(),
        "failed to evaluate Li-Firoozabadi two-cell converged local copies");

    TwoCellDecompositionSnapshot snapshot;
    snapshot.phase_identity_keys =
        two_cell_phase_identity_keys();
    snapshot.final_residual_l2 =
        static_cast<double>(
            residual_norm);
    snapshot.spatial_cross_block_abs_sum =
        global_cross_abs_sum;

    std::array<double, 6>
        local_component_inventory{};
    double local_energy = 0.0;

    for (std::size_t local = 0U;
         local < 2U;
         ++local) {
        const auto* typed =
            std::get_if<
                fdp::
                    FixedThreePhaseCurrentCellLinearization3D>(
                        &*current[local]);
        require_comm(
            comm,
            typed != nullptr,
            "Li-Firoozabadi two-cell solve changed frozen three-phase cardinality");

        snapshot.natural_variables[local] =
            natural_variables_from_state(
                typed->state);

        const auto cell =
            mesh::LocalIndex{
                static_cast<
                    mesh::LocalIndex::value_type>(
                        local)};
        if (!numbering
                 ->is_owned_cell(
                     cell)) {
            continue;
        }

        const auto component =
            flow::
                build_pore_volume_component_accumulation(
                    typed->state,
                    ref::porosity);
        const auto energy =
            flow::
                build_pore_volume_energy_accumulation_snapshot(
                    typed->state,
                    ref::porosity,
                    typed->transport,
                    typed->caloric,
                    typed->rock);
        for (std::size_t index = 0U;
             index < 6U;
             ++index) {
            local_component_inventory[index] +=
                component
                    .component_accumulation_mol_per_bulk_m3[
                        index] *
                ref::bulk_volume_m3;
        }
        local_energy +=
            energy
                .total_internal_energy_j_per_bulk_m3 *
            ref::bulk_volume_m3;
    }

    const int component_reduce_error =
        MPI_Allreduce(
            local_component_inventory.data(),
            snapshot
                .global_component_inventory_mol
                .data(),
            6,
            MPI_DOUBLE,
            MPI_SUM,
            comm);
    const int energy_reduce_error =
        MPI_Allreduce(
            &local_energy,
            &snapshot
                .global_total_internal_energy_j,
            1,
            MPI_DOUBLE,
            MPI_SUM,
            comm);
    require_comm(
        comm,
        component_reduce_error ==
                MPI_SUCCESS &&
            energy_reduce_error ==
                MPI_SUCCESS,
        "failed to reduce Li-Firoozabadi two-cell conserved totals");

    const auto& owner =
        std::get<
            fdp::
                FixedThreePhaseCurrentCellLinearization3D>(
                    *current[0]);
    const auto& neighbour =
        std::get<
            fdp::
                FixedThreePhaseCurrentCellLinearization3D>(
                    *current[1]);
    const auto owner_mobility =
        flow::
            build_local_phase_mobility_linearization(
                owner.state,
                owner.transport,
                owner.saturation_constitutive);
    const auto neighbour_mobility =
        flow::
            build_local_phase_mobility_linearization(
                neighbour.state,
                neighbour.transport,
                neighbour.saturation_constitutive);
    const auto potential =
        flow::
            build_two_cell_phase_potential_upwind_linearization(
                owner_mobility,
                neighbour_mobility,
                flow::GravityVector3D{
                    0.0, 0.0, 0.0},
                flow::OwnerToNeighbourDisplacement3D{
                    ref::dx_m, 0.0, 0.0});
    const auto phase_flux =
        fd::
            build_materialized_tpfa_internal_face_phase_darcy_flux(
                two_cell_materialized_transmissibility(),
                potential);
    const auto component_flux =
        fd::
            build_materialized_tpfa_internal_face_component_molar_flux(
                phase_flux,
                owner.state,
                owner.molar_density,
                neighbour.state,
                neighbour.molar_density);
    const auto energy_rate =
        fd::
            build_internal_energy_face_rate(
                phase_flux,
                owner.transport,
                owner.caloric,
                neighbour.transport,
                neighbour.caloric,
                fd::
                    StaticThermalFaceConductance3D{
                        0.0});

    for (std::size_t phase = 0U;
         phase < 3U;
         ++phase) {
        snapshot
            .face_phase_volumetric_rate_m3_per_s[
                phase] =
            phase_flux
                .phase[phase]
                .volumetric_flux_m3_per_s;
    }
    for (std::size_t component = 0U;
         component < 6U;
         ++component) {
        snapshot
            .face_component_molar_rate_mol_per_s[
                component] =
            component_flux
                .component_molar_flux_mol_per_s[
                    component];
    }
    snapshot.face_energy_rate_w =
        energy_rate
            .total_energy_rate_w;

    std::array<double, 3>
        local_well_phase_rate{};
    std::array<double, 6>
        local_well_component_source{};
    double local_well_energy_source = 0.0;
    if (numbering
            ->is_owned_cell(
                mesh::LocalIndex{0U})) {
        const auto well =
            wdp::
                build_fixed_bhp_peaceman_well_source_3d(
                    well_context,
                    *current[0]);
        for (std::size_t phase = 0U;
             phase < 3U;
             ++phase) {
            local_well_phase_rate[phase] =
                well
                    .phase_volumetric_rate_m3_per_s[
                        phase];
        }
        for (std::size_t component = 0U;
             component < 6U;
             ++component) {
            local_well_component_source[
                component] =
                well
                    .cell_source
                    .component_molar_rate_mol_per_s[
                        component];
        }
        local_well_energy_source =
            well
                .cell_source
                .energy_rate_w;
    }

    const int well_phase_reduce_error =
        MPI_Allreduce(
            local_well_phase_rate.data(),
            snapshot
                .well_phase_volumetric_rate_m3_per_s
                .data(),
            3,
            MPI_DOUBLE,
            MPI_SUM,
            comm);
    const int well_component_reduce_error =
        MPI_Allreduce(
            local_well_component_source.data(),
            snapshot
                .well_component_source_mol_per_s
                .data(),
            6,
            MPI_DOUBLE,
            MPI_SUM,
            comm);
    const int well_energy_reduce_error =
        MPI_Allreduce(
            &local_well_energy_source,
            &snapshot
                .well_energy_source_w,
            1,
            MPI_DOUBLE,
            MPI_SUM,
            comm);
    require_comm(
        comm,
        well_phase_reduce_error ==
                MPI_SUCCESS &&
            well_component_reduce_error ==
                MPI_SUCCESS &&
            well_energy_reduce_error ==
                MPI_SUCCESS,
        "failed to reduce Li-Firoozabadi two-cell well rates");

    require_comm(
        comm,
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
        "Li-Firoozabadi two-cell PETSc cleanup failed");

    return snapshot;
}

template <typename Closure>
void run_two_cell_decomposition_invariance(
    fdp::
        Pr76ThreePhaseProductionCellEvaluatorContext3D<
            Closure>* evaluator_context,
    const flow::NaturalVariableLayout3P&
        layout,
    std::span<const double> q,
    const std::vector<std::string>&
        component_ids) {
    int world_rank = -1;
    int world_size = -1;
    require(
        MPI_Comm_rank(
            PETSC_COMM_WORLD,
            &world_rank) ==
                MPI_SUCCESS &&
            MPI_Comm_size(
                PETSC_COMM_WORLD,
                &world_size) ==
                MPI_SUCCESS &&
            world_rank >= 0 &&
            world_rank < 2 &&
            world_size == 2,
        "Li-Firoozabadi decomposition-invariance regression requires the existing two-rank PETSc gate");

    const auto serial =
        solve_two_cell_decomposition_case(
            PETSC_COMM_SELF,
            false,
            evaluator_context,
            layout,
            q,
            component_ids);

    const auto distributed =
        solve_two_cell_decomposition_case(
            PETSC_COMM_WORLD,
            true,
            evaluator_context,
            layout,
            q,
            component_ids);

    require_comm(
        PETSC_COMM_WORLD,
        serial.phase_identity_keys ==
            distributed
                .phase_identity_keys &&
            serial.phase_identity_keys ==
                std::array<std::string, 3>{
                    "equilibrium-phase-0",
                    "equilibrium-phase-1",
                    "equilibrium-phase-2"},
        "serial/MPI Li-Firoozabadi stable phase identity changed");

    for (std::size_t cell = 0U;
         cell < 2U;
         ++cell) {
        require_comm(
            PETSC_COMM_WORLD,
            serial
                    .natural_variables[cell]
                    .size() ==
                q.size() &&
            distributed
                    .natural_variables[cell]
                    .size() ==
                q.size(),
            "serial/MPI Li-Firoozabadi natural-variable width changed");

        for (std::size_t slot = 0U;
             slot < q.size();
             ++slot) {
            const double absolute =
                slot == 0U
                    ? 2.0e-4
                    : 2.0e-8;
            near_comm(
                PETSC_COMM_WORLD,
                serial
                    .natural_variables[cell][slot],
                q[slot],
                3.0e-7,
                absolute,
                "serial two-cell state left the independent stationary reference");
            near_comm(
                PETSC_COMM_WORLD,
                distributed
                    .natural_variables[cell][slot],
                q[slot],
                3.0e-7,
                absolute,
                "2-rank two-cell state left the independent stationary reference");
            near_comm(
                PETSC_COMM_WORLD,
                distributed
                    .natural_variables[cell][slot],
                serial
                    .natural_variables[cell][slot],
                5.0e-9,
                absolute,
                "serial/MPI Li-Firoozabadi stable-cell state is decomposition-dependent");
        }
    }

    for (std::size_t slot = 0U;
         slot < q.size();
         ++slot) {
        const double absolute =
            slot == 0U
                ? 2.0e-4
                : 2.0e-8;
        near_comm(
            PETSC_COMM_WORLD,
            serial
                .natural_variables[0][slot],
            serial
                .natural_variables[1][slot],
            5.0e-9,
            absolute,
            "serial Li-Firoozabadi stationary symmetry broke across cells");
        near_comm(
            PETSC_COMM_WORLD,
            distributed
                .natural_variables[0][slot],
            distributed
                .natural_variables[1][slot],
            5.0e-9,
            absolute,
            "2-rank Li-Firoozabadi stationary symmetry broke across cells");
    }

    for (std::size_t phase = 0U;
         phase < 3U;
         ++phase) {
        near_comm(
            PETSC_COMM_WORLD,
            serial
                .face_phase_volumetric_rate_m3_per_s[
                    phase],
            0.0,
            0.0,
            2.0e-10,
            "serial Li-Firoozabadi internal phase flux is not zero");
        near_comm(
            PETSC_COMM_WORLD,
            distributed
                .face_phase_volumetric_rate_m3_per_s[
                    phase],
            0.0,
            0.0,
            2.0e-10,
            "2-rank Li-Firoozabadi internal phase flux is not zero");
        near_comm(
            PETSC_COMM_WORLD,
            distributed
                .face_phase_volumetric_rate_m3_per_s[
                    phase],
            serial
                .face_phase_volumetric_rate_m3_per_s[
                    phase],
            0.0,
            2.0e-10,
            "serial/MPI Li-Firoozabadi phase flux is decomposition-dependent");

        near_comm(
            PETSC_COMM_WORLD,
            serial
                .well_phase_volumetric_rate_m3_per_s[
                    phase],
            0.0,
            0.0,
            2.0e-10,
            "serial Li-Firoozabadi fixed-BHP phase rate is not zero");
        near_comm(
            PETSC_COMM_WORLD,
            distributed
                .well_phase_volumetric_rate_m3_per_s[
                    phase],
            0.0,
            0.0,
            2.0e-10,
            "2-rank Li-Firoozabadi fixed-BHP phase rate is not zero");
        near_comm(
            PETSC_COMM_WORLD,
            distributed
                .well_phase_volumetric_rate_m3_per_s[
                    phase],
            serial
                .well_phase_volumetric_rate_m3_per_s[
                    phase],
            0.0,
            2.0e-10,
            "serial/MPI Li-Firoozabadi fixed-BHP phase rate is decomposition-dependent");
    }

    for (std::size_t component = 0U;
         component < 6U;
         ++component) {
        near_comm(
            PETSC_COMM_WORLD,
            serial
                .face_component_molar_rate_mol_per_s[
                    component],
            0.0,
            0.0,
            2.0e-5,
            "serial Li-Firoozabadi internal component flux is not zero");
        near_comm(
            PETSC_COMM_WORLD,
            distributed
                .face_component_molar_rate_mol_per_s[
                    component],
            0.0,
            0.0,
            2.0e-5,
            "2-rank Li-Firoozabadi internal component flux is not zero");
        near_comm(
            PETSC_COMM_WORLD,
            distributed
                .face_component_molar_rate_mol_per_s[
                    component],
            serial
                .face_component_molar_rate_mol_per_s[
                    component],
            0.0,
            2.0e-5,
            "serial/MPI Li-Firoozabadi component flux is decomposition-dependent");

        near_comm(
            PETSC_COMM_WORLD,
            serial
                .well_component_source_mol_per_s[
                    component],
            0.0,
            0.0,
            2.0e-5,
            "serial Li-Firoozabadi well component source is not zero");
        near_comm(
            PETSC_COMM_WORLD,
            distributed
                .well_component_source_mol_per_s[
                    component],
            0.0,
            0.0,
            2.0e-5,
            "2-rank Li-Firoozabadi well component source is not zero");
        near_comm(
            PETSC_COMM_WORLD,
            distributed
                .well_component_source_mol_per_s[
                    component],
            serial
                .well_component_source_mol_per_s[
                    component],
            0.0,
            2.0e-5,
            "serial/MPI Li-Firoozabadi well component source is decomposition-dependent");

        const double expected =
            2.0 *
            ref::
                cell_component_inventory_mol[
                    component];
        near_comm(
            PETSC_COMM_WORLD,
            serial
                .global_component_inventory_mol[
                    component],
            expected,
            8.0e-7,
            1.0e-1,
            "serial Li-Firoozabadi two-cell component inventory left the independent reference");
        near_comm(
            PETSC_COMM_WORLD,
            distributed
                .global_component_inventory_mol[
                    component],
            expected,
            8.0e-7,
            1.0e-1,
            "2-rank Li-Firoozabadi two-cell component inventory left the independent reference");
        near_comm(
            PETSC_COMM_WORLD,
            distributed
                .global_component_inventory_mol[
                    component],
            serial
                .global_component_inventory_mol[
                    component],
            2.0e-8,
            1.0e-1,
            "serial/MPI Li-Firoozabadi global component inventory is decomposition-dependent");
    }

    near_comm(
        PETSC_COMM_WORLD,
        serial.face_energy_rate_w,
        0.0,
        0.0,
        20.0,
        "serial Li-Firoozabadi internal energy rate is not zero");
    near_comm(
        PETSC_COMM_WORLD,
        distributed.face_energy_rate_w,
        0.0,
        0.0,
        20.0,
        "2-rank Li-Firoozabadi internal energy rate is not zero");
    near_comm(
        PETSC_COMM_WORLD,
        distributed.face_energy_rate_w,
        serial.face_energy_rate_w,
        0.0,
        20.0,
        "serial/MPI Li-Firoozabadi energy face rate is decomposition-dependent");

    near_comm(
        PETSC_COMM_WORLD,
        serial.well_energy_source_w,
        0.0,
        0.0,
        20.0,
        "serial Li-Firoozabadi well energy source is not zero");
    near_comm(
        PETSC_COMM_WORLD,
        distributed.well_energy_source_w,
        0.0,
        0.0,
        20.0,
        "2-rank Li-Firoozabadi well energy source is not zero");
    near_comm(
        PETSC_COMM_WORLD,
        distributed.well_energy_source_w,
        serial.well_energy_source_w,
        0.0,
        20.0,
        "serial/MPI Li-Firoozabadi well energy source is decomposition-dependent");

    const double expected_energy =
        2.0 *
        ref::
            cell_total_internal_energy_j;
    near_comm(
        PETSC_COMM_WORLD,
        serial
            .global_total_internal_energy_j,
        expected_energy,
        8.0e-7,
        4.0e7,
        "serial Li-Firoozabadi two-cell energy inventory left the independent reference");
    near_comm(
        PETSC_COMM_WORLD,
        distributed
            .global_total_internal_energy_j,
        expected_energy,
        8.0e-7,
        4.0e7,
        "2-rank Li-Firoozabadi two-cell energy inventory left the independent reference");
    near_comm(
        PETSC_COMM_WORLD,
        distributed
            .global_total_internal_energy_j,
        serial
            .global_total_internal_energy_j,
        2.0e-8,
        4.0e7,
        "serial/MPI Li-Firoozabadi global energy inventory is decomposition-dependent");

    require_comm(
        PETSC_COMM_WORLD,
        serial.final_residual_l2 <=
                3.0e-8 &&
            distributed.final_residual_l2 <=
                3.0e-8,
        "serial/MPI Li-Firoozabadi final residual closure changed");
    near_comm(
        PETSC_COMM_WORLD,
        distributed.final_residual_l2,
        serial.final_residual_l2,
        0.0,
        3.0e-8,
        "serial/MPI Li-Firoozabadi final residual norm is decomposition-dependent");

    require_comm(
        PETSC_COMM_WORLD,
        serial.spatial_cross_block_abs_sum >
                0.0 &&
            distributed
                    .spatial_cross_block_abs_sum >
                0.0,
        "serial/MPI Li-Firoozabadi internal face lost spatial Jacobian coupling");
    near_comm(
        PETSC_COMM_WORLD,
        distributed
            .spatial_cross_block_abs_sum,
        serial
            .spatial_cross_block_abs_sum,
        2.0e-10,
        1.0e-12,
        "serial/MPI Li-Firoozabadi cross-cell Jacobian coupling is decomposition-dependent");
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

    run_two_cell_decomposition_invariance(
        &evaluator_context,
        layout,
        q,
        component_ids);
}
