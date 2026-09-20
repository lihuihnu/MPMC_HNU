#include <mpmc/flow_discretization_petsc/single_phase_snes_assembly.hpp>
#include <mpmc/flow_discretization_petsc/complete_natural_variable_petsc_materialization.hpp>

#include <petscsys.h>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <optional>
#include <stdexcept>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

namespace {

namespace mesh = mpmc::mesh;
namespace flow = mpmc::flow;
namespace disc = mpmc::discretization;
namespace dp = mpmc::discretization_petsc;
namespace fdp = mpmc::flow_discretization_petsc;

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
            PETSC_COMM_WORLD) != MPI_SUCCESS) {
        throw std::runtime_error(
            "MPI_Allreduce failed in single-phase SNES regression");
    }
    if (global == 0) {
        throw std::runtime_error(
            std::string{message});
    }
}

void near_collective(
    double actual,
    double expected,
    double relative = 3.0e-8,
    double absolute = 3.0e-9) {
    const bool local =
        std::isfinite(actual) &&
        std::isfinite(expected) &&
        std::abs(actual - expected) <=
            absolute +
            relative *
                std::max(
                    std::abs(actual),
                    std::abs(expected));
    require_collective(
        local,
        "single-phase distributed numeric mismatch");
}

flow::NaturalVariableLayout1P
target_layout(std::uint64_t stable_cell) {
    if (stable_cell != UINT64_C(10) &&
        stable_cell != UINT64_C(20)) {
        throw std::invalid_argument(
            "unknown single-phase controlled cell");
    }
    return flow::NaturalVariableLayout1P{
        flow::NaturalVariableCompositionPivot1P::
            from_dependent_component(
                3U,
                1U)};
}

std::vector<double>
target_state(std::uint64_t stable_cell) {
    auto layout =
        target_layout(stable_cell);
    std::vector<double> q(
        layout.unknown_count(),
        0.0);
    q[layout.pressure_unknown_index()] =
        stable_cell == UINT64_C(10)
            ? 10.0
            : 12.0;
    q[layout.temperature_unknown_index()] =
        stable_cell == UINT64_C(10)
            ? 8.0
            : 9.0;
    q[*layout.independent_composition_unknown_index(
        0U)] =
        stable_cell == UINT64_C(10)
            ? 0.20
            : 0.25;
    q[*layout.independent_composition_unknown_index(
        2U)] =
        0.30;
    return q;
}

std::vector<double>
initial_state(std::uint64_t stable_cell) {
    auto q =
        target_state(stable_cell);
    q[0] *= 1.03;
    q[1] *= 1.02;
    q[2] += 0.01;
    q[3] -= 0.005;
    return q;
}

struct ControlledAudit {
    std::uint64_t calls{};
    std::vector<double> cell10;
    std::vector<double> cell20;
};

std::vector<double>& audit_cell(
    ControlledAudit& audit,
    std::uint64_t stable) {
    if (stable == UINT64_C(10)) {
        return audit.cell10;
    }
    if (stable == UINT64_C(20)) {
        return audit.cell20;
    }
    throw std::invalid_argument(
        "unexpected single-phase stable cell");
}

PetscErrorCode evaluate_cell(
    mesh::LocalIndex,
    mesh::GlobalEntityId cell_global,
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
        static_cast<ControlledAudit*>(
            raw_context);
    ++audit->calls;

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
            state_input;
        state_input.component_ids.assign(
            component_ids.begin(),
            component_ids.end());
        state_input.reference_pressure_pa = p;
        state_input.temperature_k =
            temperature;
        state_input.composition_pivot =
            frozen_layout.composition_pivot();

        for (std::size_t rank = 0U;
             rank < 2U;
             ++rank) {
            const std::size_t component =
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
            state_input.independent_composition.push_back(
                natural_variables[*column]);
        }

        const auto target =
            target_state(cell_global.value());
        const double target_p = target[0];
        const double base_density =
            cell_global.value() == UINT64_C(10)
                ? 6.0
                : 7.0;
        const double density_slope =
            cell_global.value() == UINT64_C(10)
                ? 0.40
                : 0.35;
        const double molar_density =
            base_density +
            density_slope *
                (p - target_p);
        if (!std::isfinite(molar_density) ||
            !(molar_density > 0.0)) {
            *status =
                fdp::NaturalVariableSnesEvaluationStatus3D::
                    domain_error;
            return PETSC_SUCCESS;
        }

        state_input.phase_properties =
            flow::PhasePropertyPrerequisiteInput{
                molar_density,
                2.0,
                1.0,
                temperature,
                0.5 * temperature};

        auto state =
            flow::NaturalVariableCellState1P::
                create(std::move(state_input));

        const std::size_t q =
            frozen_layout.unknown_count();
        std::vector<double> dc(q, 0.0);
        dc[frozen_layout.pressure_unknown_index()] =
            density_slope;
        auto molar =
            flow::make_single_phase_molar_density_linearization(
                state,
                dc);

        std::vector<double> zero(q, 0.0);
        const flow::TransportPropertyProvenance
            provenance{
                "controlled-single-phase-snes",
                "software-regression",
                "v1"};
        auto transport =
            flow::make_single_phase_transport_linearization(
                state,
                zero,
                zero,
                0.0,
                zero,
                provenance,
                provenance);

        std::vector<double> dh(q, 0.0);
        std::vector<double> du(q, 0.0);
        dh[frozen_layout.temperature_unknown_index()] =
            1.0;
        du[frozen_layout.temperature_unknown_index()] =
            0.5;
        auto caloric =
            flow::make_single_phase_caloric_linearization(
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
            flow::make_single_phase_rock_thermal_storage_linearization(
                state,
                2.0 * temperature,
                rock_gradient,
                provenance);

        audit_cell(
            *audit,
            cell_global.value()) =
            std::vector<double>{
                natural_variables.begin(),
                natural_variables.end()};

        output->emplace(
            fdp::SinglePhaseCurrentCellLinearization3D{
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

mesh::Topology make_topology(int rank) {
    mesh::Topology::EntityIds ids;
    ids.faces = {
        mesh::GlobalEntityId{
            UINT64_C(100)}};
    ids.cells =
        rank == 0
            ? std::vector<mesh::GlobalEntityId>{
                  mesh::GlobalEntityId{UINT64_C(10)},
                  mesh::GlobalEntityId{UINT64_C(20)}}
            : std::vector<mesh::GlobalEntityId>{
                  mesh::GlobalEntityId{UINT64_C(20)},
                  mesh::GlobalEntityId{UINT64_C(10)}};
    return {
        std::move(ids),
        {}};
}

mesh::PartitionSnapshot
make_partition(int rank) {
    auto topology =
        make_topology(rank);
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
make_schedule(int rank) {
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
            2.0e-12};
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

mesh::DofLayout
make_dof_layout(int rank) {
    return mesh::DofLayout::create(
        make_topology(rank),
        {
            mesh::DofVariable{
                "natural_state_1p",
                mesh::EntityKind::cell,
                4U}
        });
}

mesh::DofNumberingSnapshot
make_numbering(
    int rank,
    const mesh::DofLayout& layout,
    const mesh::PartitionSnapshot& partition) {
    mesh::GlobalEntityNumberingInput input;
    input.global_cell_count = 2U;
    input.global_face_count = 1U;
    input.faces = {
        {
            mesh::GlobalEntityId{UINT64_C(100)},
            mesh::GlobalEntityOrdinal{UINT64_C(0)}}
    };
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
    return mesh::DofNumberingSnapshot::
        create_local(
            layout,
            partition,
            std::move(input));
}

dp::PetscMpiAijSymbolicPreallocation3D
make_cell_bridge(int rank) {
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
make_cell_pattern(int rank) {
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
direct_admissibility() {
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
materialized_entry() {
    return {
        mesh::LocalIndex{0U},
        disc::TpfaInternalFaceTransmissibilityDisposition3D::
            materialized,
        direct_admissibility(),
        disc::TpfaStaticFaceTransmissibility3D{
            disc::TpfaStaticFaceTransmissibilityDisposition3D::
                positive_harmonic_combination,
            1.5,
            2.0e-12}};
}

std::optional<
    fdp::SinglePhaseCurrentCellLinearization3D>
evaluate_direct(
    std::uint64_t stable_cell,
    const std::vector<double>& q,
    ControlledAudit* audit) {
    std::optional<
        fdp::SinglePhaseCurrentCellLinearization3D>
        output;
    auto layout =
        target_layout(stable_cell);
    const std::vector<std::string>
        ids{"A", "B", "C"};
    fdp::NaturalVariableSnesEvaluationStatus3D
        status =
            fdp::NaturalVariableSnesEvaluationStatus3D::
                success;
    const PetscErrorCode error =
        evaluate_cell(
            mesh::LocalIndex{0U},
            mesh::GlobalEntityId{stable_cell},
            q,
            layout,
            ids,
            audit,
            &output,
            &status);
    if (error != PETSC_SUCCESS ||
        status !=
            fdp::NaturalVariableSnesEvaluationStatus3D::
                success ||
        !output.has_value()) {
        throw std::runtime_error(
            "single-phase controlled direct evaluation failed");
    }
    return output;
}

std::vector<fdp::SinglePhaseSnesCellInput3D>
make_cells(
    int rank,
    ControlledAudit* audit) {
    std::vector<fdp::SinglePhaseSnesCellInput3D>
        result;
    result.reserve(2U);

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
        auto target =
            target_state(stable);
        auto evaluation =
            evaluate_direct(
                stable,
                target,
                audit);

        fdp::SinglePhaseSnesCellInput3D input;
        input.cell =
            mesh::LocalIndex{
                static_cast<
                    mesh::LocalIndex::value_type>(
                        local)};
        input.cell_global =
            mesh::GlobalEntityId{stable};
        input.bulk_volume_m3 =
            stable == UINT64_C(10)
                ? 2.0
                : 5.0;
        input.porosity =
            stable == UINT64_C(10)
                ? 0.25
                : 0.30;
        input.frozen_layout =
            target_layout(stable);
        input.component_ids =
            {"A", "B", "C"};

        if (owned) {
            input.previous_component_accumulation =
                flow::build_single_phase_component_accumulation(
                    evaluation->state,
                    input.porosity);
            input.previous_energy_accumulation =
                flow::build_single_phase_energy_accumulation_snapshot(
                    evaluation->state,
                    input.porosity,
                    evaluation->transport,
                    evaluation->caloric,
                    evaluation->rock);
        }
        result.push_back(
            std::move(input));
    }
    return result;
}

std::vector<
    fdp::SinglePhaseSnesAuthoritativeFaceInput3D>
make_faces(int rank) {
    if (rank != 0) {
        return {};
    }
    return {
        {
            mesh::LocalIndex{0U},
            mesh::GlobalEntityId{UINT64_C(100)},
            materialized_entry(),
            flow::GravityVector3D{
                0.0, 0.0, -9.81},
            flow::OwnerToNeighbourDisplacement3D{
                0.0, 0.0, 1.0},
            mpmc::flow_discretization::
                StaticThermalFaceConductance3D{
                    0.0}}
    };
}

void set_owned_state(
    Vec state,
    int rank,
    bool use_target) {
    const std::uint64_t stable =
        rank == 0
            ? UINT64_C(10)
            : UINT64_C(20);
    const auto values =
        use_target
            ? target_state(stable)
            : initial_state(stable);
    const PetscInt start =
        static_cast<PetscInt>(
            rank * 4);
    for (std::size_t slot = 0U;
         slot < values.size();
         ++slot) {
        const PetscInt index =
            start +
            static_cast<PetscInt>(slot);
        const PetscScalar value =
            static_cast<PetscScalar>(
                values[slot]);
        require_collective(
            VecSetValues(
                state,
                1,
                &index,
                &value,
                INSERT_VALUES) ==
                PETSC_SUCCESS,
            "failed to insert single-phase owned state");
    }
    require_collective(
        VecAssemblyBegin(state) ==
                PETSC_SUCCESS &&
            VecAssemblyEnd(state) ==
                PETSC_SUCCESS,
        "failed to assemble single-phase owned state");
}

} // namespace

void single_phase_snes_physical_assembly_test() {
    int rank = -1;
    int size = -1;
    if (MPI_Comm_rank(
            PETSC_COMM_WORLD,
            &rank) != MPI_SUCCESS ||
        MPI_Comm_size(
            PETSC_COMM_WORLD,
            &size) != MPI_SUCCESS) {
        throw std::runtime_error(
            "MPI rank/size query failed");
    }
    require_collective(
        size == 2,
        "single-phase SNES test requires two ranks");

    auto partition =
        make_partition(rank);
    auto schedule =
        make_schedule(rank);
    auto dof_layout =
        make_dof_layout(rank);
    auto dof_numbering =
        make_numbering(
            rank,
            dof_layout,
            partition);
    auto cell_bridge =
        make_cell_bridge(rank);
    auto cell_pattern =
        make_cell_pattern(rank);

    ControlledAudit audit;
    auto cells =
        make_cells(
            rank,
            &audit);
    auto faces =
        make_faces(rank);

    std::optional<
        fdp::SinglePhaseSnesAssemblyContext3D>
        context;
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
                "natural_state_1p",
                1.0,
                std::move(cells),
                std::move(faces),
                {evaluate_cell, &audit},
                &context);
    require_collective(
        error == PETSC_SUCCESS &&
            context.has_value(),
        "failed to create single-phase SNES assembly context");

    Vec initial = nullptr;
    require_collective(
        VecCreateMPI(
            PETSC_COMM_WORLD,
            4,
            8,
            &initial) ==
            PETSC_SUCCESS,
        "failed to create single-phase initial Vec");
    set_owned_state(
        initial,
        rank,
        false);

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
    require_collective(
        error == PETSC_SUCCESS &&
            initial_status ==
                fdp::NaturalVariableSnesEvaluationStatus3D::
                    success &&
            initial_assembly.has_value() &&
            initial_assembly->component_count() ==
                3U &&
            initial_assembly
                    ->natural_variable_count() ==
                4U &&
            initial_assembly->residual_entries()
                    .size() ==
                4U,
        "single-phase complete assembly did not reduce to Nc+1");

    for (const auto& entry :
         initial_assembly->residual_entries()) {
        require_collective(
            entry.equation_kind !=
                fdp::NaturalVariableEquationKind3D::
                    fugacity_equilibrium &&
                entry.equation_slot <= 3U,
            "single-phase assembly contains a fugacity or inactive-phase row");
    }

    const auto remote =
        initial_state(
            rank == 0
                ? UINT64_C(20)
                : UINT64_C(10));
    const auto& observed_remote =
        rank == 0
            ? audit.cell20
            : audit.cell10;
    require_collective(
        observed_remote.size() ==
            remote.size(),
        "single-phase PetscSF ghost state was not evaluated");
    for (std::size_t slot = 0U;
         slot < remote.size();
         ++slot) {
        near_collective(
            observed_remote[slot],
            remote[slot],
            0.0,
            0.0);
    }

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
    require_collective(
        error == PETSC_SUCCESS &&
            unused_residual != nullptr &&
            jacobian_template != nullptr,
        "failed to materialize single-phase PETSc system");

    const std::uint64_t calls_before =
        audit.calls;
    Vec solution = nullptr;
    std::optional<
        fdp::NaturalVariableSnesSolveReport3D>
        report;
    error =
        fdp::solve_natural_variable_snes_3d(
            PETSC_COMM_WORLD,
            *initial_assembly,
            initial,
            jacobian_template,
            context->snes_evaluator(),
            &solution,
            &report);
    require_collective(
        error == PETSC_SUCCESS &&
            solution != nullptr &&
            report.has_value() &&
            report->natural_variable_count() ==
                4U &&
            static_cast<int>(
                report->converged_reason()) >
                0 &&
            report->function_domain_errors() ==
                0 &&
            report->jacobian_domain_errors() ==
                0 &&
            audit.calls > calls_before,
        "PETSc SNES failed on reduced single-phase production assembly");

    const auto target =
        target_state(
            rank == 0
                ? UINT64_C(10)
                : UINT64_C(20));
    const PetscInt start =
        static_cast<PetscInt>(
            rank * 4);
    for (std::size_t slot = 0U;
         slot < target.size();
         ++slot) {
        const PetscInt index =
            start +
            static_cast<PetscInt>(slot);
        PetscScalar value{};
        require_collective(
            VecGetValues(
                solution,
                1,
                &index,
                &value) ==
                PETSC_SUCCESS,
            "failed to read converged single-phase state");
        near_collective(
            static_cast<double>(
                PetscRealPart(value)),
            target[slot]);
    }

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
    require_collective(
        error == PETSC_SUCCESS &&
            final_status ==
                fdp::NaturalVariableSnesEvaluationStatus3D::
                    success &&
            final_assembly.has_value() &&
            final_assembly
                    ->natural_variable_count() ==
                4U,
        "independent final single-phase assembly failed");

    double local_squared = 0.0;
    for (const auto& entry :
         final_assembly->residual_entries()) {
        local_squared +=
            entry.native_value *
            entry.native_value;
    }
    double global_squared = 0.0;
    if (MPI_Allreduce(
            &local_squared,
            &global_squared,
            1,
            MPI_DOUBLE,
            MPI_SUM,
            PETSC_COMM_WORLD) !=
        MPI_SUCCESS) {
        throw std::runtime_error(
            "failed to reduce single-phase final residual");
    }
    require_collective(
        std::sqrt(global_squared) <=
            3.0e-8 &&
            report->final_function_l2_norm() <=
                3.0e-8,
        "single-phase production assembly did not converge to target");

    for (const auto& entry :
         final_assembly->jacobian_entries()) {
        require_collective(
            std::isfinite(entry.value) &&
                entry.natural_variable_column <
                    4U,
            "single-phase final Jacobian is invalid");
    }

    const PetscErrorCode solution_destroy =
        VecDestroy(&solution);
    const PetscErrorCode residual_destroy =
        VecDestroy(&unused_residual);
    const PetscErrorCode matrix_destroy =
        MatDestroy(&jacobian_template);
    const PetscErrorCode initial_destroy =
        VecDestroy(&initial);
    require_collective(
        solution_destroy == PETSC_SUCCESS &&
            residual_destroy == PETSC_SUCCESS &&
            matrix_destroy == PETSC_SUCCESS &&
            initial_destroy == PETSC_SUCCESS,
        "single-phase SNES fixture cleanup failed");
}
