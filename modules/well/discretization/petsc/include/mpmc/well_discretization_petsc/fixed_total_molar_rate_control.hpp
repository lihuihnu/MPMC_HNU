#ifndef MPMC_WELL_DISCRETIZATION_PETSC_FIXED_TOTAL_MOLAR_RATE_CONTROL_HPP
#define MPMC_WELL_DISCRETIZATION_PETSC_FIXED_TOTAL_MOLAR_RATE_CONTROL_HPP

#include <mpmc/flow_discretization_petsc/natural_variable_snes_solver.hpp>
#include <mpmc/flow_discretization_petsc/phase_transition_outer_rebuild.hpp>
#include <mpmc/well_discretization_petsc/fixed_bhp_well_source_evaluator.hpp>

#include <petscksp.h>
#include <petscmat.h>
#include <petscsnes.h>
#include <petscvec.h>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <memory>
#include <optional>
#include <span>
#include <stdexcept>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

namespace mpmc::well_discretization_petsc {

inline constexpr std::string_view
    fixed_total_molar_rate_control_convention =
        "well-discretization-petsc/single-well/fixed-total-molar-rate-control/v1";

struct FixedTotalMolarRateAuthoritativeConnectionEvaluation3D {
    mpmc::mesh::LocalIndex cell{
        mpmc::mesh::LocalIndex::value_type{0}};
    mpmc::mesh::GlobalEntityId cell_global{
        mpmc::mesh::GlobalEntityId::value_type{0}};
    mpmc::well_discretization::
        FixedBhpConnectionCellSourceLinearization3D
            linearization;
};

/// Dynamic source context used only while one rate-controlled SNES evaluation
/// is in progress.
///
/// The logical completion set/geometry is immutable. Only the scalar BHP value
/// changes between Newton evaluations. The mixed-cardinality reservoir
/// assembly already calls cell-source evaluators only on authoritative owners,
/// so this cache contains each physical connection at most once per evaluation.
class FixedTotalMolarRateWellSourceEvaluatorContext3D {
public:
    static constexpr std::string_view convention =
        fixed_total_molar_rate_control_convention;

    [[nodiscard]] static
    FixedTotalMolarRateWellSourceEvaluatorContext3D
    create(
        FixedBhpMultiConnectionWellSourceEvaluatorContext3D
            well,
        double initial_bottom_hole_pressure_pa) {
        if (!std::isfinite(
                initial_bottom_hole_pressure_pa) ||
            !(initial_bottom_hole_pressure_pa > 0.0)) {
            throw std::invalid_argument(
                "mpmc::well_discretization_petsc: rate-controlled initial BHP must be finite and positive");
        }
        return FixedTotalMolarRateWellSourceEvaluatorContext3D{
            std::move(well),
            initial_bottom_hole_pressure_pa};
    }

    [[nodiscard]] const
    FixedBhpMultiConnectionWellSourceEvaluatorContext3D&
    well() const noexcept {
        return well_;
    }

    [[nodiscard]] double
    current_bottom_hole_pressure_pa() const noexcept {
        return current_bottom_hole_pressure_pa_;
    }

    void begin_evaluation(
        double bottom_hole_pressure_pa) {
        if (!std::isfinite(bottom_hole_pressure_pa) ||
            !(bottom_hole_pressure_pa > 0.0)) {
            throw std::invalid_argument(
                "mpmc::well_discretization_petsc: rate-controlled BHP must be finite and positive");
        }
        current_bottom_hole_pressure_pa_ =
            bottom_hole_pressure_pa;
        local_authoritative_evaluations_.clear();
    }

    [[nodiscard]] std::span<
        const FixedTotalMolarRateAuthoritativeConnectionEvaluation3D>
    local_authoritative_evaluations() const noexcept {
        return local_authoritative_evaluations_;
    }

private:
    FixedTotalMolarRateWellSourceEvaluatorContext3D(
        FixedBhpMultiConnectionWellSourceEvaluatorContext3D
            well,
        double initial_bottom_hole_pressure_pa)
        : well_(std::move(well)),
          current_bottom_hole_pressure_pa_(
              initial_bottom_hole_pressure_pa) {}

    friend PetscErrorCode
    evaluate_fixed_total_molar_rate_well_source_3d(
        mpmc::mesh::LocalIndex,
        mpmc::mesh::GlobalEntityId,
        std::span<const double>,
        const mpmc::flow_discretization_petsc::
            MixedCardinalityPhysicalCurrentCellLinearization3D&,
        void*,
        std::optional<
            mpmc::flow_discretization::
                CellSourceLinearization3D>*,
        mpmc::flow_discretization_petsc::
            NaturalVariableSnesEvaluationStatus3D*);

    FixedBhpMultiConnectionWellSourceEvaluatorContext3D
        well_;
    double current_bottom_hole_pressure_pa_{};
    std::vector<
        FixedTotalMolarRateAuthoritativeConnectionEvaluation3D>
        local_authoritative_evaluations_;
};

inline PetscErrorCode
evaluate_fixed_total_molar_rate_well_source_3d(
    mpmc::mesh::LocalIndex cell,
    mpmc::mesh::GlobalEntityId cell_global,
    std::span<const double> natural_variables,
    const mpmc::flow_discretization_petsc::
        MixedCardinalityPhysicalCurrentCellLinearization3D&
            current,
    void* raw_context,
    std::optional<
        mpmc::flow_discretization::
            CellSourceLinearization3D>* output,
    mpmc::flow_discretization_petsc::
        NaturalVariableSnesEvaluationStatus3D*
            status) {
    using namespace
        mpmc::flow_discretization_petsc;

    if (raw_context == nullptr ||
        output == nullptr ||
        status == nullptr) {
        return PETSC_ERR_ARG_NULL;
    }
    output->reset();
    *status =
        NaturalVariableSnesEvaluationStatus3D::
            success;

    auto* context =
        static_cast<
            FixedTotalMolarRateWellSourceEvaluatorContext3D*>(
                raw_context);
    const auto* frozen_connection =
        context->well_.find_connection(
            cell_global);
    if (frozen_connection == nullptr) {
        return PETSC_SUCCESS;
    }

    try {
        auto dynamic_connection =
            frozen_connection
                ->with_bottom_hole_pressure(
                    context
                        ->current_bottom_hole_pressure_pa_);
        auto evaluated =
            build_fixed_bhp_peaceman_well_source_3d(
                dynamic_connection,
                current);

        const auto duplicate =
            std::find_if(
                context
                    ->local_authoritative_evaluations_
                    .begin(),
                context
                    ->local_authoritative_evaluations_
                    .end(),
                [&](const auto& entry) {
                    return entry.cell_global ==
                        cell_global;
                });
        if (duplicate !=
            context
                ->local_authoritative_evaluations_
                .end()) {
            return PETSC_ERR_PLIB;
        }

        output->emplace(
            evaluated.cell_source);
        context
            ->local_authoritative_evaluations_
            .push_back(
                {
                    cell,
                    cell_global,
                    std::move(evaluated)});
        return PETSC_SUCCESS;
    } catch (const std::invalid_argument&) {
        *status =
            NaturalVariableSnesEvaluationStatus3D::
                domain_error;
        return PETSC_SUCCESS;
    } catch (const std::range_error&) {
        *status =
            NaturalVariableSnesEvaluationStatus3D::
                domain_error;
        return PETSC_SUCCESS;
    } catch (...) {
        return PETSC_ERR_ARG_INCOMP;
    }
}

[[nodiscard]] inline
mpmc::flow_discretization_petsc::
    MixedCardinalityPhysicalCellSourceEvaluatorBinding3D
fixed_total_molar_rate_well_source_binding_3d(
    FixedTotalMolarRateWellSourceEvaluatorContext3D*
        context) noexcept {
    return {
        &evaluate_fixed_total_molar_rate_well_source_3d,
        context};
}

struct FixedTotalMolarRateWellControlSolveReport3D {
    SNESConvergedReason converged_reason{
        SNES_CONVERGED_ITERATING};
    PetscInt nonlinear_iterations{};
    PetscInt function_evaluations{};
    PetscInt jacobian_evaluations{};
    PetscInt global_scalar_count{};
    PetscInt well_global_scalar{-1};
    int well_owner_rank{-1};
    double bottom_hole_pressure_pa{};
    double target_total_molar_rate_mol_per_s{};
    double achieved_total_molar_rate_mol_per_s{};
    double final_function_l2_norm{};
    std::string snes_type;
    std::string ksp_type;
    std::string pc_type;

    [[nodiscard]] bool converged() const noexcept {
        return static_cast<int>(
                   converged_reason) >
            0;
    }

    [[nodiscard]] double
    total_molar_rate_residual_mol_per_s() const noexcept {
        return
            achieved_total_molar_rate_mol_per_s -
            target_total_molar_rate_mol_per_s;
    }
};

/// One-scalar augmentation of an existing frozen reservoir nonlinear system.
///
/// Reservoir scalar indices [0,N) remain unchanged. The unique BHP scalar is
/// appended at global index N and owned by the last MPI rank. This preserves
/// every existing reservoir row/column mapping while adding:
///
///   Jrr : existing reservoir Jacobian
///   Jrw : reservoir source residual derivative with respect to BHP
///   Jwr : total-molar well-rate derivative with respect to reservoir variables
///   Jww : total-molar well-rate derivative with respect to BHP
///
/// The control equation is production-positive:
///
///   R_w = sum_connection sum_component n_dot_i - target = 0.
///
/// Phase cardinality is frozen for the lifetime of this augmented system.
class FixedTotalMolarRateWellControlSystem3D {
public:
    static constexpr std::string_view convention =
        fixed_total_molar_rate_control_convention;

    FixedTotalMolarRateWellControlSystem3D(
        const FixedTotalMolarRateWellControlSystem3D&) =
        delete;
    FixedTotalMolarRateWellControlSystem3D& operator=(
        const FixedTotalMolarRateWellControlSystem3D&) =
        delete;
    FixedTotalMolarRateWellControlSystem3D(
        FixedTotalMolarRateWellControlSystem3D&&) =
        delete;
    FixedTotalMolarRateWellControlSystem3D& operator=(
        FixedTotalMolarRateWellControlSystem3D&&) =
        delete;

    ~FixedTotalMolarRateWellControlSystem3D() {
        if (row_scaling_ != nullptr) {
            (void)VecDestroy(&row_scaling_);
        }
        if (jacobian_ != nullptr) {
            (void)MatDestroy(&jacobian_);
        }
        if (initial_state_ != nullptr) {
            (void)VecDestroy(&initial_state_);
        }
        if (reservoir_direction_ != nullptr) {
            (void)VecDestroy(&reservoir_direction_);
        }
        if (reservoir_residual_ != nullptr) {
            (void)VecDestroy(&reservoir_residual_);
        }
        if (reservoir_state_ != nullptr) {
            (void)VecDestroy(&reservoir_state_);
        }
        if (reservoir_jacobian_ != nullptr) {
            (void)MatDestroy(&reservoir_jacobian_);
        }
    }

    [[nodiscard]] static PetscErrorCode
    create(
        MPI_Comm comm,
        mpmc::flow_discretization_petsc::
            PhaseTransitionRebuiltNaturalVariableSystem3D*
                reservoir_system,
        FixedTotalMolarRateWellSourceEvaluatorContext3D*
            source_context,
        double initial_bottom_hole_pressure_pa,
        double target_total_molar_rate_mol_per_s,
        std::unique_ptr<
            FixedTotalMolarRateWellControlSystem3D>*
                output) {
        using namespace
            mpmc::flow_discretization_petsc;

        if (reservoir_system == nullptr ||
            source_context == nullptr ||
            output == nullptr) {
            return PETSC_ERR_ARG_NULL;
        }
        output->reset();
        if (!std::isfinite(
                initial_bottom_hole_pressure_pa) ||
            !(initial_bottom_hole_pressure_pa > 0.0) ||
            !std::isfinite(
                target_total_molar_rate_mol_per_s) ||
            !(target_total_molar_rate_mol_per_s > 0.0) ||
            source_context->well()
                    .connection_count() <
                2U) {
            return PETSC_ERR_ARG_INCOMP;
        }

        int mpi_rank = -1;
        int mpi_size = -1;
        if (MPI_Comm_rank(
                comm,
                &mpi_rank) != MPI_SUCCESS ||
            MPI_Comm_size(
                comm,
                &mpi_size) != MPI_SUCCESS ||
            mpi_rank < 0 ||
            mpi_size <= 0) {
            return PETSC_ERR_MPI;
        }

        const auto& numbering =
            reservoir_system->numbering();
        if (numbering.local_rank().value() !=
                static_cast<std::uint32_t>(
                    mpi_rank) ||
            numbering.rank_count() !=
                static_cast<std::uint32_t>(
                    mpi_size) ||
            numbering.petsc_global_scalar_count() <=
                0) {
            return PETSC_ERR_ARG_INCOMP;
        }

        auto system =
            std::unique_ptr<
                FixedTotalMolarRateWellControlSystem3D>(
                    new FixedTotalMolarRateWellControlSystem3D(
                        comm,
                        reservoir_system,
                        source_context,
                        target_total_molar_rate_mol_per_s,
                        mpi_rank,
                        mpi_size));

        PetscErrorCode error =
            system->initialize(
                initial_bottom_hole_pressure_pa);
        if (error != PETSC_SUCCESS) {
            return error;
        }

        *output =
            std::move(system);
        return PETSC_SUCCESS;
    }

    [[nodiscard]] PetscInt
    reservoir_global_scalar_count() const noexcept {
        return reservoir_global_scalar_count_;
    }

    [[nodiscard]] PetscInt
    global_scalar_count() const noexcept {
        return reservoir_global_scalar_count_ +
            PetscInt{1};
    }

    [[nodiscard]] PetscInt
    well_global_scalar() const noexcept {
        return reservoir_global_scalar_count_;
    }

    [[nodiscard]] int
    well_owner_rank() const noexcept {
        return well_owner_rank_;
    }

    [[nodiscard]] Vec
    initial_state() const noexcept {
        return initial_state_;
    }

    [[nodiscard]] Mat
    jacobian_structure() const noexcept {
        return jacobian_;
    }

    [[nodiscard]]
    mpmc::flow_discretization_petsc::
        NaturalVariableSnesEvaluator3D
    snes_evaluator() noexcept {
        return {
            &snes_function,
            &snes_jacobian,
            &snes_precheck,
            this};
    }

    [[nodiscard]] PetscErrorCode
    copy_reservoir_state(
        Vec augmented_state,
        Vec* output) {
        if (augmented_state == nullptr ||
            output == nullptr) {
            return PETSC_ERR_ARG_NULL;
        }
        if (*output != nullptr) {
            return PETSC_ERR_ARG_WRONGSTATE;
        }
        PetscErrorCode error =
            sync_augmented_state(
                augmented_state,
                nullptr);
        if (error != PETSC_SUCCESS) {
            return error;
        }
        error =
            VecDuplicate(
                reservoir_state_,
                output);
        if (error == PETSC_SUCCESS) {
            error =
                VecCopy(
                    reservoir_state_,
                    *output);
        }
        return error;
    }

    [[nodiscard]] PetscErrorCode
    evaluate_whole_well_total_molar_rate(
        Vec augmented_state,
        double* total_molar_rate_mol_per_s) {
        using namespace
            mpmc::flow_discretization_petsc;
        if (augmented_state == nullptr ||
            total_molar_rate_mol_per_s ==
                nullptr) {
            return PETSC_ERR_ARG_NULL;
        }
        *total_molar_rate_mol_per_s =
            0.0;

        double bottom_hole_pressure_pa =
            0.0;
        PetscErrorCode error =
            sync_augmented_state(
                augmented_state,
                &bottom_hole_pressure_pa);
        if (error != PETSC_SUCCESS) {
            return error;
        }

        source_context_->begin_evaluation(
            bottom_hole_pressure_pa);
        NaturalVariableSnesEvaluationStatus3D
            status =
                NaturalVariableSnesEvaluationStatus3D::
                    success;
        error =
            VecSet(
                reservoir_residual_,
                PetscScalar{0.0});
        if (error == PETSC_SUCCESS) {
            error =
                reservoir_evaluator_.function(
                    reservoir_state_,
                    reservoir_residual_,
                    reservoir_evaluator_
                        .user_context,
                    &status);
        }
        if (error != PETSC_SUCCESS) {
            return error;
        }
        if (status !=
            NaturalVariableSnesEvaluationStatus3D::
                success) {
            return PETSC_ERR_ARG_OUTOFRANGE;
        }

        double local_total =
            local_total_molar_production_rate();
        double global_total = 0.0;
        if (MPI_Allreduce(
                &local_total,
                &global_total,
                1,
                MPI_DOUBLE,
                MPI_SUM,
                comm_) != MPI_SUCCESS) {
            return PETSC_ERR_MPI;
        }
        *total_molar_rate_mol_per_s =
            global_total;
        return PETSC_SUCCESS;
    }

    [[nodiscard]] PetscErrorCode
    solve(
        Vec* solution,
        std::optional<
            FixedTotalMolarRateWellControlSolveReport3D>*
                report) {
        using namespace
            mpmc::flow_discretization_petsc;
        using namespace
            mpmc::flow_discretization_petsc::
                natural_variable_snes_detail;

        if (solution == nullptr ||
            report == nullptr) {
            return PETSC_ERR_ARG_NULL;
        }
        if (*solution != nullptr) {
            return PETSC_ERR_ARG_WRONGSTATE;
        }
        report->reset();

        Vec solved_state = nullptr;
        Vec residual = nullptr;
        Mat jacobian = nullptr;
        SNES snes = nullptr;

        PetscErrorCode error =
            VecDuplicate(
                initial_state_,
                &solved_state);
        if (error == PETSC_SUCCESS) {
            error =
                VecCopy(
                    initial_state_,
                    solved_state);
        }
        if (error == PETSC_SUCCESS) {
            error =
                VecDuplicate(
                    initial_state_,
                    &residual);
        }
        if (error == PETSC_SUCCESS) {
            error =
                MatDuplicate(
                    jacobian_,
                    MAT_DO_NOT_COPY_VALUES,
                    &jacobian);
        }
        if (error == PETSC_SUCCESS) {
            error =
                SNESCreate(
                    comm_,
                    &snes);
        }
        if (error != PETSC_SUCCESS) {
            cleanup_solve(
                &solved_state,
                &residual,
                &jacobian,
                &snes);
            return error;
        }

        NaturalVariableSnesEvaluator3D evaluator =
            snes_evaluator();
        CallbackContext callback_context{
            evaluator,
            row_scaling_};

        error =
            SNESSetFunction(
                snes,
                residual,
                form_function,
                &callback_context);
        if (error == PETSC_SUCCESS) {
            error =
                SNESSetJacobian(
                    snes,
                    jacobian,
                    jacobian,
                    form_jacobian,
                    &callback_context);
        }
        if (error == PETSC_SUCCESS) {
            error =
                SNESSetType(
                    snes,
                    SNESNEWTONLS);
        }

        SNESLineSearch line_search = nullptr;
        if (error == PETSC_SUCCESS) {
            error =
                SNESGetLineSearch(
                    snes,
                    &line_search);
        }
        if (error == PETSC_SUCCESS) {
            error =
                SNESLineSearchSetType(
                    line_search,
                    SNESLINESEARCHBT);
        }
        if (error == PETSC_SUCCESS) {
            error =
                SNESLineSearchSetPreCheck(
                    line_search,
                    line_search_precheck,
                    &callback_context);
        }

        KSP ksp = nullptr;
        PC pc = nullptr;
        if (error == PETSC_SUCCESS) {
            error =
                SNESGetKSP(
                    snes,
                    &ksp);
        }
        if (error == PETSC_SUCCESS) {
            error =
                KSPSetType(
                    ksp,
                    KSPGMRES);
        }
        if (error == PETSC_SUCCESS) {
            error =
                KSPGetPC(
                    ksp,
                    &pc);
        }
        if (error == PETSC_SUCCESS) {
            error =
                PCSetType(
                    pc,
                    PCASM);
        }
        if (error == PETSC_SUCCESS) {
            error =
                PCASMSetType(
                    pc,
                    PC_ASM_RESTRICT);
        }
        if (error == PETSC_SUCCESS) {
            error =
                PCASMSetOverlap(
                    pc,
                    PetscInt{1});
        }

        constexpr const char* prefix =
            "mpmc_rate_control_";
        constexpr const char* sub_pc =
            "-mpmc_rate_control_sub_pc_type";
        constexpr const char* sub_reorder =
            "-mpmc_rate_control_sub_pc_factor_nonzeros_along_diagonal";
        constexpr const char* sub_shift =
            "-mpmc_rate_control_sub_pc_factor_shift_type";
        bool sub_pc_set = false;
        bool sub_reorder_set = false;
        bool sub_shift_set = false;

        if (error == PETSC_SUCCESS) {
            error =
                PCSetOptionsPrefix(
                    pc,
                    prefix);
        }
        if (error == PETSC_SUCCESS) {
            error =
                PetscOptionsSetValue(
                    nullptr,
                    sub_pc,
                    "lu");
            sub_pc_set =
                error == PETSC_SUCCESS;
        }
        if (error == PETSC_SUCCESS) {
            error =
                PetscOptionsSetValue(
                    nullptr,
                    sub_reorder,
                    "1.0e-10");
            sub_reorder_set =
                error == PETSC_SUCCESS;
        }
        if (error == PETSC_SUCCESS) {
            error =
                PetscOptionsSetValue(
                    nullptr,
                    sub_shift,
                    "nonzero");
            sub_shift_set =
                error == PETSC_SUCCESS;
        }
        if (error == PETSC_SUCCESS) {
            error =
                KSPSetTolerances(
                    ksp,
                    PetscReal{1.0e-12},
                    PetscReal{1.0e-14},
                    PETSC_DEFAULT,
                    PetscInt{200});
        }
        if (error == PETSC_SUCCESS) {
            error =
                SNESSetTolerances(
                    snes,
                    PetscReal{1.0e-10},
                    PetscReal{1.0e-10},
                    PetscReal{1.0e-12},
                    PetscInt{30},
                    PetscInt{2000});
        }
        if (error == PETSC_SUCCESS) {
            error =
                SNESSetErrorIfNotConverged(
                    snes,
                    PETSC_FALSE);
        }
        if (error == PETSC_SUCCESS) {
            error =
                SNESSetCheckJacobianDomainError(
                    snes,
                    PETSC_TRUE);
        }

        if (error == PETSC_SUCCESS) {
            error =
                SNESSolve(
                    snes,
                    nullptr,
                    solved_state);
        }

        PetscErrorCode clear_error =
            PETSC_SUCCESS;
        if (sub_shift_set) {
            clear_error =
                PetscOptionsClearValue(
                    nullptr,
                    sub_shift);
        }
        if (sub_reorder_set) {
            const PetscErrorCode candidate =
                PetscOptionsClearValue(
                    nullptr,
                    sub_reorder);
            if (clear_error ==
                    PETSC_SUCCESS &&
                candidate !=
                    PETSC_SUCCESS) {
                clear_error =
                    candidate;
            }
        }
        if (sub_pc_set) {
            const PetscErrorCode candidate =
                PetscOptionsClearValue(
                    nullptr,
                    sub_pc);
            if (clear_error ==
                    PETSC_SUCCESS &&
                candidate !=
                    PETSC_SUCCESS) {
                clear_error =
                    candidate;
            }
        }
        if (error == PETSC_SUCCESS &&
            clear_error !=
                PETSC_SUCCESS) {
            error =
                clear_error;
        }
        if (error != PETSC_SUCCESS) {
            cleanup_solve(
                &solved_state,
                &residual,
                &jacobian,
                &snes);
            return error;
        }

        SNESConvergedReason reason{
            SNES_CONVERGED_ITERATING};
        PetscInt nonlinear_iterations = -1;
        PetscReal function_norm = 0.0;
        error =
            SNESGetConvergedReason(
                snes,
                &reason);
        if (error == PETSC_SUCCESS) {
            error =
                SNESGetIterationNumber(
                    snes,
                    &nonlinear_iterations);
        }
        if (error == PETSC_SUCCESS) {
            error =
                SNESGetFunctionNorm(
                    snes,
                    &function_norm);
        }
        if (error != PETSC_SUCCESS ||
            static_cast<int>(reason) <= 0) {
            cleanup_solve(
                &solved_state,
                &residual,
                &jacobian,
                &snes);
            return error != PETSC_SUCCESS
                ? error
                : PETSC_ERR_NOT_CONVERGED;
        }

        double bottom_hole_pressure_pa =
            0.0;
        error =
            extract_bottom_hole_pressure(
                solved_state,
                &bottom_hole_pressure_pa);
        double achieved_rate = 0.0;
        if (error == PETSC_SUCCESS) {
            error =
                evaluate_whole_well_total_molar_rate(
                    solved_state,
                    &achieved_rate);
        }
        if (error != PETSC_SUCCESS) {
            cleanup_solve(
                &solved_state,
                &residual,
                &jacobian,
                &snes);
            return error;
        }

        report->emplace(
            FixedTotalMolarRateWellControlSolveReport3D{
                reason,
                nonlinear_iterations,
                callback_context
                    .function_evaluations,
                callback_context
                    .jacobian_evaluations,
                global_scalar_count(),
                well_global_scalar(),
                well_owner_rank_,
                bottom_hole_pressure_pa,
                target_total_molar_rate_mol_per_s_,
                achieved_rate,
                static_cast<double>(
                    function_norm),
                SNESNEWTONLS,
                KSPGMRES,
                PCASM});

        *solution =
            solved_state;
        solved_state = nullptr;
        cleanup_solve(
            &solved_state,
            &residual,
            &jacobian,
            &snes);
        return PETSC_SUCCESS;
    }

private:
    FixedTotalMolarRateWellControlSystem3D(
        MPI_Comm comm,
        mpmc::flow_discretization_petsc::
            PhaseTransitionRebuiltNaturalVariableSystem3D*
                reservoir_system,
        FixedTotalMolarRateWellSourceEvaluatorContext3D*
            source_context,
        double target_total_molar_rate_mol_per_s,
        int mpi_rank,
        int mpi_size)
        : comm_(comm),
          reservoir_system_(reservoir_system),
          source_context_(source_context),
          target_total_molar_rate_mol_per_s_(
              target_total_molar_rate_mol_per_s),
          mpi_rank_(mpi_rank),
          mpi_size_(mpi_size),
          well_owner_rank_(
              mpi_size - 1),
          reservoir_evaluator_(
              reservoir_system
                  ->snes_evaluator()),
          reservoir_global_scalar_count_(
              reservoir_system
                  ->numbering()
                  .petsc_global_scalar_count()) {}

    [[nodiscard]] PetscErrorCode
    initialize(
        double initial_bottom_hole_pressure_pa) {
        PetscErrorCode error =
            VecDuplicate(
                reservoir_system_
                    ->initial_state(),
                &reservoir_state_);
        if (error == PETSC_SUCCESS) {
            error =
                VecDuplicate(
                    reservoir_system_
                        ->initial_state(),
                    &reservoir_residual_);
        }
        if (error == PETSC_SUCCESS) {
            error =
                VecDuplicate(
                    reservoir_system_
                        ->initial_state(),
                    &reservoir_direction_);
        }
        if (error == PETSC_SUCCESS) {
            error =
                MatDuplicate(
                    reservoir_system_
                        ->jacobian_structure(),
                    MAT_DO_NOT_COPY_VALUES,
                    &reservoir_jacobian_);
        }
        if (error != PETSC_SUCCESS) {
            return error;
        }

        const PetscInt reservoir_local =
            reservoir_system_
                ->numbering()
                .petsc_local_owned_scalar_count();
        const PetscInt augmented_local =
            reservoir_local +
            (mpi_rank_ ==
                     well_owner_rank_
                 ? PetscInt{1}
                 : PetscInt{0});

        error =
            VecCreateMPI(
                comm_,
                augmented_local,
                global_scalar_count(),
                &initial_state_);
        if (error != PETSC_SUCCESS) {
            return error;
        }

        const PetscScalar*
            reservoir_values = nullptr;
        PetscScalar*
            augmented_values = nullptr;
        error =
            VecGetArrayRead(
                reservoir_system_
                    ->initial_state(),
                &reservoir_values);
        if (error == PETSC_SUCCESS) {
            error =
                VecGetArray(
                    initial_state_,
                    &augmented_values);
        }
        if (error == PETSC_SUCCESS) {
            for (PetscInt index = 0;
                 index < reservoir_local;
                 ++index) {
                augmented_values[index] =
                    reservoir_values[index];
            }
            if (mpi_rank_ ==
                well_owner_rank_) {
                augmented_values[
                    reservoir_local] =
                    static_cast<PetscScalar>(
                        initial_bottom_hole_pressure_pa);
            }
        }
        const PetscErrorCode restore_augmented =
            augmented_values != nullptr
                ? VecRestoreArray(
                      initial_state_,
                      &augmented_values)
                : PETSC_SUCCESS;
        const PetscErrorCode restore_reservoir =
            reservoir_values != nullptr
                ? VecRestoreArrayRead(
                      reservoir_system_
                          ->initial_state(),
                      &reservoir_values)
                : PETSC_SUCCESS;
        if (error == PETSC_SUCCESS &&
            restore_augmented !=
                PETSC_SUCCESS) {
            error = restore_augmented;
        }
        if (error == PETSC_SUCCESS &&
            restore_reservoir !=
                PETSC_SUCCESS) {
            error = restore_reservoir;
        }
        if (error != PETSC_SUCCESS) {
            return error;
        }

        error =
            create_augmented_jacobian_structure(
                initial_bottom_hole_pressure_pa);
        if (error != PETSC_SUCCESS) {
            return error;
        }
        return create_initial_row_scaling();
    }

    [[nodiscard]] PetscErrorCode
    sync_augmented_state(
        Vec augmented_state,
        double* bottom_hole_pressure_pa) {
        if (augmented_state == nullptr) {
            return PETSC_ERR_ARG_NULL;
        }

        PetscInt local_size = -1;
        PetscInt global_size = -1;
        PetscErrorCode error =
            VecGetLocalSize(
                augmented_state,
                &local_size);
        if (error == PETSC_SUCCESS) {
            error =
                VecGetSize(
                    augmented_state,
                    &global_size);
        }
        const PetscInt reservoir_local =
            reservoir_system_
                ->numbering()
                .petsc_local_owned_scalar_count();
        const PetscInt expected_local =
            reservoir_local +
            (mpi_rank_ ==
                     well_owner_rank_
                 ? PetscInt{1}
                 : PetscInt{0});
        if (error != PETSC_SUCCESS) {
            return error;
        }
        if (local_size != expected_local ||
            global_size !=
                global_scalar_count()) {
            return PETSC_ERR_ARG_SIZ;
        }

        const PetscScalar*
            augmented_values = nullptr;
        PetscScalar*
            reservoir_values = nullptr;
        error =
            VecGetArrayRead(
                augmented_state,
                &augmented_values);
        if (error == PETSC_SUCCESS) {
            error =
                VecGetArray(
                    reservoir_state_,
                    &reservoir_values);
        }

        double local_bhp = 0.0;
        if (error == PETSC_SUCCESS) {
            for (PetscInt index = 0;
                 index < reservoir_local;
                 ++index) {
                reservoir_values[index] =
                    augmented_values[index];
            }
            if (mpi_rank_ ==
                well_owner_rank_) {
                local_bhp =
                    static_cast<double>(
                        PetscRealPart(
                            augmented_values[
                                reservoir_local]));
            }
        }

        const PetscErrorCode restore_reservoir =
            reservoir_values != nullptr
                ? VecRestoreArray(
                      reservoir_state_,
                      &reservoir_values)
                : PETSC_SUCCESS;
        const PetscErrorCode restore_augmented =
            augmented_values != nullptr
                ? VecRestoreArrayRead(
                      augmented_state,
                      &augmented_values)
                : PETSC_SUCCESS;
        if (error == PETSC_SUCCESS &&
            restore_reservoir !=
                PETSC_SUCCESS) {
            error = restore_reservoir;
        }
        if (error == PETSC_SUCCESS &&
            restore_augmented !=
                PETSC_SUCCESS) {
            error = restore_augmented;
        }
        if (error != PETSC_SUCCESS) {
            return error;
        }

        double bhp =
            local_bhp;
        if (MPI_Bcast(
                &bhp,
                1,
                MPI_DOUBLE,
                well_owner_rank_,
                comm_) != MPI_SUCCESS) {
            return PETSC_ERR_MPI;
        }
        if (!std::isfinite(bhp) ||
            !(bhp > 0.0)) {
            if (bottom_hole_pressure_pa !=
                nullptr) {
                *bottom_hole_pressure_pa =
                    bhp;
            }
            return PETSC_SUCCESS;
        }
        if (bottom_hole_pressure_pa !=
            nullptr) {
            *bottom_hole_pressure_pa =
                bhp;
        }
        return PETSC_SUCCESS;
    }

    [[nodiscard]] PetscErrorCode
    extract_bottom_hole_pressure(
        Vec augmented_state,
        double* bottom_hole_pressure_pa) {
        if (bottom_hole_pressure_pa ==
            nullptr) {
            return PETSC_ERR_ARG_NULL;
        }
        return sync_augmented_state(
            augmented_state,
            bottom_hole_pressure_pa);
    }

    [[nodiscard]] double
    local_total_molar_production_rate() const {
        double total = 0.0;
        for (const auto& entry :
             source_context_
                 ->local_authoritative_evaluations()) {
            for (double source_rate :
                 entry
                     .linearization
                     .cell_source
                     .component_molar_rate_mol_per_s) {
                total -= source_rate;
            }
        }
        return total;
    }

    [[nodiscard]] PetscErrorCode
    copy_reservoir_residual_to_augmented(
        Vec augmented_residual) {
        const auto& numbering =
            reservoir_system_->numbering();
        const PetscInt start =
            numbering
                .petsc_owned_scalar_start();
        const PetscInt end =
            numbering
                .petsc_owned_scalar_end();
        const PetscInt count =
            end - start;
        const PetscScalar* values =
            nullptr;
        PetscErrorCode error =
            VecGetArrayRead(
                reservoir_residual_,
                &values);
        if (error != PETSC_SUCCESS) {
            return error;
        }

        std::vector<PetscInt> rows(
            static_cast<std::size_t>(
                count));
        for (PetscInt local = 0;
             local < count;
             ++local) {
            rows[
                static_cast<std::size_t>(
                    local)] =
                start + local;
        }

        error =
            VecSetValues(
                augmented_residual,
                count,
                rows.data(),
                values,
                ADD_VALUES);
        const PetscErrorCode restore =
            VecRestoreArrayRead(
                reservoir_residual_,
                &values);
        return error != PETSC_SUCCESS
            ? error
            : restore;
    }

    [[nodiscard]] PetscErrorCode
    insert_jrr(
        Mat augmented_jacobian) {
        const auto& numbering =
            reservoir_system_->numbering();
        const PetscInt start =
            numbering
                .petsc_owned_scalar_start();
        const PetscInt end =
            numbering
                .petsc_owned_scalar_end();

        for (PetscInt row = start;
             row < end;
             ++row) {
            PetscInt count = 0;
            const PetscInt* columns = nullptr;
            const PetscScalar* values = nullptr;
            PetscErrorCode error =
                MatGetRow(
                    reservoir_jacobian_,
                    row,
                    &count,
                    &columns,
                    &values);
            if (error != PETSC_SUCCESS) {
                return error;
            }
            error =
                MatSetValues(
                    augmented_jacobian,
                    1,
                    &row,
                    count,
                    columns,
                    values,
                    ADD_VALUES);
            const PetscErrorCode restore =
                MatRestoreRow(
                    reservoir_jacobian_,
                    row,
                    &count,
                    &columns,
                    &values);
            if (error != PETSC_SUCCESS) {
                return error;
            }
            if (restore != PETSC_SUCCESS) {
                return restore;
            }
        }
        return PETSC_SUCCESS;
    }

    [[nodiscard]] PetscErrorCode
    insert_well_coupling_blocks(
        Mat augmented_jacobian) {
        const auto& numbering =
            reservoir_system_->numbering();
        const std::size_t nc =
            numbering.component_count();
        const PetscInt well =
            well_global_scalar();

        for (const auto& entry :
             source_context_
                 ->local_authoritative_evaluations()) {
            const auto& record =
                numbering.cell(entry.cell);
            const auto& linearization =
                entry.linearization;
            const std::size_t q =
                record.scalar_count;
            if (linearization
                    .cell_source
                    .input_count != q ||
                linearization
                    .cell_source
                    .component_count() != nc ||
                linearization
                    .component_source_bhp_derivative_mol_per_pa_s
                    .size() != nc) {
                return PETSC_ERR_ARG_INCOMP;
            }

            double bulk_volume_m3 = 0.0;
            try {
                bulk_volume_m3 =
                    reservoir_system_
                        ->cell_bulk_volume_m3(
                            entry.cell);
            } catch (...) {
                return PETSC_ERR_ARG_INCOMP;
            }
            if (!std::isfinite(bulk_volume_m3) ||
                !(bulk_volume_m3 > 0.0)) {
                return PETSC_ERR_ARG_INCOMP;
            }

            for (std::size_t component = 0U;
                 component < nc;
                 ++component) {
                const PetscInt row =
                    record
                        .petsc_global_scalar_start +
                    static_cast<PetscInt>(
                        component);
                const PetscScalar value =
                    static_cast<PetscScalar>(
                        -linearization
                             .component_source_bhp_derivative_mol_per_pa_s[
                                 component] /
                        bulk_volume_m3);
                PetscErrorCode error =
                    MatSetValue(
                        augmented_jacobian,
                        row,
                        well,
                        value,
                        ADD_VALUES);
                if (error != PETSC_SUCCESS) {
                    return error;
                }
            }

            const PetscInt energy_row =
                record
                    .petsc_global_scalar_start +
                static_cast<PetscInt>(nc);
            PetscErrorCode error =
                MatSetValue(
                    augmented_jacobian,
                    energy_row,
                    well,
                    static_cast<PetscScalar>(
                        -linearization
                             .energy_source_bhp_derivative_w_per_pa /
                        bulk_volume_m3),
                    ADD_VALUES);
            if (error != PETSC_SUCCESS) {
                return error;
            }

            std::vector<PetscInt>
                columns(q);
            std::vector<PetscScalar>
                values(q, PetscScalar{0.0});
            for (std::size_t column = 0U;
                 column < q;
                 ++column) {
                columns[column] =
                    record
                        .petsc_global_scalar_start +
                    static_cast<PetscInt>(
                        column);
                double derivative = 0.0;
                for (std::size_t component = 0U;
                     component < nc;
                     ++component) {
                    derivative -=
                        linearization
                            .cell_source
                            .d_component_rate(
                                component,
                                column);
                }
                values[column] =
                    static_cast<PetscScalar>(
                        derivative);
            }

            error =
                MatSetValues(
                    augmented_jacobian,
                    1,
                    &well,
                    static_cast<PetscInt>(q),
                    columns.data(),
                    values.data(),
                    ADD_VALUES);
            if (error != PETSC_SUCCESS) {
                return error;
            }

            double bhp_derivative = 0.0;
            for (double value :
                 linearization
                     .component_source_bhp_derivative_mol_per_pa_s) {
                bhp_derivative -=
                    value;
            }
            error =
                MatSetValue(
                    augmented_jacobian,
                    well,
                    well,
                    static_cast<PetscScalar>(
                        bhp_derivative),
                    ADD_VALUES);
            if (error != PETSC_SUCCESS) {
                return error;
            }
        }
        return PETSC_SUCCESS;
    }

    [[nodiscard]] PetscErrorCode
    create_augmented_jacobian_structure(
        double initial_bottom_hole_pressure_pa) {
        using namespace
            mpmc::flow_discretization_petsc;

        source_context_->begin_evaluation(
            initial_bottom_hole_pressure_pa);
        PetscErrorCode error =
            VecCopy(
                reservoir_system_
                    ->initial_state(),
                reservoir_state_);
        if (error == PETSC_SUCCESS) {
            error =
                MatZeroEntries(
                    reservoir_jacobian_);
        }
        NaturalVariableSnesEvaluationStatus3D
            status =
                NaturalVariableSnesEvaluationStatus3D::
                    success;
        if (error == PETSC_SUCCESS) {
            error =
                reservoir_evaluator_.jacobian(
                    reservoir_state_,
                    reservoir_jacobian_,
                    reservoir_evaluator_
                        .user_context,
                    &status);
        }
        if (error == PETSC_SUCCESS &&
            status !=
                NaturalVariableSnesEvaluationStatus3D::
                    success) {
            error = PETSC_ERR_ARG_OUTOFRANGE;
        }
        if (error == PETSC_SUCCESS) {
            error =
                MatAssemblyBegin(
                    reservoir_jacobian_,
                    MAT_FINAL_ASSEMBLY);
        }
        if (error == PETSC_SUCCESS) {
            error =
                MatAssemblyEnd(
                    reservoir_jacobian_,
                    MAT_FINAL_ASSEMBLY);
        }
        if (error != PETSC_SUCCESS) {
            return error;
        }

        const auto& numbering =
            reservoir_system_->numbering();
        const PetscInt reservoir_local =
            numbering
                .petsc_local_owned_scalar_count();
        const PetscInt augmented_local =
            reservoir_local +
            (mpi_rank_ ==
                     well_owner_rank_
                 ? PetscInt{1}
                 : PetscInt{0});
        const PetscInt augmented_start =
            numbering
                .petsc_owned_scalar_start();
        const PetscInt augmented_end =
            numbering
                .petsc_owned_scalar_end() +
            (mpi_rank_ ==
                     well_owner_rank_
                 ? PetscInt{1}
                 : PetscInt{0});

        std::vector<PetscInt>
            diagonal_nnz(
                static_cast<std::size_t>(
                    augmented_local),
                PetscInt{0});
        std::vector<PetscInt>
            off_diagonal_nnz(
                diagonal_nnz.size(),
                PetscInt{0});

        const PetscInt start =
            numbering
                .petsc_owned_scalar_start();
        const PetscInt end =
            numbering
                .petsc_owned_scalar_end();
        for (PetscInt row = start;
             row < end;
             ++row) {
            PetscInt count = 0;
            const PetscInt* columns = nullptr;
            const PetscScalar* values = nullptr;
            error =
                MatGetRow(
                    reservoir_jacobian_,
                    row,
                    &count,
                    &columns,
                    &values);
            if (error != PETSC_SUCCESS) {
                return error;
            }

            PetscInt diagonal = 0;
            PetscInt off_diagonal = 0;
            for (PetscInt index = 0;
                 index < count;
                 ++index) {
                if (columns[index] >=
                        augmented_start &&
                    columns[index] <
                        augmented_end) {
                    ++diagonal;
                } else {
                    ++off_diagonal;
                }
            }
            const PetscInt well =
                well_global_scalar();
            if (well >= augmented_start &&
                well < augmented_end) {
                ++diagonal;
            } else {
                ++off_diagonal;
            }

            diagonal_nnz[
                static_cast<std::size_t>(
                    row - start)] =
                diagonal;
            off_diagonal_nnz[
                static_cast<std::size_t>(
                    row - start)] =
                off_diagonal;

            error =
                MatRestoreRow(
                    reservoir_jacobian_,
                    row,
                    &count,
                    &columns,
                    &values);
            if (error != PETSC_SUCCESS) {
                return error;
            }
        }

        std::size_t local_completion_width =
            0U;
        std::size_t global_completion_width =
            0U;
        for (const auto& record :
             numbering.cells()) {
            if (record.owner_rank ==
                    numbering.local_rank() &&
                source_context_->well()
                        .find_connection(
                            record.cell_global) !=
                    nullptr) {
                local_completion_width +=
                    record.scalar_count;
            }
        }
        unsigned long long local_width =
            static_cast<unsigned long long>(
                local_completion_width);
        unsigned long long global_width = 0ULL;
        if (MPI_Allreduce(
                &local_width,
                &global_width,
                1,
                MPI_UNSIGNED_LONG_LONG,
                MPI_SUM,
                comm_) != MPI_SUCCESS) {
            return PETSC_ERR_MPI;
        }
        global_completion_width =
            static_cast<std::size_t>(
                global_width);

        if (mpi_rank_ ==
            well_owner_rank_) {
            if (local_completion_width >
                    static_cast<std::size_t>(
                        std::numeric_limits<
                            PetscInt>::max()) ||
                global_completion_width <
                    local_completion_width ||
                global_completion_width -
                        local_completion_width >
                    static_cast<std::size_t>(
                        std::numeric_limits<
                            PetscInt>::max())) {
                return PETSC_ERR_ARG_OUTOFRANGE;
            }
            const std::size_t well_local_index =
                diagonal_nnz.size() - 1U;
            diagonal_nnz[well_local_index] =
                static_cast<PetscInt>(
                    local_completion_width) +
                PetscInt{1};
            off_diagonal_nnz[well_local_index] =
                static_cast<PetscInt>(
                    global_completion_width -
                    local_completion_width);
        }

        error =
            MatCreateAIJ(
                comm_,
                augmented_local,
                augmented_local,
                global_scalar_count(),
                global_scalar_count(),
                0,
                diagonal_nnz.data(),
                0,
                off_diagonal_nnz.data(),
                &jacobian_);
        if (error == PETSC_SUCCESS) {
            error =
                MatSetOption(
                    jacobian_,
                    MAT_NEW_NONZERO_ALLOCATION_ERR,
                    PETSC_TRUE);
        }
        return error;
    }

    [[nodiscard]] PetscErrorCode
    create_initial_row_scaling() {
        using namespace
            mpmc::flow_discretization_petsc;

        PetscErrorCode error =
            MatZeroEntries(
                jacobian_);
        NaturalVariableSnesEvaluationStatus3D
            status =
                NaturalVariableSnesEvaluationStatus3D::
                    success;
        if (error == PETSC_SUCCESS) {
            error =
                evaluate_jacobian(
                    initial_state_,
                    jacobian_,
                    &status);
        }
        if (error == PETSC_SUCCESS &&
            status !=
                NaturalVariableSnesEvaluationStatus3D::
                    success) {
            error = PETSC_ERR_ARG_OUTOFRANGE;
        }
        if (error == PETSC_SUCCESS) {
            error =
                MatAssemblyBegin(
                    jacobian_,
                    MAT_FINAL_ASSEMBLY);
        }
        if (error == PETSC_SUCCESS) {
            error =
                MatAssemblyEnd(
                    jacobian_,
                    MAT_FINAL_ASSEMBLY);
        }
        if (error != PETSC_SUCCESS) {
            return error;
        }

        PetscInt start = -1;
        PetscInt end = -1;
        error =
            MatGetOwnershipRange(
                jacobian_,
                &start,
                &end);
        if (error != PETSC_SUCCESS ||
            start < 0 ||
            end < start) {
            return error != PETSC_SUCCESS
                ? error
                : PETSC_ERR_ARG_INCOMP;
        }

        error =
            VecDuplicate(
                initial_state_,
                &row_scaling_);
        if (error != PETSC_SUCCESS) {
            return error;
        }

        PetscScalar* scaling = nullptr;
        error =
            VecGetArray(
                row_scaling_,
                &scaling);
        if (error != PETSC_SUCCESS) {
            return error;
        }

        for (PetscInt row = start;
             row < end;
             ++row) {
            PetscInt count = 0;
            const PetscInt* columns = nullptr;
            const PetscScalar* values = nullptr;
            error =
                MatGetRow(
                    jacobian_,
                    row,
                    &count,
                    &columns,
                    &values);
            if (error != PETSC_SUCCESS) {
                break;
            }
            double maximum = 0.0;
            for (PetscInt index = 0;
                 index < count;
                 ++index) {
                maximum =
                    std::max(
                        maximum,
                        std::abs(
                            static_cast<double>(
                                PetscRealPart(
                                    values[index]))));
            }
            const PetscErrorCode restore =
                MatRestoreRow(
                    jacobian_,
                    row,
                    &count,
                    &columns,
                    &values);
            if (error == PETSC_SUCCESS &&
                restore != PETSC_SUCCESS) {
                error = restore;
            }
            if (error != PETSC_SUCCESS) {
                break;
            }
            if (!std::isfinite(maximum) ||
                !(maximum > 0.0)) {
                error =
                    PETSC_ERR_ARG_WRONGSTATE;
                break;
            }
            scaling[
                static_cast<std::size_t>(
                    row - start)] =
                static_cast<PetscScalar>(
                    1.0 / maximum);
        }

        const PetscErrorCode restore =
            VecRestoreArray(
                row_scaling_,
                &scaling);
        return error != PETSC_SUCCESS
            ? error
            : restore;
    }

    [[nodiscard]] PetscErrorCode
    evaluate_function(
        Vec state,
        Vec residual,
        mpmc::flow_discretization_petsc::
            NaturalVariableSnesEvaluationStatus3D*
                status) {
        using namespace
            mpmc::flow_discretization_petsc;
        if (status == nullptr) {
            return PETSC_ERR_ARG_NULL;
        }

        double bhp = 0.0;
        PetscErrorCode error =
            sync_augmented_state(
                state,
                &bhp);
        if (error != PETSC_SUCCESS) {
            return error;
        }
        if (!std::isfinite(bhp) ||
            !(bhp > 0.0)) {
            *status =
                NaturalVariableSnesEvaluationStatus3D::
                    domain_error;
            return PETSC_SUCCESS;
        }

        try {
            source_context_
                ->begin_evaluation(
                    bhp);
        } catch (...) {
            *status =
                NaturalVariableSnesEvaluationStatus3D::
                    domain_error;
            return PETSC_SUCCESS;
        }

        error =
            VecSet(
                reservoir_residual_,
                PetscScalar{0.0});
        if (error == PETSC_SUCCESS) {
            error =
                reservoir_evaluator_.function(
                    reservoir_state_,
                    reservoir_residual_,
                    reservoir_evaluator_
                        .user_context,
                    status);
        }
        if (error != PETSC_SUCCESS ||
            *status !=
                NaturalVariableSnesEvaluationStatus3D::
                    success) {
            return error;
        }

        error =
            VecAssemblyBegin(
                reservoir_residual_);
        if (error == PETSC_SUCCESS) {
            error =
                VecAssemblyEnd(
                    reservoir_residual_);
        }
        if (error == PETSC_SUCCESS) {
            error =
                copy_reservoir_residual_to_augmented(
                    residual);
        }
        if (error != PETSC_SUCCESS) {
            return error;
        }

        const double local_rate =
            local_total_molar_production_rate();
        error =
            VecSetValue(
                residual,
                well_global_scalar(),
                static_cast<PetscScalar>(
                    local_rate),
                ADD_VALUES);
        if (error == PETSC_SUCCESS &&
            mpi_rank_ ==
                well_owner_rank_) {
            error =
                VecSetValue(
                    residual,
                    well_global_scalar(),
                    static_cast<PetscScalar>(
                        -target_total_molar_rate_mol_per_s_),
                    ADD_VALUES);
        }
        return error;
    }

    [[nodiscard]] PetscErrorCode
    evaluate_jacobian(
        Vec state,
        Mat jacobian,
        mpmc::flow_discretization_petsc::
            NaturalVariableSnesEvaluationStatus3D*
                status) {
        using namespace
            mpmc::flow_discretization_petsc;
        if (status == nullptr) {
            return PETSC_ERR_ARG_NULL;
        }

        double bhp = 0.0;
        PetscErrorCode error =
            sync_augmented_state(
                state,
                &bhp);
        if (error != PETSC_SUCCESS) {
            return error;
        }
        if (!std::isfinite(bhp) ||
            !(bhp > 0.0)) {
            *status =
                NaturalVariableSnesEvaluationStatus3D::
                    domain_error;
            return PETSC_SUCCESS;
        }

        try {
            source_context_
                ->begin_evaluation(
                    bhp);
        } catch (...) {
            *status =
                NaturalVariableSnesEvaluationStatus3D::
                    domain_error;
            return PETSC_SUCCESS;
        }

        error =
            MatZeroEntries(
                reservoir_jacobian_);
        if (error == PETSC_SUCCESS) {
            error =
                reservoir_evaluator_.jacobian(
                    reservoir_state_,
                    reservoir_jacobian_,
                    reservoir_evaluator_
                        .user_context,
                    status);
        }
        if (error != PETSC_SUCCESS ||
            *status !=
                NaturalVariableSnesEvaluationStatus3D::
                    success) {
            return error;
        }
        error =
            MatAssemblyBegin(
                reservoir_jacobian_,
                MAT_FINAL_ASSEMBLY);
        if (error == PETSC_SUCCESS) {
            error =
                MatAssemblyEnd(
                    reservoir_jacobian_,
                    MAT_FINAL_ASSEMBLY);
        }
        if (error == PETSC_SUCCESS) {
            error =
                insert_jrr(
                    jacobian);
        }
        if (error == PETSC_SUCCESS) {
            error =
                insert_well_coupling_blocks(
                    jacobian);
        }
        return error;
    }

    [[nodiscard]] PetscErrorCode
    precheck(
        Vec state,
        Vec search_direction,
        PetscBool* changed_direction) {
        if (state == nullptr ||
            search_direction == nullptr ||
            changed_direction == nullptr) {
            return PETSC_ERR_ARG_NULL;
        }
        *changed_direction =
            PETSC_FALSE;

        PetscInt state_local = -1;
        PetscInt direction_local = -1;
        PetscErrorCode error =
            VecGetLocalSize(
                state,
                &state_local);
        if (error == PETSC_SUCCESS) {
            error =
                VecGetLocalSize(
                    search_direction,
                    &direction_local);
        }
        if (error != PETSC_SUCCESS) {
            return error;
        }
        if (state_local !=
            direction_local) {
            return PETSC_ERR_ARG_SIZ;
        }

        const PetscScalar*
            state_values = nullptr;
        const PetscScalar*
            direction_values = nullptr;
        error =
            VecGetArrayRead(
                state,
                &state_values);
        if (error == PETSC_SUCCESS) {
            error =
                VecGetArrayRead(
                    search_direction,
                    &direction_values);
        }
        if (error != PETSC_SUCCESS) {
            if (state_values != nullptr) {
                (void)VecRestoreArrayRead(
                    state,
                    &state_values);
            }
            return error;
        }

        double local_scale = 1.0;
        bool invalid = false;
        for (PetscInt index = 0;
             index < state_local;
             ++index) {
            const double x =
                static_cast<double>(
                    PetscRealPart(
                        state_values[index]));
            const double y =
                static_cast<double>(
                    PetscRealPart(
                        direction_values[index]));
            if (!std::isfinite(x) ||
                !std::isfinite(y) ||
                !(x > 0.0)) {
                invalid = true;
                break;
            }
            if (y > 0.0 &&
                x - y <= 0.0) {
                local_scale =
                    std::min(
                        local_scale,
                        0.8 * x / y);
            }
        }

        const PetscErrorCode restore_direction =
            VecRestoreArrayRead(
                search_direction,
                &direction_values);
        const PetscErrorCode restore_state =
            VecRestoreArrayRead(
                state,
                &state_values);
        if (restore_direction !=
            PETSC_SUCCESS) {
            return restore_direction;
        }
        if (restore_state !=
            PETSC_SUCCESS) {
            return restore_state;
        }
        if (invalid) {
            return PETSC_ERR_ARG_OUTOFRANGE;
        }

        double global_scale = 1.0;
        if (MPI_Allreduce(
                &local_scale,
                &global_scale,
                1,
                MPI_DOUBLE,
                MPI_MIN,
                comm_) != MPI_SUCCESS) {
            return PETSC_ERR_MPI;
        }
        if (!std::isfinite(global_scale) ||
            !(global_scale > 0.0)) {
            return PETSC_ERR_FP;
        }
        if (global_scale < 1.0) {
            error =
                VecScale(
                    search_direction,
                    static_cast<PetscScalar>(
                        global_scale));
            if (error == PETSC_SUCCESS) {
                *changed_direction =
                    PETSC_TRUE;
            }
        }
        return error;
    }

    static PetscErrorCode
    snes_function(
        Vec state,
        Vec residual,
        void* raw_context,
        mpmc::flow_discretization_petsc::
            NaturalVariableSnesEvaluationStatus3D*
                status) {
        if (raw_context == nullptr) {
            return PETSC_ERR_ARG_NULL;
        }
        return static_cast<
            FixedTotalMolarRateWellControlSystem3D*>(
                raw_context)
            ->evaluate_function(
                state,
                residual,
                status);
    }

    static PetscErrorCode
    snes_jacobian(
        Vec state,
        Mat jacobian,
        void* raw_context,
        mpmc::flow_discretization_petsc::
            NaturalVariableSnesEvaluationStatus3D*
                status) {
        if (raw_context == nullptr) {
            return PETSC_ERR_ARG_NULL;
        }
        return static_cast<
            FixedTotalMolarRateWellControlSystem3D*>(
                raw_context)
            ->evaluate_jacobian(
                state,
                jacobian,
                status);
    }

    static PetscErrorCode
    snes_precheck(
        Vec state,
        Vec search_direction,
        void* raw_context,
        PetscBool* changed_direction) {
        if (raw_context == nullptr) {
            return PETSC_ERR_ARG_NULL;
        }
        return static_cast<
            FixedTotalMolarRateWellControlSystem3D*>(
                raw_context)
            ->precheck(
                state,
                search_direction,
                changed_direction);
    }

    static void cleanup_solve(
        Vec* state,
        Vec* residual,
        Mat* jacobian,
        SNES* snes) noexcept {
        if (snes != nullptr &&
            *snes != nullptr) {
            (void)SNESDestroy(snes);
        }
        if (jacobian != nullptr &&
            *jacobian != nullptr) {
            (void)MatDestroy(jacobian);
        }
        if (residual != nullptr &&
            *residual != nullptr) {
            (void)VecDestroy(residual);
        }
        if (state != nullptr &&
            *state != nullptr) {
            (void)VecDestroy(state);
        }
    }

    MPI_Comm comm_;
    mpmc::flow_discretization_petsc::
        PhaseTransitionRebuiltNaturalVariableSystem3D*
            reservoir_system_{};
    FixedTotalMolarRateWellSourceEvaluatorContext3D*
        source_context_{};
    double target_total_molar_rate_mol_per_s_{};
    int mpi_rank_{-1};
    int mpi_size_{};
    int well_owner_rank_{-1};
    mpmc::flow_discretization_petsc::
        NaturalVariableSnesEvaluator3D
            reservoir_evaluator_{};
    PetscInt reservoir_global_scalar_count_{};
    Vec reservoir_state_{};
    Vec reservoir_residual_{};
    Vec reservoir_direction_{};
    Mat reservoir_jacobian_{};
    Vec initial_state_{};
    Mat jacobian_{};
    Vec row_scaling_{};
};

} // namespace mpmc::well_discretization_petsc

#endif // MPMC_WELL_DISCRETIZATION_PETSC_FIXED_TOTAL_MOLAR_RATE_CONTROL_HPP
