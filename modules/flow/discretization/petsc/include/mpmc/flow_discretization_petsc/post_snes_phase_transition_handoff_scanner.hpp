#ifndef MPMC_FLOW_DISCRETIZATION_PETSC_POST_SNES_PHASE_TRANSITION_HANDOFF_SCANNER_HPP
#define MPMC_FLOW_DISCRETIZATION_PETSC_POST_SNES_PHASE_TRANSITION_HANDOFF_SCANNER_HPP

#include <mpmc/flow_discretization_petsc/post_snes_phase_transition_controller.hpp>

#include <petscsys.h>

#include <string_view>
#include <utility>
#include <vector>

namespace mpmc::flow_discretization_petsc {

inline constexpr std::string_view
    post_snes_phase_transition_handoff_scanner_convention =
        "flow_discretization_petsc/post-snes-phase-transition-handoff-scanner/v1";

struct PostSnesPhaseTransitionHandoffScannerContext3D {
    std::vector<PostSnesPhaseTransitionProposal3D>
        initial_local_owned_proposals;
    PostSnesPhaseTransitionScanner3D
        subsequent_scanner{};
    void* subsequent_scanner_context{};
    bool initial_batch_consumed{};
};

inline PetscErrorCode
scan_post_snes_phase_transition_handoff_3d(
    const PhaseTransitionRebuiltNaturalVariableSystem3D&
        system,
    Vec converged_state,
    const VariableCardinalityNaturalVariableSnesSolveReport3D&
        solve_report,
    void* raw_context,
    PostSnesPhaseTransitionScanStatus3D*
        scan_status,
    std::vector<
        PostSnesPhaseTransitionProposal3D>*
            local_owned_proposals) {
    if (raw_context == nullptr ||
        scan_status == nullptr ||
        local_owned_proposals == nullptr) {
        return PETSC_ERR_ARG_NULL;
    }

    auto* context =
        static_cast<
            PostSnesPhaseTransitionHandoffScannerContext3D*>(
                raw_context);
    if (!context->initial_batch_consumed) {
        *scan_status =
            PostSnesPhaseTransitionScanStatus3D::
                complete;
        *local_owned_proposals =
            std::move(
                context->initial_local_owned_proposals);
        context->initial_local_owned_proposals.clear();
        context->initial_batch_consumed = true;
        return PETSC_SUCCESS;
    }

    if (context->subsequent_scanner == nullptr) {
        return PETSC_ERR_ARG_WRONGSTATE;
    }
    return context->subsequent_scanner(
        system,
        converged_state,
        solve_report,
        context->subsequent_scanner_context,
        scan_status,
        local_owned_proposals);
}

} // namespace mpmc::flow_discretization_petsc

#endif // MPMC_FLOW_DISCRETIZATION_PETSC_POST_SNES_PHASE_TRANSITION_HANDOFF_SCANNER_HPP
