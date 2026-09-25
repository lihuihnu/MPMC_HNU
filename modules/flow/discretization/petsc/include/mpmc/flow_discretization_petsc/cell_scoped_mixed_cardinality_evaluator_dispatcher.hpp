#ifndef MPMC_FLOW_DISCRETIZATION_PETSC_CELL_SCOPED_MIXED_CARDINALITY_EVALUATOR_DISPATCHER_HPP
#define MPMC_FLOW_DISCRETIZATION_PETSC_CELL_SCOPED_MIXED_CARDINALITY_EVALUATOR_DISPATCHER_HPP

#include <mpmc/flow_discretization_petsc/mixed_cardinality_physical_snes_assembly.hpp>

#include <petscsys.h>

#include <algorithm>
#include <cstddef>
#include <optional>
#include <span>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace mpmc::flow_discretization_petsc {

inline constexpr std::string_view
    cell_scoped_mixed_cardinality_evaluator_dispatcher_convention =
        "flow_discretization_petsc/cell-scoped-mixed-cardinality-evaluator-dispatcher/v1";

template <typename Binding>
struct CellScopedCurrentCellEvaluatorEntry3D {
    mpmc::mesh::GlobalEntityId cell_global{
        mpmc::mesh::GlobalEntityId::value_type{0}};
    Binding binding;
};

class CellScopedMixedCardinalityEvaluatorDispatcher3D {
public:
    using SingleEntry =
        CellScopedCurrentCellEvaluatorEntry3D<
            SinglePhaseCurrentCellEvaluatorBinding3D>;
    using TwoEntry =
        CellScopedCurrentCellEvaluatorEntry3D<
            TwoPhaseCurrentCellEvaluatorBinding3D>;
    using ThreeEntry =
        CellScopedCurrentCellEvaluatorEntry3D<
            FixedThreePhaseCurrentCellEvaluatorBinding3D>;

    CellScopedMixedCardinalityEvaluatorDispatcher3D(
        std::vector<SingleEntry> single_phase,
        std::vector<TwoEntry> two_phase,
        std::vector<ThreeEntry> three_phase)
        : single_phase_(std::move(single_phase)),
          two_phase_(std::move(two_phase)),
          three_phase_(std::move(three_phase)) {
        validate_and_sort(&single_phase_);
        validate_and_sort(&two_phase_);
        validate_and_sort(&three_phase_);
    }

    [[nodiscard]]
    MixedCardinalityPhysicalCellEvaluatorBindings3D
    bindings() noexcept {
        return {
            {&evaluate_single_phase, this},
            {&evaluate_two_phase, this},
            {&evaluate_three_phase, this}};
    }

private:
    template <class Entry>
    static void validate_and_sort(
        std::vector<Entry>* entries) {
        std::sort(
            entries->begin(),
            entries->end(),
            [](const auto& first,
               const auto& second) {
                return first.cell_global <
                    second.cell_global;
            });
        for (std::size_t index = 0U;
             index < entries->size();
             ++index) {
            const auto& entry =
                (*entries)[index];
            if (entry.binding.evaluator ==
                    nullptr ||
                (index > 0U &&
                 (*entries)[index - 1U]
                         .cell_global ==
                     entry.cell_global)) {
                throw std::invalid_argument(
                    "cell-scoped mixed-cardinality evaluator dispatcher has invalid/duplicate binding");
            }
        }
    }

    template <class Entry>
    [[nodiscard]] static const Entry*
    find(
        const std::vector<Entry>& entries,
        mpmc::mesh::GlobalEntityId
            cell_global) noexcept {
        const auto found =
            std::lower_bound(
                entries.begin(),
                entries.end(),
                cell_global,
                [](const auto& entry,
                   const auto& id) {
                    return entry.cell_global <
                        id;
                });
        return found != entries.end() &&
                found->cell_global ==
                    cell_global
            ? &*found
            : nullptr;
    }

    static PetscErrorCode
    evaluate_single_phase(
        mpmc::mesh::LocalIndex cell,
        mpmc::mesh::GlobalEntityId cell_global,
        std::span<const double> natural_variables,
        const mpmc::flow::NaturalVariableLayout1P& layout,
        std::span<const std::string> component_ids,
        void* raw_context,
        std::optional<
            SinglePhaseCurrentCellLinearization3D>* output,
        NaturalVariableSnesEvaluationStatus3D* status) {
        if (raw_context == nullptr) {
            return PETSC_ERR_ARG_NULL;
        }
        auto* self =
            static_cast<
                CellScopedMixedCardinalityEvaluatorDispatcher3D*>(
                    raw_context);
        const auto* entry =
            find(self->single_phase_, cell_global);
        if (entry == nullptr) {
            return PETSC_ERR_ARG_INCOMP;
        }
        return entry->binding.evaluator(
            cell,
            cell_global,
            natural_variables,
            layout,
            component_ids,
            entry->binding.user_context,
            output,
            status);
    }

    static PetscErrorCode
    evaluate_two_phase(
        mpmc::mesh::LocalIndex cell,
        mpmc::mesh::GlobalEntityId cell_global,
        std::span<const double> natural_variables,
        const mpmc::flow::NaturalVariableLayout2P& layout,
        std::span<const std::string> component_ids,
        void* raw_context,
        std::optional<
            TwoPhaseCurrentCellLinearization3D>* output,
        NaturalVariableSnesEvaluationStatus3D* status) {
        if (raw_context == nullptr) {
            return PETSC_ERR_ARG_NULL;
        }
        auto* self =
            static_cast<
                CellScopedMixedCardinalityEvaluatorDispatcher3D*>(
                    raw_context);
        const auto* entry =
            find(self->two_phase_, cell_global);
        if (entry == nullptr) {
            return PETSC_ERR_ARG_INCOMP;
        }
        return entry->binding.evaluator(
            cell,
            cell_global,
            natural_variables,
            layout,
            component_ids,
            entry->binding.user_context,
            output,
            status);
    }

    static PetscErrorCode
    evaluate_three_phase(
        mpmc::mesh::LocalIndex cell,
        mpmc::mesh::GlobalEntityId cell_global,
        std::span<const double> natural_variables,
        const mpmc::flow::NaturalVariableLayout3P& layout,
        std::span<const std::string> component_ids,
        void* raw_context,
        std::optional<
            FixedThreePhaseCurrentCellLinearization3D>* output,
        NaturalVariableSnesEvaluationStatus3D* status) {
        if (raw_context == nullptr) {
            return PETSC_ERR_ARG_NULL;
        }
        auto* self =
            static_cast<
                CellScopedMixedCardinalityEvaluatorDispatcher3D*>(
                    raw_context);
        const auto* entry =
            find(self->three_phase_, cell_global);
        if (entry == nullptr) {
            return PETSC_ERR_ARG_INCOMP;
        }
        return entry->binding.evaluator(
            cell,
            cell_global,
            natural_variables,
            layout,
            component_ids,
            entry->binding.user_context,
            output,
            status);
    }

    std::vector<SingleEntry> single_phase_;
    std::vector<TwoEntry> two_phase_;
    std::vector<ThreeEntry> three_phase_;
};

} // namespace mpmc::flow_discretization_petsc

#endif // MPMC_FLOW_DISCRETIZATION_PETSC_CELL_SCOPED_MIXED_CARDINALITY_EVALUATOR_DISPATCHER_HPP
