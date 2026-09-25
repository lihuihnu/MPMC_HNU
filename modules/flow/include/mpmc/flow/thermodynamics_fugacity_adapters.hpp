#ifndef MPMC_FLOW_THERMODYNAMICS_FUGACITY_ADAPTERS_HPP
#define MPMC_FLOW_THERMODYNAMICS_FUGACITY_ADAPTERS_HPP

#include <mpmc/flow/fugacity_equilibrium_residual.hpp>
#include <mpmc/thermodynamics/selected_phase_fugacity.hpp>

#include <array>
#include <concepts>
#include <cstddef>
#include <span>
#include <stdexcept>
#include <utility>

namespace mpmc::flow {

template <std::floating_point T = double>
class Pr76SelectedPhaseFugacityEvaluator3P {
public:
    using Selection = thermodynamics::Pr76SelectedPhase;

    Pr76SelectedPhaseFugacityEvaluator3P(
        const thermodynamics::Pr76Phase<T>& model,
        std::array<Selection, 3> selections)
        : model_(&model),
          selections_(std::move(selections)) {}

    template <typename Number>
    [[nodiscard]]
    PhaseLnFugacityCoefficientEvaluation<Number>
    operator()(
        PhaseSlot3 slot,
        const Number& pressure_pa,
        const Number& temperature_k,
        std::span<const Number> composition) const {
        const std::size_t phase =
            natural_variable_detail::checked_phase_index(slot);
        thermodynamics::Pr76PhaseWorkspace<Number> workspace;
        auto selected =
            thermodynamics::evaluate_selected_phase_fugacity(
                *model_,
                pressure_pa,
                temperature_k,
                composition,
                selections_[phase],
                workspace);
        return {std::move(selected.ln_phi)};
    }

    [[nodiscard]] const Selection&
    selection(PhaseSlot3 slot) const {
        return selections_[
            natural_variable_detail::checked_phase_index(slot)];
    }

private:
    const thermodynamics::Pr76Phase<T>* model_;
    std::array<Selection, 3> selections_;
};

template <std::floating_point T = double>
class Sw92SelectedPhaseFugacityEvaluator3P {
public:
    using Selection =
        thermodynamics::Sw92SelectedPhase<T>;

    Sw92SelectedPhaseFugacityEvaluator3P(
        const thermodynamics::Sw92Phase<T>& model,
        std::array<Selection, 3> selections)
        : model_(&model),
          selections_(std::move(selections)) {}

    template <typename Number>
    [[nodiscard]]
    PhaseLnFugacityCoefficientEvaluation<Number>
    operator()(
        PhaseSlot3 slot,
        const Number& pressure_pa,
        const Number& temperature_k,
        std::span<const Number> composition) const {
        const std::size_t phase =
            natural_variable_detail::checked_phase_index(slot);
        thermodynamics::Sw92PhaseWorkspace<Number> workspace;
        auto selected =
            thermodynamics::evaluate_selected_phase_fugacity(
                *model_,
                pressure_pa,
                temperature_k,
                composition,
                selections_[phase],
                workspace);
        return {std::move(selected.ln_phi)};
    }

    [[nodiscard]] const Selection&
    selection(PhaseSlot3 slot) const {
        return selections_[
            natural_variable_detail::checked_phase_index(slot)];
    }

private:
    const thermodynamics::Sw92Phase<T>* model_;
    std::array<Selection, 3> selections_;
};

class CpaSelectedPhaseFugacityEvaluator3P {
public:
    using Selection =
        thermodynamics::CpaSelectedPhase;

    CpaSelectedPhaseFugacityEvaluator3P(
        const thermodynamics::CpaPtPhase& model,
        std::array<Selection, 3> selections)
        : model_(&model),
          selections_(std::move(selections)) {}

    template <typename Number>
    [[nodiscard]]
    PhaseLnFugacityCoefficientEvaluation<Number>
    operator()(
        PhaseSlot3 slot,
        const Number& pressure_pa,
        const Number& temperature_k,
        std::span<const Number> composition) const {
        const std::size_t phase =
            natural_variable_detail::checked_phase_index(slot);
        auto selected =
            thermodynamics::evaluate_selected_phase_fugacity(
                *model_,
                pressure_pa,
                temperature_k,
                composition,
                selections_[phase]);
        return {std::move(selected.ln_phi)};
    }

    [[nodiscard]] const Selection&
    selection(PhaseSlot3 slot) const {
        return selections_[
            natural_variable_detail::checked_phase_index(slot)];
    }

private:
    const thermodynamics::CpaPtPhase* model_;
    std::array<Selection, 3> selections_;
};

static_assert(
    PhaseLnFugacityCoefficientEvaluator3P<
        Pr76SelectedPhaseFugacityEvaluator3P<double>,
        double>);
static_assert(
    PhaseLnFugacityCoefficientEvaluator3P<
        Sw92SelectedPhaseFugacityEvaluator3P<double>,
        double>);
static_assert(
    PhaseLnFugacityCoefficientEvaluator3P<
        CpaSelectedPhaseFugacityEvaluator3P,
        double>);

} // namespace mpmc::flow

#endif // MPMC_FLOW_THERMODYNAMICS_FUGACITY_ADAPTERS_HPP
