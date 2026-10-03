#ifndef MPMC_FLOW_TWO_PHASE_CELL_STATE_HPP
#define MPMC_FLOW_TWO_PHASE_CELL_STATE_HPP

#include <mpmc/flow/natural_variable_cell_state.hpp>

#include <array>
#include <cmath>
#include <cstddef>
#include <limits>
#include <optional>
#include <span>
#include <stdexcept>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

// Immutable state, frozen composition pivot and natural-variable indexing.
namespace mpmc::flow {

inline constexpr std::string_view
    two_phase_natural_variable_convention =
        "flow/natural-variable/two-phase-reduction/v1";

class NaturalVariableCompositionPivot2P {
public:
    [[nodiscard]] static NaturalVariableCompositionPivot2P
    fixed_last(std::size_t component_count) {
        validate_component_count(component_count);
        return {
            component_count,
            {component_count - 1U,
             component_count - 1U}};
    }

    [[nodiscard]] static NaturalVariableCompositionPivot2P
    from_dependent_components(
        std::size_t component_count,
        std::array<std::size_t, 2>
            dependent_components) {
        validate_component_count(component_count);
        for (const auto dependent :
             dependent_components) {
            if (dependent >= component_count) {
                throw std::invalid_argument(
                    "mpmc::flow::NaturalVariableCompositionPivot2P: dependent component out of range");
            }
        }
        return {
            component_count,
            dependent_components};
    }

    [[nodiscard]] std::size_t
    component_count() const noexcept {
        return component_count_;
    }

    [[nodiscard]] const std::array<std::size_t, 2>&
    dependent_components() const noexcept {
        return dependent_components_;
    }

    [[nodiscard]] std::size_t
    dependent_component(
        std::size_t phase) const {
        if (phase >= 2U) {
            throw std::out_of_range(
                "mpmc::flow::NaturalVariableCompositionPivot2P: phase out of range");
        }
        return dependent_components_[phase];
    }

private:
    NaturalVariableCompositionPivot2P(
        std::size_t component_count,
        std::array<std::size_t, 2>
            dependent_components)
        : component_count_(component_count),
          dependent_components_(
              dependent_components) {}

    static void validate_component_count(
        std::size_t component_count) {
        if (component_count < 2U) {
            throw std::invalid_argument(
                "mpmc::flow::NaturalVariableCompositionPivot2P: at least two components are required");
        }
    }

    std::size_t component_count_{};
    std::array<std::size_t, 2>
        dependent_components_{};
};

class NaturalVariableLayout2P {
public:
    explicit NaturalVariableLayout2P(
        std::size_t component_count)
        : NaturalVariableLayout2P(
              NaturalVariableCompositionPivot2P::
                  fixed_last(component_count)) {}

    explicit NaturalVariableLayout2P(
        NaturalVariableCompositionPivot2P pivot)
        : component_count_(
              pivot.component_count()),
          composition_pivot_(
              std::move(pivot)) {}

    [[nodiscard]] std::size_t
    component_count() const noexcept {
        return component_count_;
    }

    [[nodiscard]] static constexpr std::size_t
    phase_count() noexcept {
        return 2U;
    }

    [[nodiscard]] std::size_t
    unknown_count() const noexcept {
        return 2U * component_count_ + 1U;
    }

    [[nodiscard]] std::size_t
    equation_count() const noexcept {
        return unknown_count();
    }

    [[nodiscard]] static constexpr std::size_t
    pressure_unknown_index() noexcept {
        return 0U;
    }

    [[nodiscard]] static constexpr std::size_t
    temperature_unknown_index() noexcept {
        return 1U;
    }

    [[nodiscard]] static constexpr std::size_t
    independent_saturation_unknown_index() noexcept {
        return 2U;
    }

    [[nodiscard]] const NaturalVariableCompositionPivot2P&
    composition_pivot() const noexcept {
        return composition_pivot_;
    }

    [[nodiscard]] std::size_t
    dependent_composition_component(
        std::size_t phase) const {
        return composition_pivot_.
            dependent_component(phase);
    }

    [[nodiscard]] std::size_t
    independent_composition_component(
        std::size_t phase,
        std::size_t independent_rank) const {
        if (phase >= 2U ||
            independent_rank >=
                component_count_ - 1U) {
            throw std::out_of_range(
                "mpmc::flow::NaturalVariableLayout2P: independent composition index out of range");
        }
        const std::size_t dependent =
            dependent_composition_component(
                phase);
        return independent_rank < dependent
            ? independent_rank
            : independent_rank + 1U;
    }

    [[nodiscard]] std::optional<std::size_t>
    independent_composition_unknown_index(
        std::size_t phase,
        std::size_t component) const {
        if (phase >= 2U ||
            component >= component_count_) {
            throw std::out_of_range(
                "mpmc::flow::NaturalVariableLayout2P: phase/component index out of range");
        }
        const std::size_t dependent =
            dependent_composition_component(
                phase);
        if (component == dependent) {
            return std::nullopt;
        }
        const std::size_t rank =
            component < dependent
                ? component
                : component - 1U;
        return 3U +
            phase *
                (component_count_ - 1U) +
            rank;
    }

    [[nodiscard]] std::size_t
    component_conservation_equation_index(
        std::size_t component) const {
        if (component >= component_count_) {
            throw std::out_of_range(
                "mpmc::flow::NaturalVariableLayout2P: component equation out of range");
        }
        return component;
    }

    [[nodiscard]] std::size_t
    energy_equation_index() const noexcept {
        return component_count_;
    }

    [[nodiscard]] std::size_t
    fugacity_equilibrium_equation_index(
        std::size_t component) const {
        if (component >= component_count_) {
            throw std::out_of_range(
                "mpmc::flow::NaturalVariableLayout2P: fugacity component out of range");
        }
        return component_count_ + 1U +
            component;
    }

    [[nodiscard]] NaturalVariableLayoutDescriptor
    descriptor() const {
        return {
            component_count_,
            2U,
            {
                dependent_composition_component(0U),
                dependent_composition_component(1U)}};
    }

private:
    std::size_t component_count_{};
    NaturalVariableCompositionPivot2P
        composition_pivot_;
};

struct NaturalVariableCellStateInput2P {
    std::vector<std::string> component_ids;
    double reference_pressure_pa{};
    double temperature_k{};
    double independent_saturation{};
    std::array<std::vector<double>, 2>
        independent_phase_compositions;
    std::optional<
        NaturalVariableCompositionPivot2P>
        composition_pivot;
    std::array<PhasePropertyPrerequisiteInput, 2>
        phase_properties;
};

class NaturalVariableCellState2P {
public:
    [[nodiscard]] static NaturalVariableCellState2P
    create(NaturalVariableCellStateInput2P input) {
        NaturalVariableLayout2P layout{
            input.composition_pivot
                ? *input.composition_pivot
                : NaturalVariableCompositionPivot2P::
                      fixed_last(
                          input.component_ids.size())};
        if (layout.component_count() !=
            input.component_ids.size()) {
            throw std::invalid_argument(
                "mpmc::flow::NaturalVariableCellState2P: pivot/component count mismatch");
        }
        for (std::size_t component = 0U;
             component < input.component_ids.size();
             ++component) {
            if (input.component_ids[component].empty()) {
                throw std::invalid_argument(
                    "mpmc::flow::NaturalVariableCellState2P: component id must not be empty");
            }
            for (std::size_t previous = 0U;
                 previous < component;
                 ++previous) {
                if (input.component_ids[previous] ==
                    input.component_ids[component]) {
                    throw std::invalid_argument(
                        "mpmc::flow::NaturalVariableCellState2P: component ids must be unique and ordered");
                }
            }
        }

        natural_variable_detail::
            require_finite_positive(
                input.reference_pressure_pa,
                "reference pressure [Pa]");
        natural_variable_detail::
            require_finite_positive(
                input.temperature_k,
                "temperature [K]");
        if (!std::isfinite(
                input.independent_saturation) ||
            !(input.independent_saturation > 0.0) ||
            !(input.independent_saturation < 1.0)) {
            throw std::invalid_argument(
                "mpmc::flow::NaturalVariableCellState2P: S0 must be finite and strictly inside (0,1)");
        }

        std::array<double, 2> saturation{
            input.independent_saturation,
            1.0 -
                input.independent_saturation};
        std::array<std::vector<double>, 2>
            composition;
        std::array<PhasePropertyPayload, 2>
            properties;

        for (std::size_t phase = 0U;
             phase < 2U;
             ++phase) {
            composition[phase] =
                natural_variable_detail::
                    reconstruct_positive_composition(
                        input.independent_phase_compositions[
                            phase],
                        layout.component_count() - 1U,
                        phase,
                        layout
                            .dependent_composition_component(
                                phase));
            properties[phase] =
                natural_variable_detail::
                    validate_phase_properties(
                        input.phase_properties[phase],
                        phase);
        }

        return NaturalVariableCellState2P{
            std::move(layout),
            std::move(input.component_ids),
            input.reference_pressure_pa,
            input.temperature_k,
            saturation,
            std::move(composition),
            std::move(properties)};
    }

    [[nodiscard]] const NaturalVariableLayout2P&
    layout() const noexcept {
        return layout_;
    }

    [[nodiscard]] std::span<const std::string>
    component_ids() const noexcept {
        return component_ids_;
    }

    [[nodiscard]] double
    reference_pressure_pa() const noexcept {
        return reference_pressure_pa_;
    }

    [[nodiscard]] double
    temperature_k() const noexcept {
        return temperature_k_;
    }

    [[nodiscard]] double
    phase_saturation(
        std::size_t phase) const {
        return saturation_.at(phase);
    }

    [[nodiscard]] std::span<const double>
    phase_composition(
        std::size_t phase) const {
        return composition_.at(phase);
    }

    [[nodiscard]] const PhasePropertyPayload&
    phase_properties(
        std::size_t phase) const {
        return properties_.at(phase);
    }

private:
    NaturalVariableCellState2P(
        NaturalVariableLayout2P layout,
        std::vector<std::string> component_ids,
        double reference_pressure_pa,
        double temperature_k,
        std::array<double, 2> saturation,
        std::array<std::vector<double>, 2>
            composition,
        std::array<PhasePropertyPayload, 2>
            properties)
        : layout_(std::move(layout)),
          component_ids_(
              std::move(component_ids)),
          reference_pressure_pa_(
              reference_pressure_pa),
          temperature_k_(temperature_k),
          saturation_(saturation),
          composition_(
              std::move(composition)),
          properties_(
              std::move(properties)) {}

    NaturalVariableLayout2P layout_;
    std::vector<std::string> component_ids_;
    double reference_pressure_pa_{};
    double temperature_k_{};
    std::array<double, 2> saturation_{};
    std::array<std::vector<double>, 2>
        composition_;
    std::array<PhasePropertyPayload, 2>
        properties_{};
};

} // namespace mpmc::flow

#endif // MPMC_FLOW_TWO_PHASE_CELL_STATE_HPP
