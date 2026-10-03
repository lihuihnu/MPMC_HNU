#ifndef MPMC_FLOW_SINGLE_PHASE_CELL_STATE_HPP
#define MPMC_FLOW_SINGLE_PHASE_CELL_STATE_HPP

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
    single_phase_natural_variable_convention =
        "flow/natural-variable/single-phase-reduction/v1";

class NaturalVariableCompositionPivot1P {
public:
    [[nodiscard]] static NaturalVariableCompositionPivot1P
    fixed_last(std::size_t component_count) {
        validate_component_count(component_count);
        return {
            component_count,
            component_count - 1U};
    }

    [[nodiscard]] static NaturalVariableCompositionPivot1P
    from_dependent_component(
        std::size_t component_count,
        std::size_t dependent_component) {
        validate_component_count(component_count);
        if (dependent_component >= component_count) {
            throw std::invalid_argument(
                "mpmc::flow::NaturalVariableCompositionPivot1P: dependent component out of range");
        }
        return {
            component_count,
            dependent_component};
    }

    [[nodiscard]] static NaturalVariableCompositionPivot1P
    select(std::span<const double> composition) {
        validate_component_count(
            composition.size());
        long double sum = 0.0L;
        std::size_t best = 0U;
        double best_value = -1.0;
        for (std::size_t component = 0U;
             component < composition.size();
             ++component) {
            const double value =
                composition[component];
            if (!std::isfinite(value) ||
                !(value > 0.0)) {
                throw std::invalid_argument(
                    "mpmc::flow::NaturalVariableCompositionPivot1P: positive finite composition required");
            }
            sum += static_cast<long double>(
                value);
            if (value > best_value) {
                best = component;
                best_value = value;
            }
        }
        const long double tolerance =
            4096.0L *
            static_cast<long double>(
                std::numeric_limits<double>::epsilon()) *
            static_cast<long double>(
                composition.size());
        if (!std::isfinite(sum) ||
            std::abs(sum - 1.0L) >
                tolerance) {
            throw std::invalid_argument(
                "mpmc::flow::NaturalVariableCompositionPivot1P: composition must be normalized");
        }
        return {
            composition.size(),
            best};
    }

    [[nodiscard]] std::size_t
    component_count() const noexcept {
        return component_count_;
    }

    [[nodiscard]] std::size_t
    dependent_component() const noexcept {
        return dependent_component_;
    }

private:
    NaturalVariableCompositionPivot1P(
        std::size_t component_count,
        std::size_t dependent_component)
        : component_count_(component_count),
          dependent_component_(
              dependent_component) {}

    static void validate_component_count(
        std::size_t component_count) {
        if (component_count < 2U) {
            throw std::invalid_argument(
                "mpmc::flow::NaturalVariableCompositionPivot1P: at least two components are required");
        }
    }

    std::size_t component_count_{};
    std::size_t dependent_component_{};
};

class NaturalVariableLayout1P {
public:
    explicit NaturalVariableLayout1P(
        std::size_t component_count)
        : NaturalVariableLayout1P(
              NaturalVariableCompositionPivot1P::
                  fixed_last(component_count)) {}

    explicit NaturalVariableLayout1P(
        NaturalVariableCompositionPivot1P pivot)
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
        return 1U;
    }

    [[nodiscard]] std::size_t
    unknown_count() const noexcept {
        return component_count_ + 1U;
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

    [[nodiscard]] std::size_t
    dependent_composition_component() const noexcept {
        return composition_pivot_.
            dependent_component();
    }

    [[nodiscard]] const NaturalVariableCompositionPivot1P&
    composition_pivot() const noexcept {
        return composition_pivot_;
    }

    [[nodiscard]] std::size_t
    independent_composition_component(
        std::size_t independent_rank) const {
        if (independent_rank >=
            component_count_ - 1U) {
            throw std::out_of_range(
                "mpmc::flow::NaturalVariableLayout1P: independent composition rank out of range");
        }
        const std::size_t dependent =
            dependent_composition_component();
        return independent_rank < dependent
            ? independent_rank
            : independent_rank + 1U;
    }

    [[nodiscard]] std::optional<std::size_t>
    independent_composition_unknown_index(
        std::size_t component) const {
        if (component >= component_count_) {
            throw std::out_of_range(
                "mpmc::flow::NaturalVariableLayout1P: component index out of range");
        }
        const std::size_t dependent =
            dependent_composition_component();
        if (component == dependent) {
            return std::nullopt;
        }
        const std::size_t rank =
            component < dependent
                ? component
                : component - 1U;
        return 2U + rank;
    }

    [[nodiscard]] std::size_t
    component_conservation_equation_index(
        std::size_t component) const {
        if (component >= component_count_) {
            throw std::out_of_range(
                "mpmc::flow::NaturalVariableLayout1P: component equation out of range");
        }
        return component;
    }

    [[nodiscard]] std::size_t
    energy_equation_index() const noexcept {
        return component_count_;
    }

    [[nodiscard]] NaturalVariableLayoutDescriptor
    descriptor() const {
        return {
            component_count_,
            1U,
            {dependent_composition_component()}};
    }

private:
    std::size_t component_count_{};
    NaturalVariableCompositionPivot1P
        composition_pivot_;
};

struct NaturalVariableCellStateInput1P {
    std::vector<std::string> component_ids;
    double reference_pressure_pa{};
    double temperature_k{};
    std::vector<double>
        independent_composition;
    std::optional<
        NaturalVariableCompositionPivot1P>
        composition_pivot;
    PhasePropertyPrerequisiteInput
        phase_properties;
};

class NaturalVariableCellState1P {
public:
    [[nodiscard]] static NaturalVariableCellState1P
    create(NaturalVariableCellStateInput1P input) {
        NaturalVariableLayout1P layout{
            input.composition_pivot
                ? *input.composition_pivot
                : NaturalVariableCompositionPivot1P::
                      fixed_last(
                          input.component_ids.size())};
        if (layout.component_count() !=
            input.component_ids.size()) {
            throw std::invalid_argument(
                "mpmc::flow::NaturalVariableCellState1P: pivot/component count mismatch");
        }
        for (std::size_t component = 0U;
             component < input.component_ids.size();
             ++component) {
            if (input.component_ids[component].empty()) {
                throw std::invalid_argument(
                    "mpmc::flow::NaturalVariableCellState1P: component id must not be empty");
            }
            for (std::size_t previous = 0U;
                 previous < component;
                 ++previous) {
                if (input.component_ids[previous] ==
                    input.component_ids[component]) {
                    throw std::invalid_argument(
                        "mpmc::flow::NaturalVariableCellState1P: component ids must be unique and ordered");
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

        auto composition =
            natural_variable_detail::
                reconstruct_positive_composition(
                    input.independent_composition,
                    layout.component_count() - 1U,
                    0U,
                    layout
                        .dependent_composition_component());
        auto properties =
            natural_variable_detail::
                validate_phase_properties(
                    input.phase_properties,
                    0U);

        return NaturalVariableCellState1P{
            std::move(layout),
            std::move(input.component_ids),
            input.reference_pressure_pa,
            input.temperature_k,
            std::move(composition),
            std::move(properties)};
    }

    [[nodiscard]] const NaturalVariableLayout1P&
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

    [[nodiscard]] static constexpr double
    phase_saturation() noexcept {
        return 1.0;
    }

    [[nodiscard]] std::span<const double>
    phase_composition() const noexcept {
        return composition_;
    }

    [[nodiscard]] const PhasePropertyPayload&
    phase_properties() const noexcept {
        return properties_;
    }

private:
    NaturalVariableCellState1P(
        NaturalVariableLayout1P layout,
        std::vector<std::string> component_ids,
        double reference_pressure_pa,
        double temperature_k,
        std::vector<double> composition,
        PhasePropertyPayload properties)
        : layout_(std::move(layout)),
          component_ids_(
              std::move(component_ids)),
          reference_pressure_pa_(
              reference_pressure_pa),
          temperature_k_(temperature_k),
          composition_(
              std::move(composition)),
          properties_(
              std::move(properties)) {}

    NaturalVariableLayout1P layout_;
    std::vector<std::string> component_ids_;
    double reference_pressure_pa_{};
    double temperature_k_{};
    std::vector<double> composition_;
    PhasePropertyPayload properties_{};
};

} // namespace mpmc::flow

#endif // MPMC_FLOW_SINGLE_PHASE_CELL_STATE_HPP
