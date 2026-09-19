#ifndef MPMC_FLOW_NATURAL_VARIABLE_CELL_STATE_HPP
#define MPMC_FLOW_NATURAL_VARIABLE_CELL_STATE_HPP

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
#include <vector>

namespace mpmc::flow {

inline constexpr std::size_t fixed_three_phase_count = 3U;

/// Generic numerical phase slots for the fixed-three-phase interior contract.
///
/// These slots intentionally carry no oil/gas/water or liquid/vapor identity.
enum class PhaseSlot3 : std::uint8_t {
    phase0 = 0,
    phase1 = 1,
    phase2 = 2,
};

/// Input form for phase properties required by the first non-isothermal
/// natural-variable flow contract.
///
/// Optional values distinguish an unsupported/missing property from a supplied
/// but numerically invalid property. The validated cell state never stores
/// optionals: all five properties are mandatory before a state is accepted.
struct PhasePropertyPrerequisiteInput {
    std::optional<double> molar_density_mol_per_m3;
    std::optional<double> mass_density_kg_per_m3;
    std::optional<double> dynamic_viscosity_pa_s;
    std::optional<double> specific_enthalpy_j_per_kg;
    std::optional<double> specific_internal_energy_j_per_kg;
};

/// Fully validated phase-property payload.
struct PhasePropertyPayload {
    double molar_density_mol_per_m3;
    double mass_density_kg_per_m3;
    double dynamic_viscosity_pa_s;
    double specific_enthalpy_j_per_kg;
    double specific_internal_energy_j_per_kg;
};

namespace natural_variable_detail {

[[nodiscard]] inline std::size_t checked_phase_index(
    PhaseSlot3 slot) {
    const auto index =
        static_cast<std::size_t>(slot);
    if (index >= fixed_three_phase_count) {
        throw std::invalid_argument(
            "mpmc::flow: invalid fixed-three-phase slot");
    }
    return index;
}

inline void require_finite_positive(
    double value,
    std::string_view name) {
    if (!std::isfinite(value) || value <= 0.0) {
        throw std::invalid_argument(
            std::string{"mpmc::flow: "} +
            std::string{name} +
            " must be finite and strictly positive");
    }
}

[[nodiscard]] inline double require_property(
    const std::optional<double>& value,
    std::string_view name,
    std::size_t phase_index,
    bool strictly_positive) {
    if (!value.has_value()) {
        throw std::invalid_argument(
            std::string{
                "mpmc::flow: unsupported required phase property '"} +
            std::string{name} +
            "' for phase" +
            std::to_string(phase_index));
    }
    if (!std::isfinite(*value) ||
        (strictly_positive && *value <= 0.0)) {
        throw std::invalid_argument(
            std::string{"mpmc::flow: phase"} +
            std::to_string(phase_index) +
            " property '" +
            std::string{name} +
            (strictly_positive
                 ? "' must be finite and strictly positive"
                 : "' must be finite"));
    }
    return *value;
}

[[nodiscard]] inline PhasePropertyPayload
validate_phase_properties(
    const PhasePropertyPrerequisiteInput& input,
    std::size_t phase_index) {
    return PhasePropertyPayload{
        require_property(
            input.molar_density_mol_per_m3,
            "molar_density_mol_per_m3",
            phase_index,
            true),
        require_property(
            input.mass_density_kg_per_m3,
            "mass_density_kg_per_m3",
            phase_index,
            true),
        require_property(
            input.dynamic_viscosity_pa_s,
            "dynamic_viscosity_pa_s",
            phase_index,
            true),
        require_property(
            input.specific_enthalpy_j_per_kg,
            "specific_enthalpy_j_per_kg",
            phase_index,
            false),
        require_property(
            input.specific_internal_energy_j_per_kg,
            "specific_internal_energy_j_per_kg",
            phase_index,
            false)};
}

[[nodiscard]] inline std::vector<double>
reconstruct_positive_composition(
    std::span<const double> independent,
    std::size_t expected_independent_count,
    std::size_t phase_index) {
    if (independent.size() !=
        expected_independent_count) {
        throw std::invalid_argument(
            std::string{"mpmc::flow: phase"} +
            std::to_string(phase_index) +
            " independent composition size mismatch");
    }

    std::vector<double> full;
    full.reserve(expected_independent_count + 1U);
    double sum = 0.0;
    for (const double value : independent) {
        if (!std::isfinite(value) ||
            value <= 0.0) {
            throw std::invalid_argument(
                std::string{"mpmc::flow: phase"} +
                std::to_string(phase_index) +
                " positive-support composition entries must be finite and strictly positive");
        }
        sum += value;
        if (!std::isfinite(sum)) {
            throw std::invalid_argument(
                "mpmc::flow: composition sum is non-finite");
        }
        full.push_back(value);
    }

    const double final_value = 1.0 - sum;
    if (!std::isfinite(final_value) ||
        final_value <= 0.0) {
        throw std::invalid_argument(
            std::string{"mpmc::flow: phase"} +
            std::to_string(phase_index) +
            " reconstructed composition must be strictly positive");
    }
    full.push_back(final_value);
    return full;
}

} // namespace natural_variable_detail

/// Deterministic unknown/equation indexing for the fixed-three-phase
/// non-isothermal natural-variable formulation.
///
/// For Nc components:
///   unknowns  = p_ref + T + 2 independent saturations
///             + 3*(Nc-1) independent phase compositions
///             = 3*Nc + 1
///   equations = Nc component balances + 1 energy balance
///             + 2*Nc fugacity-equality rows
///             = 3*Nc + 1
///
/// phase0 is the local fugacity-reference phase. phase2 saturation and the
/// final component fraction of every phase are dependent quantities.
class NaturalVariableLayout3P {
public:
    explicit NaturalVariableLayout3P(
        std::size_t component_count)
        : component_count_(component_count) {
        if (component_count_ < 2U) {
            throw std::invalid_argument(
                "mpmc::flow::NaturalVariableLayout3P: at least two components are required");
        }
        if (component_count_ >
            (std::numeric_limits<std::size_t>::max() - 1U) /
                fixed_three_phase_count) {
            throw std::length_error(
                "mpmc::flow::NaturalVariableLayout3P: component count overflows layout size");
        }
    }

    [[nodiscard]] std::size_t
    component_count() const noexcept {
        return component_count_;
    }

    [[nodiscard]] static constexpr std::size_t
    phase_count() noexcept {
        return fixed_three_phase_count;
    }

    [[nodiscard]] static constexpr PhaseSlot3
    fugacity_reference_phase() noexcept {
        return PhaseSlot3::phase0;
    }

    [[nodiscard]] std::size_t
    unknown_count() const noexcept {
        return fixed_three_phase_count *
                   component_count_ +
               1U;
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

    [[nodiscard]] std::optional<std::size_t>
    independent_saturation_unknown_index(
        PhaseSlot3 slot) const {
        const auto phase =
            natural_variable_detail::
                checked_phase_index(slot);
        if (phase == 2U) {
            return std::nullopt;
        }
        return 2U + phase;
    }

    [[nodiscard]] std::optional<std::size_t>
    independent_composition_unknown_index(
        PhaseSlot3 slot,
        std::size_t component) const {
        const auto phase =
            natural_variable_detail::
                checked_phase_index(slot);
        if (component >= component_count_) {
            throw std::out_of_range(
                "mpmc::flow::NaturalVariableLayout3P: component index out of range");
        }
        if (component + 1U ==
            component_count_) {
            return std::nullopt;
        }
        return 4U +
               phase * (component_count_ - 1U) +
               component;
    }

    [[nodiscard]] std::size_t
    component_conservation_equation_index(
        std::size_t component) const {
        if (component >= component_count_) {
            throw std::out_of_range(
                "mpmc::flow::NaturalVariableLayout3P: component equation index out of range");
        }
        return component;
    }

    [[nodiscard]] std::size_t
    energy_equation_index() const noexcept {
        return component_count_;
    }

    [[nodiscard]] std::size_t
    fugacity_equilibrium_equation_index(
        PhaseSlot3 non_reference_phase,
        std::size_t component) const {
        const auto phase =
            natural_variable_detail::
                checked_phase_index(
                    non_reference_phase);
        if (phase == 0U) {
            throw std::invalid_argument(
                "mpmc::flow::NaturalVariableLayout3P: fugacity reference phase has no equilibrium row against itself");
        }
        if (component >= component_count_) {
            throw std::out_of_range(
                "mpmc::flow::NaturalVariableLayout3P: fugacity component index out of range");
        }
        return component_count_ + 1U +
               (phase - 1U) * component_count_ +
               component;
    }

private:
    std::size_t component_count_;
};

struct NaturalVariableCellStateInput3P {
    std::vector<std::string> component_ids;
    double reference_pressure_pa = 0.0;
    double temperature_k = 0.0;
    std::array<double, 2> independent_saturations{};
    std::array<std::vector<double>, 3>
        independent_phase_compositions;
    std::array<PhasePropertyPrerequisiteInput, 3>
        phase_properties;
};

/// Immutable validated local state for the first flow implementation slice.
///
/// This type is deliberately local and mesh-free. It does not evaluate an EOS,
/// residual, flux, time step, phase switch, Newton update, PETSc object or well
/// model. With capillary pressure not yet implemented, all three phase
/// pressures equal reference_pressure_pa().
class NaturalVariableCellState3P {
public:
    [[nodiscard]] static NaturalVariableCellState3P
    create(NaturalVariableCellStateInput3P input) {
        NaturalVariableLayout3P layout{
            input.component_ids.size()};

        for (std::size_t index = 0U;
             index < input.component_ids.size();
             ++index) {
            if (input.component_ids[index].empty()) {
                throw std::invalid_argument(
                    "mpmc::flow::NaturalVariableCellState3P: component id must not be empty");
            }
            for (std::size_t other = 0U;
                 other < index;
                 ++other) {
                if (input.component_ids[other] ==
                    input.component_ids[index]) {
                    throw std::invalid_argument(
                        "mpmc::flow::NaturalVariableCellState3P: component ids must be unique and ordered");
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

        const double saturation0 =
            input.independent_saturations[0];
        const double saturation1 =
            input.independent_saturations[1];
        if (!std::isfinite(saturation0) ||
            !std::isfinite(saturation1) ||
            saturation0 <= 0.0 ||
            saturation1 <= 0.0) {
            throw std::invalid_argument(
                "mpmc::flow::NaturalVariableCellState3P: positive-support independent saturations must be finite and strictly positive");
        }
        const double saturation2 =
            1.0 - saturation0 - saturation1;
        if (!std::isfinite(saturation2) ||
            saturation2 <= 0.0) {
            throw std::invalid_argument(
                "mpmc::flow::NaturalVariableCellState3P: reconstructed phase2 saturation must be strictly positive");
        }

        std::array<std::vector<double>, 3>
            compositions;
        std::array<PhasePropertyPayload, 3>
            properties{};
        for (std::size_t phase = 0U;
             phase < fixed_three_phase_count;
             ++phase) {
            compositions[phase] =
                natural_variable_detail::
                    reconstruct_positive_composition(
                        input
                            .independent_phase_compositions
                            [phase],
                        layout.component_count() - 1U,
                        phase);
            properties[phase] =
                natural_variable_detail::
                    validate_phase_properties(
                        input.phase_properties[phase],
                        phase);
        }

        return NaturalVariableCellState3P{
            std::move(layout),
            std::move(input.component_ids),
            input.reference_pressure_pa,
            input.temperature_k,
            std::array<double, 3>{
                saturation0,
                saturation1,
                saturation2},
            std::move(compositions),
            std::move(properties)};
    }

    [[nodiscard]] const NaturalVariableLayout3P&
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
    phase_pressure_pa(
        PhaseSlot3 slot) const {
        (void)natural_variable_detail::
            checked_phase_index(slot);
        return reference_pressure_pa_;
    }

    [[nodiscard]] std::array<double, 3>
    phase_pressures_pa() const noexcept {
        return {
            reference_pressure_pa_,
            reference_pressure_pa_,
            reference_pressure_pa_};
    }

    [[nodiscard]] double
    phase_saturation(
        PhaseSlot3 slot) const {
        return saturations_[
            natural_variable_detail::
                checked_phase_index(slot)];
    }

    [[nodiscard]] std::span<const double>
    phase_composition(
        PhaseSlot3 slot) const {
        return compositions_[
            natural_variable_detail::
                checked_phase_index(slot)];
    }

    [[nodiscard]] const PhasePropertyPayload&
    phase_properties(
        PhaseSlot3 slot) const {
        return properties_[
            natural_variable_detail::
                checked_phase_index(slot)];
    }

private:
    NaturalVariableCellState3P(
        NaturalVariableLayout3P layout,
        std::vector<std::string> component_ids,
        double reference_pressure_pa,
        double temperature_k,
        std::array<double, 3> saturations,
        std::array<std::vector<double>, 3>
            compositions,
        std::array<PhasePropertyPayload, 3>
            properties)
        : layout_(std::move(layout)),
          component_ids_(std::move(component_ids)),
          reference_pressure_pa_(
              reference_pressure_pa),
          temperature_k_(temperature_k),
          saturations_(std::move(saturations)),
          compositions_(std::move(compositions)),
          properties_(std::move(properties)) {}

    NaturalVariableLayout3P layout_;
    std::vector<std::string> component_ids_;
    double reference_pressure_pa_;
    double temperature_k_;
    std::array<double, 3> saturations_;
    std::array<std::vector<double>, 3>
        compositions_;
    std::array<PhasePropertyPayload, 3>
        properties_;
};

} // namespace mpmc::flow

#endif // MPMC_FLOW_NATURAL_VARIABLE_CELL_STATE_HPP
