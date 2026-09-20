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


class NaturalVariableLayout3P;

class NaturalVariableCompositionPivotDescriptor {
public:
    NaturalVariableCompositionPivotDescriptor() = default;

    explicit NaturalVariableCompositionPivotDescriptor(
        std::vector<std::size_t> dependent_components)
        : dependent_components_(
              std::move(dependent_components)) {}

    [[nodiscard]] const std::vector<std::size_t>&
    dependent_components() const noexcept {
        return dependent_components_;
    }

private:
    std::vector<std::size_t> dependent_components_;
};

/// Phase-cardinality-neutral local natural-variable layout descriptor.
///
/// This is metadata only. It preserves block width/equation positions across
/// 1/2/3-phase reduced systems without exposing phase-specific coordinates.
class NaturalVariableLayoutDescriptor {
public:
    NaturalVariableLayoutDescriptor(
        std::size_t component_count,
        std::size_t phase_count,
        std::vector<std::size_t> dependent_components)
        : component_count_(component_count),
          phase_count_(phase_count),
          composition_pivot_(
              std::move(dependent_components)) {
        if (component_count_ < 2U ||
            phase_count_ == 0U ||
            phase_count_ > fixed_three_phase_count ||
            composition_pivot_.dependent_components().size() !=
                phase_count_) {
            throw std::invalid_argument(
                "mpmc::flow::NaturalVariableLayoutDescriptor: invalid component/phase cardinality");
        }
        for (const auto dependent :
             composition_pivot_.dependent_components()) {
            if (dependent >= component_count_) {
                throw std::invalid_argument(
                    "mpmc::flow::NaturalVariableLayoutDescriptor: dependent component out of range");
            }
        }
        if (component_count_ >
            (std::numeric_limits<std::size_t>::max() - 1U) /
                phase_count_) {
            throw std::length_error(
                "mpmc::flow::NaturalVariableLayoutDescriptor: layout size overflow");
        }
    }

    NaturalVariableLayoutDescriptor(
        const NaturalVariableLayout3P& layout);

    [[nodiscard]] std::size_t
    component_count() const noexcept {
        return component_count_;
    }

    [[nodiscard]] std::size_t
    phase_count() const noexcept {
        return phase_count_;
    }

    [[nodiscard]] std::size_t
    unknown_count() const noexcept {
        return phase_count_ *
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

    [[nodiscard]] std::size_t
    component_conservation_equation_index(
        std::size_t component) const {
        if (component >= component_count_) {
            throw std::out_of_range(
                "mpmc::flow::NaturalVariableLayoutDescriptor: component equation out of range");
        }
        return component;
    }

    [[nodiscard]] std::size_t
    energy_equation_index() const noexcept {
        return component_count_;
    }

    [[nodiscard]] const NaturalVariableCompositionPivotDescriptor&
    composition_pivot() const noexcept {
        return composition_pivot_;
    }

private:
    std::size_t component_count_{};
    std::size_t phase_count_{};
    NaturalVariableCompositionPivotDescriptor
        composition_pivot_;
};


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
    std::size_t phase_index,
    std::size_t dependent_component) {
    if (independent.size() !=
        expected_independent_count) {
        throw std::invalid_argument(
            std::string{"mpmc::flow: phase"} +
            std::to_string(phase_index) +
            " independent composition size mismatch");
    }
    const std::size_t component_count =
        expected_independent_count + 1U;
    if (dependent_component >= component_count) {
        throw std::invalid_argument(
            "mpmc::flow: dependent composition component is out of range");
    }

    std::vector<double> full(
        component_count,
        0.0);
    double sum = 0.0;
    std::size_t independent_index = 0U;
    for (std::size_t component = 0U;
         component < component_count;
         ++component) {
        if (component == dependent_component) {
            continue;
        }
        const double value =
            independent[independent_index++];
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
        full[component] = value;
    }

    const double dependent_value =
        1.0 - sum;
    if (!std::isfinite(dependent_value) ||
        dependent_value <= 0.0) {
        throw std::invalid_argument(
            std::string{"mpmc::flow: phase"} +
            std::to_string(phase_index) +
            " reconstructed dependent composition must be strictly positive");
    }
    full[dependent_component] =
        dependent_value;
    return full;
}

} // namespace natural_variable_detail

/// Frozen composition-coordinate pivot for one local three-phase chart.
///
/// Selection uses the largest current mole fraction in each phase. Exact ties
/// select the smallest component index. The selected indices are immutable:
/// callers must explicitly build a new pivot before changing coordinates.
class NaturalVariableCompositionPivot3P {
public:
    [[nodiscard]] static NaturalVariableCompositionPivot3P
    fixed_last(std::size_t component_count) {
        validate_component_count(component_count);
        return NaturalVariableCompositionPivot3P{
            component_count,
            {component_count - 1U,
             component_count - 1U,
             component_count - 1U}};
    }

    [[nodiscard]] static NaturalVariableCompositionPivot3P
    from_dependent_components(
        std::size_t component_count,
        std::array<std::size_t, 3>
            dependent_components) {
        validate_component_count(component_count);
        for (const auto component :
             dependent_components) {
            if (component >= component_count) {
                throw std::invalid_argument(
                    "mpmc::flow::NaturalVariableCompositionPivot3P: dependent component out of range");
            }
        }
        return NaturalVariableCompositionPivot3P{
            component_count,
            dependent_components};
    }

    [[nodiscard]] static NaturalVariableCompositionPivot3P
    select(
        const std::array<
            std::span<const double>,
            3>& phase_compositions) {
        const std::size_t component_count =
            phase_compositions[0].size();
        validate_component_count(component_count);
        std::array<std::size_t, 3>
            dependent_components{};

        const long double tolerance =
            4096.0L *
            static_cast<long double>(
                std::numeric_limits<double>::epsilon()) *
            static_cast<long double>(
                component_count);

        for (std::size_t phase = 0U;
             phase < fixed_three_phase_count;
             ++phase) {
            const auto composition =
                phase_compositions[phase];
            if (composition.size() !=
                component_count) {
                throw std::invalid_argument(
                    "mpmc::flow::NaturalVariableCompositionPivot3P: phase composition sizes must match");
            }

            long double sum = 0.0L;
            std::size_t best = 0U;
            double best_value = -1.0;
            for (std::size_t component = 0U;
                 component < component_count;
                 ++component) {
                const double value =
                    composition[component];
                if (!std::isfinite(value) ||
                    value <= 0.0) {
                    throw std::invalid_argument(
                        "mpmc::flow::NaturalVariableCompositionPivot3P: positive-support finite compositions required");
                }
                sum +=
                    static_cast<long double>(
                        value);
                if (value > best_value) {
                    best = component;
                    best_value = value;
                }
            }
            if (!std::isfinite(sum) ||
                std::abs(sum - 1.0L) >
                    tolerance) {
                throw std::invalid_argument(
                    "mpmc::flow::NaturalVariableCompositionPivot3P: phase composition must be normalized");
            }
            dependent_components[phase] =
                best;
        }

        return NaturalVariableCompositionPivot3P{
            component_count,
            dependent_components};
    }

    [[nodiscard]] static NaturalVariableCompositionPivot3P
    select(
        const std::array<
            std::vector<double>,
            3>& phase_compositions) {
        return select({
            std::span<const double>{
                phase_compositions[0]},
            std::span<const double>{
                phase_compositions[1]},
            std::span<const double>{
                phase_compositions[2]}});
    }

    [[nodiscard]] std::size_t
    component_count() const noexcept {
        return component_count_;
    }

    [[nodiscard]] std::size_t
    dependent_component(
        PhaseSlot3 slot) const {
        return dependent_components_[
            natural_variable_detail::
                checked_phase_index(slot)];
    }

    [[nodiscard]] const std::array<
        std::size_t,
        3>&
    dependent_components() const noexcept {
        return dependent_components_;
    }

private:
    NaturalVariableCompositionPivot3P(
        std::size_t component_count,
        std::array<std::size_t, 3>
            dependent_components)
        : component_count_(component_count),
          dependent_components_(
              dependent_components) {}

    static void validate_component_count(
        std::size_t component_count) {
        if (component_count < 2U) {
            throw std::invalid_argument(
                "mpmc::flow::NaturalVariableCompositionPivot3P: at least two components are required");
        }
    }

    std::size_t component_count_;
    std::array<std::size_t, 3>
        dependent_components_;
};

struct NaturalVariableCompositionUnknownIdentity3P {
    PhaseSlot3 phase{PhaseSlot3::phase0};
    std::size_t component{};
    friend bool operator==(
        const NaturalVariableCompositionUnknownIdentity3P&,
        const NaturalVariableCompositionUnknownIdentity3P&) =
        default;
};

struct NaturalVariableFugacityRowIdentity3P {
    PhaseSlot3 non_reference_phase{
        PhaseSlot3::phase1};
    std::size_t component{};
    friend bool operator==(
        const NaturalVariableFugacityRowIdentity3P&,
        const NaturalVariableFugacityRowIdentity3P&) =
        default;
};

/// Deterministic unknown/equation indexing for the fixed-three-phase
/// non-isothermal natural-variable formulation.
///
/// The component order is never changed by pivoting. Within each phase block,
/// independent composition columns remain in canonical component-index order
/// with the frozen dependent component omitted. Fugacity-equilibrium rows
/// remain in canonical component order and are independent of the pivot.
class NaturalVariableLayout3P {
public:
    explicit NaturalVariableLayout3P(
        std::size_t component_count)
        : NaturalVariableLayout3P(
              NaturalVariableCompositionPivot3P::
                  fixed_last(component_count)) {}

    explicit NaturalVariableLayout3P(
        NaturalVariableCompositionPivot3P pivot)
        : component_count_(
              pivot.component_count()),
          composition_pivot_(
              std::move(pivot)) {
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

    [[nodiscard]] const NaturalVariableCompositionPivot3P&
    composition_pivot() const noexcept {
        return composition_pivot_;
    }

    [[nodiscard]] std::size_t
    dependent_composition_component(
        PhaseSlot3 slot) const {
        return composition_pivot_.
            dependent_component(slot);
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

    [[nodiscard]] std::size_t
    independent_composition_component(
        PhaseSlot3 slot,
        std::size_t independent_rank) const {
        if (independent_rank >=
            component_count_ - 1U) {
            throw std::out_of_range(
                "mpmc::flow::NaturalVariableLayout3P: independent composition rank out of range");
        }
        const std::size_t dependent =
            dependent_composition_component(
                slot);
        return independent_rank < dependent
            ? independent_rank
            : independent_rank + 1U;
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
        const std::size_t dependent =
            dependent_composition_component(
                slot);
        if (component == dependent) {
            return std::nullopt;
        }
        const std::size_t independent_rank =
            component < dependent
                ? component
                : component - 1U;
        return 4U +
               phase *
                   (component_count_ - 1U) +
               independent_rank;
    }

    [[nodiscard]] std::optional<
        NaturalVariableCompositionUnknownIdentity3P>
    composition_unknown_identity(
        std::size_t unknown_index) const {
        if (unknown_index >=
            unknown_count()) {
            throw std::out_of_range(
                "mpmc::flow::NaturalVariableLayout3P: unknown index out of range");
        }
        if (unknown_index < 4U) {
            return std::nullopt;
        }
        const std::size_t local =
            unknown_index - 4U;
        const std::size_t block =
            component_count_ - 1U;
        const std::size_t phase =
            local / block;
        const std::size_t rank =
            local % block;
        return NaturalVariableCompositionUnknownIdentity3P{
            static_cast<PhaseSlot3>(phase),
            independent_composition_component(
                static_cast<PhaseSlot3>(phase),
                rank)};
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
               (phase - 1U) *
                   component_count_ +
               component;
    }

    [[nodiscard]] std::optional<
        NaturalVariableFugacityRowIdentity3P>
    fugacity_equilibrium_row_identity(
        std::size_t equation_index) const {
        if (equation_index >=
            equation_count()) {
            throw std::out_of_range(
                "mpmc::flow::NaturalVariableLayout3P: equation index out of range");
        }
        const std::size_t first =
            component_count_ + 1U;
        if (equation_index < first) {
            return std::nullopt;
        }
        const std::size_t local =
            equation_index - first;
        const std::size_t phase_offset =
            local / component_count_;
        const std::size_t component =
            local % component_count_;
        return NaturalVariableFugacityRowIdentity3P{
            static_cast<PhaseSlot3>(
                phase_offset + 1U),
            component};
    }

private:
    std::size_t component_count_;
    NaturalVariableCompositionPivot3P
        composition_pivot_;
};


inline NaturalVariableLayoutDescriptor::
NaturalVariableLayoutDescriptor(
    const NaturalVariableLayout3P& layout)
    : NaturalVariableLayoutDescriptor(
          layout.component_count(),
          NaturalVariableLayout3P::phase_count(),
          std::vector<std::size_t>{
              layout.composition_pivot()
                  .dependent_components()
                  .begin(),
              layout.composition_pivot()
                  .dependent_components()
                  .end()}) {}

struct NaturalVariableCellStateInput3P {
    std::vector<std::string> component_ids;
    double reference_pressure_pa = 0.0;
    double temperature_k = 0.0;
    std::array<double, 2> independent_saturations{};
    std::array<std::vector<double>, 3>
        independent_phase_compositions;
    std::optional<
        NaturalVariableCompositionPivot3P>
        composition_pivot;
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
            input.composition_pivot
                ? *input.composition_pivot
                : NaturalVariableCompositionPivot3P::
                      fixed_last(
                          input.component_ids.size())};
        if (layout.component_count() !=
            input.component_ids.size()) {
            throw std::invalid_argument(
                "mpmc::flow::NaturalVariableCellState3P: pivot/component count mismatch");
        }

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
                        phase,
                        layout.dependent_composition_component(
                            static_cast<PhaseSlot3>(
                                phase)));
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
