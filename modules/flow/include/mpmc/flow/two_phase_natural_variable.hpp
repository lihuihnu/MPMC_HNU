#ifndef MPMC_FLOW_TWO_PHASE_NATURAL_VARIABLE_HPP
#define MPMC_FLOW_TWO_PHASE_NATURAL_VARIABLE_HPP

#include <mpmc/flow/component_accumulation_time.hpp>
#include <mpmc/flow/energy_accumulation.hpp>
#include <mpmc/flow/phase_potential_upwind.hpp>

#include <algorithm>
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

namespace two_phase_detail {

[[nodiscard]] inline NaturalVariableStateIdentity
make_state_identity(
    const NaturalVariableCellState2P& state) {
    return {
        state.layout().descriptor(),
        std::vector<std::string>{
            state.component_ids().begin(),
            state.component_ids().end()},
        state.reference_pressure_pa(),
        state.temperature_k(),
        {
            state.phase_saturation(0U),
            state.phase_saturation(1U),
            0.0},
        {
            std::vector<double>{
                state.phase_composition(0U).begin(),
                state.phase_composition(0U).end()},
            std::vector<double>{
                state.phase_composition(1U).begin(),
                state.phase_composition(1U).end()},
            std::vector<double>{}}};
}

[[nodiscard]] inline bool
same_layout(
    const NaturalVariableLayoutDescriptor& first,
    const NaturalVariableLayoutDescriptor& second) {
    return first.component_count() ==
               second.component_count() &&
        first.phase_count() ==
            second.phase_count() &&
        first.unknown_count() ==
            second.unknown_count() &&
        first.composition_pivot()
                .dependent_components() ==
            second.composition_pivot()
                .dependent_components();
}

inline void validate_gradient(
    std::span<const double> gradient,
    std::size_t q,
    const char* name) {
    if (gradient.size() != q) {
        throw std::invalid_argument(
            std::string{"mpmc::flow: two-phase "} +
            name +
            " gradient shape mismatch");
    }
    for (double value : gradient) {
        if (!std::isfinite(value)) {
            throw std::invalid_argument(
                std::string{"mpmc::flow: two-phase "} +
                name +
                " gradient contains non-finite derivative");
        }
    }
}

[[nodiscard]] inline double
d_saturation(
    std::size_t phase,
    std::size_t column) {
    if (phase >= 2U) {
        throw std::out_of_range(
            "mpmc::flow: two-phase saturation phase out of range");
    }
    if (column !=
        NaturalVariableLayout2P::
            independent_saturation_unknown_index()) {
        return 0.0;
    }
    return phase == 0U
        ? 1.0
        : -1.0;
}

[[nodiscard]] inline double
d_composition(
    const NaturalVariableLayout2P& layout,
    std::size_t phase,
    std::size_t component,
    std::size_t column) {
    if (phase >= 2U ||
        component >= layout.component_count()) {
        throw std::out_of_range(
            "mpmc::flow: two-phase composition derivative index out of range");
    }
    for (std::size_t candidate = 0U;
         candidate < layout.component_count();
         ++candidate) {
        const auto candidate_column =
            layout.independent_composition_unknown_index(
                phase,
                candidate);
        if (!candidate_column ||
            *candidate_column != column) {
            continue;
        }
        if (candidate == component) {
            return 1.0;
        }
        if (component ==
            layout.dependent_composition_component(
                phase)) {
            return -1.0;
        }
        return 0.0;
    }
    return 0.0;
}

} // namespace two_phase_detail

struct TwoPhaseMolarDensityNaturalVariableLinearization {
    NaturalVariableLayoutDescriptor layout{
        std::size_t{2U},
        std::size_t{2U},
        std::vector<std::size_t>{1U, 1U}};
    std::array<double, 2>
        molar_density_mol_per_m3{};
    std::array<std::vector<double>, 2>
        gradient;
};

struct TwoPhaseTransportNaturalVariableLinearization {
    NaturalVariableStateIdentity state_identity;
    TransportPropertyProvenance mass_density_provenance;
    TransportPropertyProvenance viscosity_provenance;
    std::array<double, 2>
        mass_density_kg_per_m3{};
    std::array<double, 2>
        dynamic_viscosity_pa_s{};
    std::array<double, 2>
        relative_permeability{};
    std::array<std::vector<double>, 2>
        mass_density_gradient;
    std::array<std::vector<double>, 2>
        dynamic_viscosity_gradient;
    std::array<std::vector<double>, 2>
        relative_permeability_gradient;
};

struct TwoPhaseCaloricNaturalVariableLinearization {
    NaturalVariableStateIdentity state_identity;
    TransportPropertyProvenance enthalpy_provenance;
    TransportPropertyProvenance internal_energy_provenance;
    std::array<double, 2>
        specific_enthalpy_j_per_kg{};
    std::array<double, 2>
        specific_internal_energy_j_per_kg{};
    std::array<std::vector<double>, 2>
        specific_enthalpy_gradient;
    std::array<std::vector<double>, 2>
        specific_internal_energy_gradient;
};

struct TwoPhaseRockThermalStorageLinearization {
    NaturalVariableStateIdentity state_identity;
    TransportPropertyProvenance provenance;
    double volumetric_internal_energy_j_per_rock_m3{};
    std::vector<double>
        volumetric_internal_energy_gradient;
};

struct TwoPhaseMobilityLinearization {
    NaturalVariableStateIdentity state_identity;
    std::array<double, 2>
        phase_pressure_pa{};
    std::array<double, 2>
        mass_density_kg_per_m3{};
    std::array<double, 2>
        dynamic_viscosity_pa_s{};
    std::array<double, 2>
        relative_permeability{};
    std::array<double, 2>
        mobility_per_pa_s{};
    std::array<std::vector<double>, 2>
        phase_pressure_gradient;
    std::array<std::vector<double>, 2>
        mass_density_gradient;
    std::array<std::vector<double>, 2>
        mobility_gradient;
};

struct TwoPhaseFugacityEquilibriumLinearization {
    NaturalVariableLayoutDescriptor layout{
        std::size_t{2U},
        std::size_t{2U},
        std::vector<std::size_t>{1U, 1U}};
    std::vector<std::string> component_ids;
    std::vector<double> residual;
    std::size_t input_count{};
    std::vector<double> jacobian;

    [[nodiscard]] std::size_t
    component_count() const noexcept {
        return component_ids.size();
    }

    [[nodiscard]] std::size_t
    residual_count() const noexcept {
        return residual.size();
    }

    [[nodiscard]] double d_residual(
        std::size_t component,
        std::size_t column) const {
        if (component >= component_count() ||
            column >= input_count) {
            throw std::out_of_range(
                "mpmc::flow: two-phase fugacity Jacobian index out of range");
        }
        return jacobian.at(
            component * input_count +
            column);
    }
};

[[nodiscard]] inline
TwoPhaseMolarDensityNaturalVariableLinearization
make_two_phase_molar_density_linearization(
    const NaturalVariableCellState2P& state,
    std::array<std::vector<double>, 2>
        gradient) {
    const std::size_t q =
        state.layout().unknown_count();
    TwoPhaseMolarDensityNaturalVariableLinearization
        result;
    result.layout =
        state.layout().descriptor();
    for (std::size_t phase = 0U;
         phase < 2U;
         ++phase) {
        two_phase_detail::validate_gradient(
            gradient[phase],
            q,
            "molar-density");
        result.molar_density_mol_per_m3[phase] =
            state.phase_properties(phase)
                .molar_density_mol_per_m3;
        result.gradient[phase] =
            std::move(gradient[phase]);
    }
    return result;
}

[[nodiscard]] inline
TwoPhaseTransportNaturalVariableLinearization
make_two_phase_transport_linearization(
    const NaturalVariableCellState2P& state,
    std::array<std::vector<double>, 2>
        mass_density_gradient,
    std::array<std::vector<double>, 2>
        viscosity_gradient,
    std::array<double, 2>
        relative_permeability,
    std::array<std::vector<double>, 2>
        relative_permeability_gradient,
    TransportPropertyProvenance
        mass_density_provenance,
    TransportPropertyProvenance
        viscosity_provenance) {
    phase_transport_detail::validate_provenance(
        mass_density_provenance,
        "two-phase-mass-density");
    phase_transport_detail::validate_provenance(
        viscosity_provenance,
        "two-phase-viscosity");

    TwoPhaseTransportNaturalVariableLinearization
        result;
    result.state_identity =
        two_phase_detail::make_state_identity(
            state);
    result.mass_density_provenance =
        std::move(mass_density_provenance);
    result.viscosity_provenance =
        std::move(viscosity_provenance);

    const std::size_t q =
        state.layout().unknown_count();
    for (std::size_t phase = 0U;
         phase < 2U;
         ++phase) {
        two_phase_detail::validate_gradient(
            mass_density_gradient[phase],
            q,
            "mass-density");
        two_phase_detail::validate_gradient(
            viscosity_gradient[phase],
            q,
            "viscosity");
        two_phase_detail::validate_gradient(
            relative_permeability_gradient[phase],
            q,
            "relative-permeability");
        if (!std::isfinite(
                relative_permeability[phase]) ||
            relative_permeability[phase] < 0.0) {
            throw std::invalid_argument(
                "mpmc::flow: two-phase relative permeability must be finite and nonnegative");
        }
        const auto& property =
            state.phase_properties(phase);
        result.mass_density_kg_per_m3[phase] =
            property.mass_density_kg_per_m3;
        result.dynamic_viscosity_pa_s[phase] =
            property.dynamic_viscosity_pa_s;
        result.relative_permeability[phase] =
            relative_permeability[phase];
        result.mass_density_gradient[phase] =
            std::move(
                mass_density_gradient[phase]);
        result.dynamic_viscosity_gradient[phase] =
            std::move(
                viscosity_gradient[phase]);
        result.relative_permeability_gradient[phase] =
            std::move(
                relative_permeability_gradient[phase]);
    }
    return result;
}

[[nodiscard]] inline
TwoPhaseCaloricNaturalVariableLinearization
make_two_phase_caloric_linearization(
    const NaturalVariableCellState2P& state,
    std::array<std::vector<double>, 2>
        enthalpy_gradient,
    std::array<std::vector<double>, 2>
        internal_energy_gradient,
    TransportPropertyProvenance
        enthalpy_provenance,
    TransportPropertyProvenance
        internal_energy_provenance) {
    phase_transport_detail::validate_provenance(
        enthalpy_provenance,
        "two-phase-enthalpy");
    phase_transport_detail::validate_provenance(
        internal_energy_provenance,
        "two-phase-internal-energy");

    TwoPhaseCaloricNaturalVariableLinearization
        result;
    result.state_identity =
        two_phase_detail::make_state_identity(
            state);
    result.enthalpy_provenance =
        std::move(enthalpy_provenance);
    result.internal_energy_provenance =
        std::move(internal_energy_provenance);

    const std::size_t q =
        state.layout().unknown_count();
    for (std::size_t phase = 0U;
         phase < 2U;
         ++phase) {
        two_phase_detail::validate_gradient(
            enthalpy_gradient[phase],
            q,
            "enthalpy");
        two_phase_detail::validate_gradient(
            internal_energy_gradient[phase],
            q,
            "internal-energy");
        const auto& property =
            state.phase_properties(phase);
        result.specific_enthalpy_j_per_kg[phase] =
            property.specific_enthalpy_j_per_kg;
        result.specific_internal_energy_j_per_kg[phase] =
            property.specific_internal_energy_j_per_kg;
        result.specific_enthalpy_gradient[phase] =
            std::move(
                enthalpy_gradient[phase]);
        result.specific_internal_energy_gradient[phase] =
            std::move(
                internal_energy_gradient[phase]);
    }
    return result;
}

[[nodiscard]] inline
TwoPhaseRockThermalStorageLinearization
make_two_phase_rock_thermal_storage_linearization(
    const NaturalVariableCellState2P& state,
    double volumetric_internal_energy_j_per_rock_m3,
    std::vector<double> gradient,
    TransportPropertyProvenance provenance) {
    two_phase_detail::validate_gradient(
        gradient,
        state.layout().unknown_count(),
        "rock-internal-energy");
    phase_transport_detail::validate_provenance(
        provenance,
        "two-phase-rock-internal-energy");
    if (!std::isfinite(
            volumetric_internal_energy_j_per_rock_m3)) {
        throw std::invalid_argument(
            "mpmc::flow: two-phase rock internal energy must be finite");
    }
    return {
        two_phase_detail::make_state_identity(
            state),
        std::move(provenance),
        volumetric_internal_energy_j_per_rock_m3,
        std::move(gradient)};
}

[[nodiscard]] inline
TwoPhaseFugacityEquilibriumLinearization
make_two_phase_fugacity_equilibrium_linearization(
    const NaturalVariableCellState2P& state,
    std::vector<double> residual,
    std::vector<double> jacobian) {
    const std::size_t n =
        state.layout().component_count();
    const std::size_t q =
        state.layout().unknown_count();
    if (residual.size() != n ||
        n > std::numeric_limits<std::size_t>::max() /
                q ||
        jacobian.size() != n * q) {
        throw std::invalid_argument(
            "mpmc::flow: two-phase fugacity residual/Jacobian shape mismatch");
    }
    for (double value : residual) {
        if (!std::isfinite(value)) {
            throw std::invalid_argument(
                "mpmc::flow: two-phase fugacity residual contains non-finite value");
        }
    }
    for (double value : jacobian) {
        if (!std::isfinite(value)) {
            throw std::invalid_argument(
                "mpmc::flow: two-phase fugacity Jacobian contains non-finite value");
        }
    }
    return {
        state.layout().descriptor(),
        std::vector<std::string>{
            state.component_ids().begin(),
            state.component_ids().end()},
        std::move(residual),
        q,
        std::move(jacobian)};
}

[[nodiscard]] inline
PoreVolumeComponentAccumulationSnapshot3P
build_two_phase_component_accumulation(
    const NaturalVariableCellState2P& state,
    double porosity) {
    component_accumulation_detail::
        validate_porosity(porosity);

    const std::size_t n =
        state.layout().component_count();
    PoreVolumeComponentAccumulationSnapshot3P
        result;
    result.porosity = porosity;
    result.component_ids.assign(
        state.component_ids().begin(),
        state.component_ids().end());
    result.component_accumulation_mol_per_bulk_m3
        .assign(n, 0.0);

    for (std::size_t phase = 0U;
         phase < 2U;
         ++phase) {
        const double saturation =
            state.phase_saturation(phase);
        const double density =
            state.phase_properties(phase)
                .molar_density_mol_per_m3;
        const auto composition =
            state.phase_composition(phase);
        result.total_accumulation_mol_per_bulk_m3 +=
            porosity *
            saturation *
            density;
        for (std::size_t component = 0U;
             component < n;
             ++component) {
            result
                .component_accumulation_mol_per_bulk_m3[
                    component] +=
                porosity *
                saturation *
                density *
                composition[component];
        }
    }
    component_accumulation_detail::
        validate_snapshot_identity(result);
    return result;
}

[[nodiscard]] inline
PoreVolumeComponentAccumulationLinearization3P
build_two_phase_component_accumulation_linearization(
    const NaturalVariableCellState2P& state,
    double porosity,
    const TwoPhaseMolarDensityNaturalVariableLinearization&
        density) {
    component_accumulation_detail::
        validate_porosity(porosity);
    if (!two_phase_detail::same_layout(
            state.layout().descriptor(),
            density.layout)) {
        throw std::invalid_argument(
            "mpmc::flow: two-phase molar-density layout mismatch");
    }

    const std::size_t n =
        state.layout().component_count();
    const std::size_t q =
        state.layout().unknown_count();
    for (std::size_t phase = 0U;
         phase < 2U;
         ++phase) {
        two_phase_detail::validate_gradient(
            density.gradient[phase],
            q,
            "molar-density");
        if (!component_accumulation_detail::
                near_roundoff(
                    density
                        .molar_density_mol_per_m3[
                            phase],
                    state.phase_properties(phase)
                        .molar_density_mol_per_m3)) {
            throw std::invalid_argument(
                "mpmc::flow: two-phase molar-density primal mismatch");
        }
    }

    PoreVolumeComponentAccumulationLinearization3P
        result{
            state.layout().descriptor(),
            porosity,
            std::vector<std::string>{
                state.component_ids().begin(),
                state.component_ids().end()},
            q,
            std::vector<double>(
                n * q,
                0.0),
            std::vector<double>(
                q,
                0.0)};

    for (std::size_t column = 0U;
         column < q;
         ++column) {
        for (std::size_t phase = 0U;
             phase < 2U;
             ++phase) {
            const double saturation =
                state.phase_saturation(phase);
            const double ds =
                two_phase_detail::d_saturation(
                    phase,
                    column);
            const double c =
                density
                    .molar_density_mol_per_m3[
                        phase];
            const double dc =
                density.gradient[phase][column];
            const auto x =
                state.phase_composition(phase);

            result.total_accumulation_gradient[
                column] +=
                porosity *
                (ds * c +
                 saturation * dc);

            for (std::size_t component = 0U;
                 component < n;
                 ++component) {
                const double dx =
                    two_phase_detail::
                        d_composition(
                            state.layout(),
                            phase,
                            component,
                            column);
                result.component_jacobian[
                    component * q +
                    column] +=
                    porosity *
                    (ds * c *
                         x[component] +
                     saturation *
                         dc *
                         x[component] +
                     saturation *
                         c *
                         dx);
            }
        }
    }

    component_accumulation_time_detail::
        validate_current_linearization(
            build_two_phase_component_accumulation(
                state,
                porosity),
            result);
    return result;
}

[[nodiscard]] inline
PoreVolumeEnergyAccumulationSnapshot3P
build_two_phase_energy_accumulation_snapshot(
    const NaturalVariableCellState2P& state,
    double porosity,
    const TwoPhaseTransportNaturalVariableLinearization&
        transport,
    const TwoPhaseCaloricNaturalVariableLinearization&
        caloric,
    const TwoPhaseRockThermalStorageLinearization&
        rock) {
    energy_accumulation_detail::
        validate_porosity(porosity);
    const auto identity =
        two_phase_detail::make_state_identity(
            state);
    if (!energy_accumulation_detail::
            same_state_identity(
                identity,
                transport.state_identity) ||
        !energy_accumulation_detail::
            same_state_identity(
                identity,
                caloric.state_identity) ||
        !energy_accumulation_detail::
            same_state_identity(
                identity,
                rock.state_identity)) {
        throw std::invalid_argument(
            "mpmc::flow: two-phase energy property identity mismatch");
    }

    double fluid = 0.0;
    for (std::size_t phase = 0U;
         phase < 2U;
         ++phase) {
        fluid +=
            porosity *
            state.phase_saturation(phase) *
            transport.mass_density_kg_per_m3[
                phase] *
            caloric.specific_internal_energy_j_per_kg[
                phase];
    }
    const double rock_bulk =
        (1.0 - porosity) *
        rock.volumetric_internal_energy_j_per_rock_m3;
    const double total =
        fluid + rock_bulk;
    if (!std::isfinite(fluid) ||
        !std::isfinite(rock_bulk) ||
        !std::isfinite(total)) {
        throw std::range_error(
            "mpmc::flow: two-phase energy accumulation is non-finite");
    }
    return {
        identity,
        porosity,
        fluid,
        rock_bulk,
        total};
}

[[nodiscard]] inline
PoreVolumeEnergyAccumulationLinearization3P
build_two_phase_energy_accumulation_linearization(
    const NaturalVariableCellState2P& state,
    double porosity,
    const TwoPhaseTransportNaturalVariableLinearization&
        transport,
    const TwoPhaseCaloricNaturalVariableLinearization&
        caloric,
    const TwoPhaseRockThermalStorageLinearization&
        rock) {
    const auto primal =
        build_two_phase_energy_accumulation_snapshot(
            state,
            porosity,
            transport,
            caloric,
            rock);
    const std::size_t q =
        state.layout().unknown_count();

    std::vector<double> gradient(
        q,
        0.0);
    for (std::size_t column = 0U;
         column < q;
         ++column) {
        gradient[column] =
            (1.0 - porosity) *
            rock.volumetric_internal_energy_gradient[
                column];

        for (std::size_t phase = 0U;
             phase < 2U;
             ++phase) {
            two_phase_detail::validate_gradient(
                transport
                    .mass_density_gradient[phase],
                q,
                "mass-density");
            two_phase_detail::validate_gradient(
                caloric
                    .specific_internal_energy_gradient[
                        phase],
                q,
                "internal-energy");
            const double saturation =
                state.phase_saturation(phase);
            const double ds =
                two_phase_detail::d_saturation(
                    phase,
                    column);
            const double rho =
                transport.mass_density_kg_per_m3[
                    phase];
            const double drho =
                transport
                    .mass_density_gradient[phase][
                        column];
            const double u =
                caloric
                    .specific_internal_energy_j_per_kg[
                        phase];
            const double du =
                caloric
                    .specific_internal_energy_gradient[
                        phase][column];

            gradient[column] +=
                porosity *
                (ds * rho * u +
                 saturation * drho * u +
                 saturation * rho * du);
        }
        if (!std::isfinite(
                gradient[column])) {
            throw std::range_error(
                "mpmc::flow: two-phase energy accumulation derivative is non-finite");
        }
    }

    return {
        primal.state_identity,
        porosity,
        primal.total_internal_energy_j_per_bulk_m3,
        q,
        std::move(gradient)};
}

[[nodiscard]] inline
TwoPhaseMobilityLinearization
build_two_phase_mobility_linearization(
    const NaturalVariableCellState2P& state,
    const TwoPhaseTransportNaturalVariableLinearization&
        transport) {
    const auto identity =
        two_phase_detail::make_state_identity(
            state);
    if (!energy_accumulation_detail::
            same_state_identity(
                identity,
                transport.state_identity)) {
        throw std::invalid_argument(
            "mpmc::flow: two-phase mobility state identity mismatch");
    }

    TwoPhaseMobilityLinearization result;
    result.state_identity = identity;
    const std::size_t q =
        state.layout().unknown_count();
    for (std::size_t phase = 0U;
         phase < 2U;
         ++phase) {
        two_phase_detail::validate_gradient(
            transport.mass_density_gradient[
                phase],
            q,
            "mass-density");
        two_phase_detail::validate_gradient(
            transport.dynamic_viscosity_gradient[
                phase],
            q,
            "viscosity");
        two_phase_detail::validate_gradient(
            transport.relative_permeability_gradient[
                phase],
            q,
            "relative-permeability");

        const double mu =
            transport.dynamic_viscosity_pa_s[
                phase];
        const double kr =
            transport.relative_permeability[
                phase];
        if (!std::isfinite(mu) ||
            !(mu > 0.0) ||
            !std::isfinite(kr) ||
            kr < 0.0) {
            throw std::invalid_argument(
                "mpmc::flow: two-phase mobility input is invalid");
        }

        result.phase_pressure_pa[phase] =
            state.reference_pressure_pa();
        result.mass_density_kg_per_m3[phase] =
            transport.mass_density_kg_per_m3[
                phase];
        result.dynamic_viscosity_pa_s[phase] =
            mu;
        result.relative_permeability[phase] =
            kr;
        result.mobility_per_pa_s[phase] =
            kr / mu;
        result.phase_pressure_gradient[phase]
            .assign(q, 0.0);
        result.phase_pressure_gradient[phase][
            state.layout()
                .pressure_unknown_index()] =
            1.0;
        result.mass_density_gradient[phase] =
            transport.mass_density_gradient[
                phase];
        result.mobility_gradient[phase]
            .assign(q, 0.0);
        const double mu2 = mu * mu;
        for (std::size_t column = 0U;
             column < q;
             ++column) {
            result.mobility_gradient[phase][
                column] =
                (transport
                     .relative_permeability_gradient[
                         phase][column] *
                     mu -
                 kr *
                     transport
                         .dynamic_viscosity_gradient[
                             phase][column]) /
                mu2;
            if (!std::isfinite(
                    result.mobility_gradient[phase][
                        column])) {
                throw std::range_error(
                    "mpmc::flow: two-phase mobility derivative is non-finite");
            }
        }
    }
    return result;
}

enum class TwoPhaseUpwindCellSelection {
    owner_negative_phase_potential,
    neighbour_positive_phase_potential,
    owner_exact_zero_tie
};

struct TwoPhasePotentialUpwindLinearization3D {
    NaturalVariableStateIdentity
        owner_state_identity;
    NaturalVariableStateIdentity
        neighbour_state_identity;
    std::array<double, 2>
        face_mass_density_kg_per_m3{};
    double gravity_projection_m2_per_s2{};
    std::array<double, 2>
        phase_potential_difference_pa{};
    std::array<TwoPhaseUpwindCellSelection, 2>
        upwind_selection{
            TwoPhaseUpwindCellSelection::
                owner_exact_zero_tie,
            TwoPhaseUpwindCellSelection::
                owner_exact_zero_tie};
    std::array<double, 2>
        upwind_mobility_per_pa_s{};
    std::array<std::vector<double>, 2>
        owner_phase_potential_gradient;
    std::array<std::vector<double>, 2>
        neighbour_phase_potential_gradient;
    std::array<std::vector<double>, 2>
        owner_upwind_mobility_gradient;
    std::array<std::vector<double>, 2>
        neighbour_upwind_mobility_gradient;
};

[[nodiscard]] inline
TwoPhasePotentialUpwindLinearization3D
build_two_phase_potential_upwind_linearization(
    const TwoPhaseMobilityLinearization&
        owner,
    const TwoPhaseMobilityLinearization&
        neighbour,
    GravityVector3D gravity,
    OwnerToNeighbourDisplacement3D
        displacement) {
    if (owner.state_identity.component_ids !=
            neighbour.state_identity.component_ids ||
        owner.state_identity.layout.phase_count() !=
            2U ||
        neighbour.state_identity.layout.phase_count() !=
            2U ||
        owner.state_identity.layout.component_count() !=
            neighbour.state_identity.layout.component_count()) {
        throw std::invalid_argument(
            "mpmc::flow: two-phase owner/neighbour state identity mismatch");
    }
    phase_potential_upwind_detail::
        require_finite_geometry(
            gravity,
            displacement);

    TwoPhasePotentialUpwindLinearization3D result;
    result.owner_state_identity =
        owner.state_identity;
    result.neighbour_state_identity =
        neighbour.state_identity;
    result.gravity_projection_m2_per_s2 =
        phase_potential_upwind_detail::dot(
            gravity,
            displacement);

    const std::size_t owner_q =
        owner.state_identity.layout.unknown_count();
    const std::size_t neighbour_q =
        neighbour.state_identity.layout.unknown_count();

    for (std::size_t phase = 0U;
         phase < 2U;
         ++phase) {
        const double face_density =
            0.5 *
            (owner.mass_density_kg_per_m3[
                 phase] +
             neighbour.mass_density_kg_per_m3[
                 phase]);
        const double delta =
            (neighbour.phase_pressure_pa[
                 phase] -
             owner.phase_pressure_pa[phase]) -
            face_density *
                result.gravity_projection_m2_per_s2;
        if (!std::isfinite(face_density) ||
            !(face_density > 0.0) ||
            !std::isfinite(delta)) {
            throw std::range_error(
                "mpmc::flow: two-phase face density/potential is invalid");
        }

        result.face_mass_density_kg_per_m3[
            phase] =
            face_density;
        result.phase_potential_difference_pa[
            phase] =
            delta;
        result.owner_phase_potential_gradient[
            phase].assign(owner_q, 0.0);
        result.neighbour_phase_potential_gradient[
            phase].assign(neighbour_q, 0.0);

        for (std::size_t column = 0U;
             column < owner_q;
             ++column) {
            result.owner_phase_potential_gradient[
                phase][column] =
                -owner
                     .phase_pressure_gradient[phase][
                         column] -
                0.5 *
                    owner
                        .mass_density_gradient[phase][
                            column] *
                    result
                        .gravity_projection_m2_per_s2;
        }
        for (std::size_t column = 0U;
             column < neighbour_q;
             ++column) {
            result.neighbour_phase_potential_gradient[
                phase][column] =
                neighbour
                    .phase_pressure_gradient[phase][
                        column] -
                0.5 *
                    neighbour
                        .mass_density_gradient[phase][
                            column] *
                    result
                        .gravity_projection_m2_per_s2;
        }

        auto selection =
            TwoPhaseUpwindCellSelection::
                owner_exact_zero_tie;
        if (delta < 0.0) {
            selection =
                TwoPhaseUpwindCellSelection::
                    owner_negative_phase_potential;
        } else if (delta > 0.0) {
            selection =
                TwoPhaseUpwindCellSelection::
                    neighbour_positive_phase_potential;
        }
        result.upwind_selection[phase] =
            selection;

        result.owner_upwind_mobility_gradient[
            phase].assign(owner_q, 0.0);
        result.neighbour_upwind_mobility_gradient[
            phase].assign(neighbour_q, 0.0);

        if (selection ==
            TwoPhaseUpwindCellSelection::
                neighbour_positive_phase_potential) {
            result.upwind_mobility_per_pa_s[
                phase] =
                neighbour.mobility_per_pa_s[
                    phase];
            result.neighbour_upwind_mobility_gradient[
                phase] =
                neighbour.mobility_gradient[
                    phase];
        } else {
            result.upwind_mobility_per_pa_s[
                phase] =
                owner.mobility_per_pa_s[
                    phase];
            result.owner_upwind_mobility_gradient[
                phase] =
                owner.mobility_gradient[
                    phase];
        }
    }
    return result;
}

} // namespace mpmc::flow

#endif // MPMC_FLOW_TWO_PHASE_NATURAL_VARIABLE_HPP
