#ifndef MPMC_FLOW_SINGLE_PHASE_NATURAL_VARIABLE_HPP
#define MPMC_FLOW_SINGLE_PHASE_NATURAL_VARIABLE_HPP

#include <mpmc/flow/component_accumulation_time.hpp>
#include <mpmc/flow/energy_accumulation.hpp>
#include <mpmc/flow/phase_potential_upwind.hpp>

#include <algorithm>
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

namespace single_phase_detail {

[[nodiscard]] inline NaturalVariableStateIdentity
make_state_identity(
    const NaturalVariableCellState1P& state) {
    return {
        state.layout().descriptor(),
        std::vector<std::string>{
            state.component_ids().begin(),
            state.component_ids().end()},
        state.reference_pressure_pa(),
        state.temperature_k(),
        {1.0, 0.0, 0.0},
        {
            std::vector<double>{
                state.phase_composition().begin(),
                state.phase_composition().end()},
            std::vector<double>{},
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
            std::string{"mpmc::flow: single-phase "} +
            name +
            " gradient shape mismatch");
    }
    for (double value : gradient) {
        if (!std::isfinite(value)) {
            throw std::invalid_argument(
                std::string{"mpmc::flow: single-phase "} +
                name +
                " gradient contains non-finite derivative");
        }
    }
}

[[nodiscard]] inline double
d_composition(
    const NaturalVariableLayout1P& layout,
    std::size_t component,
    std::size_t column) {
    if (column < 2U) {
        return 0.0;
    }
    for (std::size_t candidate = 0U;
         candidate < layout.component_count();
         ++candidate) {
        const auto candidate_column =
            layout.independent_composition_unknown_index(
                candidate);
        if (!candidate_column ||
            *candidate_column != column) {
            continue;
        }
        if (candidate == component) {
            return 1.0;
        }
        if (component ==
            layout.dependent_composition_component()) {
            return -1.0;
        }
        return 0.0;
    }
    return 0.0;
}

} // namespace single_phase_detail

struct SinglePhaseMolarDensityNaturalVariableLinearization {
    NaturalVariableLayoutDescriptor layout{
        std::size_t{2U},
        std::size_t{1U},
        std::vector<std::size_t>{1U}};
    double molar_density_mol_per_m3{};
    std::vector<double> gradient;
};

struct SinglePhaseTransportNaturalVariableLinearization {
    NaturalVariableStateIdentity state_identity;
    TransportPropertyProvenance mass_density_provenance;
    TransportPropertyProvenance viscosity_provenance;
    double mass_density_kg_per_m3{};
    double dynamic_viscosity_pa_s{};
    double relative_permeability{};
    std::vector<double> mass_density_gradient;
    std::vector<double> dynamic_viscosity_gradient;
    std::vector<double> relative_permeability_gradient;
};

struct SinglePhaseCaloricNaturalVariableLinearization {
    NaturalVariableStateIdentity state_identity;
    TransportPropertyProvenance enthalpy_provenance;
    TransportPropertyProvenance internal_energy_provenance;
    double specific_enthalpy_j_per_kg{};
    double specific_internal_energy_j_per_kg{};
    std::vector<double> specific_enthalpy_gradient;
    std::vector<double> specific_internal_energy_gradient;
};

struct SinglePhaseRockThermalStorageLinearization {
    NaturalVariableStateIdentity state_identity;
    TransportPropertyProvenance provenance;
    double volumetric_internal_energy_j_per_rock_m3{};
    std::vector<double>
        volumetric_internal_energy_gradient;
};

struct SinglePhaseMobilityLinearization {
    NaturalVariableStateIdentity state_identity;
    double phase_pressure_pa{};
    double mass_density_kg_per_m3{};
    double dynamic_viscosity_pa_s{};
    double relative_permeability{};
    double mobility_per_pa_s{};
    std::vector<double> phase_pressure_gradient;
    std::vector<double> mass_density_gradient;
    std::vector<double> mobility_gradient;
};

[[nodiscard]] inline
SinglePhaseMolarDensityNaturalVariableLinearization
make_single_phase_molar_density_linearization(
    const NaturalVariableCellState1P& state,
    std::vector<double> gradient) {
    single_phase_detail::validate_gradient(
        gradient,
        state.layout().unknown_count(),
        "molar-density");
    return {
        state.layout().descriptor(),
        state.phase_properties()
            .molar_density_mol_per_m3,
        std::move(gradient)};
}

[[nodiscard]] inline
SinglePhaseTransportNaturalVariableLinearization
make_single_phase_transport_linearization(
    const NaturalVariableCellState1P& state,
    std::vector<double> mass_density_gradient,
    std::vector<double> viscosity_gradient,
    double relative_permeability,
    std::vector<double>
        relative_permeability_gradient,
    TransportPropertyProvenance
        mass_density_provenance,
    TransportPropertyProvenance
        viscosity_provenance) {
    const std::size_t q =
        state.layout().unknown_count();
    single_phase_detail::validate_gradient(
        mass_density_gradient,
        q,
        "mass-density");
    single_phase_detail::validate_gradient(
        viscosity_gradient,
        q,
        "viscosity");
    single_phase_detail::validate_gradient(
        relative_permeability_gradient,
        q,
        "relative-permeability");
    phase_transport_detail::validate_provenance(
        mass_density_provenance,
        "single-phase-mass-density");
    phase_transport_detail::validate_provenance(
        viscosity_provenance,
        "single-phase-viscosity");
    if (!std::isfinite(relative_permeability) ||
        relative_permeability < 0.0) {
        throw std::invalid_argument(
            "mpmc::flow: single-phase relative permeability must be finite and nonnegative");
    }

    const auto& property =
        state.phase_properties();
    return {
        single_phase_detail::make_state_identity(
            state),
        std::move(mass_density_provenance),
        std::move(viscosity_provenance),
        property.mass_density_kg_per_m3,
        property.dynamic_viscosity_pa_s,
        relative_permeability,
        std::move(mass_density_gradient),
        std::move(viscosity_gradient),
        std::move(
            relative_permeability_gradient)};
}

[[nodiscard]] inline
SinglePhaseCaloricNaturalVariableLinearization
make_single_phase_caloric_linearization(
    const NaturalVariableCellState1P& state,
    std::vector<double> enthalpy_gradient,
    std::vector<double> internal_energy_gradient,
    TransportPropertyProvenance
        enthalpy_provenance,
    TransportPropertyProvenance
        internal_energy_provenance) {
    const std::size_t q =
        state.layout().unknown_count();
    single_phase_detail::validate_gradient(
        enthalpy_gradient,
        q,
        "enthalpy");
    single_phase_detail::validate_gradient(
        internal_energy_gradient,
        q,
        "internal-energy");
    phase_transport_detail::validate_provenance(
        enthalpy_provenance,
        "single-phase-enthalpy");
    phase_transport_detail::validate_provenance(
        internal_energy_provenance,
        "single-phase-internal-energy");
    const auto& property =
        state.phase_properties();
    return {
        single_phase_detail::make_state_identity(
            state),
        std::move(enthalpy_provenance),
        std::move(internal_energy_provenance),
        property.specific_enthalpy_j_per_kg,
        property.specific_internal_energy_j_per_kg,
        std::move(enthalpy_gradient),
        std::move(
            internal_energy_gradient)};
}

[[nodiscard]] inline
SinglePhaseRockThermalStorageLinearization
make_single_phase_rock_thermal_storage_linearization(
    const NaturalVariableCellState1P& state,
    double volumetric_internal_energy_j_per_rock_m3,
    std::vector<double> gradient,
    TransportPropertyProvenance provenance) {
    single_phase_detail::validate_gradient(
        gradient,
        state.layout().unknown_count(),
        "rock-internal-energy");
    phase_transport_detail::validate_provenance(
        provenance,
        "single-phase-rock-internal-energy");
    if (!std::isfinite(
            volumetric_internal_energy_j_per_rock_m3)) {
        throw std::invalid_argument(
            "mpmc::flow: single-phase rock internal energy must be finite");
    }
    return {
        single_phase_detail::make_state_identity(
            state),
        std::move(provenance),
        volumetric_internal_energy_j_per_rock_m3,
        std::move(gradient)};
}

[[nodiscard]] inline
PoreVolumeComponentAccumulationSnapshot3P
build_single_phase_component_accumulation(
    const NaturalVariableCellState1P& state,
    double porosity) {
    component_accumulation_detail::
        validate_porosity(porosity);

    PoreVolumeComponentAccumulationSnapshot3P
        result;
    result.porosity = porosity;
    result.component_ids.assign(
        state.component_ids().begin(),
        state.component_ids().end());
    result.component_accumulation_mol_per_bulk_m3
        .assign(
            result.component_ids.size(),
            0.0);

    const double density =
        state.phase_properties()
            .molar_density_mol_per_m3;
    const auto composition =
        state.phase_composition();
    for (std::size_t component = 0U;
         component < composition.size();
         ++component) {
        result
            .component_accumulation_mol_per_bulk_m3[
                component] =
            porosity *
            density *
            composition[component];
    }
    result.total_accumulation_mol_per_bulk_m3 =
        porosity * density;
    component_accumulation_detail::
        validate_snapshot_identity(result);
    return result;
}

[[nodiscard]] inline
PoreVolumeComponentAccumulationLinearization3P
build_single_phase_component_accumulation_linearization(
    const NaturalVariableCellState1P& state,
    double porosity,
    const SinglePhaseMolarDensityNaturalVariableLinearization&
        density) {
    component_accumulation_detail::
        validate_porosity(porosity);
    const auto descriptor =
        state.layout().descriptor();
    if (!single_phase_detail::same_layout(
            descriptor,
            density.layout) ||
        !component_accumulation_detail::
            near_roundoff(
                density.molar_density_mol_per_m3,
                state.phase_properties()
                    .molar_density_mol_per_m3)) {
        throw std::invalid_argument(
            "mpmc::flow: single-phase molar-density linearization identity mismatch");
    }
    const std::size_t n =
        state.layout().component_count();
    const std::size_t q =
        state.layout().unknown_count();
    single_phase_detail::validate_gradient(
        density.gradient,
        q,
        "molar-density");

    PoreVolumeComponentAccumulationLinearization3P
        result{
            descriptor,
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

    const double c =
        density.molar_density_mol_per_m3;
    const auto x =
        state.phase_composition();
    for (std::size_t column = 0U;
         column < q;
         ++column) {
        result.total_accumulation_gradient[
            column] =
            porosity *
            density.gradient[column];
        for (std::size_t component = 0U;
             component < n;
             ++component) {
            const double dx =
                single_phase_detail::
                    d_composition(
                        state.layout(),
                        component,
                        column);
            result.component_jacobian[
                component * q + column] =
                porosity *
                (density.gradient[column] *
                     x[component] +
                 c * dx);
        }
    }
    component_accumulation_time_detail::
        validate_current_linearization(
            build_single_phase_component_accumulation(
                state,
                porosity),
            result);
    return result;
}

[[nodiscard]] inline
PoreVolumeEnergyAccumulationSnapshot3P
build_single_phase_energy_accumulation_snapshot(
    const NaturalVariableCellState1P& state,
    double porosity,
    const SinglePhaseTransportNaturalVariableLinearization&
        transport,
    const SinglePhaseCaloricNaturalVariableLinearization&
        caloric,
    const SinglePhaseRockThermalStorageLinearization&
        rock) {
    energy_accumulation_detail::
        validate_porosity(porosity);
    const auto identity =
        single_phase_detail::
            make_state_identity(state);
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
            "mpmc::flow: single-phase energy property identity mismatch");
    }
    const double fluid =
        porosity *
        transport.mass_density_kg_per_m3 *
        caloric.specific_internal_energy_j_per_kg;
    const double rock_bulk =
        (1.0 - porosity) *
        rock.volumetric_internal_energy_j_per_rock_m3;
    const double total =
        fluid + rock_bulk;
    if (!std::isfinite(fluid) ||
        !std::isfinite(rock_bulk) ||
        !std::isfinite(total)) {
        throw std::range_error(
            "mpmc::flow: single-phase energy accumulation is non-finite");
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
build_single_phase_energy_accumulation_linearization(
    const NaturalVariableCellState1P& state,
    double porosity,
    const SinglePhaseTransportNaturalVariableLinearization&
        transport,
    const SinglePhaseCaloricNaturalVariableLinearization&
        caloric,
    const SinglePhaseRockThermalStorageLinearization&
        rock) {
    const auto primal =
        build_single_phase_energy_accumulation_snapshot(
            state,
            porosity,
            transport,
            caloric,
            rock);
    const std::size_t q =
        state.layout().unknown_count();
    single_phase_detail::validate_gradient(
        transport.mass_density_gradient,
        q,
        "mass-density");
    single_phase_detail::validate_gradient(
        caloric.specific_internal_energy_gradient,
        q,
        "internal-energy");
    single_phase_detail::validate_gradient(
        rock.volumetric_internal_energy_gradient,
        q,
        "rock-internal-energy");

    std::vector<double> gradient(
        q,
        0.0);
    const double rho =
        transport.mass_density_kg_per_m3;
    const double u =
        caloric.specific_internal_energy_j_per_kg;
    for (std::size_t column = 0U;
         column < q;
         ++column) {
        gradient[column] =
            porosity *
                (transport
                     .mass_density_gradient[column] *
                     u +
                 rho *
                     caloric
                         .specific_internal_energy_gradient[
                             column]) +
            (1.0 - porosity) *
                rock.volumetric_internal_energy_gradient[
                    column];
        if (!std::isfinite(
                gradient[column])) {
            throw std::range_error(
                "mpmc::flow: single-phase energy accumulation derivative is non-finite");
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
SinglePhaseMobilityLinearization
build_single_phase_mobility_linearization(
    const NaturalVariableCellState1P& state,
    const SinglePhaseTransportNaturalVariableLinearization&
        transport) {
    const auto identity =
        single_phase_detail::
            make_state_identity(state);
    if (!energy_accumulation_detail::
            same_state_identity(
                identity,
                transport.state_identity)) {
        throw std::invalid_argument(
            "mpmc::flow: single-phase mobility state identity mismatch");
    }
    const std::size_t q =
        state.layout().unknown_count();
    single_phase_detail::validate_gradient(
        transport.mass_density_gradient,
        q,
        "mass-density");
    single_phase_detail::validate_gradient(
        transport.dynamic_viscosity_gradient,
        q,
        "viscosity");
    single_phase_detail::validate_gradient(
        transport.relative_permeability_gradient,
        q,
        "relative-permeability");

    const double mu =
        transport.dynamic_viscosity_pa_s;
    const double kr =
        transport.relative_permeability;
    if (!std::isfinite(mu) ||
        !(mu > 0.0) ||
        !std::isfinite(kr) ||
        kr < 0.0) {
        throw std::invalid_argument(
            "mpmc::flow: single-phase mobility inputs are invalid");
    }
    const double mobility =
        kr / mu;
    std::vector<double> mobility_gradient(
        q,
        0.0);
    const double mu2 =
        mu * mu;
    for (std::size_t column = 0U;
         column < q;
         ++column) {
        mobility_gradient[column] =
            (transport
                 .relative_permeability_gradient[column] *
                 mu -
             kr *
                 transport
                     .dynamic_viscosity_gradient[column]) /
            mu2;
        if (!std::isfinite(
                mobility_gradient[column])) {
            throw std::range_error(
                "mpmc::flow: single-phase mobility derivative is non-finite");
        }
    }

    std::vector<double>
        phase_pressure_gradient(
            q,
            0.0);
    phase_pressure_gradient[
        state.layout().pressure_unknown_index()] =
        1.0;

    return {
        identity,
        state.reference_pressure_pa(),
        transport.mass_density_kg_per_m3,
        transport.dynamic_viscosity_pa_s,
        kr,
        mobility,
        std::move(
            phase_pressure_gradient),
        transport.mass_density_gradient,
        std::move(mobility_gradient)};
}

enum class SinglePhaseUpwindCellSelection {
    owner_negative_phase_potential,
    neighbour_positive_phase_potential,
    owner_exact_zero_tie
};

struct SinglePhasePotentialUpwindLinearization3D {
    NaturalVariableStateIdentity
        owner_state_identity;
    NaturalVariableStateIdentity
        neighbour_state_identity;
    double face_mass_density_kg_per_m3{};
    double gravity_projection_m2_per_s2{};
    double phase_potential_difference_pa{};
    SinglePhaseUpwindCellSelection
        upwind_selection{
            SinglePhaseUpwindCellSelection::
                owner_exact_zero_tie};
    double upwind_mobility_per_pa_s{};
    std::vector<double>
        owner_phase_potential_gradient;
    std::vector<double>
        neighbour_phase_potential_gradient;
    std::vector<double>
        owner_upwind_mobility_gradient;
    std::vector<double>
        neighbour_upwind_mobility_gradient;
};

[[nodiscard]] inline
SinglePhasePotentialUpwindLinearization3D
build_single_phase_potential_upwind_linearization(
    const SinglePhaseMobilityLinearization&
        owner,
    const SinglePhaseMobilityLinearization&
        neighbour,
    GravityVector3D gravity,
    OwnerToNeighbourDisplacement3D
        displacement) {
    if (owner.state_identity.component_ids !=
            neighbour.state_identity.component_ids ||
        owner.state_identity.layout.phase_count() !=
            1U ||
        neighbour.state_identity.layout.phase_count() !=
            1U ||
        owner.state_identity.layout.component_count() !=
            neighbour.state_identity.layout.component_count()) {
        throw std::invalid_argument(
            "mpmc::flow: single-phase owner/neighbour state identity mismatch");
    }
    phase_potential_upwind_detail::
        require_finite_geometry(
            gravity,
            displacement);
    const double gravity_projection =
        phase_potential_upwind_detail::dot(
            gravity,
            displacement);
    const double face_density =
        0.5 *
        (owner.mass_density_kg_per_m3 +
         neighbour.mass_density_kg_per_m3);
    const double delta =
        (neighbour.phase_pressure_pa -
         owner.phase_pressure_pa) -
        face_density *
            gravity_projection;
    if (!std::isfinite(face_density) ||
        !(face_density > 0.0) ||
        !std::isfinite(delta)) {
        throw std::range_error(
            "mpmc::flow: single-phase face density/potential is invalid");
    }

    const std::size_t owner_q =
        owner.state_identity.layout.unknown_count();
    const std::size_t neighbour_q =
        neighbour.state_identity.layout.unknown_count();
    std::vector<double> owner_potential(
        owner_q,
        0.0);
    std::vector<double> neighbour_potential(
        neighbour_q,
        0.0);
    for (std::size_t column = 0U;
         column < owner_q;
         ++column) {
        owner_potential[column] =
            -owner.phase_pressure_gradient[column] -
            0.5 *
                owner.mass_density_gradient[column] *
                gravity_projection;
    }
    for (std::size_t column = 0U;
         column < neighbour_q;
         ++column) {
        neighbour_potential[column] =
            neighbour.phase_pressure_gradient[column] -
            0.5 *
                neighbour.mass_density_gradient[column] *
                gravity_projection;
    }

    SinglePhaseUpwindCellSelection selection =
        SinglePhaseUpwindCellSelection::
            owner_exact_zero_tie;
    if (delta < 0.0) {
        selection =
            SinglePhaseUpwindCellSelection::
                owner_negative_phase_potential;
    } else if (delta > 0.0) {
        selection =
            SinglePhaseUpwindCellSelection::
                neighbour_positive_phase_potential;
    }

    std::vector<double> owner_upwind(
        owner_q,
        0.0);
    std::vector<double> neighbour_upwind(
        neighbour_q,
        0.0);
    double mobility = 0.0;
    if (selection ==
            SinglePhaseUpwindCellSelection::
                neighbour_positive_phase_potential) {
        mobility =
            neighbour.mobility_per_pa_s;
        neighbour_upwind =
            neighbour.mobility_gradient;
    } else {
        mobility =
            owner.mobility_per_pa_s;
        owner_upwind =
            owner.mobility_gradient;
    }
    if (!std::isfinite(mobility) ||
        mobility < 0.0) {
        throw std::range_error(
            "mpmc::flow: single-phase upwind mobility is invalid");
    }

    return {
        owner.state_identity,
        neighbour.state_identity,
        face_density,
        gravity_projection,
        delta,
        selection,
        mobility,
        std::move(owner_potential),
        std::move(neighbour_potential),
        std::move(owner_upwind),
        std::move(neighbour_upwind)};
}

} // namespace mpmc::flow

#endif // MPMC_FLOW_SINGLE_PHASE_NATURAL_VARIABLE_HPP
