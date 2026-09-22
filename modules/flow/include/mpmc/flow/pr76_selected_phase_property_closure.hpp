#ifndef MPMC_FLOW_PR76_SELECTED_PHASE_PROPERTY_CLOSURE_HPP
#define MPMC_FLOW_PR76_SELECTED_PHASE_PROPERTY_CLOSURE_HPP

#include <mpmc/ad/runtime_differentiate.hpp>
#include <mpmc/flow/energy_accumulation.hpp>
#include <mpmc/flow/fugacity_equilibrium_linearization.hpp>
#include <mpmc/thermodynamics/selected_phase_density.hpp>
#include <mpmc/thermodynamics/selected_phase_fugacity.hpp>

#include <algorithm>
#include <array>
#include <cmath>
#include <concepts>
#include <cstddef>
#include <functional>
#include <limits>
#include <span>
#include <stdexcept>
#include <string>
#include <string_view>
#include <type_traits>
#include <utility>
#include <vector>

namespace mpmc::flow {

inline constexpr std::string_view
    pr76_selected_phase_property_closure_convention =
        "flow/pr76-selected-phase-property/pc-none/v1";

/// PR76 supplies the selected EOS branch, fugacity coefficient and molar density.
/// It does not define a transport-viscosity model or an absolute caloric
/// reference. Those two properties therefore remain an explicit caller-owned
/// model boundary. The provider must preserve the scalar type so AD derivatives
/// propagate through viscosity and enthalpy without a finite-difference fallback.
template <typename Number>
struct SelectedPhaseTransportCaloricValues {
    Number dynamic_viscosity_pa_s;
    Number specific_enthalpy_j_per_kg;
};

struct SelectedPhasePropertyProvenance {
    TransportPropertyProvenance mass_density;
    TransportPropertyProvenance viscosity;
    TransportPropertyProvenance enthalpy;
    TransportPropertyProvenance internal_energy;
};

template <typename Number>
struct Pr76SelectedPhasePropertyValues {
    Number molar_density_mol_per_m3;
    Number mass_density_kg_per_m3;
    Number dynamic_viscosity_pa_s;
    Number specific_enthalpy_j_per_kg;
    Number specific_internal_energy_j_per_kg;
    Number mixture_molar_mass_kg_per_mol;
    std::vector<Number> ln_phi;
};

namespace pr76_selected_phase_property_detail {

template <typename Number>
[[nodiscard]] inline double primal_value(
    const Number& value) {
    if constexpr (requires { value.value(); }) {
        return static_cast<double>(value.value());
    } else {
        return static_cast<double>(value);
    }
}

template <typename Number>
inline void require_finite(
    const Number& value,
    const char* name,
    bool positive) {
    const double primal = primal_value(value);
    if (!std::isfinite(primal) ||
        (positive && !(primal > 0.0))) {
        throw std::range_error(
            std::string{"mpmc::flow: PR76 selected-phase "} +
            name +
            (positive
                 ? " must be finite and strictly positive"
                 : " must be finite"));
    }
}

inline void validate_provenance(
    const SelectedPhasePropertyProvenance& provenance) {
    phase_transport_detail::validate_provenance(
        provenance.mass_density,
        "PR76 selected-phase mass density");
    phase_transport_detail::validate_provenance(
        provenance.viscosity,
        "PR76 selected-phase viscosity");
    phase_transport_detail::validate_provenance(
        provenance.enthalpy,
        "PR76 selected-phase enthalpy");
    phase_transport_detail::validate_provenance(
        provenance.internal_energy,
        "PR76 selected-phase internal energy");
}

template <typename Number>
[[nodiscard]] inline std::vector<Number>
reconstruct_composition(
    const NaturalVariableLayoutDescriptor& layout,
    std::size_t phase,
    std::span<const Number> q) {
    const std::size_t n = layout.component_count();
    std::vector<Number> composition(n, Number{0.0});
    const auto slot = static_cast<PhaseSlot3>(phase);
    const std::size_t dependent =
        layout.dependent_composition_component(slot);
    Number sum{0.0};

    for (std::size_t component = 0U;
         component < n;
         ++component) {
        const auto column =
            layout.independent_composition_unknown_index(
                slot,
                component);
        if (!column) {
            if (component != dependent) {
                throw std::logic_error(
                    "mpmc::flow: PR76 selected-phase composition pivot identity mismatch");
            }
            continue;
        }
        composition[component] = q[*column];
        require_finite(
            composition[component],
            "independent composition",
            true);
        sum += composition[component];
    }

    composition[dependent] =
        Number{1.0} - sum;
    require_finite(
        composition[dependent],
        "dependent composition",
        true);
    return composition;
}

inline void validate_component_ids(
    std::span<const std::string> expected,
    std::span<const std::string> actual) {
    if (expected.size() != actual.size() ||
        !std::equal(
            expected.begin(),
            expected.end(),
            actual.begin())) {
        throw std::invalid_argument(
            "mpmc::flow: PR76 selected-phase component identity/order mismatch");
    }
}

[[nodiscard]] inline std::vector<double>
independent_composition_values(
    const NaturalVariableLayoutDescriptor& layout,
    std::size_t phase,
    std::span<const double> q) {
    std::vector<double> result;
    result.reserve(
        layout.component_count() - 1U);
    const auto slot =
        static_cast<PhaseSlot3>(phase);
    for (std::size_t component = 0U;
         component < layout.component_count();
         ++component) {
        const auto column =
            layout.independent_composition_unknown_index(
                slot,
                component);
        if (column) {
            result.push_back(q[*column]);
        }
    }
    return result;
}

} // namespace pr76_selected_phase_property_detail

/// Fixed-selected-branch PR76 property closure.
///
/// The auxiliary provider is called as:
///
/// provider(selection, p, T, x, mixture_molar_mass, molar_density, mass_density)
///
/// and must return SelectedPhaseTransportCaloricValues<Number>. The selected
/// branch is never changed by this class. No viscosity/enthalpy correlation is
/// guessed. Component molar masses must be present in the ordered component
/// snapshot; otherwise construction fails.
///
/// Specific internal energy is derived from the total specific enthalpy using
/// the thermodynamic identity u = h - p/rho_mass. This v1 closure is explicitly
/// pc=none: every active phase is evaluated at the reference pressure.
template <std::floating_point T, typename TransportCaloricProvider>
class Pr76SelectedPhasePropertyClosure {
public:
    Pr76SelectedPhasePropertyClosure(
        const thermodynamics::Pr76Phase<T>& model,
        std::vector<thermodynamics::Pr76SelectedPhase>
            selections,
        TransportCaloricProvider provider,
        SelectedPhasePropertyProvenance provenance)
        : model_(&model),
          selections_(std::move(selections)),
          provider_(std::move(provider)),
          provenance_(std::move(provenance)) {
        if (selections_.empty() ||
            selections_.size() >
                fixed_three_phase_count) {
            throw std::invalid_argument(
                "mpmc::flow: PR76 selected-phase closure requires 1..3 frozen selections");
        }
        pr76_selected_phase_property_detail::
            validate_provenance(provenance_);

        const auto& components =
            model.parameters().components();
        component_ids_.reserve(components.size());
        molar_mass_kg_per_mol_.reserve(
            components.size());
        for (const auto& component :
             components.items()) {
            component_ids_.push_back(component.id);
            if (!component.molar_mass.has_value() ||
                component.molar_mass->unit !=
                    thermodynamics::Unit::
                        kilogram_per_mole ||
                !std::isfinite(
                    component.molar_mass->value) ||
                !(component.molar_mass->value >
                  0.0)) {
                throw std::invalid_argument(
                    "mpmc::flow: PR76 selected-phase mass density requires an explicit positive kg/mol molar mass for every ordered component");
            }
            molar_mass_kg_per_mol_.push_back(
                component.molar_mass->value);
        }
        if (component_ids_.size() < 2U) {
            throw std::invalid_argument(
                "mpmc::flow: PR76 selected-phase closure requires at least two components");
        }
    }

    [[nodiscard]] std::size_t
    phase_count() const noexcept {
        return selections_.size();
    }

    [[nodiscard]] std::span<const std::string>
    component_ids() const noexcept {
        return component_ids_;
    }

    [[nodiscard]] const SelectedPhasePropertyProvenance&
    provenance() const noexcept {
        return provenance_;
    }

    template <typename Number>
    [[nodiscard]]
    Pr76SelectedPhasePropertyValues<Number>
    evaluate(
        std::size_t phase,
        const Number& pressure_pa,
        const Number& temperature_k,
        std::span<const Number> composition) const {
        using namespace
            pr76_selected_phase_property_detail;

        if (phase >= selections_.size()) {
            throw std::out_of_range(
                "mpmc::flow: PR76 selected-phase selection index out of range");
        }
        if (composition.size() !=
            component_ids_.size()) {
            throw std::invalid_argument(
                "mpmc::flow: PR76 selected-phase composition size mismatch");
        }
        require_finite(
            pressure_pa,
            "pressure [Pa]",
            true);
        require_finite(
            temperature_k,
            "temperature [K]",
            true);

        Number mixture_molar_mass{0.0};
        for (std::size_t component = 0U;
             component < composition.size();
             ++component) {
            require_finite(
                composition[component],
                "composition",
                true);
            mixture_molar_mass +=
                composition[component] *
                molar_mass_kg_per_mol_[component];
        }
        require_finite(
            mixture_molar_mass,
            "mixture molar mass [kg/mol]",
            true);

        thermodynamics::Pr76PhaseWorkspace<Number>
            density_workspace;
        const auto density =
            thermodynamics::
                evaluate_selected_phase_molar_density(
                    *model_,
                    pressure_pa,
                    temperature_k,
                    composition,
                    selections_[phase],
                    density_workspace);

        thermodynamics::Pr76PhaseWorkspace<Number>
            fugacity_workspace;
        auto fugacity =
            thermodynamics::
                evaluate_selected_phase_fugacity(
                    *model_,
                    pressure_pa,
                    temperature_k,
                    composition,
                    selections_[phase],
                    fugacity_workspace);

        const Number mass_density =
            density.molar_density_mol_per_m3 *
            mixture_molar_mass;
        require_finite(
            density.molar_density_mol_per_m3,
            "molar density [mol/m3]",
            true);
        require_finite(
            mass_density,
            "mass density [kg/m3]",
            true);

        auto auxiliary =
            std::invoke(
                provider_,
                selections_[phase],
                pressure_pa,
                temperature_k,
                composition,
                mixture_molar_mass,
                density.molar_density_mol_per_m3,
                mass_density);
        require_finite(
            auxiliary.dynamic_viscosity_pa_s,
            "dynamic viscosity [Pa s]",
            true);
        require_finite(
            auxiliary.specific_enthalpy_j_per_kg,
            "specific enthalpy [J/kg]",
            false);

        const Number internal_energy =
            auxiliary.specific_enthalpy_j_per_kg -
            pressure_pa / mass_density;
        require_finite(
            internal_energy,
            "specific internal energy [J/kg]",
            false);

        if (fugacity.ln_phi.size() !=
            component_ids_.size()) {
            throw std::invalid_argument(
                "mpmc::flow: PR76 selected-phase fugacity component count mismatch");
        }
        for (const auto& value : fugacity.ln_phi) {
            require_finite(
                value,
                "ln(phi)",
                false);
        }

        return {
            density.molar_density_mol_per_m3,
            mass_density,
            auxiliary.dynamic_viscosity_pa_s,
            auxiliary.specific_enthalpy_j_per_kg,
            internal_energy,
            mixture_molar_mass,
            std::move(fugacity.ln_phi)};
    }

private:
    const thermodynamics::Pr76Phase<T>* model_;
    std::vector<thermodynamics::Pr76SelectedPhase>
        selections_;
    TransportCaloricProvider provider_;
    SelectedPhasePropertyProvenance provenance_;
    std::vector<std::string> component_ids_;
    std::vector<double> molar_mass_kg_per_mol_;
};

template <std::floating_point T, typename Provider>
[[nodiscard]] inline auto
make_pr76_selected_phase_property_closure(
    const thermodynamics::Pr76Phase<T>& model,
    std::vector<thermodynamics::Pr76SelectedPhase>
        selections,
    Provider provider,
    SelectedPhasePropertyProvenance provenance) {
    return Pr76SelectedPhasePropertyClosure<
        T,
        Provider>{
        model,
        std::move(selections),
        std::move(provider),
        std::move(provenance)};
}

/// Topology-neutral P=1/2/3 property/Jacobian publication on one frozen
/// natural-variable chart. Gradients are with respect to q=P*Nc+1.
///
/// Equilibrium rows use phase0 as reference and are ordered by phase then
/// canonical component. P=1 therefore publishes zero equilibrium rows.
struct Pr76SelectedPhasePropertyChartLinearization {
    NaturalVariableLayoutDescriptor layout{
        std::size_t{2U},
        std::size_t{1U},
        std::vector<std::size_t>{1U}};
    std::vector<std::string> component_ids;
    SelectedPhasePropertyProvenance provenance;

    std::vector<PhasePropertyPayload>
        phase_properties;
    std::vector<std::vector<double>>
        phase_compositions;

    std::vector<std::vector<double>>
        molar_density_gradient;
    std::vector<std::vector<double>>
        mass_density_gradient;
    std::vector<std::vector<double>>
        viscosity_gradient;
    std::vector<std::vector<double>>
        enthalpy_gradient;
    std::vector<std::vector<double>>
        internal_energy_gradient;

    std::vector<std::vector<double>> ln_phi;
    /// ln_phi_jacobian[phase][component*q + column].
    std::vector<std::vector<double>>
        ln_phi_jacobian;

    std::vector<double> equilibrium_residual;
    /// Row-major [(P-1)*Nc] x q.
    std::vector<double> equilibrium_jacobian;
};

template <
    std::size_t K = 4U,
    std::floating_point T,
    typename Provider>
[[nodiscard]] inline
Pr76SelectedPhasePropertyChartLinearization
evaluate_pr76_selected_phase_property_chart(
    const NaturalVariableLayoutDescriptor& layout,
    std::span<const std::string> component_ids,
    std::span<const double> natural_variables,
    const Pr76SelectedPhasePropertyClosure<
        T,
        Provider>& closure,
    mpmc::ad::RuntimeJacobianWorkspace<
        double,
        K>& workspace) {
    using namespace
        pr76_selected_phase_property_detail;

    const std::size_t p_count =
        layout.phase_count();
    const std::size_t n =
        layout.component_count();
    const std::size_t q =
        layout.unknown_count();
    if (p_count != closure.phase_count() ||
        component_ids.size() != n ||
        natural_variables.size() != q) {
        throw std::invalid_argument(
            "mpmc::flow: PR76 selected-phase chart shape/cardinality mismatch");
    }
    validate_component_ids(
        closure.component_ids(),
        component_ids);

    const std::size_t property_stride =
        5U + n;
    if (p_count >
            std::numeric_limits<std::size_t>::
                    max() /
                property_stride) {
        throw std::length_error(
            "mpmc::flow: PR76 selected-phase output shape overflow");
    }
    const std::size_t property_outputs =
        p_count * property_stride;
    const std::size_t equilibrium_rows =
        (p_count - 1U) * n;
    if (property_outputs >
        std::numeric_limits<std::size_t>::
                max() -
            equilibrium_rows) {
        throw std::length_error(
            "mpmc::flow: PR76 selected-phase output count overflow");
    }
    const std::size_t output_count =
        property_outputs +
        equilibrium_rows;

    using D = mpmc::ad::Dual<double, K>;
    const auto callback =
        [&](std::span<const D> active,
            std::span<D> outputs) {
            if (active.size() != q ||
                outputs.size() !=
                    output_count) {
                throw std::logic_error(
                    "mpmc::flow: PR76 selected-phase AD callback shape mismatch");
            }
            const D& pressure =
                active[
                    layout
                        .pressure_unknown_index()];
            const D& temperature =
                active[
                    layout
                        .temperature_unknown_index()];

            std::vector<std::vector<D>>
                compositions;
            std::vector<
                Pr76SelectedPhasePropertyValues<D>>
                properties;
            compositions.reserve(p_count);
            properties.reserve(p_count);

            for (std::size_t phase = 0U;
                 phase < p_count;
                 ++phase) {
                compositions.push_back(
                    reconstruct_composition(
                        layout,
                        phase,
                        active));
                properties.push_back(
                    closure.evaluate(
                        phase,
                        pressure,
                        temperature,
                        std::span<const D>{
                            compositions.back()}));

                const std::size_t base =
                    phase * property_stride;
                outputs[base + 0U] =
                    properties.back()
                        .molar_density_mol_per_m3;
                outputs[base + 1U] =
                    properties.back()
                        .mass_density_kg_per_m3;
                outputs[base + 2U] =
                    properties.back()
                        .dynamic_viscosity_pa_s;
                outputs[base + 3U] =
                    properties.back()
                        .specific_enthalpy_j_per_kg;
                outputs[base + 4U] =
                    properties.back()
                        .specific_internal_energy_j_per_kg;
                for (std::size_t component = 0U;
                     component < n;
                     ++component) {
                    outputs[
                        base + 5U +
                        component] =
                        properties.back()
                            .ln_phi[component];
                }
            }

            using std::log;
            for (std::size_t phase = 1U;
                 phase < p_count;
                 ++phase) {
                for (std::size_t component = 0U;
                     component < n;
                     ++component) {
                    const std::size_t row =
                        (phase - 1U) * n +
                        component;
                    outputs[
                        property_outputs +
                        row] =
                        log(
                            compositions[0][component] /
                            compositions[phase][component]) +
                        properties[0]
                            .ln_phi[component] -
                        properties[phase]
                            .ln_phi[component] +
                        log(pressure / pressure);
                }
            }
        };

    const auto derived =
        workspace.evaluate(
            callback,
            natural_variables,
            output_count,
            {
                q,
                output_count,
                output_count >
                        std::numeric_limits<
                            std::size_t>::max() /
                            q
                    ? 0U
                    : output_count * q});
    if (derived.output_count !=
            output_count ||
        derived.input_count != q) {
        throw std::logic_error(
            "mpmc::flow: PR76 selected-phase AD driver returned unexpected shape");
    }

    Pr76SelectedPhasePropertyChartLinearization
        result{
            layout,
            std::vector<std::string>{
                component_ids.begin(),
                component_ids.end()},
            closure.provenance(),
            {},
            {},
            {},
            {},
            {},
            {},
            {},
            {},
            {},
            {},
            {}};

    result.phase_properties.reserve(p_count);
    result.phase_compositions.reserve(p_count);
    result.molar_density_gradient.reserve(
        p_count);
    result.mass_density_gradient.reserve(
        p_count);
    result.viscosity_gradient.reserve(p_count);
    result.enthalpy_gradient.reserve(p_count);
    result.internal_energy_gradient.reserve(
        p_count);
    result.ln_phi.reserve(p_count);
    result.ln_phi_jacobian.reserve(p_count);

    for (std::size_t phase = 0U;
         phase < p_count;
         ++phase) {
        const std::size_t base =
            phase * property_stride;
        result.phase_properties.push_back(
            PhasePropertyPayload{
                derived.values[base + 0U],
                derived.values[base + 1U],
                derived.values[base + 2U],
                derived.values[base + 3U],
                derived.values[base + 4U]});

        result.phase_compositions.push_back(
            reconstruct_composition(
                layout,
                phase,
                natural_variables));

        const auto gradient_for =
            [&](std::size_t local_output) {
                const std::size_t row =
                    base + local_output;
                return std::vector<double>{
                    derived.jacobian.begin() +
                        static_cast<
                            std::ptrdiff_t>(
                            row * q),
                    derived.jacobian.begin() +
                        static_cast<
                            std::ptrdiff_t>(
                            (row + 1U) * q)};
            };
        result.molar_density_gradient.push_back(
            gradient_for(0U));
        result.mass_density_gradient.push_back(
            gradient_for(1U));
        result.viscosity_gradient.push_back(
            gradient_for(2U));
        result.enthalpy_gradient.push_back(
            gradient_for(3U));
        result.internal_energy_gradient.push_back(
            gradient_for(4U));

        std::vector<double> ln_phi(n);
        std::vector<double> ln_phi_jacobian(
            n * q);
        for (std::size_t component = 0U;
             component < n;
             ++component) {
            const std::size_t row =
                base + 5U +
                component;
            ln_phi[component] =
                derived.values[row];
            std::copy_n(
                derived.jacobian.begin() +
                    static_cast<std::ptrdiff_t>(
                        row * q),
                q,
                ln_phi_jacobian.begin() +
                    static_cast<
                        std::ptrdiff_t>(
                        component * q));
        }
        result.ln_phi.push_back(
            std::move(ln_phi));
        result.ln_phi_jacobian.push_back(
            std::move(ln_phi_jacobian));
    }

    result.equilibrium_residual.assign(
        derived.values.begin() +
            static_cast<std::ptrdiff_t>(
                property_outputs),
        derived.values.end());
    result.equilibrium_jacobian.assign(
        derived.jacobian.begin() +
            static_cast<std::ptrdiff_t>(
                property_outputs * q),
        derived.jacobian.end());
    return result;
}

/// Direct bridge from the PR76 P=3 property publication into the existing
/// fixed-three-phase flow carriers. Saturation constitutive laws and rock
/// thermal storage remain separate caller-owned models.
struct Pr76SelectedPhaseFlowLinearization3P {
    NaturalVariableCellState3P state;
    PhaseMolarDensityNaturalVariableLinearization3P
        molar_density;
    PhaseTransportPropertyNaturalVariableLinearization3P
        transport;
    PhaseCaloricPropertyNaturalVariableLinearization3P
        caloric;
    FugacityEquilibriumResidualLinearization3P
        fugacity;
};

[[nodiscard]] inline
Pr76SelectedPhaseFlowLinearization3P
make_pr76_selected_phase_flow_linearization_3p(
    const Pr76SelectedPhasePropertyChartLinearization&
        properties,
    std::span<const double> natural_variables) {
    using namespace
        pr76_selected_phase_property_detail;

    if (properties.layout.phase_count() !=
            fixed_three_phase_count ||
        properties.phase_properties.size() !=
            fixed_three_phase_count ||
        properties.phase_compositions.size() !=
            fixed_three_phase_count ||
        natural_variables.size() !=
            properties.layout.unknown_count()) {
        throw std::invalid_argument(
            "mpmc::flow: PR76 three-phase property bridge shape mismatch");
    }

    const auto layout3 =
        static_cast<NaturalVariableLayout3P>(
            properties.layout);
    const std::size_t q =
        layout3.unknown_count();
    const std::size_t n =
        layout3.component_count();

    const auto s0 =
        layout3.independent_saturation_unknown_index(
            PhaseSlot3::phase0);
    const auto s1 =
        layout3.independent_saturation_unknown_index(
            PhaseSlot3::phase1);
    if (!s0 || !s1) {
        throw std::logic_error(
            "mpmc::flow: PR76 three-phase saturation chart is malformed");
    }

    NaturalVariableCellStateInput3P state_input;
    state_input.component_ids =
        properties.component_ids;
    state_input.reference_pressure_pa =
        natural_variables[
            layout3.pressure_unknown_index()];
    state_input.temperature_k =
        natural_variables[
            layout3.temperature_unknown_index()];
    state_input.independent_saturations = {
        natural_variables[*s0],
        natural_variables[*s1]};
    state_input.composition_pivot =
        NaturalVariableCompositionPivot3P::
            from_dependent_components(
                n,
                layout3
                    .composition_pivot()
                    .dependent_components());

    for (std::size_t phase = 0U;
         phase < fixed_three_phase_count;
         ++phase) {
        state_input
            .independent_phase_compositions[phase] =
            independent_composition_values(
                properties.layout,
                phase,
                natural_variables);
        const auto& property =
            properties.phase_properties[phase];
        state_input.phase_properties[phase] = {
            property.molar_density_mol_per_m3,
            property.mass_density_kg_per_m3,
            property.dynamic_viscosity_pa_s,
            property.specific_enthalpy_j_per_kg,
            property.specific_internal_energy_j_per_kg};
    }

    auto state =
        NaturalVariableCellState3P::create(
            std::move(state_input));

    std::array<double, 3>
        molar_density{};
    std::array<std::vector<double>, 3>
        molar_density_gradient;
    std::array<std::vector<double>, 3>
        mass_density_gradient;
    std::array<std::vector<double>, 3>
        viscosity_gradient;
    std::array<std::vector<double>, 3>
        enthalpy_gradient;
    std::array<std::vector<double>, 3>
        internal_energy_gradient;

    for (std::size_t phase = 0U;
         phase < fixed_three_phase_count;
         ++phase) {
        molar_density[phase] =
            properties
                .phase_properties[phase]
                .molar_density_mol_per_m3;
        molar_density_gradient[phase] =
            properties
                .molar_density_gradient
                .at(phase);
        mass_density_gradient[phase] =
            properties
                .mass_density_gradient
                .at(phase);
        viscosity_gradient[phase] =
            properties
                .viscosity_gradient
                .at(phase);
        enthalpy_gradient[phase] =
            properties
                .enthalpy_gradient
                .at(phase);
        internal_energy_gradient[phase] =
            properties
                .internal_energy_gradient
                .at(phase);
    }

    PhaseMolarDensityNaturalVariableLinearization3P
        density{
            layout3,
            molar_density,
            std::move(
                molar_density_gradient)};

    auto transport =
        make_phase_transport_property_linearization(
            state,
            std::move(
                mass_density_gradient),
            std::move(
                viscosity_gradient),
            properties.provenance.mass_density,
            properties.provenance.viscosity);

    auto caloric =
        make_phase_caloric_property_linearization(
            state,
            std::move(
                enthalpy_gradient),
            std::move(
                internal_energy_gradient),
            properties.provenance.enthalpy,
            properties.provenance.internal_energy);

    auto residual =
        FugacityEquilibriumResidual3P<double>{
            n,
            properties.equilibrium_residual};
    auto fugacity =
        FugacityEquilibriumResidualLinearization3P{
            layout3,
            properties.component_ids,
            std::move(residual),
            q,
            properties.equilibrium_jacobian};

    return {
        std::move(state),
        std::move(density),
        std::move(transport),
        std::move(caloric),
        std::move(fugacity)};
}

} // namespace mpmc::flow

#endif // MPMC_FLOW_PR76_SELECTED_PHASE_PROPERTY_CLOSURE_HPP
