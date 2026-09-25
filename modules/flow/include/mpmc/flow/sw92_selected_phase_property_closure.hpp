#ifndef MPMC_FLOW_SW92_SELECTED_PHASE_PROPERTY_CLOSURE_HPP
#define MPMC_FLOW_SW92_SELECTED_PHASE_PROPERTY_CLOSURE_HPP

#include <mpmc/flow/pr76_selected_phase_property_closure.hpp>
#include <mpmc/thermodynamics/selected_phase_density.hpp>
#include <mpmc/thermodynamics/selected_phase_fugacity.hpp>

#include <array>
#include <cmath>
#include <concepts>
#include <cstddef>
#include <functional>
#include <span>
#include <stdexcept>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

namespace mpmc::flow {

inline constexpr std::string_view
    sw92_selected_phase_property_closure_convention =
        "flow/sw92-selected-phase-property/fixed-family-root/pc-none/v1";

template <typename Number>
using Sw92SelectedPhasePropertyValues =
    Pr76SelectedPhasePropertyValues<Number>;

using Sw92SelectedPhasePropertyChartLinearization =
    Pr76SelectedPhasePropertyChartLinearization;
using Sw92SelectedPhaseFlowLinearization1P =
    Pr76SelectedPhaseFlowLinearization1P;
using Sw92SelectedPhaseFlowLinearization2P =
    Pr76SelectedPhaseFlowLinearization2P;
using Sw92SelectedPhaseFlowLinearization3P =
    Pr76SelectedPhaseFlowLinearization3P;

class Sw92SelectedPhasePcNoneCapabilityError final
    : public std::runtime_error {
public:
    using std::runtime_error::runtime_error;
};

namespace sw92_selected_phase_property_detail {

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
            std::string{"mpmc::flow: SW92 selected-phase "} +
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
        "SW92 selected-phase mass density");
    phase_transport_detail::validate_provenance(
        provenance.viscosity,
        "SW92 selected-phase viscosity");
    phase_transport_detail::validate_provenance(
        provenance.enthalpy,
        "SW92 selected-phase enthalpy");
    phase_transport_detail::validate_provenance(
        provenance.internal_energy,
        "SW92 selected-phase internal energy");
}

template <std::floating_point T>
inline void validate_selection(
    const thermodynamics::Sw92SelectedPhase<T>& selection) {
    if (!std::isfinite(
            static_cast<double>(
                selection.nacl_molality_mol_per_kg_water)) ||
        selection.nacl_molality_mol_per_kg_water < T{0}) {
        throw std::invalid_argument(
            "mpmc::flow: SW92 selected-phase selection requires finite nonnegative NaCl molality");
    }
    if (selection.family !=
            thermodynamics::SwPhaseFamily::aqueous &&
        selection.family !=
            thermodynamics::SwPhaseFamily::nonaqueous) {
        throw std::invalid_argument(
            "mpmc::flow: SW92 selected-phase selection has invalid phase family");
    }
}

} // namespace sw92_selected_phase_property_detail

/// Fixed-family/fixed-root SW92 property closure for flow natural variables.
///
/// SW92 owns the selected EOS branch, family-dependent fugacity coefficient and
/// molar density. Transport viscosity and the absolute caloric reference remain
/// an explicit caller-owned provider:
///
/// provider(selection, p, T, x, mixture_molar_mass, molar_density, mass_density)
///
/// The provider must preserve Number so forward AD reaches viscosity/enthalpy
/// without a finite-difference production fallback. NaCl molality, AQ/NA family
/// and root index are frozen in each Sw92SelectedPhase and are never changed by
/// this closure. Specific internal energy is u=h-p/rho_mass.
///
/// This v1 bridge is pc=none: every active selected phase is evaluated at the
/// single natural-variable reference pressure.
template <std::floating_point T, typename TransportCaloricProvider>
class Sw92SelectedPhasePropertyClosure {
public:
    Sw92SelectedPhasePropertyClosure(
        const thermodynamics::Sw92Phase<T>& model,
        std::vector<thermodynamics::Sw92SelectedPhase<T>>
            selections,
        TransportCaloricProvider provider,
        SelectedPhasePropertyProvenance provenance)
        : model_(&model),
          selections_(std::move(selections)),
          provider_(std::move(provider)),
          provenance_(std::move(provenance)) {
        if (selections_.empty() ||
            selections_.size() > fixed_three_phase_count) {
            throw std::invalid_argument(
                "mpmc::flow: SW92 selected-phase closure requires 1..3 frozen selections");
        }
        sw92_selected_phase_property_detail::
            validate_provenance(provenance_);
        for (const auto& selection : selections_) {
            sw92_selected_phase_property_detail::
                validate_selection(selection);
        }

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
                    "mpmc::flow: SW92 selected-phase mass density requires an explicit positive kg/mol molar mass for every ordered component");
            }
            molar_mass_kg_per_mol_.push_back(
                component.molar_mass->value);
        }
        if (component_ids_.size() < 2U) {
            throw std::invalid_argument(
                "mpmc::flow: SW92 selected-phase closure requires at least two components");
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

    [[nodiscard]]
    const thermodynamics::Sw92SelectedPhase<T>&
    selection(std::size_t phase) const {
        return selections_.at(phase);
    }

    template <typename Number>
    [[nodiscard]]
    Sw92SelectedPhasePropertyValues<Number>
    evaluate(
        std::size_t phase,
        const Number& pressure_pa,
        const Number& temperature_k,
        std::span<const Number> composition) const {
        using namespace
            sw92_selected_phase_property_detail;

        if (phase >= selections_.size()) {
            throw std::out_of_range(
                "mpmc::flow: SW92 selected-phase selection index out of range");
        }
        if (composition.size() !=
            component_ids_.size()) {
            throw std::invalid_argument(
                "mpmc::flow: SW92 selected-phase composition size mismatch");
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

        const auto& selected =
            selections_[phase];

        thermodynamics::Sw92PhaseWorkspace<Number>
            density_workspace;
        const auto density =
            thermodynamics::
                evaluate_selected_phase_molar_density(
                    *model_,
                    pressure_pa,
                    temperature_k,
                    composition,
                    selected,
                    density_workspace);

        thermodynamics::Sw92PhaseWorkspace<Number>
            fugacity_workspace;
        auto fugacity =
            thermodynamics::
                evaluate_selected_phase_fugacity(
                    *model_,
                    pressure_pa,
                    temperature_k,
                    composition,
                    selected,
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
                selected,
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
                "mpmc::flow: SW92 selected-phase fugacity component count mismatch");
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
    const thermodynamics::Sw92Phase<T>* model_;
    std::vector<thermodynamics::Sw92SelectedPhase<T>>
        selections_;
    TransportCaloricProvider provider_;
    SelectedPhasePropertyProvenance provenance_;
    std::vector<std::string> component_ids_;
    std::vector<double> molar_mass_kg_per_mol_;
};

template <std::floating_point T, typename Provider>
[[nodiscard]] inline auto
make_sw92_selected_phase_property_closure(
    const thermodynamics::Sw92Phase<T>& model,
    std::vector<thermodynamics::Sw92SelectedPhase<T>>
        selections,
    Provider provider,
    SelectedPhasePropertyProvenance provenance) {
    return Sw92SelectedPhasePropertyClosure<
        T,
        Provider>{
        model,
        std::move(selections),
        std::move(provider),
        std::move(provenance)};
}

/// Reuse the already audited topology-neutral property/Jacobian chart driver.
/// The driver is closure-generic; this wrapper keeps SW92 model identity in the
/// public flow boundary.
template <
    std::size_t K = 4U,
    std::floating_point T,
    typename Provider>
[[nodiscard]] inline
Sw92SelectedPhasePropertyChartLinearization
evaluate_sw92_selected_phase_property_chart(
    const NaturalVariableLayoutDescriptor& layout,
    std::span<const std::string> component_ids,
    std::span<const double> natural_variables,
    const Sw92SelectedPhasePropertyClosure<
        T,
        Provider>& closure,
    mpmc::ad::RuntimeJacobianWorkspace<
        double,
        K>& workspace) {
    return evaluate_pr76_selected_phase_property_chart<K>(
        layout,
        component_ids,
        natural_variables,
        closure,
        workspace);
}

[[nodiscard]] inline
Sw92SelectedPhaseFlowLinearization1P
make_sw92_selected_phase_flow_linearization_1p(
    const Sw92SelectedPhasePropertyChartLinearization&
        properties,
    std::span<const double> natural_variables,
    double relative_permeability,
    std::vector<double>
        relative_permeability_gradient) {
    return make_pr76_selected_phase_flow_linearization_1p(
        properties,
        natural_variables,
        relative_permeability,
        std::move(relative_permeability_gradient));
}

[[nodiscard]] inline
Sw92SelectedPhaseFlowLinearization2P
make_sw92_selected_phase_flow_linearization_2p(
    const Sw92SelectedPhasePropertyChartLinearization&
        properties,
    std::span<const double> natural_variables,
    std::array<double, 2> relative_permeability,
    std::array<std::vector<double>, 2>
        relative_permeability_gradient) {
    return make_pr76_selected_phase_flow_linearization_2p(
        properties,
        natural_variables,
        relative_permeability,
        std::move(relative_permeability_gradient));
}

[[nodiscard]] inline
Sw92SelectedPhaseFlowLinearization3P
make_sw92_selected_phase_flow_linearization_3p(
    const Sw92SelectedPhasePropertyChartLinearization&
        properties,
    std::span<const double> natural_variables) {
    return make_pr76_selected_phase_flow_linearization_3p(
        properties,
        natural_variables);
}

/// The SW92 flow property bridge currently evaluates every active selected
/// family/root at p_ref. Any nonzero capillary/phase-pressure value or any
/// noncanonical phase-pressure Jacobian is an unsupported capability request,
/// not a nonlinear domain excursion.
inline void
require_sw92_selected_phase_pc_none_capability(
    const NaturalVariableCellState3P& state,
    const ThreePhaseSaturationConstitutiveNaturalVariableLinearization3P&
        constitutive) {
    phase_transport_detail::
        validate_saturation_linearization(
            state,
            constitutive);

    const auto& layout = state.layout();
    const std::size_t q =
        layout.unknown_count();
    const std::size_t pressure_column =
        layout.pressure_unknown_index();

    for (std::size_t phase = 0U;
         phase < fixed_three_phase_count;
         ++phase) {
        if (constitutive
                    .capillary_pressure_offset_pa[phase] !=
                0.0 ||
            constitutive.phase_pressure_pa[phase] !=
                state.reference_pressure_pa()) {
            throw Sw92SelectedPhasePcNoneCapabilityError{
                "mpmc::flow: SW92 selected-phase pc=none bridge does not support nonzero capillary/phase-pressure offset"};
        }

        const auto& gradient =
            constitutive.phase_pressure_gradient[phase];
        if (gradient.size() != q) {
            throw std::invalid_argument(
                "mpmc::flow: SW92 selected-phase pc=none phase-pressure Jacobian shape mismatch");
        }
        for (std::size_t column = 0U;
             column < q;
             ++column) {
            const double expected =
                column == pressure_column
                    ? 1.0
                    : 0.0;
            if (gradient[column] != expected) {
                throw Sw92SelectedPhasePcNoneCapabilityError{
                    "mpmc::flow: SW92 selected-phase pc=none bridge does not support nonzero capillary/phase-pressure Jacobian"};
            }
        }
    }
}

} // namespace mpmc::flow

#endif // MPMC_FLOW_SW92_SELECTED_PHASE_PROPERTY_CLOSURE_HPP
