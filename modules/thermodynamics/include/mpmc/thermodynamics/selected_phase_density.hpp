#ifndef MPMC_THERMODYNAMICS_SELECTED_PHASE_DENSITY_HPP
#define MPMC_THERMODYNAMICS_SELECTED_PHASE_DENSITY_HPP

#include <mpmc/thermodynamics/selected_phase_fugacity.hpp>

#include <cmath>
#include <concepts>
#include <cstddef>
#include <span>
#include <stdexcept>
#include <string_view>

namespace mpmc::thermodynamics {

inline constexpr std::string_view
    selected_phase_density_convention =
        "thermodynamics/selected-phase-molar-density/v1";

template <typename Number>
struct SelectedPhaseMolarDensityValues {
    Number molar_density_mol_per_m3;
};

template <std::floating_point T, typename Number>
[[nodiscard]] inline SelectedPhaseMolarDensityValues<Number>
evaluate_selected_phase_molar_density(
    const Pr76Phase<T>& model,
    const Number& pressure_pa,
    const Number& temperature_k,
    std::span<const Number> composition,
    const Pr76SelectedPhase& selection,
    Pr76PhaseWorkspace<Number>& workspace) {
    const auto phase = model.evaluate_full(
        pressure_pa,
        temperature_k,
        composition,
        selection.root_index,
        workspace,
        selection.root_options);
    const Number denominator =
        phase.z *
        Pr76Pure<T>::gas_constant() *
        temperature_k;
    const Number density =
        pressure_pa / denominator;
    if constexpr (std::floating_point<Number>) {
        if (!std::isfinite(density) ||
            !(density > Number{0})) {
            throw std::range_error(
                "PR76 selected phase density is non-finite or non-positive");
        }
    }
    return {density};
}

template <std::floating_point T, typename Number>
[[nodiscard]] inline SelectedPhaseMolarDensityValues<Number>
evaluate_selected_phase_molar_density(
    const Sw92Phase<T>& model,
    const Number& pressure_pa,
    const Number& temperature_k,
    std::span<const Number> composition,
    const Sw92SelectedPhase<T>& selection,
    Sw92PhaseWorkspace<Number>& workspace) {
    const auto phase = model.evaluate_full(
        pressure_pa,
        temperature_k,
        composition,
        selection.nacl_molality_mol_per_kg_water,
        selection.family,
        selection.root_index,
        workspace,
        selection.root_options);
    const Number denominator =
        phase.z *
        Sw92Pure<T>::gas_constant() *
        temperature_k;
    const Number density =
        pressure_pa / denominator;
    if constexpr (std::floating_point<Number>) {
        if (!std::isfinite(density) ||
            !(density > Number{0})) {
            throw std::range_error(
                "SW92 selected phase density is non-finite or non-positive");
        }
    }
    return {density};
}

[[nodiscard]] inline SelectedPhaseMolarDensityValues<double>
evaluate_selected_phase_molar_density(
    const CpaPtPhase& model,
    double pressure_pa,
    double temperature_k,
    std::span<const double> composition,
    const CpaSelectedPhase& selection) {
    const auto roots = model.roots(
        pressure_pa,
        temperature_k,
        composition,
        selection.options);
    if (roots.status != CpaPtRootStatus::success) {
        if (roots.status == CpaPtRootStatus::near_multiple) {
            throw std::domain_error(
                "CPA selected phase density: near-multiple/tangent PT root is not a fixed smooth branch");
        }
        throw std::runtime_error(
            "CPA selected phase density: PT root set unavailable: " +
            roots.diagnostic);
    }
    if (selection.root_index >= roots.roots.size()) {
        throw std::out_of_range(
            "CPA selected phase density: selected root index out of range");
    }
    const double density =
        roots.roots[selection.root_index]
            .molar_density_mol_per_m3;
    if (!std::isfinite(density) ||
        !(density > 0.0)) {
        throw std::range_error(
            "CPA selected phase density is non-finite or non-positive");
    }
    return {density};
}

template <std::size_t K>
[[nodiscard]] inline
SelectedPhaseMolarDensityValues<
    mpmc::ad::Dual<double, K>>
evaluate_selected_phase_molar_density(
    const CpaPtPhase& model,
    const mpmc::ad::Dual<double, K>& pressure_pa,
    const mpmc::ad::Dual<double, K>& temperature_k,
    std::span<const mpmc::ad::Dual<double, K>> composition,
    const CpaSelectedPhase& selection) {
    return {
        evaluate_cpa_selected_pt_molar_density_first_order(
            model,
            pressure_pa,
            temperature_k,
            composition,
            selection.root_index,
            selection.options)};
}

} // namespace mpmc::thermodynamics

#endif // MPMC_THERMODYNAMICS_SELECTED_PHASE_DENSITY_HPP
