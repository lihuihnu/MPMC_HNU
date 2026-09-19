#ifndef MPMC_THERMODYNAMICS_SELECTED_PHASE_FUGACITY_HPP
#define MPMC_THERMODYNAMICS_SELECTED_PHASE_FUGACITY_HPP

#include <mpmc/thermodynamics/cpa_pt_phase.hpp>
#include <mpmc/thermodynamics/pr76_phase.hpp>
#include <mpmc/thermodynamics/sw92_phase.hpp>

#include <concepts>
#include <cstddef>
#include <span>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace mpmc::thermodynamics {

enum class SelectedPhaseFugacityDerivativeSupport {
    value_only,
    scalar_generic_first_order
};

template <typename Model>
struct SelectedPhaseFugacityCapabilities;

template <std::floating_point T>
struct SelectedPhaseFugacityCapabilities<Pr76Phase<T>> {
    static constexpr SelectedPhaseFugacityDerivativeSupport derivative_support =
        SelectedPhaseFugacityDerivativeSupport::scalar_generic_first_order;
};

template <std::floating_point T>
struct SelectedPhaseFugacityCapabilities<Sw92Phase<T>> {
    static constexpr SelectedPhaseFugacityDerivativeSupport derivative_support =
        SelectedPhaseFugacityDerivativeSupport::scalar_generic_first_order;
};

template <>
struct SelectedPhaseFugacityCapabilities<CpaPtPhase> {
    static constexpr SelectedPhaseFugacityDerivativeSupport derivative_support =
        SelectedPhaseFugacityDerivativeSupport::value_only;
};

template <typename Number>
struct SelectedPhaseFugacityValues {
    std::vector<Number> ln_phi;
};

struct Pr76SelectedPhase {
    std::size_t root_index{};
    Pr76RootOptions root_options{};
};

template <std::floating_point T>
struct Sw92SelectedPhase {
    T nacl_molality_mol_per_kg_water{};
    SwPhaseFamily family{SwPhaseFamily::nonaqueous};
    std::size_t root_index{};
    Sw92RootOptions root_options{};
};

struct CpaSelectedPhase {
    std::size_t root_index{};
    CpaPtOptions options{};
};

template <std::floating_point T, typename Number>
[[nodiscard]] inline SelectedPhaseFugacityValues<Number>
evaluate_selected_phase_fugacity(
    const Pr76Phase<T>& model,
    const Number& pressure_pa,
    const Number& temperature_k,
    std::span<const Number> composition,
    const Pr76SelectedPhase& selection,
    Pr76PhaseWorkspace<Number>& workspace) {
    auto values = model.evaluate_full(
        pressure_pa,
        temperature_k,
        composition,
        selection.root_index,
        workspace,
        selection.root_options);
    return {std::move(values.ln_phi)};
}

template <std::floating_point T, typename Number>
[[nodiscard]] inline SelectedPhaseFugacityValues<Number>
evaluate_selected_phase_fugacity(
    const Sw92Phase<T>& model,
    const Number& pressure_pa,
    const Number& temperature_k,
    std::span<const Number> composition,
    const Sw92SelectedPhase<T>& selection,
    Sw92PhaseWorkspace<Number>& workspace) {
    auto values = model.evaluate_full(
        pressure_pa,
        temperature_k,
        composition,
        selection.nacl_molality_mol_per_kg_water,
        selection.family,
        selection.root_index,
        workspace,
        selection.root_options);
    return {std::move(values.ln_phi)};
}

[[nodiscard]] inline SelectedPhaseFugacityValues<double>
evaluate_selected_phase_fugacity(
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
                "CPA selected phase fugacity: near-multiple/tangent PT root is not a fixed smooth branch");
        }
        throw std::runtime_error(
            "CPA selected phase fugacity: PT root set unavailable: " +
            roots.diagnostic);
    }
    if (selection.root_index >= roots.roots.size()) {
        throw std::out_of_range(
            "CPA selected phase fugacity: selected root index out of range");
    }
    return {roots.roots[selection.root_index].ln_phi};
}

} // namespace mpmc::thermodynamics

#endif // MPMC_THERMODYNAMICS_SELECTED_PHASE_FUGACITY_HPP
