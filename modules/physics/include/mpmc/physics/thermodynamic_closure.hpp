#ifndef MPMC_PHYSICS_THERMODYNAMIC_CLOSURE_HPP
#define MPMC_PHYSICS_THERMODYNAMIC_CLOSURE_HPP

#include <cstddef>
#include <optional>
#include <stdexcept>
#include <string>
#include <string_view>
#include <vector>

namespace mpmc::physics {

// Two orthogonal decisions are retained. A valid primal may be used for a
// residual/merit evaluation even when a derivative cannot seed the next Newton
// base point. This layer contains no grid, flux or nonlinear-solver policy.
enum class ThermodynamicClosurePrimalStatus {
    valid,
    phase_regime_not_implemented,
    indeterminate
};

enum class ThermodynamicClosureLinearizationStatus { available, unavailable };

enum class ThermodynamicClosureLinearizationReason {
    none,
    phase_boundary,
    ill_conditioned_equilibrium,
    unsupported_feed_support,
    property_failure,
    arithmetic_failure,
    solution_not_accepted
};

struct ThermodynamicPhaseState {
    double mole_phase_fraction{}; // Mole phase fraction, NOT pore-volume saturation.
    std::vector<double> composition; // Mole fractions in the ordered component snapshot.
    double compressibility_factor{}; // Z = p*v_m/(R*T), dimensionless.
    double molar_density_mol_per_m3{}; // c = 1/v_m = p/(Z*R*T).
};

struct PtVleThermodynamicState {
    ThermodynamicPhaseState liquid;
    ThermodynamicPhaseState vapor;
};

struct PtVleThermodynamicLinearization {
    std::size_t component_count{};
    std::size_t input_count{}; // N+1: p, T, then N-1 independent overall mole fractions.

    std::vector<double> vapor_fraction_gradient;
    std::vector<double> liquid_composition_jacobian;
    std::vector<double> vapor_composition_jacobian;
    std::vector<double> liquid_compressibility_gradient;
    std::vector<double> vapor_compressibility_gradient;
    std::vector<double> liquid_molar_density_gradient;
    std::vector<double> vapor_molar_density_gradient;

    [[nodiscard]] double d_vapor_fraction(std::size_t column) const {
        return vapor_fraction_gradient.at(column);
    }
    [[nodiscard]] double d_liquid_composition(
        std::size_t component, std::size_t column) const {
        return liquid_composition_jacobian.at(component * input_count + column);
    }
    [[nodiscard]] double d_vapor_composition(
        std::size_t component, std::size_t column) const {
        return vapor_composition_jacobian.at(component * input_count + column);
    }
    [[nodiscard]] double d_liquid_compressibility(std::size_t column) const {
        return liquid_compressibility_gradient.at(column);
    }
    [[nodiscard]] double d_vapor_compressibility(std::size_t column) const {
        return vapor_compressibility_gradient.at(column);
    }
    [[nodiscard]] double d_liquid_molar_density(std::size_t column) const {
        return liquid_molar_density_gradient.at(column);
    }
    [[nodiscard]] double d_vapor_molar_density(std::size_t column) const {
        return vapor_molar_density_gradient.at(column);
    }
};

// Atomic snapshot: every build creates a fresh optional linearization. A caller
// replacing this object cannot accidentally retain the derivative from a prior
// thermodynamic state inside the same snapshot.
struct ThermodynamicClosureSnapshot {
    static constexpr std::string_view convention =
        "PT-VLE/thermodynamic-closure/reduced-feed-v1";

    ThermodynamicClosurePrimalStatus primal_status{
        ThermodynamicClosurePrimalStatus::indeterminate};
    ThermodynamicClosureLinearizationStatus linearization_status{
        ThermodynamicClosureLinearizationStatus::unavailable};
    ThermodynamicClosureLinearizationReason linearization_reason{
        ThermodynamicClosureLinearizationReason::none};

    double pressure_pa{};
    double temperature_k{};
    std::vector<double> feed;
    std::vector<std::string> component_ids;
    std::string thermodynamic_model;
    std::string dataset_id;
    std::string revision;

    std::optional<PtVleThermodynamicState> primal;
    std::optional<PtVleThermodynamicLinearization> linearization;
    std::string diagnostic;

    [[nodiscard]] bool residual_available() const noexcept {
        return primal_status == ThermodynamicClosurePrimalStatus::valid &&
               primal.has_value();
    }
    [[nodiscard]] bool can_seed_newton() const noexcept {
        return residual_available() &&
               linearization_status ==
                   ThermodynamicClosureLinearizationStatus::available &&
               linearization.has_value();
    }
};

} // namespace mpmc::physics

#endif // MPMC_PHYSICS_THERMODYNAMIC_CLOSURE_HPP
