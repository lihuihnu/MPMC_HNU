#ifndef MPMC_PHYSICS_THERMODYNAMIC_CLOSURE_HPP
#define MPMC_PHYSICS_THERMODYNAMIC_CLOSURE_HPP

#include <cmath>
#include <cstddef>
#include <limits>
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
    not_implemented,
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

// Variable-cardinality PT phase-set primal. Phase order is representation order
// supplied by the accepted flash result; this generic payload deliberately does
// not invent liquid/vapor/aqueous labels. Model-specific role/family metadata
// belongs in a model adapter sidecar.
struct PtPhaseSetThermodynamicState {
    std::vector<ThermodynamicPhaseState> phases;
};

// Generic reduced-feed linearization for one fixed accepted phase-set
// representation. Columns are always
//     q = (p_Pa, T_K, z_0, ..., z_{N-2})
// with z_{N-1}=1-sum(z_0..z_{N-2}). Phase order is the exact representation
// order of the accompanying primal and has no generic liquid/vapor meaning.
struct PtPhaseSetThermodynamicLinearization {
    std::size_t component_count{};
    std::size_t phase_count{};
    std::size_t input_count{}; // N+1 reduced-feed coordinates.

    std::vector<double> phase_fraction_jacobian; // [phase * input_count + q]
    std::vector<double> composition_jacobian; // [(phase*N+i)*input_count + q]
    std::vector<double> compressibility_jacobian; // [phase * input_count + q]
    std::vector<double> molar_density_jacobian; // [phase * input_count + q]

    // Numerical diagnostics belong to the local implicit solve that produced
    // this linearization. They are not global-stability certificates.
    double equilibrium_jacobian_rcond{};
    double linear_solve_backward_error{};
    double equilibrium_residual_norm{};

    [[nodiscard]] static constexpr std::size_t pressure_column() noexcept {
        return 0U;
    }
    [[nodiscard]] static constexpr std::size_t temperature_column() noexcept {
        return 1U;
    }
    [[nodiscard]] std::size_t dependent_feed_component() const {
        if (component_count < 2U) {
            throw std::logic_error(
                "phase-set linearization: no reduced-feed coordinate chart");
        }
        return component_count - 1U;
    }
    [[nodiscard]] std::size_t feed_column(
        std::size_t independent_component) const {
        if (component_count < 2U ||
            independent_component + 1U >= component_count) {
            throw std::out_of_range(
                "phase-set linearization: feed component outside reduced chart");
        }
        return 2U + independent_component;
    }
    [[nodiscard]] double d_phase_fraction(
        std::size_t phase, std::size_t column) const {
        return phase_fraction_jacobian.at(phase * input_count + column);
    }
    [[nodiscard]] double d_composition(
        std::size_t phase, std::size_t component,
        std::size_t column) const {
        return composition_jacobian.at(
            (phase * component_count + component) * input_count + column);
    }
    [[nodiscard]] double d_compressibility(
        std::size_t phase, std::size_t column) const {
        return compressibility_jacobian.at(phase * input_count + column);
    }
    [[nodiscard]] double d_molar_density(
        std::size_t phase, std::size_t column) const {
        return molar_density_jacobian.at(phase * input_count + column);
    }
};

// Variable-cardinality v2 contract. Primal and local linearization availability
// are independent, but when a linearization is published it is owned by the
// same snapshot so stale derivatives cannot survive a primal replacement.
struct PtPhaseSetThermodynamicClosureSnapshot {
    static constexpr std::string_view convention =
        "PT/phase-set/thermodynamic-closure/reduced-feed-v2";

    ThermodynamicClosurePrimalStatus primal_status{
        ThermodynamicClosurePrimalStatus::indeterminate};
    ThermodynamicClosureLinearizationStatus linearization_status{
        ThermodynamicClosureLinearizationStatus::unavailable};
    ThermodynamicClosureLinearizationReason linearization_reason{
        ThermodynamicClosureLinearizationReason::not_implemented};

    double pressure_pa{};
    double temperature_k{};
    std::vector<double> feed;
    std::vector<std::string> component_ids;
    std::string thermodynamic_model;
    std::string dataset_id;
    std::string revision;

    std::optional<PtPhaseSetThermodynamicState> primal;
    std::optional<PtPhaseSetThermodynamicLinearization> linearization;
    std::string diagnostic;

    [[nodiscard]] bool residual_available() const noexcept {
        return primal_status == ThermodynamicClosurePrimalStatus::valid &&
               primal.has_value() && !primal->phases.empty();
    }

    [[nodiscard]] bool can_seed_newton() const noexcept {
        if (!residual_available() ||
            linearization_status !=
                ThermodynamicClosureLinearizationStatus::available ||
            linearization_reason != ThermodynamicClosureLinearizationReason::none ||
            !linearization.has_value()) {
            return false;
        }
        const auto& value = *linearization;
        if (value.component_count < 2U ||
            value.component_count == std::numeric_limits<std::size_t>::max() ||
            value.component_count != component_ids.size() ||
            value.phase_count != primal->phases.size() ||
            value.input_count != value.component_count + 1U) {
            return false;
        }
        for (const auto& phase : primal->phases) {
            if (phase.composition.size() != value.component_count) { return false; }
        }
        const auto product_matches = [](std::size_t lhs, std::size_t rhs,
                                        std::size_t actual) noexcept {
            if (lhs != 0U && rhs > std::numeric_limits<std::size_t>::max() / lhs) {
                return false;
            }
            return lhs * rhs == actual;
        };
        if (!product_matches(value.phase_count, value.input_count,
                             value.phase_fraction_jacobian.size()) ||
            !product_matches(value.phase_count, value.input_count,
                             value.compressibility_jacobian.size()) ||
            !product_matches(value.phase_count, value.input_count,
                             value.molar_density_jacobian.size())) {
            return false;
        }
        if (value.phase_count != 0U &&
            value.component_count >
                std::numeric_limits<std::size_t>::max() / value.phase_count) {
            return false;
        }
        const std::size_t phase_component_count =
            value.phase_count * value.component_count;
        if (!product_matches(phase_component_count, value.input_count,
                             value.composition_jacobian.size())) {
            return false;
        }
        const auto finite_vector = [](const std::vector<double>& values) noexcept {
            for (double entry : values) {
                if (!std::isfinite(entry)) { return false; }
            }
            return true;
        };
        return finite_vector(value.phase_fraction_jacobian) &&
               finite_vector(value.composition_jacobian) &&
               finite_vector(value.compressibility_jacobian) &&
               finite_vector(value.molar_density_jacobian) &&
               std::isfinite(value.equilibrium_jacobian_rcond) &&
               value.equilibrium_jacobian_rcond > 0.0 &&
               std::isfinite(value.linear_solve_backward_error) &&
               value.linear_solve_backward_error >= 0.0 &&
               std::isfinite(value.equilibrium_residual_norm) &&
               value.equilibrium_residual_norm >= 0.0;
    }
};

} // namespace mpmc::physics

#endif // MPMC_PHYSICS_THERMODYNAMIC_CLOSURE_HPP
