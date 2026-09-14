#ifndef MPMC_MODEL_CONFIGURATION_PT_SOLVER_SETTINGS_HPP
#define MPMC_MODEL_CONFIGURATION_PT_SOLVER_SETTINGS_HPP

#include <mpmc/model_configuration/model_definition.hpp>

#include <cstdint>
#include <optional>
#include <string>
#include <string_view>

namespace mpmc::model_configuration {

inline constexpr std::string_view pt_solver_settings_v1 = "pt-solver-settings/v1";
inline constexpr std::string_view mpmc_balanced_default_v1 = "mpmc-balanced-default/v1";
inline constexpr std::string_view mpmc_balanced_default_v1_label = "MPMC balanced default v1";
enum class PtSolverSettingsKind { unspecified, preset, custom };

// Every numeric/boolean field has presence. Zero/false are explicit values, not
// requests to inherit a backend default. Counts are checked before native casts.
struct PtEosRootSettings {
    std::optional<std::int64_t> max_iterations;
    bool operator==(const PtEosRootSettings&) const = default;
};
struct PtStabilitySettings {
    std::optional<double> tpd_tolerance; // Absolute, dimensionless D/(RT).
    std::optional<double> stationarity_tolerance; // Max chemical-potential residual/(RT).
    std::optional<double> line_search_armijo_coefficient;
    std::optional<double> max_log_composition_step;
    std::optional<std::int64_t> max_iterations;
    std::optional<std::int64_t> max_backtracks;
    std::optional<std::int64_t> max_property_evaluations;
    std::optional<bool> automatic_multistart;
    std::optional<std::int64_t> max_starts; // Includes automatically generated starts.
    bool operator==(const PtStabilitySettings&) const = default;
};
struct PtTwoPhaseSettings {
    std::optional<double> fugacity_equilibrium_tolerance; // Max |ln(f1/f2)|.
    std::optional<double> absolute_mass_balance_tolerance;
    std::optional<double> relative_mass_balance_tolerance; // Per positive feed component.
    std::optional<double> minimum_phase_fraction;
    std::optional<double> minimum_log_composition_separation; // Max |ln(y_i/x_i)|.
    std::optional<double> minimum_relative_z_separation;
    std::optional<double> max_log_equilibrium_ratio_step;
    std::optional<double> residual_progress_coefficient;
    std::optional<double> gibbs_progress_coefficient;
    std::optional<std::int64_t> max_iterations;
    std::optional<std::int64_t> max_backtracks;
    std::optional<std::int64_t> max_property_evaluations;
    std::optional<std::int64_t> rachford_rice_max_iterations;
    std::optional<std::int64_t> max_split_attempts;
    bool operator==(const PtTwoPhaseSettings&) const = default;
};
struct PtThreePhaseSettings {
    std::optional<double> chemical_potential_tolerance;
    std::optional<double> absolute_mass_balance_tolerance;
    std::optional<double> relative_mass_balance_tolerance;
    std::optional<double> generalized_rr_balance_tolerance;
    std::optional<double> minimum_phase_fraction;
    std::optional<double> minimum_log_composition_separation;
    std::optional<double> max_log_step;
    std::optional<double> residual_progress_coefficient;
    std::optional<std::int64_t> max_iterations;
    std::optional<std::int64_t> max_line_search_backtracks;
    std::optional<std::int64_t> max_balance_iterations;
    std::optional<std::int64_t> max_balance_backtracks;
    std::optional<std::int64_t> max_property_evaluations;
    std::optional<std::int64_t> max_three_phase_attempts; // Independent budget per seed source class.
    std::optional<double> new_phase_seed_fraction;
    bool operator==(const PtThreePhaseSettings&) const = default;
};

// Owning scientific/numerical draft. Public settings contain neither workspace,
// provider/root branch selectors nor backend implementation types.
struct PtSolverSettings {
    std::string version;
    PtSolverSettingsKind kind = PtSolverSettingsKind::unspecified;
    std::string preset_id; // Exact ID for preset; must be empty for a custom snapshot.
    PtEosRootSettings eos_root;
    PtStabilitySettings initial_stability;
    PtTwoPhaseSettings two_phase;
    PtStabilitySettings final_two_phase_stability;
    PtThreePhaseSettings three_phase;
    PtStabilitySettings final_three_phase_stability;
    bool operator==(const PtSolverSettings&) const = default;
};

namespace detail {
// Frozen literal values, audited against the PR76 defaults at 9d856b8. Never
// construct a backend Options object here to discover its current defaults.
inline PtStabilitySettings balanced_stability_v1() {
    return {.tpd_tolerance = 1e-10,
            .stationarity_tolerance = 1e-8,
            .line_search_armijo_coefficient = 1e-4,
            .max_log_composition_step = 4.0,
            .max_iterations = 512,
            .max_backtracks = 32,
            .max_property_evaluations = 100000,
            .automatic_multistart = true,
            .max_starts = 1024};
}
} // namespace detail

// Query a complete resolved preset. The return value owns three independent
// stability groups. Edits must be relabeled as custom before preparation.
[[nodiscard]] inline PtSolverSettings resolve_pt_solver_preset(std::string_view id) {
    if (id != mpmc_balanced_default_v1) {
        throw ModelConfigurationError(ModelConfigurationErrorCode::unsupported_preset,
                                      "preset_id", "unknown solver preset; no fallback");
    }
    return {
        .version = std::string(pt_solver_settings_v1),
        .kind = PtSolverSettingsKind::preset,
        .preset_id = std::string(mpmc_balanced_default_v1),
        .eos_root = {.max_iterations = 2048},
        .initial_stability = detail::balanced_stability_v1(),
        .two_phase = {
            .fugacity_equilibrium_tolerance = 1e-11,
            .absolute_mass_balance_tolerance = 1e-12,
            .relative_mass_balance_tolerance = 1e-10,
            .minimum_phase_fraction = 1e-10,
            .minimum_log_composition_separation = 1e-7,
            .minimum_relative_z_separation = 1e-8,
            .max_log_equilibrium_ratio_step = 2.0,
            .residual_progress_coefficient = 1e-4,
            .gibbs_progress_coefficient = 1e-4,
            .max_iterations = 512,
            .max_backtracks = 24,
            .max_property_evaluations = 20000,
            .rachford_rice_max_iterations = 192,
            .max_split_attempts = 16},
        .final_two_phase_stability = detail::balanced_stability_v1(),
        .three_phase = {
            .chemical_potential_tolerance = 1e-11,
            .absolute_mass_balance_tolerance = 1e-12,
            .relative_mass_balance_tolerance = 1e-10,
            .generalized_rr_balance_tolerance = 2e-13,
            .minimum_phase_fraction = 1e-10,
            .minimum_log_composition_separation = 1e-7,
            .max_log_step = 2.0,
            .residual_progress_coefficient = 1e-4,
            .max_iterations = 512,
            .max_line_search_backtracks = 32,
            .max_balance_iterations = 192,
            .max_balance_backtracks = 48,
            .max_property_evaluations = 30000,
            .max_three_phase_attempts = 16,
            .new_phase_seed_fraction = 0.1},
        .final_three_phase_stability = detail::balanced_stability_v1()};
}

// Host policy, not solver-preset values and never a client configuration field.
// These bound nested operations individually; they are not a wall-clock deadline
// or one aggregate property-call budget for the entire flash.
struct PtSolverSafetyLimits {
    std::size_t max_root_iterations{8192};
    std::size_t max_iterations{4096};
    std::size_t max_backtracks{256};
    std::size_t max_property_evaluations{1000000};
    std::size_t max_stability_starts{1024};
    std::size_t max_split_attempts{64};
    std::size_t max_three_phase_attempts{64};
    std::size_t max_components{256};
    std::size_t max_stability_start_entries{262144};
    std::size_t max_three_phase_starts{16};
    std::size_t max_three_phase_start_entries{12288};
    bool operator==(const PtSolverSafetyLimits&) const = default;
};

} // namespace mpmc::model_configuration
#endif // MPMC_MODEL_CONFIGURATION_PT_SOLVER_SETTINGS_HPP
