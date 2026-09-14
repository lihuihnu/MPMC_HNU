#ifndef MPMC_MODEL_CONFIGURATION_PR76_SOLVER_SETTINGS_HPP
#define MPMC_MODEL_CONFIGURATION_PR76_SOLVER_SETTINGS_HPP

#include <mpmc/model_configuration/pt_solver_settings.hpp>
#include <mpmc/flash/pr76_pt_flash_backend.hpp>

#include <cmath>
#include <limits>
#include <utility>

namespace mpmc::model_configuration {
namespace detail {

template <typename T>
inline T required_solver_value(const std::optional<T>& value, const std::string& field) {
    if (!value) {
        throw ModelConfigurationError(ModelConfigurationErrorCode::missing_field,
                                      field, "explicit solver setting required");
    }
    return *value;
}
inline double solver_number(const std::optional<double>& value, const std::string& field) {
    const double number = required_solver_value(value, field);
    if (!std::isfinite(number)) {
        throw ModelConfigurationError(ModelConfigurationErrorCode::invalid_value,
                                      field, "finite solver setting required");
    }
    return number;
}
template <typename T>
inline T solver_count(const std::optional<std::int64_t>& value, std::size_t ceiling,
                      const std::string& field) {
    const auto count = required_solver_value(value, field);
    if (count < 0 || !std::in_range<T>(count)) {
        throw ModelConfigurationError(ModelConfigurationErrorCode::invalid_value,
                                      field, "count is negative or cannot be represented");
    }
    if (static_cast<std::uint64_t>(count) > ceiling) {
        throw ModelConfigurationError(ModelConfigurationErrorCode::resource_limit,
                                      field, "server solver ceiling exceeded; no clipping");
    }
    return static_cast<T>(count);
}
inline void validate_solver_limits(const PtSolverSafetyLimits& limits) {
    const auto positive = [](std::size_t value, const char* field) {
        if (value == 0) {
            throw ModelConfigurationError(ModelConfigurationErrorCode::invalid_value,
                                          std::string("safety_limits.") + field,
                                          "server safety ceiling must be positive");
        }
    };
    positive(limits.max_root_iterations, "max_root_iterations");
    positive(limits.max_iterations, "max_iterations");
    positive(limits.max_backtracks, "max_backtracks");
    positive(limits.max_property_evaluations, "max_property_evaluations");
    positive(limits.max_stability_starts, "max_stability_starts");
    positive(limits.max_split_attempts, "max_split_attempts");
    positive(limits.max_three_phase_attempts, "max_three_phase_attempts");
    positive(limits.max_components, "max_components");
    positive(limits.max_stability_start_entries, "max_stability_start_entries");
    positive(limits.max_three_phase_starts, "max_three_phase_starts");
    positive(limits.max_three_phase_start_entries, "max_three_phase_start_entries");
}
template <class Check>
inline void check_solver_domain(const std::string& field, Check&& check) {
    // Catch only the documented bad-options exception from the existing option
    // validators. Allocation failures/programming errors are not swallowed.
    try { check(); }
    catch (const std::invalid_argument& error) {
        throw ModelConfigurationError(ModelConfigurationErrorCode::invalid_settings,
                                      field, error.what());
    }
}
inline flash::StabilityOptions map_stability_settings(
    const PtStabilitySettings& s, const PtSolverSafetyLimits& limits, const std::string& field) {
    flash::StabilityOptions result{
        .tpd_tolerance = solver_number(s.tpd_tolerance, field + ".tpd_tolerance"),
        .stationarity_tolerance = solver_number(s.stationarity_tolerance, field + ".stationarity_tolerance"),
        .armijo = solver_number(s.line_search_armijo_coefficient, field + ".line_search_armijo_coefficient"),
        .log_step_limit = solver_number(s.max_log_composition_step, field + ".max_log_composition_step"),
        .max_iterations = solver_count<int>(s.max_iterations, limits.max_iterations, field + ".max_iterations"),
        .max_backtracks = solver_count<int>(s.max_backtracks, limits.max_backtracks, field + ".max_backtracks"),
        .max_evaluations = solver_count<std::size_t>(s.max_property_evaluations, limits.max_property_evaluations,
                                                   field + ".max_property_evaluations"),
        .max_components = limits.max_components,
        .max_starts = solver_count<std::size_t>(s.max_starts, limits.max_stability_starts, field + ".max_starts"),
        .max_start_entries = limits.max_stability_start_entries,
        .automatic_starts = required_solver_value(s.automatic_multistart, field + ".automatic_multistart")};
    check_solver_domain(field, [&] { flash::detail::stability_check_options(result); });
    return result;
}
inline flash::PtSplitOptions map_two_phase_settings(
    const PtSolverSettings& settings, const PtSolverSafetyLimits& limits) {
    const auto& s = settings.two_phase;
    flash::PtSplitIterationOptions iteration{
        .fugacity_tolerance = solver_number(s.fugacity_equilibrium_tolerance, "two_phase.fugacity_equilibrium_tolerance"),
        .mass_absolute_tolerance = solver_number(s.absolute_mass_balance_tolerance, "two_phase.absolute_mass_balance_tolerance"),
        .mass_relative_tolerance = solver_number(s.relative_mass_balance_tolerance, "two_phase.relative_mass_balance_tolerance"),
        .minimum_phase_fraction = solver_number(s.minimum_phase_fraction, "two_phase.minimum_phase_fraction"),
        .log_k_separation = solver_number(s.minimum_log_composition_separation, "two_phase.minimum_log_composition_separation"),
        .relative_z_separation = solver_number(s.minimum_relative_z_separation, "two_phase.minimum_relative_z_separation"),
        .max_log_step = solver_number(s.max_log_equilibrium_ratio_step, "two_phase.max_log_equilibrium_ratio_step"),
        .residual_decrease = solver_number(s.residual_progress_coefficient, "two_phase.residual_progress_coefficient"),
        .max_iterations = solver_count<int>(s.max_iterations, limits.max_iterations, "two_phase.max_iterations"),
        .max_backtracks = solver_count<int>(s.max_backtracks, limits.max_backtracks, "two_phase.max_backtracks"),
        .max_evaluations = solver_count<std::size_t>(s.max_property_evaluations, limits.max_property_evaluations,
                                                   "two_phase.max_property_evaluations"),
        .rr = {.max_iterations = solver_count<int>(s.rachford_rice_max_iterations, limits.max_iterations,
                                                  "two_phase.rachford_rice_max_iterations"),
               .max_components = limits.max_components},
        .gibbs_armijo = solver_number(s.gibbs_progress_coefficient, "two_phase.gibbs_progress_coefficient")};
    check_solver_domain("two_phase", [&] { flash::detail::split_check_options(iteration); });
    return {.iteration = iteration,
            .initial_stability = map_stability_settings(settings.initial_stability, limits, "initial_stability"),
            .final_stability = map_stability_settings(settings.final_two_phase_stability, limits, "final_two_phase_stability"),
            .max_split_attempts = solver_count<std::size_t>(s.max_split_attempts, limits.max_split_attempts,
                                                          "two_phase.max_split_attempts")};
}
inline flash::PtThreePhaseOptions map_three_phase_settings(
    const PtThreePhaseSettings& s, const PtSolverSafetyLimits& limits) {
    flash::PtThreePhaseOptions result{
        .chemical_potential_tolerance = solver_number(s.chemical_potential_tolerance, "three_phase.chemical_potential_tolerance"),
        .mass_absolute_tolerance = solver_number(s.absolute_mass_balance_tolerance, "three_phase.absolute_mass_balance_tolerance"),
        .mass_relative_tolerance = solver_number(s.relative_mass_balance_tolerance, "three_phase.relative_mass_balance_tolerance"),
        .balance_tolerance = solver_number(s.generalized_rr_balance_tolerance, "three_phase.generalized_rr_balance_tolerance"),
        .minimum_phase_fraction = solver_number(s.minimum_phase_fraction, "three_phase.minimum_phase_fraction"),
        .log_composition_separation = solver_number(s.minimum_log_composition_separation, "three_phase.minimum_log_composition_separation"),
        .max_log_step = solver_number(s.max_log_step, "three_phase.max_log_step"),
        .residual_decrease = solver_number(s.residual_progress_coefficient, "three_phase.residual_progress_coefficient"),
        .max_iterations = solver_count<int>(s.max_iterations, limits.max_iterations, "three_phase.max_iterations"),
        .max_backtracks = solver_count<int>(s.max_line_search_backtracks, limits.max_backtracks, "three_phase.max_line_search_backtracks"),
        .max_balance_iterations = solver_count<int>(s.max_balance_iterations, limits.max_iterations, "three_phase.max_balance_iterations"),
        .max_balance_backtracks = solver_count<int>(s.max_balance_backtracks, limits.max_backtracks, "three_phase.max_balance_backtracks"),
        .max_evaluations = solver_count<std::size_t>(s.max_property_evaluations, limits.max_property_evaluations,
                                                   "three_phase.max_property_evaluations"),
        .max_components = limits.max_components};
    check_solver_domain("three_phase", [&] { flash::detail::pt_three_phase_check_options(result); });
    return result;
}
} // namespace detail

// Immutable prepared NUMERICAL configuration, not a backend or model registry.
// settings() retains actual public values and preset/custom identity; host policy
// is queried separately. Returned native options explicitly contain empty hints.
class Pr76SolverConfiguration {
public:
    Pr76SolverConfiguration(const Pr76SolverConfiguration&) = default;
    Pr76SolverConfiguration(Pr76SolverConfiguration&&) noexcept = default;
    Pr76SolverConfiguration& operator=(const Pr76SolverConfiguration&) = delete;
    Pr76SolverConfiguration& operator=(Pr76SolverConfiguration&&) = delete;

    [[nodiscard]] static Pr76SolverConfiguration create(
        const PtSolverSettings& settings, PtSolverSafetyLimits limits = {}) {
        if (settings.version != pt_solver_settings_v1) {
            throw ModelConfigurationError(ModelConfigurationErrorCode::unsupported_version,
                                          "version", "unsupported solver-settings version");
        }
        if (settings.kind == PtSolverSettingsKind::preset) {
            const auto preset = resolve_pt_solver_preset(settings.preset_id);
            if (settings != preset) {
                throw ModelConfigurationError(ModelConfigurationErrorCode::invalid_settings,
                                              "preset_id", "preset was edited; declare complete custom settings");
            }
        } else if (settings.kind != PtSolverSettingsKind::custom || !settings.preset_id.empty()) {
            throw ModelConfigurationError(ModelConfigurationErrorCode::invalid_settings,
                                          "kind", "explicit preset or custom settings identity required");
        }
        detail::validate_solver_limits(limits);
        thermodynamics::Pr76RootOptions roots{
            .max_iterations = detail::solver_count<int>(settings.eos_root.max_iterations,
                limits.max_root_iterations, "eos_root.max_iterations")};
        // Same constructor-domain rule as the existing PR evaluators; no roots
        // are evaluated merely to validate configuration.
        if (roots.max_iterations == 0) {
            throw ModelConfigurationError(ModelConfigurationErrorCode::invalid_settings,
                                          "eos_root.max_iterations", "positive root iteration limit required");
        }
        flash::Pr76PtFlashBackendOptions options{
            .split = detail::map_two_phase_settings(settings, limits),
            .three_phase = detail::map_three_phase_settings(settings.three_phase, limits),
            .final_three_phase_stability = detail::map_stability_settings(
                settings.final_three_phase_stability, limits, "final_three_phase_stability"),
            .max_three_phase_attempts = detail::solver_count<std::size_t>(
                settings.three_phase.max_three_phase_attempts, limits.max_three_phase_attempts,
                "three_phase.max_three_phase_attempts"),
            .new_phase_seed_fraction = detail::solver_number(
                settings.three_phase.new_phase_seed_fraction, "three_phase.new_phase_seed_fraction"),
            .three_phase_starts = {},
            .max_three_phase_starts = limits.max_three_phase_starts,
            .max_three_phase_start_entries = limits.max_three_phase_start_entries,
            .initial_starts = {},
            .final_starts = {}};
        // Existing validator owns seed-fraction/attempt semantics. Do not infer
        // admissibility or phase count here. Feed/hint-dependent checks stay at solve.
        const flash::Pr76PtMax3Options max3{
            .two_phase = options.split,
            .three_phase = options.three_phase,
            .final_three_phase_stability = options.final_three_phase_stability,
            .max_three_phase_attempts = options.max_three_phase_attempts,
            .new_phase_seed_fraction = options.new_phase_seed_fraction,
            .three_phase_starts = {},
            .max_three_phase_starts = options.max_three_phase_starts,
            .max_three_phase_start_entries = options.max_three_phase_start_entries};
        detail::check_solver_domain("three_phase", [&] { flash::detail::pr76_pt_max3_check_options(max3); });
        return Pr76SolverConfiguration(settings, limits, roots, std::move(options));
    }

    [[nodiscard]] const PtSolverSettings& settings() const & noexcept { return settings_; }
    const PtSolverSettings& settings() const && = delete;
    [[nodiscard]] const PtSolverSafetyLimits& safety_limits() const & noexcept { return limits_; }
    const PtSolverSafetyLimits& safety_limits() const && = delete;
    [[nodiscard]] thermodynamics::Pr76RootOptions root_options() const noexcept { return roots_; }
    [[nodiscard]] const flash::Pr76PtFlashBackendOptions& backend_options() const & noexcept { return options_; }
    const flash::Pr76PtFlashBackendOptions& backend_options() const && = delete;
private:
    Pr76SolverConfiguration(PtSolverSettings settings, PtSolverSafetyLimits limits,
                           thermodynamics::Pr76RootOptions roots, flash::Pr76PtFlashBackendOptions options)
        : settings_(std::move(settings)), limits_(limits), roots_(roots), options_(std::move(options)) {}
    PtSolverSettings settings_;
    PtSolverSafetyLimits limits_;
    thermodynamics::Pr76RootOptions roots_;
    flash::Pr76PtFlashBackendOptions options_;
};

} // namespace mpmc::model_configuration
#endif // MPMC_MODEL_CONFIGURATION_PR76_SOLVER_SETTINGS_HPP
