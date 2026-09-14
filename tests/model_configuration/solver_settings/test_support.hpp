#ifndef MPMC_TEST_SOLVER_SETTINGS_SUPPORT_HPP
#define MPMC_TEST_SOLVER_SETTINGS_SUPPORT_HPP

#include <mpmc/model_configuration/pr76_solver_settings.hpp>

#include <source_location>
#include <stdexcept>
#include <string>
#include <string_view>
#include <tuple>

namespace solver_test {
namespace mc = mpmc::model_configuration;
namespace fl = mpmc::flash;
namespace th = mpmc::thermodynamics;
using Code = mc::ModelConfigurationErrorCode;

inline void require(bool condition, std::string_view message,
                    std::source_location where = std::source_location::current()) {
    if (!condition) {
        throw std::runtime_error(std::string(where.file_name()) + ":" +
                                 std::to_string(where.line()) + ": " + std::string(message));
    }
}
template <class F>
void expect_error(Code code, F&& action, std::string_view field = {}) {
    try { action(); }
    catch (const mc::ModelConfigurationError& error) {
        require(error.code() == code, error.what());
        require(!error.field().empty(), "missing error context");
        if (!field.empty()) { require(error.field() == field, error.field()); }
        return;
    }
    throw std::runtime_error("invalid solver configuration accepted");
}
inline mc::PtSolverSettings custom() {
    auto settings = mc::resolve_pt_solver_preset(mc::mpmc_balanced_default_v1);
    settings.kind = mc::PtSolverSettingsKind::custom;
    settings.preset_id.clear();
    return settings;
}
inline auto stability_values(const fl::StabilityOptions& o) {
    return std::tuple{o.tpd_tolerance, o.stationarity_tolerance, o.armijo, o.log_step_limit,
                      o.max_iterations, o.max_backtracks, o.max_evaluations, o.max_components,
                      o.max_starts, o.max_start_entries, o.automatic_starts};
}
inline auto two_phase_values(const fl::PtSplitIterationOptions& o) {
    return std::tuple{o.fugacity_tolerance, o.mass_absolute_tolerance, o.mass_relative_tolerance,
                      o.minimum_phase_fraction, o.log_k_separation, o.relative_z_separation,
                      o.max_log_step, o.residual_decrease, o.max_iterations, o.max_backtracks,
                      o.max_evaluations, o.rr.max_iterations, o.rr.max_components, o.gibbs_armijo};
}
inline auto three_phase_values(const fl::PtThreePhaseOptions& o) {
    return std::tuple{o.chemical_potential_tolerance, o.mass_absolute_tolerance, o.mass_relative_tolerance,
                      o.balance_tolerance, o.minimum_phase_fraction, o.log_composition_separation,
                      o.max_log_step, o.residual_decrease, o.max_iterations, o.max_backtracks,
                      o.max_balance_iterations, o.max_balance_backtracks, o.max_evaluations, o.max_components};
}
inline void require_options_equal(const fl::Pr76PtFlashBackendOptions& a,
                                  const fl::Pr76PtFlashBackendOptions& b) {
    require(stability_values(a.split.initial_stability) == stability_values(b.split.initial_stability), "initial stability mapping");
    require(stability_values(a.split.final_stability) == stability_values(b.split.final_stability), "final two-phase stability mapping");
    require(stability_values(a.final_three_phase_stability) == stability_values(b.final_three_phase_stability), "final three-phase stability mapping");
    require(two_phase_values(a.split.iteration) == two_phase_values(b.split.iteration), "two-phase mapping");
    require(three_phase_values(a.three_phase) == three_phase_values(b.three_phase), "three-phase mapping");
    require(a.split.max_split_attempts == b.split.max_split_attempts &&
            a.max_three_phase_attempts == b.max_three_phase_attempts &&
            a.new_phase_seed_fraction == b.new_phase_seed_fraction, "attempt/seed mapping");
    require(a.max_three_phase_starts == b.max_three_phase_starts &&
            a.max_three_phase_start_entries == b.max_three_phase_start_entries, "host quota mapping");
    require(a.initial_starts == b.initial_starts && a.final_starts == b.final_starts &&
            a.three_phase_starts.empty() && b.three_phase_starts.empty(), "numeric mapping invented hints");
}
void solve_parity(std::string_view name);
} // namespace solver_test
#endif
