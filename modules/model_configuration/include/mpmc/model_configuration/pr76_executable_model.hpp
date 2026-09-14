#ifndef MPMC_MODEL_CONFIGURATION_PR76_EXECUTABLE_MODEL_HPP
#define MPMC_MODEL_CONFIGURATION_PR76_EXECUTABLE_MODEL_HPP

#include <mpmc/model_configuration/pr76_parameters.hpp>
#include <mpmc/model_configuration/pr76_solver_settings.hpp>

#include <atomic>
#include <memory>
#include <stdexcept>

namespace mpmc::model_configuration {

// Admission failure, not a thermodynamic outcome. A caller may retry after the
// active call completes; the model does not queue work or silently share scratch.
class Pr76ModelBusyError : public std::runtime_error {
public:
    Pr76ModelBusyError() : std::runtime_error("PR76 model already has an active solve") {}
};

namespace detail {
class Pr76SolveGuard {
public:
    explicit Pr76SolveGuard(std::atomic_flag& active) : active_(active) {
        if (active_.test_and_set(std::memory_order_acquire)) {
            throw Pr76ModelBusyError();
        }
    }
    ~Pr76SolveGuard() { active_.clear(std::memory_order_release); }
    Pr76SolveGuard(const Pr76SolveGuard&) = delete;
    Pr76SolveGuard& operator=(const Pr76SolveGuard&) = delete;
private:
    std::atomic_flag& active_;
};
} // namespace detail

// Owns both immutable public snapshots and one prepared evaluator/backend graph.
// Stable address is required: backend_ borrows evaluator_. Move the unique_ptr,
// never the object. Drafts may be edited/destroyed after construction. Callers
// must keep the model alive until every call (including a rejected call) ends.
class Pr76ExecutableModel final : public flash::PtFlashBackend {
public:
    Pr76ExecutableModel(const ThermodynamicModelDefinition& definition,
                        const PtSolverSettings& settings,
                        ModelConfigurationLimits parameter_limits = {},
                        PtSolverSafetyLimits solver_limits = {},
                        ModelDataPolicy policy = ModelDataPolicy::ordinary)
        : parameter_limits_(parameter_limits),
          solver_configuration_(Pr76SolverConfiguration::create(settings, solver_limits)),
          parameter_snapshot_(prepare_parameters(definition, parameter_limits, solver_limits, policy)),
          evaluator_(thermodynamics::Pr76Phase<double>::from_parameters(
                         parameter_snapshot_.parameters(),
                         {parameter_limits.max_components, parameter_limits.max_pair_records,
                          parameter_limits.max_matrix_entries}),
                     solver_configuration_.root_options()),
          backend_(evaluator_, solver_configuration_.backend_options()) {}

    Pr76ExecutableModel(const Pr76ExecutableModel&) = delete;
    Pr76ExecutableModel(Pr76ExecutableModel&&) = delete;
    Pr76ExecutableModel& operator=(const Pr76ExecutableModel&) = delete;
    Pr76ExecutableModel& operator=(Pr76ExecutableModel&&) = delete;

    [[nodiscard]] const Pr76ModelParameters& parameter_snapshot() const & noexcept {
        return parameter_snapshot_;
    }
    const Pr76ModelParameters& parameter_snapshot() const && = delete;
    [[nodiscard]] const Pr76SolverConfiguration& solver_configuration() const & noexcept {
        return solver_configuration_;
    }
    const Pr76SolverConfiguration& solver_configuration() const && = delete;
    [[nodiscard]] const ModelConfigurationLimits& parameter_limits() const & noexcept {
        return parameter_limits_;
    }
    const ModelConfigurationLimits& parameter_limits() const && = delete;

    // Read-only snapshots/capability may be inspected while solving. Each solve
    // returns the complete, owning native envelope, with no result translation,
    // normalization, status promotion, injected hints or parameter reconstruction.
    [[nodiscard]] const flash::PtFlashBackendCapability& capability() const noexcept override {
        return backend_.capability();
    }
    [[nodiscard]] flash::PtFlashBackendResult solve(const flash::PtFlashRequest& request) override {
        const detail::Pr76SolveGuard guard(active_);
        return backend_.solve(request);
    }

private:
    static Pr76ModelParameters prepare_parameters(
        const ThermodynamicModelDefinition& definition, ModelConfigurationLimits parameter_limits,
        PtSolverSafetyLimits solver_limits, ModelDataPolicy policy) {
        // Reject the cross-layer mismatch before copying the draft/building its
        // dense matrix, even if parameter policy permits a larger component set.
        if (definition.components.size() > solver_limits.max_components) {
            throw ModelConfigurationError(ModelConfigurationErrorCode::resource_limit,
                                          "components", "model exceeds solver component ceiling");
        }
        return Pr76ModelParameters::create(definition, parameter_limits, policy);
    }

    const ModelConfigurationLimits parameter_limits_;
    const Pr76SolverConfiguration solver_configuration_;
    const Pr76ModelParameters parameter_snapshot_;
    // The evaluator itself owns a phase-model copy and mutable workspaces.
    // Reverse destruction releases the borrowing backend before its evaluator.
    flash::Pr76VleEvaluator evaluator_;
    flash::Pr76PtFlashBackend backend_;
    std::atomic_flag active_ = ATOMIC_FLAG_INIT;
};

[[nodiscard]] inline std::unique_ptr<Pr76ExecutableModel> make_pr76_executable_model(
    const ThermodynamicModelDefinition& definition, const PtSolverSettings& settings,
    ModelConfigurationLimits parameter_limits = {}, PtSolverSafetyLimits solver_limits = {},
    ModelDataPolicy policy = ModelDataPolicy::ordinary) {
    return std::make_unique<Pr76ExecutableModel>(definition, settings, parameter_limits,
                                                solver_limits, policy);
}

} // namespace mpmc::model_configuration
#endif // MPMC_MODEL_CONFIGURATION_PR76_EXECUTABLE_MODEL_HPP
