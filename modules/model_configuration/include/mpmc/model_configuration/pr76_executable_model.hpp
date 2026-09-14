#ifndef MPMC_MODEL_CONFIGURATION_PR76_EXECUTABLE_MODEL_HPP
#define MPMC_MODEL_CONFIGURATION_PR76_EXECUTABLE_MODEL_HPP

#include <mpmc/model_configuration/pr76_parameters.hpp>
#include <mpmc/model_configuration/pr76_solver_settings.hpp>

#include <atomic>
#include <cmath>
#include <memory>
#include <optional>
#include <stdexcept>
#include <string>
#include <utility>

namespace mpmc::model_configuration {

// Admission failure, not a thermodynamic outcome. A caller may retry after the
// active call completes; the model does not queue work or silently share scratch.
class Pr76ModelBusyError : public std::runtime_error {
public:
    Pr76ModelBusyError() : std::runtime_error("PR76 model already has an active solve") {}
};

// Additive location interface: the concrete exception also retains the original
// std::invalid_argument/domain_error/length_error base and its what() text.
// Only already-rejected native requests acquire this metadata.
class Pr76SolveRequestError {
public:
    virtual ~Pr76SolveRequestError() = default;
    [[nodiscard]] const std::string& field() const noexcept { return field_; }
protected:
    explicit Pr76SolveRequestError(std::string field) : field_(std::move(field)) {}
private:
    std::string field_;
};

namespace detail {
template <class Exception>
class Pr76LocatedRequestError final : public Exception, public Pr76SolveRequestError {
public:
    Pr76LocatedRequestError(const Exception& original, std::string field)
        : Exception(original), Pr76SolveRequestError(std::move(field)) {}
};

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
        try { return backend_.solve(request); }
        catch (const std::invalid_argument& error) { locate_rejection(request, error); }
        catch (const std::domain_error& error) { locate_rejection(request, error); }
        catch (const std::length_error& error) { locate_rejection(request, error); }
    }

private:
    // Diagnosis runs only after the native backend rejects. It neither admits a
    // request nor changes numerical validation/normalization. Unknown failures
    // are rethrown unchanged instead of attributing internal errors to a field.
    [[nodiscard]] std::optional<std::string> rejected_request_field(
        const flash::PtFlashRequest& request) const {
        if (request.feed.size() != evaluator_.model().size()) { return "feed"; }
        try { flash::detail::split_check_pt(request.pressure_pa, request.temperature_k); }
        catch (const std::domain_error&) {
            return !std::isfinite(request.pressure_pa) || request.pressure_pa <= 0
                ? "pressure_pa" : "temperature_k";
        }
        try { (void)flash::detail::stability_check_composition(request.feed); }
        catch (const std::domain_error&) {
            for (std::size_t i = 0; i < request.feed.size(); ++i) {
                const double value = request.feed[i];
                if (!std::isfinite(value) || value < 0 || value > 1) {
                    return "feed[" + std::to_string(i) + "]";
                }
            }
            return "feed"; // Aggregate normalization failure has no single culprit.
        }
        const auto& bounds = parameter_snapshot_.parameters().applicability();
        if (bounds.pressure_pa &&
            (static_cast<long double>(request.pressure_pa) < bounds.pressure_pa->lower ||
             static_cast<long double>(request.pressure_pa) > bounds.pressure_pa->upper)) {
            return "pressure_pa";
        }
        if (bounds.temperature_k &&
            (static_cast<long double>(request.temperature_k) < bounds.temperature_k->lower ||
             static_cast<long double>(request.temperature_k) > bounds.temperature_k->upper)) {
            return "temperature_k";
        }
        return std::nullopt;
    }
    template <class Exception>
    [[noreturn]] void locate_rejection(const flash::PtFlashRequest& request, const Exception& error) const {
        if (auto field = rejected_request_field(request)) {
            throw detail::Pr76LocatedRequestError<Exception>(error, std::move(*field));
        }
        throw; // Preserve the original exception, including its dynamic type.
    }
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
