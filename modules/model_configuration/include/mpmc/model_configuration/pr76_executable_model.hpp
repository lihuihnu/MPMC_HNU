#ifndef MPMC_MODEL_CONFIGURATION_PR76_EXECUTABLE_MODEL_HPP
#define MPMC_MODEL_CONFIGURATION_PR76_EXECUTABLE_MODEL_HPP

#include <mpmc/model_configuration/pr76_parameters.hpp>
#include <mpmc/model_configuration/pr76_solver_settings.hpp>
#include <mpmc/model_configuration/pt_solve_hints.hpp>

#include <atomic>
#include <cmath>
#include <cstddef>
#include <limits>
#include <memory>
#include <optional>
#include <span>
#include <stdexcept>
#include <string>
#include <string_view>
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

inline void pr76_check_hint_block_storage(
    const std::vector<std::vector<double>>& starts,
    std::size_t max_count, std::size_t max_entries, std::size_t max_components,
    std::string_view field) {
    if (starts.size() > max_count) {
        throw ModelConfigurationError(
            ModelConfigurationErrorCode::resource_limit, std::string(field),
            "solve hint count exceeds host ceiling");
    }
    std::size_t entries = 0U;
    for (std::size_t i = 0; i < starts.size(); ++i) {
        if (starts[i].size() > max_components) {
            throw ModelConfigurationError(
                ModelConfigurationErrorCode::resource_limit,
                std::string(field) + "[" + std::to_string(i) + "]",
                "solve hint component count exceeds host ceiling");
        }
        if (starts[i].size() > max_entries - entries) {
            throw ModelConfigurationError(
                ModelConfigurationErrorCode::resource_limit, std::string(field),
                "solve hint storage exceeds host ceiling");
        }
        entries += starts[i].size();
    }
}

inline void pr76_check_public_hint_storage(
    const PtSolveHints& hints, const PtSolverSafetyLimits& limits) {
    if (hints.version != pt_solve_hints_v1) {
        throw ModelConfigurationError(
            ModelConfigurationErrorCode::unsupported_version, "hints.version",
            "unsupported PT solve-hint version");
    }
    pr76_check_hint_block_storage(
        hints.initial_stability_starts, limits.max_stability_starts,
        limits.max_stability_start_entries, limits.max_components,
        "hints.initial_stability_starts");
    pr76_check_hint_block_storage(
        hints.final_two_phase_stability_starts, limits.max_stability_starts,
        limits.max_stability_start_entries, limits.max_components,
        "hints.final_two_phase_stability_starts");
    if (hints.three_phase_continuation_starts.size() > limits.max_three_phase_starts) {
        throw ModelConfigurationError(
            ModelConfigurationErrorCode::resource_limit,
            "hints.three_phase_continuation_starts",
            "three-phase continuation hint count exceeds host ceiling");
    }
    std::size_t entries = 0U;
    for (std::size_t h = 0; h < hints.three_phase_continuation_starts.size(); ++h) {
        const auto& hint = hints.three_phase_continuation_starts[h];
        for (std::size_t phase = 0; phase < hint.compositions.size(); ++phase) {
            const auto& composition = hint.compositions[phase];
            if (composition.size() > limits.max_components) {
                throw ModelConfigurationError(
                    ModelConfigurationErrorCode::resource_limit,
                    "hints.three_phase_continuation_starts[" + std::to_string(h) +
                        "].compositions[" + std::to_string(phase) + "]",
                    "solve hint component count exceeds host ceiling");
            }
            if (composition.size() > limits.max_three_phase_start_entries - entries) {
                throw ModelConfigurationError(
                    ModelConfigurationErrorCode::resource_limit,
                    "hints.three_phase_continuation_starts",
                    "three-phase continuation hint storage exceeds host ceiling");
            }
            entries += composition.size();
        }
    }
}

[[nodiscard]] inline flash::Pr76PtFlashBackendOptions pr76_options_with_public_hints(
    const Pr76SolverConfiguration& configuration, const PtSolveHints& hints) {
    pr76_check_public_hint_storage(hints, configuration.safety_limits());
    auto options = configuration.backend_options();
    options.initial_starts = hints.initial_stability_starts;
    options.final_starts = hints.final_two_phase_stability_starts;
    options.three_phase_starts.reserve(hints.three_phase_continuation_starts.size());
    for (const auto& public_hint : hints.three_phase_continuation_starts) {
        flash::Pr76PtThreePhaseStart native;
        native.compositions = public_hint.compositions;
        native.phase_fraction_seed = public_hint.phase_fraction_seed;
        options.three_phase_starts.push_back(std::move(native));
    }
    return options;
}

[[nodiscard]] inline std::optional<std::string> pr76_hint_composition_field(
    std::span<const double> composition, std::span<const double> feed,
    std::string prefix) {
    if (composition.size() != feed.size()) { return prefix; }
    for (std::size_t i = 0; i < composition.size(); ++i) {
        const double value = composition[i];
        if (!std::isfinite(value) || value < 0.0 || value > 1.0) {
            return prefix + "[" + std::to_string(i) + "]";
        }
    }
    try { (void)flash::detail::stability_check_composition(composition); }
    catch (const std::domain_error&) { return prefix; }
    for (std::size_t i = 0; i < composition.size(); ++i) {
        if ((feed[i] > 0.0) != (composition[i] > 0.0)) {
            return prefix + "[" + std::to_string(i) + "]";
        }
    }
    return std::nullopt;
}
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

    // Read-only snapshots/capability may be inspected while solving. The coarse
    // interface retains the original no-public-hint behavior for registry/service
    // compatibility and returns the complete owning native envelope unchanged.
    [[nodiscard]] const flash::PtFlashBackendCapability& capability() const noexcept override {
        return backend_.capability();
    }
    [[nodiscard]] flash::PtFlashBackendResult solve(const flash::PtFlashRequest& request) override {
        const detail::Pr76SolveGuard guard(active_);
        enforce_public_applicability_gap(request);
        try { return backend_.solve(request); }
        catch (const std::invalid_argument& error) { locate_rejection(request, error); }
        catch (const std::domain_error& error) { locate_rejection(request, error); }
        catch (const std::length_error& error) { locate_rejection(request, error); }
    }

    // Public per-solve initialization/continuation hints. They are copied into an
    // ephemeral backend options snapshot; the immutable model/settings snapshots
    // and the coarse PtFlashBackend path remain unchanged. Native validation and
    // numerical decisions remain authoritative for the current P/T/feed state.
    [[nodiscard]] flash::PtFlashBackendResult solve(
        const flash::PtFlashRequest& request, const PtSolveHints& hints) {
        const detail::Pr76SolveGuard guard(active_);
        auto options = detail::pr76_options_with_public_hints(solver_configuration_, hints);
        enforce_public_applicability_gap(request);
        flash::Pr76PtFlashBackend hinted_backend(evaluator_, std::move(options));
        try { return hinted_backend.solve(request); }
        catch (const std::invalid_argument& error) { locate_rejection(request, hints, error); }
        catch (const std::domain_error& error) { locate_rejection(request, hints, error); }
        catch (const std::length_error& error) { locate_rejection(request, hints, error); }
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

    // Native applicability stores only complete closed intervals. Diagnose only
    // the public information it cannot represent: one-sided violations and exact
    // equality at an exclusive endpoint. Strict violations of complete intervals
    // continue through the unchanged native property path.
    [[nodiscard]] std::optional<std::string> public_applicability_gap_field(
        const flash::PtFlashRequest& request) const {
        if (request.feed.size() != evaluator_.model().size()) { return std::nullopt; }
        try { flash::detail::split_check_pt(request.pressure_pa, request.temperature_k); }
        catch (const std::domain_error&) { return std::nullopt; }
        try { (void)flash::detail::stability_check_composition(request.feed); }
        catch (const std::domain_error&) { return std::nullopt; }

        const auto& bounds = parameter_snapshot_.definition().applicability;
        const auto gap = [](double value, const std::optional<double>& lower,
                            bool lower_exclusive, const std::optional<double>& upper,
                            bool upper_exclusive) {
            if (lower && !upper &&
                (value < *lower || (lower_exclusive && value == *lower))) {
                return true;
            }
            if (!lower && upper &&
                (value > *upper || (upper_exclusive && value == *upper))) {
                return true;
            }
            if (lower && upper &&
                ((lower_exclusive && value == *lower) ||
                 (upper_exclusive && value == *upper))) {
                return true;
            }
            return false;
        };
        if (gap(request.pressure_pa, bounds.pressure_lower_pa,
                bounds.pressure_lower_exclusive, bounds.pressure_upper_pa,
                bounds.pressure_upper_exclusive)) {
            return "pressure_pa";
        }
        if (gap(request.temperature_k, bounds.temperature_lower_k,
                bounds.temperature_lower_exclusive, bounds.temperature_upper_k,
                bounds.temperature_upper_exclusive)) {
            return "temperature_k";
        }
        return std::nullopt;
    }

    void enforce_public_applicability_gap(const flash::PtFlashRequest& request) const {
        if (auto field = public_applicability_gap_field(request)) {
            const std::domain_error error(
                "Pr76ExecutableModel: state outside declared public applicability");
            throw detail::Pr76LocatedRequestError<std::domain_error>(error, std::move(*field));
        }
    }

    [[nodiscard]] std::optional<std::string> rejected_hint_field(
        const flash::PtFlashRequest& request, const PtSolveHints& hints) const {
        const auto& options = solver_configuration_.backend_options();
        const std::size_t n = request.feed.size();
        if (n == 0U) { return std::nullopt; }

        const auto& three = hints.three_phase_continuation_starts;
        if (three.size() > options.max_three_phase_starts ||
            three.size() > std::numeric_limits<std::size_t>::max() / 3U ||
            three.size() * 3U > std::numeric_limits<std::size_t>::max() / n ||
            three.size() * 3U * n > options.max_three_phase_start_entries) {
            return "hints.three_phase_continuation_starts";
        }
        for (std::size_t h = 0; h < three.size(); ++h) {
            const auto& hint = three[h];
            const std::string prefix =
                "hints.three_phase_continuation_starts[" + std::to_string(h) + "]";
            if (!flash::detail::rr3_fractions_feasible(hint.phase_fraction_seed)) {
                return prefix + ".phase_fraction_seed";
            }
            for (std::size_t phase = 0; phase < hint.compositions.size(); ++phase) {
                if (auto field = detail::pr76_hint_composition_field(
                        hint.compositions[phase], request.feed,
                        prefix + ".compositions[" + std::to_string(phase) + "]")) {
                    return field;
                }
            }
        }

        std::size_t active = 0U;
        for (double value : request.feed) { if (value > 0.0) { ++active; } }
        const auto generated_count = [active](const flash::StabilityOptions& stability) {
            return stability.automatic_starts ? (active == 1U ? 1U : active + 2U) : 0U;
        };

        const auto& final_options = options.split.final_stability;
        const std::size_t final_generated = generated_count(final_options);
        if (!hints.final_two_phase_stability_starts.empty() &&
            final_generated <= std::numeric_limits<std::size_t>::max() - 2U) {
            const std::size_t required = final_generated + 2U;
            const bool count_caused =
                required <= final_options.max_starts &&
                hints.final_two_phase_stability_starts.size() >
                    final_options.max_starts - required;
            const std::size_t entry_capacity = final_options.max_start_entries / n;
            const bool entries_caused =
                required <= entry_capacity &&
                hints.final_two_phase_stability_starts.size() > entry_capacity - required;
            if (count_caused || entries_caused) {
                return "hints.final_two_phase_stability_starts";
            }
        }
        for (std::size_t i = 0; i < hints.final_two_phase_stability_starts.size(); ++i) {
            if (auto field = detail::pr76_hint_composition_field(
                    hints.final_two_phase_stability_starts[i], request.feed,
                    "hints.final_two_phase_stability_starts[" + std::to_string(i) + "]")) {
                return field;
            }
        }

        const auto& initial_options = options.split.initial_stability;
        const std::size_t initial_generated = generated_count(initial_options);
        if (!hints.initial_stability_starts.empty()) {
            const bool count_caused =
                initial_generated <= initial_options.max_starts &&
                hints.initial_stability_starts.size() >
                    initial_options.max_starts - initial_generated;
            const std::size_t entry_capacity = initial_options.max_start_entries / n;
            const bool entries_caused =
                initial_generated <= entry_capacity &&
                hints.initial_stability_starts.size() > entry_capacity - initial_generated;
            if (count_caused || entries_caused) {
                return "hints.initial_stability_starts";
            }
        }
        for (std::size_t i = 0; i < hints.initial_stability_starts.size(); ++i) {
            if (auto field = detail::pr76_hint_composition_field(
                    hints.initial_stability_starts[i], request.feed,
                    "hints.initial_stability_starts[" + std::to_string(i) + "]")) {
                return field;
            }
        }
        return std::nullopt;
    }

    template <class Exception>
    [[noreturn]] void locate_rejection(
        const flash::PtFlashRequest& request, const Exception& error) const {
        if (auto field = rejected_request_field(request)) {
            throw detail::Pr76LocatedRequestError<Exception>(error, std::move(*field));
        }
        throw; // Preserve the original exception, including its dynamic type.
    }

    template <class Exception>
    [[noreturn]] void locate_rejection(
        const flash::PtFlashRequest& request, const PtSolveHints& hints,
        const Exception& error) const {
        if (auto field = rejected_request_field(request)) {
            throw detail::Pr76LocatedRequestError<Exception>(error, std::move(*field));
        }
        if (auto field = rejected_hint_field(request, hints)) {
            throw detail::Pr76LocatedRequestError<Exception>(error, std::move(*field));
        }
        throw;
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

// Convert a previous accepted PR76 publication into numerical continuation hints.
// An unresolved/non-accepted point deliberately yields an empty v1 hint set,
// matching the native continuation rule that stale state is cleared. The output
// remains only a start for a future fresh solve and does not carry phase evidence.
[[nodiscard]] inline PtSolveHints make_pr76_continuation_hints(
    const flash::PtFlashBackendResult& previous) {
    PtSolveHints hints = make_pt_solve_hints_v1();
    if (previous.capability.model_profile != thermodynamics::pr76_profile) {
        throw std::invalid_argument("PR76 continuation hints: source model profile mismatch");
    }
    const auto* accepted = previous.solution.accepted_phase_set();
    if (accepted == nullptr) { return hints; }
    if (accepted->phases.empty() || accepted->phases.size() > 3U) {
        throw std::invalid_argument("PR76 continuation hints: unsupported accepted phase count");
    }
    for (const auto& phase : accepted->phases) {
        hints.initial_stability_starts.push_back(phase.composition);
        hints.final_two_phase_stability_starts.push_back(phase.composition);
    }
    if (accepted->phases.size() == 3U) {
        PtThreePhaseContinuationHint three;
        for (std::size_t phase = 0; phase < 3U; ++phase) {
            three.compositions[phase] = accepted->phases[phase].composition;
        }
        three.phase_fraction_seed = {
            accepted->phases[1].mole_phase_fraction,
            accepted->phases[2].mole_phase_fraction};
        hints.three_phase_continuation_starts.push_back(std::move(three));
    }
    return hints;
}

[[nodiscard]] inline std::unique_ptr<Pr76ExecutableModel> make_pr76_executable_model(
    const ThermodynamicModelDefinition& definition, const PtSolverSettings& settings,
    ModelConfigurationLimits parameter_limits = {}, PtSolverSafetyLimits solver_limits = {},
    ModelDataPolicy policy = ModelDataPolicy::ordinary) {
    return std::make_unique<Pr76ExecutableModel>(definition, settings, parameter_limits,
                                                solver_limits, policy);
}

} // namespace mpmc::model_configuration
#endif // MPMC_MODEL_CONFIGURATION_PR76_EXECUTABLE_MODEL_HPP
