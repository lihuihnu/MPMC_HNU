#ifndef MPMC_MODEL_CONFIGURATION_PR76_MODEL_REGISTRY_HPP
#define MPMC_MODEL_CONFIGURATION_PR76_MODEL_REGISTRY_HPP

#include <mpmc/model_configuration/pr76_executable_model.hpp>

#include <array>
#include <functional>
#include <memory>
#include <string>
#include <string_view>

namespace mpmc::model_configuration {

struct Pr76ModelRegistryLimits {
    // Counts reservations + registered models + released models retained by a
    // solve/inspection. A released-but-running model does NOT free a slot yet.
    std::size_t max_models{64};
    ModelConfigurationLimits parameters;
    PtSolverSafetyLimits solver;
};

enum class ModelRegistryErrorCode {
    invalid_limits, capacity_exceeded, invalid_handle, model_not_found,
    registry_closed, invalid_lease, entropy_unavailable, handle_sequence_exhausted
};
class ModelRegistryError : public std::runtime_error {
public:
    ModelRegistryError(ModelRegistryErrorCode code, const char* message)
        : std::runtime_error(message), code_(code) {}
    [[nodiscard]] ModelRegistryErrorCode code() const noexcept { return code_; }
private:
    ModelRegistryErrorCode code_;
};

using ModelHandleEntropy = std::array<unsigned char, 32>;
// Host composition only, never a field accepted from an untrusted draft. A
// supplied source must provide fresh cryptographic entropy and be thread-safe.
// Deterministic sources are suitable ONLY for software tests. Failure propagates.
using ModelHandleEntropySource = std::function<ModelHandleEntropy()>;
[[nodiscard]] ModelHandleEntropy system_model_handle_entropy();

struct Pr76ModelRegistryStatus {
    std::size_t registered_models{};
    std::size_t resident_models{};
    bool closed{false};
};
// Owning description: mutations to it cannot affect registry execution. It is
// not a mutable model pointer, nor a handle that admits future work.
struct Pr76RegisteredModelSnapshot {
    ThermodynamicModelDefinition definition;
    PtSolverSettings settings;
    ModelConfigurationLimits parameter_limits;
    PtSolverSafetyLimits solver_limits;
    flash::PtFlashBackendCapability capability;
};

namespace detail {
struct Pr76RegistryEntry;
struct Pr76RegistryState;
} // namespace detail

// A single already-admitted solve. Acquisition pins ownership and excludes
// another solve on the same model. Move to solve exactly once; destruction of
// an unused lease cancels that admission. It can outlive release/close/registry.
class Pr76ModelSolveLease {
public:
    Pr76ModelSolveLease() noexcept = default;
    ~Pr76ModelSolveLease();
    Pr76ModelSolveLease(Pr76ModelSolveLease&& other) noexcept;
    Pr76ModelSolveLease& operator=(Pr76ModelSolveLease&& other) noexcept;
    Pr76ModelSolveLease(const Pr76ModelSolveLease&) = delete;
    Pr76ModelSolveLease& operator=(const Pr76ModelSolveLease&) = delete;
    [[nodiscard]] flash::PtFlashBackendResult solve(const flash::PtFlashRequest& request) &&;
    [[nodiscard]] flash::PtFlashBackendResult solve(
        const flash::PtFlashRequest& request, const PtSolveHints& hints) &&;
private:
    friend class Pr76ModelRegistry;
    explicit Pr76ModelSolveLease(std::shared_ptr<detail::Pr76RegistryEntry> entry) noexcept;
    void reset() noexcept;
    std::shared_ptr<detail::Pr76RegistryEntry> entry_;
};

// Session owner, with no global cache/TTL/background thread. Public operations
// may run concurrently; use close() to race with operations, not destruction of
// the C++ registry object. Detached solve leases remain valid after destruction.
class Pr76ModelRegistry {
public:
    explicit Pr76ModelRegistry(Pr76ModelRegistryLimits limits = {},
        ModelDataPolicy policy = ModelDataPolicy::ordinary,
        ModelHandleEntropySource entropy = system_model_handle_entropy);
    ~Pr76ModelRegistry();
    Pr76ModelRegistry(const Pr76ModelRegistry&) = delete;
    Pr76ModelRegistry(Pr76ModelRegistry&&) = delete;
    Pr76ModelRegistry& operator=(const Pr76ModelRegistry&) = delete;
    Pr76ModelRegistry& operator=(Pr76ModelRegistry&&) = delete;

    [[nodiscard]] std::string create(const ThermodynamicModelDefinition& definition,
                                     const PtSolverSettings& settings);
    [[nodiscard]] Pr76RegisteredModelSnapshot describe(std::string_view handle) const;
    [[nodiscard]] Pr76ModelSolveLease acquire_solve(std::string_view handle);
    [[nodiscard]] flash::PtFlashBackendResult solve(std::string_view handle,
                                                   const flash::PtFlashRequest& request);
    [[nodiscard]] flash::PtFlashBackendResult solve(
        std::string_view handle, const flash::PtFlashRequest& request,
        const PtSolveHints& hints);
    // Removes future lookup/admission immediately; duplicate/unknown/released
    // handles all fail with model_not_found. Already admitted work may finish.
    void release(std::string_view handle);
    // Idempotent. Invalidate every handle, reject later operations and prevent
    // pending creates from publishing; already admitted leases may finish.
    void close();
    [[nodiscard]] Pr76ModelRegistryStatus status() const;
    [[nodiscard]] const Pr76ModelRegistryLimits& limits() const & noexcept;
    const Pr76ModelRegistryLimits& limits() const && = delete;
private:
    const std::shared_ptr<detail::Pr76RegistryState> state_;
};

} // namespace mpmc::model_configuration
#endif // MPMC_MODEL_CONFIGURATION_PR76_MODEL_REGISTRY_HPP
