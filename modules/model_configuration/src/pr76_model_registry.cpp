#include <mpmc/model_configuration/pr76_model_registry.hpp>

#include <atomic>
#include <cstdint>
#include <limits>
#include <map>
#include <mutex>
#include <utility>

namespace mpmc::model_configuration {
namespace detail {
struct Pr76RegistryCapacity {
    std::atomic<std::size_t> resident{0};
};
// Separate from the registry map to avoid an ownership cycle. The reservation
// is the FIRST entry member, so its destructor runs AFTER the model is freed.
struct Pr76RegistryReservation {
    explicit Pr76RegistryReservation(std::shared_ptr<Pr76RegistryCapacity> capacity) noexcept
        : capacity_(std::move(capacity)) {
        capacity_->resident.fetch_add(1, std::memory_order_relaxed);
    }
    Pr76RegistryReservation(Pr76RegistryReservation&&) noexcept = default;
    Pr76RegistryReservation(const Pr76RegistryReservation&) = delete;
    Pr76RegistryReservation& operator=(const Pr76RegistryReservation&) = delete;
    Pr76RegistryReservation& operator=(Pr76RegistryReservation&&) = delete;
    ~Pr76RegistryReservation() {
        if (capacity_) { capacity_->resident.fetch_sub(1, std::memory_order_release); }
    }
    std::shared_ptr<Pr76RegistryCapacity> capacity_;
};
struct Pr76RegistryEntry {
    Pr76RegistryEntry(Pr76RegistryReservation reservation,
        const ThermodynamicModelDefinition& definition, const PtSolverSettings& settings,
        const Pr76ModelRegistryLimits& limits, ModelDataPolicy policy)
        : reservation_(std::move(reservation)),
          model_(definition, settings, limits.parameters, limits.solver, policy) {}
    Pr76RegistryReservation reservation_;
    Pr76ExecutableModel model_;
    std::atomic_flag admitted_ = ATOMIC_FLAG_INIT;
};
struct Pr76RegistryState {
    Pr76RegistryState(Pr76ModelRegistryLimits limits, ModelDataPolicy policy,
                      ModelHandleEntropySource entropy)
        : limits_(limits), policy_(policy), entropy_(std::move(entropy)) {
        if (limits_.max_models == 0 || !entropy_) {
            throw ModelRegistryError(ModelRegistryErrorCode::invalid_limits,
                                     "model registry: positive capacity and entropy source required");
        }
    }
    const Pr76ModelRegistryLimits limits_;
    const ModelDataPolicy policy_;
    const ModelHandleEntropySource entropy_;
    const std::shared_ptr<Pr76RegistryCapacity> capacity_ = std::make_shared<Pr76RegistryCapacity>();
    std::mutex mutex_;
    std::map<std::string, std::shared_ptr<Pr76RegistryEntry>, std::less<>> entries_;
    std::uint64_t next_sequence_{1}; // Zero permanently denotes exhaustion.
    bool closed_{false};
};
} // namespace detail
namespace {
void require_open(const detail::Pr76RegistryState& state) {
    if (state.closed_) {
        throw ModelRegistryError(ModelRegistryErrorCode::registry_closed, "model registry is closed");
    }
}
void check_handle(std::string_view handle) {
    // Fixed bound before any allocation/lookup. Generic version prefix does not
    // encode the EOS, memory address, component set or parameter identity.
    if (handle.size() != 84 || !handle.starts_with("mh1_")) {
        throw ModelRegistryError(ModelRegistryErrorCode::invalid_handle, "invalid model handle");
    }
    for (char c : handle.substr(4)) {
        if (!((c >= '0' && c <= '9') || (c >= 'a' && c <= 'f'))) {
            throw ModelRegistryError(ModelRegistryErrorCode::invalid_handle, "invalid model handle");
        }
    }
}
auto find_entry(detail::Pr76RegistryState& state, std::string_view handle) {
    require_open(state);
    check_handle(handle);
    const auto found = state.entries_.find(handle);
    if (found == state.entries_.end()) {
        // No echo of bearer tokens and no unbounded released-handle tombstones.
        throw ModelRegistryError(ModelRegistryErrorCode::model_not_found, "model handle not found");
    }
    return found;
}
std::string make_handle(const ModelHandleEntropy& entropy, std::uint64_t sequence) {
    constexpr char hex[] = "0123456789abcdef";
    std::string handle = "mh1_";
    handle.reserve(84);
    for (unsigned char byte : entropy) {
        handle.push_back(hex[byte >> 4]);
        handle.push_back(hex[byte & 15]);
    }
    // Fresh entropy prevents prediction even after observing another handle.
    // The monotonic suffix guarantees no reuse in this registry even if a host
    // entropy source repeats. Wrap is forbidden, so no tombstone cache is needed.
    for (unsigned int shift = 64; shift != 0; shift -= 4) {
        handle.push_back(hex[(sequence >> (shift - 4)) & 15U]);
    }
    return handle;
}
} // namespace

Pr76ModelSolveLease::Pr76ModelSolveLease(std::shared_ptr<detail::Pr76RegistryEntry> entry) noexcept
    : entry_(std::move(entry)) {}
Pr76ModelSolveLease::~Pr76ModelSolveLease() { reset(); }
Pr76ModelSolveLease::Pr76ModelSolveLease(Pr76ModelSolveLease&& other) noexcept = default;
Pr76ModelSolveLease& Pr76ModelSolveLease::operator=(Pr76ModelSolveLease&& other) noexcept {
    if (this != &other) {
        reset();
        entry_ = std::move(other.entry_);
    }
    return *this;
}
void Pr76ModelSolveLease::reset() noexcept {
    if (entry_) {
        entry_->admitted_.clear(std::memory_order_release);
        entry_.reset();
    }
}
flash::PtFlashBackendResult Pr76ModelSolveLease::solve(const flash::PtFlashRequest& request) && {
    // Consume before entering the solver; return and exception both destroy the
    // local lease after the solver stops using its workspace. No registry lock.
    Pr76ModelSolveLease admitted(std::move(*this));
    if (!admitted.entry_) {
        throw ModelRegistryError(ModelRegistryErrorCode::invalid_lease, "model solve lease is empty or consumed");
    }
    return admitted.entry_->model_.solve(request);
}
flash::PtFlashBackendResult Pr76ModelSolveLease::solve(
    const flash::PtFlashRequest& request, const PtSolveHints& hints) && {
    // Hints are per-call data owned by the caller. The retained registry entry
    // pins only the model/workspace; release/close cannot turn hints into state.
    Pr76ModelSolveLease admitted(std::move(*this));
    if (!admitted.entry_) {
        throw ModelRegistryError(ModelRegistryErrorCode::invalid_lease, "model solve lease is empty or consumed");
    }
    return admitted.entry_->model_.solve(request, hints);
}

Pr76ModelRegistry::Pr76ModelRegistry(Pr76ModelRegistryLimits limits, ModelDataPolicy policy,
                                   ModelHandleEntropySource entropy)
    : state_(std::make_shared<detail::Pr76RegistryState>(limits, policy, std::move(entropy))) {}
Pr76ModelRegistry::~Pr76ModelRegistry() = default;
const Pr76ModelRegistryLimits& Pr76ModelRegistry::limits() const & noexcept { return state_->limits_; }

std::string Pr76ModelRegistry::create(const ThermodynamicModelDefinition& definition,
                                    const PtSolverSettings& settings) {
    const auto state = state_;
    std::uint64_t sequence = 0;
    // Reserve BEFORE entropy, validation, copies, dense matrices and model
    // construction. Concurrent creators cannot overbook an unregistered slot.
    auto reservation = [&] {
        const std::lock_guard lock(state->mutex_);
        require_open(*state);
        if (state->capacity_->resident.load(std::memory_order_acquire) >= state->limits_.max_models) {
            throw ModelRegistryError(ModelRegistryErrorCode::capacity_exceeded, "model registry capacity exceeded");
        }
        if (state->next_sequence_ == 0) {
            throw ModelRegistryError(ModelRegistryErrorCode::handle_sequence_exhausted,
                                     "model registry handle sequence exhausted");
        }
        sequence = state->next_sequence_;
        state->next_sequence_ = sequence == std::numeric_limits<std::uint64_t>::max() ? 0 : sequence + 1;
        return detail::Pr76RegistryReservation(state->capacity_);
    }();
    auto handle = make_handle(state->entropy_(), sequence);
    auto entry = std::make_shared<detail::Pr76RegistryEntry>(
        std::move(reservation), definition, settings, state->limits_, state->policy_);
    {
        const std::lock_guard lock(state->mutex_);
        require_open(*state); // A concurrent close must prevent publication.
        state->entries_.emplace(handle, std::move(entry));
    }
    return handle;
}
Pr76RegisteredModelSnapshot Pr76ModelRegistry::describe(std::string_view handle) const {
    const auto state = state_;
    auto entry = [&] {
        const std::lock_guard lock(state->mutex_);
        return find_entry(*state, handle)->second;
    }();
    // Owning entry also covers copying a description concurrently with release.
    // Copy outside the mutex; no borrowed model pointer escapes to the caller.
    const auto& model = entry->model_;
    return {model.parameter_snapshot().definition(), model.solver_configuration().settings(),
            model.parameter_limits(), model.solver_configuration().safety_limits(), model.capability()};
}
Pr76ModelSolveLease Pr76ModelRegistry::acquire_solve(std::string_view handle) {
    const auto state = state_;
    const std::lock_guard lock(state->mutex_);
    auto entry = find_entry(*state, handle)->second;
    if (entry->admitted_.test_and_set(std::memory_order_acquire)) { throw Pr76ModelBusyError(); }
    return Pr76ModelSolveLease(std::move(entry));
}
flash::PtFlashBackendResult Pr76ModelRegistry::solve(std::string_view handle,
                                                   const flash::PtFlashRequest& request) {
    return acquire_solve(handle).solve(request);
}
flash::PtFlashBackendResult Pr76ModelRegistry::solve(
    std::string_view handle, const flash::PtFlashRequest& request,
    const PtSolveHints& hints) {
    return acquire_solve(handle).solve(request, hints);
}
void Pr76ModelRegistry::release(std::string_view handle) {
    const auto state = state_;
    std::shared_ptr<detail::Pr76RegistryEntry> retired;
    {
        const std::lock_guard lock(state->mutex_);
        const auto found = find_entry(*state, handle);
        retired = std::move(found->second);
        state->entries_.erase(found); // Linearization point for invalidation.
    }
    // Retire outside the mutex. Any admitted lease still pins the model/slot.
}
void Pr76ModelRegistry::close() {
    const auto state = state_;
    decltype(state->entries_) retired;
    {
        const std::lock_guard lock(state->mutex_);
        state->closed_ = true;
        retired.swap(state->entries_);
    }
}
Pr76ModelRegistryStatus Pr76ModelRegistry::status() const {
    const auto state = state_;
    const std::lock_guard lock(state->mutex_);
    return {state->entries_.size(), state->capacity_->resident.load(std::memory_order_acquire), state->closed_};
}
} // namespace mpmc::model_configuration
