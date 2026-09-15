#include <mpmc/model_configuration/pr76_model_registry.hpp>
#include "../executable_model/test_support.hpp"

#include <atomic>
#include <barrier>
#include <future>
#include <iostream>
#include <latch>
#include <optional>
#include <typeinfo>

namespace registry_test {
using namespace model_test;
using Code = mc::ModelRegistryErrorCode;

mc::Pr76ModelRegistryLimits limits(std::size_t count = 2) {
    mc::Pr76ModelRegistryLimits result;
    result.max_models = count;
    return result;
}
auto draft() { return definition(pr76_max3_test::model()); }
constexpr auto synthetic = mc::ModelDataPolicy::allow_synthetic_tests;

template <class F>
void fails(Code expected, F&& f) {
    try { f(); }
    catch (const mc::ModelRegistryError& e) {
        require(e.code() == expected, "wrong registry error code");
        require(std::string_view(e.what()).find("mh1_") == std::string_view::npos,
                "error disclosed a bearer handle");
        return;
    }
    throw std::runtime_error("expected registry rejection");
}
template <class F>
void busy(F&& f) {
    try { f(); }
    catch (const mc::Pr76ModelBusyError&) { return; }
    throw std::runtime_error("overlapping solve admitted");
}
void status(const mc::Pr76ModelRegistry& registry, std::size_t registered,
            std::size_t resident, bool closed = false) {
    const auto value = registry.status();
    require(value.registered_models == registered && value.resident_models == resident &&
            value.closed == closed, "registry capacity/lifetime accounting");
}

// Pause at an ownership boundary, never by sleeps or timing the EOS. Tests
// resume workers on assertion failure too, so a failure does not mask deadlock.
struct Gate {
    std::latch entered{1};
    std::latch resume{1};
    std::atomic<bool> opened{false};
    void pause() { entered.count_down(); resume.wait(); }
    void open() { if (!opened.exchange(true)) { resume.count_down(); } }
};
struct ResumeOnExit {
    Gate& gate;
    ~ResumeOnExit() { gate.open(); }
};

void creation_parity() {
    mc::Pr76ModelRegistry registry(limits()); // Real OS entropy and ordinary data policy.
    const auto native = binary_model();
    auto d = definition(native);
    auto s = preset();
    const auto handle = registry.create(d, s);
    status(registry, 1, 1);
    auto snapshot = registry.describe(handle);
    require(snapshot.definition.dataset_id == d.dataset_id && snapshot.settings == s,
            "registry did not retain complete public snapshots");
    require(snapshot.solver_limits == registry.limits().solver &&
            snapshot.parameter_limits.max_components == registry.limits().parameters.max_components,
            "host policy snapshot differs");
    compare_capability(snapshot.capability, direct(native, binary).capability);
    snapshot.definition.components.clear();
    snapshot.settings = {};
    snapshot.capability.component_ids.clear();
    d = {};
    s = {};
    const auto result = registry.solve(handle, binary);
    compare(result, direct(native, binary));
    registry.release(handle);
    status(registry, 0, 0);
    compare(result, direct(native, binary)); // Returned envelope survives release.
}
void handles() {
    // Repeated test entropy exercises deterministic stale-handle non-reuse. It
    // is not a production random source or a cryptographic-quality assertion.
    mc::Pr76ModelRegistry registry(limits(1), synthetic, [] { return mc::ModelHandleEntropy{}; });
    std::vector<std::string> retired;
    for (int i = 0; i < 16; ++i) {
        const auto handle = registry.create(draft(), preset());
        for (const auto& old : retired) {
            require(old != handle, "released handle was reused");
            fails(Code::model_not_found, [&] { (void)registry.acquire_solve(old); });
        }
        registry.release(handle);
        fails(Code::model_not_found, [&] { registry.release(handle); });
        fails(Code::model_not_found, [&] { (void)registry.describe(handle); });
        fails(Code::model_not_found, [&] { (void)registry.solve(handle, single); });
        retired.push_back(handle);
        status(registry, 0, 0);
    }
    for (const auto& invalid : std::vector<std::string>{"", "0x123", "model-1", std::string(100000, 'a'),
                                                       "mh1_" + std::string(80, 'z')}) {
        fails(Code::invalid_handle, [&] { registry.release(invalid); });
        fails(Code::invalid_handle, [&] { (void)registry.describe(invalid); });
        fails(Code::invalid_handle, [&] { (void)registry.acquire_solve(invalid); });
    }
    fails(Code::model_not_found, [&] { registry.release("mh1_" + std::string(80, 'f')); });
}
void registry_scope() {
    mc::Pr76ModelRegistry a(limits(), synthetic);
    mc::Pr76ModelRegistry b(limits(), synthetic);
    const auto ha = a.create(draft(), preset());
    const auto hb = b.create(draft(), preset());
    require(ha != hb, "independent OS-random handles collided");
    fails(Code::model_not_found, [&] { (void)b.solve(ha, single); });
    fails(Code::model_not_found, [&] { a.release(hb); });
    compare(a.solve(ha, single), b.solve(hb, single));
}
void capacity() {
    mc::Pr76ModelRegistry registry(limits(), synthetic);
    const auto a = registry.create(draft(), preset());
    const auto b = registry.create(draft(), preset());
    auto lease = registry.acquire_solve(a);
    registry.release(a);
    status(registry, 1, 2);
    fails(Code::capacity_exceeded, [&] { (void)registry.create(draft(), preset()); });
    lease = {}; // Cancelling the retained admission finally frees its model.
    status(registry, 1, 1);
    const auto c = registry.create(draft(), preset());
    require(c != a, "slot reuse revived a stale handle");
    registry.release(b);
    registry.release(c);
    status(registry, 0, 0);
    fails(Code::invalid_limits, [] { mc::Pr76ModelRegistry invalid(limits(0)); });
    fails(Code::invalid_limits, [] { mc::Pr76ModelRegistry invalid(limits(), synthetic, {}); });
}
void configuration_failure() {
    mc::Pr76ModelRegistry registry(limits(1), synthetic);
    auto d = draft();
    std::get<mc::Pr76ParameterDefinition>(d.parameters).binary.clear();
    try { (void)registry.create(d, preset()); throw std::runtime_error("missing pairs accepted"); }
    catch (const mc::ModelConfigurationError& e) {
        require(e.code() == mc::ModelConfigurationErrorCode::missing_parameter, "wrong missing pair error");
    }
    status(registry, 0, 0);
    auto s = preset();
    s.eos_root.max_iterations = 1;
    try { (void)registry.create(draft(), s); throw std::runtime_error("edited preset accepted"); }
    catch (const mc::ModelConfigurationError& e) {
        require(e.code() == mc::ModelConfigurationErrorCode::invalid_settings, "settings validation bypassed");
    }
    status(registry, 0, 0);
    const auto h = registry.create(draft(), preset());
    compare(registry.solve(h, single), direct(pr76_max3_test::model(), single));
    mc::Pr76ModelRegistry ordinary(limits(1));
    try { (void)ordinary.create(draft(), preset()); throw std::runtime_error("synthetic policy bypassed"); }
    catch (const mc::ModelConfigurationError& e) {
        require(e.code() == mc::ModelConfigurationErrorCode::invalid_source, "wrong data policy error");
    }
    status(ordinary, 0, 0);
}
void parameter_quotas() {
    std::vector<mc::Pr76ModelRegistryLimits> cases(7, limits(1));
    cases[0].parameters.max_components = 2;
    cases[1].parameters.max_pair_records = 2;
    cases[2].parameters.max_identifier_bytes = 3;
    cases[3].parameters.max_total_text_bytes = 1;
    cases[4].parameters.max_matrix_entries = 8;
    cases[5].solver.max_components = 2;
    cases[6].solver.max_root_iterations = 1;
    for (const auto& policy : cases) {
        mc::Pr76ModelRegistry registry(policy, synthetic);
        try { (void)registry.create(draft(), preset()); throw std::runtime_error("host quota bypassed"); }
        catch (const mc::ModelConfigurationError& e) {
            require(e.code() == mc::ModelConfigurationErrorCode::resource_limit, "wrong model quota error");
        }
        status(registry, 0, 0);
    }
}
void entropy_failure() {
    std::atomic<int> calls{0};
    mc::Pr76ModelRegistry registry(limits(1), synthetic, [&] {
        if (calls.fetch_add(1) == 0) {
            throw mc::ModelRegistryError(Code::entropy_unavailable, "test entropy failure");
        }
        return mc::system_model_handle_entropy();
    });
    fails(Code::entropy_unavailable, [&] { (void)registry.create(draft(), preset()); });
    status(registry, 0, 0);
    const auto h = registry.create(draft(), preset());
    compare(registry.solve(h, single), direct(pr76_max3_test::model(), single));
    require(calls.load() == 2, "unexpected entropy retries/fallback");
}
void reserved_capacity(bool close_during_create) {
    Gate gate;
    mc::Pr76ModelRegistry registry(limits(1), synthetic, [&] {
        gate.pause();
        return mc::system_model_handle_entropy();
    });
    auto task = std::async(std::launch::async, [&] { return registry.create(draft(), preset()); });
    const ResumeOnExit cleanup{gate};
    gate.entered.wait();
    status(registry, 0, 1);
    fails(Code::capacity_exceeded, [&] { (void)registry.create(draft(), preset()); });
    if (close_during_create) {
        registry.close();
        status(registry, 0, 1, true);
        gate.open();
        fails(Code::registry_closed, [&] { (void)task.get(); });
        status(registry, 0, 0, true);
    } else {
        gate.open();
        const auto h = task.get();
        status(registry, 1, 1);
        registry.release(h);
        status(registry, 0, 0);
    }
}
void release_during_solve(bool destroy_registry) {
    Gate gate;
    auto registry = std::make_unique<mc::Pr76ModelRegistry>(limits(1), synthetic);
    const auto handle = registry->create(draft(), preset());
    auto admitted = registry->acquire_solve(handle);
    // Admission precedes release; numerical work resumes AFTER invalidation.
    // This covers the difficult ownership gap without relying on solver timing.
    auto task = std::async(std::launch::async, [&, lease = std::move(admitted)]() mutable {
        gate.pause();
        return std::move(lease).solve(single);
    });
    const ResumeOnExit cleanup{gate};
    gate.entered.wait();
    registry->release(handle);
    status(*registry, 0, 1);
    fails(Code::model_not_found, [&] { (void)registry->acquire_solve(handle); });
    fails(Code::capacity_exceeded, [&] { (void)registry->create(draft(), preset()); });
    if (destroy_registry) { registry.reset(); }
    gate.open();
    compare(task.get(), direct(pr76_max3_test::model(), single));
    if (registry) {
        status(*registry, 0, 0);
        const auto next = registry->create(draft(), preset());
        compare(registry->solve(next, single), direct(pr76_max3_test::model(), single));
    }
}
void close_lifetime() {
    mc::Pr76ModelRegistry registry(limits(), synthetic);
    const auto a = registry.create(draft(), preset());
    (void)registry.create(draft(), preset());
    auto lease = registry.acquire_solve(a);
    registry.close();
    registry.close();
    status(registry, 0, 1, true);
    fails(Code::registry_closed, [&] { (void)registry.create(draft(), preset()); });
    fails(Code::registry_closed, [&] { (void)registry.describe(a); });
    fails(Code::registry_closed, [&] { (void)registry.acquire_solve(a); });
    fails(Code::registry_closed, [&] { registry.release(a); });
    compare(std::move(lease).solve(single), direct(pr76_max3_test::model(), single));
    status(registry, 0, 0, true);
}
void lease_admission() {
    mc::Pr76ModelRegistry registry(limits(), synthetic);
    const auto a = registry.create(draft(), preset());
    const auto b = registry.create(draft(), preset());
    auto first = registry.acquire_solve(a);
    auto second = registry.acquire_solve(b);
    busy([&] { (void)registry.solve(a, single); });
    busy([&] { (void)registry.acquire_solve(a); });
    second = std::move(first); // Cancel B's admission while transferring A's.
    compare(registry.solve(b, single), direct(pr76_max3_test::model(), single));
    busy([&] { (void)registry.acquire_solve(a); });
    fails(Code::invalid_lease, [&] { (void)std::move(first).solve(single); });
    require(registry.describe(a).settings == preset(), "inspection blocked by solve admission");
    compare(std::move(second).solve(single), direct(pr76_max3_test::model(), single));
    fails(Code::invalid_lease, [&] { (void)std::move(second).solve(single); });
    compare(registry.solve(a, single), direct(pr76_max3_test::model(), single));
    { auto cancelled = registry.acquire_solve(a); }
    compare(registry.solve(a, single), direct(pr76_max3_test::model(), single));
}
void solve_exception() {
    mc::Pr76ModelRegistry registry(limits(1), synthetic);
    const auto handle = registry.create(draft(), preset());
    auto lease = registry.acquire_solve(handle);
    registry.release(handle);
    const fl::PtFlashRequest invalid{0.0, single.temperature_k, single.feed};
    const auto exception = [](auto&& call) {
        try { call(); }
        catch (const std::domain_error& e) {
            const auto* located = dynamic_cast<const mc::Pr76SolveRequestError*>(&e);
            require(located != nullptr && located->field() == "pressure_pa", "lease lost request location/category");
            return std::pair{std::string(typeid(e).name()), std::string(e.what())};
        }
        throw std::runtime_error("invalid request did not throw");
    };
    // The registry must preserve the exact exception emitted by its owned model.
    // Bare-backend standard-category/message parity is checked by executable tests.
    auto standalone = model_test::create(draft(), preset());
    const auto expected = exception([&] { (void)standalone->solve(invalid); });
    const auto actual = exception([&] { (void)std::move(lease).solve(invalid); });
    require(actual == expected, "registry changed native solve exception");
    status(registry, 0, 0);
    const auto next = registry.create(draft(), preset());
    compare(registry.solve(next, single), direct(pr76_max3_test::model(), single));
}
void parallel_models() {
    mc::Pr76ModelRegistry registry(limits(), synthetic);
    const auto native_a = pr76_max3_test::model();
    const auto native_b = changed_model();
    const auto a = registry.create(definition(native_a), preset());
    const auto b = registry.create(definition(native_b), preset());
    std::barrier start(2);
    const auto run = [&](const std::string& handle, const fl::PtFlashBackendResult& expected) {
        start.arrive_and_wait();
        for (int i = 0; i < 8; ++i) { compare(registry.solve(handle, single), expected); }
    };
    auto task = std::async(std::launch::async, [&] { run(a, direct(native_a, single)); });
    run(b, direct(native_b, single));
    task.get();
    registry.release(b);
    compare(registry.solve(a, single), direct(native_a, single));
}
void release_acquire_race() {
    mc::Pr76ModelRegistry registry(limits(1), synthetic);
    for (int i = 0; i < 16; ++i) {
        const auto handle = registry.create(draft(), preset());
        std::barrier start(3);
        auto acquire = std::async(std::launch::async, [&]() -> std::optional<mc::Pr76ModelSolveLease> {
            start.arrive_and_wait();
            try { return registry.acquire_solve(handle); }
            catch (const mc::ModelRegistryError& e) {
                require(e.code() == Code::model_not_found, "invalid release/acquire race outcome");
                return std::nullopt;
            }
        });
        auto release = std::async(std::launch::async, [&] { start.arrive_and_wait(); registry.release(handle); });
        start.arrive_and_wait();
        auto lease = acquire.get();
        release.get();
        fails(Code::model_not_found, [&] { (void)registry.acquire_solve(handle); });
        status(registry, 0, lease ? 1U : 0U);
        if (lease) { compare(std::move(*lease).solve(single), direct(pr76_max3_test::model(), single)); }
        status(registry, 0, 0);
    }
}
void concurrent_creation() {
    mc::Pr76ModelRegistry registry(limits(3), synthetic);
    const auto d = draft();
    const auto s = preset();
    std::barrier start(9);
    std::vector<std::future<std::optional<std::string>>> tasks;
    for (int i = 0; i < 8; ++i) {
        tasks.push_back(std::async(std::launch::async, [&]() -> std::optional<std::string> {
            start.arrive_and_wait();
            try { return registry.create(d, s); }
            catch (const mc::ModelRegistryError& e) {
                require(e.code() == Code::capacity_exceeded, "concurrent creation failed unexpectedly");
                return std::nullopt;
            }
        }));
    }
    start.arrive_and_wait();
    std::size_t created = 0;
    for (auto& task : tasks) { if (task.get()) { ++created; } }
    require(created == 3, "concurrent creators overbooked/underfilled capacity");
    status(registry, 3, 3);
    registry.close();
    status(registry, 0, 0, true);
}
} // namespace registry_test

int main(int argc, char** argv) {
    using namespace registry_test;
    try {
        require(argc == 2, "expected one test case");
        const std::string_view name = argv[1];
        if (name == "creation_parity") { creation_parity(); }
        else if (name == "handles") { handles(); }
        else if (name == "registry_scope") { registry_scope(); }
        else if (name == "capacity") { capacity(); }
        else if (name == "configuration_failure") { configuration_failure(); }
        else if (name == "parameter_quotas") { parameter_quotas(); }
        else if (name == "entropy_failure") { entropy_failure(); }
        else if (name == "reserved_capacity") { reserved_capacity(false); }
        else if (name == "close_during_creation") { reserved_capacity(true); }
        else if (name == "release_during_solve") { release_during_solve(false); }
        else if (name == "registry_destruction") { release_during_solve(true); }
        else if (name == "close_lifetime") { close_lifetime(); }
        else if (name == "lease_admission") { lease_admission(); }
        else if (name == "solve_exception") { solve_exception(); }
        else if (name == "parallel_models") { parallel_models(); }
        else if (name == "release_acquire_race") { release_acquire_race(); }
        else if (name == "concurrent_creation") { concurrent_creation(); }
        else { throw std::runtime_error("unknown registry test case"); }
        std::cout << "PASS " << name << '\n';
        return 0;
    } catch (const std::exception& e) {
        std::cerr << "FAIL: " << e.what() << '\n';
        return 1;
    }
}
