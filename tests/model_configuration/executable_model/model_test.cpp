#include "test_support.hpp"

#include <barrier>
#include <future>
#include <iostream>
#include <limits>
#include <typeinfo>

namespace model_test {
namespace {
void phase_parity(std::string_view name) {
    const bool two = name == "two_phase";
    const auto native = two ? binary_model() : pr76_max3_test::model();
    const auto& request = two ? binary : (name == "single_phase" ? single : ternary);
    auto model = create(definition(native), preset());
    // Dispatch also works through the existing coarse, model-neutral interface.
    fl::PtFlashBackend& backend = *model;
    const auto actual = backend.solve(request);
    compare(actual, direct(native, request));
    compare_capability(backend.capability(), actual.capability);
    require(!actual.solution.global_stability_proven && !actual.morphology_resolved,
            "factory strengthened scientific claims");
    if (name != "ternary_cold_start") {
        require(actual.solution.accepted_phase_count() == (two ? 2U : 1U), "expected topology did not close");
    }
    // No public hints exist yet: cold ternary search must preserve its native
    // outcome, including indeterminacy, rather than inject the known 3-phase seed.
    std::cout << "accepted phases=" << actual.solution.accepted_phase_count()
              << "; status=" << static_cast<int>(actual.solution.status)
              << "; evidence=" << actual.transition_report.evidence.size()
              << "; diagnostic=" << actual.solution.diagnostic << '\n';
}
void parameter_isolation() {
    const auto native_a = pr76_max3_test::model();
    const auto native_b = changed_model();
    auto draft = definition(native_a);
    auto a = create(draft, preset());
    const auto before = a->solve(single);
    draft = definition(native_b);
    auto b = create(draft, preset());
    draft.components.clear();
    draft.parameters = std::monostate{};
    const auto other = b->solve(single);
    compare(other, direct(native_b, single));
    compare(a->solve(single), before);
    require(a->parameter_snapshot().parameters().kij(0, 1) == 0.05 &&
            b->parameter_snapshot().parameters().kij(0, 1) == 0.15, "kij snapshots alias");
    require(before.solution.accepted_phase_count() == 1 && other.solution.accepted_phase_count() == 1,
            "isolation anchor did not close");
    require(before.solution.candidate_phase_set->phases[0].activity.ln_phi !=
            other.solution.candidate_phase_set->phases[0].activity.ln_phi,
            "different kij did not reach evaluator");
    compare(a->solve(ternary), direct(native_a, ternary));
    compare(b->solve(ternary), direct(native_b, ternary));
    b.reset();
    compare(a->solve(single), before);
    require(a->parameter_snapshot().definition().revision == "v1", "public source snapshot mutated");
}
void settings_isolation(std::string_view name) {
    const auto native = pr76_max3_test::model();
    auto draft = preset();
    auto a = create(definition(native), draft);
    draft = custom();
    th::Pr76RootOptions roots;
    fl::Pr76PtFlashBackendOptions options;
    if (name == "root_budget") {
        draft.eos_root.max_iterations = 1;
        roots.max_iterations = 1;
    } else {
        draft.initial_stability.max_property_evaluations = 1;
        options.split.initial_stability.max_evaluations = 1;
    }
    auto b = create(definition(native), draft);
    const auto snapshot = draft;
    draft = {};
    const auto before = a->solve(single);
    const auto limited = b->solve(single);
    compare(limited, direct(native, single, roots, options));
    require(limited.solution.status == fl::PtPhaseSetStatus::indeterminate &&
            limited.solution.accepted_phase_count() == 0 && !limited.solution.diagnostic.empty(),
            "budget exhaustion lost");
    compare(a->solve(single), before);
    compare(b->solve(single), limited);
    require(a->solver_configuration().settings() == preset() &&
            b->solver_configuration().settings() == snapshot, "settings snapshots alias");
}
void lifetime() {
    const auto expected = direct(pr76_max3_test::model(), single);
    const auto retained = [] {
        auto d = definition(pr76_max3_test::model());
        auto s = preset();
        auto model = create(d, s);
        const auto* original_address = model.get();
        auto moved_owner = std::move(model);
        require(!model && moved_owner.get() == original_address, "owner move changed graph address");
        d = {};
        s = {};
        return moved_owner->solve(single);
    }();
    compare(retained, expected); // All drafts, graph and workspaces destroyed.
    auto edited_result = retained;
    edited_result.capability.component_ids[0] = "unrelated";
    edited_result.solution.feed[0] = 0.0;
    compare(retained, expected);
}
void creation_errors() {
    auto d = definition(pr76_max3_test::model());
    const auto good = d;
    const auto fails = [&](mc::ModelConfigurationErrorCode code) {
        try { auto unused = create(d, preset()); }
        catch (const mc::ModelConfigurationError& e) {
            require(e.code() == code && !e.field().empty(), "wrong construction error");
            return;
        }
        throw std::runtime_error("invalid construction accepted");
    };
    d.family = mc::ThermodynamicModelFamily::soreide_whitson_1992;
    fails(mc::ModelConfigurationErrorCode::unsupported_family);
    d.family = mc::ThermodynamicModelFamily::cubic_plus_association;
    fails(mc::ModelConfigurationErrorCode::unsupported_family);
    d = good;
    std::get<mc::Pr76ParameterDefinition>(d.parameters).binary.pop_back();
    fails(mc::ModelConfigurationErrorCode::missing_parameter);
    try { auto unused = mc::make_pr76_executable_model(good, preset()); }
    catch (const mc::ModelConfigurationError& e) {
        require(e.code() == mc::ModelConfigurationErrorCode::invalid_source, "synthetic policy bypassed");
        d = good;
        auto invalid = preset();
        invalid.eos_root.max_iterations = 1; // Edited preset must fail.
        try { auto unused = create(d, invalid); }
        catch (const mc::ModelConfigurationError& settings_error) {
            require(settings_error.code() == mc::ModelConfigurationErrorCode::invalid_settings,
                    "factory bypassed settings validation");
            compare(create(d, preset())->solve(single), direct(pr76_max3_test::model(), single));
            return;
        }
        throw std::runtime_error("edited preset accepted");
    }
    throw std::runtime_error("synthetic fixture accepted by ordinary policy");
}
void resource_limits() {
    const auto d = definition(pr76_max3_test::model());
    const auto check = [&](mc::ModelConfigurationLimits p, mc::PtSolverSafetyLimits s) {
        try { auto unused = mc::make_pr76_executable_model(d, preset(), p, s,
                    mc::ModelDataPolicy::allow_synthetic_tests); }
        catch (const mc::ModelConfigurationError& e) {
            require(e.code() == mc::ModelConfigurationErrorCode::resource_limit, "wrong quota error");
            return;
        }
        throw std::runtime_error("model quota mismatch accepted");
    };
    mc::PtSolverSafetyLimits s;
    s.max_components = 2;
    check({}, s);
    mc::ModelConfigurationLimits p;
    p.max_components = 2;
    check(p, {});
    p = {};
    p.max_matrix_entries = 8;
    check(p, {});
    p = {};
    p.max_pair_records = 2;
    check(p, {});
    p.max_pair_records = 3;
    p.max_components = 3;
    p.max_matrix_entries = 9;
    s.max_components = 3;
    auto model = mc::make_pr76_executable_model(d, preset(), p, s,
                                               mc::ModelDataPolicy::allow_synthetic_tests);
    require(model->parameter_limits().max_matrix_entries == 9 &&
            model->solver_configuration().safety_limits() == s, "host limits not retained");
    compare(model->solve(single), direct(pr76_max3_test::model(), single));
}
// Preserve native exception type/message OR its returned failure envelope. The
// underlying solver determines which bad-domain cases throw and which diagnose.
template <class F>
std::pair<std::string, std::string> error(F&& f) {
    try { f(); }
    catch (const std::exception& e) { return {typeid(e).name(), e.what()}; }
    throw std::runtime_error("invalid request unexpectedly returned");
}
void request_errors() {
    const auto native = pr76_max3_test::model();
    auto model = create(definition(native), preset());
    const std::vector<fl::PtFlashRequest> invalid{
        {0.0, 250.0, single.feed}, {1e6, -1.0, single.feed},
        {1e6, 250.0, {1.0, 0.0}}, {1e6, 250.0, {-0.1, 0.5, 0.6}},
        {1e6, 250.0, {0.5, 0.5, 0.5}},
        {std::numeric_limits<double>::quiet_NaN(), 250.0, single.feed}};
    for (const auto& r : invalid) {
        const auto expected = error([&] { (void)direct(native, r); });
        const auto actual = error([&] { (void)model->solve(r); });
        require(actual == expected, "request exception changed");
        compare(model->solve(single), direct(native, single)); // Admission released.
    }
}
void applicability() {
    const auto native = changed_model(true);
    auto model = create(definition(native), preset());
    compare(model->solve(single), direct(native, single));
    for (const auto& r : std::vector<fl::PtFlashRequest>{{2e6, 250.0, single.feed},
                                                       {1e6, 300.0, single.feed}}) {
        // Property-domain failures become diagnostic results in the native search.
        const auto expected = direct(native, r);
        const auto actual = model->solve(r);
        compare(actual, expected);
        require(actual.solution.accepted_phase_count() == 0, "declared bounds bypassed");
    }
    compare(model->solve(single), direct(native, single));
}
void admission() {
    // Deterministic overlap/exception checks on the exact production guard; no
    // timing assumptions, public test hooks or long-running blocking evaluator.
    std::atomic_flag active = ATOMIC_FLAG_INIT;
    try {
        const mc::detail::Pr76SolveGuard first(active);
        for (int i = 0; i < 2; ++i) {
            bool busy = false;
            try { const mc::detail::Pr76SolveGuard second(active); }
            catch (const mc::Pr76ModelBusyError&) { busy = true; }
            require(busy, "overlap admitted or rejected guard cleared active ownership");
        }
        throw std::runtime_error("test unwind");
    } catch (const std::runtime_error& e) {
        require(std::string_view(e.what()) == "test unwind", "admission regression");
    }
    const mc::detail::Pr76SolveGuard recovered(active);
}
void parallel_models() {
    const auto native_a = pr76_max3_test::model();
    const auto native_b = changed_model();
    auto a = create(definition(native_a), preset());
    auto b = create(definition(native_b), preset());
    const auto expected_a = direct(native_a, single);
    const auto expected_b = direct(native_b, single);
    std::barrier start(2);
    const auto run = [&](mc::Pr76ExecutableModel& model, const fl::PtFlashBackendResult& expected) {
        start.arrive_and_wait();
        for (int i = 0; i < 8; ++i) { compare(model.solve(single), expected); }
    };
    auto task = std::async(std::launch::async, [&] { run(*a, expected_a); });
    run(*b, expected_b);
    task.get();
    compare(a->solve(single), expected_a);
    compare(b->solve(single), expected_b);
}
} // namespace
} // namespace model_test

int main(int argc, char** argv) {
    using namespace model_test;
    try {
        require(argc == 2, "expected one test case");
        const std::string_view name = argv[1];
        if (name == "single_phase" || name == "two_phase" || name == "ternary_cold_start") { phase_parity(name); }
        else if (name == "parameter_isolation") { parameter_isolation(); }
        else if (name == "root_budget" || name == "settings_isolation") { settings_isolation(name); }
        else if (name == "lifetime") { lifetime(); }
        else if (name == "creation_errors") { creation_errors(); }
        else if (name == "resource_limits") { resource_limits(); }
        else if (name == "request_errors") { request_errors(); }
        else if (name == "applicability") { applicability(); }
        else if (name == "admission") { admission(); }
        else if (name == "parallel_models") { parallel_models(); }
        else { throw std::runtime_error("unknown test case"); }
        std::cout << "PASS " << name << '\n';
        return 0;
    } catch (const std::exception& e) {
        std::cerr << "FAIL: " << e.what() << '\n';
        return 1;
    }
}
