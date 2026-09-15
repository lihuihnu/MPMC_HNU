#include "test_support.hpp"

#include <iostream>
#include <limits>
#include <type_traits>

bool solver_public_header_compiles();

namespace solver_test {
namespace {
void preset() {
    const auto s = mc::resolve_pt_solver_preset("mpmc-balanced-default/v1");
    require(s.version == "pt-solver-settings/v1" && s.kind == mc::PtSolverSettingsKind::preset,
            "preset identity changed");
    const auto prepared = mc::Pr76SolverConfiguration::create(s);
    require(prepared.settings() == s, "resolved snapshot changed");
    require(prepared.root_options().max_iterations == 2048, "frozen root default");
    // Independently frozen audited constants; do not derive these from the public resolver.
    const fl::StabilityOptions stability{1e-10, 1e-8, 1e-4, 4.0, 512, 32, 100000, 256, 1024, 262144, true};
    const fl::PtSplitIterationOptions split{1e-11, 1e-12, 1e-10, 1e-10, 1e-7, 1e-8,
                                          2.0, 1e-4, 512, 24, 20000, {192, 256}, 1e-4};
    const fl::PtThreePhaseOptions three{1e-11, 1e-12, 1e-10, 2e-13, 1e-10, 1e-7,
                                      2.0, 1e-4, 512, 32, 192, 48, 30000, 256};
    const fl::Pr76PtFlashBackendOptions frozen{
        {split, stability, stability, 16}, three, stability, 16, 0.1, {}, 16, 12288, {}, {}};
    require_options_equal(prepared.backend_options(), frozen);
    // A drift in upstream defaults is a review signal, not permission to mutate preset v1.
    require_options_equal(prepared.backend_options(), fl::Pr76PtFlashBackendOptions{});
    expect_error(Code::unsupported_preset, [] { (void)mc::resolve_pt_solver_preset("latest"); });
    auto bad = s;
    bad.version = "pt-solver-settings/v2";
    expect_error(Code::unsupported_version, [&] { (void)mc::Pr76SolverConfiguration::create(bad); });
    bad = s;
    bad.initial_stability.max_iterations = 1;
    expect_error(Code::invalid_settings, [&] { (void)mc::Pr76SolverConfiguration::create(bad); }, "preset_id");
    bad = custom();
    bad.preset_id = "mpmc-balanced-default/v1";
    expect_error(Code::invalid_settings, [&] { (void)mc::Pr76SolverConfiguration::create(bad); });
    bad = s;
    bad.kind = static_cast<mc::PtSolverSettingsKind>(99);
    expect_error(Code::invalid_settings, [&] { (void)mc::Pr76SolverConfiguration::create(bad); });
}
void custom_mapping() {
    auto s = custom();
    s.eos_root.max_iterations = 37;
    s.initial_stability = {2e-10, 3e-8, 0.002, 3.0, 17, 19, 23, false, 29};
    s.final_two_phase_stability = {4e-10, 5e-8, 0.003, 5.0, 31, 37, 41, true, 43};
    s.final_three_phase_stability = {6e-10, 7e-8, 0.004, 6.0, 47, 53, 59, false, 61};
    s.two_phase = {2e-11, 3e-12, 4e-10, 5e-9, 6e-7, 7e-8, 1.5, 0.008, 0.009,
                   67, 71, 73, 79, 11};
    s.three_phase = {8e-11, 9e-12, 1e-9, 3e-13, 2e-9, 3e-7, 1.8, 0.006,
                     83, 89, 97, 101, 103, 13, 0.07};
    mc::PtSolverSafetyLimits limits;
    limits.max_components = 11;
    limits.max_stability_start_entries = 1001;
    limits.max_three_phase_starts = 13;
    limits.max_three_phase_start_entries = 3003;
    const auto prepared = mc::Pr76SolverConfiguration::create(s, limits);
    // Independent direct C++ construction, all user fields deliberately distinct.
    const fl::Pr76PtFlashBackendOptions expected{
        {{2e-11, 3e-12, 4e-10, 5e-9, 6e-7, 7e-8, 1.5, 0.008, 67, 71, 73, {79, 11}, 0.009},
         {2e-10, 3e-8, 0.002, 3.0, 17, 19, 23, 11, 29, 1001, false},
         {4e-10, 5e-8, 0.003, 5.0, 31, 37, 41, 11, 43, 1001, true}, 11},
        {8e-11, 9e-12, 1e-9, 3e-13, 2e-9, 3e-7, 1.8, 0.006, 83, 89, 97, 101, 103, 11},
        {6e-10, 7e-8, 0.004, 6.0, 47, 53, 59, 11, 61, 1001, false},
        13, 0.07, {}, 13, 3003, {}, {}};
    require_options_equal(prepared.backend_options(), expected);
    require(prepared.settings() == s && prepared.safety_limits() == limits &&
            prepared.root_options().max_iterations == 37, "snapshot or root mapping");
}

// Exhaustively remove each of the 57 numeric/boolean fields. This catches
// accidental default inheritance, including when the omitted value is false/zero.
template <class Selector, class... Member>
void omit_each(Selector select, Member... members) {
    const auto omit = [&](auto member) {
        auto s = custom();
        (select(s).*member).reset();
        expect_error(Code::missing_field, [&] { (void)mc::Pr76SolverConfiguration::create(s); });
    };
    (omit(members), ...);
}
void required_fields() {
    using S = mc::PtStabilitySettings;
    const auto stability = [](auto select) {
        omit_each(select, &S::tpd_tolerance, &S::stationarity_tolerance,
                  &S::line_search_armijo_coefficient, &S::max_log_composition_step,
                  &S::max_iterations, &S::max_backtracks, &S::max_property_evaluations,
                  &S::automatic_multistart, &S::max_starts);
    };
    stability([](auto& s) -> auto& { return s.initial_stability; });
    stability([](auto& s) -> auto& { return s.final_two_phase_stability; });
    stability([](auto& s) -> auto& { return s.final_three_phase_stability; });
    using T = mc::PtTwoPhaseSettings;
    omit_each([](auto& s) -> auto& { return s.two_phase; },
        &T::fugacity_equilibrium_tolerance, &T::absolute_mass_balance_tolerance,
        &T::relative_mass_balance_tolerance, &T::minimum_phase_fraction,
        &T::minimum_log_composition_separation, &T::minimum_relative_z_separation,
        &T::max_log_equilibrium_ratio_step, &T::residual_progress_coefficient,
        &T::gibbs_progress_coefficient, &T::max_iterations, &T::max_backtracks,
        &T::max_property_evaluations, &T::rachford_rice_max_iterations, &T::max_split_attempts);
    using H = mc::PtThreePhaseSettings;
    omit_each([](auto& s) -> auto& { return s.three_phase; },
        &H::chemical_potential_tolerance, &H::absolute_mass_balance_tolerance,
        &H::relative_mass_balance_tolerance, &H::generalized_rr_balance_tolerance,
        &H::minimum_phase_fraction, &H::minimum_log_composition_separation,
        &H::max_log_step, &H::residual_progress_coefficient, &H::max_iterations,
        &H::max_line_search_backtracks, &H::max_balance_iterations,
        &H::max_balance_backtracks, &H::max_property_evaluations,
        &H::max_three_phase_attempts, &H::new_phase_seed_fraction);
    omit_each([](auto& s) -> auto& { return s.eos_root; }, &mc::PtEosRootSettings::max_iterations);
}
void invalid_domains() {
    for (double bad : {std::numeric_limits<double>::quiet_NaN(), std::numeric_limits<double>::infinity()}) {
        auto s = custom();
        s.initial_stability.tpd_tolerance = bad;
        expect_error(Code::invalid_value, [&] { (void)mc::Pr76SolverConfiguration::create(s); },
                     "initial_stability.tpd_tolerance");
    }
    const auto invalid = [](auto mutate) {
        auto s = custom();
        mutate(s);
        expect_error(Code::invalid_settings, [&] { (void)mc::Pr76SolverConfiguration::create(s); });
    };
    invalid([](auto& s) { s.eos_root.max_iterations = 0; });
    invalid([](auto& s) { s.initial_stability.tpd_tolerance = -1; });
    invalid([](auto& s) { s.initial_stability.line_search_armijo_coefficient = 1; });
    invalid([](auto& s) { s.final_two_phase_stability.max_property_evaluations = 0; });
    invalid([](auto& s) { s.final_three_phase_stability.stationarity_tolerance = 0; });
    invalid([](auto& s) { s.two_phase.minimum_phase_fraction = 0.5; });
    invalid([](auto& s) { s.two_phase.minimum_relative_z_separation = 1; });
    invalid([](auto& s) { s.two_phase.gibbs_progress_coefficient = 0; });
    invalid([](auto& s) { s.three_phase.minimum_phase_fraction = 1.0 / 3.0; });
    invalid([](auto& s) { s.three_phase.residual_progress_coefficient = 1; });
    invalid([](auto& s) { s.three_phase.max_property_evaluations = 0; });
    invalid([](auto& s) { s.three_phase.max_three_phase_attempts = 0; });
    invalid([](auto& s) { s.three_phase.new_phase_seed_fraction = 1; });
    auto s = custom();
    s.two_phase.max_property_evaluations = -1;
    expect_error(Code::invalid_value, [&] { (void)mc::Pr76SolverConfiguration::create(s); });
    s = custom();
    s.eos_root.max_iterations = std::numeric_limits<std::int64_t>::max();
    mc::PtSolverSafetyLimits large;
    large.max_root_iterations = std::numeric_limits<std::size_t>::max();
    expect_error(Code::invalid_value, [&] { (void)mc::Pr76SolverConfiguration::create(s, large); });
}
void resource_limits() {
    const auto exceeds = [](auto mutate) {
        auto s = custom();
        mutate(s);
        expect_error(Code::resource_limit, [&] { (void)mc::Pr76SolverConfiguration::create(s); });
    };
    exceeds([](auto& s) { s.eos_root.max_iterations = 8193; });
    exceeds([](auto& s) { s.initial_stability.max_iterations = 4097; });
    exceeds([](auto& s) { s.final_two_phase_stability.max_backtracks = 257; });
    exceeds([](auto& s) { s.final_three_phase_stability.max_property_evaluations = 1000001; });
    exceeds([](auto& s) { s.initial_stability.max_starts = 1025; });
    exceeds([](auto& s) { s.two_phase.max_split_attempts = 65; });
    exceeds([](auto& s) { s.three_phase.max_three_phase_attempts = 65; });
    const auto s = mc::resolve_pt_solver_preset(mc::mpmc_balanced_default_v1);
    mc::PtSolverSafetyLimits strict;
    strict.max_property_evaluations = 99999;
    expect_error(Code::resource_limit, [&] { (void)mc::Pr76SolverConfiguration::create(s, strict); });
    strict = {};
    strict.max_components = 0;
    expect_error(Code::invalid_value, [&] { (void)mc::Pr76SolverConfiguration::create(s, strict); },
                 "safety_limits.max_components");
    auto boundary = custom();
    boundary.eos_root.max_iterations = 8192;
    const auto prepared = mc::Pr76SolverConfiguration::create(boundary);
    require(prepared.root_options().max_iterations == 8192, "ceiling equality must pass unchanged");
}
void zero_false_semantics() {
    auto s = custom();
    s.initial_stability.tpd_tolerance = 0;
    s.initial_stability.max_iterations = 0;
    s.initial_stability.automatic_multistart = false;
    s.two_phase.max_iterations = 0;
    s.two_phase.rachford_rice_max_iterations = 0;
    s.two_phase.max_property_evaluations = 0; // Existing split permits this.
    s.two_phase.max_split_attempts = 0;       // Existing split permits this too.
    s.three_phase.max_iterations = 0;
    s.three_phase.max_balance_iterations = 0;
    const auto prepared = mc::Pr76SolverConfiguration::create(s);
    const auto& o = prepared.backend_options();
    require(o.split.initial_stability.tpd_tolerance == 0 && !o.split.initial_stability.automatic_starts &&
            o.split.initial_stability.max_iterations == 0 && o.split.iteration.max_iterations == 0 &&
            o.split.iteration.rr.max_iterations == 0 && o.split.iteration.max_evaluations == 0 &&
            o.split.max_split_attempts == 0 && o.three_phase.max_iterations == 0 &&
            o.three_phase.max_balance_iterations == 0, "zero/false repaired to defaults");
}
void snapshot_isolation() {
    static_assert(!std::is_copy_assignable_v<mc::Pr76SolverConfiguration>);
    static_assert(!std::is_move_assignable_v<mc::Pr76SolverConfiguration>);
    auto s = custom();
    const auto a = mc::Pr76SolverConfiguration::create(s);
    s.initial_stability.max_property_evaluations = 1;
    s.final_two_phase_stability.max_property_evaluations = 2;
    s.final_three_phase_stability.max_property_evaluations = 3;
    const auto b = mc::Pr76SolverConfiguration::create(s);
    s = {};
    require(a.settings().initial_stability.max_property_evaluations == 100000 &&
            a.backend_options().split.initial_stability.max_evaluations == 100000, "A was mutated");
    require(b.backend_options().split.initial_stability.max_evaluations == 1 &&
            b.backend_options().split.final_stability.max_evaluations == 2 &&
            b.backend_options().final_three_phase_stability.max_evaluations == 3, "stability groups alias");
    const auto fresh = mc::resolve_pt_solver_preset(mc::mpmc_balanced_default_v1);
    require(fresh.initial_stability.max_property_evaluations == 100000, "preset is mutable global state");
}
} // namespace
} // namespace solver_test

int main(int argc, char** argv) {
    using namespace solver_test;
    try {
        require(argc == 2, "one test name required");
        const std::string_view name(argv[1]);
        if (name == "preset") { preset(); }
        else if (name == "custom_mapping") { custom_mapping(); }
        else if (name == "required_fields") { required_fields(); }
        else if (name == "invalid_domains") { invalid_domains(); }
        else if (name == "resource_limits") { resource_limits(); }
        else if (name == "zero_false_semantics") { zero_false_semantics(); }
        else if (name == "snapshot_isolation") { snapshot_isolation(); }
        else if (name == "headers") { require(solver_public_header_compiles(), "public header"); }
        else { solve_parity(name); }
        std::cout << "PASS " << name << '\n';
        return 0;
    } catch (const std::exception& error) {
        std::cerr << error.what() << '\n';
        return 1;
    }
}
