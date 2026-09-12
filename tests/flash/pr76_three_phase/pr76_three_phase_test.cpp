#include <mpmc/flash/pr76_max3_phase_set.hpp>
#include <mpmc/flash/pr76_pt_flash_backend.hpp>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <iostream>
#include <optional>
#include <source_location>
#include <stdexcept>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

bool pr76_three_phase_headers();

namespace {
namespace fl = mpmc::flash;
namespace th = mpmc::thermodynamics;
using Vec = std::vector<double>;

void require(bool value, std::string_view message,
             std::source_location where = std::source_location::current()) {
    if (!value) {
        throw std::runtime_error(
            std::string(where.file_name()) + ":" + std::to_string(where.line()) +
            ": " + std::string(message));
    }
}

void near(double actual, double expected, double tolerance,
          std::string_view message) {
    require(std::isfinite(actual) &&
                std::abs(actual - expected) <=
                    tolerance * (1.0 + std::abs(expected)),
            message);
}

th::Provenance synthetic_source(std::string locator) {
    return {th::SourceKind::synthetic_test,
            "MPMC_HNU symmetric PR76 max3 structural fixture",
            "v1",
            std::move(locator),
            "Artificial PR parameters used only to exercise generic LLV-like topology; not experimental validation",
            "tests/flash/pr76_three_phase/pr76_three_phase_test.cpp",
            "Repository structural regression"};
}

th::SourcedScalar synthetic_scalar(double value, th::Unit unit,
                                    const th::Provenance& source) {
    return {value, unit, source, "synthetic SI or dimensionless", "identity"};
}

th::Pr76Phase<double> synthetic_pr_model(bool swap_heavy = false) {
    const auto source = synthetic_source("three-component symmetric max3 fixture");
    const std::array<std::string, 3> ids{"light", "heavy-b", "heavy-c"};
    const std::array<double, 3> tc{190.6, 500.0, 500.0};
    const std::array<double, 3> pc{4.6e6, 5.0e6, 5.0e6};
    const std::array<double, 3> omega{0.01, 0.10, 0.10};

    std::vector<th::Component> catalog;
    th::PrParameterInput input;
    input.model_id = std::string(th::pr76_profile);
    input.dataset_id = "synthetic-PR76-max3-symmetric";
    input.revision = "v1";
    input.applicability = {std::nullopt, std::nullopt, source};
    for (std::size_t i = 0; i < ids.size(); ++i) {
        catalog.push_back({ids[i], ids[i], th::ComponentKind::pure, source, {}});
        input.pure.push_back({
            ids[i],
            synthetic_scalar(tc[i], th::Unit::kelvin, source),
            synthetic_scalar(pc[i], th::Unit::pascal, source),
            synthetic_scalar(omega[i], th::Unit::dimensionless, source)});
    }
    input.binary.push_back({
        ids[0], ids[1], synthetic_scalar(0.05, th::Unit::dimensionless, source)});
    input.binary.push_back({
        ids[0], ids[2], synthetic_scalar(0.05, th::Unit::dimensionless, source)});
    input.binary.push_back({
        ids[1], ids[2], synthetic_scalar(0.20, th::Unit::dimensionless, source)});

    const std::vector<std::string> order = swap_heavy
        ? std::vector<std::string>{ids[0], ids[2], ids[1]}
        : std::vector<std::string>{ids.begin(), ids.end()};
    return th::Pr76Phase<double>::from_parameters(
        th::PrParameterSet::create(
            catalog, order, input, th::DataPolicy::allow_synthetic_tests));
}

const std::array<Vec, 3>& synthetic_reference_phases() {
    // Independent offline solve of the same synthetic PR76 equations at
    // p=1 MPa, T=250 K. These values are numerical structural anchors only;
    // they are not experimental or a physical parameter validation.
    static const std::array<Vec, 3> value{{
        {0.04844395156186141, 0.9076712804842868, 0.0438847679538518},
        {0.9696550397025219, 0.015172480148739028, 0.015172480148739028},
        {0.048443951561861175, 0.04388476795385164, 0.9076712804842871}}};
    return value;
}

Vec log_ratio(const Vec& numerator, const Vec& denominator) {
    require(numerator.size() == denominator.size(), "log-ratio shape mismatch");
    Vec value(numerator.size(), 0.0);
    for (std::size_t i = 0; i < value.size(); ++i) {
        value[i] = std::log(numerator[i]) - std::log(denominator[i]);
    }
    return value;
}

void generic_exact_three_phase() {
    const std::array<Vec, 3> phases{{
        {0.8, 0.1, 0.1}, {0.1, 0.8, 0.1}, {0.1, 0.1, 0.8}}};
    const Vec feed{1.0 / 3.0, 1.0 / 3.0, 1.0 / 3.0};
    const auto provider = [&](double, double, std::span<const double>, std::size_t slot) {
        fl::PtThreePhaseProperty property;
        property.activity.ln_phi.resize(3U);
        for (std::size_t i = 0; i < 3U; ++i) {
            property.activity.ln_phi[i] = -std::log(phases[slot][i]);
        }
        property.activity.branch = slot;
        property.activity.smooth = true;
        property.z = 0.1 + 0.4 * static_cast<double>(slot);
        return property;
    };
    const auto result = fl::iterate_pt_three_phase(
        1.0e6, 300.0, feed,
        log_ratio(phases[1], phases[0]),
        log_ratio(phases[2], phases[0]),
        {1.0 / 3.0, 1.0 / 3.0}, provider);
    require(result.status == fl::PtThreePhaseStatus::converged_candidate &&
                result.candidate() != nullptr && result.iterations == 0,
            "exact manufactured three-phase state did not converge immediately");
    require(result.candidate()->chemical_potential_norm <=
                result.options.chemical_potential_tolerance &&
                result.candidate()->mass_absolute <=
                    result.options.mass_absolute_tolerance &&
                result.candidate()->mass_relative <=
                    result.options.mass_relative_tolerance,
            "generic exact three-phase invariants failed");
}

void generic_disappearance_boundary() {
    const std::array<Vec, 3> phases{{
        {0.8, 0.1, 0.1}, {0.1, 0.8, 0.1}, {0.1, 0.1, 0.8}}};
    constexpr double beta1 = 0.4;
    constexpr double beta2 = 0.0;
    Vec feed(3U, 0.0);
    for (std::size_t i = 0; i < 3U; ++i) {
        feed[i] = (1.0 - beta1 - beta2) * phases[0][i] +
                  beta1 * phases[1][i] + beta2 * phases[2][i];
    }
    const auto provider = [&](double, double, std::span<const double>, std::size_t slot) {
        fl::PtThreePhaseProperty property;
        property.activity.ln_phi.resize(3U);
        for (std::size_t i = 0; i < 3U; ++i) {
            property.activity.ln_phi[i] = -std::log(phases[slot][i]);
        }
        property.activity.branch = slot;
        property.activity.smooth = true;
        property.z = 0.1 + 0.4 * static_cast<double>(slot);
        return property;
    };
    const auto result = fl::iterate_pt_three_phase(
        1.0e6, 300.0, feed,
        log_ratio(phases[1], phases[0]),
        log_ratio(phases[2], phases[0]),
        {beta1, beta2}, provider);
    require(result.status == fl::PtThreePhaseStatus::phase_disappearance &&
                result.equations_converged() && result.disappearing_phase &&
                *result.disappearing_phase == 2U,
            "zero-third-phase boundary was not retained as disappearance evidence");
}

void require_unordered_reference_match(const fl::PtThreePhaseState& state) {
    const auto& reference = synthetic_reference_phases();
    std::array<bool, 3> used{false, false, false};
    for (const auto& phase : state.phases) {
        std::optional<std::size_t> best;
        double best_error = 1.0e100;
        for (std::size_t candidate = 0; candidate < 3U; ++candidate) {
            if (used[candidate]) { continue; }
            double error = 0.0;
            for (std::size_t i = 0; i < 3U; ++i) {
                error = std::max(
                    error,
                    std::abs(phase.composition[i] - reference[candidate][i]));
            }
            if (error < best_error) { best_error = error; best = candidate; }
        }
        require(best.has_value() && best_error < 3.0e-6,
                "PR76 synthetic three-phase composition left independent structural anchor");
        used[*best] = true;
    }
}

void pr76_synthetic_max3() {
    const auto model = synthetic_pr_model();
    fl::Pr76VleEvaluator evaluator(model);
    const Vec feed{1.0 / 3.0, 1.0 / 3.0, 1.0 / 3.0};
    auto result = fl::solve_pr76_pt_max3(1.0e6, 250.0, feed, evaluator);
    require(result.status == fl::Pr76PtMax3Status::three_phase &&
                result.three_phase_candidate() != nullptr &&
                result.selected_attempt.has_value(),
            "synthetic PR76 LLV-like fixture did not close at three phases");
    const auto& candidate = *result.three_phase_candidate();
    require_unordered_reference_match(candidate);
    double beta_sum = 0.0;
    for (const auto& phase : candidate.phases) {
        beta_sum += phase.mole_phase_fraction;
        require(phase.mole_phase_fraction > result.options.three_phase.minimum_phase_fraction,
                "accepted synthetic three-phase state contains disappearing phase");
    }
    near(beta_sum, 1.0, 2e-12, "three-phase fractions do not sum to one");
    const auto& attempt = result.attempts[*result.selected_attempt];
    require(attempt.final_stability &&
                attempt.final_stability->status == fl::StabilityStatus::no_instability_found,
            "accepted PR76 three-phase state lacks final common-tangent review");

    const auto published = fl::project_pr76_pt_max3_phase_set(result);
    require(published.solution.status == fl::PtPhaseSetStatus::accepted &&
                published.solution.accepted_phase_count() == 3U &&
                published.solution.capability.maximum_phase_count == 3U,
            "accepted PR76 three-phase state was not published as generic phase set");
}

void pr76_fresh_two_phase_neighbor() {
    const auto model = synthetic_pr_model();
    fl::Pr76VleEvaluator evaluator(model);
    const auto& phases = synthetic_reference_phases();
    constexpr double vapor_fraction = 0.4;
    Vec feed(3U, 0.0);
    for (std::size_t i = 0; i < 3U; ++i) {
        feed[i] = (1.0 - vapor_fraction) * phases[0][i] +
                  vapor_fraction * phases[1][i];
    }

    fl::Pr76ThreePhaseEvaluator phase_evaluator(
        model,
        {fl::Pr76RootSide::lower_admissible,
         fl::Pr76RootSide::upper_admissible,
         fl::Pr76RootSide::lower_admissible});
    const auto boundary = fl::iterate_pt_three_phase(
        1.0e6, 250.0, feed,
        log_ratio(phases[1], phases[0]),
        log_ratio(phases[2], phases[0]),
        {vapor_fraction, 0.0}, phase_evaluator);
    require(boundary.status == fl::PtThreePhaseStatus::phase_disappearance &&
                boundary.equations_converged() && boundary.disappearing_phase &&
                *boundary.disappearing_phase == 2U,
            "PR76 exact three-phase edge did not expose the zero-share phase");

    const std::vector<Vec> starts{
        boundary.point->phases[0].composition,
        boundary.point->phases[1].composition};
    const auto neighbor = fl::solve_pr76_pt_vle(
        1.0e6, 250.0, feed, evaluator, {}, starts, starts);
    require(neighbor.solution.status ==
                fl::PtSplitStatus::two_phase_no_instability_found &&
                neighbor.solution.candidate() != nullptr &&
                neighbor.solution.final_stability &&
                neighbor.solution.final_stability->status ==
                    fl::StabilityStatus::no_instability_found,
            "surviving PR76 pair did not fresh-resolve as a closed two-phase neighbor");
}

void pr76_component_permutation() {
    const Vec feed{1.0 / 3.0, 1.0 / 3.0, 1.0 / 3.0};
    const auto first_model = synthetic_pr_model(false);
    const auto second_model = synthetic_pr_model(true);
    fl::Pr76VleEvaluator first_eval(first_model);
    fl::Pr76VleEvaluator second_eval(second_model);
    const auto first = fl::solve_pr76_pt_max3(1.0e6, 250.0, feed, first_eval);
    const auto second = fl::solve_pr76_pt_max3(1.0e6, 250.0, feed, second_eval);
    require(first.three_phase_candidate() && second.three_phase_candidate(),
            "component permutation lost synthetic three-phase solution");

    std::array<Vec, 3> remapped;
    for (std::size_t p = 0; p < 3U; ++p) {
        const auto& x = second.three_phase_candidate()->phases[p].composition;
        remapped[p] = {x[0], x[2], x[1]};
    }
    std::array<bool, 3> used{false, false, false};
    for (const auto& phase : first.three_phase_candidate()->phases) {
        std::optional<std::size_t> match;
        for (std::size_t q = 0; q < 3U; ++q) {
            if (used[q]) { continue; }
            double error = 0.0;
            for (std::size_t i = 0; i < 3U; ++i) {
                error = std::max(error, std::abs(phase.composition[i] - remapped[q][i]));
            }
            if (error < 3e-9) { match = q; break; }
        }
        require(match.has_value(), "three-phase solution changed under component permutation");
        used[*match] = true;
    }
}

void headers() {
    require(pr76_three_phase_headers(), "public header probe failed");
}

using Test = std::pair<std::string_view, void (*)()>;
constexpr Test tests[]{
    {"generic_exact", generic_exact_three_phase},
    {"generic_disappearance", generic_disappearance_boundary},
    {"pr76_synthetic_max3", pr76_synthetic_max3},
    {"pr76_fresh_two_phase_neighbor", pr76_fresh_two_phase_neighbor},
    {"pr76_component_permutation", pr76_component_permutation},
    {"headers", headers}};

} // namespace

int main(int argc, char** argv) {
    try {
        if (argc != 2) { throw std::invalid_argument("one test name required"); }
        for (const auto& [name, run] : tests) {
            if (name == argv[1]) {
                run();
                std::cout << "[PASS] " << name << '\n';
                return 0;
            }
        }
        throw std::invalid_argument("unknown test");
    } catch (const std::exception& error) {
        std::cerr << "[FAIL] " << error.what() << '\n';
        return 1;
    }
}
