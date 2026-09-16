#include <mpmc/flash/pr76_three_phase.hpp>

#include "synthetic_fixture.hpp"

#include <array>
#include <cmath>
#include <iostream>
#include <source_location>
#include <stdexcept>
#include <string>
#include <string_view>
#include <vector>

bool pr76_three_phase_headers();

namespace {
namespace fl = mpmc::flash;
using Vec = std::vector<double>;

void require(bool value, std::string_view message,
             std::source_location where = std::source_location::current()) {
    if (!value) {
        throw std::runtime_error(
            std::string(where.file_name()) + ":" + std::to_string(where.line()) +
            ": " + std::string(message));
    }
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

void generic_refine_with_backtracking() {
    const std::array<Vec, 3> phases{{
        {0.8, 0.1, 0.1}, {0.1, 0.8, 0.1}, {0.1, 0.1, 0.8}}};
    const Vec feed{1.0 / 3.0, 1.0 / 3.0, 1.0 / 3.0};
    std::size_t injected_failures = 0U;
    const auto provider = [&](double, double, std::span<const double> composition,
                              std::size_t slot) {
        const double major = std::fmax(
            composition[0], std::fmax(composition[1], composition[2]));
        // With zero initial logK the first scaled full SSI step lands in this
        // narrow manufactured hole (~0.787 major fraction).  The half step is
        // valid and the exact target (0.8) is valid, so success requires the
        // outer logK loop to reject that trial, backtrack, and continue.
        if (major > 0.78 && major < 0.795) {
            ++injected_failures;
            throw fl::StabilityPropertyError(
                fl::StabilityPropertyIssue::nonfinite_properties,
                "manufactured recoverable three-phase trial rejection");
        }
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

    const Vec zero_log_k(3U, 0.0);
    const auto result = fl::iterate_pt_three_phase(
        1.0e6, 300.0, feed,
        zero_log_k, zero_log_k, {0.25, 0.25}, provider);
    require(result.status == fl::PtThreePhaseStatus::converged_candidate &&
                result.candidate() != nullptr,
            "perturbed manufactured three-phase state did not converge");
    require(result.iterations > 0 && result.backtracks > 0 &&
                result.rejected_evaluations > 0 && injected_failures > 0,
            "nontrivial logK update/backtracking path was not exercised");
    const auto& candidate = *result.candidate();
    require(candidate.chemical_potential_norm <=
                result.options.chemical_potential_tolerance &&
                candidate.generalized_rr_residual <=
                    result.options.balance_tolerance &&
                candidate.mass_absolute <= result.options.mass_absolute_tolerance &&
                candidate.mass_relative <= result.options.mass_relative_tolerance,
            "refined manufactured three-phase equations did not close");
    for (std::size_t phase = 0; phase < 3U; ++phase) {
        require(std::abs(candidate.phases[phase].mole_phase_fraction - 1.0 / 3.0) < 1e-12,
                "refined manufactured phase fraction changed");
        for (std::size_t component = 0; component < 3U; ++component) {
            require(std::abs(candidate.phases[phase].composition[component] -
                             phases[phase][component]) < 1e-12,
                    "refined manufactured phase composition changed");
        }
    }
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

void pr76_fresh_two_phase_neighbor() {
    const auto model = pr76_max3_test::model();
    fl::Pr76VleEvaluator evaluator(model);
    const auto& phases = pr76_max3_test::reference_phases();
    constexpr double second_fraction = 0.4;
    Vec feed(3U, 0.0);
    for (std::size_t i = 0; i < 3U; ++i) {
        feed[i] = (1.0 - second_fraction) * phases[0][i] +
                  second_fraction * phases[1][i];
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
        {second_fraction, 0.0}, phase_evaluator);
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

void headers() {
    require(pr76_three_phase_headers(), "public header probe failed");
}

using Test = std::pair<std::string_view, void (*)()>;
constexpr Test tests[]{
    {"generic_exact", generic_exact_three_phase},
    {"generic_refine", generic_refine_with_backtracking},
    {"generic_disappearance", generic_disappearance_boundary},
    {"pr76_fresh_two_phase_neighbor", pr76_fresh_two_phase_neighbor},
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
