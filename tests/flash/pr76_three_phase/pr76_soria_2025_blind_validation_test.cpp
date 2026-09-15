#include <mpmc/flash/pr76_three_phase.hpp>

#include "soria_2025_ternary_blind_fixture.hpp"

#include <cmath>
#include <iomanip>
#include <iostream>
#include <stdexcept>

namespace {
namespace fl = mpmc::flash;
namespace fx = pr76_soria_2025_ternary_blind_test;

void require(bool condition, const char* message) {
    if (!condition) { throw std::runtime_error(message); }
}

const char* max3_name(fl::Pr76PtMax3Status status) {
    switch (status) {
    case fl::Pr76PtMax3Status::single_phase: return "single";
    case fl::Pr76PtMax3Status::two_phase: return "two";
    case fl::Pr76PtMax3Status::three_phase: return "three";
    case fl::Pr76PtMax3Status::phase_boundary_unresolved: return "boundary_unresolved";
    case fl::Pr76PtMax3Status::higher_phase_count_or_wrong_candidate: return "higher_or_wrong";
    case fl::Pr76PtMax3Status::indeterminate: return "indeterminate";
    }
    return "unknown";
}

const char* split_name(fl::PtSplitStatus status) {
    switch (status) {
    case fl::PtSplitStatus::single_phase_no_instability_found: return "single_stable";
    case fl::PtSplitStatus::two_phase_no_instability_found: return "two_stable";
    case fl::PtSplitStatus::phase_set_unstable: return "phase_set_unstable";
    case fl::PtSplitStatus::indeterminate: return "indeterminate";
    }
    return "unknown";
}

fl::Pr76PtMax3Result solve_fresh(double pressure_pa) {
    const auto model = fx::model();
    fl::Pr76VleEvaluator evaluator(model);
    return fl::solve_pr76_pt_max3(
        pressure_pa, fx::experimental_temperature_k, fx::m40_feed, evaluator);
}

void print_point(double pressure_pa) {
    const auto result = solve_fresh(pressure_pa);
    std::cout << std::setprecision(10)
              << "P_MPa=" << pressure_pa / 1.0e6
              << " max3=" << max3_name(result.status)
              << " base=" << split_name(result.base.solution.status)
              << " init_stability=" << static_cast<int>(result.base.solution.initial_stability.status)
              << " attempts=" << result.attempts.size();
    if (const auto* pair = result.base.solution.candidate()) {
        std::cout << " beta_v=" << pair->fractions.vapor_fraction
                  << " fugacity=" << pair->fugacity_norm
                  << " mass_abs=" << pair->fractions.mass_absolute
                  << " mass_rel=" << pair->fractions.mass_relative;
    }
    if (const auto* three = result.three_phase_candidate()) {
        std::cout << " beta3="
                  << three->phases[0].mole_phase_fraction << ','
                  << three->phases[1].mole_phase_fraction << ','
                  << three->phases[2].mole_phase_fraction
                  << " mu=" << three->chemical_potential_norm
                  << " mass3=" << three->mass_absolute;
    }
    std::cout << " diag=" << result.diagnostic << '\n';
}

void audit_pressure_topology() {
    const auto model = fx::model();
    const auto& parameters = model.parameters();
    require(parameters.dataset_id() == "Soria-UFC-2025-M40-strict-PR76-blind",
            "blind fixture dataset identity drifted");
    require(std::abs(parameters.kij(0U, 1U) - 0.05226578047) <= 1.0e-14 &&
                std::abs(parameters.kij(0U, 2U) - 0.07884923875) <= 1.0e-14 &&
                parameters.kij(1U, 2U) == 0.0,
            "blind fixture kij drifted from PR #101");

    std::cout << "M40 strict-PR76 blind diagnostic pressure scan; "
              << "T_K=" << fx::experimental_temperature_k
              << " Pexp_MPa=" << fx::experimental_transition_pressure_pa / 1.0e6
              << " Uexp_MPa=" << fx::expanded_pressure_uncertainty_pa / 1.0e6 << '\n';

    for (int pressure_mpa = 2; pressure_mpa <= 20; ++pressure_mpa) {
        print_point(static_cast<double>(pressure_mpa) * 1.0e6);
    }
    print_point(fx::experimental_transition_pressure_pa - fx::expanded_pressure_uncertainty_pa);
    print_point(fx::experimental_transition_pressure_pa);
    print_point(fx::experimental_transition_pressure_pa + fx::expanded_pressure_uncertainty_pa);

    throw std::runtime_error(
        "diagnostic-only M-40 pressure topology scan complete; inspect hosted output before defining the blind transition gate");
}

} // namespace

int main() {
    try {
        audit_pressure_topology();
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "[AUDIT] " << error.what() << '\n';
        return 1;
    }
}
