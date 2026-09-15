#include <mpmc/flash/pr76_three_phase.hpp>

#include "heringer_2026_sour_gas_fixture.hpp"

#include <array>
#include <cmath>
#include <cstddef>
#include <iostream>
#include <numeric>
#include <stdexcept>
#include <vector>

namespace {
namespace fl = mpmc::flash;
namespace fx = pr76_heringer_2026_sour_gas_test;
using Vec = std::vector<double>;

constexpr double temperature_k = 178.8;
constexpr double co2_fraction = 0.73;
constexpr double pressure_three_phase_pa = 30.2e5;
constexpr double pressure_two_phase_pa = 35.0e5;

void require(bool condition, const char* message) {
    if (!condition) { throw std::runtime_error(message); }
}

Vec feed_at_co2_fraction(double target_co2) {
    require(std::isfinite(target_co2) && target_co2 > 0.0 && target_co2 < 1.0,
            "Heringer P-Z feed requires 0 < zCO2 < 1");

    // Heringer et al. vary the amount of CO2 mixed with the sour-gas system.
    // Build that P-Z path only from Appendix B Table B1: CO2 is prescribed and
    // every non-CO2 component retains its Table-B1 relative proportion.
    double source_non_co2 = 0.0;
    for (std::size_t i = 1U; i < fx::table_b1_mole_percent.size(); ++i) {
        source_non_co2 += fx::table_b1_mole_percent[i];
    }
    require(source_non_co2 > 0.0,
            "Heringer Table B1 lost non-CO2 composition support");

    Vec feed(fx::table_b1_mole_percent.size(), 0.0);
    feed[0] = target_co2;
    for (std::size_t i = 1U; i < feed.size(); ++i) {
        feed[i] = (1.0 - target_co2) *
                  fx::table_b1_mole_percent[i] / source_non_co2;
    }
    const double sum = std::accumulate(feed.begin(), feed.end(), 0.0);
    require(std::abs(sum - 1.0) <= 1.0e-14,
            "Heringer P-Z feed is not normalized");
    require(std::abs(feed[0] - target_co2) <= 1.0e-15,
            "Heringer P-Z feed changed requested zCO2");
    return feed;
}

void print_result(double pressure_pa, const fl::Pr76PtMax3Result& result) {
    std::cerr << "Pbar=" << pressure_pa / 1.0e5
              << " max3_status=" << static_cast<int>(result.status)
              << " base_status=" << static_cast<int>(result.base.solution.status)
              << " base_attempts=" << result.base.solution.attempts.size()
              << " max3_attempts=" << result.attempts.size()
              << " selected=" << (result.selected_attempt ? 1 : 0)
              << " diag=" << result.diagnostic << '\n';
    if (const auto* state = result.three_phase_candidate()) {
        std::cerr << "  three_phase fractions="
                  << state->phases[0].mole_phase_fraction << ','
                  << state->phases[1].mole_phase_fraction << ','
                  << state->phases[2].mole_phase_fraction
                  << " mu_norm=" << state->chemical_potential_norm
                  << " rr=" << state->generalized_rr_residual
                  << " mass_abs=" << state->mass_absolute
                  << " mass_rel=" << state->mass_relative << '\n';
    }
    if (const auto* pair = result.base.solution.candidate()) {
        std::cerr << "  base_pair beta_v=" << pair->fractions.vapor_fraction
                  << " fugacity_norm=" << pair->fugacity_norm
                  << " mass_abs=" << pair->fractions.mass_absolute
                  << " mass_rel=" << pair->fractions.mass_relative << '\n';
    }
}

fl::Pr76PtMax3Result solve_fresh(double pressure_pa, const Vec& feed) {
    const auto model = fx::model();
    fl::Pr76VleEvaluator evaluator(model);
    // Deliberately provide no literature phase compositions or continuation
    // hints. The finite TPD search and any 2->3 seed are generated internally.
    return fl::solve_pr76_pt_max3(
        pressure_pa, temperature_k, feed, evaluator);
}

void require_three_phase_closure(const fl::Pr76PtMax3Result& result) {
    require(result.base.dataset_id == "Heringer-et-al-2026-sour-gas-PR76",
            "30.2-bar solve lost Heringer-2026 dataset identity");
    require(result.status == fl::Pr76PtMax3Status::three_phase &&
                result.three_phase_candidate() != nullptr &&
                result.selected_attempt.has_value(),
            "30.2 bar, zCO2=0.73 did not resolve as three phase");
    const auto& attempt = result.attempts[*result.selected_attempt];
    require(!attempt.supplied_start.has_value() && attempt.witness_trial.has_value(),
            "30.2-bar three-phase solve did not use the automatic TPD-witness route");
    require(attempt.accepted_three_phase && attempt.final_stability.has_value() &&
                attempt.final_stability->status == fl::StabilityStatus::no_instability_found,
            "30.2-bar three-phase state lacks final finite common-tangent closure");

    const auto& state = *result.three_phase_candidate();
    require(state.chemical_potential_norm <= 1.0e-11,
            "30.2-bar three-phase fugacity/chemical-potential residual exceeds gate");
    require(state.generalized_rr_residual <= 2.0e-13,
            "30.2-bar three-phase generalized RR residual exceeds gate");
    require(state.mass_absolute <= 1.0e-12 && state.mass_relative <= 1.0e-10,
            "30.2-bar three-phase material balance exceeds gate");
    double fraction_sum = 0.0;
    for (const auto& phase : state.phases) {
        require(phase.mole_phase_fraction > 1.0e-10,
                "30.2-bar accepted phase fell below disappearance gate");
        fraction_sum += phase.mole_phase_fraction;
    }
    require(std::abs(fraction_sum - 1.0) <= 1.0e-13,
            "30.2-bar phase fractions do not sum to one");
    require(!result.global_stability_proven,
            "finite 30.2-bar stability search was promoted to global proof");
}

const fl::Pr76PtSplitResult* stable_two_phase_result(
    const fl::Pr76PtMax3Result& result) {
    if (result.base.solution.status == fl::PtSplitStatus::two_phase_no_instability_found) {
        return &result.base;
    }
    return result.two_phase_neighbor();
}

void require_two_phase_closure(const fl::Pr76PtMax3Result& result) {
    require(result.base.dataset_id == "Heringer-et-al-2026-sour-gas-PR76",
            "35-bar solve lost Heringer-2026 dataset identity");
    require(result.status == fl::Pr76PtMax3Status::two_phase,
            "35 bar, zCO2=0.73 did not resolve as two phase");
    const auto* pair_result = stable_two_phase_result(result);
    require(pair_result != nullptr &&
                pair_result->solution.status ==
                    fl::PtSplitStatus::two_phase_no_instability_found &&
                pair_result->solution.candidate() != nullptr &&
                pair_result->solution.final_stability.has_value() &&
                pair_result->solution.final_stability->status ==
                    fl::StabilityStatus::no_instability_found,
            "35-bar pair lacks fresh two-phase equilibrium/stability closure");

    const auto& state = *pair_result->solution.candidate();
    require(state.fugacity_norm <= 1.0e-11,
            "35-bar two-phase fugacity residual exceeds gate");
    require(state.fractions.mass_absolute <= 1.0e-12 &&
                state.fractions.mass_relative <= 1.0e-10,
            "35-bar two-phase material balance exceeds gate");
    require(state.fractions.vapor_fraction > 1.0e-10 &&
                state.fractions.vapor_fraction < 1.0 - 1.0e-10,
            "35-bar two-phase state lies on the phase-fraction disappearance gate");
    require(!result.global_stability_proven,
            "finite 35-bar stability search was promoted to global proof");
}

void heringer_source_consistent_three_to_two_boundary() {
    const auto feed = feed_at_co2_fraction(co2_fraction);

    // Heringer et al., FPE 604 (2026) 114653, Secs. 4.1.3-4.1.4:
    // at 30.2 bar the L1-L2 interval starts at 73.6 mol% CO2, so 73% is
    // immediately on the lower-CO2 three-phase side; at 35 bar the L1-L2
    // interval expands down to about 72 mol%, so the same 73% feed is two-phase.
    const auto lower = solve_fresh(pressure_three_phase_pa, feed);
    if (lower.status != fl::Pr76PtMax3Status::three_phase) {
        print_result(pressure_three_phase_pa, lower);
    }
    require_three_phase_closure(lower);

    const auto upper = solve_fresh(pressure_two_phase_pa, feed);
    if (upper.status != fl::Pr76PtMax3Status::two_phase) {
        print_result(pressure_two_phase_pa, upper);
    }
    require_two_phase_closure(upper);
}

} // namespace

int main() {
    try {
        heringer_source_consistent_three_to_two_boundary();
        std::cout << "[PASS] heringer_2026_source_consistent_3_to_2\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "[FAIL] " << error.what() << '\n';
        return 1;
    }
}
