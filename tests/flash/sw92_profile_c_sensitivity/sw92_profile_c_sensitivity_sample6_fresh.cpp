#include <mpmc/flash/sw92_profile_c_sensitivity.hpp>

#include "../sw92_phase_assigned_pt/physical_sample6.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <iostream>
#include <source_location>
#include <stdexcept>
#include <string>
#include <string_view>

namespace {
namespace fl = mpmc::flash;
namespace th = mpmc::thermodynamics;
namespace sample6 = sw92_profile_c_sample6;

void require(bool condition, std::string_view message,
             std::source_location where = std::source_location::current()) {
    if (!condition) {
        throw std::runtime_error(std::string(where.file_name()) + ":" +
                                 std::to_string(where.line()) + ": " +
                                 std::string(message));
    }
}

void near_increment(double analytic_increment, double resolved_half_difference,
                    double relative = 8e-3, double absolute = 4e-10,
                    std::source_location where = std::source_location::current()) {
    if (!std::isfinite(analytic_increment) ||
        !std::isfinite(resolved_half_difference) ||
        std::abs(analytic_increment - resolved_half_difference) >
            absolute + relative * std::max(
                std::abs(analytic_increment), std::abs(resolved_half_difference))) {
        std::cerr << "analytic_increment=" << analytic_increment
                  << " fresh_half_difference=" << resolved_half_difference << '\n';
        require(false, "Sample-6 fresh-resolve derivative mismatch", where);
    }
}

double density(double pressure_pa, double temperature_k, double z) {
    return pressure_pa /
        (z * th::Sw92Pure<double>::gas_constant() * temperature_k);
}

fl::Sw92ProfileCPtPhaseSetResult solve(double pressure_pa, double temperature_k,
                                       const th::Sw92Phase<double>& model) {
    return fl::solve_sw92_profile_c_pt_phase_set(
        pressure_pa, temperature_k, sample6::feed(), model, 0.0);
}

void require_three_phase(const fl::Sw92ProfileCPtPhaseSetResult& result) {
    require(result.accepted_phase_set_published() &&
                result.solution.accepted_phase_count() == 3U &&
                result.solution.accepted_phase_set() != nullptr,
            "Sample-6 perturbation left the fixed authoritative three-phase topology");
}

} // namespace

int main() {
    try {
        constexpr double pressure = 1.0e7;
        constexpr double temperature = 350.0;
        // Small enough to remain in the same physical interior topology, large
        // enough that fresh-solve roundoff is below the centered increment.
        const std::array<double,2> step{1000.0, 1.0e-2};

        const auto model = sample6::model();
        const auto base = solve(pressure, temperature, model);
        require_three_phase(base);
        const auto sensitivity =
            fl::differentiate_sw92_profile_c_phase_set(base, model);
        require(sensitivity.status == fl::PtSensitivityStatus::success,
                sensitivity.diagnostic);

        for (std::size_t column = 0; column < 2U; ++column) {
            double p_plus = pressure;
            double p_minus = pressure;
            double t_plus = temperature;
            double t_minus = temperature;
            if (column == 0U) {
                p_plus += step[column];
                p_minus -= step[column];
            } else {
                t_plus += step[column];
                t_minus -= step[column];
            }

            const auto plus = solve(p_plus, t_plus, model);
            const auto minus = solve(p_minus, t_minus, model);
            require_three_phase(plus);
            require_three_phase(minus);
            const auto& a = plus.solution.accepted_phase_set()->phases;
            const auto& b = minus.solution.accepted_phase_set()->phases;

            for (std::size_t phase = 0; phase < 3U; ++phase) {
                near_increment(
                    sensitivity.d_phase_fraction(phase,column) * step[column],
                    (a[phase].mole_phase_fraction -
                     b[phase].mole_phase_fraction) / 2.0);

                for (std::size_t component = 0; component < 8U; ++component) {
                    near_increment(
                        sensitivity.d_composition(phase,component,column) *
                            step[column],
                        (a[phase].composition[component] -
                         b[phase].composition[component]) / 2.0,
                        1.2e-2, 8e-11);
                }

                require(a[phase].compressibility_factor.has_value() &&
                            b[phase].compressibility_factor.has_value(),
                        "Sample-6 fresh perturbation lost accepted Z");
                near_increment(
                    sensitivity.d_compressibility(phase,column) * step[column],
                    (*a[phase].compressibility_factor -
                     *b[phase].compressibility_factor) / 2.0,
                    8e-3, 2e-10);

                const double rho_plus = density(
                    p_plus, t_plus, *a[phase].compressibility_factor);
                const double rho_minus = density(
                    p_minus, t_minus, *b[phase].compressibility_factor);
                near_increment(
                    sensitivity.d_molar_density(phase,column) * step[column],
                    (rho_plus-rho_minus) / 2.0,
                    1.0e-2, 2e-6);
            }
        }

        std::cout << "[PASS] Sample-6 fresh PT re-solve sensitivity cross-check\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "[FAIL] " << error.what() << '\n';
        return 1;
    }
}
