#include <mpmc/ad/dual.hpp>
#include <mpmc/thermodynamics/pr76_pure.hpp>
#include <mpmc/thermodynamics/pr76_roots.hpp>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <iomanip>
#include <iostream>
#include <optional>
#include <stdexcept>
#include <string_view>
#include <utility>

namespace {
namespace ad = mpmc::ad;
namespace th = mpmc::thermodynamics;

constexpr double co2_tc_k = 304.21;
constexpr double co2_pc_pa = 7.383e6;
constexpr double co2_omega = 0.2236;
constexpr double decane_tc_k = 617.7;
constexpr double decane_pc_pa = 2.110e6;
constexpr double decane_omega = 0.4923;
constexpr double strict_kij_from_pr101 = 0.05226578047;

// Independent Wong-Sandler/NRTL parameters reported by Arenas-Quevedo et al.
// Fluid Phase Equilibria 338 (2013) 30-36, DOI 10.1016/j.fluid.2012.10.012,
// Eqs. (10)-(15), Table 6. Their binary regression source is independent of
// the UFC 2025 four-point dataset audited here.
constexpr double ws_k12 = 0.7155;
constexpr double nrtl_delta12_j_per_mol = 11.8841e3;
constexpr double nrtl_delta21_j_per_mol = -1.9705e3;
constexpr double nrtl_nonrandomness = 0.3;

struct ExperimentalPoint {
    double temperature_k;
    double x_co2;
    double pressure_pa;
};

// UFC 2025 Table 7 CO2+n-decane binary bubble points used by PR #101/#104.
constexpr std::array<ExperimentalPoint, 4> ufc_points{{
    {323.08, 0.328, 3.14e6},
    {322.96, 0.495, 5.51e6},
    {323.21, 0.777, 8.57e6},
    {323.01, 0.911, 9.47e6},
}};

void require(bool condition, const char* message) {
    if (!condition) {
        throw std::runtime_error(message);
    }
}

enum class MixingMode { classical_vdw1f, wong_sandler_nrtl };
enum class RootRole { liquid, vapor };

template <typename Number>
struct MixedCoefficients {
    Number a;
    Number b;
};

struct PurePair {
    th::Pr76Pure<double> co2;
    th::Pr76Pure<double> decane;
};

struct GenericPhase {
    double z{};
    std::array<double, 2> ln_phi{};
};

struct YResidual {
    double y_co2{};
    double difference{};
    double common_residual{};
    GenericPhase vapor{};
};

struct BranchState {
    double common_residual{};
    double component_residual_gap{};
    double y_co2{};
    double liquid_z{};
    double vapor_z{};
};

struct BubbleResult {
    double pressure_pa{};
    double bracket_width_pa{};
    BranchState state{};
};

struct ModelMetrics {
    std::array<BubbleResult, 4> bubbles{};
    double aard{};
    double relative_sse{};
};

PurePair make_pures() {
    return {th::Pr76Pure<double>{co2_tc_k, co2_pc_pa, co2_omega},
            th::Pr76Pure<double>{decane_tc_k, decane_pc_pa, decane_omega}};
}

template <typename Number>
MixedCoefficients<Number> mix_from_moles(
    double temperature_k,
    const std::array<Number, 2>& mole_numbers,
    const PurePair& pures,
    MixingMode mode) {
    const Number total = mole_numbers[0] + mole_numbers[1];
    const std::array<Number, 2> x{mole_numbers[0] / total, mole_numbers[1] / total};
    const auto co2_value = pures.co2.evaluate(temperature_k);
    const auto decane_value = pures.decane.evaluate(temperature_k);
    const std::array<double, 2> a_i{co2_value.a, decane_value.a};
    const std::array<double, 2> b_i{co2_value.b, decane_value.b};
    const double rt = th::Pr76Pure<double>::gas_constant() * temperature_k;

    if (mode == MixingMode::classical_vdw1f) {
        const double a12 = std::sqrt(a_i[0] * a_i[1]) * (1.0 - strict_kij_from_pr101);
        const Number a_mix = x[0] * x[0] * a_i[0] +
                             Number{2.0 * a12} * x[0] * x[1] +
                             x[1] * x[1] * a_i[1];
        const Number b_mix = x[0] * b_i[0] + x[1] * b_i[1];
        return {a_mix, b_mix};
    }

    const double tau12 = nrtl_delta12_j_per_mol / rt;
    const double tau21 = nrtl_delta21_j_per_mol / rt;
    const double g12 = std::exp(-nrtl_nonrandomness * tau12);
    const double g21 = std::exp(-nrtl_nonrandomness * tau21);
    const std::array<std::array<double, 2>, 2> tau{{{0.0, tau12}, {tau21, 0.0}}};
    const std::array<std::array<double, 2>, 2> g{{{1.0, g12}, {g21, 1.0}}};

    Number gex_over_rt{0.0};
    for (std::size_t i = 0; i < 2; ++i) {
        Number numerator{0.0};
        Number denominator{0.0};
        for (std::size_t j = 0; j < 2; ++j) {
            numerator += x[j] * (tau[j][i] * g[j][i]);
            denominator += x[j] * g[j][i];
        }
        gex_over_rt += x[i] * (numerator / denominator);
    }

    Number q{0.0};
    for (std::size_t i = 0; i < 2; ++i) {
        for (std::size_t j = 0; j < 2; ++j) {
            const double kij = i == j ? 0.0 : ws_k12;
            const double cross = 0.5 *
                ((b_i[i] - a_i[i] / rt) + (b_i[j] - a_i[j] / rt)) *
                (1.0 - kij);
            q += x[i] * x[j] * cross;
        }
    }

    const double c = std::log(std::sqrt(2.0) - 1.0) / std::sqrt(2.0);
    Number d = gex_over_rt / c;
    for (std::size_t i = 0; i < 2; ++i) {
        d += x[i] * (a_i[i] / (b_i[i] * rt));
    }
    const Number b_mix = q / (Number{1.0} - d);
    const Number a_mix = Number{rt} * q * d / (Number{1.0} - d);
    return {a_mix, b_mix};
}

std::optional<std::size_t> stable_root_index(const th::Pr76RootSet<double>& roots,
                                             RootRole role) {
    if (role == RootRole::liquid) {
        for (std::size_t i = 0; i < roots.count; ++i) {
            if (roots.roots[i].slope_sign > 0) {
                return i;
            }
        }
    } else {
        for (std::size_t i = roots.count; i > 0; --i) {
            if (roots.roots[i - 1].slope_sign > 0) {
                return i - 1;
            }
        }
    }
    return std::nullopt;
}

std::optional<GenericPhase> evaluate_phase(double pressure_pa,
                                           double temperature_k,
                                           const std::array<double, 2>& fractions,
                                           RootRole role,
                                           const PurePair& pures,
                                           MixingMode mode) {
    const auto mixed = mix_from_moles<double>(temperature_k, fractions, pures, mode);
    const double rt = th::Pr76Pure<double>::gas_constant() * temperature_k;
    const double p_over_rt = pressure_pa / rt;
    const double eos_a = (mixed.a / rt) * p_over_rt;
    const double eos_b = mixed.b * p_over_rt;
    if (!std::isfinite(eos_a) || !std::isfinite(eos_b) || !(eos_b > 0.0)) {
        return std::nullopt;
    }

    const auto roots = th::pr76_roots(eos_a, eos_b);
    if (roots.status != th::Pr76RootStatus::success) {
        return std::nullopt;
    }
    const auto root_index = stable_root_index(roots, role);
    if (!root_index) {
        return std::nullopt;
    }
    const double z = roots.roots[*root_index].z;
    const double molar_volume = z * rt / pressure_pa;
    if (!(molar_volume > mixed.b)) {
        return std::nullopt;
    }

    using Dual = ad::Dual<double, 2>;
    const std::array<Dual, 2> moles{
        Dual::variable(fractions[0], 0), Dual::variable(fractions[1], 1)};
    const Dual n = moles[0] + moles[1];
    const auto mixed_dual = mix_from_moles<Dual>(temperature_k, moles, pures, mode);
    const Dual nb = n * mixed_dual.b;
    const Dual n2a = (n * n) * mixed_dual.a;

    const double sqrt2 = std::sqrt(2.0);
    const double log_ratio = std::log(
        (molar_volume + mixed.b * (1.0 - sqrt2)) /
        (molar_volume + mixed.b * (1.0 + sqrt2)));
    if (!std::isfinite(log_ratio)) {
        return std::nullopt;
    }
    const double first = -std::log(pressure_pa * (molar_volume - mixed.b) / rt);
    GenericPhase result;
    result.z = z;
    for (std::size_t i = 0; i < 2; ++i) {
        const double d_nb = nb.derivative(i);
        const double d_n2a_over_n = n2a.derivative(i) / n.value();
        const double bracket = d_n2a_over_n / mixed.a - d_nb / mixed.b;
        result.ln_phi[i] = first + (d_nb / mixed.b) * (z - 1.0) +
            (1.0 / (2.0 * sqrt2)) * (mixed.a / (mixed.b * rt)) * bracket * log_ratio;
        if (!std::isfinite(result.ln_phi[i])) {
            return std::nullopt;
        }
    }
    return result;
}

std::optional<YResidual> y_residual(double pressure_pa,
                                    double temperature_k,
                                    const std::array<double, 2>& liquid_x,
                                    const GenericPhase& liquid,
                                    double y_co2,
                                    const PurePair& pures,
                                    MixingMode mode) {
    if (!(y_co2 > liquid_x[0] + 1.0e-8) || !(y_co2 < 1.0 - 1.0e-10)) {
        return std::nullopt;
    }
    const std::array<double, 2> vapor_y{y_co2, 1.0 - y_co2};
    const auto vapor = evaluate_phase(
        pressure_pa, temperature_k, vapor_y, RootRole::vapor, pures, mode);
    if (!vapor) {
        return std::nullopt;
    }
    const double r0 = std::log(vapor_y[0]) + vapor->ln_phi[0] -
                      std::log(liquid_x[0]) - liquid.ln_phi[0];
    const double r1 = std::log(vapor_y[1]) + vapor->ln_phi[1] -
                      std::log(liquid_x[1]) - liquid.ln_phi[1];
    if (!std::isfinite(r0) || !std::isfinite(r1)) {
        return std::nullopt;
    }
    return YResidual{y_co2, r0 - r1, 0.5 * (r0 + r1), *vapor};
}

bool brackets_zero(double a, double b) {
    return a == 0.0 || b == 0.0 || std::signbit(a) != std::signbit(b);
}

std::optional<BranchState> stationary_vapor(double pressure_pa,
                                            const ExperimentalPoint& point,
                                            const PurePair& pures,
                                            MixingMode mode) {
    const std::array<double, 2> liquid_x{point.x_co2, 1.0 - point.x_co2};
    const auto liquid = evaluate_phase(
        pressure_pa, point.temperature_k, liquid_x, RootRole::liquid, pures, mode);
    if (!liquid) {
        return std::nullopt;
    }

    constexpr int intervals = 480;
    const double lower = point.x_co2 + 1.0e-5;
    const double upper = 1.0 - 1.0e-8;
    std::optional<YResidual> previous;
    std::optional<std::pair<YResidual, YResidual>> best;
    for (int index = 0; index <= intervals; ++index) {
        const double f = static_cast<double>(index) / static_cast<double>(intervals);
        const double y = lower + (upper - lower) * f;
        const auto current = y_residual(
            pressure_pa, point.temperature_k, liquid_x, *liquid, y, pures, mode);
        if (!current) {
            previous.reset();
            continue;
        }
        if (previous && brackets_zero(previous->difference, current->difference)) {
            if (!best || 0.5 * (previous->y_co2 + current->y_co2) >
                             0.5 * (best->first.y_co2 + best->second.y_co2)) {
                best = std::make_pair(*previous, *current);
            }
        }
        previous = *current;
    }
    if (!best) {
        return std::nullopt;
    }

    YResidual left = best->first;
    YResidual right = best->second;
    for (int iteration = 0; iteration < 72 && right.y_co2 - left.y_co2 > 1.0e-13; ++iteration) {
        const double middle_y = 0.5 * (left.y_co2 + right.y_co2);
        const auto middle = y_residual(
            pressure_pa, point.temperature_k, liquid_x, *liquid, middle_y, pures, mode);
        if (!middle) {
            return std::nullopt;
        }
        if (brackets_zero(left.difference, middle->difference)) {
            right = *middle;
        } else {
            left = *middle;
        }
    }

    const double y = 0.5 * (left.y_co2 + right.y_co2);
    const auto final = y_residual(
        pressure_pa, point.temperature_k, liquid_x, *liquid, y, pures, mode);
    if (!final) {
        return std::nullopt;
    }
    const double relative_z_gap = std::abs(final->vapor.z - liquid->z) /
                                  std::max(final->vapor.z, liquid->z);
    if (!(y - point.x_co2 > 1.0e-4) || !(relative_z_gap > 1.0e-5)) {
        return std::nullopt;
    }
    return BranchState{final->common_residual,
                       std::abs(final->difference),
                       y,
                       liquid->z,
                       final->vapor.z};
}

std::optional<BubbleResult> solve_bubble(const ExperimentalPoint& point,
                                         const PurePair& pures,
                                         MixingMode mode) {
    constexpr double step_pa = 25.0e3;
    constexpr int max_steps = 360;
    std::optional<double> previous_pressure;
    std::optional<BranchState> previous_state;
    double left_pa = 0.0;
    double right_pa = 0.0;
    BranchState left_state{};
    BranchState right_state{};
    bool found = false;

    for (int offset = -max_steps; offset <= max_steps; ++offset) {
        const double pressure = point.pressure_pa + static_cast<double>(offset) * step_pa;
        if (!(pressure > 0.1e6)) {
            continue;
        }
        const auto state = stationary_vapor(pressure, point, pures, mode);
        if (!state || !std::isfinite(state->common_residual)) {
            previous_pressure.reset();
            previous_state.reset();
            continue;
        }
        if (previous_pressure && previous_state &&
            brackets_zero(previous_state->common_residual, state->common_residual)) {
            left_pa = *previous_pressure;
            right_pa = pressure;
            left_state = *previous_state;
            right_state = *state;
            found = true;
            break;
        }
        previous_pressure = pressure;
        previous_state = *state;
    }
    if (!found) {
        return std::nullopt;
    }

    for (int iteration = 0; iteration < 72 && right_pa - left_pa > 0.05; ++iteration) {
        const double middle_pa = 0.5 * (left_pa + right_pa);
        const auto middle = stationary_vapor(middle_pa, point, pures, mode);
        if (!middle || !std::isfinite(middle->common_residual)) {
            return std::nullopt;
        }
        if (brackets_zero(left_state.common_residual, middle->common_residual)) {
            right_pa = middle_pa;
            right_state = *middle;
        } else {
            left_pa = middle_pa;
            left_state = *middle;
        }
    }

    const double pressure = 0.5 * (left_pa + right_pa);
    const auto state = stationary_vapor(pressure, point, pures, mode);
    if (!state) {
        return std::nullopt;
    }
    return BubbleResult{pressure, right_pa - left_pa, *state};
}

ModelMetrics evaluate_model(const PurePair& pures, MixingMode mode, std::string_view label) {
    ModelMetrics metrics;
    double absolute_relative_sum = 0.0;
    double relative_sse = 0.0;

    std::cout << label << '\n';
    for (std::size_t index = 0; index < ufc_points.size(); ++index) {
        const auto bubble = solve_bubble(ufc_points[index], pures, mode);
        require(bubble.has_value(), "four-point bubble solve failed on a UFC state");
        require(bubble->bracket_width_pa <= 0.05,
                "four-point bubble-pressure bracket did not converge tightly enough");
        require(std::abs(bubble->state.common_residual) <= 1.0e-8,
                "four-point bubble root did not close common fugacity residual");
        require(bubble->state.component_residual_gap <= 1.0e-10,
                "four-point vapor composition did not close component residual difference");
        require(bubble->state.y_co2 > ufc_points[index].x_co2,
                "four-point vapor is not CO2-richer than the liquid");
        require(bubble->state.vapor_z > bubble->state.liquid_z,
                "four-point liquid/vapor root ordering is not distinct");

        metrics.bubbles[index] = *bubble;
        const double relative =
            (bubble->pressure_pa - ufc_points[index].pressure_pa) / ufc_points[index].pressure_pa;
        absolute_relative_sum += std::abs(relative);
        relative_sse += relative * relative;

        std::cout << "  T_K=" << ufc_points[index].temperature_k
                  << " xCO2=" << ufc_points[index].x_co2
                  << " Pexp_MPa=" << ufc_points[index].pressure_pa / 1.0e6
                  << " Pcalc_MPa=" << bubble->pressure_pa / 1.0e6
                  << " signed_error_pct=" << 100.0 * relative
                  << " yCO2=" << bubble->state.y_co2
                  << " ZL=" << bubble->state.liquid_z
                  << " ZV=" << bubble->state.vapor_z << '\n';
    }
    metrics.aard = absolute_relative_sum / static_cast<double>(ufc_points.size());
    metrics.relative_sse = relative_sse;
    std::cout << "  AARD_pct=" << 100.0 * metrics.aard
              << " relative_SSE=" << metrics.relative_sse << '\n';
    return metrics;
}

void audit_four_point_ab() {
    const PurePair pures = make_pures();
    std::cout << std::setprecision(14)
              << "CO2+n-C10 UFC-2025 four-point A/B; independent WS/NRTL params fixed\n";

    const auto classical = evaluate_model(pures, MixingMode::classical_vdw1f, "classical-vdW1f");
    const auto ws = evaluate_model(pures, MixingMode::wong_sandler_nrtl, "Wong-Sandler/NRTL");

    // Guard the A arm against silently drifting away from the merged PR #101
    // strict-PR76 calibration baseline. These values are not WS acceptance targets.
    require(std::abs(classical.aard - 0.13249929377) <= 2.0e-4,
            "classical four-point AARD no longer matches the merged strict-PR76 baseline");
    require(std::abs(classical.bubbles[3].pressure_pa / 1.0e6 - 8.92879) <= 2.0e-3,
            "classical high-CO2 bubble no longer matches the merged strict-PR76 branch");
    require(std::isfinite(ws.aard) && std::isfinite(ws.relative_sse),
            "WS four-point aggregate metrics are not finite");

    std::cout << "A/B summary: classical_AARD_pct=" << 100.0 * classical.aard
              << " WS_AARD_pct=" << 100.0 * ws.aard
              << " AARD_delta_pct_points=" << 100.0 * (ws.aard - classical.aard)
              << " classical_relative_SSE=" << classical.relative_sse
              << " WS_relative_SSE=" << ws.relative_sse << '\n';
}

} // namespace

int main() {
    try {
        audit_four_point_ab();
        std::cout << "[PASS] pr76_c10_ws_four_point_ab_audit\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "[FAIL] " << error.what() << '\n';
        return 1;
    }
}
