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
#include <string>
#include <string_view>
#include <utility>
#include <vector>

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

// Frozen independent Wong-Sandler/NRTL parameters from Arenas-Quevedo et al.
// Fluid Phase Equilibria 338 (2013) 30-36, DOI 10.1016/j.fluid.2012.10.012,
// Table 6. Their CO2+n-decane regression used 29 literature VLE points over
// 319.11-372.94 K. Nothing below is refit to any validation dataset.
constexpr double ws_k12 = 0.7155;
constexpr double nrtl_delta12_j_per_mol = 11.8841e3;
constexpr double nrtl_delta21_j_per_mol = -1.9705e3;
constexpr double nrtl_nonrandomness = 0.3;
constexpr double arenas_fit_temperature_max_k = 372.94;

struct ExperimentalPoint {
    double temperature_k;
    double x_co2;
    double pressure_pa;
};

struct ExternalVlePoint {
    double temperature_k;
    double x_co2;
    double pressure_pa;
    double y_co2;
};

constexpr double atmosphere_pa = 101325.0;
constexpr double atm_to_pa(double pressure_atm) {
    return pressure_atm * atmosphere_pa;
}

// UFC 2025 Table 7 points retained as the already-merged near-323 K A/B guard.
constexpr std::array<ExperimentalPoint, 4> ufc_points{{
    {323.08, 0.328, 3.14e6},
    {322.96, 0.495, 5.51e6},
    {323.21, 0.777, 8.57e6},
    {323.01, 0.911, 9.47e6},
}};

// Independent middle-temperature transfer oracle: Inomata, Tuchiya, Arai & Saito,
// J. Chem. Eng. Japan 19 (1986) 386-391, DOI 10.1252/jcej.19.386,
// Table 4, CO2+n-decane at 411.2 K. The paper reports overall pressure accuracy
// of 0.05 MPa and equilibrium-composition error below 0.5 mol%.
// This source is not the Jimenez-Gallegos et al. 2006 dataset used to regress
// the frozen Arenas-Quevedo WS/NRTL parameters.
constexpr std::array<ExternalVlePoint, 6> inomata_411k_points{{
    {411.2, 0.254, 4.81e6, 0.986},
    {411.2, 0.467, 9.82e6, 0.989},
    {411.2, 0.543, 11.52e6, 0.980},
    {411.2, 0.659, 14.80e6, 0.969},
    {411.2, 0.739, 16.82e6, 0.953},
    {411.2, 0.804, 17.99e6, 0.925},
}};

// Independent far-temperature extrapolation oracle: Sebastian, Simnick, Lin & Chao,
// J. Chem. Eng. Data 25 (1980) 138-140, DOI 10.1021/je60085a012,
// Table I, CO2+n-decane at 189.4 degC = 462.55 K.
// This source predates, and is not the Jimenez-Gallegos et al. 2006 dataset
// used to regress the frozen Arenas-Quevedo WS/NRTL parameters.
constexpr std::array<ExternalVlePoint, 4> sebastian_462k_points{{
    {462.55, 0.0913, atm_to_pa(19.36), 0.9075},
    {462.55, 0.1472, atm_to_pa(30.38), 0.9306},
    {462.55, 0.1883, atm_to_pa(40.10), 0.9410},
    {462.55, 0.2358, atm_to_pa(50.70), 0.9478},
}};

static_assert(inomata_411k_points[0].temperature_k - arenas_fit_temperature_max_k > 30.0);
static_assert(sebastian_462k_points[0].temperature_k - arenas_fit_temperature_max_k > 80.0);

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

struct ExternalMetrics {
    double pressure_aard{};
    double pressure_relative_sse{};
    double mean_abs_y_error{};
    double max_abs_y_error{};
};

th::SourcedScalar scalar(double value,
                         th::Unit unit,
                         const th::Provenance& source,
                         std::string original_unit = "SI or dimensionless",
                         std::string conversion = "identity") {
    return {value, unit, source, std::move(original_unit), std::move(conversion)};
}

PurePair make_pures() {
    const th::Provenance source{
        th::SourceKind::literature,
        "https://repositorio.ufc.br/handle/riufc/80725",
        "E. C. Q. Soria, UFC dissertation, 2025",
        "Table 6: CO2 and n-decane Tc, Pc and acentric factor",
        "Same strict-PR76 pure inputs used by PR #101/#104/#106/#107/#108",
        "Transcribed from UFC repository dissertation",
        "Frozen pure inputs for mixing-rule validation; no parameter fitting"};

    std::vector<th::Component> catalog{
        {"carbon-dioxide", "carbon dioxide", th::ComponentKind::pure, source, {}},
        {"n-decane", "n-decane", th::ComponentKind::pure, source, {}}};
    const std::vector<std::string> order{"carbon-dioxide", "n-decane"};

    th::PrParameterInput input;
    input.model_id = std::string(th::pr76_profile);
    input.dataset_id = "C10-frozen-pure-WS-validation";
    input.revision = "multi-temperature-transfer-v2";
    input.applicability = {std::nullopt, std::nullopt, source};
    input.pure = {
        {"carbon-dioxide",
         scalar(co2_tc_k, th::Unit::kelvin, source, "K", "identity"),
         scalar(co2_pc_pa, th::Unit::pascal, source, "MPa", "MPa * 1e6 -> Pa"),
         scalar(co2_omega, th::Unit::dimensionless, source)},
        {"n-decane",
         scalar(decane_tc_k, th::Unit::kelvin, source, "K", "identity"),
         scalar(decane_pc_pa, th::Unit::pascal, source, "MPa", "MPa * 1e6 -> Pa"),
         scalar(decane_omega, th::Unit::dimensionless, source)}};
    input.binary = {{"carbon-dioxide", "n-decane",
                     scalar(strict_kij_from_pr101, th::Unit::dimensionless, source)}};

    const auto parameters = th::PrParameterSet::create(catalog, order, input);
    return {th::Pr76Pure<double>::from_parameters(parameters, 0U),
            th::Pr76Pure<double>::from_parameters(parameters, 1U)};
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
            const double cross =
                0.5 * ((b_i[i] - a_i[i] / rt) + (b_i[j] - a_i[j] / rt)) * (1.0 - kij);
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
        result.ln_phi[i] =
            first + (d_nb / mixed.b) * (z - 1.0) +
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
    const auto vapor =
        evaluate_phase(pressure_pa, temperature_k, vapor_y, RootRole::vapor, pures, mode);
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
    const auto liquid =
        evaluate_phase(pressure_pa, point.temperature_k, liquid_x, RootRole::liquid, pures, mode);
    if (!liquid) {
        return std::nullopt;
    }

    constexpr int intervals = 480;
    const double lower = point.x_co2 + 1.0e-5;
    const double upper = 1.0 - 1.0e-8;
    std::optional<YResidual> previous;
    std::optional<std::pair<YResidual, YResidual>> best;

    for (int index = 0; index <= intervals; ++index) {
        const double fraction = static_cast<double>(index) / static_cast<double>(intervals);
        const double y = lower + (upper - lower) * fraction;
        const auto current =
            y_residual(pressure_pa, point.temperature_k, liquid_x, *liquid, y, pures, mode);
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
        const auto middle =
            y_residual(pressure_pa, point.temperature_k, liquid_x, *liquid,
                       middle_y, pures, mode);
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
    const auto final =
        y_residual(pressure_pa, point.temperature_k, liquid_x, *liquid, y, pures, mode);
    if (!final) {
        return std::nullopt;
    }

    const double relative_z_gap =
        std::abs(final->vapor.z - liquid->z) / std::max(final->vapor.z, liquid->z);
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

void require_closed_bubble(const BubbleResult& bubble,
                           const ExperimentalPoint& point,
                           const char* context) {
    (void)context;
    require(bubble.pressure_pa > 0.0 && std::isfinite(bubble.pressure_pa),
            "bubble pressure is not finite and positive");
    require(bubble.bracket_width_pa <= 0.05,
            "bubble-pressure bracket did not converge tightly enough");
    require(std::abs(bubble.state.common_residual) <= 1.0e-8,
            "bubble root did not close common fugacity residual");
    require(bubble.state.component_residual_gap <= 1.0e-10,
            "vapor composition did not close component residual difference");
    require(bubble.state.y_co2 > point.x_co2,
            "bubble vapor is not CO2-richer than the liquid");
    require(bubble.state.vapor_z > bubble.state.liquid_z,
            "liquid/vapor root ordering is not distinct");
}

ModelMetrics evaluate_ufc_model(const PurePair& pures,
                                MixingMode mode,
                                std::string_view label) {
    ModelMetrics metrics;
    double absolute_relative_sum = 0.0;
    double relative_sse = 0.0;

    std::cout << label << '\n';
    for (std::size_t index = 0; index < ufc_points.size(); ++index) {
        const auto bubble = solve_bubble(ufc_points[index], pures, mode);
        require(bubble.has_value(), "four-point bubble solve failed on a UFC state");
        require_closed_bubble(*bubble, ufc_points[index], "UFC");

        metrics.bubbles[index] = *bubble;
        const double relative =
            (bubble->pressure_pa - ufc_points[index].pressure_pa) /
            ufc_points[index].pressure_pa;
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

void audit_four_point_ab(const PurePair& pures) {
    std::cout << std::setprecision(14)
              << "CO2+n-C10 UFC-2025 four-point A/B; independent WS/NRTL params fixed\n";

    const auto classical =
        evaluate_ufc_model(pures, MixingMode::classical_vdw1f, "classical-vdW1f");
    const auto ws =
        evaluate_ufc_model(pures, MixingMode::wong_sandler_nrtl, "Wong-Sandler/NRTL");

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

ExternalMetrics audit_inomata_midrange_transfer(const PurePair& pures) {
    ExternalMetrics metrics;
    double pressure_abs_relative_sum = 0.0;
    double pressure_relative_sse = 0.0;
    double y_abs_error_sum = 0.0;

    std::cout << "Inomata-1986 independent WS/NRTL middle-temperature transfer"
              << " source_DOI=10.1252/jcej.19.386"
              << " fit_Tmax_K=" << arenas_fit_temperature_max_k
              << " validation_T_K=" << inomata_411k_points[0].temperature_k
              << " extrapolation_delta_K="
              << inomata_411k_points[0].temperature_k - arenas_fit_temperature_max_k
              << '\n';

    for (const auto& external : inomata_411k_points) {
        const ExperimentalPoint point{
            external.temperature_k, external.x_co2, external.pressure_pa};
        const auto bubble =
            solve_bubble(point, pures, MixingMode::wong_sandler_nrtl);
        require(bubble.has_value(),
                "frozen WS/NRTL failed to produce an Inomata-1986 bubble state");
        require_closed_bubble(*bubble, point, "Inomata-1986");

        const double pressure_relative =
            (bubble->pressure_pa - external.pressure_pa) / external.pressure_pa;
        const double y_error = bubble->state.y_co2 - external.y_co2;
        pressure_abs_relative_sum += std::abs(pressure_relative);
        pressure_relative_sse += pressure_relative * pressure_relative;
        y_abs_error_sum += std::abs(y_error);
        metrics.max_abs_y_error = std::max(metrics.max_abs_y_error, std::abs(y_error));

        std::cout << "  xCO2=" << external.x_co2
                  << " Pexp_MPa=" << external.pressure_pa / 1.0e6
                  << " Pcalc_MPa=" << bubble->pressure_pa / 1.0e6
                  << " signed_pressure_error_pct=" << 100.0 * pressure_relative
                  << " yexp_CO2=" << external.y_co2
                  << " ycalc_CO2=" << bubble->state.y_co2
                  << " y_error=" << y_error
                  << " ZL=" << bubble->state.liquid_z
                  << " ZV=" << bubble->state.vapor_z << '\n';
    }

    const double count = static_cast<double>(inomata_411k_points.size());
    metrics.pressure_aard = pressure_abs_relative_sum / count;
    metrics.pressure_relative_sse = pressure_relative_sse;
    metrics.mean_abs_y_error = y_abs_error_sum / count;

    require(std::isfinite(metrics.pressure_aard) &&
                std::isfinite(metrics.pressure_relative_sse) &&
                std::isfinite(metrics.mean_abs_y_error) &&
                std::isfinite(metrics.max_abs_y_error),
            "Inomata transfer aggregate metrics are not finite");

    std::cout << "Inomata transfer summary:"
              << " pressure_AARD_pct=" << 100.0 * metrics.pressure_aard
              << " pressure_relative_SSE=" << metrics.pressure_relative_sse
              << " mean_abs_yCO2_error=" << metrics.mean_abs_y_error
              << " max_abs_yCO2_error=" << metrics.max_abs_y_error << '\n';
    return metrics;
}

ExternalMetrics audit_sebastian_extrapolation(const PurePair& pures) {
    ExternalMetrics metrics;
    double pressure_abs_relative_sum = 0.0;
    double pressure_relative_sse = 0.0;
    double y_abs_error_sum = 0.0;

    std::cout << "Sebastian-1980 independent WS/NRTL extrapolation"
              << " source_DOI=10.1021/je60085a012"
              << " fit_Tmax_K=" << arenas_fit_temperature_max_k
              << " validation_T_K=" << sebastian_462k_points[0].temperature_k
              << " extrapolation_delta_K="
              << sebastian_462k_points[0].temperature_k - arenas_fit_temperature_max_k
              << '\n';

    for (const auto& external : sebastian_462k_points) {
        const ExperimentalPoint point{
            external.temperature_k, external.x_co2, external.pressure_pa};
        const auto bubble =
            solve_bubble(point, pures, MixingMode::wong_sandler_nrtl);
        require(bubble.has_value(),
                "frozen WS/NRTL failed to produce a Sebastian-1980 bubble state");
        require_closed_bubble(*bubble, point, "Sebastian-1980");

        const double pressure_relative =
            (bubble->pressure_pa - external.pressure_pa) / external.pressure_pa;
        const double y_error = bubble->state.y_co2 - external.y_co2;
        pressure_abs_relative_sum += std::abs(pressure_relative);
        pressure_relative_sse += pressure_relative * pressure_relative;
        y_abs_error_sum += std::abs(y_error);
        metrics.max_abs_y_error = std::max(metrics.max_abs_y_error, std::abs(y_error));

        std::cout << "  xCO2=" << external.x_co2
                  << " Pexp_atm=" << external.pressure_pa / atmosphere_pa
                  << " Pexp_MPa=" << external.pressure_pa / 1.0e6
                  << " Pcalc_MPa=" << bubble->pressure_pa / 1.0e6
                  << " signed_pressure_error_pct=" << 100.0 * pressure_relative
                  << " yexp_CO2=" << external.y_co2
                  << " ycalc_CO2=" << bubble->state.y_co2
                  << " y_error=" << y_error
                  << " ZL=" << bubble->state.liquid_z
                  << " ZV=" << bubble->state.vapor_z << '\n';
    }

    const double count = static_cast<double>(sebastian_462k_points.size());
    metrics.pressure_aard = pressure_abs_relative_sum / count;
    metrics.pressure_relative_sse = pressure_relative_sse;
    metrics.mean_abs_y_error = y_abs_error_sum / count;

    require(std::isfinite(metrics.pressure_aard) &&
                std::isfinite(metrics.pressure_relative_sse) &&
                std::isfinite(metrics.mean_abs_y_error) &&
                std::isfinite(metrics.max_abs_y_error),
            "Sebastian extrapolation aggregate metrics are not finite");

    std::cout << "Sebastian extrapolation summary:"
              << " pressure_AARD_pct=" << 100.0 * metrics.pressure_aard
              << " pressure_relative_SSE=" << metrics.pressure_relative_sse
              << " mean_abs_yCO2_error=" << metrics.mean_abs_y_error
              << " max_abs_yCO2_error=" << metrics.max_abs_y_error << '\n';
    return metrics;
}

} // namespace

int main() {
    try {
        const PurePair pures = make_pures();
        audit_four_point_ab(pures);
        (void)audit_inomata_midrange_transfer(pures);
        (void)audit_sebastian_extrapolation(pures);
        std::cout << "[PASS] pr76_c10_wong_sandler_three_temperature_transfer_audit\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "[FAIL] " << error.what() << '\n';
        return 1;
    }
}
