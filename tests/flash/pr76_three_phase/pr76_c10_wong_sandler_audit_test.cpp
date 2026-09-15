#include <mpmc/ad/dual.hpp>
#include <mpmc/flash/pr76_split.hpp>
#include <mpmc/thermodynamics/pr76_phase.hpp>
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
#include <utility>
#include <vector>

namespace {
namespace ad = mpmc::ad;
namespace fl = mpmc::flash;
namespace th = mpmc::thermodynamics;

constexpr double co2_tc_k = 304.21;
constexpr double co2_pc_pa = 7.383e6;
constexpr double co2_omega = 0.2236;
constexpr double decane_tc_k = 617.7;
constexpr double decane_pc_pa = 2.110e6;
constexpr double decane_omega = 0.4923;
constexpr double strict_kij_from_pr101 = 0.05226578047;

// Independent Wong-Sandler/NRTL parameters reported by Arenas-Quevedo et al.
// (Fluid Phase Equilibria 338 (2013) 30-36, Table 6) for CO2+n-decane.
// Their binary parameter regression used literature dataset [22]:
// Jimenez-Gallegos, Galicia-Luna, Elizalde-Solis, JCED 51 (2006) 1624-1628.
constexpr double ws_k12 = 0.7155;
constexpr double nrtl_delta12_j_per_mol = 11.8841e3;
constexpr double nrtl_delta21_j_per_mol = -1.9705e3;
constexpr double nrtl_nonrandomness = 0.3;

constexpr double target_temperature_k = 323.01;
constexpr double target_pressure_pa = 9.47e6;
constexpr std::array<double, 2> target_liquid_x{0.911, 0.089};

void require(bool condition, const char* message) {
    if (!condition) {
        throw std::runtime_error(message);
    }
}

th::Provenance ufc_source(std::string locator, std::string note = {}) {
    return {th::SourceKind::literature,
            "https://repositorio.ufc.br/handle/riufc/80725",
            "E. C. Q. Soria, UFC dissertation, 2025",
            std::move(locator),
            std::move(note),
            "Transcribed from the UFC repository dissertation",
            "Pure PR76 input and blind branch target only; not used to fit WS parameters"};
}

th::Provenance ws_source(std::string locator, std::string note = {}) {
    return {th::SourceKind::literature,
            "https://doi.org/10.1016/j.fluid.2012.10.012",
            "M. G. Arenas-Quevedo et al., Fluid Phase Equilibria 338 (2013) 30-36",
            std::move(locator),
            std::move(note),
            "WS/NRTL parameters transcribed from Table 6; equations from Eqs. (10)-(15)",
            "Independent binary parameter source; no UFC-2025 target datum used in regression"};
}

th::SourcedScalar scalar(double value,
                         th::Unit unit,
                         const th::Provenance& source,
                         std::string original_unit = "SI or dimensionless",
                         std::string conversion = "identity") {
    return {value, unit, source, std::move(original_unit), std::move(conversion)};
}

th::PrParameterSet make_parameter_set() {
    const auto pure_source = ufc_source("Table 6: CO2 and n-decane Tc, Pc and acentric factor");
    const auto kij_source = ufc_source(
        "PR101 strict-PR76 binary calibration result",
        "The production-classical A arm retains the independently calibrated scalar kij from merged PR #101.");

    std::vector<th::Component> catalog{
        {"carbon-dioxide", "carbon dioxide", th::ComponentKind::pure, pure_source, {}},
        {"n-decane", "n-decane", th::ComponentKind::pure, pure_source, {}}};
    const std::vector<std::string> order{"carbon-dioxide", "n-decane"};

    th::PrParameterInput input;
    input.model_id = std::string(th::pr76_profile);
    input.dataset_id = "UFC-2025-C10-pure-plus-PR101-kij";
    input.revision = "WS-mixing-rule-exact-bubble-audit/v3";
    input.applicability = {std::nullopt, std::nullopt, pure_source};
    input.pure = {
        {"carbon-dioxide",
         scalar(co2_tc_k, th::Unit::kelvin, pure_source, "K", "identity"),
         scalar(co2_pc_pa, th::Unit::pascal, pure_source, "MPa", "MPa * 1e6 -> Pa"),
         scalar(co2_omega, th::Unit::dimensionless, pure_source)},
        {"n-decane",
         scalar(decane_tc_k, th::Unit::kelvin, pure_source, "K", "identity"),
         scalar(decane_pc_pa, th::Unit::pascal, pure_source, "MPa", "MPa * 1e6 -> Pa"),
         scalar(decane_omega, th::Unit::dimensionless, pure_source)}};
    input.binary = {{"carbon-dioxide", "n-decane",
                     scalar(strict_kij_from_pr101, th::Unit::dimensionless, kij_source)}};
    return th::PrParameterSet::create(catalog, order, input);
}

enum class MixingMode { classical_vdw1f, wong_sandler_nrtl };
enum class RootRole { liquid, vapor };

template <typename Number>
struct MixedCoefficients {
    Number a;
    Number b;
};

struct GenericPhase {
    double z{};
    std::array<double, 2> ln_phi{};
    std::size_t root_count{};
};

struct BranchState {
    double log_sum{};
    std::array<double, 2> vapor_y{};
    double liquid_z{};
    double vapor_z{};
    double component_residual_gap{};
    bool distinct{};
};

struct YResidual {
    double y_co2{};
    double difference{};
    double common_residual{};
    GenericPhase vapor{};
};

struct YBracket {
    YResidual left{};
    YResidual right{};
};

struct PressureBracket {
    double left_pa{};
    double right_pa{};
    BranchState left_state{};
    BranchState right_state{};
};

struct BubblePressureResult {
    double pressure_pa{};
    double bracket_width_pa{};
    BranchState state{};
};

template <typename Number>
MixedCoefficients<Number> mix_from_moles(
    double temperature_k,
    const std::array<Number, 2>& mole_numbers,
    const th::Pr76Pure<double>& co2,
    const th::Pr76Pure<double>& decane,
    MixingMode mode) {
    const Number total = mole_numbers[0] + mole_numbers[1];
    const std::array<Number, 2> x{mole_numbers[0] / total, mole_numbers[1] / total};
    const auto co2_value = co2.evaluate(temperature_k);
    const auto decane_value = decane.evaluate(temperature_k);
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

std::optional<GenericPhase> evaluate_generic_phase(
    double pressure_pa,
    double temperature_k,
    const std::array<double, 2>& fractions,
    RootRole role,
    const th::Pr76Pure<double>& co2,
    const th::Pr76Pure<double>& decane,
    MixingMode mode) {
    const auto mixed = mix_from_moles<double>(temperature_k, fractions, co2, decane, mode);
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
    const auto mixed_dual = mix_from_moles<Dual>(temperature_k, moles, co2, decane, mode);
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
    result.root_count = roots.count;
    for (std::size_t i = 0; i < 2; ++i) {
        const double d_nb = nb.derivative(i);
        const double d_n2a_over_n = n2a.derivative(i) / n.value();
        const double bracket = d_n2a_over_n / mixed.a - d_nb / mixed.b;
        result.ln_phi[i] = first + (d_nb / mixed.b) * (z - 1.0) +
            (1.0 / (2.0 * sqrt2)) * (mixed.a / (mixed.b * rt)) *
            bracket * log_ratio;
        if (!std::isfinite(result.ln_phi[i])) {
            return std::nullopt;
        }
    }
    return result;
}

std::optional<YResidual> evaluate_y_residual(
    double pressure_pa,
    double temperature_k,
    const std::array<double, 2>& liquid_x,
    const GenericPhase& liquid,
    double y_co2,
    const th::Pr76Pure<double>& co2,
    const th::Pr76Pure<double>& decane,
    MixingMode mode) {
    if (!(y_co2 > liquid_x[0]) || !(y_co2 < 1.0)) {
        return std::nullopt;
    }
    const std::array<double, 2> vapor_y{y_co2, 1.0 - y_co2};
    const auto vapor = evaluate_generic_phase(
        pressure_pa, temperature_k, vapor_y, RootRole::vapor, co2, decane, mode);
    if (!vapor) {
        return std::nullopt;
    }
    const double residual0 = std::log(vapor_y[0]) + vapor->ln_phi[0] -
                             std::log(liquid_x[0]) - liquid.ln_phi[0];
    const double residual1 = std::log(vapor_y[1]) + vapor->ln_phi[1] -
                             std::log(liquid_x[1]) - liquid.ln_phi[1];
    if (!std::isfinite(residual0) || !std::isfinite(residual1)) {
        return std::nullopt;
    }
    return YResidual{y_co2,
                     residual0 - residual1,
                     0.5 * (residual0 + residual1),
                     *vapor};
}

bool brackets_zero(double a, double b) {
    return a == 0.0 || b == 0.0 || std::signbit(a) != std::signbit(b);
}

std::optional<YBracket> locate_y_bracket(
    double pressure_pa,
    double temperature_k,
    const std::array<double, 2>& liquid_x,
    const GenericPhase& liquid,
    const th::Pr76Pure<double>& co2,
    const th::Pr76Pure<double>& decane,
    MixingMode mode,
    std::optional<double> preferred_y) {
    constexpr int intervals = 360;
    const double lower = liquid_x[0] + 1.0e-5;
    const double upper = 1.0 - 1.0e-8;
    if (!(lower < upper)) {
        return std::nullopt;
    }

    std::optional<YResidual> previous;
    std::optional<YBracket> best;
    double best_metric = 0.0;
    for (int index = 0; index <= intervals; ++index) {
        const double fraction = static_cast<double>(index) / static_cast<double>(intervals);
        const double y = lower + (upper - lower) * fraction;
        const auto current = evaluate_y_residual(
            pressure_pa, temperature_k, liquid_x, liquid, y, co2, decane, mode);
        if (!current) {
            previous.reset();
            continue;
        }
        if (previous && brackets_zero(previous->difference, current->difference)) {
            const YBracket candidate{*previous, *current};
            const double midpoint = 0.5 * (candidate.left.y_co2 + candidate.right.y_co2);
            const double metric = preferred_y ? std::abs(midpoint - *preferred_y) : -midpoint;
            if (!best || metric < best_metric) {
                best = candidate;
                best_metric = metric;
            }
        }
        previous = *current;
    }
    return best;
}

std::optional<BranchState> incipient_vapor_branch(
    double pressure_pa,
    double temperature_k,
    const std::array<double, 2>& liquid_x,
    const th::Pr76Pure<double>& co2,
    const th::Pr76Pure<double>& decane,
    MixingMode mode,
    std::optional<double> preferred_y = std::nullopt) {
    const auto liquid = evaluate_generic_phase(
        pressure_pa, temperature_k, liquid_x, RootRole::liquid, co2, decane, mode);
    if (!liquid) {
        return std::nullopt;
    }
    const auto bracket = locate_y_bracket(
        pressure_pa, temperature_k, liquid_x, *liquid, co2, decane, mode, preferred_y);
    if (!bracket) {
        return std::nullopt;
    }

    YResidual left = bracket->left;
    YResidual right = bracket->right;
    for (int iteration = 0; iteration < 64 && right.y_co2 - left.y_co2 > 2.0e-13; ++iteration) {
        const double middle_y = 0.5 * (left.y_co2 + right.y_co2);
        const auto middle = evaluate_y_residual(
            pressure_pa, temperature_k, liquid_x, *liquid, middle_y, co2, decane, mode);
        if (!middle) {
            return std::nullopt;
        }
        if (brackets_zero(left.difference, middle->difference)) {
            right = *middle;
        } else {
            left = *middle;
        }
    }

    const double y_co2 = 0.5 * (left.y_co2 + right.y_co2);
    const auto final = evaluate_y_residual(
        pressure_pa, temperature_k, liquid_x, *liquid, y_co2, co2, decane, mode);
    if (!final) {
        return std::nullopt;
    }
    const double composition_gap = std::abs(y_co2 - liquid_x[0]);
    const double relative_z_gap = std::abs(final->vapor.z - liquid->z) /
                                  std::max(final->vapor.z, liquid->z);
    return BranchState{-final->common_residual,
                       {y_co2, 1.0 - y_co2},
                       liquid->z,
                       final->vapor.z,
                       std::abs(final->difference),
                       composition_gap > 1.0e-4 && relative_z_gap > 1.0e-5};
}

std::optional<PressureBracket> locate_ws_bubble_bracket(
    const th::Pr76Pure<double>& co2,
    const th::Pr76Pure<double>& decane) {
    const auto center = incipient_vapor_branch(target_pressure_pa,
                                               target_temperature_k,
                                               target_liquid_x,
                                               co2,
                                               decane,
                                               MixingMode::wong_sandler_nrtl);
    if (!center || !center->distinct || !std::isfinite(center->log_sum)) {
        return std::nullopt;
    }

    constexpr double step_pa = 10.0e3;
    constexpr int max_steps_each_direction = 250;
    std::optional<PressureBracket> best;
    double best_midpoint_distance = 0.0;

    for (const int direction : {-1, 1}) {
        double previous_pressure = target_pressure_pa;
        BranchState previous_state = *center;
        for (int step = 1; step <= max_steps_each_direction; ++step) {
            const double pressure = target_pressure_pa +
                                    static_cast<double>(direction * step) * step_pa;
            if (!(pressure > 0.0)) {
                break;
            }
            const auto current = incipient_vapor_branch(
                pressure,
                target_temperature_k,
                target_liquid_x,
                co2,
                decane,
                MixingMode::wong_sandler_nrtl,
                previous_state.vapor_y[0]);
            if (!current || !current->distinct || !std::isfinite(current->log_sum)) {
                break;
            }
            if (brackets_zero(previous_state.log_sum, current->log_sum)) {
                PressureBracket candidate;
                if (previous_pressure < pressure) {
                    candidate = {previous_pressure, pressure, previous_state, *current};
                } else {
                    candidate = {pressure, previous_pressure, *current, previous_state};
                }
                const double midpoint = 0.5 * (candidate.left_pa + candidate.right_pa);
                const double distance = std::abs(midpoint - target_pressure_pa);
                if (!best || distance < best_midpoint_distance) {
                    best = candidate;
                    best_midpoint_distance = distance;
                }
                break;
            }
            previous_pressure = pressure;
            previous_state = *current;
        }
    }
    return best;
}

std::optional<BubblePressureResult> solve_ws_bubble_pressure(
    const th::Pr76Pure<double>& co2,
    const th::Pr76Pure<double>& decane) {
    const auto bracket = locate_ws_bubble_bracket(co2, decane);
    if (!bracket) {
        return std::nullopt;
    }

    double left = bracket->left_pa;
    double right = bracket->right_pa;
    BranchState left_state = bracket->left_state;
    BranchState right_state = bracket->right_state;
    for (int iteration = 0; iteration < 64 && right - left > 0.05; ++iteration) {
        const double middle = 0.5 * (left + right);
        const double seed_y = 0.5 * (left_state.vapor_y[0] + right_state.vapor_y[0]);
        const auto middle_state = incipient_vapor_branch(
            middle,
            target_temperature_k,
            target_liquid_x,
            co2,
            decane,
            MixingMode::wong_sandler_nrtl,
            seed_y);
        if (!middle_state || !middle_state->distinct || !std::isfinite(middle_state->log_sum)) {
            return std::nullopt;
        }
        if (brackets_zero(left_state.log_sum, middle_state->log_sum)) {
            right = middle;
            right_state = *middle_state;
        } else {
            left = middle;
            left_state = *middle_state;
        }
    }

    const double pressure = 0.5 * (left + right);
    const double seed_y = 0.5 * (left_state.vapor_y[0] + right_state.vapor_y[0]);
    const auto state = incipient_vapor_branch(pressure,
                                              target_temperature_k,
                                              target_liquid_x,
                                              co2,
                                              decane,
                                              MixingMode::wong_sandler_nrtl,
                                              seed_y);
    if (!state || !state->distinct || !std::isfinite(state->log_sum)) {
        return std::nullopt;
    }
    return BubblePressureResult{pressure, right - left, *state};
}

void validate_generic_classical_fugacity(
    const th::Pr76Phase<double>& production,
    const th::Pr76Pure<double>& co2,
    const th::Pr76Pure<double>& decane) {
    fl::Pr76VleEvaluator production_evaluator(production);
    constexpr double pressure_pa = 8.0e6;
    constexpr double temperature_k = 323.01;
    const std::array<std::array<double, 2>, 2> compositions{{
        {0.911, 0.089},
        {0.995, 0.005}}};
    const std::array<fl::PtPhaseRole, 2> production_roles{
        fl::PtPhaseRole::liquid_candidate,
        fl::PtPhaseRole::vapor_candidate};
    const std::array<RootRole, 2> generic_roles{RootRole::liquid, RootRole::vapor};

    for (std::size_t state = 0; state < compositions.size(); ++state) {
        const auto production_phase = production_evaluator(
            pressure_pa, temperature_k, compositions[state], production_roles[state]);
        const auto generic_phase = evaluate_generic_phase(
            pressure_pa,
            temperature_k,
            compositions[state],
            generic_roles[state],
            co2,
            decane,
            MixingMode::classical_vdw1f);
        require(generic_phase.has_value(),
                "generic arbitrary-mixing fugacity path could not reproduce a classical phase");
        require(std::abs(generic_phase->z - production_phase.z) <= 2.0e-12,
                "generic arbitrary-mixing root differs from production classical PR76");
        for (std::size_t i = 0; i < 2; ++i) {
            require(std::abs(generic_phase->ln_phi[i] - production_phase.activity.ln_phi[i]) <=
                        2.0e-10,
                    "generic arbitrary-mixing fugacity differs from production classical PR76");
        }
    }
}

void audit_wong_sandler_exact_bubble_pressure() {
    const auto parameters = make_parameter_set();
    const auto production = th::Pr76Phase<double>::from_parameters(parameters);
    const auto co2 = th::Pr76Pure<double>::from_parameters(parameters, 0U);
    const auto decane = th::Pr76Pure<double>::from_parameters(parameters, 1U);

    validate_generic_classical_fugacity(production, co2, decane);

    const auto source = ws_source(
        "Eqs. (10)-(15), Table 6: PR + Wong-Sandler + NRTL for CO2+n-decane",
        "Table 6 reports k12=0.7155, delta12=11.8841 kJ/mol, delta21=-1.9705 kJ/mol; NRTL alpha=0.3 is fixed in the text. Binary source [22] is JCED 51 (2006) 1624-1628.");
    require(source.kind == th::SourceKind::literature &&
                source.reference.find("fluid.2012.10.012") != std::string::npos,
            "Wong-Sandler audit lost its literature provenance");

    const auto target_state = incipient_vapor_branch(target_pressure_pa,
                                                      target_temperature_k,
                                                      target_liquid_x,
                                                      co2,
                                                      decane,
                                                      MixingMode::wong_sandler_nrtl,
                                                      0.99254);
    require(target_state.has_value() && target_state->distinct,
            "direct binary fugacity solve no longer recovers the PR105 target branch");
    require(std::abs(target_state->component_residual_gap) <= 1.0e-10,
            "direct binary vapor-composition solve did not close component fugacity difference");

    const auto bubble = solve_ws_bubble_pressure(co2, decane);
    require(bubble.has_value(),
            "Wong-Sandler exact bubble-pressure solve found no continuous distinct-branch root");

    const double signed_error_pa = bubble->pressure_pa - target_pressure_pa;
    const double absolute_error_pa = std::abs(signed_error_pa);
    const double relative_error = absolute_error_pa / target_pressure_pa;

    std::cout << std::setprecision(14)
              << "CO2+n-C10 WS/NRTL exact bubble: T_K=" << target_temperature_k
              << " xCO2=" << target_liquid_x[0] << '\n'
              << "P_bubble_MPa=" << bubble->pressure_pa / 1.0e6
              << " P_exp_MPa=" << target_pressure_pa / 1.0e6
              << " signed_error_MPa=" << signed_error_pa / 1.0e6
              << " abs_error_MPa=" << absolute_error_pa / 1.0e6
              << " relative_error_pct=" << 100.0 * relative_error << '\n'
              << "yCO2=" << bubble->state.vapor_y[0]
              << " ZL=" << bubble->state.liquid_z
              << " ZV=" << bubble->state.vapor_z
              << " log_sum=" << bubble->state.log_sum
              << " sum_xK_minus_1=" << std::expm1(bubble->state.log_sum)
              << " component_residual_gap=" << bubble->state.component_residual_gap
              << " bracket_width_Pa=" << bubble->bracket_width_pa << '\n';

    require(bubble->pressure_pa > 0.0 && std::isfinite(bubble->pressure_pa),
            "Wong-Sandler bubble pressure is not finite and positive");
    require(bubble->bracket_width_pa <= 0.05,
            "Wong-Sandler bubble-pressure bracket did not converge tightly enough");
    require(std::abs(bubble->state.log_sum) <= 1.0e-8,
            "Wong-Sandler bubble-pressure root does not satisfy common fugacity closure");
    require(bubble->state.component_residual_gap <= 1.0e-10,
            "Wong-Sandler bubble vapor does not satisfy equal component fugacity residuals");
    require(bubble->state.vapor_y[0] > target_liquid_x[0],
            "Wong-Sandler bubble vapor is not CO2-richer than the specified liquid");
    require(bubble->state.vapor_z > bubble->state.liquid_z,
            "Wong-Sandler bubble root lost distinct liquid/vapor ordering");
}

} // namespace

int main() {
    try {
        audit_wong_sandler_exact_bubble_pressure();
        std::cout << "[PASS] pr76_c10_wong_sandler_exact_bubble_audit\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "[FAIL] " << error.what() << '\n';
        return 1;
    }
}
