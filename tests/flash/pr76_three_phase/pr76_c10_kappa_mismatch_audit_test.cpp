#include <mpmc/flash/pr76_split.hpp>
#include <mpmc/thermodynamics/pr76_phase.hpp>
#include <mpmc/thermodynamics/pr76_pure.hpp>

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
namespace fl = mpmc::flash;
namespace th = mpmc::thermodynamics;

enum class KappaMode { strict_pr76, soria_high_omega };

struct BinaryDatum {
    double x_co2;
    double temperature_k;
    double pressure_pa;
};
struct FitEvaluation {
    double kij{};
    double relative_sse{};
    double aard{};
    std::vector<double> predicted_pressures_pa;
};
struct BubbleResidual {
    double log_sum{};
};
struct PressureBracket {
    double left_pa{};
    double right_pa{};
    double left_residual{};
    double right_residual{};
};

constexpr double co2_tc_k = 304.21;
constexpr double co2_pc_pa = 7.383e6;
constexpr double co2_omega = 0.2236;
constexpr double decane_tc_k = 617.7;
constexpr double decane_pc_pa = 2.110e6;
constexpr double decane_omega = 0.4923;
constexpr double strict_kij_from_pr101 = 0.05226578047;
constexpr double soria_table12_kij = 0.09670;

constexpr std::array<BinaryDatum, 4> decane_323_data{{
    {0.3280, 323.08, 3.14e6},
    {0.4950, 322.96, 5.51e6},
    {0.7770, 323.21, 8.57e6},
    {0.9110, 323.01, 9.47e6},
}};

void require(bool condition, const char* message) {
    if (!condition) { throw std::runtime_error(message); }
}

th::Provenance ufc_source(std::string locator, std::string note = {}) {
    return {th::SourceKind::literature,
            "https://repositorio.ufc.br/handle/riufc/80725",
            "E. C. Q. Soria, UFC dissertation, 2025", std::move(locator),
            std::move(note),
            "Transcribed from the UFC repository dissertation",
            "Test-only mismatch-source audit; no production model change"};
}

th::Provenance assumption_source(std::string locator, std::string note) {
    return {th::SourceKind::assumption,
            "test-only/high-omega-kappa-emulation",
            "CO2+n-C10 kappa mismatch audit", std::move(locator), std::move(note),
            "Generated deterministically inside this audit",
            "The effective acentric factor is not physical input; it only injects the target kappa into existing PR76 code"};
}

th::SourcedScalar scalar(double value, th::Unit unit, const th::Provenance& source,
                         std::string original_unit = "SI or dimensionless",
                         std::string conversion = "identity") {
    return {value, unit, source, std::move(original_unit), std::move(conversion)};
}

constexpr double strict_kappa(double omega) noexcept {
    return 0.37464 + omega * (1.54226 - 0.26992 * omega);
}

constexpr double soria_kappa(double omega) noexcept {
    if (omega < 0.49) { return strict_kappa(omega); }
    return 0.379642 + omega *
        (1.48503 + omega * (-0.164423 + 0.01666 * omega));
}

double effective_omega_for_strict_kappa(double target_kappa) {
    double left = 0.0;
    double right = 1.0;
    require(target_kappa > strict_kappa(left) && target_kappa < strict_kappa(right),
            "target high-omega kappa lies outside the monotone audit inversion interval");
    for (int iteration = 0; iteration < 100; ++iteration) {
        const double middle = 0.5 * (left + right);
        if (strict_kappa(middle) < target_kappa) {
            left = middle;
        } else {
            right = middle;
        }
    }
    return 0.5 * (left + right);
}

double model_decane_omega(KappaMode mode) {
    if (mode == KappaMode::strict_pr76) { return decane_omega; }
    return effective_omega_for_strict_kappa(soria_kappa(decane_omega));
}

th::Pr76Phase<double> make_binary_model(double kij, KappaMode mode) {
    const auto pure_source = ufc_source(
        "Table 6: CO2 and n-decane Tc, Pc and acentric factor",
        "Physical Tc/Pc/omega are frozen to the dissertation values; the high-omega audit uses a documented test-only equivalent omega solely to inject Eq.21 kappa.");
    const auto kij_source = assumption_source(
        "CO2/n-decane kij trial",
        "The same scalar kij search is used for strict PR76 and the high-omega audit branch.");
    const double decane_model_omega = model_decane_omega(mode);
    const auto decane_omega_source = mode == KappaMode::strict_pr76
        ? pure_source
        : assumption_source(
            "Eq.21 high-omega kappa -> equivalent strict-PR76 omega",
            "Tc and Pc remain physical; only the internal PR76 kappa is matched to the dissertation Eq.21 value.");

    std::vector<th::Component> catalog{
        {"carbon-dioxide", "carbon dioxide", th::ComponentKind::pure, pure_source, {}},
        {"n-decane", "n-decane", th::ComponentKind::pure, pure_source, {}}};
    const std::vector<std::string> order{"carbon-dioxide", "n-decane"};

    th::PrParameterInput input;
    input.model_id = std::string(th::pr76_profile);
    input.dataset_id = mode == KappaMode::strict_pr76
        ? "audit-C10-strict-PR76"
        : "audit-C10-Soria-Eq21-kappa-emulation";
    input.revision = "UFC-Table6-Table7/same-bubble-objective/v1";
    input.applicability = {std::nullopt, std::nullopt, pure_source};
    input.pure = {
        {"carbon-dioxide",
         scalar(co2_tc_k, th::Unit::kelvin, pure_source, "K", "identity"),
         scalar(co2_pc_pa, th::Unit::pascal, pure_source, "MPa", "MPa * 1e6 -> Pa"),
         scalar(co2_omega, th::Unit::dimensionless, pure_source)},
        {"n-decane",
         scalar(decane_tc_k, th::Unit::kelvin, pure_source, "K", "identity"),
         scalar(decane_pc_pa, th::Unit::pascal, pure_source, "MPa", "MPa * 1e6 -> Pa"),
         scalar(decane_model_omega, th::Unit::dimensionless, decane_omega_source)}};
    input.binary = {{"carbon-dioxide", "n-decane",
                     scalar(kij, th::Unit::dimensionless, kij_source)}};

    auto model = th::Pr76Phase<double>::from_parameters(
        th::PrParameterSet::create(catalog, order, input));
    const auto decane_pure = th::Pr76Pure<double>::from_parameters(model.parameters(), 1U);
    const double expected_kappa = mode == KappaMode::strict_pr76
        ? strict_kappa(decane_omega) : soria_kappa(decane_omega);
    require(std::abs(decane_pure.kappa() - expected_kappa) <= 2.0e-14,
            "test-only kappa injection does not reproduce the requested alpha convention");
    return model;
}

double wilson_k(double pressure_pa, double temperature_k,
                double tc, double pc, double omega) {
    return (pc / pressure_pa) *
        std::exp(5.373 * (1.0 + omega) * (1.0 - tc / temperature_k));
}

std::optional<BubbleResidual> bubble_residual(
    double pressure_pa, const BinaryDatum& datum, fl::Pr76VleEvaluator& evaluator) {
    const std::array<double, 2> liquid_x{datum.x_co2, 1.0 - datum.x_co2};
    fl::PtSplitPhase liquid;
    try {
        liquid = evaluator(pressure_pa, datum.temperature_k, liquid_x,
                           fl::PtPhaseRole::liquid_candidate);
    } catch (const fl::StabilityPropertyError&) {
        return std::nullopt;
    }

    const double k_co2 = wilson_k(pressure_pa, datum.temperature_k,
                                  co2_tc_k, co2_pc_pa, co2_omega);
    const double k_decane = wilson_k(pressure_pa, datum.temperature_k,
                                     decane_tc_k, decane_pc_pa, decane_omega);
    const double raw0 = liquid_x[0] * k_co2;
    const double raw1 = liquid_x[1] * k_decane;
    if (!(raw0 > 0.0) || !(raw1 > 0.0) || !std::isfinite(raw0 + raw1)) {
        return std::nullopt;
    }
    std::array<double, 2> vapor_y{raw0 / (raw0 + raw1), raw1 / (raw0 + raw1)};

    for (int iteration = 0; iteration < 100; ++iteration) {
        fl::PtSplitPhase vapor;
        try {
            vapor = evaluator(pressure_pa, datum.temperature_k, vapor_y,
                              fl::PtPhaseRole::vapor_candidate);
        } catch (const fl::StabilityPropertyError&) {
            return std::nullopt;
        }
        const std::array<double, 2> log_raw{
            std::log(liquid_x[0]) + liquid.activity.ln_phi[0] - vapor.activity.ln_phi[0],
            std::log(liquid_x[1]) + liquid.activity.ln_phi[1] - vapor.activity.ln_phi[1]};
        const double largest = std::max(log_raw[0], log_raw[1]);
        const double scaled_sum = std::exp(log_raw[0] - largest) +
                                  std::exp(log_raw[1] - largest);
        if (!(scaled_sum > 0.0) || !std::isfinite(scaled_sum)) { return std::nullopt; }
        const double log_sum = largest + std::log(scaled_sum);
        const std::array<double, 2> next_y{
            std::exp(log_raw[0] - log_sum),
            std::exp(log_raw[1] - log_sum)};
        const double change = std::max(std::abs(next_y[0] - vapor_y[0]),
                                       std::abs(next_y[1] - vapor_y[1]));
        if (change <= 2.0e-12) {
            try {
                vapor = evaluator(pressure_pa, datum.temperature_k, next_y,
                                  fl::PtPhaseRole::vapor_candidate);
            } catch (const fl::StabilityPropertyError&) {
                return std::nullopt;
            }
            const std::array<double, 2> final_log_raw{
                std::log(liquid_x[0]) + liquid.activity.ln_phi[0] - vapor.activity.ln_phi[0],
                std::log(liquid_x[1]) + liquid.activity.ln_phi[1] - vapor.activity.ln_phi[1]};
            const double final_largest = std::max(final_log_raw[0], final_log_raw[1]);
            const double final_log_sum = final_largest + std::log(
                std::exp(final_log_raw[0] - final_largest) +
                std::exp(final_log_raw[1] - final_largest));
            const double relative_z = std::abs(vapor.z - liquid.z) /
                                      std::max(vapor.z, liquid.z);
            if (!(next_y[0] > liquid_x[0]) || relative_z <= 1.0e-7 ||
                !std::isfinite(final_log_sum)) {
                return std::nullopt;
            }
            return BubbleResidual{final_log_sum};
        }
        vapor_y[0] = 0.5 * vapor_y[0] + 0.5 * next_y[0];
        vapor_y[1] = 1.0 - vapor_y[0];
    }
    return std::nullopt;
}

std::optional<PressureBracket> locate_bubble_bracket(
    const BinaryDatum& datum, fl::Pr76VleEvaluator& evaluator) {
    constexpr int scan_intervals = 28;
    constexpr double lower_factor = 0.55;
    constexpr double upper_factor = 1.55;
    std::optional<double> previous_pressure;
    std::optional<double> previous_residual;
    std::optional<PressureBracket> best;
    double best_distance = datum.pressure_pa;
    for (int index = 0; index <= scan_intervals; ++index) {
        const double fraction = static_cast<double>(index) /
                                static_cast<double>(scan_intervals);
        const double pressure = datum.pressure_pa *
            (lower_factor + (upper_factor - lower_factor) * fraction);
        const auto state = bubble_residual(pressure, datum, evaluator);
        if (!state) {
            previous_pressure.reset();
            previous_residual.reset();
            continue;
        }
        if (previous_pressure && previous_residual &&
            ((*previous_residual <= 0.0 && state->log_sum >= 0.0) ||
             (*previous_residual >= 0.0 && state->log_sum <= 0.0))) {
            const double middle = 0.5 * (*previous_pressure + pressure);
            const double distance = std::abs(middle - datum.pressure_pa);
            if (!best || distance < best_distance) {
                best = PressureBracket{*previous_pressure, pressure,
                                       *previous_residual, state->log_sum};
                best_distance = distance;
            }
        }
        previous_pressure = pressure;
        previous_residual = state->log_sum;
    }
    return best;
}

std::optional<double> bubble_pressure(
    const BinaryDatum& datum, fl::Pr76VleEvaluator& evaluator) {
    const auto bracket = locate_bubble_bracket(datum, evaluator);
    if (!bracket) { return std::nullopt; }
    double left = bracket->left_pa;
    double right = bracket->right_pa;
    double f_left = bracket->left_residual;
    double f_right = bracket->right_residual;
    if (f_left == 0.0) { return left; }
    if (f_right == 0.0) { return right; }
    for (int iteration = 0; iteration < 22; ++iteration) {
        const double middle = 0.5 * (left + right);
        const auto state = bubble_residual(middle, datum, evaluator);
        if (!state) { return std::nullopt; }
        const double f_middle = state->log_sum;
        if ((f_left <= 0.0 && f_middle >= 0.0) ||
            (f_left >= 0.0 && f_middle <= 0.0)) {
            right = middle;
            f_right = f_middle;
        } else {
            left = middle;
            f_left = f_middle;
        }
    }
    (void)f_right;
    return 0.5 * (left + right);
}

std::optional<FitEvaluation> evaluate_kij(double kij, KappaMode mode) {
    auto model = make_binary_model(kij, mode);
    fl::Pr76VleEvaluator evaluator(model);
    FitEvaluation evaluation;
    evaluation.kij = kij;
    evaluation.predicted_pressures_pa.reserve(decane_323_data.size());
    double abs_relative_sum = 0.0;
    for (const auto& datum : decane_323_data) {
        const auto predicted = bubble_pressure(datum, evaluator);
        if (!predicted || !std::isfinite(*predicted)) { return std::nullopt; }
        evaluation.predicted_pressures_pa.push_back(*predicted);
        const double relative = (*predicted - datum.pressure_pa) / datum.pressure_pa;
        evaluation.relative_sse += relative * relative;
        abs_relative_sum += std::abs(relative);
    }
    evaluation.aard = abs_relative_sum /
                      static_cast<double>(decane_323_data.size());
    return evaluation;
}

FitEvaluation fit_kij(KappaMode mode) {
    constexpr double lower_bound = 0.0;
    constexpr double upper_bound = 0.16;
    constexpr int coarse_intervals = 16;
    const double coarse_step = (upper_bound - lower_bound) /
                               static_cast<double>(coarse_intervals);
    std::optional<FitEvaluation> best;
    int best_index = -1;
    for (int index = 0; index <= coarse_intervals; ++index) {
        const double kij = lower_bound + coarse_step * static_cast<double>(index);
        const auto candidate = evaluate_kij(kij, mode);
        if (candidate && (!best || candidate->relative_sse < best->relative_sse)) {
            best = *candidate;
            best_index = index;
        }
    }
    if (!best || best_index < 0) {
        throw std::runtime_error("C10 kappa audit found no valid coarse kij candidate");
    }
    double left = std::max(lower_bound,
        lower_bound + coarse_step * static_cast<double>(best_index - 1));
    double right = std::min(upper_bound,
        lower_bound + coarse_step * static_cast<double>(best_index + 1));
    constexpr double inverse_phi = 0.6180339887498948482;
    const auto cost = [mode](double kij) {
        const auto value = evaluate_kij(kij, mode);
        return value ? value->relative_sse : 1.0e12;
    };
    double c = right - inverse_phi * (right - left);
    double d = left + inverse_phi * (right - left);
    double fc = cost(c);
    double fd = cost(d);
    for (int iteration = 0; iteration < 14; ++iteration) {
        if (fc <= fd) {
            right = d;
            d = c;
            fd = fc;
            c = right - inverse_phi * (right - left);
            fc = cost(c);
        } else {
            left = c;
            c = d;
            fc = fd;
            d = left + inverse_phi * (right - left);
            fd = cost(d);
        }
    }
    const std::array<double, 6> final_candidates{
        best->kij, left, right, c, d, 0.5 * (left + right)};
    for (const double kij : final_candidates) {
        const auto candidate = evaluate_kij(kij, mode);
        if (candidate && candidate->relative_sse < best->relative_sse) {
            best = *candidate;
        }
    }
    return *best;
}

void print_fit(std::string_view name, const FitEvaluation& fit) {
    std::cout << std::setprecision(12) << name << " kij=" << fit.kij
              << " AARD=" << fit.aard
              << " relative_SSE=" << fit.relative_sse << '\n';
    for (std::size_t i = 0; i < decane_323_data.size(); ++i) {
        std::cout << "  T=" << decane_323_data[i].temperature_k
                  << " xCO2=" << decane_323_data[i].x_co2
                  << " Pexp_MPa=" << decane_323_data[i].pressure_pa / 1.0e6
                  << " Pcalc_MPa=" << fit.predicted_pressures_pa[i] / 1.0e6 << '\n';
    }
}

double alpha(double temperature_k, double tc_k, double kappa) {
    const double bracket = 1.0 + kappa * (1.0 - std::sqrt(temperature_k / tc_k));
    return bracket * bracket;
}

void audit_c10_kappa_mismatch_source() {
    const auto formula_source = ufc_source(
        "Section 4.4.1, Eqs.20-21: piecewise Peng-Robinson form factor m(omega)",
        "Eq.20 is the original quadratic for omega < 0.49; Eq.21 is the cubic branch for omega >= 0.49.");
    require(formula_source.kind == th::SourceKind::literature &&
                formula_source.reference.find("riufc/80725") != std::string::npos,
            "high-omega kappa audit lost its UFC literature provenance");
    require(decane_omega >= 0.49,
            "n-decane no longer falls on the dissertation high-omega branch");

    const double kappa_strict = strict_kappa(decane_omega);
    const double kappa_high = soria_kappa(decane_omega);
    const double equivalent_omega = model_decane_omega(KappaMode::soria_high_omega);
    const double alpha_strict_323 = alpha(323.15, decane_tc_k, kappa_strict);
    const double alpha_high_323 = alpha(323.15, decane_tc_k, kappa_high);

    const auto strict_fit = fit_kij(KappaMode::strict_pr76);
    const auto high_same_kij = evaluate_kij(strict_kij_from_pr101,
                                             KappaMode::soria_high_omega);
    const auto high_fit = fit_kij(KappaMode::soria_high_omega);
    const auto high_published_kij = evaluate_kij(soria_table12_kij,
                                                  KappaMode::soria_high_omega);
    require(high_same_kij && high_published_kij,
            "high-omega audit could not evaluate the frozen comparison kij values");

    std::cout << std::setprecision(12)
              << "C10 kappa audit: omega_physical=" << decane_omega
              << " kappa_strict=" << kappa_strict
              << " kappa_Soria_Eq21=" << kappa_high
              << " relative_kappa_shift=" << (kappa_high / kappa_strict - 1.0)
              << " equivalent_test_omega=" << equivalent_omega << '\n'
              << "alpha_323.15_strict=" << alpha_strict_323
              << " alpha_323.15_high_omega=" << alpha_high_323
              << " relative_alpha_shift=" << (alpha_high_323 / alpha_strict_323 - 1.0)
              << '\n';
    print_fit("strict-PR76 refit", strict_fit);
    print_fit("high-omega with PR101 kij", *high_same_kij);
    print_fit("high-omega refit", high_fit);
    print_fit("high-omega with Soria Table12 kij", *high_published_kij);
    std::cout << "AARD deltas vs strict refit: same_kij="
              << high_same_kij->aard - strict_fit.aard
              << " refit=" << high_fit.aard - strict_fit.aard << '\n';

    require(std::abs(strict_fit.kij - strict_kij_from_pr101) <= 5.0e-7,
            "audit no longer reproduces the merged PR101 strict C10 kij baseline");
    require(std::abs(strict_fit.aard - 0.1324992938) <= 5.0e-6,
            "audit no longer reproduces the merged PR101 strict C10 AARD baseline");
    require(std::isfinite(high_same_kij->aard) && std::isfinite(high_fit.aard) &&
                std::isfinite(high_published_kij->aard),
            "high-omega audit produced nonfinite diagnostics");
    require(high_fit.relative_sse <= high_same_kij->relative_sse,
            "high-omega refit is worse than keeping the independently fitted strict kij");
    require(decane_323_data.size() == 4U,
            "C10 mismatch audit dataset shape changed unexpectedly");
}

} // namespace

int main() {
    try {
        audit_c10_kappa_mismatch_source();
        std::cout << "[PASS] pr76_c10_high_omega_kappa_mismatch_audit\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "[FAIL] " << error.what() << '\n';
        return 1;
    }
}
