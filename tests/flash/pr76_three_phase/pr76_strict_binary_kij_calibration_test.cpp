#include <mpmc/flash/pr76_split.hpp>
#include <mpmc/thermodynamics/pr76_phase.hpp>

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

struct BinaryDatum {
    double x_co2;
    double temperature_k;
    double pressure_pa;
};
struct HeavyComponent {
    std::string id;
    double critical_temperature_k;
    double critical_pressure_pa;
    double acentric_factor;
};
struct FitEvaluation {
    double kij{};
    double relative_sse{};
    double aard{};
    std::vector<double> predicted_pressures_pa;
};
struct BubbleResidual {
    double log_sum{};
    double vapor_co2{};
    double liquid_z{};
    double vapor_z{};
};
struct PressureBracket {
    double left_pa{};
    double right_pa{};
    double left_residual{};
    double right_residual{};
};

void require(bool condition, const char* message) {
    if (!condition) { throw std::runtime_error(message); }
}

th::Provenance ufc_source(std::string locator) {
    return {th::SourceKind::literature,
            "https://repositorio.ufc.br/handle/riufc/80725",
            "E. C. Q. Soria, UFC dissertation, 2025", std::move(locator),
            "Independent binary calibration input; the M-40 ternary transition is excluded",
            "Transcribed from the UFC repository dissertation; pressures converted MPa -> Pa",
            "Test-only numeric citation and provenance metadata"};
}
th::Provenance leal_source(std::string locator) {
    return {th::SourceKind::literature,
            "https://doi.org/10.1016/j.fluid.2018.08.019",
            "M. F. Leal et al., Fluid Phase Equilibria 478 (2018)", std::move(locator),
            "Independent CO2+n-hexadecane binary bubble-pressure calibration input",
            "Transcribed from Table 5 of the public author manuscript; pressures converted MPa -> Pa",
            "Test-only numeric citation and provenance metadata"};
}
th::Provenance trial_kij_source(const std::string& pair_id) {
    return {th::SourceKind::assumption, "strict-PR76-binary-kij-calibration",
            "binary-bubble-fit-v4", "optimization variable for " + pair_id,
            "Trial kij is fitted only to independent binary bubble-pressure data; it is not literature input",
            "Generated deterministically inside this test", "Test-only calibration variable"};
}
th::SourcedScalar scalar(double value, th::Unit unit, const th::Provenance& source,
                         std::string original_unit = "SI or dimensionless",
                         std::string conversion = "identity") {
    return {value, unit, source, std::move(original_unit), std::move(conversion)};
}

constexpr HeavyComponent decane{"n-decane", 617.7, 2.110e6, 0.4923};
constexpr HeavyComponent hexadecane{"n-hexadecane", 717.0, 1.489e6, 0.742};
constexpr double co2_tc_k = 304.21;
constexpr double co2_pc_pa = 7.383e6;
constexpr double co2_omega = 0.2236;

// UFC 2025 Table 7. The two highest-x near-critical observations are excluded
// from this local one-parameter fit; no ternary observation appears here.
constexpr std::array<BinaryDatum, 4> decane_323_data{{
    {0.3280, 323.08, 3.14e6},
    {0.4950, 322.96, 5.51e6},
    {0.7770, 323.21, 8.57e6},
    {0.9110, 323.01, 9.47e6},
}};

// Leal et al. 2018 Table 5. Reported uncertainties are U(p)=0.01 MPa,
// u(T)=0.1 K and u(x_CO2)=0.0003.
constexpr std::array<BinaryDatum, 3> hexadecane_313_data{{
    {0.3081, 313.2, 2.42e6},
    {0.5415, 313.2, 5.14e6},
    {0.7461, 313.2, 8.12e6},
}};
constexpr std::array<BinaryDatum, 3> hexadecane_333_data{{
    {0.3081, 333.2, 2.98e6},
    {0.5415, 333.2, 6.46e6},
    {0.7461, 333.2, 11.14e6},
}};

void require_data_provenance() {
    const auto c10 = ufc_source("Table 7: CO2+n-decane binary bubble pressures");
    const auto c16 = leal_source("Table 5: CO2+n-hexadecane binary bubble pressures");
    require(c10.kind == th::SourceKind::literature &&
                c10.reference.find("riufc/80725") != std::string::npos,
            "CO2+n-C10 calibration provenance is not the UFC binary dataset");
    require(c16.kind == th::SourceKind::literature &&
                c16.reference.find("10.1016/j.fluid.2018.08.019") != std::string::npos,
            "CO2+n-C16 calibration provenance is not the independent Leal dataset");
}

th::Pr76Phase<double> make_binary_model(const HeavyComponent& heavy, double kij) {
    const auto pure_source = ufc_source("Table 6: Tc, Pc and acentric factor");
    const auto kij_source = trial_kij_source("carbon-dioxide/" + heavy.id);
    std::vector<th::Component> catalog{
        {"carbon-dioxide", "carbon dioxide", th::ComponentKind::pure, pure_source, {}},
        {heavy.id, heavy.id, th::ComponentKind::pure, pure_source, {}}};
    const std::vector<std::string> order{"carbon-dioxide", heavy.id};

    th::PrParameterInput input;
    input.model_id = std::string(th::pr76_profile);
    input.dataset_id = "strict-PR76-independent-binary-calibration-" + heavy.id;
    input.revision = "ufc-table6-pure/binary-bubble-observations-v4";
    input.applicability = {std::nullopt, std::nullopt, pure_source};
    input.pure = {
        {"carbon-dioxide",
         scalar(co2_tc_k, th::Unit::kelvin, pure_source, "K", "identity"),
         scalar(co2_pc_pa, th::Unit::pascal, pure_source, "MPa", "MPa * 1e6 -> Pa"),
         scalar(co2_omega, th::Unit::dimensionless, pure_source)},
        {heavy.id,
         scalar(heavy.critical_temperature_k, th::Unit::kelvin, pure_source, "K", "identity"),
         scalar(heavy.critical_pressure_pa, th::Unit::pascal, pure_source,
                "MPa", "MPa * 1e6 -> Pa"),
         scalar(heavy.acentric_factor, th::Unit::dimensionless, pure_source)}};
    input.binary = {{"carbon-dioxide", heavy.id,
                     scalar(kij, th::Unit::dimensionless, kij_source)}};
    return th::Pr76Phase<double>::from_parameters(
        th::PrParameterSet::create(catalog, order, input));
}

// The experimental inputs are bubble pressures at known saturated-liquid
// compositions. Fit the incipient VLE branch directly; global TPD onset is a
// different observable when an LL instability exists.
double wilson_k(double pressure_pa, double temperature_k,
                double tc, double pc, double omega) {
    return (pc / pressure_pa) *
        std::exp(5.373 * (1.0 + omega) * (1.0 - tc / temperature_k));
}

std::optional<BubbleResidual> bubble_residual(
    double pressure_pa, const BinaryDatum& datum, const HeavyComponent& heavy,
    fl::Pr76VleEvaluator& evaluator) {
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
    const double k_heavy = wilson_k(pressure_pa, datum.temperature_k,
                                    heavy.critical_temperature_k,
                                    heavy.critical_pressure_pa,
                                    heavy.acentric_factor);
    const double raw0 = liquid_x[0] * k_co2;
    const double raw1 = liquid_x[1] * k_heavy;
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
            return BubbleResidual{final_log_sum, next_y[0], liquid.z, vapor.z};
        }
        vapor_y[0] = 0.5 * vapor_y[0] + 0.5 * next_y[0];
        vapor_y[1] = 1.0 - vapor_y[0];
    }
    return std::nullopt;
}

std::optional<PressureBracket> locate_bubble_bracket(
    const BinaryDatum& datum, const HeavyComponent& heavy,
    fl::Pr76VleEvaluator& evaluator) {
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
        const auto state = bubble_residual(pressure, datum, heavy, evaluator);
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
    const BinaryDatum& datum, const HeavyComponent& heavy,
    fl::Pr76VleEvaluator& evaluator) {
    const auto bracket = locate_bubble_bracket(datum, heavy, evaluator);
    if (!bracket) { return std::nullopt; }
    double left = bracket->left_pa;
    double right = bracket->right_pa;
    double f_left = bracket->left_residual;
    double f_right = bracket->right_residual;
    if (f_left == 0.0) { return left; }
    if (f_right == 0.0) { return right; }
    for (int iteration = 0; iteration < 22; ++iteration) {
        const double middle = 0.5 * (left + right);
        const auto state = bubble_residual(middle, datum, heavy, evaluator);
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

template <std::size_t N>
std::optional<FitEvaluation> evaluate_kij(
    const HeavyComponent& heavy, double kij, const std::array<BinaryDatum, N>& data) {
    auto model = make_binary_model(heavy, kij);
    fl::Pr76VleEvaluator evaluator(model);
    FitEvaluation evaluation;
    evaluation.kij = kij;
    evaluation.predicted_pressures_pa.reserve(N);
    double abs_relative_sum = 0.0;
    for (const auto& datum : data) {
        const auto predicted = bubble_pressure(datum, heavy, evaluator);
        if (!predicted || !std::isfinite(*predicted)) { return std::nullopt; }
        evaluation.predicted_pressures_pa.push_back(*predicted);
        const double relative = (*predicted - datum.pressure_pa) / datum.pressure_pa;
        evaluation.relative_sse += relative * relative;
        abs_relative_sum += std::abs(relative);
    }
    evaluation.aard = abs_relative_sum / static_cast<double>(N);
    return evaluation;
}

template <std::size_t N>
FitEvaluation fit_kij(const HeavyComponent& heavy, const std::array<BinaryDatum, N>& data) {
    constexpr double lower_bound = 0.0;
    constexpr double upper_bound = 0.16;
    constexpr int coarse_intervals = 16;
    const double coarse_step = (upper_bound - lower_bound) /
                               static_cast<double>(coarse_intervals);
    std::optional<FitEvaluation> best;
    int best_index = -1;
    for (int index = 0; index <= coarse_intervals; ++index) {
        const double kij = lower_bound + coarse_step * static_cast<double>(index);
        const auto candidate = evaluate_kij(heavy, kij, data);
        if (candidate && (!best || candidate->relative_sse < best->relative_sse)) {
            best = *candidate;
            best_index = index;
        }
    }
    if (!best || best_index < 0) {
        throw std::runtime_error("strict PR76 binary bubble fit found no valid coarse candidate");
    }
    double left = std::max(lower_bound,
        lower_bound + coarse_step * static_cast<double>(best_index - 1));
    double right = std::min(upper_bound,
        lower_bound + coarse_step * static_cast<double>(best_index + 1));
    constexpr double inverse_phi = 0.6180339887498948482;
    const auto cost = [&](double kij) {
        const auto value = evaluate_kij(heavy, kij, data);
        return value ? value->relative_sse : 1.0e12;
    };
    double c = right - inverse_phi * (right - left);
    double d = left + inverse_phi * (right - left);
    double fc = cost(c);
    double fd = cost(d);
    for (int iteration = 0; iteration < 14; ++iteration) {
        if (fc <= fd) {
            right = d; d = c; fd = fc;
            c = right - inverse_phi * (right - left); fc = cost(c);
        } else {
            left = c; c = d; fc = fd;
            d = left + inverse_phi * (right - left); fd = cost(d);
        }
    }
    const std::array<double, 6> final_candidates{
        best->kij, left, right, c, d, 0.5 * (left + right)};
    for (const double kij : final_candidates) {
        const auto candidate = evaluate_kij(heavy, kij, data);
        if (candidate && candidate->relative_sse < best->relative_sse) {
            best = *candidate;
        }
    }
    return *best;
}

template <std::size_t N>
void print_fit(std::string_view name, const FitEvaluation& fit,
               const std::array<BinaryDatum, N>& data) {
    std::cout << std::setprecision(10) << name << " kij=" << fit.kij
              << " AARD=" << fit.aard << " relative_SSE=" << fit.relative_sse << '\n';
    for (std::size_t i = 0; i < N; ++i) {
        std::cout << "  T=" << data[i].temperature_k << " xCO2=" << data[i].x_co2
                  << " Pexp_MPa=" << data[i].pressure_pa / 1.0e6
                  << " Pcalc_MPa=" << fit.predicted_pressures_pa[i] / 1.0e6 << '\n';
    }
}

void calibrate_independent_binary_kij() {
    require_data_provenance();
    const auto c10 = fit_kij(decane, decane_323_data);
    const auto c16_313 = fit_kij(hexadecane, hexadecane_313_data);
    const auto c16_333 = fit_kij(hexadecane, hexadecane_333_data);
    const auto c10_zero = evaluate_kij(decane, 0.0, decane_323_data);
    const auto c16_313_zero = evaluate_kij(hexadecane, 0.0, hexadecane_313_data);
    const auto c16_333_zero = evaluate_kij(hexadecane, 0.0, hexadecane_333_data);
    require(c10_zero && c16_313_zero && c16_333_zero,
            "strict PR76 zero-kij baseline could not reproduce the selected VLE branches");

    constexpr double target_temperature_k = 323.15;
    constexpr double low_temperature_k = 313.2;
    constexpr double high_temperature_k = 333.2;
    const double weight = (target_temperature_k - low_temperature_k) /
                          (high_temperature_k - low_temperature_k);
    const double c16_323 = c16_313.kij + weight * (c16_333.kij - c16_313.kij);

    print_fit("CO2+n-C10 local-323K", c10, decane_323_data);
    print_fit("CO2+n-C16 313.2K", c16_313, hexadecane_313_data);
    print_fit("CO2+n-C16 333.2K", c16_333, hexadecane_333_data);
    std::cout << std::setprecision(10)
              << "CO2+n-C16 interpolated-323.15K kij=" << c16_323 << '\n'
              << "zero-kij relative_SSE: C10=" << c10_zero->relative_sse
              << " C16_313=" << c16_313_zero->relative_sse
              << " C16_333=" << c16_333_zero->relative_sse << '\n';

    require(c10.kij > 0.0 && c10.kij < 0.16,
            "CO2+n-C10 strict-PR76 fitted kij left the search interval");
    require(c16_313.kij > 0.0 && c16_313.kij < 0.16 &&
                c16_333.kij > 0.0 && c16_333.kij < 0.16 &&
                c16_323 > 0.0 && c16_323 < 0.16,
            "CO2+n-C16 strict-PR76 fitted/interpolated kij left the search interval");
    require(c10.relative_sse < c10_zero->relative_sse &&
                c16_313.relative_sse < c16_313_zero->relative_sse &&
                c16_333.relative_sse < c16_333_zero->relative_sse,
            "strict PR76 fitted kij did not improve the independent binary pressure objective");
    require(std::isfinite(c10.aard) && std::isfinite(c16_313.aard) &&
                std::isfinite(c16_333.aard),
            "strict PR76 calibration produced nonfinite fit diagnostics");
    require(decane_323_data.size() == 4U && hexadecane_313_data.size() == 3U &&
                hexadecane_333_data.size() == 3U,
            "strict PR76 calibration dataset shape changed unexpectedly");
}

} // namespace

int main() {
    try {
        calibrate_independent_binary_kij();
        std::cout << "[PASS] strict_pr76_independent_binary_kij_calibration\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "[FAIL] " << error.what() << '\n';
        return 1;
    }
}
