#include <mpmc/flash/pr76_split.hpp>
#include <mpmc/thermodynamics/pr76_phase.hpp>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <iomanip>
#include <iostream>
#include <limits>
#include <optional>
#include <stdexcept>
#include <string>
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

struct BubbleResidual {
    double log_sum{};
};

struct PressureBracket {
    double left_pa{};
    double right_pa{};
    double left_residual{};
    double right_residual{};
};

struct PointFit {
    double kij{};
    double predicted_pressure_pa{};
    double relative_error{};
};

struct DatasetEvaluation {
    double aard{};
    double relative_sse{};
    std::array<double, 4> predicted_pressures_pa{};
};

struct LinearRegression {
    double intercept{};
    double slope{};
    double r_squared{};
};

constexpr double co2_tc_k = 304.21;
constexpr double co2_pc_pa = 7.383e6;
constexpr double co2_omega = 0.2236;
constexpr double decane_tc_k = 617.7;
constexpr double decane_pc_pa = 2.110e6;
constexpr double decane_omega = 0.4923;
constexpr double strict_kij_from_pr101 = 0.05226578047;

// UFC 2025 Table 7: the same four binary bubble-pressure observations used by PR #101.
// Their 0.25 K temperature span makes them a direct local test of whether ordinary
// kij(T) can explain the much larger composition-dependent pressure residual.
constexpr std::array<BinaryDatum, 4> decane_323_data{{
    {0.3280, 323.08, 3.14e6},
    {0.4950, 322.96, 5.51e6},
    {0.7770, 323.21, 8.57e6},
    {0.9110, 323.01, 9.47e6},
}};

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
            "Transcribed from the UFC repository dissertation; pressures converted MPa -> Pa",
            "Test-only effective-kij structure audit; no production model change"};
}

th::Provenance trial_kij_source() {
    return {th::SourceKind::assumption,
            "strict-PR76-C10-effective-kij-structure-audit",
            "pointwise binary bubble-pressure inversion",
            "CO2/n-decane kij",
            "Each trial kij is an audit variable inferred only from one UFC binary bubble-pressure datum",
            "Generated deterministically inside this test",
            "Diagnostic effective parameter; not a production or transferable model parameter"};
}

th::SourcedScalar scalar(double value,
                         th::Unit unit,
                         const th::Provenance& source,
                         std::string original_unit = "SI or dimensionless",
                         std::string conversion = "identity") {
    return {value, unit, source, std::move(original_unit), std::move(conversion)};
}

th::Pr76Phase<double> make_binary_model(double kij) {
    const auto pure_source = ufc_source("Table 6: CO2 and n-decane Tc, Pc and acentric factor");
    const auto kij_source = trial_kij_source();
    std::vector<th::Component> catalog{
        {"carbon-dioxide", "carbon dioxide", th::ComponentKind::pure, pure_source, {}},
        {"n-decane", "n-decane", th::ComponentKind::pure, pure_source, {}}};
    const std::vector<std::string> order{"carbon-dioxide", "n-decane"};

    th::PrParameterInput input;
    input.model_id = std::string(th::pr76_profile);
    input.dataset_id = "strict-PR76-C10-effective-kij-structure-audit";
    input.revision = "UFC-Table6-Table7/pointwise-kij-diagnostic-v1";
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
                     scalar(kij, th::Unit::dimensionless, kij_source)}};

    return th::Pr76Phase<double>::from_parameters(
        th::PrParameterSet::create(catalog, order, input));
}

double wilson_k(double pressure_pa,
                double temperature_k,
                double tc,
                double pc,
                double omega) {
    return (pc / pressure_pa) *
           std::exp(5.373 * (1.0 + omega) * (1.0 - tc / temperature_k));
}

std::optional<BubbleResidual> bubble_residual(
    double pressure_pa,
    const BinaryDatum& datum,
    fl::Pr76VleEvaluator& evaluator) {
    const std::array<double, 2> liquid_x{datum.x_co2, 1.0 - datum.x_co2};
    fl::PtSplitPhase liquid;
    try {
        liquid = evaluator(pressure_pa,
                           datum.temperature_k,
                           liquid_x,
                           fl::PtPhaseRole::liquid_candidate);
    } catch (const fl::StabilityPropertyError&) {
        return std::nullopt;
    }

    const double k_co2 = wilson_k(pressure_pa,
                                  datum.temperature_k,
                                  co2_tc_k,
                                  co2_pc_pa,
                                  co2_omega);
    const double k_decane = wilson_k(pressure_pa,
                                     datum.temperature_k,
                                     decane_tc_k,
                                     decane_pc_pa,
                                     decane_omega);
    const double raw0 = liquid_x[0] * k_co2;
    const double raw1 = liquid_x[1] * k_decane;
    if (!(raw0 > 0.0) || !(raw1 > 0.0) || !std::isfinite(raw0 + raw1)) {
        return std::nullopt;
    }
    std::array<double, 2> vapor_y{raw0 / (raw0 + raw1), raw1 / (raw0 + raw1)};

    for (int iteration = 0; iteration < 100; ++iteration) {
        fl::PtSplitPhase vapor;
        try {
            vapor = evaluator(pressure_pa,
                              datum.temperature_k,
                              vapor_y,
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
        if (!(scaled_sum > 0.0) || !std::isfinite(scaled_sum)) {
            return std::nullopt;
        }
        const double log_sum = largest + std::log(scaled_sum);
        const std::array<double, 2> next_y{
            std::exp(log_raw[0] - log_sum),
            std::exp(log_raw[1] - log_sum)};
        const double change = std::max(std::abs(next_y[0] - vapor_y[0]),
                                       std::abs(next_y[1] - vapor_y[1]));
        if (change <= 2.0e-12) {
            try {
                vapor = evaluator(pressure_pa,
                                  datum.temperature_k,
                                  next_y,
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
    const BinaryDatum& datum,
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
                best = PressureBracket{*previous_pressure,
                                       pressure,
                                       *previous_residual,
                                       state->log_sum};
                best_distance = distance;
            }
        }
        previous_pressure = pressure;
        previous_residual = state->log_sum;
    }
    return best;
}

std::optional<double> bubble_pressure(const BinaryDatum& datum,
                                      fl::Pr76VleEvaluator& evaluator) {
    const auto bracket = locate_bubble_bracket(datum, evaluator);
    if (!bracket) {
        return std::nullopt;
    }
    double left = bracket->left_pa;
    double right = bracket->right_pa;
    double f_left = bracket->left_residual;
    double f_right = bracket->right_residual;
    if (f_left == 0.0) {
        return left;
    }
    if (f_right == 0.0) {
        return right;
    }
    for (int iteration = 0; iteration < 22; ++iteration) {
        const double middle = 0.5 * (left + right);
        const auto state = bubble_residual(middle, datum, evaluator);
        if (!state) {
            return std::nullopt;
        }
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

std::optional<PointFit> evaluate_point(double kij, const BinaryDatum& datum) {
    auto model = make_binary_model(kij);
    fl::Pr76VleEvaluator evaluator(model);
    const auto predicted = bubble_pressure(datum, evaluator);
    if (!predicted || !std::isfinite(*predicted)) {
        return std::nullopt;
    }
    const double relative = (*predicted - datum.pressure_pa) / datum.pressure_pa;
    return PointFit{kij, *predicted, relative};
}

PointFit fit_point(const BinaryDatum& datum) {
    constexpr double lower_bound = 0.0;
    constexpr double upper_bound = 0.16;
    constexpr int coarse_intervals = 32;
    const double coarse_step = (upper_bound - lower_bound) /
                               static_cast<double>(coarse_intervals);

    std::optional<PointFit> best;
    int best_index = -1;
    const auto cost = [](const PointFit& fit) {
        return fit.relative_error * fit.relative_error;
    };
    for (int index = 0; index <= coarse_intervals; ++index) {
        const double kij = lower_bound + coarse_step * static_cast<double>(index);
        const auto candidate = evaluate_point(kij, datum);
        if (candidate && (!best || cost(*candidate) < cost(*best))) {
            best = *candidate;
            best_index = index;
        }
    }
    if (!best || best_index < 0) {
        throw std::runtime_error("C10 pointwise kij audit found no valid coarse candidate");
    }

    double left = std::max(lower_bound,
                           lower_bound + coarse_step * static_cast<double>(best_index - 1));
    double right = std::min(upper_bound,
                            lower_bound + coarse_step * static_cast<double>(best_index + 1));
    constexpr double inverse_phi = 0.6180339887498948482;
    const auto scalar_cost = [&](double kij) {
        const auto value = evaluate_point(kij, datum);
        return value ? cost(*value) : 1.0e12;
    };

    double c = right - inverse_phi * (right - left);
    double d = left + inverse_phi * (right - left);
    double fc = scalar_cost(c);
    double fd = scalar_cost(d);
    for (int iteration = 0; iteration < 18; ++iteration) {
        if (fc <= fd) {
            right = d;
            d = c;
            fd = fc;
            c = right - inverse_phi * (right - left);
            fc = scalar_cost(c);
        } else {
            left = c;
            c = d;
            fc = fd;
            d = left + inverse_phi * (right - left);
            fd = scalar_cost(d);
        }
    }

    const std::array<double, 6> final_candidates{
        best->kij, left, right, c, d, 0.5 * (left + right)};
    for (const double kij : final_candidates) {
        const auto candidate = evaluate_point(kij, datum);
        if (candidate && cost(*candidate) < cost(*best)) {
            best = *candidate;
        }
    }
    return *best;
}

DatasetEvaluation evaluate_kij_assignments(const std::array<double, 4>& kijs) {
    DatasetEvaluation evaluation;
    double abs_relative_sum = 0.0;
    for (std::size_t i = 0; i < decane_323_data.size(); ++i) {
        const auto point = evaluate_point(kijs[i], decane_323_data[i]);
        if (!point) {
            throw std::runtime_error("C10 kij assignment left the tracked incipient-VLE branch");
        }
        evaluation.predicted_pressures_pa[i] = point->predicted_pressure_pa;
        evaluation.relative_sse += point->relative_error * point->relative_error;
        abs_relative_sum += std::abs(point->relative_error);
    }
    evaluation.aard = abs_relative_sum /
                      static_cast<double>(decane_323_data.size());
    return evaluation;
}

LinearRegression regress(const std::array<double, 4>& coordinates,
                         const std::array<double, 4>& responses) {
    double x_mean = 0.0;
    double y_mean = 0.0;
    for (std::size_t i = 0; i < coordinates.size(); ++i) {
        x_mean += coordinates[i];
        y_mean += responses[i];
    }
    x_mean /= static_cast<double>(coordinates.size());
    y_mean /= static_cast<double>(responses.size());

    double sxx = 0.0;
    double sxy = 0.0;
    double syy = 0.0;
    for (std::size_t i = 0; i < coordinates.size(); ++i) {
        const double dx = coordinates[i] - x_mean;
        const double dy = responses[i] - y_mean;
        sxx += dx * dx;
        sxy += dx * dy;
        syy += dy * dy;
    }
    require(sxx > 0.0 && syy > 0.0,
            "C10 effective-kij regression lost coordinate or response variation");
    const double slope = sxy / sxx;
    const double intercept = y_mean - slope * x_mean;

    double residual_sum = 0.0;
    for (std::size_t i = 0; i < coordinates.size(); ++i) {
        const double residual = responses[i] - (intercept + slope * coordinates[i]);
        residual_sum += residual * residual;
    }
    return LinearRegression{intercept, slope, 1.0 - residual_sum / syy};
}

std::array<double, 4> project_kij(const LinearRegression& regression,
                                  const std::array<double, 4>& coordinates) {
    std::array<double, 4> result{};
    for (std::size_t i = 0; i < result.size(); ++i) {
        result[i] = regression.intercept + regression.slope * coordinates[i];
        require(result[i] >= 0.0 && result[i] <= 0.16 && std::isfinite(result[i]),
                "C10 projected effective kij left the audit interval");
    }
    return result;
}

void print_dataset(std::string_view name,
                   const std::array<double, 4>& kijs,
                   const DatasetEvaluation& evaluation) {
    std::cout << std::setprecision(12)
              << name << " AARD=" << evaluation.aard
              << " relative_SSE=" << evaluation.relative_sse << '\n';
    for (std::size_t i = 0; i < decane_323_data.size(); ++i) {
        std::cout << "  T=" << decane_323_data[i].temperature_k
                  << " xCO2=" << decane_323_data[i].x_co2
                  << " kij=" << kijs[i]
                  << " Pexp_MPa=" << decane_323_data[i].pressure_pa / 1.0e6
                  << " Pcalc_MPa=" << evaluation.predicted_pressures_pa[i] / 1.0e6
                  << '\n';
    }
}

void audit_c10_kij_structure() {
    const auto source = ufc_source(
        "Table 7: CO2+n-decane binary bubble pressures",
        "Same four near-323 K observations used in the independently merged PR101 strict-PR76 calibration");
    require(source.kind == th::SourceKind::literature &&
                source.reference.find("riufc/80725") != std::string::npos,
            "C10 kij-structure audit lost its UFC provenance");

    std::array<double, 4> temperatures{};
    std::array<double, 4> compositions{};
    std::array<double, 4> pointwise_kijs{};
    std::array<PointFit, 4> pointwise_fits{};
    for (std::size_t i = 0; i < decane_323_data.size(); ++i) {
        temperatures[i] = decane_323_data[i].temperature_k;
        compositions[i] = decane_323_data[i].x_co2;
        pointwise_fits[i] = fit_point(decane_323_data[i]);
        pointwise_kijs[i] = pointwise_fits[i].kij;
    }

    const auto [min_t, max_t] = std::minmax_element(temperatures.begin(), temperatures.end());
    const auto [min_x, max_x] = std::minmax_element(compositions.begin(), compositions.end());
    const auto [min_k, max_k] = std::minmax_element(pointwise_kijs.begin(), pointwise_kijs.end());
    const double temperature_span = *max_t - *min_t;
    const double composition_span = *max_x - *min_x;
    const double pointwise_kij_span = *max_k - *min_k;

    std::array<double, 4> constant_kijs{};
    constant_kijs.fill(strict_kij_from_pr101);
    const auto constant = evaluate_kij_assignments(constant_kijs);
    const auto temperature_regression = regress(temperatures, pointwise_kijs);
    const auto composition_regression = regress(compositions, pointwise_kijs);
    const auto temperature_kijs = project_kij(temperature_regression, temperatures);
    const auto composition_kijs = project_kij(composition_regression, compositions);
    const auto temperature_profile = evaluate_kij_assignments(temperature_kijs);
    const auto composition_profile = evaluate_kij_assignments(composition_kijs);
    const auto pointwise_profile = evaluate_kij_assignments(pointwise_kijs);

    std::cout << std::setprecision(12)
              << "C10 effective-kij structure audit: T_span_K=" << temperature_span
              << " xCO2_span=" << composition_span
              << " pointwise_kij_span=" << pointwise_kij_span << '\n';
    print_dataset("constant PR101 kij", constant_kijs, constant);
    print_dataset("pointwise effective kij", pointwise_kijs, pointwise_profile);
    std::cout << "linear kij(T): intercept=" << temperature_regression.intercept
              << " slope_per_K=" << temperature_regression.slope
              << " R2_on_pointwise_kij=" << temperature_regression.r_squared << '\n';
    print_dataset("projected linear kij(T)", temperature_kijs, temperature_profile);
    std::cout << "linear kij(xCO2): intercept=" << composition_regression.intercept
              << " slope_per_mole_fraction=" << composition_regression.slope
              << " R2_on_pointwise_kij=" << composition_regression.r_squared << '\n';
    print_dataset("projected linear kij(xCO2)", composition_kijs, composition_profile);

    require(temperature_span <= 0.30,
            "C10 structure audit is no longer a near-isothermal test");
    require(composition_span >= 0.50,
            "C10 structure audit lost the wide composition span needed for diagnosis");
    require(std::abs(constant.aard - 0.1324992938) <= 5.0e-6,
            "C10 structure audit no longer reproduces the merged PR101 constant-kij AARD");
    require(std::isfinite(temperature_regression.r_squared) &&
                std::isfinite(composition_regression.r_squared) &&
                std::isfinite(temperature_profile.aard) &&
                std::isfinite(composition_profile.aard),
            "C10 structure audit produced nonfinite diagnostics");
    require(pointwise_profile.relative_sse <= constant.relative_sse,
            "pointwise effective kij unexpectedly worsened the constant-kij objective");
}

} // namespace

int main() {
    try {
        audit_c10_kij_structure();
        std::cout << "[PASS] pr76_c10_effective_kij_structure_audit\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "[FAIL] " << error.what() << '\n';
        return 1;
    }
}
