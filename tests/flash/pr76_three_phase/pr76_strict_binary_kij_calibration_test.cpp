#include <mpmc/flash/pr76_stability.hpp>
#include <mpmc/thermodynamics/pr76_phase.hpp>

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
            "bounded-1D-fit-v1", "optimization variable for " + pair_id,
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
    input.revision = "ufc-table6-pure/binary-observations-v1";
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

enum class HomogeneousStatus { stable, unstable };
std::optional<HomogeneousStatus> classify_feed(
    double pressure_pa, const BinaryDatum& datum, fl::Pr76StabilityEvaluator& evaluator) {
    const std::array<double, 2> feed{datum.x_co2, 1.0 - datum.x_co2};
    const auto result = fl::test_pr76_pt_stability(
        pressure_pa, datum.temperature_k, feed, evaluator);
    if (result.search.status == fl::StabilityStatus::unstable) {
        return HomogeneousStatus::unstable;
    }
    if (result.search.status == fl::StabilityStatus::no_instability_found) {
        return HomogeneousStatus::stable;
    }
    return std::nullopt;
}

std::optional<double> transition_pressure(
    const BinaryDatum& datum, fl::Pr76StabilityEvaluator& evaluator) {
    double lower = 0.82 * datum.pressure_pa;
    double upper = 1.22 * datum.pressure_pa;
    const auto lower_status = classify_feed(lower, datum, evaluator);
    const auto upper_status = classify_feed(upper, datum, evaluator);
    if (!lower_status || !upper_status ||
        *lower_status != HomogeneousStatus::unstable ||
        *upper_status != HomogeneousStatus::stable) {
        return std::nullopt;
    }
    for (int iteration = 0; iteration < 15; ++iteration) {
        const double middle = 0.5 * (lower + upper);
        const auto status = classify_feed(middle, datum, evaluator);
        if (!status) { return std::nullopt; }
        if (*status == HomogeneousStatus::unstable) {
            lower = middle;
        } else {
            upper = middle;
        }
    }
    return 0.5 * (lower + upper);
}

template <std::size_t N>
std::optional<FitEvaluation> evaluate_kij(
    const HeavyComponent& heavy, double kij, const std::array<BinaryDatum, N>& data) {
    auto model = make_binary_model(heavy, kij);
    fl::Pr76StabilityEvaluator evaluator(model);
    FitEvaluation evaluation;
    evaluation.kij = kij;
    evaluation.predicted_pressures_pa.reserve(N);
    double abs_relative_sum = 0.0;
    for (const auto& datum : data) {
        const auto predicted = transition_pressure(datum, evaluator);
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
    double left = 0.0;
    double right = 0.16;
    constexpr double inverse_phi = 0.6180339887498948482;
    const auto cost = [&](double kij) {
        const auto value = evaluate_kij(heavy, kij, data);
        return value ? value->relative_sse : 1.0e12;
    };
    double c = right - inverse_phi * (right - left);
    double d = left + inverse_phi * (right - left);
    double fc = cost(c);
    double fd = cost(d);
    for (int iteration = 0; iteration < 18; ++iteration) {
        if (fc <= fd) {
            right = d; d = c; fd = fc;
            c = right - inverse_phi * (right - left); fc = cost(c);
        } else {
            left = c; c = d; fc = fd;
            d = left + inverse_phi * (right - left); fd = cost(d);
        }
    }
    const double fitted = 0.5 * (left + right);
    const auto result = evaluate_kij(heavy, fitted, data);
    if (!result) {
        throw std::runtime_error("strict PR76 binary kij fit ended on an invalid phase bracket");
    }
    return *result;
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
              << "CO2+n-C16 interpolated-323.15K kij=" << c16_323 << '\n';

    require(c10.kij > 0.0 && c10.kij < 0.16,
            "CO2+n-C10 strict-PR76 fitted kij left the search interval");
    require(c16_313.kij > 0.0 && c16_313.kij < 0.16 &&
                c16_333.kij > 0.0 && c16_333.kij < 0.16 &&
                c16_323 > 0.0 && c16_323 < 0.16,
            "CO2+n-C16 strict-PR76 fitted/interpolated kij left the search interval");
    require(c10.aard < 0.08,
            "CO2+n-C10 strict-PR76 binary fit exceeds 8% pressure AARD");
    require(c16_313.aard < 0.08 && c16_333.aard < 0.08,
            "CO2+n-C16 strict-PR76 binary fit exceeds 8% pressure AARD");
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
