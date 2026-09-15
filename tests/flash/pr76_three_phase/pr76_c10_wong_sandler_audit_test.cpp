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
    input.revision = "WS-mixing-rule-A-B-audit/v1";
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
    bool distinct{};
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

double wilson_k(double pressure_pa,
                double temperature_k,
                double tc,
                double pc,
                double omega) {
    return (pc / pressure_pa) *
           std::exp(5.373 * (1.0 + omega) * (1.0 - tc / temperature_k));
}

std::optional<BranchState> incipient_vapor_branch(
    double pressure_pa,
    double temperature_k,
    const std::array<double, 2>& liquid_x,
    const th::Pr76Pure<double>& co2,
    const th::Pr76Pure<double>& decane,
    MixingMode mode) {
    const auto liquid = evaluate_generic_phase(
        pressure_pa, temperature_k, liquid_x, RootRole::liquid, co2, decane, mode);
    if (!liquid) {
        return std::nullopt;
    }

    const std::array<double, 2> initial_raw{
        liquid_x[0] * wilson_k(pressure_pa, temperature_k, co2_tc_k, co2_pc_pa, co2_omega),
        liquid_x[1] * wilson_k(pressure_pa, temperature_k, decane_tc_k, decane_pc_pa, decane_omega)};
    const double initial_sum = initial_raw[0] + initial_raw[1];
    if (!(initial_sum > 0.0) || !std::isfinite(initial_sum)) {
        return std::nullopt;
    }
    std::array<double, 2> vapor_y{
        initial_raw[0] / initial_sum, initial_raw[1] / initial_sum};

    for (int iteration = 0; iteration < 160; ++iteration) {
        const auto vapor = evaluate_generic_phase(
            pressure_pa, temperature_k, vapor_y, RootRole::vapor, co2, decane, mode);
        if (!vapor) {
            return std::nullopt;
        }

        const std::array<double, 2> log_raw{
            std::log(liquid_x[0]) + liquid->ln_phi[0] - vapor->ln_phi[0],
            std::log(liquid_x[1]) + liquid->ln_phi[1] - vapor->ln_phi[1]};
        const double largest = std::max(log_raw[0], log_raw[1]);
        const double scaled_sum = std::exp(log_raw[0] - largest) +
                                  std::exp(log_raw[1] - largest);
        if (!(scaled_sum > 0.0) || !std::isfinite(scaled_sum)) {
            return std::nullopt;
        }
        const double log_sum = largest + std::log(scaled_sum);
        const std::array<double, 2> next_y{
            std::exp(log_raw[0] - log_sum), std::exp(log_raw[1] - log_sum)};
        const double change = std::max(std::abs(next_y[0] - vapor_y[0]),
                                       std::abs(next_y[1] - vapor_y[1]));
        if (change <= 2.0e-12) {
            const auto final_vapor = evaluate_generic_phase(
                pressure_pa, temperature_k, next_y, RootRole::vapor, co2, decane, mode);
            if (!final_vapor) {
                return std::nullopt;
            }
            const double composition_gap = std::max(std::abs(next_y[0] - liquid_x[0]),
                                                    std::abs(next_y[1] - liquid_x[1]));
            const double relative_z_gap = std::abs(final_vapor->z - liquid->z) /
                                          std::max(final_vapor->z, liquid->z);
            return BranchState{log_sum,
                               next_y,
                               liquid->z,
                               final_vapor->z,
                               composition_gap > 1.0e-4 && relative_z_gap > 1.0e-5};
        }
        vapor_y[0] = 0.5 * vapor_y[0] + 0.5 * next_y[0];
        vapor_y[1] = 1.0 - vapor_y[0];
    }
    return std::nullopt;
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

void audit_wong_sandler_branch_recovery() {
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

    const auto classical = incipient_vapor_branch(target_pressure_pa,
                                                   target_temperature_k,
                                                   target_liquid_x,
                                                   co2,
                                                   decane,
                                                   MixingMode::classical_vdw1f);
    const auto ws = incipient_vapor_branch(target_pressure_pa,
                                            target_temperature_k,
                                            target_liquid_x,
                                            co2,
                                            decane,
                                            MixingMode::wong_sandler_nrtl);
    require(classical.has_value(),
            "classical A arm could not resolve the target state for topology comparison");
    require(ws.has_value(),
            "Wong-Sandler B arm could not resolve the target state");

    std::cout << std::setprecision(12)
              << "CO2+n-C10 target: T_K=" << target_temperature_k
              << " xCO2=" << target_liquid_x[0]
              << " P_MPa=" << target_pressure_pa / 1.0e6 << '\n'
              << "classical-vdW1f: distinct=" << classical->distinct
              << " yCO2=" << classical->vapor_y[0]
              << " ZL=" << classical->liquid_z
              << " ZV=" << classical->vapor_z
              << " log_sum=" << classical->log_sum << '\n'
              << "Wong-Sandler/NRTL independent params: k12=" << ws_k12
              << " delta12_kJmol=" << nrtl_delta12_j_per_mol / 1.0e3
              << " delta21_kJmol=" << nrtl_delta21_j_per_mol / 1.0e3
              << " alpha=" << nrtl_nonrandomness << '\n'
              << "Wong-Sandler: distinct=" << ws->distinct
              << " yCO2=" << ws->vapor_y[0]
              << " ZL=" << ws->liquid_z
              << " ZV=" << ws->vapor_z
              << " log_sum=" << ws->log_sum
              << " exp(log_sum)-1=" << std::expm1(ws->log_sum) << '\n';

    require(!classical->distinct,
            "classical vdW1f unexpectedly retained a distinct target-pressure VLE branch");
    require(ws->distinct,
            "Wong-Sandler did not recover a distinct VLE branch at the 9.47 MPa target");
    require(ws->vapor_y[0] > target_liquid_x[0],
            "recovered Wong-Sandler vapor branch is not CO2-richer than the liquid");
    require(std::abs(ws->log_sum) < 1.0e-2,
            "recovered Wong-Sandler branch is not close to incipient fugacity closure");
}

} // namespace

int main() {
    try {
        audit_wong_sandler_branch_recovery();
        std::cout << "[PASS] pr76_c10_wong_sandler_branch_audit\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "[FAIL] " << error.what() << '\n';
        return 1;
    }
}
