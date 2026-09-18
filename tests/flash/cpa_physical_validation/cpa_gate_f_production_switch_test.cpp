#include <mpmc/thermodynamics/cpa_pt_phase.hpp>

#include "cpa_thermopack_parameter_snapshot.hpp"
#include "test_support.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <iomanip>
#include <iostream>
#include <limits>
#include <span>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace {
namespace th = mpmc::thermodynamics;
namespace external = cpa_thermopack_phase_kernel;

void require(bool condition, const std::string& message) {
    if (!condition) { throw std::runtime_error(message); }
}

void require_abs(
    double actual,
    double expected,
    double tolerance,
    const std::string& message) {
    const double delta = std::abs(actual - expected);
    if (!(delta <= tolerance)) {
        std::ostringstream stream;
        stream << std::setprecision(17)
               << message
               << " actual=" << actual
               << " expected=" << expected
               << " delta=" << delta
               << " tolerance=" << tolerance;
        throw std::runtime_error(stream.str());
    }
}

double roundoff(double scale) {
    return 4096.0 * std::numeric_limits<double>::epsilon() *
        std::max(1.0, std::abs(scale));
}

struct LegacyOracle {
    double pressure_physical_pa{};
    double pressure_association_pa{};
    double pressure_pa{};
    std::vector<double> mu_cubic_over_rt;
    std::vector<double> mu_association_over_rt;
    std::vector<double> mu_total_over_rt;
    std::vector<double> ln_phi;
};

LegacyOracle legacy_analytic_oracle(
    double target_pressure_pa,
    double temperature_k,
    double rho,
    std::span<const double> composition,
    const th::CpaParameterSet& parameters,
    const th::CpaAssociationResult& association) {
    require(composition.size() == parameters.size(),
            "Gate-F legacy oracle composition dimension mismatch");

    const double rt =
        th::cpa_gas_constant_j_per_mol_k * temperature_k;

    std::vector<double> pure_a(parameters.size(), 0.0);
    for (std::size_t i = 0U; i < parameters.size(); ++i) {
        const auto& pure = parameters.pure(i);
        const double tr =
            temperature_k / pure.critical_temperature_k;
        const double alpha_base =
            1.0 + pure.c1_dimensionless *
                (1.0 - std::sqrt(tr));
        pure_a[i] =
            pure.a0_pa_m6_per_mol2 *
            alpha_base * alpha_base;
    }

    double b_mix = 0.0;
    double a_mix = 0.0;
    std::vector<double> a_sums(parameters.size(), 0.0);
    for (std::size_t i = 0U; i < parameters.size(); ++i) {
        b_mix +=
            composition[i] *
            parameters.pure(i).b_m3_per_mol;
        for (std::size_t j = 0U; j < parameters.size(); ++j) {
            const double aij =
                std::sqrt(pure_a[i] * pure_a[j]) *
                (1.0 - parameters.kij(i, j));
            a_mix +=
                composition[i] * composition[j] * aij;
            a_sums[i] += composition[j] * aij;
        }
    }

    const double b_rho = b_mix * rho;
    require(
        std::isfinite(b_rho) &&
        b_rho > 0.0 && b_rho < 1.0,
        "Gate-F legacy oracle invalid reduced density");

    LegacyOracle oracle;
    oracle.pressure_physical_pa =
        rt * rho / (1.0 - b_rho) -
        a_mix * rho * rho / (1.0 + b_rho);

    std::vector<double> association_log_x(
        parameters.size(), 0.0);
    double association_sum = 0.0;
    for (const auto& site : association.sites) {
        require(
            site.unbonded_fraction > 0.0 &&
            std::isfinite(site.unbonded_fraction),
            "Gate-F legacy oracle invalid association site fraction");
        const double multiplicity =
            static_cast<double>(site.multiplicity);
        association_log_x[site.component_index] +=
            multiplicity *
            std::log(site.unbonded_fraction);
        association_sum +=
            composition[site.component_index] *
            multiplicity *
            (1.0 - site.unbonded_fraction);
    }

    oracle.pressure_association_pa =
        -0.5 * rt * rho *
        (1.0 + association.rho_dln_g_drho) *
        association_sum;
    oracle.pressure_pa =
        oracle.pressure_physical_pa +
        oracle.pressure_association_pa;

    const double z =
        target_pressure_pa / (rho * rt);
    require(
        std::isfinite(z) && z > 0.0,
        "Gate-F legacy oracle invalid target Z");

    const double z_physical =
        oracle.pressure_physical_pa / (rho * rt);
    const double a_over_brt =
        a_mix / (b_mix * rt);
    const double log_free_volume =
        std::log1p(-b_rho);
    const double log_attraction_volume =
        std::log1p(b_rho);
    const double g = association.radial_distribution;

    oracle.mu_cubic_over_rt.resize(parameters.size());
    oracle.mu_association_over_rt.resize(parameters.size());
    oracle.mu_total_over_rt.resize(parameters.size());
    oracle.ln_phi.resize(parameters.size());

    for (std::size_t i = 0U; i < parameters.size(); ++i) {
        const double b_i =
            parameters.pure(i).b_m3_per_mol;
        const double b_ratio = b_i / b_mix;
        const double attraction_ratio =
            2.0 * a_sums[i] / a_mix - b_ratio;

        const double mu_cubic =
            b_ratio * (z_physical - 1.0) -
            log_free_volume -
            a_over_brt * attraction_ratio *
                log_attraction_volume;
        const double mu_association =
            association_log_x[i] -
            (1.9 / 8.0) * rho * b_i * g *
                association_sum;
        const double total =
            mu_cubic + mu_association;

        oracle.mu_cubic_over_rt[i] = mu_cubic;
        oracle.mu_association_over_rt[i] =
            mu_association;
        oracle.mu_total_over_rt[i] = total;
        oracle.ln_phi[i] = total - std::log(z);
    }

    return oracle;
}

void check_state(
    double target_pressure_pa,
    double temperature_k,
    double rho,
    std::span<const double> composition,
    const th::CpaParameterSet& parameters,
    const std::string& label,
    double& max_pressure_delta,
    double& max_mu_delta,
    double& max_lnphi_delta) {
    const auto production =
        th::evaluate_cpa_phase_at_density(
            temperature_k,
            rho,
            composition,
            parameters);

    const auto legacy = legacy_analytic_oracle(
        target_pressure_pa,
        temperature_k,
        rho,
        composition,
        parameters,
        production.association);

    const double d_p_physical =
        std::abs(
            production.pressure_physical_pa -
            legacy.pressure_physical_pa);
    const double d_p_association =
        std::abs(
            production.pressure_association_pa -
            legacy.pressure_association_pa);
    const double d_p_total =
        std::abs(
            production.pressure_pa -
            legacy.pressure_pa);
    max_pressure_delta = std::max(
        {max_pressure_delta,
         d_p_physical,
         d_p_association,
         d_p_total});

    require_abs(
        production.pressure_physical_pa,
        legacy.pressure_physical_pa,
        roundoff(legacy.pressure_physical_pa),
        label + ": production/legacy physical pressure mismatch");
    require_abs(
        production.pressure_association_pa,
        legacy.pressure_association_pa,
        roundoff(legacy.pressure_association_pa),
        label + ": production/legacy association pressure mismatch");
    require_abs(
        production.pressure_pa,
        legacy.pressure_pa,
        roundoff(legacy.pressure_physical_pa) +
            roundoff(legacy.pressure_association_pa) +
            roundoff(legacy.pressure_pa),
        label + ": production/legacy total pressure mismatch");

    const auto chemical =
        th::cpa_helmholtz_residual_chemical_potentials(
            temperature_k,
            rho,
            composition,
            parameters,
            production.association);

    th::CpaPtRoot root;
    root.molar_density_mol_per_m3 = rho;
    th::cpa_detail::cpa_fill_ln_phi(
        target_pressure_pa,
        temperature_k,
        composition,
        parameters,
        production,
        root);

    require(
        root.ln_phi.size() == parameters.size(),
        label + ": production ln(phi) dimension mismatch");

    for (std::size_t i = 0U; i < parameters.size(); ++i) {
        const double d_mu_cubic =
            std::abs(
                chemical.cubic_over_rt[i] -
                legacy.mu_cubic_over_rt[i]);
        const double d_mu_association =
            std::abs(
                chemical.association_over_rt[i] -
                legacy.mu_association_over_rt[i]);
        const double d_mu_total =
            std::abs(
                chemical.total_over_rt[i] -
                legacy.mu_total_over_rt[i]);
        const double d_lnphi =
            std::abs(
                root.ln_phi[i] -
                legacy.ln_phi[i]);

        max_mu_delta = std::max(
            {max_mu_delta,
             d_mu_cubic,
             d_mu_association,
             d_mu_total});
        max_lnphi_delta =
            std::max(max_lnphi_delta, d_lnphi);

        require_abs(
            chemical.cubic_over_rt[i],
            legacy.mu_cubic_over_rt[i],
            roundoff(legacy.mu_cubic_over_rt[i]),
            label + ": production/legacy cubic mu mismatch");
        require_abs(
            chemical.association_over_rt[i],
            legacy.mu_association_over_rt[i],
            roundoff(legacy.mu_association_over_rt[i]),
            label + ": production/legacy association mu mismatch");
        require_abs(
            chemical.total_over_rt[i],
            legacy.mu_total_over_rt[i],
            roundoff(legacy.mu_cubic_over_rt[i]) +
                roundoff(legacy.mu_association_over_rt[i]) +
                roundoff(legacy.mu_total_over_rt[i]),
            label + ": production/legacy total mu mismatch");
        require_abs(
            root.ln_phi[i],
            legacy.ln_phi[i],
            roundoff(legacy.mu_cubic_over_rt[i]) +
                roundoff(legacy.mu_association_over_rt[i]) +
                roundoff(legacy.mu_total_over_rt[i]) +
                roundoff(std::log(
                    target_pressure_pa /
                    (rho *
                     th::cpa_gas_constant_j_per_mol_k *
                     temperature_k))),
            label + ": production/legacy ln(phi) mismatch");
    }
}

void run_gate_f_production_switch() {
    double max_pressure_delta = 0.0;
    double max_mu_delta = 0.0;
    double max_lnphi_delta = 0.0;
    std::size_t comparisons = 0U;

    for (const bool swapped : {false, true}) {
        const auto literature =
            cpa_physical_test::parameters(swapped);
        const auto parity =
            cpa_thermopack_snapshot::parameters(swapped);

        for (std::size_t state_index = 0U;
             state_index < external::states.size();
             ++state_index) {
            const auto& state = external::states[state_index];
            for (const auto* phase :
                 {&state.liquid, &state.vapor}) {
                const auto composition =
                    cpa_physical_test::composition(
                        phase->composition_methanol,
                        swapped);
                const std::string suffix =
                    " state=" +
                    std::to_string(state_index) +
                    (phase == &state.liquid
                         ? " phase=liquid"
                         : " phase=vapor") +
                    (swapped
                         ? " order=swapped"
                         : " order=normal");

                check_state(
                    state.pressure_pa,
                    state.temperature_k,
                    phase->molar_density_mol_per_m3,
                    composition,
                    literature,
                    "literature" + suffix,
                    max_pressure_delta,
                    max_mu_delta,
                    max_lnphi_delta);
                ++comparisons;

                check_state(
                    state.pressure_pa,
                    state.temperature_k,
                    phase->molar_density_mol_per_m3,
                    composition,
                    parity,
                    "parity" + suffix,
                    max_pressure_delta,
                    max_mu_delta,
                    max_lnphi_delta);
                ++comparisons;
            }
        }
    }

    require(
        comparisons == 40U,
        "Gate-F production/legacy comparison count changed");

    std::cout
        << std::setprecision(17)
        << "CPA_HELMHOLTZ_GATE_F_PRODUCTION_SWITCH_OK"
        << " comparisons=" << comparisons
        << " snapshots=literature,thermopack-parity"
        << " component_orders=normal,swapped"
        << " max_abs_dP_pa=" << max_pressure_delta
        << " max_abs_dmu_over_rt=" << max_mu_delta
        << " max_abs_dlnphi=" << max_lnphi_delta
        << '\n';
}

} // namespace

int main() {
    try {
        run_gate_f_production_switch();
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "[FAIL] "
                  << error.what() << '\n';
        return 1;
    }
}
