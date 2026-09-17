#include <mpmc/ad/math.hpp>
#include <mpmc/ad/runtime_differentiate.hpp>
#include <mpmc/thermodynamics/cpa_pt_phase.hpp>
#include <mpmc/thermodynamics/cpa_residual_helmholtz.hpp>

#include "cpa_thermopack_parameter_snapshot.hpp"
#include "cpa_thermopack_parity_thresholds.hpp"
#include "test_support.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <iomanip>
#include <iostream>
#include <limits>
#include <numeric>
#include <span>
#include <stdexcept>
#include <string>
#include <string_view>
#include <vector>

namespace {
namespace ad = mpmc::ad;
namespace th = mpmc::thermodynamics;
namespace external = cpa_thermopack_phase_kernel;
namespace thresholds = cpa_thermopack_parity_thresholds;

void require(bool condition, const std::string& message) {
    if (!condition) { throw std::runtime_error(message); }
}

void require_close(double actual, double expected,
                   double absolute_tolerance,
                   double relative_tolerance,
                   const std::string& message) {
    const double scale = std::max(std::abs(actual), std::abs(expected));
    const double tolerance = absolute_tolerance + relative_tolerance * scale;
    if (!(std::abs(actual - expected) <= tolerance)) {
        throw std::runtime_error(
            message + " actual=" + std::to_string(actual) +
            " expected=" + std::to_string(expected) +
            " tolerance=" + std::to_string(tolerance));
    }
}

double roundoff_tolerance(double scale) {
    constexpr double multiplier = 4096.0;
    return multiplier * std::numeric_limits<double>::epsilon() *
        std::max(1.0, std::abs(scale));
}

struct HelmholtzDerivatives {
    double cubic_value{};
    double association_q_value{};
    double total_value{};
    double pressure_physical_pa{};
    double pressure_association_pa{};
    double pressure_total_pa{};
    std::vector<double> mu_cubic_over_rt;
    std::vector<double> mu_association_over_rt;
    std::vector<double> mu_total_over_rt;
    std::vector<double> ln_phi;
};

HelmholtzDerivatives derive_from_helmholtz(
    double target_pressure_pa,
    double temperature_k,
    double molar_density_mol_per_m3,
    std::span<const double> composition,
    const th::CpaParameterSet& parameters,
    const th::CpaAssociationResult& association) {
    require(composition.size() == parameters.size(),
            "Helmholtz AD adapter composition dimension mismatch");
    require(target_pressure_pa > 0.0 && std::isfinite(target_pressure_pa),
            "Helmholtz AD adapter requires positive target pressure");

    std::vector<double> inputs;
    inputs.reserve(composition.size() + 1U);
    inputs.push_back(1.0 / molar_density_mol_per_m3);
    inputs.insert(inputs.end(), composition.begin(), composition.end());

    using Number = ad::Dual<double, 4>;
    ad::RuntimeJacobianWorkspace<double, 4> workspace;
    const auto callback = [&](std::span<const Number> variables,
                              std::span<Number> outputs) {
        const Number temperature{temperature_k};
        const Number& volume = variables.front();
        const auto mole_numbers = variables.subspan(1U);
        outputs[0] = th::cpa_cubic_residual_helmholtz_reduced(
            temperature, volume, mole_numbers, parameters);
        outputs[1] = th::cpa_association_q_reduced(
            temperature, volume, mole_numbers, parameters, association);
        outputs[2] = th::cpa_residual_helmholtz_reduced(
            temperature, volume, mole_numbers, parameters, association);
    };

    const auto differentiated = ad::value_and_jacobian_runtime<4>(
        callback, std::span<const double>{inputs}, 3U, workspace,
        {1024U, 8U, 8192U});
    require(differentiated.input_count == inputs.size() &&
                differentiated.output_count == 3U,
            "Helmholtz AD adapter returned unexpected Jacobian shape");

    HelmholtzDerivatives result;
    result.cubic_value = differentiated.values[0];
    result.association_q_value = differentiated.values[1];
    result.total_value = differentiated.values[2];

    const double total_moles = std::accumulate(
        composition.begin(), composition.end(), 0.0);
    const double volume = inputs.front();
    const double rt = th::cpa_gas_constant_j_per_mol_k * temperature_k;
    const double ideal_pressure = total_moles * rt / volume;
    const std::size_t input_count = differentiated.input_count;
    const auto derivative = [&](std::size_t output, std::size_t input) {
        return differentiated.jacobian[output * input_count + input];
    };

    result.pressure_physical_pa = ideal_pressure - rt * derivative(0U, 0U);
    result.pressure_association_pa = -rt * derivative(1U, 0U);
    result.pressure_total_pa = ideal_pressure - rt * derivative(2U, 0U);

    result.mu_cubic_over_rt.resize(parameters.size());
    result.mu_association_over_rt.resize(parameters.size());
    result.mu_total_over_rt.resize(parameters.size());
    result.ln_phi.resize(parameters.size());

    const double z = target_pressure_pa /
        (molar_density_mol_per_m3 * rt);
    require(z > 0.0 && std::isfinite(z),
            "Helmholtz AD adapter produced invalid target Z");
    const double log_z = std::log(z);
    for (std::size_t i = 0U; i < parameters.size(); ++i) {
        result.mu_cubic_over_rt[i] = derivative(0U, i + 1U);
        result.mu_association_over_rt[i] = derivative(1U, i + 1U);
        result.mu_total_over_rt[i] = derivative(2U, i + 1U);
        result.ln_phi[i] = result.mu_total_over_rt[i] - log_z;
    }
    return result;
}

double direct_association_helmholtz_value(
    std::span<const double> mole_numbers,
    const th::CpaAssociationResult& association) {
    double value = 0.0;
    for (const auto& site : association.sites) {
        const double x_site = site.unbonded_fraction;
        value += mole_numbers[site.component_index] *
            static_cast<double>(site.multiplicity) *
            (std::log(x_site) - 0.5 * x_site + 0.5);
    }
    return value;
}

std::vector<double> analytic_association_mu(
    double molar_density_mol_per_m3,
    std::span<const double> composition,
    const th::CpaParameterSet& parameters,
    const th::CpaAssociationResult& association) {
    std::vector<double> association_log_x(parameters.size(), 0.0);
    double association_sum = 0.0;
    for (const auto& site : association.sites) {
        const double multiplicity = static_cast<double>(site.multiplicity);
        association_log_x[site.component_index] +=
            multiplicity * std::log(site.unbonded_fraction);
        association_sum += composition[site.component_index] * multiplicity *
            (1.0 - site.unbonded_fraction);
    }

    std::vector<double> result(parameters.size(), 0.0);
    for (std::size_t i = 0U; i < parameters.size(); ++i) {
        result[i] = association_log_x[i] -
            (1.9 / 8.0) * molar_density_mol_per_m3 *
            parameters.pure(i).b_m3_per_mol *
            association.radial_distribution * association_sum;
    }
    return result;
}

void check_internal_state(
    double temperature_k,
    double molar_density_mol_per_m3,
    std::span<const double> composition,
    const th::CpaParameterSet& parameters,
    std::string_view label,
    double& max_q_value_delta,
    double& max_pressure_delta_pa,
    double& max_mu_delta,
    double& max_ln_phi_delta) {
    const auto state = th::evaluate_cpa_phase_at_density(
        temperature_k, molar_density_mol_per_m3,
        composition, parameters);
    require(state.pressure_pa > 0.0,
            std::string(label) + ": internal state pressure must be positive");

    th::CpaPtRoot analytic_root;
    analytic_root.molar_density_mol_per_m3 = molar_density_mol_per_m3;
    th::cpa_detail::cpa_fill_ln_phi(
        state.pressure_pa, temperature_k, composition,
        parameters, state, analytic_root);

    const auto derived = derive_from_helmholtz(
        state.pressure_pa, temperature_k, molar_density_mol_per_m3,
        composition, parameters, state.association);
    const double direct_assoc = direct_association_helmholtz_value(
        composition, state.association);
    const auto analytic_assoc_mu = analytic_association_mu(
        molar_density_mol_per_m3, composition, parameters, state.association);
    const double rt = th::cpa_gas_constant_j_per_mol_k * temperature_k;
    const double z = state.pressure_pa /
        (molar_density_mol_per_m3 * rt);

    const double q_delta = std::abs(derived.association_q_value - direct_assoc);
    max_q_value_delta = std::max(max_q_value_delta, q_delta);
    require_close(
        derived.association_q_value, direct_assoc,
        1.0e-10, 1.0e-12,
        std::string(label) + ": stationary Q did not reproduce association Helmholtz value");
    require_close(
        derived.total_value,
        derived.cubic_value + derived.association_q_value,
        roundoff_tolerance(derived.total_value), 0.0,
        std::string(label) + ": total Helmholtz value did not equal cubic+association");

    const std::array<std::pair<double, double>, 3> pressures{{
        {derived.pressure_physical_pa, state.pressure_physical_pa},
        {derived.pressure_association_pa, state.pressure_association_pa},
        {derived.pressure_total_pa, state.pressure_pa}}};
    for (std::size_t i = 0U; i < pressures.size(); ++i) {
        const double delta = std::abs(pressures[i].first - pressures[i].second);
        max_pressure_delta_pa = std::max(max_pressure_delta_pa, delta);
        require_close(
            pressures[i].first, pressures[i].second,
            roundoff_tolerance(pressures[i].second), 0.0,
            std::string(label) + ": AD pressure decomposition mismatch index=" +
                std::to_string(i));
    }

    for (std::size_t i = 0U; i < parameters.size(); ++i) {
        const double analytic_total_mu = analytic_root.ln_phi[i] + std::log(z);
        const double analytic_cubic_mu = analytic_total_mu - analytic_assoc_mu[i];
        const std::array<std::pair<double, double>, 3> chemical_potentials{{
            {derived.mu_cubic_over_rt[i], analytic_cubic_mu},
            {derived.mu_association_over_rt[i], analytic_assoc_mu[i]},
            {derived.mu_total_over_rt[i], analytic_total_mu}}};
        for (std::size_t term = 0U; term < chemical_potentials.size(); ++term) {
            const double delta = std::abs(
                chemical_potentials[term].first - chemical_potentials[term].second);
            max_mu_delta = std::max(max_mu_delta, delta);
            require_close(
                chemical_potentials[term].first,
                chemical_potentials[term].second,
                roundoff_tolerance(chemical_potentials[term].second), 0.0,
                std::string(label) + ": AD residual chemical potential mismatch component=" +
                    std::to_string(i) + " term=" + std::to_string(term));
        }

        const double ln_phi_delta = std::abs(
            derived.ln_phi[i] - analytic_root.ln_phi[i]);
        max_ln_phi_delta = std::max(max_ln_phi_delta, ln_phi_delta);
        require_close(
            derived.ln_phi[i], analytic_root.ln_phi[i],
            roundoff_tolerance(analytic_root.ln_phi[i]), 0.0,
            std::string(label) + ": AD ln(phi) mismatch component=" +
                std::to_string(i));
    }
}

void run_internal_equivalence() {
    double max_q_value_delta = 0.0;
    double max_pressure_delta_pa = 0.0;
    double max_mu_delta = 0.0;
    double max_ln_phi_delta = 0.0;

    for (const bool swapped : {false, true}) {
        const auto literature = cpa_physical_test::parameters(swapped);
        const auto parity = cpa_thermopack_snapshot::parameters(swapped);
        for (std::size_t state_index = 0U;
             state_index < external::states.size(); ++state_index) {
            const auto& state = external::states[state_index];
            for (const auto* phase : {&state.liquid, &state.vapor}) {
                const auto composition = cpa_physical_test::composition(
                    phase->composition_methanol, swapped);
                const std::string suffix =
                    " state=" + std::to_string(state_index) +
                    (phase == &state.liquid ? " phase=liquid" : " phase=vapor") +
                    (swapped ? " order=swapped" : " order=normal");
                check_internal_state(
                    state.temperature_k, phase->molar_density_mol_per_m3,
                    composition, literature, "literature" + suffix,
                    max_q_value_delta, max_pressure_delta_pa,
                    max_mu_delta, max_ln_phi_delta);
                check_internal_state(
                    state.temperature_k, phase->molar_density_mol_per_m3,
                    composition, parity, "parity" + suffix,
                    max_q_value_delta, max_pressure_delta_pa,
                    max_mu_delta, max_ln_phi_delta);
            }
        }
    }

    std::cout << std::setprecision(17)
              << "CPA_HELMHOLTZ_INTERNAL_OK"
              << " max_abs_dQ_value=" << max_q_value_delta
              << " max_abs_dP_pa=" << max_pressure_delta_pa
              << " max_abs_dmu_over_rt=" << max_mu_delta
              << " max_abs_dlnphi=" << max_ln_phi_delta
              << " snapshots=literature,thermopack-parity"
              << " component_orders=normal,swapped"
              << '\n';
}

void check_external_phase(
    const external::ThermoPackPhaseKernelState& state,
    const external::PhaseReference& phase,
    const th::CpaParameterSet& parameters,
    double& max_pressure_physical,
    double& max_pressure_association,
    double& max_mu_cubic,
    double& max_mu_association,
    double& max_ln_phi) {
    const auto composition = cpa_physical_test::composition(
        phase.composition_methanol, false);
    const auto primal = th::evaluate_cpa_phase_at_density(
        state.temperature_k, phase.molar_density_mol_per_m3,
        composition, parameters);
    const auto derived = derive_from_helmholtz(
        state.pressure_pa, state.temperature_k,
        phase.molar_density_mol_per_m3,
        composition, parameters, primal.association);

    const double d_p_physical = std::abs(
        derived.pressure_physical_pa - phase.pressure_physical_pa);
    const double d_p_association = std::abs(
        derived.pressure_association_pa - phase.pressure_association_pa);
    max_pressure_physical = std::max(max_pressure_physical, d_p_physical);
    max_pressure_association = std::max(max_pressure_association, d_p_association);
    require(d_p_physical <= thresholds::phase_max_abs_pressure_physical_pa,
            "Helmholtz AD physical pressure exceeded frozen ThermoPack parity gate");
    require(d_p_association <= thresholds::phase_max_abs_pressure_association_pa,
            "Helmholtz AD association pressure exceeded frozen ThermoPack parity gate");

    const std::array<double, 2> reference_mu_cubic{{
        phase.mu_cubic_over_rt_methanol,
        phase.mu_cubic_over_rt_water}};
    const std::array<double, 2> reference_mu_assoc{{
        phase.mu_association_over_rt_methanol,
        phase.mu_association_over_rt_water}};
    const std::array<double, 2> reference_ln_phi{{
        phase.ln_phi_methanol,
        phase.ln_phi_water}};
    for (std::size_t i = 0U; i < 2U; ++i) {
        const double d_mu_cubic = std::abs(
            derived.mu_cubic_over_rt[i] - reference_mu_cubic[i]);
        const double d_mu_assoc = std::abs(
            derived.mu_association_over_rt[i] - reference_mu_assoc[i]);
        const double d_ln_phi = std::abs(
            derived.ln_phi[i] - reference_ln_phi[i]);
        max_mu_cubic = std::max(max_mu_cubic, d_mu_cubic);
        max_mu_association = std::max(max_mu_association, d_mu_assoc);
        max_ln_phi = std::max(max_ln_phi, d_ln_phi);
        require(d_mu_cubic <= thresholds::phase_max_abs_mu_cubic_over_rt,
                "Helmholtz AD cubic mu exceeded frozen ThermoPack parity gate");
        require(d_mu_assoc <= thresholds::phase_max_abs_mu_association_over_rt,
                "Helmholtz AD association mu exceeded frozen ThermoPack parity gate");
        require(d_ln_phi <= thresholds::phase_max_abs_ln_phi,
                "Helmholtz AD ln(phi) exceeded frozen ThermoPack parity gate");
    }
}

void run_external_parity() {
    const auto parameters = cpa_thermopack_snapshot::parameters(false);
    double max_pressure_physical = 0.0;
    double max_pressure_association = 0.0;
    double max_mu_cubic = 0.0;
    double max_mu_association = 0.0;
    double max_ln_phi = 0.0;

    for (const auto& state : external::states) {
        check_external_phase(
            state, state.liquid, parameters,
            max_pressure_physical, max_pressure_association,
            max_mu_cubic, max_mu_association, max_ln_phi);
        check_external_phase(
            state, state.vapor, parameters,
            max_pressure_physical, max_pressure_association,
            max_mu_cubic, max_mu_association, max_ln_phi);
    }

    std::cout << std::setprecision(17)
              << "CPA_HELMHOLTZ_THERMOPACK_OK"
              << " threshold_contract=" << thresholds::contract
              << " max_abs_dPphysical_pa=" << max_pressure_physical
              << " max_abs_dPassociation_pa=" << max_pressure_association
              << " max_abs_dmu_cubic_over_rt=" << max_mu_cubic
              << " max_abs_dmu_association_over_rt=" << max_mu_association
              << " max_abs_dlnphi=" << max_ln_phi
              << " states=10"
              << '\n';
}

} // namespace

int main() {
    try {
        run_internal_equivalence();
        run_external_parity();
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "[FAIL] " << error.what() << '\n';
        return 1;
    }
}
