#include <mpmc/flash/cpa_split.hpp>

#include "cpa_thermopack_phase_kernel_generated.hpp"
#include "test_support.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <iomanip>
#include <iostream>
#include <limits>
#include <span>
#include <stdexcept>
#include <string_view>
#include <vector>

namespace {
namespace fl = mpmc::flash;
namespace th = mpmc::thermodynamics;
namespace ref = cpa_thermopack_phase_kernel;

void require(bool value, const char* message) {
    if (!value) { throw std::runtime_error(message); }
}

struct MuDecomposition {
    std::vector<double> cubic_over_rt;
    std::vector<double> association_over_rt;
};

MuDecomposition decompose_mpmc_mu(
    double temperature_k,
    std::span<const double> composition,
    const th::CpaParameterSet& parameters,
    const th::CpaPhaseState& state) {
    const double rho = state.molar_density_mol_per_m3;
    const double rt = th::cpa_gas_constant_j_per_mol_k * temperature_k;
    const double b = state.b_mix_m3_per_mol;
    const double a = state.a_mix_pa_m6_per_mol2;
    const double b_rho = b * rho;
    require(std::isfinite(b_rho) && b_rho > 0.0 && b_rho < 1.0,
            "MPMC phase-kernel audit reached invalid reduced density");

    const auto pure_a = th::cpa_detail::cpa_pure_a(temperature_k, parameters);
    const auto sums = th::cpa_detail::cpa_a_sums(composition, pure_a, parameters);
    const double z_physical = state.pressure_physical_pa / (rho * rt);
    const double a_over_brt = a / (b * rt);
    const double log_free_volume = std::log1p(-b_rho);
    const double log_attraction_volume = std::log1p(b_rho);

    std::vector<double> association_log_x(parameters.size(), 0.0);
    double association_sum = 0.0;
    for (const auto& site : state.association.sites) {
        require(site.unbonded_fraction > 0.0 &&
                    std::isfinite(site.unbonded_fraction),
                "MPMC phase-kernel audit found invalid association site fraction");
        const double multiplicity = static_cast<double>(site.multiplicity);
        association_log_x[site.component_index] +=
            multiplicity * std::log(site.unbonded_fraction);
        association_sum += composition[site.component_index] * multiplicity *
            (1.0 - site.unbonded_fraction);
    }

    const double g = state.association.radial_distribution;
    require(std::isfinite(g) && g > 0.0 &&
                std::isfinite(association_sum) && association_sum >= 0.0,
            "MPMC phase-kernel audit found invalid association decomposition state");

    MuDecomposition result;
    result.cubic_over_rt.resize(parameters.size(), 0.0);
    result.association_over_rt.resize(parameters.size(), 0.0);
    for (std::size_t i = 0; i < parameters.size(); ++i) {
        const double b_i = parameters.pure(i).b_m3_per_mol;
        const double b_ratio = b_i / b;
        const double attraction_ratio = 2.0 * sums[i] / a - b_ratio;
        result.cubic_over_rt[i] =
            b_ratio * (z_physical - 1.0) - log_free_volume -
            a_over_brt * attraction_ratio * log_attraction_volume;
        result.association_over_rt[i] = association_log_x[i] -
            (1.9 / 8.0) * rho * b_i * g * association_sum;
        require(std::isfinite(result.cubic_over_rt[i]) &&
                    std::isfinite(result.association_over_rt[i]),
                "MPMC phase-kernel audit produced nonfinite chemical-potential contribution");
    }
    return result;
}

const th::CpaPtRoot& select_root(
    const th::CpaPtRootSet& roots, bool liquid_side) {
    require(roots.status == th::CpaPtRootStatus::success,
            "MPMC phase-kernel audit root search did not succeed");
    const th::CpaPtRoot* selected = nullptr;
    for (const auto& root : roots.roots) {
        if (root.pressure_slope_sign <= 0) { continue; }
        if (selected == nullptr ||
            (liquid_side &&
             root.molar_density_mol_per_m3 > selected->molar_density_mol_per_m3) ||
            (!liquid_side &&
             root.molar_density_mol_per_m3 < selected->molar_density_mol_per_m3)) {
            selected = &root;
        }
    }
    require(selected != nullptr,
            "MPMC phase-kernel audit found no mechanically stable density root");
    return *selected;
}

double scaled_difference(double actual, double expected) {
    return std::abs(actual - expected) /
           std::max({1.0, std::abs(actual), std::abs(expected)});
}

struct Summary {
    double max_root_rho_relative{};
    double max_root_z_absolute{};
    double max_root_ln_phi_absolute{};
    double max_common_tv_physical_pressure_relative{};
    double max_common_tv_association_pressure_relative{};
    double max_common_tv_total_pressure_absolute{};
    double max_common_tv_mu_cubic_absolute{};
    double max_common_tv_mu_association_absolute{};
};

void audit_phase(
    const th::CpaPtPhase& model,
    const th::CpaPtOptions& pt_options,
    double pressure_pa,
    double temperature_k,
    const ref::PhaseReference& external,
    bool liquid_side,
    std::string_view phase_name,
    Summary& summary) {
    const std::vector<double> composition{
        external.composition_methanol, 1.0 - external.composition_methanol};

    // First isolate the EOS kernel at the exact ThermoPack density.  This removes
    // density-root search and flash orchestration from the comparison.
    const auto state = th::evaluate_cpa_phase_at_density(
        temperature_k, external.molar_density_mol_per_m3,
        composition, model.parameters(), pt_options.phase);
    require(state.association.converged(),
            "MPMC common-TV association state did not converge");
    const auto mu = decompose_mpmc_mu(
        temperature_k, composition, model.parameters(), state);

    const double d_p_physical =
        state.pressure_physical_pa - external.pressure_physical_pa;
    const double d_p_association =
        state.pressure_association_pa - external.pressure_association_pa;
    const double d_p_total = state.pressure_pa - pressure_pa;
    const double d_mu_cubic_meoh =
        mu.cubic_over_rt[0] - external.mu_cubic_over_rt_methanol;
    const double d_mu_cubic_h2o =
        mu.cubic_over_rt[1] - external.mu_cubic_over_rt_water;
    const double d_mu_assoc_meoh =
        mu.association_over_rt[0] - external.mu_association_over_rt_methanol;
    const double d_mu_assoc_h2o =
        mu.association_over_rt[1] - external.mu_association_over_rt_water;

    summary.max_common_tv_physical_pressure_relative = std::max(
        summary.max_common_tv_physical_pressure_relative,
        scaled_difference(state.pressure_physical_pa, external.pressure_physical_pa));
    summary.max_common_tv_association_pressure_relative = std::max(
        summary.max_common_tv_association_pressure_relative,
        scaled_difference(state.pressure_association_pa, external.pressure_association_pa));
    summary.max_common_tv_total_pressure_absolute = std::max(
        summary.max_common_tv_total_pressure_absolute, std::abs(d_p_total));
    summary.max_common_tv_mu_cubic_absolute = std::max(
        summary.max_common_tv_mu_cubic_absolute,
        std::max(std::abs(d_mu_cubic_meoh), std::abs(d_mu_cubic_h2o)));
    summary.max_common_tv_mu_association_absolute = std::max(
        summary.max_common_tv_mu_association_absolute,
        std::max(std::abs(d_mu_assoc_meoh), std::abs(d_mu_assoc_h2o)));

    // Then include the independent p,T,x density-root search.  If common-TV
    // properties agree but these values do not, the discrepancy is root solving
    // rather than the EOS pressure/chemical-potential formulation itself.
    const auto roots = model.roots(
        pressure_pa, temperature_k, composition, pt_options);
    const auto& root = select_root(roots, liquid_side);
    require(root.ln_phi.size() == 2U,
            "MPMC phase-kernel audit root lost component fugacity coefficients");
    const double rho_relative =
        std::abs(root.molar_density_mol_per_m3 - external.molar_density_mol_per_m3) /
        external.molar_density_mol_per_m3;
    const double d_z = root.compressibility_factor - external.compressibility_factor;
    const double d_ln_phi_meoh = root.ln_phi[0] - external.ln_phi_methanol;
    const double d_ln_phi_h2o = root.ln_phi[1] - external.ln_phi_water;

    summary.max_root_rho_relative = std::max(
        summary.max_root_rho_relative, rho_relative);
    summary.max_root_z_absolute = std::max(
        summary.max_root_z_absolute, std::abs(d_z));
    summary.max_root_ln_phi_absolute = std::max(
        summary.max_root_ln_phi_absolute,
        std::max(std::abs(d_ln_phi_meoh), std::abs(d_ln_phi_h2o)));

    std::cout << "CPA_THERMOPACK_PHASE_KERNEL"
              << " PPa=" << pressure_pa
              << " phase=" << phase_name
              << " xMeOH=" << external.composition_methanol
              << " rho_tp=" << external.molar_density_mol_per_m3
              << " rho_mpmc=" << root.molar_density_mol_per_m3
              << " d_rho_rel=" << rho_relative
              << " Z_tp=" << external.compressibility_factor
              << " Z_mpmc=" << root.compressibility_factor
              << " d_Z=" << d_z
              << " d_lnphi_MeOH=" << d_ln_phi_meoh
              << " d_lnphi_H2O=" << d_ln_phi_h2o
              << " d_Pphysical=" << d_p_physical
              << " d_Passoc=" << d_p_association
              << " d_Ptotal_at_tp_rho=" << d_p_total
              << " d_mu_cubic_MeOH=" << d_mu_cubic_meoh
              << " d_mu_cubic_H2O=" << d_mu_cubic_h2o
              << " d_mu_assoc_MeOH=" << d_mu_assoc_meoh
              << " d_mu_assoc_H2O=" << d_mu_assoc_h2o
              << '\n';
}

} // namespace

int main() {
    try {
        const auto parameters = cpa_physical_test::parameters(false);
        const auto model = th::CpaPtPhase::from_parameters(parameters);
        const auto pt_options = fl::cpa_pt_vle_default_phase_options();
        require(std::abs(ref::gas_constant_j_per_mol_k -
                         th::cpa_gas_constant_j_per_mol_k) <=
                    8.0 * std::numeric_limits<double>::epsilon() *
                        th::cpa_gas_constant_j_per_mol_k,
                "MPMC/ThermoPack gas constants differ in phase-kernel audit");

        Summary summary;
        std::cout << std::setprecision(17);
        for (const auto& external : ref::states) {
            audit_phase(model, pt_options, external.pressure_pa,
                        external.temperature_k, external.liquid,
                        true, "liquid", summary);
            audit_phase(model, pt_options, external.pressure_pa,
                        external.temperature_k, external.vapor,
                        false, "vapor", summary);
        }

        std::cout << "CPA_THERMOPACK_PHASE_KERNEL_SUMMARY"
                  << " states=" << ref::states.size()
                  << " phase_states=" << 2U * ref::states.size()
                  << " max_root_d_rho_rel=" << summary.max_root_rho_relative
                  << " max_root_abs_d_Z=" << summary.max_root_z_absolute
                  << " max_root_abs_d_lnphi=" << summary.max_root_ln_phi_absolute
                  << " max_common_tv_rel_d_Pphysical="
                  << summary.max_common_tv_physical_pressure_relative
                  << " max_common_tv_rel_d_Passoc="
                  << summary.max_common_tv_association_pressure_relative
                  << " max_common_tv_abs_d_Ptotal="
                  << summary.max_common_tv_total_pressure_absolute
                  << " max_common_tv_abs_d_mu_cubic="
                  << summary.max_common_tv_mu_cubic_absolute
                  << " max_common_tv_abs_d_mu_assoc="
                  << summary.max_common_tv_mu_association_absolute
                  << " magnitude_gate=none_diagnostic"
                  << '\n';
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "[FAIL] " << error.what() << '\n';
        return 1;
    }
}
