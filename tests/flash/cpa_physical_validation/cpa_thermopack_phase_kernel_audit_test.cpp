#include <mpmc/flash/cpa_split.hpp>

#include "cpa_thermopack_parameter_snapshot.hpp"
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
namespace parity = cpa_thermopack_snapshot;

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

struct ParitySummary {
    double max_common_tv_total_pressure_absolute{};
    double max_common_tv_physical_pressure_absolute{};
    double max_common_tv_association_pressure_absolute{};
    double max_common_tv_mu_cubic_absolute{};
    double max_common_tv_mu_association_absolute{};
    double max_root_rho_relative{};
    double max_root_z_absolute{};
    double max_root_ln_phi_absolute{};
};

void audit_phase(
    const th::CpaPtPhase& literature_model,
    const th::CpaPtPhase& parity_model,
    const th::CpaPtOptions& pt_options,
    double pressure_pa,
    double temperature_k,
    const ref::PhaseReference& external,
    bool liquid_side,
    std::string_view phase_name,
    Summary& literature_summary,
    ParitySummary& parity_summary) {
    const std::vector<double> composition{
        external.composition_methanol, 1.0 - external.composition_methanol};

    // Literature snapshot: retain the diagnostic showing why a scientifically
    // distinct parameter dataset must not be mislabeled as external parity.
    const auto literature_state = th::evaluate_cpa_phase_at_density(
        temperature_k, external.molar_density_mol_per_m3,
        composition, literature_model.parameters(), pt_options.phase);
    require(literature_state.association.converged(),
            "MPMC literature common-TV association state did not converge");
    const auto literature_mu = decompose_mpmc_mu(
        temperature_k, composition, literature_model.parameters(), literature_state);

    const double literature_d_p_physical =
        literature_state.pressure_physical_pa - external.pressure_physical_pa;
    const double literature_d_p_association =
        literature_state.pressure_association_pa - external.pressure_association_pa;
    const double literature_d_p_total = literature_state.pressure_pa - pressure_pa;
    const double literature_d_mu_cubic_meoh =
        literature_mu.cubic_over_rt[0] - external.mu_cubic_over_rt_methanol;
    const double literature_d_mu_cubic_h2o =
        literature_mu.cubic_over_rt[1] - external.mu_cubic_over_rt_water;
    const double literature_d_mu_assoc_meoh =
        literature_mu.association_over_rt[0] - external.mu_association_over_rt_methanol;
    const double literature_d_mu_assoc_h2o =
        literature_mu.association_over_rt[1] - external.mu_association_over_rt_water;

    literature_summary.max_common_tv_physical_pressure_relative = std::max(
        literature_summary.max_common_tv_physical_pressure_relative,
        scaled_difference(
            literature_state.pressure_physical_pa, external.pressure_physical_pa));
    literature_summary.max_common_tv_association_pressure_relative = std::max(
        literature_summary.max_common_tv_association_pressure_relative,
        scaled_difference(
            literature_state.pressure_association_pa, external.pressure_association_pa));
    literature_summary.max_common_tv_total_pressure_absolute = std::max(
        literature_summary.max_common_tv_total_pressure_absolute,
        std::abs(literature_d_p_total));
    literature_summary.max_common_tv_mu_cubic_absolute = std::max(
        literature_summary.max_common_tv_mu_cubic_absolute,
        std::max(std::abs(literature_d_mu_cubic_meoh),
                 std::abs(literature_d_mu_cubic_h2o)));
    literature_summary.max_common_tv_mu_association_absolute = std::max(
        literature_summary.max_common_tv_mu_association_absolute,
        std::max(std::abs(literature_d_mu_assoc_meoh),
                 std::abs(literature_d_mu_assoc_h2o)));

    const auto literature_roots = literature_model.roots(
        pressure_pa, temperature_k, composition, pt_options);
    const auto& literature_root = select_root(literature_roots, liquid_side);
    require(literature_root.ln_phi.size() == 2U,
            "MPMC literature phase-kernel root lost component fugacity coefficients");
    const double literature_rho_relative =
        std::abs(literature_root.molar_density_mol_per_m3 -
                 external.molar_density_mol_per_m3) /
        external.molar_density_mol_per_m3;
    const double literature_d_z =
        literature_root.compressibility_factor - external.compressibility_factor;
    const double literature_d_ln_phi_meoh =
        literature_root.ln_phi[0] - external.ln_phi_methanol;
    const double literature_d_ln_phi_h2o =
        literature_root.ln_phi[1] - external.ln_phi_water;

    literature_summary.max_root_rho_relative = std::max(
        literature_summary.max_root_rho_relative, literature_rho_relative);
    literature_summary.max_root_z_absolute = std::max(
        literature_summary.max_root_z_absolute, std::abs(literature_d_z));
    literature_summary.max_root_ln_phi_absolute = std::max(
        literature_summary.max_root_ln_phi_absolute,
        std::max(std::abs(literature_d_ln_phi_meoh),
                 std::abs(literature_d_ln_phi_h2o)));

    // Parity snapshot: this is the only parameter source used for the actual
    // MPMC_HNU <-> ThermoPack phase-kernel parity path.
    const auto parity_state = th::evaluate_cpa_phase_at_density(
        temperature_k, external.molar_density_mol_per_m3,
        composition, parity_model.parameters(), pt_options.phase);
    require(parity_state.association.converged(),
            "named parity common-TV association state did not converge");
    const auto parity_mu = decompose_mpmc_mu(
        temperature_k, composition, parity_model.parameters(), parity_state);

    const double parity_d_p_physical =
        parity_state.pressure_physical_pa - external.pressure_physical_pa;
    const double parity_d_p_association =
        parity_state.pressure_association_pa - external.pressure_association_pa;
    const double parity_d_p_total = parity_state.pressure_pa - pressure_pa;
    const double parity_d_mu_cubic_meoh =
        parity_mu.cubic_over_rt[0] - external.mu_cubic_over_rt_methanol;
    const double parity_d_mu_cubic_h2o =
        parity_mu.cubic_over_rt[1] - external.mu_cubic_over_rt_water;
    const double parity_d_mu_assoc_meoh =
        parity_mu.association_over_rt[0] - external.mu_association_over_rt_methanol;
    const double parity_d_mu_assoc_h2o =
        parity_mu.association_over_rt[1] - external.mu_association_over_rt_water;

    parity_summary.max_common_tv_total_pressure_absolute = std::max(
        parity_summary.max_common_tv_total_pressure_absolute,
        std::abs(parity_d_p_total));
    parity_summary.max_common_tv_physical_pressure_absolute = std::max(
        parity_summary.max_common_tv_physical_pressure_absolute,
        std::abs(parity_d_p_physical));
    parity_summary.max_common_tv_association_pressure_absolute = std::max(
        parity_summary.max_common_tv_association_pressure_absolute,
        std::abs(parity_d_p_association));
    parity_summary.max_common_tv_mu_cubic_absolute = std::max(
        parity_summary.max_common_tv_mu_cubic_absolute,
        std::max(std::abs(parity_d_mu_cubic_meoh),
                 std::abs(parity_d_mu_cubic_h2o)));
    parity_summary.max_common_tv_mu_association_absolute = std::max(
        parity_summary.max_common_tv_mu_association_absolute,
        std::max(std::abs(parity_d_mu_assoc_meoh),
                 std::abs(parity_d_mu_assoc_h2o)));

    const auto parity_roots = parity_model.roots(
        pressure_pa, temperature_k, composition, pt_options);
    const auto& parity_root = select_root(parity_roots, liquid_side);
    require(parity_root.ln_phi.size() == 2U,
            "named parity root lost component fugacity coefficients");
    const double parity_rho_relative =
        std::abs(parity_root.molar_density_mol_per_m3 -
                 external.molar_density_mol_per_m3) /
        external.molar_density_mol_per_m3;
    const double parity_d_z =
        parity_root.compressibility_factor - external.compressibility_factor;
    const double parity_d_ln_phi_meoh =
        parity_root.ln_phi[0] - external.ln_phi_methanol;
    const double parity_d_ln_phi_h2o =
        parity_root.ln_phi[1] - external.ln_phi_water;

    parity_summary.max_root_rho_relative = std::max(
        parity_summary.max_root_rho_relative, parity_rho_relative);
    parity_summary.max_root_z_absolute = std::max(
        parity_summary.max_root_z_absolute, std::abs(parity_d_z));
    parity_summary.max_root_ln_phi_absolute = std::max(
        parity_summary.max_root_ln_phi_absolute,
        std::max(std::abs(parity_d_ln_phi_meoh),
                 std::abs(parity_d_ln_phi_h2o)));

    std::cout << "CPA_THERMOPACK_LITERATURE_PHASE_KERNEL"
              << " PPa=" << pressure_pa
              << " phase=" << phase_name
              << " xMeOH=" << external.composition_methanol
              << " d_rho_rel=" << literature_rho_relative
              << " d_Z=" << literature_d_z
              << " d_lnphi_MeOH=" << literature_d_ln_phi_meoh
              << " d_lnphi_H2O=" << literature_d_ln_phi_h2o
              << " d_Pphysical=" << literature_d_p_physical
              << " d_Passoc=" << literature_d_p_association
              << " d_Ptotal_at_tp_rho=" << literature_d_p_total
              << " d_mu_cubic_MeOH=" << literature_d_mu_cubic_meoh
              << " d_mu_cubic_H2O=" << literature_d_mu_cubic_h2o
              << " d_mu_assoc_MeOH=" << literature_d_mu_assoc_meoh
              << " d_mu_assoc_H2O=" << literature_d_mu_assoc_h2o
              << '\n';

    std::cout << "CPA_THERMOPACK_PARITY_PHASE_KERNEL"
              << " PPa=" << pressure_pa
              << " phase=" << phase_name
              << " dataset=" << parity_model.parameters().dataset_id()
              << " d_Pphysical=" << parity_d_p_physical
              << " d_Passoc=" << parity_d_p_association
              << " d_Ptotal_at_tp_rho=" << parity_d_p_total
              << " d_mu_cubic_MeOH=" << parity_d_mu_cubic_meoh
              << " d_mu_cubic_H2O=" << parity_d_mu_cubic_h2o
              << " d_mu_assoc_MeOH=" << parity_d_mu_assoc_meoh
              << " d_mu_assoc_H2O=" << parity_d_mu_assoc_h2o
              << " d_rho_rel=" << parity_rho_relative
              << " d_Z=" << parity_d_z
              << " d_lnphi_MeOH=" << parity_d_ln_phi_meoh
              << " d_lnphi_H2O=" << parity_d_ln_phi_h2o
              << '\n';
}

} // namespace

int main() {
    try {
        const auto literature_parameters = cpa_physical_test::parameters(false);
        const auto parity_parameters = parity::parameters(false);
        require(parity_parameters.dataset_id() == parity::parity_dataset_id &&
                    parity_parameters.revision() == parity::parity_revision,
                "phase-kernel parity audit lost named parameter-snapshot identity");

        const auto literature_model =
            th::CpaPtPhase::from_parameters(literature_parameters);
        const auto parity_model = th::CpaPtPhase::from_parameters(parity_parameters);
        const auto pt_options = fl::cpa_pt_vle_default_phase_options();
        require(std::abs(ref::gas_constant_j_per_mol_k -
                         th::cpa_gas_constant_j_per_mol_k) <=
                    8.0 * std::numeric_limits<double>::epsilon() *
                        th::cpa_gas_constant_j_per_mol_k,
                "MPMC/ThermoPack gas constants differ in phase-kernel audit");

        require(std::abs(literature_parameters.pure(0).critical_temperature_k -
                         parity_parameters.pure(0).critical_temperature_k) > 0.0 &&
                    std::abs(literature_parameters.pure(1).critical_temperature_k -
                             parity_parameters.pure(1).critical_temperature_k) > 0.0,
                "literature and parity parameter snapshots unexpectedly collapsed");
        require(parity_parameters.pure(0).critical_temperature_k ==
                    ref::methanol_alpha_critical_temperature_k &&
                    parity_parameters.pure(1).critical_temperature_k ==
                    ref::water_alpha_critical_temperature_k,
                "named parity snapshot no longer matches frozen ThermoPack alpha Tc");

        Summary literature_summary;
        ParitySummary parity_summary;
        std::cout << std::setprecision(17);
        for (const auto& external : ref::states) {
            audit_phase(literature_model, parity_model, pt_options,
                        external.pressure_pa, external.temperature_k,
                        external.liquid, true, "liquid",
                        literature_summary, parity_summary);
            audit_phase(literature_model, parity_model, pt_options,
                        external.pressure_pa, external.temperature_k,
                        external.vapor, false, "vapor",
                        literature_summary, parity_summary);
        }

        std::cout << "CPA_THERMOPACK_LITERATURE_PHASE_KERNEL_SUMMARY"
                  << " states=" << ref::states.size()
                  << " phase_states=" << 2U * ref::states.size()
                  << " max_root_d_rho_rel=" << literature_summary.max_root_rho_relative
                  << " max_root_abs_d_Z=" << literature_summary.max_root_z_absolute
                  << " max_root_abs_d_lnphi=" << literature_summary.max_root_ln_phi_absolute
                  << " max_common_tv_rel_d_Pphysical="
                  << literature_summary.max_common_tv_physical_pressure_relative
                  << " max_common_tv_rel_d_Passoc="
                  << literature_summary.max_common_tv_association_pressure_relative
                  << " max_common_tv_abs_d_Ptotal="
                  << literature_summary.max_common_tv_total_pressure_absolute
                  << " max_common_tv_abs_d_mu_cubic="
                  << literature_summary.max_common_tv_mu_cubic_absolute
                  << " max_common_tv_abs_d_mu_assoc="
                  << literature_summary.max_common_tv_mu_association_absolute
                  << " magnitude_gate=none_diagnostic"
                  << '\n';

        std::cout << "CPA_THERMOPACK_PARITY_PHASE_KERNEL_SUMMARY"
                  << " dataset=" << parity_parameters.dataset_id()
                  << " revision=" << parity_parameters.revision()
                  << " states=" << ref::states.size()
                  << " phase_states=" << 2U * ref::states.size()
                  << " max_common_tv_abs_d_Pphysical="
                  << parity_summary.max_common_tv_physical_pressure_absolute
                  << " max_common_tv_abs_d_Passoc="
                  << parity_summary.max_common_tv_association_pressure_absolute
                  << " max_common_tv_abs_d_Ptotal="
                  << parity_summary.max_common_tv_total_pressure_absolute
                  << " max_common_tv_abs_d_mu_cubic="
                  << parity_summary.max_common_tv_mu_cubic_absolute
                  << " max_common_tv_abs_d_mu_assoc="
                  << parity_summary.max_common_tv_mu_association_absolute
                  << " max_root_d_rho_rel=" << parity_summary.max_root_rho_relative
                  << " max_root_abs_d_Z=" << parity_summary.max_root_z_absolute
                  << " max_root_abs_d_lnphi=" << parity_summary.max_root_ln_phi_absolute
                  << " magnitude_gate=none_parity_threshold_not_frozen"
                  << '\n';
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "[FAIL] " << error.what() << '\n';
        return 1;
    }
}
