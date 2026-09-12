#include <mpmc/flash/cpa_split.hpp>

#include "test_support.hpp"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <stdexcept>
#include <string_view>
#include <vector>

namespace {
namespace fl = mpmc::flash;
namespace th = mpmc::thermodynamics;

void require(bool condition, const char* message) {
    if (!condition) { throw std::runtime_error(message); }
}

th::CpaPtOptions physical_pt_options() {
    th::CpaPtOptions options;
    options.scan_intervals = 256U;
    options.max_evaluations = 4096U;
    // The outer flash gate is 1e-11 in log fugacity. Keep the density-root
    // residual comfortably below that scale without changing any flash/TPD tolerance.
    options.pressure_absolute_tolerance_pa = 1.0e-7;
    options.pressure_relative_tolerance = 1.0e-12;
    return options;
}

fl::PtSplitOptions physical_split_options() {
    fl::PtSplitOptions options;
    options.initial_stability.automatic_starts = false;
    options.final_stability.automatic_starts = false;
    return options;
}

void parameter_provenance() {
    const auto parameters = cpa_physical_test::parameters(false);
    require(parameters.size() == 2U, "physical fixture lost two-component snapshot");
    require(parameters.components().at(0).id == "METHANOL" &&
                parameters.components().at(1).id == "WATER",
            "physical fixture component order changed");
    require(std::abs(parameters.kij(0U, 1U) -
                     cpa_physical_test::kij_methanol_water) < 1.0e-15,
            "literature CR-1 kij changed");
    require(parameters.pure(0U).sites.size() == 2U &&
                parameters.pure(0U).sites[0].multiplicity == 1U &&
                parameters.pure(0U).sites[1].multiplicity == 1U,
            "methanol 2B scheme changed");
    require(parameters.pure(1U).sites.size() == 2U &&
                parameters.pure(1U).sites[0].multiplicity == 2U &&
                parameters.pure(1U).sites[1].multiplicity == 2U,
            "water 4C scheme changed");
    require(parameters.association_records().size() == 4U,
            "explicit self/cross association records changed");
    for (const auto& record : parameters.pure_records()) {
        require(record.critical_temperature_k.source.kind == th::SourceKind::literature &&
                    record.a0_pa_m6_per_mol2.source.kind == th::SourceKind::literature &&
                    record.b_m3_per_mol.source.kind == th::SourceKind::literature &&
                    record.c1_dimensionless.source.kind == th::SourceKind::literature,
                "pure parameter lost literature provenance");
    }
    for (const auto& record : parameters.association_records()) {
        require(record.epsilon_j_per_mol.source.kind == th::SourceKind::literature &&
                    record.beta_dimensionless.source.kind == th::SourceKind::literature,
                "association pair lost literature provenance");
    }
}

void require_active_association(
    const th::CpaPtPhase& model, const th::CpaPtOptions& options,
    double pressure_pa, const fl::PtSplitState& point) {
    const auto roots = model.roots(
        pressure_pa, cpa_physical_test::temperature_k,
        point.fractions.liquid, options);
    require(roots.status == th::CpaPtRootStatus::success &&
                point.liquid.activity.branch < roots.roots.size(),
            "physical liquid root became unavailable");
    const auto& root = roots.roots[point.liquid.activity.branch];
    const auto state = th::evaluate_cpa_phase_at_density(
        cpa_physical_test::temperature_k, root.molar_density_mol_per_m3,
        point.fractions.liquid, model.parameters(), options.phase);
    require(state.association.status == th::CpaAssociationStatus::success &&
                !state.association.sites.empty(),
            "physical VLE did not traverse converged association state");
    require(std::any_of(
                state.association.sites.begin(), state.association.sites.end(),
                [](const th::CpaAssociationSiteState& site) {
                    return site.unbonded_fraction < 0.999;
                }),
            "physical VLE association sites are effectively inactive");
}

fl::CpaPtSplitResult solve_point(
    const cpa_physical_test::ExperimentalVlePoint& experimental,
    bool swapped = false) {
    const auto parameters = cpa_physical_test::parameters(swapped);
    const auto model = th::CpaPtPhase::from_parameters(parameters);
    fl::CpaVleEvaluator evaluator(model, physical_pt_options());
    const auto feed = cpa_physical_test::feed(experimental, swapped);
    const auto starts = cpa_physical_test::starts(experimental, swapped);
    return fl::solve_cpa_pt_vle(
        experimental.pressure_pa, cpa_physical_test::temperature_k,
        feed, evaluator, physical_split_options(), starts, starts);
}

double experimental_pair_residual(
    fl::CpaVleEvaluator& evaluator,
    const cpa_physical_test::ExperimentalVlePoint& experimental) {
    const auto x = cpa_physical_test::composition(experimental.liquid_methanol);
    const auto y = cpa_physical_test::composition(experimental.vapor_methanol);
    const auto liquid = evaluator(
        experimental.pressure_pa, cpa_physical_test::temperature_k,
        x, fl::PtPhaseRole::liquid_candidate);
    const auto vapor = evaluator(
        experimental.pressure_pa, cpa_physical_test::temperature_k,
        y, fl::PtPhaseRole::vapor_candidate);
    double norm = 0.0;
    for (std::size_t i = 0; i < x.size(); ++i) {
        norm = std::max(norm, std::abs(
            std::log(x[i]) + liquid.activity.ln_phi[i] -
            std::log(y[i]) - vapor.activity.ln_phi[i]));
    }
    return norm;
}

void print_solve_diagnostic(
    const cpa_physical_test::ExperimentalVlePoint& experimental,
    double direct_residual, const fl::CpaPtSplitResult& result) {
    std::cout << "CPA_PHYSICAL_AUDIT PkPa=" << experimental.pressure_pa / 1000.0
              << " direct_mu=" << direct_residual
              << " initial=" << static_cast<int>(result.solution.initial_stability.status)
              << " status=" << static_cast<int>(result.solution.status)
              << " attempts=" << result.solution.attempts.size();
    if (result.solution.candidate()) {
        const auto& point = *result.solution.candidate();
        std::cout << " x_calc=" << point.fractions.liquid[0]
                  << " y_calc=" << point.fractions.vapor[0]
                  << " beta=" << point.fractions.vapor_fraction
                  << " fug=" << point.fugacity_norm;
    }
    if (result.solution.final_stability) {
        const auto& final = *result.solution.final_stability;
        std::cout << " final=" << static_cast<int>(final.status);
        if (final.lowest_sampled) {
            std::cout << " tpd_min=" << final.lowest_sampled->value
                      << " tpd_stat=" << final.lowest_sampled->stationarity
                      << " tpd_branch=" << final.lowest_sampled->branch;
        }
        for (const auto& trial : final.trials) {
            if (trial.status == fl::StabilityTrialStatus::negative_tpd && trial.point) {
                std::cout << " neg_tpd=" << trial.point->value
                          << ",neg_stat=" << trial.point->stationarity;
            }
        }
    } else {
        std::cout << " final=none";
    }
    std::cout << " diag=" << result.solution.diagnostic << '\n';
}

void literature_vle_points() {
    const auto parameters = cpa_physical_test::parameters(false);
    const auto model = th::CpaPtPhase::from_parameters(parameters);
    const auto pt_options = physical_pt_options();
    fl::CpaVleEvaluator evaluator(model, pt_options);
    const auto split_options = physical_split_options();

    double sum_liquid_error = 0.0;
    double sum_vapor_error = 0.0;
    double max_liquid_error = 0.0;
    double max_vapor_error = 0.0;
    std::size_t accepted = 0U;

    for (const auto& experimental : cpa_physical_test::points()) {
        const double direct_residual = experimental_pair_residual(evaluator, experimental);
        const auto feed = cpa_physical_test::feed(experimental, false);
        const auto starts = cpa_physical_test::starts(experimental, false);
        const auto result = fl::solve_cpa_pt_vle(
            experimental.pressure_pa, cpa_physical_test::temperature_k,
            feed, evaluator, split_options, starts, starts);
        print_solve_diagnostic(experimental, direct_residual, result);
        if (result.solution.status != fl::PtSplitStatus::two_phase_no_instability_found ||
            result.solution.candidate() == nullptr || !result.solution.final_stability ||
            result.solution.final_stability->status != fl::StabilityStatus::no_instability_found) {
            continue;
        }
        const auto& point = *result.solution.candidate();
        const double dx = std::abs(
            point.fractions.liquid[0] - experimental.liquid_methanol);
        const double dy = std::abs(
            point.fractions.vapor[0] - experimental.vapor_methanol);
        sum_liquid_error += dx;
        sum_vapor_error += dy;
        max_liquid_error = std::max(max_liquid_error, dx);
        max_vapor_error = std::max(max_vapor_error, dy);
        ++accepted;
        require(point.fugacity_norm <= result.solution.options.iteration.fugacity_tolerance,
                "accepted physical point lost fugacity equality");
        require(point.fractions.mass_absolute <=
                    result.solution.options.iteration.mass_absolute_tolerance &&
                    point.fractions.mass_relative <=
                    result.solution.options.iteration.mass_relative_tolerance,
                "accepted physical point lost material balance");
        require_active_association(model, pt_options, experimental.pressure_pa, point);
    }

    require(accepted == cpa_physical_test::points().size(),
            "CPA physical VLE did not accept every declared literature point");
    const double mean_dx = sum_liquid_error / static_cast<double>(accepted);
    const double mean_dy = sum_vapor_error / static_cast<double>(accepted);
    std::cout << "CPA_PHYSICAL_VLE_SUMMARY mean_dx=" << mean_dx
              << " mean_dy=" << mean_dy
              << " max_dx=" << max_liquid_error
              << " max_dy=" << max_vapor_error << '\n';
    require(mean_dx <= 0.025 && mean_dy <= 0.025 &&
                max_liquid_error <= 0.05 && max_vapor_error <= 0.05,
            "CPA literature VLE deviations exceed physical-validation envelope");
}

void component_permutation() {
    const auto& experimental = cpa_physical_test::points()[1];
    const auto first = solve_point(experimental, false);
    const auto second = solve_point(experimental, true);
    require(first.solution.status == fl::PtSplitStatus::two_phase_no_instability_found &&
                second.solution.status == fl::PtSplitStatus::two_phase_no_instability_found &&
                first.solution.candidate() != nullptr && second.solution.candidate() != nullptr,
            "CPA physical component permutation changed accepted topology");
    const auto& a = *first.solution.candidate();
    const auto& b = *second.solution.candidate();
    require(std::abs(a.fractions.vapor_fraction - b.fractions.vapor_fraction) < 1.0e-9 &&
                std::abs(a.fractions.liquid[0] - b.fractions.liquid[1]) < 1.0e-8 &&
                std::abs(a.fractions.vapor[0] - b.fractions.vapor[1]) < 1.0e-8,
            "CPA physical VLE changed under ordered-component permutation");
}

using Test = std::pair<std::string_view, void (*)()>;
constexpr Test tests[]{
    {"parameter_provenance", parameter_provenance},
    {"literature_vle_points", literature_vle_points},
    {"component_permutation", component_permutation}};

} // namespace

int main(int argc, char** argv) {
    try {
        if (argc != 2) { throw std::invalid_argument("one test name required"); }
        for (const auto& [name, run] : tests) {
            if (name == argv[1]) {
                run();
                std::cout << "[PASS] " << name << '\n';
                return 0;
            }
        }
        throw std::invalid_argument("unknown test name");
    } catch (const std::exception& error) {
        std::cerr << "[FAIL] " << error.what() << '\n';
        return 1;
    }
}
