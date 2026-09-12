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
using Vec = std::vector<double>;

void require(bool condition, const char* message) {
    if (!condition) { throw std::runtime_error(message); }
}

th::CpaPtOptions physical_pt_options() {
    th::CpaPtOptions options;
    options.scan_intervals = 256U;
    options.max_evaluations = 4096U;
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
    require(parameters.size() == 2U,
            "CPA physical fixture lost its two-component snapshot");
    require(parameters.components().at(0).id == "METHANOL" &&
                parameters.components().at(1).id == "WATER",
            "CPA physical fixture component order changed");
    require(std::abs(parameters.kij(0U, 1U) -
                     cpa_physical_test::kij_methanol_water) < 1.0e-15,
            "CPA physical fixture lost the literature CR-1 kij");
    require(parameters.pure(0U).sites.size() == 2U &&
                parameters.pure(0U).sites[0].multiplicity == 1U &&
                parameters.pure(0U).sites[1].multiplicity == 1U,
            "methanol is no longer represented as the literature 2B scheme");
    require(parameters.pure(1U).sites.size() == 2U &&
                parameters.pure(1U).sites[0].multiplicity == 2U &&
                parameters.pure(1U).sites[1].multiplicity == 2U,
            "water is no longer represented as the literature 4C scheme");
    require(parameters.association_records().size() == 4U,
            "CPA physical fixture lost explicit self/cross association pairs");
    for (const auto& record : parameters.pure_records()) {
        require(record.critical_temperature_k.source.kind == th::SourceKind::literature &&
                    record.a0_pa_m6_per_mol2.source.kind == th::SourceKind::literature &&
                    record.b_m3_per_mol.source.kind == th::SourceKind::literature &&
                    record.c1_dimensionless.source.kind == th::SourceKind::literature,
                "CPA physical pure parameter lost literature provenance");
    }
    for (const auto& record : parameters.association_records()) {
        require(record.epsilon_j_per_mol.source.kind == th::SourceKind::literature &&
                    record.beta_dimensionless.source.kind == th::SourceKind::literature,
                "CPA physical association pair lost literature provenance");
    }
}

void require_active_association(
    const th::CpaPtPhase& model,
    const th::CpaPtOptions& options,
    double pressure_pa,
    const fl::PtSplitState& point) {
    const auto roots = model.roots(
        pressure_pa, cpa_physical_test::temperature_k,
        point.fractions.liquid, options);
    require(roots.status == th::CpaPtRootStatus::success &&
                point.liquid.activity.branch < roots.roots.size(),
            "CPA physical liquid root became unavailable");
    const auto& root = roots.roots[point.liquid.activity.branch];
    const auto state = th::evaluate_cpa_phase_at_density(
        cpa_physical_test::temperature_k,
        root.molar_density_mol_per_m3,
        point.fractions.liquid,
        model.parameters(), options.phase);
    require(state.association.status == th::CpaAssociationStatus::success &&
                !state.association.sites.empty(),
            "CPA physical VLE did not traverse a converged association state");
    const bool bonded = std::any_of(
        state.association.sites.begin(), state.association.sites.end(),
        [](const th::CpaAssociationSiteState& site) {
            return site.unbonded_fraction < 0.999;
        });
    require(bonded,
            "CPA physical VLE association sites are effectively inactive");
}

fl::CpaPtSplitResult solve_point(
    const cpa_physical_test::ExperimentalVlePoint& experimental,
    bool swapped = false) {
    const auto parameters = cpa_physical_test::parameters(swapped);
    const auto model = th::CpaPtPhase::from_parameters(parameters);
    const auto pt_options = physical_pt_options();
    fl::CpaVleEvaluator evaluator(model, pt_options);
    const auto feed = cpa_physical_test::feed(experimental, swapped);
    const auto starts = cpa_physical_test::starts(experimental, swapped);
    return fl::solve_cpa_pt_vle(
        experimental.pressure_pa, cpa_physical_test::temperature_k,
        feed, evaluator, physical_split_options(), starts, starts);
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
        const auto feed = cpa_physical_test::feed(experimental, false);
        const auto starts = cpa_physical_test::starts(experimental, false);
        const auto result = fl::solve_cpa_pt_vle(
            experimental.pressure_pa, cpa_physical_test::temperature_k,
            feed, evaluator, split_options, starts, starts);
        require(result.solution.status ==
                    fl::PtSplitStatus::two_phase_no_instability_found &&
                    result.solution.candidate() != nullptr &&
                    result.solution.final_stability.has_value() &&
                    result.solution.final_stability->status ==
                        fl::StabilityStatus::no_instability_found,
                "CPA literature methanol-water VLE point did not close as an accepted two-phase state");

        const auto& point = *result.solution.candidate();
        const double liquid = point.fractions.liquid[0];
        const double vapor = point.fractions.vapor[0];
        const double liquid_error = std::abs(liquid - experimental.liquid_methanol);
        const double vapor_error = std::abs(vapor - experimental.vapor_methanol);
        sum_liquid_error += liquid_error;
        sum_vapor_error += vapor_error;
        max_liquid_error = std::max(max_liquid_error, liquid_error);
        max_vapor_error = std::max(max_vapor_error, vapor_error);
        ++accepted;

        require(point.fugacity_norm <= result.solution.options.iteration.fugacity_tolerance,
                "CPA physical VLE accepted point lost fugacity equality");
        require(point.fractions.mass_absolute <=
                    result.solution.options.iteration.mass_absolute_tolerance &&
                    point.fractions.mass_relative <=
                    result.solution.options.iteration.mass_relative_tolerance,
                "CPA physical VLE accepted point lost material balance");
        require_active_association(model, pt_options, experimental.pressure_pa, point);

        std::cout << "CPA_PHYSICAL_VLE PkPa=" << experimental.pressure_pa / 1000.0
                  << " x_exp=" << experimental.liquid_methanol
                  << " x_calc=" << liquid
                  << " y_exp=" << experimental.vapor_methanol
                  << " y_calc=" << vapor
                  << " beta=" << point.fractions.vapor_fraction
                  << '\n';
    }

    require(accepted == cpa_physical_test::points().size(),
            "CPA physical VLE did not evaluate the declared literature points");
    const double mean_liquid_error = sum_liquid_error / static_cast<double>(accepted);
    const double mean_vapor_error = sum_vapor_error / static_cast<double>(accepted);
    std::cout << "CPA_PHYSICAL_VLE_SUMMARY mean_dx=" << mean_liquid_error
              << " mean_dy=" << mean_vapor_error
              << " max_dx=" << max_liquid_error
              << " max_dy=" << max_vapor_error << '\n';

    // Folas/Derawi report a CR-1 correlation at 333.15 K with about 0.6% AAD
    // in bubble pressure and 0.8 percentage points in vapor mole fraction.
    // This inverse PT-flash check is deliberately looser because experimental
    // pressure is imposed and both phase compositions are solved simultaneously.
    require(mean_liquid_error <= 0.025 && mean_vapor_error <= 0.025 &&
                max_liquid_error <= 0.05 && max_vapor_error <= 0.05,
            "CPA literature VLE deviations exceed the declared physical-validation envelope");
}

void component_permutation() {
    const auto& experimental = cpa_physical_test::points()[2];
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
