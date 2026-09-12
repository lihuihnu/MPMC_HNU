#include <mpmc/flash/cpa_split.hpp>

#include "test_support.hpp"
#include "../../thermodynamics/cpa_pt_phase/test_support.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <string_view>
#include <vector>

namespace {
namespace fl = mpmc::flash;
namespace th = mpmc::thermodynamics;
using Vec = std::vector<double>;

void require(bool value, const char* message) {
    if (!value) { throw std::runtime_error(message); }
}

void explicit_density_sides() {
    const auto parameters = cpa_stability_test::pure();
    const auto model = th::CpaPtPhase::from_parameters(parameters);
    fl::CpaVleEvaluator evaluator(model);
    const Vec x{1.0};
    const auto roots = model.roots(1.0e5, 200.0, x);
    require(roots.status == th::CpaPtRootStatus::success && roots.roots.size() == 3U,
            "CPA explicit-side fixture lost its three density roots");

    const auto liquid = evaluator(
        1.0e5, 200.0, x, fl::PtPhaseRole::liquid_candidate);
    const auto vapor = evaluator(
        1.0e5, 200.0, x, fl::PtPhaseRole::vapor_candidate);
    require(liquid.activity.branch != vapor.activity.branch,
            "CPA explicit split sides collapsed onto the same root in a three-root state");
    require(vapor.z > liquid.z,
            "CPA lower-density candidate does not have larger compressibility factor");
    require(roots.roots[liquid.activity.branch].molar_density_mol_per_m3 >
                roots.roots[vapor.activity.branch].molar_density_mol_per_m3,
            "CPA split role-to-density-side mapping changed");
}

fl::CpaPtSplitResult run_structural_vle(bool swapped) {
    const auto parameters = cpa_stability_test::binary(swapped);
    const auto model = th::CpaPtPhase::from_parameters(parameters);
    fl::CpaVleEvaluator evaluator(model, cpa_stability_test::fast_pt_options());
    const auto feed = cpa_split_test::feed(swapped);
    const auto starts = cpa_split_test::starts(swapped);
    fl::PtSplitOptions options;
    options.initial_stability.automatic_starts = false;
    options.final_stability.automatic_starts = false;
    return fl::solve_cpa_pt_vle(
        cpa_split_test::pressure_pa, cpa_split_test::temperature_k,
        feed, evaluator, options, starts, starts);
}

void structural_anchor_residual() {
    const auto parameters = cpa_stability_test::binary(false);
    const auto model = th::CpaPtPhase::from_parameters(parameters);
    fl::CpaVleEvaluator evaluator(model, cpa_stability_test::fast_pt_options());
    const auto& phases = cpa_split_test::coexistence_phases();
    const auto liquid = evaluator(
        cpa_split_test::pressure_pa, cpa_split_test::temperature_k,
        phases[0], fl::PtPhaseRole::liquid_candidate);
    const auto vapor = evaluator(
        cpa_split_test::pressure_pa, cpa_split_test::temperature_k,
        phases[1], fl::PtPhaseRole::vapor_candidate);
    double norm = 0.0;
    for (std::size_t i = 0; i < phases[0].size(); ++i) {
        const double residual =
            std::log(phases[0][i]) + liquid.activity.ln_phi[i] -
            std::log(phases[1][i]) - vapor.activity.ln_phi[i];
        norm = std::max(norm, std::abs(residual));
    }
    require(norm < 1.0e-9 && vapor.z > liquid.z,
            "CPA structural VLE anchors no longer form the intended common-tangent pair");
}

void binary_vle_baseline() {
    const auto result = run_structural_vle(false);
    require(result.solution.initial_stability.status == fl::StabilityStatus::unstable,
            "CPA VLE structural feed no longer has instability evidence");
    require(result.solution.status == fl::PtSplitStatus::two_phase_no_instability_found &&
                result.solution.candidate() != nullptr &&
                result.solution.final_stability.has_value() &&
                result.solution.final_stability->status ==
                    fl::StabilityStatus::no_instability_found,
            "CPA VLE baseline did not close through split and final common-tangent review");
    const auto& point = *result.solution.candidate();
    require(point.fugacity_norm <= result.solution.options.iteration.fugacity_tolerance &&
                point.fractions.mass_absolute <=
                    result.solution.options.iteration.mass_absolute_tolerance &&
                point.fractions.mass_relative <=
                    result.solution.options.iteration.mass_relative_tolerance,
            "CPA VLE baseline lost fugacity equality or material balance");
    require(point.vapor.z > point.liquid.z &&
                std::abs(point.fractions.vapor_fraction -
                         cpa_split_test::lower_density_phase_fraction) < 2.0e-8,
            "CPA VLE baseline lost density ordering or structural phase fraction");
    const auto& expected = cpa_split_test::coexistence_phases();
    require(std::abs(point.fractions.liquid[0] - expected[0][0]) < 2.0e-8 &&
                std::abs(point.fractions.vapor[0] - expected[1][0]) < 2.0e-8,
            "CPA VLE baseline left the structural coexistence anchors");
    require(result.dataset_id == "synthetic-cpa-stability-binary" &&
                result.component_ids == std::vector<std::string>({"A", "B"}) &&
                result.split_convention == std::string(fl::cpa_pt_vle_convention) &&
                result.stability_convention == std::string(fl::cpa_pt_stability_convention),
            "CPA VLE result lost model/provenance identity");
}

void component_permutation() {
    const auto first = run_structural_vle(false);
    const auto second = run_structural_vle(true);
    require(first.solution.status == fl::PtSplitStatus::two_phase_no_instability_found &&
                second.solution.status == fl::PtSplitStatus::two_phase_no_instability_found &&
                first.solution.candidate() && second.solution.candidate(),
            "CPA VLE component permutation changed accepted topology");
    const auto& a = *first.solution.candidate();
    const auto& b = *second.solution.candidate();
    require(std::abs(a.fractions.vapor_fraction - b.fractions.vapor_fraction) < 2.0e-10,
            "CPA VLE phase fraction changed under component permutation");
    require(std::abs(a.fractions.liquid[0] - b.fractions.liquid[1]) < 2.0e-9 &&
                std::abs(a.fractions.liquid[1] - b.fractions.liquid[0]) < 2.0e-9 &&
                std::abs(a.fractions.vapor[0] - b.fractions.vapor[1]) < 2.0e-9 &&
                std::abs(a.fractions.vapor[1] - b.fractions.vapor[0]) < 2.0e-9,
            "CPA VLE phase compositions changed under runtime component permutation");
}

void associating_phase_provider_smoke() {
    const auto parameters = cpa_pt_test::associating_binary(false);
    require(parameters.association_pair(
                parameters.components().index_of("A"), "H",
                parameters.components().index_of("A"), "H") != nullptr,
            "CPA associating split fixture lost its association contract");
    const auto model = th::CpaPtPhase::from_parameters(parameters);
    fl::CpaVleEvaluator evaluator(model);
    const Vec x{0.7, 0.3};
    const auto roots = model.roots(1.0e6, 330.0, x);
    require(roots.status == th::CpaPtRootStatus::success && !roots.roots.empty(),
            "associating CPA split smoke state has no resolved root");
    const auto low = evaluator(
        1.0e6, 330.0, x, fl::PtPhaseRole::vapor_candidate);
    const auto high = evaluator(
        1.0e6, 330.0, x, fl::PtPhaseRole::liquid_candidate);
    require(low.activity.ln_phi.size() == x.size() &&
                high.activity.ln_phi.size() == x.size(),
            "associating CPA split evaluator lost fugacity dimension");
    for (double value : low.activity.ln_phi) {
        require(std::isfinite(value), "nonfinite associating lower-density ln(phi)");
    }
    for (double value : high.activity.ln_phi) {
        require(std::isfinite(value), "nonfinite associating upper-density ln(phi)");
    }
}

using Test = std::pair<std::string_view, void (*)()>;
constexpr Test tests[]{
    {"explicit_density_sides", explicit_density_sides},
    {"structural_anchor_residual", structural_anchor_residual},
    {"binary_vle", binary_vle_baseline},
    {"component_permutation", component_permutation},
    {"associating_provider", associating_phase_provider_smoke}};

} // namespace

int main(int argc, char** argv) {
    try {
        if (argc != 2) { throw std::invalid_argument("one test name required"); }
        for (const auto& [name, run] : tests) {
            if (name == argv[1]) {
                run();
                return 0;
            }
        }
        throw std::invalid_argument("unknown test name");
    } catch (const std::exception&) {
        return 1;
    }
}
