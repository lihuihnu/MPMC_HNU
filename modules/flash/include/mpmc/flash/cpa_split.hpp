#ifndef MPMC_FLASH_CPA_SPLIT_HPP
#define MPMC_FLASH_CPA_SPLIT_HPP

#include <mpmc/flash/cpa_stability.hpp>
#include <mpmc/flash/pt_split.hpp>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <optional>
#include <span>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace mpmc::flash {

inline constexpr const char* cpa_pt_vle_convention =
    "CPA/PT/two-phase/density-side-logK-RR-common-tangent/v1";

// The generic PT split requires 1e-11 log-fugacity closure. CPA density roots
// therefore need materially tighter pressure closure than the standalone
// phase-property default; otherwise root error becomes the outer residual floor.
[[nodiscard]] inline thermodynamics::CpaPtOptions
cpa_pt_vle_default_phase_options() {
    thermodynamics::CpaPtOptions options;
    options.pressure_absolute_tolerance_pa = 1.0e-7;
    options.pressure_relative_tolerance = 1.0e-12;
    return options;
}

// Keep the generic global split budget unchanged, but cap one numerical role
// assignment so a wrong density-side initialization cannot consume all provider
// calls before the alternate assignment is attempted.
[[nodiscard]] inline PtSplitOptions cpa_pt_vle_default_split_options() {
    PtSplitOptions options;
    options.max_evaluations_per_attempt = 8192U;
    return options;
}

// The three-argument call is the all-admissible-root minimum-Gibbs stability
// provider. The four-argument call requests one numerical density side for the
// split iteration. These candidate roles are not physical morphology labels.
class CpaVleEvaluator {
public:
    explicit CpaVleEvaluator(
        const thermodynamics::CpaPtPhase& model,
        thermodynamics::CpaPtOptions pt_options =
            cpa_pt_vle_default_phase_options())
        : stability_(model, std::move(pt_options)) {}

    CpaVleEvaluator(const CpaVleEvaluator&) = delete;
    CpaVleEvaluator& operator=(const CpaVleEvaluator&) = delete;
    CpaVleEvaluator(CpaVleEvaluator&&) = delete;
    CpaVleEvaluator& operator=(CpaVleEvaluator&&) = delete;

    [[nodiscard]] const thermodynamics::CpaPtPhase& model() const & noexcept {
        return stability_.model();
    }
    const thermodynamics::CpaPtPhase& model() const && = delete;

    [[nodiscard]] const thermodynamics::CpaPtOptions& pt_options() const & noexcept {
        return stability_.pt_options();
    }
    const thermodynamics::CpaPtOptions& pt_options() const && = delete;

    [[nodiscard]] StabilityPhase operator()(
        double pressure_pa, double temperature_k,
        std::span<const double> composition) {
        return stability_(pressure_pa, temperature_k, composition);
    }

    [[nodiscard]] PtSplitPhase operator()(
        double pressure_pa, double temperature_k,
        std::span<const double> composition,
        PtPhaseRole role) {
        namespace th = thermodynamics;
        if (role != PtPhaseRole::liquid_candidate &&
            role != PtPhaseRole::vapor_candidate) {
            throw std::invalid_argument("CPA VLE: unknown candidate role");
        }

        const auto roots = model().roots(
            pressure_pa, temperature_k, composition, pt_options());
        switch (roots.status) {
        case th::CpaPtRootStatus::success:
            break;
        case th::CpaPtRootStatus::near_multiple:
            throw StabilityPropertyError(
                StabilityPropertyIssue::root_topology,
                "CPA VLE: unresolved tangent/near-multiple density-root topology");
        case th::CpaPtRootStatus::no_root:
            throw StabilityPropertyError(
                StabilityPropertyIssue::root_range,
                "CPA VLE: finite density search found no representable root");
        case th::CpaPtRootStatus::property_failure:
            throw StabilityPropertyError(
                StabilityPropertyIssue::root_range,
                "CPA VLE: density-state property failure: " + roots.diagnostic);
        case th::CpaPtRootStatus::iteration_limit:
        case th::CpaPtRootStatus::evaluation_limit:
            throw StabilityPropertyError(
                StabilityPropertyIssue::root_iteration_limit,
                "CPA VLE: bounded density-root search did not complete: " +
                    roots.diagnostic);
        case th::CpaPtRootStatus::root_limit:
            throw StabilityPropertyError(
                StabilityPropertyIssue::root_topology,
                "CPA VLE: root-count quota prevented topology resolution");
        }

        std::optional<std::size_t> selected;
        for (std::size_t k = 0; k < roots.roots.size(); ++k) {
            if (roots.roots[k].pressure_slope_sign <= 0) { continue; }
            if (!selected) {
                selected = k;
                continue;
            }
            const bool denser =
                roots.roots[k].molar_density_mol_per_m3 >
                roots.roots[*selected].molar_density_mol_per_m3;
            const bool lighter =
                roots.roots[k].molar_density_mol_per_m3 <
                roots.roots[*selected].molar_density_mol_per_m3;
            if ((role == PtPhaseRole::liquid_candidate && denser) ||
                (role == PtPhaseRole::vapor_candidate && lighter)) {
                selected = k;
            }
        }
        if (!selected) {
            throw StabilityPropertyError(
                StabilityPropertyIssue::no_admissible_branch,
                "CPA VLE: no mechanically admissible simple density root");
        }

        const auto& root = roots.roots[*selected];
        if (!std::isfinite(root.compressibility_factor) ||
            !(root.compressibility_factor > 0.0) ||
            root.ln_phi.size() != composition.size()) {
            throw StabilityPropertyError(
                StabilityPropertyIssue::nonfinite_properties,
                "CPA VLE: selected root has invalid phase-property payload");
        }
        return {{root.ln_phi, *selected, true}, root.compressibility_factor};
    }

private:
    CpaStabilityEvaluator stability_;
};

struct CpaPtSplitResult {
    PtSplitResult solution;
    std::string dataset_id;
    std::string revision;
    std::vector<std::string> component_ids;
    std::string model_profile{std::string(thermodynamics::cpa_profile)};
    std::string phase_convention{std::string(thermodynamics::cpa_pt_convention)};
    std::string stability_convention{cpa_pt_stability_convention};
    std::string split_convention{cpa_pt_vle_convention};
    thermodynamics::CpaPtOptions pt_options;
};

[[nodiscard]] inline CpaPtSplitResult solve_cpa_pt_vle(
    double pressure_pa, double temperature_k,
    std::span<const double> feed,
    CpaVleEvaluator& evaluator,
    PtSplitOptions options = cpa_pt_vle_default_split_options(),
    std::span<const std::vector<double>> initial_starts = {},
    std::span<const std::vector<double>> final_starts = {}) {
    if (feed.size() != evaluator.model().size()) {
        throw std::invalid_argument(
            "CPA VLE: feed does not match ordered parameter snapshot");
    }

    CpaPtSplitResult result;
    const auto& parameters = evaluator.model().parameters();
    result.dataset_id = parameters.dataset_id();
    result.revision = parameters.revision();
    result.pt_options = evaluator.pt_options();
    result.component_ids.reserve(parameters.size());
    for (const auto& component : parameters.components().items()) {
        result.component_ids.push_back(component.id);
    }
    result.solution = solve_pt_vle(
        pressure_pa, temperature_k, feed,
        evaluator, evaluator,
        options, initial_starts, final_starts);
    return result;
}

} // namespace mpmc::flash

#endif // MPMC_FLASH_CPA_SPLIT_HPP
