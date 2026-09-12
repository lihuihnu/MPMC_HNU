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
#include <string_view>
#include <utility>
#include <vector>

namespace mpmc::flash {

inline constexpr std::string_view cpa_pt_vle_convention =
    "CPA/PT/VLE/minGibbs-stability-explicit-density-sides/logK-SSI-RR/v1";

// Numerical density-root sides only. They are not liquid/vapor morphology labels
// and root vector indices are never persisted as phase identity across states.
enum class CpaRootSide {
    lower_density_admissible,
    upper_density_admissible
};

namespace detail {

inline void cpa_flash_require_root_success(
    const thermodynamics::CpaPtRootSet& roots,
    std::string_view context) {
    namespace th = thermodynamics;
    switch (roots.status) {
    case th::CpaPtRootStatus::success:
        return;
    case th::CpaPtRootStatus::near_multiple:
        throw StabilityPropertyError(
            StabilityPropertyIssue::root_topology,
            std::string(context) + ": unresolved near-multiple density-root topology");
    case th::CpaPtRootStatus::no_root:
        throw StabilityPropertyError(
            StabilityPropertyIssue::root_range,
            std::string(context) + ": no representable density root");
    case th::CpaPtRootStatus::property_failure:
        throw StabilityPropertyError(
            StabilityPropertyIssue::root_range,
            std::string(context) + ": density-state property failure: " +
                roots.diagnostic);
    case th::CpaPtRootStatus::iteration_limit:
    case th::CpaPtRootStatus::evaluation_limit:
        throw StabilityPropertyError(
            StabilityPropertyIssue::root_iteration_limit,
            std::string(context) + ": bounded density-root search did not complete: " +
                roots.diagnostic);
    case th::CpaPtRootStatus::root_limit:
        throw StabilityPropertyError(
            StabilityPropertyIssue::root_topology,
            std::string(context) + ": root-count quota prevented topology resolution");
    }
    throw StabilityPropertyError(
        StabilityPropertyIssue::root_topology,
        std::string(context) + ": unknown CPA root status");
}

[[nodiscard]] inline std::vector<std::size_t> cpa_flash_admissible_roots(
    const thermodynamics::CpaPtRootSet& roots) {
    std::vector<std::size_t> admissible;
    admissible.reserve(roots.roots.size());
    for (std::size_t k = 0; k < roots.roots.size(); ++k) {
        if (roots.roots[k].pressure_slope_sign > 0) {
            admissible.push_back(k);
        }
    }
    if (admissible.empty()) {
        throw StabilityPropertyError(
            StabilityPropertyIssue::no_admissible_branch,
            "CPA flash: no mechanically admissible simple density root");
    }
    return admissible;
}

[[nodiscard]] inline std::size_t cpa_flash_select_root_side(
    const thermodynamics::CpaPtRootSet& roots,
    std::span<const std::size_t> admissible,
    CpaRootSide side) {
    if (admissible.empty()) {
        throw StabilityPropertyError(
            StabilityPropertyIssue::no_admissible_branch,
            "CPA flash: empty admissible-root set");
    }
    const auto less_density = [&](std::size_t first, std::size_t second) {
        return roots.roots[first].molar_density_mol_per_m3 <
               roots.roots[second].molar_density_mol_per_m3;
    };
    if (side == CpaRootSide::lower_density_admissible) {
        return *std::min_element(admissible.begin(), admissible.end(), less_density);
    }
    return *std::max_element(admissible.begin(), admissible.end(), less_density);
}

[[nodiscard]] inline PtSplitPhase cpa_flash_evaluate_side(
    const thermodynamics::CpaPtPhase& model,
    const thermodynamics::CpaPtOptions& pt_options,
    double pressure_pa,
    double temperature_k,
    std::span<const double> composition,
    CpaRootSide side) {
    const auto roots = model.roots(
        pressure_pa, temperature_k, composition, pt_options);
    cpa_flash_require_root_success(roots, "CPA split");
    const auto admissible = cpa_flash_admissible_roots(roots);
    const std::size_t selected = cpa_flash_select_root_side(
        roots, admissible, side);
    const auto& root = roots.roots[selected];
    if (!std::isfinite(root.compressibility_factor) ||
        !(root.compressibility_factor > 0.0) ||
        root.ln_phi.size() != composition.size()) {
        throw StabilityPropertyError(
            StabilityPropertyIssue::nonfinite_properties,
            "CPA split: selected density root has invalid phase properties");
    }
    for (const double value : root.ln_phi) {
        if (!std::isfinite(value)) {
            throw StabilityPropertyError(
                StabilityPropertyIssue::nonfinite_properties,
                "CPA split: selected density root has nonfinite ln(phi)");
        }
    }
    return {{root.ln_phi, selected, true}, root.compressibility_factor};
}

} // namespace detail

// Stability uses the minimum-Gibbs mechanically admissible root at every trial
// composition. Fixed two-phase iteration instead requests opposite density-root
// sides so two phase slots do not collapse onto the same minimum-Gibbs envelope.
// The requested sides are numerical continuation choices only; final all-root
// common-tangent stability remains the topology acceptance gate.
class CpaVleEvaluator {
public:
    explicit CpaVleEvaluator(
        const thermodynamics::CpaPtPhase& model,
        thermodynamics::CpaPtOptions pt_options = {})
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
        if (role == PtPhaseRole::liquid_candidate) {
            return detail::cpa_flash_evaluate_side(
                model(), pt_options(), pressure_pa, temperature_k, composition,
                CpaRootSide::upper_density_admissible);
        }
        if (role == PtPhaseRole::vapor_candidate) {
            return detail::cpa_flash_evaluate_side(
                model(), pt_options(), pressure_pa, temperature_k, composition,
                CpaRootSide::lower_density_admissible);
        }
        throw std::invalid_argument("CPA VLE: unknown candidate role");
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
    std::string split_convention{std::string(cpa_pt_vle_convention)};
    thermodynamics::CpaPtOptions pt_options;
};

[[nodiscard]] inline CpaPtSplitResult solve_cpa_pt_vle(
    double pressure_pa, double temperature_k,
    std::span<const double> feed,
    CpaVleEvaluator& evaluator,
    PtSplitOptions options = {},
    std::span<const std::vector<double>> initial_starts = {},
    std::span<const std::vector<double>> final_starts = {}) {
    if (feed.size() != evaluator.model().size()) {
        throw std::invalid_argument(
            "CPA VLE: feed does not match ordered model snapshot");
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
        evaluator, evaluator, options, initial_starts, final_starts);
    return result;
}

} // namespace mpmc::flash

#endif // MPMC_FLASH_CPA_SPLIT_HPP
