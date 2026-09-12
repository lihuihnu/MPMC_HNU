#ifndef MPMC_FLASH_CPA_STABILITY_HPP
#define MPMC_FLASH_CPA_STABILITY_HPP

#include <mpmc/flash/pt_stability.hpp>
#include <mpmc/thermodynamics/cpa_pt_phase.hpp>

#include <cmath>
#include <cstddef>
#include <span>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace mpmc::flash {

inline constexpr const char* cpa_pt_stability_convention =
    "CPA/PT/stability/minimum-Gibbs-mechanically-admissible-root/v1";

// Owns the configured CPA PT model. Each call performs a fresh bounded density
// root search at the supplied composition; no phase/root identity is persisted
// across calls. Instances are intended for sequential reuse, not concurrent use.
class CpaStabilityEvaluator {
public:
    explicit CpaStabilityEvaluator(
        const thermodynamics::CpaPtPhase& model,
        thermodynamics::CpaPtOptions pt_options = {})
        : model_(model), pt_options_(pt_options) {
        if (model_.size() == 0U) {
            throw std::invalid_argument("CPA stability: empty configured model");
        }
        // Validate the numerical contract at construction rather than allowing
        // an invalid option set to fail only after a TPD search has started.
        thermodynamics::cpa_detail::validate_cpa_pt_options(pt_options_);
    }

    CpaStabilityEvaluator(const CpaStabilityEvaluator&) = delete;
    CpaStabilityEvaluator& operator=(const CpaStabilityEvaluator&) = delete;
    CpaStabilityEvaluator(CpaStabilityEvaluator&&) = delete;
    CpaStabilityEvaluator& operator=(CpaStabilityEvaluator&&) = delete;

    [[nodiscard]] const thermodynamics::CpaPtPhase& model() const & noexcept {
        return model_;
    }
    const thermodynamics::CpaPtPhase& model() const && = delete;

    [[nodiscard]] const thermodynamics::CpaPtOptions& pt_options() const & noexcept {
        return pt_options_;
    }
    const thermodynamics::CpaPtOptions& pt_options() const && = delete;

    [[nodiscard]] StabilityPhase operator()(
        double pressure_pa, double temperature_k,
        std::span<const double> composition) {
        namespace th = thermodynamics;
        const auto roots = model_.roots(
            pressure_pa, temperature_k, composition, pt_options_);
        switch (roots.status) {
        case th::CpaPtRootStatus::success:
            break;
        case th::CpaPtRootStatus::near_multiple:
            throw StabilityPropertyError(
                StabilityPropertyIssue::root_topology,
                "CPA stability: finite density search retained unresolved near-multiple topology");
        case th::CpaPtRootStatus::no_root:
            throw StabilityPropertyError(
                StabilityPropertyIssue::root_range,
                "CPA stability: finite density search found no representable pressure root");
        case th::CpaPtRootStatus::property_failure:
            throw StabilityPropertyError(
                StabilityPropertyIssue::root_range,
                "CPA stability: density-state property evaluation failed: " +
                    roots.diagnostic);
        case th::CpaPtRootStatus::iteration_limit:
        case th::CpaPtRootStatus::evaluation_limit:
            throw StabilityPropertyError(
                StabilityPropertyIssue::root_iteration_limit,
                "CPA stability: bounded density-root search did not complete: " +
                    roots.diagnostic);
        case th::CpaPtRootStatus::root_limit:
            throw StabilityPropertyError(
                StabilityPropertyIssue::root_topology,
                "CPA stability: root-count quota prevented topology resolution");
        }

        std::vector<std::size_t> admissible;
        admissible.reserve(roots.roots.size());
        for (std::size_t k = 0; k < roots.roots.size(); ++k) {
            // Positive dP/drho is equivalent to negative dP/dv at fixed n,T.
            // This is only mechanical admissibility; TPD stability is separate.
            if (roots.roots[k].pressure_slope_sign > 0) {
                admissible.push_back(k);
            }
        }
        if (admissible.empty()) {
            throw StabilityPropertyError(
                StabilityPropertyIssue::no_admissible_branch,
                "CPA stability: no mechanically admissible simple density root");
        }

        const auto gibbs_offset = [&](std::size_t root_index) {
            const auto& root = roots.roots[root_index];
            if (root.ln_phi.size() != composition.size()) {
                throw StabilityPropertyError(
                    StabilityPropertyIssue::nonfinite_properties,
                    "CPA stability: root fugacity dimension mismatch");
            }
            double value = 0.0;
            double correction = 0.0;
            for (std::size_t i = 0; i < composition.size(); ++i) {
                detail::stability_add(
                    composition[i] * root.ln_phi[i], value, correction);
            }
            if (!std::isfinite(value)) {
                throw StabilityPropertyError(
                    StabilityPropertyIssue::nonfinite_properties,
                    "CPA stability: nonrepresentable Gibbs root comparison");
            }
            return value;
        };

        std::size_t best_position = 0U;
        double best_gibbs = gibbs_offset(admissible[0]);
        for (std::size_t position = 1U; position < admissible.size(); ++position) {
            const double value = gibbs_offset(admissible[position]);
            if (value < best_gibbs) {
                best_position = position;
                best_gibbs = value;
            }
        }

        const std::size_t best_root = admissible[best_position];
        StabilityPhase phase{
            roots.roots[best_root].ln_phi,
            best_root,
            true};

        // The same-composition ideal/reference terms cancel. If two admissible
        // CPA roots have indistinguishable Gibbs offsets, the lower envelope is
        // nonsmooth even though each density root is individually simple.
        for (const std::size_t candidate : admissible) {
            if (candidate == best_root) { continue; }
            const double candidate_gibbs = gibbs_offset(candidate);
            double magnitude = 1.0;
            double correction = 0.0;
            for (std::size_t i = 0; i < composition.size(); ++i) {
                detail::stability_add(
                    composition[i] *
                        (std::abs(roots.roots[candidate].ln_phi[i]) +
                         std::abs(roots.roots[best_root].ln_phi[i])),
                    magnitude, correction);
            }
            if (!std::isfinite(magnitude)) {
                throw StabilityPropertyError(
                    StabilityPropertyIssue::nonfinite_properties,
                    "CPA stability: nonrepresentable Gibbs tie scale");
            }
            if (std::abs(candidate_gibbs - best_gibbs) <=
                256.0 * detail::stability_eps * magnitude) {
                phase.smooth = false;
            }
        }
        return phase;
    }

private:
    thermodynamics::CpaPtPhase model_;
    thermodynamics::CpaPtOptions pt_options_;
};

struct CpaStabilityResult {
    StabilityResult search;
    std::string dataset_id;
    std::string revision;
    std::vector<std::string> component_ids;
    std::string model_profile{std::string(thermodynamics::cpa_profile)};
    std::string phase_convention{std::string(thermodynamics::cpa_pt_convention)};
    std::string stability_convention{cpa_pt_stability_convention};
    thermodynamics::CpaPtOptions pt_options;
};

[[nodiscard]] inline CpaStabilityResult test_cpa_pt_stability(
    double pressure_pa, double temperature_k,
    std::span<const double> feed,
    CpaStabilityEvaluator& evaluator,
    StabilityOptions options = {},
    std::span<const std::vector<double>> extra_starts = {}) {
    if (feed.size() != evaluator.model().size()) {
        throw std::invalid_argument(
            "CPA stability: feed does not match ordered model snapshot");
    }

    CpaStabilityResult result;
    const auto& parameters = evaluator.model().parameters();
    result.dataset_id = parameters.dataset_id();
    result.revision = parameters.revision();
    result.pt_options = evaluator.pt_options();
    result.component_ids.reserve(parameters.size());
    for (const auto& component : parameters.components().items()) {
        result.component_ids.push_back(component.id);
    }
    result.search = test_pt_stability(
        pressure_pa, temperature_k, feed,
        evaluator, options, extra_starts);
    return result;
}

} // namespace mpmc::flash

#endif // MPMC_FLASH_CPA_STABILITY_HPP
