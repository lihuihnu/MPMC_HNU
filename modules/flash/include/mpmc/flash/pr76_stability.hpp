#ifndef MPMC_FLASH_PR76_STABILITY_HPP
#define MPMC_FLASH_PR76_STABILITY_HPP

#include <mpmc/flash/pt_stability.hpp>
#include <mpmc/thermodynamics/pr76_phase.hpp>

#include <cmath>
#include <cstddef>
#include <span>
#include <string>
#include <utility>
#include <vector>

namespace mpmc::flash {

// Owns the model snapshot; scratch is reused sequentially, not concurrently.
// Neither this adapter nor the generic search duplicates any PR76 EOS formula.
class Pr76StabilityEvaluator {
public:
    explicit Pr76StabilityEvaluator(const thermodynamics::Pr76Phase<double>& model,
                                   thermodynamics::Pr76RootOptions root_options = {})
        : model_(model), root_options_(root_options) {
        if (model_.size() == 0 || root_options_.max_iterations <= 0) {
            throw std::invalid_argument("PR76 stability: invalid model or root iteration limit");
        }
    }
    Pr76StabilityEvaluator(const Pr76StabilityEvaluator&) = delete;
    Pr76StabilityEvaluator& operator=(const Pr76StabilityEvaluator&) = delete;
    Pr76StabilityEvaluator(Pr76StabilityEvaluator&&) = delete;
    Pr76StabilityEvaluator& operator=(Pr76StabilityEvaluator&&) = delete;

    [[nodiscard]] const thermodynamics::Pr76Phase<double>& model() const & noexcept { return model_; }
    const thermodynamics::Pr76Phase<double>& model() const && = delete;
    [[nodiscard]] thermodynamics::Pr76RootOptions root_options() const noexcept { return root_options_; }

    [[nodiscard]] StabilityPhase operator()(double pressure_pa, double temperature_k,
                                           std::span<const double> composition) {
        namespace th = thermodynamics;
        try {
            const auto roots = model_.roots_full(pressure_pa, temperature_k, composition,
                                                 workspace_, root_options_);
            switch (roots.status) {
            case th::Pr76RootStatus::near_multiple:
                throw StabilityPropertyError(StabilityPropertyIssue::root_topology,
                                             "PR76 stability: unresolved root topology");
            case th::Pr76RootStatus::iteration_limit:
                throw StabilityPropertyError(StabilityPropertyIssue::root_iteration_limit,
                                             "PR76 stability: root iteration limit");
            case th::Pr76RootStatus::unrepresentable:
                throw StabilityPropertyError(StabilityPropertyIssue::root_range,
                                             "PR76 stability: nonrepresentable root");
            case th::Pr76RootStatus::success: break;
            }
            std::vector<StabilityPhase> candidates;
            candidates.reserve(roots.count);
            for (std::size_t k = 0; k < roots.count; ++k) {
                // H'>0 <=> dp/dv<0. This is mechanical admissibility, NOT TPD stability.
                if (roots.roots[k].slope_sign <= 0) { continue; }
                auto phase = model_.evaluate_full(pressure_pa, temperature_k, composition,
                                                   k, workspace_, root_options_);
                candidates.push_back({std::move(phase.ln_phi), k, true});
            }
            if (candidates.empty()) {
                throw StabilityPropertyError(StabilityPropertyIssue::no_admissible_branch,
                                             "PR76 stability: no mechanically admissible fluid root");
            }
            const auto difference = [&](const StabilityPhase& a, const StabilityPhase& b) {
                double value = 0, correction = 0;
                for (std::size_t i = 0; i < composition.size(); ++i) {
                    detail::stability_add(composition[i] * (a.ln_phi[i] - b.ln_phi[i]),
                                          value, correction);
                }
                if (!std::isfinite(value)) {
                    throw StabilityPropertyError(StabilityPropertyIssue::nonfinite_properties,
                                                 "PR76 stability: nonrepresentable Gibbs comparison");
                }
                return value;
            };
            std::size_t best = 0;
            for (std::size_t k = 1; k < candidates.size(); ++k) {
                if (difference(candidates[k], candidates[best]) < 0) { best = k; }
            }
            if (!roots.roots[candidates[best].branch].derivative_valid) {
                throw StabilityPropertyError(StabilityPropertyIssue::ill_conditioned_root,
                    "PR76 stability: selected root is too ill-conditioned for a reliable search");
            }
            // Gibbs comparison at fixed composition cancels the ideal/reference terms.
            // A near-tie makes the minimum envelope nonsmooth; do not assert stationarity.
            for (std::size_t k = 0; k < candidates.size(); ++k) {
                if (k == best) { continue; }
                double magnitude = 1;
                for (std::size_t i = 0; i < composition.size(); ++i) {
                    magnitude += composition[i] * (std::abs(candidates[k].ln_phi[i]) +
                                                   std::abs(candidates[best].ln_phi[i]));
                }
                if (!std::isfinite(magnitude)) {
                    throw StabilityPropertyError(StabilityPropertyIssue::nonfinite_properties,
                                                 "PR76 stability: nonrepresentable Gibbs comparison scale");
                }
                if (std::abs(difference(candidates[k], candidates[best])) <=
                    256.0 * detail::stability_eps * magnitude) {
                    candidates[best].smooth = false;
                }
            }
            return std::move(candidates[best]);
        } catch (const th::Pr76PhaseError& error) {
            StabilityPropertyIssue issue = StabilityPropertyIssue::root_range;
            switch (error.code()) {
            case th::Pr76PhaseErrorCode::near_multiple: issue = StabilityPropertyIssue::root_topology; break;
            case th::Pr76PhaseErrorCode::iteration_limit: issue = StabilityPropertyIssue::root_iteration_limit; break;
            case th::Pr76PhaseErrorCode::ill_conditioned_derivative:
                issue = StabilityPropertyIssue::ill_conditioned_root; break;
            case th::Pr76PhaseErrorCode::unrepresentable_root: break;
            }
            throw StabilityPropertyError(issue, error.what());
        } catch (const std::range_error& error) {
            throw StabilityPropertyError(StabilityPropertyIssue::root_range, error.what());
        }
    }
private:
    thermodynamics::Pr76Phase<double> model_;
    thermodynamics::Pr76PhaseWorkspace<double> workspace_;
    thermodynamics::Pr76RootOptions root_options_;
};

struct Pr76StabilityResult {
    StabilityResult search;
    std::string dataset_id, revision;
    std::vector<std::string> component_ids;
    std::string model_profile{thermodynamics::pr76_profile};
    std::string phase_convention{thermodynamics::pr76_pt_convention};
    thermodynamics::Pr76RootOptions root_options;
};

[[nodiscard]] inline Pr76StabilityResult test_pr76_pt_stability(
    double pressure_pa, double temperature_k, std::span<const double> feed,
    Pr76StabilityEvaluator& evaluator, StabilityOptions options = {},
    std::span<const std::vector<double>> extra_starts = {}) {
    if (feed.size() != evaluator.model().size()) {
        throw std::invalid_argument("PR76 stability: feed does not match ordered model snapshot");
    }
    Pr76StabilityResult result;
    const auto& parameters = evaluator.model().parameters();
    result.dataset_id = parameters.dataset_id();
    result.revision = parameters.revision();
    result.root_options = evaluator.root_options();
    for (const auto& component : parameters.components().items()) {
        result.component_ids.push_back(component.id);
    }
    result.search = test_pt_stability(pressure_pa, temperature_k, feed, evaluator, options, extra_starts);
    return result;
}

} // namespace mpmc::flash
#endif // MPMC_FLASH_PR76_STABILITY_HPP
