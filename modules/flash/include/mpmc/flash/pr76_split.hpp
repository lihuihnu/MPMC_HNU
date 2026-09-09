#ifndef MPMC_FLASH_PR76_SPLIT_HPP
#define MPMC_FLASH_PR76_SPLIT_HPP

#include <mpmc/flash/pr76_stability.hpp>
#include <mpmc/flash/pt_split.hpp>
#include <cstddef>
#include <optional>
#include <span>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace mpmc::flash {

// One owned parameter snapshot. Reusable sequential scratch; no simultaneous calls
// on the same instance. The three-argument call is the all-allowed-root stability
// provider; the four-argument call evaluates an explicit candidate branch.
class Pr76VleEvaluator {
public:
    explicit Pr76VleEvaluator(const thermodynamics::Pr76Phase<double>& model,
                              thermodynamics::Pr76RootOptions roots = {})
        : stability_(model, roots) {}
    Pr76VleEvaluator(const Pr76VleEvaluator&) = delete;
    Pr76VleEvaluator& operator=(const Pr76VleEvaluator&) = delete;
    Pr76VleEvaluator(Pr76VleEvaluator&&) = delete;
    Pr76VleEvaluator& operator=(Pr76VleEvaluator&&) = delete;
    [[nodiscard]] const thermodynamics::Pr76Phase<double>& model() const & noexcept {
        return stability_.model();
    }
    const thermodynamics::Pr76Phase<double>& model() const && = delete;
    [[nodiscard]] thermodynamics::Pr76RootOptions root_options() const noexcept {
        return stability_.root_options();
    }
    [[nodiscard]] StabilityPhase operator()(double p, double t, std::span<const double> w) {
        return stability_(p, t, w);
    }
    [[nodiscard]] PtSplitPhase operator()(double p, double t, std::span<const double> w, PtPhaseRole role) {
        namespace th = thermodynamics;
        if (role != PtPhaseRole::liquid_candidate && role != PtPhaseRole::vapor_candidate) {
            throw std::invalid_argument("PR76 VLE: unknown candidate role");
        }
        try {
            const auto roots = model().roots_full(p, t, w, workspace_, root_options());
            switch (roots.status) {
            case th::Pr76RootStatus::near_multiple:
                throw StabilityPropertyError(StabilityPropertyIssue::root_topology, "PR76 VLE: unresolved root topology");
            case th::Pr76RootStatus::iteration_limit:
                throw StabilityPropertyError(StabilityPropertyIssue::root_iteration_limit, "PR76 VLE: root iteration limit");
            case th::Pr76RootStatus::unrepresentable:
                throw StabilityPropertyError(StabilityPropertyIssue::root_range, "PR76 VLE: unrepresentable root");
            case th::Pr76RootStatus::success: break;
            }
            std::optional<std::size_t> selected;
            for (std::size_t k = 0; k < roots.count; ++k) {
                if (roots.roots[k].slope_sign <= 0) { continue; }
                selected = k;
                if (role == PtPhaseRole::liquid_candidate) { break; }
            }
            if (!selected) {
                throw StabilityPropertyError(StabilityPropertyIssue::no_admissible_branch, "PR76 VLE: no admissible root");
            }
            if (!roots.roots[*selected].derivative_valid) {
                throw StabilityPropertyError(StabilityPropertyIssue::ill_conditioned_root, "PR76 VLE: ill-conditioned candidate root");
            }
            auto values = model().evaluate_full(p, t, w, *selected, workspace_, root_options());
            // The smallest/largest admissible root is a REQUESTED candidate role,
            // not proof of phase identity. Distinct density, equations and the
            // all-root common-tangent review are separate acceptance conditions.
            return {{std::move(values.ln_phi), *selected, true}, values.z};
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
    Pr76StabilityEvaluator stability_;
    thermodynamics::Pr76PhaseWorkspace<double> workspace_;
};

struct Pr76PtSplitResult {
    PtSplitResult solution;
    std::string dataset_id, revision;
    std::vector<std::string> component_ids;
    std::string model_profile{thermodynamics::pr76_profile};
    std::string phase_convention{thermodynamics::pr76_pt_convention};
    thermodynamics::Pr76RootOptions root_options;
};

[[nodiscard]] inline Pr76PtSplitResult solve_pr76_pt_vle(
    double pressure_pa, double temperature_k, std::span<const double> feed,
    Pr76VleEvaluator& evaluator, PtSplitOptions options = {},
    std::span<const std::vector<double>> initial_starts = {},
    std::span<const std::vector<double>> final_starts = {}) {
    if (feed.size() != evaluator.model().size()) {
        throw std::invalid_argument("PR76 VLE: feed does not match ordered parameter snapshot");
    }
    Pr76PtSplitResult result;
    const auto& parameters = evaluator.model().parameters();
    result.dataset_id = parameters.dataset_id();
    result.revision = parameters.revision();
    result.root_options = evaluator.root_options();
    for (const auto& component : parameters.components().items()) {
        result.component_ids.push_back(component.id);
    }
    result.solution = solve_pt_vle(pressure_pa, temperature_k, feed, evaluator, evaluator,
                                   options, initial_starts, final_starts);
    return result;
}

} // namespace mpmc::flash
#endif // MPMC_FLASH_PR76_SPLIT_HPP
