#ifndef MPMC_FLASH_SW92_STABILITY_HPP
#define MPMC_FLASH_SW92_STABILITY_HPP

#include <mpmc/flash/pt_stability.hpp>
#include <mpmc/thermodynamics/sw92_phase.hpp>

#include <cmath>
#include <cstddef>
#include <span>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace mpmc::flash {

/// Fixed-family SW92 property provider for the generic PT stability search.
/// The phase family and NaCl molality are constructor state and never change
/// during one search. `StabilityPhase::branch` is only the increasing-Z cubic
/// root index inside that family; it does not encode aqueous/nonaqueous identity.
class Sw92FamilyStabilityEvaluator {
public:
    Sw92FamilyStabilityEvaluator(
        const thermodynamics::Sw92Phase<double>& model,
        double nacl_molality_mol_per_kg_water,
        thermodynamics::SwPhaseFamily family,
        thermodynamics::Sw92RootOptions root_options = {})
        : model_(model), molality_(nacl_molality_mol_per_kg_water),
          family_(family), root_options_(root_options) {
        if (model_.size() == 0 || root_options_.max_iterations <= 0) {
            throw std::invalid_argument("SW92 stability: invalid model or root iteration limit");
        }
        switch (family_) {
        case thermodynamics::SwPhaseFamily::aqueous:
        case thermodynamics::SwPhaseFamily::nonaqueous:
            break;
        default:
            throw std::invalid_argument("SW92 stability: unknown phase family");
        }
        if (!std::isfinite(molality_) || molality_ < 0.0) {
            throw std::domain_error("SW92 stability: finite NaCl molality >=0 mol/kg H2O required");
        }
        if (const auto& bounds = model_.parameters().applicability().nacl_molality_mol_per_kg_water;
            bounds && (static_cast<long double>(molality_) < bounds->lower ||
                       static_cast<long double>(molality_) > bounds->upper)) {
            throw std::domain_error("SW92 stability: NaCl molality outside declared dataset interval");
        }
    }

    Sw92FamilyStabilityEvaluator(const Sw92FamilyStabilityEvaluator&) = delete;
    Sw92FamilyStabilityEvaluator& operator=(const Sw92FamilyStabilityEvaluator&) = delete;
    Sw92FamilyStabilityEvaluator(Sw92FamilyStabilityEvaluator&&) = delete;
    Sw92FamilyStabilityEvaluator& operator=(Sw92FamilyStabilityEvaluator&&) = delete;

    [[nodiscard]] const thermodynamics::Sw92Phase<double>& model() const & noexcept {
        return model_;
    }
    const thermodynamics::Sw92Phase<double>& model() const && = delete;
    [[nodiscard]] double nacl_molality_mol_per_kg_water() const noexcept { return molality_; }
    [[nodiscard]] thermodynamics::SwPhaseFamily family() const noexcept { return family_; }
    [[nodiscard]] thermodynamics::Sw92RootOptions root_options() const noexcept {
        return root_options_;
    }

    [[nodiscard]] StabilityPhase operator()(double pressure_pa, double temperature_k,
                                           std::span<const double> composition) {
        namespace th = thermodynamics;
        try {
            const auto roots = model_.roots(pressure_pa, temperature_k, composition,
                                            molality_, family_, workspace_, root_options_);
            switch (roots.status) {
            case th::Sw92RootStatus::near_multiple:
                throw StabilityPropertyError(StabilityPropertyIssue::root_topology,
                    "SW92 stability: unresolved root topology");
            case th::Sw92RootStatus::iteration_limit:
                throw StabilityPropertyError(StabilityPropertyIssue::root_iteration_limit,
                    "SW92 stability: root iteration limit");
            case th::Sw92RootStatus::unrepresentable:
                throw StabilityPropertyError(StabilityPropertyIssue::root_range,
                    "SW92 stability: nonrepresentable root");
            case th::Sw92RootStatus::success:
                break;
            }

            std::vector<StabilityPhase> candidates;
            candidates.reserve(roots.count);
            for (std::size_t root_index = 0; root_index < roots.count; ++root_index) {
                // H'>0 <=> dp/dv<0. This mechanical filter is not a TPD test.
                if (roots.roots[root_index].slope_sign <= 0) { continue; }
                auto values = model_.evaluate(pressure_pa, temperature_k, composition,
                                              molality_, family_, root_index,
                                              workspace_, root_options_);
                candidates.push_back({std::move(values.ln_phi), root_index, true});
            }
            if (candidates.empty()) {
                throw StabilityPropertyError(StabilityPropertyIssue::no_admissible_branch,
                    "SW92 stability: no mechanically admissible fluid root in selected family");
            }

            // At fixed p,T,x,family, ideal/reference terms cancel between cubic
            // roots, so sum_i x_i ln(phi_i) ranks their molar Gibbs energies.
            const auto difference = [&](const StabilityPhase& a, const StabilityPhase& b) {
                double value = 0.0;
                double correction = 0.0;
                for (std::size_t i = 0; i < composition.size(); ++i) {
                    detail::stability_add(composition[i] * (a.ln_phi[i] - b.ln_phi[i]),
                                          value, correction);
                }
                if (!std::isfinite(value)) {
                    throw StabilityPropertyError(StabilityPropertyIssue::nonfinite_properties,
                        "SW92 stability: nonrepresentable same-family Gibbs comparison");
                }
                return value;
            };

            std::size_t best = 0;
            for (std::size_t candidate = 1; candidate < candidates.size(); ++candidate) {
                if (difference(candidates[candidate], candidates[best]) < 0.0) {
                    best = candidate;
                }
            }
            if (!roots.roots[candidates[best].branch].derivative_valid) {
                throw StabilityPropertyError(StabilityPropertyIssue::ill_conditioned_root,
                    "SW92 stability: selected root is too ill-conditioned for a reliable search");
            }

            // A near-tie makes the same-family minimum envelope nonsmooth.
            for (std::size_t candidate = 0; candidate < candidates.size(); ++candidate) {
                if (candidate == best) { continue; }
                double magnitude = 1.0;
                for (std::size_t i = 0; i < composition.size(); ++i) {
                    magnitude += composition[i] *
                        (std::abs(candidates[candidate].ln_phi[i]) +
                         std::abs(candidates[best].ln_phi[i]));
                }
                if (!std::isfinite(magnitude)) {
                    throw StabilityPropertyError(StabilityPropertyIssue::nonfinite_properties,
                        "SW92 stability: nonrepresentable Gibbs comparison scale");
                }
                if (std::abs(difference(candidates[candidate], candidates[best])) <=
                    256.0 * detail::stability_eps * magnitude) {
                    candidates[best].smooth = false;
                }
            }
            return std::move(candidates[best]);
        } catch (const th::Sw92PhaseError& error) {
            StabilityPropertyIssue issue = StabilityPropertyIssue::root_range;
            switch (error.code()) {
            case th::Sw92PhaseErrorCode::near_multiple:
                issue = StabilityPropertyIssue::root_topology;
                break;
            case th::Sw92PhaseErrorCode::iteration_limit:
                issue = StabilityPropertyIssue::root_iteration_limit;
                break;
            case th::Sw92PhaseErrorCode::unrepresentable_root:
                break;
            }
            throw StabilityPropertyError(issue, error.what());
        } catch (const std::range_error& error) {
            throw StabilityPropertyError(StabilityPropertyIssue::root_range, error.what());
        }
    }

private:
    thermodynamics::Sw92Phase<double> model_;
    thermodynamics::Sw92PhaseWorkspace<double> workspace_;
    double molality_{};
    thermodynamics::SwPhaseFamily family_{thermodynamics::SwPhaseFamily::aqueous};
    thermodynamics::Sw92RootOptions root_options_;
};

struct Sw92FamilyStabilityResult {
    StabilityResult search;
    std::string dataset_id;
    std::string revision;
    std::vector<std::string> component_ids;
    std::string model_profile{thermodynamics::sw92_corrected_profile};
    std::string phase_convention{thermodynamics::sw92_pt_convention};
    double nacl_molality_mol_per_kg_water{};
    thermodynamics::SwPhaseFamily family{thermodynamics::SwPhaseFamily::aqueous};
    thermodynamics::Sw92RootOptions root_options;
};

/// Fixed-family feed stability. This does not select or compare aqueous and
/// nonaqueous reference families. Cross-family orchestration is a separate layer.
[[nodiscard]] inline Sw92FamilyStabilityResult test_sw92_pt_family_stability(
    double pressure_pa, double temperature_k, std::span<const double> feed,
    Sw92FamilyStabilityEvaluator& evaluator, StabilityOptions options = {},
    std::span<const std::vector<double>> extra_starts = {}) {
    if (feed.size() != evaluator.model().size()) {
        throw std::invalid_argument("SW92 stability: feed does not match ordered model snapshot");
    }
    Sw92FamilyStabilityResult result;
    const auto& parameters = evaluator.model().parameters();
    result.dataset_id = parameters.dataset_id();
    result.revision = parameters.revision();
    result.nacl_molality_mol_per_kg_water = evaluator.nacl_molality_mol_per_kg_water();
    result.family = evaluator.family();
    result.root_options = evaluator.root_options();
    for (const auto& component : parameters.components().items()) {
        result.component_ids.push_back(component.id);
    }
    result.search = test_pt_stability(pressure_pa, temperature_k, feed,
                                      evaluator, options, extra_starts);
    return result;
}

} // namespace mpmc::flash

#endif // MPMC_FLASH_SW92_STABILITY_HPP
