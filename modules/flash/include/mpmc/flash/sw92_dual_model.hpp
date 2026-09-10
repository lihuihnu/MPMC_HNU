#ifndef MPMC_FLASH_SW92_DUAL_MODEL_HPP
#define MPMC_FLASH_SW92_DUAL_MODEL_HPP

#include <mpmc/flash/sw92_split.hpp>

#include <cmath>
#include <cstddef>
#include <limits>
#include <optional>
#include <span>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

namespace mpmc::flash {

/// Evidence-aligned compatibility profile: two complete, independent fixed-family
/// equilibrium calculations. This identifier does not denote a joint AQ/NA phase set.
inline constexpr std::string_view sw92_whitson_dual_model_algorithm =
    "SW92-equilibrium/whitson-dual-model-observables/v1";

struct Sw92DualModelPassOptions {
    thermodynamics::Sw92RootOptions root_options;
    PtSplitOptions split_options;
};

struct Sw92DualModelOptions {
    Sw92DualModelPassOptions aqueous;
    Sw92DualModelPassOptions nonaqueous;
};

/// Non-owning start views used only for the duration of one orchestration call.
struct Sw92DualModelStarts {
    std::span<const std::vector<double>> aqueous_initial{};
    std::span<const std::vector<double>> aqueous_final{};
    std::span<const std::vector<double>> nonaqueous_initial{};
    std::span<const std::vector<double>> nonaqueous_final{};
};

enum class Sw92DualModelStatus {
    observables_available,
    both_family_observables_unavailable,
    aqueous_observable_unavailable,
    nonaqueous_observable_unavailable,
    multicomponent_phase_label_not_implemented,
    binary_phase_label_indeterminate,
    cross_model_ratio_unrepresentable
};

/// One extracted phase from one independently accepted fixed-family run.
/// `source_role` records which numerical low/high-Z candidate supplied the
/// water-rich/poor binary observable; it is not an AQ/NA family identity.
struct Sw92DualModelTargetPhase {
    std::vector<double> composition;
    PtPhaseRole source_role{PtPhaseRole::liquid_candidate};
};

/// Binary-only compatibility observables currently authorized by the algorithm
/// contract. The ratio is explicitly cross-model: y_NA / x_AQ componentwise.
struct Sw92DualModelObservables {
    std::size_t water_index{};
    Sw92DualModelTargetPhase aqueous;
    Sw92DualModelTargetPhase nonaqueous;
    std::vector<double> cross_model_equilibrium_ratio;
};

/// Profile-A result. The two family runs remain complete and independent.
/// No method projects their extracted phases into PtPhaseSetResult or claims one
/// material balance/common tangent across AQ and NA models.
struct Sw92WhitsonDualModelResult {
    Sw92FamilyPtSplitResult aqueous_run;
    Sw92FamilyPtSplitResult nonaqueous_run;
    std::optional<Sw92DualModelObservables> observables;
    Sw92DualModelStatus status{Sw92DualModelStatus::both_family_observables_unavailable};
    std::string diagnostic;

    std::string dataset_id;
    std::string revision;
    std::vector<std::string> component_ids;
    std::string model_profile{thermodynamics::sw92_corrected_profile};
    std::string phase_convention{thermodynamics::sw92_pt_convention};
    std::string equilibrium_algorithm{sw92_whitson_dual_model_algorithm};
    double nacl_molality_mol_per_kg_water{};

    [[nodiscard]] bool compatibility_observables_available() const noexcept {
        return status == Sw92DualModelStatus::observables_available && observables.has_value();
    }
};

namespace detail {

inline bool sw92_family_two_phase_accepted(const Sw92FamilyPtSplitResult& run) noexcept {
    return run.solution.status == PtSplitStatus::two_phase_no_instability_found &&
           run.solution.candidate() != nullptr;
}

inline std::optional<Sw92DualModelTargetPhase> sw92_select_binary_water_target(
    const Sw92FamilyPtSplitResult& run, std::size_t water_index, bool water_richer,
    bool& unresolved_tie) {
    unresolved_tie = false;
    const auto* point = run.solution.candidate();
    if (point == nullptr || water_index >= point->fractions.liquid.size() ||
        water_index >= point->fractions.vapor.size()) {
        return std::nullopt;
    }
    const double liquid_water = point->fractions.liquid[water_index];
    const double vapor_water = point->fractions.vapor[water_index];
    if (!std::isfinite(liquid_water) || !std::isfinite(vapor_water)) {
        return std::nullopt;
    }
    const double scale = 1.0 + std::abs(liquid_water) + std::abs(vapor_water);
    const double guard = 256.0 * std::numeric_limits<double>::epsilon() * scale;
    const double difference = liquid_water - vapor_water;
    if (std::abs(difference) <= guard) {
        unresolved_tie = true;
        return std::nullopt;
    }
    const bool choose_liquid = water_richer ? difference > 0.0 : difference < 0.0;
    Sw92DualModelTargetPhase selected;
    selected.composition = choose_liquid ? point->fractions.liquid : point->fractions.vapor;
    selected.source_role = choose_liquid ? PtPhaseRole::liquid_candidate
                                         : PtPhaseRole::vapor_candidate;
    return selected;
}

} // namespace detail

/// Run the evidence-aligned Whitson dual-model compatibility profile.
///
/// Both passes are constructed from the SAME ordered thermodynamic snapshot and
/// fixed molality. AQ and NA use independent numerical options/resources and each
/// performs its own feed stability, material balance, fugacity iteration and final
/// same-family stability review. The two phase fractions are never combined.
///
/// Automatic physical target-phase extraction is intentionally limited to the
/// audited binary case. Multicomponent calls still execute and retain both complete
/// family runs, but report that physical phase labeling is not yet implemented.
[[nodiscard]] inline Sw92WhitsonDualModelResult solve_sw92_whitson_dual_model_observables(
    double pressure_pa, double temperature_k, std::span<const double> feed,
    const thermodynamics::Sw92Phase<double>& model,
    double nacl_molality_mol_per_kg_water,
    Sw92DualModelOptions options = {}, Sw92DualModelStarts starts = {}) {
    Sw92FamilyVleEvaluator aqueous_evaluator(
        model, nacl_molality_mol_per_kg_water,
        thermodynamics::SwPhaseFamily::aqueous, options.aqueous.root_options);
    Sw92FamilyVleEvaluator nonaqueous_evaluator(
        model, nacl_molality_mol_per_kg_water,
        thermodynamics::SwPhaseFamily::nonaqueous, options.nonaqueous.root_options);

    Sw92WhitsonDualModelResult result;
    const auto& parameters = model.parameters();
    result.dataset_id = parameters.dataset_id();
    result.revision = parameters.revision();
    result.nacl_molality_mol_per_kg_water = nacl_molality_mol_per_kg_water;
    for (const auto& component : parameters.components().items()) {
        result.component_ids.push_back(component.id);
    }

    result.aqueous_run = solve_sw92_pt_family_vle(
        pressure_pa, temperature_k, feed, aqueous_evaluator,
        options.aqueous.split_options, starts.aqueous_initial, starts.aqueous_final);
    result.nonaqueous_run = solve_sw92_pt_family_vle(
        pressure_pa, temperature_k, feed, nonaqueous_evaluator,
        options.nonaqueous.split_options, starts.nonaqueous_initial, starts.nonaqueous_final);

    const bool aqueous_accepted = detail::sw92_family_two_phase_accepted(result.aqueous_run);
    const bool nonaqueous_accepted = detail::sw92_family_two_phase_accepted(result.nonaqueous_run);
    if (!aqueous_accepted || !nonaqueous_accepted) {
        if (!aqueous_accepted && !nonaqueous_accepted) {
            result.status = Sw92DualModelStatus::both_family_observables_unavailable;
            result.diagnostic = "AQ and NA fixed-family runs did not both provide accepted two-phase targets";
        } else if (!aqueous_accepted) {
            result.status = Sw92DualModelStatus::aqueous_observable_unavailable;
            result.diagnostic = "AQ fixed-family run did not provide an accepted two-phase target";
        } else {
            result.status = Sw92DualModelStatus::nonaqueous_observable_unavailable;
            result.diagnostic = "NA fixed-family run did not provide an accepted two-phase target";
        }
        return result;
    }

    if (model.size() != 2) {
        result.status = Sw92DualModelStatus::multicomponent_phase_label_not_implemented;
        result.diagnostic =
            "both family runs retained; automatic physical target-phase labeling is binary-only";
        return result;
    }

    const std::size_t water_index = parameters.water_index();
    bool aqueous_tie = false;
    bool nonaqueous_tie = false;
    auto aqueous_target = detail::sw92_select_binary_water_target(
        result.aqueous_run, water_index, true, aqueous_tie);
    auto nonaqueous_target = detail::sw92_select_binary_water_target(
        result.nonaqueous_run, water_index, false, nonaqueous_tie);
    if (!aqueous_target || !nonaqueous_target) {
        result.status = Sw92DualModelStatus::binary_phase_label_indeterminate;
        result.diagnostic = (aqueous_tie || nonaqueous_tie)
            ? "binary water-rich/water-poor ordering is unresolved at numerical roundoff scale"
            : "accepted family run lacks a representable binary target phase";
        return result;
    }

    Sw92DualModelObservables observables;
    observables.water_index = water_index;
    observables.aqueous = std::move(*aqueous_target);
    observables.nonaqueous = std::move(*nonaqueous_target);
    observables.cross_model_equilibrium_ratio.resize(model.size());
    for (std::size_t i = 0; i < model.size(); ++i) {
        const double denominator = observables.aqueous.composition[i];
        const double numerator = observables.nonaqueous.composition[i];
        if (!std::isfinite(denominator) || !std::isfinite(numerator) ||
            !(denominator > 0.0) || numerator < 0.0) {
            result.status = Sw92DualModelStatus::cross_model_ratio_unrepresentable;
            result.diagnostic = "cross-model equilibrium ratio has a nonrepresentable component";
            return result;
        }
        const double ratio = numerator / denominator;
        if (!std::isfinite(ratio)) {
            result.status = Sw92DualModelStatus::cross_model_ratio_unrepresentable;
            result.diagnostic = "cross-model equilibrium ratio overflowed";
            return result;
        }
        observables.cross_model_equilibrium_ratio[i] = ratio;
    }

    result.observables = std::move(observables);
    result.status = Sw92DualModelStatus::observables_available;
    return result;
}

} // namespace mpmc::flash

#endif // MPMC_FLASH_SW92_DUAL_MODEL_HPP
