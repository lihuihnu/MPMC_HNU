#ifndef MPMC_FLOW_PHASE_SET_TRANSITION_HPP
#define MPMC_FLOW_PHASE_SET_TRANSITION_HPP

#include <mpmc/flow/component_accumulation.hpp>
#include <mpmc/flow/energy_accumulation.hpp>
#include <mpmc/flow/natural_variable_cell_state.hpp>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <optional>
#include <span>
#include <stdexcept>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

namespace mpmc::flow {

inline constexpr std::string_view
    phase_set_transition_convention =
        "flow/natural-variable/phase-set-transition/v1";

enum class PhaseSetTransitionTrigger {
    stability_witness,
    final_phase_set_instability,
    phase_disappearance,
    provider_topology_witness,
    provider_boundary_route
};

enum class PhaseSetTransitionCandidateStatus {
    target_resolve_required,
    target_resolved
};

struct PhaseSetTransitionPhaseCandidate {
    PhaseSetTransitionPhaseCandidate() = default;

    PhaseSetTransitionPhaseCandidate(
        double mole_phase_fraction_value,
        std::vector<double> composition_value,
        double molar_density_value,
        std::optional<std::size_t>
            provider_branch = std::nullopt,
        bool provider_smooth = false)
        : mole_phase_fraction(
              mole_phase_fraction_value),
          composition(
              std::move(composition_value)),
          molar_density_mol_per_m3(
              molar_density_value),
          provider_activity_branch(
              provider_branch),
          provider_activity_smooth(
              provider_smooth) {}

    double mole_phase_fraction{};
    std::vector<double> composition;
    double molar_density_mol_per_m3{};

    /// Opaque provider activity branch from the accepted PT phase. This is
    /// thermodynamic branch provenance only; it is never a physical phase
    /// identity and must not be used to infer oil/gas/water or morphology.
    std::optional<std::size_t>
        provider_activity_branch;
    bool provider_activity_smooth{false};
};

struct PhaseSetTransitionCandidate {
    std::size_t source_phase_count{};
    std::size_t target_phase_count{};
    PhaseSetTransitionTrigger trigger{
        PhaseSetTransitionTrigger::
            provider_topology_witness};
    PhaseSetTransitionCandidateStatus status{
        PhaseSetTransitionCandidateStatus::
            target_resolve_required};

    double pressure_pa{};
    double temperature_k{};
    std::vector<std::string> component_ids;
    std::vector<
        PhaseSetTransitionPhaseCandidate>
        target_phases;

    std::string evidence_profile;
    std::string diagnostic;
};

struct PhaseSetTransitionProjectionOptions {
    double material_balance_tolerance{
        1.0e-10};
};

class PhaseSetTransitionProjection {
public:
    PhaseSetTransitionProjection(
        std::size_t source_phase_count,
        std::size_t target_phase_count,
        NaturalVariableLayoutDescriptor layout,
        std::vector<double> natural_variables,
        std::vector<double> target_saturations,
        std::vector<double> overall_composition,
        double mixture_molar_density_mol_per_m3)
        : source_phase_count_(
              source_phase_count),
          target_phase_count_(
              target_phase_count),
          layout_(std::move(layout)),
          natural_variables_(
              std::move(natural_variables)),
          target_saturations_(
              std::move(target_saturations)),
          overall_composition_(
              std::move(overall_composition)),
          mixture_molar_density_mol_per_m3_(
              mixture_molar_density_mol_per_m3) {
        validate();
    }

    [[nodiscard]] std::size_t
    source_phase_count() const noexcept {
        return source_phase_count_;
    }

    [[nodiscard]] std::size_t
    target_phase_count() const noexcept {
        return target_phase_count_;
    }

    [[nodiscard]] const NaturalVariableLayoutDescriptor&
    layout() const noexcept {
        return layout_;
    }

    [[nodiscard]] std::span<const double>
    natural_variables() const noexcept {
        return natural_variables_;
    }

    [[nodiscard]] std::span<const double>
    target_saturations() const noexcept {
        return target_saturations_;
    }

    [[nodiscard]] std::span<const double>
    overall_composition() const noexcept {
        return overall_composition_;
    }

    [[nodiscard]] double
    mixture_molar_density_mol_per_m3()
        const noexcept {
        return mixture_molar_density_mol_per_m3_;
    }

    [[nodiscard]] static constexpr bool
    requires_target_natural_variable_resolve()
        noexcept {
        return true;
    }

private:
    void validate() const {
        if (source_phase_count_ == 0U ||
            source_phase_count_ >
                fixed_three_phase_count ||
            target_phase_count_ == 0U ||
            target_phase_count_ >
                fixed_three_phase_count ||
            source_phase_count_ ==
                target_phase_count_ ||
            layout_.phase_count() !=
                target_phase_count_ ||
            natural_variables_.size() !=
                layout_.unknown_count() ||
            target_saturations_.size() !=
                target_phase_count_ ||
            overall_composition_.size() !=
                layout_.component_count() ||
            !std::isfinite(
                mixture_molar_density_mol_per_m3_) ||
            !(mixture_molar_density_mol_per_m3_ >
              0.0)) {
            throw std::invalid_argument(
                "mpmc::flow: malformed phase-set transition projection");
        }

        double saturation_sum = 0.0;
        for (double value :
             target_saturations_) {
            const bool valid =
                target_phase_count_ == 1U
                    ? std::isfinite(value) &&
                          value == 1.0
                    : std::isfinite(value) &&
                          value > 0.0 &&
                          value < 1.0;
            if (!valid) {
                throw std::invalid_argument(
                    "mpmc::flow: projected saturation is incompatible with target phase cardinality");
            }
            saturation_sum += value;
        }
        const double saturation_tolerance =
            4096.0 *
            std::numeric_limits<double>::
                epsilon() *
            static_cast<double>(
                target_phase_count_);
        if (!std::isfinite(
                saturation_sum) ||
            std::abs(
                saturation_sum - 1.0) >
                saturation_tolerance) {
            throw std::invalid_argument(
                "mpmc::flow: projected phase saturations are not normalized");
        }

        double composition_sum = 0.0;
        for (double value :
             overall_composition_) {
            if (!std::isfinite(value) ||
                !(value > 0.0)) {
                throw std::invalid_argument(
                    "mpmc::flow: projected overall composition must remain on strict positive support");
            }
            composition_sum += value;
        }
        const double composition_tolerance =
            4096.0 *
            std::numeric_limits<double>::
                epsilon() *
            static_cast<double>(
                overall_composition_.size());
        if (!std::isfinite(
                composition_sum) ||
            std::abs(
                composition_sum - 1.0) >
                composition_tolerance) {
            throw std::invalid_argument(
                "mpmc::flow: projected overall composition is not normalized");
        }
    }

    std::size_t source_phase_count_{};
    std::size_t target_phase_count_{};
    NaturalVariableLayoutDescriptor layout_;
    std::vector<double> natural_variables_;
    std::vector<double> target_saturations_;
    std::vector<double> overall_composition_;
    double mixture_molar_density_mol_per_m3_{};
};

struct PhaseSetTransitionConservationHistory {
    std::size_t source_phase_count{};
    std::size_t target_phase_count{};
    PoreVolumeComponentAccumulationSnapshot3P
        previous_component_accumulation;
    PoreVolumeEnergyAccumulationSnapshot3P
        previous_energy_accumulation;
};

struct NaturalVariableActiveSetCell {
    std::size_t phase_count{};
    std::size_t scalar_offset{};
    std::size_t scalar_count{};
};

class NaturalVariableActiveSetLayout {
public:
    [[nodiscard]] static
    NaturalVariableActiveSetLayout
    create(
        std::size_t component_count,
        std::span<const std::size_t>
            phase_counts) {
        if (component_count < 2U ||
            phase_counts.empty()) {
            throw std::invalid_argument(
                "mpmc::flow: active-set layout requires Nc>=2 and at least one cell");
        }

        std::vector<
            NaturalVariableActiveSetCell>
            cells;
        cells.reserve(
            phase_counts.size());

        std::size_t total = 0U;
        for (const std::size_t phase_count :
             phase_counts) {
            if (phase_count == 0U ||
                phase_count >
                    fixed_three_phase_count ||
                component_count >
                    (std::numeric_limits<
                         std::size_t>::max() -
                     1U) /
                        phase_count) {
                throw std::invalid_argument(
                    "mpmc::flow: active-set cell phase count is invalid");
            }
            const std::size_t q =
                phase_count *
                    component_count +
                1U;
            if (q >
                std::numeric_limits<
                    std::size_t>::max() -
                    total) {
                throw std::length_error(
                    "mpmc::flow: active-set scalar count overflow");
            }
            cells.push_back(
                {
                    phase_count,
                    total,
                    q});
            total += q;
        }

        return {
            component_count,
            std::move(cells),
            total};
    }

    [[nodiscard]] std::size_t
    component_count() const noexcept {
        return component_count_;
    }

    [[nodiscard]] std::size_t
    cell_count() const noexcept {
        return cells_.size();
    }

    [[nodiscard]] const
    NaturalVariableActiveSetCell&
    cell(std::size_t index) const {
        return cells_.at(index);
    }

    [[nodiscard]] std::size_t
    total_scalar_count() const noexcept {
        return total_scalar_count_;
    }

private:
    NaturalVariableActiveSetLayout(
        std::size_t component_count,
        std::vector<
            NaturalVariableActiveSetCell>
            cells,
        std::size_t total_scalar_count)
        : component_count_(
              component_count),
          cells_(std::move(cells)),
          total_scalar_count_(
              total_scalar_count) {}

    std::size_t component_count_{};
    std::vector<
        NaturalVariableActiveSetCell>
        cells_;
    std::size_t total_scalar_count_{};
};

namespace phase_set_transition_detail {

[[nodiscard]] inline double
sum_values(
    std::span<const double> values) {
    double sum = 0.0;
    double correction = 0.0;
    for (double value : values) {
        const double increment =
            value - correction;
        const double next =
            sum + increment;
        correction =
            (next - sum) -
            increment;
        sum = next;
    }
    return sum;
}

inline void
validate_component_ids(
    std::span<const std::string>
        component_ids) {
    if (component_ids.size() < 2U) {
        throw std::invalid_argument(
            "mpmc::flow: phase-set transition requires at least two components");
    }
    for (std::size_t i = 0U;
         i < component_ids.size();
         ++i) {
        if (component_ids[i].empty()) {
            throw std::invalid_argument(
                "mpmc::flow: phase-set transition component id is empty");
        }
        for (std::size_t j = i + 1U;
             j < component_ids.size();
             ++j) {
            if (component_ids[i] ==
                component_ids[j]) {
                throw std::invalid_argument(
                    "mpmc::flow: phase-set transition component ids must be unique and ordered");
            }
        }
    }
}

[[nodiscard]] inline std::size_t
select_dependent_component(
    std::span<const double> composition) {
    if (composition.size() < 2U) {
        throw std::invalid_argument(
            "mpmc::flow: phase-set transition composition is too short");
    }
    std::size_t best = 0U;
    double best_value = -1.0;
    double sum = 0.0;
    for (std::size_t component = 0U;
         component < composition.size();
         ++component) {
        const double value =
            composition[component];
        if (!std::isfinite(value) ||
            !(value > 0.0)) {
            throw std::invalid_argument(
                "mpmc::flow: transition projection requires strict positive-support phase compositions");
        }
        sum += value;
        if (value > best_value) {
            best = component;
            best_value = value;
        }
    }

    const double tolerance =
        4096.0 *
        std::numeric_limits<double>::
            epsilon() *
        static_cast<double>(
            composition.size());
    if (!std::isfinite(sum) ||
        std::abs(sum - 1.0) >
            tolerance) {
        throw std::invalid_argument(
            "mpmc::flow: transition phase composition is not normalized");
    }
    return best;
}

inline void
validate_candidate_direction(
    const PhaseSetTransitionCandidate&
        candidate) {
    if (candidate.source_phase_count == 0U ||
        candidate.source_phase_count >
            fixed_three_phase_count ||
        candidate.target_phase_count == 0U ||
        candidate.target_phase_count >
            fixed_three_phase_count ||
        candidate.source_phase_count ==
            candidate.target_phase_count) {
        throw std::invalid_argument(
            "mpmc::flow: phase-set transition must change between 1/2/3 active phases");
    }

    if (candidate.target_phase_count >
            candidate.source_phase_count &&
        candidate.trigger ==
            PhaseSetTransitionTrigger::
                phase_disappearance) {
        throw std::invalid_argument(
            "mpmc::flow: phase appearance cannot be triggered by disappearance evidence");
    }

    if (candidate.target_phase_count <
            candidate.source_phase_count &&
        candidate.trigger !=
            PhaseSetTransitionTrigger::
                phase_disappearance &&
        candidate.trigger !=
            PhaseSetTransitionTrigger::
                provider_boundary_route) {
        throw std::invalid_argument(
            "mpmc::flow: phase removal requires disappearance/boundary evidence");
    }
}

[[nodiscard]] inline std::vector<double>
normalized_inventory_composition(
    const PoreVolumeComponentAccumulationSnapshot3P&
        component_history) {
    if (component_history.component_ids.size() < 2U ||
        component_history
                .component_accumulation_mol_per_bulk_m3
                .size() !=
            component_history.component_ids.size()) {
        throw std::invalid_argument(
            "mpmc::flow: malformed component inventory for phase transition");
    }

    double total = 0.0;
    for (double value :
         component_history
             .component_accumulation_mol_per_bulk_m3) {
        if (!std::isfinite(value) ||
            !(value > 0.0)) {
            throw std::invalid_argument(
                "mpmc::flow: transition inventory must remain on strict positive support");
        }
        total += value;
    }
    if (!std::isfinite(total) ||
        !(total > 0.0)) {
        throw std::invalid_argument(
            "mpmc::flow: transition inventory total is invalid");
    }

    std::vector<double> result;
    result.reserve(
        component_history.component_ids.size());
    for (double value :
         component_history
             .component_accumulation_mol_per_bulk_m3) {
        result.push_back(
            value / total);
    }
    return result;
}

} // namespace phase_set_transition_detail

[[nodiscard]] inline
std::vector<double>
phase_set_transition_overall_composition(
    const PhaseSetTransitionCandidate&
        candidate) {
    using namespace
        phase_set_transition_detail;

    validate_candidate_direction(
        candidate);
    validate_component_ids(
        candidate.component_ids);

    if (candidate.status !=
            PhaseSetTransitionCandidateStatus::
                target_resolved ||
        candidate.target_phases.size() !=
            candidate.target_phase_count) {
        throw std::invalid_argument(
            "mpmc::flow: target topology must be freshly resolved before transition projection");
    }

    const std::size_t n =
        candidate.component_ids.size();
    std::vector<double> overall(
        n,
        0.0);
    double beta_sum = 0.0;

    for (const auto& phase :
         candidate.target_phases) {
        if (!std::isfinite(
                phase.mole_phase_fraction) ||
            !(phase.mole_phase_fraction >
              0.0) ||
            phase.composition.size() != n ||
            !std::isfinite(
                phase.molar_density_mol_per_m3) ||
            !(phase.molar_density_mol_per_m3 >
              0.0)) {
            throw std::invalid_argument(
                "mpmc::flow: resolved transition phase is malformed");
        }

        (void)select_dependent_component(
            phase.composition);

        beta_sum +=
            phase.mole_phase_fraction;
        for (std::size_t component = 0U;
             component < n;
             ++component) {
            overall[component] +=
                phase.mole_phase_fraction *
                phase.composition[component];
        }
    }

    const double tolerance =
        8192.0 *
        std::numeric_limits<double>::
            epsilon() *
        static_cast<double>(
            candidate.target_phase_count);
    if (!std::isfinite(beta_sum) ||
        std::abs(beta_sum - 1.0) >
            tolerance) {
        throw std::invalid_argument(
            "mpmc::flow: target mole phase fractions are not normalized");
    }

    const double overall_sum =
        sum_values(overall);
    if (!std::isfinite(overall_sum) ||
        std::abs(overall_sum - 1.0) >
            tolerance) {
        throw std::invalid_argument(
            "mpmc::flow: transition target phase set violates material balance");
    }

    return overall;
}

inline void
validate_phase_set_transition_material_balance(
    const PhaseSetTransitionCandidate&
        candidate,
    const PoreVolumeComponentAccumulationSnapshot3P&
        source_component_inventory,
    PhaseSetTransitionProjectionOptions
        options = {}) {
    if (!std::isfinite(
            options.material_balance_tolerance) ||
        options.material_balance_tolerance <
            0.0) {
        throw std::invalid_argument(
            "mpmc::flow: invalid phase-transition material-balance tolerance");
    }
    if (source_component_inventory.component_ids !=
        candidate.component_ids) {
        throw std::invalid_argument(
            "mpmc::flow: transition candidate component order differs from source inventory");
    }

    const auto source =
        phase_set_transition_detail::
            normalized_inventory_composition(
                source_component_inventory);
    const auto target =
        phase_set_transition_overall_composition(
            candidate);

    double maximum = 0.0;
    for (std::size_t component = 0U;
         component < source.size();
         ++component) {
        maximum = std::max(
            maximum,
            std::abs(
                source[component] -
                target[component]));
    }

    if (maximum >
        options.material_balance_tolerance) {
        throw std::invalid_argument(
            "mpmc::flow: transition target phase set does not preserve source overall composition");
    }
}

[[nodiscard]] inline
PhaseSetTransitionProjection
project_phase_set_transition_candidate(
    const PhaseSetTransitionCandidate&
        candidate,
    const PoreVolumeComponentAccumulationSnapshot3P&
        source_component_inventory,
    PhaseSetTransitionProjectionOptions
        options = {}) {
    using namespace
        phase_set_transition_detail;

    validate_phase_set_transition_material_balance(
        candidate,
        source_component_inventory,
        options);

    if (!std::isfinite(
            candidate.pressure_pa) ||
        !(candidate.pressure_pa > 0.0) ||
        !std::isfinite(
            candidate.temperature_k) ||
        !(candidate.temperature_k > 0.0) ||
        candidate.evidence_profile.empty()) {
        throw std::invalid_argument(
            "mpmc::flow: transition candidate lacks valid p/T/evidence provenance");
    }

    const std::size_t p =
        candidate.target_phase_count;
    const std::size_t n =
        candidate.component_ids.size();

    std::vector<std::size_t>
        dependent_components;
    dependent_components.reserve(p);

    std::vector<double>
        volume_weights;
    volume_weights.reserve(p);
    double total_volume_weight = 0.0;

    for (const auto& phase :
         candidate.target_phases) {
        dependent_components.push_back(
            select_dependent_component(
                phase.composition));
        const double weight =
            phase.mole_phase_fraction /
            phase.molar_density_mol_per_m3;
        if (!std::isfinite(weight) ||
            !(weight > 0.0)) {
            throw std::range_error(
                "mpmc::flow: transition phase volume weight is invalid");
        }
        volume_weights.push_back(
            weight);
        total_volume_weight += weight;
    }
    if (!std::isfinite(
            total_volume_weight) ||
        !(total_volume_weight > 0.0)) {
        throw std::range_error(
            "mpmc::flow: transition total phase volume is invalid");
    }

    NaturalVariableLayoutDescriptor layout{
        n,
        p,
        dependent_components};
    std::vector<double> q(
        layout.unknown_count(),
        0.0);
    q[layout.pressure_unknown_index()] =
        candidate.pressure_pa;
    q[layout.temperature_unknown_index()] =
        candidate.temperature_k;

    std::vector<double> saturations(
        p,
        0.0);
    for (std::size_t phase = 0U;
         phase < p;
         ++phase) {
        saturations[phase] =
            volume_weights[phase] /
            total_volume_weight;
        const auto slot =
            static_cast<PhaseSlot3>(
                phase);
        if (const auto saturation_column =
                layout
                    .independent_saturation_unknown_index(
                        slot)) {
            q[*saturation_column] =
                saturations[phase];
        }

        for (std::size_t component = 0U;
             component < n;
             ++component) {
            const auto column =
                layout
                    .independent_composition_unknown_index(
                        slot,
                        component);
            if (column) {
                q[*column] =
                    candidate.target_phases[
                        phase]
                        .composition[
                            component];
            }
        }
    }

    const auto overall =
        phase_set_transition_overall_composition(
            candidate);
    const double mixture_density =
        1.0 /
        total_volume_weight;

    return {
        candidate.source_phase_count,
        candidate.target_phase_count,
        std::move(layout),
        std::move(q),
        std::move(saturations),
        overall,
        mixture_density};
}

[[nodiscard]] inline
PhaseSetTransitionConservationHistory
migrate_phase_set_transition_history(
    std::size_t source_phase_count,
    std::size_t target_phase_count,
    std::span<const std::string>
        target_component_ids,
    double target_porosity,
    const PoreVolumeComponentAccumulationSnapshot3P&
        previous_component_accumulation,
    const PoreVolumeEnergyAccumulationSnapshot3P&
        previous_energy_accumulation) {
    using namespace
        phase_set_transition_detail;

    if (source_phase_count == 0U ||
        source_phase_count >
            fixed_three_phase_count ||
        target_phase_count == 0U ||
        target_phase_count >
            fixed_three_phase_count ||
        source_phase_count ==
            target_phase_count) {
        throw std::invalid_argument(
            "mpmc::flow: conservation history migration requires a real 1/2/3 phase-count change");
    }

    validate_component_ids(
        target_component_ids);
    component_accumulation_detail::
        validate_snapshot_identity(
            previous_component_accumulation);
    energy_accumulation_detail::
        validate_snapshot(
            previous_energy_accumulation);

    if (previous_component_accumulation
            .component_ids.size() !=
            target_component_ids.size() ||
        !std::equal(
            previous_component_accumulation
                .component_ids.begin(),
            previous_component_accumulation
                .component_ids.end(),
            target_component_ids.begin()) ||
        previous_energy_accumulation
                .state_identity
                .component_ids.size() !=
            target_component_ids.size() ||
        !std::equal(
            previous_energy_accumulation
                .state_identity
                .component_ids.begin(),
            previous_energy_accumulation
                .state_identity
                .component_ids.end(),
            target_component_ids.begin()) ||
        previous_energy_accumulation
                .state_identity
                .layout.phase_count() !=
            source_phase_count ||
        !std::isfinite(target_porosity) ||
        !(target_porosity > 0.0) ||
        !(target_porosity < 1.0)) {
        throw std::invalid_argument(
            "mpmc::flow: phase transition history identity mismatch");
    }

    const double porosity_tolerance =
        8192.0 *
        std::numeric_limits<double>::
            epsilon() *
        std::max(
            {1.0,
             std::abs(target_porosity),
             std::abs(
                 previous_component_accumulation
                     .porosity),
             std::abs(
                 previous_energy_accumulation
                     .porosity)});

    if (std::abs(
            previous_component_accumulation
                    .porosity -
            target_porosity) >
            porosity_tolerance ||
        std::abs(
            previous_energy_accumulation
                    .porosity -
            target_porosity) >
            porosity_tolerance) {
        throw std::invalid_argument(
            "mpmc::flow: phase transition cannot change frozen porosity history");
    }

    return {
        source_phase_count,
        target_phase_count,
        previous_component_accumulation,
        previous_energy_accumulation};
}

} // namespace mpmc::flow

#endif // MPMC_FLOW_PHASE_SET_TRANSITION_HPP
