bool phase_set_transition_header();

#include <mpmc/flow/phase_set_transition.hpp>
#include <mpmc/flow/phase_set_transition_flash_adapter.hpp>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <iostream>
#include <optional>
#include <stdexcept>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

namespace flow = mpmc::flow;
namespace flash = mpmc::flash;

namespace {

void require_transition(
    bool condition,
    std::string_view message) {
    if (!condition) {
        throw std::runtime_error(
            std::string{message});
    }
}

void near_transition(
    double actual,
    double expected,
    double tolerance = 1.0e-12) {
    if (!std::isfinite(actual) ||
        !std::isfinite(expected) ||
        std::abs(actual - expected) >
            tolerance *
                std::max(
                    1.0,
                    std::max(
                        std::abs(actual),
                        std::abs(expected)))) {
        throw std::runtime_error(
            "phase-transition numeric mismatch");
    }
}

template <class Function>
void expect_invalid_transition(
    Function&& function) {
    try {
        function();
    } catch (const std::invalid_argument&) {
        return;
    }
    throw std::runtime_error(
        "expected invalid phase transition");
}

flow::PoreVolumeComponentAccumulationSnapshot3P
source_inventory() {
    return {
        0.25,
        {"A", "B", "C"},
        {0.30, 0.40, 0.30},
        1.0};
}

flow::PoreVolumeEnergyAccumulationSnapshot3P
source_energy_history() {
    flow::NaturalVariableStateIdentity identity;
    identity.layout =
        flow::NaturalVariableLayoutDescriptor{
            3U,
            1U,
            {1U}};
    identity.component_ids =
        {"A", "B", "C"};
    identity.reference_pressure_pa =
        10.0;
    identity.temperature_k =
        8.0;
    identity.saturation =
        {1.0, 0.0, 0.0};
    identity.phase_composition[0] =
        {0.30, 0.40, 0.30};
    identity.phase_composition[1].clear();
    identity.phase_composition[2].clear();

    return {
        std::move(identity),
        0.25,
        70.0,
        30.0,
        100.0};
}

flow::PhaseSetTransitionCandidate
appearance_candidate() {
    flow::PhaseSetTransitionCandidate candidate;
    candidate.source_phase_count = 1U;
    candidate.target_phase_count = 2U;
    candidate.trigger =
        flow::PhaseSetTransitionTrigger::
            stability_witness;
    candidate.status =
        flow::PhaseSetTransitionCandidateStatus::
            target_resolved;
    candidate.pressure_pa = 10.0;
    candidate.temperature_k = 8.0;
    candidate.component_ids =
        {"A", "B", "C"};
    candidate.evidence_profile =
        "controlled/stability/v1";
    candidate.diagnostic =
        "fresh two-phase target";

    candidate.target_phases = {
        {
            0.40,
            {0.20, 0.50, 0.30},
            10.0},
        {
            0.60,
            {
                11.0 / 30.0,
                1.0 / 3.0,
                0.30},
            5.0}
    };
    return candidate;
}

flash::PtFlashBackendResult
accepted_flash_two_phase() {
    flash::PtFlashBackendResult result;
    result.capability.backend_id =
        "controlled";
    result.capability.model_profile =
        "controlled-model";
    result.capability.algorithm_profile =
        "controlled-algorithm";
    result.capability.publication_profile =
        "controlled-publication";
    result.capability.configuration_profile =
        "controlled-configuration";
    result.capability.dataset_id =
        "controlled-dataset";
    result.capability.revision =
        "v1";
    result.capability.component_ids =
        {"A", "B", "C"};
    result.capability.supported_phase_counts =
        {1U, 2U, 3U};
    result.capability.transition_capability.edges =
        {
            {
                1U,
                2U,
                flash::
                    PtPhaseTransitionSupport::
                        fresh_target_resolve,
                true},
            {
                2U,
                1U,
                flash::
                    PtPhaseTransitionSupport::
                        fresh_target_resolve,
                true},
            {
                2U,
                3U,
                flash::
                    PtPhaseTransitionSupport::
                        fresh_target_resolve,
                true},
            {
                3U,
                2U,
                flash::
                    PtPhaseTransitionSupport::
                        fresh_target_resolve,
                true}
        };
    result.capability
        .performs_initial_stability_search =
        true;
    result.capability
        .performs_final_phase_set_review =
        true;

    result.solution.capability
        .maximum_phase_count = 3U;
    result.solution.status =
        flash::PtPhaseSetStatus::accepted;
    result.solution.pressure_pa = 10.0;
    result.solution.temperature_k = 8.0;
    result.solution.feed =
        {0.30, 0.40, 0.30};

    flash::PtCandidatePhaseSet set;
    const auto candidate =
        appearance_candidate();
    for (const auto& phase :
         candidate.target_phases) {
        flash::PtCandidatePhase projected;
        projected.mole_phase_fraction =
            phase.mole_phase_fraction;
        projected.composition =
            phase.composition;
        projected.activity.ln_phi =
            {0.0, 0.0, 0.0};
        projected.activity.branch = 0U;
        projected.activity.smooth = true;
        set.phases.push_back(
            std::move(projected));
    }
    result.solution.candidate_phase_set =
        std::move(set);
    result.provider_result_convention =
        "controlled-result/v1";
    result.morphology_resolved = false;

    result.transition_report.evidence.push_back(
        {
            1U,
            std::optional<std::size_t>{2U},
            flash::PtPhaseTransitionTrigger::
                initial_stability_witness,
            flash::PtPhaseTransitionResolution::
                accepted_target,
            true,
            true,
            "controlled/stability/v1",
            "fresh two-phase target"});
    return result;
}

} // namespace

void phase_set_transition_contract() {
    const auto inventory =
        source_inventory();

    {
        const auto candidate =
            appearance_candidate();
        const auto projection =
            flow::
                project_phase_set_transition_candidate(
                    candidate,
                    inventory);

        require_transition(
            projection.source_phase_count() ==
                    1U &&
                projection.target_phase_count() ==
                    2U &&
                projection.layout().phase_count() ==
                    2U &&
                projection.layout().unknown_count() ==
                    7U &&
                projection.natural_variables()
                        .size() ==
                    7U,
            "1->2 projection cardinality mismatch");

        near_transition(
            projection.target_saturations()[0],
            0.25);
        near_transition(
            projection.target_saturations()[1],
            0.75);
        near_transition(
            projection
                .mixture_molar_density_mol_per_m3(),
            6.25);

        const auto q =
            projection.natural_variables();
        near_transition(q[0], 10.0);
        near_transition(q[1], 8.0);
        near_transition(q[2], 0.25);

        const auto overall =
            projection.overall_composition();
        near_transition(overall[0], 0.30);
        near_transition(overall[1], 0.40);
        near_transition(overall[2], 0.30);

        require_transition(
            projection
                    .layout()
                    .composition_pivot()
                    .dependent_components()[0] ==
                1U &&
                projection
                    .layout()
                    .composition_pivot()
                    .dependent_components()[1] ==
                0U,
            "transition pivot selection changed");
    }

    {
        auto unresolved =
            appearance_candidate();
        unresolved.status =
            flow::
                PhaseSetTransitionCandidateStatus::
                    target_resolve_required;
        unresolved.target_phases.clear();

        expect_invalid_transition(
            [&] {
                (void)flow::
                    project_phase_set_transition_candidate(
                        unresolved,
                        inventory);
            });
    }

    {
        auto mismatch =
            appearance_candidate();
        mismatch.target_phases[1]
            .composition =
            {0.40, 0.30, 0.30};

        expect_invalid_transition(
            [&] {
                (void)flow::
                    project_phase_set_transition_candidate(
                        mismatch,
                        inventory);
            });
    }

    {
        flow::PhaseSetTransitionCandidate
            disappearance;
        disappearance.source_phase_count =
            2U;
        disappearance.target_phase_count =
            1U;
        disappearance.trigger =
            flow::PhaseSetTransitionTrigger::
                phase_disappearance;
        disappearance.status =
            flow::
                PhaseSetTransitionCandidateStatus::
                    target_resolved;
        disappearance.pressure_pa = 10.0;
        disappearance.temperature_k = 8.0;
        disappearance.component_ids =
            {"A", "B", "C"};
        disappearance.evidence_profile =
            "controlled/disappearance/v1";
        disappearance.target_phases = {
            {
                1.0,
                {0.30, 0.40, 0.30},
                8.0}
        };

        const auto projection =
            flow::
                project_phase_set_transition_candidate(
                    disappearance,
                    inventory);
        require_transition(
            projection.layout().phase_count() ==
                    1U &&
                projection.layout().unknown_count() ==
                    4U,
            "2->1 projection cardinality mismatch");
        near_transition(
            projection.target_saturations()[0],
            1.0);

        auto invalid_direction =
            disappearance;
        invalid_direction.trigger =
            flow::PhaseSetTransitionTrigger::
                stability_witness;
        expect_invalid_transition(
            [&] {
                (void)flow::
                    project_phase_set_transition_candidate(
                        invalid_direction,
                        inventory);
            });
    }

    {
        const auto energy =
            source_energy_history();
        const auto migrated =
            flow::
                migrate_phase_set_transition_history(
                    1U,
                    2U,
                    std::array<std::string, 3>{
                        "A", "B", "C"},
                    0.25,
                    inventory,
                    energy);

        require_transition(
            migrated.source_phase_count ==
                    1U &&
                migrated.target_phase_count ==
                    2U &&
                migrated
                        .previous_component_accumulation
                        .component_accumulation_mol_per_bulk_m3 ==
                    inventory
                        .component_accumulation_mol_per_bulk_m3 &&
                migrated
                        .previous_component_accumulation
                        .total_accumulation_mol_per_bulk_m3 ==
                    inventory
                        .total_accumulation_mol_per_bulk_m3 &&
                migrated
                        .previous_energy_accumulation
                        .total_internal_energy_j_per_bulk_m3 ==
                    energy
                        .total_internal_energy_j_per_bulk_m3 &&
                migrated
                        .previous_energy_accumulation
                        .state_identity
                        .layout.phase_count() ==
                    1U,
            "transition conservation history was recomputed instead of preserved");
    }

    {
        const std::array<std::size_t, 3>
            phase_counts{1U, 2U, 3U};
        const auto active =
            flow::NaturalVariableActiveSetLayout::
                create(
                    3U,
                    phase_counts);

        require_transition(
            active.cell_count() == 3U &&
                active.cell(0U).scalar_offset ==
                    0U &&
                active.cell(0U).scalar_count ==
                    4U &&
                active.cell(1U).scalar_offset ==
                    4U &&
                active.cell(1U).scalar_count ==
                    7U &&
                active.cell(2U).scalar_offset ==
                    11U &&
                active.cell(2U).scalar_count ==
                    10U &&
                active.total_scalar_count() ==
                    21U,
            "per-cell 1/2/3 active-set scalar layout mismatch");
    }

    {
        auto flash_result =
            accepted_flash_two_phase();
        const std::array<double, 2>
            density{10.0, 5.0};

        const auto adapted =
            flow::
                make_phase_set_transition_candidate_from_flash(
                    1U,
                    flash_result,
                    density);
        require_transition(
            adapted.has_value() &&
                adapted->target_phase_count ==
                    2U &&
                adapted->trigger ==
                    flow::
                        PhaseSetTransitionTrigger::
                            stability_witness &&
                adapted->status ==
                    flow::
                        PhaseSetTransitionCandidateStatus::
                            target_resolved,
            "accepted flash phase-set was not adapted into a fresh flow transition candidate");

        const auto projection =
            flow::
                project_phase_set_transition_candidate(
                    *adapted,
                    inventory);
        near_transition(
            projection.target_saturations()[0],
            0.25);

        const auto same =
            flow::
                make_phase_set_transition_candidate_from_flash(
                    2U,
                    flash_result,
                    density);
        require_transition(
            !same.has_value(),
            "same-cardinality flash result must not trigger a transition");
    }
}

int main() {
    try {
        phase_set_transition_contract();
        if (!phase_set_transition_header()) {
            throw std::runtime_error(
                "phase-set transition public-header probe failed");
        }
        std::cout
            << "[PASS] phase-set transition contract\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr
            << "[FAIL] "
            << error.what()
            << '\n';
        return 1;
    }
}
