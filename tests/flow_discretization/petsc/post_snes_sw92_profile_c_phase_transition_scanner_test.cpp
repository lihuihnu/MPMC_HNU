#include <mpmc/flow_discretization_petsc/post_snes_sw92_profile_c_phase_transition_scanner.hpp>

#include "../../flash/sw92_phase_assigned_pt/physical_sample6.hpp"
#include "test_support.hpp"

#include <petscsys.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <optional>
#include <span>
#include <stdexcept>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

bool post_snes_sw92_profile_c_phase_transition_scanner_header();
void sw92_transactional_phase_transition_restart_test();

namespace {

namespace fdp = mpmc::flow_discretization_petsc;
namespace fl = mpmc::flash;
namespace flow = mpmc::flow;
namespace mesh = mpmc::mesh;
namespace sample6 = sw92_profile_c_sample6;
namespace st = sw92_test;
namespace th = mpmc::thermodynamics;

void require_sw92_scanner(
    bool condition,
    std::string_view message) {
    if (!condition) {
        throw std::runtime_error(
            std::string{message});
    }
}

void near_sw92_scanner(
    double actual,
    double expected,
    double relative = 5.0e-8,
    double absolute = 5.0e-11) {
    const double scale =
        std::max(
            {1.0,
             std::abs(actual),
             std::abs(expected)});
    require_sw92_scanner(
        std::isfinite(actual) &&
            std::isfinite(expected) &&
            std::abs(actual - expected) <=
                absolute +
                    relative * scale,
        "SW92 authoritative scanner numeric mismatch");
}

th::Sw92Phase<double>
binary_model() {
    return th::Sw92Phase<double>::
        from_parameters(
            st::binary_parameters(
                st::co2));
}

void same_cardinality_is_noop() {
    auto model = binary_model();

    fdp::PostSnesPtFlashSourceCellSnapshot3D
        source;
    source.cell_global =
        mesh::GlobalEntityId{UINT64_C(8101)};
    source.source_phase_count = 2U;
    source.pressure_pa = 3.0e6;
    source.temperature_k = 340.0;
    source.component_ids = {
        st::co2.id,
        st::water.id};
    source.overall_composition = {
        0.70,
        0.30};

    std::optional<
        fdp::PostSnesSw92ProfileCTransitionScanCellResult3D>
        transition;
    fdp::PostSnesPhaseTransitionScanStatus3D
        status =
            fdp::PostSnesPhaseTransitionScanStatus3D::
                indeterminate;
    const PetscErrorCode error =
        fdp::
            scan_post_snes_sw92_profile_c_source_cell_3d(
                source,
                model,
                {},
                &transition,
                &status);
    require_sw92_scanner(
        error == PETSC_SUCCESS &&
            status ==
                fdp::
                    PostSnesPhaseTransitionScanStatus3D::
                        complete &&
            !transition.has_value(),
        "authoritative 2P SW92 state generated a false transition");
}

flow::PoreVolumeComponentAccumulationSnapshot3P
source_component_history(
    const fdp::PostSnesPtFlashSourceCellSnapshot3D&
        source) {
    constexpr double porosity = 0.25;
    constexpr double total = 40.0;
    std::vector<double> amounts;
    amounts.reserve(
        source.overall_composition.size());
    for (double fraction :
         source.overall_composition) {
        amounts.push_back(
            total * fraction);
    }
    return {
        porosity,
        source.component_ids,
        std::move(amounts),
        total};
}

flow::PoreVolumeEnergyAccumulationSnapshot3P
source_energy_history(
    const fdp::PostSnesPtFlashSourceCellSnapshot3D&
        source) {
    flow::NaturalVariableStateIdentity3P
        identity;
    identity.layout =
        flow::NaturalVariableLayoutDescriptor{
            source.component_ids.size(),
            2U,
            std::vector<std::size_t>{
                0U,
                0U}};
    identity.component_ids =
        source.component_ids;
    identity.reference_pressure_pa =
        source.pressure_pa;
    identity.temperature_k =
        source.temperature_k;
    identity.saturation =
        {0.5, 0.5, 0.0};
    identity.phase_composition[0] =
        source.overall_composition;
    identity.phase_composition[1] =
        source.overall_composition;

    return {
        std::move(identity),
        0.25,
        2.0e6,
        1.0e6,
        3.0e6};
}

void sample6_two_to_three_authoritative_target() {
    const auto model =
        sample6::model();

    fdp::PostSnesPtFlashSourceCellSnapshot3D
        source;
    source.cell_global =
        mesh::GlobalEntityId{UINT64_C(8201)};
    source.source_phase_count = 2U;
    source.pressure_pa = 1.0e7;
    source.temperature_k = 350.0;
    source.component_ids =
        sample6::normal_order();
    source.overall_composition =
        sample6::feed();

    std::optional<
        fdp::PostSnesSw92ProfileCTransitionScanCellResult3D>
        transition;
    fdp::PostSnesPhaseTransitionScanStatus3D
        status =
            fdp::PostSnesPhaseTransitionScanStatus3D::
                indeterminate;
    const PetscErrorCode error =
        fdp::
            scan_post_snes_sw92_profile_c_source_cell_3d(
                source,
                model,
                {},
                &transition,
                &status);
    require_sw92_scanner(
        error == PETSC_SUCCESS &&
            status ==
                fdp::
                    PostSnesPhaseTransitionScanStatus3D::
                        complete &&
            transition.has_value(),
        "Sample-6 SW92 2->3 authoritative scanner did not publish a transition");

    const auto& proposal =
        transition->proposal;
    const auto& target =
        transition->target;
    require_sw92_scanner(
        proposal.cell_global ==
                source.cell_global &&
            proposal.candidate
                    .source_phase_count ==
                2U &&
            proposal.candidate
                    .target_phase_count ==
                3U &&
            proposal.candidate.status ==
                flow::
                    PhaseSetTransitionCandidateStatus::
                        target_resolved &&
            proposal.candidate
                    .component_ids ==
                source.component_ids &&
            !proposal.candidate
                 .evidence_profile.empty() &&
            target.cell_global ==
                source.cell_global &&
            target.component_ids ==
                source.component_ids &&
            target.phases.size() == 3U &&
            target.projection
                    .source_phase_count() ==
                2U &&
            target.projection
                    .target_phase_count() ==
                3U &&
            target.projection
                    .layout()
                    .phase_count() ==
                3U &&
            target.projection
                    .natural_variables()
                    .size() ==
                3U *
                    source.component_ids.size() +
                    1U &&
            target.transition_evidence_profile ==
                proposal.candidate
                    .evidence_profile,
        "Sample-6 authoritative target materialization shape/provenance mismatch");

    require_sw92_scanner(
        target.phases[0]
                    .metadata
                    .thermodynamic_family ==
                th::SwPhaseFamily::aqueous &&
            target.phases[1]
                    .metadata
                    .thermodynamic_family ==
                th::SwPhaseFamily::nonaqueous &&
            target.phases[2]
                    .metadata
                    .thermodynamic_family ==
                th::SwPhaseFamily::nonaqueous,
        "Sample-6 authoritative target lost AQ/NA family metadata");

    double saturation_sum = 0.0;
    for (std::size_t phase = 0U;
         phase < 3U;
         ++phase) {
        const auto& materialized =
            target.phases[phase];
        const auto& candidate =
            proposal.candidate
                .target_phases[phase];
        require_sw92_scanner(
            materialized.selection
                    .nacl_molality_mol_per_kg_water ==
                0.0 &&
            materialized.selection
                    .family ==
                materialized.metadata
                    .thermodynamic_family &&
            candidate
                    .provider_activity_branch ==
                std::optional<std::size_t>{
                    materialized.selection
                        .root_index} &&
            candidate
                    .provider_activity_smooth &&
            materialized
                    .molar_density_mol_per_m3 >
                0.0,
            "Sample-6 target family/root/density provenance changed");
        near_sw92_scanner(
            candidate
                .molar_density_mol_per_m3,
            materialized
                .molar_density_mol_per_m3,
            2.0e-11,
            2.0e-12);
        saturation_sum +=
            target.projection
                .target_saturations()[phase];
    }
    near_sw92_scanner(
        saturation_sum,
        1.0,
        0.0,
        2.0e-12);

    const auto projected_overall =
        target.projection
            .overall_composition();
    require_sw92_scanner(
        projected_overall.size() ==
            source.overall_composition.size(),
        "Sample-6 target overall composition size changed");
    for (std::size_t component = 0U;
         component <
             projected_overall.size();
         ++component) {
        near_sw92_scanner(
            projected_overall[component],
            source.overall_composition[
                component],
            2.0e-8,
            2.0e-11);
    }

    const auto component_history =
        source_component_history(source);
    flow::
        validate_phase_set_transition_material_balance(
            proposal.candidate,
            component_history,
            {2.0e-8});

    const auto energy_history =
        source_energy_history(source);
    const auto migrated =
        flow::
            migrate_phase_set_transition_history(
                2U,
                3U,
                target.component_ids,
                0.25,
                component_history,
                energy_history);
    require_sw92_scanner(
        migrated.source_phase_count ==
                2U &&
            migrated.target_phase_count ==
                3U &&
            migrated
                    .previous_component_accumulation
                    .component_ids ==
                component_history.component_ids &&
            migrated
                    .previous_component_accumulation
                    .component_accumulation_mol_per_bulk_m3 ==
                component_history
                    .component_accumulation_mol_per_bulk_m3 &&
            migrated
                    .previous_energy_accumulation
                    .total_internal_energy_j_per_bulk_m3 ==
                energy_history
                    .total_internal_energy_j_per_bulk_m3,
        "SW92 transition proposal changed frozen component/energy history");

    fl::Sw92ProfileCPtFlashBackend
        role_neutral_backend{model};
    require_sw92_scanner(
        role_neutral_backend
            .capability()
            .phase_metadata_namespace.empty(),
        "SW92 generic PT backend unexpectedly acquired provider role/family metadata");
}

} // namespace

void post_snes_sw92_profile_c_phase_transition_scanner_test() {
    require_sw92_scanner(
        post_snes_sw92_profile_c_phase_transition_scanner_header(),
        "SW92 authoritative transition scanner public header probe failed");
    same_cardinality_is_noop();
    sample6_two_to_three_authoritative_target();
    sw92_transactional_phase_transition_restart_test();
}
