#include <mpmc/flash/sw92_profile_c_phase_set.hpp>

#include "../sw92_phase_assigned_pt/physical_sample6.hpp"
#include <fugacity_adapter_regression.hpp>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <iostream>
#include <source_location>
#include <stdexcept>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

bool sw92_profile_c_phase_set_header();

namespace {
namespace fl = mpmc::flash;
namespace fo = mpmc::flow;
namespace fa = mpmc::test::flow_adapter;
namespace th = mpmc::thermodynamics;
namespace sample6 = sw92_profile_c_sample6;
using Vec = std::vector<double>;

void require(bool condition, std::string_view message,
             std::source_location where = std::source_location::current()) {
    if (!condition) {
        throw std::runtime_error(std::string(where.file_name()) + ":" +
                                 std::to_string(where.line()) + ": " +
                                 std::string(message));
    }
}

void near(double actual, double expected, double relative = 4e-8,
          double absolute = 4e-11,
          std::source_location where = std::source_location::current()) {
    if (!std::isfinite(actual) || !std::isfinite(expected) ||
        std::abs(actual - expected) >
            absolute + relative * std::abs(expected)) {
        std::cerr << "actual=" << actual << " expected=" << expected << '\n';
        require(false, "numeric mismatch", where);
    }
}

th::Sw92Phase<double> binary_model(bool reverse = false) {
    return th::Sw92Phase<double>::from_parameters(
        sw92_test::binary_parameters(sw92_test::co2, reverse));
}

Vec edge_feed(const Vec& first, const Vec& second, double second_fraction) {
    require(first.size() == second.size() && second_fraction > 0.0 &&
                second_fraction < 1.0,
            "invalid edge-feed construction");
    Vec feed(first.size(), 0.0);
    for (std::size_t i = 0; i < feed.size(); ++i) {
        feed[i] = (1.0 - second_fraction) * first[i] +
                  second_fraction * second[i];
    }
    return feed;
}

void require_publication_basics(
    const fl::Sw92ProfileCPtPhaseSetResult& result,
    std::size_t phase_count) {
    require(result.solution.status == fl::PtPhaseSetStatus::accepted &&
                result.accepted_phase_set_published() &&
                result.solution.accepted_phase_count() == phase_count &&
                result.solution.capability.maximum_phase_count == 3U,
            "authoritative publication status/count/capability mismatch");
    require(!result.solution.global_stability_proven &&
                !result.global_stability_proven && !result.morphology_resolved,
            "publication invented global proof or morphology");
    require(result.solution.accepted_phase_set() != nullptr &&
                result.phase_metadata.size() == phase_count,
            "publication payload/metadata size mismatch");
    require(result.publication_convention ==
                fl::sw92_profile_c_phase_set_publication_convention &&
                result.orchestration_convention ==
                    fl::sw92_phase_assigned_pt_convention &&
                result.boundary_convention ==
                    fl::sw92_phase_assigned_boundary_convention,
            "publication provenance conventions changed");
}

void publish_one_phase() {
    const auto result = fl::solve_sw92_profile_c_pt_phase_set(
        3.0e6, 340.0, Vec{0.995, 0.005}, binary_model(), 0.0);
    require_publication_basics(result, 1U);
    const auto& phase = result.solution.accepted_phase_set()->phases.front();
    require(result.phase_metadata.front().physical_role ==
                fl::Sw92PhaseAssignedPtPhysicalRole::nonaqueous_unclassified &&
                result.phase_metadata.front().thermodynamic_family ==
                    th::SwPhaseFamily::nonaqueous,
            "single-H publication invented hydrocarbon morphology/family");
    require(phase.mole_phase_fraction == 1.0 &&
                phase.composition == Vec({0.995, 0.005}) &&
                phase.activity.ln_phi.size() == 2U && phase.activity.smooth,
            "single-H publication changed the accepted fixed-family state");
    require(!phase.compressibility_factor,
            "single-H publication recomputed/fabricated a Z value");
}

void publish_two_phase() {
    const auto result = fl::solve_sw92_profile_c_pt_phase_set(
        3.0e6, 340.0, Vec{0.7, 0.3}, binary_model(), 0.0);
    require_publication_basics(result, 2U);
    const auto& phases = result.solution.accepted_phase_set()->phases;
    require(result.phase_metadata[0].physical_role ==
                fl::Sw92PhaseAssignedPtPhysicalRole::aqueous &&
                result.phase_metadata[0].thermodynamic_family ==
                    th::SwPhaseFamily::aqueous &&
                result.phase_metadata[1].physical_role ==
                    fl::Sw92PhaseAssignedPtPhysicalRole::nonaqueous_unclassified &&
                result.phase_metadata[1].thermodynamic_family ==
                    th::SwPhaseFamily::nonaqueous,
            "W+H publication changed role/family semantics");
    require(phases[0].compressibility_factor &&
                phases[1].compressibility_factor &&
                phases[0].activity.smooth && phases[1].activity.smooth,
            "W+H publication lost owned phase properties");
    near(phases[0].composition[0], 0.0059450651661861014);
    near(phases[1].composition[0], 0.98878773430905531);
    near(phases[1].mole_phase_fraction, 0.70617094335057182);
}

void publish_sample6_three_phase() {
    const auto result = fl::solve_sw92_profile_c_pt_phase_set(
        1.0e7, 350.0, sample6::feed(), sample6::model(), 0.0);
    require_publication_basics(result, 3U);
    require(result.dataset_id ==
                "MR2017-Sample6-SW92-corrected-original-PR76-base" &&
                result.component_ids == sample6::normal_order(),
            "Sample-6 publication lost ordered model provenance");

    const auto& phases = result.solution.accepted_phase_set()->phases;
    require(result.phase_metadata[0].physical_role ==
                fl::Sw92PhaseAssignedPtPhysicalRole::aqueous &&
                result.phase_metadata[0].thermodynamic_family ==
                    th::SwPhaseFamily::aqueous,
            "Sample-6 W metadata changed");
    for (std::size_t phase = 1U; phase < 3U; ++phase) {
        require(result.phase_metadata[phase].physical_role ==
                    fl::Sw92PhaseAssignedPtPhysicalRole::nonaqueous_unclassified &&
                    result.phase_metadata[phase].thermodynamic_family ==
                        th::SwPhaseFamily::nonaqueous,
                "Sample-6 H publication invented L/V morphology");
    }

    const auto expected = std::array<const std::array<long double, 8>*, 3>{
        &sample6::golden.w, &sample6::golden.h0, &sample6::golden.h1};
    for (std::size_t phase = 0; phase < 3U; ++phase) {
        require(phases[phase].compressibility_factor &&
                    phases[phase].activity.smooth &&
                    phases[phase].activity.ln_phi.size() == 8U,
                "Sample-6 publication lost phase activity/Z");
        near(phases[phase].mole_phase_fraction,
             static_cast<double>(sample6::golden.fractions[phase]), 8e-8, 3e-12);
        near(*phases[phase].compressibility_factor,
             static_cast<double>(sample6::golden.z[phase]), 8e-8, 3e-12);
        for (std::size_t i = 0; i < 8U; ++i) {
            near(phases[phase].composition[i],
                 static_cast<double>((*expected[phase])[i]), 8e-8, 3e-12);
            const double chemical_potential =
                std::log(phases[phase].composition[i]) +
                phases[phase].activity.ln_phi[i];
            near(chemical_potential,
                 static_cast<double>(sample6::golden.common[i]), 8e-8, 4e-10);
        }
    }
}


void flow_fugacity_adapter_on_sample6_three_phase() {
    // Use the already-regressed component-permuted representation of the same
    // accepted Sample-6 state. The current natural-variable baseline makes the
    // final component dependent; normal Sample-6 ordering puts a ~7.45e-12
    // aqueous trace component in that subtractive 1-sum coordinate, which is
    // numerically ill-conditioned even though the physical state is unchanged.
    const auto model = sample6::model(true);
    const auto result = fl::solve_sw92_profile_c_pt_phase_set(
        1.0e7, 350.0, sample6::feed(true), model, 0.0);
    require_publication_basics(result, 3U);

    const auto& phases =
        result.solution.accepted_phase_set()->phases;
    std::array<
        fo::Sw92SelectedPhaseFugacityEvaluator3P<double>::Selection,
        3>
        selections{};
    std::array<double, 3> fractions{};
    std::array<Vec, 3> compositions;

    for (std::size_t phase = 0U; phase < 3U; ++phase) {
        selections[phase].nacl_molality_mol_per_kg_water =
            result.nacl_molality_mol_per_kg_water;
        selections[phase].family =
            result.phase_metadata[phase].thermodynamic_family;
        selections[phase].root_index =
            phases[phase].activity.branch;
        fractions[phase] =
            phases[phase].mole_phase_fraction;
        compositions[phase] =
            phases[phase].composition;
    }

    fo::Sw92SelectedPhaseFugacityEvaluator3P<double>
        adapter{model, selections};
    fa::verify_accepted_state_and_jacobian(
        1.0e7,
        350.0,
        fractions,
        compositions,
        adapter,
        2.0e-7,
        {2.0e-9, 1.0e-4, 5.0e-4, 1.0e-2});
}

void fresh_disappearance_neighbors() {
    const auto model = sample6::model();
    const Vec w = sample6::w();
    const Vec h0 = sample6::h0();
    const Vec h1 = sample6::h1();

    const Vec w_h_feed = edge_feed(w, h0, 0.45);
    const std::vector<Vec> h_diagnostic{h1};
    const auto w_h_neighbor = fl::detail::sw92_phase_assigned_resolve_w_h_neighbor(
        1.0e7, 350.0, w_h_feed, w, h0, model, 0.0,
        fl::Sw92PhaseAssignedPtOptions{}, h_diagnostic);
    require(w_h_neighbor.status ==
                fl::Sw92PhaseAssignedBoundaryStatus::resolved_to_w_h &&
                w_h_neighbor.neighbor_locally_closed() &&
                w_h_neighbor.re_solve_attempted && w_h_neighbor.phases.size() == 2U,
            "fresh physical W+H disappearance neighbor no longer closes");
    const auto published_w_h = fl::solve_sw92_profile_c_pt_phase_set(
        1.0e7, 350.0, w_h_feed, model, 0.0);
    require_publication_basics(published_w_h, 2U);

    // On the H0+H1 edge, a direct fresh no-W re-solve must *not* authorize
    // removal of W: it finds a robust AQ/W appearance witness. The full
    // boundary-aware topology graph is nevertheless allowed to continue from
    // that witness into the W-containing C1 topology. If C1 closes and C2a1
    // finds no additional-H evidence, authoritative publication is W+H -- not
    // the rejected all-NA H0+H1 neighbor. This distinction is the point of the
    // regression: preserve the incipient-W evidence without falsely requiring
    // every rival-topology rejection to make the complete flash indeterminate.
    const Vec no_w_feed = edge_feed(h0, h1, 0.90);
    const std::vector<Vec> starts{h0, h1};
    const auto no_w_neighbor = fl::detail::sw92_phase_assigned_resolve_no_w_neighbor(
        1.0e7, 350.0, no_w_feed, starts, model, 0.0,
        fl::Sw92PhaseAssignedPtOptions{});
    require(no_w_neighbor.status ==
                fl::Sw92PhaseAssignedBoundaryStatus::neighbor_topology_not_closed &&
                no_w_neighbor.re_solve_attempted &&
                !no_w_neighbor.neighbor_locally_closed() &&
                no_w_neighbor.neighbor_no_w &&
                no_w_neighbor.neighbor_no_w->status ==
                    fl::Sw92PhaseAssignedNoWStatus::aqueous_phase_witness_found &&
                no_w_neighbor.neighbor_no_w->selected_water_witness() != nullptr,
            "fresh H0+H1 neighbor lost its incipient-W blocking evidence");

    const auto routed = fl::solve_sw92_profile_c_pt_phase_set(
        1.0e7, 350.0, no_w_feed, model, 0.0);
    require_publication_basics(routed, 2U);
    require(routed.phase_metadata[0].physical_role ==
                fl::Sw92PhaseAssignedPtPhysicalRole::aqueous &&
                routed.phase_metadata[0].thermodynamic_family ==
                    th::SwPhaseFamily::aqueous &&
                routed.phase_metadata[1].physical_role ==
                    fl::Sw92PhaseAssignedPtPhysicalRole::nonaqueous_unclassified &&
                routed.phase_metadata[1].thermodynamic_family ==
                    th::SwPhaseFamily::nonaqueous,
            "incipient-W evidence was bypassed by an authoritative all-NA publication");
}

void provenance_guard() {
    auto source = fl::solve_sw92_phase_assigned_pt_boundary_aware(
        1.0e7, 350.0, sample6::feed(), sample6::model(), 0.0);
    require(source.base.c2b1 && source.base.c2b2,
            "Sample-6 provenance fixture lost C2 chain");

    auto c2_mismatch = source;
    c2_mismatch.base.c2b1->revision += "-mismatch";
    const auto rejected_c2 = fl::project_sw92_profile_c_pt_phase_set(c2_mismatch);
    require(rejected_c2.solution.status == fl::PtPhaseSetStatus::indeterminate &&
                !rejected_c2.accepted_phase_set_published() &&
                !rejected_c2.solution.candidate_phase_set,
            "publication accepted mismatched C2b.1 provenance");

    auto no_w_mismatch = std::move(source);
    no_w_mismatch.base.no_w.hydrocarbon_flash.revision += "-mismatch";
    const auto rejected_no_w = fl::project_sw92_profile_c_pt_phase_set(no_w_mismatch);
    require(rejected_no_w.solution.status == fl::PtPhaseSetStatus::indeterminate &&
                !rejected_no_w.accepted_phase_set_published() &&
                !rejected_no_w.solution.candidate_phase_set,
            "publication accepted mismatched fixed-family provenance");
}

void component_permutation() {
    const auto normal = fl::solve_sw92_profile_c_pt_phase_set(
        1.0e7, 350.0, sample6::feed(false), sample6::model(false), 0.0);
    const auto reverse = fl::solve_sw92_profile_c_pt_phase_set(
        1.0e7, 350.0, sample6::feed(true), sample6::model(true), 0.0);
    require_publication_basics(normal, 3U);
    require_publication_basics(reverse, 3U);
    require(normal.component_ids == sample6::normal_order() &&
                reverse.component_ids == sample6::reversed_order(),
            "published ordered-component provenance did not follow permutation");

    const auto& a = normal.solution.accepted_phase_set()->phases;
    const auto& b = reverse.solution.accepted_phase_set()->phases;
    for (std::size_t phase = 0; phase < 3U; ++phase) {
        near(a[phase].mole_phase_fraction, b[phase].mole_phase_fraction);
        require(a[phase].compressibility_factor && b[phase].compressibility_factor,
                "permuted Sample-6 publication lost Z");
        near(*a[phase].compressibility_factor, *b[phase].compressibility_factor);
        require(normal.phase_metadata[phase].physical_role ==
                    reverse.phase_metadata[phase].physical_role &&
                    normal.phase_metadata[phase].thermodynamic_family ==
                        reverse.phase_metadata[phase].thermodynamic_family,
                "component permutation changed phase role/family metadata");
        for (std::size_t i = 0; i < a[phase].composition.size(); ++i) {
            const std::size_t ri = a[phase].composition.size() - 1U - i;
            near(a[phase].composition[i], b[phase].composition[ri]);
            near(a[phase].activity.ln_phi[i], b[phase].activity.ln_phi[ri]);
        }
    }
}

void h_slot_symmetry() {
    auto source = fl::solve_sw92_phase_assigned_pt_boundary_aware(
        1.0e7, 350.0, sample6::feed(), sample6::model(), 0.0);
    const auto original = fl::project_sw92_profile_c_pt_phase_set(source);
    require_publication_basics(original, 3U);
    require(source.base.c2b1 && source.base.c2b1->point && source.base.c2b2,
            "Sample-6 H-slot fixture lost final three-phase state");

    auto& state = *source.base.c2b1->point;
    std::swap(state.hydrocarbon0_phase, state.hydrocarbon1_phase);
    std::swap(state.log_k_h0, state.log_k_h1);
    std::swap(state.chemical_potential_residual_h0,
              state.chemical_potential_residual_h1);
    std::swap(state.raw_hydrocarbon0_sum, state.raw_hydrocarbon1_sum);
    std::swap(state.aqueous_h0_log_distance, state.aqueous_h1_log_distance);
    std::swap(state.aqueous_minus_h0_water_fraction,
              state.aqueous_minus_h1_water_fraction);
    std::swap(source.phases[1], source.phases[2]);
    std::swap(source.base.c2b2->retained_h0_trivial_point,
              source.base.c2b2->retained_h1_trivial_point);

    const auto swapped = fl::project_sw92_profile_c_pt_phase_set(source);
    require_publication_basics(swapped, 3U);
    const auto& a = original.solution.accepted_phase_set()->phases;
    const auto& b = swapped.solution.accepted_phase_set()->phases;
    require(a[0].composition == b[0].composition &&
                a[1].composition == b[2].composition &&
                a[2].composition == b[1].composition,
            "H-slot exchange changed the physical published phase set");
    for (std::size_t phase = 1U; phase < 3U; ++phase) {
        require(swapped.phase_metadata[phase].physical_role ==
                    fl::Sw92PhaseAssignedPtPhysicalRole::nonaqueous_unclassified &&
                    swapped.phase_metadata[phase].thermodynamic_family ==
                        th::SwPhaseFamily::nonaqueous,
                "H-slot exchange introduced asymmetric morphology metadata");
    }
}

void headers() {
    require(sw92_profile_c_phase_set_header(),
            "Profile-C phase-set public header self-containment probe failed");
}

} // namespace

int main(int argc, char** argv) {
    if (argc != 2) {
        std::cerr << "expected one test-case name\n";
        return 2;
    }
    const std::string_view name{argv[1]};
    try {
        if (name == "publish_one_phase") publish_one_phase();
        else if (name == "publish_two_phase") publish_two_phase();
        else if (name == "publish_sample6_three_phase") publish_sample6_three_phase();
        else if (name == "flow_fugacity_adapter") flow_fugacity_adapter_on_sample6_three_phase();
        else if (name == "fresh_disappearance_neighbors") fresh_disappearance_neighbors();
        else if (name == "provenance_guard") provenance_guard();
        else if (name == "component_permutation") component_permutation();
        else if (name == "h_slot_symmetry") h_slot_symmetry();
        else if (name == "headers") headers();
        else throw std::invalid_argument("unknown test-case name");
    } catch (const std::exception& error) {
        std::cerr << error.what() << '\n';
        return 1;
    }
    return 0;
}
