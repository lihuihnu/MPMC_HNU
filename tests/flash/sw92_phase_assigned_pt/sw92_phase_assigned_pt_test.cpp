#include <mpmc/flash/sw92_phase_assigned_pt.hpp>

#include "../sw92_phase_assigned_three_phase_closure/fixture.hpp"
#include "physical_sample6.hpp"
#include "test_support.hpp"

#include <array>
#include <cmath>
#include <cstddef>
#include <iostream>
#include <source_location>
#include <stdexcept>
#include <string>
#include <string_view>
#include <vector>

bool sw92_phase_assigned_pt_header();

namespace {
namespace fl = mpmc::flash;
namespace th = mpmc::thermodynamics;
namespace sample6 = sw92_profile_c_sample6;
using Vec = std::vector<double>;

void require(bool condition, std::string_view message,
             std::source_location where = std::source_location::current()) {
    if (!condition) {
        throw std::runtime_error(std::string(where.file_name()) + ":" +
                                 std::to_string(where.line()) + ": " + std::string(message));
    }
}

void near(double actual, long double expected, long double relative = 3e-8L,
          long double absolute = 3e-11L,
          std::source_location where = std::source_location::current()) {
    if (!std::isfinite(actual) ||
        std::abs(static_cast<long double>(actual) - expected) >
            absolute + relative * std::abs(expected)) {
        std::cerr << "actual=" << actual << " expected=" << expected << '\n';
        require(false, "reference mismatch", where);
    }
}

Vec log_ratio(const Vec& numerator, const Vec& denominator) {
    require(numerator.size() == denominator.size(), "log-ratio dimension mismatch");
    Vec values(numerator.size());
    for (std::size_t i = 0; i < values.size(); ++i) {
        require(numerator[i] > 0.0 && denominator[i] > 0.0,
                "physical reference requires positive log-ratio support");
        values[i] = std::log(numerator[i]) - std::log(denominator[i]);
    }
    return values;
}

th::Sw92Phase<double> binary_model(bool reverse = false) {
    return th::Sw92Phase<double>::from_parameters(
        sw92_test::binary_parameters(sw92_test::co2, reverse));
}

void wet_binary_routes_to_w_h() {
    const auto result = fl::solve_sw92_phase_assigned_pt(
        3.0e6, 340.0, Vec{0.7, 0.3}, binary_model(), 0.0);
    require(result.status == fl::Sw92PhaseAssignedPtStatus::w_h_locally_closed &&
                result.locally_closed_phase_candidate() &&
                result.phases.size() == 2U && result.c1 && result.c2a1,
            "traceable wet binary did not route to locally closed W+H");
    require(result.phases[0].physical_role ==
                fl::Sw92PhaseAssignedPtPhysicalRole::aqueous &&
                result.phases[1].physical_role ==
                fl::Sw92PhaseAssignedPtPhysicalRole::nonaqueous_unclassified &&
                result.phases[0].thermodynamic_family == th::SwPhaseFamily::aqueous &&
                result.phases[1].thermodynamic_family == th::SwPhaseFamily::nonaqueous,
            "top-level W/H role-family mapping changed");
    near(result.phases[0].composition[0],
         0.0059450651661861014407375682315186536L);
    near(result.phases[1].composition[0],
         0.98878773430905531065855598854496999L);
    near(result.phases[1].mole_phase_fraction,
         0.70617094335057182413614298038232832L);
    require(!result.global_stability_proven &&
                !result.accepted_phase_set_published && !result.morphology_resolved,
            "top-level driver exceeded finite/morphology-neutral semantics");
}

void dry_binary_routes_to_no_w_single() {
    const auto result = fl::solve_sw92_phase_assigned_pt(
        3.0e6, 340.0, Vec{0.995, 0.005}, binary_model(), 0.0);
    require(result.status ==
                fl::Sw92PhaseAssignedPtStatus::no_w_single_h_locally_closed &&
                result.locally_closed_phase_candidate() &&
                result.phases.size() == 1U && !result.c1,
            "dry binary did not close at no-W single NA/H candidate");
    require(result.phases.front().physical_role ==
                fl::Sw92PhaseAssignedPtPhysicalRole::nonaqueous_unclassified &&
                result.phases.front().thermodynamic_family ==
                    th::SwPhaseFamily::nonaqueous,
            "dry top-level phase mapping changed");
    near(result.phases.front().composition[0], 0.995L, 0.0L, 2e-15L);
}

void zero_water_routes_to_no_w_single() {
    const auto result = fl::solve_sw92_phase_assigned_pt(
        3.0e6, 340.0, Vec{1.0, 0.0}, binary_model(), 0.0);
    require(result.status ==
                fl::Sw92PhaseAssignedPtStatus::no_w_single_h_locally_closed &&
                result.phases.size() == 1U &&
                !result.no_w.aqueous_appearance_search,
            "zero-water inventory did not close directly in no-W topology");
}

void synthetic_ternary_routes_to_three_phase() {
    const auto model = c2b2_test::model(false);
    const auto feed = c2b2_test::feed_from_c1_beta(0.73L, false);
    const auto result = fl::solve_sw92_phase_assigned_pt(
        3.0e6, 260.0, feed, model, 0.0);
    if (!(result.status ==
              fl::Sw92PhaseAssignedPtStatus::w_h0_h1_locally_closed &&
          result.locally_closed_phase_candidate() &&
          result.phases.size() == 3U && result.c1 && result.c2a1 &&
          result.c2b1 && result.c2b2)) {
        std::cerr << "top-level status=" << static_cast<int>(result.status)
                  << " noW=" << static_cast<int>(result.no_w.status)
                  << " c1_attempts=" << result.c1_attempts
                  << " c2b1_attempts=" << result.c2b1_attempts
                  << " diagnostic=" << result.diagnostic << '\n';
        require(false,
                "synthetic structural state did not traverse top-level three-phase path");
    }
    require(result.c2b2->w_present_locally_closed() &&
                result.phases[0].physical_role ==
                    fl::Sw92PhaseAssignedPtPhysicalRole::aqueous &&
                result.phases[1].physical_role ==
                    fl::Sw92PhaseAssignedPtPhysicalRole::nonaqueous_unclassified &&
                result.phases[2].physical_role ==
                    fl::Sw92PhaseAssignedPtPhysicalRole::nonaqueous_unclassified,
            "three-phase top-level publication semantics changed");
    for (const auto& phase : result.phases) {
        require(phase.mole_phase_fraction > 0.0 &&
                    phase.mole_phase_fraction < 1.0 &&
                    phase.composition.size() == 3U &&
                    phase.compressibility_factor.has_value(),
                "top-level three-phase instance is incomplete");
    }
    near(result.phases[0].composition[0], c2b2_test::w[0]);
}

void check_physical_three_phase_state(
    const fl::Sw92PhaseAssignedThreePhaseState& state) {
    const auto phases = std::array<const fl::Sw92PhaseAssignedThreePhasePhase*, 3>{
        &state.aqueous_phase, &state.hydrocarbon0_phase, &state.hydrocarbon1_phase};
    const auto expected_composition =
        std::array<const std::array<long double, 8>*, 3>{
            &sample6::golden.w, &sample6::golden.h0, &sample6::golden.h1};
    for (std::size_t phase = 0; phase < phases.size(); ++phase) {
        require(phases[phase]->composition.size() == sample6::golden.feed.size(),
                "physical Sample-6 component dimension changed");
        for (std::size_t i = 0; i < sample6::golden.feed.size(); ++i) {
            near(phases[phase]->composition[i], (*expected_composition[phase])[i],
                 8e-8L, 3e-12L);
        }
        near(phases[phase]->mole_phase_fraction, sample6::golden.fractions[phase],
             8e-8L, 3e-12L);
        near(phases[phase]->compressibility_factor, sample6::golden.z[phase],
             8e-8L, 3e-12L);
    }
    for (std::size_t i = 0; i < sample6::golden.common.size(); ++i) {
        near(state.common_log_activity[i], sample6::golden.common[i],
             8e-8L, 3e-11L);
    }
    require(state.aqueous_phase.thermodynamic_family == th::SwPhaseFamily::aqueous &&
                state.hydrocarbon0_phase.thermodynamic_family ==
                    th::SwPhaseFamily::nonaqueous &&
                state.hydrocarbon1_phase.thermodynamic_family ==
                    th::SwPhaseFamily::nonaqueous,
            "physical Sample-6 family mapping changed");
    require(state.aqueous_phase.composition[0] >
                std::max(state.hydrocarbon0_phase.composition[0],
                         state.hydrocarbon1_phase.composition[0]),
            "physical Sample-6 W phase lost water-rich role ordering");
    require(state.hydrocarbon0_phase.compressibility_factor <
                state.hydrocarbon1_phase.compressibility_factor,
            "physical Sample-6 H slots lost deterministic lower-Z canonical order");
    require(state.chemical_potential_norm <= 1.0e-11 &&
                state.mass_absolute <= 1.0e-12 &&
                state.mass_relative <= 1.0e-10 &&
                state.generalized_rr_residual <= 2.0e-13,
            "physical Sample-6 equilibrium/balance gate changed");
}

void physical_sample6_direct_three_phase() {
    const auto model = sample6::model();
    const Vec feed = sample6::feed();
    const Vec w = sample6::w();
    const Vec h0 = sample6::h0();
    const Vec h1 = sample6::h1();
    const auto result = fl::iterate_sw92_phase_assigned_three_phase_candidate(
        1.0e7, 350.0, feed, log_ratio(h0, w), log_ratio(h1, w),
        {static_cast<double>(sample6::golden.fractions[1]),
         static_cast<double>(sample6::golden.fractions[2])},
        model, 0.0);
    require(result.status == fl::Sw92PhaseAssignedThreePhaseStatus::converged_candidate &&
                result.candidate() != nullptr,
            "physical Sample-6 direct three-phase kernel did not converge");
    check_physical_three_phase_state(*result.candidate());
    require(result.dataset_id ==
                "MR2017-Sample6-SW92-corrected-original-PR76-base" &&
                !result.global_stability_proven &&
                !result.accepted_phase_set_published,
            "physical Sample-6 provenance/candidate semantics changed");
}

void physical_sample6_top_level_three_phase() {
    const auto result = fl::solve_sw92_phase_assigned_pt(
        1.0e7, 350.0, sample6::feed(), sample6::model(), 0.0);
    if (!(result.status == fl::Sw92PhaseAssignedPtStatus::w_h0_h1_locally_closed &&
          result.phases.size() == 3U && result.c1 && result.c2a1 &&
          result.c2b1 && result.c2b2 && result.c2b1->candidate())) {
        std::cerr << "physical Sample-6 top-level status=" <<
            static_cast<int>(result.status)
                  << " noW=" << static_cast<int>(result.no_w.status)
                  << " c1_attempts=" << result.c1_attempts
                  << " c2b1_attempts=" << result.c2b1_attempts
                  << " diagnostic=" << result.diagnostic << '\n';
        if (result.c1) {
            std::cerr << "c1=" << static_cast<int>(result.c1->status)
                      << " diagnostic=" << result.c1->diagnostic << '\n';
        }
        if (result.c2a1) {
            std::cerr << "c2a1=" << static_cast<int>(result.c2a1->status)
                      << " negative=" << result.c2a1->negative_witnesses.size()
                      << " diagnostic=" << result.c2a1->diagnostic << '\n';
        }
        if (result.c2b1) {
            std::cerr << "c2b1=" << static_cast<int>(result.c2b1->status)
                      << " diagnostic=" << result.c2b1->diagnostic << '\n';
        }
        if (result.c2b2) {
            std::cerr << "c2b2=" << static_cast<int>(result.c2b2->status)
                      << " diagnostic=" << result.c2b2->diagnostic << '\n';
        }
        require(false,
                "physical Mortezazadeh-Rasaei Sample-6 did not traverse the top-level three-phase flash path");
    }
    require(result.c2b2->w_present_locally_closed(),
            "physical Sample-6 did not pass the final W-present H-multiplicity review");
    check_physical_three_phase_state(*result.c2b1->candidate());
    for (std::size_t phase = 0; phase < result.phases.size(); ++phase) {
        near(result.phases[phase].mole_phase_fraction,
             sample6::golden.fractions[phase], 8e-8L, 3e-12L);
    }
}

void physical_sample6_direct_permutation() {
    const auto model = sample6::model(true);
    const Vec feed = sample6::feed(true);
    const Vec w = sample6::w(true);
    const Vec h0 = sample6::h0(true);
    const Vec h1 = sample6::h1(true);
    const auto result = fl::iterate_sw92_phase_assigned_three_phase_candidate(
        1.0e7, 350.0, feed, log_ratio(h0, w), log_ratio(h1, w),
        {static_cast<double>(sample6::golden.fractions[1]),
         static_cast<double>(sample6::golden.fractions[2])},
        model, 0.0);
    require(result.candidate() != nullptr,
            "physical Sample-6 reversed component order did not converge");
    const auto& state = *result.candidate();
    for (std::size_t i = 0; i < sample6::golden.feed.size(); ++i) {
        const std::size_t reversed = sample6::golden.feed.size() - 1U - i;
        near(state.aqueous_phase.composition[reversed], sample6::golden.w[i],
             8e-8L, 3e-12L);
        near(state.hydrocarbon0_phase.composition[reversed], sample6::golden.h0[i],
             8e-8L, 3e-12L);
        near(state.hydrocarbon1_phase.composition[reversed], sample6::golden.h1[i],
             8e-8L, 3e-12L);
    }
}

void binary_component_permutation() {
    const auto normal = fl::solve_sw92_phase_assigned_pt(
        3.0e6, 340.0, Vec{0.7, 0.3}, binary_model(false), 0.0);
    const auto reverse = fl::solve_sw92_phase_assigned_pt(
        3.0e6, 340.0, Vec{0.3, 0.7}, binary_model(true), 0.0);
    require(normal.status == reverse.status &&
                normal.status == fl::Sw92PhaseAssignedPtStatus::w_h_locally_closed &&
                normal.phases.size() == reverse.phases.size(),
            "component permutation changed top-level topology status");
    near(normal.phases[0].mole_phase_fraction,
         reverse.phases[0].mole_phase_fraction);
    near(normal.phases[0].composition[0], reverse.phases[0].composition[1]);
    near(normal.phases[1].composition[0], reverse.phases[1].composition[1]);
}

void attempt_quota_is_explicit() {
    fl::Sw92PhaseAssignedPtOptions options;
    options.max_c1_attempts = 1U;
    const auto result = fl::solve_sw92_phase_assigned_pt(
        3.0e6, 340.0, Vec{0.7, 0.3}, binary_model(), 0.0, options);
    require(result.status == fl::Sw92PhaseAssignedPtStatus::w_h_locally_closed ||
                (result.status == fl::Sw92PhaseAssignedPtStatus::topology_unresolved &&
                 result.c1_attempt_limit_reached),
            "bounded C1 attempt quota was hidden or misreported");
}

void headers() {
    require(sw92_phase_assigned_pt_header(),
            "top-level Profile-C PT public header probe failed");
}

} // namespace

int main(int argc, char** argv) {
    if (argc != 2) {
        std::cerr << "expected one test-case name\n";
        return 2;
    }
    const std::string_view name{argv[1]};
    try {
        if (name == "wet_w_h") wet_binary_routes_to_w_h();
        else if (name == "dry_no_w") dry_binary_routes_to_no_w_single();
        else if (name == "zero_water") zero_water_routes_to_no_w_single();
        else if (name == "synthetic_three_phase") synthetic_ternary_routes_to_three_phase();
        else if (name == "physical_three_phase_direct") physical_sample6_direct_three_phase();
        else if (name == "physical_three_phase_top_level") physical_sample6_top_level_three_phase();
        else if (name == "physical_three_phase_permutation") physical_sample6_direct_permutation();
        else if (name == "permutation") binary_component_permutation();
        else if (name == "attempt_quota") attempt_quota_is_explicit();
        else if (name == "headers") headers();
        else throw std::invalid_argument("unknown test-case name");
    } catch (const std::exception& error) {
        std::cerr << error.what() << '\n';
        return 1;
    }
    return 0;
}
