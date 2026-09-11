#include <mpmc/flash/sw92_phase_assigned_boundary.hpp>

#include "../sw92_phase_assigned_pt/physical_sample6.hpp"
#include "../sw92_phase_assigned_three_phase_closure/fixture.hpp"

#include <algorithm>
#include <cstddef>
#include <iostream>
#include <source_location>
#include <stdexcept>
#include <string>
#include <string_view>
#include <vector>

bool sw92_phase_assigned_boundary_header();

namespace {
namespace fl = mpmc::flash;
namespace sample6 = sw92_profile_c_sample6;
using Vec = std::vector<double>;

void require(bool condition, std::string_view message,
             std::source_location where = std::source_location::current()) {
    if (!condition) {
        throw std::runtime_error(std::string(where.file_name()) + ":" +
                                 std::to_string(where.line()) + ": " + std::string(message));
    }
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

void physical_w_h_edge_re_solve() {
    const auto model = sample6::model();
    const Vec w = sample6::w();
    const Vec h0 = sample6::h0();
    const Vec h1 = sample6::h1();
    const Vec feed = edge_feed(w, h0, 0.45);
    const std::vector<Vec> diagnostics{h1};

    const auto result = fl::detail::sw92_phase_assigned_resolve_w_h_neighbor(
        1.0e7, 350.0, feed, w, h0, model, 0.0,
        fl::Sw92PhaseAssignedPtOptions{}, diagnostics);
    require(result.status == fl::Sw92PhaseAssignedBoundaryStatus::resolved_to_w_h &&
                result.neighbor_locally_closed() && result.re_solve_attempted &&
                result.neighbor_c1 && result.neighbor_c1->candidate() &&
                result.neighbor_c2a1 && result.phases.size() == 2U,
            "physical W+H edge did not re-solve as a locally closed two-phase neighbor");
    require(result.neighbor_c2a1->status ==
                fl::Sw92PhaseAssignedHSideWitnessStatus::
                    no_additional_nonaqueous_witness_found,
            "zero-amount third H remained a robust instability on the physical W+H edge");
}

void physical_no_w_edge_re_solve() {
    const auto model = sample6::model();
    const Vec h0 = sample6::h0();
    const Vec h1 = sample6::h1();
    const Vec feed = edge_feed(h0, h1, 0.90);
    const std::vector<Vec> starts{h0, h1};

    const auto result = fl::detail::sw92_phase_assigned_resolve_no_w_neighbor(
        1.0e7, 350.0, feed, starts, model, 0.0,
        fl::Sw92PhaseAssignedPtOptions{});
    if (result.status !=
        fl::Sw92PhaseAssignedBoundaryStatus::resolved_to_no_w_two_h) {
        std::cerr << "physical no-W edge status=" << static_cast<int>(result.status)
                  << " diagnostic=" << result.diagnostic << '\n';
        if (result.neighbor_no_w) {
            std::cerr << "neighbor no-W status="
                      << static_cast<int>(result.neighbor_no_w->status)
                      << " diagnostic=" << result.neighbor_no_w->diagnostic << '\n';
        }
    }
    require(result.status ==
                fl::Sw92PhaseAssignedBoundaryStatus::resolved_to_no_w_two_h &&
                result.neighbor_locally_closed() && result.re_solve_attempted &&
                result.neighbor_no_w && result.phases.size() == 2U,
            "physical H0+H1 edge did not re-solve as a locally closed no-W neighbor");
}

void synthetic_h_drop_guard() {
    const auto chain = c2b2_test::solve_chain(0.60L, false, 0.25);
    const auto review = c2b2_test::review_chain(chain);
    require(review.status == fl::Sw92PhaseAssignedC2b2Status::route_to_w_h,
            "synthetic H-disappearance route changed");

    const auto result = fl::resolve_sw92_phase_assigned_c2b2_boundary(
        chain.c1, chain.c2a1, chain.c2b1, review,
        c2b2_test::model(false));
    require(result.re_solve_attempted && !result.neighbor_locally_closed() &&
                result.status ==
                    fl::Sw92PhaseAssignedBoundaryStatus::neighbor_topology_not_closed,
            "artificially large H-disappearance threshold silently deleted a still-required phase");
}

void synthetic_w_drop_guard() {
    const auto chain = c2b2_test::solve_chain(0.80L, false, 0.25);
    const auto review = c2b2_test::review_chain(chain);
    require(review.status ==
                fl::Sw92PhaseAssignedC2b2Status::route_to_no_w_h0_h1,
            "synthetic W-disappearance route changed");

    const auto result = fl::resolve_sw92_phase_assigned_c2b2_boundary(
        chain.c1, chain.c2a1, chain.c2b1, review,
        c2b2_test::model(false));
    require(result.re_solve_attempted && !result.neighbor_locally_closed(),
            "artificially large W-disappearance threshold silently deleted a still-required W phase");
    require(result.status ==
                fl::Sw92PhaseAssignedBoundaryStatus::neighbor_topology_not_closed ||
            result.status ==
                fl::Sw92PhaseAssignedBoundaryStatus::higher_phase_count_or_wrong_candidate,
            "W-disappearance neighbor failure was not retained explicitly");
}

void boundary_aware_interior_unchanged() {
    const auto result = fl::solve_sw92_phase_assigned_pt_boundary_aware(
        1.0e7, 350.0, sample6::feed(), sample6::model(), 0.0);
    require(result.status == fl::Sw92PhaseAssignedPtStatus::w_h0_h1_locally_closed &&
                result.locally_closed_phase_candidate() && result.phases.size() == 3U &&
                !result.boundary && !result.boundary_c2b1 && !result.boundary_c2b2,
            "boundary-aware wrapper changed the validated physical interior three-phase solution");
}

void boundary_aware_guard_path() {
    fl::Sw92PhaseAssignedPtOptions options;
    options.c2b1.minimum_phase_fraction = 0.28;
    const auto feed = c2b2_test::feed_from_c1_beta(0.73L, false);
    const auto result = fl::solve_sw92_phase_assigned_pt_boundary_aware(
        3.0e6, 260.0, feed, c2b2_test::model(false), 0.0, options);
    require(result.base.status == fl::Sw92PhaseAssignedPtStatus::topology_unresolved &&
                result.disappearance_attempts > 0U &&
                result.boundary_c2b1 && result.boundary_c2b2 && result.boundary &&
                result.boundary->re_solve_attempted &&
                !result.boundary->neighbor_locally_closed() &&
                result.status == fl::Sw92PhaseAssignedPtStatus::topology_unresolved &&
                result.phases.empty(),
            "boundary-aware top-level path did not re-solve and explicitly reject an artificial phase drop");
}

void headers() {
    require(sw92_phase_assigned_boundary_header(),
            "boundary public header self-containment probe failed");
}

} // namespace

int main(int argc, char** argv) {
    if (argc != 2) {
        std::cerr << "expected one test-case name\n";
        return 2;
    }
    const std::string_view name{argv[1]};
    try {
        if (name == "physical_w_h_edge") physical_w_h_edge_re_solve();
        else if (name == "physical_no_w_edge") physical_no_w_edge_re_solve();
        else if (name == "synthetic_h_drop_guard") synthetic_h_drop_guard();
        else if (name == "synthetic_w_drop_guard") synthetic_w_drop_guard();
        else if (name == "boundary_aware_interior") boundary_aware_interior_unchanged();
        else if (name == "boundary_aware_guard") boundary_aware_guard_path();
        else if (name == "headers") headers();
        else throw std::invalid_argument("unknown test-case name");
    } catch (const std::exception& error) {
        std::cerr << error.what() << '\n';
        return 1;
    }
    return 0;
}
