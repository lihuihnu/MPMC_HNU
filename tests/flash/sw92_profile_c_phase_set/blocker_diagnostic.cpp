#include <mpmc/flash/sw92_profile_c_phase_set.hpp>

#include "../sw92_phase_assigned_pt/physical_sample6.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <iostream>
#include <vector>

namespace {
namespace fl = mpmc::flash;
namespace sample6 = sw92_profile_c_sample6;
using Vec = std::vector<double>;

Vec edge_feed(const Vec& first, const Vec& second, double second_fraction) {
    Vec feed(first.size(), 0.0);
    for (std::size_t i = 0; i < feed.size(); ++i) {
        feed[i] = (1.0 - second_fraction) * first[i] +
                  second_fraction * second[i];
    }
    return feed;
}

double log_distance(const Vec& a, const Vec& b) {
    double distance = 0.0;
    for (std::size_t i = 0; i < a.size(); ++i) {
        if (a[i] <= 0.0 || b[i] <= 0.0) { continue; }
        distance = std::max(distance, std::abs(std::log(a[i]) - std::log(b[i])));
    }
    return distance;
}

void print_trials(const fl::Sw92PhaseAssignedHSideWitnessResult& result,
                  const Vec& retained_h, const Vec& h0, const Vec& h1) {
    if (!result.nonaqueous_search) { return; }
    std::cout << "effective_tpd=" << result.effective_tpd_tolerance
              << " common_allowance=" << result.common_reference_allowance << '\n';
    for (std::size_t i = 0; i < result.nonaqueous_search->trials.size(); ++i) {
        const auto& trial = result.nonaqueous_search->trials[i];
        std::cout << "trial[" << i << "] status=" << static_cast<int>(trial.status);
        if (trial.point) {
            std::cout << " tpd=" << trial.point->value
                      << " guard=" << trial.point->roundoff_guard
                      << " stationarity=" << trial.point->stationarity
                      << " d_retained=" << log_distance(trial.point->composition, retained_h)
                      << " d_h0=" << log_distance(trial.point->composition, h0)
                      << " d_h1=" << log_distance(trial.point->composition, h1)
                      << " water=" << trial.point->composition.front();
        }
        std::cout << '\n';
    }
}

} // namespace

int main() {
    const auto model = sample6::model();
    const Vec h0 = sample6::h0();
    const Vec h1 = sample6::h1();
    const Vec feed = edge_feed(h0, h1, 0.90);
    const auto source = fl::solve_sw92_phase_assigned_pt_boundary_aware(
        1.0e7, 350.0, feed, model, 0.0);
    const auto published = fl::project_sw92_profile_c_pt_phase_set(source);

    std::cout << "source.status=" << static_cast<int>(source.status)
              << " base.status=" << static_cast<int>(source.base.status)
              << " no_w.status=" << static_cast<int>(source.base.no_w.status)
              << " no_w.count=" << source.base.no_w.retained_na_candidate_count()
              << " source.phases=" << source.phases.size()
              << " boundary=" << static_cast<bool>(source.boundary)
              << " c1=" << static_cast<bool>(source.base.c1)
              << " c2a1=" << static_cast<bool>(source.base.c2a1)
              << " c2b1=" << static_cast<bool>(source.base.c2b1)
              << " c2b2=" << static_cast<bool>(source.base.c2b2)
              << " published.status=" << static_cast<int>(published.solution.status)
              << " published.count=" << published.solution.accepted_phase_count()
              << '\n';

    for (std::size_t i = 0;
         i < source.base.no_w.retained_na_candidate_compositions.size(); ++i) {
        const auto& candidate = source.base.no_w.retained_na_candidate_compositions[i];
        std::cout << "no_w.candidate[" << i << "] fraction="
                  << source.base.no_w.retained_na_candidate_fractions[i]
                  << " d_h0=" << log_distance(candidate, h0)
                  << " d_h1=" << log_distance(candidate, h1)
                  << " water=" << candidate.front() << '\n';
    }

    if (source.base.c1 && source.base.c1->candidate()) {
        const auto& retained_h = source.base.c1->candidate()->nonaqueous_phase.composition;
        std::cout << "c1.H d_h0=" << log_distance(retained_h, h0)
                  << " d_h1=" << log_distance(retained_h, h1)
                  << " water=" << retained_h.front() << '\n';
        fl::Sw92PhaseAssignedHSideWitnessOptions options;
        const std::vector<Vec> physical_starts{h0, h1};
        const auto challenged = fl::test_sw92_phase_assigned_h_side_na_witness(
            *source.base.c1, model, options, physical_starts);
        std::cout << "c2a1.physical-start.status="
                  << static_cast<int>(challenged.status)
                  << " usable=" << challenged.usable_witness_count()
                  << " negatives=" << challenged.negative_witnesses.size() << '\n';
        print_trials(challenged, retained_h, h0, h1);
    }

    if (source.base.c2a1) {
        std::cout << "c2a1.status=" << static_cast<int>(source.base.c2a1->status)
                  << " usable=" << source.base.c2a1->usable_witness_count()
                  << " negatives=" << source.base.c2a1->negative_witnesses.size()
                  << '\n';
    }
    if (source.boundary) {
        std::cout << "boundary.status=" << static_cast<int>(source.boundary->status)
                  << " neighbor_closed=" << source.boundary->neighbor_locally_closed()
                  << '\n';
    }
    std::cout << "source.diagnostic=" << source.diagnostic << '\n'
              << "publication.diagnostic=" << published.solution.diagnostic << '\n';
    return 0;
}
