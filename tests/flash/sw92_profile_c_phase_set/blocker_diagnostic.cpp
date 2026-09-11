#include <mpmc/flash/sw92_profile_c_phase_set.hpp>

#include "../sw92_phase_assigned_pt/physical_sample6.hpp"

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

} // namespace

int main() {
    const auto model = sample6::model();
    const Vec feed = edge_feed(sample6::h0(), sample6::h1(), 0.90);
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
    if (source.base.c2a1) {
        std::cout << "c2a1.status=" << static_cast<int>(source.base.c2a1->status)
                  << " usable=" << source.base.c2a1->usable_witness_count() << '\n';
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
