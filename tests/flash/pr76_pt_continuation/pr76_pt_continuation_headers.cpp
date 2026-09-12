#include <mpmc/flash/pr76_pt_continuation.hpp>

#include <type_traits>

static_assert(
    mpmc::flash::Pr76PtContinuationResult::convention ==
    mpmc::flash::pr76_pt_continuation_convention);
static_assert(!mpmc::flash::Pr76PtContinuationResult::global_stability_proven);
static_assert(std::is_default_constructible_v<mpmc::flash::Pr76PtPathState>);

int main() {
    mpmc::flash::Pr76PtTransitionBracket bracket;
    return bracket.exact_boundary_resolved ? 1 : 0;
}
