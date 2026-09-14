#ifndef MPMC_MODEL_CONFIGURATION_PT_SOLVE_HINTS_HPP
#define MPMC_MODEL_CONFIGURATION_PT_SOLVE_HINTS_HPP

#include <array>
#include <string>
#include <string_view>
#include <vector>

namespace mpmc::model_configuration {

// Versioned, transport-neutral numerical initialization hints for one PT solve.
// They are never phase-count evidence and never alter the immutable model/settings
// snapshots. Every composition is in the model snapshot's ordered component list.
inline constexpr std::string_view pt_solve_hints_v1 = "pt-solve-hints/v1";

// Three-phase continuation seed. phase_fraction_seed stores the independent
// fractions of phases 1 and 2; phase 0 is 1-beta1-beta2, matching the existing
// generalized-RR native contract. The current solve still decides the phase set.
struct PtThreePhaseContinuationHint {
    std::array<std::vector<double>, 3> compositions;
    std::array<double, 2> phase_fraction_seed{1.0 / 3.0, 1.0 / 3.0};

    bool operator==(const PtThreePhaseContinuationHint&) const = default;
};

struct PtSolveHints {
    // Required exact version when this DTO is supplied. Empty does not mean v1.
    std::string version;

    // Extra multistart compositions for the feed-stability search. Native
    // automatic starts, when enabled by settings, remain independent additions.
    std::vector<std::vector<double>> initial_stability_starts;

    // Extra compositions for the final two-phase common-tangent review. The
    // solver still adds the two freshly computed phase compositions itself.
    std::vector<std::vector<double>> final_two_phase_stability_starts;

    // Optional numerical seeds considered only after the established two-phase
    // final review has found additional-phase evidence. They cannot force 3 phases.
    std::vector<PtThreePhaseContinuationHint> three_phase_continuation_starts;

    bool operator==(const PtSolveHints&) const = default;
};

[[nodiscard]] inline PtSolveHints make_pt_solve_hints_v1() {
    PtSolveHints hints;
    hints.version = std::string(pt_solve_hints_v1);
    return hints;
}

} // namespace mpmc::model_configuration

#endif // MPMC_MODEL_CONFIGURATION_PT_SOLVE_HINTS_HPP
