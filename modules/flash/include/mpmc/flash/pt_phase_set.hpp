#ifndef MPMC_FLASH_PT_PHASE_SET_HPP
#define MPMC_FLASH_PT_PHASE_SET_HPP

#include <mpmc/flash/pt_stability.hpp>

#include <cstddef>
#include <optional>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

namespace mpmc::flash {

// Capability of one flash algorithm/entry point, not a thermodynamic proof that
// the underlying model can never admit more phases.
struct PtPhaseSetCapability {
    std::size_t maximum_phase_count{1};
};

// One candidate phase instance. Vector position and provider branch are
// diagnostics only; neither is a universal liquid/vapor/aqueous phase identity.
struct PtCandidatePhase {
    double mole_phase_fraction{};
    std::vector<double> composition;
    StabilityPhase activity;
    std::optional<double> compressibility_factor;
};

struct PtCandidatePhaseSet {
    std::vector<PtCandidatePhase> phases;
};

// Acceptance is deliberately independent of phase count. Callers obtain the
// actual accepted phase count from accepted_phase_set()->phases.size().
enum class PtPhaseSetStatus {
    accepted,
    phase_set_unstable,
    indeterminate
};

struct PtPhaseSetResult {
    static constexpr std::string_view convention = "PT/phase-set-v1";

    PtPhaseSetCapability capability;
    PtPhaseSetStatus status{PtPhaseSetStatus::indeterminate};
    double pressure_pa{};
    double temperature_k{};
    std::vector<double> feed;
    std::optional<PtCandidatePhaseSet> candidate_phase_set;
    bool global_stability_proven{false};
    std::string diagnostic;

    // Structural publication guard only. Thermodynamic validity remains the
    // producer's responsibility; this prevents an obviously malformed empty or
    // over-capability set from being exposed through the accepted helper.
    [[nodiscard]] const PtCandidatePhaseSet* accepted_phase_set() const & noexcept {
        if (status != PtPhaseSetStatus::accepted || !candidate_phase_set) {
            return nullptr;
        }
        const std::size_t count = candidate_phase_set->phases.size();
        if (count == 0 || count > capability.maximum_phase_count) {
            return nullptr;
        }
        return &*candidate_phase_set;
    }
    const PtCandidatePhaseSet* accepted_phase_set() const && = delete;

    [[nodiscard]] std::size_t accepted_phase_count() const noexcept {
        const auto* accepted = accepted_phase_set();
        return accepted == nullptr ? 0 : accepted->phases.size();
    }
};

namespace detail {

// Pure copies of already-reviewed max3 payloads. Model adapters retain their
// public types, provenance, two-phase projector and diagnostic wording. These
// helpers never select roots, solve equilibrium or repeat scientific acceptance.
template <class Phase>
[[nodiscard]] PtCandidatePhase project_max3_owned_phase(const Phase& source) {
    PtCandidatePhase phase;
    phase.mole_phase_fraction = source.mole_phase_fraction;
    phase.composition = source.composition;
    phase.activity = source.activity;
    phase.compressibility_factor = source.z;
    return phase;
}

template <class Source, class ProjectTwoPhase>
[[nodiscard]] PtPhaseSetResult project_max3_solution(
    const Source& source, ProjectTwoPhase project_two_phase,
    std::string_view missing_candidate_diagnostic) {
    using Status = decltype(source.status);
    PtPhaseSetResult result;
    const auto copy_feed = [&] {
        result.capability.maximum_phase_count = 3U;
        result.pressure_pa = source.base.solution.initial_stability.pressure_pa;
        result.temperature_k = source.base.solution.initial_stability.temperature_k;
        result.feed = source.base.solution.initial_stability.feed;
        result.global_stability_proven = false;
    };
    const auto project_two = [&](const auto& two_phase) {
        auto projected = project_two_phase(two_phase);
        projected.capability.maximum_phase_count = 3U;
        return projected;
    };
    switch (source.status) {
    case Status::single_phase:
        result = project_two(source.base);
        break;
    case Status::two_phase:
        if (const auto* neighbor = source.two_phase_neighbor()) {
            result = project_two(*neighbor);
        } else {
            result = project_two(source.base);
        }
        break;
    case Status::three_phase: {
        copy_feed();
        const auto* candidate = source.three_phase_candidate();
        if (candidate == nullptr) {
            result.status = PtPhaseSetStatus::indeterminate;
            result.diagnostic = missing_candidate_diagnostic;
            break;
        }
        PtCandidatePhaseSet set;
        set.phases.reserve(3U);
        for (const auto& phase : candidate->phases) {
            set.phases.push_back(project_max3_owned_phase(phase));
        }
        result.candidate_phase_set = std::move(set);
        result.status = PtPhaseSetStatus::accepted;
        result.diagnostic = source.diagnostic;
        break;
    }
    case Status::higher_phase_count_or_wrong_candidate:
        result = project_two(source.base);
        result.status = PtPhaseSetStatus::phase_set_unstable;
        result.diagnostic = source.diagnostic;
        break;
    case Status::phase_boundary_unresolved:
    case Status::indeterminate:
        copy_feed();
        result.status = PtPhaseSetStatus::indeterminate;
        result.diagnostic = source.diagnostic;
        break;
    }
    return result;
}

} // namespace detail
} // namespace mpmc::flash

#endif // MPMC_FLASH_PT_PHASE_SET_HPP
