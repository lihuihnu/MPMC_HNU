#ifndef MPMC_FLASH_PT_PHASE_SET_HPP
#define MPMC_FLASH_PT_PHASE_SET_HPP

#include <mpmc/flash/pt_stability.hpp>

#include <cstddef>
#include <optional>
#include <string>
#include <string_view>
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

} // namespace mpmc::flash

#endif // MPMC_FLASH_PT_PHASE_SET_HPP
