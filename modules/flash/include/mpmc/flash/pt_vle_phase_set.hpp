#ifndef MPMC_FLASH_PT_VLE_PHASE_SET_HPP
#define MPMC_FLASH_PT_VLE_PHASE_SET_HPP

#include <mpmc/flash/pt_phase_set.hpp>
#include <mpmc/flash/pt_split.hpp>

#include <optional>
#include <string>
#include <utility>

namespace mpmc::flash {
namespace detail {

[[nodiscard]] inline PtCandidatePhase project_vle_phase(
    double mole_phase_fraction, const std::vector<double>& composition,
    const PtSplitPhase& phase) {
    PtCandidatePhase projected;
    projected.mole_phase_fraction = mole_phase_fraction;
    projected.composition = composition;
    projected.activity = phase.activity;
    projected.compressibility_factor = phase.z;
    return projected;
}

[[nodiscard]] inline std::optional<PtCandidatePhaseSet> project_vle_pair(
    const PtSplitResult& source) {
    const auto* candidate = source.candidate();
    if (candidate == nullptr) {
        return std::nullopt;
    }
    PtCandidatePhaseSet projected;
    projected.phases.reserve(2);
    const double vapor_fraction = candidate->fractions.vapor_fraction;
    // Legacy VLE stores beta_V explicitly and defines beta_L = 1-beta_V.
    // The projection performs only that contract-prescribed complement; it does
    // not rerun thermodynamic properties, RR, stability, or equilibrium.
    projected.phases.push_back(project_vle_phase(
        1.0 - vapor_fraction, candidate->fractions.liquid, candidate->liquid));
    projected.phases.push_back(project_vle_phase(
        vapor_fraction, candidate->fractions.vapor, candidate->vapor));
    return projected;
}

} // namespace detail

// Structural projection of the existing VLE result into the generic phase-set
// contract. It performs no provider calls, root solves, normalization, phase
// relabeling, Rachford-Rice solve, stability search, or equilibrium iteration.
[[nodiscard]] inline PtPhaseSetResult project_pt_vle_phase_set(
    const PtSplitResult& source) {
    PtPhaseSetResult result;
    result.capability.maximum_phase_count = 2;
    result.pressure_pa = source.initial_stability.pressure_pa;
    result.temperature_k = source.initial_stability.temperature_k;
    result.feed = source.initial_stability.feed;
    result.global_stability_proven = source.global_stability_proven;
    result.diagnostic = source.diagnostic;

    switch (source.status) {
    case PtSplitStatus::single_phase_no_instability_found: {
        if (!source.initial_stability.reference) {
            result.status = PtPhaseSetStatus::indeterminate;
            result.diagnostic =
                "PT phase-set projection: accepted single-phase status lacks reference phase";
            return result;
        }
        PtCandidatePhase phase;
        phase.mole_phase_fraction = 1.0;
        phase.composition = source.initial_stability.feed;
        phase.activity = *source.initial_stability.reference;
        // Generic stability providers do not expose a compressibility factor.
        phase.compressibility_factor.reset();
        PtCandidatePhaseSet set;
        set.phases.push_back(std::move(phase));
        result.candidate_phase_set = std::move(set);
        result.status = PtPhaseSetStatus::accepted;
        return result;
    }
    case PtSplitStatus::two_phase_no_instability_found:
        result.candidate_phase_set = detail::project_vle_pair(source);
        if (!result.candidate_phase_set) {
            result.status = PtPhaseSetStatus::indeterminate;
            result.diagnostic =
                "PT phase-set projection: accepted two-phase status lacks converged pair";
            return result;
        }
        result.status = PtPhaseSetStatus::accepted;
        return result;
    case PtSplitStatus::phase_set_unstable:
        result.candidate_phase_set = detail::project_vle_pair(source);
        result.status = PtPhaseSetStatus::phase_set_unstable;
        return result;
    case PtSplitStatus::indeterminate:
        // A converged pair may exist while final stability is unresolved. Retain
        // it as a candidate, but never expose it through accepted_phase_set().
        result.candidate_phase_set = detail::project_vle_pair(source);
        result.status = PtPhaseSetStatus::indeterminate;
        return result;
    }
    result.status = PtPhaseSetStatus::indeterminate;
    result.diagnostic = "PT phase-set projection: unknown legacy VLE status";
    return result;
}

} // namespace mpmc::flash

#endif // MPMC_FLASH_PT_VLE_PHASE_SET_HPP
