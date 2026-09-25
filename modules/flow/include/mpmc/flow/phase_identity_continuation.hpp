#ifndef MPMC_FLOW_PHASE_IDENTITY_CONTINUATION_HPP
#define MPMC_FLOW_PHASE_IDENTITY_CONTINUATION_HPP

#include <mpmc/flow/cross_cardinality_phase_identity.hpp>
#include <mpmc/flow/phase_set_transition.hpp>

#include <algorithm>
#include <cstddef>
#include <span>
#include <stdexcept>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

namespace mpmc::flow {

inline constexpr std::string_view
    phase_identity_continuation_convention =
        "flow/phase-set-transition/physical-phase-identity-continuation/v1";

class PhaseIdentityContinuationSnapshot {
public:
    PhaseIdentityContinuationSnapshot(
        std::size_t source_phase_count,
        std::size_t target_phase_count,
        FrozenActivePhaseIdentityMap source,
        FrozenActivePhaseIdentityMap target,
        PhaseSetTransitionTrigger trigger,
        std::string evidence_profile,
        std::string diagnostic)
        : source_phase_count_(source_phase_count),
          target_phase_count_(target_phase_count),
          source_(std::move(source)),
          target_(std::move(target)),
          trigger_(trigger),
          evidence_profile_(std::move(evidence_profile)),
          diagnostic_(std::move(diagnostic)) {
        validate();
    }

    [[nodiscard]] std::size_t source_phase_count() const noexcept {
        return source_phase_count_;
    }
    [[nodiscard]] std::size_t target_phase_count() const noexcept {
        return target_phase_count_;
    }
    [[nodiscard]] const FrozenActivePhaseIdentityMap&
    source() const noexcept { return source_; }
    [[nodiscard]] const FrozenActivePhaseIdentityMap&
    target() const noexcept { return target_; }
    [[nodiscard]] PhaseSetTransitionTrigger trigger() const noexcept {
        return trigger_;
    }
    [[nodiscard]] std::string_view evidence_profile() const noexcept {
        return evidence_profile_;
    }
    [[nodiscard]] std::string_view diagnostic() const noexcept {
        return diagnostic_;
    }

    [[nodiscard]] std::vector<FrozenPhysicalPhaseIdentity>
    phase_universe() const {
        std::vector<FrozenPhysicalPhaseIdentity> result;
        result.reserve(source_.phase_count() + target_.phase_count());
        result.insert(
            result.end(),
            source_.identities().begin(),
            source_.identities().end());
        result.insert(
            result.end(),
            target_.identities().begin(),
            target_.identities().end());
        std::sort(
            result.begin(),
            result.end(),
            frozen_phase_identity_less);
        result.erase(
            std::unique(result.begin(), result.end()),
            result.end());
        return result;
    }

private:
    void validate() const {
        if (source_phase_count_ == 0U ||
            source_phase_count_ > fixed_three_phase_count ||
            target_phase_count_ == 0U ||
            target_phase_count_ > fixed_three_phase_count ||
            source_phase_count_ == target_phase_count_ ||
            source_.phase_count() != source_phase_count_ ||
            target_.phase_count() != target_phase_count_ ||
            evidence_profile_.empty()) {
            throw std::invalid_argument(
                "mpmc::flow::PhaseIdentityContinuationSnapshot: malformed transition identity continuation");
        }

        if (target_phase_count_ > source_phase_count_) {
            for (const auto& identity : source_.identities()) {
                if (!target_.find(identity).has_value()) {
                    throw std::invalid_argument(
                        "mpmc::flow::PhaseIdentityContinuationSnapshot: phase appearance cannot silently rename an existing physical phase");
                }
            }
        } else {
            for (const auto& identity : target_.identities()) {
                if (!source_.find(identity).has_value()) {
                    throw std::invalid_argument(
                        "mpmc::flow::PhaseIdentityContinuationSnapshot: phase disappearance cannot introduce a new physical phase identity");
                }
            }
        }
    }

    std::size_t source_phase_count_{};
    std::size_t target_phase_count_{};
    FrozenActivePhaseIdentityMap source_;
    FrozenActivePhaseIdentityMap target_;
    PhaseSetTransitionTrigger trigger_{
        PhaseSetTransitionTrigger::provider_topology_witness};
    std::string evidence_profile_;
    std::string diagnostic_;
};

[[nodiscard]] inline PhaseIdentityContinuationSnapshot
make_phase_identity_continuation_snapshot(
    const PhaseSetTransitionCandidate& candidate,
    FrozenActivePhaseIdentityMap source,
    FrozenActivePhaseIdentityMap target) {
    if (candidate.status !=
            PhaseSetTransitionCandidateStatus::target_resolved ||
        candidate.source_phase_count != source.phase_count() ||
        candidate.target_phase_count != target.phase_count() ||
        candidate.evidence_profile.empty()) {
        throw std::invalid_argument(
            "mpmc::flow: phase identity continuation requires a freshly resolved transition candidate and matching explicit identity maps");
    }
    return {
        candidate.source_phase_count,
        candidate.target_phase_count,
        std::move(source),
        std::move(target),
        candidate.trigger,
        candidate.evidence_profile,
        candidate.diagnostic};
}

} // namespace mpmc::flow

#endif // MPMC_FLOW_PHASE_IDENTITY_CONTINUATION_HPP
