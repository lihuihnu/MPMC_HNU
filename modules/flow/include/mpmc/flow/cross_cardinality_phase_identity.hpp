#ifndef MPMC_FLOW_CROSS_CARDINALITY_PHASE_IDENTITY_HPP
#define MPMC_FLOW_CROSS_CARDINALITY_PHASE_IDENTITY_HPP

#include <mpmc/flow/natural_variable_cell_state.hpp>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <optional>
#include <span>
#include <stdexcept>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

namespace mpmc::flow {

inline constexpr std::string_view
    cross_cardinality_phase_identity_convention =
        "flow/cross-cardinality-phase-identity/v1";

/// Opaque physical phase identity frozen by the caller for one nonlinear solve.
///
/// Neither natural-variable slot, provider vector position, compressibility
/// factor, density nor a liquid/vapor heuristic is a physical identity here.
/// The producer must preserve this key explicitly through continuation,
/// phase-set transition or another audited orchestration layer.
struct FrozenPhysicalPhaseIdentity {
    std::string provenance_scope;
    std::string opaque_phase_key;

    friend bool operator==(
        const FrozenPhysicalPhaseIdentity&,
        const FrozenPhysicalPhaseIdentity&) =
        default;
};

[[nodiscard]] inline bool
frozen_phase_identity_less(
    const FrozenPhysicalPhaseIdentity& left,
    const FrozenPhysicalPhaseIdentity& right) noexcept {
    return left.provenance_scope <
               right.provenance_scope ||
        (left.provenance_scope ==
             right.provenance_scope &&
         left.opaque_phase_key <
             right.opaque_phase_key);
}

/// Local slot -> explicit physical phase identity map.
///
/// The map is immutable for one frozen active-set solve. It deliberately does
/// not infer identity from PhaseSlot3 or phase-array position.
class FrozenActivePhaseIdentityMap {
public:
    explicit FrozenActivePhaseIdentityMap(
        std::vector<FrozenPhysicalPhaseIdentity>
            identities)
        : identities_(
              std::move(identities)) {
        validate();
    }

    [[nodiscard]] std::size_t
    phase_count() const noexcept {
        return identities_.size();
    }

    [[nodiscard]] std::span<
        const FrozenPhysicalPhaseIdentity>
    identities() const noexcept {
        return identities_;
    }

    [[nodiscard]] const FrozenPhysicalPhaseIdentity&
    identity(
        std::size_t active_phase_index) const {
        return identities_.at(
            active_phase_index);
    }

    [[nodiscard]] std::optional<std::size_t>
    find(
        const FrozenPhysicalPhaseIdentity&
            identity) const noexcept {
        const auto found =
            std::find(
                identities_.begin(),
                identities_.end(),
                identity);
        if (found == identities_.end()) {
            return std::nullopt;
        }
        return static_cast<std::size_t>(
            std::distance(
                identities_.begin(),
                found));
    }

private:
    void validate() const {
        if (identities_.empty() ||
            identities_.size() >
                fixed_three_phase_count) {
            throw std::invalid_argument(
                "mpmc::flow::FrozenActivePhaseIdentityMap: active phase count must be in [1,3]");
        }
        for (std::size_t index = 0U;
             index < identities_.size();
             ++index) {
            const auto& identity =
                identities_[index];
            if (identity.provenance_scope.empty() ||
                identity.opaque_phase_key.empty()) {
                throw std::invalid_argument(
                    "mpmc::flow::FrozenActivePhaseIdentityMap: phase identity requires nonempty provenance scope and opaque key");
            }
            for (std::size_t previous = 0U;
                 previous < index;
                 ++previous) {
                if (identities_[previous] ==
                    identity) {
                    throw std::invalid_argument(
                        "mpmc::flow::FrozenActivePhaseIdentityMap: duplicate physical phase identity");
                }
            }
        }
    }

    std::vector<FrozenPhysicalPhaseIdentity>
        identities_;
};

enum class CrossCardinalityPhasePresence {
    active_on_both_sides,
    owner_only_active,
    neighbour_only_active
};

struct CrossCardinalityFacePhaseIdentityBinding {
    FrozenPhysicalPhaseIdentity identity;
    std::optional<std::size_t>
        owner_active_phase_index;
    std::optional<std::size_t>
        neighbour_active_phase_index;
    CrossCardinalityPhasePresence presence{
        CrossCardinalityPhasePresence::
            active_on_both_sides};
};

/// Model-neutral phase pairing plan for one owner->neighbour face.
///
/// Binding order is canonical lexicographic order of the explicit opaque
/// identities, not local slot order. Reversing face orientation therefore does
/// not silently redefine which physical phases are compared.
class CrossCardinalityFacePhaseIdentityPlan {
public:
    CrossCardinalityFacePhaseIdentityPlan(
        std::size_t owner_phase_count,
        std::size_t neighbour_phase_count,
        std::vector<
            CrossCardinalityFacePhaseIdentityBinding>
            bindings,
        bool slot_aligned_same_active_set)
        : owner_phase_count_(
              owner_phase_count),
          neighbour_phase_count_(
              neighbour_phase_count),
          bindings_(
              std::move(bindings)),
          slot_aligned_same_active_set_(
              slot_aligned_same_active_set) {
        validate();
    }

    [[nodiscard]] std::size_t
    owner_phase_count() const noexcept {
        return owner_phase_count_;
    }

    [[nodiscard]] std::size_t
    neighbour_phase_count() const noexcept {
        return neighbour_phase_count_;
    }

    [[nodiscard]] std::span<
        const CrossCardinalityFacePhaseIdentityBinding>
    bindings() const noexcept {
        return bindings_;
    }

    [[nodiscard]] bool
    same_cardinality() const noexcept {
        return owner_phase_count_ ==
            neighbour_phase_count_;
    }

    [[nodiscard]] bool
    same_active_identity_set() const noexcept {
        return std::all_of(
            bindings_.begin(),
            bindings_.end(),
            [](const auto& binding) {
                return binding.presence ==
                    CrossCardinalityPhasePresence::
                        active_on_both_sides;
            });
    }

    /// True only when both sides have the same active identities and each
    /// identity occupies the same local slot. Existing fixed-cardinality TPFA
    /// kernels may be used directly only under this stronger condition.
    [[nodiscard]] bool
    slot_aligned_same_active_set() const noexcept {
        return slot_aligned_same_active_set_;
    }

    [[nodiscard]] bool
    requires_absent_phase_extension() const noexcept {
        return std::any_of(
            bindings_.begin(),
            bindings_.end(),
            [](const auto& binding) {
                return binding.presence !=
                    CrossCardinalityPhasePresence::
                        active_on_both_sides;
            });
    }

    [[nodiscard]]
    const CrossCardinalityFacePhaseIdentityBinding*
    find(
        const FrozenPhysicalPhaseIdentity&
            identity) const noexcept {
        const auto found =
            std::lower_bound(
                bindings_.begin(),
                bindings_.end(),
                identity,
                [](const auto& binding,
                   const auto& value) {
                    return frozen_phase_identity_less(
                        binding.identity,
                        value);
                });
        if (found == bindings_.end() ||
            found->identity != identity) {
            return nullptr;
        }
        return &*found;
    }

private:
    void validate() const {
        if (owner_phase_count_ == 0U ||
            owner_phase_count_ >
                fixed_three_phase_count ||
            neighbour_phase_count_ == 0U ||
            neighbour_phase_count_ >
                fixed_three_phase_count ||
            bindings_.empty()) {
            throw std::invalid_argument(
                "mpmc::flow::CrossCardinalityFacePhaseIdentityPlan: invalid phase cardinality");
        }

        std::size_t owner_seen = 0U;
        std::size_t neighbour_seen = 0U;
        for (std::size_t index = 0U;
             index < bindings_.size();
             ++index) {
            const auto& binding =
                bindings_[index];
            if (binding.identity.provenance_scope.empty() ||
                binding.identity.opaque_phase_key.empty() ||
                (index > 0U &&
                 !frozen_phase_identity_less(
                     bindings_[index - 1U].identity,
                     binding.identity))) {
                throw std::invalid_argument(
                    "mpmc::flow::CrossCardinalityFacePhaseIdentityPlan: bindings must contain unique canonical identities");
            }

            const bool owner =
                binding.owner_active_phase_index
                    .has_value();
            const bool neighbour =
                binding.neighbour_active_phase_index
                    .has_value();
            if (!owner && !neighbour) {
                throw std::invalid_argument(
                    "mpmc::flow::CrossCardinalityFacePhaseIdentityPlan: phase is absent on both sides");
            }
            if (owner &&
                *binding.owner_active_phase_index >=
                    owner_phase_count_) {
                throw std::invalid_argument(
                    "mpmc::flow::CrossCardinalityFacePhaseIdentityPlan: owner phase slot is out of range");
            }
            if (neighbour &&
                *binding.neighbour_active_phase_index >=
                    neighbour_phase_count_) {
                throw std::invalid_argument(
                    "mpmc::flow::CrossCardinalityFacePhaseIdentityPlan: neighbour phase slot is out of range");
            }

            const auto expected_presence =
                owner && neighbour
                    ? CrossCardinalityPhasePresence::
                          active_on_both_sides
                    : (owner
                           ? CrossCardinalityPhasePresence::
                                 owner_only_active
                           : CrossCardinalityPhasePresence::
                                 neighbour_only_active);
            if (binding.presence !=
                expected_presence) {
                throw std::invalid_argument(
                    "mpmc::flow::CrossCardinalityFacePhaseIdentityPlan: presence metadata disagrees with phase slots");
            }
            owner_seen +=
                owner ? 1U : 0U;
            neighbour_seen +=
                neighbour ? 1U : 0U;
        }

        if (owner_seen !=
                owner_phase_count_ ||
            neighbour_seen !=
                neighbour_phase_count_) {
            throw std::invalid_argument(
                "mpmc::flow::CrossCardinalityFacePhaseIdentityPlan: active phase coverage is incomplete");
        }

        const bool actually_slot_aligned =
            same_cardinality() &&
            same_active_identity_set() &&
            std::all_of(
                bindings_.begin(),
                bindings_.end(),
                [](const auto& binding) {
                    return binding
                               .owner_active_phase_index ==
                        binding
                            .neighbour_active_phase_index;
                });
        if (actually_slot_aligned !=
            slot_aligned_same_active_set_) {
            throw std::invalid_argument(
                "mpmc::flow::CrossCardinalityFacePhaseIdentityPlan: slot-alignment flag is inconsistent");
        }
    }

    std::size_t owner_phase_count_{};
    std::size_t neighbour_phase_count_{};
    std::vector<
        CrossCardinalityFacePhaseIdentityBinding>
        bindings_;
    bool slot_aligned_same_active_set_{};
};

[[nodiscard]] inline
CrossCardinalityFacePhaseIdentityPlan
make_cross_cardinality_face_phase_identity_plan(
    const FrozenActivePhaseIdentityMap& owner,
    const FrozenActivePhaseIdentityMap& neighbour) {
    std::vector<FrozenPhysicalPhaseIdentity>
        union_identities;
    union_identities.reserve(
        owner.phase_count() +
        neighbour.phase_count());
    union_identities.insert(
        union_identities.end(),
        owner.identities().begin(),
        owner.identities().end());
    union_identities.insert(
        union_identities.end(),
        neighbour.identities().begin(),
        neighbour.identities().end());
    std::sort(
        union_identities.begin(),
        union_identities.end(),
        frozen_phase_identity_less);
    union_identities.erase(
        std::unique(
            union_identities.begin(),
            union_identities.end()),
        union_identities.end());

    std::vector<
        CrossCardinalityFacePhaseIdentityBinding>
        bindings;
    bindings.reserve(
        union_identities.size());
    for (const auto& identity :
         union_identities) {
        const auto owner_slot =
            owner.find(identity);
        const auto neighbour_slot =
            neighbour.find(identity);
        bindings.push_back(
            {
                identity,
                owner_slot,
                neighbour_slot,
                owner_slot && neighbour_slot
                    ? CrossCardinalityPhasePresence::
                          active_on_both_sides
                    : (owner_slot
                           ? CrossCardinalityPhasePresence::
                                 owner_only_active
                           : CrossCardinalityPhasePresence::
                                 neighbour_only_active)});
    }

    bool slot_aligned =
        owner.phase_count() ==
        neighbour.phase_count();
    if (slot_aligned) {
        for (std::size_t slot = 0U;
             slot < owner.phase_count();
             ++slot) {
            if (owner.identity(slot) !=
                neighbour.identity(slot)) {
                slot_aligned = false;
                break;
            }
        }
    }

    return {
        owner.phase_count(),
        neighbour.phase_count(),
        std::move(bindings),
        slot_aligned};
}

/// Frozen semantics for a phase that is inactive on one side of a face.
///
/// - mobility is structurally zero and has an all-zero derivative on the
///   inactive side;
/// - phase pressure and mass density are *not* zero, copied, or inferred;
///   they require an explicit smooth extension before a two-sided potential
///   can be formed;
/// - composition/enthalpy of an inactive phase are not defined. A nonzero
///   advective payload must therefore come from an actually active upwind
///   phase. If the inactive side is selected as upwind, structural-zero
///   mobility makes that phase's advective flux zero.
struct CrossCardinalityAbsentPhaseSemantics {
    static constexpr double
        mobility_per_pa_s = 0.0;
    static constexpr bool
        mobility_gradient_is_structural_zero =
            true;
    static constexpr bool
        phase_pressure_requires_explicit_extension =
            true;
    static constexpr bool
        mass_density_requires_explicit_extension =
            true;
    static constexpr bool
        absent_composition_is_undefined =
            true;
    static constexpr bool
        absent_enthalpy_is_undefined =
            true;
    static constexpr bool
        nonzero_advective_payload_requires_active_upwind_phase =
            true;
};

/// Explicit hypothetical phase-potential state on the side where the phase is
/// inactive. This is *not* an active-phase state and carries no mobility,
/// composition, viscosity or caloric payload.
class AbsentPhasePotentialExtensionLinearization {
public:
    AbsentPhasePotentialExtensionLinearization(
        FrozenPhysicalPhaseIdentity identity,
        double phase_pressure_pa,
        double mass_density_kg_per_m3,
        std::size_t input_count,
        std::vector<double>
            phase_pressure_gradient,
        std::vector<double>
            mass_density_gradient,
        std::string provenance)
        : identity_(
              std::move(identity)),
          phase_pressure_pa_(
              phase_pressure_pa),
          mass_density_kg_per_m3_(
              mass_density_kg_per_m3),
          input_count_(
              input_count),
          phase_pressure_gradient_(
              std::move(
                  phase_pressure_gradient)),
          mass_density_gradient_(
              std::move(
                  mass_density_gradient)),
          provenance_(
              std::move(provenance)) {
        validate();
    }

    [[nodiscard]] const FrozenPhysicalPhaseIdentity&
    identity() const noexcept {
        return identity_;
    }

    [[nodiscard]] double
    phase_pressure_pa() const noexcept {
        return phase_pressure_pa_;
    }

    [[nodiscard]] double
    mass_density_kg_per_m3() const noexcept {
        return mass_density_kg_per_m3_;
    }

    [[nodiscard]] std::size_t
    input_count() const noexcept {
        return input_count_;
    }

    [[nodiscard]] std::span<const double>
    phase_pressure_gradient() const noexcept {
        return phase_pressure_gradient_;
    }

    [[nodiscard]] std::span<const double>
    mass_density_gradient() const noexcept {
        return mass_density_gradient_;
    }

    [[nodiscard]] std::string_view
    provenance() const noexcept {
        return provenance_;
    }

private:
    void validate() const {
        if (identity_.provenance_scope.empty() ||
            identity_.opaque_phase_key.empty() ||
            !std::isfinite(
                phase_pressure_pa_) ||
            !(phase_pressure_pa_ > 0.0) ||
            !std::isfinite(
                mass_density_kg_per_m3_) ||
            !(mass_density_kg_per_m3_ > 0.0) ||
            input_count_ == 0U ||
            phase_pressure_gradient_.size() !=
                input_count_ ||
            mass_density_gradient_.size() !=
                input_count_ ||
            provenance_.empty()) {
            throw std::invalid_argument(
                "mpmc::flow::AbsentPhasePotentialExtensionLinearization: invalid explicit phase-potential extension");
        }
        for (const double value :
             phase_pressure_gradient_) {
            if (!std::isfinite(value)) {
                throw std::invalid_argument(
                    "mpmc::flow::AbsentPhasePotentialExtensionLinearization: pressure gradient contains non-finite derivative");
            }
        }
        for (const double value :
             mass_density_gradient_) {
            if (!std::isfinite(value)) {
                throw std::invalid_argument(
                    "mpmc::flow::AbsentPhasePotentialExtensionLinearization: density gradient contains non-finite derivative");
            }
        }
    }

    FrozenPhysicalPhaseIdentity identity_;
    double phase_pressure_pa_{};
    double mass_density_kg_per_m3_{};
    std::size_t input_count_{};
    std::vector<double>
        phase_pressure_gradient_;
    std::vector<double>
        mass_density_gradient_;
    std::string provenance_;
};

} // namespace mpmc::flow

#endif // MPMC_FLOW_CROSS_CARDINALITY_PHASE_IDENTITY_HPP
