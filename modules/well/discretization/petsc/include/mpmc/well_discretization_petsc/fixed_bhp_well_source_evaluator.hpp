#ifndef MPMC_WELL_DISCRETIZATION_PETSC_FIXED_BHP_WELL_SOURCE_EVALUATOR_HPP
#define MPMC_WELL_DISCRETIZATION_PETSC_FIXED_BHP_WELL_SOURCE_EVALUATOR_HPP

#include <mpmc/flow/cross_cardinality_phase_identity.hpp>
#include <mpmc/flow_discretization_petsc/mixed_cardinality_physical_snes_assembly.hpp>
#include <mpmc/well_discretization/energy_rate.hpp>
#include <mpmc/well_discretization/fixed_bhp_connection_source.hpp>

#include <algorithm>
#include <cmath>
#include <optional>
#include <span>
#include <stdexcept>
#include <string>
#include <string_view>
#include <utility>
#include <variant>
#include <vector>

namespace mpmc::well_discretization_petsc {

inline constexpr std::string_view
    fixed_bhp_peaceman_well_source_evaluator_convention =
        "well-discretization-petsc/fixed-bhp-single-connection/phase-identity-rebindable-peaceman/v3";

/// Backward-compatible name retained for callers of the first 3P-only bridge.
inline constexpr std::string_view
    fixed_bhp_three_phase_peaceman_well_source_evaluator_convention =
        fixed_bhp_peaceman_well_source_evaluator_convention;

/// One well-side injection-enthalpy value keyed by stable physical phase
/// identity rather than by the current local active-phase slot.
struct FixedBhpPhaseIdentitySpecificEnthalpy3D {
    mpmc::flow::FrozenPhysicalPhaseIdentity
        identity;
    double specific_enthalpy_j_per_kg{};
};

/// Stable phase-identity registry retained across phase-set rebuilds.
///
/// The registry may include currently inactive phases. Rebinding resolves the
/// current FrozenActivePhaseIdentityMap to an exact slot-ordered
/// FixedBhpInjectionEnthalpy3D without inferring identity from slot position.
struct FixedBhpPhaseIdentityInjectionEnthalpy3D {
    std::string provenance;
    std::vector<
        FixedBhpPhaseIdentitySpecificEnthalpy3D>
        phases;
};

namespace fixed_bhp_phase_identity_detail {

inline void validate_registry(
    const FixedBhpPhaseIdentityInjectionEnthalpy3D&
        registry) {
    if (registry.provenance.empty() ||
        registry.phases.empty() ||
        registry.phases.size() > 3U) {
        throw std::invalid_argument(
            "mpmc::well_discretization_petsc: invalid fixed-BHP phase-identity injection-enthalpy registry");
    }
    for (std::size_t index = 0U;
         index < registry.phases.size();
         ++index) {
        const auto& entry =
            registry.phases[index];
        if (entry.identity.provenance_scope.empty() ||
            entry.identity.opaque_phase_key.empty() ||
            !std::isfinite(
                entry.specific_enthalpy_j_per_kg)) {
            throw std::invalid_argument(
                "mpmc::well_discretization_petsc: malformed fixed-BHP phase-identity enthalpy binding");
        }
        for (std::size_t previous = 0U;
             previous < index;
             ++previous) {
            if (registry.phases[previous].identity ==
                entry.identity) {
                throw std::invalid_argument(
                    "mpmc::well_discretization_petsc: duplicate fixed-BHP phase identity in injection-enthalpy registry");
            }
        }
    }
}

[[nodiscard]] inline
mpmc::well_discretization::
    FixedBhpInjectionEnthalpy3D
resolve_injection_enthalpy(
    const FixedBhpPhaseIdentityInjectionEnthalpy3D&
        registry,
    const mpmc::flow::
        FrozenActivePhaseIdentityMap&
            active_phases) {
    validate_registry(registry);

    mpmc::well_discretization::
        FixedBhpInjectionEnthalpy3D
        result;
    result.provenance =
        registry.provenance;
    result.specific_enthalpy_j_per_kg.reserve(
        active_phases.phase_count());

    for (const auto& identity :
         active_phases.identities()) {
        const auto found =
            std::find_if(
                registry.phases.begin(),
                registry.phases.end(),
                [&](const auto& entry) {
                    return entry.identity ==
                        identity;
                });
        if (found == registry.phases.end()) {
            throw std::invalid_argument(
                "mpmc::well_discretization_petsc: fixed-BHP phase-identity registry does not cover one active phase");
        }
        result.specific_enthalpy_j_per_kg.push_back(
            found->specific_enthalpy_j_per_kg);
    }
    return result;
}

} // namespace fixed_bhp_phase_identity_detail

[[nodiscard]] inline
mpmc::well_discretization::
    FixedBhpInjectionEnthalpy3D
resolve_fixed_bhp_injection_enthalpy_by_phase_identity_3d(
    const FixedBhpPhaseIdentityInjectionEnthalpy3D&
        registry,
    const mpmc::flow::
        FrozenActivePhaseIdentityMap&
            active_phases) {
    return fixed_bhp_phase_identity_detail::
        resolve_injection_enthalpy(
            registry,
            active_phases);
}

/// Immutable configuration for one fixed-BHP Peaceman connection.
///
/// Legacy construction may still provide slot-ordered enthalpy for one frozen
/// chart. Phase-transition-aware construction instead retains a stable
/// physical-phase registry plus the currently bound active phase map. Rebinding
/// keeps stable cell identity, Peaceman geometry/properties, BHP and source
/// provenance unchanged while resolving a new active-slot enthalpy vector.
class FixedBhpPeacemanWellSourceEvaluatorContext3D {
public:
    static constexpr std::string_view convention =
        fixed_bhp_peaceman_well_source_evaluator_convention;

    [[nodiscard]] static
    FixedBhpPeacemanWellSourceEvaluatorContext3D
    create(
        mpmc::mesh::GlobalEntityId target_cell_global,
        const mpmc::well::PeacemanWellIndex3D&
            connection,
        double bottom_hole_pressure_pa,
        mpmc::well_discretization::
            FixedBhpInjectionEnthalpy3D
                injection_enthalpy,
        std::string source_provenance) {
        validate_common(
            bottom_hole_pressure_pa,
            injection_enthalpy,
            source_provenance);
        auto normalized_connection =
            normalize_connection(connection);

        return FixedBhpPeacemanWellSourceEvaluatorContext3D{
            target_cell_global,
            std::move(normalized_connection),
            bottom_hole_pressure_pa,
            std::move(injection_enthalpy),
            std::move(source_provenance),
            std::nullopt,
            std::nullopt};
    }

    [[nodiscard]] static
    FixedBhpPeacemanWellSourceEvaluatorContext3D
    create_phase_identity_bound(
        mpmc::mesh::GlobalEntityId target_cell_global,
        const mpmc::well::PeacemanWellIndex3D&
            connection,
        double bottom_hole_pressure_pa,
        FixedBhpPhaseIdentityInjectionEnthalpy3D
            phase_identity_injection_enthalpy,
        mpmc::flow::FrozenActivePhaseIdentityMap
            active_phases,
        std::string source_provenance) {
        auto injection_enthalpy =
            resolve_fixed_bhp_injection_enthalpy_by_phase_identity_3d(
                phase_identity_injection_enthalpy,
                active_phases);
        validate_common(
            bottom_hole_pressure_pa,
            injection_enthalpy,
            source_provenance);
        auto normalized_connection =
            normalize_connection(connection);

        return FixedBhpPeacemanWellSourceEvaluatorContext3D{
            target_cell_global,
            std::move(normalized_connection),
            bottom_hole_pressure_pa,
            std::move(injection_enthalpy),
            std::move(source_provenance),
            std::move(
                phase_identity_injection_enthalpy),
            std::move(active_phases)};
    }

    /// Compatibility factory for the former 3P-only injection payload.
    [[nodiscard]] static
    FixedBhpPeacemanWellSourceEvaluatorContext3D
    create(
        mpmc::mesh::GlobalEntityId target_cell_global,
        const mpmc::well::PeacemanWellIndex3D&
            connection,
        double bottom_hole_pressure_pa,
        mpmc::well_discretization::
            InjectionPhaseSpecificEnthalpy3P
                injection_enthalpy,
        std::string source_provenance) {
        return create(
            target_cell_global,
            connection,
            bottom_hole_pressure_pa,
            {
                std::move(
                    injection_enthalpy
                        .provenance),
                std::vector<double>{
                    injection_enthalpy
                        .specific_enthalpy_j_per_kg
                        .begin(),
                    injection_enthalpy
                        .specific_enthalpy_j_per_kg
                        .end()}},
            std::move(source_provenance));
    }

    [[nodiscard]]
    FixedBhpPeacemanWellSourceEvaluatorContext3D
    rebind(
        mpmc::mesh::GlobalEntityId target_cell_global,
        const mpmc::flow::
            FrozenActivePhaseIdentityMap&
                target_active_phases) const {
        if (target_cell_global !=
                target_cell_global_ ||
            !phase_identity_injection_enthalpy_
                 .has_value()) {
            throw std::invalid_argument(
                "mpmc::well_discretization_petsc: fixed-BHP completion rebinding requires the same stable cell and a phase-identity registry");
        }
        return create_phase_identity_bound(
            target_cell_global_,
            connection_,
            bottom_hole_pressure_pa_,
            *phase_identity_injection_enthalpy_,
            target_active_phases,
            source_provenance_);
    }

    [[nodiscard]] mpmc::mesh::GlobalEntityId
    target_cell_global() const noexcept {
        return target_cell_global_;
    }

    [[nodiscard]] const
    mpmc::well::PeacemanWellIndex3D&
    connection() const noexcept {
        return connection_;
    }

    [[nodiscard]] double
    bottom_hole_pressure_pa() const noexcept {
        return bottom_hole_pressure_pa_;
    }

    [[nodiscard]] const
    mpmc::well_discretization::
        FixedBhpInjectionEnthalpy3D&
    injection_enthalpy() const noexcept {
        return injection_enthalpy_;
    }

    [[nodiscard]] bool
    phase_identity_rebindable() const noexcept {
        return phase_identity_injection_enthalpy_
            .has_value();
    }

    [[nodiscard]] const std::optional<
        FixedBhpPhaseIdentityInjectionEnthalpy3D>&
    phase_identity_injection_enthalpy() const noexcept {
        return phase_identity_injection_enthalpy_;
    }

    [[nodiscard]] const std::optional<
        mpmc::flow::
            FrozenActivePhaseIdentityMap>&
    active_phase_identities() const noexcept {
        return active_phase_identities_;
    }

    [[nodiscard]] const std::string&
    source_provenance() const noexcept {
        return source_provenance_;
    }

private:
    static void validate_common(
        double bottom_hole_pressure_pa,
        const mpmc::well_discretization::
            FixedBhpInjectionEnthalpy3D&
                injection_enthalpy,
        const std::string& source_provenance) {
        if (!std::isfinite(bottom_hole_pressure_pa) ||
            !(bottom_hole_pressure_pa > 0.0) ||
            source_provenance.empty() ||
            injection_enthalpy.provenance.empty() ||
            injection_enthalpy
                    .specific_enthalpy_j_per_kg
                    .empty() ||
            injection_enthalpy
                    .specific_enthalpy_j_per_kg
                    .size() >
                3U) {
            throw std::invalid_argument(
                "mpmc::well_discretization_petsc: invalid fixed-BHP well-source configuration");
        }
        for (double value :
             injection_enthalpy
                 .specific_enthalpy_j_per_kg) {
            if (!std::isfinite(value)) {
                throw std::invalid_argument(
                    "mpmc::well_discretization_petsc: injection enthalpy must be finite [J/kg]");
            }
        }
    }

    [[nodiscard]] static
    mpmc::well::PeacemanWellIndex3D
    normalize_connection(
        const mpmc::well::PeacemanWellIndex3D&
            connection) {
        return mpmc::well::
            make_peaceman_well_index_3d(
                connection.cell_dimensions_m,
                connection.permeability_m2,
                connection.well_direction,
                connection.wellbore_radius_m,
                connection.skin_factor);
    }

    FixedBhpPeacemanWellSourceEvaluatorContext3D(
        mpmc::mesh::GlobalEntityId target_cell_global,
        mpmc::well::PeacemanWellIndex3D connection,
        double bottom_hole_pressure_pa,
        mpmc::well_discretization::
            FixedBhpInjectionEnthalpy3D
                injection_enthalpy,
        std::string source_provenance,
        std::optional<
            FixedBhpPhaseIdentityInjectionEnthalpy3D>
            phase_identity_injection_enthalpy,
        std::optional<
            mpmc::flow::
                FrozenActivePhaseIdentityMap>
            active_phase_identities)
        : target_cell_global_(
              target_cell_global),
          connection_(
              std::move(connection)),
          bottom_hole_pressure_pa_(
              bottom_hole_pressure_pa),
          injection_enthalpy_(
              std::move(injection_enthalpy)),
          source_provenance_(
              std::move(source_provenance)),
          phase_identity_injection_enthalpy_(
              std::move(
                  phase_identity_injection_enthalpy)),
          active_phase_identities_(
              std::move(
                  active_phase_identities)) {}

    mpmc::mesh::GlobalEntityId
        target_cell_global_;
    mpmc::well::PeacemanWellIndex3D
        connection_;
    double bottom_hole_pressure_pa_{};
    mpmc::well_discretization::
        FixedBhpInjectionEnthalpy3D
            injection_enthalpy_;
    std::string source_provenance_;
    std::optional<
        FixedBhpPhaseIdentityInjectionEnthalpy3D>
        phase_identity_injection_enthalpy_;
    std::optional<
        mpmc::flow::FrozenActivePhaseIdentityMap>
        active_phase_identities_;
};

[[nodiscard]] inline
FixedBhpPeacemanWellSourceEvaluatorContext3D
rebind_fixed_bhp_peaceman_well_source_context_3d(
    const FixedBhpPeacemanWellSourceEvaluatorContext3D&
        context,
    mpmc::mesh::GlobalEntityId target_cell_global,
    const mpmc::flow::
        FrozenActivePhaseIdentityMap&
            target_active_phases) {
    return context.rebind(
        target_cell_global,
        target_active_phases);
}

using FixedBhpThreePhasePeacemanWellSourceEvaluatorContext3D =
    FixedBhpPeacemanWellSourceEvaluatorContext3D;

inline constexpr std::string_view
    fixed_bhp_multi_connection_well_source_evaluator_convention =
        "well-discretization-petsc/fixed-bhp-single-well/multi-connection-owner-only/v1";

/// One logical fixed-BHP well spanning multiple stable reservoir cells.
///
/// The context is immutable after construction. Connections are sorted by
/// stable cell GlobalEntityId for deterministic lookup. A stable cell may
/// appear at most once, and every connection must carry exactly the same BHP.
/// This prevents one physical completion from being assembled twice and makes
/// the shared control parameter explicit without introducing a well unknown.
class FixedBhpMultiConnectionWellSourceEvaluatorContext3D {
public:
    static constexpr std::string_view convention =
        fixed_bhp_multi_connection_well_source_evaluator_convention;

    [[nodiscard]] static
    FixedBhpMultiConnectionWellSourceEvaluatorContext3D
    create(
        std::string well_id,
        std::vector<
            FixedBhpPeacemanWellSourceEvaluatorContext3D>
            connections) {
        if (well_id.empty() ||
            connections.size() < 2U) {
            throw std::invalid_argument(
                "mpmc::well_discretization_petsc: multi-connection fixed-BHP well requires nonempty well id and at least two connections");
        }

        const double bottom_hole_pressure_pa =
            connections.front()
                .bottom_hole_pressure_pa();
        if (!std::isfinite(bottom_hole_pressure_pa) ||
            !(bottom_hole_pressure_pa > 0.0)) {
            throw std::invalid_argument(
                "mpmc::well_discretization_petsc: multi-connection fixed-BHP well has invalid shared BHP");
        }

        std::sort(
            connections.begin(),
            connections.end(),
            [](const auto& left,
               const auto& right) {
                return left
                           .target_cell_global()
                           .value() <
                    right
                        .target_cell_global()
                        .value();
            });

        for (std::size_t index = 0U;
             index < connections.size();
             ++index) {
            if (connections[index]
                    .bottom_hole_pressure_pa() !=
                bottom_hole_pressure_pa) {
                throw std::invalid_argument(
                    "mpmc::well_discretization_petsc: multi-connection fixed-BHP well connections do not share one frozen BHP");
            }
            if (index > 0U &&
                connections[index - 1U]
                        .target_cell_global() ==
                    connections[index]
                        .target_cell_global()) {
                throw std::invalid_argument(
                    "mpmc::well_discretization_petsc: duplicate stable-cell completion in multi-connection fixed-BHP well");
            }
        }

        return FixedBhpMultiConnectionWellSourceEvaluatorContext3D{
            std::move(well_id),
            bottom_hole_pressure_pa,
            std::move(connections)};
    }

    [[nodiscard]] std::string_view
    well_id() const noexcept {
        return well_id_;
    }

    [[nodiscard]] double
    bottom_hole_pressure_pa() const noexcept {
        return bottom_hole_pressure_pa_;
    }

    [[nodiscard]] std::size_t
    connection_count() const noexcept {
        return connections_.size();
    }

    [[nodiscard]] std::span<
        const FixedBhpPeacemanWellSourceEvaluatorContext3D>
    connections() const noexcept {
        return connections_;
    }

    [[nodiscard]]
    const FixedBhpPeacemanWellSourceEvaluatorContext3D*
    find_connection(
        mpmc::mesh::GlobalEntityId
            cell_global) const noexcept {
        const auto found =
            std::lower_bound(
                connections_.begin(),
                connections_.end(),
                cell_global.value(),
                [](const auto& connection,
                   std::uint64_t value) {
                    return connection
                               .target_cell_global()
                               .value() <
                        value;
                });
        return found != connections_.end() &&
                       found->target_cell_global() ==
                           cell_global
            ? &*found
            : nullptr;
    }

    [[nodiscard]]
    FixedBhpPeacemanWellSourceEvaluatorContext3D*
    find_connection(
        mpmc::mesh::GlobalEntityId
            cell_global) noexcept {
        const auto found =
            std::lower_bound(
                connections_.begin(),
                connections_.end(),
                cell_global.value(),
                [](const auto& connection,
                   std::uint64_t value) {
                    return connection
                               .target_cell_global()
                               .value() <
                        value;
                });
        return found != connections_.end() &&
                       found->target_cell_global() ==
                           cell_global
            ? &*found
            : nullptr;
    }

private:
    FixedBhpMultiConnectionWellSourceEvaluatorContext3D(
        std::string well_id,
        double bottom_hole_pressure_pa,
        std::vector<
            FixedBhpPeacemanWellSourceEvaluatorContext3D>
            connections)
        : well_id_(std::move(well_id)),
          bottom_hole_pressure_pa_(
              bottom_hole_pressure_pa),
          connections_(
              std::move(connections)) {}

    std::string well_id_;
    double bottom_hole_pressure_pa_{};
    std::vector<
        FixedBhpPeacemanWellSourceEvaluatorContext3D>
        connections_;
};

namespace fixed_bhp_well_source_evaluator_detail {

[[nodiscard]] inline
mpmc::well_discretization::
    FixedBhpActivePhaseSourceInput3D
make_source_input(
    const mpmc::flow_discretization_petsc::
        SinglePhaseCurrentCellLinearization3D&
            current) {
    const auto mobility =
        mpmc::flow::
            build_single_phase_mobility_linearization(
                current.state,
                current.transport);

    return {
        mobility.state_identity,
        {mobility.phase_pressure_pa},
        {mobility.phase_pressure_gradient},
        {mobility.mobility_per_pa_s},
        {mobility.mobility_gradient},
        {current.molar_density
             .molar_density_mol_per_m3},
        {current.molar_density.gradient},
        {current.transport
             .mass_density_kg_per_m3},
        {current.transport
             .mass_density_gradient},
        {current.caloric
             .specific_enthalpy_j_per_kg},
        {current.caloric
             .specific_enthalpy_gradient}};
}

[[nodiscard]] inline
mpmc::well_discretization::
    FixedBhpActivePhaseSourceInput3D
make_source_input(
    const mpmc::flow_discretization_petsc::
        TwoPhaseCurrentCellLinearization3D&
            current) {
    const auto mobility =
        mpmc::flow::
            build_two_phase_mobility_linearization(
                current.state,
                current.transport);

    return {
        mobility.state_identity,
        std::vector<double>{
            mobility.phase_pressure_pa.begin(),
            mobility.phase_pressure_pa.end()},
        std::vector<std::vector<double>>{
            mobility.phase_pressure_gradient.begin(),
            mobility.phase_pressure_gradient.end()},
        std::vector<double>{
            mobility.mobility_per_pa_s.begin(),
            mobility.mobility_per_pa_s.end()},
        std::vector<std::vector<double>>{
            mobility.mobility_gradient.begin(),
            mobility.mobility_gradient.end()},
        std::vector<double>{
            current.molar_density
                .molar_density_mol_per_m3.begin(),
            current.molar_density
                .molar_density_mol_per_m3.end()},
        std::vector<std::vector<double>>{
            current.molar_density.gradient.begin(),
            current.molar_density.gradient.end()},
        std::vector<double>{
            current.transport
                .mass_density_kg_per_m3.begin(),
            current.transport
                .mass_density_kg_per_m3.end()},
        std::vector<std::vector<double>>{
            current.transport
                .mass_density_gradient.begin(),
            current.transport
                .mass_density_gradient.end()},
        std::vector<double>{
            current.caloric
                .specific_enthalpy_j_per_kg.begin(),
            current.caloric
                .specific_enthalpy_j_per_kg.end()},
        std::vector<std::vector<double>>{
            current.caloric
                .specific_enthalpy_gradient.begin(),
            current.caloric
                .specific_enthalpy_gradient.end()}};
}

[[nodiscard]] inline
mpmc::well_discretization::
    FixedBhpActivePhaseSourceInput3D
make_source_input(
    const mpmc::flow_discretization_petsc::
        FixedThreePhaseCurrentCellLinearization3D&
            current) {
    const auto mobility =
        mpmc::flow::
            build_local_phase_mobility_linearization(
                current.state,
                current.transport,
                current.saturation_constitutive);

    return {
        mobility.state_identity,
        std::vector<double>{
            mobility.phase_pressure_pa.begin(),
            mobility.phase_pressure_pa.end()},
        std::vector<std::vector<double>>{
            mobility.phase_pressure_gradient.begin(),
            mobility.phase_pressure_gradient.end()},
        std::vector<double>{
            mobility.mobility_per_pa_s.begin(),
            mobility.mobility_per_pa_s.end()},
        std::vector<std::vector<double>>{
            mobility.mobility_gradient.begin(),
            mobility.mobility_gradient.end()},
        std::vector<double>{
            current.molar_density
                .molar_density_mol_per_m3.begin(),
            current.molar_density
                .molar_density_mol_per_m3.end()},
        std::vector<std::vector<double>>{
            current.molar_density.gradient.begin(),
            current.molar_density.gradient.end()},
        std::vector<double>{
            current.transport
                .mass_density_kg_per_m3.begin(),
            current.transport
                .mass_density_kg_per_m3.end()},
        std::vector<std::vector<double>>{
            current.transport
                .mass_density_gradient.begin(),
            current.transport
                .mass_density_gradient.end()},
        std::vector<double>{
            current.caloric
                .specific_enthalpy_j_per_kg.begin(),
            current.caloric
                .specific_enthalpy_j_per_kg.end()},
        std::vector<std::vector<double>>{
            current.caloric
                .specific_enthalpy_gradient.begin(),
            current.caloric
                .specific_enthalpy_gradient.end()}};
}

} // namespace fixed_bhp_well_source_evaluator_detail

[[nodiscard]] inline
mpmc::well_discretization::
    FixedBhpConnectionCellSourceLinearization3D
build_fixed_bhp_peaceman_well_source_3d(
    const FixedBhpPeacemanWellSourceEvaluatorContext3D&
        context,
    const mpmc::flow_discretization_petsc::
        MixedCardinalityPhysicalCurrentCellLinearization3D&
            current) {
    auto source_input =
        std::visit(
            [](const auto& typed) {
                return
                    fixed_bhp_well_source_evaluator_detail::
                        make_source_input(typed);
            },
            current);

    const std::size_t phase_count =
        source_input.state_identity.layout
            .phase_count();
    if (context.injection_enthalpy()
            .specific_enthalpy_j_per_kg
            .size() !=
            phase_count ||
        (context.active_phase_identities()
             .has_value() &&
         context.active_phase_identities()
                 ->phase_count() !=
             phase_count)) {
        throw std::invalid_argument(
            "mpmc::well_discretization_petsc: fixed-BHP phase-identity/enthalpy binding does not match the frozen target chart");
    }

    return mpmc::well_discretization::
        make_fixed_bhp_connection_cell_source_3d(
            context.connection(),
            source_input,
            context.injection_enthalpy(),
            context.bottom_hole_pressure_pa(),
            context.source_provenance());
}

[[nodiscard]] inline
mpmc::well_discretization::
    FixedBhpConnectionCellSourceLinearization3D
build_fixed_bhp_three_phase_peaceman_well_source_3d(
    const FixedBhpThreePhasePeacemanWellSourceEvaluatorContext3D&
        context,
    const mpmc::flow_discretization_petsc::
        FixedThreePhaseCurrentCellLinearization3D&
            current) {
    return build_fixed_bhp_peaceman_well_source_3d(
        context,
        mpmc::flow_discretization_petsc::
            MixedCardinalityPhysicalCurrentCellLinearization3D{
                current});
}

/// Owner-only mixed-cardinality source callback.
///
/// Non-target owned cells publish no source. The target consumes exactly its
/// current frozen 1P/2P/3P chart and never pads inactive phases.
inline PetscErrorCode
evaluate_fixed_bhp_peaceman_well_source_3d(
    mpmc::mesh::LocalIndex,
    mpmc::mesh::GlobalEntityId cell_global,
    std::span<const double> natural_variables,
    const mpmc::flow_discretization_petsc::
        MixedCardinalityPhysicalCurrentCellLinearization3D&
            current,
    void* raw_context,
    std::optional<
        mpmc::flow_discretization::
            CellSourceLinearization3D>* output,
    mpmc::flow_discretization_petsc::
        NaturalVariableSnesEvaluationStatus3D*
            status) {
    using namespace
        mpmc::flow_discretization_petsc;

    if (raw_context == nullptr ||
        output == nullptr ||
        status == nullptr) {
        return PETSC_ERR_ARG_NULL;
    }
    output->reset();
    *status =
        NaturalVariableSnesEvaluationStatus3D::
            success;

    auto* context =
        static_cast<
            FixedBhpPeacemanWellSourceEvaluatorContext3D*>(
                raw_context);
    if (cell_global !=
        context->target_cell_global()) {
        return PETSC_SUCCESS;
    }

    const std::size_t q =
        std::visit(
            [](const auto& typed) {
                return typed.transport
                    .state_identity.layout
                    .unknown_count();
            },
            current);
    if (natural_variables.size() != q) {
        return PETSC_ERR_ARG_SIZ;
    }

    const std::size_t phase_count =
        std::visit(
            [](const auto& typed) {
                return typed.transport
                    .state_identity.layout
                    .phase_count();
            },
            current);
    if (context->injection_enthalpy()
            .specific_enthalpy_j_per_kg
            .size() !=
        phase_count) {
        return PETSC_ERR_SUP;
    }

    try {
        auto adapted =
            build_fixed_bhp_peaceman_well_source_3d(
                *context,
                current);
        if (adapted.cell_source.input_count !=
            q) {
            return PETSC_ERR_ARG_SIZ;
        }
        output->emplace(
            std::move(
                adapted.cell_source));
        return PETSC_SUCCESS;
    } catch (const std::invalid_argument&) {
        *status =
            NaturalVariableSnesEvaluationStatus3D::
                domain_error;
        return PETSC_SUCCESS;
    } catch (const std::range_error&) {
        *status =
            NaturalVariableSnesEvaluationStatus3D::
                domain_error;
        return PETSC_SUCCESS;
    } catch (...) {
        return PETSC_ERR_ARG_INCOMP;
    }
}

/// Owner-only multi-connection callback for one logical fixed-BHP well.
///
/// The surrounding mixed-cardinality assembly calls this evaluator only for
/// locally owned cells. This callback performs a deterministic stable-cell
/// lookup and delegates exactly one matching completion to the validated
/// single-connection evaluator. Non-completion owned cells publish no source.
inline PetscErrorCode
evaluate_fixed_bhp_multi_connection_well_source_3d(
    mpmc::mesh::LocalIndex cell,
    mpmc::mesh::GlobalEntityId cell_global,
    std::span<const double> natural_variables,
    const mpmc::flow_discretization_petsc::
        MixedCardinalityPhysicalCurrentCellLinearization3D&
            current,
    void* raw_context,
    std::optional<
        mpmc::flow_discretization::
            CellSourceLinearization3D>* output,
    mpmc::flow_discretization_petsc::
        NaturalVariableSnesEvaluationStatus3D*
            status) {
    if (raw_context == nullptr ||
        output == nullptr ||
        status == nullptr) {
        return PETSC_ERR_ARG_NULL;
    }
    output->reset();
    *status =
        mpmc::flow_discretization_petsc::
            NaturalVariableSnesEvaluationStatus3D::
                success;

    auto* context =
        static_cast<
            FixedBhpMultiConnectionWellSourceEvaluatorContext3D*>(
                raw_context);
    auto* connection =
        context->find_connection(
            cell_global);
    if (connection == nullptr) {
        return PETSC_SUCCESS;
    }

    return evaluate_fixed_bhp_peaceman_well_source_3d(
        cell,
        cell_global,
        natural_variables,
        current,
        connection,
        output,
        status);
}

/// Compatibility callback retaining the former strict-3P behavior.
inline PetscErrorCode
evaluate_fixed_bhp_three_phase_peaceman_well_source_3d(
    mpmc::mesh::LocalIndex cell,
    mpmc::mesh::GlobalEntityId cell_global,
    std::span<const double> natural_variables,
    const mpmc::flow_discretization_petsc::
        MixedCardinalityPhysicalCurrentCellLinearization3D&
            current,
    void* raw_context,
    std::optional<
        mpmc::flow_discretization::
            CellSourceLinearization3D>* output,
    mpmc::flow_discretization_petsc::
        NaturalVariableSnesEvaluationStatus3D*
            status) {
    if (raw_context == nullptr ||
        output == nullptr ||
        status == nullptr) {
        return evaluate_fixed_bhp_peaceman_well_source_3d(
            cell,
            cell_global,
            natural_variables,
            current,
            raw_context,
            output,
            status);
    }
    auto* context =
        static_cast<
            FixedBhpPeacemanWellSourceEvaluatorContext3D*>(
                raw_context);
    if (cell_global !=
        context->target_cell_global()) {
        return evaluate_fixed_bhp_peaceman_well_source_3d(
            cell,
            cell_global,
            natural_variables,
            current,
            raw_context,
            output,
            status);
    }
    if (!std::holds_alternative<
            mpmc::flow_discretization_petsc::
                FixedThreePhaseCurrentCellLinearization3D>(
                    current)) {
        output->reset();
        *status =
            mpmc::flow_discretization_petsc::
                NaturalVariableSnesEvaluationStatus3D::
                    success;
        return PETSC_ERR_SUP;
    }
    return evaluate_fixed_bhp_peaceman_well_source_3d(
        cell,
        cell_global,
        natural_variables,
        current,
        raw_context,
        output,
        status);
}

[[nodiscard]] inline
mpmc::flow_discretization_petsc::
    MixedCardinalityPhysicalCellSourceEvaluatorBinding3D
fixed_bhp_peaceman_well_source_binding_3d(
    FixedBhpPeacemanWellSourceEvaluatorContext3D*
        context) noexcept {
    return {
        &evaluate_fixed_bhp_peaceman_well_source_3d,
        context};
}

[[nodiscard]] inline
mpmc::flow_discretization_petsc::
    MixedCardinalityPhysicalCellSourceEvaluatorBinding3D
fixed_bhp_multi_connection_well_source_binding_3d(
    FixedBhpMultiConnectionWellSourceEvaluatorContext3D*
        context) noexcept {
    return {
        &evaluate_fixed_bhp_multi_connection_well_source_3d,
        context};
}

[[nodiscard]] inline
mpmc::flow_discretization_petsc::
    MixedCardinalityPhysicalCellSourceEvaluatorBinding3D
fixed_bhp_three_phase_peaceman_well_source_binding_3d(
    FixedBhpThreePhasePeacemanWellSourceEvaluatorContext3D*
        context) noexcept {
    return {
        &evaluate_fixed_bhp_three_phase_peaceman_well_source_3d,
        context};
}

} // namespace mpmc::well_discretization_petsc

#endif // MPMC_WELL_DISCRETIZATION_PETSC_FIXED_BHP_WELL_SOURCE_EVALUATOR_HPP
