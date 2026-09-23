#ifndef MPMC_WELL_DISCRETIZATION_PETSC_FIXED_BHP_WELL_SOURCE_EVALUATOR_HPP
#define MPMC_WELL_DISCRETIZATION_PETSC_FIXED_BHP_WELL_SOURCE_EVALUATOR_HPP

#include <mpmc/flow_discretization_petsc/mixed_cardinality_physical_snes_assembly.hpp>
#include <mpmc/well_discretization/cell_source_adapter.hpp>

#include <cmath>
#include <optional>
#include <span>
#include <stdexcept>
#include <string>
#include <string_view>
#include <utility>
#include <variant>

namespace mpmc::well_discretization_petsc {

inline constexpr std::string_view
    fixed_bhp_three_phase_peaceman_well_source_evaluator_convention =
        "well-discretization-petsc/fixed-bhp-single-connection/three-phase-peaceman/v1";

/// Immutable configuration for the first production well/source bridge.
///
/// This slice supports exactly one stable target cell and a frozen BHP. The
/// mixed-cardinality reservoir may contain 1P/2P/3P cells, but the configured
/// connection is intentionally valid only while its target cell is frozen 3P.
/// Phase-transition-aware completion remapping is a separate future contract.
class FixedBhpThreePhasePeacemanWellSourceEvaluatorContext3D {
public:
    static constexpr std::string_view convention =
        fixed_bhp_three_phase_peaceman_well_source_evaluator_convention;

    [[nodiscard]] static
    FixedBhpThreePhasePeacemanWellSourceEvaluatorContext3D
    create(
        mpmc::mesh::GlobalEntityId target_cell_global,
        const mpmc::well::PeacemanWellIndex3D&
            connection,
        double bottom_hole_pressure_pa,
        mpmc::well_discretization::
            InjectionPhaseSpecificEnthalpy3P
                injection_enthalpy,
        std::string source_provenance) {
        if (!std::isfinite(bottom_hole_pressure_pa) ||
            !(bottom_hole_pressure_pa > 0.0) ||
            source_provenance.empty() ||
            injection_enthalpy.provenance.empty()) {
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

        // Rebuild the Peaceman connection from its primitive inputs so derived
        // fields cannot be corrupted before entering the production bridge.
        auto normalized_connection =
            mpmc::well::
                make_peaceman_well_index_3d(
                    connection.cell_dimensions_m,
                    connection.permeability_m2,
                    connection.well_direction,
                    connection.wellbore_radius_m,
                    connection.skin_factor);

        return FixedBhpThreePhasePeacemanWellSourceEvaluatorContext3D{
            target_cell_global,
            std::move(normalized_connection),
            bottom_hole_pressure_pa,
            std::move(injection_enthalpy),
            std::move(source_provenance)};
    }

    [[nodiscard]] mpmc::mesh::GlobalEntityId
    target_cell_global() const noexcept {
        return target_cell_global_;
    }

    [[nodiscard]] const mpmc::well::PeacemanWellIndex3D&
    connection() const noexcept {
        return connection_;
    }

    [[nodiscard]] double
    bottom_hole_pressure_pa() const noexcept {
        return bottom_hole_pressure_pa_;
    }

    [[nodiscard]] const mpmc::well_discretization::
        InjectionPhaseSpecificEnthalpy3P&
    injection_enthalpy() const noexcept {
        return injection_enthalpy_;
    }

    [[nodiscard]] const std::string&
    source_provenance() const noexcept {
        return source_provenance_;
    }

private:
    FixedBhpThreePhasePeacemanWellSourceEvaluatorContext3D(
        mpmc::mesh::GlobalEntityId target_cell_global,
        mpmc::well::PeacemanWellIndex3D connection,
        double bottom_hole_pressure_pa,
        mpmc::well_discretization::
            InjectionPhaseSpecificEnthalpy3P
                injection_enthalpy,
        std::string source_provenance)
        : target_cell_global_(
              target_cell_global),
          connection_(
              std::move(connection)),
          bottom_hole_pressure_pa_(
              bottom_hole_pressure_pa),
          injection_enthalpy_(
              std::move(injection_enthalpy)),
          source_provenance_(
              std::move(source_provenance)) {}

    mpmc::mesh::GlobalEntityId target_cell_global_;
    mpmc::well::PeacemanWellIndex3D connection_;
    double bottom_hole_pressure_pa_{};
    mpmc::well_discretization::
        InjectionPhaseSpecificEnthalpy3P
            injection_enthalpy_;
    std::string source_provenance_;
};

/// Evaluate the already-validated well physics chain on one frozen 3P current
/// cell and retain the BHP derivative sidecar for future well-unknown work.
[[nodiscard]] inline
mpmc::well_discretization::
    WellConnectionCellSourceAdapterResult3P
build_fixed_bhp_three_phase_peaceman_well_source_3d(
    const FixedBhpThreePhasePeacemanWellSourceEvaluatorContext3D&
        context,
    const mpmc::flow_discretization_petsc::
        FixedThreePhaseCurrentCellLinearization3D&
            current) {
    const auto mobility =
        mpmc::flow::
            build_local_phase_mobility_linearization(
                current.state,
                current.transport,
                current.saturation_constitutive);

    const auto component_rate =
        mpmc::well_discretization::
            make_connection_component_molar_rate_linearization_3p(
                context.connection(),
                mobility,
                current.molar_density,
                context.bottom_hole_pressure_pa());

    const auto energy_rate =
        mpmc::well_discretization::
            make_connection_advective_energy_rate_linearization_3p(
                context.connection(),
                mobility,
                current.transport,
                current.caloric,
                context.injection_enthalpy(),
                context.bottom_hole_pressure_pa());

    return mpmc::well_discretization::
        make_connection_cell_source_adapter_3p(
            component_rate,
            energy_rate,
            context.source_provenance());
}

/// Mixed-cardinality source callback.
///
/// The surrounding MixedCardinalityPhysicalSnesAssemblyContext3D already calls
/// this callback only for locally owned cells. Non-target owned cells publish
/// no source. The configured target must be frozen 3P; 1P/2P target states are
/// explicitly unsupported rather than padded with fictitious inactive phases.
inline PetscErrorCode
evaluate_fixed_bhp_three_phase_peaceman_well_source_3d(
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
    using namespace mpmc::flow_discretization_petsc;

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
            FixedBhpThreePhasePeacemanWellSourceEvaluatorContext3D*>(
                raw_context);
    if (cell_global !=
        context->target_cell_global()) {
        return PETSC_SUCCESS;
    }

    const auto* three_phase =
        std::get_if<
            FixedThreePhaseCurrentCellLinearization3D>(
                &current);
    if (three_phase == nullptr) {
        return PETSC_ERR_SUP;
    }

    const std::size_t q =
        three_phase->transport.state_identity.layout
            .unknown_count();
    if (natural_variables.size() != q) {
        return PETSC_ERR_ARG_SIZ;
    }

    try {
        auto adapted =
            build_fixed_bhp_three_phase_peaceman_well_source_3d(
                *context,
                *three_phase);
        if (adapted.cell_source.input_count != q) {
            return PETSC_ERR_ARG_SIZ;
        }
        output->emplace(
            std::move(adapted.cell_source));
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
