#ifndef MPMC_WELL_DISCRETIZATION_PETSC_FIXED_BHP_WELL_SOURCE_EVALUATOR_HPP
#define MPMC_WELL_DISCRETIZATION_PETSC_FIXED_BHP_WELL_SOURCE_EVALUATOR_HPP

#include <mpmc/flow_discretization_petsc/mixed_cardinality_physical_snes_assembly.hpp>
#include <mpmc/well_discretization/energy_rate.hpp>
#include <mpmc/well_discretization/fixed_bhp_connection_source.hpp>

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
        "well-discretization-petsc/fixed-bhp-single-connection/variable-cardinality-peaceman/v2";

/// Backward-compatible name retained for callers of the first 3P-only bridge.
inline constexpr std::string_view
    fixed_bhp_three_phase_peaceman_well_source_evaluator_convention =
        fixed_bhp_peaceman_well_source_evaluator_convention;

/// Immutable configuration for one fixed-BHP Peaceman connection.
///
/// The injection enthalpy vector is cardinality-specific for the currently
/// frozen target chart and therefore has exactly 1, 2 or 3 entries. This
/// evaluator does not remap completion phase identity across a topology change.
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

        auto normalized_connection =
            mpmc::well::
                make_peaceman_well_index_3d(
                    connection.cell_dimensions_m,
                    connection.permeability_m2,
                    connection.well_direction,
                    connection.wellbore_radius_m,
                    connection.skin_factor);

        return FixedBhpPeacemanWellSourceEvaluatorContext3D{
            target_cell_global,
            std::move(normalized_connection),
            bottom_hole_pressure_pa,
            std::move(injection_enthalpy),
            std::move(source_provenance)};
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

    [[nodiscard]] const std::string&
    source_provenance() const noexcept {
        return source_provenance_;
    }

private:
    FixedBhpPeacemanWellSourceEvaluatorContext3D(
        mpmc::mesh::GlobalEntityId target_cell_global,
        mpmc::well::PeacemanWellIndex3D connection,
        double bottom_hole_pressure_pa,
        mpmc::well_discretization::
            FixedBhpInjectionEnthalpy3D
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

    mpmc::mesh::GlobalEntityId
        target_cell_global_;
    mpmc::well::PeacemanWellIndex3D
        connection_;
    double bottom_hole_pressure_pa_{};
    mpmc::well_discretization::
        FixedBhpInjectionEnthalpy3D
            injection_enthalpy_;
    std::string source_provenance_;
};

using FixedBhpThreePhasePeacemanWellSourceEvaluatorContext3D =
    FixedBhpPeacemanWellSourceEvaluatorContext3D;

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

    if (context.injection_enthalpy()
            .specific_enthalpy_j_per_kg
            .size() !=
        source_input.state_identity.layout
            .phase_count()) {
        throw std::invalid_argument(
            "mpmc::well_discretization_petsc: fixed-BHP injection enthalpy cardinality does not match the frozen target chart");
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
    if (cell_global !=
            static_cast<
                FixedBhpPeacemanWellSourceEvaluatorContext3D*>(
                    raw_context)
                ->target_cell_global()) {
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
        if (output != nullptr) {
            output->reset();
        }
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
fixed_bhp_three_phase_peaceman_well_source_binding_3d(
    FixedBhpThreePhasePeacemanWellSourceEvaluatorContext3D*
        context) noexcept {
    return {
        &evaluate_fixed_bhp_three_phase_peaceman_well_source_3d,
        context};
}

} // namespace mpmc::well_discretization_petsc

#endif // MPMC_WELL_DISCRETIZATION_PETSC_FIXED_BHP_WELL_SOURCE_EVALUATOR_HPP
