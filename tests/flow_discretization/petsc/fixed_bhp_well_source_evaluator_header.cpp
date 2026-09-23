#include <mpmc/well_discretization_petsc/fixed_bhp_well_source_evaluator.hpp>
#include <mpmc/well_discretization_petsc/fixed_total_molar_rate_control.hpp>

#include <string_view>
#include <type_traits>

bool fixed_bhp_well_source_evaluator_header() {
    using Callback =
        mpmc::flow_discretization_petsc::
            MixedCardinalityPhysicalCellSourceEvaluator3D;
    static_assert(
        std::is_same_v<
            decltype(
                &mpmc::well_discretization_petsc::
                    evaluate_fixed_bhp_peaceman_well_source_3d),
            Callback>);
    static_assert(
        std::is_same_v<
            decltype(
                &mpmc::well_discretization_petsc::
                    evaluate_fixed_bhp_multi_connection_well_source_3d),
            Callback>);
    static_assert(
        mpmc::well_discretization_petsc::
            FixedBhpPeacemanWellSourceEvaluatorContext3D::
                convention ==
        mpmc::well_discretization_petsc::
            fixed_bhp_peaceman_well_source_evaluator_convention);
    static_assert(
        std::is_same_v<
            mpmc::well_discretization_petsc::
                FixedBhpThreePhasePeacemanWellSourceEvaluatorContext3D,
            mpmc::well_discretization_petsc::
                FixedBhpPeacemanWellSourceEvaluatorContext3D>);
    return
        mpmc::well_discretization_petsc::
            fixed_bhp_peaceman_well_source_evaluator_convention ==
            std::string_view{
                "well-discretization-petsc/fixed-bhp-single-connection/phase-identity-rebindable-peaceman/v3"} &&
        mpmc::well_discretization_petsc::
            fixed_bhp_multi_connection_well_source_evaluator_convention ==
            std::string_view{
                "well-discretization-petsc/fixed-bhp-single-well/multi-connection-owner-only/v1"} &&
        mpmc::well_discretization_petsc::
            FixedTotalMolarRateWellControlSystem3D::
                convention ==
            std::string_view{
                "well-discretization-petsc/single-well/fixed-total-molar-rate-control/v1"};
}
