#include <mpmc/well_discretization_petsc/fixed_bhp_well_source_evaluator.hpp>

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
                    evaluate_fixed_bhp_three_phase_peaceman_well_source_3d),
            Callback>);
    static_assert(
        mpmc::well_discretization_petsc::
            FixedBhpThreePhasePeacemanWellSourceEvaluatorContext3D::
                convention ==
        mpmc::well_discretization_petsc::
            fixed_bhp_three_phase_peaceman_well_source_evaluator_convention);
    return
        mpmc::well_discretization_petsc::
            fixed_bhp_three_phase_peaceman_well_source_evaluator_convention ==
        std::string_view{
            "well-discretization-petsc/fixed-bhp-single-connection/three-phase-peaceman/v1"};
}
