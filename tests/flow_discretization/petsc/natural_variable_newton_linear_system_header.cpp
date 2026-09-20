#include <mpmc/flow_discretization_petsc/natural_variable_newton_linear_system.hpp>

#include <string_view>
#include <type_traits>
#include <utility>

static_assert(
    mpmc::flow_discretization_petsc::
        FrozenNaturalVariableNewtonLinearCorrection3D::
            convention ==
    mpmc::flow_discretization_petsc::
        frozen_natural_variable_newton_linear_system_convention);

static_assert(
    std::is_same_v<
        decltype(
            (std::declval<
                const mpmc::flow_discretization_petsc::
                    NaturalVariableLinearCorrectionEntry3D&>()
                 .petsc_global_scalar)),
        const PetscInt&>);

bool natural_variable_newton_linear_system_header() {
    return true;
}
