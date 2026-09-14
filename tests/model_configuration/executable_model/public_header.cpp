// Must be the first include: verify the new public header is self-contained.
#include <mpmc/model_configuration/pr76_executable_model.hpp>
#include <type_traits>
#include <utility>

namespace mc = mpmc::model_configuration;
static_assert(!std::is_copy_constructible_v<mc::Pr76ExecutableModel>);
static_assert(!std::is_move_constructible_v<mc::Pr76ExecutableModel>);
static_assert(!std::is_copy_assignable_v<mc::Pr76ExecutableModel>);
static_assert(!std::is_move_assignable_v<mc::Pr76ExecutableModel>);
static_assert(std::is_base_of_v<mpmc::flash::PtFlashBackend, mc::Pr76ExecutableModel>);
static_assert(std::is_same_v<decltype(std::declval<mc::Pr76ExecutableModel&>().parameter_snapshot()),
                             const mc::Pr76ModelParameters&>);
static_assert(std::is_same_v<decltype(std::declval<mc::Pr76ExecutableModel&>().solver_configuration()),
                             const mc::Pr76SolverConfiguration&>);
static_assert(std::is_same_v<decltype(mc::make_pr76_executable_model(
    std::declval<const mc::ThermodynamicModelDefinition&>(), std::declval<const mc::PtSolverSettings&>())),
    std::unique_ptr<mc::Pr76ExecutableModel>>);
