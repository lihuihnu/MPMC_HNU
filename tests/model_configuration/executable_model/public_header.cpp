// Must be the first include: verify the transport-neutral hint DTO is self-contained.
#include <mpmc/model_configuration/pt_solve_hints.hpp>
#include <mpmc/model_configuration/pr76_executable_model.hpp>
#include <type_traits>
#include <utility>

namespace mc = mpmc::model_configuration;
static_assert(!std::is_copy_constructible_v<mc::Pr76ExecutableModel>);
static_assert(!std::is_move_constructible_v<mc::Pr76ExecutableModel>);
static_assert(!std::is_copy_assignable_v<mc::Pr76ExecutableModel>);
static_assert(!std::is_move_assignable_v<mc::Pr76ExecutableModel>);
static_assert(std::is_base_of_v<mpmc::flash::PtFlashBackend, mc::Pr76ExecutableModel>);
static_assert(std::is_copy_constructible_v<mc::PtSolveHints>);
static_assert(std::is_move_constructible_v<mc::PtSolveHints>);
static_assert(std::is_same_v<decltype(std::declval<mc::Pr76ExecutableModel&>().parameter_snapshot()),
                             const mc::Pr76ModelParameters&>);
static_assert(std::is_same_v<decltype(std::declval<mc::Pr76ExecutableModel&>().solver_configuration()),
                             const mc::Pr76SolverConfiguration&>);
static_assert(std::is_same_v<decltype(std::declval<mc::Pr76ExecutableModel&>().solve(
                                 std::declval<const mpmc::flash::PtFlashRequest&>(),
                                 std::declval<const mc::PtSolveHints&>())),
                             mpmc::flash::PtFlashBackendResult>);
static_assert(std::is_same_v<decltype(mc::make_pr76_executable_model(
    std::declval<const mc::ThermodynamicModelDefinition&>(), std::declval<const mc::PtSolverSettings&>())),
    std::unique_ptr<mc::Pr76ExecutableModel>>);
