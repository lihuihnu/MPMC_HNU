#include <mpmc/model_configuration/pr76_model_registry.hpp>
#include <type_traits>
#include <utility>

namespace mc = mpmc::model_configuration;
static_assert(!std::is_copy_constructible_v<mc::Pr76ModelRegistry>);
static_assert(!std::is_move_constructible_v<mc::Pr76ModelRegistry>);
static_assert(!std::is_copy_constructible_v<mc::Pr76ModelSolveLease>);
static_assert(std::is_nothrow_move_constructible_v<mc::Pr76ModelSolveLease>);
static_assert(std::is_nothrow_move_assignable_v<mc::Pr76ModelSolveLease>);
static_assert(std::is_same_v<decltype(std::declval<const mc::Pr76ModelRegistry&>().describe("")),
                             mc::Pr76RegisteredModelSnapshot>);
template <class T>
concept LvalueSolve = requires(T& model) { model.solve(mpmc::flash::PtFlashRequest{}); };
static_assert(!LvalueSolve<mc::Pr76ModelSolveLease>);
static_assert(std::is_same_v<decltype(std::declval<mc::Pr76ModelSolveLease&&>().solve(
    mpmc::flash::PtFlashRequest{})), mpmc::flash::PtFlashBackendResult>);
