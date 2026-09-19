#include <mpmc/flow/fugacity_equilibrium_residual.hpp>

#include <cstddef>
#include <span>
#include <type_traits>
#include <utility>
#include <vector>

namespace {

struct ProbeEvaluator {
    template <typename Number>
    mpmc::flow::
        PhaseLnFugacityCoefficientEvaluation<
            Number>
    operator()(
        mpmc::flow::PhaseSlot3,
        const Number&,
        const Number&,
        std::span<const Number>
            composition) const {
        return {
            std::vector<Number>(
                composition.size(),
                Number{})};
    }
};

static_assert(
    mpmc::flow::
        PhaseLnFugacityCoefficientEvaluator3P<
            ProbeEvaluator,
            double>);

static_assert(
    std::is_same_v<
        decltype(
            std::declval<
                const mpmc::flow::
                    FugacityEquilibriumResidual3P<
                        double>&>()
                .values()),
        std::span<const double>>);

} // namespace

bool fugacity_equilibrium_residual_header() {
    return true;
}
