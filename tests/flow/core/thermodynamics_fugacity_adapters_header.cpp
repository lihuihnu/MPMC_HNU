#include <mpmc/flow/thermodynamics_fugacity_adapters.hpp>

#include <type_traits>

static_assert(std::is_class_v<
    mpmc::flow::CpaSelectedPhaseFugacityEvaluator3P>);
static_assert(std::is_class_v<
    mpmc::flow::Pr76SelectedPhaseFugacityEvaluator3P<double>>);
static_assert(std::is_class_v<
    mpmc::flow::Sw92SelectedPhaseFugacityEvaluator3P<double>>);

bool thermodynamics_fugacity_adapters_header() {
    return true;
}
