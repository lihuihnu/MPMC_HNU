#include <mpmc/flow/pr76_selected_phase_property_closure.hpp>

#include <type_traits>

static_assert(std::is_class_v<
    mpmc::flow::Pr76SelectedPhasePropertyChartLinearization>);
static_assert(std::is_class_v<
    mpmc::flow::Pr76SelectedPhaseFlowLinearization3P>);

bool pr76_selected_phase_property_closure_header() {
    return true;
}
