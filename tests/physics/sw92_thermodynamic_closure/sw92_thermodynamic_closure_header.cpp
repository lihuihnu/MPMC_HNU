#include <mpmc/physics/sw92_thermodynamic_closure.hpp>

#include <string_view>

bool sw92_thermodynamic_closure_header() {
    mpmc::physics::PtPhaseSetThermodynamicClosureSnapshot snapshot;
    mpmc::physics::PtPhaseSetThermodynamicLinearization linearization;
    linearization.component_count = 2U;
    linearization.phase_count = 1U;
    linearization.input_count = 3U;
    return mpmc::physics::PtPhaseSetThermodynamicClosureSnapshot::convention ==
               std::string_view{
                   "PT/phase-set/thermodynamic-closure/reduced-feed-v2"} &&
           mpmc::physics::sw92_profile_c_closure_convention ==
               std::string_view{
                   "SW92/phase-assigned/Profile-C/physics-thermodynamic-closure-reduced-feed/v2"} &&
           !snapshot.can_seed_newton() &&
           linearization.dependent_feed_component() == 1U &&
           linearization.feed_column(0U) == 2U;
}
