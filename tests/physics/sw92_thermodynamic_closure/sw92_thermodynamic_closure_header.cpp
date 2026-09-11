#include <mpmc/physics/sw92_thermodynamic_closure.hpp>

#include <string_view>

bool sw92_thermodynamic_closure_header() {
    mpmc::physics::PtPhaseSetThermodynamicClosureSnapshot snapshot;
    mpmc::physics::PtPhaseSetThermodynamicLinearization linearization;
    linearization.component_count = 2U;
    linearization.phase_count = 1U;
    linearization.input_count = 3U;
    linearization.phase_fraction_jacobian.assign(3U, 0.0);
    linearization.composition_jacobian.assign(5U, 0.0); // Intentionally malformed: expected 6.
    linearization.compressibility_jacobian.assign(3U, 0.0);
    linearization.molar_density_jacobian.assign(3U, 0.0);
    linearization.equilibrium_jacobian_rcond = 1.0;
    linearization.linear_solve_backward_error = 0.0;
    linearization.equilibrium_residual_norm = 0.0;

    snapshot.component_ids = {"A", "B"};
    snapshot.primal_status = mpmc::physics::ThermodynamicClosurePrimalStatus::valid;
    snapshot.linearization_status =
        mpmc::physics::ThermodynamicClosureLinearizationStatus::available;
    snapshot.linearization_reason =
        mpmc::physics::ThermodynamicClosureLinearizationReason::none;
    snapshot.primal.emplace();
    snapshot.primal->phases.resize(1U);
    snapshot.linearization = linearization;

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
