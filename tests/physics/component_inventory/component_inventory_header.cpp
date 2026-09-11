#include <mpmc/physics/component_inventory.hpp>

#include <limits>
#include <string_view>

bool component_inventory_header() {
    mpmc::physics::PtComponentInventorySnapshot snapshot;
    mpmc::physics::PtComponentInventoryLinearization linearization;
    linearization.component_count = 2U;
    linearization.input_count = 3U;
    linearization.total_molar_density_gradient.assign(3U, 0.0);
    linearization.component_molar_density_jacobian.assign(6U, 0.0);

    // A public fixed-VLE snapshot can be externally malformed. In particular,
    // the legacy availability helper does not validate derivative dimensions.
    // The inventory adapter must reject input_count==0 without dividing by zero
    // while retaining the independently valid primal inventory.
    mpmc::physics::ThermodynamicClosureSnapshot malformed_source;
    malformed_source.primal_status =
        mpmc::physics::ThermodynamicClosurePrimalStatus::valid;
    malformed_source.linearization_status =
        mpmc::physics::ThermodynamicClosureLinearizationStatus::available;
    malformed_source.linearization_reason =
        mpmc::physics::ThermodynamicClosureLinearizationReason::none;
    malformed_source.pressure_pa = 1.0e6;
    malformed_source.temperature_k = 300.0;
    malformed_source.feed = {0.4, 0.6};
    malformed_source.component_ids = {"A", "B"};
    malformed_source.primal.emplace();
    malformed_source.primal->liquid = {
        0.5, {0.4, 0.6}, 1.0, 1000.0};
    malformed_source.primal->vapor = {
        0.5, {0.4, 0.6}, 1.0, 2000.0};
    malformed_source.linearization.emplace();
    malformed_source.linearization->component_count = 2U;
    malformed_source.linearization->input_count = 0U;

    const auto guarded =
        mpmc::physics::build_pt_component_inventory(malformed_source);
    if (!guarded.inventory_available() || guarded.linearization_available() ||
        guarded.linearization_reason !=
            mpmc::physics::ThermodynamicClosureLinearizationReason::solution_not_accepted) {
        return false;
    }

    auto nonfinite_primal = guarded;
    nonfinite_primal.primal->component_molar_density_mol_per_m3[0] =
        std::numeric_limits<double>::quiet_NaN();

    return mpmc::physics::PtComponentInventorySnapshot::convention ==
               std::string_view{
                   "PT/component-inventory/phase-resolved-reduced-feed/v1"} &&
           !snapshot.inventory_available() &&
           !snapshot.linearization_available() &&
           !nonfinite_primal.inventory_available() &&
           linearization.dependent_feed_component() == 1U &&
           linearization.feed_column(0U) == 2U;
}
