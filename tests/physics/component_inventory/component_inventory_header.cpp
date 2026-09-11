#include <mpmc/physics/component_inventory.hpp>

#include <string_view>

bool component_inventory_header() {
    mpmc::physics::PtComponentInventorySnapshot snapshot;
    mpmc::physics::PtComponentInventoryLinearization linearization;
    linearization.component_count = 2U;
    linearization.input_count = 3U;
    linearization.total_molar_density_gradient.assign(3U, 0.0);
    linearization.component_molar_density_jacobian.assign(6U, 0.0);
    return mpmc::physics::PtComponentInventorySnapshot::convention ==
               std::string_view{
                   "PT/component-inventory/phase-resolved-reduced-feed/v1"} &&
           !snapshot.inventory_available() &&
           !snapshot.linearization_available() &&
           linearization.dependent_feed_component() == 1U &&
           linearization.feed_column(0U) == 2U;
}
