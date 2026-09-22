#include <mpmc/flow/pr76_methane_ethane_propane_properties.hpp>

#include <type_traits>

static_assert(std::is_class_v<
    mpmc::flow::Pr76MethaneEthanePropanePropertyProvider>);

bool pr76_methane_ethane_propane_properties_header() {
    return
        mpmc::flow::
            pr76_methane_ethane_propane_property_convention ==
        "flow/pr76/methane-ethane-propane/NIST-SRD69-LBC/v1";
}
