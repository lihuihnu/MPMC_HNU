#include <mpmc/flow/pr76_li_firoozabadi_sour_gas_properties.hpp>

#include <type_traits>

static_assert(std::is_class_v<
    mpmc::flow::Pr76LiFiroozabadiSourGasPropertyProvider>);

bool pr76_li_firoozabadi_sour_gas_properties_header() {
    return true;
}
