#include <mpmc/flow/sw92_co2_water_properties.hpp>

bool sw92_co2_water_properties_header() {
    return mpmc::flow::
               sw92_co2_water_property_convention ==
           "flow/sw92/co2-water/zero-salinity/Chung1988-NISTCO2-IAPWS95H2O-SW92-departure/v2";
}
