#include <mpmc/well_discretization_petsc/fixed_total_molar_rate_timestep_driver.hpp>

#include <string_view>

bool fixed_total_molar_rate_timestep_driver_header() {
    return
        mpmc::well_discretization_petsc::
            fixed_total_molar_rate_timestep_driver_convention ==
        std::string_view{
            "well-discretization-petsc/single-well/fixed-total-molar-rate-physical-timestep-driver/v1"};
}
