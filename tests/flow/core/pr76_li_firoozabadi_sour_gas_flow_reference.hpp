// Generated only by reference_li_firoozabadi_sour_gas_flow.py.
// Independent of MPMC production C++ output.
#ifndef MPMC_TEST_PR76_LI_FIROOZABADI_SOUR_GAS_FLOW_REFERENCE_HPP
#define MPMC_TEST_PR76_LI_FIROOZABADI_SOUR_GAS_FLOW_REFERENCE_HPP

#include <array>

namespace pr76_li_firoozabadi_sour_gas_flow_reference {
inline constexpr double pressure_pa = 2000000.0;
inline constexpr double temperature_k = 178.8;
inline constexpr double porosity = 0.3;
inline constexpr double timestep_seconds = 1.0;
inline constexpr double bottom_hole_pressure_pa = pressure_pa;
inline constexpr double dx_m = 304.8;
inline constexpr double dy_m = 304.8;
inline constexpr double dz_m = 6.096;
inline constexpr double bulk_volume_m3 = 566336.93184;
inline constexpr double permeability_m2 = 4.9346165e-13;
inline constexpr double wellbore_radius_m = 0.0762;
inline constexpr std::array<double, 3> phase_saturation{
    0.69062384541647425,
    0.18515603298995780,
    0.12422012159356794};
inline constexpr std::array<double, 3> molar_density_mol_per_m3{
    1549.6059815988284,
    23296.967613285993,
    28407.099160395410};
inline constexpr std::array<double, 3> mass_density_kg_per_m3{
    39.451620206140094,
    812.67677379555095,
    1152.1242703461601};
inline constexpr std::array<double, 3> dynamic_viscosity_pa_s{
    0.000010951380704796704,
    0.00028297191613066484,
    0.00036514189669459766};
inline constexpr std::array<double, 3> specific_enthalpy_j_per_kg{
    -169779.39041907499,
    -518174.30519072997,
    -477763.37551774654};
inline constexpr std::array<double, 3> specific_internal_energy_j_per_kg{
    -220474.39329981164,
    -520635.30821738120,
    -479499.29936856438};
inline constexpr std::array<double, 3> mixture_molar_mass_kg_per_mol{
    0.025459130046358826,
    0.034883371402039926,
    0.040557617792682890};
inline constexpr std::array<double, 3> ideal_gas_molar_enthalpy_j_per_mol{
    -3731.2336342079718,
    -4796.5558190130165,
    -4239.7987244640812};
inline constexpr std::array<double, 3> residual_molar_enthalpy_j_per_mol{
    -591.20194566278620,
    -13279.110919949202,
    -15137.145655126713};
inline constexpr std::array<double, 6> component_accumulation_mol_per_bulk_m3{
    1336.8753386532859,
    319.72517289733769,
    89.464800728176187,
    312.17117649811223,
    480.49787939410598,
    135.01630913555379};
inline constexpr double total_internal_energy_j_per_bulk_m3 = -45891793.599850742;
inline constexpr double quartz_specific_heat_j_per_kg_k = 494.81656050955414;
inline constexpr double quartz_volumetric_heat_capacity_j_per_m3_k = 1305227.1233121019;
inline constexpr std::array<double, 6> cell_component_inventory_mol{
    757121877.54546288,
    181072173.45069175,
    50667220.752072300,
    176794066.30682400,
    272123694.77168434,
    76464722.264190495};
inline constexpr double cell_total_internal_energy_j = -25990217583974.018;
inline constexpr std::array<double, 3> well_phase_volumetric_rate_m3_per_s{
    0.0, 0.0, 0.0};
inline constexpr std::array<double, 6> well_component_molar_rate_mol_per_s{
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
inline constexpr double well_energy_rate_w = 0.0;
} // namespace pr76_li_firoozabadi_sour_gas_flow_reference

#endif
