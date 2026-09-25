#ifndef MPMC_TEST_FLOW_CORE_SPE3_KENYON_BEHIE_SATURATION_REFERENCE_HPP
#define MPMC_TEST_FLOW_CORE_SPE3_KENYON_BEHIE_SATURATION_REFERENCE_HPP

#include <string_view>

namespace mpmc::flow::test_reference::spe3 {

/// Physical reference used only by the flow saturation-constitutive regression.
///
/// Primary publication:
/// D. E. Kenyon and G. A. Behie,
/// "Third SPE Comparative Solution Project: Gas Cycling of Retrograde
/// Condensate Reservoirs", Journal of Petroleum Technology 39(8), 981-997,
/// 1987, DOI 10.2118/12278-PA.
///
/// Public machine-readable source used to transcribe the minimal brackets below:
/// OPM/opm-data, spe3/SPE3CASE2.DATA at
/// b65b5209d9135e66592fce07d00bf6f8d6d1a651.  The deck states that SGOF/SWOF
/// values are taken from Kenyon & Behie's table 2 and is distributed under
/// ODbL-1.0 / DBCL-1.0.  Only the six rows required to bracket the tested
/// state are reproduced here.
inline constexpr std::string_view publication_doi =
    "10.2118/12278-PA";
inline constexpr std::string_view source_revision =
    "OPM/opm-data@b65b5209d9135e66592fce07d00bf6f8d6d1a651:"
    "spe3/SPE3CASE2.DATA";
inline constexpr std::string_view source_license =
    "ODbL-1.0 / DBCL-1.0";

inline constexpr double psi_to_pa =
    6894.757293168;

/// SWOF source rows: Sw, Krw, Krow, Pcow[psi], with Pcow = Po - Pw.
struct WaterOilRow {
    double water_saturation;
    double water_relative_permeability;
    double oil_relative_permeability;
    double pcow_psi;
};

inline constexpr WaterOilRow water_lower{
    0.28, 0.020, 0.400, 15.5};
inline constexpr WaterOilRow water_upper{
    0.32, 0.033, 0.315, 12.0};

/// These two SWOF rows bracket So=0.55 after So=1-Sw.
inline constexpr WaterOilRow oil_low_saturation_source{
    0.48, 0.119, 0.112, 4.2}; // So=0.52
inline constexpr WaterOilRow oil_high_saturation_source{
    0.44, 0.090, 0.150, 5.3}; // So=0.56

/// SGOF source rows: Sg, Krg, Krog, Pcgo[psi], with Pcgo = Pg - Po.
struct GasOilRow {
    double gas_saturation;
    double gas_relative_permeability;
    double oil_relative_permeability;
    double pcgo_psi;
};

inline constexpr GasOilRow gas_lower{
    0.12, 0.026, 0.400, 0.0};
inline constexpr GasOilRow gas_upper{
    0.16, 0.040, 0.315, 0.0};

/// Interior three-phase point.  Physical role mapping is explicit:
/// phase0=oil, phase1=water, phase2=gas.
inline constexpr double oil_saturation = 0.55;
inline constexpr double water_saturation = 0.30;
inline constexpr double gas_saturation = 0.15;

/// Independent hand-calculated piecewise-linear oracle at the point above.
inline constexpr double oil_relative_permeability = 0.1405;
inline constexpr double water_relative_permeability = 0.0265;
inline constexpr double gas_relative_permeability = 0.0365;
inline constexpr double pcow_psi = 13.75;
inline constexpr double water_minus_oil_pressure_pa =
    -94802.91278106;
inline constexpr double gas_minus_oil_pressure_pa = 0.0;

/// Exact slopes of the source-table segments used by the oracle.
inline constexpr double d_kro_d_so = 0.95;
inline constexpr double d_krw_d_sw = 0.325;
inline constexpr double d_krg_d_sg = 0.35;
inline constexpr double d_water_minus_oil_pressure_d_sw_pa =
    603291.2631522;

} // namespace mpmc::flow::test_reference::spe3

#endif
