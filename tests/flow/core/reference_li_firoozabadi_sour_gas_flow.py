"""Independent low-temperature flow-property oracle for the Li-Firoozabadi sour-gas benchmark.

The three-phase equilibrium is imported from the existing independent Decimal
PR76 oracle.  This script does not import MPMC production C++ and does not fit
C++ output.

Sources/model boundary:
- Li & Firoozabadi, SPE Journal 17(4), 2012, DOI 10.2118/129844-PA:
  PR76 sour-gas EOS state/parameters.
- NIST Chemistry WebBook SRD 69: molar masses, low-temperature ideal-gas Cp,
  and critical molar volumes/densities.
- Lohrenz, Bray & Clark, JPT 1964, DOI 10.2118/915-PA:
  dense-fluid viscosity correlation already used by MPMC.
- Peaceman, SPEJ 1983, DOI 10.2118/10528-PA:
  Cartesian well index.
- NIST SRD 30, Anderson (1936), citation Z00788:
  99.93% quartz density and low-temperature heat capacity.
- Odeh/SPE1 engineering benchmark input as distributed by OPM opm-data:
  1000 x 1000 x 20 ft cell, phi=0.3, K=500 mD, rw=0.25 ft.

The stationary short step uses p_bhp=p_cell and pc=none, so the independent
Peaceman well-rate oracle is exactly zero.  No relative-permeability value is
part of the external numerical oracle because there is no face flux and zero
well drawdown.
"""

from __future__ import annotations

import argparse
import importlib.util
from decimal import Decimal as D, localcontext
from pathlib import Path

HERE = Path(__file__).resolve().parent
FLASH_REFERENCE = (
    HERE.parent.parent
    / "flash"
    / "pr76_three_phase"
    / "reference_sour_gas_decimal.py"
)

COMPONENTS = (
    "carbon-dioxide",
    "nitrogen",
    "hydrogen-sulfide",
    "methane",
    "ethane",
    "propane",
)
TC = tuple(map(D, ("304.211", "126.2", "373.2", "190.564", "305.322", "369.825")))
PC_PA = tuple(D(v) * D(100000) for v in ("73.819", "33.9", "89.4", "45.992", "48.718", "42.462"))
OMEGA = tuple(map(D, ("0.225", "0.039", "0.081", "0.01141", "0.10574", "0.15813")))
MW = tuple(map(D, ("0.0440095", "0.0280134", "0.034081", "0.0160425", "0.0300690", "0.0440956")))
VC_M3_PER_MOL = (
    D("0.0919e-3"),
    D(1) / D("11.18") * D("1e-3"),
    D(1) / D("10.2") * D("1e-3"),
    D("0.09860e-3"),
    D("0.147e-3"),
    D("0.200e-3"),
)

CP_T = tuple(map(D, ("100", "200", "298.15")))
CP = (
    tuple(map(D, ("29.208", "32.359", "37.129"))),
    tuple(map(D, ("29.10366858", "29.10731520", "29.12379083308665"))),
    tuple(map(D, ("33.259", "33.380", "34.192"))),
    tuple(map(D, ("33.28", "33.51", "35.69"))),
    tuple(map(D, ("35.70", "42.30", "52.49"))),
    tuple(map(D, ("41.30", "56.07", "73.60"))),
)

NONZERO_KIJ = {
    (0, 1): D("-0.020"),
    (0, 2): D("0.120"),
    (1, 2): D("0.200"),
    (0, 3): D("0.125"),
    (1, 3): D("0.031"),
    (2, 3): D("0.100"),
    (0, 4): D("0.135"),
    (1, 4): D("0.042"),
    (2, 4): D("0.080"),
    (0, 5): D("0.150"),
    (1, 5): D("0.091"),
    (2, 5): D("0.080"),
}

R = D("8.31446261815324")
PRESSURE_PA = D("2000000")
TEMPERATURE_K = D("178.8")
POROSITY = D("0.3")
DX_M = D("1000") * D("0.3048")
DY_M = D("1000") * D("0.3048")
DZ_M = D("20") * D("0.3048")
BULK_VOLUME_M3 = DX_M * DY_M * DZ_M
PERMEABILITY_M2 = D("500") * D("9.869233e-16")
WELLBORE_RADIUS_M = D("0.25") * D("0.3048")
TIMESTEP_SECONDS = D("1")

QUARTZ_T0_K = D("169.1")
QUARTZ_T1_K = D("184.8")
QUARTZ_CP0_J_PER_KG_K = D("469.3")
QUARTZ_CP1_J_PER_KG_K = D("510.6")
QUARTZ_DENSITY_KG_PER_M3 = D("2.6378") * D("1000")


def load_flash_reference():
    spec = importlib.util.spec_from_file_location(
        "mpmc_li_firoozabadi_decimal_reference",
        FLASH_REFERENCE,
    )
    if spec is None or spec.loader is None:
        raise RuntimeError("failed to load independent PR76 Decimal reference")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def kij_matrix():
    result = [[D(0) for _ in COMPONENTS] for _ in COMPONENTS]
    for (i, j), value in NONZERO_KIJ.items():
        result[i][j] = value
        result[j][i] = value
    return result


def integrated_cp(temperature, component):
    if temperature < CP_T[0] or temperature > CP_T[-1]:
        raise ValueError("low-temperature Cp oracle is valid only on [100,298.15] K")
    magnitude = D(0)
    values = CP[component]
    for segment in range(len(CP_T) - 1):
        t0 = CP_T[segment]
        t1 = CP_T[segment + 1]
        if temperature >= t1:
            continue
        slope = (values[segment + 1] - values[segment]) / (t1 - t0)
        if temperature <= t0:
            width = t1 - t0
            magnitude += values[segment] * width + D("0.5") * slope * width * width
        else:
            delta = t1 - temperature
            magnitude += values[segment + 1] * delta - D("0.5") * slope * delta * delta
    return -magnitude


def pr_parameters():
    kij = kij_matrix()
    pure_a = []
    pure_da = []
    pure_b = []
    for tc, pc, omega in zip(TC, PC_PA, OMEGA):
        kappa = D("0.37464") + D("1.54226") * omega - D("0.26992") * omega * omega
        ac = D("0.45724") * (R * tc) * (R * tc) / pc
        b = D("0.07780") * R * tc / pc
        sqrt_tr = (TEMPERATURE_K / tc).sqrt()
        factor = D(1) + kappa * (D(1) - sqrt_tr)
        a = ac * factor * factor
        da_dt = -ac * factor * kappa / (TEMPERATURE_K * tc).sqrt()
        pure_a.append(a)
        pure_da.append(da_dt)
        pure_b.append(b)
    return kij, pure_a, pure_da, pure_b


def phase_property(composition, zfactor, parameters):
    kij, pure_a, pure_da, pure_b = parameters
    n = len(COMPONENTS)
    attraction = D(0)
    d_attraction_dt = D(0)
    covolume = D(0)
    for i in range(n):
        covolume += composition[i] * pure_b[i]
        for j in range(n):
            aij = (pure_a[i] * pure_a[j]).sqrt() * (D(1) - kij[i][j])
            daij_dt = D("0.5") * aij * (
                pure_da[i] / pure_a[i] + pure_da[j] / pure_a[j]
            )
            weight = composition[i] * composition[j]
            attraction += weight * aij
            d_attraction_dt += weight * daij_dt

    b_reduced = covolume * PRESSURE_PA / (R * TEMPERATURE_K)
    sqrt2 = D(2).sqrt()
    logarithm = (
        (zfactor + (D(1) + sqrt2) * b_reduced)
        / (zfactor + (D(1) - sqrt2) * b_reduced)
    ).ln()
    departure = (
        R * TEMPERATURE_K * (zfactor - D(1))
        + (TEMPERATURE_K * d_attraction_dt - attraction)
        / (D(2) * sqrt2 * covolume)
        * logarithm
    )

    ideal = sum(
        composition[i] * integrated_cp(TEMPERATURE_K, i)
        for i in range(n)
    )
    mixture_mw = sum(composition[i] * MW[i] for i in range(n))
    molar_density = PRESSURE_PA / (zfactor * R * TEMPERATURE_K)
    mass_density = molar_density * mixture_mw

    pa_per_psia = D("6894.757293168")
    dilute_numerator = D(0)
    dilute_denominator = D(0)
    pseudo_tc_rankine = D(0)
    pseudo_pc_psia = D(0)
    pseudo_vc = D(0)
    for i in range(n):
        tc_rankine = TC[i] * D("1.8")
        pc_psia = PC_PA[i] / pa_per_psia
        mw_g_per_mol = MW[i] * D(1000)
        xi = (
            D("5.4402")
            * tc_rankine ** (D(1) / D(6))
            / (mw_g_per_mol.sqrt() * pc_psia ** (D(2) / D(3)))
        )
        tr = TEMPERATURE_K / TC[i]
        if tr <= D("1.5"):
            pure_viscosity_cp = D("34.0e-5") * tr ** D("0.94") / xi
        else:
            pure_viscosity_cp = (
                D("17.78e-5")
                * (D("4.58") * tr - D("1.67")) ** D("0.625")
                / xi
            )
        sqrt_mw = mw_g_per_mol.sqrt()
        dilute_numerator += composition[i] * pure_viscosity_cp * sqrt_mw
        dilute_denominator += composition[i] * sqrt_mw
        pseudo_tc_rankine += composition[i] * tc_rankine
        pseudo_pc_psia += composition[i] * pc_psia
        pseudo_vc += composition[i] * VC_M3_PER_MOL[i]

    dilute_cp = dilute_numerator / dilute_denominator
    xi_m = (
        D("5.4402")
        * pseudo_tc_rankine ** (D(1) / D(6))
        / ((mixture_mw * D(1000)).sqrt() * pseudo_pc_psia ** (D(2) / D(3)))
    )
    reduced_density = molar_density * pseudo_vc
    polynomial = (
        D("0.1023")
        + D("0.023364") * reduced_density
        + D("0.058533") * reduced_density ** 2
        - D("0.040758") * reduced_density ** 3
        + D("0.0093324") * reduced_density ** 4
    )
    viscosity_pa_s = (dilute_cp + (polynomial ** 4 - D("1e-4")) / xi_m) * D("1e-3")

    specific_enthalpy = (ideal + departure) / mixture_mw
    specific_internal_energy = specific_enthalpy - PRESSURE_PA / mass_density

    return {
        "mixture_mw": mixture_mw,
        "molar_density": molar_density,
        "mass_density": mass_density,
        "viscosity": viscosity_pa_s,
        "ideal_molar_enthalpy": ideal,
        "residual_molar_enthalpy": departure,
        "specific_enthalpy": specific_enthalpy,
        "specific_internal_energy": specific_internal_energy,
    }


def quartz_heat_capacity(temperature):
    if temperature < QUARTZ_T0_K or temperature > QUARTZ_T1_K:
        raise ValueError("quartz benchmark interpolation is valid only on [169.1,184.8] K")
    slope = (
        QUARTZ_CP1_J_PER_KG_K - QUARTZ_CP0_J_PER_KG_K
    ) / (QUARTZ_T1_K - QUARTZ_T0_K)
    return QUARTZ_CP0_J_PER_KG_K + slope * (temperature - QUARTZ_T0_K)


def solve_reference(precision):
    with localcontext() as context:
        context.prec = precision
        flash = load_flash_reference()
        _, phases, fractions, zfactors, _ = flash.solve_reference(precision)
        parameters = pr_parameters()
        properties = [
            phase_property(phase, zfactor, parameters)
            for phase, zfactor in zip(phases, zfactors)
        ]

        phase_pore_volumes = [
            fractions[i] / properties[i]["molar_density"]
            for i in range(3)
        ]
        pore_volume_sum = sum(phase_pore_volumes)
        saturations = [
            value / pore_volume_sum
            for value in phase_pore_volumes
        ]

        component_accumulation = [
            POROSITY
            * sum(
                saturations[phase]
                * properties[phase]["molar_density"]
                * phases[phase][component]
                for phase in range(3)
            )
            for component in range(len(COMPONENTS))
        ]

        fluid_energy = POROSITY * sum(
            saturations[phase]
            * properties[phase]["mass_density"]
            * properties[phase]["specific_internal_energy"]
            for phase in range(3)
        )

        quartz_cp = quartz_heat_capacity(TEMPERATURE_K)
        quartz_volumetric_heat_capacity = QUARTZ_DENSITY_KG_PER_M3 * quartz_cp

        # Rock internal energy is referenced to zero at the benchmark T=178.8 K.
        # Therefore the reference-state rock energy is exactly zero while the
        # non-isothermal Jacobian retains the sourced rho*Cp temperature slope.
        total_energy = fluid_energy

        return {
            "phases": phases,
            "fractions": fractions,
            "zfactors": zfactors,
            "properties": properties,
            "saturations": saturations,
            "component_accumulation": component_accumulation,
            "energy_accumulation": total_energy,
            "quartz_cp": quartz_cp,
            "quartz_volumetric_heat_capacity": quartz_volumetric_heat_capacity,
        }


def literal(value):
    return format(value, ".17g")


def array(name, values):
    return (
        f"inline constexpr std::array<double, {len(values)}> {name}{{\n    "
        + ",\n    ".join(literal(value) for value in values)
        + "};\n"
    )


def render(result):
    props = result["properties"]
    text = """// Generated only by reference_li_firoozabadi_sour_gas_flow.py.
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
"""
    text += array("phase_saturation", result["saturations"])
    text += array("molar_density_mol_per_m3", [p["molar_density"] for p in props])
    text += array("mass_density_kg_per_m3", [p["mass_density"] for p in props])
    text += array("dynamic_viscosity_pa_s", [p["viscosity"] for p in props])
    text += array("specific_enthalpy_j_per_kg", [p["specific_enthalpy"] for p in props])
    text += array("specific_internal_energy_j_per_kg", [p["specific_internal_energy"] for p in props])
    text += array("mixture_molar_mass_kg_per_mol", [p["mixture_mw"] for p in props])
    text += array("ideal_gas_molar_enthalpy_j_per_mol", [p["ideal_molar_enthalpy"] for p in props])
    text += array("residual_molar_enthalpy_j_per_mol", [p["residual_molar_enthalpy"] for p in props])
    text += array("component_accumulation_mol_per_bulk_m3", result["component_accumulation"])
    text += f"inline constexpr double total_internal_energy_j_per_bulk_m3 = {literal(result['energy_accumulation'])};\n"
    text += f"inline constexpr double quartz_specific_heat_j_per_kg_k = {literal(result['quartz_cp'])};\n"
    text += f"inline constexpr double quartz_volumetric_heat_capacity_j_per_m3_k = {literal(result['quartz_volumetric_heat_capacity'])};\n"
    text += array("cell_component_inventory_mol", [
        value * BULK_VOLUME_M3
        for value in result["component_accumulation"]
    ])
    text += f"inline constexpr double cell_total_internal_energy_j = {literal(result['energy_accumulation'] * BULK_VOLUME_M3)};\n"
    text += """inline constexpr std::array<double, 3> well_phase_volumetric_rate_m3_per_s{
    0.0, 0.0, 0.0};
inline constexpr std::array<double, 6> well_component_molar_rate_mol_per_s{
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
inline constexpr double well_energy_rate_w = 0.0;
} // namespace pr76_li_firoozabadi_sour_gas_flow_reference

#endif
"""
    return text


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--precision", type=int, choices=(80, 96), default=80)
    parser.add_argument("--write", action="store_true")
    args = parser.parse_args()

    result = solve_reference(args.precision)
    output = render(result)
    path = HERE / "pr76_li_firoozabadi_sour_gas_flow_reference.hpp"

    if args.write:
        path.write_text(output, encoding="utf-8")
    elif output != path.read_text(encoding="utf-8"):
        raise AssertionError(
            "Li-Firoozabadi sour-gas flow reference mismatch; do not fit production output"
        )

    print(
        f"Decimal({args.precision}) Li-Firoozabadi sour-gas flow reference checked"
    )
    print(
        "saturations:",
        *(literal(value) for value in result["saturations"]),
    )
    print(
        "viscosity [Pa s]:",
        *(literal(p["viscosity"]) for p in result["properties"]),
    )


if __name__ == "__main__":
    main()
