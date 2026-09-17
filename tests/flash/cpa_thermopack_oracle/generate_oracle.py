#!/usr/bin/env python3
"""Generate an independent ThermoPack CPA TP-flash oracle.

This script intentionally imports only ThermoPack and Python standard-library
modules.  It does not import MPMC_HNU production or test code, so the generated
values cannot be circularly derived from the implementation under test.

The caller must provide ThermoPack built from the pinned revision recorded in
PR #115 / cpa_thermopack_mapping.md on PYTHONPATH / the active environment.
"""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

from thermopack.cpa import SRK_CPA


THERMOPACK_COMMIT = "d68c794c7342bfc6938eb424a1fbb88b7780b738"
TEMPERATURE_K = 333.15
Kij_METHANOL_WATER = -0.055

# ThermoPack set_pure_params uses conventional CPA units:
#   a0  [Pa L^2 mol^-2]
#   b   [L mol^-1]
#   eps [J mol^-1]
#   beta, c1 dimensionless.
# These are exactly the SI-converted MPMC_HNU fixture numbers expressed in the
# units required by the external program's public API.
PURE = {
    1: {
        "component": "MEOH",
        "params": [4.0531e5, 0.030978, 24591.0, 0.0161, 0.43102],
        "association_scheme": "2B",
    },
    2: {
        "component": "H2O",
        "params": [1.2277e5, 0.014515, 16655.0, 0.0692, 0.67359],
        "association_scheme": "4C",
    },
}

# Existing traceable 333.15 K Kurihara et al. points used by the repository's
# physical CPA regression. Experimental x/y are used ONLY to construct the same
# material-balanced overall feed as the MPMC_HNU validation fixture. They are not
# supplied to ThermoPack as equilibrium phase guesses or expected solutions.
EXPERIMENTAL_CARRIERS = [
    {"pressure_pa": 3.9223e4, "x_meoh_exp": 0.1686, "y_meoh_exp": 0.5714},
    {"pressure_pa": 4.8852e4, "x_meoh_exp": 0.3039, "y_meoh_exp": 0.6943},
    {"pressure_pa": 5.6652e4, "x_meoh_exp": 0.4461, "y_meoh_exp": 0.7742},
    {"pressure_pa": 6.3998e4, "x_meoh_exp": 0.6044, "y_meoh_exp": 0.8383},
    {"pressure_pa": 7.2832e4, "x_meoh_exp": 0.7776, "y_meoh_exp": 0.9141},
]


def _require_close(actual: float, expected: float, message: str) -> None:
    scale = max(1.0, abs(actual), abs(expected))
    if not math.isfinite(actual) or abs(actual - expected) > 5.0e-13 * scale:
        raise RuntimeError(f"{message}: expected {expected!r}, got {actual!r}")


def configure_matched_model() -> SRK_CPA:
    # Start from the pinned ThermoPack component records because they also own
    # the 2B/4C association topology and CR-1 combining-rule structure.
    eos = SRK_CPA("MEOH,H2O", mixing="vdW", alpha="Classic",
                  parameter_reference="Default")

    # Freeze the same formulation as MPMC_HNU: simplified RDF and STANDARD
    # association mixing (CR-1), not Elliott's Delta combining rule.
    eos.set_cpa_formulation(True, False)

    # Do not rely on database defaults even when the current values match.
    # Explicitly install and then read back the exact pure parameter vector.
    for index, record in PURE.items():
        eos.set_pure_params(index, record["params"])
        observed = eos.get_pure_params(index)
        if len(observed) != 5:
            raise RuntimeError("ThermoPack returned an unexpected CPA pure-parameter shape")
        for actual, expected in zip(observed, record["params"]):
            _require_close(float(actual), float(expected),
                           f"pure parameter mismatch for {record['component']}")

    # ThermoPack's DEFAULT MEOH/H2O binary record uses kij_a=-0.09.  The MPMC
    # literature fixture at 333.15 K uses Folas CR-1 kij_a=-0.055.  Override the
    # cubic interaction explicitly; kij_eps=0 preserves the CR-1 arithmetic
    # epsilon combination. The database beta combining rule remains GEOMETRIC.
    eos.set_kij(1, 2, Kij_METHANOL_WATER, 0.0)
    observed_kij = eos.get_kij(1, 2)
    _require_close(float(observed_kij[0]), Kij_METHANOL_WATER,
                   "ThermoPack cubic kij override failed")
    _require_close(float(observed_kij[1]), 0.0,
                   "ThermoPack association-energy kij override failed")
    return eos


def generate() -> dict:
    eos = configure_matched_model()
    states = []

    for point in EXPERIMENTAL_CARRIERS:
        z_meoh = 0.9 * point["x_meoh_exp"] + 0.1 * point["y_meoh_exp"]
        z = [z_meoh, 1.0 - z_meoh]
        flash = eos.two_phase_tpflash(TEMPERATURE_K, point["pressure_pa"], z)

        x = [float(v) for v in flash.x]
        y = [float(v) for v in flash.y]
        beta_v = float(flash.betaV)
        beta_l = float(flash.betaL)
        phase = int(flash.phase)

        if phase != int(eos.TWOPH):
            raise RuntimeError(
                f"expected a ThermoPack two-phase result at {point['pressure_pa']} Pa, "
                f"got phase code {phase}")
        if len(x) != 2 or len(y) != 2:
            raise RuntimeError("ThermoPack returned an unexpected binary phase-composition shape")
        if not (0.0 < beta_v < 1.0 and 0.0 < beta_l < 1.0):
            raise RuntimeError("ThermoPack returned a non-interior phase fraction")
        _require_close(beta_v + beta_l, 1.0, "ThermoPack phase fractions do not normalize")

        reconstructed = [beta_l * x[i] + beta_v * y[i] for i in range(2)]
        mb_max_abs = max(abs(reconstructed[i] - z[i]) for i in range(2))
        if mb_max_abs > 2.0e-10:
            raise RuntimeError(
                f"ThermoPack oracle extraction failed material balance: {mb_max_abs}")

        states.append({
            "temperature_k": TEMPERATURE_K,
            "pressure_pa": point["pressure_pa"],
            "feed": {"MEOH": z[0], "H2O": z[1]},
            "phase_code": phase,
            "beta_vapor": beta_v,
            "beta_liquid": beta_l,
            "liquid": {"MEOH": x[0], "H2O": x[1]},
            "vapor": {"MEOH": y[0], "H2O": y[1]},
            "extraction_material_balance_max_abs": mb_max_abs,
        })

    pure_readback = {
        PURE[i]["component"]: [float(v) for v in eos.get_pure_params(i)]
        for i in sorted(PURE)
    }
    kij_readback = [float(v) for v in eos.get_kij(1, 2)]

    return {
        "schema": "MPMC_HNU/CPA/ThermoPack-two-phase-TP-oracle/v1",
        "generator": {
            "software": "thermotools/thermopack",
            "commit": THERMOPACK_COMMIT,
            "api": "SRK_CPA + two_phase_tpflash",
        },
        "model": {
            "components": ["MEOH", "H2O"],
            "physical_eos": "SRK",
            "mixing": "vdW",
            "alpha": "Classic",
            "simplified_cpa": True,
            "association_delta_rule": "STANDARD/CR-1",
            "association_schemes": {"MEOH": "2B", "H2O": "4C"},
            "pure_parameter_units": [
                "Pa L^2 mol^-2", "L mol^-1", "J mol^-1", "1", "1"
            ],
            "pure_parameter_readback": pure_readback,
            "kij_readback": {"kij_a": kij_readback[0], "kij_eps": kij_readback[1]},
            "cross_beta_rule": "GEOMETRIC from pinned ThermoPack CPA binary record",
        },
        "states": states,
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--out", type=Path, default=None)
    args = parser.parse_args()

    result = generate()
    rendered = json.dumps(result, indent=2, sort_keys=True)
    if args.out is not None:
        args.out.parent.mkdir(parents=True, exist_ok=True)
        args.out.write_text(rendered + "\n", encoding="utf-8")
    print("THERMOPACK_CPA_ORACLE_BEGIN")
    print(rendered)
    print("THERMOPACK_CPA_ORACLE_END")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
