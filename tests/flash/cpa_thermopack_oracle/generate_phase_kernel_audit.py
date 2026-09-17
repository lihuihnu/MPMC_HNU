#!/usr/bin/env python3
"""Generate an independent ThermoPack CPA phase-kernel diagnostic oracle.

The frozen TP-flash oracle owns the five equilibrium phase compositions.  This
script evaluates ThermoPack single-phase properties at those same T/P/x states
and constructs an association-off *shadow* CPA model with identical SRK
parameters.  At the full-model phase volume, full minus shadow residual
properties isolate the association contribution without importing MPMC_HNU code
or copying ThermoPack's association implementation.
"""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
import sys

import numpy as np
from thermopack.cpa import SRK_CPA

# Import only the sibling external-oracle configuration.  It imports ThermoPack
# and standard-library modules, never MPMC_HNU production/test C++ code.
sys.path.insert(0, str(Path(__file__).resolve().parent))
from generate_oracle import (  # noqa: E402
    Kij_METHANOL_WATER,
    PURE,
    TEMPERATURE_K,
    THERMOPACK_COMMIT,
    _require_close,
    _set_pinned_cpa_formulation,
    configure_matched_model,
)

BASE_ORACLE = Path(__file__).resolve().parent / "thermopack_d68c794_meoh_h2o_33315k.json"


def _value(obj):
    """Unwrap ThermoPack v3 Property while remaining compatible with scalars."""
    return obj.val if hasattr(obj, "val") else obj


def _scalar(obj) -> float:
    value = _value(obj)
    array = np.asarray(value)
    if array.size != 1:
        raise RuntimeError(f"expected scalar ThermoPack property, got shape {array.shape}")
    result = float(array.reshape(-1)[0])
    if not math.isfinite(result):
        raise RuntimeError("ThermoPack scalar property became nonfinite")
    return result


def _vector(obj, size: int) -> list[float]:
    value = _value(obj)
    array = np.asarray(value, dtype=float).reshape(-1)
    if array.size != size:
        raise RuntimeError(
            f"expected ThermoPack vector of length {size}, got shape {array.shape}")
    result = [float(v) for v in array]
    if not all(math.isfinite(v) for v in result):
        raise RuntimeError("ThermoPack vector property became nonfinite")
    return result


def configure_nonassociating_shadow() -> SRK_CPA:
    """Build an SRK-identical CPA model with association Delta forced to zero."""
    eos = SRK_CPA("MEOH,H2O", mixing="vdW", alpha="Classic",
                  parameter_reference="Default")
    _set_pinned_cpa_formulation(eos)

    for index, record in PURE.items():
        params = list(record["params"])
        # beta multiplies Delta directly. beta=0 leaves a0/b/c1/epsilon and the
        # cubic mixing model untouched while making every association Delta zero.
        params[3] = 0.0
        eos.set_pure_params(index, params)
        observed = [float(v) for v in eos.get_pure_params(index)]
        for position, (actual, expected) in enumerate(zip(observed, params)):
            _require_close(
                actual, float(expected),
                f"nonassociating shadow pure parameter mismatch component={index} field={position}")

    eos.set_kij(1, 2, Kij_METHANOL_WATER, 0.0)
    observed_kij = [float(v) for v in eos.get_kij(1, 2)]
    _require_close(observed_kij[0], Kij_METHANOL_WATER,
                   "shadow cubic kij mismatch")
    _require_close(observed_kij[1], 0.0,
                   "shadow association-energy kij mismatch")
    return eos


def phase_properties(full: SRK_CPA, shadow: SRK_CPA,
                     pressure_pa: float, composition: list[float],
                     phase_flag: int, phase_name: str) -> dict:
    volume = _scalar(full.specific_volume(TEMPERATURE_K, pressure_pa,
                                          composition, phase_flag))
    if not volume > 0.0:
        raise RuntimeError("ThermoPack returned non-positive phase volume")
    density = 1.0 / volume
    z = _scalar(full.zfac(TEMPERATURE_K, pressure_pa, composition, phase_flag))
    ln_phi = _vector(full.thermo(TEMPERATURE_K, pressure_pa,
                                 composition, phase_flag), 2)

    # Re-evaluate at exactly the returned full-model TV state so the
    # association decomposition is independent of root-finding differences.
    pressure_total_tv = _scalar(full.pressure_tv(
        TEMPERATURE_K, volume, composition, property_flag="IR"))
    pressure_shadow_tv = _scalar(shadow.pressure_tv(
        TEMPERATURE_K, volume, composition, property_flag="IR"))

    mu_res_total = _vector(full.chemical_potential_tv(
        TEMPERATURE_K, volume, composition, property_flag="R"), 2)
    mu_res_shadow = _vector(shadow.chemical_potential_tv(
        TEMPERATURE_K, volume, composition, property_flag="R"), 2)

    rgas = float(full.Rgas)
    if not math.isfinite(rgas) or not rgas > 0.0:
        raise RuntimeError("ThermoPack returned invalid gas constant")
    rt = rgas * TEMPERATURE_K
    association_mu_over_rt = [
        (mu_res_total[i] - mu_res_shadow[i]) / rt for i in range(2)
    ]
    association_pressure = pressure_total_tv - pressure_shadow_tv

    _require_close(pressure_total_tv, pressure_pa,
                   f"ThermoPack {phase_name} TP/TV pressure mismatch")
    z_from_volume = pressure_pa * volume / rt
    _require_close(z, z_from_volume,
                   f"ThermoPack {phase_name} Z/volume inconsistency")

    return {
        "phase": phase_name,
        "composition": {"MEOH": composition[0], "H2O": composition[1]},
        "molar_volume_m3_per_mol": volume,
        "molar_density_mol_per_m3": density,
        "compressibility_factor": z,
        "ln_phi": {"MEOH": ln_phi[0], "H2O": ln_phi[1]},
        "pressure_total_tv_pa": pressure_total_tv,
        "pressure_cubic_shadow_tv_pa": pressure_shadow_tv,
        "pressure_association_effect_pa": association_pressure,
        "mu_residual_total_j_per_mol": {
            "MEOH": mu_res_total[0], "H2O": mu_res_total[1]
        },
        "mu_residual_cubic_shadow_j_per_mol": {
            "MEOH": mu_res_shadow[0], "H2O": mu_res_shadow[1]
        },
        "mu_association_over_rt": {
            "MEOH": association_mu_over_rt[0],
            "H2O": association_mu_over_rt[1],
        },
    }


def generate() -> dict:
    frozen = json.loads(BASE_ORACLE.read_text(encoding="utf-8"))
    if frozen["schema"] != "MPMC_HNU/CPA/ThermoPack-two-phase-TP-oracle/v1":
        raise RuntimeError("unexpected base ThermoPack oracle schema")
    if frozen["generator"]["commit"] != THERMOPACK_COMMIT:
        raise RuntimeError("base ThermoPack oracle revision drifted")

    full = configure_matched_model()
    shadow = configure_nonassociating_shadow()

    states = []
    for source in frozen["states"]:
        pressure_pa = float(source["pressure_pa"])
        liquid = [float(source["liquid"]["MEOH"]),
                  float(source["liquid"]["H2O"])]
        vapor = [float(source["vapor"]["MEOH"]),
                 float(source["vapor"]["H2O"])]
        states.append({
            "temperature_k": float(source["temperature_k"]),
            "pressure_pa": pressure_pa,
            "liquid": phase_properties(full, shadow, pressure_pa, liquid,
                                        int(full.LIQPH), "liquid"),
            "vapor": phase_properties(full, shadow, pressure_pa, vapor,
                                       int(full.VAPPH), "vapor"),
        })

    return {
        "schema": "MPMC_HNU/CPA/ThermoPack-phase-kernel-audit/v1",
        "generator": {
            "software": "thermotools/thermopack",
            "commit": THERMOPACK_COMMIT,
            "phase_property_api": [
                "specific_volume", "zfac", "thermo", "pressure_tv",
                "chemical_potential_tv"
            ],
            "decomposition": (
                "full SRK-CPA minus association-Delta-zero shadow at identical T,V,n"
            ),
        },
        "model": frozen["model"],
        "gas_constant_j_per_mol_k": float(full.Rgas),
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
    print("THERMOPACK_CPA_PHASE_KERNEL_AUDIT_BEGIN")
    print(rendered)
    print("THERMOPACK_CPA_PHASE_KERNEL_AUDIT_END")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
