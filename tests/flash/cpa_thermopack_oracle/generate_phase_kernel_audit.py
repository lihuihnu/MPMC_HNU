#!/usr/bin/env python3
"""Generate an independent ThermoPack CPA phase-kernel diagnostic oracle.

The frozen TP-flash oracle owns the five equilibrium phase compositions. This
script evaluates ThermoPack single-phase properties at those same T/P/x states.
The cubic contribution is reconstructed from the pinned ThermoPack component
source Tc plus its public CPA a0/b/c1/kij read-back using the documented classic
alpha SRK equations. Association is then full ThermoPack residual minus that
cubic contribution at identical T,V,n.

Why source Tc instead of get_critical_parameters(): the latter solves for the
pure-fluid critical point of the active EOS. It is not the component Tc stored in
cbeos%single(i)%Tc and used by cbCalcAlphaTerm. Pinned SelectCubicEOS initializes
that alpha Tc from comp(i)%p_comp%tc before CPA overwrites a0/b/c1.
"""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
import sys

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from generate_oracle import (  # noqa: E402
    TEMPERATURE_K,
    THERMOPACK_COMMIT,
    configure_matched_model,
)

BASE_ORACLE = Path(__file__).resolve().parent / "thermopack_d68c794_meoh_h2o_33315k.json"
COMPONENTS = ("MEOH", "H2O")


def _value(obj):
    return obj.val if hasattr(obj, "val") else obj


def _scalar(obj) -> float:
    array = np.asarray(_value(obj))
    if array.size != 1:
        raise RuntimeError(f"expected scalar ThermoPack property, got shape {array.shape}")
    result = float(array.reshape(-1)[0])
    if not math.isfinite(result):
        raise RuntimeError("ThermoPack scalar property became nonfinite")
    return result


def _vector(obj, size: int) -> list[float]:
    array = np.asarray(_value(obj), dtype=float).reshape(-1)
    if array.size != size:
        raise RuntimeError(
            f"expected ThermoPack vector of length {size}, got shape {array.shape}")
    result = [float(v) for v in array]
    if not all(math.isfinite(v) for v in result):
        raise RuntimeError("ThermoPack vector property became nonfinite")
    return result


def _require_tp_root_consistency(actual: float, expected: float, message: str) -> None:
    scale = max(1.0, abs(actual), abs(expected))
    if not math.isfinite(actual) or abs(actual - expected) > 3.0e-8 * scale:
        raise RuntimeError(f"{message}: expected {expected!r}, got {actual!r}")


def _source_alpha_tc(thermopack_source: Path) -> tuple[list[float], dict]:
    files = {
        "MEOH": thermopack_source / "fluids" / "Methanol.json",
        "H2O": thermopack_source / "fluids" / "Water.json",
    }
    values = []
    provenance = {}
    for component in COMPONENTS:
        path = files[component]
        record = json.loads(path.read_text(encoding="utf-8"))
        if record["ident"] != component:
            raise RuntimeError(f"unexpected ThermoPack component file identity: {path}")
        critical = record["critical"]
        tc = float(critical["temperature"])
        if not math.isfinite(tc) or not tc > 0.0:
            raise RuntimeError(f"invalid source Tc for {component}")
        values.append(tc)
        provenance[component] = {
            "temperature_k": tc,
            "ref": str(critical.get("ref", "")),
            "source_file": str(path.relative_to(thermopack_source)),
        }
    return values, provenance


def _cubic_srkcpa_at_tv(eos, temperature_k: float, density: float,
                        composition: list[float], critical_t: list[float]) -> dict:
    """Documented SRK-CPA cubic pressure and residual chemical potential / RT."""
    pure = [[float(v) for v in eos.get_pure_params(i)] for i in (1, 2)]
    kij = float(eos.get_kij(1, 2)[0])
    rgas = float(eos.Rgas)
    rt = rgas * temperature_k

    a0 = [record[0] * 1.0e-6 for record in pure]  # Pa L^2 -> Pa m^6
    b = [record[1] * 1.0e-3 for record in pure]   # L/mol -> m^3/mol
    c1 = [record[4] for record in pure]
    ai = [
        a0[i] * (1.0 + c1[i] *
                 (1.0 - math.sqrt(temperature_k / critical_t[i]))) ** 2
        for i in range(2)
    ]
    aij = [[0.0, 0.0], [0.0, 0.0]]
    for i in range(2):
        for j in range(2):
            kij_ij = kij if i != j else 0.0
            aij[i][j] = math.sqrt(ai[i] * ai[j]) * (1.0 - kij_ij)

    a_mix = sum(
        composition[i] * composition[j] * aij[i][j]
        for i in range(2) for j in range(2)
    )
    b_mix = sum(composition[i] * b[i] for i in range(2))
    b_rho = b_mix * density
    if not 0.0 < b_rho < 1.0:
        raise RuntimeError("ThermoPack audit cubic state crossed SRK covolume singularity")

    p_cubic = (
        rt * density / (1.0 - b_rho)
        - a_mix * density * density / (1.0 + b_rho)
    )
    z_cubic = p_cubic / (density * rt)
    sums = [
        sum(composition[j] * aij[i][j] for j in range(2))
        for i in range(2)
    ]
    a_over_brt = a_mix / (b_mix * rt)
    log_free = math.log1p(-b_rho)
    log_attr = math.log1p(b_rho)
    mu_cubic_over_rt = []
    for i in range(2):
        b_ratio = b[i] / b_mix
        attraction_ratio = 2.0 * sums[i] / a_mix - b_ratio
        mu_cubic_over_rt.append(
            b_ratio * (z_cubic - 1.0)
            - log_free
            - a_over_brt * attraction_ratio * log_attr
        )

    return {
        "a_i_pa_m6_per_mol2": {COMPONENTS[i]: ai[i] for i in range(2)},
        "a_mix_pa_m6_per_mol2": a_mix,
        "b_mix_m3_per_mol": b_mix,
        "pressure_pa": p_cubic,
        "mu_residual_over_rt": {
            COMPONENTS[i]: mu_cubic_over_rt[i] for i in range(2)
        },
    }


def phase_properties(full, pressure_pa: float, composition: list[float],
                     phase_flag: int, phase_name: str,
                     critical_t: list[float]) -> dict:
    volume = _scalar(full.specific_volume(
        TEMPERATURE_K, pressure_pa, composition, phase_flag))
    if not volume > 0.0:
        raise RuntimeError("ThermoPack returned non-positive phase volume")
    density = 1.0 / volume
    z = _scalar(full.zfac(TEMPERATURE_K, pressure_pa, composition, phase_flag))
    ln_phi = _vector(full.thermo(
        TEMPERATURE_K, pressure_pa, composition, phase_flag), 2)
    pressure_total_tv = _scalar(full.pressure_tv(
        TEMPERATURE_K, volume, composition, property_flag="IR"))
    mu_res_total = _vector(full.chemical_potential_tv(
        TEMPERATURE_K, volume, composition, property_flag="R"), 2)

    rgas = float(full.Rgas)
    rt = rgas * TEMPERATURE_K
    cubic = _cubic_srkcpa_at_tv(
        full, TEMPERATURE_K, density, composition, critical_t)
    association_pressure = pressure_total_tv - cubic["pressure_pa"]
    total_mu_over_rt = [value / rt for value in mu_res_total]
    association_mu_over_rt = [
        total_mu_over_rt[i] - cubic["mu_residual_over_rt"][COMPONENTS[i]]
        for i in range(2)
    ]

    _require_tp_root_consistency(
        pressure_total_tv, pressure_pa,
        f"ThermoPack {phase_name} TP/TV pressure mismatch")
    _require_tp_root_consistency(
        z, pressure_pa * volume / rt,
        f"ThermoPack {phase_name} Z/volume inconsistency")

    return {
        "phase": phase_name,
        "composition": {"MEOH": composition[0], "H2O": composition[1]},
        "molar_volume_m3_per_mol": volume,
        "molar_density_mol_per_m3": density,
        "compressibility_factor": z,
        "ln_phi": {"MEOH": ln_phi[0], "H2O": ln_phi[1]},
        "pressure_total_tv_pa": pressure_total_tv,
        "pressure_physical_pa": cubic["pressure_pa"],
        "pressure_association_pa": association_pressure,
        "mu_residual_total_over_rt": {
            COMPONENTS[i]: total_mu_over_rt[i] for i in range(2)
        },
        "mu_cubic_over_rt": cubic["mu_residual_over_rt"],
        "mu_association_over_rt": {
            COMPONENTS[i]: association_mu_over_rt[i] for i in range(2)
        },
        "a_i_pa_m6_per_mol2": cubic["a_i_pa_m6_per_mol2"],
        "a_mix_pa_m6_per_mol2": cubic["a_mix_pa_m6_per_mol2"],
        "b_mix_m3_per_mol": cubic["b_mix_m3_per_mol"],
    }


def generate(thermopack_source: Path) -> dict:
    frozen = json.loads(BASE_ORACLE.read_text(encoding="utf-8"))
    if frozen["schema"] != "MPMC_HNU/CPA/ThermoPack-two-phase-TP-oracle/v1":
        raise RuntimeError("unexpected base ThermoPack oracle schema")
    if frozen["generator"]["commit"] != THERMOPACK_COMMIT:
        raise RuntimeError("base ThermoPack oracle revision drifted")

    full = configure_matched_model()
    critical_t, critical_provenance = _source_alpha_tc(thermopack_source)

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
            "liquid": phase_properties(
                full, pressure_pa, liquid, int(full.LIQPH), "liquid", critical_t),
            "vapor": phase_properties(
                full, pressure_pa, vapor, int(full.VAPPH), "vapor", critical_t),
        })

    return {
        "schema": "MPMC_HNU/CPA/ThermoPack-phase-kernel-audit/v2",
        "generator": {
            "software": "thermotools/thermopack",
            "commit": THERMOPACK_COMMIT,
            "phase_property_api": [
                "get_pure_params", "get_kij", "specific_volume", "zfac",
                "thermo", "pressure_tv", "chemical_potential_tv"
            ],
            "alpha_tc_source": "pinned component JSON used by SelectCubicEOS/initCubicTcPcAcf",
            "decomposition": (
                "full ThermoPack residual minus documented SRK cubic evaluated "
                "from pinned component Tc and ThermoPack read-back a0/b/c1/kij "
                "at identical T,V,n"
            ),
            "tp_root_consistency_relative_tolerance": 3.0e-8,
            "tp_root_consistency_basis": (
                "pinned saft_volume_solver pressure convergence: 1e8*machine_prec*P"
            ),
        },
        "model": frozen["model"],
        "gas_constant_j_per_mol_k": float(full.Rgas),
        "alpha_critical_temperature_source": critical_provenance,
        "states": states,
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--thermopack-source", type=Path, required=True)
    parser.add_argument("--out", type=Path, default=None)
    args = parser.parse_args()
    source = args.thermopack_source.resolve()
    if not (source / ".git").exists():
        raise RuntimeError("--thermopack-source must be the pinned ThermoPack checkout")
    result = generate(source)
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
